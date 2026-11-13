"""Unit tests against hand-calculated reference values (pytest).
Each expected value is computed by hand from the published equation, independent of the code."""
import os, sys, math, json, subprocess
import numpy as np, pytest
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import egfr
from doac_rules import RULES, RULES_US

# (age, weight, height, female, scr µmol/L) -> hand values
# CG man 60 y, 80 kg, Scr 88.4 (1.0 mg/dL): (140-60)*80/72 = 88.889
# CG woman 80 y, 60 kg, Scr 88.4: 60*60/72*0.85 = 42.5
# CKD-EPI 2021 man 50 y, Scr 79.56 (0.9 mg/dL, ratio 1): 142*0.9938^50
# CKD-EPI 2009 woman 70 y, Scr 61.88 (0.7 mg/dL, ratio 1): 141*0.993^70*1.018
# EKFC man 60 y, Scr 80 (=Q): 107.3*0.990^20 ; EKFC woman 30 y, Scr 124 (2×Q): 107.3*2^-1.132
# DuBois 70 kg 175 cm: 0.007184*70^0.425*175^0.725 ; Mosteller: sqrt(70*175/3600)
CASES = [
    ("cg", lambda: egfr.cockcroft_gault(60, 80, 88.4, False), (140 - 60) * 80 / 72),
    ("cg_f", lambda: egfr.cockcroft_gault(80, 60, 88.4, True), 60 * 60 / 72 * 0.85),
    ("ckdepi21", lambda: egfr.ckd_epi_2021(50, 0.9 * 88.4, False), 142 * 0.9938 ** 50),
    ("ckdepi09_f", lambda: egfr.ckd_epi_2009(70, 0.7 * 88.4, True), 141 * 0.993 ** 70 * 1.018),
    ("ekfc", lambda: egfr.ekfc(60, 80, False), 107.3 * 0.990 ** 20),
    ("ekfc_f_high", lambda: egfr.ekfc(30, 124, True), 107.3 * 2 ** -1.132),
    ("dubois", lambda: egfr.bsa_dubois(70, 175), 0.007184 * 70 ** 0.425 * 175 ** 0.725),
    ("mosteller", lambda: egfr.bsa_mosteller(70, 175), math.sqrt(70 * 175 / 3600)),
]

@pytest.mark.parametrize("name,fn,expected", CASES)
def test_python_reference(name, fn, expected):
    assert float(np.asarray(fn())) == pytest.approx(expected, rel=1e-9)

def _js(expr):
    core = os.path.join(HERE, "doac_core.js")
    out = subprocess.run(["node", "-e", f"const D=require({json.dumps(core)});console.log(JSON.stringify({expr}))"],
                         capture_output=True, text=True, check=True).stdout
    return json.loads(out)

JS_CASES = [
    ("D.cg(60,80,88.4,false)", (140 - 60) * 80 / 72),
    ("D.cg(80,60,88.4,true)", 60 * 60 / 72 * 0.85),
    ("D.ckdepi2021(50,0.9*88.4,false)", 142 * 0.9938 ** 50),
    ("D.ckdepi2009(70,0.7*88.4,true)", 141 * 0.993 ** 70 * 1.018),
    ("D.ekfc(60,80,false)", 107.3 * 0.990 ** 20),
    ("D.ekfc(30,124,true)", 107.3 * 2 ** -1.132),
    ("D.bsaDuBois(70,175)", 0.007184 * 70 ** 0.425 * 175 ** 0.725),
]

@pytest.mark.parametrize("expr,expected", JS_CASES)
def test_js_reference(expr, expected):
    assert _js(expr) == pytest.approx(expected, rel=1e-9)

# Dose-rule boundaries (EU SmPC / US PI wording)
P = dict(age=70, weight_kg=np.array([80.0]), scr_umol=np.array([90.0]))
@pytest.mark.parametrize("drug,rf,expected", [
    ("rivaroxaban", 14.99, 0), ("rivaroxaban", 15, 1), ("rivaroxaban", 49.99, 1), ("rivaroxaban", 50, 2),
    ("edoxaban", 50, 1), ("edoxaban", 50.01, 2), ("dabigatran", 29.99, 0), ("dabigatran", 30, 1.5),
    ("dabigatran", 50.01, 2), ("apixaban", 29.99, 1), ("apixaban", 30, 2), ("apixaban", 14.99, 0)])
def test_eu_boundaries(drug, rf, expected):
    kw = dict(age=np.array([70]), weight_kg=np.array([80.0]), scr_umol=np.array([90.0]))
    assert RULES[drug](np.array([rf]), **kw)[0] == expected

def test_eu_special_rules():
    a = np.array
    assert RULES["edoxaban"](a([80.0]), age=a([70]), weight_kg=a([60.0]), scr_umol=a([90.0]))[0] == 1   # ≤60 kg
    assert RULES["edoxaban"](a([80.0]), age=a([70]), weight_kg=a([80.0]), scr_umol=a([90.0]), pgp=a([True]))[0] == 1
    assert RULES["dabigatran"](a([80.0]), age=a([80]), weight_kg=a([80.0]), scr_umol=a([90.0]))[0] == 1  # ≥80 y
    assert RULES["dabigatran"](a([80.0]), age=a([76]), weight_kg=a([80.0]), scr_umol=a([90.0]))[0] == 1.5
    assert RULES["apixaban"](a([80.0]), age=a([80]), weight_kg=a([60.0]), scr_umol=a([90.0]))[0] == 1   # 2 of 3
    assert RULES["apixaban"](a([80.0]), age=a([80]), weight_kg=a([61.0]), scr_umol=a([132.9]))[0] == 2  # 1 of 3

@pytest.mark.parametrize("drug,rf,expected", [
    ("rivaroxaban", 50, 1), ("rivaroxaban", 50.01, 2), ("edoxaban", 95, 2), ("edoxaban", 95.01, 3),
    ("dabigatran", 30, 1), ("dabigatran", 30.01, 2), ("dabigatran", 14.99, 0)])
def test_us_boundaries(drug, rf, expected):
    kw = dict(age=np.array([70]), weight_kg=np.array([80.0]), scr_umol=np.array([90.0]))
    assert RULES_US[drug](np.array([rf]), **kw)[0] == expected
