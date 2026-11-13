"""Software verification: Python reference (egfr.py + doac_rules.py, used for the paper)
vs independent JavaScript implementation (doac_core.js, used by the web calculator).

1) 50 000 random virtual patients (continuous values, wide ranges, random P-gp/verapamil)
2) boundary cases: Scr solved so that CG lands exactly on every dose threshold (±0.01 mL/min)
Compares every renal estimate (max relative difference) and every dose class (EU + US).
Usage: python cross_validate.py   (needs node on PATH)
"""
import json, os, subprocess, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import egfr
from doac_rules import RULES, RULES_US

rng = np.random.default_rng(20261112)
N = 50_000
p = pd.DataFrame({
    "age": rng.uniform(18, 100, N), "weight": rng.uniform(35, 160, N),
    "height": rng.uniform(140, 205, N), "scr": np.exp(rng.uniform(np.log(30), np.log(700), N)),
    "female": rng.integers(0, 2, N).astype(bool), "pgp": rng.random(N) < 0.15,
    "verapamil": rng.random(N) < 0.10})
# boundary cases: CG exactly at thresholds ± 0.01
rows = []
for thr in [15, 30, 49.99, 50, 50.01, 95, 95.01]:
    for d in [-0.01, 0, 0.01]:
        for age in [45, 76, 82]:
            for wt in [55, 60, 61, 90]:
                for fem in [False, True]:
                    target = thr + d
                    scr = (140 - age) * wt / (72 * target) * (0.85 if fem else 1) * 88.4
                    rows.append(dict(age=age, weight=wt, height=170, scr=scr, female=fem, pgp=False, verapamil=False))
p = pd.concat([p, pd.DataFrame(rows)], ignore_index=True)
# apixaban creatinine boundary
p = pd.concat([p, pd.DataFrame([dict(age=80, weight=60, height=165, scr=s, female=True, pgp=False, verapamil=False)
                                for s in [132.99, 133.0, 133.01]])], ignore_index=True)

# ---- Python reference ----
f = p.female.values
py = pd.DataFrame({"CG": egfr.cockcroft_gault(p.age, p.weight, p.scr, f)})
bD = egfr.bsa_dubois(p.weight, p.height); bM = egfr.bsa_mosteller(p.weight, p.height)
for n, fn in [("CKDEPI09", egfr.ckd_epi_2009), ("CKDEPI21", egfr.ckd_epi_2021), ("EKFC", egfr.ekfc)]:
    py[n] = fn(p.age.values, p.scr.values, f)
    py[n + "_abs"] = egfr.deindex(py[n], bD); py[n + "_absM"] = egfr.deindex(py[n], bM)
EST = list(py.columns)
kw = dict(age=p.age.values, weight_kg=p.weight.values, scr_umol=p.scr.values)
pydose = {}
for region, rules in [("EU", RULES), ("US", RULES_US)]:
    for drug, rule in rules.items():
        for e in EST:
            extra = {}
            if region == "EU" and drug == "edoxaban": extra = dict(pgp=p.pgp.values)
            if region == "EU" and drug == "dabigatran": extra = dict(pgp=p.verapamil.values)
            pydose[(region, drug, e)] = rule(py[e].values, **kw, **extra)

# ---- JavaScript ----
inp = os.path.join(HERE, "_cases.json"); out = os.path.join(HERE, "_js.json")
p.to_json(inp, orient="records", double_precision=15)  # full precision: boundary cases sit exactly on thresholds
js = f"""
const D=require({json.dumps(os.path.join(HERE,'doac_core.js'))});const fs=require('fs');
const cases=JSON.parse(fs.readFileSync({json.dumps(inp)}));
const r=cases.map(c=>({{EU:D.doses(c,'EU'),US:D.doses(c,'US')}}));
fs.writeFileSync({json.dumps(out)},JSON.stringify(r));"""
subprocess.run(["node", "-e", js], check=True)
jr = json.load(open(out))

# ---- compare ----
report = []
for e in EST:
    jv = np.array([r["EU"]["estimates"][e] for r in jr])
    rel = np.max(np.abs(jv - py[e].values) / np.abs(py[e].values))
    report.append(dict(item=f"estimate {e}", n=len(p), max_rel_diff=f"{rel:.1e}", agreement="—"))
mism_total = 0
for (region, drug, e), pv in pydose.items():
    jv = np.array([r[region]["doses"][drug][e] for r in jr], dtype=float)
    mism = int(np.sum(jv != pv)); mism_total += mism
    report.append(dict(item=f"dose {region} {drug} [{e}]", n=len(p), max_rel_diff="—",
                       agreement=f"{100*(1-mism/len(p)):.3f}% ({mism} mismatches)"))
rep = pd.DataFrame(report)
rep.to_csv(os.path.join(os.path.dirname(HERE), "results", "software_verification.csv"), index=False)
print(f"cases: {len(p)} (random {N}, boundary {len(p)-N})")
print(rep[rep.item.str.startswith("estimate")].to_string(index=False))
print(f"dose comparisons: {len(pydose)} rules×estimators × {len(p)} cases = {len(pydose)*len(p):,}; total mismatches: {mism_total}")
os.remove(inp); os.remove(out)
