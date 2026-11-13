"""DOAC dose classes for stroke prevention in NVAF.

`rf` = renal function value plugged into the rule (CG CrCl, or an eGFR used
in its place). Ordinal dose level:
    0 = not recommended / contraindicated
    1 = reduced dose
    1.5 = dabigatran EU "150 or 110 by individual assessment"
    2 = standard dose
US edoxaban adds level 3 = "do not use, CrCl > 95" (not ordinal; handled
separately in sim.py).

EU rules (SmPC 4.2):
- Xarelto: CrCl 15-49 -> 15 mg od; <15 not recommended; else 20 mg od.
- Lixiana: CrCl 15-50 OR weight <=60 kg OR P-gp inhibitor (ciclosporin,
  dronedarone, erythromycin, ketoconazole) -> 30 mg od; <15 NR; else 60 mg.
- Pradaxa: CrCl <30 contraindicated; age >=80 or verapamil -> 110 bid;
  age 75-80 or CrCl 30-50 -> 150 or 110 individually.
- Eliquis: >=2 of {age>=80, weight<=60, Scr>=133 µmol/L} -> 2.5 bid;
  CrCl 15-29 -> 2.5 bid; <15 NR.
US rules (USPI):
- Xarelto: CrCl <=50 -> 15 mg; >50 -> 20 mg.
- Savaysa (NVAF): >95 do not use; >50-95 -> 60 mg; 15-50 -> 30 mg; <15 NR.
  (weight / P-gp reduction applies to DVT/PE only.)
- Pradaxa: >30 -> 150 bid; 15-30 -> 75 bid; <15 no recommendation.
- Eliquis: >=2 of 3 criteria only; no CrCl-based reduction.
Not modelled: bleeding-risk discretion.
"""
import numpy as np

LABELS = {
    "rivaroxaban": {0: "NR", 1: "15 mg", 2: "20 mg"},
    "edoxaban": {0: "NR", 1: "30 mg", 2: "60 mg", 3: "NR (>95)"},
    "dabigatran": {0: "CI", 1: "110 mg x2", 1.5: "150/110 x2", 2: "150 mg x2"},
    "apixaban": {0: "NR", 1: "2.5 mg x2", 2: "5 mg x2"},
}
LABELS_US = dict(LABELS, dabigatran={0: "NR", 1: "75 mg x2", 2: "150 mg x2"})


# ---------------- EU ----------------
def rivaroxaban(rf, **_):
    return np.select([rf < 15, rf < 50], [0, 1], default=2).astype(float)


def edoxaban(rf, weight_kg, pgp=False, **_):
    lvl = np.select([rf < 15, rf <= 50], [0, 1], default=2).astype(float)
    return np.where((lvl == 2) & ((weight_kg <= 60) | pgp), 1.0, lvl)


def dabigatran(rf, age, pgp=False, **_):          # pgp here = verapamil
    lvl = np.full(np.shape(rf), 2.0)
    lvl = np.where((rf <= 50) | ((age >= 75) & (age < 80)), 1.5, lvl)
    lvl = np.where((age >= 80) | pgp, 1.0, lvl)
    return np.where(rf < 30, 0.0, lvl)


def apixaban(rf, age, weight_kg, scr_umol, **_):
    n_crit = (age >= 80).astype(int) + (weight_kg <= 60).astype(int) \
        + (scr_umol >= 133).astype(int)
    lvl = np.where(n_crit >= 2, 1.0, 2.0)
    lvl = np.where(rf < 30, 1.0, lvl)
    return np.where(rf < 15, 0.0, lvl)


# ---------------- US ----------------
def rivaroxaban_us(rf, **_):
    return np.where(rf <= 50, 1.0, 2.0)


def edoxaban_us(rf, **_):
    return np.select([rf < 15, rf <= 50, rf <= 95], [0, 1, 2], default=3).astype(float)


def dabigatran_us(rf, **_):
    return np.select([rf < 15, rf <= 30], [0, 1], default=2).astype(float)


def apixaban_us(rf, age, weight_kg, scr_umol, **_):
    n_crit = (age >= 80).astype(int) + (weight_kg <= 60).astype(int) \
        + (scr_umol >= 133).astype(int)
    return np.where(n_crit >= 2, 1.0, 2.0)


RULES = {"rivaroxaban": rivaroxaban, "edoxaban": edoxaban,
         "dabigatran": dabigatran, "apixaban": apixaban}
RULES_US = {"rivaroxaban": rivaroxaban_us, "edoxaban": edoxaban_us,
            "dabigatran": dabigatran_us, "apixaban": apixaban_us}
