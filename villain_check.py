"""Additional check: Villain et al. 2023 (Clin Kidney J 16:2661), Berlin Initiative Study, n = 224 AF patients >=70 y.
Published: age median 87 (IQR 81-91), 51% women, BMI 27 (25-30), CG median 45 (36-60) mL/min,
CKD-EPI cr 56 (41-68), EKFC cr 50 (37-60) mL/min/1.73 m2; drug mix rivaroxaban 44%, apixaban 38%,
edoxaban 9%, dabigatran 8%; EMA (EU) rules; eGFR de-indexed with DuBois BSA.
Published discrepancy in dose status vs CG: 7% with de-indexed EKFC (creatinine).
Inputs: age ~ N(86, 7.4) truncated to 70-100; BMI ~ N(27.5, 3.7); heights assumed (EU); Scr women = 0.8 x men,
men's median and SD of log calibrated to the published CG median and IQR. Discordance is computed as a
different EU dose class (not dose status of the prescribed dose), weighted by the drug mix."""
import numpy as np
import contsim as cs
from doac_rules import RULES
N = cs.N; rng = cs.rng
def pop(m, s):
    rng.bit_generator.state = st
    fem = rng.random(N) < 0.51
    age = np.clip(rng.normal(86, 7.4, N), 70, 100)
    ht = np.where(fem, rng.normal(163, 6, N), rng.normal(176, 7, N))
    bmi = np.clip(rng.normal(27.5, 3.7, N), 16, 50); wt = bmi * (ht / 100) ** 2
    scr = np.exp(np.log(np.where(fem, 0.8 * m, m)) + s * rng.normal(size=N))
    return cs.finish(age, fem, wt, ht, scr)
st = rng.bit_generator.state
best = None
for m in np.arange(80, 120, 1.0):
    for s in np.arange(0.20, 0.45, 0.02):
        g = pop(m, s); q = np.percentile(g.CG, [25, 50, 75])
        e = np.sum((q - [36, 45, 60]) ** 2)
        if best is None or e < best[0]: best = (e, m, s, q)
g = pop(best[1], best[2])
print(f"calibration: men Scr {best[1]}, sdlog {best[2]:.2f}; CG quartiles {np.round(best[3],1)} (published 36/45/60)")
print(f"check (not calibrated): CKD-EPI 2009 median {np.median(g.CKDEPI09):.1f} (published CKD-EPI cr 56), "
      f"EKFC median {np.median(g.EKFC):.1f} (published 50)")
mix = {"rivaroxaban": .44, "apixaban": .38, "edoxaban": .09, "dabigatran": .08}
kw = dict(age=g.age.values, weight_kg=g.weight.values, scr_umol=g.scr.values)
for est in ["EKFC_abs", "CKDEPI09_abs", "CKDEPI21_abs"]:
    tot = sum(w * np.mean(RULES[d](g[est].values, **kw) != RULES[d](g.CG.values, **kw)) for d, w in mix.items()) / sum(mix.values())
    print(f"{est}: discordance {100*tot:.1f}%")
