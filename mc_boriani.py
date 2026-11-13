"""Probabilistic sensitivity analysis of the European scenario (Boriani-like).
Each draw: sample the ASSUMED parameters (uniform ranges below), re-calibrate the Scr distribution
to the published CKD-EPI stage distribution (I 15.3, II 50.0, IIIa 20.3, IIIb 10.6, IV/V 3.8 %)
on an (age, sex, Scr) grid, then simulate 200,000 virtual patients and compute discordance.
Ranges: age SD 9–12 y; BMI mean 27.3–28.3, SD 4.0–5.5 kg/m2; height M 173–179 / F 160–166 cm, SD 6–8;
female:male Scr median ratio 0.75–0.85; female share 0.40–0.43 (reported 41.3%). Age mean fixed at 70.
Output: median and 2.5–97.5th percentiles over 300 draws."""
import numpy as np, pandas as pd
from scipy.stats import norm
import egfr
from doac_rules import RULES, RULES_US
rng = np.random.default_rng(42)
# calibration grid
A, S, C = np.meshgrid(np.arange(18, 101), [0, 1], np.arange(30, 801, 2.5), indexing="ij")
A, S, C = A.ravel(), S.ravel().astype(bool), C.ravel()
st = np.digitize(egfr.ckd_epi_2009(A, C, S), [30, 45, 60, 90])  # 0 IV/V, 1 IIIb, 2 IIIa, 3 II, 4 I
target = np.array([3.8, 10.6, 20.3, 50.0, 15.3])
EST = ["CKDEPI21", "EKFC", "CKDEPI21_abs", "EKFC_abs"]
N, M = 300, 200_000
rows, params = [], []
for i in range(N):
    p = dict(age_sd=rng.uniform(9, 12), bmi_m=rng.uniform(27.3, 28.3), bmi_sd=rng.uniform(4.0, 5.5),
             hm=rng.uniform(173, 179), hf=rng.uniform(160, 166), hsd=rng.uniform(6, 8),
             fr=rng.uniform(0.75, 0.85), pf=rng.uniform(0.40, 0.43))
    wa = norm.pdf(A, 70, p["age_sd"]) * np.where(S, p["pf"], 1 - p["pf"])
    best = None
    for m in np.arange(70, 110, 1.0):
        for s in np.arange(0.16, 0.56, 0.02):
            w = wa * norm.pdf(np.log(C), np.log(np.where(S, m * p["fr"], m)), s) / C
            d = np.bincount(st, weights=w, minlength=5); d = 100 * d / d.sum()
            e = np.sum((d - target) ** 2)
            if best is None or e < best[0]: best = (e, m, s)
    p.update(scr_m=best[1], scr_sdlog=best[2], calib_rmse=np.sqrt(best[0] / 5))
    # simulate patients
    fem = rng.random(M) < p["pf"]
    age = np.clip(rng.normal(70, p["age_sd"], M), 18, 100)
    ht = np.where(fem, rng.normal(p["hf"], p["hsd"], M), rng.normal(p["hm"], p["hsd"], M))
    bmi = np.clip(rng.normal(p["bmi_m"], p["bmi_sd"], M), 16, 50)
    wt = bmi * (ht / 100) ** 2
    scr = np.exp(rng.normal(np.log(np.where(fem, p["scr_m"] * p["fr"], p["scr_m"])), p["scr_sdlog"]))
    est = {"CG": egfr.cockcroft_gault(age, wt, scr, fem)}
    bsa = egfr.bsa_dubois(wt, ht)
    for n, fn in [("CKDEPI21", egfr.ckd_epi_2021), ("EKFC", egfr.ekfc)]:
        est[n] = fn(age, scr, fem); est[n + "_abs"] = egfr.deindex(est[n], bsa)
    kw = dict(age=age, weight_kg=wt, scr_umol=scr)
    r = {"draw": i}
    for reg, rules in [("EU", RULES), ("US", RULES_US)]:
        for drug, rule in rules.items():
            ref = rule(est["CG"], **kw)
            for e in EST:
                alt = rule(est[e], **kw)
                if reg == "US" and drug == "edoxaban":
                    ch = (ref == 3) != (alt == 3); d = np.where((ref == 3) | (alt == 3), 0, alt - ref)
                else:
                    ch = np.zeros(M, bool); d = alt - ref
                r[f"{reg}|{drug}|{e}|any"] = 100 * np.mean((d != 0) | ch)
                r[f"{reg}|{drug}|{e}|hi"] = 100 * np.mean(d > 0)
                r[f"{reg}|{drug}|{e}|lo"] = 100 * np.mean(d < 0)
    for drug in RULES:
        eu = RULES[drug](est["CG"], **kw); us = RULES_US[drug](est["CG"], **kw)
        if drug == "edoxaban": us = np.where(us == 3, 0, us)
        diff = eu != us; poss = diff & (eu == 1.5)
        r[f"label|{drug}|definite"] = 100 * np.mean(diff & ~poss); r[f"label|{drug}|possible"] = 100 * np.mean(poss)
        if drug == "edoxaban": r["label|edoxaban|CG>95"] = 100 * np.mean(est["CG"] > 95); r["label|edoxaban|wt<=60_only"] = 100*np.mean(diff & (est["CG"] <= 95))
        if drug == "dabigatran": r["label|dabigatran|age>=80"] = 100 * np.mean(diff & (eu == 1))
    rows.append(r); params.append(p)
R = pd.DataFrame(rows); P = pd.DataFrame(params)
R.to_csv("results/mc_boriani_draws.csv", index=False); P.to_csv("results/mc_boriani_params.csv", index=False)
q = R.drop(columns="draw").quantile([0.5, 0.025, 0.975]).T.round(1); q.columns = ["median", "p2.5", "p97.5"]
q.to_csv("results/mc_boriani_summary.csv")
print(P.describe().round(2).loc[["mean", "min", "max"]].to_string()); pd.set_option("display.max_rows", 200); print(q.to_string())
