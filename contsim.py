"""Continuous-sampling version of all cohort analyses (replaces the re-weighted grid for frequencies).

Same sampler as mc_boriani.py: patients are drawn from continuous distributions (no 5 kg / 5 µmol/L steps),
age clipped to 18-100, BMI clipped to 16-50 where weight is derived from BMI and height, weight clipped to 30-200 kg.
For each cohort the distributional specification is the one used for the grid weights (Table S1).
Where a cohort's grid model had no height density (Yao, Chan, Ammar, Cheng), height is drawn from sex-specific
normals (US: M 176±7, F 162±7; Taiwan: M 168±6, F 156±6 cm); height only affects non-indexed estimates.

Outputs: results/cont_validation_points.csv, cont_validation_cheng.csv, cont_estimator_by_population.csv,
         cont_vs_grid_boriani.csv
"""
import numpy as np, pandas as pd
import egfr
from doac_rules import RULES, rivaroxaban, apixaban

N = 1_000_000
rng = np.random.default_rng(2026)


def lognorm(med, sdlog, n):
    return np.exp(rng.normal(np.log(med), sdlog, n))


def finish(age, fem, wt, ht, scr):
    age = np.clip(age, 18, 100); wt = np.clip(wt, 30, 200)
    d = dict(age=age, female=fem, weight=wt, height=ht, scr=scr)
    d["CG"] = egfr.cockcroft_gault(age, wt, scr, fem)
    bsa = egfr.bsa_dubois(wt, ht)
    for n, fn in [("CKDEPI09", egfr.ckd_epi_2009), ("CKDEPI21", egfr.ckd_epi_2021), ("EKFC", egfr.ekfc)]:
        d[n] = fn(age, scr, fem); d[n + "_abs"] = egfr.deindex(d[n], bsa)
    return pd.DataFrame(d)


def bmi_based(age_m, age_sd, pf, hm, hsd_m, hf, hsd_f, bmi_m, bmi_sd, scr_m, scr_f, sdlog):
    fem = rng.random(N) < pf
    age = rng.normal(age_m, age_sd, N)
    ht = np.where(fem, rng.normal(hf, hsd_f, N), rng.normal(hm, hsd_m, N))
    bmi = np.clip(rng.normal(bmi_m, bmi_sd, N), 16, 50)
    wt = bmi * (ht / 100) ** 2
    scr = lognorm(np.where(fem, scr_f, scr_m), sdlog, N)
    return finish(age, fem, wt, ht, scr)


def weight_based(age_m, age_sd, pf, wt_fn, hts, scr_m, scr_f, sdlog):
    fem = rng.random(N) < pf
    age = rng.normal(age_m, age_sd, N)
    wt = wt_fn(fem)
    ht = np.where(fem, rng.normal(hts[2], hts[3], N), rng.normal(hts[0], hts[1], N))
    scr = lognorm(np.where(fem, scr_f, scr_m), sdlog, N)
    return finish(age, fem, wt, ht, scr)


US_H, TW_H = (176, 7, 162, 7), (168, 6, 156, 6)
POP = {
    "Manzano-Fernández (Spain)": lambda: bmi_based(74, 9, 0.5, 172, 7, 158, 6, 30.1, 5.1, 92, 76, np.log(1.16 / 0.79) / 1.349),
    "Yao, ORBIT-AF II (USA)": lambda: weight_based(71, 10.4, 0.414,
        lambda f: lognorm(np.where(f, 80, 95), np.log(106 / 75) / 1.349, N), US_H, YAO_M, YAO_M * 0.82, np.log(1.2 / 0.8) / 1.349),
    "Chan (Taiwan)": lambda: weight_based(71.1, 12.8, 0.43,
        lambda f: rng.normal(np.where(f, 56.0, 68.5), 12.5, N), TW_H, CHAN_M, CHAN_M * 0.8, np.log(1.45 / 0.8) / 1.349),
    # Ammar 2026 full text (Table 1): age 69.7 (11.8), 42.9% women, weight 95.6 (26.3) kg, height 1.7 (0.1) m,
    # Scr 0.96 (0.27) mg/dL; CG < 15 excluded. Sex split of weight (men 101.3, women 88.0, SD 25.5) ASSUMED;
    # Scr women = 0.8 x men ASSUMED; men's median and SD of log calibrated to the published Scr mean and SD.
    "Ammar (USA)": lambda: weight_based(69.7, 11.8, 0.429,
        lambda f: rng.normal(np.where(f, 88.0, 101.3), 25.5, N), US_H, AM_M, AM_M * 0.8, AM_S),
    "Boriani (Europe)": lambda: bmi_based(70, 10.4, 0.413, 176, 7, 163, 6, 27.8, 4.7, BOR_M, BOR_M * 0.8, BOR_S),
}


def calibrate():
    """Re-calibrate, on continuous samples, the Scr inputs that were calibrated in the original analyses:
    Yao: men's Scr median so that median CG = 82.6 mL/min; Chan: men's median so that overall Scr median =
    1.03 mg/dL (91 µmol/L); Boriani: median and sdlog so that the CKD-EPI 2009 stage distribution matches
    I 15.3, II 50.0, IIIa 20.3, IIIb 10.6, IV/V 3.8 %. Common random numbers keep the search smooth."""
    global YAO_M, CHAN_M, BOR_M, BOR_S, AM_M, AM_S
    YAO_M = CHAN_M = BOR_M = AM_M = 90.0; BOR_S = AM_S = 0.28
    st0 = rng.bit_generator.state
    best = None
    for m in np.arange(70, 110, 0.5):
        YAO_M = m; rng.bit_generator.state = st0; g = POP["Yao, ORBIT-AF II (USA)"]()
        e = abs(np.median(g.CG) - 82.6)
        if best is None or e < best[0]: best = (e, m)
    YAO_M = best[1]; best = None
    for m in np.arange(80, 120, 0.5):
        CHAN_M = m; rng.bit_generator.state = st0; g = POP["Chan (Taiwan)"]()
        e = abs(np.median(g.scr) - 91.06)
        if best is None or e < best[0]: best = (e, m)
    CHAN_M = best[1]; best = None
    tgt = np.array([3.8, 10.6, 20.3, 50.0, 15.3])
    for m in np.arange(84, 106, 1.0):
        for s in np.arange(0.20, 0.42, 0.02):
            BOR_M, BOR_S = m, s; rng.bit_generator.state = st0; g = POP["Boriani (Europe)"]()
            d = 100 * np.bincount(np.digitize(g.CKDEPI09, [30, 45, 60, 90]), minlength=5) / N
            e = np.sum((d - tgt) ** 2)
            if best is None or e < best[0]: best = (e, m, s, d)
    BOR_M, BOR_S = best[1], best[2]; bor_d = best[3]; best = None
    for m in np.arange(80, 100, 0.5):
        for sd in np.arange(0.18, 0.36, 0.01):
            AM_M, AM_S = m, sd; rng.bit_generator.state = st0; g = POP["Ammar (USA)"]()
            e = (g.scr.mean() - 0.96 * 88.4) ** 2 + (g.scr.std() - 0.27 * 88.4) ** 2
            if best is None or e < best[0]: best = (e, m, sd, g.scr.mean(), g.scr.std())
    AM_M, AM_S = best[1], best[2]
    print(f"Ammar calibration: men Scr {AM_M}, sdlog {AM_S:.2f}; Scr mean {best[3]:.1f} SD {best[4]:.1f} (target 84.9, 23.9)")
    best = (None, None, None, bor_d)
    print(f"Calibrated: Yao men Scr {YAO_M}; Chan men Scr {CHAN_M}; Boriani men Scr {BOR_M}, sdlog {BOR_S:.2f}, "
          f"stages {np.round(best[3], 1)} (target {tgt})")
    rng.bit_generator.state = st0

if __name__ == "__main__":
    calibrate()
    P = {k: f() for k, f in POP.items()}
    m = P["Manzano-Fernández (Spain)"]; print(f"Manzano overall Scr median {np.median(m.scr):.1f} (published 84)")
    kw = lambda g: dict(age=g.age.values, weight_kg=g.weight.values, scr_umol=g.scr.values)

    # ---------- validation points (same rules / subsets as validate_*.py) ----------
    pts = []
    g = P["Manzano-Fernández (Spain)"]
    def dabi_m(rf, **_): return np.select([rf < 30, rf < 50], [0, 1], default=2).astype(float)
    subs = {"all": np.ones(N, bool), "CG<60": (g.CG < 60).values, "≥75 y": (g.age >= 75).values}
    for drug, rule in [("dabigatran", dabi_m), ("rivaroxaban", rivaroxaban), ("apixaban", apixaban)]:
        d = rule(g.CKDEPI09_abs.values, **kw(g)) != rule(g.CG.values, **kw(g))
        for sn, m in subs.items():
            pts.append(("Manzano-Fernández 2015", f"{drug}, {sn}", 100 * d[m].mean()))
    g = P["Yao, ORBIT-AF II (USA)"]
    def riva_y(rf): return np.select([rf < 30, rf <= 50], [0, 1], default=2)
    d = riva_y(g.CKDEPI09.values) - riva_y(g.CG.values)
    pts += [("Yao 2023 (ORBIT-AF II)", "rivaroxaban, any", 100 * (d != 0).mean()),
            ("Yao 2023 (ORBIT-AF II)", "rivaroxaban, lower dose", 100 * (d < 0).mean()),
            ("Yao 2023 (ORBIT-AF II)", "rivaroxaban, higher dose", 100 * (d > 0).mean())]
    g = P["Chan (Taiwan)"]
    def cat(rf): return np.select([rf < 15, rf <= 50], [0, 1], default=2)
    d = cat(g.CKDEPI09.values) != cat(g.CG.values)
    for lab, m in [("all", np.ones(N, bool)), ("≥50 kg & <75 y", (g.weight >= 50) & (g.age < 75)),
                   ("<50 kg", g.weight < 50), ("≥75 y", g.age >= 75), ("<50 kg & ≥75 y", (g.weight < 50) & (g.age >= 75))]:
        pts.append(("Chan 2022", lab, 100 * d[np.asarray(m)].mean()))
    g = P["Ammar (USA)"]; g = g[g.CG >= 15]
    from doac_rules import RULES_US
    rv = RULES_US["rivaroxaban"]; ref = rv(g.CG.values)
    print(f"Ammar check (not calibrated): CG mean {g.CG.mean():.1f} SD {g.CG.std():.1f} (published 98.4, 47.7); "
          f"CKD-EPI 2021 {g.CKDEPI21.mean():.1f} (78.2), NI {g.CKDEPI21_abs.mean():.1f} (94.0); "
          f"CKD-EPI 2009 {g.CKDEPI09.mean():.1f} (74.8), NI {g.CKDEPI09_abs.mean():.1f} (89.9)")
    d21a = rv(g.CKDEPI21_abs.values) - ref
    am = [("rivaroxaban, CKD-EPI 2021 NI", 4.8, 100 * (d21a != 0).mean()),
          ("rivaroxaban, CKD-EPI 2021 NI, higher dose", 100 * 20 / 476, 100 * (d21a > 0).mean()),
          ("rivaroxaban, CKD-EPI 2021 indexed", 7.6, 100 * (rv(g.CKDEPI21.values) != ref).mean()),
          ("rivaroxaban, CKD-EPI 2009 NI", 3.8, 100 * (rv(g.CKDEPI09_abs.values) != ref).mean())]
    old = pd.read_csv("results/validation_points.csv")
    old = old[old.cohort != "Ammar 2026"]
    old = pd.concat([old, pd.DataFrame([dict(cohort="Ammar 2026", label=l, pub=round(p, 1), sim=np.nan) for l, p, _ in am])])
    pts += [("Ammar 2026", l, v) for l, _, v in am]
    V = old.merge(pd.DataFrame(pts, columns=["cohort", "label", "sim_cont"]), on=["cohort", "label"])
    V = V.rename(columns={"sim": "sim_grid"})
    V.to_csv("results/cont_validation_points.csv", index=False)
    for col in ["sim_grid", "sim_cont"]:
        dd = V[col] - V.pub; sd = dd.std(ddof=1)
        print(f"{col}: mean diff {dd.mean():.2f}, SD {sd:.2f}, LoA {dd.mean()-1.96*sd:.1f} to {dd.mean()+1.96*sd:.1f}, "
              f"mean |d| {dd.abs().mean():.2f}")
    print(V.assign(dc=V.sim_cont - V.pub).groupby("cohort").dc.mean().round(2).to_string())
    print(V.round(1).to_string(index=False))

    # ---------- Cheng: calibrate Scr to CG mean/SD, compare means within CG stage ----------
    best = None
    base = dict(fem=rng.random(N) < 0.41, age=rng.normal(75.9, 11.0, N), wt=rng.normal(66.5, 13.0, N), z=rng.normal(0, 1, N))
    ht = np.where(base["fem"], rng.normal(156, 6, N), rng.normal(168, 6, N))
    for m in np.arange(80, 131, 2.5):
        for s in np.arange(0.25, 0.50, 0.025):
            scr = np.exp(np.log(np.where(base["fem"], m * 0.78, m)) + s * base["z"])
            cg = egfr.cockcroft_gault(np.clip(base["age"], 18, 100), np.clip(base["wt"], 30, 200), scr, base["fem"])
            e = (cg.mean() - 56.4) ** 2 + (cg.std() - 25.9) ** 2
            if best is None or e < best[0]: best = (e, m, s)
    _, m, s = best
    scr = np.exp(np.log(np.where(base["fem"], m * 0.78, m)) + s * base["z"])
    g = finish(base["age"], base["fem"], base["wt"], ht, scr)
    print(f"Cheng calibration: Scr men {m}, sdlog {s:.3f}; CG {g.CG.mean():.1f} ± {g.CG.std():.1f}")
    st = pd.cut(g.CG, [-1, 15, 30, 45, 60, 1e9], right=False, labels=["<15", "15–29", "30–44", "45–59", "≥60"])
    C = g.groupby(st, observed=True)[["CKDEPI21", "EKFC"]].mean()
    oldB = pd.read_csv("results/validation_points_cheng.csv")
    oldB["sim_cont"] = [C.loc[r.stage, "CKDEPI21" if r.equation == "CKD-EPI" else "EKFC"] for r in oldB.itertuples()]
    oldB = oldB.rename(columns={"sim": "sim_grid"}); oldB.to_csv("results/cont_validation_cheng.csv", index=False)
    dd = oldB.sim_cont - oldB.pub; sd = dd.std(ddof=1)
    print(oldB.round(1).to_string(index=False))
    print(f"Cheng cont: mean diff {dd.mean():.2f}, LoA {dd.mean()-1.96*sd:.1f} to {dd.mean()+1.96*sd:.1f}, mean |d| {dd.abs().mean():.2f}")

    # ---------- Table 4 with continuous sampling ----------
    rows = []
    for pop, g in P.items():
        for drug in ["rivaroxaban", "edoxaban", "dabigatran"]:
            ref = RULES[drug](g.CG.values, **kw(g)); r = dict(population=pop, drug=drug)
            for e in ["CKDEPI09", "CKDEPI21", "EKFC", "CKDEPI21_abs", "EKFC_abs"]:
                dd = RULES[drug](g[e].values, **kw(g)) - ref
                r[e] = f"{100*(dd!=0).mean():.1f} ({100*(dd>0).mean():.1f})"
            rows.append(r)
    T = pd.DataFrame(rows); T.to_csv("results/cont_estimator_by_population.csv", index=False)
    pd.set_option("display.width", 220); print(T.to_string(index=False))

    # ---------- Boriani: grid vs continuous with identical parameters ----------
    g = P["Boriani (Europe)"]
    print(f"Boriani continuous: weight ≤60 kg {100*(g.weight<=60).mean():.1f}%, <50 kg {100*(g.weight<50).mean():.1f}%, "
          f"Scr <40 {100*(g.scr<40).mean():.1f}%, CG>95 {100*(g.CG>95).mean():.1f}%")

    # ---------- Chan: age–weight correlation sensitivity ----------
    rows = []
    for rho in [0, -0.2, -0.3, -0.4]:
        fem = rng.random(N) < 0.43
        z1, z2 = rng.normal(size=N), rng.normal(size=N)
        age = np.clip(71.1 + 12.8 * z1, 18, 100)
        wt = np.where(fem, 56.0, 68.5) + 12.5 * (rho * z1 + np.sqrt(1 - rho ** 2) * z2)
        ht = np.where(fem, rng.normal(156, 6, N), rng.normal(168, 6, N))
        g = finish(age, fem, wt, ht, lognorm(np.where(fem, CHAN_M * 0.8, CHAN_M), np.log(1.45 / 0.8) / 1.349, N))
        for est in ["CKDEPI09", "EKFC_abs"]:
            d = cat(g[est].values) != cat(g.CG.values)
            for lab, msk in [("all", np.ones(N, bool)), ("≥50 kg & <75 y", (g.weight >= 50) & (g.age < 75)),
                             ("<50 kg", g.weight < 50), ("≥75 y", g.age >= 75), ("<50 kg & ≥75 y", (g.weight < 50) & (g.age >= 75))]:
                rows.append(dict(rho=rho, est=est, subgroup=lab, sim=round(100 * d[np.asarray(msk)].mean(), 1)))
    R = pd.DataFrame(rows); R.to_csv("results/cont_chan_corr.csv", index=False)
    print(R.pivot_table(index=["est", "subgroup"], columns="rho", values="sim").to_string())
