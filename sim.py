"""Grid simulation: DOAC dose by Cockcroft-Gault (reference) vs by eGFR.
Usage: python sim.py  -> results/*.csv, figures/*

Weighting: indexed eGFR does not depend on height, and the number of valid
heights (BMI 16-50) differs by weight. Each grid row is weighted 1/n_heights
(of its weight) so every (age, weight, sex, Scr) combination counts once.
"""
import itertools
import os
import numpy as np
import pandas as pd
import egfr
from doac_rules import RULES, RULES_US, LABELS

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

EQ = {"CKDEPI09": egfr.ckd_epi_2009, "CKDEPI21": egfr.ckd_epi_2021, "EKFC": egfr.ekfc}
BSA = {"DuBois": egfr.bsa_dubois, "Mosteller": egfr.bsa_mosteller}


def build_grid():
    age = np.arange(18, 101, 1)
    wt = np.arange(40, 151, 5)
    ht = np.arange(145, 201, 5)
    sex = np.array([0, 1])                     # 1 = female
    scr = np.arange(40, 401, 5)                # µmol/L
    g = pd.DataFrame(list(itertools.product(age, wt, ht, sex, scr)),
                     columns=["age", "weight", "height", "female", "scr"])
    bmi = g.weight / (g.height / 100) ** 2
    g = g[(bmi >= 16) & (bmi <= 50)].reset_index(drop=True)
    g["w"] = 1.0 / g.groupby(["age", "weight", "female", "scr"]).height.transform("size")
    return g


def add_renal(g):
    f = g.female.astype(bool).values
    g["CG"] = egfr.cockcroft_gault(g.age, g.weight, g.scr, f)
    for name, fn in EQ.items():
        g[name] = fn(g.age.values, g.scr.values, f)
        for b, bfn in BSA.items():
            suffix = "_abs" if b == "DuBois" else "_absM"
            g[name + suffix] = egfr.deindex(g[name], bfn(g.weight, g.height))
    return g


EST = [e + s for s in ["", "_abs", "_absM"] for e in EQ]


def wmean(mask, w):
    return float(np.sum(w * mask) / np.sum(w) * 100)


def compare(g, rules, scenario, pgp=False):
    w = g.w.values
    kw = dict(age=g.age.values, weight_kg=g.weight.values, scr_umol=g.scr.values, pgp=pgp)
    rows = []
    for drug, rule in rules.items():
        ref = rule(g.CG.values, **kw)
        for est in EST:
            alt = rule(g[est].values, **kw)
            row = dict(scenario=scenario, drug=drug, estimator=est)
            if scenario.startswith("US") and drug == "edoxaban":
                # level 3 (CrCl>95, do not use) is not ordinal
                r3, a3 = ref == 3, alt == 3
                row.update(concordant=wmean(alt == ref, w),
                           used_though_CG_gt95=wmean(r3 & ~a3, w),
                           excluded_though_CG_le95=wmean(~r3 & a3, w))
                both = ~r3 & ~a3
                d = np.where(both, alt - ref, 0)
            else:
                d = alt - ref
                row.update(concordant=wmean(d == 0, w))
            row.update(higher_dose=wmean(d > 0, w), lower_dose=wmean(d < 0, w),
                       CI_breached=wmean((ref == 0) & (alt > 0) & (alt != 3), w))
            rows.append(row)
            if scenario == "EU" and not pgp:
                g[f"{drug}_{est}_d"] = d
        if scenario == "EU" and not pgp:
            g[f"{drug}_CG"] = ref
    return pd.DataFrame(rows)


def bands(g):
    g["wband"] = pd.cut(g.weight, [39, 60, 80, 100, 150], labels=["40-60", "65-80", "85-100", ">100"])
    g["aband"] = pd.cut(g.age, [17, 64, 79, 100], labels=["18-64", "65-79", ">=80"])
    rows = []
    for drug in RULES:
        for est in ["CKDEPI21", "EKFC", "CKDEPI21_abs", "EKFC_abs"]:
            d = g[f"{drug}_{est}_d"]
            for sexv, sexn in [(0, "male"), (1, "female"), (None, "both")]:
                m = g if sexv is None else g[g.female == sexv]
                dd = d.loc[m.index]
                t = pd.DataFrame({"a": m.aband, "w_": m.wband, "w": m.w,
                                  "hi": (dd > 0) * m.w, "lo": (dd < 0) * m.w})
                s = t.groupby(["a", "w_"], observed=True)[["hi", "lo", "w"]].sum()
                s["higher_dose"] = (s.hi / s.w * 100).round(1)
                s["lower_dose"] = (s.lo / s.w * 100).round(1)
                s = s[["higher_dose", "lower_dose"]].reset_index()
                s.columns = ["age_band", "weight_band", "higher_dose", "lower_dose"]
                s["drug"], s["estimator"], s["sex"] = drug, est, sexn
                rows.append(s)
    return pd.concat(rows)


def profiles():
    cases = [  # age, weight, height, female, scr
        (85, 55, 160, 1, 90), (85, 55, 160, 1, 110), (80, 50, 155, 1, 80),
        (60, 120, 180, 0, 110), (45, 140, 180, 0, 100), (75, 70, 170, 0, 130)]
    df = add_renal(pd.DataFrame(cases, columns=["age", "weight", "height", "female", "scr"]))
    kw = dict(age=df.age.values, weight_kg=df.weight.values, scr_umol=df.scr.values)
    for drug, rule in RULES.items():
        for est in ["CG", "CKDEPI21", "EKFC", "CKDEPI21_abs", "EKFC_abs"]:
            df[f"{drug}_{est}"] = [LABELS[drug][v] for v in rule(df[est].values, **kw)]
    return df


def figure1(g):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    cmap = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#e34948"])
    ests = [("CKDEPI21", "CKD-EPI 2021"), ("EKFC", "EKFC")]
    drugs = list(RULES)
    fig, axes = plt.subplots(len(drugs), len(ests), figsize=(7.2, 9.6), sharex=True, sharey=True)
    for i, drug in enumerate(drugs):
        for j, (est, title) in enumerate(ests):
            d = g[f"{drug}_{est}_d"]
            t = pd.DataFrame({"age": g.age, "weight": g.weight, "w": g.w,
                              "net": ((d > 0).astype(float) - (d < 0).astype(float)) * g.w})
            s = t.groupby(["weight", "age"])[["net", "w"]].sum()
            mat = (s.net / s.w * 100).unstack("age")
            ax = axes[i, j]
            im = ax.imshow(mat.values, origin="lower", aspect="auto", cmap=cmap,
                           norm=TwoSlopeNorm(0, -100, 100),
                           extent=[17.5, 100.5, 37.5, 152.5], interpolation="nearest")
            if i == 0:
                ax.set_title(title, fontsize=10)
            if j == 0:
                ax.set_ylabel(f"{drug.capitalize()}\nWeight (kg)", fontsize=9)
            if i == len(drugs) - 1:
                ax.set_xlabel("Age (years)", fontsize=9)
            ax.tick_params(labelsize=8, color="#888", length=2)
            for sp in ax.spines.values():
                sp.set_color("#bbb")
    cb = fig.colorbar(im, ax=axes, shrink=0.5, pad=0.03)
    cb.set_label("Net discordance (% of Scr × sex grid)\n← eGFR lower dose   eGFR higher dose →", fontsize=8)
    cb.ax.tick_params(labelsize=8)
    fig.savefig(os.path.join(FIG, "Figure_1.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "Figure_1.pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    g = add_renal(build_grid())
    print("grid points:", len(g), "| unique (age,wt,sex,scr):", round(g.w.sum()))
    res = pd.concat([compare(g, RULES, "EU"),
                     compare(g, RULES, "EU_Pgp_or_verapamil", pgp=True)[lambda x: x.drug.isin(["edoxaban", "dabigatran"])],
                     compare(g, RULES_US, "US")])
    res.round(1).to_csv(os.path.join(OUT, "discordance_summary.csv"), index=False)
    b = bands(g)
    b.to_csv(os.path.join(OUT, "discordance_by_age_weight_sex.csv"), index=False)
    profiles().round(1).to_csv(os.path.join(OUT, "profiles.csv"), index=False)
    figure1(g)
    pd.set_option("display.width", 200)
    print(res.round(1).to_string(index=False))
