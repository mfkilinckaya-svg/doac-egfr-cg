"""Figure 2: external validation – simulated vs published values (main model, independent margins)."""
import matplotlib as _m; _m.rcParams["axes.unicode_minus"] = True
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
A = pd.read_csv("results/cont_validation_points.csv").rename(columns={"sim_cont": "sim"})
B = pd.read_csv("results/cont_validation_cheng.csv").rename(columns={"sim_cont": "sim"})
def ba(x, y):
    d = y - x; sd = d.std(ddof=1); return d.mean(), d.mean() - 1.96 * sd, d.mean() + 1.96 * sd
def ccc(x, y):
    # Lin's CCC: all moments with denominator n (ddof=0). The earlier version mixed pandas .var() (ddof=1)
    # with np.cov(bias=True) (ddof=0), which shrank CCC by (n-1)/n.
    x, y = np.asarray(x, float), np.asarray(y, float)
    return 2*np.cov(x, y, bias=True)[0,1] / (x.var() + y.var() + (x.mean()-y.mean())**2)
STY = {"Manzano-Fernández 2015":("#2a78d6","o"),"Yao 2023 (ORBIT-AF II)":("#eb6834","s"),
       "Chan 2022":("#1baf7a","^"),"Ammar 2026":("#eda100","D")}
ink, muted, grid = "#1c2328", "#5d6872", "#e3e7ea"
fig, ax = plt.subplots(1, 2, figsize=(9.6, 4.9))
for a in ax:
    for s in ["top","right"]: a.spines[s].set_visible(False)
    for s in ["left","bottom"]: a.spines[s].set_color("#b9c1c8")
    a.tick_params(colors=muted, labelsize=8); a.grid(color=grid, lw=0.6); a.set_axisbelow(True)
a = ax[0]; a.plot([0,55],[0,55], color="#9aa6b0", lw=1, ls="--", zorder=1)
for c,(col,m) in STY.items():
    d = A[A.cohort==c]; a.scatter(d.pub, d.sim, s=42, marker=m, color=col, edgecolor="white", lw=1.2, label=c, zorder=3)
a.set_xlim(0,55); a.set_ylim(0,55); a.set_aspect("equal")
a.set_xlabel("Published dose discordance (%)", color=ink, fontsize=9); a.set_ylabel("Simulated dose discordance (%)", color=ink, fontsize=9)
mad = np.mean(np.abs(A.sim-A.pub)); m_, lo_, hi_ = ba(A.pub, A.sim)
a.set_title(f"A  Dose discordance vs Cockcroft–Gault\n{len(A)} estimates, 4 cohorts; mean difference {m_:.1f}\n(95% LoA {lo_:.1f} to {hi_:.1f}) points; CCC {ccc(A.pub,A.sim):.2f}",
                                                  loc="left", fontsize=9, color=ink)
a.legend(frameon=False, fontsize=7.5, loc="upper left", labelcolor=ink)
b = ax[1]; b.plot([20,95],[20,95], color="#9aa6b0", lw=1, ls="--", zorder=1)
for eq,(col,m) in {"CKD-EPI":("#2a78d6","o"),"EKFC":("#eb6834","s")}.items():
    d = B[B.equation==eq]; b.scatter(d.pub, d.sim, s=42, marker=m, color=col, edgecolor="white", lw=1.2, label=eq, zorder=3)
    for _,r in d.iterrows(): b.annotate("CG "+r.stage, (r.pub, r.sim), xytext=(5,-9), textcoords="offset points", fontsize=6.5, color=muted)
b.set_xlim(20,95); b.set_ylim(20,95); b.set_aspect("equal")
b.set_xlabel("Published mean eGFR (mL/min/1.73 m²)", color=ink, fontsize=9); b.set_ylabel("Simulated mean eGFR (mL/min/1.73 m²)", color=ink, fontsize=9)
madb = np.mean(np.abs(B.sim-B.pub)); m_, lo_, hi_ = ba(B.pub, B.sim)
b.set_title(f"B  Mean eGFR within CG stage, Cheng 2026\n8 estimates; mean difference {m_:.1f}\n(95% LoA {lo_:.1f} to {hi_:.1f}) mL/min/1.73 m²; CCC {ccc(B.pub,B.sim):.2f}",
                                                   loc="left", fontsize=9, color=ink)
b.legend(frameon=False, fontsize=7.5, loc="upper left", labelcolor=ink)
fig.tight_layout(); fig.savefig("figures/Figure_2_validation.png", dpi=300, bbox_inches="tight", pad_inches=0.15); fig.savefig("figures/Figure_2_validation.pdf", bbox_inches="tight", pad_inches=0.15)
print(f"A: CCC {ccc(A.pub,A.sim):.3f}, mean |diff| {mad:.2f}; without apixaban CCC {ccc(A[~A.label.str.startswith('apix')].pub, A[~A.label.str.startswith('apix')].sim):.3f}")
print(f"B: CCC {ccc(B.pub,B.sim):.3f}, mean |diff| {madb:.2f}")
