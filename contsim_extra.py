"""Extra continuous-sampling numbers for the text: Boriani CG stages, interacting-drug agreement,
Mosteller vs DuBois, Ammar direction, edoxaban weight <=60 kg share."""
import numpy as np, pandas as pd
import contsim as cs, egfr
from doac_rules import RULES
cs.calibrate()
g = cs.POP["Boriani (Europe)"]()
st = np.digitize(g.CG, [30, 45, 60, 90])
print("Boriani CG stages I..IV/V:", np.round(100*np.bincount(st, minlength=5)[::-1]/len(g), 1), "published I 30.2, II 37.6")
kw = dict(age=g.age.values, weight_kg=g.weight.values, scr_umol=g.scr.values)
for drug in ["edoxaban", "dabigatran"]:
    for pgp in [False, True]:
        ref = RULES[drug](g.CG.values, pgp=pgp, **kw)
        agr = [100*np.mean(RULES[drug](g[e].values, pgp=pgp, **kw) == ref) for e in ["CKDEPI21","EKFC","CKDEPI21_abs","EKFC_abs"]]
        print(drug, "pgp" if pgp else "none", "agreement %:", np.round(agr, 1))
# Mosteller
bsaM = egfr.bsa_mosteller(g.weight.values, g.height.values)
for drug in ["rivaroxaban", "edoxaban", "dabigatran"]:
    ref = RULES[drug](g.CG.values, **kw)
    for e in ["CKDEPI21", "EKFC"]:
        dD = 100*np.mean(RULES[drug](g[e+"_abs"].values, **kw) != ref)
        dM = 100*np.mean(RULES[drug](egfr.deindex(g[e].values, bsaM), **kw) != ref)
        print(f"{drug} {e}: DuBois {dD:.1f} Mosteller {dM:.1f} diff {dM-dD:+.1f}")
a = cs.POP["Ammar (USA)"]()
from doac_rules import rivaroxaban
d = rivaroxaban(a.CKDEPI21_abs.values) - rivaroxaban(a.CG.values)
print(f"Ammar: discordance {100*np.mean(d!=0):.1f}, higher share {100*np.mean(d>0)/np.mean(d!=0):.0f}%")
