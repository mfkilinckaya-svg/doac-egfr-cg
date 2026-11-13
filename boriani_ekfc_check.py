"""Uncalibrated check of EKFC in the European scenario (Boriani et al. 2025, Eur J Intern Med 136:86).
Scr was calibrated only to the CKD-EPI 2009 stage distribution (men's median 95 umol/L, SD of log 0.28; see contsim.py).
Boriani et al. report that EKFC placed 52.21% of patients in stage II (60-89) and CG placed 37.58% in stage II
and 30.19% in stage I; these were not used for calibration and are compared here."""
import numpy as np
import contsim as cs
cs.YAO_M = cs.CHAN_M = cs.AM_M = 90.0; cs.AM_S = 0.25
cs.BOR_M, cs.BOR_S = 95.0, 0.28
g = cs.POP["Boriani (Europe)"]()
for e in ["CKDEPI09", "CKDEPI21", "EKFC", "CG"]:
    st = np.digitize(g[e], [30, 45, 60, 90])
    d = 100 * np.bincount(st, minlength=5)[::-1] / len(g)
    print(f"{e}: stage I {d[0]:.1f}, II {d[1]:.1f}, IIIa {d[2]:.1f}, IIIb {d[3]:.1f}, IV/V {d[4]:.1f}")
print("published: EKFC stage II 52.21; CG stage II 37.58, stage I 30.19; CKD-EPI (calibration target) I 15.3 II 50.0 IIIa 20.3 IIIb 10.6 IV/V 3.8")
