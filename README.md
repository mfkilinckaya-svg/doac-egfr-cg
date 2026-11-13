# DOAC dosing by eGFR versus Cockcroft–Gault
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24191222.svg)](https://doi.org/10.5281/zenodo.24191222)


Code, results and web calculator for the manuscript

> Kilinckaya MF. *Why DOAC dosing by eGFR departs from the label: a simulation of Cockcroft–Gault, CKD-EPI and EKFC under EU and US labels, checked against published cohorts.* (submitted)

The analysis asks, for apixaban, dabigatran, edoxaban and rivaroxaban, how often and in which direction the dose changes when the laboratory's eGFR (CKD-EPI 2009, CKD-EPI 2021 or EKFC; indexed or de-indexed to mL/min) is used instead of the Cockcroft–Gault (CG) creatinine clearance on which the EU (SmPC) and US (USPI) labels are based. No patient data are used.

## Contents

| File | What it does | Manuscript output |
| --- | --- | --- |
| `egfr.py` | CG, CKD-EPI 2009 (without race), CKD-EPI 2021, EKFC; DuBois and Mosteller BSA; de-indexing | all |
| `doac_rules.py` | EU and US label dose rules for the four DOACs (atrial fibrillation) | all |
| `sim.py` | Full-factorial grid (age 18–100 y, weight 40–150 kg, height 145–200 cm, sex, Scr 40–400 µmol/L) | Fig. 1, Tables S2, S3 |
| `mc_boriani.py` | Monte Carlo analysis of the European atrial fibrillation scenario (300 draws × 200,000 patients) | Tables 2 and 4, Table S4 |
| `contsim.py` | Continuous-sampling populations matched to the published cohorts; comparison with published discordance; estimator table; age–weight correlation sensitivity | Fig. 2 data, Table 3, Table S5 |
| `contsim_extra.py` | Additional numbers quoted in the text (CG stage check, interacting drugs, Mosteller vs DuBois) | Results text |
| `figure_validation.py` | Draws Fig. 2 (Bland–Altman mean difference, limits of agreement, Lin's CCC) | Fig. 2 |
| `calculator/doac_core.js` | Independent JavaScript implementation of all equations and dose rules (v1.0.1) | web calculator |
| `calculator/template.html`, `build_calculator.py`, `index.html` | Single-file web calculator (`index.html` is built from the template) | web calculator |
| `calculator/test_reference_values.py` | 35 pytest checks against hand-calculated values and every label threshold, in both implementations | Table S6 |
| `calculator/cross_validate.py` | Python vs JavaScript in 50,507 virtual patients, including cases exactly on every threshold | Table S6 |
| `results/` | All outputs used in the manuscript (CSV and logs) | — |
| `figures/` | Fig. 1 and Fig. 2 (PNG and vector PDF) | — |

`results/validation_points.csv` and `results/validation_points_cheng.csv` hold the published values from the comparison cohorts (column `pub`); the column `sim` there is from an earlier grid-based version and is not used. The values used in the manuscript are in `results/cont_validation_points.csv` and `results/cont_validation_cheng.csv` (column `sim_cont`).

## Reproduce

Python 3.11 or later and Node.js 18 or later (for the JavaScript tests).

```bash
pip install -r requirements.txt
python sim.py                 # grid: Fig. 1, Tables S2–S3
python mc_boriani.py          # Monte Carlo: Tables 2, 4, S4
python contsim.py             # cohort comparisons: Table 3, S5, Fig. 2 data
python contsim_extra.py       # numbers quoted in the Results
python figure_validation.py   # Fig. 2
python -m pytest calculator/test_reference_values.py
python calculator/cross_validate.py
python calculator/build_calculator.py
```

Random seeds are fixed, so the outputs in `results/` are reproduced exactly.

## Web calculator

Online: **https://mfkilinckaya-svg.github.io/doac-egfr-cg/calculator/** — or open `calculator/index.html` in any browser (no installation, no data leave the browser). It shows, for one patient, CG and each eGFR estimate and the resulting dose under the EU and US labels. It is a research tool, not a medical device: dosing follows the product label and clinical judgement.

## Conventions

- Scr in µmol/L (1 mg/dL = 88.4 µmol/L); CG with actual body weight.
- Renal function is not rounded before thresholds are applied; label ranges are coded as written (for example rivaroxaban SmPC 15–49 mL/min as ≥ 15 and < 50).
- Non-indexed eGFR = eGFR × BSA / 1.73 (DuBois in the main analysis, Mosteller as sensitivity).

## Citation

See `CITATION.cff`. Archived release (v1.0.0): https://doi.org/10.5281/zenodo.24191222

## Licence

MIT (see `LICENSE`).
