"""Renal function equations for the DOAC dosing simulation.

All functions are vectorised (numpy). Creatinine input in µmol/L.
Conversion: 1 mg/dL = 88.4 µmol/L.

NOTE: CKD-EPI 2009 / 2021 / EKFC coefficients must be cross-checked against
Table 1 of the EKFC manuscript (egfr.py there) before submission.
"""
import numpy as np

UMOL_PER_MGDL = 88.4


def cockcroft_gault(age, weight_kg, scr_umol, female):
    """Cockcroft-Gault CrCl (mL/min), actual body weight, not BSA-indexed.
    CrCl = (140 - age) * weight / (72 * Scr[mg/dL]) * 0.85 if female
    (Cockcroft & Gault, Nephron 1976; formula as printed in EU SmPCs of
    edoxaban and dabigatran)."""
    scr_mg = np.asarray(scr_umol) / UMOL_PER_MGDL
    crcl = (140 - np.asarray(age)) * np.asarray(weight_kg) / (72 * scr_mg)
    return np.where(female, crcl * 0.85, crcl)


def ckd_epi_2009(age, scr_umol, female):
    """CKD-EPI 2009 creatinine (Levey 2009), race coefficient omitted
    (non-Black). mL/min/1.73 m²."""
    scr = np.asarray(scr_umol) / UMOL_PER_MGDL
    female = np.asarray(female, dtype=bool)
    kappa = np.where(female, 0.7, 0.9)
    alpha = np.where(female, -0.329, -0.411)
    r = scr / kappa
    egfr = 141 * np.minimum(r, 1) ** alpha * np.maximum(r, 1) ** -1.209 \
        * 0.993 ** np.asarray(age)
    return np.where(female, egfr * 1.018, egfr)


def ckd_epi_2021(age, scr_umol, female):
    """CKD-EPI 2021 race-free creatinine (Inker 2021). mL/min/1.73 m²."""
    scr = np.asarray(scr_umol) / UMOL_PER_MGDL
    female = np.asarray(female, dtype=bool)
    kappa = np.where(female, 0.7, 0.9)
    alpha = np.where(female, -0.241, -0.302)
    r = scr / kappa
    egfr = 142 * np.minimum(r, 1) ** alpha * np.maximum(r, 1) ** -1.200 \
        * 0.9938 ** np.asarray(age)
    return np.where(female, egfr * 1.012, egfr)


def ekfc_q(age, female):
    """EKFC Q value (µmol/L): age-dependent below 25 y, fixed 80 (M) /
    62 (F) from 25 y (Pottel 2021)."""
    age = np.asarray(age, dtype=float)
    female = np.asarray(female, dtype=bool)
    a = np.clip(age, 2, 25)
    lnq_m = 3.200 + 0.259 * a - 0.543 * np.log(a) - 0.00763 * a**2 + 0.0000790 * a**3
    lnq_f = 3.080 + 0.177 * a - 0.223 * np.log(a) - 0.00596 * a**2 + 0.0000686 * a**3
    q_young = np.where(female, np.exp(lnq_f), np.exp(lnq_m))
    q_adult = np.where(female, 62.0, 80.0)
    return np.where(age <= 25, q_young, q_adult)  # EKFC makalesi kodu (egfr_kfre.py) ile uyumlu


def ekfc(age, scr_umol, female):
    """EKFC creatinine equation (Pottel, Ann Intern Med 2021).
    eGFR = 107.3 * (Scr/Q)^-0.322 (Scr/Q<1) or ^-1.132 (>=1),
    * 0.990^(age-40) if age > 40. mL/min/1.73 m²."""
    age = np.asarray(age, dtype=float)
    r = np.asarray(scr_umol) / ekfc_q(age, female)
    egfr = 107.3 * np.where(r < 1, r ** -0.322, r ** -1.132)
    return np.where(age > 40, egfr * 0.990 ** (age - 40), egfr)


def bsa_dubois(weight_kg, height_cm):
    """DuBois & DuBois 1916: 0.007184 * W^0.425 * H^0.725 (m²)."""
    return 0.007184 * np.asarray(weight_kg) ** 0.425 * np.asarray(height_cm) ** 0.725


def bsa_mosteller(weight_kg, height_cm):
    """Mosteller 1987: sqrt(W*H/3600) (m²)."""
    return np.sqrt(np.asarray(weight_kg) * np.asarray(height_cm) / 3600)


def deindex(egfr_indexed, bsa):
    """Indexed eGFR (mL/min/1.73 m²) -> absolute eGFR (mL/min)."""
    return egfr_indexed * bsa / 1.73
