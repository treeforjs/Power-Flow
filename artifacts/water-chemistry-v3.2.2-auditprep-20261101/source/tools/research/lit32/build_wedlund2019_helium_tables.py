"""Write recommended He+/He2+ + H2O charge-changing and ionization cross sections and polynomial fit
coefficients transcribed from Simon Wedlund et al., A&A 630, A35 (2019), arXiv:1902.05048 (A&A proofs
version): Tables A.1 (PDF page 17), B.1 (PDF page 19), 1 (PDF page 5) and 3 (PDF page 11).

Numbers were transcribed from the pdftotext -layout text layer and checked a second time against the
page images rendered at 250 dpi (A.1, B.1) and 300 dpi (Tables 1 and 3). The extrapolation asterisks were
read from the page images only (their placement in the text layer is ambiguous).

Usage: python3 -I build_wedlund2019_helium_tables.py <extracted_dir>
"""
import json
import math
import os
import sys

SRC = "ion-wedlund-2019-h-he-h2o-cx"
DOI = "10.1051/0004-6361/201834848"
CIT = ("C. Simon Wedlund, D. Bodewits, M. Alho, R. Hoekstra, E. Behar, G. Gronoff, H. Gunell, H. Nilsson, E. Kallio, "
       "A. Beth, Solar wind charge exchange in cometary atmospheres. I. Charge-changing and ionization cross sections "
       "for He and H particles in H2O, Astron. Astrophys. 630, A35 (2019); arXiv:1902.05048")
SPEED = [100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000,
         2500, 3000, 3500, 4000, 4500, 5000]
E_A1 = [0.052, 0.117, 0.209, 0.326, 0.470, 0.639, 0.835, 1.06, 1.30, 1.58, 1.88, 2.56, 3.34, 4.23, 5.22, 8.16,
        11.7, 16.0, 20.9, 32.6, 47.0, 63.9, 83.5, 105.7, 130.5]
E_B1 = [0.052, 0.117, 0.209, 0.326, 0.470, 0.639, 0.835, 1.057, 1.305, 1.579, 1.879, 2.558, 3.341, 4.228, 5.220,
        8.156, 11.740, 15.990, 20.880, 32.620, 46.980, 63.940, 83.520, 105.700, 130.500]

# Table A.1, Helium-H2O cross sections [m2]
A1 = {
    "s21": [1.62e-20, 2.74e-20, 3.92e-20, 5.12e-20, 6.28e-20, 7.39e-20, 8.43e-20, 9.39e-20, 1.03e-19, 1.10e-19,
            1.17e-19, 1.29e-19, 1.37e-19, 1.42e-19, 1.45e-19, 1.44e-19, 1.36e-19, 1.24e-19, 1.11e-19, 8.44e-20,
            6.22e-20, 4.50e-20, 3.23e-20, 2.31e-20, 1.66e-20],
    "s20": [1.03e-19, 8.46e-20, 6.65e-20, 5.40e-20, 4.58e-20, 4.02e-20, 3.65e-20, 3.39e-20, 3.21e-20, 3.08e-20,
            3.00e-20, 2.90e-20, 2.88e-20, 2.88e-20, 2.90e-20, 2.95e-20, 2.92e-20, 2.77e-20, 2.52e-20, 1.84e-20,
            1.15e-20, 6.31e-21, 3.11e-21, 1.40e-21, 5.85e-22],
    "s12": [2.80e-28, 2.34e-27, 1.06e-26, 3.47e-26, 9.12e-26, 2.07e-25, 4.19e-25, 7.81e-25, 1.36e-24, 2.24e-24,
            3.52e-24, 7.80e-24, 1.54e-23, 2.77e-23, 4.64e-23, 1.34e-22, 3.05e-22, 5.89e-22, 1.01e-21, 2.29e-21,
            4.12e-21, 6.31e-21, 8.65e-21, 1.09e-20, 1.28e-20],
    "s10": [1.08e-19, 7.91e-20, 6.37e-20, 5.66e-20, 5.36e-20, 5.29e-20, 5.33e-20, 5.45e-20, 5.61e-20, 5.80e-20,
            5.99e-20, 6.38e-20, 6.72e-20, 6.99e-20, 7.17e-20, 7.27e-20, 6.94e-20, 6.32e-20, 5.56e-20, 4.03e-20,
            2.76e-20, 1.84e-20, 1.22e-20, 8.05e-21, 5.34e-21],
}
A1_EXTRAP = {"s21": [], "s20": [], "s12": [100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 800],
             "s10": [100]}
# Table B.1, Helium-H2O columns [m2]
B1 = {
    "s22": [1.41e-21, 2.46e-21, 3.63e-21, 4.95e-21, 6.40e-21, 7.99e-21, 9.69e-21, 1.15e-20, 1.35e-20, 1.55e-20,
            1.76e-20, 2.22e-20, 2.70e-20, 3.20e-20, 3.72e-20, 5.08e-20, 6.44e-20, 7.76e-20, 9.00e-20, 1.11e-19,
            1.27e-19, 1.38e-19, 1.44e-19, 1.46e-19, 1.44e-19],
    "s11": [5.22e-21, 7.77e-21, 1.03e-20, 1.27e-20, 1.51e-20, 1.76e-20, 2.00e-20, 2.23e-20, 2.47e-20, 2.70e-20,
            2.93e-20, 3.38e-20, 3.81e-20, 4.22e-20, 4.61e-20, 5.50e-20, 6.25e-20, 6.88e-20, 7.38e-20, 8.07e-20,
            8.42e-20, 8.49e-20, 8.38e-20, 8.12e-20, 7.78e-20],
}
B1_EXTRAP = {"s22": [100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000, 1250],
             "s11": [100, 150, 200, 250, 300, 350, 400]}

# Fit coefficients as printed (Table 1 for charge changing, Table 3 for ionization)
FITS = {
    "s21": dict(table="1", page=5, n=4, p=["−129.6349", "85.3069", "−25.7100", "3.50927", "−0.18016"],
                vkm=[75, 5350], ekev=[0.03, 150], conf="high"),
    "s20": dict(table="1", page=5, n=5, p=["3327.3456", "−3277.1313", "1272.0445", "−244.7593", "23.35760", "−0.88492"],
                vkm=[100, 5350], ekev=[0.05, 150], conf="high"),
    "s12": dict(table="1", page=5, n=4, p=["−314.9414", "205.4565", "−57.4465", "7.3111", "−0.34819"],
                vkm=[820, 10000], ekev=[3.50, 520], conf="high"),
    "s10": dict(table="1", page=5, n=5, p=["−5450.8180", "4667.4695", "−1592.4968", "269.5820", "−22.63262", "0.75346"],
                vkm=[120, 5000], ekev=[0.08, 130], conf="high"),
    "s22": dict(table="3", page=11, n=4, p=["−265.3697", "176.7016", "−48.3734", "5.91427", "−0.27030"],
                vkm=[1400, 10000], ekev=[10.0, 520], conf="high"),
    "s11": dict(table="3", page=11, n=4, p=["−169.8607", "108.1412", "−29.7492", "3.66917", "−0.16967"],
                vkm=[450, 10000], ekev=[1.06, 520], conf="high"),
}

TEXT = {
    "s21": dict(
        slug="he2plus-h2o-electron-capture-sigma21", proj="he2plus",
        proc="He2+ + H2O -> He+ + [H2O]+ (one-electron capture by the projectile, all target product channels summed)",
        validity_text=("Sect. 3.1.1 (PDF page 4): 'A least-squares polynomial fit of degree 3 in log10 of the He2+ speed v_i "
                       "was performed. Expected validity range v_i = 75 - 5350 km/s (E_i = 0.03 - 150 keV/u). Confidence: high.' "
                       "Data selection: Seredyuk et al. (2005) 0.025-2 keV/u, Greenwood et al. (2004) 2-5.19 keV/u, Greenwood "
                       "et al. (2000) up to 6.67 keV/u, Rudd et al. (1985b) above 10 keV/u. Uncertainties ~10% at low energies "
                       "(15-25% below 0.3 keV/u, 95% confidence) and 12% at high energies (Rudd et al. 1985b)."),
        disc=("Discrepancy: the Sect. 3.1.1 text states a polynomial of degree 3, Table 1 prints degree n = 4 with five "
              "coefficients p0-p4. The Table 1 coefficients are recorded; the fit-vs-table check below shows that the "
              "five-coefficient polynomial reproduces Table A.1.")),
    "s20": dict(
        slug="he2plus-h2o-double-capture-sigma20", proj="he2plus",
        proc="He2+ + H2O -> He + [H2O]2+ (two-electron capture by the projectile, all target product channels summed)",
        validity_text=("Sect. 3.1.1 (PDF page 5): 'A polynomial fit of order 5 best represents the datasets. Expected validity "
                       "range: v_i = 100 - 5350 km/s (E_i = 0.05 - 150 keV/u). Confidence: high.' Data selection: Greenwood "
                       "et al. (2004) recommendation for 0.052-5.19 keV/u and Rudd et al. (1985b) above 5 keV/u. Uncertainties "
                       "~20% below 5 keV/u (30-40% below 0.3 keV/u, Greenwood et al. 2004) and 16% above (Rudd et al. 1985b)."),
        disc=""),
    "s12": dict(
        slug="heplus-h2o-electron-loss-sigma12", proj="heplus",
        proc="He+ + H2O -> He2+ + e + [H2O] (projectile one-electron loss / stripping)",
        validity_text=("Sect. 3.1.2 (PDF page 5): 'A polynomial fit in log10 of order 4 was used. Validity range: v_i = 820 - "
                       "10 000 km/s (E_i = 3.5 - 520 keV/u). Confidence: high.' Data: Rudd et al. (1985d) measurements in H2O, "
                       "3.50-112.5 keV/u (820-4640 km/s), uncertainties 21-33%; no measurements below or above. High-energy "
                       "peak defined by the additive rule sigma(H2) + sigma(O2)/2 with Sataka et al. (1990) H2 and O2 data at "
                       "75-450 keV/u, multiplied by a constant 1.64 to match Rudd et al. (1985d) at 75 keV/u."),
        disc=("Note: the source's 'high' confidence class (< 25%) is printed although the quoted measurement uncertainties "
              "are 21-33% and the part of the curve above 112.5 keV/u is an additive-rule reconstruction scaled by 1.64.")),
    "s10": dict(
        slug="heplus-h2o-electron-capture-sigma10", proj="heplus",
        proc="He+ + H2O -> He + [H2O]+ (one-electron capture by the projectile, all target product channels summed)",
        validity_text=("Sect. 3.1.2 (PDF page 6): 'Polynomial fit of order 5 was performed. Validity range: v_i = 120 - 5000 "
                       "km/s (E_i = 0.08 - 130 keV/u). Confidence: high.' Data selection: Greenwood et al. (2000) below "
                       "1.7 keV/u (uncertainty < 7%), Rudd et al. (1985d) for E_i >= 2.5 keV/u (14-20%), plus one Koopman "
                       "(1968) point at 73 eV/u (118 km/s) multiplied by a calibration factor 7.15 (sigma10 ~ 9e-20 m2) with "
                       "reduced weight. Koopman (1968) values are a factor 7.15 below Greenwood et al. (2000) at 0.33 keV/u."),
        disc=("Note: the low-energy rise below ~350 km/s rests on the single rescaled Koopman (1968) point; the source asks for "
              "measurements below 0.5 keV/u.")),
    "s22": dict(
        slug="he2plus-h2o-ionization-sigma22", proj="he2plus",
        proc="He2+ + H2O -> He2+ + [H2O]q+ + q e (target ionization without change of the projectile charge, sigma22; all dissociation paths of the target)",
        validity_text=("Sect. 4.1 (PDF page 10): 'A polynomial fit of order 4 of the cross section as a function of the "
                       "logarithm of the impact speed was performed. Expected validity range is v_i = 1400 - 10 000 km/s "
                       "(E_i = 10 - 520 keV/u). Confidence: high.' Data: Rudd et al. (1985b) total electron production "
                       "10-300 keV/u (1400-7500 km/s), uncertainty 13%; Toburen et al. (1980) 75-500 keV/u, uncertainty 20%."),
        disc=""),
    "s11": dict(
        slug="heplus-h2o-ionization-sigma11", proj="heplus",
        proc="He+ + H2O -> He+ + [H2O]q+ + q e (target ionization without change of the projectile charge, sigma11; all dissociation paths of the target)",
        validity_text=("Sect. 4.1 (PDF page 10): 'A polynomial fit of order 4 in log10 v_i was performed. Expected validity is "
                       "450 - 10 000 km/s (E_i = 1 - 520 keV/u). Confidence: high.' Data: Rudd et al. (1985d) 1.25-112.5 keV/u "
                       "(490-4650 km/s) and Toburen et al. (1980) 75-500 keV/u. Uncertainties 20% below 30 keV/u, < 8% at "
                       "30-450 keV/u (Rudd et al. 1985d), ~20% above 450 keV/u (Toburen et al. 1980)."),
        disc=("Note: Table 3 gives E_i = 1.06 - 520 keV/u; the Sect. 4.1 text gives 1 - 520 keV/u. The Table 3 values are "
              "recorded.")),
}

IONIZ_DEF = ("Measured quantities behind the ionization recommendations are total electron production or positive-ion "
             "production cross sections (Sect. 4, PDF page 9-10), which can contain transfer-ionization contributions.")
CONF_DEF = "high (< 25%), medium (25-75%), and low (> 75%) uncertainty (Table 1/3 Notes; Sect. 2.1)"
EXTRAP_RULES = ("Sect. 2.3 (PDF page 2-3): the validity range of each fit is confined to the velocity range of available "
                "measurements (or reconstructed datasets); where needed, smooth extrapolations were made in power laws of "
                "the velocity down to 100 km/s and for very high energies; these extrapolations have large uncertainties and "
                "are only given for reference in the appendix tables (A.1, B.1). The extrapolation exponents are not printed.")
AMU_KEV = 0.5 * 1.66053906660e-27 / 1.602176634e-16  # keV per (m/s)^2 per u


def parse(s):
    return float(s.replace("−", "-"))


def polyval(p, v_m_s):
    lv = math.log10(v_m_s)
    return 10.0 ** sum(c * lv ** j for j, c in enumerate(p))


def sig_round_equal(a, b):
    for n in (2, 3, 4):
        if float("%.*e" % (n - 1, b)) == a:
            return True
    return False


def main():
    out = sys.argv[1]
    for key in ("s21", "s20", "s12", "s10", "s22", "s11"):
        tx, fit = TEXT[key], FITS[key]
        is_a1 = key in A1
        y = A1[key] if is_a1 else B1[key]
        ext = A1_EXTRAP[key] if is_a1 else B1_EXTRAP[key]
        sym = "sigma" + key[1:]
        tab_name = "A.1" if is_a1 else "B.1"
        tab_page = 17 if is_a1 else 19
        conf_cls = "a (<25%)"
        mass_u = 4.002602
        t = {
            "table_id": "ion-wedlund2019-" + tx["slug"],
            "source_id": SRC, "citation": CIT, "doi": DOI,
            "locator": "Table %s (arXiv:1902.05048 PDF page %d), column s%s; fit coefficients in Table %s (PDF page %d)" % (
                tab_name, tab_page, key[1:], fit["table"], fit["page"]),
            "process": tx["proc"],
            "quantity": "cross_section_m2",
            "x_name": "energy_eV_per_amu_lab",
            "x": [round(e * 1000.0, 3) for e in (E_A1 if is_a1 else E_B1)],
            "x_impact_speed_km_s": SPEED,
            "y": y,
            "y_uncertainty": None,
            "y_uncertainty_rel": None,
            "confidence_class": conf_cls,
            "extrapolated_rows_speed_km_s": ext,
            "original_units": "m2; x printed as keV/amu with impact speed in km/s (speed is the primary variable of the source fits)",
            "conversion_to_SI": ("none for y (already m2, column header 'Helium-H2O cross sections [m2]'); x[eV/amu] = printed "
                                 "keV/amu * 1000. For a 4He projectile of mass 4.002602 u, E_lab[eV] = x * 4.002602 (not applied)."),
            "tag": "recommended",
            "extraction_method": ("transcribed from the text layer (pdftotext -layout) and checked against the page image "
                                  "rendered at 250 dpi; extrapolation asterisks read from the page image"),
            "notes": ("Recommended values are evaluated from the fitted polynomials log10(sigma) = sum p_j (log10 v)^j (v in m/s), "
                      "i.e. smoothed recommendations, not raw measurements. Source confidence classes (superscript in the "
                      "column header): a <25%, b 25-75%, c >75%; this column is class a. Rows marked with an asterisk in the "
                      "source ('smoothly extrapolated from the fits') are listed in extrapolated_rows_speed_km_s. Second-pass "
                      "check: all 25 values and the asterisk marks were read a second time from the page image and agree with "
                      "the text layer. Source text is the A&A proofs version on arXiv. Basis: " + tx["validity_text"] +
                      (" " + IONIZ_DEF if key in B1 else "") + (" " + tx["disc"] if tx["disc"] else "")),
        }
        with open(os.path.join(out, t["table_id"] + ".json"), "w") as f:
            json.dump(t, f, indent=1)
        print("wrote", t["table_id"], len(y))

        p = [parse(s) for s in fit["p"]]
        rows = []
        for v, ytab in zip(SPEED, y):
            fv = polyval(p, v * 1000.0)
            rows.append({"speed_km_s": v, "table_m2": ytab, "fit_m2": float("%.4g" % fv),
                         "rel_diff": round(fv / ytab - 1.0, 4),
                         "inside_validity": fit["vkm"][0] <= v <= fit["vkm"][1],
                         "table_equals_rounded_fit": sig_round_equal(ytab, fv)})
        ins = [abs(r["rel_diff"]) for r in rows if r["inside_validity"]]
        outs = [abs(r["rel_diff"]) for r in rows if not r["inside_validity"]]
        neq = sum(r["table_equals_rounded_fit"] for r in rows)
        check = ("Fit evaluated at the impact speeds of the recommended-value table %s (Table %s): %d rows inside the stated "
                 "validity range, max |fit/table - 1| = %.3g%%; %d rows outside, max |fit/table - 1| = %s; %d of 25 table "
                 "values equal the fit rounded to 2, 3 or 4 significant figures.") % (
            t["table_id"] + ".json", tab_name, len(ins), 100 * max(ins) if ins else float("nan"), len(outs),
            ("%.3g%%" % (100 * max(outs))) if outs else "n/a", neq)
        vk = [round(AMU_KEV * (fit["vkm"][0] * 1e3) ** 2, 4), round(AMU_KEV * (fit["vkm"][1] * 1e3) ** 2, 4)]
        c = {
            "table_id": "ion-wedlund2019-fitcoeffs-%s-h2o-%s" % (tx["proj"], sym),
            "source_id": SRC, "citation": CIT, "doi": DOI,
            "locator": ("Table %s '%s', arXiv:1902.05048 PDF page %d (A&A proofs article page %d of 20), row %s" % (
                fit["table"],
                ("Recommended charge-changing cross section polynomial fits for (He2+, He+, and He0) projectiles colliding "
                 "with H2O vapor") if fit["table"] == "1" else
                ("Recommended ionization cross-section polynomial fits for (He2+, He+, He0) projectiles colliding with H2O "
                 "vapor"), fit["page"], fit["page"], sym)),
            "process": tx["proc"],
            "quantity": "fit_coefficient_log10_sigma_m2_vs_log10_v_m_per_s",
            "x_name": "polynomial_order_j",
            "x": list(range(fit["n"] + 1)),
            "y": p,
            "y_uncertainty": None,
            "y_uncertainty_rel": None,
            "degree_n": fit["n"],
            "fit_equation": ("log10(sigma) = sum_{j=0}^{n} p_j (log10 v_i)^j, with v_i the impactor speed in m s^-1 and sigma "
                             "in m^2 (Eq. (3), arXiv PDF page 2, Sect. 2.1 'Fit and validity'; restated in the Notes of "
                             "Tables 1 and 3)"),
            "validity_speed_km_s": fit["vkm"],
            "validity_keV_per_u": fit["ekev"],
            "validity_text": tx["validity_text"],
            "confidence_level": fit["conf"],
            "confidence_definition": CONF_DEF,
            "extrapolation_rules": EXTRAP_RULES,
            "original_units": "dimensionless coefficients for sigma in m2 and v_i in m/s",
            "conversion_to_SI": "none",
            "tag": "recommended",
            "extraction_method": ("transcribed from the text layer (pdftotext -layout) and checked against the page image "
                                  "rendered at 300 dpi"),
            "fit_vs_recommended_table_check": check,
            "fit_vs_recommended_table_rows": rows,
            "notes": ("Coefficients p_j are dimensionless as printed (for sigma in m2 and v_i in m/s); no unit conversion "
                      "applied. Second-pass check: every coefficient digit, sign, the degree, the validity range and the "
                      "confidence entry were read a second time from the 300 dpi render of the table and agree with the text "
                      "layer. Sanity check: " + check + " Speed-energy consistency (E/u = 0.5 m_u v^2): printed speed and "
                      "keV/u limits agree within 10% (see validity_keV_per_u_recomputed_from_speed)." +
                      (" " + tx["disc"] if tx["disc"] else "")),
            "validity_keV_per_u_recomputed_from_speed": vk,
            "y_as_printed": fit["p"],
        }
        with open(os.path.join(out, c["table_id"] + ".json"), "w") as f:
            json.dump(c, f, indent=1)
        print("wrote", c["table_id"], check)


if __name__ == "__main__":
    main()
