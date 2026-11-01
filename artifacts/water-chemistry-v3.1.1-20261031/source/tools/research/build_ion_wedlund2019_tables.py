"""Write recommended H+/H/H- + H2O charge-changing and ionization cross sections transcribed from
Simon Wedlund et al., A&A 630, A35 (2019), arXiv:1902.05048 (A&A proofs version), Tables A.2 and B.1.
Values were transcribed from the rendered table images and agree with the PDF text layer.

Usage: python3 -I build_ion_wedlund2019_tables.py <extracted_dir>
"""
import json
import os
import sys

SRC = "ion-wedlund-2019-h-he-h2o-cx"
CIT = ("C. Simon Wedlund, D. Bodewits, M. Alho, R. Hoekstra, E. Behar, G. Gronoff, H. Gunell, H. Nilsson, E. Kallio, "
       "A. Beth, Solar wind charge exchange in cometary atmospheres. I. Charge-changing and ionization cross sections "
       "for He and H particles in H2O, Astron. Astrophys. 630, A35 (2019); arXiv:1902.05048")
SPEED = [100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000,
         2500, 3000, 3500, 4000, 4500, 5000]
E_A2 = [0.052, 0.117, 0.209, 0.326, 0.470, 0.639, 0.835, 1.06, 1.30, 1.58, 1.88, 2.56, 3.34, 4.23, 5.22, 8.16,
        11.7, 16.0, 20.9, 32.6, 47.0, 63.9, 83.5, 105.7, 130.5]
E_B1 = [0.052, 0.117, 0.209, 0.326, 0.470, 0.639, 0.835, 1.057, 1.305, 1.579, 1.879, 2.558, 3.341, 4.228, 5.220,
        8.156, 11.740, 15.990, 20.880, 32.620, 46.980, 63.940, 83.520, 105.700, 130.500]

A2 = {
    "s10": [1.34e-19, 2.03e-19, 2.24e-19, 2.19e-19, 2.06e-19, 1.92e-19, 1.79e-19, 1.68e-19, 1.59e-19, 1.52e-19,
            1.46e-19, 1.37e-19, 1.30e-19, 1.25e-19, 1.21e-19, 1.11e-19, 1.02e-19, 9.06e-20, 7.91e-20, 5.64e-20,
            3.76e-20, 2.39e-20, 1.47e-20, 8.84e-21, 5.27e-21],
    "s1m1": [4.8e-25, 1.9e-24, 4.3e-24, 8.2e-24, 1.4e-23, 2.3e-23, 3.5e-23, 5.1e-23, 7.2e-23, 9.9e-23, 1.3e-22,
             2.2e-22, 3.3e-22, 4.5e-22, 5.9e-22, 9.2e-22, 1.1e-21, 1.1e-21, 9.8e-22, 5.7e-22, 2.6e-22, 1.1e-22,
             4.0e-23, 1.5e-23, 5.3e-24],
    "s01": [5.17e-23, 1.63e-22, 3.80e-22, 7.26e-22, 1.22e-21, 1.86e-21, 2.64e-21, 3.56e-21, 4.60e-21, 5.74e-21,
            6.97e-21, 9.63e-21, 1.24e-20, 1.53e-20, 1.82e-20, 2.47e-20, 3.02e-20, 3.44e-20, 3.75e-20, 4.08e-20,
            4.15e-20, 4.05e-20, 3.87e-20, 3.64e-20, 3.39e-20],
    "s0m1": [9.4e-22, 1.8e-21, 2.7e-21, 3.6e-21, 4.5e-21, 5.3e-21, 6.0e-21, 6.6e-21, 7.1e-21, 7.4e-21, 7.7e-21,
             8.1e-21, 8.2e-21, 8.0e-21, 7.8e-21, 6.8e-21, 5.7e-21, 4.7e-21, 3.8e-21, 2.5e-21, 1.6e-21, 1.1e-21,
             7.4e-22, 5.2e-22, 3.7e-22],
    "sm11": [4.6e-22, 7.5e-22, 1.1e-21, 1.4e-21, 1.8e-21, 2.2e-21, 2.7e-21, 3.1e-21, 3.5e-21, 3.9e-21, 4.4e-21,
             5.2e-21, 6.0e-21, 6.8e-21, 7.5e-21, 9.1e-21, 1.0e-20, 1.1e-20, 1.1e-20, 1.2e-20, 1.1e-20, 1.0e-20,
             9.3e-21, 8.2e-21, 7.1e-21],
    "sm10": [5.7e-20, 6.8e-20, 7.4e-20, 7.8e-20, 8.0e-20, 8.1e-20, 8.2e-20, 8.2e-20, 8.1e-20, 8.1e-20, 8.0e-20,
             7.8e-20, 7.7e-20, 7.5e-20, 7.3e-20, 6.9e-20, 6.5e-20, 6.1e-20, 5.8e-20, 5.3e-20, 4.9e-20, 4.5e-20,
             4.1e-20, 3.8e-20, 3.5e-20],
}
A2_EXTRAP = {"s01": [100], "sm11": [100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600]}
B1 = {
    "s11H": [3.35e-23, 1.05e-22, 2.34e-22, 4.35e-22, 7.21e-22, 1.10e-21, 1.58e-21, 2.16e-21, 2.85e-21, 3.64e-21,
             4.53e-21, 6.59e-21, 8.96e-21, 1.16e-20, 1.44e-20, 2.19e-20, 2.93e-20, 3.60e-20, 4.18e-20, 4.99e-20,
             5.40e-20, 5.50e-20, 5.39e-20, 5.16e-20, 4.85e-20],
    "s00H": [4.6e-22, 1.5e-21, 3.2e-21, 5.4e-21, 7.8e-21, 1.1e-20, 1.3e-20, 1.6e-20, 1.9e-20, 2.2e-20, 2.5e-20,
             3.0e-20, 3.5e-20, 3.9e-20, 4.2e-20, 4.9e-20, 5.4e-20, 5.6e-20, 5.7e-20, 5.5e-20, 5.1e-20, 4.6e-20,
             4.1e-20, 3.6e-20, 3.2e-20],
}
B1_EXTRAP = {"s11H": [100, 150, 200, 250, 300, 350], "s00H": [100]}

COLS = [
    ("s10", "A2", "hplus-h2o-electron-capture-sigma10",
     "H+ + H2O -> H + [H2O]+ (one-electron capture by the projectile, all target product channels summed)",
     "recommended", "a (<25%)",
     "Basis (Sec. 3.2.1): measurements of Lindsay et al. 1997 (0.5-1.5 keV) and Greenwood et al. 2000 (1.5-7 keV) "
     "and higher-energy data; BELOW 500 eV the curve is the theoretical total charge-transfer cross section of Mada et "
     "al. 2007 multiplied by 2.2 to match the 500 eV data, so the 0.052-0.47 keV/u rows are theory-shaped, not measured. "
     "The source cites Baribaud et al. sigma10 = (14 +/- 8)e-20 m2 at 5 keV and notes Koopman 1968 low by ~10x."),
    ("s1m1", "A2", "hplus-h2o-double-capture-sigma1m1",
     "H+ + H2O -> H- + [H2O]2+ (two-electron capture)", "model_derived", "c (>75%)",
     "Basis: additive rule sigma(H2)+sigma(O2)/2 at low energy; only H2O measurement is Toburen & Nakai 1969 (75-250 keV)."),
    ("s01", "A2", "h-h2o-electron-loss-sigma01",
     "H + H2O -> H+ + e + [H2O] (projectile one-electron loss / stripping)", "model_derived", "b (25-75%)",
     "Basis: Dagnac et al. 1970 and Luna et al. 2007 (1.5-90 keV); below 1.5 keV the additive rule (H on H2 and on O) "
     "scaled down by the source (factor 3.8 noted at 2 keV). Value at 100 km/s is marked as extrapolated."),
    ("s0m1", "A2", "h-h2o-electron-capture-sigma0m1",
     "H + H2O -> H- + [H2O]+ (one-electron capture by neutral H)", "model_derived", "c (>75%) in Table A.2 footnote; text: medium overall, low below 8 keV",
     "Basis: Luna et al. 2007 (8-100 keV); below 8 keV additive rule data multiplied by a constant 4.8."),
    ("sm11", "A2", "hminus-h2o-two-electron-loss-sigmam11",
     "H- + H2O -> H+ + 2e + [H2O] (two-electron loss)", "model_derived", "c (>75%)",
     "Basis: single H2O point (Baribaud et al. 1971, 0.7e-20 m2 at 5 keV) and additive rule; rows <= 600 km/s are marked "
     "as extrapolated from the fit."),
    ("sm10", "A2", "hminus-h2o-one-electron-loss-sigmam10",
     "H- + H2O -> H + e + [H2O] (one-electron loss, collisional detachment of H-)", "model_derived", "c (>75%)",
     "Basis (Sec. 3.2.3): additive rule sigma(H2)+sigma(O2)/2 scaled by 1/2.3 to the single H2O measurement of "
     "Baribaud et al. 1971 (7.5e-20 m2 at 5 keV, 13% error). The source states the extrapolation to other energies has "
     "uncertainties that cannot be assessed. Relevant to channel 3 (H- detachment in H2O) above ~50 eV only."),
    ("s11H", "B1", "hplus-h2o-ionization-sigma11",
     "H+ + H2O -> H+ + [H2O]+ + e (target ionization without charge change)", "recommended", "a (<25%)",
     "Rows <= 350 km/s (<= 0.639 keV/u) are marked as extrapolated. PRINTED UNIT INCONSISTENCY retained: the Hydrogen-H2O "
     "column header reads [x1e-22 m2] while the entries are printed in E-notation (e.g. 4.85E-20); the entries are "
     "consistent with m2 (cf. the helium columns headed [m2] and the Table 4 fit, which returns m2) and are stored as m2."),
    ("s00H", "B1", "h-h2o-ionization-sigma00",
     "H + H2O -> H + [H2O]+ + e (target ionization by fast neutral H)", "model_derived", "c (>75%)",
     "Value at 100 km/s marked as extrapolated. Same printed header inconsistency as sigma11 (see that table); stored as m2."),
]


def main():
    out = sys.argv[1]
    for key, tab, slug, proc, tag, conf, basis in COLS:
        y = A2[key] if tab == "A2" else B1[key]
        ext = (A2_EXTRAP if tab == "A2" else B1_EXTRAP).get(key, [])
        t = {
            "table_id": "ion-wedlund2019-" + slug,
            "source_id": SRC, "citation": CIT, "doi": "10.1051/0004-6361/201834848",
            "locator": ("Table A.2 (arXiv:1902.05048 PDF page 18)" if tab == "A2" else "Table B.1 (arXiv:1902.05048 PDF page 19)")
                       + ", column " + key.replace("H", "") + "; fit coefficients in Table " + ("2" if tab == "A2" else "4"),
            "process": proc,
            "quantity": "cross_section_m2",
            "x_name": "energy_eV_per_amu_lab",
            "x": [round(e * 1000.0, 3) for e in (E_A2 if tab == "A2" else E_B1)],
            "x_impact_speed_km_s": SPEED,
            "y": y,
            "y_uncertainty": None,
            "confidence_class": conf,
            "extrapolated_rows_speed_km_s": ext,
            "original_units": "m2; x printed as keV/amu with impact speed in km/s (speed is the primary variable of the source fits)",
            "conversion_to_SI": "none for y (already m2); x[eV/amu] = printed keV/amu * 1000. For an H projectile of mass m_u = 1.00794 u, E_lab[eV] = x * m_u (not applied).",
            "tag": tag,
            "extraction_method": "manual transcription of the printed recommended-value table (rendered at 200 dpi, checked against PDF text layer)",
            "notes": ("Recommended values are evaluated from the fitted polynomials log10(sigma) = sum p_j (log10 v)^j (v in m/s), "
                      "i.e. smoothed recommendations, not raw measurements. Source confidence classes: a <25%, b 25-75%, c >75%. "
                      "Source text is the A&A proofs version on arXiv. " + basis),
        }
        with open(os.path.join(out, t["table_id"] + ".json"), "w") as f:
            json.dump(t, f, indent=1)
        print("wrote", t["table_id"], len(y))


if __name__ == "__main__":
    main()
