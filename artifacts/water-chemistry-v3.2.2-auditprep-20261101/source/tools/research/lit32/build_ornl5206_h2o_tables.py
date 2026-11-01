"""Write the H2O-target tables of the 1977 ORNL "Red Book", C. F. Barnett et al., Atomic data for controlled
fusion research, ORNL-5206 Vol. 1 (1977), OSTI 7112109, DOI 10.2172/7112109.

Pages (report page / PDF page of the OSTI scan): A.3.26 / 137, A.4.16 / 159, A.5.38 / 235, A.6.18 / 291,
A.6.20 / 293. The scan has an OCR text layer with frequent character errors; every number below was read
from the page image rendered at 150 dpi, compared with the OCR text layer, and read a second time from the
image where the two disagreed.

Usage: python3 -I build_ornl5206_h2o_tables.py <extracted_dir>
"""
import json
import os
import sys

SRC = "ion-barnett-1977-ornl5206-vol1-redbook"
DOI = "10.2172/7112109"
CIT = ("C. F. Barnett, J. A. Ray, E. Ricci, M. I. Wilker, E. W. McDaniel, E. W. Thomas, H. B. Gilbody, Atomic data "
       "for controlled fusion research, Vol. 1, ORNL-5206 (Oak Ridge National Laboratory, February 1977); OSTI ID 7112109")
CM2 = 1e-4
EXM = ("transcribed from page image (150 dpi render of the OSTI scan) and compared with the OCR text layer; "
       "OCR character errors were resolved from the image")


def tab(table_id, report_page, pdf_page, title, process, quantity, x_name, x, y_src, unit_cm2, unc_rel, tag, notes,
        extra=None):
    y = [None if v is None else (float("%.6g" % (v * CM2)) if unit_cm2 else v) for v in y_src]
    pairs = [(a, b) for a, b in zip(x, y) if b is not None]
    t = {
        "table_id": table_id,
        "source_id": SRC, "citation": CIT, "doi": DOI,
        "locator": "Report page %s '%s' (PDF page %d of the OSTI scan)" % (report_page, title, pdf_page),
        "process": process,
        "quantity": quantity,
        "x_name": x_name,
        "x": [p[0] for p in pairs],
        "y": [p[1] for p in pairs],
        "y_as_printed": [v for v in y_src if v is not None],
        "y_uncertainty": unc_rel,
        "y_uncertainty_rel": unc_rel,
        "original_units": ("energy keV (laboratory energy of the projectile); cross section cm^2" if unit_cm2 else
                           "energy keV (laboratory energy of the projectile); yield in atoms per incident ion (dimensionless)"),
        "conversion_to_SI": ("cross section converted with 1 cm^2 = 1e-4 m^2; energy kept in keV" if unit_cm2 else
                             "none (dimensionless yield); energy kept in keV"),
        "tag": tag,
        "extraction_method": EXM,
        "notes": notes + (" Source units: cm^2, converted to m^2 with 1 cm^2 = 1e-4 m^2 (y_as_printed keeps the printed "
                          "cm^2 values)." if unit_cm2 else "") +
                 " Second-pass check: every number was read a second time from the page image and agrees with the "
                 "transcription.",
    }
    if extra:
        t.update(extra)
    return t


def main():
    out = sys.argv[1]
    T = []
    # A.3.26 H3+ dissociation in water vapor (Dagnac et al. 1971), +-30%
    E = [4.5, 6.0, 8.0, 10.0, 20.0, 40.0, 50.0, 60.0]
    title = "Cross Sections for the Formation of H+, H2+, and (H2^0+H^0) from the Dissociation of H3+ in Water Vapor"
    ref = ("Single reference: R. Dagnac, J. Angles, D. Blanc, D. Molina and R. Pradal, C. R. Acad. Sci. Paris B 273, 460 "
           "(1971). Stated accuracy +-30%. Energy is the H3+ laboratory energy (keV per ion, not per amu).")
    cols = {
        "neutral": ([None, 4.5e-16, 5.4e-16, 6.1e-16, 7.4e-16, 6.4e-16, 5.6e-16, None],
                    "H3+ + H2O -> (H^0 + H2^0) + ... (formation of fast neutral fragments, column sigma(H0+H2^0))"),
        "hplus": ([1.3e-16, 1.4e-16, 1.7e-16, 2.0e-16, 3.2e-16, 6.2e-16, 8.2e-16, None],
                  "H3+ + H2O -> H+ + ... (formation of fast H+ from dissociation of H3+, column sigma(H+))"),
        "h2plus": ([1.1e-16, 1.2e-16, 1.4e-16, 1.5e-16, 1.9e-16, 2.3e-16, 2.4e-16, 2.5e-16],
                   "H3+ + H2O -> H2+ + ... (formation of fast H2+ from dissociation of H3+, column sigma(H2+))"),
    }
    for k, (ys, proc) in cols.items():
        T.append(tab("ion-ornl5206-h3plus-h2o-dissociation-%s-formation" % k, "A.3.26", 137, title, proc,
                     "cross_section_m2", "energy_keV_lab", E, ys, True, 0.30, "measured",
                     ref + " Collision-induced dissociation of a molecular ion on H2O; the column header for the neutral "
                     "channel prints sigma(H0+H2^0) and the title prints (H2^0+H^0)."))
    # A.4.16 H+ + H2O electron capture
    E = [0.05, 0.07, 0.10, 0.20, 0.50, 0.75, 1.0, 2.0, 5.0, 7.5, 10.0, 20.0, 50.0, 75.0, 100.0, 200.0, 500.0, 750.0,
         1000.0, 2500.0]
    ys = [3.0e-15, 2.9e-15, 2.8e-15, 2.6e-15, 2.3e-15, 2.2e-15, 2.1e-15, 1.8e-15, 1.4e-15, 1.2e-15, 1.1e-15, 7.6e-16,
          3.5e-16, 1.8e-16, 9.0e-17, 1.1e-17, 4.3e-19, 9.2e-20, 3.4e-20, 1.6e-21]
    T.append(tab("ion-ornl5206-hplus-h2o-electron-capture-sigma10", "A.4.16", 159,
                 "Electron Capture Cross Sections for H+ in H2O and CH4",
                 "H+ + H2O -> H^0 + [H2O]+ (one-electron capture, sigma10)", "cross_section_m2", "energy_keV_lab", E, ys,
                 True, 0.25, "recommended",
                 "Compiled recommended curve from K. H. Berkner, R. V. Pyle, J. W. Stearns, Nucl. Fusion 10, 145 (1970); "
                 "E. S. Chambers, UCRL-14214 (1965); M. A. Coplan, K. W. Ogilvie, J. Chem. Phys. 52, 4154 (1970); "
                 "R. Dagnac, D. Blanc, D. Molina, J. Phys. B 3, 1239 (1970); D. W. Koopman, Phys. Rev. 166, 57 (1968); "
                 "L. H. Toburen, M. Y. Nakai, R. A. Langley, Phys. Rev. 171, 114 (1968). Stated accuracy +-25%. "
                 "Independent cross-check of Simon Wedlund et al. 2019 sigma10 (ion-wedlund2019-hplus-h2o-electron-capture-"
                 "sigma10.json). The neighbouring CH4 column was not transcribed."))
    # A.5.38 H + H2O electron loss
    E = [100.0, 200.0, 500.0, 1000.0, 1500.0, 2000.0]
    ys = [3.18e-16, 2.70e-16, 1.55e-16, 9.81e-17, 7.10e-17, 5.70e-17]
    T.append(tab("ion-ornl5206-h-h2o-electron-loss-sigma01", "A.5.38", 235,
                 "Cross Sections for Electron Loss of H Atoms in N2, O2, and H2O",
                 "H + H2O -> H+ + e + [H2O] (projectile electron loss, sigma01)", "cross_section_m2", "energy_keV_lab", E,
                 ys, True, 0.25, "measured",
                 "H2O column from a single reference: L. H. Toburen, M. Y. Nakai, R. A. Langley, Phys. Rev. 171, 114 "
                 "(1968). Stated accuracy +-25%. The N2 and O2 columns on the same page were not transcribed."))
    # A.6.18 / A.6.20 yields
    ref63 = ("H0 yield: C. F. Barnett, M. Rankin, J. A. Ray, Proc. 6th Intl. Conf. on Ionization Phenomena in Gases "
             "(SERMA, Paris 1963) Vol. 1 p. 63; H- yield: N. B. Brooks, M. Roos, R. P. Bastide, P. H. Rose, "
             "A. B. Wittkower, Phys. Lett. 6, 169 (1963). Accuracy: random error < +-20%. These are yields (atoms per "
             "incident ion) observed after the beam traversed the gas, not cross sections; the target thickness is "
             "not given on the page.")
    E = [10.0, 15.0, 20.0, 30.0]
    T.append(tab("ion-ornl5206-h2plus-h2o-hminus-yield", "A.6.18", 291,
                 "Yield of H0 and H- Observed for an H2+ Beam Traversing H2 and H2O",
                 "H2+ + H2O -> H- + ... (H- atoms per incident H2+ ion)", "yield_atoms_per_incident_ion",
                 "energy_keV_lab", E, [6.15e-2, 8.17e-2, 7.78e-2, 6.27e-2], False, 0.20, "measured", ref63))
    E = [60.0, 80.0, 100.0, 150.0, 200.0, 300.0, 400.0]
    T.append(tab("ion-ornl5206-h2plus-h2o-h0-yield", "A.6.18", 291,
                 "Yield of H0 and H- Observed for an H2+ Beam Traversing H2 and H2O",
                 "H2+ + H2O -> H0 + ... (fast H0 atoms per incident H2+ ion)", "yield_atoms_per_incident_ion",
                 "energy_keV_lab", E, [1.16, 0.940, 0.782, 0.547, 0.398, 0.235, 0.151], False, 0.20, "measured",
                 ref63 + " The 40 keV row is blank for all columns."))
    E = [10.0, 15.0, 20.0, 30.0]
    T.append(tab("ion-ornl5206-h3plus-h2o-hminus-yield", "A.6.20", 293,
                 "Yield of H0 and H- Observed for an H3+ Beam Traversing H2 and H2O",
                 "H3+ + H2O -> H- + ... (H- atoms per incident H3+ ion)", "yield_atoms_per_incident_ion",
                 "energy_keV_lab", E, [5.50e-2, 1.36e-1, 1.66e-1, 9.30e-1], False, 0.20, "measured",
                 ref63 + " PRINTED VALUE SUSPECT: the 30 keV entry is printed 9.30 E-01, a factor ~5.6 above the "
                 "20 keV entry; a misprint for 9.30 E-02 is likely but the printed value is kept."))
    E = [40.0, 60.0, 80.0, 100.0, 150.0, 200.0, 300.0, 400.0, 500.0]
    T.append(tab("ion-ornl5206-h3plus-h2o-h0-yield", "A.6.20", 293,
                 "Yield of H0 and H- Observed for an H3+ Beam Traversing H2 and H2O",
                 "H3+ + H2O -> H0 + ... (fast H0 atoms per incident H3+ ion)", "yield_atoms_per_incident_ion",
                 "energy_keV_lab", E, [2.53, 2.24, 1.98, 1.73, 1.30, 1.01, 6.92, 5.05, 3.90], False, 0.20, "measured",
                 ref63 + " PRINTED VALUES INCONSISTENT: the 300, 400 and 500 keV entries are printed with exponent E 00 "
                 "(6.92, 5.05, 3.90), which exceeds the maximum of 3 H atoms per H3+ ion and breaks the monotonic "
                 "decrease; exponent E-01 (0.692, 0.505, 0.390) is likely intended. The printed values are kept here; "
                 "do not use those three rows without correction."))
    for t in T:
        with open(os.path.join(out, t["table_id"] + ".json"), "w") as f:
            json.dump(t, f, indent=1)
        print("wrote", t["table_id"], len(t["y"]))


if __name__ == "__main__":
    main()
