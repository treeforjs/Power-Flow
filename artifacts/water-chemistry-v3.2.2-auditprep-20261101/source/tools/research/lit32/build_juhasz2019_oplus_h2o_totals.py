"""Write the two total fragment-ion production cross sections for 6.6-keV 16O+ + H2O given in the text of
Z. Juhasz et al., Phys. Rev. A 100, 032713 (2019), arXiv:1910.00657v2, Sect. III.D 'Total cross sections'
(arXiv PDF page 8). Transcribed from the text layer and checked against the page image rendered at 200 dpi.

Usage: python3 -I build_juhasz2019_oplus_h2o_totals.py <extracted_dir>
"""
import json
import os
import sys

SRC = "ion-juhasz-2019-oplus-h2o-fragment-emission"
DOI = "10.1103/PhysRevA.100.032713"
CIT = ("Z. Juhasz, B. Sulik, E. Lattouf, E. Bene, B. A. Huber, P. Herczku, S. T. S. Kovacs, A. Mery, J.-C. Poully, "
       "J. Rangama, J. A. Tanis, V. Vizcaino, J.-Y. Chesnel, Anion and cation emission from water molecules after "
       "collisions with 6.6-keV 16O+ ions, Phys. Rev. A 100, 032713 (2019); arXiv:1910.00657")
COMMON = ("Single energy, 6.6 keV 16O+ (0.4125 keV/u). Obtained by the authors by integrating the measured absolute "
          "double-differential fragment-ion cross sections over energy and solid angle; the authors note the evaluation "
          "requires the contribution of the slowest fragments, which is hard to measure, and that the observed angular "
          "anisotropy makes the 40-140 deg range carry ~90% of the total. About two-thirds of the detected fragments are "
          "hydrogen ions (TOF); the authors estimate ~4e-16 cm^2 for H+ and ~3e-18 cm^2 for H- emission (same page). "
          "These are fragment-ion production cross sections of the target (dissociative ionization / capture and "
          "related processes summed); they are not electron-capture or projectile-loss cross sections. Source units: "
          "cm^2, converted with 1 cm^2 = 1e-4 m^2. Second-pass check: value and uncertainty read a second time from the "
          "page image and agree with the text layer.")


def main():
    out = sys.argv[1]
    rows = [
        ("ion-juhasz2019-oplus-h2o-cation-fragment-total-6p6kev",
         "16O+ + H2O -> positive fragment ions (H+, O+, OH+, ... summed) + ...", 6e-20, 3e-20, "(6 +- 3) x 10^-16 cm2"),
        ("ion-juhasz2019-oplus-h2o-anion-fragment-total-6p6kev",
         "16O+ + H2O -> negative fragment ions (H-, O-, OH- summed) + ...", 5e-22, 3e-22, "(5 +- 3) x 10^-18 cm2"),
    ]
    for tid, proc, y, dy, printed in rows:
        t = {
            "table_id": tid, "source_id": SRC, "citation": CIT, "doi": DOI,
            "locator": "Sect. III.D 'Total cross sections', arXiv:1910.00657v2 PDF page 8 (left column, text)",
            "process": proc,
            "quantity": "cross_section_m2",
            "x_name": "energy_keV_lab",
            "x": [6.6],
            "y": [y],
            "y_as_printed": printed,
            "y_uncertainty": [round(dy / y, 3)],
            "y_uncertainty_rel": [round(dy / y, 3)],
            "y_uncertainty_abs_m2": [dy],
            "original_units": "keV (lab, total projectile energy); cm^2",
            "conversion_to_SI": "1 cm^2 = 1e-4 m^2; energy kept in keV",
            "tag": "measured",
            "extraction_method": "transcribed from the text layer",
            "notes": COMMON,
        }
        with open(os.path.join(out, tid + ".json"), "w") as f:
            json.dump(t, f, indent=1)
        print("wrote", tid)


if __name__ == "__main__":
    main()
