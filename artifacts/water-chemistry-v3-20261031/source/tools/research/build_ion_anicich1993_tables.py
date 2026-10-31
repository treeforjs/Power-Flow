"""Write recommended thermal (300 K) rate coefficients and product distributions transcribed from
Anicich, J. Phys. Chem. Ref. Data 22, 1469 (1993), Table of Reactions (NIST reprint PDF).
Values were read from the rendered page images (the PDF text layer is OCR and contains errors).

Usage: python3 -I build_ion_anicich1993_tables.py <extracted_dir>
"""
import json
import os
import sys

SRC = "ion-anicich-1993-jpcrd-evaluated"
CIT = ("V. G. Anicich, Evaluated bimolecular ion-molecule gas phase kinetics of positive ions for use in modeling "
       "planetary atmospheres, cometary comae, and interstellar clouds, J. Phys. Chem. Ref. Data 22, 1469-1569 (1993)")
BASE_NOTES = ("Evaluated recommendation; the source states that all data are referenced to 300 K unless a footnote "
              "says otherwise (no temperature footnote on this entry). Neutral products are inferred by the evaluator "
              "from atom balance, not observed (source Sec. 2). Uncertainty is the evaluator's stated relative "
              "uncertainty. Thermal value only: it is not an energy-dependent cross section and does not apply to "
              "suprathermal ions without separate data.")

REACTIONS = [
    # (key, page, pdf_page, reactants, channels [(products, fraction)], k_cm3_s, rel_unc, refs, extra_note)
    ("h2oplus-h2o", 1512, 45, "H2O+ + H2O", [("H3O+ + OH", 1.00)], 1.85e-9, 0.15,
     ["7304: Huntress & Pinizzotto, J. Chem. Phys. 59, 4742 (1973)",
      "7202: Bolden & Twiddy, Faraday Disc. 53, 192 (1972)"],
     "Symmetric charge transfer (H2O+ + H2O -> H2O + H2O+) is invisible in thermal rate measurements and is not part of this entry."),
    ("ohplus-h2o", 1511, 44, "OH+ + H2O", [("H2O+ + OH", 0.55), ("H3O+ + O", 0.45)], 2.89e-9, 0.10,
     ["7304: Huntress & Pinizzotto, J. Chem. Phys. 59, 4742 (1973)"], ""),
    ("hplus-h2o", 1473, 6, "H+ + H2O", [("H2O+ + H", 1.00)], 8.20e-9, 0.10,
     ["7401: Huntress, Kim & Theard, Chem. Phys. Lett. 29, 189 (1974); Huntress, Planet. Space Sci. 23, 377 (1975)"], ""),
    ("oplus-h2o", 1509, 42, "O+ + H2O", [("H2O+ + O", 1.00)], 2.60e-9, 0.15,
     ["8925: Heninger, Fenistein, Mauclaire, Marx & Murad, Geophys. Res. Lett. 16, 139 (1989)",
      "8010: Adams, Smith & Paulson, J. Chem. Phys. 72, 288 (1980); Smith, Adams & Miller, J. Chem. Phys. 69, 308 (1978)",
      "7202: Bolden & Twiddy, Faraday Disc. 53, 192 (1972)",
      "7006: Howard, Rundle & Kaufman, J. Chem. Phys. 53, 3745 (1970)"],
     "O+ electronic state not specified in the table entry (ground-state O+(4S) is implied for thermal flow-tube data)."),
    ("h2plus-h2o", 1474, 7, "H2+ + H2O", [("H2O+ + H2", 0.53), ("H3O+ + H", 0.47)], 7.30e-9, 0.10,
     ["7503: Kim & Huntress, J. Chem. Phys. 62, 2820 (1975)"], ""),
]


def main():
    out = sys.argv[1]
    for key, page, pdfp, reac, chans, k, ru, refs, extra in REACTIONS:
        loc = "J. Phys. Chem. Ref. Data 22, p. %d, Table of Reactions, entry '%s' (NIST reprint PDF page %d)" % (page, reac, pdfp)
        proc = reac + " -> " + " | ".join(p for p, _ in chans)
        t = {
            "table_id": "ion-anicich1993-%s-k300" % key,
            "source_id": SRC, "citation": CIT, "doi": "10.1063/1.555940",
            "locator": loc,
            "process": proc + " (total rate coefficient)",
            "quantity": "rate_coefficient_m3_s",
            "x_name": "temperature_K",
            "x": [300.0],
            "y": [k * 1e-6],
            "y_uncertainty": [k * 1e-6 * ru],
            "y_uncertainty_relative": ru,
            "original_units": "cm3/s; printed value %.3g cm3/s +/- %d%%" % (k, round(ru * 100)),
            "conversion_to_SI": "k[m3/s] = k[cm3/s] * 1e-6",
            "tag": "recommended",
            "extraction_method": "manual transcription from rendered page image of the evaluated Table of Reactions",
            "evaluated_from": refs,
            "notes": BASE_NOTES + (" " + extra if extra else ""),
        }
        with open(os.path.join(out, t["table_id"] + ".json"), "w") as f:
            json.dump(t, f, indent=1)
        print("wrote", t["table_id"])
        if len(chans) > 1:
            b = {
                "table_id": "ion-anicich1993-%s-branching300" % key,
                "source_id": SRC, "citation": CIT, "doi": "10.1063/1.555940",
                "locator": loc + ", column 'Prod. Dist.'",
                "process": proc,
                "quantity": "branching_fraction",
                "x_name": "temperature_K",
                "x": [300.0] * len(chans),
                "y": [f for _, f in chans],
                "y_labels": [reac + " -> " + p for p, _ in chans],
                "y_uncertainty": None,
                "original_units": "dimensionless product distribution",
                "conversion_to_SI": "none (dimensionless)",
                "tag": "recommended",
                "extraction_method": "manual transcription from rendered page image of the evaluated Table of Reactions",
                "evaluated_from": refs,
                "notes": "One y entry per product channel (see y_labels); all at the same nominal 300 K. No uncertainty "
                         "is printed for product distributions. " + BASE_NOTES,
            }
            with open(os.path.join(out, b["table_id"] + ".json"), "w") as f:
                json.dump(b, f, indent=1)
            print("wrote", b["table_id"])


if __name__ == "__main__":
    main()
