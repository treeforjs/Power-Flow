"""Write normalized tables transcribed from Lishawa et al., J. Chem. Phys. 93, 3196 (1990)
(DTIC reprint PL-TR-91-2018, ADA232204). Numbers were transcribed manually from the printed
Tables I-III of the reprint (rendered at 250-400 dpi and checked against the text layer).

Usage: python3 -I build_ion_lishawa_tables.py <extracted_dir>
"""
import json
import os
import sys

SRC = "ion-lishawa-1990-h2oplus-h2o"
CIT = "C. R. Lishawa, R. A. Dressler, J. A. Gardner, R. H. Salter, E. Murad, J. Chem. Phys. 93, 3196-3206 (1990)"
NA = None

COMMON_NOTES = (
    "Product ion laboratory energies are TOF band centres (secondary ion laboratory kinetic energy), not centre-of-mass "
    "kinetic energy release; the source converts TOF to energy assuming production near the collision chamber exit and 0 deg "
    "scattering. Number in parentheses in the source (normalized statistical weight of the band) is kept as *_weight. "
    "Primary ions were produced by electron impact at 20-30 eV; A-state ions are argued to have decayed, vibrational "
    "excitation of the X state cannot be excluded (source Sec. II). Isotopic targets/projectiles (D2O, D2-18O, H2-18O) "
    "were used to separate channels; masses therefore differ from H2O+ + H2O."
)


def table(table_id, locator, process, x, y, yu, extra, notes):
    return {
        "table_id": table_id,
        "source_id": SRC,
        "citation": CIT,
        "doi": "10.1063/1.458852",
        "locator": locator,
        "process": process,
        "quantity": "product_ion_laboratory_energy_eV",
        "x_name": "energy_eV_cm",
        "x": x,
        "y": y,
        "y_uncertainty": yu,
        "original_units": "eV (laboratory frame product ion energy); x in eV centre-of-mass",
        "conversion_to_SI": "1 eV = 1.602176634e-19 J (no conversion applied; values kept in eV)",
        "tag": "measured",
        "extraction_method": "manual transcription from printed table in the DTIC reprint scan (visual check at 250-400 dpi against PDF text layer)",
        "additional_columns": extra,
        "notes": COMMON_NOTES + " " + notes,
    }


def main():
    out = sys.argv[1]
    tabs = []
    # Table I: atom pickup (2a) H2-18O+ + D2O -> H2D-18O+ + OD
    tabs.append(table(
        "ion-lishawa1990-t1-atom-pickup-2a-product-lab-energy",
        "J. Chem. Phys. 93, 3198, Table I (DTIC reprint PDF page 5)",
        "H2(18)O+ + D2O -> H2D(18)O+ + OD (atom pickup, reaction 2a)",
        [2.15, 4.65, 6.95, 9.75, 14.8, 19.8],
        [2.7, 7.1, 10.3, 13.8, 23.7, 32.9],
        [1.2, 1.5, 1.7, 2.2, 2.3, 3.0],
        {"E_prim_lab_eV": [4.3, 9.3, 13.9, 19.5, 29.5, 39.5]},
        "Single forward band; source states the atom pickup cross section is at least an order of magnitude below proton transfer.",
    ))
    # Table II (3a)
    tabs.append(table(
        "ion-lishawa1990-t2-proton-transfer-3a-product-lab-energy",
        "J. Chem. Phys. 93, 3199, Table II upper block (DTIC reprint PDF page 6)",
        "H2O+ + D2(18)O -> OH + HD2(18)O+ (proton transfer, reaction 3a)",
        [2.8, 5.2, 7.9, 10.5, 16.1, 20.0, 24.1],
        [0.04, 0.06, 0.13, 0.10, 0.11, 0.07, 0.15],
        [0.05, 0.07, 0.07, 0.05, 0.05, 0.07, 0.07],
        {
            "E_prim_lab_eV": [5.1, 9.5, 14.3, 19.0, 29.1, 36.4, 43.8],
            "band1_weight": [1.00, 0.99, 0.99, 0.97, 0.97, 0.89, 0.86],
            "band2_eV": [NA, NA, NA, 3.3, 10.0, 12.8, 14.6],
            "band2_uncertainty_eV": [NA, NA, NA, 1.2, 1.8, 2.0, 2.3],
            "band2_weight": [NA, NA, NA, 0.03, 0.03, 0.10, 0.14],
            "band3_eV": [NA, 6.4, 8.1, 14.4, NA, 32.0, NA],
            "band3_uncertainty_eV": [NA, 1.2, 1.8, 2.5, NA, 3.0, NA],
            "band3_weight": [NA, 0.01, 0.01, "<0.01", NA, 0.01, NA],
            "E_sec_zero_KER_calc_eV": [1.3, 2.5, 3.7, 4.9, 7.5, 9.4, 11.3],
        },
        "y is the 1st (near-thermal) band. Weight of band 1 was multiplied by 6 in the source to correct the TOF-mode "
        "extraction efficiency. Below E_cm ~12 eV the 2nd band weights are stated to be estimates. Source estimates the "
        "cross sections of 2nd and 3rd band ions to be less than 5e-18 and 1e-18 cm2 respectively (text p. 3204). "
        "Null = '...' (not observed/determined) in the source.",
    ))
    tabs.append(table(
        "ion-lishawa1990-t2-deuteron-transfer-3b-product-lab-energy",
        "J. Chem. Phys. 93, 3199, Table II lower block (DTIC reprint PDF page 6)",
        "D2O+ + D2(18)O -> OD + D3(18)O+ (deuteron transfer, reaction 3b)",
        [1.3, 2.1, 4.7, 7.3, 10.1, 12.6, 15.2, 17.9, 20.6],
        [0.23, 0.07, 0.12, 0.13, 0.18, 0.19, 0.21, 0.09, 0.17],
        [0.10, 0.07, 0.05, 0.07, 0.08, 0.10, 0.10, 0.15, 0.15],
        {
            "E_prim_lab_eV": [2.5, 4.1, 9.0, 14.0, 19.2, 24.1, 29.1, 34.1, 39.4],
            "band1_weight": [1.00, 1.00, 0.98, 0.99, 0.98, 0.95, 0.94, 0.95, 0.81],
            "band2_eV": [NA, NA, NA, NA, 5.5, 6.2, 8.6, 9.3, 15.1],
            "band2_uncertainty_eV": [NA, NA, NA, NA, 1.2, 1.2, 1.5, 1.5, 2.0],
            "band2_weight": [NA, NA, NA, NA, 0.01, 0.04, 0.05, 0.05, 0.19],
            "band3_eV": [NA, NA, 6.7, 9.5, 16.1, 18.2, 22.7, NA, NA],
            "band3_uncertainty_eV": [NA, NA, 1.3, 1.3, 1.8, 2.2, 2.3, NA, NA],
            "band3_weight": [NA, NA, 0.02, 0.01, 0.01, 0.01, 0.01, NA, NA],
            "E_sec_zero_KER_calc_eV": [0.7, 1.1, 2.5, 3.8, 5.2, 6.6, 7.9, 9.3, 10.7],
        },
        "y is the 1st (near-thermal) band. Null = '...' in the source. The PDF text layer misreads two cells "
        "(0.18 +/- 0.08 and 18.2 +/- 2.2); the printed image values are used.",
    ))
    # Table III (4a), (4b): symmetric charge exchange
    tabs.append(table(
        "ion-lishawa1990-t3-charge-exchange-4a-product-lab-energy",
        "J. Chem. Phys. 93, 3201, Table III upper block (DTIC reprint PDF page 8)",
        "H2O+ + D2(18)O -> H2O + D2(18)O+ (symmetric charge exchange, reaction 4a)",
        [2.81, 6.05, 10.6, 16.1, 21.8, 27.1],
        [0.03, 0.03, 0.16, 0.03, 0.05, 0.05],
        [0.05, 0.05, 0.08, 0.05, 0.05, 0.05],
        {
            "E_prim_lab_eV": [5.1, 11.0, 19.3, 29.2, 39.6, 49.2],
            "band1_weight": [NA, NA, 0.90, 0.95, 0.96, 0.95],
            "band2_eV": [NA, NA, 4.9, 11.8, 18.0, 21.0],
            "band2_uncertainty_eV": [NA, NA, 1.2, 1.5, 2.2, 2.2],
            "band2_weight": [NA, NA, 0.10, 0.05, 0.04, 0.05],
            "E_sec_zero_KER_calc_eV": [1.3, 2.7, 4.8, 7.2, 9.8, 12.2],
        },
        "y is the 1st (near-thermal) band: charge exchange product ions are formed essentially at thermal laboratory energy. "
        "Source: below E_cm = 10 eV the 2nd band of 4a is poorly resolved (no weights given). Source estimates 2nd-band "
        "cross section below 1e-17 cm2 (text p. 3204) and attributes 2nd-band ions at least partly to dissociating oxonium ions.",
    ))
    tabs.append(table(
        "ion-lishawa1990-t3-charge-exchange-4b-product-lab-energy",
        "J. Chem. Phys. 93, 3201, Table III lower block (DTIC reprint PDF page 8)",
        "D2(18)O+ + H2O -> D2(18)O + H2O+ (symmetric charge exchange, reaction 4b)",
        [3.96, 8.64, 12.7, 17.1, 21.7],
        [0.06, 0.05, 0.07, 0.06, 0.06],
        [0.05, 0.05, 0.08, 0.05, 0.05],
        {
            "E_prim_lab_eV": [8.8, 19.2, 28.3, 38.3, 48.3],
            "band1_weight": [0.99, 0.98, 0.97, 0.98, 0.93],
            "band2_eV": [5.0, 7.8, 9.9, 16.4, 21.8],
            "band2_uncertainty_eV": [1.2, 1.3, 1.5, 2.0, 2.2],
            "band2_weight": [0.01, 0.02, 0.03, 0.02, 0.07],
            "E_sec_zero_KER_calc_eV": [2.2, 4.8, 7.0, 9.4, 9.8],
        },
        "y is the 1st (near-thermal) band. Printed inconsistency retained, not corrected: the last E_sec(E_kin=0) entry "
        "is printed as 9.8 eV at E_prim = 48.3 eV, whereas the source's Eq. (8) with m3=18, m1=22, m2=18 gives about 12.0 eV "
        "(the other rows of this block agree with Eq. (8) to rounding).",
    ))
    for t in tabs:
        with open(os.path.join(out, t["table_id"] + ".json"), "w") as f:
            json.dump(t, f, indent=1)
        print("wrote", t["table_id"])


if __name__ == "__main__":
    main()
