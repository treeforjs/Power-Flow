"""Write negative-ion tables transcribed from retained primary sources.

Usage: python3 -I build_ion_neg_tables.py <extracted_dir>
"""
import json
import os
import sys


def main():
    out = sys.argv[1]
    t = {
        "table_id": "ion-neg-miller2011-hminus-h-ad-k-16meV",
        "source_id": "ion-neg-miller-2011-hminus-h-ad",
        "citation": "K. A. Miller, H. Bruhns, J. Eliasek, M. Cizek, H. Kreckel, X. Urbain, D. W. Savin, Associative detachment of H- + H -> H2 + e-, Phys. Rev. A 84, 052709 (2011); arXiv:1110.0683",
        "doi": "10.1103/PhysRevA.84.052709",
        "locator": "arXiv:1110.0683v1, Table II (PDF page 26); systematic uncertainty Table I (PDF page 25)",
        "process": "H- + H -> H2 + e- (associative detachment), merged beams",
        "quantity": "rate_coefficient_m3_s",
        "x_name": "energy_eV_cm",
        "x": [0.016, 0.016, 0.016],
        "y": [5.6e-15, 5.2e-15, 5.2e-15],
        "y_uncertainty": [0.7e-15, 0.6e-15, 0.6e-15],
        "additional_columns": {"helium_gas_cell_pressure_1e-4_Torr": [1.0, 2.0, 3.0],
                               "statistical_uncertainty_m3_s": [0.5e-15, 0.4e-15, 0.4e-15],
                               "theory_same_energy_m3_s": 5.0e-15},
        "original_units": "1e-9 cm3/s (printed 5.6, 5.2, 5.2; statistical +/-0.5, 0.4, 0.4; 'relative' +/-0.7, 0.6, 0.6)",
        "conversion_to_SI": "k[m3/s] = printed * 1e-9 * 1e-6",
        "tag": "measured",
        "extraction_method": "transcribed from the arXiv preprint text layer (born-digital)",
        "notes": ("Merged-beams rate coefficient <sigma v_r> at mean relative energy <E_r> = 16 meV, measured at three "
                  "helium stripping-cell pressures as a linearity test (not three energies). y_uncertainty is the source's "
                  "'relative uncertainty' column; Table I gives a total 1-sigma systematic uncertainty of 24%. The energy-"
                  "resolved data (4 meV to 4.83 eV) are only in Figs. 6-7 (not digitized). The source recommends the "
                  "thermal rate coefficient of Kreckel et al., Science 329, 69 (2010) (its ref. [3]) for T < 1e4 K; that "
                  "expression was not retrieved. The source notes that its results lie 2.2 +/- 0.9 times above the 300 K "
                  "flowing-afterglow value of Martinez et al., Astrophys. J. 705 (2009) (its ref. [4]), attributing the "
                  "discrepancy to the normalization data of Howard, Fehsenfeld & McFarland 1974 (its ref. [7]). Relevance "
                  "to a water discharge: H- and H are both minor fragments; H- + H AD only matters where both are abundant."),
    }
    with open(os.path.join(out, t["table_id"] + ".json"), "w") as f:
        json.dump(t, f, indent=1)
    print("wrote", t["table_id"])


if __name__ == "__main__":
    main()
