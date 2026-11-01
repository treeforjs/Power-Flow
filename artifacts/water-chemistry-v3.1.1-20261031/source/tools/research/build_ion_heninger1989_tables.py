"""Write O+(4S) + H2O rate coefficient and cross section vs collision energy transcribed from
Heninger, Fenistein, Mauclaire, Marx, Murad, Geophys. Res. Lett. 16, 139 (1989), Table 1
(DTIC reprint ADA208592, Internet Archive copy).

Usage: python3 -I build_ion_heninger1989_tables.py <extracted_dir>
"""
import json
import os
import sys

SRC = "ion-heninger-1989-oplus-h2o-dtic-ada208592"
CIT = ("M. Heninger, S. Fenistein, G. Mauclaire, R. Marx, E. Murad, Review of the reaction of O+ with H2O and its "
       "bearing on composition measurements from the space shuttle, Geophys. Res. Lett. 16, 139-141 (1989)")
E = [0.0327, 0.23, 0.46, 0.58, 0.92, 1.46, 1.84, 2.92, 3.68, 4.63, 5.83, 7.3, 9.2, 11.6]
K = [2.5, 2.9, 2.8, 2.8, 2.6, 2.5, 2.4, 2.0, 2.0, 2.1, 1.8, 1.9, 2.6, 2.2]          # 1e-9 cm3/s
SIG = [230, 125, 85.5, 77.7, 57.6, 42.8, 37, 24.9, 21.9, 20.2, 15.4, 14.8, 17.7, 13.5]  # 1e-16 cm2
NOTES = ("Triple-cell ICR; ground-state O+(4S) from 30 eV electron impact on CO2; O+ kinetic energy set by resonant "
         "cyclotron excitation, so the energy distribution width is stated by the source to be difficult to ascertain. "
         "Rate obtained from the decay of O+ versus time with H2O at ~1e-6 Torr (total O+ removal by H2O; the product is "
         "written as H2O+ + O in reaction (1), product ions were not separately quantified in Table 1). H2O pressure "
         "calibrated via Kr+ + H2O at thermal energy. The source states rate constants are accurate to +/- 0.5e-9 cm3/s; "
         "no separate uncertainty is given for the cross sections. The source notes its cross sections are ~3x (0.2 eV) and "
         "~2x (10 eV) higher than the earlier beam data of Murad and Lai and considers the ICR data more reliable.")


def write(out, t):
    with open(os.path.join(out, t["table_id"] + ".json"), "w") as f:
        json.dump(t, f, indent=1)
    print("wrote", t["table_id"])


def main():
    out = sys.argv[1]
    base = {"source_id": SRC, "citation": CIT, "doi": "10.1029/GL016i002p00139",
            "locator": "Geophys. Res. Lett. 16, p. 140, Table 1 (DTIC ADA208592 reprint PDF page 2)",
            "process": "O+(4S) + H2O -> H2O+ + O (reaction (1); total O+ loss rate)",
            "x_name": "energy_eV_cm", "x": E, "tag": "measured",
            "extraction_method": "manual transcription from printed Table 1 (rendered at 300 dpi; text layer agrees)"}
    k = dict(base)
    k.update({"table_id": "ion-heninger1989-oplus-h2o-k-vs-ecm", "quantity": "rate_coefficient_m3_s",
              "y": [round(v * 1e-15, 20) for v in K], "y_uncertainty": [0.5e-15] * len(K),
              "original_units": "1e-9 cm3 molecule-1 s-1 (printed values: " + ", ".join(str(v) for v in K) + ")",
              "conversion_to_SI": "k[m3/s] = printed value * 1e-9 * 1e-6",
              "notes": NOTES + " k at a given E_cm is a rate constant for a non-thermal ion energy distribution centred at that E_cm, not a Maxwellian k(T)."})
    write(out, k)
    s = dict(base)
    s.update({"table_id": "ion-heninger1989-oplus-h2o-sigma-vs-ecm", "quantity": "cross_section_m2",
              "y": [round(v * 1e-20, 26) for v in SIG], "y_uncertainty": None,
              "original_units": "1e-16 cm2 (printed values: " + ", ".join(str(v) for v in SIG) + ")",
              "conversion_to_SI": "sigma[m2] = printed value * 1e-16 * 1e-4",
              "notes": NOTES + " Cross sections are as printed by the source (derived from k and the relative velocity)."})
    write(out, s)


if __name__ == "__main__":
    main()
