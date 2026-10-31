"""Write extracted tables for channel 5 (rec-ii-) and channel 6 (neutral-).

All numbers below are transcribed by hand from the retained source PDFs and
checked against rendered page images; nothing is computed except unit
conversion (cm3 -> m3: x1e-6; cm6 -> m6: x1e-12; cm2 -> m2: x1e-4).
Rate expressions are stored as parameters only; they are not evaluated onto
temperature grids.

Usage: python3 -I build_neutral_rec_ii_tables.py <extracted_dir>
"""
import json
import os
import sys

CM3 = 1e-6
CM6 = 1e-12


def r(v):
    """Round away binary noise from unit conversion (13 significant digits)."""
    if v is None:
        return None
    return float("%.13g" % v)


def write(out_dir, t):
    path = os.path.join(out_dir, t["table_id"] + ".json")
    with open(path, "w") as f:
        json.dump(t, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print("wrote", path)


IUPAC_SRC = "neutral-iupac-atkinson-2004-vol1"
IUPAC_CIT = ("Atkinson et al., Atmos. Chem. Phys. 4, 1461-1738 (2004), "
             "doi:10.5194/acp-4-1461-2004")
B05_SRC = "neutral-baulch-2005-combustion-suppl2"
B05_CIT = ("Baulch et al., J. Phys. Chem. Ref. Data 34, 757-1397 (2005), "
           "doi:10.1063/1.1748524")
TH_SRC = "neutral-tsang-hampson-1986-combustion"
TH_CIT = ("Tsang and Hampson, J. Phys. Chem. Ref. Data 15, 1087-1279 (1986), "
          "doi:10.1063/1.555759")

CONV_K = "cm3 molecule-1 s-1 multiplied by 1e-6 gives m3 s-1"
CONV_K3 = "cm6 molecule-2 s-1 multiplied by 1e-12 gives m6 s-1"


def expr(form, params_orig, units_orig, tmin, tmax, reliability, si_params):
    return {
        "form": form,
        "parameters_original_units": params_orig,
        "original_units": units_orig,
        "parameters_SI": si_params,
        "validity_T_K": [tmin, tmax],
        "reliability": reliability,
        "note": "Parameters as printed by the source; not evaluated onto a grid here.",
    }


def tables():
    T = []

    # ---------------- Channel 6: IUPAC 2004 ----------------
    T.append({
        "table_id": "neutral-iupac2004-o1d-h2o-total-k",
        "source_id": IUPAC_SRC,
        "locator": "Data sheet I.A2.19, printed p. 1510 (PDF page 50); also Table 1 summary, reaction 19, printed p. 1463",
        "process": "O(1D) + H2O -> products (total k = k1 + k2 + k3; k1: 2 OH, k2: H2 + O2, k3: O(3P) + H2O)",
        "quantity": "rate_coefficient_m3_s",
        "x_name": "temperature_K",
        "x": [298],
        "y": [r(2.2e-10 * CM3)],
        "y_uncertainty": [{"type": "delta_log10_k", "value": 0.1}],
        "original_units": "cm3 molecule-1 s-1",
        "conversion_to_SI": CONV_K,
        "tag": "recommended",
        "extraction_method": "manual transcription from text layer, verified against rendered page image",
        "expression": expr("k = A (temperature independent)", {"A": 2.2e-10}, "cm3 molecule-1 s-1",
                           200, 350, "delta log10 k = +/-0.1 at 298 K; delta(E/R) = +/-100 K",
                           {"A_m3_s": r(2.2e-10 * CM3)}),
        "notes": ("IUPAC preferred value: k = 2.2e-10 cm3 molecule-1 s-1 independent of temperature over 200-350 K; "
                  "k1 (2 OH) = 2.2e-10 at 298 K; k2 (H2 + O2) < 2.2e-12 at 298 K; k3 (O(3P) + H2O) < 1.2e-11 at 298 K. "
                  "Preferred k is the mean of Streit et al. 1976, Amimoto et al. 1979, Lee and Slanger 1979, Gericke and Comes 1981, "
                  "Wine and Ravishankara 1981 (all thermal, laser-photolysis). Thermal evaluation only; no translational-energy "
                  "dependence for suprathermal O(1D) and no OH product state/recoil distribution is given in this data sheet."),
    })
    T.append({
        "table_id": "neutral-iupac2004-o1d-h2o-branching",
        "source_id": IUPAC_SRC,
        "locator": "Data sheet I.A2.19, 'Branching Ratios' rows of the rate-coefficient table, printed p. 1510 (PDF page 50)",
        "process": "O(1D) + H2O: minor-channel branching ratios k2/k (H2 + O2) and k3/k (O(3P) + H2O)",
        "quantity": "branching_fraction",
        "x_name": "temperature_K",
        "x": [298, 298, 298],
        "y": [0.01, 0.049, 0.006],
        "y_labels": ["k2/k, H2 + O2 (Zellner et al. 1980, FP-RA)",
                     "k3/k, O(3P) + H2O (Wine and Ravishankara 1982, PLP-RF; temperature 298 +/- 2 K)",
                     "k2/k, H2 + O2 (Glinski and Birks 1985)"],
        "y_uncertainty": [{"plus": 0.005, "minus": 0.01}, {"plus": 0.032, "minus": 0.032},
                          {"plus": 0.007, "minus": 0.006}],
        "original_units": "dimensionless",
        "conversion_to_SI": "none",
        "tag": "measured",
        "extraction_method": "manual transcription from text layer, verified against rendered page image",
        "notes": ("Measured values as quoted by the IUPAC compilation (primary papers not retained). IUPAC preferred upper "
                  "limits at 298 K: k2 < 2.2e-12 and k3 < 1.2e-11 cm3 molecule-1 s-1, i.e. the 2 OH channel dominates. "
                  "The remaining fraction (2 OH) is not a separately stated measured value and is not filled here."),
    })
    T.append({
        "table_id": "neutral-iupac2004-oh-oh-h2o-o",
        "source_id": IUPAC_SRC,
        "locator": "Data sheet I.A2.21, printed p. 1514 (PDF page 54); Table 1 summary reaction 21, printed p. 1463",
        "process": "OH + OH -> H2O + O(3P)",
        "quantity": "rate_coefficient_m3_s",
        "x_name": "temperature_K",
        "x": [298],
        "y": [r(1.48e-12 * CM3)],
        "y_uncertainty": [{"type": "delta_log10_k", "value": 0.15}],
        "original_units": "cm3 molecule-1 s-1",
        "conversion_to_SI": CONV_K,
        "tag": "recommended",
        "extraction_method": "manual transcription from text layer",
        "expression": expr("k = A (T/298)^n exp(-E_over_R/T)", {"A": 6.2e-14, "n": 2.6, "E_over_R_K": -945},
                           "cm3 molecule-1 s-1", 200, 350, "delta log10 k = +/-0.15 at 298 K; delta(E/R) = +/-250 K",
                           {"A_m3_s": r(6.2e-14 * CM3), "n": 2.6, "T_ref_K": 298, "E_over_R_K": -945}),
        "notes": "Printed form: 6.2e-14 (T/298)^2.6 exp(945/T); E_over_R stored as -945 K so that exp(-E_over_R/T) = exp(+945/T).",
    })
    T.append({
        "table_id": "neutral-iupac2004-o-oh-o2-h",
        "source_id": IUPAC_SRC,
        "locator": "Data sheet I.A2.15, printed p. 1505 (PDF page 45); Table 1 summary reaction 15, printed p. 1463",
        "process": "O(3P) + OH -> O2 + H",
        "quantity": "rate_coefficient_m3_s",
        "x_name": "temperature_K",
        "x": [298],
        "y": [r(3.5e-11 * CM3)],
        "y_uncertainty": [{"type": "delta_log10_k", "value": 0.1}],
        "original_units": "cm3 molecule-1 s-1",
        "conversion_to_SI": CONV_K,
        "tag": "recommended",
        "extraction_method": "manual transcription from text layer",
        "expression": expr("k = A exp(-E_over_R/T)", {"A": 2.4e-11, "E_over_R_K": -110}, "cm3 molecule-1 s-1",
                           150, 500, "delta log10 k = +/-0.1 at 298 K; delta(E/R) = +/-100 K",
                           {"A_m3_s": r(2.4e-11 * CM3), "E_over_R_K": -110}),
        "notes": "Printed form: 2.4e-11 exp(110/T) over 150-500 K. Fragment-fragment reaction (O + OH), not a reaction with H2O.",
    })
    T.append({
        "table_id": "neutral-iupac2004-oh-h2-h2o-h",
        "source_id": IUPAC_SRC,
        "locator": "Data sheet I.A2.20, printed p. 1512 (PDF page 52); Table 1 summary reaction 20, printed p. 1463",
        "process": "OH + H2 -> H2O + H",
        "quantity": "rate_coefficient_m3_s",
        "x_name": "temperature_K",
        "x": [298],
        "y": [r(6.7e-15 * CM3)],
        "y_uncertainty": [{"type": "delta_log10_k", "value": 0.1}],
        "original_units": "cm3 molecule-1 s-1",
        "conversion_to_SI": CONV_K,
        "tag": "recommended",
        "extraction_method": "manual transcription from text layer",
        "expression": expr("k = A exp(-E_over_R/T)", {"A": 7.7e-12, "E_over_R_K": 2100}, "cm3 molecule-1 s-1",
                           200, 450, "delta log10 k = +/-0.1 at 298 K; delta(E/R) = +/-200 K",
                           {"A_m3_s": r(7.7e-12 * CM3), "E_over_R_K": 2100}),
        "notes": "Reverse of H + H2O -> OH + H2. High-temperature form in neutral-baulch2005-oh-h2-h2o-h.",
    })
    T.append({
        "table_id": "neutral-iupac2004-o1d-h2-oh-h",
        "source_id": IUPAC_SRC,
        "locator": "Data sheet I.A2.18, printed p. 1508 (PDF page 48); Table 1 summary reaction 18, printed p. 1463",
        "process": "O(1D) + H2 -> OH + H (total k including O(3P) + H2 quenching)",
        "quantity": "rate_coefficient_m3_s",
        "x_name": "temperature_K",
        "x": [298],
        "y": [r(1.1e-10 * CM3)],
        "y_uncertainty": [{"type": "delta_log10_k", "value": 0.1}],
        "original_units": "cm3 molecule-1 s-1",
        "conversion_to_SI": CONV_K,
        "tag": "recommended",
        "extraction_method": "manual transcription from text layer",
        "expression": expr("k = A (temperature independent)", {"A": 1.1e-10}, "cm3 molecule-1 s-1", 200, 350,
                           "delta log10 k = +/-0.1 at 298 K; delta(E/R) = +/-100 K", {"A_m3_s": r(1.1e-10 * CM3)}),
        "notes": ("Total k = k1 (OH + H) + k2 (O(3P) + H2). Source comment: channel (1) appears to be the dominant pathway "
                  "(>95%) (Wine and Ravishankara, 1982). Fragment reaction (O(1D) + H2)."),
    })

    # ---------------- Channel 6: Baulch 2005 ----------------
    T.append({
        "table_id": "neutral-baulch2005-h-h2o-oh-h2",
        "source_id": B05_SRC,
        "locator": "Data sheet 'H + H2O -> OH + H2', Preferred Values, printed p. 922 (PDF page 166)",
        "process": "H + H2O -> OH + H2 (endothermic, dH298 = +60.96 kJ/mol)",
        "quantity": "rate_coefficient_m3_s",
        "x_name": "temperature_K",
        "x": [],
        "y": [],
        "y_uncertainty": None,
        "original_units": "cm3 molecule-1 s-1",
        "conversion_to_SI": CONV_K,
        "tag": "recommended",
        "extraction_method": "manual transcription, verified against rendered page image",
        "expression": expr("k = A T^n exp(-E_over_R/T), T in K", {"A": 7.5e-16, "n": 1.6, "E_over_R_K": 9030},
                           "cm3 molecule-1 s-1", 800, 2500, "delta log10 k = +/-0.2 over 800-2500 K",
                           {"A_m3_s": r(7.5e-16 * CM3), "n": 1.6, "T_ref_K": 1, "E_over_R_K": 9030}),
        "notes": ("Source states no point values; x/y left empty and the preferred expression is stored as parameters. "
                  "Direct measurements exist only above 1160 K (Madronich and Felder 1984; Michael and Sutherland 1988; "
                  "Davidson et al. 1988); the preferred expression is derived from k(OH + H2) and the equilibrium constant. "
                  "Thermal Maxwellian coefficient only; not an energy-resolved cross section for fast H atoms."),
    })
    T.append({
        "table_id": "neutral-baulch2005-oh-h2-h2o-h",
        "source_id": B05_SRC,
        "locator": "Data sheet 'OH + H2 -> H2O + H', Preferred Values, printed p. 1029 (PDF page 273)",
        "process": "OH + H2 -> H2O + H",
        "quantity": "rate_coefficient_m3_s",
        "x_name": "temperature_K",
        "x": [],
        "y": [],
        "y_uncertainty": None,
        "original_units": "cm3 molecule-1 s-1",
        "conversion_to_SI": CONV_K,
        "tag": "recommended",
        "extraction_method": "manual transcription, verified against rendered page image",
        "expression": expr("k = A T^n exp(-E_over_R/T), T in K", {"A": 3.6e-16, "n": 1.52, "E_over_R_K": 1740},
                           "cm3 molecule-1 s-1", 250, 2500,
                           "delta log10 k = +/-0.1 at 250 K, rising to +/-0.3 at 2500 K",
                           {"A_m3_s": r(3.6e-16 * CM3), "n": 1.52, "T_ref_K": 1, "E_over_R_K": 1740}),
        "notes": "Source states no point values; expression stored as parameters only.",
    })
    T.append({
        "table_id": "neutral-baulch2005-oh-oh-h2o-o",
        "source_id": B05_SRC,
        "locator": "Data sheet 'OH + OH -> H2O + O', Preferred Values, printed p. 1032 (PDF page 276)",
        "process": "OH + OH -> H2O + O(3P)",
        "quantity": "rate_coefficient_m3_s",
        "x_name": "temperature_K",
        "x": [],
        "y": [],
        "y_uncertainty": None,
        "original_units": "cm3 molecule-1 s-1",
        "conversion_to_SI": CONV_K,
        "tag": "recommended",
        "extraction_method": "manual transcription, verified against rendered page image",
        "expression": expr("k = A T^n exp(-E_over_R/T), T in K", {"A": 5.56e-20, "n": 2.42, "E_over_R_K": -970},
                           "cm3 molecule-1 s-1", 250, 2400, "delta log10 k = +/-0.15 over 250-2400 K",
                           {"A_m3_s": r(5.56e-20 * CM3), "n": 2.42, "T_ref_K": 1, "E_over_R_K": -970}),
        "notes": ("Printed form: 5.56e-20 T^2.42 exp(970/T). Source adopts the NASA/IUPAC 298 K value "
                  "(comment j: based on 1.5e-12 at 298 K). See also neutral-iupac2004-oh-oh-h2o-o (200-350 K)."),
    })
    T.append({
        "table_id": "neutral-baulch2005-h-oh-m-h2o",
        "source_id": B05_SRC,
        "locator": "Data sheet 'H + OH (+M) -> H2O (+M)', Preferred Values, printed p. 913 (PDF page 157)",
        "process": "H + OH + M -> H2O + M, low-pressure limit k0 with M = H2O (three-body)",
        "quantity": "termolecular_rate_coefficient_m6_s",
        "x_name": "temperature_K",
        "x": [],
        "y": [],
        "y_uncertainty": None,
        "original_units": "cm6 molecule-2 s-1",
        "conversion_to_SI": CONV_K3,
        "tag": "recommended",
        "extraction_method": "manual transcription, verified against rendered page image",
        "expression": expr("k0 = A T^n, T in K", {"A": 3.9e-25, "n": -2.0}, "cm6 molecule-2 s-1", 300, 3000,
                           "delta log10 k0 = +/-0.5 for M = H2O over 300-3000 K",
                           {"A_m6_s": r(3.9e-25 * CM6), "n": -2.0, "T_ref_K": 1}),
        "notes": ("Same sheet gives k0 = 2.3e-26 T^-2.0 (M = Ar, +/-0.3) and 6.1e-26 T^-2.0 (M = N2, +/-0.5), both 300-3000 K. "
                  "M = H2O expression accepted from Baulch et al. (1972) evaluation; measured points are flame/shock-tube "
                  "data at 1400-1912 K. Three-body; listed for completeness."),
    })

    # ---------------- Channel 6: Tsang and Hampson 1986 ----------------
    T.append({
        "table_id": "neutral-tsanghampson1986-oh-h-o-h2",
        "source_id": TH_SRC,
        "locator": "Data sheet 6,4 (a) 'OH + H -> O + H2', printed p. 1114 (PDF page 29); summary index printed p. 1091",
        "process": "OH + H -> O(3P) + H2",
        "quantity": "rate_coefficient_m3_s",
        "x_name": "temperature_K",
        "x": [],
        "y": [],
        "y_uncertainty": None,
        "original_units": "cm3 molecule-1 s-1",
        "conversion_to_SI": CONV_K,
        "tag": "recommended",
        "extraction_method": "manual transcription, verified against rendered page image",
        "expression": expr("k = A T^n exp(-E_over_R/T), T in K", {"A": 8.1e-21, "n": 2.8, "E_over_R_K": 1950},
                           "cm3 molecule-1 s-1", 300, 2500, "uncertainty factor 2",
                           {"A_m3_s": r(8.1e-21 * CM3), "n": 2.8, "T_ref_K": 1, "E_over_R_K": 1950}),
        "notes": ("Source comment: there has been no direct measurement of k_a; the expression (Cohen and Westberg 1983 review) "
                  "is based on the reverse process (O + H2). Model-derived from reverse-rate data and thermochemistry "
                  "within the evaluation; tagged recommended as an evaluated value."),
    })

    # ---------------- Channel 5 ----------------
    T.append({
        "table_id": "rec-ii-smith1978-h3oplus-h2o3-mn-300k",
        "source_id": "rec-ii-smith-1978-dtic-ada056649",
        "locator": ("Smith, Church and Miller, J. Chem. Phys. 68, 1224 (1978), Table II (last four rows), reprinted in "
                    "DTIC ADA056649, PDF page 16 (report page 10)"),
        "process": "H3O+(H2O)3 + X- -> neutrals (binary mutual neutralization), X- = Cl-, NO3-, NO3-.HNO3",
        "quantity": "rate_coefficient_m3_s",
        "x_name": "temperature_K",
        "x": [300, 300, 300],
        "y": [r(4.8e-8 * CM3), r(5.5e-8 * CM3), r(5.7e-8 * CM3)],
        "y_labels": ["H3O+(H2O)3 + Cl- (reference: present work)",
                     "H3O+(H2O)3 + NO3- (reference 1: Smith and Church 1976)",
                     "H3O+(H2O)3 + NO3-.HNO3 (reference 1: Smith and Church 1976)"],
        "y_uncertainty": [r(0.6e-8 * CM3), r(1.0e-8 * CM3), r(1.0e-8 * CM3)],
        "original_units": "alpha x 1e8 cm3 s-1",
        "conversion_to_SI": "value x 1e-8 cm3 s-1 multiplied by 1e-6 gives m3 s-1",
        "tag": "measured",
        "extraction_method": "manual transcription from rendered page image of the scanned reprint (text layer is OCR-garbled)",
        "notes": ("Flowing afterglow (helium carrier), 300 K, Langmuir-probe/mass-spectrometer analysis. The anions are NOT water "
                  "fragments (no OH-, O-, H- data in this source); these values are the closest measured hydrated-hydronium "
                  "mutual-neutralization coefficients found and are recorded as analog data only, not as H3O+ + OH- values. "
                  "Average cross sections printed in the same table: 0.93, 1.27, 1.54 (x1e-12 cm2) for the three rows. "
                  "The paper notes that all molecular and cluster-ion pairs studied fall within (4-10)e-8 cm3 s-1 and that "
                  "absorbing-sphere theory predicts roughly a factor 2 larger alpha for low-electron-affinity anions (e.g. O-, O2-); "
                  "no measurement for such anions is reported."),
    })
    T.append({
        "table_id": "rec-ii-bogot2025-h3oplus-ohminus-product-ratios",
        "source_id": "rec-ii-bogot-2025-h3oplus-ohminus-mn",
        "locator": ("Research Square preprint rs-4777257/v1 of Bogot et al., Nat. Chem. 17, 541 (2025): main text p. 4 "
                    "(ratios quoted from Bogot et al., Science 383, 285 (2024), ref. 13) and supplementary materials "
                    "section 'Channel assignment by PCA analysis' / Fig. S4 (re-analysis confirming the 3:1 ratio) and "
                    "two-body analysis / Fig. S2 (2H2O share)"),
        "process": "H3O+ + OH- mutual neutralization (merged beams, DESIREE, collision energy below 0.1 eV): neutral product channels",
        "quantity": "product_channel_ratio_approximate",
        "x_name": "energy_eV_cm",
        "x": [None, None],
        "y": [3.0, 0.02],
        "y_labels": ["ratio (H2O + OH + H) : (2 OH + H2), stated as '~3:1'",
                     "fraction of detected MN events in 2 H2O (proton-transfer) channel, stated as '~2%' (main text) / 'over ~2%' (SM)"],
        "y_uncertainty": None,
        "original_units": "dimensionless",
        "conversion_to_SI": "none",
        "tag": "measured",
        "extraction_method": "manual transcription of approximate values stated in text; no figure digitization",
        "notes": ("Values are approximate as printed ('~'); no uncertainty or energy dependence given; x is null because the "
                  "collision energy is stated only as 'sub 100 meV'. No absolute cross section or rate coefficient is given. "
                  "Deuterated analogue (same preprint): D vs D2 channel ratio ~2:1 and 2D2O < ~1%; the authors attribute the "
                  "difference to suppression of a ~30% proton-transfer contribution to the H2O + OH + H channel. KER of the "
                  "3-body channels is reported as predominantly released as kinetic energy (distributions only in figures; "
                  "not digitized). Preprint (CC BY 4.0) is retained; the Nature Chemistry version of record and the Science 2024 "
                  "paper were not retrieved (paywalled)."),
    })
    return T


def main():
    out = sys.argv[1]
    for t in tables():
        write(out, t)


if __name__ == "__main__":
    main()
