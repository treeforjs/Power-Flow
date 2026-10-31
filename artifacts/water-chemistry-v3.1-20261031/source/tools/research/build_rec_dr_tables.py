"""Write dissociative-recombination tables (branching at nominal zero collision energy, thermal rate
coefficients, fit parameters) transcribed from retained primary sources. Each value was read from the rendered
page image and/or the PDF text layer; locators are given per table.

Usage: python3 -I build_rec_dr_tables.py <extracted_dir>
"""
import json
import os
import sys


def write(out, t):
    with open(os.path.join(out, t["table_id"] + ".json"), "w") as f:
        json.dump(t, f, indent=1, ensure_ascii=False)
    print("wrote", t["table_id"])


def branching(table_id, source_id, citation, doi, locator, process, labels, y, u, method, notes, tag="measured"):
    return {"table_id": table_id, "source_id": source_id, "citation": citation, "doi": doi, "locator": locator,
            "process": process, "quantity": "branching_fraction", "x_name": "energy_eV_cm",
            "x": [0.0] * len(y), "y": y, "y_labels": labels, "y_uncertainty": u,
            "original_units": "dimensionless fraction of DR events", "conversion_to_SI": "none (dimensionless)",
            "tag": tag, "extraction_method": method, "notes": notes}


def rate300(table_id, source_id, citation, doi, locator, process, k_cm3, u_cm3, method, notes, tag="measured"):
    return {"table_id": table_id, "source_id": source_id, "citation": citation, "doi": doi, "locator": locator,
            "process": process, "quantity": "rate_coefficient_m3_s", "x_name": "temperature_K", "x": [300.0],
            "y": [k_cm3 * 1e-6], "y_uncertainty": None if u_cm3 is None else [u_cm3 * 1e-6],
            "original_units": "cm3/s (printed %.3g%s)" % (k_cm3, "" if u_cm3 is None else " +/- %.2g" % u_cm3),
            "conversion_to_SI": "alpha[m3/s] = alpha[cm3/s] * 1e-6", "tag": tag,
            "extraction_method": method, "notes": notes}


def main():
    out = sys.argv[1]

    # Rosen et al. 2000, H2O+ (CRYRING)
    s, c, d = ("rec-dr-rosen-2000-h2oplus-cryring",
               "S. Rosen et al., Recombination of simple molecular ions studied in storage ring: dissociative recombination of H2O+, Faraday Discuss. 115, 295-302 (2000)",
               "10.1039/a909314a")
    write(out, branching(
        "rec-dr-rosen2000-h2oplus-branching-e0", s, c, d,
        "Faraday Discuss. 115, p. 299 Eq. (4) (channel definitions) and p. 301 Eq. (5) (values)",
        "H2O+ + e -> neutral products (vibrationally relaxed H2O+, 5 s storage; nominal 0 eV)",
        ["a: OH + H", "b: O + H2 as defined in Eq. (4); Eq. (5) prints 'OH + H2' (atom-imbalanced printing inconsistency, retained, not repaired)", "c: O + H + H"],
        [0.20, 0.09, 0.71], [0.05, 0.04, 0.06],
        "PDF text layer checked against rendered pages (see prior assessment ASSESSMENT.txt of 2026-10-30)",
        "Electron beam kT_perp = 1 meV, kT_par = 0.01 meV. Ground-state-product energy releases quoted by the source "
        "(calculated, not measured KER distributions): a 7.5 eV, b 7.6 eV, c 3.1 eV. The source also quotes ASTRID values "
        "from Jensen et al. 1999 (a 0.30, b 0.13, c 0.57) and flowing-afterglow values for vibrationally excited ions from "
        "Rowe et al. 1988 (a 0.55, b 0.24, c >0.24); those are quotations, not separately verified here."))
    write(out, rate300(
        "rec-dr-rosen2000-h2oplus-alpha-300k", s, c, d, "Faraday Discuss. 115, p. 301, first paragraph",
        "H2O+ + e -> neutral products (total DR), isotropic Maxwellian at 300 K", 4.3e-7, 0.6e-7,
        "PDF text layer checked against rendered page",
        "Obtained by the source by convolving the measured cross section (0-30 eV, derived as <sigma v>/|v_e - v_i|, "
        "not fully deconvolved; ~20% total 1-sigma uncertainty) with an isotropic Maxwellian. The source states "
        "sigma ~ E^-1.24 below 0.1 eV and ~E^-1.55 between 0.1 and 2 eV, with resonances near 4 and 11 eV; the cross "
        "section itself is only in Fig. 4 (not tabulated, not digitized here)."))

    # Vejby-Christensen et al. 1997 (ASTRID)
    s, c, d = ("rec-dr-vejby-christensen-1997-h2oplus-h3oplus",
               "L. Vejby-Christensen et al., Complete branching ratios for the dissociative recombination of H2O+, H3O+, and CH3+, Astrophys. J. 483, 531-540 (1997)",
               "10.1086/304242")
    write(out, branching(
        "rec-dr-vejby1997-h2oplus-branching-e0", s, c, d, "Astrophys. J. 483, p. 540, Eq. (9) (Sec. 3.3.3)",
        "H2O+ + e -> neutral products (fully relaxed H2O+, zero relative energy)",
        ["a: OH + H", "b: O + H2", "c: O + H + H"], [0.22, 0.10, 0.68], [0.10, 0.10, 0.20],
        "PDF text layer checked against rendered page image",
        "Superseded in precision by later storage-ring work (Jensen et al. 1999; Rosen et al. 2000) but retained as an "
        "independent measurement."))
    write(out, branching(
        "rec-dr-vejby1997-h3oplus-branching-e0", s, c, d, "Astrophys. J. 483, p. 539, Eq. (7) (Sec. 3.3.1)",
        "H3O+ + e -> neutral products (zero relative energy)",
        ["a: H2O + H", "b: OH + H2", "c: OH + H + H", "d: O + H2 + H"], [0.33, 0.18, 0.48, 0.01],
        [0.08, 0.07, 0.08, 0.04], "PDF text layer checked against rendered page image",
        "Jensen et al. 2000 (same ASTRID group) states that the earlier work did not subtract background in the 0.6-5 eV "
        "valley of the cross section and reports improved branching ratios (see rec-dr-jensen2000-h3oplus-branching-e0)."))

    # Jensen et al. 2000 (ASTRID)
    s, c, d = ("rec-dr-jensen-2000-h3oplus-astrid",
               "M. J. Jensen et al., Dissociative recombination of H3O+, HD2O+, and D3O+, Astrophys. J. 543, 764-774 (2000)",
               "10.1086/317137")
    write(out, branching(
        "rec-dr-jensen2000-h3oplus-branching-e0", s, c, d, "Astrophys. J. 543, p. 771, Table 1, column 'H3O+ (This Work)'",
        "H3O+ + e -> neutral products (E = 0)",
        ["a: H2O + H", "b: OH + H2", "c: OH + H + H", "d: O + H2 + H"], [0.25, 0.14, 0.60, 0.013],
        [0.01, 0.01, 0.02, 0.005], "PDF text layer checked against rendered page image",
        "Ground-state-product energy releases quoted by the source (Eq. 2): a 6.4 eV, b 5.8 eV, c 1.3 eV, d 1.5 eV. "
        "The source derives f(H) = 1.46 +/- 0.03 and f(OH) = 0.74 +/- 0.02 and states that branching ratios are assumed "
        "not to change drastically below ~0.3 eV (assumption, not measurement). Table 1 also lists D3O+ (0.25/0.15/0.57/0.03), "
        "HD2O+, previous ASTRID H3O+ (0.33/0.18/0.48/0.01) and flowing-afterglow H3O+ (0.05/0.36/0.29/0.30, Williams et al. 1996) "
        "for comparison; the flowing-afterglow set is stated by the source to disagree completely."))
    write(out, rate300(
        "rec-dr-jensen2000-h3oplus-alpha-300k", s, c, d, "Astrophys. J. 543, p. 768, Eq. (7)",
        "H3O+ + e -> neutral products (total DR), Maxwellian electrons at 300 K", 4.3e-7, 0.6e-7,
        "PDF text layer (born-digital PDF)",
        "Obtained by integrating the measured cross section (corrected for finite energy resolution) with a Maxwellian; "
        "below the measured range the cross section was assumed to follow sigma ~ E^-1. The source reports sigma ~ E^-1.15 "
        "at low energy, structure 0.2-1 eV, a valley 0.6-5 eV and a peak near 10 eV; cross sections only in Figs. 3-4 "
        "(not tabulated, not digitized here). Electron beam kT_perp ~ 25 meV (source Sec. 2)."))

    # Neau et al. 2000 (CRYRING), annual-report summary + JCP abstract
    s, c, d = ("rec-dr-neau-2000-h3oplus-cryring",
               "A. Neau et al., Dissociative recombination of D3O+ and H3O+: Absolute cross sections and branching ratios, J. Chem. Phys. 113, 1762-1770 (2000); numbers here from the Manne Siegbahn Laboratory Annual Report 1999 contribution 4.8 (pp. 67-68) by the same authors and from the JCP abstract",
               "10.1063/1.481979")
    write(out, branching(
        "rec-dr-neau-msl1999-h3oplus-branching-e0", s, c, d,
        "Manne Siegbahn Laboratory Annual Report 1999, contribution 4.8, p. 68, Sec. 4.8.3.1 (retained PDF MSL99-67-neau-annual-report.pdf)",
        "H3O+ + e -> neutral products (0 eV)",
        ["a: H2O + H", "b: OH + H2", "c: OH + H + H", "d: O + H2 + H"], [0.18, 0.11, 0.67, 0.04],
        [0.05, 0.05, 0.06, 0.06], "PDF text layer of the retained annual-report contribution",
        "Non-refereed laboratory annual-report summary; the archival JCP paper (not retained, paywalled) states in its "
        "abstract 67-70% three-body OX+X+X, 17-18% water, negligible O, which is consistent with these values. Differences "
        "from the ASTRID values (Jensen 2000: 0.25/0.14/0.60/0.013) are acknowledged by the source abstract; not reconciled here."))
    write(out, rate300(
        "rec-dr-neau2000-h3oplus-alpha-300k", s, c, d, "J. Chem. Phys. 113, 1762 (2000), abstract (Crossref-held copy in the source directory)",
        "H3O+ + e -> neutral products (total DR), thermal 300 K", 7.6e-7, None,
        "transcribed from the Crossref-held abstract; full text not retained",
        "The abstract states thermal rates were deduced from the measured cross sections (0.001-28 eV) for electron "
        "temperatures 50 K to 30000 K; only the 300 K value is available here (the alpha(Te) table/fit in the full text "
        "was not retrieved). No uncertainty in the abstract. Differs from ASTRID (4.3e-7) and FALP (1.0e-6) values; not reconciled."))

    # Herd, Adams & Smith 1990 (FALP)
    s, c, d = ("rec-dr-herd-1990-h3oplus-oh-yield",
               "C. R. Herd, N. G. Adams, D. Smith, OH production in the dissociative recombination of H3O+, HCO2+, and N2OH+: comparison with theory and interstellar implications, Astrophys. J. 349, 388-392 (1990)",
               "10.1086/168322")
    write(out, {"table_id": "rec-dr-herd1990-h3oplus-oh-yield-300k", "source_id": s, "citation": c, "doi": d,
                "locator": "Astrophys. J. 349, p. 390, Table 2, column 'H3O+ + e' (ADS scan page 3); abstract p. 388",
                "process": "H3O+ + e -> OH + (H2 or 2H) yield, flowing afterglow at 300 K",
                "quantity": "branching_fraction", "x_name": "temperature_K", "x": [300.0, 300.0, 300.0, 300.0],
                "y": [0.46, 0.19, 0.65, 0.35],
                "y_labels": ["OH(v''=0) + H2(2H)", "OH(v''>0) + H2(2H)", "total OH yield", "other products"],
                "y_uncertainty": [None, None, 0.15, None],
                "original_units": "percent (46%, 19%, 65%, 35%); total OH (65 +/- 15)% from the abstract",
                "conversion_to_SI": "fraction = percent / 100", "tag": "measured",
                "extraction_method": "manual transcription from ADS page scan (rendered at 250 dpi); abstract from text layer",
                "notes": "FALP with LIF detection of OH; ions thermalised in He carrier at 300 K (internal state distribution "
                         "differs from storage-ring experiments). OH(v''=0) vs OH(v''>0) separated by NO quenching. The source "
                         "does not resolve OH + H2 from OH + 2H."})
    write(out, rate300(
        "rec-dr-herd1990-h3oplus-alpha-300k", s, c, d, "Astrophys. J. 349, p. 390, Table 2, row 'alpha_e' (footnote a: Smith and Adams 1984; Adams and Smith 1988b)",
        "H3O+ + e -> neutral products (total DR), FALP at 300 K", 1.0e-6, None,
        "manual transcription from ADS page scan",
        "Quoted by Herd et al. from earlier FALP work of the same group (footnote a); no uncertainty printed. "
        "Larger than storage-ring thermal values (4.3e-7 ASTRID, 7.6e-7 CRYRING); not reconciled.",
        tag="measured"))

    # Kalosi et al. 2023 (CSR), OH+ fit parameters
    s, c, d = ("rec-dr-kalosi-2023-ohplus-csr",
               "A. Kalosi et al., Dissociative recombination of rotationally cold OH+ and its implications for the cosmic ray ionization rate in diffuse clouds, Astrophys. J. Lett. 955, L26 (2023); arXiv:2309.17015",
               "10.3847/2041-8213/acf71d")
    write(out, {"table_id": "rec-dr-kalosi2023-ohplus-alpha-kinetic-T-fit", "source_id": s, "citation": c, "doi": d,
                "locator": "arXiv:2309.17015v1, Appendix E, Eqs. (E5)-(E7), Tables 1 and 2 (PDF pp. 17-18)",
                "process": "OH+ (X, v=0, N=0/1) + e -> O + H (total DR), kinetic-temperature (Maxwellian) rate coefficient",
                "quantity": "rate_coefficient_m3_s", "x_name": "temperature_K", "x": [], "y": [], "y_uncertainty": None,
                "expression": {
                    "form_E5": "alpha(T)[cm3/s] = A*(300/T)^n + T^-1.5 * sum_{r=1..4} c_r*exp(-T_r/T)",
                    "E5_parameters_central": {"A": 1.10e-7, "n": 0.767, "c1": 3.46e-4, "c2": -9.16e-4, "c3": -1.85e-3, "c4": 3.33e-2,
                                              "T1": 129, "T2": 1220, "T3": 5900, "T4": 43400},
                    "E5_parameters_lower_error_band": {"A": 9.76e-8, "n": 0.655, "c1": 1.88e-4, "c2": -9.18e-4, "c3": -2.31e-3, "c4": -6.83e-3,
                                                       "T1": 114, "T2": 868, "T3": 2980, "T4": 10100},
                    "E5_parameters_upper_error_band": {"A": 1.32e-7, "n": 0.812, "c1": 4.76e-4, "c2": -8.19e-4, "c3": -7.32e-4, "c4": 4.06e-2,
                                                       "T1": 137, "T2": 1290, "T3": 6990, "T4": 37500},
                    "E5_fit_max_relative_deviation": "0.5% below 6000 K, 1.4% from 6000 to 20000 K (source)",
                    "form_E6_E7": "alpha(T)[cm3/s] = A*(T/300)^beta*exp(-gamma/T); F(T) = F0*exp(g*(1/T - 1/300)) (log-normal uncertainty factor)",
                    "E6_E7_parameters_by_interval": [
                        {"T_range_K": [10, 30], "A": 1.58e-7, "beta": -0.595, "gamma": -2.30, "F0": 1.32, "g": 1.32},
                        {"T_range_K": [30, 100], "A": 1.98e-7, "beta": -0.519, "gamma": -0.792, "F0": 1.27, "g": 2.68},
                        {"T_range_K": [100, 1000], "A": 1.78e-7, "beta": -1.099, "gamma": 52.4, "F0": 1.28, "g": 0.764},
                        {"T_range_K": [1000, 4000], "A": 1.02e-7, "beta": -0.845, "gamma": -198, "F0": 1.31, "g": 8.30},
                        {"T_range_K": [4000, 11000], "A": 1.48e-8, "beta": -0.301, "gamma": -2285, "F0": 1.28, "g": 2.26},
                        {"T_range_K": [11000, 20000], "A": 1.98e-9, "beta": 0.142, "gamma": -6885, "F0": 1.40, "g": 29.7}],
                    "validity_K": [10, 20000],
                    "units": "cm3/s; T in K"},
                "original_units": "cm3/s (fit-function parameters)",
                "conversion_to_SI": "alpha[m3/s] = alpha[cm3/s] * 1e-6 (apply after evaluating the source expression)",
                "tag": "measured",
                "extraction_method": "transcribed from Tables 1 and 2 of the arXiv preprint text layer (born-digital); fit parameters only, no grid evaluated",
                "notes": "Derived by the source from merged-beam rate coefficients measured at the Cryogenic Storage Ring for "
                         "OH+ with ~80% N=0 and ~20% N=1; applies to OH+ in its ground electronic, vibrational and rotational "
                         "level and to a Maxwellian electron (kinetic) temperature. x/y are empty because the source gives only "
                         "fit functions; no values were computed here. The text states that below ~1.6 eV collision energy the "
                         "O(1D) + H(n=1) channel dominates and above it O(3P) + H(n=2) opens (citing Stromholm et al. 1997); no "
                         "branching fractions are tabulated. The arXiv version is the preprint (draft of 2023-10-02); the "
                         "published ApJL version was not compared."})


if __name__ == "__main__":
    main()
