#!/usr/bin/env python3
"""Build normalized metal (Cr, Fe, Ni) data tables under data/extracted/.

Usage:
    metal_build_extracted.py <runtime_root>

Inputs (read as data only, never executed):
  data/sources/metal-nist-asd-ionization-energies/nist_asd_ie_CrFeNi_I.csv
  data/sources/metal-nist-isotopic-compositions/nist_composition_{Cr,Fe,Ni}.html
  data/sources/metal-nist-codata-constants/codata_allascii.txt
  data/sources/metal-chianti-diparams/{cr,fe,ni}_1.diparams
Manually transcribed constants (visually verified against page images of the
retained scans) are held below with their locators:
  Lennon et al. 1988 JPCRD 17, 1285, Table 1 (journal pp. 1333-1335)
  Lotz, IPP 1/76 (1968), Table 1 (report p. 6) and Tables 3-5
Outputs:
  data/extracted/metal-nist-ionization-energies.json
  data/extracted/metal-nist-masses.json
  data/extracted/metal-lennon1988-ionization.json
  data/extracted/metal-lotz-ipp1-76-ionization.json
  data/extracted/metal-chianti-dere2007-ionization.json
"""
import csv
import hashlib
import html
import json
import math
import os
import re
import sys

CM2_TO_M2 = 1.0e-4
EV_TO_J = 1.602176634e-19  # exact (SI 2019)

# Common evaluation grid for model-derived tables (eV). Points below the
# respective threshold are dropped per table; the threshold itself is added.
GRID_EV = [8, 9, 10, 11, 12, 13, 14, 15, 17, 20, 25, 30, 35, 40, 50, 60, 70, 80,
           100, 120, 150, 200, 250, 300, 400, 500, 700, 1000, 1500, 2000, 3000,
           5000, 7000, 10000]

ELEMENTS = {"Cr": 24, "Fe": 26, "Ni": 28}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(root, path):
    return os.path.relpath(path, root)


def sig(x, n=6):
    """Round to n significant figures for model-derived output."""
    if x == 0:
        return 0.0
    return float(f"{x:.{n}g}")


# ---------------------------------------------------------------- NIST ASD IE
def build_ie(root):
    src = os.path.join(root, "data/sources/metal-nist-asd-ionization-energies/nist_asd_ie_CrFeNi_I.csv")
    rows = []
    with open(src, newline="", encoding="utf-8") as f:
        rd = csv.reader(f)
        header = next(rd)
        for r in rd:
            if not r or not r[0].startswith("="):
                continue
            clean = [c.strip()[2:-1] if c.strip().startswith('="') else c.strip() for c in r]
            rows.append(dict(zip(header, clean)))
    species = [r["Sp. Name"] for r in rows]
    ie = [float(r["Ionization Energy (eV)"]) for r in rows]
    unc = [float(r["Uncertainty (eV)"]) for r in rows]
    detail = [
        {
            "species": r["Sp. Name"],
            "ground_shells": r["Ground Shells (a)"],
            "ground_level": r["Ground Level"],
            "ionized_level": r["Ionized Level"],
            "asd_reference_codes": r["References"],
        }
        for r in rows
    ]
    table = {
        "table_id": "metal-nist-asd-ie-CrFeNi-I",
        "source_id": "metal-nist-asd-ionization-energies",
        "locator": "NIST ASD v5.12 Ionization Energies output (CSV, format=2), rows Cr I, Fe I, Ni I; retrieved file "
        + rel(root, src),
        "process": "e + M(ground) -> M+(ground) + 2e; first ionization energy, M = Cr, Fe, Ni",
        "quantity": "ionization_energy_eV",
        "x_name": "species",
        "x": species,
        "y": ie,
        "y_uncertainty": unc,
        "threshold_eV": None,
        "original_units": "eV",
        "conversion_to_SI": "J = eV * 1.602176634e-19 (exact)",
        "tag": "recommended",
        "extraction_method": "parsed from downloaded NIST ASD CSV; values copied without modification",
        "notes": "Matches requested check values: Fe I 7.9024681(12) eV, Cr I 6.76651(4) eV, Ni I 7.639878(17) eV. "
        "Ground/ionized levels and ASD bibliography codes in 'level_detail'. Raw sha256 " + sha256(src),
        "level_detail": detail,
    }
    return {"description": "First ionization energies of Cr, Fe, Ni from NIST ASD", "tables": [table]}


# -------------------------------------------------------- NIST masses, CODATA
def parse_composition(path):
    t = open(path, encoding="utf-8", errors="replace").read()
    m = re.search(r"<pre>(.*?)</pre>", t, flags=re.S | re.I)
    body = html.unescape(re.sub(r"<[^>]+>", "", m.group(1)))
    recs = []
    for block in re.split(r"\n\s*\n", body):
        d = {}
        for line in block.splitlines():
            if "=" in line:
                k, v = line.split("=", 1)
                d[k.strip()] = v.strip()
        if "Mass Number" in d:
            recs.append(d)
    return recs


def split_unc(s):
    """'55.93493633(49)' -> (value, uncertainty, estimated_flag)."""
    est = "#" in s
    s2 = s.replace("#", "")
    m = re.match(r"^([0-9.]+)\(([0-9]+)\)$", s2)
    if not m:
        return float(s2) if s2 else None, None, est
    val, u = m.group(1), m.group(2)
    dec = len(val.split(".")[1]) if "." in val else 0
    return float(val), int(u) * 10 ** (-dec), est


def parse_codata(path):
    out = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            name = line[:60].strip()
            if name in ("electron mass in u", "atomic mass constant", "atomic mass constant energy equivalent in MeV"):
                val = line[60:85].replace(" ", "").replace("...", "")
                unc = line[85:110].replace(" ", "")
                out[name] = (float(val), None if "exact" in unc else float(unc), line[110:].strip())
    return out


def build_masses(root):
    sdir = os.path.join(root, "data/sources/metal-nist-isotopic-compositions")
    cpath = os.path.join(root, "data/sources/metal-nist-codata-constants/codata_allascii.txt")
    codata = parse_codata(cpath)
    me_u, me_u_unc, _ = codata["electron mass in u"]
    u_kg, u_kg_unc, _ = codata["atomic mass constant"]
    u_mev, _, _ = codata["atomic mass constant energy equivalent in MeV"]
    ie = {t["x"][i].split()[0]: t["y"][i] for t in build_ie(root)["tables"] for i in range(len(t["x"]))}
    tables = []
    for el in ("Cr", "Fe", "Ni"):
        p = os.path.join(sdir, f"nist_composition_{el}.html")
        recs = parse_composition(p)
        stable = [r for r in recs if r.get("Isotopic Composition")]
        a_vals = [int(r["Mass Number"]) for r in stable]
        mass, munc = [], []
        frac, fracu = [], []
        for r in stable:
            v, u, _ = split_unc(r["Relative Atomic Mass"])
            mass.append(v)
            munc.append(u)
            fv, fu, _ = split_unc(r["Isotopic Composition"])
            frac.append(fv)
            fracu.append(fu)
        saw = stable[0]["Standard Atomic Weight"]
        saw_v, saw_u, _ = split_unc(saw)
        tables.append({
            "table_id": f"metal-nist-isotopes-{el}",
            "source_id": "metal-nist-isotopic-compositions",
            "locator": f"NIST Atomic Weights and Isotopic Compositions, ele={el}, ascii2 output; isotopes with a listed "
            "representative isotopic composition; file " + rel(root, p),
            "process": f"atomic mass data for {el}",
            "quantity": "relative_atomic_mass_u",
            "x_name": "mass_number",
            "x": a_vals,
            "y": mass,
            "y_uncertainty": munc,
            "threshold_eV": None,
            "original_units": "u (unified atomic mass unit, Da)",
            "conversion_to_SI": f"kg = u * {u_kg} (CODATA 2022 atomic mass constant)",
            "tag": "recommended",
            "extraction_method": "parsed from downloaded NIST HTML <pre> block; values copied without modification",
            "notes": f"Standard atomic weight {saw} (as listed). Representative isotopic mole fractions "
            f"(same order as x): {frac}, uncertainties {fracu}. NIST landing page states isotope masses from AME2012, "
            "atomic weights from Meija et al. 2013, compositions from Berglund and Wieser 2009 (compiled, not "
            "critically evaluated by NIST). Raw sha256 " + sha256(p),
            "standard_atomic_weight": {"value": saw_v, "uncertainty": saw_u, "as_listed": saw},
            "isotopic_mole_fraction": frac,
            "isotopic_mole_fraction_uncertainty": fracu,
        })
        # model-derived singly charged ion masses (ground state ion)
        x_lab = [f"{el}-{a}" for a in a_vals] + [f"{el} (standard atomic weight)"]
        neutral = mass + [saw_v]
        ie_u = ie[el] / (u_mev * 1.0e6)
        ion_u = [sig(m - me_u + ie_u, 12) for m in neutral]
        tables.append({
            "table_id": f"metal-ion-mass-{el}+",
            "source_id": "metal-nist-isotopic-compositions",
            "locator": f"derived from table metal-nist-isotopes-{el}, CODATA 2022 'electron mass in u' and "
            "'atomic mass constant energy equivalent in MeV' (codata_allascii.txt), NIST ASD first ionization energy",
            "process": f"{el}+ ion mass (singly charged, ground state)",
            "quantity": "ion_mass_u",
            "x_name": "isotope_or_mixture",
            "x": x_lab,
            "y": ion_u,
            "y_uncertainty": None,
            "threshold_eV": None,
            "original_units": "u",
            "conversion_to_SI": f"kg = u * {u_kg}",
            "tag": "model_derived",
            "extraction_method": "computed: m(M+) = m(M) - m_e + IE/c^2, with m_e = "
            f"{me_u} u, IE/c^2 = IE[eV]/({u_mev}e6 eV/u) = {ie_u:.3e} u",
            "notes": "Neutral masses are relative atomic masses of the neutral atoms (NIST). The binding-energy term "
            "is below 1e-8 u and is included only for completeness. Last entry uses the standard atomic weight, "
            "i.e. the natural-abundance mean mass. Values rounded to 12 significant figures.",
        })
    return {"description": "Relative atomic masses (NIST) and derived singly charged ion masses for Cr, Fe, Ni",
            "codata_2022_used": {"electron_mass_u": [me_u, me_u_unc], "atomic_mass_constant_kg": [u_kg, u_kg_unc],
                                 "atomic_mass_constant_MeV": u_mev, "file_sha256": sha256(cpath)},
            "tables": tables}


# ------------------------------------------------------------- Lennon 1988
LENNON = {
    # species: (I_eV, A, [B1..BN], estimated_error_percent, journal page)
    "Cr": (6.8, 1.029, [-0.429, 0.052, 0.074], 40, 1333),
    "Fe": (7.9, 1.142, [-0.920, 1.782, -1.694], 40, 1334),
    "Ni": (7.6, 0.832, [-0.766, 1.258, -0.679], 200, 1335),
}
LENNON_BASIS = {
    "Cr": "Sec. 6.24 (p. 1360): 'The McGuire calculation is recommended in the absence of any experimental data' "
          "(McGuire, Phys. Rev. A 16, 62 (1977), Born approximation). Fig. 134 shows McGuire and the recommended curve.",
    "Fe": "Sec. 6.26 (p. 1361): 'Our recommended cross sections are derived from a fit to the calculation of "
          "McGuire' (Phys. Rev. A 16, 62 (1977)). Fig. 140 shows McGuire (Born approx.) and the recommended curve.",
    "Ni": "No Ni I entry in the Sec. 6 text was found in the scan; Fig. 143 (p. 1324) shows McGuire (Born approx.) "
          "and the recommended curve, i.e. the recommendation follows the McGuire calculation.",
}


def lennon_sigma_cm2(E, I, A, B):
    if E < I:
        return 0.0
    s = A * math.log(E / I) + sum(b * (1.0 - I / E) ** (i + 1) for i, b in enumerate(B))
    return 1.0e-13 * s / (I * E)


def build_lennon(root):
    pdf = os.path.join(root, "data/sources/metal-lennon1988-jpcrd/jpcrd_reprint_1.555809.pdf")
    h = sha256(pdf)
    tables = []
    for el, (I, A, B, err, page) in LENNON.items():
        tables.append({
            "table_id": f"metal-lennon1988-{el}I-fit-parameters",
            "source_id": "metal-lennon1988-jpcrd",
            "locator": f"Lennon et al., J. Phys. Chem. Ref. Data 17, 1285 (1988), Table 1 (continued), p. {page}, row "
            f"'{el} I'; Eq. (1) on p. 1325",
            "process": f"e + {el}(ground) -> {el}+ + 2e (single ionization, recommended)",
            "quantity": "fit_parameters_Lennon1988_eq1",
            "x_name": "parameter",
            "x": ["I_eV", "A"] + [f"B{i+1}" for i in range(len(B))],
            "y": [I, A] + B,
            "y_uncertainty": None,
            "threshold_eV": I,
            "original_units": "I in eV; A and B_i in 1e-13 eV^2 cm^2",
            "conversion_to_SI": "sigma[m^2] = 1e-4 * sigma[cm^2]; sigma[cm^2] = 1e-13/(I*E) * [A ln(E/I) + "
            "sum_i B_i (1 - I/E)^i], E and I in eV",
            "tag": "recommended",
            "extraction_method": "manual transcription from page image of retained scan (OCR text layer cross-checked); "
            "each value visually verified",
            "notes": f"Source-estimated error {err}% (Table 1 last column). Basis: {LENNON_BASIS[el]} "
            "Eq. (1) as printed in the scan shows B_i (1 - I/E) without a visible exponent; the indexed sum is only "
            "meaningful with exponent i. The exponent-i form reproduces the plotted recommended peaks read from the "
            "scan: Cr I about 6.4e-16 cm^2 near 20 eV (Fig. 134), Fe I about 5.0e-16 cm^2 near 25 eV (Fig. 140), "
            "Ni I about 4.4e-16 cm^2 near 25 eV (Fig. 143); for Fe I the no-exponent reading would give about "
            "3.8e-16 cm^2 at 25 eV. "
            f"I = {I} eV is the fit threshold used by the source (NIST ASD: see metal-nist-ionization-energies). "
            f"Plotted range of recommended curve: threshold to 1e4 eV. PDF sha256 {h}",
        })
        xs = [I] + [e for e in GRID_EV if e > I]
        ys = [sig(lennon_sigma_cm2(e, I, A, B) * CM2_TO_M2) for e in xs]
        tables.append({
            "table_id": f"metal-lennon1988-{el}I-sigma-evaluated",
            "source_id": "metal-lennon1988-jpcrd",
            "locator": f"evaluated from table metal-lennon1988-{el}I-fit-parameters (Table 1, p. {page}; Eq. (1), p. 1325)",
            "process": f"e + {el}(ground) -> {el}+ + 2e (single ionization)",
            "quantity": "cross_section_m2",
            "x_name": "energy_eV",
            "x": xs,
            "y": ys,
            "y_uncertainty": None,
            "threshold_eV": I,
            "original_units": "cm^2 (formula output), converted to m^2",
            "conversion_to_SI": "m^2 = cm^2 * 1e-4 (applied)",
            "tag": "model_derived",
            "extraction_method": "analytic evaluation of the source fit, Eq. (1) with exponent i on (1 - I/E), at a "
            "fixed energy grid chosen here; no smoothing or interpolation of data points",
            "notes": f"Underlying data are the McGuire (1977) Born-approximation calculation, not measurements. "
            f"Source-estimated error {err}%. Grid points are evaluation points only, not source energies. "
            "Valid range per source figure: threshold to 1e4 eV. Rounded to 6 significant figures.",
        })
    return {"description": "Lennon et al. (1988) recommended single-ionization cross sections for Cr I, Fe I, Ni I: "
            "fit parameters (transcribed) and evaluated tables (model-derived)", "tables": tables}


# ------------------------------------------------------------------- Lotz
# (subshell, q, a [1e-14 cm^2 eV^2], b, c, P [eV]); c = None where printed '-'
LOTZ = {
    "Cr": {"config": "3d5 4s (listed as '4s 3d5')", "shells": [
        ("4s", 1, 4.0, 0.0, None, 6.77),
        ("3d", 5, 2.7, 0.9, 0.2, 8.25),
        ("3p", 6, 4.0, 0.6, 0.4, 48.0)]},
    "Fe": {"config": "3d6 4s2", "shells": [
        ("4s", 2, 4.0, 0.4, 0.6, 7.87),
        ("3d", 6, 2.6, 0.92, 0.19, 9.0),
        ("3p", 6, 4.0, 0.6, 0.4, 59.0)]},
    "Ni": {"config": "3d8 4s2", "shells": [
        ("4s", 2, 4.0, 0.4, 0.6, 8.68),
        ("3d", 8, 2.4, 0.94, 0.17, 10.0),
        ("3p", 6, 4.0, 0.6, 0.4, 73.0)]},
}


def lotz_sigma_cm2(E, shells):
    s = 0.0
    for _, q, a, b, c, P in shells:
        if E < P:
            continue
        damp = 1.0 if b == 0.0 else (1.0 - b * math.exp(-c * (E / P - 1.0)))
        s += a * 1.0e-14 * q * math.log(E / P) / (E * P) * damp
    return s


def build_lotz(root):
    pdf = os.path.join(root, "data/sources/metal-lotz1969-zphys/IPP_1_76_Lotz1968_Sc_to_Zn.pdf")
    h = sha256(pdf)
    tables = []
    for el, d in LOTZ.items():
        sh = d["shells"]
        x, y = [], []
        for name, q, a, b, c, P in sh:
            for k, v in (("q", q), ("a", a), ("b", b), ("c", c), ("P_eV", P)):
                x.append(f"{name}:{k}")
                y.append(v)
        tables.append({
            "table_id": f"metal-lotz-{el}I-subshell-constants",
            "source_id": "metal-lotz1969-zphys",
            "locator": f"Lotz, IPP 1/76 (Feb. 1968), report p. 6 'Table 1' row '{el} I' (q_i, a_i, b_i, c_i); "
            "P_1 from 'Table 3' (4s), P_2 from 'Table 4' (3d), P_3 from 'Table 5' (3p), column I (neutral atom); "
            "formula on report p. 1",
            "process": f"e + {el}(ground) -> {el}+ + 2e (single ionization, empirical formula)",
            "quantity": "fit_parameters_Lotz",
            "x_name": "subshell:parameter",
            "x": x,
            "y": y,
            "y_uncertainty": None,
            "threshold_eV": sh[0][5],
            "original_units": "a_i in 1e-14 cm^2 eV^2; P_i in eV; q_i, b_i, c_i dimensionless",
            "conversion_to_SI": "sigma[cm^2] = sum_i a_i q_i ln(E/P_i)/(E P_i) {1 - b_i exp[-c_i (E/P_i - 1)]}, "
            "term i included for E >= P_i; sigma[m^2] = 1e-4 sigma[cm^2]",
            "tag": "model_derived",
            "extraction_method": "manual transcription from page images of the retained scanned report (no text "
            "layer); each value visually verified",
            "notes": f"Ground configuration {d['config']}. Report text: 'From scandium to zinc no cross-section curve is "
            "known experimentally'; 3d constants were set equal to those of 2p electrons in the same ionization "
            "stage (extrapolated beyond six electrons), i.e. the constants are reasoned estimates, not fits to data. "
            "Probable error stated as not higher than +40/-30 %. c_1 printed as '-' for Cr I (b_1 = 0, term "
            "undamped) and stored as null. For Ni I, P_1 = 8.68 eV is marked '*': energy to remove one 4s electron "
            "while the other remains 4s (differs from the ionization potential). The report is the IPP version "
            "linked by the MPG repository to Z. Phys. 220, 466 (1969); the journal version was not checked for "
            f"changes. PDF sha256 {h}",
        })
        P1 = sh[0][5]
        xs = [P1] + [e for e in GRID_EV if e > P1]
        ys = [sig(lotz_sigma_cm2(e, sh) * CM2_TO_M2) for e in xs]
        tables.append({
            "table_id": f"metal-lotz-{el}I-sigma-evaluated",
            "source_id": "metal-lotz1969-zphys",
            "locator": f"evaluated from table metal-lotz-{el}I-subshell-constants",
            "process": f"e + {el}(ground) -> {el}+ + 2e (single ionization, empirical formula)",
            "quantity": "cross_section_m2",
            "x_name": "energy_eV",
            "x": xs,
            "y": ys,
            "y_uncertainty": None,
            "threshold_eV": P1,
            "original_units": "cm^2 (formula output), converted to m^2",
            "conversion_to_SI": "m^2 = cm^2 * 1e-4 (applied)",
            "tag": "model_derived",
            "extraction_method": "analytic evaluation of the Lotz formula with the transcribed constants at a fixed "
            "energy grid chosen here (subshells 4s, 3d, 3p)",
            "notes": "Empirical-formula estimate without experimental anchoring for Sc-Zn. Stated probable error "
            "+40/-30 %. Grid points are evaluation points only. Rounded to 6 significant figures.",
        })
    return {"description": "Lotz empirical ionization formula constants for Cr I, Fe I, Ni I (IPP report 1/76) "
            "and evaluated cross sections", "tables": tables}


# --------------------------------------------------------- CHIANTI / Dere 2007
def read_diparams(path):
    with open(path, encoding="utf-8") as f:
        lines = f.read().splitlines()
    z, ion, nspl, nfac, neaev = (int(v) for v in lines[0].split())
    l2 = [float(v) for v in lines[1].split()]
    l3 = [float(v) for v in lines[2].split()]
    comments = [ln for ln in lines[3:] if ln.strip() and ln.strip() != "-1"]
    return {"Z": z, "ion": ion, "nspl": nspl, "nfac": nfac, "neaev": neaev,
            "btf": l2[0], "x": l2[1:], "ev1": l3[0], "y": l3[1:], "comments": comments}


def build_chianti(root):
    sdir = os.path.join(root, "data/sources/metal-chianti-diparams")
    ver = open(os.path.join(sdir, "VERSION"), encoding="utf-8").read().strip()
    tables = []
    for el in ("cr", "fe", "ni"):
        p = os.path.join(sdir, f"{el}_1.diparams")
        d = read_diparams(p)
        assert d["nfac"] == 1 and d["neaev"] == 0 and d["nspl"] == len(d["x"]) == len(d["y"])
        El = el.capitalize()
        tables.append({
            "table_id": f"metal-chianti-{El}I-diparams-knots",
            "source_id": "metal-chianti-diparams",
            "locator": f"CHIANTI {ver} dbase/{el}/{el}_1/{el}_1.diparams lines 1-3 (file header cites Dere 2007, "
            "A&A 466, 771)",
            "process": f"e + {El}(ground) -> {El}+ + 2e (direct single ionization)",
            "quantity": "scaled_cross_section_BT_knots",
            "x_name": "scaled_energy_x",
            "x": d["x"],
            "y": d["y"],
            "y_uncertainty": None,
            "threshold_eV": d["ev1"],
            "original_units": "x dimensionless; y = scaled cross section Sigma in 1e-14 cm^2 eV^2",
            "conversion_to_SI": "u = 1 - f + exp(ln f/(1 - x)), E = u*I; sigma[cm^2] = 1e-14 (ln u + 1) Sigma/(u I^2) "
            f"(ChiantiPy util.descale_bti and Ion.diCross); f = {d['btf']}, I = {d['ev1']} eV; m^2 = 1e-4 cm^2",
            "tag": "recommended",
            "extraction_method": "parsed from downloaded CHIANTI file; values copied without modification",
            "notes": "File comment: '" + " / ".join(c.strip() for c in d["comments"]) + "'. CHIANTI interpolates the "
            "knots with a cubic spline (scipy splrep, s=0) in the scaled variables. nfac=1, no excitation-"
            "autoionization term (neaev=0). Which experimental/theoretical data Dere (2007) fitted for this neutral "
            "atom could not be confirmed because the A&A article was not retrievable. Raw sha256 " + sha256(p),
        })
        f, I = d["btf"], d["ev1"]
        xs, ys = [], []
        for x, y in zip(d["x"], d["y"]):
            if x >= 1.0:
                continue
            u = 1.0 - f + math.exp(math.log(f) / (1.0 - x))
            E = u * I
            s = 1.0e-14 * (math.log(u) + 1.0) * y / (u * I * I)
            xs.append(sig(E, 6))
            ys.append(sig(s * CM2_TO_M2, 6))
        tables.append({
            "table_id": f"metal-chianti-{El}I-sigma-at-knots",
            "source_id": "metal-chianti-diparams",
            "locator": f"exact inverse scaling of table metal-chianti-{El}I-diparams-knots",
            "process": f"e + {El}(ground) -> {El}+ + 2e (direct single ionization)",
            "quantity": "cross_section_m2",
            "x_name": "energy_eV",
            "x": xs,
            "y": ys,
            "y_uncertainty": None,
            "threshold_eV": I,
            "original_units": "1e-14 cm^2 eV^2 (scaled), converted",
            "conversion_to_SI": "applied: sigma[m^2] = 1e-4 * 1e-14 (ln u + 1) Sigma/(u I^2)",
            "tag": "recommended",
            "extraction_method": "closed-form inverse of the CHIANTI/Dere Burgess-Tully-type scaling applied to each "
            "tabulated knot; no interpolation",
            "notes": "Only knots x < 1 map to finite energy; the x = 1 knot (Sigma = %s) is the E -> infinity limit "
            "and is omitted here. Values between knots require the CHIANTI cubic spline in scaled space, which "
            "is not evaluated in this file. Rounded to 6 significant figures." % d["y"][-1],
        })
    return {"description": f"CHIANTI {ver} direct-ionization fits (Dere 2007) for Cr I, Fe I, Ni I",
            "tables": tables}


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    root = os.path.abspath(sys.argv[1])
    out = os.path.join(root, "data/extracted")
    builds = {
        "metal-nist-ionization-energies.json": build_ie,
        "metal-nist-masses.json": build_masses,
        "metal-lennon1988-ionization.json": build_lennon,
        "metal-lotz-ipp1-76-ionization.json": build_lotz,
        "metal-chianti-dere2007-ionization.json": build_chianti,
    }
    for name, fn in builds.items():
        doc = fn(root)
        doc["generated_by"] = "tools/research/metal_build_extracted.py"
        path = os.path.join(out, name)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(doc, f, indent=1, ensure_ascii=False)
            f.write("\n")
        print(f"wrote {rel(root, path)} ({len(doc['tables'])} tables) sha256={sha256(path)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
