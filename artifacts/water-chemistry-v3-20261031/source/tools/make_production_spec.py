#!/usr/bin/env python3
"""Write data/tables/chemistry-v3-spec.json, the production chemistry v3 spec.

Every choice below is documented in SCOPE.md; numbers come from retained source
files (LXCat downloads, the frozen water datasets, extracted tables) through the
builder's source types, never typed in here except where marked as a citation
constant (orbital constants of the NIST BEB database, Tan/Montenegro orbital
branching weights, Lennon 1988 fit parameters copied from the extracted table).

usage: make_production_spec.py OUT.json [ion_groups.json]
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_tables import parse_lxcat  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MIT = 'data/sources/lxcat-muroranit-h2o-cs-20261031/Cross section.txt'
RK = '/home/treefor/projects/runtimes/h2o-collision-data-20261029/data/collisions/h2o-gas-20261029.json'
RK_EL = '/home/treefor/projects/runtimes/h2o-collision-data-20261029/data/collisions/h2o-gas-20261029-elastic.json'
K = {"bohr_m": 5.2918e-11, "rydberg_eV": 13.6057, "mc2_eV": 510998.95069}
# NIST BEB orbital constants of H2O (B, U in eV; N; Q = 1), physics.nist.gov Ionization database,
# as bound in the frozen datasets (channel H2O:nist_BEB_valence).
ORB = {"1b1": (12.61, 61.91), "3a1": (15.57, 59.52), "1b2": (19.83, 48.36), "2a1": (36.88, 70.71)}


def orb(label, w):
    B, U = ORB[label]
    return {"label": label, "B": B, "U": U, "N": 2, "w": w}


def mit(keyword, process, occurrence=0):
    return {"type": "lxcat", "file": MIT, "keyword": keyword, "process": process, "occurrence": occurrence}


def rk(channel, f=RK):
    return {"type": "raikou", "file": f, "channel": channel}


def rbeb_total(lo, hi):
    return {"type": "beb", "orbitals": [orb(l, 1.0) for l in ORB], **K, "relativistic_from_eV": 0.0, "lo": lo, "hi": hi, "per_decade": 96,
            "locator": "RBEB (Kim, Santos & Parente, PRA 62, 052710 (2000)) total over the four NIST valence orbitals of H2O; energy dependence only"}


def high_energy_tail(base_at_1keV_src, mid_shape, mid_hi):
    """[1 keV, mid_hi]: base scaled with the shape of mid_shape; [mid_hi, 1e7]: RBEB shape; both continuous."""
    mid = {"type": "scaled_shape", "shape": mid_shape, "anchor": base_at_1keV_src, "E0": 1000.0, "lo": 1000.0, "hi": mid_hi}
    top = {"type": "scaled_shape", "shape": rbeb_total(mid_hi, 1e7), "anchor": mid, "E0": mid_hi, "lo": mid_hi, "hi": 1e7}
    return mid, top


def inelastic_bound():
    # any single inelastic process <= the total inelastic cross section (Munoz 2007 Table II recommended,
    # 7.5 eV - 10 keV; held constant above, an upper bound of a decreasing function)
    return {"source": {"type": "extend_constant", "base": rk("H2O:munoz2007:recommended_inelastic"), "hi": 1e7},
            "desc": "Munoz 2007 recommended total inelastic (held constant above 10 keV)"}


def first_positive(blocks, keyword, process):
    b = next(b for b in blocks if b['keyword'] == keyword and b['process'] == process)
    return next((E, v) for E, v in b['table'] if v > 0)


def inelastic_bound_low():
    # total inelastic (Munoz 2007, from 7.5 eV) as the bound of any unresolved inelastic interval
    return {"source": {"type": "extend_constant", "base": rk("H2O:munoz2007:recommended_inelastic"), "hi": 1e7},
            "desc": "Munoz 2007 recommended total inelastic (held constant above 10 keV)"}


def onset_bound(E1, v1, what):
    # monotone onset: the first positive tabulated value bounds the unresolved interval below it;
    # above the data the total inelastic cross section bounds the process
    return {"source": {"type": "piecewise_bound", "parts": [
                {"type": "inline", "E": [0.0, E1], "sigma": [v1, v1], "locator": "%s: first positive tabulated value %g m2 at %g eV" % (what, v1, E1)},
                {"type": "extend_constant", "base": rk("H2O:munoz2007:recommended_inelastic"), "hi": 1e7}]},
            "desc": "%s onset bound %g m2 below %g eV; total inelastic above" % (what, v1, E1)}


def electron_h2o():
    procs = []
    blocks = parse_lxcat(os.path.join(ROOT, MIT))
    # ---- elastic (momentum transfer preserved by the screened angular law)
    el_mid, _ = high_energy_tail(rk("H2O:elastic_momentum"), {"type": "iam_elastic", "file": RK_EL, "quantity": "momentum"}, 1e4)
    procs.append({
        "id": "e + H2O -> e + H2O (elastic)", "kind": "ELASTIC", "angular": "SCREENED", "sigma_is_momentum_transfer": True,
        "tag": "model_derived", "products": [["e", "SURVIVOR"], ["H2O", "GAS"]], "known_zero_below_eV": 0.0,
        "sigma": {"type": "piecewise", "segments": [
            {"source": mit("ELASTIC", "E + H2O, Elastic"), "hi": 100.0, "class": "MODEL", "tag": "model_derived"},
            {"source": rk("H2O:elastic_momentum"), "lo": 100.0, "hi": 1000.0, "class": "DATA", "tag": "recommended"},
            {"source": el_mid, "class": "MODEL", "tag": "model_derived"}]},
        "unresolved_bound": {"source": {"type": "extend_constant", "base": {"type": "iam_elastic", "file": RK_EL, "quantity": "integral"}, "hi": 1e7},
                             "desc": "IAM integral elastic at 10 keV held constant (decreasing function)"},
        "notes": "Momentum-transfer cross section of the MuroranIT complete set (Kawaguchi et al., Jpn. J. Appl. Phys. 55, 07LD03 (2016)) to 100 eV "
                 "(swarm-validated; above ~150 eV it departs from the recommended values, x5.6 at 1 keV); 100 eV - 1 keV: Song 2021 Table 4 "
                 "(recommended; MuroranIT/Song = 0.98 at 100 eV); 1-10 keV: IAM (NIST SRD 64) momentum-transfer energy dependence normalized to "
                 "Song at 1 keV. Sampled as an integral "
                 "elastic cross section sigma_m/(1-<mu>) with the screened law xi(E) from sigma_m/sigma_el ratios (Song 2021 Tables 3/4 0.1-100 eV; "
                 "Song Table 4 / Munoz 2007 elastic integral 100-1000 eV; IAM 1-10 keV); isotropic below 0.1 eV."})
    rot_bound = {"source": {"type": "extend_constant", "base": rk("H2O:munoz2007:rotational_integral_300K"), "hi": 1e7},
                 "desc": "Munoz 2007 thermal rotational integral (any single rotational transition is below it), constant above 10 keV"}
    for b in blocks:
        if b['keyword'] != 'EXCITATION':
            continue
        loss = float(b['header'][1].split()[0])
        src = mit("EXCITATION", b['process'])
        state = b['process'].split('->')[1].strip().replace(', Excitation', '')
        if loss < 0.1:
            procs.append({"id": "e + H2O -> %s (rotational, %g eV)" % (state, loss), "kind": "EXCITE", "loss_eV": loss, "sink": "ROT",
                          "angular": "ISO", "tag": "model_derived", "products": [["e", "SURVIVOR"], ["H2O", "GAS"]],
                          "sigma": {"type": "piecewise", "segments": [{"source": src, "class": "MODEL", "tag": "model_derived"}]},
                          "unresolved_bound": rot_bound})
        elif loss < 1.0:
            procs.append({"id": "e + H2O -> %s (vibrational, %g eV)" % (state, loss), "kind": "EXCITE", "loss_eV": loss, "sink": "VIB",
                          "angular": "ISO", "tag": "model_derived", "products": [["e", "SURVIVOR"], ["H2O", "GAS"]],
                          "sigma": {"type": "piecewise", "segments": [{"source": src, "class": "MODEL", "tag": "model_derived"}]},
                          "unresolved_bound": inelastic_bound()})
        else:
            mid, top = high_energy_tail(src, rk("H2O:munoz2007:recommended_inelastic"), 1e4)
            procs.append({"id": "e + H2O -> e + OH + H via %s (%g eV)" % (state, loss), "kind": "DISSOC", "loss_eV": loss,
                          "kinetic_fraction": 1.0, "sink": "FRAGMENT", "angular": "ISO", "tag": "model_derived",
                          "products": [["e", "SURVIVOR"], ["OH", "RESERVOIR"], ["H", "RESERVOIR"]],
                          "sigma": {"type": "piecewise", "segments": [
                              {"source": src, "hi": 1000.0, "class": "MODEL", "tag": "model_derived"},
                              {"source": mid, "class": "MODEL", "tag": "model_derived"},
                              {"source": top, "class": "MODEL", "tag": "model_derived"}]},
                          "unresolved_bound": inelastic_bound(),
                          "notes": "Electronic excitation of the MuroranIT set treated as dissociation to ground OH(X) + H(1s) with the excess "
                                   "(state energy - D0) as fragment kinetic energy; 1-10 keV energy dependence of the Munoz 2007 recommended total "
                                   "inelastic cross section, above 10 keV of the RBEB valence total, both normalized at the junction."})
    # ---- ionization channels: MuroranIT swarm-fitted onsets (model) + Song 2021 Table 12 (recommended,
    # Straub et al. 1998 lineage) to 1 keV + measured gross-ionization energy dependence (Munoz 2007) to 5 keV
    # + RBEB energy dependence above, each normalized at the junction (continuous; fixed 1 keV fractions)
    ion_channels = [
        ("H2O+", "E + H2O -> E + E + H2O+, Ionization", "H2O:ion_yield:H2O+", 13.5, [["H2O+", "BIRTH"]], [orb("1b1", 1.0), orb("3a1", 1.0), orb("1b2", 0.08)], "ION_INTERNAL"),
        ("OH+", "E + H2O -> E + E + OH+ + H, Ionization", "H2O:ion_yield:OH+", 17.5, [["OH+", "BIRTH"], ["H", "RESERVOIR"]], [orb("1b2", 0.70), orb("2a1", 0.14)], "FRAGMENT"),
        ("H+", "E + H2O -> E + E + H+ + OH, Ionization", "H2O:ion_yield:H+", 20.0, [["H+", "BIRTH"], ["OH", "RESERVOIR"]], [orb("1b2", 0.22), orb("2a1", 0.58)], "FRAGMENT"),
        ("O+", "E + H2O -> E + E + O+ + H2, Ionization", "H2O:ion_yield:O+", 25.0, [["O+", "BIRTH"], ["H2", "RESERVOIR"]], [orb("2a1", 0.28)], "FRAGMENT"),
        ("H2+", "E + H2O -> E + E + H_2^{+} + O, Ionization", "H2O:ion_yield:H2+", 30.0, [["H2+", "BIRTH"], ["O", "RESERVOIR"]], [orb("2a1", 1.0)], "FRAGMENT"),
        ("O++", "E + H2O -> E + E + E + O^{2+} + H2, Ionization", "H2O:ion_yield:O_atomic_2plus", 90.0, [["O++", "BIRTH"], ["e", "BIRTH"], ["H2", "RESERVOIR"]], [orb("2a1", 1.0)], "FRAGMENT"),
    ]
    for ion, mproc, song, song_lo, prods, orbs, sink in ion_channels:
        song_src = rk(song)
        mid, top = high_energy_tail(song_src, rk("H2O:munoz2007:measured_ionization"), 5000.0)
        segs = []
        mit_b = next(b for b in blocks if b['keyword'] == 'IONIZATION' and b['process'] == mproc)
        mit_lo = mit_b['table'][0][0]
        if ion == 'OH+':
            pass   # Song starts below MuroranIT
        elif mit_lo < song_lo:
            segs.append({"source": {"type": "restrict", "base": mit("IONIZATION", mproc), "lo": mit_lo, "hi": song_lo}, "class": "MODEL", "tag": "model_derived"})
        segs.append({"source": song_src, "class": "DATA", "tag": "recommended"})
        segs.append({"source": mid, "class": "MODEL", "tag": "model_derived"})
        segs.append({"source": top, "class": "MODEL", "tag": "model_derived"})
        p = {"id": "e + H2O -> 2e + %s" % (' + '.join(n for n, _ in prods)) if ion != 'O++' else "e + H2O -> 3e + O++ + H2",
             "kind": "IONIZE", "tag": "recommended", "sink": sink, "kinetic_fraction": 0.0, "orbitals": orbs,
             "products": [["e", "SURVIVOR"], ["e", "BIRTH"]] + prods, "sigma": {"type": "piecewise", "segments": segs},
             "unresolved_bound": inelastic_bound_low()}
        if ion == 'H2O+':
            p["unresolved_bound"] = {"source": {"type": "inline", "E": [12.6, 12.65], "sigma": [2.8343e-24, 2.8343e-24],
                                                "locator": "MuroranIT H2O+ value at 12.677 eV (first positive point of a monotone onset) as the bound of the 12.6178-12.65 eV sliver"},
                                     "desc": "onset sliver bound 2.8343e-24 m2 (MuroranIT value at 12.677 eV)"}
        procs.append(p)
    # ---- ion-pair formation: the parts of the H+ and OH+ yields below the dissociative-ionization
    # thresholds (H+ + OH: 18.70 eV, OH+ + H: 18.12 eV) can only be ion pairs (H+ + OH-: 16.87 eV, OH+ + H-: 17.37 eV)
    procs.append({"id": "e + H2O -> e + H+ + OH- (ion pair)", "kind": "DISSOC", "loss_eV": "threshold", "kinetic_fraction": 1.0,
                  "sink": "FRAGMENT", "angular": "ISO", "tag": "model_derived",
                  "products": [["e", "SURVIVOR"], ["H+", "BIRTH"], ["OH-", "BIRTH"]],
                  "sigma": {"type": "piecewise", "segments": [{"source": {"type": "restrict", "base": mit("IONIZATION", "E + H2O -> E + E + H+ + OH, Ionization"),
                                                                           "lo": 16.95, "hi": 18.6995}, "class": "MODEL", "tag": "model_derived"}]},
                  "unresolved_bound": inelastic_bound(),
                  "notes": "H+ yield of the MuroranIT set below the H+ + OH threshold (18.70 eV, ATcT) assigned to H+ + OH- (16.87 eV); above 18.70 eV the "
                           "ion-pair share of the H+ yield is not separable (open)."})
    procs.append({"id": "e + H2O -> e + OH+ + H- (ion pair)", "kind": "DISSOC", "loss_eV": "threshold", "kinetic_fraction": 1.0,
                  "sink": "FRAGMENT", "angular": "ISO", "tag": "model_derived",
                  "products": [["e", "SURVIVOR"], ["OH+", "BIRTH"], ["H-", "BIRTH"]],
                  "sigma": {"type": "piecewise", "segments": [{"source": {"type": "restrict", "base": rk("H2O:ion_yield:OH+"), "lo": 17.5, "hi": 18.1195},
                                                               "class": "DATA", "tag": "recommended"}]},
                  "unresolved_bound": inelastic_bound(),
                  "notes": "OH+ yield of Song 2021 Table 12 below the OH+ + H threshold (18.12 eV) assigned to OH+ + H- (17.37 eV)."})
    # ---- dissociative attachment (MuroranIT)
    for proc, prods in (("E + H2O -> H- + OH, Attachment", [["H-", "BIRTH"], ["OH", "RESERVOIR"]]),
                        ("E + H2O -> OH- + H, Attachment", [["OH-", "BIRTH"], ["H", "RESERVOIR"]]),
                        ("E + H2O -> O- + H2, Attachment", [["O-", "BIRTH"], ["H2", "RESERVOIR"]])):
        E1, v1 = first_positive(blocks, "ATTACHMENT", proc)
        procs.append({"id": "e + H2O -> %s" % ' + '.join(n for n, _ in prods), "kind": "CAPTURE", "kinetic_fraction": 1.0, "sink": "FRAGMENT",
                      "tag": "model_derived", "products": prods,
                      "sigma": {"type": "piecewise", "segments": [{"source": mit("ATTACHMENT", proc), "class": "MODEL", "tag": "model_derived"}]},
                      "unresolved_bound": onset_bound(E1, v1, "MuroranIT " + proc)})
    xi = {"segments": [
        {"lo": 0.1, "hi": 100.0, "sigma_m": rk("H2O:elastic_momentum"), "sigma_el": rk("H2O:elastic_integral_song2021")},
        {"lo": 100.0, "hi": 1000.0, "sigma_m": rk("H2O:elastic_momentum"), "sigma_el": rk("H2O:munoz2007:elastic_integral")},
        {"lo": 1000.0, "hi": 10000.0, "sigma_m": {"type": "iam_elastic", "file": RK_EL, "quantity": "momentum"},
         "sigma_el": {"type": "iam_elastic", "file": RK_EL, "quantity": "integral"}}]}
    return {"projectile": "e", "target": "H2O", "processes": procs, "xi": xi}


def electron_h():
    # NIST BEB for H(1s): B = U = 13.6057 eV, N = 1, Q = 0.5668 (frozen dataset H-1s-ionization-20261029)
    o = {"label": "1s", "B": 13.6057, "U": 13.6057, "N": 1, "w": 1.0, "Q": 0.5668}
    beb = {"type": "beb", "orbitals": [o], **K, "lo": 14.0, "hi": 10000.0, "per_decade": 1000,
           "locator": "NIST BEB H(1s) constants (frozen dataset H-1s-ionization-20261029, Kim & Rudd PRA 50, 3954 (1994))"}
    return {"projectile": "e", "target": "H", "processes": [{
        "id": "e + H -> 2e + H+", "kind": "IONIZE", "tag": "calculated", "sink": "ION_INTERNAL", "kinetic_fraction": 0.0,
        "orbitals": [{"label": "1s", "B": 13.6057, "U": 13.6057, "N": 1, "w": 1.0}],
        "products": [["e", "SURVIVOR"], ["e", "BIRTH"], ["H+", "BIRTH"]],
        "sigma": {"type": "piecewise", "segments": [{"source": beb, "class": "DATA", "tag": "calculated"}]},
        "unresolved_bound": {"source": {"type": "extend_constant", "base": beb, "hi": 1e7}, "desc": "BEB at 10 keV held constant"},
        "notes": "Energy sharing uses the Q = 1 BEB SDCS shape of the same orbital (declared approximation)."}]}


def electron_metals():
    groups = []
    lennon = json.load(open(os.path.join(ROOT, 'data/extracted/metal-lennon1988-ionization.json')))
    fits = {t['table_id'].split('-')[2][:2]: t for t in lennon['tables'] if 'fit-parameters' in t['table_id']}
    for el in ('Fe', 'Cr', 'Ni'):
        t = fits[el]
        par = dict(zip(t['x'], t['y']))
        model = {"type": "lennon1988", "I_eV": par['I_eV'], "A": par['A'], "B": [par['B1'], par['B2'], par['B3']], "lo": par['I_eV'],
                 "hi": 10000.0, "per_decade": 1000, "locator": t['locator']}
        groups.append({"projectile": "e", "target": el, "processes": [{
            "id": "e + %s -> 2e + %s+" % (el, el), "kind": "IONIZE", "tag": "recommended", "sink": "ION_INTERNAL", "kinetic_fraction": 0.0,
            "orbitals": [{"label": "valence", "B": par['I_eV'], "U": par['I_eV'], "N": 2, "w": 1.0}],
            "products": [["e", "SURVIVOR"], ["e", "BIRTH"], ["%s+" % el, "BIRTH"]],
            "sigma": {"type": "piecewise", "segments": [{"source": model, "class": "DATA", "tag": "recommended"}]},
            "unresolved_bound": {"source": {"type": "extend_constant", "base": model, "hi": 1e7}, "desc": "fit value at 10 keV held constant"},
            "notes": "Lennon et al. 1988 recommended (McGuire 1977 Born basis), quoted error %s; energy sharing from a one-orbital BEB SDCS "
                     "with B = U = fit threshold (model)." % ('40 %' if el != 'Ni' else '200 %')}]})
    return groups


U_EV = 931494103.72
ALPHA_H2O = {"alpha_A3": 1.501, "locator": "NIST CCCBDB SRD 101 experimental polarizability of H2O (data/sources/ion-nist-cccbdb-h2o-polarizability)"}


def ext(table_id):
    return json.load(open(os.path.join(ROOT, 'data/extracted', table_id + '.json')))


def merged_points(table_ids, cm_to_lab, merge_rel=0.03):
    """Digitized points of one or more tables (E_cm, sigma), converted to lab energy; points whose energies lie
    within merge_rel of the first point of a group are merged by arithmetic mean (scatter of single
    measurements; documented, no other smoothing)."""
    pts = []
    for t in table_ids:
        d = ext(t)
        pts += list(zip(d['x'], d['y']))
    pts.sort()
    out, grp = [], []
    for E, v in pts:
        if grp and E > grp[0][0] * (1 + merge_rel):
            out.append((sum(e for e, _ in grp) / len(grp), sum(x for _, x in grp) / len(grp)))
            grp = []
        grp.append((E, v))
    out.append((sum(e for e, _ in grp) / len(grp), sum(x for _, x in grp) / len(grp)))
    return [E * cm_to_lab for E, _ in out], [v for _, v in out]


def ion_groups(mass):
    def lab(p_, t_):
        return (mass[p_] + mass[t_]) / mass[t_]

    def langevin(p_, q, lo=1e-5, hi=1e3):
        return {"id": "%s + H2O momentum transfer (Langevin polarization model)" % p_, "kind": "ELASTIC", "angular": "ISO", "tag": "model_derived",
                "products": [[p_, "SURVIVOR"], ["H2O", "GAS"]], "known_zero_below_eV": 0.0,
                "sigma": {"type": "piecewise", "segments": [{"source": {"type": "langevin", "projectile": p_, "target": "H2O", "charge": q, **ALPHA_H2O,
                                                                        "lo": lo, "hi": hi}, "class": "MODEL", "tag": "model_derived"}]},
                "notes": "No measured ion mobility or momentum-transfer data for this ion in H2O were obtained. Classical polarization (Langevin) "
                         "momentum transfer sigma_m = 1.105 k_L / v with the measured H2O polarizability; it omits the ion-dipole term and the "
                         "short-range repulsion, so it is a lower estimate at collision energies above ~1 eV (open item)."}

    def const_k(p_, k_m3_s, frac, E_cm_hi, locator):
        return {"type": "rate_k", "projectile": p_, "target": "H2O", "k_m3_s": k_m3_s * frac, "lo": 1e-5, "hi": E_cm_hi * lab(p_, "H2O"), "per_decade": 64,
                "locator": locator + (" x branch %g" % frac if frac != 1 else "")}

    def bound_k(p_, k_m3_s):
        return {"source": {"type": "rate_k", "projectile": p_, "target": "H2O", "k_m3_s": k_m3_s, "lo": 1e-5, "hi": 1e7, "per_decade": 16,
                           "locator": "thermal rate coefficient / v (Langevin-type decreasing form) as the declared bound"},
                "desc": "k_300K / v"}

    groups = []
    # ---- H2O+ + H2O
    k_pt = ext('ion-anicich1993-h2oplus-h2o-k300')['y'][0]
    Ept, spt = merged_points(['ion-lishawa1990-fig5-proton-transfer-3a-sigma-digitized'], lab("H2O+", "H2O"))
    Ecx, scx = merged_points(['ion-lishawa1990-fig7-charge-exchange-4a-sigma-digitized', 'ion-lishawa1990-fig7-charge-exchange-4b-sigma-digitized'],
                             lab("H2O+", "H2O"))
    groups.append({"projectile": "H2O+", "target": "H2O", "processes": [
        {"id": "H2O+ + H2O -> H3O+ + OH (proton transfer)", "kind": "CAPTURE", "kinetic_fraction": 1.0, "sink": "FRAGMENT", "tag": "measured",
         "products": [["H3O+", "BIRTH"], ["OH", "RESERVOIR"]], "known_zero_below_eV": 0.0,
         "sigma": {"type": "piecewise", "segments": [
             {"source": const_k("H2O+", k_pt, 1.0, Ept[0] / lab("H2O+", "H2O"), "Anicich 1993 JPCRD k(300 K) = 1.85e-9 cm3/s held constant to the first beam point"),
              "class": "MODEL", "tag": "model_derived"},
             {"source": {"type": "inline", "E": Ept, "sigma": spt, "locator": "Lishawa et al. JCP 93, 3196 (1990) Fig. 5 reaction 3a (H2O+ + D2-18O), digitized points, E_cm -> E_lab; points within 3% in energy averaged"},
              "class": "DATA", "tag": "measured"}]},
         "unresolved_bound": {"source": {"type": "inline", "E": [Ept[-1], 1e7], "sigma": [max(spt[-3:]), max(spt[-3:])], "locator": "largest of the last three measured points held constant"},
                              "desc": "last measured magnitude held constant"},
         "notes": "Isotopic analogue data (deuterated, 18O target); exothermic by 1.00 eV (ATcT); products from an isotropic compound decay with all excess as kinetic energy (the measured product energies of Lishawa Table II indicate direct, forward-scattered dynamics; not represented)."},
        {"id": "H2O+ + H2O -> H2O(fast) + H2O+ (symmetric charge exchange)", "kind": "CHARGE_X", "angular": "FORWARD", "tag": "measured",
         "products": [["H2O", "RESERVOIR"], ["H2O+", "BIRTH"]], "known_zero_below_eV": 0.0,
         "sigma": {"type": "piecewise", "segments": [{"source": {"type": "inline", "E": Ecx, "sigma": scx,
                                                                 "locator": "Lishawa et al. 1990 Fig. 7 reactions 4a and 4b (no isotope effect stated by the source), digitized points, E_cm -> E_lab; points within 3% averaged"},
                                                      "class": "DATA", "tag": "measured"}]},
         "unresolved_bound": bound_k("H2O+", k_pt)},
        langevin("H2O+", 1)]})
    groups.append({"projectile": "H3O+", "target": "H2O", "processes": [langevin("H3O+", 1)]})
    # ---- fragment ions + H2O (thermal rate coefficients, Anicich 1993; constant k up to E_cm = 1 eV)
    for ion, tid, branches in (("OH+", "ion-anicich1993-ohplus-h2o", [([["H2O+", "BIRTH"], ["OH", "RESERVOIR"]], 0), ([["H3O+", "BIRTH"], ["O", "RESERVOIR"]], 1)]),
                               ("H2+", "ion-anicich1993-h2plus-h2o", [([["H2O+", "BIRTH"], ["H2", "RESERVOIR"]], 0), ([["H3O+", "BIRTH"], ["H", "RESERVOIR"]], 1)]),
                               ("H+", "ion-anicich1993-hplus-h2o", [([["H2O+", "BIRTH"], ["H", "RESERVOIR"]], None)])):
        k = ext(tid + '-k300')['y'][0]
        fr = ext(tid + '-branching300')['y'] if branches[0][1] is not None else [1.0]
        procs = []
        for prods, bi in branches:
            f = fr[bi] if bi is not None else 1.0
            procs.append({"id": "%s + H2O -> %s (thermal)" % (ion, ' + '.join(n for n, _ in prods)), "kind": "CAPTURE", "kinetic_fraction": 1.0,
                          "sink": "FRAGMENT", "tag": "model_derived", "products": prods, "known_zero_below_eV": 0.0,
                          "sigma": {"type": "piecewise", "segments": [{"source": const_k(ion, k, f, 1.0, "Anicich 1993 JPCRD k(300 K) %g cm3/s" % (k * 1e6)),
                                                                       "class": "MODEL", "tag": "model_derived"}]},
                          "unresolved_bound": bound_k(ion, k * f),
                          "notes": "Measured 300 K rate coefficient and branching (recommended) extended as a constant rate coefficient to E_cm = 1 eV (model); energy dependence above is unresolved."})
        if ion == "H+":
            d = ext('ion-wedlund2019-hplus-h2o-electron-capture-sigma10')
            procs.append({"id": "H+ + H2O -> H(fast) + H2O+ (charge transfer, keV)", "kind": "CHARGE_X", "angular": "FORWARD", "tag": "recommended",
                          "products": [["H", "RESERVOIR"], ["H2O+", "BIRTH"]], "known_zero_below_eV": 0.0,
                          "sigma": {"type": "piecewise", "segments": [{"source": {"type": "inline", "E": [x * mass["H+"] / U_EV for x in d['x']], "sigma": d['y'],
                                                                                  "locator": "Simon Wedlund et al. A&A 2019 Table A.2 sigma_10 (recommended fit), eV/u x m(H+)/u -> E_lab"},
                                                                       "class": "DATA", "tag": "recommended"}]},
                          "unresolved_bound": bound_k("H+", k)})
        procs.append(langevin(ion, 1))
        groups.append({"projectile": ion, "target": "H2O", "processes": procs})
    d = ext('ion-heninger1989-oplus-h2o-sigma-vs-ecm')
    k_o = ext('ion-anicich1993-oplus-h2o-k300')['y'][0]
    groups.append({"projectile": "O+", "target": "H2O", "processes": [
        {"id": "O+ + H2O -> H2O+ + O (charge transfer)", "kind": "CAPTURE", "kinetic_fraction": 1.0, "sink": "FRAGMENT", "tag": "measured",
         "products": [["H2O+", "BIRTH"], ["O", "RESERVOIR"]], "known_zero_below_eV": 0.0,
         "sigma": {"type": "piecewise", "segments": [
             {"source": const_k("O+", k_o, 1.0, d['x'][0], "Anicich 1993 k(300 K) 2.6e-9 cm3/s below the first ICR point"), "class": "MODEL", "tag": "model_derived"},
             {"source": {"type": "inline", "E": [x * lab("O+", "H2O") for x in d['x']], "sigma": d['y'],
                         "locator": "Heninger et al. GRL 1989 Table 1 (O+(4S) + H2O, ICR), E_cm -> E_lab"}, "class": "DATA", "tag": "measured"}]},
         "unresolved_bound": {"source": {"type": "inline", "E": [11.6 * lab("O+", "H2O"), 1e7], "sigma": [1.77e-19, 1.77e-19],
                                         "locator": "largest of the last three Heninger points held constant"}, "desc": "last measured magnitude held constant"}},
        langevin("O+", 1)]})
    groups.append({"projectile": "O++", "target": "H2O", "processes": [langevin("O++", 2)]})
    # ---- negative ions + H2O
    d = ext('ion-wedlund2019-hminus-h2o-one-electron-loss-sigmam10')
    nodata = lambda what, prods, proj: {"id": what + " (no data)", "kind": "CAPTURE", "kinetic_fraction": 1.0, "sink": "FRAGMENT", "tag": "model_derived",
                                        "products": prods, "known_zero_below_eV": 0.0,
                                        "sigma": {"type": "piecewise", "segments": [{"source": {"type": "none", "lo": 1e-5, "hi": 1e7, "locator": what}, "class": "UNRESOLVED", "tag": "model_derived"}]},
                                        "notes": "Placeholder: no rate or cross section obtained; never sampled; its exposure is recorded in the coverage ledger."}
    groups.append({"projectile": "H-", "target": "H2O", "processes": [
        {"id": "H- + H2O -> H + e + H2O (collisional detachment)", "kind": "CAPTURE", "kinetic_fraction": 1.0, "sink": "FRAGMENT", "tag": "model_derived",
         "products": [["H", "RESERVOIR"], ["e", "BIRTH"], ["H2O", "GAS"]], "known_zero_below_eV": 0.0,
         "sigma": {"type": "piecewise", "segments": [{"source": {"type": "inline", "E": [x * mass["H-"] / U_EV for x in d['x']], "sigma": d['y'],
                                                                 "locator": "Simon Wedlund et al. 2019 Table B.1 sigma_-10 (low confidence, H2/O2 additive rule scaled to one H2O measurement)"},
                                                      "class": "MODEL", "tag": "model_derived"}]},
         "unresolved_bound": None},
        nodata("H- + H2O -> OH- + H2", [["OH-", "BIRTH"], ["H2", "RESERVOIR"]], "H-"),
        langevin("H-", -1)]})
    groups.append({"projectile": "O-", "target": "H2O", "processes": [
        nodata("O- + H2O -> OH- + OH", [["OH-", "BIRTH"], ["OH", "RESERVOIR"]], "O-"),
        nodata("O- + H2O -> O + e + H2O (detachment)", [["O", "RESERVOIR"], ["e", "BIRTH"], ["H2O", "GAS"]], "O-"),
        langevin("O-", -1)]})
    groups.append({"projectile": "OH-", "target": "H2O", "processes": [
        nodata("OH- + H2O -> OH + e + H2O (detachment)", [["OH", "RESERVOIR"], ["e", "BIRTH"], ["H2O", "GAS"]], "OH-"),
        langevin("OH-", -1)]})
    km = sum(ext('ion-neg-miller2011-hminus-h-ad-k-16meV')['y']) / 3.0
    groups.append({"projectile": "H-", "target": "H", "processes": [
        {"id": "H- + H -> H2 + e (associative detachment)", "kind": "CAPTURE", "kinetic_fraction": 1.0, "sink": "FRAGMENT", "tag": "model_derived",
         "products": [["H2", "RESERVOIR"], ["e", "BIRTH"]], "known_zero_below_eV": 0.0,
         "sigma": {"type": "piecewise", "segments": [{"source": {"type": "rate_k", "projectile": "H-", "target": "H", "k_m3_s": km, "lo": 1e-5,
                                                                 "hi": 0.1 * (mass["H-"] + mass["H"]) / mass["H"],
                                                                 "locator": "Miller et al. 2011 merged beams, mean of the three 16 meV values %.4g m3/s, held constant to E_cm = 0.1 eV" % km},
                                                      "class": "MODEL", "tag": "model_derived"}]},
         "unresolved_bound": None}]})
    # ---- dissociative recombination
    dr = []
    br = ext('rec-dr-rosen2000-h2oplus-branching-e0')['y']
    al = ext('rec-dr-rosen2000-h2oplus-alpha-300k')['y'][0]
    for prods, f in (([["OH", "RESERVOIR"], ["H", "RESERVOIR"]], br[0]), ([["O", "RESERVOIR"], ["H2", "RESERVOIR"]], br[1]),
                     ([["O", "RESERVOIR"], ["H", "RESERVOIR"], ["H", "RESERVOIR"]], br[2])):
        dr.append({"id": "e + H2O+ -> %s (DR)" % ' + '.join(n for n, _ in prods), "kind": "CAPTURE", "kinetic_fraction": 1.0, "sink": "FRAGMENT",
                   "tag": "model_derived", "products": prods, "known_zero_below_eV": 0.0,
                   "sigma": {"type": "piecewise", "segments": [{"source": {"type": "dr_powerlaw", "alpha_m3_s": al, "T_K": 300.0, "fraction": f,
                                                                           "segments": [[0.0, 0.1, 1.24], [0.1, 2.0, 1.55]], "lo": 1e-5, "hi": 2.0,
                                                                           "locator": "Rosen et al. Faraday Discuss. 115, 295 (2000): alpha(300 K) 4.3e-7 cm3/s, sigma ~ E^-1.24 below 0.1 eV and E^-1.55 to 2 eV (text), branching at 0 eV"},
                                                                "class": "MODEL", "tag": "model_derived"}]},
                   "unresolved_bound": None,
                   "notes": "Ground-state products with all exothermicity as kinetic energy (the source quotes ground-product releases, not measured KER); branching held at its 0 eV value; resonances at 4 and 11 eV not represented (unresolved above 2 eV). The printed Eq. 5 product inconsistency of the source (channel b) is kept as the balanced Eq. 4 lineage O + H2."})
    groups.append({"projectile": "e", "target": "H2O+", "processes": dr})
    dr = []
    br = ext('rec-dr-jensen2000-h3oplus-branching-e0')['y']
    al = ext('rec-dr-jensen2000-h3oplus-alpha-300k')['y'][0]
    for prods, f in (([["H2O", "RESERVOIR"], ["H", "RESERVOIR"]], br[0]), ([["OH", "RESERVOIR"], ["H2", "RESERVOIR"]], br[1]),
                     ([["OH", "RESERVOIR"], ["H", "RESERVOIR"], ["H", "RESERVOIR"]], br[2]), ([["O", "RESERVOIR"], ["H2", "RESERVOIR"], ["H", "RESERVOIR"]], br[3])):
        dr.append({"id": "e + H3O+ -> %s (DR)" % ' + '.join(n for n, _ in prods), "kind": "CAPTURE", "kinetic_fraction": 1.0, "sink": "FRAGMENT",
                   "tag": "model_derived", "products": prods, "known_zero_below_eV": 0.0,
                   "sigma": {"type": "piecewise", "segments": [{"source": {"type": "dr_powerlaw", "alpha_m3_s": al, "T_K": 300.0, "fraction": f,
                                                                           "segments": [[0.0, 0.2, 1.15]], "lo": 1e-5, "hi": 0.2,
                                                                           "locator": "Jensen et al. ApJ 543, 764 (2000): alpha(300 K) 4.3e-7 cm3/s, sigma ~ E^-1.15 at low energy (text), branching at 0 eV"},
                                                                "class": "MODEL", "tag": "model_derived"}]},
                   "unresolved_bound": None,
                   "notes": "Other measured alpha(300 K): Neau 2000 7.6e-7, Herd 1990 1.0e-6 cm3/s (spread kept, not reconciled). Structure above 0.2 eV and the 10 eV peak are unresolved."})
    groups.append({"projectile": "e", "target": "H3O+", "processes": dr})
    groups.append({"projectile": "e", "target": "OH+", "processes": [
        nodata("e + OH+ -> O + H (DR; Kalosi 2023 alpha(T) fit retained, cross section not obtained)", [["O", "RESERVOIR"], ["H", "RESERVOIR"]], "e")]})
    groups.append({"projectile": "H3O+", "target": "OH-", "processes": [
        nodata("H3O+ + OH- -> H2O + OH + H (mutual neutralization)", [["H2O", "RESERVOIR"], ["OH", "RESERVOIR"], ["H", "RESERVOIR"]], "H3O+")]})
    return groups


def main(out, ion_path=None):
    common = json.load(open(os.path.join(ROOT, 'tests/specs/common-species.json')))
    nist = json.load(open(os.path.join(ROOT, 'data/extracted/metal-nist-ionization-energies.json')))['tables'][0]
    ie = dict(zip(nist['x'], nist['y']))
    masses = json.load(open(os.path.join(ROOT, 'data/extracted/metal-nist-masses.json')))
    std_w = {}
    for t in masses['tables']:
        if t['table_id'].startswith('metal-nist-isotopes-'):
            el = t['table_id'].split('-')[-1]
            frac = None
            import re
            m = re.search(r'Representative isotopic mole fractions \(same order as x\): \[([^\]]*)\]', t['notes'])
            frac = [float(v) for v in m.group(1).split(',')]
            std_w[el] = sum(a * b for a, b in zip(t['y'], frac)) / sum(frac)
    for el in ('Fe', 'Cr', 'Ni'):
        common['thermo']['atoms'][el] = {"mass_u": std_w[el], "ref_eV": 0.0,
                                         "note": "abundance-weighted NIST isotope masses (natural vapour); reference energy of the neutral atom set to 0"}
        common['species'].append({"name": el, "kind": "neutral", "charge": 0, "atoms": {el: 1}, "ref_eV": 0.0, "ref_source": "neutral atom reference 0"})
        common['species'].append({"name": el + "+", "kind": "cation", "charge": 1, "atoms": {el: 1}, "ref_eV": ie[el + ' I'],
                                  "ref_source": "NIST ASD v5.12 first ionization energy %s eV" % ie[el + ' I']})
    groups = [electron_h2o(), electron_h()] + electron_metals()
    import build_tables as bt
    files = {}
    sp = bt.species_table(common, files)
    mass = {x['name']: x['mass_eV'] for x in sp}
    groups += ion_groups(mass)
    spec = {"id": "pf-chemistry-v3-production", "energy_bottom_eV": 1e-5, "energy_top_eV": 1e7, **common, "groups": groups}
    json.dump(spec, open(out, 'w'), indent=1)
    print('groups', len(groups), 'processes', sum(len(g['processes']) for g in groups))


if __name__ == '__main__':
    main(*sys.argv[1:])
