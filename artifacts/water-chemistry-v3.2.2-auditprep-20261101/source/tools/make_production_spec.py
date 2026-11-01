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
            d_ = {"id": "e + H2O -> %s (rotational, %g eV)" % (state, loss), "kind": "EXCITE", "loss_eV": loss, "sink": "ROT",
                  "angular": "ISO", "tag": "model_derived", "products": [["e", "SURVIVOR"], ["H2O", "GAS"]],
                  "sigma": {"type": "piecewise", "segments": [{"source": src, "class": "MODEL", "tag": "model_derived"}]},
                  "unresolved_bound": rot_bound}
            t_end = b['table'][-1]
            if t_end[1] == 0.0 and t_end[0] < 100.0:
                # v3.2 (decision of 2026-10-31): MuroranIT ends the J=0->2 and J=0->3 tables at 5 eV; continue them to 1 keV
                # with the J=0->1 energy dependence of the same set scaled at the cut (continuous; model_derived)
                e_cut = max(E for E, v in b['table'] if v > 0)
                j1 = next(x for x in blocks if x['keyword'] == 'EXCITATION' and '(J=1)' in x['process'])
                base = {"type": "restrict", "base": src, "lo": b['table'][0][0], "hi": e_cut}
                tail = {"type": "scaled_shape", "shape": mit("EXCITATION", j1['process']), "anchor": base, "E0": e_cut, "lo": e_cut, "hi": 1000.0}
                d_["sigma"]["segments"] = [{"source": base, "class": "MODEL", "tag": "model_derived"},
                                           {"source": tail, "class": "MODEL", "tag": "model_derived"}]
                d_["documented_bound"] = {"source": rot_bound["source"], "desc": rot_bound["desc"], "energies_eV": [5.0, 10.0, 30.0, 100.0, 300.0, 1000.0],
                                          "rule": "J=0->1 energy dependence of MuroranIT scaled at %g eV" % e_cut,
                                          "uncertainty": "factor 3 above %g eV (shape of a dipole-allowed transition used for a dipole-forbidden one)" % e_cut}
                d_["notes"] = ("MuroranIT (swarm-fitted) to %g eV; above, to 1 keV, the J=0->1 energy dependence of the same set scaled at %g eV "
                               "(v3.2, model_derived, uncertainty factor 3: the dipole-allowed J=1 shape for a dipole-forbidden transition); "
                               "v3/v3.1 left 5 eV - 1 keV unresolved, bounded by the Munoz 2007 thermal rotational integral, which made "
                               "electrons above 5 eV dominate the bounded share of breakdown cells. Above 1 keV unresolved as J=1." % (e_cut, e_cut))
            procs.append(d_)
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
            if ion in ('O+', 'H2+', 'O++'):
                # v3.2 (decision of 2026-10-31): between the thermochemical threshold and the swarm-fitted appearance energy of
                # the MuroranIT table (where its value is zero) the cross section is modelled as zero (model_derived); v3/v3.1
                # left that onset sliver unresolved, bounded by the total inelastic cross section
                segs.append({"source": {"type": "zero", "lo": 10.0, "hi": mit_lo}, "class": "MODEL", "tag": "model_derived"})
            segs.append({"source": {"type": "restrict", "base": mit("IONIZATION", mproc), "lo": mit_lo, "hi": song_lo}, "class": "MODEL", "tag": "model_derived"})
        segs.append({"source": song_src, "class": "DATA", "tag": "recommended"})
        segs.append({"source": mid, "class": "MODEL", "tag": "model_derived"})
        segs.append({"source": top, "class": "MODEL", "tag": "model_derived"})
        p = {"id": "e + H2O -> 2e + %s" % (' + '.join(n for n, _ in prods)) if ion != 'O++' else "e + H2O -> 3e + O++ + H2",
             "kind": "IONIZE", "tag": "recommended", "sink": sink, "kinetic_fraction": 0.0, "orbitals": orbs,
             "products": [["e", "SURVIVOR"], ["e", "BIRTH"]] + prods, "sigma": {"type": "piecewise", "segments": segs},
             "unresolved_bound": inelastic_bound_low()}
        if ion in ('O+', 'H2+', 'O++'):
            p["notes"] = ("v3.2: zero (model_derived) between the thermochemical threshold and the MuroranIT appearance energy %g eV, where the "
                          "swarm-fitted table starts at zero (dissociative ionization appears above its thermochemical threshold); v3/v3.1 left "
                          "this onset sliver unresolved, bounded by the total inelastic cross section, which made hot breakdown cells "
                          "(T_e 20 eV) refuse." % mit_lo)
        if ion == 'H2O+':
            p["unresolved_bound"] = {"source": {"type": "inline", "E": [12.6, 12.65], "sigma": [2.8343e-24, 2.8343e-24],
                                                "locator": "MuroranIT H2O+ value at 12.677 eV (first positive point of a monotone onset) as the bound of the 12.6178-12.65 eV sliver"},
                                     "desc": "onset sliver bound 2.8343e-24 m2 (MuroranIT value at 12.677 eV)"}
        procs.append(p)
    # ---- ion-pair formation: the parts of the H+ and OH+ yields below the dissociative-ionization
    # thresholds (H+ + OH: 18.70 eV, OH+ + H: 18.12 eV) can only be ion pairs (H+ + OH-: 16.87 eV, OH+ + H-: 17.37 eV)
    def ion_pair_tail(win, e_end):
        # v3.2 (decision of 2026-10-31): above the threshold window, the ion-pair cross section continued with the
        # energy dependence of the Munoz 2007 recommended total inelastic cross section scaled at the window end (to 10 keV),
        # then with the RBEB shape (to 10 MeV); continuous, model_derived, uncertainty factor 10
        mid = {"type": "scaled_shape", "shape": rk("H2O:munoz2007:recommended_inelastic"), "anchor": win, "E0": e_end, "lo": e_end, "hi": 1e4}
        top = {"type": "scaled_shape", "shape": rbeb_total(1e4, 1e7), "anchor": mid, "E0": 1e4, "lo": 1e4, "hi": 1e7}
        return [{"source": mid, "class": "MODEL", "tag": "model_derived"}, {"source": top, "class": "MODEL", "tag": "model_derived"}]
    IP_NOTE = (" v3.2: above the window the cross section is continued with the energy dependence of the Munoz 2007 recommended total "
               "inelastic cross section scaled at the window end (to 10 keV), then the RBEB shape (model_derived, uncertainty factor 10; "
               "the Munoz total inelastic bound of v3 is documented). v3/v3.1 left it unresolved above the window, bounded by the total "
               "inelastic cross section, which made the bounded share of breakdown cells exceed 1 % for T_e above a few eV.")
    win_hp = {"type": "restrict", "base": mit("IONIZATION", "E + H2O -> E + E + H+ + OH, Ionization"), "lo": 16.95, "hi": 18.6995}
    win_ohp = {"type": "restrict", "base": rk("H2O:ion_yield:OH+"), "lo": 17.5, "hi": 18.1195}
    ip_doc = {"source": inelastic_bound()["source"], "desc": inelastic_bound()["desc"], "energies_eV": [20.0, 30.0, 100.0, 1000.0, 1e4, 1e6],
              "rule": "Munoz total inelastic energy dependence scaled at the window end; RBEB shape above 10 keV", "uncertainty": "factor 10"}
    procs.append({"id": "e + H2O -> e + H+ + OH- (ion pair)", "kind": "DISSOC", "loss_eV": "threshold", "kinetic_fraction": 1.0,
                  "sink": "FRAGMENT", "angular": "ISO", "tag": "model_derived",
                  "products": [["e", "SURVIVOR"], ["H+", "BIRTH"], ["OH-", "BIRTH"]], "documented_bound": ip_doc,
                  "sigma": {"type": "piecewise", "segments": [{"source": win_hp, "class": "MODEL", "tag": "model_derived"}] + ion_pair_tail(win_hp, 18.6995)},
                  "unresolved_bound": inelastic_bound(),
                  "notes": "H+ yield of the MuroranIT set below the H+ + OH threshold (18.70 eV, ATcT) assigned to H+ + OH- (16.87 eV); above 18.70 eV the "
                           "ion-pair share of the H+ yield is not separable (open)." + IP_NOTE})
    procs.append({"id": "e + H2O -> e + OH+ + H- (ion pair)", "kind": "DISSOC", "loss_eV": "threshold", "kinetic_fraction": 1.0,
                  "sink": "FRAGMENT", "angular": "ISO", "tag": "model_derived",
                  "products": [["e", "SURVIVOR"], ["OH+", "BIRTH"], ["H-", "BIRTH"]], "documented_bound": ip_doc,
                  "sigma": {"type": "piecewise", "segments": [{"source": win_ohp, "class": "DATA", "tag": "recommended"}] + ion_pair_tail(win_ohp, 18.1195)},
                  "unresolved_bound": inelastic_bound(),
                  "notes": "OH+ yield of Song 2021 Table 12 below the OH+ + H threshold (18.12 eV) assigned to OH+ + H- (17.37 eV)." + IP_NOTE})
    # ---- dissociative attachment (MuroranIT)
    for proc, prods in (("E + H2O -> H- + OH, Attachment", [["H-", "BIRTH"], ["OH", "RESERVOIR"]]),
                        ("E + H2O -> OH- + H, Attachment", [["OH-", "BIRTH"], ["H", "RESERVOIR"]]),
                        ("E + H2O -> O- + H2, Attachment", [["O-", "BIRTH"], ["H2", "RESERVOIR"]])):
        E1, v1 = first_positive(blocks, "ATTACHMENT", proc)
        bk = next(x for x in blocks if x['keyword'] == "ATTACHMENT" and x['process'] == proc)
        e_end = bk['table'][-1][0]
        assert bk['table'][-1][1] == 0.0
        procs.append({"id": "e + H2O -> %s" % ' + '.join(n for n, _ in prods), "kind": "CAPTURE", "kinetic_fraction": 1.0, "sink": "FRAGMENT",
                      "tag": "model_derived", "products": prods,
                      "sigma": {"type": "piecewise", "segments": [{"source": mit("ATTACHMENT", proc), "class": "MODEL", "tag": "model_derived"},
                                                                  {"source": {"type": "zero", "lo": e_end, "hi": 1e7}, "class": "MODEL", "tag": "model_derived"}]},
                      "unresolved_bound": onset_bound(E1, v1, "MuroranIT " + proc),
                      "documented_bound": {"source": inelastic_bound()["source"], "desc": inelastic_bound()["desc"], "energies_eV": [20.0, 30.0, 100.0, 1000.0],
                                           "rule": "zero above the end of the resonance region (the source table returns to zero at %g eV)" % e_end,
                                           "uncertainty": "the non-resonant attachment above the resonances is not represented; anion production there is carried by the ion-pair records"},
                      "notes": "MuroranIT dissociative attachment (resonances near 6.5, 8.6 and 11.8 eV) to %g eV, where the table returns to zero. v3.2 "
                               "(decision of 2026-10-31): above it the cross section is modelled as zero (model_derived). Reasoning: "
                               "dissociative attachment to H2O proceeds through the temporary negative-ion resonances of H2O- (the 2B1, 2A1 and 2B2 "
                               "states near 6.5, 8.6 and 11.8 eV); above the last resonance band there is no further resonant capture channel, so "
                               "the resonant attachment cross section vanishes, and anion production at higher energies proceeds through "
                               "ion-pair formation (polar dissociation), carried by the ion-pair records. v3/v3.1 left it unresolved, bounded by the "
                               "total inelastic cross section (now documented as the bound)." % e_end})
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


E_TH = 1.5 * 8.617333262e-5 * 300.0   # 300 K mean CM kinetic energy 3/2 kT (eV)
E_MT_JOIN = 2000.0   # v3.2: lower edge of the ZBL momentum-transfer segment (log-log bridge from the Langevin value at 1 keV)
ZBL_LOCATOR = ("J. F. Ziegler, J. P. Biersack, U. Littmark, The Stopping and Range of Ions in Solids (Pergamon, 1985): universal reduced nuclear "
               "stopping s_n(eps) and reduced energy eps = 32.53 M2 E[keV] / (Z1 Z2 (M1 + M2)(Z1^0.23 + Z2^0.23)); S_n = 8.462e-15 Z1 Z2 M1 s_n / "
               "((M1 + M2)(Z1^0.23 + Z2^0.23)) eV cm2")
MT_DIFF = {"H+": "7.4", "H-": "7.4", "H2+": "5.2", "O++": "2.4"}


# ======================================================================== v3.2: ions and anions + H2O above 1 keV
E_KEV = 1000.0   # v3.2 lower edge: new processes are placeholders below it (v3.1 behaviour kept), data / model / bounded above
HARTREE_EV = 27.211386245988   # CODATA 2018
C_KM_S = 299792.458
# Bondi van der Waals radii (A. Bondi, J. Phys. Chem. 68, 441 (1964)): radius of the heaviest atom of a species
R_BONDI_A = {"H": 1.20, "O": 1.52}
R_T_A = 1.52   # H2O target
RB = {"1b1": 12.61, "3a1": 14.73, "1b2": 18.55, "2a1": 32.2}   # Rudd et al. 1992 Table I(b) (extracted table checked below)
# target-ionization reference of the singly charged dressed ions (O+, OH+, H2O+, H3O+, H2+): "proton" (decision B of 2026-10-31:
# 1 x proton sigma11 at equal speed) or "heplus" (He+ sigma11 at equal speed, Simon Wedlund 2019 Table 3)
SINGLY_CHARGED_IONIZATION_REF = "heplus"   # decision of 2026-10-31
GOB = 'ion-gobet2004-table2-partial-ionization-fit-parameters'
W11 = 'ion-wedlund2019-fitcoeffs-hplus-h2o-sigma11'


def e_of_speed(v_km_s, m_eV):
    b = v_km_s / C_KM_S
    return m_eV * (1.0 / (1.0 - b * b) ** 0.5 - 1.0)


def rorbs(w):
    return [{"label": l, "B": RB[l], "N": 2, "w": w.get(l, 0.0)} for l in ("1b1", "3a1", "1b2", "2a1")]


def rudd_h2o():
    d = ext('ion-rudd1992-table5-sdcs-parameters-h2o')
    return [float(x) for x in d['y']]


def wfit(table, proj, lo, hi, per_decade=64):
    return {"type": "wedlund_fit", "table_id": table, "projectile": proj, "lo": lo, "hi": hi, "per_decade": per_decade}


def v32_ion_extend(groups, mass, refs):
    """v3.2 content above 1 keV for every ion and anion + H2O group (see SCOPE.md section 12).
    Decision of 2026-10-31 (option B): where no measurement exists, the dominant channels are sampled from equal-speed
    model estimates (model_derived; the physics upper bound of each is documented, not used); H3O+ capture and the
    rearrangements stay unresolved with bounds. New processes are placeholders below 1 keV (v3.1 behaviour kept)."""
    import math
    ob_ = ext('ion-rudd1992-table1b-h2o-orbital-binding-energies')
    assert all(abs(dict(zip(ob_['x'], ob_['y']))[k] - v) == 0 for k, v in RB.items())
    RUDD = rudd_h2o()
    a0 = K["bohr_m"]
    I_T = refs["H2O+"] - refs["H2O"]   # first ionization energy of H2O from the ATcT reference energies
    mH = mass["H+"]

    def r_of(P):   # equal-speed energy ratio E_proton / E_P
        return mH / mass[P]

    def species_atoms(sp):
        return {"H+": "H", "H2+": "H", "H-": "H"}.get(sp, "O")

    def R_sum(sp):
        return R_BONDI_A[species_atoms(sp)] + R_T_A

    def q_of(prods, proj):
        return sum(refs[n] for n, _ in prods) - refs[proj] - refs["H2O"]

    # proton reference cross sections (Simon Wedlund 2019 fits; the fits evaluated down to 100 km/s reproduce the source's
    # own appendix extrapolation within 3 %)
    v100 = e_of_speed(100, mH)
    hp_lo, hp_hi = e_of_speed(400, mH), e_of_speed(30000, mH)   # sigma11 validity (Table 4)
    e10_hi = e_of_speed(20000, mH)                               # sigma10 validity upper edge
    sig11_ref = wfit(W11, "H+", hp_lo, hp_hi)
    sig11_ext = wfit(W11, "H+", v100, hp_hi)
    sig10_ext = wfit('ion-wedlund2019-fitcoeffs-hplus-h2o-sigma10', "H+", v100, e10_hi)
    rt_shape = {"type": "rudd_total", "projectile": "H+", "target": "H2O", "orbitals": rorbs({"1b1": 1, "3a1": 1, "1b2": 1, "2a1": 1}), "rudd": RUDD,
                **{k: K[k] for k in ("bohr_m", "rydberg_eV")}, "lo": hp_hi, "hi": 1.1e7, "per_decade": 16,
                "locator": "Rudd et al. 1992 water-vapour SDCS (Table V; inner-shell set for 2a1) integrated over W, four valence orbitals"}
    sig11_hi = {"type": "scaled_shape", "shape": rt_shape, "anchor": sig11_ref, "E0": hp_hi, "lo": hp_hi, "hi": 1.1e7}

    def ob_bound(q):
        Rc = (2.0 * math.sqrt(q) + 1.0) / (I_T / HARTREE_EV) * a0
        return {"type": "constant", "value": math.pi * Rc * Rc, "lo": E_KEV, "hi": 1e7,
                "locator": "classical over-barrier capture radius R_c = (2 sqrt(q) + 1)/I_T in atomic units, q = %d, I_T = %.4f eV (ATcT); "
                           "sigma <= pi R_c^2 (H. Ryufuku, K. Sasaki, T. Watanabe, Phys. Rev. A 21, 745 (1980); A. Niehaus, J. Phys. B 19, 2925 (1986))" % (q, I_T)}

    def ion_bound(proj, q):
        return {"type": "close_distant_bound", "projectile": proj, "R_sum_A": R_sum(proj), "q": q, "ref_mode": "equal_velocity",
                "ref_projectile": "H+", "ref_parts": [sig11_ref], "hold_below": True, "lo": E_KEV, "hi": 1e7, "per_decade": 32}

    efree = {"type": "sum_of", "hold": True, "parts": [rk("H2O:munoz2007:elastic_integral"), rk("H2O:munoz2007:inelastic_integral")]}

    def loss_bound(proj, I_eV):
        return {"type": "close_distant_bound", "projectile": proj, "R_sum_A": R_sum(proj), "q": 1, "ref_mode": "electron",
                "ref_parts": [efree], "hold_below": True, "ref_min_eV": I_eV / 4.0, "lo": E_KEV, "hi": 1e7, "per_decade": 32}

    def bounded(pid, kind, prods, bsrc, desc, notes, **kw):
        d = {"id": pid, "kind": kind, "tag": "model_derived", "products": prods, "placeholder_below_eV": E_KEV,
             "sigma": {"type": "piecewise", "segments": [{"source": {"type": "none", "lo": 1e-5, "hi": 1e7, "locator": pid}, "class": "UNRESOLVED",
                                                          "tag": "model_derived"}]},
             "unresolved_bound": {"source": bsrc, "desc": desc}, "notes": notes}
        d.update(kw)
        return d

    G = {(g["projectile"], g["target"]): g for g in groups}

    def geo_doc(P, rule, unc):
        return {"source": {"type": "constant", "value": math.pi * (R_sum(P) * 1e-10) ** 2, "lo": E_KEV, "hi": 1e7,
                           "locator": "geometric close-collision cross section pi (R_P + R_T)^2, Bondi radii"},
                "desc": "pi (%.2f A)^2 (close collisions, probability <= 1)" % R_sum(P), "rule": rule, "uncertainty": unc}

    def speed_of(E, m):
        return 299792458.0 * math.sqrt(E * (E + 2 * m)) / (E + m)

    def insert_before_mt(g, procs):
        i_mt = next(i for i, p in enumerate(g["processes"]) if p["kind"] == "ELASTIC")
        g["processes"][i_mt:i_mt] = procs

    # ------------------------------------------------ target ionization channels (Gobet 2004 product shares at equal speed)
    CHANS = [("H2O+", [["H2O+", "BIRTH"]], {"1b1": 1.0, "3a1": 1.0, "1b2": 0.08}, "ION_INTERNAL"),
             ("OH+", [["OH+", "BIRTH"], ["H", "RESERVOIR"]], {"1b2": 0.70, "2a1": 0.14}, "FRAGMENT"),
             ("O+", [["O+", "BIRTH"], ["H2", "RESERVOIR"]], {"2a1": 0.28}, "FRAGMENT"),
             ("H+", [["H+", "BIRTH"], ["OH", "RESERVOIR"]], {"1b2": 0.22, "2a1": 0.58}, "FRAGMENT")]
    SHARE_NOTE = ("Product share: Gobet et al. 2004 Table II partial fits for proton impact (events with at least one free electron; O2+ row "
                  "omitted, < 2 %%) as fractions of their sum at the same projectile speed, measured for 20-150 keV protons and held at the "
                  "edge values outside. Ejected-electron energies from the Rudd et al. 1992 SDCS at the projectile speed (water-vapour Table V "
                  "set; universal inner-shell set for 2a1); orbital weights as the electron-impact channel of the same ion; isotropic ejection.")

    def ionization_channels(P, q, total_segments, tag_note, bound=None, extra_notes=""):
        """total_segments: [(lo, hi, base_source_in_P_energy, class, tag)]"""
        r = r_of(P)
        gob = [{"type": "rudd_analytic", "table_id": GOB, "row": rw, "projectile": P, "bohr_m": a0, "lo": 2e4 / r, "hi": 1.5e5 / r} for rw in ("H2O+", "OH+", "O+", "H+")]
        out = []
        for k_, (row, prods, w, sink) in enumerate(CHANS):
            fr = {"type": "branch_fraction", "parts": gob, "index": k_, "lo": 2e4 / r, "hi": 1.5e5 / r, "range_lo": E_KEV, "range_hi": 1e7}
            segs = [{"source": {"type": "product_of", "a": fr, "b": base, "lo": lo_, "hi": hi_, "per_decade": 32}, "class": cls, "tag": tg}
                    for lo_, hi_, base, cls, tg in total_segments]
            d = {"id": "%s + H2O -> %s + e + %s (target ionization)" % (P, P, ' + '.join(n for n, _ in prods)), "kind": "ION_IONIZE",
                 "tag": total_segments[0][4] if P == "H+" else "model_derived", "sink": sink, "kinetic_fraction": 0.0, "orbitals": rorbs(w),
                 "rudd": RUDD, "placeholder_below_eV": E_KEV, "products": [[P, "SURVIVOR"], ["e", "BIRTH"]] + prods,
                 "sigma": {"type": "piecewise", "segments": segs}, "notes": tag_note + " " + SHARE_NOTE + extra_notes + " Below 1 keV: placeholder (v3.1 behaviour)."}
            if bound is not None:
                d["documented_bound"] = dict(bound)
            out.append(d)
        return out

    def ev_segments(P, factor, ref_lo, ref_mid, tail):
        """total q^2 x proton sigma11 at equal speed: [1 keV, E(ref_mid)] from the fit, the Rudd-shape tail above"""
        r = r_of(P)
        e_mid = min(1e7, ref_mid / r)
        segs = [(E_KEV, e_mid, {"type": "equal_velocity", "ref": ref_lo, "ref_projectile": "H+", "projectile": P, "factor": factor,
                                "lo": E_KEV, "hi": e_mid}, "MODEL", "model_derived")]
        if e_mid < 1e7:
            segs.append((e_mid, 1e7, {"type": "equal_velocity", "ref": tail, "ref_projectile": "H+", "projectile": P, "factor": factor,
                                      "lo": e_mid, "hi": 1e7}, "MODEL", "model_derived"))
        return segs

    def ion_bound_doc(P, q, rule=None):
        return {"source": ion_bound(P, q), "desc": "pi (%.2f A)^2 + %d x proton sigma11 at equal speed (close + distant collisions)" % (R_sum(P), q * q),
                "rule": rule or "%d x proton sigma11 (Simon Wedlund 2019 fit) at the same speed" % (q * q),
                "uncertainty": "factor 3 (first-order charge scaling; slow dressed ions ionize also by quasi-molecular electron promotion)"}

    def breakup_model(P, prods, what, note_extra=""):
        I = q_of(prods, P)
        onset = I / 4.0 * mass[P] / mass["e"]
        return {"id": "%s + H2O -> %s (%s)" % (P, ' + '.join(n for n, _ in prods[:-1]) + ' + H2O', what), "kind": "BREAKUP", "tag": "model_derived",
                "angular": "FORWARD", "loss_eV": 0.0, "kinetic_fraction": 1.0, "sink": "NONE", "placeholder_below_eV": E_KEV, "products": prods,
                "sigma": {"type": "piecewise", "segments": [{"source": {"type": "free_collision", "projectile": P, "ref": efree, "I_eV": I, "n": 1.0,
                                                                         "lo": E_KEV, "hi": 1e7, "per_decade": 32}, "class": "MODEL", "tag": "model_derived"}]},
                "documented_bound": {"source": loss_bound(P, I),
                                     "desc": "pi (%.2f A)^2 + free-electron e + H2O (elastic + inelastic) at the equal-speed electron energy above I/4" % R_sum(P),
                                     "rule": "free-collision estimate, I = %.3f eV, onset %.4g eV" % (I, onset),
                                     "uncertainty": "factor 3 above the onset; molecular (close-collision) breakup below the onset is not represented"},
                "notes": "Free-collision estimate (decision B of 2026-10-31): the most weakly bound projectile electron scattered as a free electron of the "
                         "projectile speed by H2O; e + H2O integral elastic + inelastic cross section (Munoz et al. 2007, IAM-SCAR, 1 eV - 10 keV, "
                         "held outside) at E_e = (m_e/m_P) E, times the isotropic share 1 - I/(4 E_e) of collisions transferring at least the "
                         "breakup threshold I = %.3f eV (ATcT reference energies); zero below the onset %.4g eV. Checked on H-: the same estimate "
                         "is within a factor 2.5 above the Simon Wedlund 2019 H- + H2O sigma-10 at 1-130 keV (RESULTS.txt). %sBelow 1 keV: placeholder."
                         % (I, onset, note_extra)}

    def cx_model(P, prods, q=1, bound_q=1):
        r = r_of(P)
        e_hi = min(1e7, e10_hi / r)
        d = {"id": "%s + H2O -> %s (charge transfer above 1 keV)" % (P, ' + '.join(n for n, _ in prods)), "kind": "CHARGE_X", "angular": "FORWARD",
             "tag": "model_derived", "products": prods, "known_zero_below_eV": E_KEV,
             "sigma": {"type": "piecewise", "segments": [{"source": {"type": "equal_velocity", "ref": sig10_ext, "ref_projectile": "H+", "projectile": P,
                                                                     "factor": 1.0, "lo": E_KEV, "hi": e_hi}, "class": "MODEL", "tag": "model_derived"}]},
             "documented_bound": {"source": ob_bound(bound_q), "desc": "classical over-barrier pi R_c^2 (q = %d)" % bound_q,
                                  "rule": "proton sigma10 (Simon Wedlund 2019 fit) at the same speed",
                                  "uncertainty": "factor 3 (near-resonant capture like H+ + H2O; the proton values cover the same speeds)"},
             "notes": "Equal-speed estimate (decision B of 2026-10-31): the recommended proton capture cross section sigma10 of Simon Wedlund 2019 (Table 2 fit, "
                      "100-20000 km/s) at the same speed; the capture is near-resonant like H+ + H2O (ionization energies of the neutral "
                      "products 13.6 eV for O and H, 13.0 eV for OH, 15.4 eV for H2; symmetric for H2O). Forward charge-transfer kinematics. "
                      "Below 1 keV this reaction is the v3.1 record (thermal / measured / bounded gap); this record is zero there by construction."}
        if e_hi < 1e7:
            d["unresolved_bound"] = {"source": {"type": "equal_velocity", "ref": {"type": "powerlaw_from", "ref": sig10_ext, "E0": e10_hi, "p": 0.0,
                                                                                    "lo": e10_hi, "hi": 1.1e7 * r}, "ref_projectile": "H+",
                                                "projectile": P, "factor": 1.0, "lo": e_hi, "hi": 1e7},
                                     "desc": "proton sigma10 at the fit's upper validity edge held (equal speed)"}
        return d

    # ---------------------------------------------------------------- H+ + H2O (data)
    g = G[("H+", "H2O")]
    cx = next(p for p in g["processes"] if p["id"].startswith("H+ + H2O -> H + H2O+"))
    tab = cx["sigma"]["segments"][1]["source"]
    e_tab_hi = tab["E"][-1]
    w10 = wfit('ion-wedlund2019-fitcoeffs-hplus-h2o-sigma10', "H+", e_tab_hi, e10_hi)
    cx["sigma"]["segments"].append({"source": w10, "class": "DATA", "tag": "recommended"})
    kb = dict(cx["unresolved_bound"]["source"], hi=tab["E"][0])
    cx["unresolved_bound"] = {"source": {"type": "piecewise_bound", "parts": [kb, {"type": "powerlaw_from", "ref": w10, "E0": e10_hi, "p": 0.0, "lo": e10_hi, "hi": 1e7}]},
                              "desc": "k_300K / v in the gap below the keV data; sigma10 at the fit's upper validity edge (%.4g eV) held above (decreasing)" % e10_hi}
    cx["notes"] += (" v3.2: above the Table A.2 range (%.4g eV) the recommended polynomial fit of Simon Wedlund 2019 Table 2 (sigma10, high confidence, "
                    "validity 100-20000 km/s) to %.4g eV; above it unresolved, bounded by the fit value at the edge held constant." % (e_tab_hi, e10_hi))
    hsegs = [(E_KEV, 2e4, sig11_ref, "MODEL", "model_derived"), (2e4, 1.5e5, sig11_ref, "DATA", "recommended"),
             (1.5e5, hp_hi, sig11_ref, "MODEL", "model_derived"), (hp_hi, 1e7, dict(sig11_hi, hi=1e7), "MODEL", "model_derived")]
    new = ionization_channels("H+", 1, hsegs, "Total: Simon Wedlund 2019 recommended sigma11 fit (Table 4, high confidence, validity 400-30000 km/s; the "
                              "text gives 400-10000 km/s); above 30000 km/s the Rudd 1992 SDCS integral scaled to the fit at the edge (model). Tagged "
                              "data where the product share is measured (20-150 keV), model elsewhere.")
    e1m1_hi = e_of_speed(7600, mH)
    w1m1 = wfit('ion-wedlund2019-fitcoeffs-hplus-h2o-sigma1m1', "H+", E_KEV, e1m1_hi)
    new.append({"id": "H+ + H2O -> H- + OH+ + H+ (double electron capture)", "kind": "CHARGE_X", "angular": "FORWARD", "tag": "model_derived",
                "loss_eV": 0.0, "kinetic_fraction": 1.0, "sink": "NONE", "placeholder_below_eV": E_KEV,
                "products": [["H-", "BIRTH"], ["OH+", "BIRTH"], ["H+", "BIRTH"]],
                "sigma": {"type": "piecewise", "segments": [{"source": w1m1, "class": "MODEL", "tag": "model_derived"}]},
                "unresolved_bound": {"source": {"type": "powerlaw_from", "ref": w1m1, "E0": e1m1_hi, "p": 0.0, "lo": e1m1_hi, "hi": 1e7},
                                     "desc": "sigma1-1 at the fit's upper validity edge held above (decreasing)"},
                "notes": "Simon Wedlund 2019 Table 2 sigma1-1 fit (low confidence, > 75 %%; validity 100-7600 km/s) from 1 keV; the doubly charged "
                         "target ion is represented by its two-body breakup OH+ + H+ with the fragment kinetic-energy release not represented "
                         "(loss_eV 0). Below 1 keV: placeholder; above the validity edge unresolved, bounded by the edge value held."})
    insert_before_mt(g, new)

    # ---------------------------------------------------------------- heavy positive ions + H2O
    def zero_above_kev(proc):
        proc["sigma"]["segments"].append({"source": {"type": "zero", "lo": E_KEV, "hi": 1e7}, "class": "ZERO", "tag": proc["tag"]})
        proc["notes"] = proc.get("notes", "") + (" v3.2: zero from 1 keV by construction: the reaction is carried there by the record "
                                                 "'... (charge transfer above 1 keV)' with forward charge-transfer kinematics.")

    IONS = {"H2+": ("H2O+ + H2", [["H2", "RESERVOIR"], ["H2O+", "BIRTH"]], [["H+", "BIRTH"], ["H", "RESERVOIR"], ["H2O", "GAS"]], "projectile dissociation"),
            "O+": ("H2O+ + O", [["O", "RESERVOIR"], ["H2O+", "BIRTH"]], [["O++", "BIRTH"], ["e", "BIRTH"], ["H2O", "GAS"]], "projectile electron loss"),
            "OH+": ("H2O+ + OH", [["OH", "RESERVOIR"], ["H2O+", "BIRTH"]], [["O", "RESERVOIR"], ["H+", "BIRTH"], ["H2O", "GAS"]], "projectile dissociation"),
            "H2O+": ("H2O(fast) + H2O+", None, [["OH+", "BIRTH"], ["H", "RESERVOIR"], ["H2O", "GAS"]], "projectile dissociation"),
            "H3O+": (None, None, [["H2O+", "BIRTH"], ["H", "RESERVOIR"], ["H2O", "GAS"]], "projectile dissociation")}
    for ion, (cx_key, cx_prods, bu_prods, bu_what) in IONS.items():
        g = G[(ion, "H2O")]
        add = []
        if ion == "H2O+":
            pt = next(p for p in g["processes"] if p["id"].startswith("H2O+ + H2O -> H3O+ + OH"))
            lt = pt["sigma"]["segments"][2]["source"]
            k_last = lt["sigma"][-1] * speed_of(lt["E"][-1], mass["H2O+"])
            pt["sigma"]["segments"].append({"source": {"type": "rate_k", "projectile": "H2O+", "target": "H2O", "k_m3_s": k_last, "lo": E_KEV, "hi": 1e7, "per_decade": 16,
                                                       "locator": "last measured Lishawa 1990 point (%.4g m2 at %.4g eV) continued as a constant rate (1/v decline)" % (lt["sigma"][-1], lt["E"][-1])},
                                            "class": "MODEL", "tag": "model_derived"})
            pt["documented_bound"] = {"source": dict(pt["unresolved_bound"]["source"], E=[E_KEV, 1e7]), "desc": pt["unresolved_bound"]["desc"],
                                      "rule": "constant rate of the last measured point (sigma ~ 1/v)", "uncertainty": "factor 10 (no keV measurement)"}
            pt["notes"] += (" v3.2 above 1 keV (decision of 2026-10-31): the last measured point continued as a constant rate "
                            "(sigma ~ 1/v, a falling form, model_derived, uncertainty factor 10); the v3.1 bound (last measured magnitude "
                            "held) is documented. Compound kinematics as below 1 keV (keV product angles not known).")
            proc = next(p for p in g["processes"] if ("-> %s" % cx_key) in p["id"])
            m_ = cx_model(ion, [["H2O", "RESERVOIR"], ["H2O+", "BIRTH"]])
            proc["sigma"]["segments"].append(m_["sigma"]["segments"][0])
            proc["documented_bound"] = m_["documented_bound"]
            proc["notes"] += " v3.2 above 1 keV: " + m_["notes"].replace("Below 1 keV this reaction is the v3.1 record (thermal / measured / bounded gap); this record is zero there by construction.", "")
        elif cx_key is not None:
            proc = next(p for p in g["processes"] if ("-> %s" % cx_key) in p["id"])
            zero_above_kev(proc)
            add.append(cx_model(ion, cx_prods))
            if ion in ("OH+", "H2+"):
                pt = next(p for p in g["processes"] if "-> H3O+ +" in p["id"])
                kpt = pt["unresolved_bound"]["source"]["k_m3_s"]
                pt["sigma"]["segments"].append({"source": {"type": "rate_k", "projectile": ion, "target": "H2O", "k_m3_s": kpt, "lo": E_KEV, "hi": 1e7, "per_decade": 16,
                                                           "locator": "Anicich 1993 k(300 K) x branch / v continued above 1 keV (1/v decline)"},
                                                "class": "MODEL", "tag": "model_derived"})
                pt["documented_bound"] = geo_doc(ion, "measured thermal rate coefficient x branch / v (sigma ~ 1/v)", "factor 10 (no keV measurement)")
                pt["notes"] += (" v3.2 above 1 keV (decision of 2026-10-31): the measured thermal rate continued as a constant rate "
                                "(sigma ~ 1/v, model_derived, uncertainty factor 10); the geometric close-collision bound is documented.")
        else:   # H3O+: endothermic capture; energy-defect (Massey) scaled equal-speed proton estimate (decision of 2026-10-31)
            prods = [["H2O", "RESERVOIR"], ["H", "RESERVOIR"], ["H2O+", "BIRTH"]]
            dE = q_of(prods, "H3O+")
            dE_ref = (refs["H+"] - refs["H"]) - (refs["H2O+"] - refs["H2O"])   # IE(H) - IE(H2O): defect of H+ + H2O
            mI = mass["H3O+"]
            ref10 = wfit('ion-wedlund2019-fitcoeffs-hplus-h2o-sigma10', "H3O+", max(E_KEV, e_of_speed(100, mI)), min(1e7, e_of_speed(20000, mI)))
            add.append({"id": "H3O+ + H2O -> H2O + H + H2O+ (charge transfer)", "kind": "BREAKUP", "angular": "FORWARD", "tag": "model_derived",
                        "loss_eV": 0.0, "kinetic_fraction": 1.0, "sink": "NONE", "placeholder_below_eV": E_KEV, "products": prods,
                        "sigma": {"type": "piecewise", "segments": [{"source": {"type": "adiabatic_scaled", "ref": ref10, "projectile": "H3O+", "a_A": 7.0,
                                                                                 "dE_eV": dE, "dE_ref_eV": dE_ref, "lo": E_KEV, "hi": 1e7},
                                                                      "class": "MODEL", "tag": "model_derived"}]},
                        "documented_bound": {"source": ob_bound(1), "desc": "classical over-barrier pi R_c^2 (q = 1)",
                                             "rule": "proton sigma10 at the same speed x exp(-a (|dE| - |dE_ref|)/(h v)), a = 7 A, dE = %.3f eV, dE_ref = %.3f eV" % (dE, abs(dE_ref)),
                                             "uncertainty": "factor 5 above 50 keV, an order of magnitude below 10 keV (adiabatic suppression of an endothermic capture)"},
                        "notes": "Capture by H3O+ (the neutral H3O is unbound: products H2O + H, fast and forward; slow H2O+). No data. Above 1 keV "
                                 "(decision of 2026-10-31): the recommended proton capture cross section of Simon Wedlund 2019 at the same "
                                 "speed (near-resonant reference, energy defect |IE(H) - IE(H2O)| = %.3f eV) times the Massey adiabatic factor "
                                 "exp(-a (|dE| - |dE_ref|)/(h v)) for the defect of this endothermic capture, dE = %.3f eV (ATcT), with Hasted's "
                                 "adiabatic length a = 7 A (H. S. W. Massey, Rep. Prog. Phys. 12, 248 (1949); J. B. Hasted, Proc. R. Soc. A 205, 421 "
                                 "(1951)). Kinematics: projectile-side breakup record (fast H3O -> H2O + H with no extra release; the target leaves "
                                 "as H2O+). Below 1 keV: placeholder." % (abs(dE_ref), dE)})
        if SINGLY_CHARGED_IONIZATION_REF == "heplus":
            mI = mass[ion]
            eh_hi = min(1e7, e_of_speed(10000, mI))
            sh = wfit('ion-wedlund2019-fitcoeffs-heplus-h2o-sigma11', ion, e_of_speed(100, mI), eh_hi)
            isegs = [(E_KEV, eh_hi, sh, "MODEL", "model_derived")]
            if eh_hi < 1e7:
                isegs.append((eh_hi, 1e7, {"type": "scaled_shape", "shape": wfit(W11, ion, eh_hi, 1e7), "anchor": sh, "E0": eh_hi, "lo": eh_hi, "hi": 1e7},
                              "MODEL", "model_derived"))
            add += ionization_channels(ion, 1, isegs,
                                       "Equal-speed estimate: target ionization by He+ without charge change (singly charged dressed analogue), "
                                       "Simon Wedlund 2019 Table 3 sigma11 fit (high confidence, 450-10000 km/s; evaluated down to 100 km/s as the "
                                       "source's appendix extrapolation) at the same speed; above 10000 km/s the proton speed dependence scaled to it. "
                                       "Uncertainty factor 3. Validation: O+ at 6.6 keV gives 1.4e-20 m2 against the measured total positive "
                                       "fragment-ion production (6 +- 3)e-20 m2 of Juhasz et al. 2019.", ion_bound_doc(ion, 1, "He+ sigma11 (Simon Wedlund 2019 fit) at the same speed"))
        else:
            add += ionization_channels(ion, 1, ev_segments(ion, 1.0, sig11_ext, hp_hi, sig11_hi),
                                       "Equal-speed estimate (decision B of 2026-10-31): 1 x the proton sigma11 of Simon Wedlund 2019 (Table 4 fit, evaluated "
                                       "down to 100 km/s as the source's appendix extrapolation; Rudd-shape tail above 30000 km/s) at the same speed "
                                       "(first-order scaling with the net charge q = 1); uncertainty factor 3. Comparison: at 6.6 keV O+ this gives "
                                       "about 6e-22 m2, against the measured total positive fragment-ion production (6 +- 3)e-20 m2 of Juhasz et al. "
                                       "2019 (He+ sigma11 at the same speed: 1.4e-20 m2).", ion_bound_doc(ion, 1))
        add.append(breakup_model(ion, bu_prods, bu_what))
        insert_before_mt(g, add)

    # O++: He2+ + H2O at equal speed (Simon Wedlund 2019 He2+ fits, functions of the speed: evaluated at the O++ speed)
    g = G[("O++", "H2O")]
    mO2 = mass["O++"]
    e21_hi = e_of_speed(5350, mO2)
    s21 = wfit('ion-wedlund2019-fitcoeffs-he2plus-h2o-sigma21', "O++", E_KEV, e21_hi)
    cxo = cx_model("O++", [["O+", "BIRTH"], ["H2O+", "BIRTH"]], bound_q=2)
    cxo.update({"id": "O++ + H2O -> O+ + H2O+ (electron capture)", "known_zero_below_eV": 0.0, "placeholder_below_eV": E_KEV,
                "sigma": {"type": "piecewise", "segments": [
                    {"source": s21, "class": "MODEL", "tag": "model_derived"},
                    {"source": {"type": "scaled_shape", "shape": wfit('ion-wedlund2019-fitcoeffs-hplus-h2o-sigma10', "O++", e21_hi, 1e7), "anchor": s21,
                                "E0": e21_hi, "lo": e21_hi, "hi": 1e7}, "class": "MODEL", "tag": "model_derived"}]},
                "notes": "Equal-speed estimate (decision B of 2026-10-31): one-electron capture by He2+ from H2O, Simon Wedlund 2019 Table 1 sigma21 fit "
                         "(high confidence, 75-5350 km/s), evaluated at the O++ speed (doubly charged analogue; O++ carries a core, He2+ is bare): "
                         "uncertainty factor 3; above the validity edge (%.4g eV) the proton sigma10 speed dependence scaled to it (falling, model). "
                         "Below 1 keV: placeholder (v3.1 behaviour)." % e21_hi})
    cxo["documented_bound"] = dict(cxo["documented_bound"], rule="He2+ sigma21 (Simon Wedlund 2019 fit) at the same speed")
    e22_hi = e_of_speed(10000, mO2)
    s22 = wfit('ion-wedlund2019-fitcoeffs-he2plus-h2o-sigma22', "O++", e_of_speed(100, mO2), e22_hi)
    s22_tail = {"type": "scaled_shape", "shape": wfit(W11, "O++", e22_hi, 1e7), "anchor": s22, "E0": e22_hi, "lo": e22_hi, "hi": 1e7}
    osegs = [(E_KEV, e22_hi, s22, "MODEL", "model_derived"), (e22_hi, 1e7, s22_tail, "MODEL", "model_derived")]
    add = [cxo] + ionization_channels("O++", 2, osegs,
                                      "Equal-speed estimate (decision B of 2026-10-31): target ionization by He2+ without charge change, Simon Wedlund 2019 "
                                      "Table 3 sigma22 fit (high confidence, 1400-10000 km/s; evaluated down to 100 km/s as the source's appendix "
                                      "extrapolation), at the O++ speed; above 10000 km/s the proton sigma11 speed dependence scaled to it at the "
                                      "edge. Uncertainty factor 3.", ion_bound_doc("O++", 2, "He2+ sigma22 (Simon Wedlund 2019 fit) at the same speed"))
    insert_before_mt(g, add)

    # ---------------------------------------------------------------- anions + H2O
    g = G[("H-", "H2O")]
    det = g["processes"][0]
    assert det["kind"] == "BREAKUP"
    t_ = det["sigma"]["segments"][0]["source"]
    e_m10_hi = e_of_speed(20000, mass["H-"])
    wm10 = wfit('ion-wedlund2019-fitcoeffs-hminus-h2o-sigmam10', "H-", t_["E"][-1], e_m10_hi)
    det["sigma"]["segments"].append({"source": wm10, "class": "MODEL", "tag": "model_derived"})
    det["sigma"]["segments"].append({"source": {"type": "powerlaw_from", "ref": wm10, "E0": e_m10_hi, "p": -1.0, "lo": e_m10_hi, "hi": 1e7}, "class": "MODEL", "tag": "model_derived"})
    lb = dict(det["unresolved_bound"]["source"], hi=t_["E"][0])
    det["unresolved_bound"] = {"source": lb, "desc": "Langevin capture rate k_L/v between the threshold and the table start"}
    det["notes"] += (" v3.2: above the Table B.1 range the Simon Wedlund 2019 Table 2 sigma-10 fit (low confidence; validity 100-20000 km/s) to %.4g eV "
                     "(model_derived as the table); above, the edge value continued as v^-2 (first-order high-velocity scaling of electron loss; model)." % e_m10_hi)
    e_lo11, e_hi11 = e_of_speed(650, mass["H-"]), e_of_speed(7600, mass["H-"])
    wm11 = wfit('ion-wedlund2019-fitcoeffs-hminus-h2o-sigmam11', "H-", E_KEV, e_hi11)   # below 650 km/s: the source's appendix extrapolation
    dd = {"id": "H- + H2O -> H+ + e + e + H2O (double detachment)", "kind": "BREAKUP", "angular": "FORWARD", "tag": "model_derived", "loss_eV": 0.0,
          "kinetic_fraction": 1.0, "sink": "NONE", "placeholder_below_eV": E_KEV,
          "products": [["H+", "BIRTH"], ["e", "BIRTH"], ["e", "BIRTH"], ["H2O", "GAS"]],
          "sigma": {"type": "piecewise", "segments": [
              {"source": wm11, "class": "MODEL", "tag": "model_derived"},
              {"source": {"type": "scaled_shape", "shape": dict(wm10, lo=e_hi11, hi=e_m10_hi), "anchor": wm11, "E0": e_hi11, "lo": e_hi11, "hi": e_m10_hi},
               "class": "MODEL", "tag": "model_derived"},
              {"source": {"type": "powerlaw_from", "ref": {"type": "scaled_shape", "shape": dict(wm10, lo=e_hi11, hi=e_m10_hi), "anchor": wm11, "E0": e_hi11,
                                                           "lo": e_hi11, "hi": e_m10_hi}, "E0": e_m10_hi, "p": -1.0, "lo": e_m10_hi, "hi": 1e7},
               "class": "MODEL", "tag": "model_derived"}]},
          "notes": "Simon Wedlund 2019 Table 2 sigma-11 fit (low confidence; validity 650-7600 km/s, %.4g-%.4g eV; additive-rule reconstruction "
                   "anchored to the one H2O measurement of Baribaud et al. 1971 at 5 keV), evaluated from 1 keV (below 650 km/s as the source's "
                   "appendix extrapolation); above 7600 km/s to 20000 km/s the sigma-10 fit's speed dependence scaled to it at the edge (model, uncertainty "
                   "factor 3), and above 20000 km/s that value continued as v^-2. Below 1 keV: placeholder." % (e_lo11, e_hi11)}
    rearr = next(p for p in g["processes"] if p["id"].startswith("H- + H2O -> OH- + H2"))
    k_um = ext('ion-neg-cx-umist-rate22-hminus-h2o-ohminus-k300')['y'][0]
    rearr.update({"id": "H- + H2O -> OH- + H2 (proton transfer; no data below 1 keV)", "placeholder": False, "placeholder_below_eV": E_KEV,
                  "sigma": {"type": "piecewise", "segments": [{"source": {"type": "rate_k", "projectile": "H-", "target": "H2O", "k_m3_s": k_um, "lo": E_KEV, "hi": 1e7,
                                                                          "per_decade": 16, "locator": "UMIST Rate22 k(300 K) = 4.8e-9 cm3/s (measured, Martinez et al. 2010) / v"},
                                                               "class": "MODEL", "tag": "model_derived"}]},
                  "documented_bound": geo_doc("H-", "measured thermal rate coefficient / v (sigma ~ 1/v)", "factor 10 (no keV measurement)"),
                  "notes": "No cross section obtained. Below 1 keV: placeholder (v3.1). Above 1 keV (decision of 2026-10-31): the measured "
                           "thermal rate coefficient over the speed (a falling 1/v form, model_derived, uncertainty factor 10); the geometric "
                           "close-collision bound is documented. Compound kinematics."})
    rearr.pop("unresolved_bound", None)
    insert_before_mt(g, [dd] + ionization_channels("H-", 1, ev_segments("H-", 1.0, sig11_ext, hp_hi, sig11_hi),
                                                  "Equal-speed estimate (decision B of 2026-10-31): 1 x the proton sigma11 of Simon Wedlund 2019 at the same speed "
                                                  "(first-order scaling with |q| = 1); uncertainty factor 3.", ion_bound_doc("H-", 1)))

    for an, neutral, rearr_id in (("O-", "O", "O- + H2O -> OH- + OH"), ("OH-", "OH", None)):
        g = G[(an, "H2O")]
        procs = []
        for p in g["processes"]:
            if "(detachment)" in p["id"]:
                prods = [[neutral, "RESERVOIR"], ["e", "BIRTH"], ["H2O", "GAS"]]
                d = breakup_model(an, prods, "detachment",
                                  "The geometric bound covers total electron loss; Hasted 1956 and the Rahman 1986 compilation were not retrieved. ")
                d["id"] = "%s + H2O -> %s + e + H2O (detachment)" % (an, neutral)
                procs.append(d)
            elif rearr_id and p["id"].startswith(rearr_id):
                p = dict(p)
                p.update({"id": rearr_id + " (no data below 1 keV)", "placeholder": False, "placeholder_below_eV": E_KEV,
                          "sigma": {"type": "piecewise", "segments": [{"source": {"type": "langevin", "projectile": an, "target": "H2O", "charge": -1, **ALPHA_H2O,
                                                                                  "factor": 1.0, "lo": E_KEV, "hi": 1e7}, "class": "MODEL", "tag": "model_derived"}]},
                          "documented_bound": geo_doc(an, "Langevin capture rate k_L / v (capture-type, unit reaction probability)",
                                                      "factor 10 (upper side: unit reaction probability; endothermic by 0.32 eV, ATcT)"),
                          "notes": "No cross section obtained. Below 1 keV: placeholder (v3.1). Above 1 keV (decision of 2026-10-31): a "
                                   "capture-type estimate, the Langevin capture rate over the speed with unit reaction probability (model_derived, "
                                   "upper side, uncertainty factor 10); the geometric close-collision bound is documented. Compound kinematics."})
                p.pop("unresolved_bound", None)
                procs.append(p)
            elif p["kind"] == "ELASTIC":
                procs += ionization_channels(an, 1, ev_segments(an, 1.0, sig11_ext, hp_hi, sig11_hi),
                                             "Equal-speed estimate (decision B of 2026-10-31): 1 x the proton sigma11 of Simon Wedlund 2019 at the same speed "
                                             "(first-order scaling with |q| = 1); uncertainty factor 3.", ion_bound_doc(an, 1))
                procs.append(p)
            else:
                procs.append(p)
        g["processes"] = procs
    return groups

def v32_gap_bridges(groups, mass, refs):
    """v3.2 (decision of 2026-10-31): sampled log-log bridges over the v3.1 gaps between the low-energy data of an ion and
    its v3.2 value at 1 keV (or the Wedlund table start for H+ and H-); model_derived; the v3.1 gap bounds are documented."""
    G = {(g["projectile"], g["target"]): g for g in groups}
    GAP_E = [3.0, 10.0, 30.0, 100.0, 300.0, 999.0]

    def first(g, start):
        return next(p for p in g["processes"] if p["id"].startswith(start))

    def doc(proc, src, desc, rule, unc, energies=GAP_E):
        proc["documented_bound"] = {"source": src, "desc": desc, "rule": rule, "uncertainty": unc, "energies_eV": energies}

    def bridge(a, Ea, b, Eb):
        return {"source": {"type": "loglog_bridge", "a": a, "Ea": Ea, "b": b, "Eb": Eb, "per_decade": 32}, "class": "MODEL", "tag": "model_derived"}

    NOTE = (" v3.2 gap bridge (decision of 2026-10-31): between %.4g and %.4g eV a log-log interpolation from the last low-energy "
            "value to the value at the upper end (model_derived, uncertainty factor %s); v3.1 left the gap unresolved, bounded by "
            "%s (documented, with the ratio bound / sigma).")
    # OH+ and H2+: thermal charge transfer bridged to the keV charge-transfer estimate at 1 keV; proton transfer k f / v throughout
    for ion, cx_key in (("OH+", "H2O+ + OH"), ("H2+", "H2O+ + H2")):
        g = G[(ion, "H2O")]
        th = first(g, "%s + H2O -> %s (thermal)" % (ion, cx_key))
        kev = first(g, "%s + H2O -> " % ion + {"OH+": "OH + H2O+", "H2+": "H2 + H2O+"}[ion] + " (charge transfer above 1 keV)")
        a = th["sigma"]["segments"][0]["source"]
        E1 = a["hi"]
        th["sigma"]["segments"].insert(1, bridge(a, E1, kev["sigma"]["segments"][0]["source"], E_KEV))
        doc(th, th["unresolved_bound"]["source"], th["unresolved_bound"]["desc"], "log-log bridge from k(300 K) x branch / v at E_cm 1 eV to the "
            "1 keV equal-speed proton sigma10", "factor 5")
        th["notes"] += NOTE % (E1, E_KEV, "5", "k_300K / v")
        pt = first(g, "%s + H2O -> H3O+ +" % ion)
        pt["sigma"]["segments"][-1]["source"] = dict(pt["sigma"]["segments"][-1]["source"], lo=pt["sigma"]["segments"][0]["source"]["hi"])
        pt["documented_bound"] = dict(pt["documented_bound"], energies_eV=GAP_E + [1e4, 1e5])
        pt["notes"] += " The thermal rate is continued as a constant rate from E_cm 1 eV (the v3.1 gap) as above 1 keV."
    # H+: thermal constant rate bridged to the start of the Simon Wedlund 2019 table
    g = G[("H+", "H2O")]
    cx = first(g, "H+ + H2O -> H + H2O+")
    a, tab = cx["sigma"]["segments"][0]["source"], cx["sigma"]["segments"][1]["source"]
    cx["sigma"]["segments"].insert(1, bridge(a, a["hi"], tab, tab["E"][0]))
    doc(cx, dict(cx["unresolved_bound"]["source"]["parts"][0]), "k_300K / v (the v3.1 gap bound)",
        "log-log bridge from k(300 K) / v at E_cm 1 eV to the first Simon Wedlund 2019 table value", "factor 3", [3.0, 10.0, 30.0, 50.0])
    cx["notes"] += NOTE % (a["hi"], tab["E"][0], "3", "k_300K / v")
    # H2O+: proton transfer as a constant rate from the last measured point; symmetric charge exchange bridged to 1 keV
    g = G[("H2O+", "H2O")]
    pt = first(g, "H2O+ + H2O -> H3O+ + OH")
    lt = pt["sigma"]["segments"][2]["source"]
    pt["sigma"]["segments"][-1]["source"] = dict(pt["sigma"]["segments"][-1]["source"], lo=lt["E"][-1])
    pt["documented_bound"] = dict(pt["documented_bound"], source=dict(pt["documented_bound"]["source"], E=[lt["E"][-1], 1e7]), energies_eV=GAP_E + [1e4, 1e5])
    pt["notes"] += " The constant-rate continuation starts at the last measured point (%.4g eV), closing the v3.1 gap below 1 keV." % lt["E"][-1]
    cx = first(g, "H2O+ + H2O -> H2O(fast) + H2O+")
    lc = cx["sigma"]["segments"][1]["source"]
    kseg = cx["sigma"]["segments"][-1]["source"]
    cx["sigma"]["segments"].insert(len(cx["sigma"]["segments"]) - 1, bridge(lc, lc["E"][-1], kseg, E_KEV))
    k31 = cx["unresolved_bound"]["source"]
    cx["documented_bound"] = dict(cx["documented_bound"], source={"type": "piecewise_bound", "parts": [dict(k31, hi=E_KEV), cx["documented_bound"]["source"]]},
                                  desc="k_300K / v below 1 keV (the v3.1 gap bound); " + cx["documented_bound"]["desc"] + " above",
                                  energies_eV=[100.0, 300.0, 999.0, 1e4, 5e4, 1e5])
    cx["notes"] += NOTE % (lc["E"][-1], E_KEV, "5", "k_300K / v")
    # O+: Heninger table bridged to the keV charge-transfer estimate at 1 keV
    g = G[("O+", "H2O")]
    th = first(g, "O+ + H2O -> H2O+ + O (charge transfer)")
    hen = th["sigma"]["segments"][1]["source"]
    kev = first(g, "O+ + H2O -> O + H2O+ (charge transfer above 1 keV)")
    th["sigma"]["segments"].insert(2, bridge(hen, hen["E"][-1], kev["sigma"]["segments"][0]["source"], E_KEV))
    doc(th, th["unresolved_bound"]["source"], th["unresolved_bound"]["desc"], "log-log bridge from the last Heninger 1989 point to the 1 keV "
        "equal-speed proton sigma10", "factor 3", [30.0, 100.0, 300.0, 999.0])
    th["notes"] += NOTE % (hen["E"][-1], E_KEV, "3", "the last measured magnitude held")
    # H-: detachment rises from zero at the threshold to the start of the Simon Wedlund 2019 table
    g = G[("H-", "H2O")]
    det = g["processes"][0]
    t_ = det["sigma"]["segments"][0]["source"]
    q = refs["H"] - refs["H-"]
    m, M = mass["H-"], mass["H2O"]
    E_th = q * (2 * (m + M) + q) / (2 * M)
    det["sigma"]["segments"].insert(0, {"source": {"type": "threshold_powerlaw", "ref": t_, "E_th": E_th, "E0": t_["E"][0], "E1": t_["E"][1],
                                                   "lo": E_th, "per_decade": 64}, "class": "MODEL", "tag": "model_derived"})
    doc(det, det["unresolved_bound"]["source"], det["unresolved_bound"]["desc"], "power law in (E - E_th) from zero at the threshold "
        "(%.4g eV) to the first Simon Wedlund 2019 table value, slope of the table's first interval" % E_th, "factor 3", [1.0, 3.0, 10.0, 30.0])
    det["notes"] += (" v3.2 (decision of 2026-10-31): between the threshold (%.4g eV) and the table start the cross section rises from zero as a "
                     "power law in (E - E_th) continuous in value and log-log slope at the table start (model_derived, factor 3); v3.1 left "
                     "it unresolved, bounded by the Langevin capture rate (documented)." % E_th)
    return groups


def ion_groups(mass):
    def lab(p_, t_):
        return (mass[p_] + mass[t_]) / mass[t_]

    ATOMS = {"H+": ["H"], "H2+": ["H", "H"], "O+": ["O"], "O++": ["O"], "OH+": ["O", "H"], "H2O+": ["O", "H", "H"], "H3O+": ["O", "H", "H", "H"],
             "H-": ["H"], "O-": ["O"], "OH-": ["O", "H"]}

    def langevin(p_, q, lo=1e-5, hi=1e3):
        # v3.1: Langevin polarization model to 1 keV. v3.2: above 1 keV the ZBL universal nuclear-stopping momentum
        # transfer (screened Coulomb), joined by a log-log bridge from the Langevin value at 1 keV to the ZBL value at
        # 2 keV (continuous; between the two models)
        lg = {"type": "langevin", "projectile": p_, "target": "H2O", "charge": q, **ALPHA_H2O, "lo": lo, "hi": hi}
        zbl = {"type": "zbl_momentum_transfer", "projectile": p_, "target": "H2O", "projectile_atoms": ATOMS[p_], "target_atoms": ["O", "H", "H"],
               "lo": E_MT_JOIN, "hi": 1e7, "per_decade": 32, "locator": ZBL_LOCATOR}
        hi_tag = "model_derived"
        if p_ == "H+":   # ICRU 49 (NIST PSTAR) recommended nuclear stopping of protons in water vapour
            zbl = {"type": "nuclear_stopping_mt", "table_id": "nist-pstar-water-vapor-nuclear-stopping", "projectile": "H+", "target": "H2O",
                   "lo": E_MT_JOIN, "hi": 1e7}
            hi_tag = "recommended"
        return {"id": "%s + H2O momentum transfer (Langevin polarization model; ZBL nuclear stopping above 1 keV)" % p_, "kind": "ELASTIC", "angular": "ISO",
                "tag": "model_derived", "products": [[p_, "SURVIVOR"], ["H2O", "GAS"]], "known_zero_below_eV": 0.0, "constant_rate_below": True,
                "sigma": {"type": "piecewise", "segments": [
                    {"source": lg, "class": "MODEL", "tag": "model_derived"},
                    {"source": {"type": "loglog_bridge", "a": lg, "Ea": hi, "b": zbl, "Eb": E_MT_JOIN}, "class": "MODEL", "tag": "model_derived"},
                    {"source": zbl, "class": "DATA" if hi_tag == "recommended" else "MODEL", "tag": hi_tag}]},
                "notes": ("H+ above 2 keV: ICRU Report 49 nuclear stopping of protons in water vapour (NIST PSTAR, recommended; Moliere-screened classical "
                          "orbits), converted the same way. " if p_ == "H+" else "") + "No measured ion mobility or momentum-transfer data for this ion in H2O were obtained. To 1 keV (v3.1): classical polarization (Langevin) "
                         "momentum transfer sigma_m = 1.105 k_L / v with the measured H2O polarizability; it omits the ion-dipole term and the "
                         "short-range repulsion, so it is a lower estimate at collision energies above ~1 eV (open item). Above 2 keV (v3.2): the "
                         "equivalent momentum-transfer cross section of the ZBL universal screened-Coulomb potential, sigma_m = 2 S_n / (gamma E) with the "
                         "nuclear stopping S_n summed over the atoms of the projectile (each at the common speed) and of H2O (Bragg additivity); "
                         "isotropic CM scattering with this sigma_m reproduces the summed nuclear stopping exactly. 1-2 keV: log-log bridge between the "
                         "two models (continuous). Uncertainty: ZBL nuclear stopping about 10-20 %% for atoms (larger for molecular binding effects, "
                         "neglected); the bridge region lies between the two models (they differ by up to a factor %s at 1 keV)." % MT_DIFF.get(p_, "~1.2")}

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
         "products": [["H3O+", "BIRTH"], ["OH", "RESERVOIR"]], "known_zero_below_eV": 0.0, "constant_rate_below": True,
         "sigma": {"type": "piecewise", "segments": [
             {"source": const_k("H2O+", k_pt, 1.0, E_TH, "Anicich 1993 JPCRD k(300 K) = 1.85e-9 cm3/s (+-15 %) held constant up to the 300 K mean CM energy 3/2 kT"),
              "class": "MODEL", "tag": "model_derived"},
             {"source": {"type": "loglog_bridge", "a": const_k("H2O+", k_pt, 1.0, E_TH, "Anicich 1993 k(300 K)"), "Ea": E_TH * lab("H2O+", "H2O"),
                         "b": {"type": "inline", "E": Ept, "sigma": spt, "locator": "Lishawa 1990 Fig. 5 first point"}, "Eb": Ept[0]},
              "class": "MODEL", "tag": "model_derived"},
             {"source": {"type": "inline", "E": Ept, "sigma": spt, "locator": "Lishawa et al. JCP 93, 3196 (1990) Fig. 5 reaction 3a (H2O+ + D2-18O), digitized points, E_cm -> E_lab; points within 3% in energy averaged"},
              "class": "DATA", "tag": "measured"}]},
         "unresolved_bound": {"source": {"type": "inline", "E": [Ept[-1], 1e7], "sigma": [max(spt[-3:]), max(spt[-3:])], "locator": "largest of the last three measured points held constant"},
                              "desc": "last measured magnitude held constant"},
         "notes": "Isotopic analogue data (deuterated, 18O target); exothermic by 1.00 eV (ATcT); products from an isotropic compound decay with all excess as kinetic energy (the measured product energies of Lishawa Table II indicate direct, forward-scattered dynamics; not represented). "
                  "Below the first beam point (E_cm 1.23 eV): Anicich k(300 K) as a constant rate up to the 300 K mean CM energy 3/2 kT = 0.0388 eV, then log-log interpolation to the beam point (continuous); "
                  "the constant-k value at E_cm 1.23 eV would be 43 % above the beam measurement (two measurements, both kept; the difference is the uncertainty of this region); constant sigma*v continues to zero relative energy."},
        {"id": "H2O+ + H2O -> H2O(fast) + H2O+ (symmetric charge exchange)", "kind": "CHARGE_X", "angular": "FORWARD", "tag": "measured",
         "products": [["H2O", "RESERVOIR"], ["H2O+", "BIRTH"]], "known_zero_below_eV": 0.0, "constant_rate_below": True,
         "sigma": {"type": "piecewise", "segments": [
             {"source": {"type": "rate_anchor", "projectile": "H2O+", "lo": 1e-5, "E1": Ecx[0],
                         "ref": {"type": "inline", "E": Ecx, "sigma": scx, "locator": "Lishawa 1990 Fig. 7 first point"}},
              "class": "MODEL", "tag": "model_derived"},
             {"source": {"type": "inline", "E": Ecx, "sigma": scx,
                         "locator": "Lishawa et al. 1990 Fig. 7 reactions 4a and 4b (no isotope effect stated by the source), digitized points, E_cm -> E_lab; points within 3% averaged"},
              "class": "DATA", "tag": "measured"}]},
         "unresolved_bound": bound_k("H2O+", k_pt),
         "notes": "No thermal measurement of symmetric charge exchange exists (identical products); below the first beam point (E_cm 0.79 eV) a constant rate anchored at that point (continuous, model_derived)."},
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
                          "sink": "FRAGMENT", "tag": "model_derived", "products": prods, "known_zero_below_eV": 0.0, "constant_rate_below": True,
                          "sigma": {"type": "piecewise", "segments": [{"source": const_k(ion, k, f, 1.0, "Anicich 1993 JPCRD k(300 K) %g cm3/s" % (k * 1e6)),
                                                                       "class": "MODEL", "tag": "model_derived"}]},
                          "unresolved_bound": bound_k(ion, k * f),
                          "notes": "Measured 300 K rate coefficient and branching (recommended) extended as a constant rate coefficient to E_cm = 1 eV (model); energy dependence above is unresolved."})
        if ion == "H+":
            # one reaction (charge transfer H+ + H2O -> H + H2O+): thermal constant rate to E_cm = 1 eV, unresolved gap
            # bounded by k/v, Simon Wedlund 2019 from 52 eV/u; forward charge-transfer kinematics throughout
            d = ext('ion-wedlund2019-hplus-h2o-electron-capture-sigma10')
            procs = [{"id": "H+ + H2O -> H + H2O+ (charge transfer)", "kind": "CHARGE_X", "angular": "FORWARD", "tag": "recommended",
                      "products": [["H", "RESERVOIR"], ["H2O+", "BIRTH"]], "known_zero_below_eV": 0.0, "constant_rate_below": True,
                      "sigma": {"type": "piecewise", "segments": [
                          {"source": const_k(ion, k, 1.0, 1.0, "Anicich 1993 JPCRD k(300 K) %g cm3/s" % (k * 1e6)), "class": "MODEL", "tag": "model_derived"},
                          {"source": {"type": "inline", "E": [x * mass["H+"] / U_EV for x in d['x']], "sigma": d['y'],
                                      "locator": "Simon Wedlund et al. A&A 2019 Table A.2 sigma_10 (recommended fit), eV/u x m(H+)/u -> E_lab"},
                           "class": "DATA", "tag": "recommended"}]},
                      "unresolved_bound": bound_k("H+", k),
                      "notes": "Measured 300 K rate coefficient (Anicich 1993) as a constant rate up to E_cm = 1 eV and to zero relative energy; "
                               "no data between E_cm 1 eV and 52 eV/u (unresolved, bounded by k/v); keV recommended cross section above. One process for the one "
                               "reaction (v3 listed a thermal capture record and a keV record separately)."}]
        procs.append(langevin(ion, 1))
        groups.append({"projectile": ion, "target": "H2O", "processes": procs})
    d = ext('ion-heninger1989-oplus-h2o-sigma-vs-ecm')
    k_o = ext('ion-anicich1993-oplus-h2o-k300')['y'][0]
    groups.append({"projectile": "O+", "target": "H2O", "processes": [
        {"id": "O+ + H2O -> H2O+ + O (charge transfer)", "kind": "CAPTURE", "kinetic_fraction": 1.0, "sink": "FRAGMENT", "tag": "measured",
         "products": [["H2O+", "BIRTH"], ["O", "RESERVOIR"]], "known_zero_below_eV": 0.0, "constant_rate_below": True,
         "notes": "Below the first ICR point (E_cm 0.033 eV, thermal) a constant rate equal to that point's sigma*v (continuous); the Anicich k(300 K) = 2.6e-9 cm3/s is 24 % above it (uncertainty of this region).",
         "sigma": {"type": "piecewise", "segments": [
             {"source": {"type": "rate_anchor", "projectile": "O+", "lo": 1e-5, "E1": d['x'][0] * lab("O+", "H2O"),
                         "ref": {"type": "inline", "E": [x * lab("O+", "H2O") for x in d['x']], "sigma": d['y'], "locator": "Heninger 1989 Table 1 first point"}},
              "class": "MODEL", "tag": "model_derived"},
             {"source": {"type": "inline", "E": [x * lab("O+", "H2O") for x in d['x']], "sigma": d['y'],
                         "locator": "Heninger et al. GRL 1989 Table 1 (O+(4S) + H2O, ICR), E_cm -> E_lab"}, "class": "DATA", "tag": "measured"}]},
         "unresolved_bound": {"source": {"type": "inline", "E": [11.6 * lab("O+", "H2O"), 1e7], "sigma": [1.77e-19, 1.77e-19],
                                         "locator": "largest of the last three Heninger points held constant"}, "desc": "last measured magnitude held constant"}},
        langevin("O+", 1)]})
    groups.append({"projectile": "O++", "target": "H2O", "processes": [langevin("O++", 2)]})
    # ---- negative ions + H2O
    d = ext('ion-wedlund2019-hminus-h2o-one-electron-loss-sigmam10')
    nodata = lambda what, prods, proj: {"id": what + " (no data)", "kind": "CAPTURE", "kinetic_fraction": 1.0, "sink": "FRAGMENT", "tag": "model_derived",
                                        "products": prods, "known_zero_below_eV": 0.0, "placeholder": True,
                                        "sigma": {"type": "piecewise", "segments": [{"source": {"type": "none", "lo": 1e-5, "hi": 1e7, "locator": what}, "class": "UNRESOLVED", "tag": "model_derived"}]},
                                        "notes": "Placeholder: no rate or cross section obtained; never sampled; its exposure is recorded in the coverage ledger."}
    groups.append({"projectile": "H-", "target": "H2O", "processes": [
        {"id": "H- + H2O -> H + e + H2O (collisional detachment)", "kind": "BREAKUP", "angular": "FORWARD", "loss_eV": 0.0, "kinetic_fraction": 1.0,
         "sink": "NONE", "tag": "model_derived",
         "products": [["H", "RESERVOIR"], ["e", "BIRTH"], ["H2O", "GAS"]],
         "sigma": {"type": "piecewise", "segments": [{"source": {"type": "inline", "E": [x * mass["H-"] / U_EV for x in d['x']], "sigma": d['y'],
                                                                 "locator": "Simon Wedlund et al. 2019 Table B.1 sigma_-10 (low confidence, H2/O2 additive rule scaled to one H2O measurement)"},
                                                      "class": "MODEL", "tag": "model_derived"}]},
         "unresolved_bound": {"source": {"type": "langevin", "projectile": "H-", "target": "H2O", "charge": -1, **ALPHA_H2O, "factor": 1.0, "lo": 1e-5, "hi": 1e7},
                              "desc": "Langevin capture rate k_L/v (ion-molecule capture limit) between the 0.75 eV threshold and the table start"},
         "notes": "Exactly zero below the endothermic threshold (EA(H) = 0.754 eV, ATcT); no data between the threshold and 52 eV/u: unresolved, bounded by the Langevin capture rate."},
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
         "products": [["H2", "RESERVOIR"], ["e", "BIRTH"]], "known_zero_below_eV": 0.0, "constant_rate_below": True,
         "sigma": {"type": "piecewise", "segments": [{"source": {"type": "rate_k", "projectile": "H-", "target": "H", "k_m3_s": km, "lo": 1e-5,
                                                                 "hi": 0.1 * (mass["H-"] + mass["H"]) / mass["H"],
                                                                 "locator": "Miller et al. 2011 merged beams, mean of the three 16 meV values %.4g m3/s, held constant to E_cm = 0.1 eV" % km},
                                                      "class": "MODEL", "tag": "model_derived"}]},
         "unresolved_bound": {"source": {"type": "rate_k", "projectile": "H-", "target": "H", "k_m3_s": km, "lo": 1e-5, "hi": 1e7, "per_decade": 16,
                                         "locator": "measured 16 meV rate / v (the source reports a rate falling with energy; declared bound)"},
                              "desc": "measured 16 meV rate coefficient / v"}}]})
    # ---- dissociative recombination
    dr = []
    br = ext('rec-dr-rosen2000-h2oplus-branching-e0')['y']
    al = ext('rec-dr-rosen2000-h2oplus-alpha-300k')['y'][0]
    rz = ext('rec-dr-rosen2000-h2oplus-sigma-digitized')
    for prods, f in (([["OH", "RESERVOIR"], ["H", "RESERVOIR"]], br[0]), ([["O", "RESERVOIR"], ["H2", "RESERVOIR"]], br[1]),
                     ([["O", "RESERVOIR"], ["H", "RESERVOIR"], ["H", "RESERVOIR"]], br[2])):
        pl = {"type": "dr_powerlaw", "alpha_m3_s": al, "T_K": 300.0, "fraction": f, "segments": [[0.0, 0.1, 1.24], [0.1, 2.0, 1.55]], "lo": 1e-9, "hi": 2.0,
              "locator": "Rosen et al. Faraday Discuss. 115, 295 (2000): alpha(300 K) 4.3e-7 cm3/s, sigma ~ E^-1.24 below 0.1 eV and E^-1.55 to 2 eV (text), branching at 0 eV"}
        dr.append({"id": "e + H2O+ -> %s (DR)" % ' + '.join(n for n, _ in prods), "kind": "CAPTURE", "kinetic_fraction": 1.0, "sink": "FRAGMENT",
                   "tag": "model_derived", "products": prods, "known_zero_below_eV": 0.0, "constant_rate_below": True,
                   "sigma": {"type": "piecewise", "segments": [
                       {"source": pl, "class": "MODEL", "tag": "model_derived"},
                       {"source": {"type": "product", "factor": f, "a": {"type": "joined_table", "E0": 2.0,
                                   "anchor": dict(pl, fraction=1.0), "table": {"type": "extracted", "table_id": "rec-dr-rosen2000-h2oplus-sigma-digitized"}}},
                        "class": "DATA", "tag": "measured"}]},
                   "unresolved_bound": {"source": dict(pl, lo=rz['x'][-1], hi=1e7, per_decade=16, segments=[[0.0, 0.1, 1.24], [0.1, 1e7, 1.55]]),
                                        "desc": "E^-1.55 power-law continuation above the last digitized point (%.4g eV)" % rz['x'][-1]},
                   "notes": "Ground-state products with all exothermicity as kinetic energy (the source quotes ground-product releases, not measured KER); branching held at its 0 eV value at all energies (model). "
                            "2.0-%.4g eV: Rosen Fig. 4 measured cross section, figure-digitized (resonances near 4 and 11 eV included), joined continuously to the stated E^-1.55 law at 2 eV; "
                            "above %.4g eV unresolved (the scattered last bins of the figure exceed a factor-2 reading uncertainty), bounded by the E^-1.55 continuation. "
                            "The printed Eq. 5 product inconsistency of the source (channel b) is kept as the balanced Eq. 4 lineage O + H2." % (rz['x'][-1], rz['x'][-1])})
    groups.append({"projectile": "e", "target": "H2O+", "processes": dr, "energy_bottom_eV": 1e-9})
    dr = []
    br = ext('rec-dr-jensen2000-h3oplus-branching-e0')['y']
    al = ext('rec-dr-jensen2000-h3oplus-alpha-300k')['y'][0]
    jz = ext('rec-dr-jensen2000-h3oplus-sigma-digitized')
    for prods, f in (([["H2O", "RESERVOIR"], ["H", "RESERVOIR"]], br[0]), ([["OH", "RESERVOIR"], ["H2", "RESERVOIR"]], br[1]),
                     ([["OH", "RESERVOIR"], ["H", "RESERVOIR"], ["H", "RESERVOIR"]], br[2]), ([["O", "RESERVOIR"], ["H2", "RESERVOIR"], ["H", "RESERVOIR"]], br[3])):
        pl = {"type": "dr_powerlaw", "alpha_m3_s": al, "T_K": 300.0, "fraction": f, "segments": [[0.0, 0.2, 1.15]], "lo": 1e-9, "hi": 0.2,
              "locator": "Jensen et al. ApJ 543, 764 (2000): alpha(300 K) 4.3e-7 cm3/s, sigma ~ E^-1.15 at low energy (text), branching at 0 eV"}
        dr.append({"id": "e + H3O+ -> %s (DR)" % ' + '.join(n for n, _ in prods), "kind": "CAPTURE", "kinetic_fraction": 1.0, "sink": "FRAGMENT",
                   "tag": "model_derived", "products": prods, "known_zero_below_eV": 0.0, "constant_rate_below": True,
                   "sigma": {"type": "piecewise", "segments": [
                       {"source": pl, "class": "MODEL", "tag": "model_derived"},
                       {"source": {"type": "product", "factor": f, "a": {"type": "joined_table", "E0": 0.2,
                                   "anchor": dict(pl, fraction=1.0), "table": {"type": "extracted", "table_id": "rec-dr-jensen2000-h3oplus-sigma-digitized"}}},
                        "class": "DATA", "tag": "measured"}]},
                   "unresolved_bound": {"source": dict(pl, lo=jz['x'][-1], hi=1e7, per_decade=16, segments=[[0.0, 1e7, 1.15]]),
                                        "desc": "E^-1.15 power-law continuation above the last digitized point (%.4g eV)" % jz['x'][-1]},
                   "notes": "0.2-%.4g eV: Jensen Fig. 3 measured cross section (this work, divided by the plotted x10), figure-digitized, joined continuously to the stated E^-1.15 law at 0.2 eV; "
                            "above 6 eV the overlapping previous-ASTRID points are traced together with this work (the source states agreement outside the 0.6-5 eV valley); above %.4g eV unresolved, "
                            "bounded by the E^-1.15 continuation. Branching held at its 0 eV value (model). Other measured alpha(300 K): Neau 2000 7.6e-7, Herd 1990 1.0e-6 cm3/s (spread kept, not reconciled)." % (jz['x'][-1], jz['x'][-1])})
    groups.append({"projectile": "e", "target": "H3O+", "processes": dr, "energy_bottom_eV": 1e-9})
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
    refs = {x['name']: x['ref_eV'] for x in sp}
    groups += v32_gap_bridges(v32_ion_extend(ion_groups(mass), mass, refs), mass, refs)
    spec = {"id": "pf-chemistry-v3.2-production", "energy_bottom_eV": 1e-5, "energy_top_eV": 1e7, **common, "groups": groups}
    json.dump(spec, open(out, 'w'), indent=1)
    print('groups', len(groups), 'processes', sum(len(g['processes']) for g in groups))


if __name__ == '__main__':
    main(*sys.argv[1:])
