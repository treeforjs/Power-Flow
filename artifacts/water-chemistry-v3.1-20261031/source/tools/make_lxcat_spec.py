#!/usr/bin/env python3
"""Spec of one complete LXCat e + H2O cross-section set, for swarm comparisons.

usage: make_lxcat_spec.py LXCAT_TXT SET_NAME OUT_SPEC.json

Mapping (energy-loss bookkeeping only; product branching is not the subject of
these comparison specs): ELASTIC -> ELASTIC isotropic with the momentum-transfer
cross section; EFFECTIVE -> ELASTIC isotropic with sigma_eff minus the sum of
every other process of the set on the union grid (BOLSIG+ convention; negative
differences are clamped at zero and the clamped amount recorded); EXCITATION ->
EXCITE with the set's energy loss; IONIZATION -> IONIZE to the named ion
(H2O+ for totals) with NIST BEB orbitals for the energy sharing; ATTACHMENT ->
CAPTURE to H- + OH, O- + H2 or OH- + H according to the label (H- + OH for
unlabelled totals).
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_tables import parse_lxcat, nudge_duplicates  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ORBS = [{"label": "1b1", "B": 12.61, "U": 61.91, "N": 2, "w": 1.0}, {"label": "3a1", "B": 15.57, "U": 59.52, "N": 2, "w": 1.0},
        {"label": "1b2", "B": 19.83, "U": 48.36, "N": 2, "w": 1.0}, {"label": "2a1", "B": 36.88, "U": 70.71, "N": 2, "w": 1.0}]


def lerp(pts, E):
    if E < pts[0][0] or E > pts[-1][0]:
        return 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= E <= x1:
            return y0 if x1 == x0 else y0 + (y1 - y0) * (E - x0) / (x1 - x0)
    return pts[-1][1]


def main(path, name, out):
    rel = os.path.relpath(os.path.abspath(path), ROOT)
    blocks = parse_lxcat(path)
    common = json.load(open(os.path.join(ROOT, 'tests/specs/common-species.json')))
    procs = []
    occ = {}
    effective = None
    others = []
    for b in blocks:
        key = (b['keyword'], b['process'], '')
        k = occ.get(key, 0)
        occ[key] = k + 1
        src = {"type": "lxcat", "file": rel, "keyword": b['keyword'], "process": b['process'], "occurrence": k}
        hdr = [h.strip() for h in b['header']]
        third = None
        for h in hdr[1:3]:
            try:
                third = float(h.split()[0])
                break
            except (ValueError, IndexError):
                pass
        pid = '%s: %s [%s #%d]' % (name, b['process'], b['keyword'], k)
        if b['keyword'] == 'EFFECTIVE':
            effective = (b, pid)
            continue
        pts = nudge_duplicates(b['table'])
        if b['keyword'] == 'ELASTIC':
            procs.append({"id": pid, "kind": "ELASTIC", "angular": "ISO", "tag": "recommended", "products": [["e", "SURVIVOR"], ["H2O", "GAS"]],
                          "sigma": src, "known_zero_below_eV": 0.0})
        elif b['keyword'] == 'EXCITATION':
            procs.append({"id": pid, "kind": "EXCITE", "loss_eV": third, "sink": "ELEC" if third > 1 else ("VIB" if third > 0.1 else "ROT"),
                          "tag": "recommended", "products": [["e", "SURVIVOR"], ["H2O", "GAS"]], "sigma": src})
        elif b['keyword'] == 'IONIZATION':
            lab = b['process']
            if 'H+' in lab.replace('H2+', '').replace('H_2', '') and 'OH' in lab and 'OH+' not in lab:
                prods, ion = [["H+", "BIRTH"], ["OH", "RESERVOIR"]], 'H+'
            elif 'OH+' in lab:
                prods, ion = [["OH+", "BIRTH"], ["H", "RESERVOIR"]], 'OH+'
            elif 'O^{2+}' in lab or 'O++' in lab:
                prods, ion = [["e", "BIRTH"], ["O++", "BIRTH"], ["H2", "RESERVOIR"]], 'O++'
            elif 'O+' in lab and 'H2' in lab:
                prods, ion = [["O+", "BIRTH"], ["H2", "RESERVOIR"]], 'O+'
            elif 'H_2^{+}' in lab or 'H2+' in lab:
                prods, ion = [["H2+", "BIRTH"], ["O", "RESERVOIR"]], 'H2+'
            else:
                prods, ion = [["H2O+", "BIRTH"]], 'H2O+'
            if ion == 'O++':
                full = [["e", "SURVIVOR"], ["e", "BIRTH"]] + [["O++", "BIRTH"], ["e", "BIRTH"], ["H2", "RESERVOIR"]]
            else:
                full = [["e", "SURVIVOR"], ["e", "BIRTH"]] + prods
            procs.append({"id": pid, "kind": "IONIZE", "tag": "recommended", "sink": "ION_INTERNAL" if ion == 'H2O+' else "FRAGMENT",
                          "kinetic_fraction": 0.0, "orbitals": ORBS, "products": full, "sigma": src})
        elif b['keyword'] == 'ATTACHMENT':
            lab = b['process']
            if 'OH-' in lab or 'OH(-)' in lab:
                prods = [["OH-", "BIRTH"], ["H", "RESERVOIR"]]
            elif 'O-' in lab or 'O(-)' in lab:
                prods = [["O-", "BIRTH"], ["H2", "RESERVOIR"]]
            else:
                prods = [["H-", "BIRTH"], ["OH", "RESERVOIR"]]
            if not b['table']:
                continue
            procs.append({"id": pid, "kind": "CAPTURE", "tag": "recommended", "products": prods, "sigma": src})
        others.append(pts)
    notes = {}
    if effective is not None:
        b, pid = effective
        eff = nudge_duplicates(b['table'])
        grid = sorted(set([E for E, _ in eff] + [E for o in others for E, _ in o if eff[0][0] <= E <= eff[-1][0]]))
        vals, clamped = [], 0.0
        for E in grid:
            v = lerp(eff, E) - sum(lerp(o, E) for o in others)
            if v < 0:
                clamped = max(clamped, -v / max(lerp(eff, E), 1e-300))
                v = 0.0
            vals.append(v)
        procs.insert(0, {"id": pid + ' minus inelastic', "kind": "ELASTIC", "angular": "ISO", "tag": "recommended",
                         "products": [["e", "SURVIVOR"], ["H2O", "GAS"]], "known_zero_below_eV": 0.0,
                         "sigma": {"type": "inline", "E": grid, "sigma": vals,
                                   "locator": "derived: LXCat %s EFFECTIVE (%s line %d) minus the sum of the other processes of the set" % (name, rel, b['line'])}})
        notes['effective_clamped_max_relative'] = clamped
    spec = {"id": "lxcat-%s-swarm-comparison" % name, "energy_bottom_eV": 1e-4, "energy_top_eV": 1e7, **common,
            "notes": notes, "groups": [{"projectile": "e", "target": "H2O", "processes": procs}]}
    json.dump(spec, open(out, 'w'), indent=1)
    print(name, len(procs), 'processes', notes)


if __name__ == '__main__':
    main(*sys.argv[1:4])
