#!/usr/bin/env python3
"""Refusal check of the coverage policy for a representative breakdown mix (one cell, one step).

usage: breakdown_mix_share.py BLOB MANIFEST OUT.json

The production rule refuses a step iff, for some bounded region r, B_r > max(1, 0.01 A) (physical events), and iff
an unadmitted unbounded region is exposed. B_r and A both scale with cell volume x dt, so the scale-free test is
max_r B_r / A > 0.01 (the max(1, .) floor only makes the rule more lenient). For each mix this tool evaluates
A = sum of the real collision rates (blob cross sections, interpolated as the kernel does) and B_r = sum of the
bound rates of region r (blob bound rows) per unit volume, from:
  * electrons: Maxwellian at T_e (2, 5, 10, 20 eV) on H2O, and on the H2O+ / H3O+ ions (dissociative recombination);
  * thermal ions and anions (300 K relative Maxwellian on H2O);
  * a fraction f of every ion and anion species at E_b (gas motion neglected): epithermal drifting ions at 3-1000 eV
    (f = 0.1 and 1) and keV ions at 10-100 keV (f = 1e-3 .. 0.1).
The gas thermal motion is neglected for electrons and beams; the exposure energy is the relative energy, as in
the kernel's coverage accounting. Unbounded exposure is reported per region with its admission status.
"""
import json
import math
import struct
import sys

sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from ion_beam_share import read_blob  # noqa: E402

C = 299792458.0
KB = 8.617333262e-5
NONE = 0xFFFFFFFF


def main(blob, manifest, out):
    sec = read_blob(blob)
    man = json.load(open(manifest))
    ns, ng, npr, nreg = struct.unpack('<4I', sec['meta'])
    species = [struct.unpack_from('<ddiI', sec['species'], 64 * i) for i in range(ns)]
    groups = [struct.unpack_from('<16I', sec['groups'], 64 * i) for i in range(ng)]
    flags = [struct.unpack_from('<I', sec['procs'], 512 * i + 84)[0] for i in range(npr)]
    grid = struct.unpack('<%dd' % (len(sec['grid']) // 8), sec['grid'])
    sigma = struct.unpack('<%dd' % (len(sec['sigma']) // 8), sec['sigma'])
    bound = struct.unpack('<%dd' % (len(sec['bound']) // 8), sec['bound'])
    bin_region = struct.unpack('<%dI' % (len(sec['bin_region']) // 4), sec['bin_region'])
    bin_row = struct.unpack('<%dI' % (len(sec['bin_row']) // 4), sec['bin_row'])
    bin_off = struct.unpack('<%dI' % (len(sec['bin_off']) // 4), sec['bin_off'])
    names = sec['strings'].decode().split('\n')
    sp = names[:ns]
    region_names = names[ns + 2 * npr:ns + 2 * npr + nreg]
    rinfo = {r['index']: r for r in man['region_info']}
    pol_adm = None
    mass = {sp[i]: species[i][0] for i in range(ns)}
    gidx = {(sp[g[0]], sp[g[1]]): k for k, g in enumerate(groups)}

    def speed(E, m):
        return C * math.sqrt(E * (E + 2 * m)) / (E + m)

    def at(gk, E):
        """(sum sigma v, region, bound sigma v) of group gk at target-frame energy E, as the kernel evaluates them"""
        proj, targ, goff, gn, poff, pn, soff, boff, bn = groups[gk][:9]
        G = grid[goff:goff + gn]
        m = species[proj][0]
        if E < G[0]:
            i = -1
        elif E >= G[-1]:
            i = gn - 1
        else:
            lo, hi = 0, gn - 1
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if E < G[mid]:
                    hi = mid
                else:
                    lo = mid
            i = lo

        def row_eval(row):
            if i < 0:
                return 0.0
            if i >= gn - 1:
                return row[gn - 1] if E == G[gn - 1] else 0.0
            return row[i] + (row[i + 1] - row[i]) * ((E - G[i]) / (G[i + 1] - G[i]))
        v = speed(E, m)
        sv = 0.0
        if i < 0:
            v0 = speed(G[0], m)
            for j in range(pn):
                if flags[poff + j] & 1:
                    sv += sigma[soff + j * gn] * v0
        else:
            for j in range(pn):
                sv += row_eval(sigma[soff + j * gn:soff + (j + 1) * gn]) * v
        b = 0 if i < 0 else (gn if i >= gn - 1 else i + 1)
        reg = bin_region[bin_off[gk] + b]
        row = bin_row[bin_off[gk] + b]
        bsv = row_eval(bound[boff + row * gn:boff + (row + 1) * gn]) * v if (reg != NONE and row != NONE) else 0.0
        return sv, reg, bsv, v

    def maxwell_E(kT, n=600):
        """(E, weight) quadrature of a Maxwellian energy distribution f(E) = 2 sqrt(E/pi) kT^-1.5 exp(-E/kT)"""
        lo, hi = 1e-6 * kT, 40 * kT
        pts = []
        for k in range(n):
            E = lo * (hi / lo) ** ((k + 0.5) / n)
            dE = E * math.log(hi / lo) / n
            pts.append((E, 2 * math.sqrt(E / math.pi) * kT ** -1.5 * math.exp(-E / kT) * dE))
        return pts

    def accumulate(acc, gk, n_proj, n_targ, dist):
        for E, w in dist:
            sv, reg, bsv, v = at(gk, E)
            acc['A'] += n_proj * n_targ * w * sv
            if reg != NONE:
                r = acc['R'].setdefault(reg, {'B': 0.0, 'X': 0.0})
                r['B'] += n_proj * n_targ * w * bsv
                r['X'] += n_proj * n_targ * w * v

    n_gas = 2.5e25                 # H2O, 1 atm 300 K (cancels except for the electron-ion terms)
    n_e = 1e21
    pos = {"H3O+": 0.50, "H2O+": 0.30, "OH+": 0.10, "H+": 0.05, "O+": 0.04, "H2+": 0.005, "O++": 0.005}
    neg = {"OH-": 0.6, "O-": 0.2, "H-": 0.2}
    neg_frac = 0.1
    results = []
    # ion populations besides the 300 K Maxwellian: none; epithermal drifting ions (fraction f of every species at
    # 3-1000 eV, as at high E/N); keV ions (fraction f at 10-100 keV)
    scen = [(0.0, 0.0)] + [(E, f) for E in (3.0, 10.0, 30.0, 100.0, 300.0, 1000.0) for f in (0.1, 1.0)] + \
        [(E, f) for E in (1e4, 5e4, 1e5) for f in (1e-3, 1e-2, 1e-1)]
    for Te in (2.0, 5.0, 10.0, 20.0):
        for Eb, f in scen:
            if True:
                acc = {'A': 0.0, 'R': {}}
                n_neg = neg_frac * n_e
                n_pos = n_e + n_neg * 1.0
                de = maxwell_E(Te)
                accumulate(acc, gidx[("e", "H2O")], n_e, n_gas, de)
                for ion in ("H2O+", "H3O+"):
                    if ("e", ion) in gidx:
                        accumulate(acc, gidx[("e", ion)], n_e, n_pos * pos[ion], de)
                for comp, ntot in ((pos, n_pos), (neg, n_neg)):
                    for ion, share in comp.items():
                        if (ion, "H2O") not in gidx:
                            continue
                        gk = gidx[(ion, "H2O")]
                        mu = mass[ion] * mass["H2O"] / (mass[ion] + mass["H2O"])
                        # relative-speed Maxwellian at 300 K, as projectile energy in the target frame E = (m_P/mu) E_rel
                        dth = [(E * mass[ion] / mu, w) for E, w in maxwell_E(KB * 300.0, 200)]
                        accumulate(acc, gk, ntot * share * (1 - f), n_gas, dth)
                        if f > 0:
                            accumulate(acc, gk, ntot * share * f, n_gas, [(Eb, 1.0)])
                rows = []
                for reg, r in acc['R'].items():
                    ri = rinfo[reg]
                    rows.append({'region': reg, 'name': region_names[reg], 'bounded': ri['bounded'], 'placeholders_only': ri['placeholders_only'],
                                 'above_grid': ri['above_grid'], 'B_over_A': r['B'] / acc['A'] if ri['bounded'] else None,
                                 'exposure_rate_m3_s': r['X']})
                b_rows = sorted([x for x in rows if x['bounded']], key=lambda x: -x['B_over_A'])
                unb = [x for x in rows if not x['bounded'] and not x['placeholders_only'] and not x['above_grid'] and x['exposure_rate_m3_s'] > 0]
                worst = b_rows[0] if b_rows else None
                results.append({'Te_eV': Te, 'beam_eV': Eb, 'beam_fraction': f, 'A_per_m3_s': acc['A'],
                                'max_B_over_A': worst['B_over_A'] if worst else 0.0, 'worst_region': worst['name'] if worst else None,
                                'refuses': bool((worst and worst['B_over_A'] > 0.01) or unb),
                                'unadmitted_unbounded_exposed': [x['name'] for x in unb], 'top_bounded': b_rows[:5]})
                print('Te %4.1f eV  ions %7.0f eV  f %-6g  max B/A %.3g  %s  [%s]' % (Te, Eb, f, results[-1]['max_B_over_A'],
                                                                                   'REFUSES' if results[-1]['refuses'] else 'ok',
                                                                                   (worst['name'][:110] if worst else '')))
    json.dump({'blob': blob, 'mix': {'n_gas_m3': n_gas, 'n_e_m3': n_e, 'positive_ion_shares': pos, 'anion_shares': neg,
                                     'anion_density_over_n_e': neg_frac, 'ion_temperature_K': 300.0},
               'rule': 'refuse iff max_r B_r / A > 0.01 or an unadmitted unbounded region is exposed', 'cases': results},
              open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(*sys.argv[1:4])
