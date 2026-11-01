#!/usr/bin/env python3
"""Expected event share of the bounded (unresolved) coverage regions for a monoenergetic ion beam in water vapour.

usage: ion_beam_share.py BLOB MANIFEST OUT.json [E_eV ...]   (default energies 1e4 2e4 3e4 5e4 1e5 eV)

For every heavy-projectile group with target H2O and every beam energy E (projectile kinetic energy in the gas
frame; the thermal motion of the gas is neglected, |v_gas| << v_beam above 1 keV): A(E) = sum_j sigma_j(E), the
sampled (real) cross sections of the group at E, read from the blob rows exactly as the kernel interpolates them,
and for the coverage region r of the interval containing E its bound B_r(E) (blob bound row; zero rows count as
bounded by zero). The expected number of events the region could hide per real event of the beam is B_r / A, the
quantity the production rule compares with 1 % (max(1, 0.01 A) per step). Unbounded regions are listed with
share null. Placeholder-only regions (never sampled, admitted) are listed separately.
"""
import json
import math
import struct
import sys


def read_blob(path):
    b = open(path, 'rb').read()
    assert b[:8] == b'PFCHEMV3'
    _, n = struct.unpack_from('<II', b, 8)
    sec = {}
    for i in range(n):
        o = 16 + 40 * i
        name = b[o:o + 24].rstrip(b'\0').decode()
        pos, ln = struct.unpack_from('<QQ', b, o + 24)
        sec[name] = b[pos:pos + ln]
    return sec


def main(blob, manifest, out, energies):
    sec = read_blob(blob)
    man = json.load(open(manifest))
    ns, ng, npr, nreg = struct.unpack('<4I', sec['meta'])
    groups = [struct.unpack_from('<16I', sec['groups'], 64 * i) for i in range(ng)]
    grid = struct.unpack('<%dd' % (len(sec['grid']) // 8), sec['grid'])
    sigma = struct.unpack('<%dd' % (len(sec['sigma']) // 8), sec['sigma'])
    bound = struct.unpack('<%dd' % (len(sec['bound']) // 8), sec['bound'])
    bin_region = struct.unpack('<%dI' % (len(sec['bin_region']) // 4), sec['bin_region'])
    bin_row = struct.unpack('<%dI' % (len(sec['bin_row']) // 4), sec['bin_row'])
    bin_off = struct.unpack('<%dI' % (len(sec['bin_off']) // 4), sec['bin_off'])
    names = sec['strings'].decode().split('\n')
    sp_names = names[:ns]
    proc_names = names[ns:ns + npr]
    region_names = names[ns + 2 * npr:ns + 2 * npr + nreg]
    rinfo = {r['index']: r for r in man['region_info']}
    NONE = 0xFFFFFFFF
    res = []
    for gi, g in enumerate(groups):
        proj, targ, goff, gn, poff, pn, soff, boff, bn = g[:9]
        if sp_names[proj] == 'e' or sp_names[targ] != 'H2O':
            continue
        G = grid[goff:goff + gn]
        for E in energies:
            # grid_locate + row_eval of the kernel
            if E < G[0] or E > G[-1]:
                continue
            lo, hi = 0, gn - 1
            if E >= G[-1]:
                i = gn - 1
            else:
                while hi - lo > 1:
                    mid = (lo + hi) // 2
                    if E < G[mid]:
                        hi = mid
                    else:
                        lo = mid
                i = lo

            def row_eval(row):
                if i >= gn - 1:
                    return row[gn - 1] if E == G[gn - 1] else 0.0
                return row[i] + (row[i + 1] - row[i]) * ((E - G[i]) / (G[i + 1] - G[i]))
            per = []
            A = 0.0
            for j in range(pn):
                v = row_eval(sigma[soff + j * gn:soff + (j + 1) * gn])
                A += v
                if v > 0:
                    per.append({'process': proc_names[poff + j], 'sigma_m2': v})
            b = i + 1 if i < gn - 1 else gn
            reg = bin_region[bin_off[gi] + b]
            row = bin_row[bin_off[gi] + b]
            entry = {'group': '%s+%s' % (sp_names[proj], sp_names[targ]), 'E_eV': E, 'A_sigma_sampled_m2': A, 'sampled': per}
            if reg == NONE:
                entry['region'] = None
            else:
                ri = rinfo[reg]
                B = row_eval(bound[boff + row * gn:boff + (row + 1) * gn]) if row != NONE else None
                entry['region'] = {'index': reg, 'name': region_names[reg], 'processes': ri['processes'], 'bounded': ri['bounded'],
                                   'placeholders_only': ri['placeholders_only'], 'B_sigma_bound_m2': B,
                                   'share_B_over_A': (B / A if (B is not None and A > 0) else None)}
            res.append(entry)
    json.dump({'blob': blob, 'energies_eV': energies, 'rows': res}, open(out, 'w'), indent=1)
    for r in res:
        rg = r['region']
        s = 'resolved' if rg is None else ('placeholders only' if rg['placeholders_only'] else
                                           ('UNBOUNDED' if rg['share_B_over_A'] is None else 'B/A = %.3g' % rg['share_B_over_A']))
        print('%-10s %8.3g eV  A %.3e m2  %s%s' % (r['group'], r['E_eV'], r['A_sigma_sampled_m2'], s,
                                                    '' if rg is None else '  [' + '; '.join(p[:40] for p in rg['processes']) + ']'))


if __name__ == '__main__':
    E = [float(x) for x in sys.argv[4:]] or [1e4, 2e4, 3e4, 5e4, 1e5]
    main(sys.argv[1], sys.argv[2], sys.argv[3], E)
