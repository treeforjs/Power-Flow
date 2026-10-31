#!/usr/bin/env python3
"""Human-readable species, process, coverage and seam tables of a built chemistry v3 manifest.

usage: report_processes.py MANIFEST.json OUT.txt
"""
import json
import sys


def fmt_e(x):
    if x >= 1e299:
        return 'inf'
    return '%.6g' % x


def main(mp, out):
    m = json.load(open(mp))
    L = []
    L.append('Chemistry v3 tables: %s' % m['spec_id'])
    L.append('blob sha256 %s (%d bytes); %d species, %d groups, %d processes, %d coverage regions' % (
        m['blob_sha256'], m['blob_bytes'], len(m['species']), len(m['groups']), len(m['processes']), len(m['regions'])))
    L.append('')
    L.append('Species (rest energy m c^2 and 0 K reference energy R = DHf(0 K) relative to H2(g), O2(g); ions in the stationary-electron convention):')
    L.append('  %-6s %-8s %3s %22s %14s  %s' % ('name', 'kind', 'q', 'mass (eV)', 'R (eV)', 'reference source'))
    for s in m['species']:
        L.append('  %-6s %-8s %+3d %22.6f %14.6f  %s' % (s['name'], s['kind'], s['charge'], s['mass_eV'], s['ref_eV'], s.get('ref_source', '')))
    L.append('')
    L.append('Processes. Coverage classes per energy interval (projectile kinetic energy in the target rest frame, eV):')
    L.append('  ZERO = known zero (below the thermochemical/kinematic threshold or a source statement), DATA = measured/recommended/')
    L.append('  calculated tables, MODEL = model_derived values (tagged in every event), UNRESOLVED = no data: sampled as zero,')
    L.append('  exposure and bound recorded in the coverage ledger.')
    group = None
    for p in m['processes']:
        if p['group'] != group:
            group = p['group']
            L.append('')
            L.append('== %s' % group)
        L.append('- %s' % p['id'])
        L.append('    kind %s, tag %s, Q = %+.6f eV, kinematic threshold %.6g eV%s%s' % (
            p['kind'], p['tag'], p['q_eV'], p['kinematic_threshold_eV'],
            (', loss %.6g eV' % p['loss_eV']) if p.get('loss_eV') not in (None, 0.0) else '',
            (', angular %s' % p['angular']) if p['kind'] in ('ELASTIC', 'EXCITE', 'DISSOC', 'CHARGE_X') else ''))
        L.append('    products: %s; internal sink %s, kinetic fraction %g' % (
            ', '.join('%s(%s)' % (a, b) for a, b in p['products']), p['sink'], p['kinetic_fraction']))
        L.append('    coverage: %s' % '; '.join('[%s, %s) %s' % (fmt_e(c['lo']), fmt_e(c['hi']), c['class']) for c in p['coverage']))
        for s in p['segments']:
            L.append('    segment [%s, %s] %s/%s: %s' % (fmt_e(s['lo']), fmt_e(s['hi']), s['class'], s['tag'], s['source'][:400]))
        errs = ', '.join('%.3g (%s)' % (e['max_rel_midpoint_error'], e['type']) for e in p['interpolation'])
        L.append('    interpolation: max relative midpoint error %s' % errs)
        if p.get('unresolved_bound'):
            L.append('    unresolved bound: %s' % p['unresolved_bound'])
        if p.get('source_values_below_threshold_set_to_zero'):
            c = p['source_values_below_threshold_set_to_zero']
            L.append('    %d source values below the thermochemical threshold set to zero (largest %.3g m2 at %.6g eV)' % (
                len(c), max(v for _, v in c), max(c, key=lambda x: x[1])[0]))
        for s in p.get('seams', []):
            if s['relative_jump'] is not None:
                L.append('    junction at %.6g eV: relative jump %+.4g' % (s['energy_eV'], s['relative_jump']))
        if p.get('notes'):
            L.append('    notes: %s' % p['notes'])
    open(out, 'w').write('\n'.join(L) + '\n')
    print('wrote', out, len(L), 'lines')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
