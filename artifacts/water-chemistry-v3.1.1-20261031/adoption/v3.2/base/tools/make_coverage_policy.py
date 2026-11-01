#!/usr/bin/env python3
"""Recommended production coverage policy of chemistry v3.2: the admitted-region list native v8 loads.

usage: make_coverage_policy.py MANIFEST.json OUT.json
Rule (decision of 2026-10-31): bounded regions refuse iff B_r > max(event_limit, fraction * A); an unbounded
region refuses iff its exposure X_r > 0, unless it is admitted. Admitted unbounded regions: (i) regions above the
top of the tables (10 MeV), (ii) regions whose every unresolved process is a placeholder there (never sampled;
exposure recorded): a no-data process, or (v3.2) a process below its placeholder_below_eV edge (new processes
keep the v3.1 placeholder treatment below 1 keV). Every other unbounded region refuses.
"""
import hashlib
import json
import sys


def main(mp, out):
    m = json.load(open(mp))
    placeholders = sorted(p['id'] for p in m['processes'] if p.get('placeholder'))
    partial = sorted('%s (below %g eV)' % (p['id'], p['placeholder_below_eV']) for p in m['processes']
                     if not p.get('placeholder') and p.get('placeholder_below_eV', 0.0) > 0)
    admitted, refused, bounded = [], [], []
    for r in m['region_info']:
        rec = {'index': r['index'], 'name': r['name'], 'group': r['group'], 'decade_eV': r['decade_eV'], 'processes': r['processes'],
               'placeholders_here': r.get('placeholder_processes', [])}
        if r['bounded']:
            bounded.append(rec)
            continue
        if r['above_grid']:
            admitted.append(dict(rec, reason='above the top of the tables (10 MeV)'))
        elif r['placeholders_only']:
            admitted.append(dict(rec, reason='placeholder process(es) here (no data, or below the placeholder edge of a v3.2 process): '
                                             'never sampled, exposure recorded'))
        else:
            refused.append(dict(rec, reason='unbounded, not admitted (refuses when exposed)'))
    # v3.2 model estimates (decision B of 2026-10-31): sampled equal-speed / free-collision estimates and the documented
    # physics bound each sits below (the bound applies to the sum of the product channels of one rule)
    est = {}
    for p in m['processes']:
        db = p.get('documented_bound')
        if not db:
            continue
        key = (p['group'], db['desc'], db['rule'])
        e = est.setdefault(key, {'group': p['group'], 'processes': [], 'rule': db['rule'], 'uncertainty': db['uncertainty'],
                                 'documented_bound': db['desc'], 'pts': {}})
        e['processes'].append(p['id'])
        for pt in db['points']:
            q = e['pts'].setdefault(pt['E_eV'], {'E_eV': pt['E_eV'], 'estimate_m2': 0.0, 'bound_m2': pt['bound_m2']})
            q['estimate_m2'] += pt['sigma_m2']
    estimates = []
    for e in est.values():
        pts = []
        for E in sorted(e.pop('pts').items()):
            q = E[1]
            q['bound_over_estimate'] = (q['bound_m2'] / q['estimate_m2']) if (q['bound_m2'] is not None and q['estimate_m2'] > 0) else None
            pts.append(q)
        e['points'] = pts
        estimates.append(e)
    pol = {
        'schema': 'pfc-chem-v3.2-coverage-policy-1',
        'blob_sha256': m['blob_sha256'],
        'manifest_sha256': hashlib.sha256(open(mp, 'rb').read()).hexdigest(),
        'region_count': len(m['region_info']),
        'policy': {
            'bounded_regions': 'refuse the step iff B_r > max(event_limit, fraction * A)',
            'event_limit_physical_events': 1.0,
            'fraction': 0.01,
            'A': 'applied physical events of the step: sum of event weights W over all real events of both collision stages',
            'B_r': 'bound of region r in the step: sum over carriers and stages of W n v sigma_bound dt (physical events)',
            'unbounded_regions': 'refuse the step iff the exposure X_r > 0, unless the region index is in admitted_unbounded_regions',
            'region_identity': 'regions are identified by index and name (blob strings section); a mismatch of blob_sha256 refuses the policy',
        },
        'admitted_unbounded_regions': admitted,
        'refusing_unbounded_regions': refused,
        'bounded_regions': bounded,
        'placeholder_processes': placeholders,
        'placeholder_below_edge_processes': partial,
        'model_estimates_note': ('v3.2 (decision of 2026-10-31): where no measurement exists above 1 keV, the dominant channels of the heavy '
                                 'ions and anions are sampled from equal-speed or free-collision model estimates (model_derived, resolved, not '
                                 'part of the refusal rule). Each estimate is listed with its documented physics upper bound and the ratio '
                                 'bound / estimate (documentation only).'),
        'model_estimates': estimates,
    }
    json.dump(pol, open(out, 'w'), indent=1)
    print('regions', len(m['region_info']), 'bounded', len(bounded), 'admitted unbounded', len(admitted), 'refusing unbounded', len(refused))
    for r in refused:
        print('  refusing:', r['index'], r['group'], 'decade', r['decade_eV'], [p[:45] for p in r['processes']])


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
