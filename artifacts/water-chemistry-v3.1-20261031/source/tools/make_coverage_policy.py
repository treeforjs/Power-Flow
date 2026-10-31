#!/usr/bin/env python3
"""Recommended production coverage policy of chemistry v3.1: the admitted-region list native v8 loads.

usage: make_coverage_policy.py MANIFEST.json OUT.json
Rule (owner decision 2026-10-31): bounded regions refuse iff B_r > max(event_limit, fraction * A); an unbounded
region refuses iff its exposure X_r > 0, unless it is admitted. Admitted unbounded regions: (i) regions above the
top of the tables (10 MeV), (ii) regions whose every unresolved process is a no-data placeholder (never sampled;
exposure recorded). Every other unbounded region refuses.
"""
import hashlib
import json
import sys


def main(mp, out):
    m = json.load(open(mp))
    placeholders = {p['id'] for p in m['processes'] if p['id'].endswith('(no data)') or '(no data)' in p['id']
                    or all(seg['source'].startswith('no data') for seg in p['segments'])}
    admitted, refused, bounded = [], [], []
    for r in m['region_info']:
        rec = {'index': r['index'], 'name': r['name'], 'group': r['group'], 'decade_eV': r['decade_eV'], 'processes': r['processes']}
        if r['bounded']:
            bounded.append(rec)
            continue
        if r['above_grid']:
            admitted.append(dict(rec, reason='above the top of the tables (10 MeV)'))
        elif r.get('placeholders_only') or all(p in placeholders for p in r['processes']):
            admitted.append(dict(rec, reason='no-data placeholder process(es): never sampled, exposure recorded'))
        else:
            refused.append(dict(rec, reason='unbounded, not admitted (refuses when exposed)'))
    pol = {
        'schema': 'pfc-chem-v3.1-coverage-policy-1',
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
        'placeholder_processes': sorted(placeholders),
    }
    json.dump(pol, open(out, 'w'), indent=1)
    print('regions', len(m['region_info']), 'bounded', len(bounded), 'admitted unbounded', len(admitted), 'refusing unbounded', len(refused))
    for r in refused:
        print('  refusing:', r['index'], r['group'], 'decade', r['decade_eV'], [p[:45] for p in r['processes']])


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
