#!/usr/bin/python -I
"""Collects the coordinator v6.1 qualification evidence into qualification/QUALIFICATION-SUMMARY.json
(post-processing only). Inputs: identity-*.json, checks-v6-vs-v61.json and timings-*.json written by
tools/compare_v6_v61.py, the build receipt, the source v6.1 run receipt and the profile.

  write_summary_v61.py BUILD_LABEL SOURCE_RUN_DIR PROFILE_DIR
"""
import json, sys
from pathlib import Path

OUT = Path('/home/treefor/projects/.workspace/outputs/pic-water-common-commit-v6.1-20261031')
Q = OUT / 'qualification'


def load(p):
    return json.loads(Path(p).read_text())


def main(label, source_run, profile_dir):
    br = load(OUT / ('build-' + label) / 'build-receipt.json')
    s = {'schema': 'coordinator-v6.1-qualification-summary-v1',
         'build': {'label': label, 'binary_sha256': br['binary_sha256'], 'status': br['status']},
         'reference': 'coordinator v6 build q3 (binary eb1936f1...), its qualification-q3 runs and the pf-pilot-1 record'}
    ident = {}
    for f in sorted(Q.glob('identity-*.json')):
        r = load(f)
        ident[f.stem.replace('identity-', '')] = {'reference_v6': r['reference_v6'], 'steps_compared': r['steps_compared'], 'bit_identical': r['bit_identical'],
                                                  'files': r['files'], 'ledger_deterministic_columns': r['ledger_deterministic_columns'],
                                                  'ledger_columns_differing': r['ledger_columns_differing']}
    s['1_bit_identity_v6_vs_v61'] = ident
    s['1_all_bit_identical'] = all(v['bit_identical'] for v in ident.values())
    c = load(Q / 'checks-v6-vs-v61.json')
    s['2_checks'] = {k: c[k] for k in ('cases_v6', 'cases_v61', 'passed_v6', 'passed_v61', 'same_cases_and_outcomes', 'refusal_points')}
    src = load(Path(source_run) / 'run-receipt.json')
    s['3_source_v61'] = {'run': str(source_run), 'all_passed': src['all_passed'],
                         'test_source_v6_1': {k: src['normal/test_source_v6_1'][k] for k in ('passed', 'failed', 'cases')},
                         'sealed_v6_suite_on_v61': {k: src['normal/test_source_v6_on_v6_1'][k] for k in ('passed', 'failed', 'cases')}}
    t = {}
    for f in sorted(Q.glob('timings-*.json')):
        data = load(f)
        for run, rec in data.get('runs', data).items():
            t[run] = {k: round(v['median'], 4) for k, v in rec.items() if isinstance(v, dict)}
            t[run]['steps'] = rec.get('steps', 2000)
    s['4_timings_median_ms'] = t
    p = load(Path(profile_dir) / 'profile.json')
    s['5_profile'] = {k: p[k] for k in p if k != 'patches'}
    s['5_profile']['patches_bit_identical'] = p['all_patches_bit_identical']
    s['5_profile']['note'] = ('source_owner_digests_equal_v6 compares a freshly bound sealed v6 owner with the live v6.1 owner: equal on '
                              'the pilot fixture (no release yet); on the qual fixture the live owner has released neutrals and its state '
                              'differs from a fresh owner by design. Owner-state identity with v6 is shown by the step digests (accepted '
                              'source seal) and by the source v6.1 identity tests.')
    (Q / 'QUALIFICATION-SUMMARY.json').write_text(json.dumps(s, indent=1) + '\n')
    print(json.dumps({'summary': str(Q / 'QUALIFICATION-SUMMARY.json'), 'all_bit_identical': s['1_all_bit_identical'],
                      'checks': s['2_checks']['passed_v61'], 'source': s['3_source_v61']['all_passed']}))


if __name__ == '__main__':
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
