#!/usr/bin/env python3
"""Collect the v7 qualification and calibration records into one summary (post-processing only).

  write_summary_v7.py OUT_ROOT
Reads OUT_ROOT/qualification/*.json (identity, checks, receipts) and OUT_ROOT/calibration/*.json and writes
OUT_ROOT/qualification/QUALIFICATION-SUMMARY.json.
"""
import json, sys
from pathlib import Path

root = Path(sys.argv[1])
q = root / 'qualification'
summary = {'schema': 'coordinator-v7-qualification-summary-v1', 'identity_vs_v61': {}, 'checks': None, 'receipts': {}, 'calibration': {}}
for f in sorted(q.glob('identity-*.json')):
    r = json.loads(f.read_text())
    summary['identity_vs_v61'][f.stem] = {'reference': r.get('reference_v61'), 'candidate': r.get('candidate_v7'), 'steps': r['steps_compared'],
                                          'bit_identical': r['bit_identical'], 'columns_differing': r['ledger_columns_differing'],
                                          'deterministic_columns': r['ledger_deterministic_columns'],
                                          'frames': r['files'].get('frames.bin')}
for f in sorted(q.glob('checks-*.json')):
    r = json.loads(f.read_text())
    summary['checks'] = {k: r[k] for k in ('cases_v61', 'cases_v7', 'passed_v61', 'passed_v7', 'v61_cases_first_in_v7_same_names_and_outcomes')}
    summary['checks']['evidence_text_differences'] = [d['name'] for d in r['evidence_differences']]
    summary['checks']['new_v7_cases'] = r['new_v7_cases']
for f in sorted(q.glob('receipts-*.json')):
    r = json.loads(f.read_text())
    summary['receipts'][f.stem] = {'receipts': r['receipts'], 'verified': r['verified']}
for f in sorted((root / 'calibration').glob('*.json')):
    r = json.loads(f.read_text())
    if 'observables' not in r:
        summary['calibration'][f.stem] = r
        continue
    keep = {}
    for name in ('probe:W_E_J', 'probe:W_B_J', 'cum:E_field_delta_J', 'cum:E_circuit_work_J', 'cum:E_plasma_work_J', 'cum:E_joule_J',
                 'cum:E_thermal_delta_J', 'E_census_kinetic_J', 'arr_cum_expected_electrons', 'cov_applied_expected_events_cum',
                 'q_collected_cum_C', 'cum:parent_events', 'cum:sinks', 'carriers_electron', 'carriers_ion', 'carriers_h2o'):
        if name in r['observables']:
            o = r['observables'][name]
            keep[name] = {k: o[k] for k in ('kind', 'rel_error_max', 'order_p', 'order_p3_end', 'dt_eq_s', 'end_values', 'z_end') if k in o}
    summary['calibration'][f.stem] = {'window_s': r['window_s'], 'comparison_times': r['comparison_times'], 'runs': r['runs'],
                                      'aggregates': r.get('aggregates'), 'observables': keep}
(q / 'QUALIFICATION-SUMMARY.json').write_text(json.dumps(summary, indent=1) + '\n')
print(json.dumps({'identity': {k: v['bit_identical'] for k, v in summary['identity_vs_v61'].items()}, 'checks': summary['checks'] and
                  {k: summary['checks'][k] for k in ('cases_v7', 'passed_v7')}, 'receipts': summary['receipts']}, indent=1))
