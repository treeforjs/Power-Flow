#!/usr/bin/env python3
"""Summarize one qualification attempt: v4 -> v5 port table and key results.
usage: summarize.py <attempt-dir>   (writes PORT-TABLE.txt and summary.json there)"""
import json
import sys
from pathlib import Path

A = Path(sys.argv[1])
V4 = Path('/home/treefor/projects/Power-Flow/artifacts/audit-reproducer-ports-20261030/source-v4/attempt-003/normal/port-results.json')


def load(p):
    try:
        return json.loads(Path(p).read_text())
    except Exception:
        return {}


v4 = {p['id']: p for p in load(V4).get('probes', [])}
v5 = {p['id']: p for p in load(A / 'port-normal/port-results.json').get('probes', [])}
rc = load(A / 'receipt.json')
lines = ['Source-v4 audit port table re-run against owner v5', f'attempt: {A}', '']
for pid in v4:
    lines.append(f'{pid} [{v4[pid].get("audit_finding", "")}]')
    lines.append(f'  v4: {v4[pid].get("v4_verdict")}')
    lines.append(f'  v5: {v5.get(pid, {}).get("v5_verdict", "MISSING")}')
for pid in v5:
    if pid not in v4:
        lines.append(f'{pid}: v5 only: {v5[pid].get("v5_verdict")}')


def probe(pid, *keys):
    p = v5.get(pid, {})
    return {k: p.get(k) for k in keys}


summary = {
    'tests_normal': rc.get('tests_normal'), 'tests_sanitized': rc.get('tests_sanitized'),
    'byte_parity': rc.get('normal_sanitized_byte_parity'),
    'compile_reject_all_as_expected': rc.get('compile_reject_all_as_expected'),
    'precision_oracle': {k: rc.get('precision_oracle_normal', {}).get(k) for k in ('passed', 'failed', 'stats')},
    'port_fraction_oracle': {k: rc.get('port_fraction_oracle_normal', {}).get(k) for k in ('passed', 'failed')},
    'count_oracle': rc.get('count_oracle_normal'),
    'minrepro_normal': rc.get('minrepro_normal', {}).get('stdout'),
    'weight_sweep': probe('EXT-A08-weight-precision-sweep', 'admitted_weights', 'refused_weights', 'typed_scoped_refusals', 'smallest_admitted_weight',
                          'admitted_with_subnormal_weighted_operand', 'worst_admitted_relative_exchange_error',
                          'worst_admitted_current_entry_error_relative_to_gross', 'admitted_current_entries_zeroed'),
    'exposure_sweep': probe('EXT-A08-exposure-scalar-precision-sweep', 'admitted_area_fractions', 'refused_area_fractions', 'typed_scoped_refusals',
                            'admitted_with_relative_error_above_tolerance', 'worst_admitted_relative_error', 'worst_area_fraction'),
    'v4_unchanged': rc.get('v4_unchanged'),
}
(A / 'PORT-TABLE.txt').write_text('\n'.join(lines) + '\n')
(A / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
print('\n'.join(lines))
print(json.dumps(summary, indent=1)[:4000])
