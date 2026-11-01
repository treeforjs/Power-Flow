#!/usr/bin/env python3
"""Small metadata/ledger refusal checks, with no numerical simulation."""
import copy
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('history_audit', ROOT / 'audit.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)
m = a.load(ROOT / 'INPUTS-SHA256.json')
p = Path(m['paths']['pilot'])
schema = a.load(p / 'ledger-schema.json')
data = (p / 'ledger.bin').read_bytes()
rows = a.rows_from_bytes(data, schema)
genesis = a.load(p / 'genesis.json')
summary = a.load(p / 'run-summary.json')
passed = []


def refuses(name, fn, message):
    try:
        fn()
    except ValueError as e:
        assert message in str(e), (name, str(e))
        passed.append(name)
    else:
        raise AssertionError('accepted invalid case: ' + name)


r = a.summarize(rows, genesis, summary, summary['final_time_s'])
assert r['counts']['h2o_macroparticle_births']['recorded_cumulative'] == 0
assert r['counts']['expected_h2o_molecules']['recorded_cumulative'] > 174000
passed.append('zero_sampled_macros_preserve_positive_expected_source')
refuses('20ns_extension', lambda: a.summarize(rows, genesis, summary, 20e-9), 'unsupported requested endpoint')
refuses('missing_accepted_row', lambda: a.summarize(rows[:999] + rows[1000:], genesis, summary), 'noncontiguous')
bad = copy.deepcopy(rows)
bad[100]['time_s'] += 1e-16
refuses('clock_gap', lambda: a.summarize(bad, genesis, summary), 'clock discontinuity')
bad_total = copy.deepcopy(rows)
bad_total[-1]['arr_cum_expected_h2o_molecules'] = 0
refuses('discarded_expected_count', lambda: a.summarize(bad_total, genesis, summary), 'step/cumulative mismatch')
bad_summary = copy.deepcopy(summary)
bad_summary['source_owner']['surface'] = 'ss304-foil-sites-v7'
refuses('wrong_surface', lambda: a.summarize(rows, genesis, bad_summary), 'surface differs')
refuses('truncated_ledger', lambda: a.rows_from_bytes(data[:-1], schema), 'incomplete ledger row')
bad_schema = copy.deepcopy(schema)
bad_schema['columns'][-1] = bad_schema['columns'][0]
refuses('duplicate_schema_column', lambda: a.rows_from_bytes(data, bad_schema), 'duplicate ledger columns')
bad_manifest = copy.deepcopy(m)
first = next(iter(bad_manifest['files_sha256']))
bad_manifest['files_sha256'][first] = '0' * 64
refuses('changed_input_pin', lambda: a.audit(bad_manifest), 'input pin differs')
print(json.dumps({'schema': 'early-source-history-audit-checks-v1', 'passed': passed, 'count': len(passed)}, indent=2))
