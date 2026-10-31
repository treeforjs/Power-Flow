#!/usr/bin/python -I
"""Collect the multi-step (joint, native candidate003) and Hall-only loop series of one label into one table.

  series_summary.py SERIES_DIRECTORY   (OUT/series-LABEL)

Writes SERIES_DIRECTORY/series-summary.json with, per run: kind, backend, carriers, steps completed, stop
reason, step-time statistics, median per-phase times (ms) and push parity counts, and a CPU/GPU comparison
per carrier count.
"""
import json, sys
from pathlib import Path

root = Path(sys.argv[1])
rows = []
for d in sorted(p for p in root.iterdir() if p.is_dir()):
    rec = d / 'series-receipt.json'
    if not rec.exists():
        continue
    r = json.loads(rec.read_text())
    s = d / ('multistep-summary.json' if r['kind'] == 'multistep' else 'hall-loop-summary.json')
    if not s.exists():
        rows.append({'run': d.name, 'kind': r['kind'], 'backend': r['backend'], 'carriers': r['carriers'], 'error': 'no summary',
                     'returncode': r['command']['returncode']})
        continue
    j = json.loads(s.read_text())
    ledger = d / ('step-ledger.jsonl' if r['kind'] == 'multistep' else 'hall-loop-ledger.jsonl')
    totals = []
    if ledger.exists():
        for line in ledger.read_text().splitlines():
            x = json.loads(line)
            totals.append(x['step_total_ms'])
    st = sorted(totals)
    row = {'run': d.name, 'kind': r['kind'], 'backend': r['backend'], 'carriers': r['carriers'],
           'returncode': r['command']['returncode'], 'wall_s': r['command']['elapsed_s'],
           'steps_completed': j.get('steps_completed'), 'stop': j.get('stop') or j.get('genesis_refusal'),
           'genesis_refusal': j.get('genesis_refusal'),
           'step_ms': {'first': totals[0], 'median': st[len(st) // 2], 'mean': sum(totals) / len(totals), 'max': st[-1],
                       'median_after_first': sorted(totals[1:])[len(totals[1:]) // 2] if len(totals) > 1 else None} if totals else None,
           'phase_ms_median': {k: v['median'] for k, v in (j.get('phase_ms') or {}).items()},
           'push_parity': j.get('push_parity'), 'final_journal_rows': j.get('final_journal_rows'),
           'final_journal_bytes': j.get('final_journal_bytes'), 'final_time_s': j.get('final_time_s'),
           'configuration_sha256': j.get('configuration_sha256'), 'motion_provider': j.get('motion_provider')}
    rows.append(row)
comparison = {}
for r in rows:
    key = '%s/%s' % (r['kind'], r['carriers'])
    comparison.setdefault(key, {})[r['backend']] = {'steps': r.get('steps_completed'),
                                                    'median_step_ms': (r.get('step_ms') or {}).get('median'),
                                                    'stop': r.get('stop')}
out = {'schema': 'stage1-v5-series-summary-v1', 'series_directory': str(root), 'runs': rows, 'cpu_vs_gpu': comparison}
(root / 'series-summary.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(comparison, indent=1))
