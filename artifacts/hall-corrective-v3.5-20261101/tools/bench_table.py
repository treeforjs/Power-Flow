#!/usr/bin/python -I
"""Tables of the test_v34 bench runs: median over the steps after the first of every timing (ms per step), host CPU
(thread CPU time of the Hall calls: prepare + preflight + publish) per carrier, by mode, mesh and carrier count.
Usage: bench_table.py OUTDIR LABEL [LABEL...] > table.txt   (writes OUTDIR/bench-table.json as well)"""
import json, statistics, sys
from pathlib import Path

out = Path(sys.argv[1]); labels = sys.argv[2:]
rows = {}
for lab in labels:
    for d in sorted((out / lab).glob('gpu-bench-*')):
        f = d / 'bench.json'
        if not f.exists():
            continue
        j = json.loads(f.read_text()); m = j['meta']
        allsteps = j['steps']
        # provider statistics are cumulative over the run: per-step values as differences
        for k in ('push_upload_ms', 'push_admit_ms', 'push_assemble_ms', 'push_device_ms'):
            prev = 0.0
            for s_ in allsteps:
                if k in s_:
                    cur = s_[k]; s_[k] = cur - prev; prev = cur
        steps = allsteps[1:] if len(allsteps) > 1 else allsteps
        def med(path):
            vals = []
            for s in steps:
                x = s
                for k in path:
                    x = x.get(k) if isinstance(x, dict) else None
                    if x is None:
                        break
                if isinstance(x, (int, float)):
                    vals.append(float(x))
            return statistics.median(vals) if vals else None
        keys = {
            'host_cpu_ms': None, 'prepare_cpu_ms': ('prepare_cpu_ms',), 'preflight_cpu_ms': ('preflight_cpu_ms',), 'publish_cpu_ms': ('publish_cpu_ms',),
            'prepare_wall_ms': ('prepare_wall_ms',), 'preflight_wall_ms': ('preflight_wall_ms',), 'iterations': ('iterations',),
            'admit_census_ms': ('admit_census_ms',), 'motion_ms': ('motion_ms',), 'support_and_bank_ms': ('support_and_bank_ms',),
            'field_prepare_ms': ('field_prepare_ms',), 'solve_ms': ('solve_ms',), 'residual_check_ms': ('residual_check_ms',),
            'record_seal_ms': ('record_seal_ms',), 'field_seal_ms': ('field_seal_ms',), 'heat_ms': ('heat_ms',), 'trial_digest_ms': ('trial_digest_ms',),
            'push_upload_ms': ('push_upload_ms',), 'push_admit_ms': ('push_admit_ms',), 'push_assemble_ms': ('push_assemble_ms',),
            'push_device_ms': ('push_device_ms',), 'push_gpu_total_ms': ('push_gpu_total_ms',),
            'cpu_trial_digest_ms': ('cpu', 'trial_digest_ms'), 'cpu_support_and_bank_ms': ('cpu', 'support_and_bank_ms'),
            'cpu_field_prepare_ms': ('cpu', 'field_prepare_ms'), 'cpu_finish_ms': ('cpu', 'finish_ms'),
            'ds_census_cpu_ms': ('device_step', 'census_cpu_ms'), 'ds_motion_cpu_ms': ('device_step', 'motion_cpu_ms'),
            'ds_record_cpu_ms': ('device_step', 'record_cpu_ms'), 'ds_preflight_cpu_ms': ('device_step', 'preflight_cpu_ms'),
            'ds_motion_wait_ms': ('device_step', 'motion_wait_ms'), 'ds_push_gpu_total_ms': ('device_step', 'push_gpu_total_ms'),
        }
        r = {k: (med(p) if p else None) for k, p in keys.items()}
        r['host_cpu_ms'] = (r['prepare_cpu_ms'] or 0) + (r['preflight_cpu_ms'] or 0) + (r['publish_cpu_ms'] or 0)
        r['host_cpu_us_per_carrier'] = 1e3 * r['host_cpu_ms'] / m['carriers']
        r['steps_used'] = len(steps); r['run'] = str(d)
        rows[(m['mode'], m['transverse'], m['carriers'])] = r
um = {45: '200um', 90: '100um', 180: '50um'}
order = ['owner', 'v33', 'device', 'lane']
cols = sorted({(t, n) for (_, t, n) in rows}, key=lambda x: (x[1], -x[0]))
items = ['host_cpu_ms', 'host_cpu_us_per_carrier', 'prepare_cpu_ms', 'preflight_cpu_ms', 'prepare_wall_ms', 'preflight_wall_ms', 'iterations',
         'admit_census_ms', 'motion_ms', 'trial_digest_ms', 'support_and_bank_ms', 'field_prepare_ms', 'solve_ms', 'residual_check_ms', 'record_seal_ms',
         'field_seal_ms', 'heat_ms', 'push_upload_ms', 'push_admit_ms', 'push_assemble_ms', 'push_device_ms', 'push_gpu_total_ms',
         'cpu_trial_digest_ms', 'cpu_support_and_bank_ms', 'cpu_field_prepare_ms', 'cpu_finish_ms',
         'ds_census_cpu_ms', 'ds_motion_cpu_ms', 'ds_record_cpu_ms', 'ds_preflight_cpu_ms', 'ds_motion_wait_ms', 'ds_push_gpu_total_ms']
for mode in order:
    if not any(k[0] == mode for k in rows):
        continue
    print('mode', mode, '(median ms per step over the steps after the first; host = thread CPU of prepare + preflight + publish)')
    print('  %-28s' % 'item' + ''.join('%14s' % ('%.0e/%s' % (n, um.get(t, t))) for t, n in cols))
    for it in items:
        vals = [rows.get((mode, t, n), {}).get(it) for t, n in cols]
        if all(v is None for v in vals):
            continue
        print('  %-28s' % it + ''.join('%14s' % ('-' if v is None else ('%.3f' % v if abs(v) < 100 else '%.1f' % v)) for v in vals))
    print()
(out / 'bench-table.json').write_text(json.dumps({'%s/%s/%s' % k: v for k, v in rows.items()}, indent=1) + '\n')
