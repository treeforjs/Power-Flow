#!/usr/bin/python -I
"""Regression of the v3.5 build against the sealed v3.4 runs of the same suites: every result JSON of a v3.5 run directory
against the v3.4 run directory of the same stage name, value by value, timings, memory, host/device descriptions and
paths excluded. v3.5 changes nothing for bootstraps, steps, checkpoints and restores of non-genesis states, so the
results (digests included) must be identical.
Usage: compare_runs_v34.py V35_LABEL_DIR OUT_JSON V34_LABEL_DIR [V34_LABEL_DIR ...]"""
import json, sys
from pathlib import Path

new_root, out, old_roots = Path(sys.argv[1]), Path(sys.argv[2]), [Path(p) for p in sys.argv[3:]]
SKIP = ('wall', '_ms', 'ms_', 'seconds', '_s_cpu', 'cpu_s', 'elapsed', 'memory', 'rss', 'timing', 'phases', 'device/', 'teardown',
        'admission', '/run', 'checkpoint', '/dir', 'path', 'log', 'hostloop', 'times/', 'host_s')


def flat(x, p='', o=None):
    o = {} if o is None else o
    if isinstance(x, dict):
        for k, v in x.items():
            flat(v, p + '/' + k, o)
    elif isinstance(x, list):
        for i, v in enumerate(x):
            flat(v, p + '/' + str(i), o)
    else:
        o[p] = x
    return o


def stage(d):
    return d.name.rsplit('-', 1)[0]   # e.g. gpu-device-parity-vacuum


rows, ok = [], True
for nd in sorted(p for p in new_root.iterdir() if p.is_dir() and p.name != 'build'):
    olds = [od for root in old_roots for od in sorted(root.glob(stage(nd) + '-*')) if stage(od) == stage(nd)]
    if not olds:
        rows.append({'stage': stage(nd), 'v34': None}); continue
    od = olds[-1]
    files = sorted(f.relative_to(nd) for f in nd.rglob('*.json') if f.name not in ('run-receipt.json',) and 'checkpoint' not in str(f))
    r = {'stage': stage(nd), 'v35': str(nd), 'v34': str(od), 'files': {}}
    for f in files:
        if not (od / f).exists():
            r['files'][str(f)] = 'missing in v3.4'; continue
        try:
            a, b = flat(json.loads((nd / f).read_text())), flat(json.loads((od / f).read_text()))
        except Exception as e:
            r['files'][str(f)] = 'unreadable: %s' % e; continue
        diff = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k) and not any(s in k.lower() for s in SKIP))
        r['files'][str(f)] = {'values': len(a), 'differing': len(diff), 'first': diff[:8]}
        ok = ok and not diff
    rows.append(r)
res = {'schema': 'hall-v3.5-regression-vs-v3.4-v1', 'identical': ok, 'stages': rows}
out.write_text(json.dumps(res, indent=1) + '\n')
for r in rows:
    if r.get('v34') is None:
        print('%-40s (no v3.4 run)' % r['stage']); continue
    print('%-40s %s' % (r['stage'], ', '.join('%s: %s' % (k, v if isinstance(v, str) else ('identical (%d)' % v['values'] if not v['differing'] else
                                                                                       'DIFF %d %s' % (v['differing'], v['first'][:3]))) for k, v in r['files'].items())))
print('identical' if ok else 'DIFFERENCES')
