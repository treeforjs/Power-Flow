#!/usr/bin/env python3
"""Step-size refusals of the owners against the attempted omega_p dt (post-processing only).

  refusals.py --out OUT.json RUN_DIR [RUN_DIR ...]

Every attempt of every step (step-receipts.jsonl) with omega_p,max of its decision inputs times the attempted dt;
attempts are grouped in bins of omega_p dt; per bin: attempts, refused attempts by refusal text, refusal rate per
attempt with a 95 % upper bound (rule of three when none is refused, else Wilson).
"""
import json, math, sys
from collections import defaultdict
from pathlib import Path

EDGES = [0, .01, .02, .03, .04, .045, .05, .055, .06, .065, .07, .1, 1]


def wilson_upper(k, n, z=1.96):
    if n == 0:
        return None
    if k == 0:
        return 3.0 / n
    p = k / n
    return (p + z * z / (2 * n) + z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / (1 + z * z / n)


def main():
    args = sys.argv[1:]
    out = Path(args[args.index('--out') + 1]); del args[args.index('--out'):args.index('--out') + 2]
    bins = defaultdict(lambda: {'attempts': 0, 'refused': 0, 'texts': defaultdict(int)})
    runs = {}
    for d in map(Path, args):
        n_att = n_ref = 0
        for line in (d / 'step-receipts.jsonl').read_text().splitlines():
            r = json.loads(line)
            if 'decision' in r:
                wp = float.fromhex(r['decision']['inputs']['omega_p_max'])
            else:
                continue
            for a in r['attempts']:
                x = wp * float.fromhex(a['dt_s'])
                b = next(i for i in range(len(EDGES) - 1) if EDGES[i] <= x < EDGES[i + 1])
                key = '%.3f-%.3f' % (EDGES[b], EDGES[b + 1])
                bins[key]['attempts'] += 1
                n_att += 1
                if a['outcome'] == 'refused':
                    bins[key]['refused'] += 1
                    bins[key]['texts'][a['refusal'][:90]] += 1
                    n_ref += 1
        runs[str(d)] = {'attempts': n_att, 'refused': n_ref}
    res = {}
    for k in sorted(bins):
        v = bins[k]
        res[k] = {'attempts': v['attempts'], 'refused': v['refused'], 'rate': v['refused'] / v['attempts'],
                  'rate_upper95': wilson_upper(v['refused'], v['attempts']), 'refusals': dict(v['texts'])}
    o = {'runs': runs, 'omega_p_dt_bins': res}
    out.write_text(json.dumps(o, indent=1) + '\n')
    for k, v in res.items():
        print('%-12s attempts %5d refused %3d rate %.4f (<= %.4f)' % (k, v['attempts'], v['refused'], v['rate'], v['rate_upper95']))


if __name__ == '__main__':
    sys.exit(main())
