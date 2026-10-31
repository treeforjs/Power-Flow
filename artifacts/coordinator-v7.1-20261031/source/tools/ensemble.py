#!/usr/bin/env python3
"""Seed-ensemble tally of the stochastic counts at a common time (post-processing only).

  ensemble.py --time-after-genesis SECONDS --out OUT.json LABEL=RUN_DIR[,RUN_DIR...] [LABEL=...]

For every label (one step length), sums over its runs (different seeds) of ionization events, expected events
(coverage), electron and H2O births and sinks up to the given time; reports each label's totals, the ratio
events/expected with its Poisson standard error, and the difference of the per-expected rates between labels.
"""
import json, math, struct, sys
from pathlib import Path


def ledger(d):
    cols = json.loads((d / 'ledger-schema.json').read_text())['columns']
    raw = (d / 'ledger.bin').read_bytes()
    n = len(cols)
    rows = [struct.unpack_from('<%dd' % n, raw, 8 * n * r) for r in range(len(raw) // (8 * n))]
    return {c: [r[i] for r in rows] for i, c in enumerate(cols)}


def main():
    a = sys.argv[1:]
    t = float(a[a.index('--time-after-genesis') + 1]); del a[a.index('--time-after-genesis'):a.index('--time-after-genesis') + 2]
    out = Path(a[a.index('--out') + 1]); del a[a.index('--out'):a.index('--out') + 2]
    res = {}
    for spec in a:
        label, dirs = spec.split('=', 1)
        tot = {'runs': 0, 'events': 0, 'expected_events': 0., 'electron_births': 0, 'h2o_births': 0, 'sinks': 0}
        for d in map(Path, dirs.split(',')):
            L = ledger(d)
            t0 = L['time_s'][0] - L['dt_s'][0]
            k = min(range(len(L['time_s'])), key=lambda i: abs(L['time_s'][i] - (t0 + t)))
            if abs(L['time_s'][k] - (t0 + t)) > 1e-3 * L['dt_s'][k]:
                raise SystemExit('%s: no row at t0 + %g' % (d, t))
            tot['runs'] += 1
            tot['events'] += int(sum(L['parent_events'][:k + 1]))
            tot['expected_events'] += L['cov_applied_expected_events_cum'][k]
            tot['electron_births'] += int(L['arr_cum_electron_births'][k])
            tot['h2o_births'] += int(L['arr_cum_h2o_births'][k])
            tot['sinks'] += int(sum(L['sinks'][:k + 1]))
        tot['events_over_expected'] = tot['events'] / tot['expected_events']
        tot['standard_error'] = math.sqrt(tot['events']) / tot['expected_events']
        res[label] = tot
    labels = list(res)
    if len(labels) == 2:
        x, y = res[labels[0]], res[labels[1]]
        diff = x['events_over_expected'] - y['events_over_expected']
        res['difference'] = {'labels': labels, 'events_over_expected': diff,
                             'sigma': diff / math.hypot(x['standard_error'], y['standard_error'])}
    out.write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    sys.exit(main())
