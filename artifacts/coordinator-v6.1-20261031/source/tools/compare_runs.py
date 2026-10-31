#!/usr/bin/env python3
"""Compare two coordinator v6 run directories.

  compare_runs.py identity A B OUT.json   two fresh processes, same backends: byte identity of the
                                           deterministic per-step outputs (ledger without timings,
                                           owner digests per step, genesis record)
  compare_runs.py backends CPU GPU OUT.json CPU vs GPU backends: discrete per-step outcomes equal,
                                           GPU push parity, field accepted through its guard receipt
                                           (independent residual, absolute Gauss, divB, energy closure)
                                           and the CPU/GPU difference of the continuous ledger columns
"""
import hashlib, json, struct, sys
from pathlib import Path

DISCRETE = ['carriers_electron', 'carriers_ion', 'carriers_h2o', 'births_electron', 'births_h2o', 'sinks', 'parent_events',
            'arr_cum_electron_births', 'arr_cum_h2o_births', 'census_rows_step', 'newborns_step', 'h2o_total_exactly_conserved',
            'nuclei_exactly_conserved', 'hall_iterations']


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def ledger(d):
    d = Path(d)
    cols = json.loads((d / 'ledger-schema.json').read_text())['columns']
    raw = (d / 'ledger.bin').read_bytes()
    n = len(cols)
    rows = [struct.unpack_from('<%dd' % n, raw, k * 8 * n) for k in range(len(raw) // (8 * n))]
    return cols, rows


def identity(a, b):
    out = {'mode': 'identity', 'a': str(a), 'b': str(b), 'files': {}}
    ok = True
    for f in ('deterministic.bin', 'step-digests.txt', 'genesis.json', 'ledger-schema.json'):
        ha, hb = sha(Path(a) / f), sha(Path(b) / f)
        out['files'][f] = {'a_sha256': ha, 'b_sha256': hb, 'identical': ha == hb}
        ok = ok and ha == hb
    sa = json.loads((Path(a) / 'run-summary.json').read_text())
    sb = json.loads((Path(b) / 'run-summary.json').read_text())
    out['steps'] = [sa['steps_completed'], sb['steps_completed']]
    out['byte_identical'] = ok and sa['steps_completed'] == sb['steps_completed']
    return out


def backends(c, g):
    cc, rc = ledger(c)
    cg, rg = ledger(g)
    assert cc == cg, 'ledger schemas differ'
    idx = {n: i for i, n in enumerate(cc)}
    out = {'mode': 'backends', 'cpu': str(c), 'gpu': str(g), 'steps': [len(rc), len(rg)]}
    out['discrete_columns_identical'] = all(rc[k][idx[n]] == rg[k][idx[n]] for k in range(min(len(rc), len(rg))) for n in DISCRETE) and len(rc) == len(rg)
    rel = {}
    for n in cc:
        if n.startswith('t_') or n in ('alloc_swap_window', 'rss_kib') or n in DISCRETE:
            continue
        worst = 0.0
        for a, b in zip(rc, rg):
            x, y = a[idx[n]], b[idx[n]]
            if x != x and y != y:
                continue
            s = max(abs(x), abs(y))
            if s > 0:
                worst = max(worst, abs(x - y) / s)
        rel[n] = worst
    out['max_relative_difference_by_column'] = rel
    out['columns_bit_identical'] = sorted(n for n, v in rel.items() if v == 0)
    out['columns_differing'] = {n: v for n, v in sorted(rel.items(), key=lambda kv: -kv[1]) if v > 0}
    guard = {}
    for name, rows in (('cpu', rc), ('gpu', rg)):
        guard[name] = {
            'max_residual_over_target': max(r[idx['solve_residual_over_target']] for r in rows),
            'max_gauss_over_bound': max(r[idx['gauss_max_abs_C']] / r[idx['gauss_bound_C']] for r in rows),
            'max_divB_over_bound': max(r[idx['divB_max_T_m']] / r[idx['divB_bound_T_m']] for r in rows),
            'max_abs_step_energy_residual_J': max(abs(r[idx['E_step_residual_J']]) for r in rows),
            'roundoff_bound_J': max(r[idx['E_roundoff_bound_J']] for r in rows)}
    out['guard_receipt'] = guard
    sg = json.loads((Path(g) / 'run-summary.json').read_text())
    out['gpu_push_parity'] = sg.get('push_parity')
    out['gpu_field_solves'] = sg.get('gpu_field_solves')
    gr = guard['gpu']
    out['gpu_field_accepted_through_guards'] = gr['max_residual_over_target'] <= 1 and gr['max_gauss_over_bound'] <= 1 and \
        gr['max_divB_over_bound'] <= 1 and gr['max_abs_step_energy_residual_J'] <= gr['roundoff_bound_J']
    p = out['gpu_push_parity'] or {}
    out['push_bit_identical_every_iterate'] = p.get('iterates_compared', 0) > 0 and p.get('iterates_compared') == p.get('iterates_identical')
    return out


if __name__ == '__main__':
    mode, a, b, dest = sys.argv[1:5]
    r = identity(a, b) if mode == 'identity' else backends(a, b)
    Path(dest).write_text(json.dumps(r, indent=1) + '\n')
    summary = {k: v for k, v in r.items() if k not in ('max_relative_difference_by_column', 'columns_bit_identical', 'columns_differing')}
    print(json.dumps(summary, indent=1))
    if mode == 'backends':
        print('columns differing (max rel):', json.dumps(dict(list(r['columns_differing'].items())[:12]), indent=1))
