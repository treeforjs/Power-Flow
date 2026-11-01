#!/usr/bin/python -I
"""Step-size calibration of coordinator v8 runs over one physical window (post-processing only; standard library).

The procedure and output of coordinator v7's tools/calibration.py (pic-water-common-commit-v7-20261031, numpy) in plain
Python, plus the step counts and wall times of every run (speed-up against the fixed runs):

  calibration_v8.py --name NAME --out OUT.json --fixed DT=RUN_DIR [DT=RUN_DIR ...] [--controller LABEL=RUN_DIR ...] [--every S]

Three fixed runs dt, dt/2, dt/4 over the same window (the smallest is the reference); observables at common times
t0 + k*every from ledger.bin, probes.bin and cumulative sums of per-step ledger columns. Per observable: largest deviation
from the reference over the comparison times scaled by the largest reference magnitude (a probe field component by the
largest magnitude of its field vector at that probe); observed order p from e(dt)/e(dt/2) = 2^p + 1 (errors against the dt/4
reference); for a controller run the equivalent fixed step dt_eq with C (dt_eq^p - dt_ref^p) = e_ctl.
"""
import argparse, array, json, math, statistics, sys
from pathlib import Path


def table(d, binname, schema):
    cols = json.loads((d / schema).read_text())['columns']
    a = array.array('d')
    raw = (d / binname).read_bytes()
    n = len(cols)
    a.frombytes(raw[:(len(raw) // (8 * n)) * 8 * n])
    rows = [list(a[i * n:(i + 1) * n]) for i in range(len(a) // n)]
    return cols, rows


CUMULATIVE = ['E_joule_J', 'E_circuit_work_J', 'E_plasma_work_J', 'E_thermal_delta_J', 'E_wall_heat_J', 'E_field_delta_J',
              'E_outer_loss_J', 'E_collected_kinetic_J', 'E_escaped_kinetic_J']
STATE = ['E_census_kinetic_J', 'q_collected_cum_C', 'q_escaped_cum_C', 'arr_cum_expected_electrons', 'cov_applied_expected_events_cum',
         'cov_unresolved_upper_events_cum', 'q_free_gross_C']
COUNTS = ['carriers_electron', 'carriers_ion', 'carriers_h2o', 'arr_cum_electron_births', 'arr_cum_h2o_births']
STEP_COUNTS = ['parent_events', 'sinks']


def cumsum(v):
    out, s = [0.], 0.
    for x in v:
        s += x
        out.append(s)
    return out


def load(d):
    d = Path(d)
    lc, L = table(d, 'ledger.bin', 'ledger-schema.json')
    pc, Pr = table(d, 'probes.bin', 'probes-schema.json')
    cc, C = table(d, 'control-ledger.bin', 'control-schema.json')
    li = {c: i for i, c in enumerate(lc)}
    summ = json.loads((d / 'run-summary.json').read_text())
    run = {'dir': str(d), 'steps': len(L), 'time': [Pr[0][1]] + [r[li['time_s']] for r in L], 'dt': [r[li['dt_s']] for r in L],
           'wall_s': summ.get('wall_s'), 'step_ms_median': summ.get('step_ms', {}).get('median')}
    obs = {}
    for c in CUMULATIVE + STEP_COUNTS:
        obs['cum:' + c] = cumsum(r[li[c]] for r in L)
    for c in STATE + COUNTS:
        obs[c] = [math.nan] + [r[li[c]] for r in L]
    for i, c in enumerate(pc[2:], start=2):
        obs['probe:' + c] = [r[i] for r in Pr]
    ci = {c: i for i, c in enumerate(cc)}
    run['control'] = {c: [r[ci[c]] for r in C] for c in cc}
    run['criterion_codes'] = json.loads((d / 'control-schema.json').read_text()).get('criterion_codes', {})
    run['obs'] = obs
    return run


def kind(name):
    return 'stochastic' if name in COUNTS or name in ('cum:parent_events', 'cum:sinks') else 'deterministic'


def at_times(run, times, tol):
    idx = []
    for t in times:
        k = min(range(len(run['time'])), key=lambda i: abs(run['time'][i] - t))
        if abs(run['time'][k] - t) > tol:
            raise SystemExit('%s: no state within %g s of t=%r (closest %r)' % (run['dir'], tol, t, run['time'][k]))
        idx.append(k)
    return idx


def nanmax(v):
    v = [x for x in v if not math.isnan(x)]
    return max(v) if v else math.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--name', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--fixed', nargs='+', required=True)
    ap.add_argument('--controller', nargs='*', default=[])
    ap.add_argument('--every', type=float, default=0)
    a = ap.parse_args()
    fixed = sorted(((float(x.split('=', 1)[0]), load(x.split('=', 1)[1])) for x in a.fixed), key=lambda t: -t[0])
    if len(fixed) != 3 or not all(abs(fixed[i][0] / fixed[i + 1][0] - 2) < 1e-12 for i in range(2)):
        raise SystemExit('need fixed runs at dt, dt/2, dt/4')
    ctrls = [(x.split('=', 1)[0], load(x.split('=', 1)[1])) for x in a.controller]
    ref = fixed[-1][1]
    t0 = ref['time'][0]
    window = min(r['time'][-1] for _, r in fixed + ctrls) - t0
    every = a.every if a.every > 0 else window
    nk = int(math.floor(window / every * (1 + 1e-9)))
    times = [t0 + every * k for k in range(1, nk + 1)]
    tol = 1e-3 * fixed[-1][0]
    runs = [('fixed-%g' % dt, dt, r) for dt, r in fixed] + [(lab, None, r) for lab, r in ctrls]
    sel = {lab: at_times(r, times, tol) for lab, _, r in runs}
    reflab = 'fixed-%g' % fixed[-1][0]
    out = {'name': a.name, 'window_s': window, 'comparison_times': len(times), 'every_s': every, 'reference': reflab, 'runs': {}, 'observables': {}}
    f10 = fixed[0][1]
    for lab, dt, r in runs:
        c = r['control']
        names = {v: k for k, v in r['criterion_codes'].items()}
        crit = {}
        for x in c['criterion']:
            crit[names.get(int(x), str(int(x)))] = crit.get(names.get(int(x), str(int(x))), 0) + 1
        lte = [x for x, v in zip(c['lte_pred_rel'], c['lte_valid']) if v > 0]
        out['runs'][lab] = {'dir': r['dir'], 'steps': r['steps'], 'dt_min_s': min(r['dt']), 'dt_max_s': max(r['dt']), 'dt_mean_s': statistics.fmean(r['dt']),
                            'wall_s': r['wall_s'], 'step_ms_median': r['step_ms_median'], 'refused_attempts': int(sum(x - 1 for x in c['attempts'])),
                            'criteria': crit, 'steps_saved_vs_largest_fixed': f10['steps'] / r['steps'],
                            'wall_speedup_vs_largest_fixed': (f10['wall_s'] / r['wall_s']) if r['wall_s'] else None,
                            'max_omega_p_dt': max(c['omega_p_dt']), 'max_omega_c_dt': max(c['omega_c_dt']),
                            'max_measured_gyro': max(c['measured_gyro']), 'max_displacement_pred_cells': max(c['displacement_pred_cells']),
                            'max_measured_displacement_cells': max(c['measured_displacement_cells']),
                            'median_lte_pred_rel': statistics.median(lte) if lte else None}
    dts = [dt for dt, _ in fixed]
    for name in ref['obs']:
        qref = [ref['obs'][name][i] for i in sel[reflab]]
        if all(math.isnan(x) for x in qref):
            continue
        scale = nanmax([abs(x) for x in qref])
        if name.startswith('probe:') and name.count(':') == 2:
            pt, comp = name.split(':')[1], name.split(':')[2]
            vec = [math.sqrt(sum(ref['obs']['probe:%s:%s%s' % (pt, comp[0], x)][i] ** 2 for x in 'xyz')) for i in sel[reflab]]
            scale = nanmax(vec)
        rec = {'kind': kind(name), 'reference_end': qref[-1], 'scale': scale, 'abs_error_max': {}, 'rel_error_max': {}, 'end_values': {}}
        for lab, dt, r in runs:
            q = [r['obs'][name][i] for i in sel[lab]]
            e = nanmax([abs(x - y) for x, y in zip(q, qref)])
            rec['abs_error_max'][lab] = e
            rec['rel_error_max'][lab] = e / scale if scale > 0 else (0. if e == 0 else math.inf)
            rec['end_values'][lab] = q[-1]
        if rec['kind'] == 'deterministic':
            e1, e2 = rec['abs_error_max']['fixed-%g' % dts[0]], rec['abs_error_max']['fixed-%g' % dts[1]]
            if e1 > 0 and e2 > 0 and e1 / e2 > 1:
                p = math.log2(max(e1 / e2 - 1, 1e-300))
                rec['order_p'] = p
                if .05 < p < 20:
                    Cc = e2 / (2 ** p - 1)
                    rec['dt_eq_s'] = {lab: dts[2] * (rec['abs_error_max'][lab] / Cc + 1) ** (1 / p) for lab, dt, _ in runs if dt is None}
        out['observables'][name] = rec
    for group, pref in (('probe-field-all', 'probe:'), ('probe-field-gap', 'probe:gap')):
        names = [n for n in out['observables'] if n.startswith(pref) and not n.startswith('probe:W_')]
        if names:
            g = {lab: max(out['observables'][n]['rel_error_max'][lab] for n in names) for lab, _, _ in runs}
            out.setdefault('aggregates', {})[group] = {'entries': len(names), 'max_rel_error': g}
    Path(a.out).write_text(json.dumps(out, indent=1) + '\n')
    print('%s: window %.4g s, %d comparison times, reference %s' % (a.name, window, len(times), reflab))
    labs = [lab for lab, _, _ in runs]
    print('%-40s %-6s %s' % ('observable', 'p', '  '.join('%-11s' % l[-11:] for l in labs)))
    for name, rec in out['observables'].items():
        if rec['kind'] != 'deterministic' or (name.startswith('probe:') and not name.startswith('probe:W_')):
            continue
        p = rec.get('order_p')
        print('%-40s %-6s %s' % (name[:40], '%.2f' % p if p is not None else '-', '  '.join('%-11.2e' % rec['rel_error_max'][l] for l in labs)))
    for gname, g in out.get('aggregates', {}).items():
        print('%-40s %-6s %s' % (gname, '', '  '.join('%-11.2e' % g['max_rel_error'][l] for l in labs)))
    print('steps: ' + '  '.join('%s=%d' % (l, out['runs'][l]['steps']) for l in labs))
    print('dt_mean_fs: ' + '  '.join('%s=%.2f' % (l, 1e15 * out['runs'][l]['dt_mean_s']) for l in labs))
    print('wall_s: ' + '  '.join('%s=%.1f' % (l, out['runs'][l]['wall_s'] or 0) for l in labs))


if __name__ == '__main__':
    sys.exit(main())
