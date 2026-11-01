#!/usr/bin/env python3
"""Step-size calibration of coordinator v7 runs over one physical window (post-processing only).

  calibration.py --name NAME --out OUT.json --fixed DT=RUN_DIR [DT=RUN_DIR ...] [--controller LABEL=RUN_DIR ...]
                 [--every SECONDS]

The fixed-dt runs must include three step lengths dt, dt/2, dt/4 over the same window (the smallest is the
reference). Observables are taken at common times t0 + k*every (k = 1..; every defaults to the window end) from
the step ledger (ledger.bin), the probes (probes.bin) and cumulative sums of per-step ledger columns:

  deterministic   field at the probe entries (E, B), field energies W_E, W_B, census kinetic energy, cumulative
                  Joule heat, circuit work, plasma work, thermal change, wall heat, charge collected/escaped, expected
                  emission (arrival), expected and unresolved-upper parent events (coverage)
  stochastic      carrier counts, births, parent events, sinks (compared with Poisson sigma; no order)

Per observable: the largest deviation from the reference over the comparison times, scaled by the largest
magnitude of the reference over those times (relative error; a probe field component by the largest magnitude of
its field vector at that probe), for every fixed run and controller run; the
observed order p from e(dt)/e(dt/2) = (4^p - 1)/(2^p - 1) (errors against the dt/4 reference) and the
three-run order p3 = log2(|Q_dt - Q_dt/2| / |Q_dt/2 - Q_dt/4|) at the window end; for a controller run the
equivalent fixed step dt_eq with C (dt_eq^p - dt_ref^p) = e_ctl.
"""
import argparse, json, math, sys
from pathlib import Path
import numpy as np


def table(d, binname, schema):
    cols = json.loads((d / schema).read_text())['columns']
    raw = np.frombuffer((d / binname).read_bytes(), dtype='<f8')
    n = len(cols)
    return cols, raw[: (len(raw) // n) * n].reshape(-1, n)


CUMULATIVE = ['E_joule_J', 'E_circuit_work_J', 'E_plasma_work_J', 'E_thermal_delta_J', 'E_wall_heat_J', 'E_field_delta_J',
              'E_outer_loss_J', 'E_collected_kinetic_J', 'E_escaped_kinetic_J']
STATE = ['E_census_kinetic_J', 'q_collected_cum_C', 'q_escaped_cum_C', 'arr_cum_expected_electrons', 'cov_applied_expected_events_cum',
         'cov_unresolved_upper_events_cum', 'q_free_gross_C']
COUNTS = ['carriers_electron', 'carriers_ion', 'carriers_h2o', 'arr_cum_electron_births', 'arr_cum_h2o_births']
STEP_COUNTS = ['parent_events', 'sinks']


def load(d):
    d = Path(d)
    lc, L = table(d, 'ledger.bin', 'ledger-schema.json')
    pc, P = table(d, 'probes.bin', 'probes-schema.json')
    cc, C = table(d, 'control-ledger.bin', 'control-schema.json')
    li = {c: i for i, c in enumerate(lc)}
    run = {'dir': str(d), 'steps': len(L), 'time': np.concatenate([[P[0, 1]], L[:, li['time_s']]]), 'dt': L[:, li['dt_s']]}
    obs = {}
    for c in CUMULATIVE:
        obs['cum:' + c] = np.concatenate([[0.], np.cumsum(L[:, li[c]])])
    for c in STATE + COUNTS:
        obs[c] = np.concatenate([[np.nan], L[:, li[c]]])
    for c in STEP_COUNTS:
        obs['cum:' + c] = np.concatenate([[0.], np.cumsum(L[:, li[c]])])
    for i, c in enumerate(pc[2:], start=2):
        obs['probe:' + c] = P[:, i]
    ci = {c: i for i, c in enumerate(cc)}
    run['control'] = {c: C[:, ci[c]] for c in cc}
    run['obs'] = obs
    return run


def kind(name):
    if name in COUNTS or name in ('cum:parent_events', 'cum:sinks'):
        return 'stochastic'
    return 'deterministic'


def at_times(run, times, tol):
    idx = []
    for t in times:
        k = int(np.argmin(np.abs(run['time'] - t)))
        if abs(run['time'][k] - t) > tol:
            raise SystemExit('%s: no state within %g s of t=%r (closest %r)' % (run['dir'], tol, t, run['time'][k]))
        idx.append(k)
    return np.array(idx)


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
    times = t0 + every * np.arange(1, nk + 1)
    tol = 1e-3 * fixed[-1][0]
    runs = [('fixed-%g' % dt, dt, r) for dt, r in fixed] + [(lab, None, r) for lab, r in ctrls]
    sel = {lab: at_times(r, times, tol) for lab, _, r in runs}
    out = {'name': a.name, 'window_s': float(window), 'comparison_times': len(times), 'every_s': every, 'reference': 'fixed-%g' % fixed[-1][0],
           'runs': {}, 'observables': {}}
    for lab, dt, r in runs:
        c = r['control']
        out['runs'][lab] = {'dir': r['dir'], 'steps': r['steps'], 'dt_min_s': float(np.min(r['dt'])), 'dt_max_s': float(np.max(r['dt'])),
                            'dt_mean_s': float(np.mean(r['dt'])), 'refused_attempts': int(np.sum(c['attempts'] - 1)),
                            'criteria': {str(int(k)): int(v) for k, v in zip(*np.unique(c['criterion'], return_counts=True))},
                            'max_omega_p_dt': float(np.max(c['omega_p_dt'])), 'max_measured_plasma_dt': float(np.max(c['measured_plasma_dt'])),
                            'max_omega_c_dt': float(np.max(c['omega_c_dt'])), 'max_measured_gyro': float(np.max(c['measured_gyro'])),
                            'max_displacement_pred_cells': float(np.max(c['displacement_pred_cells'])),
                            'max_measured_displacement_cells': float(np.max(c['measured_displacement_cells'])),
                            'median_lte_pred_rel': float(np.median(c['lte_pred_rel'][c['lte_valid'] > 0])) if np.any(c['lte_valid'] > 0) else None}
    dts = [dt for dt, _ in fixed]
    for name in ref['obs']:
        qref = ref['obs'][name][sel[out['reference']]]
        if np.all(np.isnan(qref)):
            continue
        scale = float(np.nanmax(np.abs(qref)))
        if name.startswith('probe:') and name.count(':') == 2:
            # a field component is scaled by the largest magnitude of its field vector at that probe
            pt, comp = name.split(':')[1], name.split(':')[2]
            vec = np.sqrt(sum(ref['obs']['probe:%s:%s%s' % (pt, comp[0], x)][sel[out['reference']]] ** 2 for x in 'xyz'))
            scale = float(np.nanmax(vec))
        rec = {'kind': kind(name), 'reference_end': float(qref[-1]), 'scale': scale, 'abs_error_max': {}, 'rel_error_max': {}, 'end_values': {}}
        for lab, dt, r in runs:
            q = r['obs'][name][sel[lab]]
            e = float(np.nanmax(np.abs(q - qref)))
            rec['abs_error_max'][lab] = e
            rec['rel_error_max'][lab] = e / scale if scale > 0 else (0. if e == 0 else math.inf)
            rec['end_values'][lab] = float(q[-1])
        if rec['kind'] == 'deterministic':
            e1, e2 = rec['abs_error_max']['fixed-%g' % dts[0]], rec['abs_error_max']['fixed-%g' % dts[1]]
            if e1 > 0 and e2 > 0 and e1 / e2 > 1:
                p = math.log2(max(e1 / e2 - 1, 1e-300))
                rec['order_p'] = p
                if .05 < p < 20:
                    # e(dt) = C ((dt/dt_ref)^p - 1) with dt_ref the reference step
                    C = e2 / (2 ** p - 1)
                    rec['dt_eq_s'] = {lab: dts[2] * (rec['abs_error_max'][lab] / C + 1) ** (1 / p) for lab, dt, _ in runs if dt is None}
            v = [rec['end_values']['fixed-%g' % d] for d in dts]
            if v[1] != v[2] and v[0] != v[1]:
                rec['order_p3_end'] = math.log2(abs(v[0] - v[1]) / abs(v[1] - v[2]))
        else:
            for lab, dt, r in runs:
                q = rec['end_values'][lab]
                s = math.sqrt(max(abs(q) + abs(rec['reference_end']), 1.))
                rec.setdefault('z_end', {})[lab] = (q - rec['reference_end']) / s
        out['observables'][name] = rec
    # Aggregate of the probe field entries: largest relative error over all entries (each scaled by its own maximum).
    for group, pref in (('probe-field-all', 'probe:'), ('probe-field-blob', 'probe:blob'), ('probe-field-gap', 'probe:gap')):
        names = [n for n in out['observables'] if n.startswith(pref) and not n.startswith('probe:W_')]
        if names:
            g = {lab: max(out['observables'][n]['rel_error_max'][lab] for n in names) for lab, _, _ in runs}
            out.setdefault('aggregates', {})[group] = {'entries': len(names), 'max_rel_error': g}
    Path(a.out).write_text(json.dumps(out, indent=1) + '\n')
    # short text report
    print('%s: window %.4g s, %d comparison times, reference %s' % (a.name, window, len(times), out['reference']))
    print('%-44s %-6s %s' % ('observable', 'p', '  '.join('%-12s' % lab for lab, _, _ in runs)))
    for name, rec in out['observables'].items():
        if rec['kind'] != 'deterministic' or name.startswith('probe:') and not name.startswith('probe:W_') and 'blob' not in name and 'gap-centre' not in name:
            continue
        p = rec.get('order_p')
        print('%-44s %-6s %s' % (name[:44], '%.2f' % p if p is not None else '-', '  '.join('%-12.3e' % rec['rel_error_max'][lab] for lab, _, _ in runs)))
    for gname, g in out.get('aggregates', {}).items():
        print('%-44s %-6s %s' % (gname, '', '  '.join('%-12.3e' % g['max_rel_error'][lab] for lab, _, _ in runs)))
    print('steps: ' + '  '.join('%s=%d' % (lab, out['runs'][lab]['steps']) for lab, _, _ in runs))


if __name__ == '__main__':
    sys.exit(main())
