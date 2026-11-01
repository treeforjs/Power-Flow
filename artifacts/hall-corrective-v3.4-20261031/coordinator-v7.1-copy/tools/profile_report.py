#!/usr/bin/env python3
"""Host profile of coupled steps (post-processing only).

  profile_report.py --out OUT.json [--text OUT.txt] RUN_DIR [RUN_DIR ...]

Each run directory is a load-fixture run (run-summary.json, ledger.bin, subphase-timings.bin,
coordinator-timings.bin). Per run: the median over the steps after the first of every timing item (ms), grouped by
owner. Across runs with the same mesh and host-work mode: a least-squares fit t = a + b N of every item over the
carrier counts N (b in us per carrier); across meshes at the same N and mode: the ratio 50 um / 200 um.

Owner attribution of the phases (each phase is the wall time of the coordinator call that wraps the owner call):
  Hall        hall admit census, motion (host part), deposit support, field prepare (host), field solve, residual check,
              v2 seal, field seal, heat; preflight_final (can_publish Hall)
  native      native begin, birth credit (scope prepare + transport fork; incl. the arrival credit adapter), census
              (native census rows), motion admission, seal (native seal + chemistry seal), can_publish native, swap window
  source      source field, quotes, emission exposures, emission births, source trials, can_publish sources
  arrival     arrival plan, can_publish arrival, post-commit (arrival ledger commit and owner digests)
  chemistry   first and second exposure, can_publish chemistry
  coordinator observe (prepare and preflight), predicate, binding checks (midpoint digest), controller inputs and
              decision, ledger totals, ledger row and I/O, digest line, controller state, control row, probes, receipt, frames
  device      GPU field solve, GPU push device run (and the push's host upload/admit/assemble, reported separately)
"""
import argparse, json, statistics, struct
from pathlib import Path

GROUPS = {
    'hall': ['t_hall_admit_census_ms', 't_hall_motion_ms', 't_hall_deposit_support_ms', 't_hall_field_prepare_ms', 't_hall_field_solve_ms',
             't_hall_residual_check_ms', 't_hall_v2_seal_ms', 't_hall_field_seal_ms', 't_hall_heat_ms', 'sub:t_cp_hall_ms'],
    'native': ['t_native_begin_ms', 't_birth_credit_ms', 't_census_ms', 't_native_motion_admission_ms', 't_seal_ms', 'sub:t_cp_native_ms',
               't_swap_window_ms'],
    'source': ['t_source_field_ms', 't_source_quotes_ms', 't_emission_exposure_ms', 't_emission_births_ms', 't_source_trials_ms',
               'sub:t_cp_sources_ms'],
    'arrival': ['t_arrival_plan_ms', 'sub:t_cp_arrival_ms', 't_post_commit_ms'],
    'chemistry': ['t_first_exposure_ms', 't_second_exposure_ms', 'sub:t_cp_chemistry_ms'],
    'coordinator': ['coord:t_observe_in_prepare_ms', 'coord:t_observe_in_preflight_ms', 'sub:t_cp_predicate_ms', 't_binding_checks_ms',
                    'coord:t_control_inputs_and_decision_ms', 'coord:t_totals_ms', 'coord:t_ledger_row_ms', 'coord:t_ledger_io_ms',
                    'coord:t_digest_line_ms', 'coord:t_controller_state_ms', 'coord:t_control_row_ms', 'coord:t_probes_ms', 'coord:t_receipt_ms',
                    'coord:t_frames_ms'],
    'device': ['coord:t_device_field_solve_ms', 'coord:t_device_push_run_ms'],
    'push host parts': ['coord:t_device_push_upload_ms', 'coord:t_device_push_admit_ms', 'coord:t_device_push_assemble_ms'],
}


def table(d, b, s):
    cols = json.loads((d / s).read_text())['columns']
    raw = (d / b).read_bytes()
    n = len(cols)
    rows = [struct.unpack_from('<%dd' % n, raw, 8 * n * r) for r in range(len(raw) // (8 * n))]
    return {c: [r[i] for r in rows] for i, c in enumerate(cols)}


def load(d):
    d = Path(d)
    s = json.loads((d / 'run-summary.json').read_text())
    L = table(d, 'ledger.bin', 'ledger-schema.json')
    S = table(d, 'subphase-timings.bin', 'subphase-schema.json')
    C = table(d, 'coordinator-timings.bin', 'coordinator-timings-schema.json')
    items = {}
    for k, v in L.items():
        if k.startswith('t_'): items[k] = v
    for k, v in S.items():
        if k.startswith('t_'): items['sub:' + k] = v
    for k, v in C.items():
        if k.startswith('t_'): items['coord:' + k] = v
    # full host step: prepare + publication (t_step_total_ms) + the coordinator's work outside them (decision before
    # prepare; totals, ledger, digests, controller state, probes, receipt and frames after publication)
    outside = ['coord:t_control_inputs_and_decision_ms', 'coord:t_totals_ms', 'coord:t_ledger_row_ms', 'coord:t_ledger_io_ms',
               'coord:t_digest_line_ms', 'coord:t_controller_state_ms', 'coord:t_control_row_ms', 'coord:t_probes_ms', 'coord:t_receipt_ms',
               'coord:t_frames_ms']
    nsteps = len(L['t_step_total_ms'])
    items['step_full_ms'] = [L['t_step_total_ms'][i] + sum(items[k][i] for k in outside if k in items) for i in range(nsteps)]
    med = {k: statistics.median(v[1:] if len(v) > 1 else v) for k, v in items.items()}
    spec = (s.get('mesh') or {}).get('hall_mesh_spec') or {'transverse': 45}
    rec = {'dir': str(d), 'carriers': int(statistics.median(C['carriers'])), 'mesh_um': round(9000 / spec['transverse']),
           'host_work': s.get('host_work'), 'steps': s['steps_completed'], 'stop': s['stop'], 'rss_kib_final': s.get('rss_kib_final'),
           'genesis_s': s.get('genesis_s'), 'step_ms': med['step_full_ms'], 'items_ms': med,
           'owners_ms': {g: sum(med.get(k, 0.) for k in ks) for g, ks in GROUPS.items()}}
    rec['owners_ms']['unattributed'] = rec['step_ms'] - sum(v for g, v in rec['owners_ms'].items() if g not in ('device', 'push host parts'))
    return rec


def fit(points):
    n = len(points)
    if n < 2: return None
    mx = sum(x for x, _ in points) / n
    my = sum(y for _, y in points) / n
    sxx = sum((x - mx) ** 2 for x, _ in points)
    if sxx == 0: return None
    b = sum((x - mx) * (y - my) for x, y in points) / sxx
    return {'a_ms': my - b * mx, 'b_us_per_carrier': 1e3 * b}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--text')
    ap.add_argument('runs', nargs='+')
    a = ap.parse_args()
    runs = [load(d) for d in a.runs]
    out = {'runs': runs, 'fits': {}, 'mesh_ratio': {}}
    keys = sorted({k for r in runs for k in r['items_ms']} | {'owner:' + g for g in GROUPS} | {'step'})
    def value(r, k):
        return r['step_ms'] if k == 'step' else r['owners_ms'][k[6:]] if k.startswith('owner:') else r['items_ms'].get(k, 0.)
    for mesh in sorted({r['mesh_um'] for r in runs}):
        for mode in sorted({r['host_work'] for r in runs}):
            sel = [r for r in runs if r['mesh_um'] == mesh and r['host_work'] == mode]
            out['fits']['%dum/%s' % (mesh, mode)] = {k: fit([(r['carriers'], value(r, k)) for r in sel]) for k in keys}
    for n in sorted({r['carriers'] for r in runs}):
        for mode in sorted({r['host_work'] for r in runs}):
            x = [r for r in runs if r['carriers'] == n and r['host_work'] == mode]
            m200 = next((r for r in x if r['mesh_um'] == 200), None)
            m50 = next((r for r in x if r['mesh_um'] == 50), None)
            if m200 and m50:
                out['mesh_ratio']['%d/%s' % (n, mode)] = {k: (value(m50, k) / value(m200, k) if value(m200, k) > 0 else None) for k in keys}
    Path(a.out).write_text(json.dumps(out, indent=1) + '\n')
    lines = []
    head = '%-34s' % 'item (median ms per step)' + ''.join('%14s' % ('%g/%dum/%s' % (r['carriers'], r['mesh_um'], r['host_work'])) for r in runs)
    lines.append(head)
    lines.append('%-34s' % 'host step (full)' + ''.join('%14.1f' % r['step_ms'] for r in runs))
    for g in list(GROUPS) + ['unattributed']:
        lines.append('%-34s' % ('  ' + g) + ''.join('%14.1f' % r['owners_ms'][g] for r in runs))
    lines.append('%-34s' % 'host us per carrier (step/N)' + ''.join('%14.2f' % (1e3 * r['step_ms'] / r['carriers']) for r in runs))
    lines.append('')
    for g, ks in GROUPS.items():
        lines.append('[%s]' % g)
        for k in ks:
            lines.append('%-34s' % ('  ' + k.split(':')[-1]) + ''.join('%14.2f' % r['items_ms'].get(k, 0.) for r in runs))
    text = '\n'.join(lines) + '\n'
    if a.text: Path(a.text).write_text(text)
    print(text)


if __name__ == '__main__':
    main()
