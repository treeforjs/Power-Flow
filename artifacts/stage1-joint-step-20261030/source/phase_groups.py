#!/usr/bin/python -I
"""Per-phase wall-time groups for the Stage 1 joint step (P1, one process).

Reads t1-timings.json of each T1 run directory given and writes
phase-groups.json next to them. Measured values come from the real step
(coordinator timers, link-wrapped Hall motion calls, publication probes);
the split of the monolithic Hall proposal into field solve and field-state
hashing uses the probe replays recorded in the same file and is labeled
"derived". Usage: phase_groups.py RUN_DIRECTORY
"""
import json, statistics, sys
from pathlib import Path

run = Path(sys.argv[1])
profile = json.loads(Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/phase-breakdown-summary.json').read_text())
rows = []
for d in sorted(run.glob('t1-*')):
    f = d / 't1-timings.json'
    if not f.exists():
        continue
    t = json.loads(f.read_text())['P1']
    ph, wrap, pub, probe = t['prepare_phases'], t['prepare_wrapped_hall_calls'], t['publish'], t['phase_probes']
    it = ph['hall_iterations']
    motion = wrap['prepare_canonical_motion']['seconds']
    state_hash, admit, model_prepare = probe['model_state_sha256_s'], probe['model_admit_constraints_s'], probe['model_prepare_s']
    solve = it * max(model_prepare - 2 * state_hash - 2 * admit, 0.0)
    chemistry = ph['first_exposure_s'] + ph['second_exposure_s']
    prepare = t['prepare_s']
    can_publish = pub['can_publish_s']
    commit = pub['swap_window_s'] + pub['release_s']
    hashing = prepare - motion - solve - chemistry + can_publish
    # Field-state SHA-256 (JSON dump of E, B, both low words and the charge
    # ledger): 2 per Model::prepare, 2 in heat binding, 1 per observation
    # (prepare and can_publish), 2 in the final Hall preflight.
    field_hashes = 2 * it + 2 + 2 + 2
    rows.append({
        'run': d.name, 'hall_iterations': it, 'census_rows': probe['census_rows'],
        'measured': {'prepare_s': prepare, 'can_publish_s': can_publish, 'swap_window_s': pub['swap_window_s'],
                     'release_s': pub['release_s'], 'step_total_s': prepare + pub['total_s'],
                     'coordinator_phases_s': ph, 'wrapped_hall_calls': wrap,
                     'allocations_prepare': t['prepare_allocations'],
                     'allocations_can_publish': pub['allocations_can_publish'],
                     'allocations_between_first_and_last_swap': pub['allocations_between_first_and_last_swap']},
        'probe_replays_s': probe,
        'groups_ms': {'gather_push': 1e3 * motion, 'field_solve_derived': 1e3 * solve,
                      'hashing_and_validation': 1e3 * hashing, 'chemistry': 1e3 * chemistry, 'commit': 1e3 * commit},
        'field_state_hash': {'count_per_step': field_hashes, 'ms_each_probe': 1e3 * state_hash,
                             'ms_total_derived': 1e3 * field_hashes * state_hash},
    })
groups = list(rows[0]['groups_ms']) if rows else []
summary = {g: statistics.median(r['groups_ms'][g] for r in rows) for g in groups}
total = sum(summary.values())
out = {'schema': 'stage1-phase-groups-v1',
       'definitions': {
           'gather_push': 'link-wrapped hall_driver::prepare_canonical_motion wall time summed over all fixed-point iterations (real step; includes its internal census/mesh/field binary hashing)',
           'field_solve_derived': 'iterations x (Model::prepare - 2 state_sha256 - 2 admit), from probe replays on the post-step state',
           'hashing_and_validation': 'prepare + can_publish - gather_push - field_solve - chemistry (owner observations, JSON field-state hashing, heat binding, native typed admission, seals, final preflights)',
           'chemistry': 'first + second parent exposures (both SSA half-exposures, private native publication)',
           'commit': 'first native swap to end of reservation release (allocation-free window plus release)'},
       'runs': rows, 'median_groups_ms': summary,
       'median_percent': {g: 100 * v / total for g, v in summary.items()} if total else {},
       'electron_only_reference': {k: {'ms_per_step': v['ms_per_step'], 'groups': v['groups']} for k, v in profile['results'].items()},
       'comparability': 'electron-only profile: about 29.5k particles, electrostatic GPU Poisson, 210 steps; Stage 1: 2 carriers, Hall implicit Maxwell, one step. Group sizes are comparable only per phase type.'}
(run / 'phase-groups.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps({'median_groups_ms': summary, 'median_percent': out['median_percent']}, indent=1))
