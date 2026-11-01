#!/usr/bin/env python3
"""Checkpointed conditional mechanics driven by observed actual-pilot T bytes.

This is a prescribed-temperature comparison, not native caloric coupling. It
does not alter the pilot, infer supports, reflect a foil, or authorize fields.
Every accepted stage contains all four trials (two independent bodies, two
conditional supports); rejection preserves the last complete stage.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
BOUNDARIES = ('free', 'both_end_clamped')
BODIES = ('lower', 'upper')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def json_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def write_json(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_bytes(json_bytes(value) + b'\n')
    temp.replace(path)


def load_bound_inputs(audit_path):
    """Verify the copied bytes, certificate binding and separately named density."""
    audit = read(audit_path)
    if (audit.get('schema') != 'actual-native-material-feedback-prefix-audit-v4'
            or not audit.get('prescribed_temperature_geometry_comparison_ready')
            or not audit.get('canonical_order_verified')):
        raise ValueError('verified actual-pilot canonical prefix required')
    snapshot = Path(audit['snapshot'])
    manifest_path = snapshot / 'snapshot-manifest.json'
    if sha(manifest_path) != audit['snapshot_manifest_sha256']:
        raise ValueError('snapshot manifest hash mismatch')
    manifest = read(manifest_path)
    for name, entry in manifest['files'].items():
        path = snapshot / name
        if not path.resolve().is_relative_to(snapshot.resolve()) or sha(path) != entry['sha256']:
            raise ValueError('snapshot artifact hash mismatch')
    prefix_path = Path(audit['canonical_supported_prefix_path'])
    if sha(prefix_path) != audit['canonical_supported_prefix_sha256']:
        raise ValueError('canonical prefix hash mismatch')
    prefix = read(prefix_path)
    if (prefix.get('schema') != 'actual-pilot-prescribed-temperature-prefix-v4'
            or prefix['native_snapshot_manifest_sha256'] != sha(manifest_path)
            or prefix['material_density_scenario_sha256'] != audit['conditional_density_scenario_sha256']):
        raise ValueError('canonical prefix identity mismatch')
    scenario_path = Path(audit['conditional_density_scenario_path'])
    scenario = read(scenario_path)
    # The audit binds canonical scenario JSON, rather than pretty-print bytes.
    scenario_hash = hashlib.sha256(json_bytes(scenario)).hexdigest()
    if scenario_hash != audit['conditional_density_scenario_sha256']:
        raise ValueError('conditional scenario hash mismatch')
    baseline_path = ROOT / 'data/materials/ss304-v1.json'
    baseline = read(baseline_path)
    if (scenario.get('record_id') != 'ss304-pilot8000-density-scenario-v4'
            or scenario['conditional_density_scenario']['baseline_record_sha256'] != sha(baseline_path)):
        raise ValueError('distinct pilot-density scenario/baseline binding required')
    for key in ('thermal_expansion', 'young_modulus', 'proof_stress_screening', 'valid_temperature_K'):
        if scenario[key] != baseline[key]:
            raise ValueError('pilot density scenario may not silently replace property tables')
    if prefix['geometry']['reference_temperature_K'] != 300. or prefix['geometry']['reference_density_kg_m3'] != 8000.:
        raise ValueError('actual pilot reference convention changed')
    if len(prefix['states']) < 2:
        raise ValueError('actual accepted temperature history is incomplete')
    previous = None
    for state in prefix['states']:
        if (type(state['revision']) is not int or state['revision'] < 0
                or not isinstance(state['time_s'], (int, float)) or not math.isfinite(state['time_s']) or state['time_s'] < 0
                or (previous and (state['revision'] != previous['revision'] + 1 or state['time_s'] <= previous['time_s']))
                or state['Tmin_K'] < 298.2 or state['Tmax_K'] >= 373.15
                or state['max_phase'] != 0 or state['max_liquid_fraction'] != 0):
            raise ValueError('prefix clocks, phase, or open first coefficient segment violated')
        previous = state
    if prefix['states'][0]['revision'] != 0:
        raise ValueError('initial accepted state absent; prehistory may not be invented')
    return audit, prefix, scenario


def cell_temperature(nodal, nz):
    """Exact Q1 reference-cell average of T; no enthalpy-conservation claim."""
    import numpy as np
    t = np.asarray(nodal, dtype=float)
    if t.ndim != 2 or min(t.shape) < 2 or not np.isfinite(t).all() or t.min() < 298.2 or t.max() >= 373.15:
        raise ValueError('finite complete nodal first-segment temperature required')
    corners = np.stack((t[:-1, :-1], t[:-1, 1:], t[1:, :-1], t[1:, 1:]))
    return np.tile(corners.mean(axis=0).ravel(), nz), np.tile(corners.max(axis=0).ravel(), nz)


def resource_guard():
    if os.sched_getaffinity(0) != {7} or os.getpriority(os.PRIO_PROCESS, 0) < 15:
        raise ValueError('CPU7 and nice>=15 are required')
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        if os.environ.get(name) != '1':
            raise ValueError('single-thread frozen numerical environment required')


class ChunkBudget(Exception):
    pass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--nz', type=int, choices=(2, 4), default=2)
    parser.add_argument('--max-stages', type=int, default=8)
    parser.add_argument('--max-seconds', type=float, default=13.)
    parser.add_argument('--plot-only', action='store_true')
    args = parser.parse_args()
    if not 1 <= args.max_stages <= 8 or not 0 < args.max_seconds <= 14:
        raise ValueError('one chunk is bounded by 8 stages and <15 seconds')
    resource_guard()
    started = time.monotonic()
    def budget(signum, frame):
        raise ChunkBudget()
    signal.signal(signal.SIGALRM, budget)
    signal.setitimer(signal.ITIMER_REAL, args.max_seconds)
    audit, prefix, record = load_bound_inputs(args.audit)
    import numpy as np
    from mhdlab.solid_foil_coupling import von_mises_Pa
    from mhdlab.solid_foil_geometry import deform_hex_mesh, midplane_gap
    from mhdlab.solid_foil_mechanics import structured_hex_mesh, transient_foil, mesh_geometry_sha256
    from mhdlab.thermal_expansion import SolidExpansionMaterial
    expansion, modulus, density = (record[k] for k in ('thermal_expansion', 'young_modulus', 'density_reference'))
    material = SolidExpansionMaterial(np.array(expansion['temperature_K']), np.array(expansion['relative_length_change']),
        np.array(modulus['temperature_K']), np.array(modulus['value_Pa']), density['value_kg_m3'], density['temperature_K'],
        'ss304-pilot8000-density-scenario-v4;sha256=' + audit['conditional_density_scenario_sha256'])
    rho = material.density_at_stress_free_temperature(300.)
    if not np.isclose(rho, 8000., rtol=1e-14, atol=0.):
        raise ValueError('density scenario does not preserve actual pilot reference mass')
    g = prefix['geometry']
    z = np.linspace(-g['invariant_length_m'] / 2, g['invariant_length_m'] / 2, args.nz + 1)
    meshes = {body: structured_hex_mesh(g[body]['x_edges_m'], g[body]['y_edges_m'], z) for body in BODIES}
    states = prefix['states']
    for state in states:
        all_values = []
        for body in BODIES:
            t = np.asarray(state[body]['nodal_temperature_K'])
            if t.shape != (len(g[body]['y_edges_m']), len(g[body]['x_edges_m'])):
                raise ValueError('canonical temperature/geometry shape mismatch')
            cell_temperature(t, args.nz)
            all_values.extend(t.ravel())
        if min(all_values) != state['Tmin_K'] or max(all_values) != state['Tmax_K']:
            raise ValueError('canonical temperature extrema mismatch')
    if any(np.any(np.asarray(states[0][body]['nodal_temperature_K']) != 300.) for body in BODIES):
        raise ValueError('initial stress-free temperature state must be exactly300K')
    code_files = ('tools/run_material_feedback_geometry_v4.py', 'mhdlab/solid_foil_mechanics.py',
                  'mhdlab/solid_foil_geometry.py', 'mhdlab/thermal_expansion.py', 'mhdlab/solid_foil_coupling.py')
    config = {'schema': 'conditional-actual-pilot-geometry-config-v4', 'audit_sha256': sha(args.audit),
        'prefix_sha256': audit['canonical_supported_prefix_sha256'], 'snapshot_manifest_sha256': audit['snapshot_manifest_sha256'],
        'density_scenario_sha256': audit['conditional_density_scenario_sha256'], 'nz': args.nz,
        'mesh_sha256': {b: mesh_geometry_sha256(m) for b, m in meshes.items()},
        'source_sha256': {name: sha(ROOT / name) for name in code_files}, 'poisson_ratio': .29,
        'poisson_ratio_basis': 'explicit conditional assumption; not a measured foil property',
        'supports': list(BOUNDARIES), 'support_basis': 'unknown physical fixture; free and immobile both-end-plane clamp brackets only',
        'temperature_projection': 'exact reference Q1 corner-average T, then piecewise constant per solid cell; z invariant',
        'native_caloric_ledger_conservation_claimed': False, 'thermal_mechanical_two_way_feedback': False,
        'owner_final_build_run_seal_verified': False,
        'source_scope': 'hash-verified observed accepted checkpoint bytes with native coordinate certificate; formal source-to-executable run seal not established',
        'initial_state': 'zero u/v, stress-free300K at first accepted C19 time; no prehistory reconstructed',
        'physical_coupled_ready': False, 'production_field_handoff_authorized': False,
        'spatial_convergence_claimed': False, 'proof_screen': record['proof_stress_screening'],
        'excluded': ['Lorentz forces', 'plasticity', 'melt', 'fixture reconstruction', 'caloric heat debit', 'EM geometry feedback']}
    config_hash = hashlib.sha256(json_bytes(config)).hexdigest()
    args.output.mkdir(parents=True, exist_ok=True)
    config_path = args.output / 'config.json'
    if config_path.exists():
        if read(config_path) != config:
            raise ValueError('restart configuration/source/input identity changed')
    else:
        write_json(config_path, config)
    current = int((args.output / 'CURRENT').read_text()) if (args.output / 'CURRENT').exists() else 0
    keys = [f'{bc}__{body}' for bc in BOUNDARIES for body in BODIES]
    u = {key: np.zeros(meshes[key.split('__')[1]].nodes_m.shape) for key in keys}
    v = {key: value.copy() for key, value in u.items()}
    if not 0 <= current < len(states):
        raise ValueError('restart revision outside the bound history')
    if current:
        directory = args.output / f'step_{current:06d}'
        receipt = read(directory / 'receipt.json')
        if (receipt['config_sha256'] != config_hash or receipt['revision'] != current
                or receipt['time_s'] != states[current]['time_s'] or sha(directory / 'state.npz') != receipt['state_sha256']):
            raise ValueError('restart state identity/hash mismatch')
        with np.load(directory / 'state.npz', allow_pickle=False) as archive:
            u = {key: archive[key + '__u'].copy() for key in keys}
            v = {key: archive[key + '__v'].copy() for key in keys}
    if args.plot_only:
        plot_results(args.output)
        return
    proof = record['proof_stress_screening']
    status = 'chunk_complete'
    completed = 0
    try:
        if (args.output / 'rejection.json').exists():
            status = 'rejected_trial_preserved_last_accepted'
        else:
            for index in range(current + 1, min(len(states), current + args.max_stages + 1)):
                previous, candidate = states[index - 1], states[index]
                trial_u, trial_v, arrays, body_rows, gaps = {}, {}, {}, {}, {}
                for bc in BOUNDARIES:
                    geometries = {}
                    for body in BODIES:
                        key, mesh = f'{bc}__{body}', meshes[body]
                        oldT, oldmax = cell_temperature(previous[body]['nodal_temperature_K'], args.nz)
                        newT, newmax = cell_temperature(candidate[body]['nodal_temperature_K'], args.nz)
                        result = transient_foil(mesh, [previous['time_s'], candidate['time_s']], np.array([oldT, newT]),
                            material=material, stress_free_temperature_K=300., poisson_ratio=.29, boundary=bc,
                            initial_displacement_m=u[key], initial_velocity_m_s=v[key])
                        vm = von_mises_Pa(result['gauss_stress_Pa'])
                        screen = np.interp(np.array([oldmax, newmax]), proof['temperature_K'], proof['value_Pa'])
                        ratio = float(np.max(vm / screen[..., None]))
                        if ratio > 1:
                            raise ValueError(f'{key} nominal proof-screen exceeded: ratio={ratio:.17g}; no plastic response supported')
                        trial_u[key], trial_v[key] = result['displacement_m'][-1], result['velocity_m_s'][-1]
                        geom = deform_hex_mesh(mesh.nodes_m, mesh.cells, trial_u[key], reference_density_kg_m3=rho,
                            cell_temperature_k=newT, admissible_temperature_range_k=(298.2, 773.15))
                        geometries[body] = geom
                        geom_receipt = geom.receipt()
                        geom_receipt.pop('material_face_ids')
                        body_rows[key] = {'max_displacement_m': float(np.linalg.norm(trial_u[key], axis=1).max()),
                            'max_von_mises_Pa': float(vm[-1].max()), 'max_proof_ratio_both_endpoints': ratio,
                            'proof_temperature_basis': 'maximum native corner T per cell, conservative decreasing proof table',
                            'kinetic_energy_J': float(result['kinetic_energy_J'][-1]),
                            'stress_work_increment_J': float(result['stress_work_increment_J'][0]),
                            'mechanical_balance_residual_J': float(result['mechanical_balance_residual_J'][0]),
                            'geometry': geom_receipt, 'heat_debit_J': None}
                        arrays.update({key + '__u': trial_u[key], key + '__v': trial_v[key],
                            key + '__nodes_m': geom.nodes_m, key + '__cell_density_kg_m3': geom.cell_density_kg_m3})
                    gap = midplane_gap(geometries['lower'], geometries['upper'], z_m=0.)
                    gaps[bc] = {'minimum_gap_m': gap['minimum_gap_m'], 'maximum_gap_m': gap['maximum_gap_m'],
                        'minimum_gap_change_m': gap['minimum_gap_m'] - .005,
                        'maximum_gap_change_m': gap['maximum_gap_m'] - .005,
                        'parallel_plate_max_field_ratio_fixed_voltage': float(np.max(.005 / gap['gap_m'])),
                        'parallel_plate_only_no_EM_solve': True}
                    arrays.update({bc + '__gap_x_m': gap['x_m'], bc + '__gap_m': gap['gap_m']})
                final = args.output / f'step_{index:06d}'
                pending = args.output / f'.step_{index:06d}.pending'
                if pending.exists():
                    shutil.rmtree(pending)  # only this runner's interrupted uncommitted trial
                pending.mkdir()
                np.savez(pending / 'state.npz', **arrays)
                receipt = {'schema': 'conditional-actual-pilot-geometry-stage-v4', 'config_sha256': config_hash,
                    'revision': index, 'parent_revision': index - 1, 'time_s': candidate['time_s'],
                    'accepted_interval_s': [previous['time_s'], candidate['time_s']],
                    'native_temperature_range_K': [candidate['Tmin_K'], candidate['Tmax_K']],
                    'state_sha256': sha(pending / 'state.npz'), 'bodies': body_rows, 'gaps': gaps,
                    'physical_coupled_ready': False, 'heat_debit_J': None}
                write_json(pending / 'receipt.json', receipt)
                if final.exists():
                    raise ValueError('unreferenced completed stage exists; inspect rather than silently replace it')
                # Defer the cooperative alarm across the short commit critical
                # section. The caller's outer timeout still preserves CURRENT.
                signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGALRM})
                try:
                    pending.replace(final)
                    marker = args.output / 'CURRENT.tmp'
                    marker.write_text(str(index) + '\n')
                    marker.replace(args.output / 'CURRENT')
                    current, u, v = index, trial_u, trial_v
                    completed += 1
                finally:
                    signal.pthread_sigmask(signal.SIG_UNBLOCK, {signal.SIGALRM})
                print(json.dumps({'accepted_revision': current, 'time_s': candidate['time_s'], 'gaps': gaps}), flush=True)
    except ChunkBudget:
        status = 'chunk_budget_preserved_last_accepted'
    except ValueError as error:
        status = 'rejected_trial_preserved_last_accepted'
        write_json(args.output / 'rejection.json', {'reason': str(error), 'last_accepted_revision': current,
            'last_accepted_time_s': states[current]['time_s'], 'attempted_revision': current + 1,
            'config_sha256': config_hash, 'physical_coupled_ready': False})
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.)
    if current == len(states) - 1:
        status = 'supported_prefix_complete'
    summary = {'schema': 'conditional-actual-pilot-geometry-summary-v4', 'status': status,
        'accepted_revision': current, 'accepted_time_s': states[current]['time_s'],
        'available_smooth_prefix_last_revision': len(states) - 1, 'stages_accepted_this_chunk': completed,
        'chunk_elapsed_s': time.monotonic() - started, 'config_sha256': config_hash,
        'physical_coupled_ready': False, 'production_field_handoff_authorized': False}
    write_json(args.output / 'summary.json', summary)
    print(json.dumps(summary), flush=True)


def plot_results(output):
    """Separate bounded postprocessing command; no solver is called."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    current = int((output / 'CURRENT').read_text())
    rows = [read(output / f'step_{k:06d}/receipt.json') for k in range(1, current + 1)]
    if not rows:
        raise ValueError('no accepted conditional geometry to plot')
    fig, ax = plt.subplots(1, 2, figsize=(11, 4), facecolor='white', constrained_layout=True)
    times = np.array([row['time_s'] for row in rows]) * 1e9
    for bc, label in (('free', 'Free'), ('both_end_clamped', 'Both end planes clamped')):
        ax[0].plot(times, [row['gaps'][bc]['minimum_gap_change_m'] * 1e9 for row in rows], label=label)
        ax[1].plot(times, [max(row['bodies'][f'{bc}__{b}']['max_von_mises_Pa'] for b in BODIES) / 1e6 for row in rows], label=label)
    ax[0].set_ylabel('Minimum midplane gap change (nm)')
    ax[1].set_ylabel('Maximum sampled von Mises stress (MPa)')
    for axis in ax:
        axis.set_xlabel('Actual pilot time (ns)')
        axis.grid(alpha=.25)
        axis.legend()
    fig.suptitle('Conditional geometry from independent actual-pilot temperatures; nominal 8000 kg/m³ at300K\nNo inferred fixture, caloric feedback, magnetic load or spatial-convergence claim', fontsize=10)
    fig.savefig(output / 'conditional-geometry-history.png', dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
