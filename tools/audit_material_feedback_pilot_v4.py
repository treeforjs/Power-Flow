#!/usr/bin/env python3
"""Read-only admission audit of the actual evolving-material pilot.

No old heating trajectory is substituted and no mechanics is launched. An
optional canonical owner export is audited only when source-bound and sealed.
The export schema is emitted with --export-contract for an owner to implement;
the audit does not create missing states, mirrors, supports or density data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import struct


ROOT = Path(__file__).resolve().parents[1]
PILOT = Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-material-feedback-pilot-20261029')
OUTPUT = Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-expansion-20261029/foil-support-v4/pilot-audit')
EXPORT_SCHEMA = 'material-feedback-mechanics-export-v1'
HEAT_MAGIC_V2 = 0x454C484541543034


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f'{name} must be finite real SI data')
    return float(value)


def matrix(value, ny, nx, name):
    if not isinstance(value, list) or len(value) != ny or any(not isinstance(row, list) or len(row) != nx for row in value):
        raise ValueError(f'{name} needs complete canonical (y,x) data, shape {(ny, nx)}')
    return [[number(item, name) for item in row] for row in value]


def edges(value, name):
    if not isinstance(value, list) or len(value) < 2:
        raise ValueError(f'{name} needs complete reference node axes')
    result = [number(x, name) for x in value]
    if any(b <= a for a, b in zip(result, result[1:])):
        raise ValueError(f'{name} must strictly increase')
    return result


def pilot_density_scenario(material, baseline_sha256, *, density_kg_m3, reference_temperature_K):
    """New owner-nominal density scenario; never mutate the baseline record."""
    if density_kg_m3 != 8000. or reference_temperature_K != 300.:
        raise ValueError('only the explicitly selected pilot8000-at300K density scenario is admitted here')
    scenario = json.loads(json.dumps(material))
    law = material['thermal_expansion']
    temperatures, strain = law['temperature_K'], law['relative_length_change']
    k = next(k for k in range(len(temperatures) - 1) if temperatures[k] <= 300. <= temperatures[k + 1])
    s300 = strain[k] + (strain[k + 1] - strain[k]) * (300. - temperatures[k]) / (temperatures[k + 1] - temperatures[k])
    scenario['record_id'] = 'ss304-pilot8000-density-scenario-v4'
    scenario['density_reference'] = {
        'source': 'explicit native material-feedback pilot nominal density convention, not measured foil density',
        'temperature_K': law['reference_temperature_K'],
        'value_kg_m3': density_kg_m3 * (1 + s300) ** 3,
        'temperature_assignment': 'back-mapped through the unchanged bounded expansion law so rho_ref(300 K)=8000 kg/m3',
        'mass_conservation': 'fixed reference mass at300K; rho_current=rho_ref/det(F_total)',
    }
    scenario['conditional_density_scenario'] = {
        'baseline_record_sha256': baseline_sha256,
        'baseline_density_reference': material['density_reference'],
        'owner_reference_density_kg_m3': density_kg_m3, 'owner_reference_temperature_K': reference_temperature_K,
        'expansion_and_elastic_tables_changed': False, 'baseline_record_changed': False,
        'measured_foil_density': False, 'physical_EOS_ready': False,
        'scope': 'source-T-driven conditional geometry comparison until native caloric and work identities align',
    }
    scenario['material']['qualification'] += ' Density only is a separately identified native-pilot nominal8000-at300K scenario; not a measured or batch-qualified foil density.'
    return scenario


def read_heat_checkpoint_v2(path):
    """Decode only the inspected ELHEAT04 layout; retain native DOF ordering."""
    data = Path(path).read_bytes()
    cursor = 0

    def get(format_):
        nonlocal cursor
        count = struct.calcsize('<' + format_)
        if cursor + count > len(data):
            raise ValueError('truncated native heat checkpoint')
        value = struct.unpack_from('<' + format_, data, cursor)
        cursor += count
        return value[0] if len(value) == 1 else value

    def text():
        nonlocal cursor
        count = get('Q')
        if count > 1000000 or cursor + count > len(data):
            raise ValueError('invalid native checkpoint string')
        value = data[cursor:cursor + count].decode('utf-8')
        cursor += count
        return value

    if get('Q') != HEAT_MAGIC_V2:
        raise ValueError('unrecognized native heat checkpoint ABI; no inferred binary interpretation')
    signature = text()
    words = signature.split()
    if len(words) < 18 or words[0] != 'MASS_LUMPED_DUAL_CELL_BE_V2' or words[11] != 'nominal304_v2_q=2':
        raise ValueError('checkpoint is not the inspected native v2/q2 material-feedback pilot')
    if words[17] != 'solid-only-C19-no-particle-sources/v1' or not words[12].startswith('rho_ref='):
        raise ValueError('unexpected pilot source/density signature')
    geometry = {'nx': int(words[1]), 'ny': int(words[2]), 'x0_m': float.fromhex(words[3]),
        'y0_m': float.fromhex(words[4]), 'width_m': float.fromhex(words[5]),
        'thickness_m': float.fromhex(words[6]), 'invariant_length_m': float.fromhex(words[7]),
        'initial_temperature_K': float.fromhex(words[8]), 'mesh_revision': words[9], 'body': words[10],
        'reference_density_kg_m3': float.fromhex(words[12].split('=', 1)[1]),
        'material_headers_sha256': words[14], 'caloric_callback_contract': words[16]}
    time, revision, count = get('d'), get('Q'), get('Q')
    if count != (geometry['nx'] + 1) * (geometry['ny'] + 1) or count > 1000000:
        raise ValueError('checkpoint native DOF count differs from signed geometry')
    values = [get('d') for _ in range(count)]
    inventory_count = get('Q')
    if inventory_count != 0:
        raise ValueError('solid-only pilot unexpectedly contains a surface inventory')
    cumulative = [get('d') for _ in range(19)]
    watermark_count = get('Q')
    if watermark_count > 100000:
        raise ValueError('invalid checkpoint watermark count')
    watermarks = {}
    for _ in range(watermark_count):
        key, sequence = text(), get('Q')
        if not key or not sequence or key in watermarks:
            raise ValueError('invalid checkpoint watermark')
        watermarks[key] = sequence
    if cursor != len(data) or not all(math.isfinite(x) for x in [time, *values, *cumulative]):
        raise ValueError('trailing or nonfinite native heat checkpoint data')
    return {'geometry': geometry, 'signature_sha256': hashlib.sha256(signature.encode()).hexdigest(),
        'time_s': time, 'revision': revision, 'native_temperature_dofs_K': values,
        'cumulative': cumulative, 'watermarks': watermarks, 'file_sha256': hashlib.sha256(data).hexdigest()}


def canonical_temperature(checkpoint, certificate):
    """Apply the native MFEM vertex/DOF certificate; never assume a permutation."""
    geometry = checkpoint['geometry']
    if certificate.get('schema') != 'native-pilot-Q1-vertex-dof-v4' or '4.10' not in certificate.get('mfem_version', ''):
        raise ValueError('MFEM4.10 scalar Q1 coordinate certificate required')
    for key in ('nx', 'ny', 'x0_m', 'y0_m', 'width_m', 'thickness_m'):
        if certificate.get(key) != geometry[key]:
            raise ValueError('coordinate certificate is bound to different native geometry')
    x, y = edges(certificate['x_edges_m'], 'native x'), edges(certificate['y_edges_m'], 'native y')
    values = checkpoint['native_temperature_dofs_K']
    if len(x) != geometry['nx'] + 1 or len(y) != geometry['ny'] + 1 or certificate['node_count'] != len(values):
        raise ValueError('coordinate certificate native dimensions differ')
    result = [[None] * len(x) for _ in y]
    used = set()
    for vertex in certificate['vertices']:
        i, j, dof = vertex['i'], vertex['j'], vertex['dof']
        if (type(i) is not int or type(j) is not int or type(dof) is not int or not 0 <= i < len(x)
                or not 0 <= j < len(y) or not 0 <= dof < len(values) or dof in used or result[j][i] is not None
                or vertex['x_m'] != x[i] or vertex['y_m'] != y[j]):
            raise ValueError('coordinate certificate is not a complete one-to-one physical map')
        used.add(dof)
        result[j][i] = values[dof]
    if len(used) != len(values):
        raise ValueError('coordinate certificate omits native DOFs')
    return result


def snapshot_case(pilot_root, case, destination, *, coordinate_paths=None):
    """Copy completed CURRENT prefix and source inputs before any state analysis."""
    case = case.resolve()
    if not case.is_relative_to(pilot_root.resolve()):
        raise ValueError('case must belong to the identified actual pilot output')
    maximum = int((case / 'CURRENT').read_text().strip())
    if maximum < 1 or maximum > 100000:
        raise ValueError('actual pilot lacks a bounded nonempty accepted prefix')
    destination.mkdir(parents=True, exist_ok=False)
    files = {}

    def copy(source, target):
        target.parent.mkdir(parents=True, exist_ok=True)
        before = sha(source)
        shutil.copyfile(source, target)
        if sha(source) != before or sha(target) != before:
            raise ValueError('owner file changed during snapshot copy')
        files[str(target.relative_to(destination))] = {'source': str(source), 'sha256': before, 'bytes': target.stat().st_size}

    copy(pilot_root / 'component-contract.json', destination / 'component-contract.json')
    contract = read(destination / 'component-contract.json')
    source_root = Path(contract['code'])
    for path in sorted(source_root.rglob('*')):
        if path.is_file() and path.suffix in ('.cpp', '.hpp', '.h', '.json', '.csv', '.py', '.txt'):
            copy(path, destination / 'source/pilot' / path.relative_to(source_root))
    thermal_root = source_root.parent / 'mfem_electrode_heat_feedback_v2'
    for path in sorted(thermal_root.rglob('*')):
        if path.is_file() and path.suffix in ('.hpp', '.h', '.json'):
            copy(path, destination / 'source/thermal-v2' / path.relative_to(thermal_root))
    binary = source_root.parents[1] / 'build-solid-material-pilot/solid-material-feedback-pilot'
    copy(binary, destination / 'executable/solid-material-feedback-pilot')
    for body, path in (coordinate_paths or {}).items():
        copy(path, destination / 'certificates' / (body + '.json'))
    mfem_root = Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi')
    mfem_hashes = {str(path): sha(path) for path in (mfem_root / 'lib/libmfem.a', mfem_root / 'include/mfem/config/_config.hpp')}
    for step in range(maximum + 1):
        name = f'step_{step:06d}'
        for filename in ('upper.bin', 'lower.bin', 'fields.bin', 'row.csv'):
            copy(case / 'checkpoints' / name / filename, destination / 'checkpoints' / name / filename)
    if any(sha(item['source']) != item['sha256'] for item in files.values()):
        raise ValueError('source or accepted checkpoint changed while capturing immutable prefix')
    manifest = {'schema': 'readonly-actual-pilot-prefix-snapshot-v4', 'case': str(case),
        'captured_CURRENT_step': maximum, 'files': files, 'owner_files_modified': False,
        'frozen_MFEM_sha256': mfem_hashes,
        'source_code_to_executable_build_seal_verified': False,
        'source_binding_scope': 'copied accepted-state bytes and current source/executable observations; owner final build/run seal pending'}
    (destination / 'snapshot-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def audit_checkpoint_snapshot(snapshot, *, coordinate_paths, material, baseline_sha256):
    manifest = read(snapshot / 'snapshot-manifest.json')
    for name, item in manifest['files'].items():
        if sha(snapshot / name) != item['sha256']:
            raise ValueError('copied source/checkpoint snapshot hash mismatch')
    seal_text = (snapshot / 'source/thermal-v2/binding_seals.hpp').read_text()
    match = re.search(r'material_headers_sha256\[\]="([0-9a-f]{64})"', seal_text)
    if not match:
        raise ValueError('copied thermal-v2 material binding seal missing')
    material_headers_hash = match.group(1)
    certificates = {}
    certificate_bindings = {}
    for body, path in coordinate_paths.items():
        certificate = read(path)
        binary = Path(certificate['executable_path'])
        certificates[body] = certificate
        certificate_bindings[body] = {'path': str(path), 'sha256': sha(path),
            'helper_executable': str(binary), 'helper_executable_sha256': sha(binary),
            'helper_source_sha256': sha(ROOT / 'tools/pilot_heat_coordinates_v4.cpp')}
    stages = []
    reference_signatures = {}
    body_geometry = {}
    for step in range(manifest['captured_CURRENT_step'] + 1):
        directory = snapshot / 'checkpoints' / f'step_{step:06d}'
        row = [float(value) for value in (directory / 'row.csv').read_text().strip().split(',')]
        if len(row) != 27 or row[0] != step:
            raise ValueError('actual pilot row schema/revision differs from inspected source')
        stage = {'revision': step, 'time_s': row[1]}
        actual_min, actual_max = math.inf, -math.inf
        for body in ('upper', 'lower'):
            checkpoint = read_heat_checkpoint_v2(directory / (body + '.bin'))
            if checkpoint['revision'] != step or checkpoint['time_s'] != row[1] or checkpoint['geometry']['body'] != body:
                raise ValueError('native upper/lower/field-row accepted clock/body mismatch')
            if checkpoint['geometry']['material_headers_sha256'] != material_headers_hash:
                raise ValueError('actual native material signature differs from copied owner header binding')
            if step == 0:
                reference_signatures[body] = checkpoint['signature_sha256']
                body_geometry[body] = checkpoint['geometry']
                if checkpoint['geometry']['initial_temperature_K'] != 300. or any(value != 300. for value in checkpoint['native_temperature_dofs_K']):
                    raise ValueError('actual initial300K state absent; prehistory may not be reconstructed')
            elif checkpoint['signature_sha256'] != reference_signatures[body]:
                raise ValueError('native source/material/mesh signature changed in accepted history')
            values = checkpoint['native_temperature_dofs_K']
            actual_min, actual_max = min(actual_min, min(values)), max(actual_max, max(values))
            stage[body] = {'native_dof_count': len(values), 'native_min_K': min(values), 'native_max_K': max(values)}
            if body in certificates:
                stage[body]['nodal_temperature_K'] = canonical_temperature(checkpoint, certificates[body])
        if actual_min != row[6] or actual_max != row[7]:
            raise ValueError('native checkpoint extrema differ from accepted field/thermal row')
        stage['Tmin_K'], stage['Tmax_K'] = actual_min, actual_max
        stage['max_phase'], stage['max_liquid_fraction'] = row[8], row[9]
        stages.append(stage)
    if any(b['time_s'] <= a['time_s'] for a, b in zip(stages, stages[1:])):
        raise ValueError('actual accepted times do not increase')
    result = {'schema': 'actual-native-material-feedback-prefix-audit-v4', 'snapshot': str(snapshot),
        'snapshot_manifest_sha256': sha(snapshot / 'snapshot-manifest.json'), 'accepted_states': len(stages),
        'accepted_interval_s': [stages[0]['time_s'], stages[-1]['time_s']],
        'native_body_geometry': body_geometry, 'native_body_signature_sha256': reference_signatures,
        'canonical_order_verified': set(certificates) == {'upper', 'lower'},
        'coordinate_certificate_bindings': certificate_bindings,
        'frozen_MFEM_sha256': manifest['frozen_MFEM_sha256'],
        'native_material_headers_sha256': material_headers_hash,
        'native_cell_energy_export_available': False, 'caloric_native_energy_mapping_verified': False,
        'owner_final_build_run_seal_verified': False, 'physical_mechanics_admission': False,
        'computed_proof_screen': 'pending actual conditional mechanics computation',
        'support_condition': 'unknown; free/short-plane/both-end-clamped are conditional scenarios only'}
    range_prefix = smooth_prefix = 0
    range_open = smooth_open = True
    for stage in stages:
        good = (298.2 <= stage['Tmin_K'] and stage['Tmax_K'] <= 773.15
                and stage['max_phase'] == 0 and stage['max_liquid_fraction'] == 0)
        range_open = range_open and good
        smooth_open = smooth_open and good and stage['Tmin_K'] > 293.15 and stage['Tmax_K'] < 373.15
        range_prefix += int(range_open)
        smooth_prefix += int(smooth_open)
    result.update(range_prefix_state_count=range_prefix, smooth_prefix_state_count=smooth_prefix,
        range_prefix_last_time_s=stages[range_prefix - 1]['time_s'] if range_prefix else None,
        smooth_prefix_last_time_s=stages[smooth_prefix - 1]['time_s'] if smooth_prefix else None,
        smooth_prefix_method='conservative global-node first open293.15..373.15K coefficient segment',
        stage_temperature_summary=[{k: stage[k] for k in ('revision', 'time_s', 'Tmin_K', 'Tmax_K', 'max_phase', 'max_liquid_fraction')} for stage in stages])
    geometries = list(body_geometry.values())
    if any(g['reference_density_kg_m3'] != 8000. or g['initial_temperature_K'] != 300. for g in geometries):
        raise ValueError('actual pilot density/reference changed; selected nominal scenario no longer matches')
    scenario = pilot_density_scenario(material, baseline_sha256, density_kg_m3=8000., reference_temperature_K=300.)
    result['conditional_density_scenario'] = scenario
    result['conditional_density_scenario_sha256'] = hashlib.sha256(json.dumps(scenario, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    result['prescribed_temperature_geometry_comparison_ready'] = smooth_prefix >= 2 and result['canonical_order_verified']
    if result['prescribed_temperature_geometry_comparison_ready']:
        result['canonical_supported_prefix'] = {
            'schema': 'actual-pilot-prescribed-temperature-prefix-v4', 'native_snapshot_manifest_sha256': result['snapshot_manifest_sha256'],
            'material_density_scenario_sha256': result['conditional_density_scenario_sha256'],
            'temperature_source': 'both independent actual native heat checkpoints; source-bound copied prefix',
            'geometry': {'upper': {name: certificates['upper'][name] for name in ('x_edges_m', 'y_edges_m')},
                         'lower': {name: certificates['lower'][name] for name in ('x_edges_m', 'y_edges_m')},
                         'invariant_length_m': geometries[0]['invariant_length_m'], 'reference_density_kg_m3': 8000., 'reference_temperature_K': 300.},
            'states': stages[:smooth_prefix], 'native_caloric_ledger_conservation_claimed': False,
            'thermal_mechanical_two_way_feedback': False, 'structural_support_measured': False,
        }
    return result


def export_contract():
    return {
        'schema': EXPORT_SCHEMA,
        'meaning': 'Required owner-authored export; this file is a contract, not measured or simulated data.',
        'required_top_level': ['schema', 'status=sealed', 'pilot_contract_sha256', 'source_files',
            'executable', 'material_source_sha256', 'geometry', 'states'],
        'source_files': 'object mapping absolute source/data paths to SHA256; include material/scenario and current inputs',
        'executable': {'path': 'absolute executable path', 'sha256': 'exact run executable SHA256'},
        'geometry': {'invariant_length_m': .025, 'reference_density_kg_m3': 8000.,
            'reference_temperature_K': 300.,
            'upper': {'x_edges_m': 'full x nodes', 'y_edges_m': 'full y nodes'},
            'lower': {'x_edges_m': 'full x nodes', 'y_edges_m': 'full y nodes'}},
        'states': [{'revision': 'consecutive nonnegative accepted revision', 'time_s': 'actual accepted time',
            'upper': {'nodal_temperature_K': 'canonical (ny+1,nx+1) absolute K',
                      'cell_base_energy_J': 'canonical (ny,nx), native lumped h, full axial extent once',
                      'stored_energy_J': 'native body stored reference energy'},
            'lower': 'same fields from independent native lower owner, never inferred by this audit'}],
        'optional_per_body': ['nodal_phase (same shape; 0 means solid)', 'nodal_liquid_fraction'],
        'guards': ['hash verification before admission', 'complete both-foil accepted history',
                   'no temperature clipping or lower-body reconstruction', 'strict coefficient smooth-segment prefix',
                   'density convention reconciliation remains separate', 'computed proof screen still required'],
    }


def prefix_audit(export, material):
    """Find a conservative complete prefix; it is not mechanical proof admission."""
    geometry = export['geometry']
    height = number(geometry['invariant_length_m'], 'height')
    density = number(geometry['reference_density_kg_m3'], 'reference density')
    reference_temperature = number(geometry['reference_temperature_K'], 'reference temperature')
    if height != .025 or density <= 0:
        raise ValueError('explicit 25 mm invariant extent and positive source density required')
    axes = {}
    for body in ('upper', 'lower'):
        x = edges(geometry[body]['x_edges_m'], body + ' x')
        y = edges(geometry[body]['y_edges_m'], body + ' y')
        if not math.isclose(x[-1] - x[0], .005, rel_tol=2e-14, abs_tol=0):
            raise ValueError('pilot foil width differs from 5 mm target')
        if not math.isclose(y[-1] - y[0], .0002, rel_tol=2e-14, abs_tol=0):
            raise ValueError('pilot foil thickness differs from 0.2 mm target')
        axes[body] = (x, y)
    if not math.isclose(axes['upper'][1][0] - axes['lower'][1][-1], .005, rel_tol=2e-14, abs_tol=0):
        raise ValueError('pilot reference gap differs from 5 mm target')
    states = export['states']
    if not isinstance(states, list) or len(states) < 2:
        raise ValueError('at least two complete accepted states are required; endpoints alone are not a history')
    knots = sorted(set(material['thermal_expansion']['temperature_K'] + material['young_modulus']['temperature_K']))
    lower_bound = max(298.2, material['valid_temperature_K'][0])
    upper_bound = min(773.15, material['valid_temperature_K'][1])
    previous_revision, previous_time = None, None
    valid_count = 0
    smooth_count = 0
    first_range_failure = first_smooth_failure = None
    initial_segments = {}
    all_min, all_max = math.inf, -math.inf
    maximum_energy_residual = 0.
    energy_complete = True
    stage_rows = []
    for index, state in enumerate(states):
        revision = state['revision']
        time = number(state['time_s'], 'accepted time')
        if type(revision) is not int or revision < 0 or time < 0:
            raise ValueError('invalid accepted revision/time')
        if previous_revision is not None and (revision != previous_revision + 1 or time <= previous_time):
            raise ValueError('skipped/nonincreasing accepted revision or time; no silent interpolation')
        previous_revision, previous_time = revision, time
        stage_min, stage_max = math.inf, -math.inf
        range_ok = smooth_ok = True
        for body in ('upper', 'lower'):
            x, y = axes[body]
            body_state = state[body]  # Missing lower bulk is an error, never a symmetry assumption.
            temperature = matrix(body_state['nodal_temperature_K'], len(y), len(x), body + ' temperature')
            flat = [value for row in temperature for value in row]
            stage_min, stage_max = min(stage_min, min(flat)), max(stage_max, max(flat))
            range_ok = range_ok and all(lower_bound <= value <= upper_bound for value in flat)
            if 'nodal_phase' in body_state:
                phase = matrix(body_state['nodal_phase'], len(y), len(x), body + ' phase')
                range_ok = range_ok and all(value == 0 for row in phase for value in row)
            if 'nodal_liquid_fraction' in body_state:
                fraction = matrix(body_state['nodal_liquid_fraction'], len(y), len(x), body + ' liquid fraction')
                range_ok = range_ok and all(value == 0 for row in fraction for value in row)
            cell_segments = []
            for j in range(len(y) - 1):
                for i in range(len(x) - 1):
                    corners = [temperature[j][i], temperature[j][i + 1], temperature[j + 1][i + 1], temperature[j + 1][i]]
                    segments = [next((k for k in range(len(knots) - 1) if knots[k] < value < knots[k + 1]), None) for value in corners]
                    cell_segments.append(segments[0] if segments[0] is not None and len(set(segments)) == 1 else None)
            if index == 0:
                initial_segments[body] = cell_segments
            smooth_ok = smooth_ok and None not in cell_segments and cell_segments == initial_segments[body]
            if 'cell_base_energy_J' not in body_state or 'stored_energy_J' not in body_state:
                energy_complete = False
            else:
                energies = matrix(body_state['cell_base_energy_J'], len(y) - 1, len(x) - 1, body + ' native cell energy')
                flattened = [value for row in energies for value in row]
                stored = number(body_state['stored_energy_J'], body + ' stored energy')
                error = abs(math.fsum(flattened) - stored)
                bound = 5e-13 * math.fsum(abs(value) for value in flattened)
                if error > bound:
                    raise ValueError('native cell energy sum differs from body stored-energy ledger')
                maximum_energy_residual = max(maximum_energy_residual, error)
        all_min, all_max = min(all_min, stage_min), max(all_max, stage_max)
        if range_ok and first_range_failure is None:
            valid_count += 1
        elif first_range_failure is None:
            first_range_failure = {'index': index, 'revision': revision, 'time_s': time, 'Tmin_K': stage_min, 'Tmax_K': stage_max}
        if range_ok and smooth_ok and first_smooth_failure is None:
            smooth_count += 1
        elif first_smooth_failure is None:
            first_smooth_failure = {'index': index, 'revision': revision, 'time_s': time,
                'reason': 'temperature/phase range or coefficient-knot landing/crossing/cell span'}
        stage_rows.append({'index': index, 'revision': revision, 'time_s': time, 'Tmin_K': stage_min, 'Tmax_K': stage_max})
    # Same immutable expansion law converts the baseline table density to T_ref.
    expansion = material['thermal_expansion']
    temperatures, strains = expansion['temperature_K'], expansion['relative_length_change']
    if not temperatures[0] <= reference_temperature <= temperatures[-1]:
        raise ValueError('source reference temperature lies outside baseline mechanics record')
    k = min(next((k for k in range(len(temperatures) - 1) if reference_temperature <= temperatures[k + 1]), len(temperatures) - 2), len(temperatures) - 2)
    s = strains[k] + (strains[k + 1] - strains[k]) * (reference_temperature - temperatures[k]) / (temperatures[k + 1] - temperatures[k])
    baseline_density = material['density_reference']['value_kg_m3'] / (1 + s) ** 3
    density_compatible = math.isclose(density, baseline_density, rel_tol=1e-14, abs_tol=0)
    return {
        'accepted_states': len(states), 'accepted_interval_s': [states[0]['time_s'], states[-1]['time_s']],
        'all_history_temperature_range_K': [all_min, all_max],
        'material_range_K': [lower_bound, upper_bound], 'range_prefix_state_count': valid_count,
        'range_prefix_last_time_s': states[valid_count - 1]['time_s'] if valid_count else None,
        'smooth_prefix_state_count': smooth_count,
        'smooth_prefix_last_time_s': states[smooth_count - 1]['time_s'] if smooth_count else None,
        'first_range_failure': first_range_failure, 'first_smooth_failure': first_smooth_failure,
        'smooth_screen': 'every cell corner stays in its initial common open expansion/E segment; conservative nodal admission',
        'native_cell_energy_complete': energy_complete, 'maximum_native_cell_sum_residual_J': maximum_energy_residual,
        'native_reference_density_kg_m3': density, 'baseline_mechanics_reference_density_kg_m3': baseline_density,
        'mass_convention_compatible': density_compatible,
        'relative_source_mass_difference': density / baseline_density - 1,
        'data_ready_for_conditional_mechanics_trial': smooth_count >= 2 and energy_complete,
        'conservative_physical_mechanics_mapping_ready': False,
        'computed_proof_screen': 'not performed; temperature alone cannot establish elastic admission',
        'support_condition': 'not specified by pilot; conditional brackets only',
        'stage_temperature_summary': stage_rows,
    }


def verify_export(path, contract_path, material):
    original = sha(path)
    export = read(path)
    if export.get('schema') != EXPORT_SCHEMA or export.get('status') != 'sealed':
        raise ValueError('unrecognized or unsealed actual pilot history export')
    if export.get('pilot_contract_sha256') != sha(contract_path):
        raise ValueError('history bound to a different pilot contract')
    files = export.get('source_files', {})
    if not isinstance(files, dict) or not files:
        raise ValueError('source/data hash manifest required')
    for name, expected in files.items():
        if sha(name) != expected:
            raise ValueError('source/data file differs from sealed history: ' + name)
    binary = export['executable']
    if sha(binary['path']) != binary['sha256']:
        raise ValueError('pilot executable differs from sealed history')
    material_hash = export.get('material_source_sha256')
    if material_hash not in files.values():
        raise ValueError('material source seal must be present in verified source/data manifest')
    result = prefix_audit(export, material)
    if sha(path) != original or any(sha(name) != expected for name, expected in files.items()):
        raise ValueError('live history/source changed during read; immutable export required')
    result.update({'export_sha256': original, 'verified_source_files': files,
                   'verified_executable': binary, 'source_sealed': True})
    return result


def audit(pilot_root, history_export=None):
    contract_path = pilot_root / 'component-contract.json'
    if not contract_path.is_file():
        return {'status': 'actual_pilot_contract_not_available', 'pilot_output': str(pilot_root),
                'mechanics_data_ready': False, 'old_source_substitution': False}
    contract = read(contract_path)
    if contract.get('schema') != 'solid-only-native-material-feedback-pilot-v1':
        raise ValueError('unexpected actual pilot component contract')
    code = Path(contract['code'])
    source_manifest = code / 'copied-inputs.json'
    copied_inputs = []
    if source_manifest.is_file():
        for entry in read(source_manifest):
            target = (code / entry['path']).resolve()
            if not target.is_relative_to(code.resolve()):
                raise ValueError('copied-input path escapes actual pilot source')
            copied_inputs.append(dict(entry, actual_sha256=sha(target), verified=sha(target) == entry['sha256']))
    report = {'schema': 'ss304-material-feedback-mechanics-audit-v4', 'pilot_output': str(pilot_root),
        'pilot_source': str(code), 'contract_sha256': sha(contract_path), 'component_contract': contract,
        'status': 'actual_pilot_identified_history_not_yet_bound', 'source_sealed': False,
        'mechanics_data_ready': False, 'physical_mechanics_admission': False,
        'copied_input_verification': copied_inputs,
        'available_output_files': [{'path': str(p.relative_to(pilot_root)), 'bytes': p.stat().st_size}
            for p in sorted(pilot_root.rglob('*')) if p.is_file()],
        'density_admission_options': [
            'Separate Swiss-density source-T mechanics benchmark; no native mass/energy conservation claim.',
            'Separately versioned conditional native8000 mechanics density scenario, explicit integration-owner choice and validation.',
            'Physical two-way coupling remains blocked until a consistent mass/density convention is bound.'],
        'old_source_substitution': False, 'live_owner_files_modified': False, 'mechanics_solve_performed': False}
    if any(not entry['verified'] for entry in copied_inputs):
        raise ValueError('copied pilot input identity mismatch')
    if history_export is not None:
        material_path = ROOT / 'data/materials/ss304-v1.json'
        result = verify_export(history_export, contract_path, read(material_path))
        report.update(status='sealed_actual_history_audited', source_sealed=True,
            mechanics_data_ready=result['data_ready_for_conditional_mechanics_trial'],
            baseline_mechanics_record_sha256=sha(material_path), history=result)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pilot-output', type=Path, default=PILOT)
    parser.add_argument('--history-export', type=Path)
    parser.add_argument('--case', type=Path, help='Actual case under pilot-output; copies its completed CURRENT prefix before analysis')
    parser.add_argument('--coordinate-upper', type=Path, help='Read-only MFEM4.10 Q1 upper coordinate certificate')
    parser.add_argument('--coordinate-lower', type=Path, help='Read-only MFEM4.10 Q1 lower coordinate certificate')
    parser.add_argument('--output', type=Path, default=OUTPUT / 'audit.json')
    parser.add_argument('--export-contract', action='store_true')
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('audit destination exists; use a new report filename to preserve evidence')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.case is not None:
        if args.history_export or args.export_contract:
            raise ValueError('choose actual checkpoint audit or canonical export audit, not both')
        coordinate_paths = {body: path for body, path in (('upper', args.coordinate_upper), ('lower', args.coordinate_lower)) if path is not None}
        if len(coordinate_paths) == 1:
            raise ValueError('both foil coordinate certificates are required together')
        case = args.case if args.case.is_absolute() else args.pilot_output / args.case
        snapshot = args.output.with_name(args.output.stem + '-snapshot')
        snapshot_case(args.pilot_output, case, snapshot, coordinate_paths=coordinate_paths)
        material_path = ROOT / 'data/materials/ss304-v1.json'
        result = audit_checkpoint_snapshot(snapshot,
            coordinate_paths={body: snapshot / 'certificates' / (body + '.json') for body in coordinate_paths},
            material=read(material_path), baseline_sha256=sha(material_path))
        scenario = result.pop('conditional_density_scenario')
        scenario_path = args.output.with_name(args.output.stem + '-pilot8000-material.json')
        if scenario_path.exists():
            raise ValueError('conditional material scenario destination exists')
        scenario_path.write_bytes(json.dumps(scenario, sort_keys=True, separators=(',', ':')).encode())
        if sha(scenario_path) != result['conditional_density_scenario_sha256']:
            raise ValueError('conditional material scenario hash mismatch')
        result['conditional_density_scenario_path'] = str(scenario_path)
        prefix = result.pop('canonical_supported_prefix', None)
        if prefix is not None:
            prefix_path = args.output.with_name(args.output.stem + '-supported-prefix.json')
            if prefix_path.exists():
                raise ValueError('supported prefix destination exists')
            prefix_path.write_text(json.dumps(prefix, separators=(',', ':'), allow_nan=False) + '\n')
            result['canonical_supported_prefix_path'] = str(prefix_path)
            result['canonical_supported_prefix_sha256'] = sha(prefix_path)
        result['status'] = 'actual_copied_checkpoint_prefix_audited'
    else:
        result = export_contract() if args.export_contract else audit(args.pilot_output, args.history_export)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'status': result.get('status', 'export_contract_only'), 'output': str(args.output),
                      'mechanics_data_ready': result.get('mechanics_data_ready', result.get('prescribed_temperature_geometry_comparison_ready', False))}, indent=2))


if __name__ == '__main__':
    main()
