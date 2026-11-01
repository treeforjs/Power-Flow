"""Admission and projection checks for the actual-pilot comparison runner."""
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('pilot_geometry_v4', ROOT / 'tools/run_material_feedback_geometry_v4.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def test_reference_q1_average_preserves_independent_cells_and_z_order():
    nodal = [[300., 304., 312.], [308., 320., 340.]]
    averaged, maximum = runner.cell_temperature(nodal, 2)
    np.testing.assert_array_equal(averaged, [308., 319., 308., 319.])
    np.testing.assert_array_equal(maximum, [320., 340., 320., 340.])
    # Nonlinear h averaging is intentionally a different operation.
    assert averaged[0] ** 2 != np.mean(np.array(nodal)[:, :2] ** 2)


@pytest.mark.parametrize('values', [
    [[298.1, 300.], [300., 300.]], [[300., 373.15], [300., 300.]],
    [[300., np.nan], [300., 300.]], [300., 300.], [[300.], [300.]],
])
def test_projection_rejects_invalid_or_unsupported_nodal_data(values):
    with pytest.raises(ValueError):
        runner.cell_temperature(values, 2)


def admission_fixture(tmp_path):
    baseline_path = ROOT / 'data/materials/ss304-v1.json'
    baseline = json.loads(baseline_path.read_text())
    scenario = dict(baseline, record_id='ss304-pilot8000-density-scenario-v4',
        conditional_density_scenario={'baseline_record_sha256': runner.sha(baseline_path)})
    scenario_hash = hashlib.sha256(runner.json_bytes(scenario)).hexdigest()
    scenario_path = tmp_path / 'scenario.json'
    runner.write_json(scenario_path, scenario)
    manifest_path = tmp_path / 'snapshot-manifest.json'
    runner.write_json(manifest_path, {'files': {}})
    prefix = {'schema': 'actual-pilot-prescribed-temperature-prefix-v4',
        'native_snapshot_manifest_sha256': runner.sha(manifest_path),
        'material_density_scenario_sha256': scenario_hash,
        'geometry': {'reference_temperature_K': 300., 'reference_density_kg_m3': 8000.},
        'states': [{'revision': i, 'time_s': 1e-9 * (i + 1), 'Tmin_K': 300., 'Tmax_K': 301.,
                    'max_phase': 0, 'max_liquid_fraction': 0} for i in range(2)]}
    prefix_path = tmp_path / 'prefix.json'
    runner.write_json(prefix_path, prefix)
    audit = {'schema': 'actual-native-material-feedback-prefix-audit-v4',
        'prescribed_temperature_geometry_comparison_ready': True, 'canonical_order_verified': True,
        'snapshot': str(tmp_path), 'snapshot_manifest_sha256': runner.sha(manifest_path),
        'canonical_supported_prefix_path': str(prefix_path), 'canonical_supported_prefix_sha256': runner.sha(prefix_path),
        'conditional_density_scenario_sha256': scenario_hash, 'conditional_density_scenario_path': str(scenario_path)}
    audit_path = tmp_path / 'audit.json'
    runner.write_json(audit_path, audit)
    return audit_path, audit, prefix_path, prefix


def test_canonical_prefix_hash_binds_stage_bytes(tmp_path):
    path, audit, prefix_path, prefix = admission_fixture(tmp_path)
    runner.load_bound_inputs(path)
    prefix['states'][1]['Tmax_K'] = 302.
    runner.write_json(prefix_path, prefix)
    with pytest.raises(ValueError, match='prefix hash'):
        runner.load_bound_inputs(path)


def test_stale_or_knot_crossing_prefix_rejected_even_with_matching_hash(tmp_path):
    path, audit, prefix_path, prefix = admission_fixture(tmp_path)
    prefix['states'][1]['Tmax_K'] = 373.15
    runner.write_json(prefix_path, prefix)
    audit['canonical_supported_prefix_sha256'] = runner.sha(prefix_path)
    runner.write_json(path, audit)
    with pytest.raises(ValueError, match='coefficient segment'):
        runner.load_bound_inputs(path)
