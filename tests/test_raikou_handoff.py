"""Independent arithmetic and rejection checks for the external receiver contract."""

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from mhdlab.raikou_handoff import collected_path_weights, prepare_resolved_projection, strict_json, validate_raikou_handoff


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "configs" / "raikou_handoff_manufactured"


def bundle():
    manifest = json.loads((FIXTURE / "manifest.json").read_text())
    raw = {name: (FIXTURE / ref["relative_path"]).read_bytes() for name, ref in manifest["artifacts"].items()}
    return manifest, raw


def changed(mutate):
    """Rebind edited controls so rejection tests reach their scientific invariant."""
    manifest, raw = bundle()
    data = {name: json.loads(value) for name, value in raw.items() if name != "runtime_output"}
    mutate(data, manifest)
    links = {"receipt": {"source_sha256": "source", "runtime_output_sha256": "runtime_output"},
             "collection": {"mesh_sha256": "mesh"}, "rates": {"state_mapping_sha256": "states"},
             "populations": {"source_receipt_sha256": "receipt", "state_mapping_sha256": "states",
                             "mesh_sha256": "mesh", "rates_sha256": "rates", "electron_history_sha256": "electrons"},
             "opacity": {"populations_sha256": "populations", "collection_sha256": "collection"}}
    for name in ("states", "source", "runtime_output", "receipt", "mesh", "collection", "rates",
                 "electrons", "populations", "opacity"):
        if name != "runtime_output":
            for field, target in links.get(name, {}).items():
                data[name][field] = sha256(raw[target]).hexdigest()
            raw[name] = json.dumps(data[name], allow_nan=False).encode()
        manifest["artifacts"][name]["sha256"] = sha256(raw[name]).hexdigest()
    return manifest, raw


def test_manufactured_bundle_never_claims_physical_source_acceptance():
    report = validate_raikou_handoff(*bundle())
    assert report["receiver_receipt_structurally_valid"] is True
    assert report["physical_receiver_receipt_present"] is False
    assert report["independent_runtime_execution_verified"] is False
    assert report["thin_projection_admissible_under_owner_declarations"] is False
    assert "manufactured_evidence_only" in report["blockers"]
    assert "missing_rate_family:excitation" in report["blockers"]
    assert report["source_placement_checked_against_receiving_faces"] is True
    np.testing.assert_allclose(report["registered_photon_time_s"], [-5e-9, -3e-9, -1e-9])
    np.testing.assert_allclose(report["geometry"]["path_weights_m"], [0.001, 0.002])


def test_actual_gas_energy_does_not_have_to_equal_expected_source_translation():
    manifest, raw = bundle()
    receipt = json.loads(raw["receipt"])
    group = receipt["groups"][0]
    assert group["moments"]["translational_energy_j"] == 6.6271152e-20
    assert group["stationary_energy_ledger"]["gas_added_translational_j"] == 2e-19
    assert validate_raikou_handoff(manifest, raw)["receiver_receipt_structurally_valid"]


def test_hash_changed_bytes_reject_before_admission():
    manifest, raw = bundle()
    raw["source"] += b" "
    with pytest.raises(ValueError, match="source evidence hash"):
        validate_raikou_handoff(manifest, raw)


@pytest.mark.parametrize("field,value", [("before_parent_count", 12), ("after_parent_count", 6),
                                          ("gas_inserted_parent_count", 4)])
def test_parent_inventory_and_actual_gas_insertion_must_close(field, value):
    def mutate(data, _):
        data["receipt"]["inventory_ledgers"][0][field] = value
    with pytest.raises(ValueError, match="released_number"):
        validate_raikou_handoff(*changed(mutate))


def test_parent_inventory_debited_once_across_state_partitions():
    def split(data, _):
        source = deepcopy(data["source"]["groups"][0])
        source["group_id"] = "water1"
        source["state_id"] = "H2Op_X_000_000"
        data["source"]["groups"].append(source)
        received = deepcopy(data["receipt"]["groups"][0])
        received["group_id"] = "water1"
        received["identity"]["state_id"] = source["state_id"]
        data["receipt"]["groups"].append(received)
        pool = data["receipt"]["inventory_ledgers"][0]
        pool["before_parent_count"], pool["after_parent_count"], pool["gas_inserted_parent_count"] = 20, 10, 10
    assert validate_raikou_handoff(*changed(split))["source_group_count"] == 2

    def duplicate(data, manifest):
        split(data, manifest)
        data["receipt"]["inventory_ledgers"].append(deepcopy(data["receipt"]["inventory_ledgers"][0]))
    with pytest.raises(ValueError, match="exactly once"):
        validate_raikou_handoff(*changed(duplicate))


def test_rare_positive_component_cannot_be_dropped_by_global_tolerance():
    def mutate(data, _):
        rare = deepcopy(data["source"]["groups"][0])
        rare["group_id"] = "rare"
        rare["nodes"][0]["physical_particle_count"] = 1e-100
        rare["nodes"] = rare["nodes"][:1]
        data["source"]["groups"].append(rare)
    with pytest.raises(ValueError, match="every source component"):
        validate_raikou_handoff(*changed(mutate))


@pytest.mark.parametrize("field,value", [("temperature_first_moment_k", 2500),
                                          ("local_time_temperature_moment_k_s", 1.3e-6),
                                          ("translational_energy_j", 2e-19)])
def test_joint_moments_cannot_be_replaced_by_separate_means(field, value):
    def mutate(data, _):
        data["receipt"]["groups"][0]["moments"][field] = value
    with pytest.raises(ValueError, match=field):
        validate_raikou_handoff(*changed(mutate))


def test_covariance_readback_must_be_explicit_and_correct():
    def mutate(data, _):
        data["receipt"]["groups"][0]["covariance_local_time_temperature_k_s"] = 0
    with pytest.raises(ValueError, match="covariance"):
        validate_raikou_handoff(*changed(mutate))


@pytest.mark.parametrize("field,value", [("wall_thermal_after_j", -1e-18),
                                          ("gas_added_translational_j", 6.6271152e-20),
                                          ("wall_total_change_j", 0)])
def test_unavailable_thermal_energy_or_expected_instead_of_actual_debit_rejects(field, value):
    def mutate(data, _):
        data["receipt"]["groups"][0]["stationary_energy_ledger"][field] = value
    with pytest.raises(ValueError):
        validate_raikou_handoff(*changed(mutate))


def test_support_impulse_cannot_be_discarded():
    def mutate(data, _):
        data["receipt"]["groups"][0]["external_support_impulse_kg_m_per_s"] = [0, 0, 0]
    with pytest.raises(ValueError, match="momentum_ledger"):
        validate_raikou_handoff(*changed(mutate))


def test_atomic_n2_requires_distinct_2s_and_2p():
    def mutate(data, _):
        data["states"]["states"] = [s for s in data["states"]["states"] if s["state_id"] != "H_2s"]
        data["populations"]["density_m3"].pop("H_2s")
    with pytest.raises(ValueError, match="2s and 2p"):
        validate_raikou_handoff(*changed(mutate))


@pytest.mark.parametrize("state_id,field", [("H2Op_X_000_000", "vibrational_quantum_numbers"),
                                            ("OHm_X_v0_J0", "rotational_quantum_numbers")])
def test_molecular_cations_and_anions_require_rovibrational_resolution(state_id, field):
    def mutate(data, _):
        next(s for s in data["states"]["states"] if s["state_id"] == state_id).pop(field)
    with pytest.raises(ValueError, match=field):
        validate_raikou_handoff(*changed(mutate))


def test_populations_must_be_complete_resolved_state_arrays():
    def mutate(data, _):
        data["populations"]["density_m3"].pop("H2_B_v0_J1")
    with pytest.raises(ValueError, match="complete"):
        validate_raikou_handoff(*changed(mutate))


def test_mesh_order_cannot_be_assumed_from_array_shape():
    def mutate(data, _):
        data["populations"]["cell_ids"].reverse()
    with pytest.raises(ValueError, match="ordering"):
        validate_raikou_handoff(*changed(mutate))


@pytest.mark.parametrize("change", ["extrapolate_gate", "lost_time_precision", "missing_registration"])
def test_registered_gate_does_not_guess_or_extend_clock(change):
    def mutate(data, _):
        populations = data["populations"]
        if change == "extrapolate_gate":
            populations["photon_gate_s"][1] = 0
        elif change == "lost_time_precision":
            populations["clock_registration"]["photon_minus_source_time_s"] = 1e100
        else:
            populations.pop("clock_registration")
    with pytest.raises(ValueError):
        validate_raikou_handoff(*changed(mutate))


def test_resolved_velocity_distribution_can_replace_an_explicit_thermal_closure():
    def mutate(data, _):
        data["populations"]["velocity_closures"]["H_2p"] = {
            "model": "resolved_velocity_grid", "justification": "Manufactured two-stream distribution",
            "velocity_nodes_m_per_s": [[-1000, 0, 0], [1000, 0, 0]],
            "number_fraction": [[[0.25, 0.75], [0.75, 0.25]]] * 3}
    assert validate_raikou_handoff(*changed(mutate))["receiver_receipt_structurally_valid"]

    def bad(data, manifest):
        mutate(data, manifest)
        data["populations"]["velocity_closures"]["H_2p"]["number_fraction"][0][0] = [0.5, 0.75]
    with pytest.raises(ValueError, match="sum to one"):
        validate_raikou_handoff(*changed(bad))


def test_missing_emitting_velocity_or_opacity_support_rejects():
    def velocity(data, _):
        data["populations"]["velocity_closures"].pop("H_2p")
    with pytest.raises(ValueError, match="velocity"):
        validate_raikou_handoff(*changed(velocity))

    def opacity(data, _):
        data["opacity"]["transition_ids"].pop()
    with pytest.raises(ValueError, match="every retained transition"):
        validate_raikou_handoff(*changed(opacity))


def test_optical_depth_above_declared_limit_blocks_thin_projection():
    def mutate(data, _):
        data["opacity"]["maximum_optical_depth"] = 0.02
    report = validate_raikou_handoff(*changed(mutate))
    assert "resolved_radiative_transfer_required" in report["blockers"]


def test_status_relabeling_is_not_a_runtime_receiver_receipt():
    def mutate(data, manifest):
        manifest["status"] = "owner_supplied"
        for row in data.values():
            row["status"] = "owner_supplied"
    with pytest.raises(ValueError, match="runtime receipt"):
        validate_raikou_handoff(*changed(mutate))


def test_strict_json_rejects_duplicate_or_nonfinite_constants():
    with pytest.raises(ValueError, match="duplicate"):
        strict_json(b'{"a":1,"a":2}')
    with pytest.raises(ValueError, match="nonfinite"):
        strict_json(b'{"a":NaN}')
    with pytest.raises(ValueError, match="nonfinite"):
        strict_json(b'{"a":1e999}')


def test_finite_rays_disjoint_support_negative_direction_and_truncation():
    boxes = [[[0, 0, 0], [1, 1, 1]], [[2, 0, 0], [4, 1, 1]]]
    forward = collected_path_weights(boxes, [[-1, .5, .5]], [[5, .5, .5]], [1])
    backward = collected_path_weights(boxes, [[5, .5, .5]], [[-1, .5, .5]], [1])
    truncated = collected_path_weights(boxes, [[.5, .5, .5]], [[3, .5, .5]], [1])
    np.testing.assert_allclose(forward["path_weights_m"], [1, 2])
    np.testing.assert_allclose(backward["path_weights_m"], [1, 2])
    np.testing.assert_allclose(truncated["path_weights_m"], [.5, 1])


def test_parallel_ray_on_internal_face_has_one_half_open_owner():
    boxes = [[[0, 0, 0], [1, 1, 1]], [[0, 1, 0], [1, 2, 1]]]
    internal = collected_path_weights(boxes, [[-1, 1, .5]], [[2, 1, .5]], [1])
    outer = collected_path_weights(boxes, [[-1, 2, .5]], [[2, 2, .5]], [1])
    np.testing.assert_allclose(internal["path_weights_m"], [0, 1])
    np.testing.assert_allclose(outer["path_weights_m"], [0, 1])


def test_cell_subdivision_and_normalized_duplicate_rays_preserve_path():
    whole = [[[0, 0, 0], [2, 1, 1]]]
    split = [[[0, 0, 0], [1, 1, 1]], [[1, 0, 0], [2, 1, 1]]]
    a = collected_path_weights(whole, [[-1, .5, .5]], [[3, .5, .5]], [1])
    b = collected_path_weights(split, [[-1, .5, .5]], [[3, .5, .5]], [1])
    repeated = collected_path_weights(whole, [[-1, .5, .5]] * 2, [[3, .5, .5]] * 2, [.25, .75])
    assert sum(a["path_weights_m"]) == pytest.approx(sum(b["path_weights_m"]))
    np.testing.assert_allclose(repeated["path_weights_m"], a["path_weights_m"])
    assert a["detector_etendue_applied"] is False


@pytest.mark.parametrize("change", ["overlap", "bad_weights", "zero_ray", "nonfinite", "zero_cell"])
def test_geometry_rejects_double_paths_and_unsupported_inputs(change):
    boxes, starts, ends, weights = [[[0, 0, 0], [1, 1, 1]]], [[-1, .5, .5]], [[2, .5, .5]], [1]
    if change == "overlap":
        boxes *= 2
    elif change == "bad_weights":
        weights = [2]
    elif change == "zero_ray":
        ends = starts
    elif change == "nonfinite":
        starts[0][0] = float("nan")
    else:
        boxes[0][1][0] = 0
    with pytest.raises(ValueError):
        collected_path_weights(boxes, starts, ends, weights)


def test_cli_loads_current_source_and_makes_no_external_evolution_call():
    completed = subprocess.run([sys.executable, str(ROOT / "tools" / "validate_raikou_handoff.py"),
                                str(FIXTURE / "manifest.json")], capture_output=True, text=True, check=True)
    report = json.loads(completed.stdout)
    assert report["external_model_evolution_calls"] == 0
    assert report["validator_source_sha256"] == sha256((ROOT / "mhdlab" / "raikou_handoff.py").read_bytes()).hexdigest()
    assert report["physical_receiver_receipt_present"] is False


def test_projection_requires_explicit_manufactured_opt_in():
    instrument = json.loads((FIXTURE / "instrument.json").read_text())
    with pytest.raises(ValueError, match="allow_manufactured_control"):
        prepare_resolved_projection(*bundle(), instrument)


def test_instrument_status_and_population_reference_must_match_bundle():
    instrument = json.loads((FIXTURE / "instrument.json").read_text())
    instrument["status"] = "owner_supplied"
    with pytest.raises(ValueError, match="status must match"):
        prepare_resolved_projection(*bundle(), instrument, allow_manufactured_control=True)
    instrument["status"] = "manufactured"
    instrument["population_history_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="this population history"):
        prepare_resolved_projection(*bundle(), instrument, allow_manufactured_control=True)


def test_preparation_cannot_bypass_required_transfer_even_for_manufactured_control():
    def mutate(data, _):
        data["opacity"]["maximum_optical_depth"] = 2
    manifest, raw = changed(mutate)
    instrument = json.loads((FIXTURE / "instrument.json").read_text())
    instrument["population_history_sha256"] = manifest["artifacts"]["populations"]["sha256"]
    with pytest.raises(ValueError, match="resolved radiative transfer"):
        prepare_resolved_projection(manifest, raw, instrument, allow_manufactured_control=True)


def test_cli_projects_bound_populations_with_geometry_once_and_exact_line_units():
    completed = subprocess.run([sys.executable, str(ROOT / "tools" / "validate_raikou_handoff.py"),
                                str(FIXTURE / "manifest.json"), "--instrument-config", str(FIXTURE / "instrument.json"),
                                "--allow-manufactured-control"], capture_output=True, text=True, check=True)
    report = json.loads(completed.stdout)
    evidence = report["projection_evidence"]
    projection = evidence["projection"]
    assert evidence["velocity_projection_performed_by_adapter"] is False
    assert evidence["profile_quality_independently_verified"] is False
    assert projection["status"] == "manufactured"
    assert evidence["admission"]["physical_receiver_receipt_present"] is False
    # 1e10*1mm + 2e10*2mm = 5e7 upper emitters/m2 for each synthetic line.
    assert all(line["gate_mean_upper_column_m2"] == pytest.approx(5e7) for line in projection["lines"])
    expected = sum(6.62607015e-34 * 299792458.0 / (4 * np.pi * wavelength * 1e-9) * A * 5e7
                   for wavelength, A in [(130, 2), (150, 3)])
    assert sum(projection["bin_integrated_radiance_w_m2_sr"]) == pytest.approx(expected, rel=1e-13)
    assert np.dot(np.diff(projection["bin_edges_nm"]), projection["bin_mean_radiance_w_m2_sr_nm"]) == pytest.approx(expected, rel=1e-13)
    assert report["external_model_evolution_calls"] == 0


def test_positive_pool_cannot_be_hidden_under_absolute_count_tolerance():
    manifest, raw = bundle()
    manifest["conservation_tolerance"]["absolute"]["released_number"] = 1
    with pytest.raises(ValueError, match="rare positive pools"):
        validate_raikou_handoff(manifest, raw)


@pytest.mark.parametrize("missing", ["source_placement", "receiver_readback"])
def test_unverified_source_placement_remains_a_physical_admission_blocker(missing):
    def mutate(data, _):
        if missing == "source_placement":
            data["source"]["patches"][0].pop("placement")
        else:
            data["receipt"].pop("patch_placements")
    report = validate_raikou_handoff(*changed(mutate))
    assert report["source_placement_checked_against_receiving_faces"] is False
    assert "physical_source_placement_unverified" in report["blockers"]


@pytest.mark.parametrize("mismatch", ["readback", "wrong_receiving_cell", "area", "incomplete_face_support"])
def test_source_placement_must_match_receiving_faces_and_readback(mismatch):
    def mutate(data, _):
        source_patch = data["source"]["patches"][0]
        readback = data["receipt"]["patch_placements"][0]
        if mismatch == "readback":
            readback["placement"]["origin_m"][0] = 1e-6
        elif mismatch == "wrong_receiving_cell":
            readback["receiving_cell_ids"] = ["gas1"]
        elif mismatch == "area":
            source_patch["area_m2"] = 2e-12
        else:
            source_patch["placement"]["span_v_m"][0] = 0.005
            source_patch["area_m2"] = 5e-9
            readback["placement"] = deepcopy(source_patch["placement"])
            readback["receiving_cell_ids"] = ["gas0", "gas1"]
    with pytest.raises(ValueError, match="placement|complete source patch"):
        validate_raikou_handoff(*changed(mutate))


def test_a_permissive_owner_tau_limit_cannot_admit_moderately_thick_emission():
    def mutate(data, _):
        data["opacity"]["admission_limit"] = 0.9
        data["opacity"]["maximum_optical_depth"] = 0.5
    report = validate_raikou_handoff(*changed(mutate))
    assert "resolved_radiative_transfer_required" in report["blockers"]
    assert report["thin_emission_assessment"]["hard_optical_depth_limit"] == 0.01


def test_thin_admission_respects_a_tighter_radiance_error_budget():
    def mutate(data, _):
        data["opacity"]["maximum_optical_depth"] = 0.009
        data["opacity"]["relative_radiance_error_budget"] = 0.005
    report = validate_raikou_handoff(*changed(mutate))
    assert "resolved_radiative_transfer_required" in report["blockers"]
    assert report["thin_emission_assessment"]["fixed_emissivity_absorption_attenuation_bound"] == pytest.approx(1 - np.exp(-.009))


def test_thin_error_budget_cannot_exceed_one_percent():
    def mutate(data, _):
        data["opacity"]["relative_radiance_error_budget"] = 0.1
    with pytest.raises(ValueError, match="at most 1 percent"):
        validate_raikou_handoff(*changed(mutate))
