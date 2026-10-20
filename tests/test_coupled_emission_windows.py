"""Manufactured polynomial controls; these fixtures are not source ODE solutions."""

from dataclasses import FrozenInstanceError
from fractions import Fraction
import hashlib
import json

import pytest

from mhdlab.coupled_emission_windows import (
    CoupledEmissionCell,
    WindowQuery,
    admit_saved_width_source,
    integrate_emission_window,
)


F = Fraction
KB = F("1.380649e-23")
QE = F("1.602176634e-19")
PROVENANCE = (
    "Manufactured exact polynomial/algebra-only fixture; not a source ODE solution."
)
ALIASES = (
    "report", "assembly", "input_config", "reuse_map", "uniform_power",
    "edge_power", "middle_power", "independent_completion", "independent_report",
    "scientific_review",
)


def _ratio(record):
    return F(int(record["numerator"]), int(record["denominator"]))


def _assert_record(record, expected):
    expected = F(expected)
    assert _ratio(record["exact"]) == expected
    assert record["value"] == float(expected)
    projection = record["projection"]
    assert projection["exact_zero"] is (expected == 0)
    assert projection["positive"] is (expected > 0)
    assert projection["negative"] is (expected < 0)
    assert projection["finite"] is True
    assert projection["underflowed"] is False
    assert projection["overflowed"] is False
    assert _ratio(projection["exact_rounding_error"]) == F(record["value"]) - expected


def _assert_stat(result, name, expected):
    stat = result["statistics"][name]
    assert stat["available"] is True
    assert stat["reason"] is None
    _assert_record(stat["record"], expected)


def _cell(**changes):
    # M(u)=2+3u+u^2, T(u)=10+2u, t=4+2u.
    kwargs = dict(
        start_s=F(4), stop_s=F(6),
        emitted_controls_m2=(F(2), F(3), F(13, 3), F(6)),
        temperature_controls_k=(F(10), F(32, 3), F(34, 3), F(12)),
        temperature_min_k=F(8), temperature_max_k=F(16),
        initial_coverage_m2=F(10), state_id="synthetic_parent_state",
        source_history_id="synthetic_polynomial", parameter_status="manufactured",
        parameter_provenance=PROVENANCE,
    )
    kwargs.update(changes)
    return CoupledEmissionCell(**kwargs)


def _integrate(cells=None, **changes):
    kwargs = dict(
        start_s=F(4), stop_s=F(6), desorption_energy_ev=F(1, 2),
        internal_energy_ev=F(1, 4), areal_heat_capacity_j_m2_k=F(2),
        power_w_m2=F(3), query_id="synthetic_window",
        parameter_provenance=PROVENANCE,
    )
    kwargs.update(changes)
    return integrate_emission_window((_cell(),) if cells is None else cells, **kwargs)


def test_full_cell_exact_number_time_temperature_and_joint_moments():
    result = _integrate()
    expected = {
        "released_m2": F(4),
        "local_time_first_moment_s_m2": F(13, 3),
        "local_time_second_moment_s2_m2": F(6),
        "temperature_first_moment_k_m2": F(133, 3),
        "temperature_second_moment_k2_m2": F(1478, 3),
        "local_time_temperature_moment_k_s_m2": F(148, 3),
    }
    for key, value in expected.items():
        _assert_record(result["areal"][key], value)
    for key, value in {
        "mean_local_time_s": F(13, 12), "variance_local_time_s2": F(47, 144),
        "mean_temperature_k": F(133, 12), "variance_temperature_k2": F(47, 144),
        "mean_local_time_temperature_k_s": F(37, 3),
        "covariance_local_time_temperature_k_s": F(47, 144),
    }.items():
        _assert_stat(result, key, value)
    assert _ratio(result["local_time_anchor_exact_s"]) == 4


def test_partial_cell_restriction_has_no_second_interval_length_factor():
    result = _integrate(start_s=F(9, 2), stop_s=F(11, 2))
    for key, value in {
        "released_m2": F(2),
        "local_time_first_moment_s_m2": F(25, 24),
        "local_time_second_moment_s2_m2": F(17, 24),
        "temperature_first_moment_k_m2": F(529, 24),
        "temperature_second_moment_k2_m2": F(2917, 12),
        "local_time_temperature_moment_k_s_m2": F(559, 48),
        "translational_energy_j_m2": F(529, 12) * KB,
    }.items():
        _assert_record(result["areal"][key], value)
    for key, value in {
        "mean_local_time_s": F(25, 48), "variance_local_time_s2": F(191, 2304),
        "mean_temperature_k": F(529, 48), "variance_temperature_k2": F(191, 2304),
        "mean_local_time_temperature_k_s": F(559, 96),
        "covariance_local_time_temperature_k_s": F(191, 2304),
    }.items():
        _assert_stat(result, key, value)
    assert _ratio(result["local_time_anchor_exact_s"]) == F(9, 2)


def test_quadratic_temperature_tests_product_degree_and_covariance():
    cell = _cell(
        start_s=F(7), stop_s=F(10), emitted_controls_m2=(0, F(2, 3), F(4, 3), 2),
        temperature_controls_k=(5, 5, F(16, 3), 6),
        temperature_min_k=4, temperature_max_k=8,
    )
    result = _integrate((cell,), start_s=7, stop_s=10)
    _assert_record(result["areal"]["temperature_second_moment_k2_m2"], F(856, 15))
    _assert_record(result["areal"]["local_time_temperature_moment_k_s_m2"], F(33, 2))
    _assert_record(result["areal"]["translational_energy_j_m2"], F(64, 3) * KB)
    for key, value in {
        "mean_local_time_s": F(3, 2), "variance_local_time_s2": F(3, 4),
        "mean_temperature_k": F(16, 3), "variance_temperature_k2": F(4, 45),
        "mean_local_time_temperature_k_s": F(33, 4),
        "covariance_local_time_temperature_k_s": F(1, 4),
    }.items():
        _assert_stat(result, key, value)


def test_multi_cell_joint_moments_use_one_common_local_anchor():
    left = _cell(start_s=0, stop_s=1, initial_coverage_m2=20)
    right = _cell(
        start_s=1, stop_s=2, initial_coverage_m2=20,
        emitted_controls_m2=(6, F(23, 3), F(29, 3), 12),
        temperature_controls_k=(12, F(38, 3), F(40, 3), 14),
    )
    whole = _integrate((left, right), start_s=0, stop_s=2)
    left_window = _integrate((left,), start_s=0, stop_s=1)
    right_window = _integrate((right,), start_s=1, stop_s=2)
    for key, value in {
        "released_m2": 10, "local_time_first_moment_s_m2": F(34, 3),
        "local_time_second_moment_s2_m2": 16,
        "temperature_first_moment_k_m2": F(368, 3),
        "temperature_second_moment_k2_m2": F(4552, 3),
        "local_time_temperature_moment_k_s_m2": F(436, 3),
    }.items():
        _assert_record(whole["areal"][key], value)
    for key, value in {
        "mean_local_time_s": F(17, 15), "variance_local_time_s2": F(71, 225),
        "mean_temperature_k": F(184, 15), "variance_temperature_k2": F(284, 225),
        "mean_local_time_temperature_k_s": F(218, 15),
        "covariance_local_time_temperature_k_s": F(142, 225),
    }.items():
        _assert_stat(whole, key, value)
    def q(result, name):
        return _ratio(result["areal"][name]["exact"])
    assert q(whole, "local_time_first_moment_s_m2") == (
        q(left_window, "local_time_first_moment_s_m2")
        + q(right_window, "local_time_first_moment_s_m2")
        + q(right_window, "released_m2")
    )
    assert q(whole, "local_time_second_moment_s2_m2") == (
        q(left_window, "local_time_second_moment_s2_m2")
        + q(right_window, "local_time_second_moment_s2_m2")
        + 2 * q(right_window, "local_time_first_moment_s_m2")
        + q(right_window, "released_m2")
    )
    assert q(whole, "local_time_temperature_moment_k_s_m2") == (
        q(left_window, "local_time_temperature_moment_k_s_m2")
        + q(right_window, "local_time_temperature_moment_k_s_m2")
        + q(right_window, "temperature_first_moment_k_m2")
    )


def test_energy_channels_and_nonzero_polynomial_defect_are_not_repaired():
    result = _integrate()
    translation = F(266, 3) * KB
    wall = translation + 3 * QE
    for key, value in {
        "translational_energy_j_m2": translation, "binding_energy_j_m2": 2 * QE,
        "internal_energy_j_m2": QE, "total_wall_energy_debit_j_m2": wall,
        "input_energy_j_m2": 6, "wall_thermal_change_j_m2": 4,
        "adsorbate_reference_change_j_m2": 2 * QE,
        "gas_energy_j_m2": translation + QE,
        "total_energy_residual_j_m2": wall - 2,
    }.items():
        _assert_record(result["areal"][key], value)
    assert result["areal"]["total_energy_residual_j_m2"]["value"] < 0
    changed = _integrate(power_w_m2=7, areal_heat_capacity_j_m2_k=3)
    for key in ("released_m2", "translational_energy_j_m2", "binding_energy_j_m2",
                "internal_energy_j_m2", "total_wall_energy_debit_j_m2"):
        assert changed["areal"][key] == result["areal"][key]
    _assert_record(changed["areal"]["input_energy_j_m2"], 14)
    _assert_record(changed["areal"]["wall_thermal_change_j_m2"], 6)
    _assert_record(changed["areal"]["total_energy_residual_j_m2"], wall - 8)


def test_huge_exact_origin_does_not_erase_local_time_or_joint_moments():
    origin = F(2**60)
    translated = _integrate((_cell(start_s=origin, stop_s=origin + 2),),
                            start_s=origin, stop_s=origin + 2)
    ordinary = _integrate()
    assert translated["areal"] == ordinary["areal"]
    assert translated["statistics"] == ordinary["statistics"]
    assert _ratio(translated["local_time_anchor_exact_s"]) == origin
    assert _ratio(translated["window_stop_exact_s"]) - origin == 2


def test_exact_zero_emission_keeps_input_thermal_defect_and_unavailable_means():
    result = _integrate((_cell(emitted_controls_m2=(2, 2, 2, 2)),))
    for name, record in result["areal"].items():
        expected = {"input_energy_j_m2": 6, "wall_thermal_change_j_m2": 4,
                    "total_energy_residual_j_m2": -2}.get(name, 0)
        _assert_record(record, expected)
    assert len(result["statistics"]) == 6
    for record in result["statistics"].values():
        assert record == {"available": False, "reason": "exact_zero_emitted_mass", "record": None}


def test_cell_is_immutable_and_outputs_are_detached_strict_json():
    cell = _cell()
    with pytest.raises((FrozenInstanceError, AttributeError)):
        cell.state_id = "changed"
    result = _integrate((cell,))
    json.dumps(result, allow_nan=False)
    result["areal"]["released_m2"]["exact"]["numerator"] = "999"
    fresh = _integrate((cell,))
    _assert_record(fresh["areal"]["released_m2"], 4)
    assert cell.state_id == "synthetic_parent_state"


@pytest.mark.parametrize("changes", [
    {"emitted_controls_m2": (2, 4, 3, 6)},
    {"emitted_controls_m2": (-1, 1, 3, 6)},
    {"emitted_controls_m2": (2, 3, 6, 11)},
    {"temperature_controls_k": (10, 17, 11, 12)},
    {"emitted_controls_m2": (2, 3, 6)},
    {"temperature_controls_k": (10, True, 11, 12)},
    {"emitted_controls_m2": (2, 3, float("nan"), 6)},
    {"initial_coverage_m2": F(-1, 10**1000)},
    {"parameter_status": "sourced"},
    {"parameter_provenance": " "},
])
def test_invalid_bernstein_support_and_labels_are_rejected(changes):
    with pytest.raises(ValueError):
        _cell(**changes)


@pytest.mark.parametrize("start,stop", [(3, 5), (5, 7), (5, 5), (6, 4)])
def test_window_support_is_strict_without_padding(start, stop):
    with pytest.raises(ValueError):
        _integrate(start_s=start, stop_s=stop)


@pytest.mark.parametrize("variant", ["gap", "overlap", "mass_jump", "temperature_jump",
                                     "mass_slope", "temperature_slope", "state", "history"])
def test_cell_chain_requires_continuity_and_joint_association(variant):
    left = _cell(start_s=0, stop_s=1, initial_coverage_m2=20)
    changes = dict(start_s=1, stop_s=2, initial_coverage_m2=20,
                   emitted_controls_m2=(6, F(23, 3), F(29, 3), 12),
                   temperature_controls_k=(12, F(38, 3), F(40, 3), 14))
    changes.update({
        "gap": {"start_s": F(5, 4)}, "overlap": {"start_s": F(3, 4)},
        "mass_jump": {"emitted_controls_m2": (7, 8, 10, 12)},
        "temperature_jump": {"temperature_controls_k": (13, 13, 14, 14)},
        "mass_slope": {"emitted_controls_m2": (6, 8, 10, 12)},
        "temperature_slope": {"temperature_controls_k": (12, 13, 14, 14)},
        "state": {"state_id": "different_state"},
        "history": {"source_history_id": "different_history"},
    }[variant])
    with pytest.raises(ValueError):
        _integrate((left, _cell(**changes)), start_s=0, stop_s=2)


def test_existing_frozen_cell_is_revalidated_after_control_forgery():
    cell = _cell()
    object.__setattr__(cell, "emitted_controls_m2", (F(2), F(4), F(3), F(6)))
    with pytest.raises(ValueError):
        _integrate((cell,))


@pytest.mark.parametrize("changes", [
    {"desorption_energy_ev": F(-1, 10**1000)},
    {"internal_energy_ev": True},
    {"areal_heat_capacity_j_m2_k": 0},
    {"power_w_m2": float("inf")},
])
def test_cost_power_and_heat_capacity_domain_guards(changes):
    with pytest.raises(ValueError):
        _integrate(**changes)


@pytest.mark.parametrize("case", ["mass_underflow", "mass_overflow", "time_variance_underflow"])
def test_unsupported_nonzero_projections_are_rejected(case):
    with pytest.raises(ValueError):
        if case == "mass_underflow":
            tiny = F(1, 10**400)
            cell = _cell(emitted_controls_m2=(0, tiny, 2*tiny, 3*tiny))
            _integrate((cell,))
        elif case == "mass_overflow":
            huge = F(10**400)
            cell = _cell(emitted_controls_m2=(0, huge, 2*huge, 3*huge),
                         initial_coverage_m2=4*huge)
            _integrate((cell,))
        else:
            tiny_duration = F(1, 10**200)
            cell = _cell(start_s=0, stop_s=tiny_duration,
                         emitted_controls_m2=(0, F(1, 3), F(2, 3), 1))
            _integrate((cell,), start_s=0, stop_s=tiny_duration)


@pytest.mark.parametrize("start,stop", [(True, F(1, 2)), (F(-1, 4), F(1, 2)),
                                        (F(1, 2), F(5, 4)), (F(3, 4), F(1, 4))])
def test_query_fractions_do_not_hide_invalid_support(start, stop):
    with pytest.raises(ValueError):
        WindowQuery("synthetic", start, stop, PROVENANCE)


def _admission_inputs():
    files = {alias: b"{}" for alias in ALIASES}
    hashes = {alias: hashlib.sha256(raw).hexdigest() for alias, raw in files.items()}
    return files, hashes


@pytest.mark.parametrize("variant", ["missing_alias", "extra_alias", "wrong_hash", "duplicate_json", "nonfinite_json", "wrong_schema"])
def test_source_admission_rejects_fresh_invalid_bytes_without_old_source_calls(variant):
    files, hashes = _admission_inputs()
    if variant == "missing_alias":
        files.pop("middle_power")
    elif variant == "extra_alias":
        files["unexpected"] = b"{}"
    elif variant == "wrong_hash":
        # Invalid JSON and a mismatching raw digest: the digest must reject first.
        files["report"] = b"not-json"
    else:
        files["report"] = {
            "duplicate_json": b'{"schema":"one","schema":"two"}',
            "nonfinite_json": b'{"schema":"one","bad":NaN}',
            "wrong_schema": b'{"schema":"manufactured-invalid-source-schema"}',
        }[variant]
        hashes["report"] = hashlib.sha256(files["report"]).hexdigest()
    with pytest.raises(ValueError) as exc:
        admit_saved_width_source(files, expected_sha256=hashes, parameter_provenance=PROVENANCE)
    if variant == "wrong_hash":
        assert "sha" in str(exc.value).lower() or "digest" in str(exc.value).lower()


@pytest.mark.parametrize("token,message", [
    (b"1e-400", "nonzero JSON float underflows its representation"),
    (b"-0", "negative-zero JSON integer is not canonical"),
])
def test_strict_json_tokens_reject_before_source_schema_admission(token, message):
    files, hashes = _admission_inputs()
    files["report"] = b'{"schema":"invalid-test-schema","probe":' + token + b'}'
    hashes["report"] = hashlib.sha256(files["report"]).hexdigest()
    with pytest.raises(ValueError, match=message):
        admit_saved_width_source(files, expected_sha256=hashes, parameter_provenance=PROVENANCE)
