"""Independent unit, support and gate checks for the observation-only bridge."""

import json

import numpy as np
import pytest
from scipy.special import ndtr

from mhdlab.resolved_emission_observation import ResolvedTransition, project_resolved_emission


HC = 6.62607015e-34 * 299792458.0


def _inputs():
    return dict(
        time_s=np.array([0., 1., 2.]),
        upper_state_density_m3={"B:v0:J1": np.full((3, 1), 2e16)},
        path_weights_m=np.array([0.003]),
        transitions=[ResolvedTransition("B0J1-X0J0", "B:v0:J1", 150., 1e8)],
        bin_edges_nm=np.linspace(145., 155., 201), gate_start_s=0., gate_end_s=2.,
        fwhm_nm=0.1, time_registration="manufactured source and gate use the same clock",
        population_provenance="manufactured analytic upper-level density",
        geometry_provenance="manufactured 3-mm uniform ray", status="manufactured",
        optical_model="optically_thin",
    )


def _project(**changes):
    values = _inputs()
    values.update(changes)
    return project_resolved_emission(**values)


def test_uniform_slab_has_absolute_line_and_bin_units():
    result = _project()
    expected = 2e16 * 0.003 * 1e8 * HC / (4 * np.pi * 150e-9)
    assert result["lines"][0]["integrated_radiance_w_m2_sr"] == pytest.approx(expected)
    assert sum(result["bin_integrated_radiance_w_m2_sr"]) == pytest.approx(expected)
    widths = np.diff(result["bin_edges_nm"])
    np.testing.assert_allclose(np.asarray(result["bin_mean_radiance_w_m2_sr_nm"]) * widths,
                               result["bin_integrated_radiance_w_m2_sr"])
    assert result["lines"][0]["gate_mean_upper_column_m2"] == pytest.approx(6e13)
    assert result["population_closure_supplied"] is False
    assert result["density_temperature_fit_performed"] is False
    assert result["flash_simulation_performed"] is False
    json.dumps(result, allow_nan=False)


def test_cell_split_preserves_emission_and_collection_weight_is_applied_once():
    whole = _project()
    split = _project(
        upper_state_density_m3={"B:v0:J1": np.full((3, 3), 2e16)},
        path_weights_m=[0.0005, 0.001, 0.0015],
    )
    np.testing.assert_allclose(split["bin_integrated_radiance_w_m2_sr"],
                               whole["bin_integrated_radiance_w_m2_sr"])
    weighted = _project(path_weights_m=[0.0015])
    assert weighted["lines"][0]["integrated_radiance_w_m2_sr"] == pytest.approx(
        whole["lines"][0]["integrated_radiance_w_m2_sr"] / 2)


def test_disjoint_emitting_support_does_not_cross_multiply_cell_averages():
    result = _project(upper_state_density_m3={"B:v0:J1": [[2e16, 0], [0, 2e16], [2e16, 0]]},
                      path_weights_m=[0.003, 0.0])
    assert result["gate_mean_upper_state_density_m3"]["B:v0:J1"] == pytest.approx([1e16, 1e16])
    assert result["lines"][0]["gate_mean_upper_column_m2"] == pytest.approx(3e13)


def test_partial_gate_interpolates_both_endpoints_and_integrates_piecewise_history():
    # First cell is 2+3t; second cell changes slope at t=1. The gate is [0.25,1.5].
    result = _project(time_s=[0., 1., 2.],
                      upper_state_density_m3={"B:v0:J1": [[2., 0.], [5., 4.], [8., 2.]]},
                      path_weights_m=[2., 3.], gate_start_s=0.25, gate_end_s=1.5)
    # First mean=2+3*(0.25+1.5)/2; second trapezoids have areas 1.875+1.75.
    means = [4.625, (1.875 + 1.75) / 1.25]
    assert result["gate_mean_upper_state_density_m3"]["B:v0:J1"] == pytest.approx(means)
    assert result["lines"][0]["gate_mean_upper_column_m2"] == pytest.approx(2*means[0]+3*means[1])


def test_signed_before_breakdown_clock_has_analytic_partial_gate_average():
    result = _project(time_s=[-2., -1., 0.],
                      upper_state_density_m3={"B:v0:J1": [[2., 0.], [5., 4.], [8., 2.]]},
                      path_weights_m=[2., 3.], gate_start_s=-1.75, gate_end_s=-0.5,
                      time_registration="signed seconds relative to breakdown at zero")
    means = [4.625, (1.875 + 1.75) / 1.25]
    assert result["gate_mean_upper_state_density_m3"]["B:v0:J1"] == pytest.approx(means)
    assert result["lines"][0]["gate_mean_upper_column_m2"] == pytest.approx(2*means[0]+3*means[1])
    assert result["time_s"] == [-2., -1., 0.]
    assert result["gate"]["start_s"] == -1.75
    assert result["gate"]["end_s"] == -0.5


@pytest.mark.parametrize("shift", [-11., 7.5])
def test_emission_is_invariant_under_common_declared_clock_origin_shift(shift):
    values = dict(time_s=np.array([-2., -1., 0.]),
                  upper_state_density_m3={"B:v0:J1": [[2., 0.], [5., 4.], [8., 2.]]},
                  path_weights_m=[2., 3.], gate_start_s=-1.75, gate_end_s=-0.5)
    original = _project(**values)
    values.update(time_s=values["time_s"] + shift,
                  gate_start_s=values["gate_start_s"] + shift,
                  gate_end_s=values["gate_end_s"] + shift,
                  time_registration=f"same registered clock with declared origin shift {shift} s")
    shifted = _project(**values)
    assert shifted["time_s"] == pytest.approx(np.array(original["time_s"]) + shift)
    assert shifted["gate"]["start_s"] == pytest.approx(original["gate"]["start_s"] + shift)
    np.testing.assert_allclose(shifted["bin_integrated_radiance_w_m2_sr"],
                               original["bin_integrated_radiance_w_m2_sr"])
    np.testing.assert_allclose(shifted["bin_mean_radiance_w_m2_sr_nm"],
                               original["bin_mean_radiance_w_m2_sr_nm"])


def test_distinct_upper_states_and_shared_upper_branching_are_explicit():
    result = _project(
        upper_state_density_m3={"B:v0:J1": np.full((3, 1), 2e16),
                                "C:v1:J2": np.full((3, 1), 7e16)},
        transitions=[ResolvedTransition("first", "B:v0:J1", 149., 1e8),
                     ResolvedTransition("branch", "B:v0:J1", 150., 2e8),
                     ResolvedTransition("second", "C:v1:J2", 151., 1e8)],
        fwhm_nm=[0.1, 0.2, 0.05],
    )
    lines = result["lines"]
    assert [line["gate_mean_upper_column_m2"] for line in lines] == pytest.approx([6e13, 6e13, 2.1e14])
    assert lines[1]["integrated_radiance_w_m2_sr"] / lines[0]["integrated_radiance_w_m2_sr"] == pytest.approx(2*149/150)


@pytest.mark.parametrize("width", [0.01, 0.2, 1.0])
def test_line_area_is_independent_of_gaussian_width_and_nonuniform_bins(width):
    result = _project(bin_edges_nm=[140., 149.9, 150., 150.2, 160.], fwhm_nm=width)
    assert sum(result["bin_integrated_radiance_w_m2_sr"]) == pytest.approx(
        result["lines"][0]["integrated_radiance_w_m2_sr"])


def test_clipped_window_retains_lost_photons_without_renormalization():
    result = _project(bin_edges_nm=[150., 150.1], fwhm_nm=0.2)
    sigma = 0.2 / np.sqrt(8*np.log(2))
    fraction = ndtr(0.1/sigma) - 0.5
    assert sum(result["bin_integrated_radiance_w_m2_sr"]) == pytest.approx(
        fraction * result["lines"][0]["integrated_radiance_w_m2_sr"])
    assert sum(result["bin_integrated_radiance_w_m2_sr"]) < result["lines"][0]["integrated_radiance_w_m2_sr"] / 2


@pytest.mark.parametrize("changes", [
    {"time_s": [0.]}, {"time_s": [0., 0., 1.]}, {"time_s": [-1., 0., np.inf]},
    {"time_s": [0., np.nan, 2.]}, {"time_s": [0., 2., 1.]},
    {"gate_start_s": 1., "gate_end_s": 1.}, {"gate_start_s": -0.1},
    {"gate_end_s": 2.1}, {"gate_end_s": np.inf}, {"gate_start_s": "0 s"},
    {"path_weights_m": [-1.]}, {"path_weights_m": [1., 2.]},
    {"path_weights_m": [np.inf]}, {"path_weights_m": [True]},
    {"upper_state_density_m3": {"total_band": np.ones((3, 1))}},
    {"upper_state_density_m3": {"B:v0:J1": -np.ones((3, 1))}},
    {"upper_state_density_m3": {"B:v0:J1": np.ones((1, 3))}},
    {"upper_state_density_m3": {"B:v0:J1": np.full((3, 1), np.nan)}},
    {"bin_edges_nm": [151., 150.]}, {"bin_edges_nm": [150., 150.]},
    {"bin_edges_nm": [0., 150.]}, {"bin_edges_nm": [150.]},
    {"fwhm_nm": 0.}, {"fwhm_nm": [0.1, 0.2]}, {"fwhm_nm": "0.1 nm"},
    {"optical_model": "escape_factor"}, {"optical_model": ""},
    {"time_registration": ""}, {"population_provenance": "  "},
    {"geometry_provenance": None}, {"status": "actual_flash_result"},
    {"transitions": []},
])
def test_invalid_inputs_and_unregistered_or_unsupported_models_reject(changes):
    with pytest.raises(ValueError):
        _project(**changes)


@pytest.mark.parametrize("values", [
    ("", "B", 150., 1e8), ("a", "", 150., 1e8),
    ("a", "B", -150., 1e8), ("a", "B", "150 nm", 1e8),
    ("a", "B", 150., 0.), ("a", "B", 150., np.inf),
])
def test_transition_values_require_explicit_ids_and_physical_numeric_units(values):
    with pytest.raises(ValueError):
        ResolvedTransition(*values)


def test_duplicate_transition_ids_reject_but_shared_upper_state_does_not():
    transition = ResolvedTransition("same", "B:v0:J1", 150., 1e8)
    with pytest.raises(ValueError, match="unique"):
        _project(transitions=[transition, transition])


def test_empty_emitting_path_returns_zero_without_inventing_a_depth():
    result = _project(path_weights_m=[0.])
    assert result["lines"][0]["integrated_radiance_w_m2_sr"] == 0.
    assert np.count_nonzero(result["bin_integrated_radiance_w_m2_sr"]) == 0


def test_overflowing_column_is_rejected_instead_of_serializing_infinity():
    with pytest.raises(ValueError):
        _project(upper_state_density_m3={"B:v0:J1": np.full((3, 1), 1e300)},
                 path_weights_m=[1e300], bin_edges_nm=[149.9, 150.1])
