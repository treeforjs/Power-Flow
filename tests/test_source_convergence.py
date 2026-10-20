"""Independent manufactured-field and integrated-source convergence checks."""

import json
import math

import numpy as np
import pytest

from mhdlab.constants import KB, QE
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.source_convergence import (
    compare_observables, manufactured_planar_surface, observed_self_order,
    source_observables,
)
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_release import TemkinReservoir


STATE = "manufactured:gas"
MASS = 3e-26


def pool(**overrides):
    values = dict(name="manufactured-pool", state_id=STATE, attempt_frequency_s=1.,
                  activation_energy_ev=0., desorption_energy_ev=.3, internal_energy_ev=.2,
                  coverage_coefficient=0., reference_coverage_m2=1e20, max_coverage_m2=1e20,
                  kinetic_provenance="manufactured kinetic verification law",
                  energy_provenance="manufactured independent energy reference")
    values.update(overrides)
    return TemkinReservoir(**values)


def history(x_edges, z_edges, coverage, *, temperature=600., times=None, pools=None):
    times = np.array([0., math.log(2.)]) if times is None else np.asarray(times, dtype=float)
    x_edges, z_edges = np.asarray(x_edges, dtype=float), np.asarray(z_edges, dtype=float)
    shape = (len(times), len(z_edges) - 1, len(x_edges) - 1)
    temperature = np.broadcast_to(temperature, shape).copy()
    original = FlashSurfaceHistory(times, x_edges, temperature[:, 0, :], "synthetic", "manufactured",
                                   "manufactured analytic source-statistics fixture", "manufactured single face")
    attributes = original.attributes
    attributes["thermal_input_status"] = "manufactured"
    surface = PlanarSurfaceHistory(times, x_edges, z_edges, temperature, attributes)
    return evolve_temkin_planar_release(surface, np.asarray(coverage, dtype=float), [pool()] if pools is None else pools,
                                        parameter_provenance="manufactured source convergence verification")


def metrics(source, **kwargs):
    masses = {record.state_id: MASS for record in source.reservoirs}
    return source_observables(source, state_masses_kg=masses, **kwargs)["metrics"]


def patch_count_fixture(*, shift_x=0., shift_z=0.):
    x, z = np.array([0., 1., 3.]) + shift_x, np.array([0., 2., 5.]) + shift_z
    area = np.diff(z)[:, None] * np.diff(x)[None, :]
    emitted_counts = np.array([[1., 0.], [0., 3.]])
    # nu*dt=ln2 gives half depletion. Counts are prescribed independently of
    # area so an accidental equal-cell weighting cannot reproduce the result.
    return history(x, z, (2 * emitted_counts / area)[None])


@pytest.mark.parametrize("nx,nz,nt", [(1, 1, 1), (3, 2, 4), (4, 3, 3)])
def test_manufactured_surface_matches_independent_spatial_gaussian_quadrature(nx, nz, nt):
    width, length, duration = 4e-3, 3e-2, 8e-8
    surface = manufactured_planar_surface(width_bins=nx, current_bins=nz, snapshot_intervals=nt,
                                          duration_s=duration, width_m=width, current_length_m=length)
    nodes, weights = np.polynomial.legendre.leggauss(3)
    for t, time in enumerate(surface.time_s):
        for j, (z0, z1) in enumerate(zip(surface.z_edges_m[:-1], surface.z_edges_m[1:])):
            z_points = .5 * (z0 + z1) + .5 * (z1 - z0) * nodes
            for i, (x0, x1) in enumerate(zip(surface.x_edges_m[:-1], surface.x_edges_m[1:])):
                x_points = .5 * (x0 + x1) + .5 * (x1 - x0) * nodes
                point_temperature = 300 + (time / duration) ** 2 * (
                    250 + 1200 * (2 * x_points[None, :] / width) ** 2
                    * (.6 + .8 * z_points[:, None] / length))
                expected = np.einsum("i,j,ij->", weights, weights, point_temperature) / 4
                assert surface.temperature_k[t, j, i] == pytest.approx(expected, rel=2e-14, abs=0.)
    np.testing.assert_allclose(surface.time_s, np.linspace(0., duration, nt + 1), rtol=0., atol=0.)
    assert surface.attributes["source_kind"] == "synthetic"
    assert surface.attributes["thermal_input_status"] == "manufactured"
    assert "exact knot" in surface.attributes["manufactured_snapshot_interpolation"]


def test_manufactured_knots_follow_quadratic_time_not_temporal_bin_average_or_centroid():
    surface = manufactured_planar_surface(width_bins=1, current_bins=1, snapshot_intervals=2,
                                          duration_s=2., width_m=4., current_length_m=3.)
    # Averaging (2x/W)^2 over the full width gives 1/3, whereas the centroid
    # sample is zero. The resulting amplitude is 250+1200/3=650 K.
    np.testing.assert_allclose(surface.temperature_k[:, 0, 0], [300., 462.5, 950.], rtol=0., atol=2e-13)
    assert surface.temperature_k[1, 0, 0] != 300 + 650 / 3  # whole-interval time mean
    assert surface.temperature_k[-1, 0, 0] != 550.  # a spatial centroid sample


def test_actual_patch_areas_give_exact_number_and_energy_ledgers():
    source = patch_count_fixture()
    result = source_observables(source, state_masses_kg={STATE: MASS})
    values = result["metrics"]
    assert values["initial_number"] == pytest.approx([8.], rel=3e-14, abs=0.)
    assert values["released_number"] == pytest.approx([4.], rel=3e-14, abs=0.)
    assert values["remaining_number"] == pytest.approx([4.], rel=3e-14, abs=0.)
    energy_per_particle = {
        "translation_energy_j": 2 * KB * 600., "internal_energy_j": .2 * QE,
        "binding_energy_j": .3 * QE, "total_wall_energy_debit_j": 2 * KB * 600. + .5 * QE,
    }
    for name, expected in energy_per_particle.items():
        assert values[name][0] == pytest.approx(4 * expected, rel=4e-14, abs=0.)
        assert result["units"][name] == "J"
    assert abs(values["inventory_balance_residual_number"][0]) < 1e-14
    assert result["thermal_input_status"] == "manufactured"
    assert "no exact continuum" in result["qualifications"]["finite_reference"]
    json.dumps(result, allow_nan=False)


def test_resolved_rectangle_moments_include_within_patch_spread_and_covariance():
    values = metrics(patch_count_fixture())
    expected = {
        "resolved_mean_x_m": 13 / 8, "resolved_mean_z_m": 23 / 8,
        "resolved_second_raw_x_m2": 10 / 3, "resolved_second_raw_z_m2": 121 / 12,
        "resolved_variance_x_m2": 133 / 192, "resolved_variance_z_m2": 349 / 192,
        "resolved_covariance_xz_m2": 45 / 64,
    }
    for name, number in expected.items():
        assert values[name][0] == pytest.approx(number, rel=4e-14, abs=0.)


def test_width_uniform_control_preserves_current_distribution_and_erases_width_correlation():
    values = metrics(patch_count_fixture())
    expected = {"uniform_width_mean_x_m": 1.5, "uniform_width_variance_x_m2": .75,
                "uniform_width_second_raw_x_m2": 3., "uniform_width_covariance_xz_m2": 0.}
    for name, number in expected.items():
        assert values[name][0] == pytest.approx(number, rel=4e-14, abs=0.)
    for suffix in ("mean_z_m", "variance_z_m2", "second_raw_z_m2"):
        assert values["uniform_width_" + suffix] == values["resolved_" + suffix]
    assert values["uniform_width_mean_z_m"][0] != 2.5  # full-face redistribution
    assert values["resolved_covariance_xz_m2"][0] > 0


def test_single_rectangle_has_nonzero_physical_variance():
    values = metrics(history([2., 5.], [10., 14.], [[[1.]]]))
    assert values["resolved_mean_x_m"] == pytest.approx([3.5], rel=2e-14, abs=0.)
    assert values["resolved_mean_z_m"] == pytest.approx([12.], rel=2e-14, abs=0.)
    assert values["resolved_variance_x_m2"] == pytest.approx([9 / 12], rel=2e-14, abs=0.)
    assert values["resolved_variance_z_m2"] == pytest.approx([16 / 12], rel=2e-14, abs=0.)


def test_spatial_translation_preserves_centered_variances_and_covariance():
    original = metrics(patch_count_fixture())
    shifted = metrics(patch_count_fixture(shift_x=1e12, shift_z=2e12))
    for prefix in ("resolved", "uniform_width"):
        for suffix in ("variance_x_m2", "variance_z_m2", "covariance_xz_m2"):
            assert shifted[f"{prefix}_{suffix}"] == pytest.approx(original[f"{prefix}_{suffix}"], rel=3e-13, abs=1e-15)
    assert shifted["resolved_mean_x_m"][0] == 1e12 + 13 / 8
    assert shifted["resolved_mean_z_m"][0] == 2e12 + 23 / 8


def test_uniform_density_observables_are_invariant_to_nonuniform_partition():
    first = history([0., 1., 3.], [0., 2., 5.], np.ones((1, 2, 2)))
    second = history([0., .1, .3, 1.1, 3.], [0., .4, .6, 5.], np.ones((1, 3, 4)))
    first_values, second_values = metrics(first), metrics(second)
    for name in first_values:
        if name == "inventory_balance_residual_number":
            continue
        np.testing.assert_allclose(second_values[name], first_values[name], rtol=3e-13, atol=1e-15 if "covariance" in name else 0.)
    assert second_values["released_number"] == pytest.approx([7.5], rel=3e-14, abs=0.)


def test_thermal_velocity_and_energy_use_emitted_temperature_mixture():
    source = history([0., 1., 2.], [0., 1.], [[[2., 2.]]], temperature=[[[300., 1200.]]])
    values = metrics(source)
    base = math.sqrt(math.pi * KB * 300 / (2 * MASS))
    assert values["released_number"] == pytest.approx([2.], rel=3e-14, abs=0.)
    assert values["mean_velocity_normal_m_s"][0] == pytest.approx(1.5 * base, rel=3e-14, abs=0.)
    assert values["mean_velocity_normal_m_s"][0] != pytest.approx(math.sqrt(2.5) * base, rel=1e-3, abs=0.)
    assert values["mean_speed_squared_m2_s2"][0] == pytest.approx(4 * KB * 750 / MASS, rel=3e-14, abs=0.)
    assert values["translation_energy_j"][0] == pytest.approx(2 * KB * (300 + 1200), rel=3e-14, abs=0.)


def test_thermal_ratio_restores_a_representable_speed_after_intermediate_underflow():
    record = pool(reference_coverage_m2=1e300, max_coverage_m2=1e300)
    source = history([0., 1.], [0., 1.], [[[1e300]]], temperature=1e-310, pools=[record])
    values = source_observables(source, state_masses_kg={STATE: 1e-310})["metrics"]
    # kB*T underflows on its own, but division by this supplied mass restores
    # kB*T/m=kB, so both velocity moments remain positive and representable.
    assert KB * 1e-310 == 0.
    assert values["mean_velocity_normal_m_s"][0] == pytest.approx(math.sqrt(math.pi * KB / 2), rel=3e-13, abs=0.)
    assert values["mean_speed_squared_m2_s2"][0] == pytest.approx(4 * KB, rel=3e-13, abs=0.)


def test_unqualified_subnormal_source_probabilities_are_rejected_explicitly():
    source = history([0., 1., 2.], [0., 1.], [[[1e-300, 1e10]]])
    with pytest.raises(ValueError):
        source_observables(source, state_masses_kg={STATE: MASS})


def test_local_arrhenius_release_is_integrated_after_nonlinear_temperature_response():
    record = pool(activation_energy_ev=1., attempt_frequency_s=1e5)
    duration = .001
    source = history([0., 1., 2.], [0., 1.], [[[1., 1.]]], temperature=[[[300., 1200.]]],
                     times=[0., duration], pools=[record])
    expected = sum(-math.expm1(-1e5 * duration * math.exp(-QE / (KB * temperature)))
                   for temperature in (300., 1200.))
    averaged_first = 2 * (-math.expm1(-1e5 * duration * math.exp(-QE / (KB * 750))))
    result = metrics(source)["released_number"][0]
    assert result == pytest.approx(expected, rel=3e-13, abs=0.)
    assert result > 100 * averaged_first


def test_empty_reservoir_has_zero_extensive_budgets_and_null_normalized_moments():
    second_pool = pool(name="empty-pool", state_id="manufactured:empty")
    source = history([0., 1.], [0., 1.], [[[2.]], [[0.]]], pools=[pool(), second_pool])
    result = source_observables(source, state_masses_kg={STATE: MASS, second_pool.state_id: MASS})
    for name, values in result["metrics"].items():
        if any(token in name for token in ("mean_", "variance_", "second_raw_", "covariance_")):
            assert values[1] is None
        else:
            assert values[1] == 0.
    assert result["qualifications"]["empty_release_reservoirs"] == ["empty-pool"]
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("origin", [0., 1e9])
def test_absolute_emission_time_mean_and_centered_variance_match_global_exponential(origin):
    record = pool(attempt_frequency_s=1.3)
    source = history([0., 1.], [0., 1.], [[[2.]]], times=[origin, origin + .5, origin + 1.], pools=[record])
    result = source_observables(source, state_masses_kg={STATE: MASS}, include_time_moments=True,
                                time_quadrature_order=64)
    values = result["metrics"]
    rate, duration = record.attempt_frequency_s, 1.
    mean = 1 / rate - duration / math.expm1(rate * duration)
    variance = 1 / rate ** 2 - duration ** 2 * math.exp(rate * duration) / math.expm1(rate * duration) ** 2
    assert values["emission_mean_time_s"][0] == pytest.approx(origin + mean, rel=0., abs=2e-13 if origin == 0 else np.spacing(origin))
    assert values["emission_variance_time_s2"][0] == pytest.approx(variance, rel=3e-12, abs=0.)
    assert result["units"]["emission_mean_time_s"] == "s"
    assert result["units"]["emission_variance_time_s2"] == "s^2"
    assert "absolute" in result["qualifications"]["time_reference"]


def test_interval_splitting_preserves_emission_time_moments():
    record = pool(attempt_frequency_s=math.log(2.))
    one = history([0., 1.], [0., 1.], [[[1.]]], times=[0., 3.], pools=[record])
    split = history([0., 1.], [0., 1.], [[[1.]]], times=[0., 1., 3.], pools=[record])
    first, second = metrics(one, include_time_moments=True, time_quadrature_order=64), metrics(split, include_time_moments=True, time_quadrature_order=64)
    rate = math.log(2.)
    expected_mean, expected_variance = 1 / rate - 3 / 7, 1 / rate ** 2 - 72 / 49
    for values in (first, second):
        assert values["emission_mean_time_s"][0] == pytest.approx(expected_mean, rel=3e-12, abs=0.)
        assert values["emission_variance_time_s2"][0] == pytest.approx(expected_variance, rel=3e-12, abs=0.)


def test_comparison_is_labelled_finite_and_uses_per_component_absolute_and_relative_errors():
    result = compare_observables({"metrics": {"number": [8., 12.], "signed": [-4., 3.]}},
                                 {"metrics": {"number": [10., 10.], "signed": [-2., 0.]}})
    assert result["reference_kind"] == "finite_finest_resolution_reference"
    number = result["metrics"]["number"]
    assert number["absolute_difference"] == [2., 2.]
    assert number["relative_difference"] == [.2, .2]
    assert number["difference_l2"] == pytest.approx(math.sqrt(8), rel=2e-15, abs=0.)
    assert number["relative_l2"] == pytest.approx(.2, rel=2e-15, abs=0.)
    signed = result["metrics"]["signed"]
    assert signed["absolute_difference"] == [2., 3.]
    assert signed["relative_difference"] == [1., None]
    assert signed["warnings"]
    json.dumps(result, allow_nan=False)


def test_zero_and_null_reference_metrics_keep_undefined_ratios_explicit():
    result = compare_observables({"metrics": {"zero": [0., 2.], "missing": [None, 3.], "absent": [None]}},
                                 {"metrics": {"zero": [0., 0.], "missing": [8., None], "absent": [None]}})
    assert result["metrics"]["zero"]["absolute_difference"] == [0., 2.]
    assert result["metrics"]["zero"]["relative_difference"] == [None, None]
    assert result["metrics"]["zero"]["relative_l2"] is None
    assert result["metrics"]["missing"]["absolute_difference"] == [None, None]
    assert result["metrics"]["absent"]["difference_l2"] is None
    assert result["metrics"]["absent"]["reference_l2"] is None
    assert result["metrics"]["absent"]["relative_l2"] is None
    assert result["metrics"]["absent"]["warnings"]
    json.dumps(result, allow_nan=False)


def test_unrepresentable_relative_comparison_is_null_with_warning_not_infinity():
    result = compare_observables({"metrics": {"value": [1e300]}}, {"metrics": {"value": [1e-300]}})
    record = result["metrics"]["value"]
    assert record["absolute_difference"] == [1e300]
    assert record["relative_difference"] == [None]
    assert record["relative_l2"] is None
    assert record["warnings"]
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("scale", [1., 1e-200, 1e200])
@pytest.mark.parametrize("sign", [1., -1.])
def test_observed_order_uses_three_solution_differences_not_finest_reference_errors(scale, sign):
    exact, coefficient = 7 * scale, sign * 3 * scale
    levels = [exact + coefficient * spacing ** 2 for spacing in (1., .5, .25)]
    result = observed_self_order(*levels)
    assert result["observed_order"] == pytest.approx(2., rel=2e-13, abs=0.)
    assert not result["warnings"]
    # Relative to the finite finest level, the first two errors have ratio5,
    # which would spuriously report log2(5) instead of the three-level order2.
    assert result["observed_order"] != pytest.approx(math.log2(5), rel=1e-3, abs=0.)


def test_vector_observed_order_uses_difference_norms():
    exact, coefficient = np.array([1., 7., -2.]), np.array([3., -4., 5.])
    result = observed_self_order(*(exact + coefficient * spacing ** 2 for spacing in (1., .5, .25)))
    assert result["observed_order"] == pytest.approx(2., rel=2e-14, abs=0.)
    assert result["coarse_difference_l2"] / result["fine_difference_l2"] == pytest.approx(4., rel=2e-14, abs=0.)


def test_observed_order_uses_the_declared_nonbinary_refinement_ratio():
    ratio = 3.
    levels = [2 + 7 * spacing ** 2 for spacing in (1., 1 / ratio, 1 / ratio ** 2)]
    result = observed_self_order(*levels, refinement_ratio=ratio)
    assert result["observed_order"] == pytest.approx(2., rel=2e-14, abs=0.)
    assert result["refinement_ratio"] == ratio


@pytest.mark.parametrize("levels", [
    (0., 0., 0.), (2., 2., 2.), (1., 1., 0.), (1., 0., 0.),
    (1., 1.25, 1.1875), (0., 1., 3.), (1e16, 1e16, 1e16),
    ([2., 0.], [0., 0.], [1., 0.]),
])
def test_unresolved_oscillatory_or_nonasymptotic_sequences_do_not_claim_order(levels):
    result = observed_self_order(*levels)
    assert result["observed_order"] is None
    assert result["warnings"]
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("overrides", [
    {"width_bins": True}, {"current_bins": 0}, {"snapshot_intervals": 1.5},
    {"width_bins": "3"}, {"duration_s": 0.}, {"duration_s": True},
    {"width_m": -1.}, {"width_m": np.inf}, {"current_length_m": "1"},
])
def test_fixture_rejects_ambiguous_geometry_or_cadence(overrides):
    parameters = dict(width_bins=2, current_bins=2, snapshot_intervals=2,
                      duration_s=1., width_m=1., current_length_m=1.)
    parameters.update(overrides)
    with pytest.raises(ValueError):
        manufactured_planar_surface(**parameters)


@pytest.mark.parametrize("mass", [0., -1., True, "3e-26", 3e-26 + 0j, np.nan, np.inf])
def test_observables_require_explicit_positive_exact_state_masses(mass):
    with pytest.raises(ValueError):
        source_observables(patch_count_fixture(), state_masses_kg={STATE: mass})


def test_observables_reject_missing_mass_history_type_and_invalid_quadrature():
    source = patch_count_fixture()
    with pytest.raises(ValueError):
        source_observables(source, state_masses_kg={"manufactured": MASS})
    with pytest.raises(ValueError):
        source_observables(object(), state_masses_kg={STATE: MASS})
    with pytest.raises(ValueError):
        source_observables(source, state_masses_kg={STATE: MASS}, include_time_moments=1)
    for order in (1, True, "16", 2.5):
        with pytest.raises(ValueError):
            source_observables(source, state_masses_kg={STATE: MASS}, time_quadrature_order=order)


@pytest.mark.parametrize("ratio", [0., 1., -.5, True, "2", 2 + 0j, np.inf, np.nan])
def test_order_requires_a_finite_refinement_ratio_above_one(ratio):
    with pytest.raises(ValueError):
        observed_self_order(4., 1., .25, refinement_ratio=ratio)


@pytest.mark.parametrize("levels", [
    ([1., 2.], [1.], [1., 2.]), (True, 1., 0.), ("3", 1., 0.),
    (3 + 0j, 1., 0.), (np.nan, 1., 0.), (np.inf, 1., 0.),
    (np.ma.array([3.], mask=[True]), [1.], [0.]),
])
def test_order_rejects_invalid_or_inconsistent_level_values(levels):
    with pytest.raises(ValueError):
        observed_self_order(*levels)


@pytest.mark.parametrize("left,right", [
    ({"metrics": {"a": [1.]}}, {"metrics": {"b": [1.]}}),
    ({"metrics": {"a": [1., 2.]}}, {"metrics": {"a": [1.]}}),
    ({"metrics": {"a": [True]}}, {"metrics": {"a": [1.]}}),
    ({"metrics": {"a": [np.nan]}}, {"metrics": {"a": [1.]}}),
    ({"metrics": {"a": [1.]}, "units": {"a": "m"}}, {"metrics": {"a": [1.]}, "units": {"a": "cm"}}),
    ({"metrics": {"a": [1.]}, "reservoir_names": ["one"]}, {"metrics": {"a": [1.]}, "reservoir_names": ["two"]}),
])
def test_comparison_rejects_mismatched_or_invalid_metrics(left, right):
    with pytest.raises(ValueError):
        compare_observables(left, right)
