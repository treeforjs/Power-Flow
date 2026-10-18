"""Independent references for externally supplied coverage-response fits."""

from dataclasses import FrozenInstanceError
import json
import math

import numpy as np
import pytest
from scipy.optimize import least_squares, minimize_scalar

from mhdlab.coverage_regression import (
    BinnedObservation,
    CalibrationScale,
    CoveragePredictionGrid,
    fit_surface_coverage,
)


def grid(coverage=(0., 1., 2., 4.), predictions=None, **changes):
    coverage_input = coverage
    coverage = np.asarray(coverage)
    if predictions is None:
        predictions = np.column_stack((1. + coverage, 2. + .4 * coverage,
                                       .5 + 2. * coverage))
    nbin = np.shape(predictions)[1]
    arguments = dict(
        monolayer_areal_density_m2=2.75e18,
        monolayer_definition="Manufactured convention: one layer is 2.75e18 particles/m2",
        monolayer_provenance="Explicit manufactured layer convention; no kinetic Nref supplied",
        signal_units="manufactured_signal",
        bin_edges=np.arange(nbin + 1, dtype=float), bin_axis_units="manufactured_bin_coordinate",
        bin_support="same explicitly supplied fixed-bin support",
        status="manufactured", provenance="Manufactured external response grid",
    )
    arguments.update(changes)
    return CoveragePredictionGrid(coverage_input, predictions, **arguments)


def observation(values=(2., 3., 4.), covariance=None, background=None, **changes):
    nbin = np.size(values)
    arguments = dict(
        signal_units="manufactured_signal", bin_edges=np.arange(nbin + 1, dtype=float),
        bin_axis_units="manufactured_bin_coordinate",
        bin_support="same explicitly supplied fixed-bin support",
        status="manufactured", provenance="Manufactured external binned observations",
    )
    arguments.update(changes)
    return BinnedObservation(values, np.eye(nbin) if covariance is None else covariance,
                             np.zeros(nbin) if background is None else background, **arguments)


def calibration(mode="fixed", **arguments):
    if mode == "fixed" and not arguments:
        arguments["fixed_scale"] = 1.
    return CalibrationScale(mode, provenance="Manufactured scale convention", **arguments)


def interpolate_response(sample, coverage):
    """The external contract specifies linear interpolation between its rows."""
    return np.array([np.interp(coverage, sample.coverage_ml, column)
                     for column in sample.predictions.T])


def numerical_joint_reference(sample, data, *, prior_mean, prior_sigma):
    """Optimize both unknowns with a separate residual-based numerical solver.

    Each supplied coverage segment is considered separately, including its
    endpoints. This oracle does not use the profiled stationary-point formula.
    """
    chol = np.linalg.cholesky(data.covariance)
    target = data.values - data.known_background

    def residual(parameters):
        coverage, scale = parameters
        prediction = scale * interpolate_response(sample, coverage)
        whitened = np.linalg.solve(chol, target - prediction)
        return np.r_[whitened, (scale - prior_mean) / prior_sigma]

    candidates = []
    for left, right in zip(sample.coverage_ml[:-1], sample.coverage_ml[1:]):
        for scale_start in (.05, 1., 3.):
            answer = least_squares(residual, [float((left + right) / 2), scale_start],
                                   bounds=([left, 0.], [right, np.inf]),
                                   ftol=1e-13, xtol=1e-13, gtol=1e-13, max_nfev=1000)
            assert answer.success
            candidates.append((float(np.dot(answer.fun, answer.fun)), *answer.x))
    return min(candidates)


def gaussian_example():
    coverage = np.array([.2, .8, 1.7, 3.2, 5.])
    predictions = np.column_stack((1. + .4 * coverage, .5 + coverage**2 / 3.,
                                   2. + 1.5 / (1. + coverage)))
    sample = grid(coverage, predictions)
    background = np.array([.4, .8, .1])
    values = 1.1 * interpolate_response(sample, 2.1) + np.array([.3, -.4, .2]) + background
    covariance = np.array([[.36, .12, -.06], [.12, .25, .05], [-.06, .05, .49]])
    return sample, observation(values, covariance, background)


def test_known_coverage_off_grid_with_unscaled_background_and_explicit_layer_conversion():
    sample = grid()
    truth, scale = 2.3, 1.8
    background = np.array([4., 7., 11.])
    data = observation(scale * interpolate_response(sample, truth) + background,
                       [[.3, .1, 0.], [.1, .7, -.04], [0., -.04, .5]], background)
    answer = fit_surface_coverage(sample, data, calibration(fixed_scale=scale))
    best = answer["best_fit"]
    assert best["coverage_ml"] == pytest.approx(truth, rel=2e-12, abs=0.)
    assert best["calibration_scale"] == scale
    assert best["chi2"] < 1e-22
    assert best["areal_inventory_m2"] == pytest.approx(truth * 2.75e18, rel=2e-12, abs=0.)
    assert best["coverage_boundary"] == "interior"
    other_convention = grid(monolayer_areal_density_m2=8 * 2.75e18)
    other = fit_surface_coverage(other_convention, data, calibration(fixed_scale=scale))["best_fit"]
    assert other["coverage_ml"] == best["coverage_ml"]
    assert other["areal_inventory_m2"] == 8 * best["areal_inventory_m2"]
    json.dumps(answer, allow_nan=False)


def test_gaussian_scale_profile_matches_independent_one_dimensional_optimization():
    sample, data = gaussian_example()
    prior_mean, prior_sigma = .9, .2
    answer = fit_surface_coverage(sample, data, calibration("gaussian_prior", prior_mean=prior_mean,
                                                          prior_sigma=prior_sigma))
    target = data.values - data.known_background
    for index, response in enumerate(sample.predictions):
        precision_response = np.linalg.solve(data.covariance, response)
        expected_scale = max(0., (np.dot(precision_response, target) + prior_mean / prior_sigma**2) /
                             (np.dot(response, precision_response) + 1. / prior_sigma**2))

        def loss(scale):
            error = target - scale * response
            return np.dot(error, np.linalg.solve(data.covariance, error)) + ((scale-prior_mean)/prior_sigma)**2

        numerical = minimize_scalar(loss, bounds=(0., 5.), method="bounded",
                                    options={"xatol": 1e-13})
        assert numerical.success
        profile = answer["grid_profile"]
        assert profile["calibration_scale"][index] == pytest.approx(expected_scale, rel=2e-13, abs=0.)
        assert profile["calibration_scale"][index] == pytest.approx(numerical.x, rel=3e-8, abs=0.)
        assert profile["chi2"][index] == pytest.approx(loss(expected_scale), rel=3e-13, abs=0.)
        assert profile["chi2"][index] == pytest.approx(profile["data_chi2"][index] +
                                                    profile["prior_chi2"][index], rel=2e-13, abs=0.)


def test_gaussian_global_fit_matches_independent_two_parameter_residual_optimization():
    sample, data = gaussian_example()
    expected_loss, expected_coverage, expected_scale = numerical_joint_reference(
        sample, data, prior_mean=.9, prior_sigma=.2)
    best = fit_surface_coverage(sample, data, calibration("gaussian_prior", prior_mean=.9,
                                                         prior_sigma=.2))["best_fit"]
    assert best["coverage_ml"] == pytest.approx(expected_coverage, rel=3e-7, abs=0.)
    assert best["calibration_scale"] == pytest.approx(expected_scale, rel=3e-7, abs=0.)
    assert best["chi2"] == pytest.approx(expected_loss, rel=3e-12, abs=0.)


def test_selected_candidate_and_reported_signal_replay_the_same_objective():
    sample, data = gaussian_example()
    answer = fit_surface_coverage(sample, data, calibration("gaussian_prior", prior_mean=.9,
                                                          prior_sigma=.2))
    best = answer["best_fit"]
    assert best["chi2"] == min(candidate["chi2"] for candidate in answer["search_candidates"])
    signal = best["calibration_scale"] * interpolate_response(sample, best["coverage_ml"])
    np.testing.assert_allclose(best["prediction_signal"], signal, rtol=3e-14, atol=0.)
    np.testing.assert_allclose(best["prediction_with_background"], signal + data.known_background,
                               rtol=3e-14, atol=0.)
    error = data.values - data.known_background - signal
    objective = error @ np.linalg.solve(data.covariance, error) + ((best["calibration_scale"]-.9)/.2)**2
    assert best["chi2"] == pytest.approx(objective, rel=3e-13, abs=0.)


def test_full_correlated_covariance_changes_the_coverage_fit():
    sample = grid([0., .5, 1.], [[0., 0.], [.5, 1.], [1., 2.]])
    correlated = observation([1., 1.], [[1., .8], [.8, 1.]])
    diagonal = observation([1., 1.], np.eye(2))
    full = fit_surface_coverage(sample, correlated, calibration())["best_fit"]
    diag = fit_surface_coverage(sample, diagonal, calibration())["best_fit"]
    assert full["coverage_ml"] == pytest.approx(1. / 3., rel=2e-13, abs=0.)
    assert diag["coverage_ml"] == pytest.approx(.6, rel=2e-13, abs=0.)
    assert abs(full["coverage_ml"] - diag["coverage_ml"]) > .25


def test_free_absolute_scale_leaves_amplitude_only_coverage_unidentified():
    coverage = np.array([.25, 1., 3., 9.])
    shape = np.array([1., 2., 4.])
    sample = grid(coverage, coverage[:, None] * shape)
    answer = fit_surface_coverage(sample, observation(6. * shape), calibration("free"))
    assert np.max(answer["grid_profile"]["chi2"]) < 1e-24
    assert np.asarray(answer["grid_profile"]["calibration_scale"]) == pytest.approx(6. / coverage, rel=3e-13, abs=0.)
    assert answer["diagnostics"]["proportional_prediction_shapes"]
    assert answer["diagnostics"]["coverage_amplitude_degeneracy"]
    assert answer["diagnostics"]["flat_profile"]
    assert answer["diagnostics"]["coverage_unconstrained"]


def test_nonlinear_response_shape_identifies_coverage_with_free_absolute_scale():
    coverage = np.array([0., 1., 2., 3., 4.])
    sample = grid(coverage, np.column_stack((coverage, coverage**2, coverage**3)))
    truth, scale = 2.4, 1.7
    answer = fit_surface_coverage(sample, observation(scale * interpolate_response(sample, truth)),
                                  calibration("free"))
    best = answer["best_fit"]
    assert best["coverage_ml"] == pytest.approx(truth, rel=2e-11, abs=0.)
    assert best["calibration_scale"] == pytest.approx(scale, rel=2e-11, abs=0.)
    assert best["chi2"] < 1e-20
    assert not answer["diagnostics"]["coverage_amplitude_degeneracy"]
    assert not answer["diagnostics"]["coverage_unconstrained"]


@pytest.mark.parametrize("mode", ["free", "gaussian_prior"])
def test_weak_shape_information_survives_a_dominant_amplitude_bin(mode):
    # The exact solution is attained by the midpoint response at scale one.
    # Subtracting nearly equal O(1e40) Gram products can erase bin two.
    sample = grid([0., 1.], [[1e10, 1.], [2e10, 3.]])
    scale = calibration(mode, **({"prior_mean": 1., "prior_sigma": .2}
                                if mode == "gaussian_prior" else {}))
    best = fit_surface_coverage(sample, observation([1.5e10, 2.]), scale)["best_fit"]
    assert best["coverage_ml"] == pytest.approx(.5, rel=2e-7, abs=0.)
    assert best["calibration_scale"] == pytest.approx(1., rel=2e-7, abs=0.)
    assert best["chi2"] < 1e-12


def test_unrepresentable_interior_coverage_cannot_be_reported_as_an_endpoint_fit():
    # Midpoint prediction=1, while the midpoint coverage1e16+1 cannot be
    # represented. Reporting coverage1e16 would replay prediction0 instead.
    sample = grid([1e16, 1e16 + 2.], [[0.], [2.]])
    with pytest.raises(ValueError):
        fit_surface_coverage(sample, observation([1.]), calibration())


def test_whitening_cannot_silently_erase_a_positive_response():
    # A finite scale1e300 times these responses produces ordinary signals.
    # Whitening them by sqrt(1e300) underflows in float64 and must fail.
    sample = grid([0., 1.], [[1e-300], [2e-300]])
    data = observation([1.5], [[1e300]])
    with pytest.raises(ValueError):
        fit_surface_coverage(sample, data, calibration(fixed_scale=1e300))


def test_large_parameter_independent_chi_square_does_not_erase_profile_differences():
    # The first bin contributes1e16 at every coverage. The second bin gives
    # (1-2N)^2, so the minimum isN=.5 and both endpoint delta-chi2 values are1.
    sample = grid([0., 1.], [[0., 0.], [0., 2.]])
    answer = fit_surface_coverage(sample, observation([1e8, 1.]), calibration())
    assert answer["best_fit"]["coverage_ml"] == pytest.approx(.5, rel=2e-13, abs=0.)
    np.testing.assert_allclose(answer["grid_profile"]["delta_chi2"], [1., 1.], rtol=2e-13, atol=0.)
    assert answer["diagnostics"]["profile_span_chi2"] == pytest.approx(1., rel=2e-13, abs=0.)
    assert not answer["diagnostics"]["flat_profile"]


def test_exported_signal_and_residual_retain_the_profile_fit_before_parameter_rounding():
    # The feasible free-scale fit is exactly one signal unit. A float64
    # summary of1/49 can lose a unit in its final bit when multiplied by49;
    # exported fitted signals must retain the fit used for the objective.
    sample = grid([0., 1.], [[49.], [98.]])
    answer = fit_surface_coverage(sample, observation([1.], [[1e-32]]), calibration("free"))
    best = answer["best_fit"]
    assert best["prediction_signal"] == [1.]
    assert best["prediction_with_background"] == [1.]
    assert best["whitened_residual"] == [0.]
    assert best["chi2"] == 0.


@pytest.mark.parametrize("target,expected,boundary", [(0., 1., "lower"), (8., 3., "upper")])
def test_supported_coverage_bounds_are_not_extrapolated(target, expected, boundary):
    sample = grid([1., 2., 3.], [[1., 2.], [2., 4.], [3., 6.]])
    best = fit_surface_coverage(sample, observation([target, 2 * target]), calibration())["best_fit"]
    assert best["coverage_ml"] == expected
    assert best["coverage_boundary"] == boundary


def test_zero_response_grid_marks_free_scale_and_coverage_unidentified():
    answer = fit_surface_coverage(grid([0., 1., 3.], np.zeros((3, 2))), observation([1., 2.]),
                                  calibration("free"))
    assert answer["best_fit"]["calibration_scale"] == 0.
    assert answer["best_fit"]["chi2"] == pytest.approx(5., rel=2e-13, abs=0.)
    for key in ("zero_model_at_best", "scale_unidentified_at_best", "coverage_unconstrained"):
        assert answer["diagnostics"][key]


def test_zero_scale_constraint_is_explicit_for_negative_background_subtracted_signal():
    data = observation([.1, .2], background=[1., 1.])
    answer = fit_surface_coverage(grid([1., 2., 3.], [[1., 2.], [2., 3.], [4., 5.]]), data,
                                  calibration("free"))
    best = answer["best_fit"]
    assert best["calibration_scale"] == 0.
    assert best["zero_scale_boundary"]
    assert all(answer["grid_profile"]["zero_scale_boundary"])
    assert best["chi2"] == pytest.approx(.9**2 + .8**2, rel=2e-13, abs=0.)
    assert answer["diagnostics"]["coverage_unconstrained"]


def test_gaussian_calibration_constraint_keeps_the_nonzero_prior_penalty_at_scale_zero():
    sample = grid([1., 2., 3.], [[1., 2.], [2., 3.], [4., 5.]])
    answer = fit_surface_coverage(sample, observation([-2., -3.]),
                                  calibration("gaussian_prior", prior_mean=.2, prior_sigma=3.))
    best = answer["best_fit"]
    assert best["calibration_scale"] == 0.
    assert best["zero_scale_boundary"]
    assert best["data_chi2"] == pytest.approx(13., rel=2e-13, abs=0.)
    assert best["prior_chi2"] == pytest.approx((.2/3.)**2, rel=2e-13, abs=0.)
    assert best["chi2"] == pytest.approx(13. + (.2/3.)**2, rel=2e-13, abs=0.)


def test_signed_observation_and_background_are_supported_without_clipping():
    sample = grid([0., 1., 2.], [[0., 0.], [1., 2.], [2., 4.]])
    data = observation([-.2, .6], background=[-1.2, -1.4])
    best = fit_surface_coverage(sample, data, calibration())["best_fit"]
    assert best["coverage_ml"] == pytest.approx(1., rel=2e-13, abs=0.)
    assert best["chi2"] < 1e-24


@pytest.mark.parametrize("changes", [
    {"signal_units": "other_signal"}, {"bin_axis_units": "other_axis"},
    {"bin_support": "different gate"}, {"bin_edges": [0., 1., 2., 4.]},
])
def test_support_and_units_must_match_exactly(changes):
    with pytest.raises(ValueError):
        fit_surface_coverage(grid(), observation(**changes), calibration())


@pytest.mark.parametrize("bad", [
    [[1., .5, 0.], [.4, 1., 0.], [0., 0., 1.]],
    [[1., 2., 0.], [2., 1., 0.], [0., 0., 1.]],
    np.zeros((3, 3)), np.diag([1., 1., 0.]),
    np.diag([1., np.nan, 1.]), np.diag([1., np.inf, 1.]),
    np.eye(2),
])
def test_covariance_must_be_finite_symmetric_and_positive_definite(bad):
    with pytest.raises(ValueError):
        observation(covariance=bad)


@pytest.mark.parametrize("bad", [
    [True, False, True], ["1", "2", "3"], [1 + 0j, 2., 3.],
    np.ma.array([1., 2., 3.], mask=False), [1., np.nan, 3.], [1., np.inf, 3.],
])
def test_observation_rejects_coercions_masked_and_nonfinite_values(bad):
    with pytest.raises(ValueError):
        observation(bad)


@pytest.mark.parametrize("bad", [
    np.ones((4, 3), dtype=bool), np.full((4, 3), "1"), np.ones((4, 3), dtype=complex),
    np.ma.array(np.ones((4, 3)), mask=False), np.full((4, 3), np.nan),
    np.full((4, 3), np.inf), -np.ones((4, 3)), np.ones((3, 3)),
])
def test_predictions_reject_invalid_shape_type_range_and_masked_arrays(bad):
    with pytest.raises(ValueError):
        grid(predictions=bad)


@pytest.mark.parametrize("bad", [True, "2.75e18", 0., -1., np.nan, np.inf])
def test_monolayer_density_is_an_explicit_positive_finite_real_scalar(bad):
    with pytest.raises(ValueError):
        grid(monolayer_areal_density_m2=bad)


@pytest.mark.parametrize("bad", [[0., 1., 1., 3.], [0., 2., 1., 3.], [-1., 0., 1., 2.],
                                 [False, True, True, True], [0., 1., np.inf, 3.]])
def test_coverage_grid_requires_finite_ordered_nonnegative_real_values(bad):
    with pytest.raises(ValueError):
        grid(coverage=bad)


@pytest.mark.parametrize("mode,arguments", [
    ("unknown", {}), ("fixed", {}), ("fixed", {"fixed_scale": True}),
    ("fixed", {"fixed_scale": -1.}), ("fixed", {"fixed_scale": "1"}),
    ("fixed", {"fixed_scale": np.inf}),
    ("fixed", {"fixed_scale": 1., "prior_mean": 1., "prior_sigma": .2}),
    ("free", {"fixed_scale": 1.}), ("free", {"prior_mean": 1., "prior_sigma": .2}),
    ("gaussian_prior", {"prior_mean": 1., "prior_sigma": 0.}),
    ("gaussian_prior", {"prior_mean": -1., "prior_sigma": .2}),
    ("gaussian_prior", {"prior_mean": 1., "prior_sigma": True}),
])
def test_calibration_modes_have_strict_parameter_contracts(mode, arguments):
    with pytest.raises(ValueError):
        CalibrationScale(mode, provenance="manufactured", **arguments)


def test_array_inputs_are_detached_and_cannot_be_made_writeable():
    coverage = np.array([0., 1., 2., 4.])
    predictions = np.arange(12., dtype=float).reshape(4, 3)
    values, covariance, background = np.array([2., 3., 4.]), np.eye(3), np.array([.1, .2, .3])
    sample, data = grid(coverage, predictions), observation(values, covariance, background)
    coverage[:] = 99.
    predictions[:] = 99.
    values[:] = 99.
    covariance[:] = 99.
    background[:] = 99.
    np.testing.assert_array_equal(sample.coverage_ml, [0., 1., 2., 4.])
    np.testing.assert_array_equal(sample.predictions, np.arange(12.).reshape(4, 3))
    np.testing.assert_array_equal(data.values, [2., 3., 4.])
    np.testing.assert_array_equal(data.covariance, np.eye(3))
    np.testing.assert_array_equal(data.known_background, [.1, .2, .3])
    for array in (sample.coverage_ml, sample.predictions, sample.bin_edges,
                  data.values, data.covariance, data.known_background, data.bin_edges):
        with pytest.raises(ValueError):
            array.setflags(write=True)
    with pytest.raises((FrozenInstanceError, AttributeError)):
        sample.monolayer_areal_density_m2 = 1.
    with pytest.raises((FrozenInstanceError, AttributeError)):
        data.status = "sourced"
