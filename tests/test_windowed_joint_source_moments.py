"""Independent analytic mixtures for partial source-window joint moments."""
from dataclasses import replace
from fractions import Fraction
import json
import math

import numpy as np
import pytest

from mhdlab.constants import KB, QE
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.source_convergence import manufactured_planar_surface
from mhdlab.temkin_history import STEP_FIELDS, DIAGNOSTIC_FIELDS
from mhdlab.temkin_planar import TemkinPlanarReleaseHistory
from mhdlab.temkin_release import TemkinReservoir, release_step
from mhdlab.temkin_windowed_source import SourceWindowSelection, condition_temkin_planar_source
from mhdlab.thermal_release import SurfaceInventory
from mhdlab.windowed_joint_source_moments import (
    UNITS, WindowedJointMomentRangeError, windowed_joint_source_moments,
)


def pool(name="a", state="M:a", **changes):
    parameters = dict(name=name, state_id=state, attempt_frequency_s=3., activation_energy_ev=.08,
                      desorption_energy_ev=.3, internal_energy_ev=.1, coverage_coefficient=0.,
                      reference_coverage_m2=1e18, max_coverage_m2=1e22,
                      kinetic_provenance="Manufactured independent constant-hazard fixture",
                      energy_provenance="Manufactured supplied wall costs")
    parameters.update(changes)
    return TemkinReservoir(**parameters)


def source_record(*, pools=None, initial=None, temperature=None, x_edges=None, z_edges=None,
                  origin=10., duration=2., windows=None):
    pools = [pool()] if pools is None else pools
    x_edges = np.array([-.002, .004]) if x_edges is None else np.asarray(x_edges)
    z_edges = np.array([0., .003]) if z_edges is None else np.asarray(z_edges)
    shape = (len(z_edges)-1, len(x_edges)-1)
    initial = np.full((len(pools), *shape), 1e18) if initial is None else np.asarray(initial)
    temperature = np.full(shape, 700.) if temperature is None else np.asarray(temperature)
    base = manufactured_planar_surface(width_bins=1, current_bins=1, snapshot_intervals=1,
                                      duration_s=1., width_m=.002, current_length_m=.004)
    attributes = dict(base.attributes)
    for name in ("current_direction_representation", "current_patch_intervals_m", "patch_temperature_approximation"):
        attributes.pop(name, None)
    surface = PlanarSurfaceHistory(np.array([origin, origin+duration]), x_edges, z_edges,
                                   np.stack([temperature, temperature]), attributes)
    trial = release_step(SurfaceInventory(initial.reshape(len(pools), -1)), pools,
                         temperature.reshape(-1), duration)
    history = TemkinPlanarReleaseHistory(
        time_edges_s=surface.time_s, x_edges_m=x_edges, z_edges_m=z_edges,
        temperature_k=temperature[None], coverage_m2=np.stack([initial, trial.trial_inventory.coverage_m2.reshape(initial.shape)]),
        **{name: getattr(trial, name).reshape((1, *initial.shape)) for name in (*STEP_FIELDS, *DIAGNOSTIC_FIELDS)},
        reservoirs=tuple(pools), surface_attributes=surface.attributes,
        parameter_provenance="Manufactured conditional moment fixture, no measured physical inputs")
    windows = [(.1, .7)] if windows is None else windows
    return condition_temkin_planar_source(history, [SourceWindowSelection(0, *bounds) for bounds in windows],
                                         log_mass_tolerance=1e-7, status="manufactured",
                                         provenance="Manufactured independent window moments")


def inspect(source, **changes):
    parameters = dict(quadrature_order=64)
    parameters.update(changes)
    return windowed_joint_source_moments(source, **parameters)


def exponential_moments(rate, duration):
    tail = math.exp(-rate*duration)
    norm = -math.expm1(-rate*duration)
    mean = 1/rate-duration*tail/norm
    variance = 1/rate**2-duration**2*tail/norm**2
    return mean, variance


@pytest.fixture(scope="module")
def mixture():
    return source_record(pools=[pool(), pool("b", "M:b", attempt_frequency_s=1.1, activation_energy_ev=.03)],
                         initial=[[[1e18, 2e18], [3e18, 0.]], [[.5e18, 0.], [1e18, 4e18]]],
                         temperature=[[600., 900.], [1100., 700.]],
                         x_edges=[-.004, -.001, .005], z_edges=[0., .002, .008],
                         windows=[(.1, .6), (.8, 1.7)])


def test_single_component_closed_truncated_exponential_and_uniform_rectangles():
    source = source_record()
    report = inspect(source)
    component = source.components[0]
    window = component.source_window
    rate = window.reservoir.attempt_frequency_s*math.exp(-window.reservoir.activation_energy_ev*QE/(KB*window.temperature_k))
    mean, variance = exponential_moments(rate, window.stop_s-window.start_s)
    values = report["resolved"]
    assert values["mean_time_s"] == pytest.approx(component.interval_origin_s+window.start_s+mean, rel=1e-13)
    assert values["mean_time_offset_s"] == pytest.approx(mean, rel=3e-12)
    assert values["variance_time_s2"] == pytest.approx(variance, rel=3e-10)
    assert values["mean_x_m"] == pytest.approx(.001)
    assert values["variance_x_m2"] == pytest.approx(.006**2/12)
    assert values["mean_z_m"] == pytest.approx(.0015)
    assert values["variance_z_m2"] == pytest.approx(.003**2/12)
    assert values["mean_width_centered_square_m2"] == pytest.approx(.006**2/12)
    assert values["variance_width_centered_square_m4"] == pytest.approx(.006**4/180)
    assert values["variance_expected_translational_energy_j2"] == 0.
    assert values["variance_kinetic_energy_j2"] == pytest.approx(2*(KB*700)**2, rel=3e-15)
    assert values["second_raw_kinetic_energy_j2"] == pytest.approx(6*(KB*700)**2, rel=3e-15)
    assert report["resolved"] == report["matched_uniform_width"]


def test_analytic_mixture_counts_time_position_energy_covariances(mixture):
    report = inspect(mixture)
    rows = []
    for component in mixture.components:
        w = component.source_window
        rate = w.reservoir.attempt_frequency_s*math.exp(-w.reservoir.activation_energy_ev*QE/(KB*w.temperature_k))
        mean, variance = exponential_moments(rate, w.stop_s-w.start_s)
        count = w.initial_coverage_m2*w.patch_area_m2*math.exp(-rate*w.start_s)*-math.expm1(-rate*(w.stop_s-w.start_s))
        rows.append((count, component.interval_origin_s+w.start_s+mean, variance,
                     sum(component.x_bounds_m)/2, sum(component.z_bounds_m)/2, 2*KB*w.temperature_k))
    total = math.fsum(row[0] for row in rows)
    means = [math.fsum(row[0]*row[column] for row in rows)/total for column in (1, 3, 4, 5)]
    def covariance(a, b, ma, mb):
        return math.fsum(row[0]*(row[a]-ma)*(row[b]-mb) for row in rows)/total
    values = report["resolved"]
    assert values["released_number"] == pytest.approx(total, rel=3e-12)
    assert values["mean_time_s"] == pytest.approx(means[0], rel=3e-13)
    assert values["variance_time_s2"] == pytest.approx(covariance(1, 1, means[0], means[0])+math.fsum(row[0]*row[2] for row in rows)/total, rel=2e-10)
    assert values["cov_time_x_s_m"] == pytest.approx(covariance(1, 3, means[0], means[1]), rel=3e-10)
    assert values["cov_time_z_s_m"] == pytest.approx(covariance(1, 4, means[0], means[2]), rel=3e-10)
    assert values["cov_time_expected_translational_energy_s_j"] == pytest.approx(covariance(1, 5, means[0], means[3]), rel=3e-10)
    assert values["cov_expected_translational_energy_x_j_m"] == pytest.approx(covariance(5, 3, means[3], means[1]), rel=3e-12)
    assert values["cov_expected_translational_energy_z_j_m"] == pytest.approx(covariance(5, 4, means[3], means[2]), rel=3e-12)


def test_uniform_width_preserves_all_other_marginals_and_clears_width_covariances(mixture):
    report = inspect(mixture)
    left, right = mixture.x_edges_m[[0, -1]]
    uniform = report["matched_uniform_width"]
    for name in ("released_number", "mean_time_s", "mean_time_offset_s", "variance_time_s2", "mean_z_m", "variance_z_m2",
                 "mean_expected_translational_energy_j", "variance_expected_translational_energy_j2", "variance_kinetic_energy_j2",
                 "mean_conditional_kinetic_energy_variance_j2", "second_raw_kinetic_energy_j2",
                 "cov_time_z_s_m", "cov_time_expected_translational_energy_s_j", "cov_expected_translational_energy_z_j_m"):
        assert uniform[name] == report["resolved"][name]
    for name in ("cov_time_x_s_m", "cov_expected_translational_energy_x_j_m", "cov_time_width_centered_square_s_m2",
                 "cov_expected_translational_energy_width_centered_square_j_m2"):
        assert uniform[name] == 0.
    assert uniform["mean_x_m"] == pytest.approx((left+right)/2)
    assert uniform["variance_x_m2"] == pytest.approx((right-left)**2/12)
    assert report["resolved"]["cov_time_x_s_m"] != 0.
    assert report["resolved"]["cov_expected_translational_energy_x_j_m"] != 0.
    assert report["resolved"]["variance_x_m2"] != uniform["variance_x_m2"]


def test_groups_and_component_records_recombine_without_changing_source_labels(mixture):
    report = inspect(mixture)
    assert len(report["components"]) == len(mixture.components) == 12
    assert len(report["per_selection"]) == 2
    assert len(report["per_reservoir"]) == 2
    assert report["status"] == "manufactured"
    assert report["parameter_provenance"] == mixture.parameter_provenance
    assert report["source"] == mixture.manifest()
    total = report["resolved"]["released_number"]
    for groups in (report["per_selection"], report["per_reservoir"]):
        assert math.fsum(group["resolved"]["released_number"] for group in groups) == pytest.approx(total, rel=3e-15)
        for name in ("mean_x_m", "mean_z_m", "mean_expected_translational_energy_j"):
            assert math.fsum(group["resolved"]["released_number"]*group["resolved"][name] for group in groups)/total == pytest.approx(report["resolved"][name], rel=3e-14)
    for index, (record, component) in enumerate(zip(report["components"], mixture.components)):
        assert record["component_index"] == index
        for name in ("selection_index", "interval_index", "reservoir_index", "current_index", "width_index", "interval_origin_s"):
            assert record[name] == getattr(component, name)
        assert record["temperature_k"] == component.source_window.temperature_k
        assert record["state_id"] == component.source_window.reservoir.state_id
        assert record["released_number"] == component.source_window.manifest()["quantities"]["released_number"]["value"]
        assert record["probability"] == pytest.approx(record["released_number"]/total, rel=2e-15)
        assert record["variance_kinetic_energy_j2"] == pytest.approx(2*(KB*record["temperature_k"])**2, rel=3e-15)


def test_full_kinetic_energy_scatter_includes_within_component_flux_variance(mixture):
    report = inspect(mixture)
    values = report["resolved"]
    assert values["variance_expected_translational_energy_j2"] > 0.
    assert values["mean_conditional_kinetic_energy_variance_j2"] > 0.
    assert values["variance_kinetic_energy_j2"] == pytest.approx(
        values["variance_expected_translational_energy_j2"]+values["mean_conditional_kinetic_energy_variance_j2"], rel=3e-15)
    assert values["second_raw_kinetic_energy_j2"] == pytest.approx(
        values["variance_kinetic_energy_j2"]+values["mean_expected_translational_energy_j"]**2, rel=3e-15)


def test_huge_origin_does_not_erase_time_correlations_or_order_changes():
    parameters = dict(temperature=[[500., 1500.]], initial=[[[1e18, 2e18]]], x_edges=[0., 1., 3.], windows=[(.01, .5)])
    ordinary = inspect(source_record(origin=0., **parameters))
    far_source = source_record(origin=1e16, **parameters)
    far = inspect(far_source)
    for name in ("mean_time_offset_s", "variance_time_s2", "cov_time_x_s_m", "cov_time_expected_translational_energy_s_j"):
        assert far["resolved"][name] == ordinary["resolved"][name]
    assert far["resolved"]["cov_time_x_s_m"] != 0.
    assert not far["mean_time_anchor"]["global_mean_time_precision_qualified"]
    assert far["resolved"]["mean_time_s"] == 1e16
    coarse = inspect(far_source, quadrature_order=2)
    assert coarse["mean_time_anchor"]["interval_origin_s"] == far["mean_time_anchor"]["interval_origin_s"]
    assert coarse["mean_time_anchor"]["window_start_offset_s"] == far["mean_time_anchor"]["window_start_offset_s"]
    assert coarse["resolved"]["mean_time_offset_s"] != far["resolved"]["mean_time_offset_s"]


def test_probability_underflow_does_not_remove_finite_energy_variance():
    source = source_record(pools=[pool(activation_energy_ev=0., max_coverage_m2=1e101)],
                           initial=[[[1e-280, 1e100]]], temperature=[[1e150, 700.]],
                           x_edges=[0., 1., 2.], z_edges=[0., 1.])
    report = inspect(source)
    tiny, dominant = report["components"]
    assert tiny["probability"] == 0.
    assert tiny["probability_underflowed"]
    assert tiny["log_probability"] < -745.
    assert tiny["released_number"] > 0.
    assert tiny["fourth_raw_kinetic_energy_j4"]["overflowed"]
    n = Fraction(tiny["released_number"])
    big = Fraction(dominant["released_number"])
    a = 2*Fraction(KB)*Fraction(1e150)
    b = 2*Fraction(KB)*Fraction(700.)
    expected = n*big/(n+big)**2*(a-b)**2
    assert report["resolved"]["variance_expected_translational_energy_j2"] == float(expected)
    assert float(expected) > 1e-130


def test_extensive_count_survives_areal_release_underflow_without_floor():
    source = source_record(pools=[pool(attempt_frequency_s=1., activation_energy_ev=0.)], initial=[[[1.]]],
                           x_edges=[0., 1e75], z_edges=[0., 1e75], origin=0., duration=1000., windows=[(800., 801.)])
    quantity = source.components[0].source_window.manifest()["quantities"]
    assert quantity["released_m2"]["value"] == 0.
    assert not quantity["released_m2"]["qualified"]
    assert quantity["released_number"]["qualified"]
    report = inspect(source)
    assert report["resolved"]["released_number"] == quantity["released_number"]["value"] > 0.
    assert len(report["components"]) == 1


def test_zero_initial_groups_have_null_normalized_statistics():
    source = source_record(pools=[pool(), pool("empty", "M:empty")], initial=[[[1e18]], [[0.]]])
    report = inspect(source)
    empty = report["per_reservoir"][1]
    assert empty["resolved"]["released_number"] == 0.
    assert all(value is None for key, value in empty["resolved"].items() if key != "released_number")
    assert empty["mean_time_anchor"] is None
    assert len(report["components"]) == 1


def test_entirely_empty_source_retains_explicit_skips_and_null_groups():
    source = source_record(initial=[[[0.]]])
    report = inspect(source)
    assert report["components"] == []
    assert report["resolved"]["released_number"] == 0.
    assert report["resolved"]["variance_time_s2"] is None
    assert report["source"]["skipped_zero_initials"]


def test_detached_serializable_output_units_and_omitted_optional_components(mixture):
    report = inspect(mixture)
    serialized = json.dumps(report, allow_nan=False)
    assert json.loads(serialized) == report
    assert set(report["resolved"]) == set(report["units"]) == set(UNITS)
    assert report["units"]["cov_time_expected_translational_energy_s_j"] == "s J"
    report["source"]["components"][0]["source_window"]["temperature_k"] = -1
    report["components"][0]["timing"]["mean_local_delay_s"] = -1
    fresh = inspect(mixture, include_component_records=False)
    assert "components" not in fresh
    assert fresh["source"] == mixture.manifest()
    assert not fresh["host_budget_accepted"]


@pytest.mark.parametrize("order", [True, 1, 4097, 2.5, np.ma.array(4), float("nan")])
def test_invalid_rule_orders_reject(order):
    with pytest.raises(ValueError, match="quadrature_order"):
        inspect(None, quadrature_order=order)


@pytest.mark.parametrize("value", [1, None, np.bool_(True)])
def test_optional_component_flag_requires_explicit_boolean(value):
    with pytest.raises(ValueError, match="include_component_records"):
        inspect(None, include_component_records=value)


def test_non_source_rejects():
    with pytest.raises(ValueError, match="TemkinPlanarSourceWindows"):
        inspect({})


@pytest.mark.parametrize("field,value", [("initial_coverage_m2", 2e18), ("temperature_k", 701.)])
def test_forged_physical_source_does_not_override_canonical_baseline(field, value):
    source = source_record()
    object.__setattr__(source.components[0].source_window, field, value)
    with pytest.raises(ValueError):
        inspect(source)


def test_forged_cached_scalar_quantity_rejects_constructor_replay():
    source = source_record()
    window = source.components[0].source_window
    quantities = dict(window._quantities)
    quantities["released_number"] = replace(quantities["released_number"], value=1.)
    object.__setattr__(window, "_quantities", quantities)
    with pytest.raises(ValueError, match="modified"):
        inspect(source)


def test_positive_unrepresented_local_variance_rejects_instead_of_true_zero():
    source = source_record(pools=[pool(attempt_frequency_s=1e200, activation_energy_ev=0.)],
                           origin=0., duration=1e-190, windows=[(0., 1e-190)])
    with pytest.raises(WindowedJointMomentRangeError, match="local timing"):
        inspect(source)


def test_unrepresented_positive_extensive_release_is_not_omitted():
    source = source_record(pools=[pool(attempt_frequency_s=1., activation_energy_ev=0.)], initial=[[[1.]]],
                           x_edges=[0., 1.], z_edges=[0., 1.], origin=0., duration=1000., windows=[(800., 801.)])
    assert source.components[0].source_window.manifest()["quantities"]["released_number"]["value"] == 0.
    with pytest.raises(WindowedJointMomentRangeError, match="extensive released_number"):
        inspect(source)


def test_aggregate_only_domain_is_distinct_from_standalone_component_second_moments():
    source = source_record(pools=[pool(activation_energy_ev=0., max_coverage_m2=1e101)],
                           initial=[[[1e-280, 1e100]]], temperature=[[1e250, 700.]],
                           x_edges=[0., 1., 2.], z_edges=[0., 1.])
    aggregate = inspect(source, include_component_records=False)
    assert aggregate["resolved"]["variance_kinetic_energy_j2"] > 1e70
    assert math.isfinite(aggregate["resolved"]["variance_kinetic_energy_j2"])
    with pytest.raises(WindowedJointMomentRangeError, match="variance_kinetic_energy"):
        inspect(source)
