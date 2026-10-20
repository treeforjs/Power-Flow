"""Independent law, geometry and mixture checks of joint emitted-source moments.

Parameters are manufactured mathematical fixtures.  Analytic first-order rates
provide independent emission-time integrals; no packet sampling is performed.
"""

import copy
from decimal import Decimal, localcontext
import json
import math

import numpy as np
import pytest

from mhdlab.constants import KB, QE
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.joint_source_moments import JointSourceMomentRangeError, joint_source_moments
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_release import LAW_IDENTITY, TemkinReservoir


STATE = "manufactured:gas"
BASE_TEMPERATURE = 300.
RATE_CENTER = math.log(2.)
RATE_EDGE = 2 * RATE_CENTER
ORDER = 64


def pool(**changes):
    values = dict(
        name="manufactured-pool", state_id=STATE, attempt_frequency_s=4 * RATE_CENTER,
        activation_energy_ev=2 * RATE_CENTER * KB * BASE_TEMPERATURE / QE,
        desorption_energy_ev=.9, internal_energy_ev=.4, coverage_coefficient=0.,
        reference_coverage_m2=10., max_coverage_m2=10.,
        kinetic_provenance="manufactured two-temperature analytic rate fixture",
        energy_provenance="manufactured independently specified gas-reference costs",
    )
    values.update(changes)
    return TemkinReservoir(**values)


def make_history(x, z, temperature, coverage, *, times=None, pools=None):
    x, z = np.asarray(x, dtype=float), np.asarray(z, dtype=float)
    times = np.asarray([0., 1.] if times is None else times, dtype=float)
    shape = (len(times), len(z)-1, len(x)-1)
    temperature = np.broadcast_to(temperature, shape).copy()
    base = FlashSurfaceHistory(times, x, temperature[:, 0, :], "synthetic", "manufactured-joint",
                               "manufactured source correlation verification", "manufactured stationary face")
    attributes = base.attributes
    attributes["thermal_input_status"] = "manufactured"
    surface = PlanarSurfaceHistory(times, x, z, temperature, attributes)
    return evolve_temkin_planar_release(surface, np.asarray(coverage, dtype=float), [pool()] if pools is None else pools,
                                        parameter_provenance="manufactured joint source verification")


def symmetric_history(*, origin=0., shift_x=0., shift_z=0., split=False):
    times = np.array([0., .5, 1.]) if split else np.array([0., 1.])
    return make_history(np.array([-2., -1., 1., 2.])+shift_x, np.array([0., 1.])+shift_z,
                        [[[600., 300., 600.]]], [[[4/3, 1., 4/3]]], times=times+origin)


def result(source, **changes):
    return joint_source_moments(source, quadrature_order=ORDER, **changes)


def exponential_moments(rate):
    denominator = math.expm1(rate)
    mean = 1 / rate - 1 / denominator
    variance = 1 / rate**2 - math.exp(rate) / denominator**2
    return mean, variance


def symmetric_oracle():
    mean_c, var_c = exponential_moments(RATE_CENTER)
    mean_e, var_e = exponential_moments(RATE_EDGE)
    delta = mean_e-mean_c
    thermal_unit = KB * BASE_TEMPERATURE
    mean_time = (2 * mean_e+mean_c)/3
    return {
        "released_number": 3., "mean_time_s": mean_time,
        "variance_time_s2": (2 * var_e+var_c)/3 + 2 * delta**2/9,
        "mean_x_m": 0., "variance_x_m2": 5/3,
        "mean_z_m": .5, "variance_z_m2": 1/12,
        "mean_width_centered_square_m2": 5/3,
        "variance_width_centered_square_m4": 64/45,
        "mean_expected_translational_energy_j": 10 * thermal_unit/3,
        "variance_expected_translational_energy_j2": 8 * thermal_unit**2/9,
        "cov_time_x_s_m": 0., "cov_time_z_s_m": 0.,
        "cov_time_width_centered_square_s_m2": 4 * delta/9,
        "cov_expected_translational_energy_width_centered_square_j_m2": 8 * thermal_unit/9,
        "cov_time_expected_translational_energy_s_j": 4 * thermal_unit * delta/9,
    }


def assert_nonzero_expected(actual, expected, *, rel=3e-11):
    for key, value in expected.items():
        if value == 0:
            assert actual[key] == pytest.approx(0., rel=0., abs=2e-13)
        else:
            assert actual[key] == pytest.approx(value, rel=rel, abs=0.), key


def test_joint_time_and_rectangle_moments_match_independent_analytic_mixture():
    report = result(symmetric_history())
    assert_nonzero_expected(report["reservoirs"][0]["resolved"], symmetric_oracle())
    assert report["schema"] == "mhdlab-joint-source-moments-v1"
    assert report["law_identity"] == LAW_IDENTITY
    assert report["source_history_schema"] == "mhdlab-temkin-planar-release-v1"
    assert report["surface_attributes"]["thermal_input_status"] == "manufactured"
    json.dumps(report, allow_nan=False)


def test_symmetric_edge_heating_is_visible_in_time_width_square_covariance():
    stats = result(symmetric_history())["reservoirs"][0]["resolved"]
    assert stats["cov_time_x_s_m"] == pytest.approx(0., rel=0., abs=2e-13)
    assert stats["cov_time_width_centered_square_s_m2"] < -.01
    assert stats["cov_time_width_centered_square_s_m2"] == pytest.approx(
        symmetric_oracle()["cov_time_width_centered_square_s_m2"], rel=3e-11, abs=0.)
    # A center-only position treatment gives a different x² difference; exact
    # uniform-bin moments are 7/3 on each edge and 1/3 in the wide middle bin.
    assert stats["mean_width_centered_square_m2"] == pytest.approx(5/3, rel=2e-13, abs=0.)


def test_time_energy_covariance_uses_flux_expected_2kbt_and_separate_energy_reference():
    report = result(symmetric_history())
    stats = report["reservoirs"][0]["resolved"]
    covariance = stats["cov_time_expected_translational_energy_s_j"]
    assert covariance < 0
    assert covariance == pytest.approx(KB * BASE_TEMPERATURE * stats["cov_time_width_centered_square_s_m2"], rel=3e-11, abs=0.)
    assert stats["variance_expected_translational_energy_j2"] == pytest.approx(
        8 * (KB * BASE_TEMPERATURE)**2/9, rel=3e-11, abs=0.)
    # This is variance across component expected energies, not the larger
    # variance of sampled kinetic energies in their local Maxwellian mixtures.
    expected_mean = 10 * KB * BASE_TEMPERATURE/3
    assert stats["mean_expected_translational_energy_j"] == pytest.approx(expected_mean, rel=2e-13, abs=0.)
    assert expected_mean != pool().activation_energy_ev * QE
    assert expected_mean != pool().desorption_energy_ev * QE


def test_constant_temperature_has_zero_between_component_expected_energy_variance():
    source = make_history([-3., -1., 2.], [0., 1., 4.], 600., [[[.1, .2], [.3, .7]]],
                          times=[0., .25, 1.])
    entry = result(source)["reservoirs"][0]
    groups = [entry]
    groups.extend(entry["conditional_width"])
    groups.extend(entry["conditional_interval"])
    for group in groups:
        for hypothesis in ("resolved", "uniform_width"):
            stats = group[hypothesis]
            assert stats["mean_expected_translational_energy_j"] == pytest.approx(2*KB*600, rel=3e-13, abs=0.)
            assert stats["variance_expected_translational_energy_j2"] == 0.
            assert stats["cov_time_expected_translational_energy_s_j"] == 0.
            assert stats["cov_expected_translational_energy_width_centered_square_j_m2"] == 0.


def test_uniform_width_removes_spatial_correlations_and_keeps_time_energy_mixture():
    entry = result(symmetric_history())["reservoirs"][0]
    resolved, uniform = entry["resolved"], entry["uniform_width"]
    common = ("released_number", "mean_time_s", "variance_time_s2", "mean_z_m", "variance_z_m2",
              "mean_expected_translational_energy_j", "variance_expected_translational_energy_j2",
              "cov_time_z_s_m", "cov_time_expected_translational_energy_s_j")
    for name in common:
        assert uniform[name] == pytest.approx(resolved[name], rel=3e-13, abs=0.)
    assert uniform["mean_x_m"] == 0.
    assert uniform["variance_x_m2"] == pytest.approx(4/3, rel=2e-13, abs=0.)
    assert uniform["mean_width_centered_square_m2"] == pytest.approx(4/3, rel=2e-13, abs=0.)
    assert uniform["variance_width_centered_square_m4"] == pytest.approx(64/45, rel=2e-13, abs=0.)
    for name in ("cov_time_x_s_m", "cov_time_width_centered_square_s_m2",
                 "cov_expected_translational_energy_width_centered_square_j_m2"):
        assert uniform[name] == pytest.approx(0., rel=0., abs=0.)


def test_width_conditionals_use_actual_bin_width_and_fixed_full_face_square_reference():
    entry = result(symmetric_history())["reservoirs"][0]
    expected_means = [exponential_moments(RATE_EDGE)[0], exponential_moments(RATE_CENTER)[0], exponential_moments(RATE_EDGE)[0]]
    for index, conditional in enumerate(entry["conditional_width"]):
        resolved, uniform = conditional["resolved"], conditional["uniform_width"]
        width = conditional["x_right_m"]-conditional["x_left_m"]
        assert conditional["width_index"] == index
        assert resolved["released_number"] == pytest.approx(1., rel=3e-13, abs=0.)
        assert uniform["released_number"] == pytest.approx(3 * width/4, rel=3e-13, abs=0.)
        assert resolved["mean_time_s"] == pytest.approx(expected_means[index], rel=3e-11, abs=0.)
        assert uniform["mean_time_s"] == pytest.approx(symmetric_oracle()["mean_time_s"], rel=3e-11, abs=0.)
        expected_square, expected_square_variance = ((1/3, 4/45) if index == 1 else (7/3, 34/45))
        for stats in (resolved, uniform):
            assert stats["mean_width_centered_square_m2"] == pytest.approx(expected_square, rel=3e-13, abs=0.)
            assert stats["variance_width_centered_square_m4"] == pytest.approx(expected_square_variance, rel=3e-13, abs=0.)
    assert entry["conditional_width"][0]["resolved"]["mean_expected_translational_energy_j"] == pytest.approx(2 * KB * 600, rel=3e-13, abs=0.)
    assert entry["conditional_width"][1]["resolved"]["mean_expected_translational_energy_j"] == pytest.approx(2 * KB * 300, rel=3e-13, abs=0.)


def time_window_integrals(rate, start, stop):
    before, after = math.exp(-rate * start), math.exp(-rate * stop)
    number = before-after
    first = (start+1/rate)*before-(stop+1/rate)*after
    second = (start**2+2*start/rate+2/rate**2)*before-(stop**2+2*stop/rate+2/rate**2)*after
    return number, first, second


def test_interval_conditionals_match_independent_exponential_window_integrals():
    entry = result(symmetric_history(split=True))["reservoirs"][0]
    for interval in entry["conditional_interval"]:
        start, stop = interval["time_start_s"], interval["time_stop_s"]
        # Initial extensive populations are 4/3 in each edge and 2 in the
        # center. Integrate the original unconditioned exponential law directly.
        edge = np.array(time_window_integrals(RATE_EDGE, start, stop)) * (8/3)
        center = np.array(time_window_integrals(RATE_CENTER, start, stop)) * 2
        number, first, second = edge+center
        mean = first/number
        variance = second/number-mean**2
        energy_first = edge[1]*(2*KB*600)+center[1]*(2*KB*300)
        energy_mean = (edge[0]*(2*KB*600)+center[0]*(2*KB*300))/number
        expected_covariance = energy_first/number-mean*energy_mean
        for name in ("resolved", "uniform_width"):
            stats = interval[name]
            assert stats["released_number"] == pytest.approx(number, rel=3e-13, abs=0.)
            assert stats["mean_time_s"] == pytest.approx(mean, rel=3e-11, abs=0.)
            assert stats["variance_time_s2"] == pytest.approx(variance, rel=3e-10, abs=0.)
            assert stats["cov_time_expected_translational_energy_s_j"] == pytest.approx(expected_covariance, rel=3e-9, abs=0.)


@pytest.mark.parametrize("conditionals", ["conditional_width", "conditional_interval"])
def test_conditional_count_and_first_moment_sums_reconstruct_global_population(conditionals):
    entry = result(symmetric_history(split=True))["reservoirs"][0]
    for hypothesis in ("resolved", "uniform_width"):
        global_stats = entry[hypothesis]
        records = [item[hypothesis] for item in entry[conditionals]]
        total = math.fsum(item["released_number"] for item in records)
        assert total == pytest.approx(global_stats["released_number"], rel=3e-13, abs=0.)
        for name in ("mean_time_s", "mean_x_m", "mean_z_m", "mean_width_centered_square_m2",
                     "mean_expected_translational_energy_j"):
            first = math.fsum(item["released_number"]*item[name] for item in records)
            expected = total*global_stats[name]
            assert first == pytest.approx(expected, rel=3e-11, abs=2e-13 if expected == 0 else 0.)


def test_current_patch_area_is_used_once_in_joint_time_current_covariance():
    source = make_history([10., 12.], [0., 1., 4.], [[[300.], [600.]]], [[[1.], [2/3]]])
    entry = result(source)["reservoirs"][0]
    mean_c, _ = exponential_moments(RATE_CENTER)
    mean_e, _ = exponential_moments(RATE_EDGE)
    expected = 3/8*(mean_e-mean_c)
    for hypothesis in ("resolved", "uniform_width"):
        stats = entry[hypothesis]
        assert stats["released_number"] == pytest.approx(4., rel=3e-13, abs=0.)
        assert stats["mean_z_m"] == pytest.approx(2., rel=3e-13, abs=0.)
        assert stats["cov_time_z_s_m"] == pytest.approx(expected, rel=3e-11, abs=0.)
        assert stats["cov_time_width_centered_square_s_m2"] == pytest.approx(0., rel=0., abs=2e-13)


def test_large_time_and_spatial_origins_preserve_centered_covariances():
    baseline = result(symmetric_history(split=True))["reservoirs"][0]
    shifted_report = result(symmetric_history(origin=1e9, shift_x=1e12, shift_z=2e12, split=True))
    shifted = shifted_report["reservoirs"][0]
    for hypothesis in ("resolved", "uniform_width"):
        for name, expected in baseline[hypothesis].items():
            if name in ("mean_time_s", "mean_x_m", "mean_z_m"):
                continue
            assert shifted[hypothesis][name] == pytest.approx(expected, rel=3e-10, abs=2e-13 if expected == 0 else 0.), name
        assert shifted[hypothesis]["mean_time_s"] == pytest.approx(1e9+baseline[hypothesis]["mean_time_s"], rel=0., abs=np.spacing(1e9))
    assert shifted_report["coordinate_reference"]["x_face_center_m"] == 1e12
    assert shifted_report["coordinate_reference"]["time_origin_s"] == 1e9


def test_interval_splitting_preserves_global_joint_time_energy_position_law():
    single = result(symmetric_history())["reservoirs"][0]
    split = result(symmetric_history(split=True))["reservoirs"][0]
    for hypothesis in ("resolved", "uniform_width"):
        for name, expected in single[hypothesis].items():
            if name == "mean_x_m":
                # Symmetry gives the analytic zero. A rounded nonzero baseline
                # is not an exact nonzero reference for a relative comparison.
                assert expected == pytest.approx(0., rel=0., abs=2e-13)
                assert split[hypothesis][name] == pytest.approx(0., rel=0., abs=2e-13)
                continue
            assert split[hypothesis][name] == pytest.approx(expected, rel=3e-11, abs=2e-13 if expected == 0 else 0.), name
    assert_nonzero_expected(split["resolved"], symmetric_oracle())


def test_empty_width_bin_and_pool_have_explicitly_undefined_distributions():
    empty = pool(name="empty", state_id="manufactured:empty")
    source = make_history([-2., -1., 1., 2.], [0., 1.], [[[600., 300., 600.]]],
                          [[[4/3, 0., 4/3]], [[0., 0., 0.]]], pools=[pool(), empty])
    report = result(source)
    first, second = report["reservoirs"]
    center = first["conditional_width"][1]
    assert center["resolved"]["released_number"] == 0.
    assert center["resolved"]["mean_time_s"] is None
    assert center["uniform_width"]["released_number"] > 0
    for hypothesis in ("resolved", "uniform_width"):
        stats = second[hypothesis]
        assert stats["released_number"] == 0.
        assert all(value is None for name, value in stats.items() if name != "released_number")
    json.dumps(report, allow_nan=False)


def test_component_records_are_explicit_and_retain_actual_area_and_time_law():
    source = symmetric_history()
    compact = result(source)
    assert "component_records" not in compact
    detailed = result(source, include_component_records=True)
    records = detailed["component_records"]
    assert len(records) == detailed["diagnostics"]["positive_components"] == 3
    assert detailed["diagnostics"]["component_time_underflows"] == []
    for record in records:
        assert record["released_number"] == pytest.approx(1., rel=3e-13, abs=0.)
        assert record["expected_translational_energy_j"] == pytest.approx(2*KB*record["temperature_k"], rel=3e-13, abs=0.)
        assert record["timing"]["law_identity"] == LAW_IDENTITY
        assert record["timing"]["quadrature_order"] == ORDER


def test_covariance_units_and_source_context_are_explicit():
    report = result(symmetric_history())
    expected = {"mean_time_s": "s", "variance_time_s2": "s^2", "mean_x_m": "m",
                "mean_width_centered_square_m2": "m^2", "variance_width_centered_square_m4": "m^4",
                "mean_expected_translational_energy_j": "J", "variance_expected_translational_energy_j2": "J^2",
                "cov_time_x_s_m": "s m", "cov_time_width_centered_square_s_m2": "s m^2",
                "cov_time_expected_translational_energy_s_j": "s J",
                "cov_expected_translational_energy_width_centered_square_j_m2": "J m^2"}
    for name, unit in expected.items():
        assert report["units"][name] == unit
    assert report["surface_attributes"]["source_kind"] == "synthetic"
    assert report["parameter_provenance"].startswith("manufactured")


@pytest.mark.parametrize("order", [0, 1, True, "16", 16+0j, np.nan, np.inf, np.ma.array(16, mask=True)])
def test_quadrature_order_requires_unmasked_supported_integer(order):
    with pytest.raises(ValueError):
        joint_source_moments(symmetric_history(), quadrature_order=order)


@pytest.mark.parametrize("value", [1, "yes", None, np.array([True])])
def test_component_record_flag_requires_explicit_boolean(value):
    with pytest.raises(ValueError):
        joint_source_moments(symmetric_history(), quadrature_order=ORDER, include_component_records=value)


@pytest.mark.parametrize("field,value", [
    ("law_identity", "constant_hazard"), ("time_model", "uniform_in_interval"),
    ("schema", "unqualified"), ("coverage_model", "shared_sites_without_model"),
])
def test_history_semantics_are_revalidated(field, value):
    forged = copy.copy(symmetric_history())
    object.__setattr__(forged, field, value)
    with pytest.raises(ValueError):
        result(forged)


def test_invalid_depth_units_and_masked_source_data_cannot_bypass_revalidation():
    source = symmetric_history()
    forged = copy.copy(source)
    attributes = dict(source.surface_attributes)
    attributes["material_sample_depth_units"] = "cm"
    object.__setattr__(forged, "surface_attributes", attributes)
    with pytest.raises(ValueError):
        result(forged)
    forged = copy.copy(source)
    object.__setattr__(forged, "temperature_k", np.ma.array(source.temperature_k, mask=True))
    with pytest.raises(ValueError):
        result(forged)
    with pytest.raises(ValueError):
        joint_source_moments(object(), quadrature_order=ORDER)


def test_unrepresentable_energy_variance_rejects_instead_of_claiming_constant_energy():
    tiny_temperature = 1e-200/(2*KB)
    record = pool(activation_energy_ev=0., desorption_energy_ev=0., internal_energy_ev=0.)
    source = make_history([-1., 0., 1.], [0., 1.], [[[tiny_temperature, 2*tiny_temperature]]],
                          [[[1., 1.]]], pools=[record])
    # Both expected energies are positive/representable, while their positive
    # between-component variance is below the minimum float64 magnitude.
    assert 2*KB*tiny_temperature > 0
    assert (2*KB*tiny_temperature)**2 == 0
    with pytest.raises(JointSourceMomentRangeError):
        result(source)


def test_underflow_fourth_spatial_moment_is_explicitly_outside_supported_range():
    source = make_history([0., 1e-100], [0., 1.], 600., [[[1.]]])
    assert 1e-100**2 > 0
    assert 1e-100**4 == 0
    with pytest.raises(JointSourceMomentRangeError):
        result(source)


def test_tiny_middle_width_patch_retains_position_and_square_moment():
    source = make_history([-1., 1e-20, 2e-20, 1.], [0., 1.], 600., [[[0., 1., 0.]]])
    stats = result(source)["reservoirs"][0]["resolved"]
    assert stats["mean_x_m"] == pytest.approx(1.5e-20, rel=3e-13, abs=0.)
    assert stats["variance_x_m2"] == pytest.approx(1e-40/12, rel=3e-13, abs=0.)
    assert stats["mean_width_centered_square_m2"] == pytest.approx(7e-40/3, rel=3e-13, abs=0.)


def test_tiny_current_cluster_retains_mean_and_centered_variance():
    source = make_history([-.5, .5], [-1., 1e-20, 2e-20, 3e-20, 1.], 600.,
                          [[[0.], [1.], [1.], [0.]]])
    entry = result(source)["reservoirs"][0]
    for hypothesis in ("resolved", "uniform_width"):
        assert entry[hypothesis]["mean_z_m"] == pytest.approx(2e-20, rel=3e-13, abs=0.)
        assert entry[hypothesis]["variance_z_m2"] == pytest.approx((2e-20)**2/12, rel=3e-13, abs=0.)


def test_asymmetric_endpoint_midpoint_retains_below_span_rounding_mean():
    upper = np.nextafter(1., math.inf)
    source = make_history([-1., upper], [-1., upper], 600., [[[1.]]])
    report = result(source)
    midpoint = 2.**-53
    for hypothesis in ("resolved", "uniform_width"):
        stats = report["reservoirs"][0][hypothesis]
        assert stats["mean_x_m"] == midpoint
        assert stats["mean_z_m"] == midpoint
    assert report["coordinate_reference"]["x_face_center_m"] == midpoint


def decimal_time_energy_mixture(source, component_records):
    """Reassemble exposed component quadrature moments in high precision.

    This verifies joint-mixture arithmetic and preservation of start+offset.
    It does not independently rederive conditional-time quadrature estimates.
    """
    with localcontext() as context:
        context.prec = 96
        convert = Decimal.from_float
        counts = [convert(record["released_number"]) for record in component_records]
        total = sum(counts)
        weights = [number/total for number in counts]
        times = [convert(float(source.time_edges_s[record["interval_index"]]))
                 +convert(record["timing"]["mean_offset_s"]) for record in component_records]
        energies = [convert(record["expected_translational_energy_j"]) for record in component_records]
        mean_time = sum(weight*time for weight, time in zip(weights, times))
        mean_energy = sum(weight*energy for weight, energy in zip(weights, energies))
        covariance = sum(weight*(time-mean_time)*(energy-mean_energy)
                         for weight, time, energy in zip(weights, times, energies))
        variance = sum(weight*(convert(record["timing"]["variance_offset_s2"])+(time-mean_time)**2)
                       for weight, record, time in zip(weights, component_records, times))
    return float(covariance), float(variance), [float(weight) for weight in weights]


def test_late_local_offsets_preserve_time_energy_correlation_below_absolute_time_spacing():
    record = pool(attempt_frequency_s=1e20, activation_energy_ev=1e4,
                  reference_coverage_m2=1e18, max_coverage_m2=1e18,
                  desorption_energy_ev=0., internal_energy_ev=0.)
    source = make_history([0., .5, 1.], [0., 1.], [[[1e8, 2e8]]], [[[1e18, 1e18]]],
                          times=[1., 1.+1e-12], pools=[record])
    report = joint_source_moments(source, quadrature_order=32, include_component_records=True)
    components = report["component_records"]
    assert len(components) == 2
    assert all(1.+component["timing"]["mean_offset_s"] == 1. for component in components)
    expected_covariance, expected_variance, _ = decimal_time_energy_mixture(source, components)
    assert expected_covariance < 0 and expected_variance > 0
    for hypothesis in ("resolved", "uniform_width"):
        stats = report["reservoirs"][0][hypothesis]
        assert stats["mean_time_s"] == 1.
        assert stats["cov_time_expected_translational_energy_s_j"] == pytest.approx(expected_covariance, rel=3e-12, abs=0.)
        assert stats["variance_time_s2"] == pytest.approx(expected_variance, rel=3e-12, abs=0.)


def test_rare_early_component_preserves_dominant_late_local_covariance():
    record = pool(attempt_frequency_s=1e20, activation_energy_ev=1e4,
                  reference_coverage_m2=1e18, max_coverage_m2=1e18,
                  desorption_energy_ev=0., internal_energy_ev=0.)
    temperature = np.array([[[8e5, 8e5]], [[8e5, 8e5]], [[1e8, 2e8]]])
    source = make_history([0., .5, 1.], [0., 1.], temperature, [[[1e18, 1e18]]],
                          times=[0., 1., 1.+1e-12], pools=[record])
    report = joint_source_moments(source, quadrature_order=32, include_component_records=True)
    components = report["component_records"]
    assert len(components) == 4
    expected_covariance, expected_variance, weights = decimal_time_energy_mixture(source, components)
    early_probability = sum(weight for weight, component in zip(weights, components) if component["interval_index"] == 0)
    assert 1e-44 < early_probability < 1e-42
    assert expected_covariance < -1e-36
    for hypothesis in ("resolved", "uniform_width"):
        stats = report["reservoirs"][0][hypothesis]
        assert stats["cov_time_expected_translational_energy_s_j"] == pytest.approx(expected_covariance, rel=3e-12, abs=0.)
        assert stats["variance_time_s2"] == pytest.approx(expected_variance, rel=3e-12, abs=0.)
