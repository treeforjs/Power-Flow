"""Independent analytic and numerical checks of the manufactured periodic mode."""

from dataclasses import FrozenInstanceError, replace
from fractions import Fraction
import math

import numpy as np
import pytest
from scipy.integrate import quad

from mhdlab.conductive_heating_benchmark import (
    MAX_ABSOLUTE_PHASE_RAD,
    PeriodicConductiveHeating,
    UnsupportedConductiveBenchmarkRangeError,
)


@pytest.fixture
def mode():
    # Manufactured coefficients; these are not a stainless-steel property set.
    return PeriodicConductiveHeating(2.0, 3.0, 5.0, 7.0, 0.4, 1.2)


def direct(model, y, t):
    sigma, mu, cv, kap, k, b = (
        model.conductivity_s_m, model.permeability_h_m,
        model.volumetric_heat_capacity_j_m3_k, model.thermal_conductivity_w_m_k,
        model.wavenumber_rad_m, model.initial_field_t,
    )
    lam, delta = 2 * k * k / (mu * sigma), 4 * kap * k * k / cv
    field = b * math.exp(-lam * t / 2) * math.cos(k * y)
    current = b * k / mu * math.exp(-lam * t / 2) * math.sin(k * y)
    # Positive time-convolution oracle instead of the module's closed form.
    thermal = quad(
        lambda s: b * b * k * k / (2 * mu * mu * sigma * cv)
        * math.exp(-lam * s) * (1 - math.exp(-delta * (t - s)) * math.cos(2 * k * y)),
        0, t, epsabs=1e-13, epsrel=1e-12,
    )[0]
    return np.array([field, current, current * current / sigma, thermal])


def values(record):
    return np.array([record.magnetic_field_t, record.current_density_a_m2,
                     record.joule_power_w_m3, record.temperature_rise_k])


def test_properties_and_dimensions(mode):
    assert mode.magnetic_diffusivity_m2_s == pytest.approx(1 / 6)
    assert mode.thermal_diffusivity_m2_s == pytest.approx(7 / 5)
    assert mode.resistive_energy_decay_rate_s == pytest.approx(2 * 0.4**2 / 6)
    assert mode.thermal_mode_decay_rate_s == pytest.approx(4 * 7 * 0.4**2 / 5)
    assert mode.period_m == pytest.approx(2 * math.pi / 0.4)


@pytest.mark.parametrize("time", [0.0, 0.1, 2.0, 30.0])
@pytest.mark.parametrize("position", [-2.3, 0.0, 0.8, 6.7])
def test_point_fields_against_independent_convolution(mode, time, position):
    np.testing.assert_allclose(values(mode.fields(position, time)), direct(mode, position, time), rtol=2e-12, atol=5e-15)


def test_initial_mode_and_signed_field(mode):
    y = np.linspace(-mode.period_m, mode.period_m, 17)
    fields = mode.fields(y, 0)
    np.testing.assert_allclose(fields.magnetic_field_t, mode.initial_field_t * np.cos(mode.wavenumber_rad_m * y), atol=1e-15)
    assert np.all(fields.temperature_rise_k == 0)
    negative = replace(mode, initial_field_t=-mode.initial_field_t).fields(y, 2)
    positive = mode.fields(y, 2)
    np.testing.assert_array_equal(negative.magnetic_field_t, -positive.magnetic_field_t)
    np.testing.assert_array_equal(negative.current_density_a_m2, -positive.current_density_a_m2)
    np.testing.assert_array_equal(negative.joule_power_w_m3, positive.joule_power_w_m3)
    np.testing.assert_array_equal(negative.temperature_rise_k, positive.temperature_rise_k)


def test_induction_and_ampere_derivatives(mode):
    y, t, dy, dt = 1.7, 2.0, 2e-4, 2e-4
    center = mode.fields(y, t)
    left, right = mode.fields(y - dy, t), mode.fields(y + dy, t)
    b_y = (right.magnetic_field_t - left.magnetic_field_t) / (2 * dy)
    b_yy = (right.magnetic_field_t - 2 * center.magnetic_field_t + left.magnetic_field_t) / dy**2
    b_t = (mode.fields(y, t + dt).magnetic_field_t - mode.fields(y, t - dt).magnetic_field_t) / (2 * dt)
    assert center.current_density_a_m2 == pytest.approx(-b_y / mode.permeability_h_m, rel=2e-8)
    assert b_t == pytest.approx(mode.magnetic_diffusivity_m2_s * b_yy, rel=2e-7)


@pytest.mark.parametrize("kappa", [0.0, 5 / 12, 7.0])
def test_thermal_pde_including_resonance(mode, kappa):
    model = replace(mode, thermal_conductivity_w_m_k=kappa)
    y, t, dy, dt = 1.1, 3.0, 3e-4, 2e-4
    center = model.fields(y, t)
    temperature_t = (model.fields(y, t + dt).temperature_rise_k - model.fields(y, t - dt).temperature_rise_k) / (2 * dt)
    temperature_yy = (model.fields(y + dy, t).temperature_rise_k - 2 * center.temperature_rise_k + model.fields(y - dy, t).temperature_rise_k) / dy**2
    assert model.volumetric_heat_capacity_j_m3_k * temperature_t == pytest.approx(
        model.thermal_conductivity_w_m_k * temperature_yy + center.joule_power_w_m3, rel=2e-6, abs=1e-12,
    )


@pytest.mark.parametrize("kappa", [0.0, 5 / 12, 7.0])
def test_exact_cell_means_against_independent_quadrature(mode, kappa):
    model = replace(mode, thermal_conductivity_w_m_k=kappa)
    edges, t = np.array([-3.7, -0.1, 0.0, 0.002, 1.4, 8.1]), 4.0
    actual = values(model.cell_averages(edges, t))
    for cell, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
        expected = [quad(lambda y: direct(model, y, t)[column], left, right, epsabs=1e-13)[0] / (right - left) for column in range(4)]
        np.testing.assert_allclose(actual[:, cell], expected, rtol=5e-11, atol=1e-14)


def test_cell_ampere_and_joule_square_are_distinct(mode):
    edges = np.linspace(-mode.period_m, mode.period_m, 25)
    fields, cells = mode.fields(edges, 3), mode.cell_averages(edges, 3)
    ampere = -np.diff(fields.magnetic_field_t) / (mode.permeability_h_m * np.diff(edges))
    np.testing.assert_allclose(cells.current_density_a_m2, ampere, rtol=3e-14, atol=3e-17)
    whole = mode.cell_averages([0, mode.period_m], 3)
    assert abs(whole.current_density_a_m2[0]) < 1e-30
    assert whole.joule_power_w_m3[0] > 1e-3
    assert whole.joule_power_w_m3[0] > whole.current_density_a_m2[0]**2 / mode.conductivity_s_m


@pytest.mark.parametrize("time", [0.0, 1e-8, 0.5, 50.0, 1e5])
def test_periodic_energy_closure_and_independent_spatial_integrals(mode, time):
    budget = mode.energy_budget(time)
    period = mode.period_m
    magnetic = quad(lambda y: float(mode.fields(y, time).magnetic_field_t)**2 / (2 * mode.permeability_h_m), 0, period, epsabs=1e-12)[0] / period
    thermal = quad(lambda y: mode.volumetric_heat_capacity_j_m3_k * float(mode.fields(y, time).temperature_rise_k), 0, period, epsabs=1e-12)[0] / period
    power = quad(lambda y: float(mode.fields(y, time).joule_power_w_m3), 0, period, epsabs=1e-12)[0] / period
    assert budget.mean_magnetic_energy_j_m3 == pytest.approx(magnetic, rel=3e-13, abs=1e-15)
    assert budget.mean_thermal_energy_j_m3 == pytest.approx(thermal, rel=3e-13, abs=1e-15)
    assert budget.mean_joule_power_w_m3 == pytest.approx(power, rel=3e-13, abs=1e-15)
    assert budget.mean_total_energy_j_m3 == pytest.approx(budget.initial_mean_magnetic_energy_j_m3, rel=2e-14)
    assert abs(budget.energy_residual_j_m3) < 3e-15


def test_mean_heating_rate_and_conduction_redistribution(mode):
    t, dt = 7.0, 1e-4
    rate = (mode.energy_budget(t + dt).mean_thermal_energy_j_m3 - mode.energy_budget(t - dt).mean_thermal_energy_j_m3) / (2 * dt)
    assert rate == pytest.approx(mode.energy_budget(t).mean_joule_power_w_m3, rel=3e-10)
    adiabatic = replace(mode, thermal_conductivity_w_m_k=0)
    assert adiabatic.energy_budget(t) == mode.energy_budget(t)
    assert float(mode.fields(0, t).temperature_rise_k) > float(adiabatic.fields(0, t).temperature_rise_k)
    assert float(mode.fields(mode.period_m / 4, t).temperature_rise_k) < float(adiabatic.fields(mode.period_m / 4, t).temperature_rise_k)


def test_exact_resonance_and_nearby_values(mode):
    resonant = replace(mode, thermal_conductivity_w_m_k=mode.volumetric_heat_capacity_j_m3_k / (2 * mode.permeability_h_m * mode.conductivity_s_m))
    x, y, t = resonant.resistive_energy_decay_rate_s * 3, 1.2, 3.0
    a = mode.initial_field_t**2 / (4 * mode.permeability_h_m * mode.volumetric_heat_capacity_j_m3_k)
    expected = a * (-math.expm1(-x) - x * math.exp(-x) * math.cos(2 * mode.wavenumber_rad_m * y))
    assert float(resonant.fields(y, t).temperature_rise_k) == pytest.approx(expected, rel=2e-14)
    for factor in [1 - 1e-10, 1 + 1e-10]:
        neighbor = replace(resonant, thermal_conductivity_w_m_k=resonant.thermal_conductivity_w_m_k * factor)
        assert float(neighbor.fields(y, t).temperature_rise_k) == pytest.approx(expected, rel=2e-10)


def test_early_node_temperature_and_late_limits(mode):
    t = 1e-9
    a = mode.initial_field_t**2 / (4 * mode.permeability_h_m * mode.volumetric_heat_capacity_j_m3_k)
    expected = a * mode.resistive_energy_decay_rate_s * mode.thermal_mode_decay_rate_s * t * t / 2
    assert float(mode.fields(0, t).temperature_rise_k) == pytest.approx(expected, rel=1e-9)
    np.testing.assert_allclose(mode.fields([0, 1, 2, 3], 1e6).temperature_rise_k, a, rtol=2e-14)
    adiabatic = replace(mode, thermal_conductivity_w_m_k=0)
    positions = np.array([0, 1, 2, 3])
    np.testing.assert_allclose(adiabatic.fields(positions, 1e6).temperature_rise_k, 2 * a * np.sin(mode.wavenumber_rad_m * positions)**2, rtol=2e-14)


def test_compensated_tiny_node_temperature():
    model = PeriodicConductiveHeating(1e308, 1e-300, 1e-300, 2.5e-309, 1e-140, 1e300)
    time = 1e-300
    sigma, mu, cv, kap, k, b = map(Fraction, (model.conductivity_s_m, model.permeability_h_m, model.volumetric_heat_capacity_j_m3_k, model.thermal_conductivity_w_m_k, model.wavenumber_rad_m, model.initial_field_t))
    x, z = 2 * k * k * Fraction(time) / (mu * sigma), 4 * kap * k * k * Fraction(time) / cv
    # Both x,z round to zero individually, but the scaled O(x*z) rise is finite.
    expected = float(b * b / (4 * mu * cv) * x * z / 2)
    assert float(x) == float(z) == 0
    assert float(model.fields(0, time).temperature_rise_k) == pytest.approx(expected, rel=3e-12)


def test_compensated_tiny_cell_geometry():
    model = PeriodicConductiveHeating(1e-300, 1e-300, 1, 0, 1e-200, 1e100)
    width = 1e-200
    cells = model.cell_averages([0, width], 0)
    sigma, mu, k, b, w = map(Fraction, (model.conductivity_s_m, model.permeability_h_m, model.wavenumber_rad_m, model.initial_field_t, width))
    expected_j = float(b * k * k * w / (2 * mu))
    expected_q = float(b * b * k**4 * w * w / (3 * mu * mu * sigma))
    assert cells.current_density_a_m2[0] == pytest.approx(expected_j, rel=3e-12, abs=0)
    assert cells.joule_power_w_m3[0] == pytest.approx(expected_q, rel=3e-12, abs=0)


def test_independent_current_survives_underflowed_magnetic_amplitude():
    model = PeriodicConductiveHeating(1, 1e-300, 1e300, 1, 1e13, 0.1)
    time, y = 1e-323, math.pi / (4 * model.wavenumber_rad_m)
    fields = model.fields(y, time)
    assert float(fields.magnetic_field_t) == 0
    assert float(fields.current_density_a_m2) > 0
    assert float(fields.joule_power_w_m3) > 0
    assert float(fields.joule_power_w_m3) == pytest.approx(float(fields.current_density_a_m2)**2, rel=3e-12, abs=0)
    with pytest.raises(UnsupportedConductiveBenchmarkRangeError):
        _ = model.resistive_energy_decay_rate_s


def test_zero_field_all_diagnostics(mode):
    zero = replace(mode, initial_field_t=0)
    assert np.all(values(zero.fields([-1, 0, 1], 3)) == 0)
    assert np.all(values(zero.cell_averages([-1, 0, 1], 3)) == 0)
    assert set(vars_dict(zero.energy_budget(3)).values()) == {0.0}


def vars_dict(record):
    return {name: getattr(record, name) for name in record.__dataclass_fields__}


def test_output_immutability_and_shapes(mode):
    for result in [mode.fields([0, 1], 2), mode.cell_averages([0, 1, 2], 2)]:
        for value in vars_dict(result).values():
            assert not value.flags.writeable
            with pytest.raises(ValueError):
                value.setflags(write=True)
    assert mode.fields(0, 2).magnetic_field_t.shape == ()
    with pytest.raises(FrozenInstanceError):
        mode.initial_field_t = 4


@pytest.mark.parametrize("field", ["conductivity_s_m", "permeability_h_m", "volumetric_heat_capacity_j_m3_k", "wavenumber_rad_m"])
@pytest.mark.parametrize("bad", [0, -1, math.inf, math.nan, True, "1", 1 + 0j])
def test_positive_coefficient_validation(mode, field, bad):
    with pytest.raises(ValueError):
        replace(mode, **{field: bad})


@pytest.mark.parametrize("bad", [-1, math.inf, math.nan, True, "1", 1 + 0j])
def test_kappa_validation(mode, bad):
    with pytest.raises(ValueError):
        replace(mode, thermal_conductivity_w_m_k=bad)


@pytest.mark.parametrize("bad", [-1, math.inf, math.nan, True, "1", 1 + 0j])
def test_time_validation(mode, bad):
    for method, argument in [(mode.fields, [0]), (mode.cell_averages, [0, 1])]:
        with pytest.raises(ValueError):
            method(argument, bad)
    with pytest.raises(ValueError):
        mode.energy_budget(bad)


@pytest.mark.parametrize("bad", [[], [[0]], [True], ["1"], [1 + 0j], [math.inf], [math.nan], np.ma.array([0], mask=[True])])
def test_coordinate_validation(mode, bad):
    with pytest.raises(ValueError):
        mode.fields(bad, 1)


@pytest.mark.parametrize("bad", [[0], 1, [0, 0], [1, 0], [[0, 1]], [0, math.inf]])
def test_edge_validation(mode, bad):
    with pytest.raises(ValueError):
        mode.cell_averages(bad, 1)


def test_explicit_phase_and_output_range_guards(mode):
    with pytest.raises(UnsupportedConductiveBenchmarkRangeError, match="k.y"):
        mode.fields((MAX_ABSOLUTE_PHASE_RAD + 1) / mode.wavenumber_rad_m, 1)
    huge = replace(mode, initial_field_t=1e308, permeability_h_m=1e-300)
    with pytest.raises(UnsupportedConductiveBenchmarkRangeError, match="result"):
        huge.energy_budget(1)
    with pytest.raises(UnsupportedConductiveBenchmarkRangeError, match="result"):
        huge.fields(1, 0)
