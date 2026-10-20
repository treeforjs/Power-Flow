"""Physical orbit oracles independent of transport or reaction integration."""

import numpy as np
import pytest

from mhdlab.kinetic_orbits import (
    SPEED_OF_LIGHT_M_S,
    advance_orbits,
    kinetic_energy_j,
    velocity_m_s,
)


ELECTRON_MASS = 9.1093837139e-31
PROTON_MASS = 1.67262192595e-27
ELEMENTARY_CHARGE = 1.602176634e-19


def inputs(**overrides):
    values = dict(
        position_m=np.zeros((2, 3)),
        momentum_kg_m_s=np.array([[1e-24, 2e-24, 3e-24], [3e-24, -2e-24, 1e-24]]),
        mass_kg=np.array([ELECTRON_MASS, PROTON_MASS]),
        charge_c=np.array([-ELEMENTARY_CHARGE, ELEMENTARY_CHARGE]),
        dt_s=1e-12,
        electric_field_v_m=np.zeros((2, 3)),
        magnetic_field_t=np.zeros((2, 3)),
    )
    values.update(overrides)
    return values


@pytest.mark.parametrize("relativistic", [False, True])
def test_zero_charge_is_analytic_ballistic_motion_in_arbitrary_fields(relativistic):
    case = inputs(charge_c=np.zeros(2), electric_field_v_m=np.full((2, 3), 1e8), magnetic_field_t=np.full((2, 3), 3.0))
    initial_position = case["position_m"].copy()
    initial_momentum = case["momentum_kg_m_s"].copy()
    if relativistic:
        energy_total = np.sqrt((case["mass_kg"] * SPEED_OF_LIGHT_M_S**2)**2 + np.sum(initial_momentum**2, axis=1) * SPEED_OF_LIGHT_M_S**2)
        velocity = initial_momentum * (SPEED_OF_LIGHT_M_S**2 / energy_total)[:, None]
    else:
        velocity = initial_momentum / case["mass_kg"][:, None]
    result = advance_orbits(**case, relativistic=relativistic)
    np.testing.assert_allclose(result.position_m, initial_position + case["dt_s"] * velocity, rtol=2e-15)
    np.testing.assert_array_equal(result.momentum_kg_m_s, initial_momentum)
    np.testing.assert_array_equal(result.kinetic_delta_j, 0.0)
    np.testing.assert_array_equal(result.electric_displacement_work_j, 0.0)
    np.testing.assert_array_equal(case["position_m"], initial_position)
    np.testing.assert_array_equal(case["momentum_kg_m_s"], initial_momentum)


def test_uniform_electric_field_has_exact_nonrelativistic_acceleration_and_work():
    case = inputs(electric_field_v_m=np.array([[2e4, -1e4, 3e4], [1e5, 2e4, -1e4]]))
    dt, charge, mass = case["dt_s"], case["charge_c"], case["mass_kg"]
    initial_p = case["momentum_kg_m_s"]
    force = charge[:, None] * case["electric_field_v_m"]
    expected_p = initial_p + dt * force
    expected_dx = dt * initial_p / mass[:, None] + 0.5 * dt**2 * force / mass[:, None]
    result = advance_orbits(**case, relativistic=False)
    np.testing.assert_allclose(result.momentum_kg_m_s, expected_p, rtol=2e-15)
    np.testing.assert_allclose(result.position_m, expected_dx, rtol=2e-15)
    expected_energy_change = np.sum((expected_p - initial_p) * (expected_p + initial_p), axis=1) / (2.0 * mass)
    np.testing.assert_allclose(result.kinetic_delta_j, expected_energy_change, rtol=1e-13)
    np.testing.assert_allclose(result.electric_displacement_work_j, expected_energy_change, rtol=1e-13)


def test_pure_magnetic_field_rotates_electron_and_ion_with_opposite_signs():
    mass = np.array([ELECTRON_MASS, PROTON_MASS])
    charge = np.array([-ELEMENTARY_CHARGE, ELEMENTARY_CHARGE])
    speed, field, dt = 1e6, 0.7, 1e-12
    momentum = np.column_stack([mass * speed, np.zeros(2), np.zeros(2)])
    result = advance_orbits(
        np.zeros((2, 3)), momentum, mass, charge, dt,
        np.zeros((2, 3)), np.tile([0.0, 0.0, field], (2, 1)), relativistic=False,
    )
    angle = 2.0 * np.arctan(charge * field * dt / (2.0 * mass))
    np.testing.assert_allclose(result.momentum_kg_m_s[:, 0] / (mass * speed), np.cos(angle), rtol=2e-15)
    np.testing.assert_allclose(result.momentum_kg_m_s[:, 1] / (mass * speed), -np.sin(angle), rtol=2e-15)
    assert result.momentum_kg_m_s[0, 1] > 0.0
    assert result.momentum_kg_m_s[1, 1] < 0.0
    np.testing.assert_allclose(result.gyro_angle_rad, np.abs(charge) * field * dt / mass)
    np.testing.assert_array_equal(result.electric_displacement_work_j, 0.0)


@pytest.mark.parametrize("relativistic", [False, True])
def test_many_pure_b_steps_preserve_momentum_norm_and_kinetic_energy(relativistic):
    mass = np.array([ELECTRON_MASS, PROTON_MASS])
    charge = np.array([-ELEMENTARY_CHARGE, ELEMENTARY_CHARGE])
    factor = SPEED_OF_LIGHT_M_S * (10.0 if relativistic else 0.001)
    momentum = mass[:, None] * factor * np.array([[1.0, 0.3, -0.2], [-0.1, 0.7, 1.0]])
    initial_norm = np.linalg.norm(momentum, axis=1)
    initial_energy = kinetic_energy_j(momentum, mass, relativistic=relativistic)
    position = np.zeros((2, 3))
    for _ in range(400):
        step = advance_orbits(position, momentum, mass, charge, 1e-11, np.zeros((2, 3)), np.tile([0.2, -0.1, 1.0], (2, 1)), relativistic=relativistic)
        position, momentum = step.position_m, step.momentum_kg_m_s
        assert np.max(step.magnetic_rotation_relative_error) < 2e-15
    np.testing.assert_allclose(np.linalg.norm(momentum, axis=1), initial_norm, rtol=2e-13)
    np.testing.assert_allclose(kinetic_energy_j(momentum, mass, relativistic=relativistic), initial_energy, rtol=3e-13)
    if relativistic:
        assert np.all(np.linalg.norm(velocity_m_s(momentum, mass), axis=1) < SPEED_OF_LIGHT_M_S)


@pytest.mark.parametrize("charge_sign,mass", [(-1.0, ELECTRON_MASS), (1.0, PROTON_MASS)])
def test_nonrelativistic_exb_drift_averaged_over_integer_boris_turn(charge_sign, mass):
    # A full discrete rotation sums the oscillatory velocity to zero. The
    # remaining displacement must be the charge/mass independent E x B drift.
    count, electric, magnetic = 80, 1000.0, 0.2
    charge = charge_sign * ELEMENTARY_CHARGE
    frequency = abs(charge) * magnetic / mass
    dt = 2.0 * np.tan(np.pi / count) / frequency
    position, momentum = np.zeros((1, 3)), np.zeros((1, 3))
    for _ in range(count):
        step = advance_orbits(position, momentum, np.array([mass]), np.array([charge]), dt,
                              np.array([[electric, 0.0, 0.0]]), np.array([[0.0, 0.0, magnetic]]), relativistic=False)
        position, momentum = step.position_m, step.momentum_kg_m_s
    mean_velocity = position[0] / (count * dt)
    np.testing.assert_allclose(mean_velocity, [0.0, -electric / magnetic, 0.0], atol=2e-10, rtol=2e-13)


def test_relativistic_electric_acceleration_reports_work_quadrature_error():
    mass = np.array([ELECTRON_MASS])
    charge = np.array([-ELEMENTARY_CHARGE])
    electric = np.array([[1e8, 0.0, 0.0]])
    dt = 1e-10
    step = advance_orbits(np.zeros((1, 3)), np.zeros((1, 3)), mass, charge, dt, electric, np.zeros((1, 3)))
    expected_p = charge[:, None] * electric * dt
    expected_gamma = np.sqrt(1.0 + np.sum(expected_p**2, axis=1) / (mass * SPEED_OF_LIGHT_M_S)**2)
    expected_energy = (expected_gamma - 1.0) * mass * SPEED_OF_LIGHT_M_S**2
    np.testing.assert_allclose(step.momentum_kg_m_s, expected_p, rtol=2e-15)
    np.testing.assert_allclose(step.kinetic_delta_j, expected_energy, rtol=2e-15)
    np.testing.assert_allclose(step.work_residual_j, step.kinetic_delta_j - charge * np.sum(electric * step.position_m, axis=1))
    assert abs(step.work_residual_j[0]) > 0.1 * expected_energy[0]
    assert np.linalg.norm(velocity_m_s(step.momentum_kg_m_s, mass)[0]) < SPEED_OF_LIGHT_M_S


def test_relativistic_cold_particle_energy_avoids_gamma_minus_one_cancellation():
    mass = np.array([ELECTRON_MASS])
    speed = 1e-5
    momentum = np.array([[mass[0] * speed, 0.0, 0.0]])
    energy = kinetic_energy_j(momentum, mass)
    assert energy[0] > 0.0
    np.testing.assert_allclose(energy, 0.5 * mass * speed**2, rtol=2e-15)


def test_large_gyro_angle_is_reported_without_clipping_or_energy_loss():
    case = inputs(magnetic_field_t=np.tile([0.0, 0.0, 1.0], (2, 1)), dt_s=1.0)
    step = advance_orbits(**case, relativistic=False)
    assert step.gyro_angle_rad[0] > 1e10
    assert step.boris_rotation_angle_rad[0] < np.pi
    np.testing.assert_allclose(np.linalg.norm(step.momentum_kg_m_s, axis=1), np.linalg.norm(case["momentum_kg_m_s"], axis=1), rtol=3e-15)


@pytest.mark.parametrize("overrides", [
    {"position_m": np.zeros((2, 2))}, {"position_m": np.zeros((1, 3))},
    {"momentum_kg_m_s": np.zeros(3)}, {"momentum_kg_m_s": np.full((2, 3), np.nan)},
    {"mass_kg": ELECTRON_MASS}, {"mass_kg": [ELECTRON_MASS]}, {"mass_kg": [0.0, PROTON_MASS]},
    {"mass_kg": [-ELECTRON_MASS, PROTON_MASS]}, {"mass_kg": [np.inf, PROTON_MASS]},
    {"charge_c": [ELEMENTARY_CHARGE]}, {"charge_c": [np.nan, ELEMENTARY_CHARGE]},
    {"electric_field_v_m": np.zeros(3)}, {"electric_field_v_m": np.full((2, 3), np.inf)},
    {"magnetic_field_t": np.zeros((3, 3))}, {"dt_s": 0.0}, {"dt_s": -1.0},
    {"dt_s": np.inf}, {"dt_s": np.nan}, {"dt_s": [1e-12]}, {"relativistic": "yes"},
])
def test_invalid_or_inconsistent_inputs_are_rejected(overrides):
    with pytest.raises(ValueError):
        advance_orbits(**inputs(**overrides))


def test_empty_particle_batch_is_supported():
    step = advance_orbits(np.empty((0, 3)), np.empty((0, 3)), np.empty(0), np.empty(0), 1e-12, np.empty((0, 3)), np.empty((0, 3)))
    assert step.position_m.shape == (0, 3)
    assert step.kinetic_delta_j.shape == (0,)
