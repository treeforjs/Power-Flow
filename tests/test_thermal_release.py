"""Analytic thermal-release checks using supplied reference parameters."""

from dataclasses import FrozenInstanceError, replace
import math

import numpy as np
import pytest

from mhdlab.constants import KB, QE
from mhdlab.thermal_release import BindingReservoir, SurfaceInventory, energy_budgets, release_step


def reservoir(**overrides):
    values = dict(name="reference_pool", state_id="reference:gas_state", attempt_frequency_s=2.0,
                  activation_energy_ev=0.0, desorption_energy_ev=0.7, internal_energy_ev=0.2)
    values.update(overrides)
    return BindingReservoir(**values)


def test_exact_analytic_decay_and_nonmutating_immutable_trial():
    source = np.array([[10.0, 20.0]])
    inventory = SurfaceInventory(source)
    source[0, 0] = 999
    trial = release_step(inventory, [reservoir()], [300.0, 600.0], 0.25)
    expected_remaining = np.array([[10, 20]]) * np.exp(-0.5)
    expected_release = np.array([[10, 20]]) * (-np.expm1(-0.5))
    np.testing.assert_allclose(trial.trial_inventory.coverage_m2, expected_remaining, rtol=2e-15)
    np.testing.assert_allclose(trial.released_m2, expected_release, rtol=2e-15)
    np.testing.assert_array_equal(inventory.coverage_m2, [[10, 20]])
    np.testing.assert_allclose(trial.average_flux_m2_s, expected_release / 0.25, rtol=2e-15)
    assert trial.reservoir_names == ("reference_pool",)
    assert trial.state_ids == ("reference:gas_state",)
    assert np.max(np.abs(trial.inventory_balance_residual_m2)) < 4e-15
    for array in (inventory.coverage_m2, trial.released_m2, trial.temperature_k):
        with pytest.raises(ValueError):
            array.flags.writeable = True
    with pytest.raises(FrozenInstanceError):
        inventory.coverage_m2 = source


def test_stiff_exhaustion_avoids_attempt_times_dt_overflow():
    pool = reservoir(attempt_frequency_s=1e300)
    trial = release_step(SurfaceInventory([[1e20, 0.0]]), [pool], [300.0, 300.0], 1e100)
    np.testing.assert_array_equal(trial.released_m2, [[1e20, 0]])
    np.testing.assert_array_equal(trial.trial_inventory.coverage_m2, [[0, 0]])
    assert trial.average_flux_m2_s[0, 0] == pytest.approx(1e-80)


def test_energy_budgets_use_2kbt_and_distinguish_barrier_from_binding():
    pool = reservoir(activation_energy_ev=0.05, desorption_energy_ev=0.7, internal_energy_ev=0.2)
    temperatures = np.array([400.0, 800.0])
    trial = release_step(SurfaceInventory([[1e20, 2e20]]), [pool], temperatures, 0.2)
    np.testing.assert_allclose(trial.translational_energy_j_m2 / trial.released_m2,
                               2 * KB * temperatures[None, :], rtol=2e-15)
    np.testing.assert_allclose(trial.internal_energy_j_m2 / trial.released_m2, 0.2 * QE, rtol=2e-15)
    np.testing.assert_allclose(trial.binding_energy_j_m2 / trial.released_m2, 0.7 * QE, rtol=2e-15)
    np.testing.assert_allclose(trial.total_wall_energy_debit_j_m2,
                               trial.translational_energy_j_m2 + trial.internal_energy_j_m2
                               + trial.binding_energy_j_m2, rtol=2e-15)
    changed_binding = release_step(SurfaceInventory([[1e20, 2e20]]),
                                   [replace(pool, desorption_energy_ev=1.4)], temperatures, 0.2)
    np.testing.assert_array_equal(changed_binding.released_m2, trial.released_m2)
    np.testing.assert_allclose(changed_binding.binding_energy_j_m2, 2 * trial.binding_energy_j_m2)


def test_subdivision_invariance_at_fixed_temperature():
    pool = reservoir(attempt_frequency_s=1e9, activation_energy_ev=0.4)
    initial = SurfaceInventory([[4e20, 8e20]])
    whole = release_step(initial, [pool], [400.0, 500.0], 1e-5)
    first = release_step(initial, [pool], [400.0, 500.0], 3e-6)
    second = release_step(first.trial_inventory, [pool], [400.0, 500.0], 7e-6)
    np.testing.assert_allclose(second.trial_inventory.coverage_m2, whole.trial_inventory.coverage_m2, rtol=4e-14)
    np.testing.assert_allclose(first.released_m2 + second.released_m2, whole.released_m2, rtol=4e-14)
    np.testing.assert_allclose(first.total_wall_energy_debit_j_m2 + second.total_wall_energy_debit_j_m2,
                               whole.total_wall_energy_debit_j_m2, rtol=4e-14)


def test_independent_pools_and_width_bins_can_emit_same_state():
    pools = [reservoir(name="fast", attempt_frequency_s=4),
             reservoir(name="slow", attempt_frequency_s=1)]
    trial = release_step(SurfaceInventory([[3, 6], [5, 10]]), pools, [300, 900], 0.5)
    expected = np.array([[3, 6], [5, 10]]) * np.array([math.exp(-2), math.exp(-0.5)])[:, None]
    np.testing.assert_allclose(trial.trial_inventory.coverage_m2, expected, rtol=2e-15)
    assert trial.state_ids[0] == trial.state_ids[1]


def test_arrhenius_hazard_temperature_ratio_without_depletion_linearization():
    pool = reservoir(attempt_frequency_s=1e6, activation_energy_ev=0.2)
    temperatures = np.array([300.0, 600.0])
    trial = release_step(SurfaceInventory([[1.0, 1.0]]), [pool], temperatures, 1e-6)
    hazards = -np.log(trial.trial_inventory.coverage_m2[0])
    expected_ratio = math.exp(pool.activation_energy_ev * QE / KB * (1 / 300 - 1 / 600))
    assert hazards[1] / hazards[0] == pytest.approx(expected_ratio, rel=3e-12)


def test_no_premature_tiny_hazard_cutoff_and_balance_residual_is_explicit():
    pool = reservoir(attempt_frequency_s=1e-200)
    initial = SurfaceInventory([[1e300]])
    trial = release_step(initial, [pool], [300.0], 1e-200)
    # h=1e-400 cannot be represented alone, but N*h=1e-100 can.
    assert trial.released_m2[0, 0] == pytest.approx(1e-100, rel=2e-13, abs=0)
    assert trial.average_flux_m2_s[0, 0] == pytest.approx(1e100, rel=2e-13)
    assert trial.trial_inventory.coverage_m2[0, 0] == initial.coverage_m2[0, 0]
    assert trial.inventory_balance_residual_m2[0, 0] == trial.released_m2[0, 0]


def test_representable_survivors_remain_when_survival_factor_underflows():
    pool = reservoir(attempt_frequency_s=1000.0)
    trial = release_step(SurfaceInventory([[1e300]]), [pool], [300.0], 1.0)
    assert trial.trial_inventory.coverage_m2[0, 0] == pytest.approx(math.exp(math.log(1e300) - 1000), rel=2e-12, abs=0)
    assert trial.trial_inventory.coverage_m2[0, 0] > 0


def test_large_barrier_and_subnormal_temperature_have_no_artificial_temperature_floor():
    pool = reservoir(activation_energy_ev=1e300)
    trial = release_step(SurfaceInventory([[1.0]]), [pool], [1e-300], 1.0)
    np.testing.assert_array_equal(trial.released_m2, [[0]])
    np.testing.assert_array_equal(trial.trial_inventory.coverage_m2, [[1]])


@pytest.mark.parametrize("overrides", [
    {"name": ""}, {"name": " leading"}, {"state_id": " "},
    {"attempt_frequency_s": 0}, {"attempt_frequency_s": -1},
    {"activation_energy_ev": -1}, {"desorption_energy_ev": -1}, {"internal_energy_ev": -1},
    {"attempt_frequency_s": True}, {"activation_energy_ev": 1j},
    {"desorption_energy_ev": "0.5"}, {"internal_energy_ev": np.nan},
])
def test_reservoir_validation(overrides):
    with pytest.raises(ValueError):
        reservoir(**overrides)


@pytest.mark.parametrize("coverage", [[], [1, 2], [[]], [[-1]], [[np.inf]], [[1j]], [[True]],
                                       [[1.0, True]], [["1.0"]], [[1], [1, 2]]])
def test_inventory_rejects_invalid_or_coerced_numeric_inputs(coverage):
    with pytest.raises(ValueError):
        SurfaceInventory(coverage)


@pytest.mark.parametrize("temperature,dt", [
    ([0], 1), ([-1], 1), ([np.nan], 1), ([True], 1), ([300, True], 1),
    ([300 + 0j], 1), (["300"], 1), ([[300]], 1), ([300, 400], 1),
    ([300], 0), ([300], -1), ([300], True), ([300], 1j), ([300], np.inf), ([300], "1"),
])
def test_step_rejects_invalid_temperature_and_time(temperature, dt):
    with pytest.raises(ValueError):
        release_step(SurfaceInventory([[1.0]]), [reservoir()], temperature, dt)


def test_reservoir_order_and_count_are_explicit_and_unique():
    with pytest.raises(ValueError, match="count"):
        release_step(SurfaceInventory([[1], [2]]), [reservoir()], [300], 1)
    with pytest.raises(ValueError, match="unique"):
        release_step(SurfaceInventory([[1], [2]]), [reservoir(), reservoir()], [300], 1)
    for pools in ([], [object()], "pool"):
        with pytest.raises(ValueError, match="reservoir"):
            release_step(SurfaceInventory([[1]]), pools, [300], 1)


def test_unrepresentable_flux_or_energy_reports_overflow_instead_of_inf():
    with pytest.raises(ValueError, match="flux"):
        release_step(SurfaceInventory([[1e300]]), [reservoir(attempt_frequency_s=1e300)], [300], 1e-200)
    with pytest.raises(ValueError, match="float64"):
        release_step(SurfaceInventory([[1e300]]), [reservoir()], [1e100], 1)


def test_public_energy_budgets_match_step_and_preserve_scaled_extremes():
    pools = [reservoir()]
    trial = release_step(SurfaceInventory([[1e20, 2e20]]), pools, [300, 600], 0.2)
    budgets = energy_budgets(trial.released_m2, pools, [300, 600])
    assert set(budgets) == {"translational_energy_j_m2", "internal_energy_j_m2",
                            "binding_energy_j_m2", "total_wall_energy_debit_j_m2"}
    for name, values in budgets.items():
        np.testing.assert_array_equal(values, getattr(trial, name))
        with pytest.raises(ValueError):
            values.flags.writeable = True
    large = energy_budgets([[1e300]], [reservoir(desorption_energy_ev=1e20)], [300])
    assert large["binding_energy_j_m2"][0, 0] == pytest.approx(1e300 * (1e20 * QE), rel=2e-15)
    cold = energy_budgets([[1e300]], [reservoir(desorption_energy_ev=0, internal_energy_ev=0)], [1e-300])
    assert cold["translational_energy_j_m2"][0, 0] == pytest.approx(2 * KB, rel=2e-15)


@pytest.mark.parametrize("release,temperature", [
    ([[True]], [300]), ([[-1]], [300]), ([[1]], [True]), ([[1]], [0]),
    ([[1j]], [300]), ([[1, 2]], [300]), ([[1]], [[300]]),
])
def test_public_energy_budgets_validate_inputs(release, temperature):
    with pytest.raises(ValueError):
        energy_budgets(release, [reservoir()], temperature)
