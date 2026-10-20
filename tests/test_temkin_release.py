"""Independent frozen-temperature, coverage-dependent release verification.

All parameters below describe manufactured mathematical cases.  They are not
species or material fits.  Nonlinear oracles integrate a normalized differential
equation rather than reproducing the release implementation's inversion.
"""

from dataclasses import FrozenInstanceError, replace
import math

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from mhdlab.constants import KB, QE
from mhdlab.temkin_release import (
    LAW_IDENTITY, TemkinReservoir, TemkinReleaseTrial, UnsupportedTemkinRangeError,
    release_step,
)
from mhdlab.thermal_release import (
    BindingReservoir, SurfaceInventory, ThermalReleaseTrial,
    energy_budgets as baseline_energy_budgets, release_step as baseline_release_step,
)


KB_EV_K = KB / QE
DIAGNOSTICS = (
    "released_m2", "average_flux_m2_s", "translational_energy_j_m2",
    "internal_energy_j_m2", "binding_energy_j_m2", "total_wall_energy_debit_j_m2",
    "inventory_balance_residual_m2",
)


def reservoir(**overrides):
    parameters = dict(
        name="manufactured-pool", state_id="manufactured:ground",
        attempt_frequency_s=7.5e3, activation_energy_ev=.18,
        desorption_energy_ev=.23, internal_energy_ev=.17,
        coverage_coefficient=.65, reference_coverage_m2=8e17,
        max_coverage_m2=8e17,
        kinetic_provenance="manufactured coverage-dependent rate for verification",
        energy_provenance="manufactured independent constant release-energy cost",
    )
    parameters.update(overrides)
    return TemkinReservoir(**parameters)


def constant_barrier_pool(pool):
    return BindingReservoir(
        pool.name, pool.state_id, pool.attempt_frequency_s, pool.activation_energy_ev,
        pool.desorption_energy_ev, pool.internal_energy_ev,
    )


def independent_survivor_fraction(pool, temperature, initial, dt):
    """DOP853 oracle for v=-log(N/N0), with a dimensionless clock.

    v' = exp[c*(exp(-v)-1)] in tau=initial_rate*t is independent of
    a special-function inverse and avoids subtracting nearly equal coverages.
    """
    if initial == 0:
        return 0., 0.
    activation = pool.activation_energy_ev / (KB_EV_K * temperature)
    coefficient = activation * pool.coverage_coefficient * initial / pool.reference_coverage_m2
    initial_rate = pool.attempt_frequency_s * math.exp(-activation + coefficient)
    duration = initial_rate * dt
    solution = solve_ivp(
        lambda _, value: [math.exp(coefficient * math.expm1(-float(value[0])))],
        (0., duration), [0.], method="DOP853", rtol=2e-13, atol=2e-15,
    )
    assert solution.success
    decay = float(solution.y[0, -1])
    return math.exp(-decay), -math.expm1(-decay)


@pytest.mark.parametrize("temperature", [250., 700., 1400.])
@pytest.mark.parametrize("duration", [1e-8, 3e-4, .03])
def test_zero_coverage_coefficient_recovers_analytic_first_order_law(temperature, duration):
    pool = reservoir(coverage_coefficient=0.)
    initial = np.array([[0., 2e16, 6e17]])
    result = release_step(SurfaceInventory(initial), [pool], [temperature] * 3, duration)
    hazard = pool.attempt_frequency_s * duration * math.exp(-pool.activation_energy_ev / (KB_EV_K * temperature))
    np.testing.assert_allclose(result.trial_inventory.coverage_m2, initial * math.exp(-hazard), rtol=2e-12, atol=0.)
    np.testing.assert_allclose(result.released_m2, initial * (-math.expm1(-hazard)), rtol=2e-12, atol=0.)
    baseline = baseline_release_step(SurfaceInventory(initial), [constant_barrier_pool(pool)],
                                     [temperature] * 3, duration)
    for field in DIAGNOSTICS[:-1]:
        np.testing.assert_allclose(getattr(result, field), getattr(baseline, field), rtol=2e-12, atol=0.)


def test_zero_activation_energy_is_constant_first_order_even_for_large_alpha():
    # No rate barrier remains when Ed=0; the unused dimensionless coefficient
    # need not satisfy the positive-Ed barrier domain inequality.
    pool = reservoir(activation_energy_ev=0., coverage_coefficient=1e5,
                     attempt_frequency_s=2., reference_coverage_m2=1., max_coverage_m2=4.)
    initial = np.array([[1., 3.]])
    result = release_step(SurfaceInventory(initial), [pool], [1., 1e9], .3)
    np.testing.assert_allclose(result.trial_inventory.coverage_m2, initial * math.exp(-.6), rtol=2e-12, atol=0.)
    np.testing.assert_allclose(result.released_m2, initial * (-math.expm1(-.6)), rtol=2e-12, atol=0.)


@pytest.mark.parametrize("coverage_fraction", [.04, .4, .95])
@pytest.mark.parametrize("initial_clock", [.001, .5, 5., 30.])
def test_nonlinear_release_matches_independent_log_survivor_ode(coverage_fraction, initial_clock):
    pool = reservoir()
    temperature = 700.
    initial = coverage_fraction * pool.reference_coverage_m2
    barrier = pool.activation_energy_ev / (KB_EV_K * temperature)
    rate = pool.attempt_frequency_s * math.exp(-barrier * (1 - pool.coverage_coefficient * coverage_fraction))
    duration = initial_clock / rate
    expected_after, expected_release = independent_survivor_fraction(pool, temperature, initial, duration)
    result = release_step(SurfaceInventory([[initial]]), [pool], [temperature], duration)
    assert result.trial_inventory.coverage_m2[0, 0] / initial == pytest.approx(expected_after, rel=4e-11, abs=0.)
    assert result.released_m2[0, 0] / initial == pytest.approx(expected_release, rel=4e-11, abs=0.)


@pytest.mark.parametrize("coefficient", [3., 20., 200., 512.])
def test_large_supported_coverage_feedback_matches_independent_ode(coefficient):
    # For the upper endpoint, use an exact binary scaling of the supplied
    # constants: Ed*QE/(KB*T)=1024 exactly for these double inputs.  Forming
    # Ed=1024*(KB/QE) instead can round just outside the declared domain.
    temperature = QE if coefficient == 512 else 1.
    activation = 1024 * KB if coefficient == 512 else 2 * coefficient * KB_EV_K
    pool = reservoir(attempt_frequency_s=1e5, activation_energy_ev=activation,
                     coverage_coefficient=.5, reference_coverage_m2=1e18, max_coverage_m2=1e18)
    initial = 1e18
    duration = 20. / (pool.attempt_frequency_s * math.exp(-coefficient))
    expected_after, expected_release = independent_survivor_fraction(pool, temperature, initial, duration)
    result = release_step(SurfaceInventory([[initial]]), [pool], [temperature], duration)
    assert result.trial_inventory.coverage_m2[0, 0] / initial == pytest.approx(expected_after, rel=5e-11, abs=0.)
    assert result.released_m2[0, 0] / initial == pytest.approx(expected_release, rel=5e-11, abs=0.)


def test_unsupported_nonlinear_range_is_explicit_and_does_not_freeze_the_rate():
    pool = reservoir(attempt_frequency_s=1., activation_energy_ev=513 * KB_EV_K,
                     coverage_coefficient=1., reference_coverage_m2=1., max_coverage_m2=1.)
    with pytest.raises(UnsupportedTemkinRangeError):
        release_step(SurfaceInventory([[1.]]), [pool], [1.], 1.)


def test_multiple_pools_and_temperatures_remain_local():
    pools = [reservoir(), reservoir(name="manufactured-second", state_id="manufactured:excited",
                                    coverage_coefficient=.2, activation_energy_ev=.11,
                                    reference_coverage_m2=2e17, max_coverage_m2=3e17)]
    temperatures = np.array([350., 600., 1100.])
    initial = np.array([[7e17, 2e17, 0.], [3e16, 1e17, 2e17]])
    result = release_step(SurfaceInventory(initial), pools, temperatures, .004)
    assert result.reservoir_names == tuple(pool.name for pool in pools)
    assert result.state_ids == tuple(pool.state_id for pool in pools)
    for row, pool in enumerate(pools):
        for column, temperature in enumerate(temperatures):
            fraction_after, fraction_released = independent_survivor_fraction(pool, temperature, initial[row, column], .004)
            assert result.trial_inventory.coverage_m2[row, column] == pytest.approx(initial[row, column] * fraction_after, rel=4e-11, abs=0.)
            assert result.released_m2[row, column] == pytest.approx(initial[row, column] * fraction_released, rel=4e-11, abs=0.)


def test_fixed_temperature_evolution_obeys_semigroup_and_additive_ledgers():
    pool = reservoir()
    initial = SurfaceInventory([[.1 * pool.reference_coverage_m2, .9 * pool.reference_coverage_m2]])
    temperature = [500., 1000.]
    whole = release_step(initial, [pool], temperature, .008)
    first = release_step(initial, [pool], temperature, .003)
    second = release_step(first.trial_inventory, [pool], temperature, .005)
    np.testing.assert_allclose(second.trial_inventory.coverage_m2, whole.trial_inventory.coverage_m2, rtol=2e-11, atol=0.)
    for field in DIAGNOSTICS:
        if field in ("average_flux_m2_s", "inventory_balance_residual_m2"):
            continue
        np.testing.assert_allclose(getattr(first, field) + getattr(second, field), getattr(whole, field), rtol=2e-11, atol=0.)


def test_inversion_time_residual_matches_independent_coverage_integral():
    pool = reservoir()
    initial, temperature, duration = 6e17, 650., .003
    result = release_step(SurfaceInventory([[initial]]), [pool], [temperature], duration)
    fraction = result.trial_inventory.coverage_m2[0, 0] / initial
    activation = pool.activation_energy_ev / (KB_EV_K * temperature)
    coefficient = activation * pool.coverage_coefficient * initial / pool.reference_coverage_m2
    clock = pool.attempt_frequency_s * math.exp(-activation + coefficient) * duration
    integral, error = quad(lambda value: math.exp(coefficient * (1 - value)) / value,
                           fraction, 1., epsabs=1e-13, epsrel=1e-13)
    assert error < 1e-11 * clock
    independent_residual = math.log(integral / clock)
    assert abs(independent_residual) < 2e-11
    assert result.inversion_log_time_residual[0, 0] == pytest.approx(independent_residual, rel=0., abs=2e-13)
    assert result.inversion_iterations[0, 0] >= 0


def test_reference_coverage_rescaling_preserves_the_same_barrier_function():
    first_pool = reservoir(coverage_coefficient=.25, reference_coverage_m2=2e18,
                           max_coverage_m2=1e18)
    second_pool = replace(first_pool, coverage_coefficient=.5, reference_coverage_m2=4e18)
    inventory = SurfaceInventory([[2e17, 7e17]])
    first = release_step(inventory, [first_pool], [600., 1000.], .007)
    second = release_step(inventory, [second_pool], [600., 1000.], .007)
    np.testing.assert_allclose(first.released_m2, second.released_m2, rtol=2e-12, atol=0.)
    np.testing.assert_allclose(first.trial_inventory.coverage_m2, second.trial_inventory.coverage_m2, rtol=2e-12, atol=0.)


def test_rescaling_all_areal_densities_preserves_fractional_kinetics():
    pool = reservoir()
    scale = 1e-11
    scaled = replace(pool, reference_coverage_m2=scale * pool.reference_coverage_m2,
                     max_coverage_m2=scale * pool.max_coverage_m2)
    initial = np.array([[.3, .9]]) * pool.reference_coverage_m2
    first = release_step(SurfaceInventory(initial), [pool], [600., 900.], .004)
    second = release_step(SurfaceInventory(scale * initial), [scaled], [600., 900.], .004)
    np.testing.assert_allclose(second.released_m2, scale * first.released_m2, rtol=3e-12, atol=0.)
    np.testing.assert_allclose(second.trial_inventory.coverage_m2, scale * first.trial_inventory.coverage_m2, rtol=3e-12, atol=0.)


def test_reference_density_is_not_replaced_by_each_initial_inventory():
    pool = reservoir()
    fractions = np.array([.2, .8])
    initial = fractions * pool.reference_coverage_m2
    temperature = 700.
    activation = pool.activation_energy_ev / (KB_EV_K * temperature)
    initial_rates = pool.attempt_frequency_s * np.exp(-activation * (1 - pool.coverage_coefficient * fractions))
    duration = 1e-8 / max(initial_rates)
    result = release_step(SurfaceInventory([initial]), [pool], [temperature] * 2, duration)
    fractional_rates = result.released_m2[0] / initial / duration
    np.testing.assert_allclose(fractional_rates, initial_rates, rtol=4e-8, atol=0.)
    assert fractional_rates[1] / fractional_rates[0] == pytest.approx(
        math.exp(activation * pool.coverage_coefficient * (.8 - .2)), rel=4e-8, abs=0.)


def test_full_declared_domain_endpoint_can_have_zero_barrier_without_clipping():
    pool = reservoir(coverage_coefficient=.5, reference_coverage_m2=2.,
                     max_coverage_m2=4., attempt_frequency_s=3.)
    duration = 1e-9
    result = release_step(SurfaceInventory([[4.]]), [pool], [500.], duration)
    assert result.released_m2[0, 0] / (4 * duration) == pytest.approx(3., rel=5e-8, abs=0.)


def test_alpha_above_one_is_valid_on_an_explicit_smaller_coverage_domain():
    pool = reservoir(coverage_coefficient=3., reference_coverage_m2=1e18,
                     max_coverage_m2=2e17)
    initial = 1.5e17
    expected_after, expected_release = independent_survivor_fraction(pool, 700., initial, .002)
    result = release_step(SurfaceInventory([[initial]]), [pool], [700.], .002)
    assert result.trial_inventory.coverage_m2[0, 0] / initial == pytest.approx(expected_after, rel=4e-11, abs=0.)
    assert result.released_m2[0, 0] / initial == pytest.approx(expected_release, rel=4e-11, abs=0.)


def test_energy_cost_uses_supplied_reference_cost_and_flux_thermal_mean_once():
    pool = reservoir(activation_energy_ev=.9, coverage_coefficient=.7,
                     attempt_frequency_s=1e8)
    temperature = 400.
    result = release_step(SurfaceInventory([[.9 * pool.reference_coverage_m2]]), [pool], [temperature], 1e-5)
    amount = result.released_m2[0, 0]
    assert amount > 0
    per_particle = {
        "translational_energy_j_m2": 2 * KB * temperature,
        "binding_energy_j_m2": QE * pool.desorption_energy_ev,
        "internal_energy_j_m2": QE * pool.internal_energy_ev,
        "total_wall_energy_debit_j_m2": 2 * KB * temperature + QE * (pool.desorption_energy_ev + pool.internal_energy_ev),
    }
    for field, expected in per_particle.items():
        assert getattr(result, field)[0, 0] / amount == pytest.approx(expected, rel=3e-13, abs=0.)
    assert per_particle["binding_energy_j_m2"] != QE * pool.activation_energy_ev


def test_tiny_nonlinear_release_survives_underflow_of_the_hazard_itself():
    pool = reservoir(attempt_frequency_s=1., activation_energy_ev=1000 * KB_EV_K,
                     coverage_coefficient=.001, reference_coverage_m2=1e300,
                     max_coverage_m2=1e300, desorption_energy_ev=0., internal_energy_ev=0.)
    result = release_step(SurfaceInventory([[1e300]]), [pool], [1.], 1.)
    expected = math.exp(math.log(1e300) - 999.)
    assert result.released_m2[0, 0] > 0
    assert result.released_m2[0, 0] == pytest.approx(expected, rel=3e-11, abs=0.)
    assert result.trial_inventory.coverage_m2[0, 0] == 1e300
    assert result.inventory_balance_residual_m2[0, 0] == result.released_m2[0, 0]


def test_tiny_nonlinear_survivor_survives_underflow_of_fraction_itself():
    pool = reservoir(attempt_frequency_s=1., activation_energy_ev=KB_EV_K,
                     coverage_coefficient=.5, reference_coverage_m2=1e300,
                     max_coverage_m2=1e300, desorption_energy_ev=0., internal_energy_ev=0.)
    duration = 800 * math.e
    # For large v, v=nu*exp(-A)*dt+Ein(c), c=A*alpha*N0/Nref.
    # The neglected correction is O(c*exp(-v)), far below float64 precision.
    correction, quadrature_error = quad(lambda value: -math.expm1(-value) / value if value else 1., 0., .5,
                                        epsabs=1e-14, epsrel=1e-14)
    assert quadrature_error < 1e-13
    expected = math.exp(math.log(1e300) - 800 - correction)
    result = release_step(SurfaceInventory([[1e300]]), [pool], [1.], duration)
    assert result.trial_inventory.coverage_m2[0, 0] > 0
    assert result.trial_inventory.coverage_m2[0, 0] == pytest.approx(expected, rel=3e-11, abs=0.)
    assert result.released_m2[0, 0] == 1e300
    assert result.inventory_balance_residual_m2[0, 0] == result.trial_inventory.coverage_m2[0, 0]


def test_extreme_stiffness_depletes_positively_and_records_finite_ledgers():
    pool = reservoir(attempt_frequency_s=1e300, activation_energy_ev=KB_EV_K,
                     coverage_coefficient=.8, reference_coverage_m2=1e20,
                     max_coverage_m2=1e20)
    result = release_step(SurfaceInventory([[0., 1e20]]), [pool], [1., 1.], 1e300)
    np.testing.assert_array_equal(result.trial_inventory.coverage_m2, [[0., 0.]])
    np.testing.assert_array_equal(result.released_m2, [[0., 1e20]])
    for field in DIAGNOSTICS:
        assert np.all(np.isfinite(getattr(result, field)))
    assert result.average_flux_m2_s[0, 1] == pytest.approx(1e-280, rel=2e-13, abs=0.)


def test_trial_and_inputs_are_immutable_and_law_identity_is_explicit():
    pool = reservoir()
    original = np.array([[2e17, 6e17]])
    temperatures = np.array([500., 900.])
    inventory = SurfaceInventory(original)
    result = release_step(inventory, [pool], temperatures, .001)
    assert isinstance(result, TemkinReleaseTrial)
    assert not isinstance(result, ThermalReleaseTrial)
    assert result.law_identity == LAW_IDENTITY == "temkin_first_order_constant_attempt_frequency_v1"
    original[0, 0], temperatures[0] = 0., 1.
    np.testing.assert_array_equal(inventory.coverage_m2, [[2e17, 6e17]])
    np.testing.assert_array_equal(result.temperature_k, [500., 900.])
    for field in DIAGNOSTICS + ("temperature_k",):
        with pytest.raises(ValueError):
            getattr(result, field).setflags(write=True)
    for field in ("inversion_iterations", "inversion_regime", "inversion_log_time_residual"):
        assert getattr(result, field).shape == inventory.coverage_m2.shape
        assert np.all(np.isfinite(getattr(result, field)))
        with pytest.raises(ValueError):
            getattr(result, field).setflags(write=True)
    with pytest.raises(ValueError):
        result.trial_inventory.coverage_m2.setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        pool.coverage_coefficient = .1
    with pytest.raises(FrozenInstanceError):
        result.dt_s = 1.


@pytest.mark.parametrize("overrides", [
    {"name": ""}, {"state_id": " "}, {"kinetic_provenance": ""}, {"energy_provenance": " "},
    {"kinetic_provenance": True}, {"energy_provenance": None},
    {"attempt_frequency_s": 0.}, {"attempt_frequency_s": -1.},
    {"activation_energy_ev": -1.}, {"desorption_energy_ev": -1.}, {"internal_energy_ev": -1.},
    {"coverage_coefficient": -1.}, {"reference_coverage_m2": 0.}, {"max_coverage_m2": 0.},
    {"coverage_coefficient": True}, {"reference_coverage_m2": "8e17"},
    {"max_coverage_m2": 8e17 + 0j}, {"attempt_frequency_s": "7.5e3"},
    {"activation_energy_ev": False}, {"desorption_energy_ev": .2 + 0j},
    {"internal_energy_ev": "0"}, {"coverage_coefficient": np.nan},
    {"reference_coverage_m2": np.inf}, {"max_coverage_m2": np.nan},
    {"coverage_coefficient": 1.01},
    {"coverage_coefficient": np.nextafter(1., math.inf)},
    {"coverage_coefficient": .6, "reference_coverage_m2": 1., "max_coverage_m2": 2.},
])
def test_reservoir_rejects_ambiguous_parameters_missing_sources_and_negative_barrier_domain(overrides):
    with pytest.raises(ValueError):
        reservoir(**overrides)


@pytest.mark.parametrize("temperature,duration", [
    ([0.], .1), ([-1.], .1), ([np.inf], .1), ([np.nan], .1), ([True], .1),
    (["500"], .1), ([500. + 0j], .1), ([[500.]], .1), ([500., 600.], .1),
    ([500.], 0.), ([500.], -1.), ([500.], True), ([500.], ".1"),
    ([500.], .1 + 0j), ([500.], np.inf), ([500.], np.nan),
    (np.ma.array([500.], mask=[True]), .1),
])
def test_release_step_rejects_invalid_temperature_and_time(temperature, duration):
    with pytest.raises(ValueError):
        release_step(SurfaceInventory([[1e17]]), [reservoir()], temperature, duration)


def test_release_step_enforces_reservoir_count_unique_names_and_inventory_domain():
    pool = reservoir()
    for pools in ([], [pool, pool], [constant_barrier_pool(pool)], "pool"):
        with pytest.raises(ValueError):
            release_step(SurfaceInventory([[1e17]]), pools, [500.], .1)
    with pytest.raises(ValueError):
        release_step(SurfaceInventory([[1e17], [1e17]]), [pool, pool], [500.], .1)
    with pytest.raises(ValueError):
        release_step(SurfaceInventory([[np.nextafter(pool.max_coverage_m2, math.inf)]]), [pool], [500.], .1)
    with pytest.raises(ValueError):
        release_step([[1e17]], [pool], [500.], .1)
    assert not isinstance(pool, BindingReservoir)


def test_baseline_release_energy_and_history_interfaces_refuse_the_distinct_law():
    from mhdlab.flash_surface import FlashSurfaceHistory
    from mhdlab.surface_release_history import evolve_surface_release
    from mhdlab.thermal_boundary import sample_emission_packets

    pool = reservoir()
    inventory = SurfaceInventory([[1e17]])
    with pytest.raises(ValueError):
        baseline_release_step(inventory, [pool], [500.], .001)
    with pytest.raises(ValueError):
        baseline_energy_budgets([[1e15]], [pool], [500.])
    surface = FlashSurfaceHistory(
        [0., .001], [0., .1], [[500.], [500.]], "synthetic", "manufactured",
        "manufactured temperature for unsupported-law rejection", "single manufactured face",
    )
    with pytest.raises(ValueError):
        evolve_surface_release(surface, inventory, [pool], parameter_provenance="manufactured unsupported-law test")
    trial = release_step(inventory, [pool], [500.], .001)
    with pytest.raises(ValueError):
        sample_emission_packets(trial, {pool.state_id: 1.}, 0, np.full((1, 6), .5),
                                outward_normal_sign=1, current_length_m=.1)
