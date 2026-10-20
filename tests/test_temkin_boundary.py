"""Temkin conditional times, thermal moments and conservative packet pairing."""

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from mhdlab.constants import KB, QE
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.surface_release_history import evolve_surface_release
from mhdlab.temkin_boundary import TemkinEmissionPackets, sample_temkin_emission_packets
from mhdlab.temkin_emission_time import NoTemkinEmissionError, TIME_IDENTITY
from mhdlab.temkin_history import evolve_temkin_release
from mhdlab.temkin_release import LAW_IDENTITY, TemkinReservoir
from mhdlab.thermal_boundary import ThermalEmissionPackets
from mhdlab.thermal_release import BindingReservoir, SurfaceInventory


def history(*, temperatures=(450, 900), coverage=((1e18, 5e17),), coefficient=.4,
            barrier=.3, rate=1e8, times=(0, 1e-7), edges=(0, .001, .004)):
    surface = FlashSurfaceHistory(times, edges, [temperatures] * 2, "synthetic", "temkin-packets",
                                   "manufactured prescribed temperatures", "single test face")
    pools = [TemkinReservoir(f"pool-{i}", f"M:{i}", rate, barrier, .5 + .1 * i, .2 + .1 * i,
                             coefficient, 1e18, 1e18, "manufactured kinetic law", "manufactured energy cost")
             for i in range(len(coverage))]
    return evolve_temkin_release(surface, SurfaceInventory(coverage), pools,
                                 parameter_provenance="manufactured Temkin packet verification")


def draws(n):
    return np.full((n, 6), .5)


def sample(source=None, uniforms=None, **kwargs):
    source = history() if source is None else source
    parameters = dict(outward_normal_sign=1, current_length_m=.025)
    parameters.update(kwargs)
    masses = {pool.state_id: (index + 1) * 3e-26 for index, pool in enumerate(source.reservoirs)}
    return sample_temkin_emission_packets(source, masses, 0, draws(8) if uniforms is None else uniforms, **parameters)


def test_every_component_inventory_state_temperature_and_energy_ledger_is_retained():
    source = history(coverage=((1e18, 5e17), (2e17, 8e17)))
    uniforms = draws(13)
    uniforms[:, 0] = np.linspace(.1, .9, len(uniforms))
    uniforms[:, 2] = np.linspace(.12, .88, len(uniforms))
    result = sample(source, uniforms)
    assert isinstance(result, TemkinEmissionPackets)
    assert not isinstance(result, ThermalEmissionPackets)
    assert result.law_identity == LAW_IDENTITY
    assert result.time_model == result.time_identity == TIME_IDENTITY
    assert result.required_minimum_n == 4
    for component, (pool, width) in enumerate(zip(result.component_reservoir_index, result.component_width_index)):
        selected = (result.reservoir_index == pool) & (result.source_width_index == width)
        expected = source.released_m2[0, pool, width] * np.diff(source.x_edges_m)[width] * .025
        assert np.count_nonzero(selected) == result.packet_counts[component]
        assert math.fsum(result.particle_weight[selected]) == pytest.approx(expected, rel=3e-15)
        assert result.component_released_particles[component] == pytest.approx(expected, rel=3e-15)
        expected_energy = 2 * KB * source.temperature_k[0, width] * expected
        assert result.component_expected_translational_energy_j[component] == pytest.approx(expected_energy, rel=5e-15)
        sampled = math.fsum(.5 * result.mass_kg[i] * np.dot(result.velocity_m_s[i], result.velocity_m_s[i])
                            * result.particle_weight[i] for i in np.flatnonzero(selected))
        assert result.component_sampled_translational_energy_j[component] == pytest.approx(sampled, rel=5e-15)
        assert result.component_translational_energy_sampling_residual_j[component] == pytest.approx(sampled - expected_energy, rel=1e-13)
        np.testing.assert_array_equal(result.temperature_k[selected], source.temperature_k[0, width])
        assert all(state == source.reservoirs[pool].state_id for state, use in zip(result.state_ids, selected) if use)
        assert np.all(result.resolved_x_m[selected] > source.x_edges_m[width])
        assert np.all(result.resolved_x_m[selected] < source.x_edges_m[width + 1])
    assert math.fsum(result.particle_weight) == pytest.approx(result.total_released_particles, rel=3e-15)
    assert result.unequal_macro_weights
    np.testing.assert_array_equal(result.momentum_kg_m_s, result.mass_kg[:, None] * result.velocity_m_s)
    assert result.surface_attributes == source.surface_attributes
    assert result.parameter_provenance == source.parameter_provenance


def test_pairing_changes_width_mapping_without_resampling_times_states_weights_or_velocities():
    uniforms = draws(9)
    uniforms[:, 1] = np.linspace(.1, .9, len(uniforms))
    uniforms[:, 2] = np.linspace(.15, .85, len(uniforms))
    source = history()
    first = sample(source, uniforms)
    modified = uniforms.copy()
    modified[:, 2] = .2
    second = sample(source, modified)
    for name in ("time_offset_s", "time_s", "velocity_m_s", "momentum_kg_m_s", "particle_weight",
                 "source_width_index", "reservoir_index", "temperature_k"):
        np.testing.assert_array_equal(getattr(first, name), getattr(second, name))
    assert first.state_ids == second.state_ids
    assert np.any(first.resolved_x_m != first.uniform_x_m)
    assert np.any(first.resolved_x_m != second.resolved_x_m)
    for name in ("time_offset_s", "particle_weight", "velocity_m_s"):
        with pytest.raises(ValueError):
            getattr(first, name).flags.writeable = True
    with pytest.raises(TypeError):
        first.surface_attributes["run_id"] = "changed"


def test_conditional_times_match_independent_nonlinear_ode_and_differ_from_initial_hazard():
    source = history(temperatures=(600,), coverage=((1e18,),), coefficient=.8, edges=(0, .004))
    uniforms = draws(17)
    probabilities = np.linspace(.05, .95, len(uniforms))
    uniforms[:, 1] = probabilities
    result = sample(source, uniforms)
    pool = source.reservoirs[0]
    beta = pool.activation_energy_ev * QE / (KB * 600)
    duration = source.time_edges_s[-1] - source.time_edges_s[0]
    # Integrate normalized remaining coverage independently of the implicit
    # release/time inversion implementation.
    solution = solve_ivp(lambda t, theta: -pool.attempt_frequency_s * theta * np.exp(-beta * (1 - .8 * theta)),
                         [0, duration], [1.], rtol=2e-12, atol=1e-14, dense_output=True)
    assert solution.success
    released_fraction = 1 - solution.y[0, -1]
    actual_cdf = (1 - solution.sol(result.time_offset_s)[0]) / released_fraction
    np.testing.assert_allclose(actual_cdf, probabilities, rtol=3e-10, atol=2e-11)
    initial_rate = pool.attempt_frequency_s * math.exp(-beta * (1 - .8))
    incorrect = -np.log1p(-probabilities * (-math.expm1(-initial_rate * duration))) / initial_rate
    assert np.max(np.abs(incorrect - result.time_offset_s)) > .02 * duration


def test_outward_flux_distribution_thermal_moments_and_mass_scaling():
    source = history(temperatures=(600,), coverage=((1e18,),), coefficient=0, edges=(0, .004))
    n = 8192
    uniforms = draws(n)
    quantiles = (np.arange(n) + .5) / n
    uniforms[:, 3] = quantiles
    uniforms[:, 4] = quantiles
    uniforms[:, 5] = (np.arange(n) * .6180339887498949) % 1
    uniforms[0, 5] = .5
    result = sample(source, uniforms)
    normal_cdf = -np.expm1(-result.mass_kg * result.velocity_m_s[:, 1] ** 2 / (2 * KB * 600))
    np.testing.assert_allclose(normal_cdf, quantiles, rtol=5e-14, atol=5e-15)
    assert np.all(result.velocity_m_s[:, 1] > 0)
    kinetic = .5 * result.mass_kg * np.sum(result.velocity_m_s ** 2, axis=1)
    assert kinetic.mean() / (2 * KB * 600) == pytest.approx(1, rel=5e-4)
    assert result.component_sampled_translational_energy_j[0] / result.component_expected_translational_energy_j[0] == pytest.approx(1, rel=5e-4)
    heavier = sample_temkin_emission_packets(source, {"M:0": 12e-26}, 0, uniforms,
                                             outward_normal_sign=-1, current_length_m=.025)
    np.testing.assert_allclose(np.abs(heavier.velocity_m_s), .5 * np.abs(result.velocity_m_s), rtol=2e-14)
    assert np.all(heavier.velocity_m_s[:, 1] < 0)
    np.testing.assert_allclose(heavier.momentum_kg_m_s[:, 0], 2 * result.momentum_kg_m_s[:, 0], rtol=2e-14)
    np.testing.assert_array_equal(heavier.time_offset_s, result.time_offset_s)


def test_rare_components_are_seeded_and_zero_source_generates_no_packets():
    source = history(coverage=((1e2, 1e18),))
    result = sample(source, draws(6))
    np.testing.assert_array_equal(result.packet_counts, [1, 5])
    assert result.unequal_macro_weights
    assert np.all(result.particle_weight > 0)
    with pytest.raises(ValueError, match="at least 2"):
        sample(source, draws(1))
    with pytest.raises(NoTemkinEmissionError):
        sample(history(coverage=((0, 0),)))


def test_selected_interval_is_replayed_before_emission_and_malformed_metadata_rejects():
    source = history()
    object.__setattr__(source, "released_m2", source.released_m2 * 1.1)
    with pytest.raises(ValueError):
        sample(source)
    source = history()
    attributes = dict(source.surface_attributes)
    attributes.pop("current_averaging_interval_m_json")
    object.__setattr__(source, "surface_attributes", attributes)
    with pytest.raises(ValueError, match="missing required"):
        sample(source)


def test_later_interval_preserves_original_index_midpoint_temperature_and_inventory():
    surface = FlashSurfaceHistory([0, 1e-8, 3e-8], [0, .004], [[450], [700], [950]],
                                   "synthetic", "later-interval", "manufactured temperature ramp", "test face")
    pool = history().reservoirs[0]
    source = evolve_temkin_release(surface, SurfaceInventory([[1e18]]), [pool], parameter_provenance="manufactured ramp")
    result = sample_temkin_emission_packets(source, {pool.state_id: 3e-26}, 1, draws(5),
                                             outward_normal_sign=1, current_length_m=.025)
    assert result.interval_index == 1
    assert result.time_start_s == 1e-8
    assert result.time_end_s == 3e-8
    np.testing.assert_array_equal(result.temperature_k, 825.)
    assert result.total_released_particles == pytest.approx(source.released_m2[1, 0, 0] * .004 * .025, rel=3e-15)


def test_constant_hazard_history_cannot_be_reinterpreted_as_temkin_packets():
    surface = FlashSurfaceHistory([0, 1e-7], [0, .004], [[600], [600]],
                                   "synthetic", "baseline", "manufactured baseline source", "test face")
    baseline = evolve_surface_release(surface, SurfaceInventory([[1e18]]), [BindingReservoir("base", "M", 1e8, .3, .5)],
                                      parameter_provenance="manufactured baseline")
    with pytest.raises(ValueError, match="distinct validated TemkinReleaseHistory"):
        sample_temkin_emission_packets(baseline, {"M": 3e-26}, 0, draws(5),
                                       outward_normal_sign=1, current_length_m=.025)


@pytest.mark.parametrize("column,value", [(0, 0), (1, 1), (2, -.1), (3, np.nan), (4, True), (5, .5j)])
def test_uniforms_require_strict_open_real_variates(column, value):
    uniforms = np.full((8, 6), .5, dtype=object)
    uniforms[0, column] = value
    with pytest.raises(ValueError):
        sample(uniforms=uniforms)


def test_masked_variates_cannot_be_silently_supplied():
    uniforms = np.ma.array(draws(8), mask=False)
    uniforms.mask[0, 1] = True
    with pytest.raises(ValueError, match="masked"):
        sample(uniforms=uniforms)


def test_exact_state_masses_valid_geometry_and_nonrelativistic_range_are_required():
    source = history()
    for mapping in ({}, {"M": 3e-26}, {"M:0": 0}, {"M:0": True}, {"M:0": 1j}, {"M:0": 1e-40}):
        with pytest.raises(ValueError):
            sample_temkin_emission_packets(source, mapping, 0, draws(8), outward_normal_sign=1, current_length_m=.025)
    for arguments in ({"outward_normal_sign": 0}, {"outward_normal_sign": True},
                      {"current_length_m": 0}, {"current_length_m": "0.025"}):
        with pytest.raises(ValueError):
            sample(source, **arguments)
    for index in (-1, 1, True):
        with pytest.raises(ValueError):
            sample_temkin_emission_packets(source, {"M:0": 3e-26}, index, draws(8),
                                            outward_normal_sign=1, current_length_m=.025)
