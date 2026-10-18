"""Deterministic boundary-contract checks, without a gap or radiation solve."""

import math

import numpy as np
import pytest

from mhdlab.constants import KB
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.surface_release_history import evolve_surface_release
from mhdlab.thermal_boundary import NoThermalEmissionError, sample_emission_packets
from mhdlab.thermal_release import BindingReservoir, SurfaceInventory


def history(temperatures=(300, 900), coverage=((1e18, 1e18),), rate=2., times=(0., 1.), edges=(0., .001, .004)):
    surface = FlashSurfaceHistory(times, edges, [temperatures, temperatures], "synthetic", "boundary-verification",
                                  "explicit analytic parameters", "single planar face")
    pools = [BindingReservoir(f"pool-{index}", f"state-{index}", rate, 0., 0.5, 0.2)
             for index in range(len(coverage))]
    return evolve_surface_release(surface, SurfaceInventory(coverage), pools, parameter_provenance="analytic fixture")


def draws(number):
    return np.full((number, 6), .5)


def sample(source, uniform_draws=None, **kwargs):
    values = dict(outward_normal_sign=1, current_length_m=.025)
    values.update(kwargs)
    mass_map = {pool.state_id: (index + 1) * 3e-26 for index, pool in enumerate(source.reservoirs)}
    return sample_emission_packets(source, mass_map, 0, draws(8) if uniform_draws is None else uniform_draws, **values)


def test_exact_component_inventory_weights_with_paired_source_temperature_mixture():
    source = history()
    u = draws(8)
    u[:, 0] = [.1, .9, .1, .4, .6, .7, .8, .9]
    u[:, 2] = np.linspace(.1, .9, 8)
    packets = sample(source, u)
    assert packets.required_minimum_n == 2
    np.testing.assert_array_equal(packets.packet_counts, [2, 6])
    np.testing.assert_allclose(packets.component_probabilities, [.25, .75], rtol=2e-15)
    assert packets.unequal_macro_weights is False
    for component, (pool, width) in enumerate(zip(packets.component_reservoir_index, packets.component_width_index)):
        selected = (packets.reservoir_index == pool) & (packets.source_width_index == width)
        expected = source.released_m2[0, pool, width] * np.diff(source.x_edges_m)[width] * .025
        assert math.fsum(packets.particle_weight[selected]) == pytest.approx(expected, rel=2e-15)
        assert packets.component_released_particles[component] == pytest.approx(expected, rel=2e-15)
        np.testing.assert_array_equal(packets.temperature_k[selected], source.temperature_k[0, width])
        assert np.all(packets.resolved_x_m[selected] >= source.x_edges_m[width])
        assert np.all(packets.resolved_x_m[selected] <= source.x_edges_m[width + 1])
    assert math.fsum(packets.particle_weight) == pytest.approx(packets.total_released_particles, rel=2e-15)
    np.testing.assert_allclose(packets.uniform_x_m, .004 * u[:, 2], rtol=2e-15)
    # Both position columns refer to this one set of state/time/momentum arrays;
    # choosing the uniform mapping does not draw or alter those quantities.
    np.testing.assert_array_equal(packets.momentum_kg_m_s, packets.velocity_m_s * packets.mass_kg[:, None])
    repeated = sample(source, u.copy())
    for name in ("state_ids", "reservoir_names"):
        assert getattr(repeated, name) == getattr(packets, name)
    for name in ("time_offset_s", "velocity_m_s", "momentum_kg_m_s", "particle_weight"):
        np.testing.assert_array_equal(getattr(repeated, name), getattr(packets, name))


def test_stratification_retains_rare_components_without_overdrawing_inventory():
    source = history(coverage=((1e10, 1e18),))
    packets = sample(source, draws(6))
    np.testing.assert_array_equal(packets.packet_counts, [1, 5])
    assert packets.unequal_macro_weights
    assert np.all(packets.particle_weight > 0)
    for component, width in enumerate(packets.component_width_index):
        selected = packets.source_width_index == width
        assert math.fsum(packets.particle_weight[selected]) == pytest.approx(packets.component_released_particles[component], rel=2e-15)


def test_exponential_emission_time_inverse_and_small_hazard_conditional_limit():
    source = history()
    u = draws(8)
    u[:, 1] = np.linspace(.1, .9, 8)
    packets = sample(source, u)
    expected = -np.log1p(-u[:, 1] * (-np.expm1(-2.0))) / 2.0
    np.testing.assert_allclose(packets.time_offset_s, expected, rtol=3e-15)
    assert not np.allclose(packets.time_offset_s, u[:, 1])
    slow = history(rate=1e-200, times=(0., 1e-200), coverage=((1e300, 1e300),))
    slow_packets = sample(slow, u)
    np.testing.assert_allclose(slow_packets.time_offset_s, 1e-200 * u[:, 1], rtol=2e-15, atol=0)
    assert slow_packets.total_released_particles > 0


def test_flux_maxwellian_energy_and_mass_scaling_by_deterministic_quadrature():
    source = history(temperatures=(600,), coverage=((1e18,),), edges=(0., .004))
    number = 8192
    u = draws(number)
    quantiles = (np.arange(number) + .5) / number
    u[:, 3] = quantiles
    u[:, 4] = quantiles
    u[:, 5] = (np.arange(number) * 0.6180339887498949) % 1
    u[0, 5] = .5  # all variates are strictly open
    packets = sample(source, u)
    energies = .5 * packets.mass_kg * np.sum(packets.velocity_m_s ** 2, axis=1)
    assert energies.mean() / (2 * KB * 600) == pytest.approx(1, rel=5e-4)
    weighted_energy = math.fsum(energies * packets.particle_weight)
    assert weighted_energy == pytest.approx(packets.component_sampled_translational_energy_j[0], rel=2e-15)
    assert packets.component_sampled_translational_energy_j[0] / packets.component_expected_translational_energy_j[0] == pytest.approx(1, rel=5e-4)
    assert packets.component_translational_energy_sampling_residual_j[0] == (
        packets.component_sampled_translational_energy_j[0] - packets.component_expected_translational_energy_j[0])
    heavier = sample_emission_packets(source, {"state-0": 12e-26}, 0, u,
                                      outward_normal_sign=-1, current_length_m=.025)
    np.testing.assert_allclose(np.abs(heavier.velocity_m_s), .5 * np.abs(packets.velocity_m_s), rtol=2e-14)
    assert np.all(heavier.velocity_m_s[:, 1] < 0)
    np.testing.assert_allclose(heavier.momentum_kg_m_s[:, 0], 2 * packets.momentum_kg_m_s[:, 0], rtol=2e-14)


def test_mixed_source_temperature_energy_keeps_origin_components_and_weighted_diagnostics():
    source = history(temperatures=(300, 900), coverage=((1e18, 1e18),))
    n = 2048
    u = draws(n)
    u[:, 0] = (np.arange(n) + .5) / n
    initial = sample(source, u)
    for width in (0, 1):
        selected = np.flatnonzero(initial.source_width_index == width)
        quantiles = (np.arange(len(selected)) + .5) / len(selected)
        u[selected, 3] = quantiles
        u[selected, 4] = quantiles
    packets = sample(source, u)
    expected = packets.component_released_particles * (2 * KB * np.array([300, 900]))
    np.testing.assert_allclose(packets.component_expected_translational_energy_j, expected, rtol=2e-15)
    np.testing.assert_allclose(packets.component_sampled_translational_energy_j, expected, rtol=2e-3)
    assert np.unique(packets.temperature_k).tolist() == [300, 900]


def test_internal_state_and_mass_map_survive_multiple_reservoirs():
    source = history(coverage=((1e18, 2e18), (3e18, 4e18)))
    packets = sample(source, draws(4))
    assert packets.required_minimum_n == 4
    assert packets.state_ids == ("state-0", "state-0", "state-1", "state-1")
    np.testing.assert_array_equal(packets.mass_kg, [3e-26, 3e-26, 6e-26, 6e-26])
    np.testing.assert_array_equal(packets.packet_counts, [1, 1, 1, 1])
    assert packets.parameter_provenance == source.parameter_provenance
    assert packets.surface_attributes["run_id"] == source.surface_attributes["run_id"]
    with pytest.raises(TypeError):
        packets.surface_attributes["run_id"] = "changed"


def test_separate_offsets_retain_subtime_information_and_outputs_are_immutable():
    source = history(times=(1e8, 1e8 + 1), rate=1e8)
    packets = sample(source)
    assert np.all(packets.time_offset_s > 0)
    assert np.all(packets.time_offset_s < 1e-6)
    # Absolute best-float times can lose part/all of the offset; that does not
    # authorize replacing the local offset before handing it downstream.
    assert np.any(packets.time_s - packets.time_start_s != packets.time_offset_s)
    for name in ("time_s", "time_offset_s", "particle_weight", "packet_counts", "velocity_m_s"):
        with pytest.raises(ValueError):
            getattr(packets, name).flags.writeable = True


def test_no_emission_and_insufficient_packet_count_do_not_generate_particles():
    with pytest.raises(NoThermalEmissionError):
        sample(history(coverage=((0, 0),)))
    with pytest.raises(ValueError, match="at least 2"):
        sample(history(), draws(1))


@pytest.mark.parametrize("column,value", [(0, 0), (1, 1), (2, -.1), (3, np.nan), (4, True), (5, .5j)])
def test_uniforms_are_strict_open_real_variates(column, value):
    u = np.full((8, 6), .5, dtype=object)
    u[0, column] = value
    with pytest.raises(ValueError):
        sample(history(), u)


def test_exact_state_masses_and_geometry_inputs_are_required():
    source = history()
    for mass_map in ({}, {"state": 3e-26}, {"state-0": True}, {"state-0": 0}, {"state-0": 1j}):
        with pytest.raises(ValueError):
            sample_emission_packets(source, mass_map, 0, draws(4), outward_normal_sign=1, current_length_m=.025)
    for kwargs in ({"outward_normal_sign": 0}, {"outward_normal_sign": True},
                   {"current_length_m": 0}, {"current_length_m": "0.025"}):
        with pytest.raises(ValueError):
            sample(source, **kwargs)
    for index in (-1, 1, True):
        with pytest.raises(ValueError):
            sample_emission_packets(source, {"state-0": 3e-26}, index, draws(4),
                                      outward_normal_sign=1, current_length_m=.025)
