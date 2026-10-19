"""Conditioned packet inventories, thermal mixtures and representation guards."""

from dataclasses import FrozenInstanceError, fields, replace
from fractions import Fraction
import json
import math
from types import MappingProxyType

import numpy as np
import pytest

from mhdlab.constants import KB, QE
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.source_convergence import manufactured_planar_surface
from mhdlab.temkin_planar import TemkinPlanarReleaseHistory
from mhdlab.temkin_history import STEP_FIELDS, DIAGNOSTIC_FIELDS
from mhdlab.temkin_release import TemkinReservoir, release_step
from mhdlab.thermal_release import SurfaceInventory
from mhdlab.temkin_windowed_boundary import (
    TemkinWindowedBoundaryRangeError, TemkinWindowedEmissionPackets, sample_temkin_windowed_emission_packets,
)
from mhdlab.temkin_windowed_source import SourceWindowSelection, condition_temkin_planar_source


def pool(name="pool_a", state="M:a", **changes):
    parameters = dict(name=name, state_id=state, attempt_frequency_s=2., activation_energy_ev=0.,
                      desorption_energy_ev=.3, internal_energy_ev=.1, coverage_coefficient=.35,
                      reference_coverage_m2=2e18, max_coverage_m2=2e20,
                      kinetic_provenance="Manufactured fixed-reference rate fixture",
                      energy_provenance="Manufactured independent wall/state costs")
    parameters.update(changes)
    return TemkinReservoir(**parameters)


def source_record(*, pools=None, initial=None, temperature=None, x_edges=None, z_edges=None,
                  origin=10., duration=2., windows=None, attributes=None):
    pools = [pool()] if pools is None else pools
    x_edges = np.array([0., .002]) if x_edges is None else np.asarray(x_edges)
    z_edges = np.array([0., .004]) if z_edges is None else np.asarray(z_edges)
    shape = (len(z_edges)-1, len(x_edges)-1)
    initial = np.full((len(pools), *shape), 1e18) if initial is None else np.asarray(initial)
    temperature = np.full(shape, 700.) if temperature is None else np.asarray(temperature)
    base = manufactured_planar_surface(width_bins=1, current_bins=1, snapshot_intervals=1,
                                      duration_s=1., width_m=.002, current_length_m=.004)
    metadata = dict(base.attributes) if attributes is None else dict(attributes)
    for name in ("current_direction_representation", "current_patch_intervals_m", "patch_temperature_approximation"):
        metadata.pop(name, None)
    surface = PlanarSurfaceHistory(np.array([origin, origin+duration]), x_edges, z_edges,
                                   np.stack([temperature, temperature]), metadata)
    # This fixture is one explicitly frozen interval. Direct canonical assembly
    # avoids an unrelated midpoint interpolation of an already constant field,
    # which loses minimum-subnormal temperatures before reaching a sampler.
    trial = release_step(SurfaceInventory(initial.reshape(len(pools), -1)), pools,
                         temperature.reshape(-1), duration)
    history = TemkinPlanarReleaseHistory(
        time_edges_s=surface.time_s, x_edges_m=x_edges, z_edges_m=z_edges,
        temperature_k=temperature[None],
        coverage_m2=np.stack([initial, trial.trial_inventory.coverage_m2.reshape(initial.shape)]),
        **{name: getattr(trial, name).reshape((1, *initial.shape)) for name in (*STEP_FIELDS, *DIAGNOSTIC_FIELDS)},
        reservoirs=tuple(pools), surface_attributes=surface.attributes,
        parameter_provenance="Manufactured packet source; no measured physical parameters")
    selections = [SourceWindowSelection(0, .1, .7)] if windows is None else [SourceWindowSelection(0, *bounds) for bounds in windows]
    return condition_temkin_planar_source(history, selections, log_mass_tolerance=1e-7,
                                         status="manufactured", provenance="Manufactured conditioned packet fixture")


def draws(number, seed=713):
    return .05 + .9 * np.random.default_rng(seed).random((number, 7))


def sample(source, values=None, masses=None, **changes):
    if values is None:
        values = draws(max(32, len(source.components)))
    if masses is None:
        masses = {component.source_window.reservoir.state_id: 3e-26 for component in source.components}
    parameters = dict(face_normal_coordinate_m=.0004, outward_normal_sign=-1)
    parameters.update(changes)
    return sample_temkin_windowed_emission_packets(source, masses, values, **parameters)


@pytest.fixture(scope="module")
def mixed_source():
    return source_record(pools=[pool(), pool("pool_b", "M:b", attempt_frequency_s=.8, internal_energy_ev=.2)],
                         initial=[[[1e18, 0.], [2e18, 1.5e18]], [[.5e18, 1e18], [0., 3e18]]],
                         temperature=[[700., 900.], [1100., 600.]],
                         x_edges=[-.002, 0., .004], z_edges=[0., .003, .008],
                         windows=[(.1, .7), (1., 1.8)])


def test_mandatory_components_counts_state_temperature_and_exact_weight_closure(mixed_source):
    variates = draws(64)
    masses = {"M:a": 3e-26, "M:b": 6e-26}
    packets = sample(mixed_source, variates, masses)
    required = len(mixed_source.components)
    assert packets.required_minimum_n == required == 12
    np.testing.assert_array_equal(packets.component_index[:required], np.arange(required))
    assert np.all(packets.packet_counts >= 1)
    assert int(np.sum(packets.packet_counts)) == len(variates)
    assert packets.unequal_macro_weights
    for index, component in enumerate(mixed_source.components):
        selected = packets.component_index == index
        window = component.source_window
        count = window.manifest()["quantities"]["released_number"]["value"]
        assert packets.component_released_particles[index] == count
        assert packets.packet_counts[index] == np.count_nonzero(selected)
        assert all(state == window.reservoir.state_id for state in np.array(packets.state_ids)[selected])
        assert np.all(packets.temperature_k[selected] == window.temperature_k)
        assert np.all(packets.mass_kg[selected] == masses[window.reservoir.state_id])
        np.testing.assert_allclose(np.sum(packets.particle_weight[selected]), count, rtol=2e-15)
        exact = sum((Fraction(float(value)) for value in packets.particle_weight[selected]), Fraction()) - Fraction(count)
        assert packets.component_particle_balance_residual[index] == float(exact)
    assert packets.total_particle_balance_residual == float(
        sum((Fraction(float(value)) for value in packets.particle_weight), Fraction())
        - sum((Fraction(float(value)) for value in packets.component_released_particles), Fraction()))


def test_paired_hypotheses_change_width_only_and_keep_original_local_timing(mixed_source):
    variates = draws(48)
    packets = sample(mixed_source, variates)
    np.testing.assert_array_equal(packets.resolved_position_m[:, 1:], packets.uniform_width_position_m[:, 1:])
    np.testing.assert_array_equal(packets.resolved_position_m[:, 0], packets.resolved_x_m)
    np.testing.assert_array_equal(packets.uniform_width_position_m[:, 0], packets.uniform_x_m)
    assert np.any(packets.resolved_x_m != packets.uniform_x_m)
    assert np.all(packets.resolved_position_m[:, 1] == .0004)
    assert np.all(packets.velocity_m_s[:, 1] < 0)
    for packet, component_index in enumerate(packets.component_index):
        component = mixed_source.components[component_index]
        window = component.source_window
        assert packets.selection_index[packet] == component.selection_index
        assert packets.interval_index[packet] == component.interval_index
        assert packets.reservoir_index[packet] == component.reservoir_index
        assert packets.source_current_index[packet] == component.current_index
        assert packets.source_width_index[packet] == component.width_index
        assert packets.interval_origin_s[packet] == component.interval_origin_s
        assert packets.window_start_offset_s[packet] == window.start_s
        assert packets.window_stop_offset_s[packet] == window.stop_s
        u = variates[packet]
        duration = window.stop_s-window.start_s
        rate = window.reservoir.attempt_frequency_s
        expected_delay = -math.log1p(-u[1] * -math.expm1(-rate*duration))/rate
        assert packets.local_time_delay_s[packet] == pytest.approx(expected_delay, rel=2e-11)
        assert packets.log_local_time_delay_s[packet] == pytest.approx(math.log(expected_delay), abs=2e-11)
        assert component.x_bounds_m[0] <= packets.resolved_x_m[packet] <= component.x_bounds_m[1]
        assert packets.resolved_width_offset_m[packet] == pytest.approx(np.diff(component.x_bounds_m)[0]*u[2], rel=2e-15)
        assert packets.source_current_offset_m[packet] == pytest.approx(np.diff(component.z_bounds_m)[0]*u[6], rel=2e-15)
    assert np.all(packets.global_time_precision_qualified)


def test_three_velocity_components_energy_sampling_and_state_costs_are_independent(mixed_source):
    variates = np.full((24, 7), .5)
    variates[:, 3] = .2
    variates[:, 4] = .8
    variates[:, 5] = .125
    packets = sample(mixed_source, variates)
    sigma = np.sqrt(KB*packets.temperature_k/packets.mass_kg)
    expected_normal = -sigma*np.sqrt(-2*np.log1p(-variates[:, 3]))
    expected_radius = sigma*np.sqrt(-2*np.log(variates[:, 4]))
    np.testing.assert_allclose(packets.velocity_m_s[:, 0], expected_radius / math.sqrt(2), rtol=2e-14)
    np.testing.assert_allclose(packets.velocity_m_s[:, 1], expected_normal, rtol=2e-14)
    np.testing.assert_allclose(packets.velocity_m_s[:, 2], expected_radius / math.sqrt(2), rtol=2e-14)
    np.testing.assert_allclose(packets.momentum_kg_m_s, packets.mass_kg[:, None] * packets.velocity_m_s, rtol=2e-15)
    sampled_total = 0.
    for index, component in enumerate(mixed_source.components):
        selected = packets.component_index == index
        actual = math.fsum(.5*float(m)*float(w)*math.hypot(*velocity)**2
                           for m, w, velocity in zip(packets.mass_kg[selected], packets.particle_weight[selected], packets.velocity_m_s[selected]))
        count = packets.component_released_particles[index]
        expected = 2*KB*component.source_window.temperature_k*count
        assert packets.component_sampled_translational_energy_j[index] == pytest.approx(actual, rel=3e-14)
        assert packets.component_expected_translational_energy_j[index] == pytest.approx(expected, rel=3e-12)
        assert packets.component_translational_energy_sampling_residual_j[index] == math.fsum((packets.component_sampled_translational_energy_j[index], -packets.component_expected_translational_energy_j[index]))
        assert packets.component_sampled_translational_energy_j[index] < expected / 2
        assert packets.component_expected_internal_energy_j[index] == pytest.approx(count * QE * component.source_window.reservoir.internal_energy_ev, rel=3e-12)
        assert packets.component_expected_binding_energy_j[index] == pytest.approx(count * QE * component.source_window.reservoir.desorption_energy_ev, rel=3e-12)
        sampled_total += actual
    assert packets.sampled_translational_energy_j == pytest.approx(sampled_total, rel=3e-14)
    assert packets.translational_energy_sampling_residual_j < 0
    assert packets.wall_energy_sampling_residual_j == math.fsum((packets.sampled_total_wall_energy_debit_j, -packets.expected_total_wall_energy_debit_j))


def test_record_arrays_and_manifest_are_immutable_or_detached(mixed_source):
    packets = sample(mixed_source)
    manifest = packets.manifest()
    arrays = {field.name: getattr(packets, field.name) for field in fields(packets)
              if isinstance(getattr(packets, field.name), np.ndarray)}
    assert set(manifest["array_units"]) == set(arrays)
    for array in arrays.values():
        assert not array.flags.writeable
        with pytest.raises(ValueError):
            array.setflags(write=True)
    assert manifest["host_budget_accepted"] is False
    assert manifest["source"]["host_budget_accepted"] is False
    assert manifest["array_units"]["log_local_time_delay_s"] == "ln(delay / 1 s)"
    manifest["source"]["components"].clear()
    manifest["components"][0]["packet_weight"] = -1.
    assert packets.manifest()["source"]["components"]
    assert packets.manifest()["components"][0]["packet_weight"] > 0
    json.dumps(packets.manifest(), allow_nan=False)
    with pytest.raises((FrozenInstanceError, AttributeError)):
        packets.total_released_particles = 0.
    with pytest.raises(TypeError):
        TemkinWindowedEmissionPackets()


def test_large_global_origin_does_not_erase_retained_local_delay():
    source = source_record(origin=1e16, duration=2., windows=[(.01, .1)])
    packets = sample(source)
    assert np.all(packets.best_effort_time_s == 1e16)
    assert not np.any(packets.global_time_precision_qualified)
    assert np.all(packets.local_time_delay_s > 0)
    assert np.ptp(packets.local_time_delay_s) > 0
    assert np.all(packets.interval_origin_s == 1e16)
    assert np.all(packets.window_start_offset_s == .01)
    assert np.all(packets.global_time_absolute_rounding_error_s > 0)


def test_late_window_preserves_a_count_whose_original_release_fraction_underflows():
    source = source_record(pools=[pool(attempt_frequency_s=1., reference_coverage_m2=1e300, max_coverage_m2=1e300)],
                           initial=[[[1e300]]], duration=1000., origin=0., windows=[(800., 801.)])
    packets = sample(source)
    window = source.components[0].source_window
    assert math.exp(window.time_law.log_source_fraction) == 0
    assert packets.total_released_particles > 0
    expected_log = math.log(1e300) + math.log(window.patch_area_m2) - 800 + math.log(-math.expm1(-1.))
    assert math.log(packets.total_released_particles) == pytest.approx(expected_log, abs=3e-10)
    assert np.all(packets.window_start_offset_s == 800.)


def test_extensive_release_is_sampled_when_only_areal_quantities_underflow():
    source = source_record(pools=[pool(attempt_frequency_s=1.)], initial=[[[1.]]], duration=1000., origin=0.,
                           x_edges=[0., 1e150], z_edges=[0., 1e150], windows=[(800., 801.)])
    window = source.components[0].source_window
    assert window.manifest()["quantities"]["released_m2"]["value"] == 0
    assert not window.all_budget_values_qualified
    packets = sample(source)
    assert packets.total_released_particles > 0
    assert packets.expected_translational_energy_j > 0


def test_positive_probability_underflow_is_flagged_without_omitting_mandatory_component():
    source = source_record(pools=[pool(max_coverage_m2=1e200, reference_coverage_m2=1e200)],
                           initial=[[[1e-280, 1e100]]], x_edges=[0., 1., 2.], z_edges=[0., 1.])
    packets = sample(source, draws(8))
    assert packets.component_probabilities[0] == 0
    assert packets.component_probability_underflowed[0]
    assert not packets.component_probability_representation_qualified[0]
    assert math.isfinite(packets.component_log_probabilities[0])
    assert packets.packet_counts[0] >= 1
    assert packets.component_released_particles[0] > 0
    assert np.all(packets.particle_weight[packets.component_index == 0] > 0)


def test_repeated_packet_weights_with_coarse_subnormal_rounding_are_rejected():
    source = source_record(pools=[pool(attempt_frequency_s=1., desorption_energy_ev=0., internal_energy_ev=0.)],
                           initial=[[[2.**-500]]], temperature=[[1e300]], duration=1., origin=0.,
                           x_edges=[0., 2.**-286], z_edges=[0., 2.**-286], windows=[(0., math.log(2.))])
    mass = {"M:a": 1e300}
    one = sample(source, np.full((1, 7), .5), mass)
    assert one.total_released_particles == 2*np.nextafter(0., 1.)
    with pytest.raises(TemkinWindowedBoundaryRangeError, match="packet weight"):
        sample(source, np.full((3, 7), .5), mass)


def test_required_local_delay_rejects_unqualified_subnormal_projection():
    source = source_record(pools=[pool(attempt_frequency_s=1., desorption_energy_ev=0., internal_energy_ev=0.,
                                     max_coverage_m2=1e300, reference_coverage_m2=1e300)],
                           initial=[[[1e300]]], temperature=[[1e300]], duration=2e-323, origin=0.,
                           windows=[(0., 2e-323)])
    variates = np.full((1, 7), .5)
    variates[:, 1] = .1
    with pytest.raises(ValueError, match="quantile|quantization"):
        sample(source, variates, {"M:a": 1e300})


@pytest.mark.parametrize("temperature,mass", [(1., 1e-300), (1e300, 1e-30)])
def test_nonrelativistic_limit_is_explicit(temperature, mass):
    source = source_record(temperature=[[temperature]])
    with pytest.raises(TemkinWindowedBoundaryRangeError, match="nonrelativistic"):
        sample(source, np.full((1, 7), .5), {"M:a": mass})


@pytest.mark.parametrize("temperature", [float(np.nextafter(0., 1.)), 1.5e-301])
def test_positive_momentum_underflow_or_coarse_rounding_is_rejected(temperature):
    source = source_record(pools=[pool(max_coverage_m2=1e300, reference_coverage_m2=1e300,
                                     internal_energy_ev=0., desorption_energy_ev=0.)],
                           initial=[[[1e300]]], temperature=[[temperature]], origin=0., duration=1., windows=[(0., 1.)])
    variates = np.full((1, 7), .5)
    variates[:, 5] = .125
    with pytest.raises(TemkinWindowedBoundaryRangeError, match="momentum"):
        sample(source, variates, {"M:a": float(np.nextafter(0., 1.))})


def test_underflowed_or_coarsely_projected_thermal_speed_is_rejected():
    source = source_record(pools=[pool(max_coverage_m2=1e300, reference_coverage_m2=1e300,
                                     internal_energy_ev=0., desorption_energy_ev=0.)],
                           initial=[[[1e300]]], temperature=[[float(np.nextafter(0., 1.))]],
                           origin=0., duration=1., windows=[(0., 1.)])
    with pytest.raises(TemkinWindowedBoundaryRangeError, match="speed|radius"):
        sample(source, np.full((1, 7), .5), {"M:a": 1e301})


def test_unqualified_positive_component_energy_rejects_the_whole_source():
    source = source_record(pools=[pool(internal_energy_ev=0., desorption_energy_ev=0.)],
                           initial=[[[1e-300]]], temperature=[[1e-100]])
    with pytest.raises(TemkinWindowedBoundaryRangeError, match="every positive"):
        sample(source)


@pytest.mark.parametrize("bad", [np.zeros((4, 7)), np.ones((4, 7)), np.full((4, 7), math.nan),
                                  np.full((4, 6), .5), np.empty((0, 7)), np.full((4, 7), True),
                                  np.ma.array(np.full((4, 7), .5), mask=True)])
def test_invalid_explicit_variates_are_rejected(bad):
    with pytest.raises(ValueError):
        sample(source_record(), bad)


@pytest.mark.parametrize("changes", [dict(outward_normal_sign=True), dict(outward_normal_sign=0),
                                     dict(face_normal_coordinate_m=True), dict(face_normal_coordinate_m=math.inf)])
def test_face_inputs_are_explicit_and_finite(changes):
    with pytest.raises(ValueError):
        sample(source_record(), **changes)


def test_missing_exact_state_mass_and_insufficient_component_budget_are_rejected(mixed_source):
    with pytest.raises(ValueError, match="exact-state"):
        sample(mixed_source, masses={"M": 3e-26})
    with pytest.raises(ValueError, match="at least"):
        sample(mixed_source, draws(len(mixed_source.components)-1))


@pytest.mark.parametrize("bad", [0., -1., True, math.nan, math.inf, np.ma.masked])
def test_invalid_exact_state_mass_is_rejected(bad):
    with pytest.raises(ValueError):
        sample(source_record(), masses={"M:a": bad})


def test_forged_cached_window_numbers_are_replayed_and_rejected():
    source = source_record()
    window = source.components[0].source_window
    records = dict(window._quantities)
    records["released_number"] = replace(records["released_number"], value=1.)
    object.__setattr__(window, "_quantities", MappingProxyType(records))
    with pytest.raises(ValueError, match="modified|association|canonical|source-window"):
        sample(source)


def test_forged_component_temperature_is_rejected_against_canonical_baseline():
    source = source_record()
    object.__setattr__(source.components[0].source_window, "temperature_k", 701.)
    with pytest.raises(ValueError):
        sample(source)


def test_wrong_source_type_and_empty_conditioned_population_are_rejected():
    with pytest.raises(ValueError, match="TemkinPlanarSourceWindows"):
        sample_temkin_windowed_emission_packets(None, {}, np.full((1, 7), .5), face_normal_coordinate_m=0., outward_normal_sign=1)
    source = source_record(initial=[[[0.]]])
    with pytest.raises(ValueError, match="no positive"):
        sample(source)
