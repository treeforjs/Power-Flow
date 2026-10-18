"""Qualified Temkin histories, finite-volume moments and temporal refinement."""

from dataclasses import replace
import math

import numpy as np
import pytest

from mhdlab.constants import KB, QE
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.surface_release_history import SurfaceReleaseHistory
from mhdlab.temkin_history import TemkinReleaseHistory, evolve_temkin_release
from mhdlab.temkin_release import LAW_IDENTITY, TemkinReservoir, energy_budgets
from mhdlab.thermal_boundary import sample_emission_packets
from mhdlab.thermal_release import SurfaceInventory


def pool(**overrides):
    values = dict(name="manufactured sites", state_id="M:v0", attempt_frequency_s=1e8,
                  activation_energy_ev=.3, desorption_energy_ev=.2, internal_energy_ev=.05,
                  coverage_coefficient=.4, reference_coverage_m2=1e18, max_coverage_m2=1e18,
                  kinetic_provenance="manufactured coverage-dependent barrier",
                  energy_provenance="manufactured constant desorption/internal energy convention")
    values.update(overrides)
    return TemkinReservoir(**values)


def surface(temperatures=((450, 650), (700, 900)), times=(0, 1e-7), **metadata):
    return FlashSurfaceHistory(times, (0, .2, 1), temperatures, "synthetic", "temkin-test",
                               "manufactured prescribed temperature", "test boundary", **metadata)


def history(**kwargs):
    values = dict(surface=surface(), inventory=SurfaceInventory([[1e18, 5e17]]), reservoirs=[pool()],
                  parameter_provenance="manufactured history verification", subdivisions_per_interval=2)
    values.update(kwargs)
    return evolve_temkin_release(**values)


def test_zero_coefficient_has_analytic_first_order_limit_but_remains_a_distinct_type():
    sample = surface(temperatures=((600, 600), (600, 600)))
    p = pool(coverage_coefficient=0)
    result = history(surface=sample, reservoirs=[p], subdivisions_per_interval=4)
    hazard = p.attempt_frequency_s * 1e-7 * math.exp(-p.activation_energy_ev * QE / (KB * 600))
    initial = np.array([[1e18, 5e17]])
    np.testing.assert_allclose(result.coverage_m2[-1], initial * math.exp(-hazard), rtol=4e-15)
    np.testing.assert_allclose(result.released_m2.sum(axis=0), initial * (-math.expm1(-hazard)), rtol=5e-15)
    assert isinstance(result, TemkinReleaseHistory)
    assert not isinstance(result, SurfaceReleaseHistory)
    assert result.law_identity == LAW_IDENTITY
    assert np.all(result.inversion_iterations == 0)
    with pytest.raises(ValueError, match="SurfaceReleaseHistory"):
        sample_emission_packets(result, {p.state_id: 1e-26}, 0, np.full((3, 6), .5),
                                outward_normal_sign=1, current_length_m=1)


def test_local_number_energy_and_uniform_controls_remain_separate_per_reservoir():
    second = pool(name="second pool", state_id="M:v1", activation_energy_ev=.5,
                  coverage_coefficient=.2, desorption_energy_ev=.4)
    initial = np.array([[1e18, 5e17], [2e17, 8e17]])
    result = history(inventory=SurfaceInventory(initial), reservoirs=[pool(), second])
    np.testing.assert_allclose(result.coverage_m2[-1] + result.released_m2.sum(axis=0), initial, rtol=4e-15)
    np.testing.assert_allclose(result.translational_energy_j_m2,
                               2 * KB * result.temperature_k[:, None] * result.released_m2, rtol=4e-15)
    np.testing.assert_allclose(result.total_wall_energy_debit_j_m2,
                               result.translational_energy_j_m2 + result.internal_energy_j_m2 + result.binding_energy_j_m2,
                               rtol=4e-15)
    for name, values in result.matched_uniform().items():
        np.testing.assert_allclose(np.sum(values * np.diff(result.x_edges_m), axis=-1),
                                   result.width_integral(name), rtol=3e-15)
        np.testing.assert_array_equal(values[..., 0], values[..., 1])
    np.testing.assert_array_equal(initial, [[1e18, 5e17], [2e17, 8e17]])
    with pytest.raises(ValueError):
        result.coverage_m2.flags.writeable = True
    with pytest.raises(TypeError):
        result.surface_attributes["run_id"] = "changed"


def test_fractional_midpoint_handles_adjacent_representable_times():
    start, stop = 1., np.nextafter(1., np.inf)
    result = history(surface=surface(temperatures=((300, 400), (900, 1200)), times=(start, stop)),
                     subdivisions_per_interval=1)
    np.testing.assert_array_equal(result.temperature_k, [[600, 800]])
    with pytest.raises(ValueError, match="not representable"):
        history(surface=surface(times=(start, stop)), subdivisions_per_interval=2)


def test_midpoint_temporal_refinement_has_second_order_self_convergence():
    sample = surface(temperatures=((300, 350), (1000, 1100)))
    results = [history(surface=sample, subdivisions_per_interval=n).coverage_m2[-1] for n in (8, 16, 32)]
    coarse_difference = np.linalg.norm(results[0] - results[1])
    fine_difference = np.linalg.norm(results[1] - results[2])
    assert 3.6 < coarse_difference / fine_difference < 4.4


def test_entire_declared_coverage_domain_is_enforced():
    with pytest.raises(ValueError):
        history(inventory=SurfaceInventory([[1.01e18, 5e17]]))


@pytest.mark.parametrize("metadata", [
    {"temperature_location": "material_sample_plane", "material_sample_depth_m": 1e-6},
    {"current_averaging_interval_m": (0, .025)},
])
def test_unqualified_temperature_associations_cannot_drive_temkin(metadata):
    with pytest.raises(ValueError):
        history(surface=surface(**metadata))


def test_qualified_proxy_and_current_average_references_survive_history_creation():
    sample = surface(temperature_location="boundary_temperature_proxy", material_sample_depth_m=1e-6,
                     boundary_temperature_reference="manufactured normal-gradient check",
                     current_averaging_interval_m=(0, .025), current_uniformity_reference="manufactured invariant-current field")
    result = history(surface=sample)
    assert result.surface_attributes["material_sample_depth_m"] == 1e-6
    assert result.surface_attributes["boundary_temperature_reference"] == sample.boundary_temperature_reference
    assert result.surface_attributes["current_uniformity_reference"] == sample.current_uniformity_reference
    attributes = dict(result.surface_attributes)
    attributes.pop("current_averaging_interval_m_json")
    with pytest.raises(ValueError, match="missing required"):
        replace(result, surface_attributes=attributes)


@pytest.mark.parametrize("mutation", ["law", "time", "regime", "iterations", "kinetics", "energy"])
def test_mislabeled_or_inconsistent_interval_records_reject(mutation):
    result = history()
    if mutation == "law":
        arguments = {"law_identity": "constant_hazard_first_order"}
    elif mutation == "time":
        arguments = {"time_model": "constant_hazard_emission_time"}
    elif mutation == "regime":
        arguments = {"inversion_regime": np.zeros_like(result.inversion_regime)}
    elif mutation == "iterations":
        arguments = {"inversion_iterations": np.full_like(result.inversion_iterations, 129)}
    elif mutation == "kinetics":
        arguments = {"reservoirs": [replace(result.reservoirs[0], attempt_frequency_s=2e8)]}
    else:
        arguments = {"translational_energy_j_m2": result.translational_energy_j_m2 * 2}
    with pytest.raises(ValueError):
        replace(result, **arguments)


def test_equivalent_qualified_solve_diagnostics_and_last_bit_roots_are_portable():
    result = history(subdivisions_per_interval=1)
    iterations = result.inversion_iterations.copy()
    root = result.inversion_regime == 3
    assert np.any(root)
    iterations[root] += 1
    residual = result.inversion_log_time_residual.copy()
    residual[root] = 1e-12
    # A different qualified inversion can report a different iteration count
    # and slightly different release without changing the law or energy ledger.
    released = result.released_m2[0] * (1 + 1e-12)
    remaining = result.coverage_m2[0] - released
    budget = energy_budgets(released, result.reservoirs, result.temperature_k[0])
    balance = np.array([[math.fsum((float(remaining[p, x]), -float(result.coverage_m2[0, p, x]), float(released[p, x])))
                         for x in range(released.shape[1])] for p in range(released.shape[0])])
    equivalent = replace(result, coverage_m2=np.stack([result.coverage_m2[0], remaining]),
                         released_m2=released[None], inversion_iterations=iterations,
                         inversion_log_time_residual=residual, inventory_balance_residual_m2=balance[None],
                         **{name: value[None] for name, value in budget.items()})
    np.testing.assert_array_equal(equivalent.inversion_iterations, iterations)
    assert equivalent.law_identity == result.law_identity


@pytest.mark.parametrize("subdivisions", [0, -1, True, 1.5])
def test_invalid_temporal_partition_rejects(subdivisions):
    with pytest.raises(ValueError):
        history(subdivisions_per_interval=subdivisions)
