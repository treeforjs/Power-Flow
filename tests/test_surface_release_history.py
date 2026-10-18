import math

import numpy as np
import pytest

from mhdlab.constants import KB, QE
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.surface_release_history import evolve_surface_release
from mhdlab.thermal_release import BindingReservoir, SurfaceInventory


def surface(temperatures=((500, 500), (500, 500)), times=(0, 1e-7), edges=(-0.0025, 0, 0.0025)):
    return FlashSurfaceHistory(times, edges, temperatures, "synthetic", "verification",
                               "manufactured history", "single stationary face")


def pool(name="A", barrier=0.3):
    return BindingReservoir(name, "parent-ground", 1e8, barrier, 0.2, 0.05)


def test_constant_temperature_analytic_depletion_and_integrated_energy():
    initial = np.array([[1e18, 2e18]])
    history = evolve_surface_release(surface(), SurfaceInventory(initial), [pool()],
                                     parameter_provenance="manufactured verification parameters",
                                     subdivisions_per_interval=7)
    hazard = 1e8 * 1e-7 * math.exp(-0.3 * QE / (KB * 500))
    expected = initial * math.exp(-hazard)
    np.testing.assert_allclose(history.coverage_m2[-1], expected, rtol=2e-15, atol=0)
    released = initial * (-math.expm1(-hazard))
    np.testing.assert_allclose(history.released_m2.sum(axis=0), released, rtol=8e-15, atol=0)
    np.testing.assert_allclose(history.translational_energy_j_m2.sum(axis=0),
                               2 * KB * 500 * released, rtol=8e-15, atol=0)
    np.testing.assert_allclose(history.total_wall_energy_debit_j_m2,
                               history.translational_energy_j_m2 + history.internal_energy_j_m2
                               + history.binding_energy_j_m2, rtol=2e-15, atol=0)
    np.testing.assert_array_equal(initial, [[1e18, 2e18]])
    with pytest.raises(ValueError):
        history.coverage_m2[0, 0, 0] = 0


def test_width_dependent_heating_and_uniform_null_preserve_each_reservoir_and_interval():
    history = evolve_surface_release(surface(((900, 400), (900, 400)), edges=(-0.0025, -0.001, 0.0025)),
                                     SurfaceInventory([[1e18, 1e18], [2e17, 2e17]]),
                                     [pool(), pool("B", 0.6)], parameter_provenance="manufactured",
                                     subdivisions_per_interval=4)
    assert np.all(history.released_m2[:, :, 0] > history.released_m2[:, :, 1])
    uniform = history.matched_uniform()
    widths = np.diff(history.x_edges_m)
    for name, values in uniform.items():
        np.testing.assert_allclose(values[:, :, 0], values[:, :, 1], rtol=0, atol=0)
        np.testing.assert_allclose((values * widths).sum(axis=-1), history.width_integral(name),
                                   rtol=3e-15, atol=0)


def test_midpoint_release_converges_for_linearly_varying_temperature():
    input_surface = surface(((300, 350), (1000, 800)))
    def final_coverage(n):
        return evolve_surface_release(input_surface, SurfaceInventory([[1e18, 1e18]]), [pool()],
                                      parameter_provenance="manufactured", subdivisions_per_interval=n).coverage_m2[-1]
    reference = final_coverage(1024)
    errors = [np.linalg.norm(final_coverage(n) - reference) for n in (8, 16, 32)]
    assert 3.7 < errors[0] / errors[1] < 4.3
    assert 3.7 < errors[1] / errors[2] < 4.3


def test_every_input_temperature_knot_is_preserved():
    history = evolve_surface_release(surface(((300, 400), (700, 600), (400, 500)), times=(2e-8, 3e-8, 8e-8)),
                                     SurfaceInventory([[1e18, 1e18]]), [pool()],
                                     parameter_provenance="manufactured", subdivisions_per_interval=3)
    assert len(history.time_edges_s) == 7
    assert history.time_edges_s[0] == 2e-8
    assert history.time_edges_s[3] == 3e-8
    assert history.time_edges_s[-1] == 8e-8
    np.testing.assert_allclose(history.average_flux_m2_s * np.diff(history.time_edges_s)[:, None, None],
                               history.released_m2, rtol=2e-15, atol=0)


@pytest.mark.parametrize("subdivisions", [0, -1, 1.5, True, np.inf, 100_001])
def test_invalid_subdivision_rejected(subdivisions):
    with pytest.raises(ValueError):
        evolve_surface_release(surface(), SurfaceInventory([[1e18, 1e18]]), [pool()],
                               parameter_provenance="manufactured", subdivisions_per_interval=subdivisions)


def test_missing_parameter_source_and_wrong_inventory_shape_are_rejected():
    with pytest.raises(ValueError, match="provenance"):
        evolve_surface_release(surface(), SurfaceInventory([[1e18, 1e18]]), [pool()], parameter_provenance=" ")
    with pytest.raises(ValueError):
        evolve_surface_release(surface(), SurfaceInventory([[1e18]]), [pool()], parameter_provenance="manufactured")


def test_midpoint_temperature_is_resolved_even_when_absolute_midpoint_time_rounds():
    start = 1e-7
    stop = np.nextafter(start, np.inf)
    result = evolve_surface_release(surface(((300, 400), (900, 1200)), times=(start, stop)),
                                    SurfaceInventory([[1e18, 1e18]]), [pool()],
                                    parameter_provenance="manufactured")
    np.testing.assert_array_equal(result.temperature_k, [[600, 800]])


@pytest.mark.parametrize("temperature,binding", [(300, 1e20), (1e-300, 0)])
def test_finite_scaled_energy_products_remain_valid_in_history(temperature, binding):
    result = evolve_surface_release(surface(((temperature, temperature), (temperature, temperature)), times=(0, 1)),
                                    SurfaceInventory([[1e300, 1e300]]),
                                    [BindingReservoir("extreme", "M", 1, 0, binding)],
                                    parameter_provenance="extreme manufactured numerical case")
    assert np.all(np.isfinite(result.total_wall_energy_debit_j_m2))
    assert np.all(result.translational_energy_j_m2 > 0)
