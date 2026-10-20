"""Local analytical depletion, patch-area ledgers and comparison controls."""

from dataclasses import replace
import math

import numpy as np
import pytest

from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PLANAR_SURFACE_SCHEMA, PlanarSurfaceHistory, evolve_planar_release
from mhdlab.thermal_release import BindingReservoir


KB_EV_K = 1.380649e-23 / 1.602176634e-19


def attributes(**overrides):
    return FlashSurfaceHistory([0.0, 1.0], [0.0, 1.0], [[300.0], [300.0]],
                               "synthetic", "patch-test", "analytical patch fixture", "test face",
                               **overrides).attributes


def surface(**overrides):
    values = dict(time_s=[0.0, .7], x_edges_m=[0.0, .2, 1.0], z_edges_m=[0.0, .3, 1.0],
                  temperature_k=np.full((2, 2, 2), 600.0), source_attributes=attributes())
    values.update(overrides)
    return PlanarSurfaceHistory(**values)


def pools():
    return (BindingReservoir("fast", "state-A", 2.0, 0.0, .5, .2),
            BindingReservoir("slow", "state-B", .5, 0.0, .8, .4))


def test_constant_temperature_local_depletion_matches_analytic_solution():
    sample = surface()
    coverage = 1e18 * np.array([[[1.0, 2.0], [3.0, 4.0]], [[2.0, 1.0], [4.0, 3.0]]])
    result = evolve_planar_release(sample, coverage, pools(), parameter_provenance="analytic constant hazard", subdivisions_per_interval=4)
    expected_remaining = coverage * np.exp(-np.array([2.0, .5])[:, None, None] * .7)
    expected_released = coverage * (-np.expm1(-np.array([2.0, .5])[:, None, None] * .7))
    np.testing.assert_allclose(result.coverage_m2[-1], expected_remaining, rtol=3e-15, atol=0.0)
    np.testing.assert_allclose(np.sum(result.released_m2, axis=0), expected_released, rtol=3e-15, atol=0.0)
    assert result.released_m2.shape == (4, 2, 2, 2)
    np.testing.assert_allclose(result.patch_areas_m2, [[.06, .24], [.14, .56]], rtol=2e-15, atol=0.0)
    expected_total = np.sum(expected_released * np.array([[.06, .24], [.14, .56]]), axis=(-2, -1))
    np.testing.assert_allclose(np.sum(result.extensive_totals(), axis=0), expected_total, rtol=3e-15, atol=0.0)
    np.testing.assert_allclose(result.extensive_totals("coverage_m2")[0] - result.extensive_totals("coverage_m2")[-1],
                               expected_total, rtol=3e-15, atol=0.0)
    np.testing.assert_array_equal(coverage, 1e18 * np.array([[[1.0, 2.0], [3.0, 4.0]], [[2.0, 1.0], [4.0, 3.0]]]))


def test_local_hot_cold_release_is_not_release_at_current_averaged_temperature():
    temperature = np.array([[[300.0], [1200.0]]] * 2)
    sample = surface(time_s=[0.0, 1e-6], x_edges_m=[0.0, 1.0], z_edges_m=[0.0, .5, 1.0], temperature_k=temperature)
    reservoir = BindingReservoir("specified-barrier", "state-A", 1e3, 1.0, 1.0)
    result = evolve_planar_release(sample, np.full((1, 2, 1), 1e18), (reservoir,), parameter_provenance="manufactured high-barrier comparison")
    rates = np.exp(-1.0 / (KB_EV_K * np.array([300.0, 1200.0])))
    expected_local = 1e18 * (-np.expm1(-1e3 * 1e-6 * rates))
    np.testing.assert_allclose(result.released_m2[0, 0, :, 0], expected_local, rtol=1e-13, atol=0.0)
    wrong_uniform_release = 1e18 * (-math.expm1(-1e3 * 1e-6 * math.exp(-1.0 / (KB_EV_K * 750.0))))
    correct_total = result.extensive_totals()[0, 0]
    assert 160.0 < correct_total / wrong_uniform_release < 170.0
    assert result.temperature_k[0, 0, 0] == 300.0
    assert result.temperature_k[0, 1, 0] == 1200.0


def test_nonuniform_patch_areas_and_energy_totals_match_direct_local_sum():
    sample = surface(temperature_k=np.array([[[300., 800.], [600., 1200.]]] * 2))
    result = evolve_planar_release(sample, np.full((2, 2, 2), 1e18), pools(), parameter_provenance="specified local energy test")
    areas = np.array([[.06, .24], [.14, .56]])
    for field in ("released_m2", "translational_energy_j_m2", "internal_energy_j_m2",
                  "binding_energy_j_m2", "total_wall_energy_debit_j_m2"):
        np.testing.assert_allclose(result.patch_extensive(field), getattr(result, field) * areas, rtol=2e-15, atol=0.0)
        np.testing.assert_allclose(result.extensive_totals(field), np.sum(getattr(result, field) * areas, axis=(-2, -1)), rtol=2e-15, atol=0.0)
    independent_translation = result.released_m2 * (2 * 1.380649e-23 * result.temperature_k[:, None, :, :])
    np.testing.assert_allclose(result.translational_energy_j_m2, independent_translation, rtol=2e-15, atol=0.0)


def test_width_null_preserves_each_current_patch_time_reservoir_and_energy():
    sample = surface(temperature_k=np.array([[[300., 900.], [700., 1200.]]] * 2))
    coverage = 1e18 * np.array([[[1., 5.], [2., 3.]], [[4., 1.], [3., 2.]]])
    result = evolve_planar_release(sample, coverage, pools(), parameter_provenance="controlled width redistribution", subdivisions_per_interval=2)
    null = result.matched_uniform_width()
    full = result.matched_uniform_full_face()
    for name, values in null.items():
        before = np.sum(getattr(result, name) * np.diff(sample.x_edges_m), axis=-1)
        after = np.sum(values * np.diff(sample.x_edges_m), axis=-1)
        np.testing.assert_allclose(after, before, rtol=3e-15, atol=0.0)
        np.testing.assert_array_equal(values[..., 0], values[..., 1])
        np.testing.assert_allclose(np.sum(full[name] * result.patch_areas_m2, axis=(-2, -1)), result.extensive_totals(name), rtol=3e-15, atol=0.0)
    assert not np.array_equal(null["released_m2"][:, :, 0, 0], null["released_m2"][:, :, 1, 0])
    np.testing.assert_array_equal(full["released_m2"][:, :, 0, 0], full["released_m2"][:, :, 1, 1])


def test_depth_proxy_and_unknown_metadata_are_preserved_without_auto_conversion():
    attrs = attributes(temperature_location="material_sample_plane", material_sample_depth_m=12.5e-6)
    attrs["extraction"] = {"dataset": "lattice", "selected_blocks": [2, 4]}
    sample = surface(source_attributes=attrs)
    attrs["extraction"]["selected_blocks"][0] = 99
    assert sample.source_attributes["extraction"]["selected_blocks"] == (2, 4)
    assert sample.source_attributes["material_sample_depth_m"] == 12.5e-6
    with pytest.raises(ValueError, match="not a boundary temperature"):
        evolve_planar_release(sample, np.ones((2, 2, 2)), pools(), parameter_provenance="specified")
    proxy = attributes(temperature_location="boundary_temperature_proxy", material_sample_depth_m=12.5e-6,
                       boundary_temperature_reference="independent depth convergence record")
    result = evolve_planar_release(surface(source_attributes=proxy), np.ones((2, 2, 2)), pools(), parameter_provenance="specified")
    assert result.surface_attributes["material_sample_depth_m"] == 12.5e-6
    assert result.surface_attributes["boundary_temperature_reference"] == "independent depth convergence record"
    assert result.surface_attributes["current_patch_intervals_m"] == ((0.0, .3), (.3, 1.0))
    assert result.surface_attributes["schema"] == PLANAR_SURFACE_SCHEMA


def test_preaveraged_current_trace_cannot_be_inflated_into_resolved_patches():
    attrs = attributes(current_averaging_interval_m=(0., 1.), current_uniformity_reference="a declared average remains an average")
    with pytest.raises(ValueError, match="cannot be inflated"):
        surface(source_attributes=attrs)


def test_fractional_midpoint_preserves_adjacent_representable_time_interval():
    start, stop = 1.0, np.nextafter(1.0, np.inf)
    sample = surface(time_s=[start, stop], temperature_k=np.array([[[1000.] * 2] * 2, [[2000.] * 2] * 2]))
    result = evolve_planar_release(sample, np.ones((2, 2, 2)), pools(), parameter_provenance="adjacent-time check")
    np.testing.assert_array_equal(result.temperature_k, 1500.0)
    with pytest.raises(ValueError, match="not representable"):
        evolve_planar_release(sample, np.ones((2, 2, 2)), pools(), parameter_provenance="adjacent-time check", subdivisions_per_interval=2)


def test_arrays_and_metadata_are_immutable_and_revalidated():
    sample = surface()
    with pytest.raises(ValueError):
        sample.temperature_k.flags.writeable = True
    with pytest.raises(TypeError):
        sample.source_attributes["run_id"] = "changed"
    result = evolve_planar_release(sample, np.ones((2, 2, 2)), pools(), parameter_provenance="immutable output")
    with pytest.raises(ValueError):
        result.released_m2.flags.writeable = True
    wrong_energy = result.translational_energy_j_m2 * 2
    with pytest.raises(ValueError, match="energy convention"):
        replace(result, translational_energy_j_m2=wrong_energy)


@pytest.mark.parametrize("changes", [
    {"time_s": [0., 0.]}, {"x_edges_m": [0., 0.]}, {"z_edges_m": [1., 0.]},
    {"temperature_k": np.ones((2, 2))}, {"temperature_k": np.zeros((2, 2, 2))},
    {"temperature_k": np.full((2, 2, 2), np.nan)}, {"temperature_k": np.ones((2, 2, 2), dtype=bool)},
    {"temperature_k": np.ones((2, 2, 2), dtype=complex)},
])
def test_invalid_geometry_and_temperature_rejected(changes):
    with pytest.raises(ValueError):
        surface(**changes)


@pytest.mark.parametrize("coverage", [np.ones((2, 2)), np.full((2, 2, 2), -1.),
                                      np.ones((2, 2, 2), dtype=bool), np.ones((2, 2, 2), dtype=complex)])
def test_invalid_local_inventory_rejected(coverage):
    with pytest.raises(ValueError):
        evolve_planar_release(surface(), coverage, pools(), parameter_provenance="specified")


def test_missing_or_conflicting_source_metadata_rejected():
    for key, value in (("temperature_location", None), ("material_sample_depth_units", "cm"),
                       ("coordinate_system", "cylindrical"), ("extra", np.inf),
                       ("source_kind", {}), ("schema", []), ("material_sample_depth_m", 10**1000)):
        attrs = attributes()
        attrs[key] = value
        with pytest.raises(ValueError):
            surface(source_attributes=attrs)
    attrs = attributes()
    attrs["current_patch_intervals_m"] = [[0., 1.]]
    with pytest.raises(ValueError, match="conflicts"):
        surface(source_attributes=attrs)
