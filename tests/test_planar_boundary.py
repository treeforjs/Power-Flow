"""Local source inventories and paired width controls at a planar launch face."""

from dataclasses import replace
import math

import numpy as np
import pytest

from mhdlab.constants import KB
from mhdlab.flash_surface import SCHEMA as FLASH_SCHEMA, FlashSurfaceHistory
from mhdlab.planar_boundary import sample_planar_emission_packets
from mhdlab.planar_release import PlanarSurfaceHistory, evolve_planar_release
from mhdlab.surface_release_history import SurfaceReleaseHistory
from mhdlab.thermal_boundary import NoThermalEmissionError, sample_emission_packets
from mhdlab.thermal_release import BindingReservoir


MASSES = {"A:v0": 1e-26, "B:v1": 3e-26}


def history(*, z_edges=(0, .25, 1), coverage=None, raw_geometry=False):
    nz = len(z_edges) - 1
    attributes = FlashSurfaceHistory([0, 1], [0, 1], [[300], [300]], "synthetic",
                                    "planar-packet-test", "manufactured boundary fixture", "test face",
                                    **({"temperature_location": "boundary_temperature_proxy",
                                        "material_sample_depth_m": 1e-5,
                                        "boundary_temperature_reference": "manufactured field has no normal gradient"}
                                       if raw_geometry else {})).attributes
    if raw_geometry:
        attributes["resolved_patch_extraction"] = {
            "schema": "mhdlab-flash-planar-extraction-v1",
            "patch_extraction_manifests": [{"sampling": {
                "face_normal_coordinate_m": .2, "outward_sign": 1, "material_sample_depth_m": 1e-5,
                "current_interval_m": [lo, hi],
            }} for lo, hi in zip(z_edges[:-1], z_edges[1:])],
        }
    temperature = np.linspace(450, 900, nz * 2).reshape(nz, 2)
    surface = PlanarSurfaceHistory([0, .5], [0, .2, 1], z_edges,
                                  np.stack([temperature] * 2), attributes)
    pools = [BindingReservoir("fast", "A:v0", 2, 0, .5, .2),
             BindingReservoir("slow", "B:v1", .5, 0, .8, .4)]
    initial = np.full((2, nz, 2), 1e18) if coverage is None else coverage
    return evolve_planar_release(surface, initial, pools, parameter_provenance="manufactured packet verification")


def draws(n):
    return np.linspace(.11, .89, n * 7).reshape(n, 7)


def packets(source=None, uniforms=None, **kwargs):
    arguments = dict(face_normal_coordinate_m=.2, outward_normal_sign=-1)
    arguments.update(kwargs)
    return sample_planar_emission_packets(source or history(), MASSES, 0,
                                          draws(11) if uniforms is None else uniforms, **arguments)


def test_each_local_component_conserves_number_and_transfers_thermal_energy_diagnostics():
    source = history()
    result = packets(source)
    assert result.required_minimum_n == 8
    # Each row starts with four mandatory components. Three extra packets
    # split by row area 1:3 using largest remainder, producing row totals 5:6.
    np.testing.assert_array_equal(result.minimum_packets_per_current_row, [4, 4])
    np.testing.assert_array_equal(result.packets_per_current_row, [5, 6])
    assert result.unequal_macro_weights
    expected_particles = source.released_m2[0] * source.patch_areas_m2
    for component, (pool, current, width) in enumerate(zip(result.component_reservoir_index,
                                                          result.component_current_index,
                                                          result.component_width_index)):
        selected = ((result.reservoir_index == pool) & (result.source_current_index == current)
                    & (result.source_width_index == width))
        assert np.count_nonzero(selected) == result.packet_counts[component]
        assert math.fsum(result.particle_weight[selected]) == pytest.approx(expected_particles[pool, current, width], rel=3e-15)
        assert result.component_released_particles[component] == pytest.approx(expected_particles[pool, current, width], rel=3e-15)
        expected_energy = 2 * KB * source.temperature_k[0, current, width] * expected_particles[pool, current, width]
        assert result.component_expected_translational_energy_j[component] == pytest.approx(expected_energy, rel=4e-15)
        realized = math.fsum(.5 * result.mass_kg[i] * float(np.dot(result.velocity_m_s[i], result.velocity_m_s[i]))
                             * result.particle_weight[i] for i in np.flatnonzero(selected))
        assert result.component_sampled_translational_energy_j[component] == pytest.approx(realized, rel=4e-15)
        assert result.component_translational_energy_sampling_residual_j[component] == pytest.approx(realized - expected_energy, rel=6e-15)
    assert result.total_released_particles == pytest.approx(np.sum(expected_particles), rel=3e-15)
    assert result.expected_translational_energy_j == pytest.approx(
        np.sum(source.translational_energy_j_m2[0] * source.patch_areas_m2), rel=3e-15)
    assert result.sampled_translational_energy_j != result.expected_translational_energy_j
    assert result.translational_energy_sampling_residual_j == pytest.approx(
        result.sampled_translational_energy_j - result.expected_translational_energy_j)


def test_comparison_changes_only_width_position_and_seventh_variate_only_changes_current_position():
    source = history()
    uniforms = draws(11)
    result = packets(source, uniforms)
    np.testing.assert_array_equal(result.resolved_position_m[:, 1:], result.uniform_width_position_m[:, 1:])
    np.testing.assert_array_equal(result.resolved_position_m[:, 1], .2)
    assert np.any(result.resolved_position_m[:, 0] != result.uniform_width_position_m[:, 0])
    assert np.all(result.velocity_m_s[:, 1] < 0)
    np.testing.assert_array_equal(result.momentum_kg_m_s, result.mass_kg[:, None] * result.velocity_m_s)
    for i, (pool, current, width) in enumerate(zip(result.reservoir_index, result.source_current_index, result.source_width_index)):
        assert source.z_edges_m[current] < result.resolved_position_m[i, 2] < source.z_edges_m[current + 1]
        assert source.x_edges_m[width] < result.resolved_position_m[i, 0] < source.x_edges_m[width + 1]
        assert result.state_ids[i] == source.reservoirs[pool].state_id
        assert result.temperature_k[i] == source.temperature_k[0, current, width]
    changed = uniforms.copy()
    changed[:, 6] = .9
    second = packets(source, changed)
    np.testing.assert_array_equal(second.momentum_kg_m_s, result.momentum_kg_m_s)
    np.testing.assert_array_equal(second.time_s, result.time_s)
    np.testing.assert_array_equal(second.particle_weight, result.particle_weight)
    np.testing.assert_array_equal(second.resolved_position_m[:, :2], result.resolved_position_m[:, :2])
    assert np.any(second.resolved_position_m[:, 2] != result.resolved_position_m[:, 2])


def test_row_sampler_outputs_are_used_without_velocity_time_or_weight_changes():
    source = history(raw_geometry=True)
    uniforms = draws(11)
    result = packets(source, uniforms, outward_normal_sign=1)
    start = 0
    for slot, current in enumerate(result.row_current_indices):
        count = int(result.packets_per_current_row[current])
        row = SurfaceReleaseHistory(
            time_edges_s=source.time_edges_s, x_edges_m=source.x_edges_m,
            temperature_k=source.temperature_k[:, current], coverage_m2=source.coverage_m2[:, :, current],
            **{name: getattr(source, name)[:, :, current] for name in (
                "released_m2", "translational_energy_j_m2", "internal_energy_j_m2", "binding_energy_j_m2",
                "total_wall_energy_debit_j_m2", "inventory_balance_residual_m2")},
            reservoirs=source.reservoirs, surface_attributes=result.row_surface_attributes[slot],
            parameter_provenance=source.parameter_provenance,
        )
        underlying = sample_emission_packets(row, MASSES, 0, uniforms[start:start + count, :6],
                                             outward_normal_sign=1,
                                             current_length_m=float(np.diff(source.z_edges_m)[current]))
        for name in ("velocity_m_s", "momentum_kg_m_s", "time_s", "time_offset_s", "particle_weight"):
            np.testing.assert_array_equal(getattr(result, name)[start:start + count], getattr(underlying, name))
        attributes = result.row_surface_attributes[slot]
        assert attributes["schema"] == FLASH_SCHEMA
        assert attributes["planar_row_adapter"]["current_interval_m"] == tuple(source.z_edges_m[current:current + 2])
        assert attributes["planar_row_adapter"]["parent_surface_attributes"] is source.surface_attributes
        start += count
    assert result.surface_attributes is source.surface_attributes
    assert result.surface_attributes["temperature_location"] == "boundary_temperature_proxy"
    assert result.surface_attributes["material_sample_depth_m"] == 1e-5
    with pytest.raises(TypeError):
        result.surface_attributes["run_id"] = "changed"
    with pytest.raises(TypeError):
        result.row_surface_attributes[0]["planar_row_adapter"]["current_index"] = 99
    with pytest.raises(ValueError):
        result.resolved_position_m.flags.writeable = True


def test_zero_rows_and_components_receive_no_packets_and_all_positive_components_are_seeded():
    coverage = np.full((2, 3, 2), 1e18)
    coverage[:, 1] = 0
    coverage[1, 2, 0] = 0
    source = history(z_edges=(0, .2, .5, 1), coverage=coverage)
    result = packets(source, draws(7))
    np.testing.assert_array_equal(result.packets_per_current_row, [4, 0, 3])
    np.testing.assert_array_equal(result.packet_counts, np.ones(7, dtype=int))
    assert 1 not in result.source_current_index
    assert result.row_current_indices == (0, 2)
    with pytest.raises(ValueError, match="at least 7"):
        packets(source, draws(6))


def test_equal_population_largest_remainder_tie_uses_current_index_without_an_outer_draw():
    source = history(z_edges=(0, .5, 1))
    first = packets(source)
    np.testing.assert_array_equal(first.packets_per_current_row, [6, 5])
    uniforms = draws(11)
    uniforms[:, 0] = .99
    second = packets(source, uniforms)
    np.testing.assert_array_equal(second.packets_per_current_row, first.packets_per_current_row)


def test_empty_source_and_missing_exact_state_masses_reject():
    empty = history(coverage=np.zeros((2, 2, 2)))
    with pytest.raises(NoThermalEmissionError):
        packets(empty)
    with pytest.raises(ValueError, match="exact-state mass"):
        sample_planar_emission_packets(history(), {"A": 1e-26, "B": 3e-26}, 0, draws(11),
                                      face_normal_coordinate_m=.2, outward_normal_sign=1)


@pytest.mark.parametrize("column,value", [(0, 0), (1, 1), (6, 0), (6, 1), (6, np.nan)])
def test_uniforms_are_strictly_open_in_all_columns(column, value):
    uniforms = draws(11)
    uniforms[0, column] = value
    with pytest.raises(ValueError):
        packets(uniforms=uniforms)


@pytest.mark.parametrize("arguments", [
    {"face_normal_coordinate_m": np.inf}, {"face_normal_coordinate_m": True},
    {"face_normal_coordinate_m": 10 ** 1000},
    {"outward_normal_sign": 0}, {"outward_normal_sign": True},
])
def test_launch_geometry_requires_explicit_valid_values(arguments):
    with pytest.raises(ValueError):
        packets(**arguments)


@pytest.mark.parametrize("columns", [6, 8])
def test_exactly_seven_uniform_columns_are_required(columns):
    with pytest.raises(ValueError, match="shape \\(n,7\\)"):
        packets(uniforms=np.full((11, columns), .5))


@pytest.mark.parametrize("arguments", [{"face_normal_coordinate_m": .3}, {"outward_normal_sign": -1}])
def test_launch_geometry_must_match_retained_raw_face(arguments):
    values = dict(face_normal_coordinate_m=.2, outward_normal_sign=1)
    values.update(arguments)
    with pytest.raises(ValueError, match="conflicts with preserved raw"):
        packets(history(raw_geometry=True), **values)


def test_retained_raw_current_patch_association_cannot_be_changed():
    source = history(raw_geometry=True)
    attrs = dict(source.surface_attributes)
    attrs["resolved_patch_extraction"] = {"patch_extraction_manifests": [
        {"sampling": {"face_normal_coordinate_m": .2, "outward_sign": 1,
                      "material_sample_depth_m": 1e-5, "current_interval_m": [0, .9]}}
    ] * 2}
    malformed = replace(source, surface_attributes=attrs)
    with pytest.raises(ValueError, match="resolved current interval"):
        packets(malformed, outward_normal_sign=1)
