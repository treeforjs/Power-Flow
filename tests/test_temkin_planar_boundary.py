"""Local-area Temkin inventories and paired stationary 3D launches."""

import math

import numpy as np
import pytest

from mhdlab.constants import KB
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_boundary import PlanarEmissionPackets
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.temkin_boundary import sample_temkin_emission_packets
from mhdlab.temkin_emission_time import FrozenTemkinEmission, NoTemkinEmissionError, TIME_IDENTITY
from mhdlab.temkin_history import TemkinReleaseHistory, STEP_FIELDS, DIAGNOSTIC_FIELDS
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_planar_boundary import TemkinPlanarEmissionPackets, sample_temkin_planar_emission_packets
from mhdlab.temkin_release import LAW_IDENTITY, TemkinReservoir


MASSES = {"A:v0": 3e-26, "B:v1": 6e-26}


def history(*, z=(0, .01, .025), coverage=None, barrier=.3, raw=False, subdivisions=1):
    attributes = FlashSurfaceHistory([0, 1], [0, 1], [[600], [600]], "synthetic", "planar-temkin-test",
                                    "manufactured temperature patches", "test face",
                                    **({"temperature_location": "boundary_temperature_proxy", "material_sample_depth_m": 1e-5,
                                        "boundary_temperature_reference": "manufactured field has zero normal gradient"} if raw else {})).attributes
    x = (0, .001, .004)
    if raw:
        attributes["resolved_patch_extraction"] = {"patch_extraction_manifests": [{"sampling": {
            "face_normal_coordinate_m": .2, "outward_sign": 1, "material_sample_depth_m": 1e-5,
            "current_interval_m": [lo, hi], "x_edges_m": x,
        }} for lo, hi in zip(z[:-1], z[1:])]}
    temperature = np.linspace(450, 1000, (len(z) - 1) * 2).reshape(len(z) - 1, 2)
    surface = PlanarSurfaceHistory([0, 1e-7], x, z, [temperature, temperature], attributes)
    pools = [TemkinReservoir("A", "A:v0", 1e8, barrier, .5, .1, .4, 1e18, 1e18, "manufactured kinetics", "manufactured cost"),
             TemkinReservoir("B", "B:v1", 5e7, barrier, .7, .2, .3, 1e18, 1e18, "manufactured kinetics", "manufactured cost")]
    initial = np.full((2, len(z) - 1, 2), 1e18) if coverage is None else coverage
    return evolve_temkin_planar_release(surface, initial, pools, parameter_provenance="manufactured planar packet check",
                                        subdivisions_per_interval=subdivisions)


def draws(n):
    return np.linspace(.12, .88, n * 7).reshape(n, 7)


def sample(source=None, uniforms=None, **kwargs):
    parameters = dict(face_normal_coordinate_m=.2, outward_normal_sign=1)
    parameters.update(kwargs)
    return sample_temkin_planar_emission_packets(history() if source is None else source, MASSES, 0,
                                                 draws(13) if uniforms is None else uniforms, **parameters)


def test_each_local_inventory_uses_its_actual_area_once_and_keeps_thermal_energy_mixture():
    source = history(raw=True)
    result = sample(source)
    assert isinstance(result, TemkinPlanarEmissionPackets)
    assert not isinstance(result, PlanarEmissionPackets)
    assert result.law_identity == LAW_IDENTITY
    assert result.time_model == result.time_identity == TIME_IDENTITY
    assert result.source_history_schema == source.schema
    assert result.coverage_model == source.coverage_model
    expected = source.released_m2[0] * source.patch_areas_m2
    for component, (pool, current, width) in enumerate(zip(result.component_reservoir_index,
                                                          result.component_current_index, result.component_width_index)):
        selected = ((result.reservoir_index == pool) & (result.source_current_index == current)
                    & (result.source_width_index == width))
        assert np.count_nonzero(selected) == result.packet_counts[component]
        assert math.fsum(result.particle_weight[selected]) == pytest.approx(expected[pool, current, width], rel=4e-15)
        np.testing.assert_array_equal(result.temperature_k[selected], source.temperature_k[0, current, width])
        energy = expected[pool, current, width] * 2 * KB * source.temperature_k[0, current, width]
        assert result.component_expected_translational_energy_j[component] == pytest.approx(energy, rel=5e-15)
        realized = math.fsum(.5 * result.mass_kg[i] * np.dot(result.velocity_m_s[i], result.velocity_m_s[i])
                             * result.particle_weight[i] for i in np.flatnonzero(selected))
        assert result.component_sampled_translational_energy_j[component] == pytest.approx(realized, rel=5e-15)
    assert result.total_released_particles == pytest.approx(np.sum(expected), rel=4e-15)
    assert result.expected_translational_energy_j == pytest.approx(
        source.extensive_totals("translational_energy_j_m2").sum(), rel=4e-15)
    assert result.sampled_translational_energy_j != result.expected_translational_energy_j
    assert result.unequal_macro_weights
    assert result.current_length_m == .025


def test_pairing_keeps_current_normal_positions_time_state_velocity_and_weight():
    source, uniforms = history(), draws(13)
    result = sample(source, uniforms)
    np.testing.assert_array_equal(result.resolved_position_m[:, 1:], result.uniform_width_position_m[:, 1:])
    np.testing.assert_array_equal(result.resolved_position_m[:, 0], result.resolved_x_m)
    np.testing.assert_array_equal(result.uniform_width_position_m[:, 0], result.uniform_x_m)
    np.testing.assert_array_equal(result.resolved_position_m[:, 1], .2)
    assert np.any(result.resolved_x_m != result.uniform_x_m)
    for packet, (current, width, pool) in enumerate(zip(result.source_current_index, result.source_width_index, result.reservoir_index)):
        assert source.z_edges_m[current] < result.resolved_position_m[packet, 2] < source.z_edges_m[current + 1]
        assert source.x_edges_m[width] < result.resolved_x_m[packet] < source.x_edges_m[width + 1]
        assert result.state_ids[packet] == source.reservoirs[pool].state_id
    modified = uniforms.copy()
    modified[:, 6] = .95
    second = sample(source, modified)
    for name in ("time_s", "time_offset_s", "momentum_kg_m_s", "velocity_m_s", "particle_weight"):
        np.testing.assert_array_equal(getattr(result, name), getattr(second, name))
    np.testing.assert_array_equal(result.resolved_position_m[:, :2], second.resolved_position_m[:, :2])
    assert np.any(result.resolved_position_m[:, 2] != second.resolved_position_m[:, 2])


def test_row_allocation_is_exact_rational_with_minima_and_deterministic_ties():
    # Zero barrier makes row population proportional to area despite different T.
    source = history(barrier=0)
    result = sample(source, draws(13))
    np.testing.assert_array_equal(result.minimum_packets_per_current_row, [4, 4])
    np.testing.assert_array_equal(result.packets_per_current_row, [6, 7])
    tied = sample(history(z=(0, .0125, .025), barrier=0), draws(11))
    np.testing.assert_array_equal(tied.packets_per_current_row, [6, 5])
    initial = np.full((2, 3, 2), 1e18)
    initial[:, 1] = 0
    initial[1, 2, 0] = 0
    sparse = history(z=(0, .005, .01, .025), coverage=initial)
    seeded = sample(sparse, draws(7))
    np.testing.assert_array_equal(seeded.packets_per_current_row, [4, 0, 3])
    np.testing.assert_array_equal(seeded.packet_counts, np.ones(7, dtype=int))
    with pytest.raises(ValueError, match="at least 7"):
        sample(sparse, draws(6))


def test_rows_use_actual_temkin_sampler_and_retain_law_specific_conditional_times():
    source, uniforms = history(raw=True), draws(13)
    result = sample(source, uniforms)
    start = 0
    for slot, current in enumerate(result.row_current_indices):
        count = int(result.packets_per_current_row[current])
        row = TemkinReleaseHistory(
            time_edges_s=source.time_edges_s, x_edges_m=source.x_edges_m,
            temperature_k=source.temperature_k[:, current], coverage_m2=source.coverage_m2[:, :, current],
            **{name: getattr(source, name)[:, :, current] for name in (*STEP_FIELDS, *DIAGNOSTIC_FIELDS)},
            reservoirs=source.reservoirs, surface_attributes=result.row_surface_attributes[slot],
            parameter_provenance=source.parameter_provenance,
        )
        direct = sample_temkin_emission_packets(row, MASSES, 0, uniforms[start:start + count, :6],
                                                outward_normal_sign=1, current_length_m=float(np.diff(source.z_edges_m)[current]))
        for name in ("time_s", "time_offset_s", "velocity_m_s", "momentum_kg_m_s", "particle_weight"):
            np.testing.assert_array_equal(getattr(result, name)[start:start + count], getattr(direct, name))
        for packet in range(start, start + count):
            pool, width = result.reservoir_index[packet], result.source_width_index[packet]
            emission = FrozenTemkinEmission(source.reservoirs[pool], float(source.coverage_m2[0, pool, current, width]),
                                            float(source.temperature_k[0, current, width]), 1e-7)
            assert emission.cdf(float(result.time_offset_s[packet])) == pytest.approx(uniforms[packet, 1], rel=2e-9)
        adapter = result.row_surface_attributes[slot]["temkin_planar_row_adapter"]
        assert adapter["current_interval_m"] == tuple(source.z_edges_m[current:current + 2])
        assert "no whole-face uniformity" in adapter["qualification"]
        assert adapter["parent_surface_attributes_reference"] == "TemkinPlanarEmissionPackets.surface_attributes"
        start += count
    assert result.surface_attributes == source.surface_attributes
    assert isinstance(result.surface_attributes["resolved_patch_extraction"]["patch_extraction_manifests"], tuple)
    with pytest.raises(TypeError):
        result.row_surface_attributes[0]["temkin_planar_row_adapter"]["current_index"] = 99
    with pytest.raises(ValueError):
        result.resolved_position_m.flags.writeable = True


@pytest.mark.parametrize("arguments", [{"face_normal_coordinate_m": .3}, {"outward_normal_sign": -1}])
def test_launch_face_must_match_retained_raw_geometry(arguments):
    with pytest.raises(ValueError, match="conflicts with preserved raw"):
        sample(history(raw=True), **arguments)


def test_forged_law_arrays_metadata_and_empty_source_reject():
    source = history()
    object.__setattr__(source, "released_m2", source.released_m2 * 1.1)
    with pytest.raises(ValueError):
        sample(source)
    source = history(raw=True)
    attributes = dict(source.surface_attributes)
    attributes["material_sample_depth_m"] *= 2
    object.__setattr__(source, "surface_attributes", attributes)
    with pytest.raises(ValueError):
        sample(source)
    with pytest.raises(NoTemkinEmissionError):
        sample(history(coverage=np.zeros((2, 2, 2))))


def test_later_interval_keeps_actual_times_depleted_inventory_and_adapter_index():
    source = history(subdivisions=2)
    result = sample_temkin_planar_emission_packets(source, MASSES, 1, draws(13),
                                                   face_normal_coordinate_m=.2, outward_normal_sign=1)
    assert result.interval_index == 1
    assert result.time_start_s == 5e-8
    assert result.time_end_s == 1e-7
    assert result.total_released_particles == pytest.approx(source.extensive_totals()[1].sum(), rel=4e-15)
    for metadata in result.row_surface_attributes:
        assert metadata["temkin_planar_row_adapter"]["source_interval_index"] == 1


@pytest.mark.parametrize("columns", [6, 8])
def test_exactly_seven_uniform_columns_are_required(columns):
    with pytest.raises(ValueError, match="shape \\(n,7\\)"):
        sample(uniforms=np.full((13, columns), .5))


@pytest.mark.parametrize("column,value", [(0, 0), (1, 1), (6, 0), (6, 1), (6, np.nan)])
def test_uniforms_are_strict_open_seven_dimensional_variates(column, value):
    uniforms = draws(13)
    uniforms[0, column] = value
    with pytest.raises(ValueError):
        sample(uniforms=uniforms)


def test_masked_variates_missing_masses_and_invalid_geometry_reject():
    masked = np.ma.array(draws(13), mask=False)
    masked.mask[0, 6] = True
    with pytest.raises(ValueError, match="masked"):
        sample(uniforms=masked)
    with pytest.raises(ValueError, match="exact-state mass"):
        sample_temkin_planar_emission_packets(history(), {"A": 3e-26}, 0, draws(13),
                                              face_normal_coordinate_m=.2, outward_normal_sign=1)
    for arguments in ({"face_normal_coordinate_m": np.inf}, {"outward_normal_sign": True}):
        with pytest.raises(ValueError):
            sample(**arguments)
