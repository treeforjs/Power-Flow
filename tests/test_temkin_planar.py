"""Local nonlinear depletion, physical areas and matched planar controls."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from mhdlab.constants import KB, QE
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_boundary import sample_planar_emission_packets
from mhdlab.planar_release import PlanarReleaseHistory, PlanarSurfaceHistory
from mhdlab.temkin_history import RELEASE_FIELDS
from mhdlab.temkin_planar import COVERAGE_MODEL, SCHEMA, TemkinPlanarReleaseHistory, evolve_temkin_planar_release
from mhdlab.temkin_release import LAW_IDENTITY, TemkinReservoir
from mhdlab.thermal_boundary import sample_emission_packets


def pool(**overrides):
    values = dict(name="manufactured sites", state_id="M:v0", attempt_frequency_s=1e8,
                  activation_energy_ev=.3, desorption_energy_ev=.2, internal_energy_ev=.05,
                  coverage_coefficient=.4, reference_coverage_m2=1e18, max_coverage_m2=1e18,
                  kinetic_provenance="manufactured independent-pool barrier",
                  energy_provenance="manufactured independent energy reference")
    values.update(overrides)
    return TemkinReservoir(**values)


def surface(*, times=(0., 1e-7), temperatures=None, nested=False, **metadata):
    base = FlashSurfaceHistory([0., 1.], [0., 1.], [[300.], [300.]], "synthetic", "temkin-planar-test",
                               "manufactured local patch temperatures", "test boundary", **metadata)
    attrs = base.attributes
    if nested:
        attrs["resolved_patch_extraction"] = {
            "schema": "mhdlab-flash-planar-extraction-v1", "other_source_evidence": ["retained"],
            "patch_extraction_manifests": [
                {"sampling": {"current_interval_m": interval, "x_edges_m": [0., .2, 1.],
                              "material_sample_depth_m": base.material_sample_depth_m},
                 "snapshots": [{"path": "manufactured-local.h5", "sha256": "synthetic-source-label"}]}
                for interval in ([0., .3], [.3, 1.])],
        }
    if temperatures is None:
        temperatures = np.array([[[450., 650.], [550., 900.]]] * len(times))
    return PlanarSurfaceHistory(times, [0., .2, 1.], [0., .3, 1.], temperatures, attrs)


def history(**overrides):
    values = dict(surface=surface(), coverage_m2=np.array([[[1e18, 5e17], [2e17, 8e17]]]),
                  reservoirs=[pool()], parameter_provenance="manufactured planar qualification",
                  subdivisions_per_interval=2)
    values.update(overrides)
    return evolve_temkin_planar_release(**values)


def test_nonlinear_patch_depletion_matches_independent_fractional_ode():
    sample = surface()
    initial = np.array([[[1e18, 5e17], [2e17, 8e17]]])
    reservoir = pool()
    result = history(surface=sample, coverage_m2=initial, subdivisions_per_interval=1)
    for z, x in np.ndindex(2, 2):
        temperature = sample.temperature_k[0, z, x]
        c0 = initial[0, z, x] / reservoir.reference_coverage_m2
        def rhs(_, fraction):
            barrier = reservoir.activation_energy_ev * (1 - reservoir.coverage_coefficient * fraction)
            return -reservoir.attempt_frequency_s * fraction * np.exp(-barrier * QE / (KB * temperature))
        numerical = solve_ivp(rhs, (0., 1e-7), [c0], method="DOP853", rtol=2e-12, atol=1e-14)
        assert numerical.success
        expected = numerical.y[0, -1] * reservoir.reference_coverage_m2
        np.testing.assert_allclose(result.coverage_m2[-1, 0, z, x], expected, rtol=3e-12)
    assert result.schema == SCHEMA
    assert result.law_identity == LAW_IDENTITY
    assert result.coverage_model == COVERAGE_MODEL
    assert isinstance(result, TemkinPlanarReleaseHistory)
    assert not isinstance(result, PlanarReleaseHistory)


def test_zero_coefficient_constant_hazard_limit_is_local_and_conservative():
    second = pool(name="second", state_id="M:v1", attempt_frequency_s=3e7, coverage_coefficient=0)
    first = pool(coverage_coefficient=0)
    initial = np.array([[[1e18, 5e17], [2e17, 8e17]], [[4e17, 3e17], [1e17, 6e17]]])
    result = history(coverage_m2=initial, reservoirs=[first, second], subdivisions_per_interval=4)
    for p, reservoir in enumerate((first, second)):
        hazard = reservoir.attempt_frequency_s * 1e-7 * np.exp(-reservoir.activation_energy_ev * QE / (KB * result.temperature_k[0]))
        np.testing.assert_allclose(result.coverage_m2[-1, p], initial[p] * np.exp(-hazard), rtol=5e-15)
        np.testing.assert_allclose(result.released_m2.sum(axis=0)[p], initial[p] * (-np.expm1(-hazard)), rtol=2e-14)
    np.testing.assert_allclose(result.patch_areas_m2, [[.06, .24], [.14, .56]], rtol=2e-15)
    np.testing.assert_allclose(result.coverage_m2[-1] + result.released_m2.sum(axis=0), initial, rtol=5e-15)
    assert np.all(result.inversion_regime == 1)
    np.testing.assert_allclose(result.average_flux_m2_s, result.released_m2 / np.diff(result.time_edges_s)[:, None, None, None])


def test_current_temperature_average_before_nonlinear_release_is_incorrect():
    local = np.array([[[300., 300.], [1200., 1200.]]] * 2)
    p = pool(activation_energy_ev=1., attempt_frequency_s=1e3, coverage_coefficient=.1)
    result = history(surface=surface(times=(0., 1e-6), temperatures=local), reservoirs=[p],
                     coverage_m2=np.full((1, 2, 2), 1e18), subdivisions_per_interval=1)
    average = np.sum(local[0] * result.patch_areas_m2) / np.sum(result.patch_areas_m2)
    wrong = history(surface=surface(times=(0., 1e-6), temperatures=np.full((2, 2, 2), average)),
                    reservoirs=[p], coverage_m2=np.full((1, 2, 2), 1e18), subdivisions_per_interval=1)
    assert result.extensive_totals()[0, 0] / wrong.extensive_totals()[0, 0] > 8
    assert result.released_m2[0, 0, 1, 0] > result.released_m2[0, 0, 0, 0] * 1e10
    np.testing.assert_array_equal(result.temperature_k[0], local[0])


def test_independent_pool_coverage_does_not_sum_into_a_shared_barrier():
    first = pool()
    second = pool(name="independent second", state_id="M:v1")
    initial = np.full((1, 2, 2), 8e17)
    one = history(reservoirs=[first], coverage_m2=initial)
    two = history(reservoirs=[first, second], coverage_m2=np.concatenate((initial, initial)))
    np.testing.assert_array_equal(one.released_m2[:, 0], two.released_m2[:, 0])
    np.testing.assert_array_equal(two.released_m2[:, 0], two.released_m2[:, 1])
    # The declared domain belongs to each pool; summed coverage exceeds Nmax.
    assert np.all(np.sum(two.coverage_m2[0], axis=0) > first.max_coverage_m2)


def test_width_and_full_face_controls_preserve_their_distinct_number_energy_totals():
    result = history()
    null = result.matched_uniform_width()
    full = result.matched_uniform_full_face()
    for name in RELEASE_FIELDS:
        source = getattr(result, name)
        np.testing.assert_allclose(np.sum(null[name] * np.diff(result.x_edges_m), axis=-1),
                                   np.sum(source * np.diff(result.x_edges_m), axis=-1), rtol=4e-15)
        np.testing.assert_array_equal(null[name][..., 0], null[name][..., 1])
        np.testing.assert_allclose(result.patch_extensive(name), source * result.patch_areas_m2, rtol=2e-15)
        np.testing.assert_allclose(np.sum(full[name] * result.patch_areas_m2, axis=(-2, -1)),
                                   result.extensive_totals(name), rtol=4e-15)
        np.testing.assert_array_equal(full[name][:, :, 0, 0], full[name][:, :, 1, 1])
    assert not np.array_equal(null["released_m2"][:, :, 0], null["released_m2"][:, :, 1])
    np.testing.assert_allclose(result.translational_energy_j_m2,
                               result.released_m2 * 2 * KB * result.temperature_k[:, None], rtol=3e-15)
    np.testing.assert_allclose(result.total_wall_energy_debit_j_m2,
                               result.translational_energy_j_m2 + result.binding_energy_j_m2 + result.internal_energy_j_m2,
                               rtol=3e-15)


def test_extreme_patch_fractions_do_not_discard_representable_uniform_moments():
    base = FlashSurfaceHistory([0., 1.], [0., 1.], [[300.], [300.]], "synthetic", "extreme-patch-test",
                               "manufactured floating-point scaling check", "test boundary")
    sample = PlanarSurfaceHistory([0., 1.], [0., 1e-200, 1e200], [0., 1.], np.full((2, 1, 2), 300.), base.attributes)
    p = pool(attempt_frequency_s=1., activation_energy_ev=0., coverage_coefficient=0.,
             reference_coverage_m2=1e300, max_coverage_m2=1e300)
    result = evolve_temkin_planar_release(sample, [[[1e300, 0.]]], [p], parameter_provenance="manufactured scaling check")
    for method in (result.matched_uniform_width, result.matched_uniform_full_face):
        for name, control in method().items():
            np.testing.assert_allclose(np.sum(control * result.patch_areas_m2, axis=(-2, -1)),
                                       result.extensive_totals(name), rtol=3e-15)
            assert np.all(control > 0)


@pytest.mark.parametrize("initial,temperature", [(1., 300.), (1e77, 1e23)])
def test_unresolved_positive_uniform_density_is_rejected_instead_of_losing_moments(initial, temperature):
    base = FlashSurfaceHistory([0., 1.], [0., 1.], [[300.], [300.]], "synthetic", "unresolved-patch-test",
                               "manufactured floating-point range check", "test boundary")
    sample = PlanarSurfaceHistory([0., 1.], [0., 1e-200, 1e200], [0., 1.], np.full((2, 1, 2), temperature), base.attributes)
    p = pool(attempt_frequency_s=1., activation_energy_ev=0., coverage_coefficient=0., desorption_energy_ev=0., internal_energy_ev=0.,
             reference_coverage_m2=max(1e18, initial), max_coverage_m2=max(1e18, initial))
    result = evolve_temkin_planar_release(sample, [[[initial, 0.]]], [p], parameter_provenance="manufactured range check")
    for method in (result.matched_uniform_width, result.matched_uniform_full_face):
        with pytest.raises(ValueError, match="below the qualified normal range"):
            method()


def test_proxy_depth_and_nested_raw_extraction_metadata_are_retained_immutably():
    sample = surface(nested=True, temperature_location="boundary_temperature_proxy", material_sample_depth_m=1e-5,
                     boundary_temperature_reference="manufactured normal-gradient convergence")
    result = history(surface=sample)
    assert result.surface_attributes["material_sample_depth_m"] == 1e-5
    assert result.surface_attributes["current_averaging_interval_m_json"] == "null"
    nested = result.surface_attributes["resolved_patch_extraction"]
    assert nested["patch_extraction_manifests"][1]["sampling"]["current_interval_m"] == (.3, 1.)
    assert nested["other_source_evidence"] == ("retained",)
    with pytest.raises(TypeError):
        nested["patch_extraction_manifests"][0]["sampling"]["current_interval_m"] = (0., 1.)
    for name in ("coverage_m2", "released_m2", "inversion_iterations", "inversion_regime"):
        with pytest.raises(ValueError):
            getattr(result, name).flags.writeable = True


def test_material_plane_and_whole_current_reduction_are_not_boundary_sources():
    with pytest.raises(ValueError, match="not a boundary temperature"):
        history(surface=surface(temperature_location="material_sample_plane", material_sample_depth_m=1e-5))
    with pytest.raises(ValueError, match="cannot be inflated"):
        surface(current_averaging_interval_m=(0., 1.), current_uniformity_reference="still averaged")


def test_fractional_midpoint_and_time_refinement_remain_qualified():
    start, stop = 1., np.nextafter(1., np.inf)
    temperatures = np.array([[[300., 400.], [500., 600.]], [[900., 1000.], [1100., 1200.]]])
    result = history(surface=surface(times=(start, stop), temperatures=temperatures), subdivisions_per_interval=1)
    np.testing.assert_array_equal(result.temperature_k, (temperatures[:1] + temperatures[1:]) / 2)
    with pytest.raises(ValueError, match="not representable"):
        history(surface=surface(times=(start, stop)), subdivisions_per_interval=2)
    sample = surface(temperatures=np.array([[[300., 350.], [400., 450.]], [[1000., 1050.], [1100., 1150.]]]))
    finals = [history(surface=sample, subdivisions_per_interval=n).coverage_m2[-1] for n in (8, 16, 32)]
    assert 3.6 < np.linalg.norm(finals[0] - finals[1]) / np.linalg.norm(finals[1] - finals[2]) < 4.4


@pytest.mark.parametrize("mutation", ["law", "time", "schema", "shared_sites", "kinetics", "energy", "domain", "regime", "false_tiny", "false_saturated", "diagnostic"])
def test_mislabeled_or_inconsistent_law_and_energy_records_are_rejected(mutation):
    original = history()
    changes = {
        "law": {"law_identity": "constant_hazard"}, "time": {"time_model": "constant_hazard_emission"},
        "schema": {"schema": "mhdlab-planar-release-v1"}, "shared_sites": {"coverage_model": "shared_site_density"},
        "kinetics": {"reservoirs": [replace(pool(), attempt_frequency_s=2e8)]},
        "energy": {"binding_energy_j_m2": original.binding_energy_j_m2 * 2},
        "domain": {"reservoirs": [replace(pool(), max_coverage_m2=9e17)]},
        "regime": {"inversion_regime": np.ones_like(original.inversion_regime)},
        "false_tiny": {"inversion_regime": np.full_like(original.inversion_regime, 2),
                       "inversion_iterations": np.zeros_like(original.inversion_iterations),
                       "inversion_log_time_residual": np.zeros_like(original.inversion_log_time_residual)},
        "false_saturated": {"inversion_regime": np.full_like(original.inversion_regime, 4),
                            "inversion_iterations": np.zeros_like(original.inversion_iterations),
                            "inversion_log_time_residual": np.zeros_like(original.inversion_log_time_residual)},
        "diagnostic": {"inversion_iterations": np.full_like(original.inversion_iterations, 129)},
    }
    with pytest.raises(ValueError):
        replace(original, **changes[mutation])


@pytest.mark.parametrize("key", ["current_patch_intervals_m", "patch_temperature_approximation", "nested_depth", "nested_interval", "nested_width"])
def test_source_geometry_and_nested_extraction_must_agree(key):
    original = history(surface=surface(nested=True))
    from mhdlab.planar_release_io import plain_metadata
    attrs = plain_metadata(original.surface_attributes)
    if key == "current_patch_intervals_m":
        attrs[key] = [[0., .5], [.5, 1.]]
    elif key == "patch_temperature_approximation":
        attrs[key] = "exact point values"
    else:
        sampling = attrs["resolved_patch_extraction"]["patch_extraction_manifests"][0]["sampling"]
        sampling[{"nested_depth": "material_sample_depth_m", "nested_interval": "current_interval_m", "nested_width": "x_edges_m"}[key]] = 1
    with pytest.raises(ValueError):
        replace(original, surface_attributes=attrs)


@pytest.mark.parametrize("changes", [
    {"coverage_m2": np.ones((1, 2))}, {"coverage_m2": np.full((1, 2, 2), -1.)},
    {"coverage_m2": np.ones((1, 2, 2), dtype=bool)}, {"coverage_m2": np.ones((1, 2, 2), dtype=complex)},
    {"coverage_m2": np.full((1, 2, 2), 1.1e18)}, {"subdivisions_per_interval": True},
    {"subdivisions_per_interval": 0}, {"subdivisions_per_interval": 100001},
    {"reservoirs": [pool(), pool()]}, {"reservoirs": None}, {"parameter_provenance": ""},
])
def test_invalid_inventory_partition_and_qualification_rejected(changes):
    with pytest.raises(ValueError):
        history(**changes)


def test_constant_hazard_packet_adapter_cannot_consume_temkin_planar_history():
    with pytest.raises(ValueError, match="SurfaceReleaseHistory"):
        sample_emission_packets(history(), {"M:v0": 1e-26}, 0, np.full((4, 6), .5),
                                outward_normal_sign=1, current_length_m=1.)
    with pytest.raises(ValueError, match="PlanarReleaseHistory"):
        sample_planar_emission_packets(history(), {"M:v0": 1e-26}, 0, np.full((4, 7), .5),
                                       face_normal_coordinate_m=0., outward_normal_sign=1)
    with pytest.raises(ValueError):
        history().extensive_totals("inversion_regime")
