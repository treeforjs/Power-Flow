"""Analytic conservation, surface, Piola, and rejection controls."""

from dataclasses import replace

import numpy as np
import pytest

from mhdlab.thermal_expansion_mapping import (
    affine_foil_geometry, conservative_remap_1d, projected_gap_y_m,
    qualify_field_handoff, transform_field,
)


def foil(F=np.eye(3), **kwargs):
    return affine_foil_geometry(F, reference_density_kg_m3=7900., **kwargs)


def test_overlap_density_remap_is_analytic_and_conserves_integral():
    transfer = conservative_remap_1d([0, 1, 3], [0, .5, 2, 3], [2., 4.])
    np.testing.assert_allclose(transfer.values, [2, 10 / 3, 4])
    assert transfer.source_total == 10
    assert transfer.target_total == pytest.approx(10)
    assert transfer.relative_conservation_error < 5e-13


def test_extensive_remap_does_not_treat_amount_as_density():
    transfer = conservative_remap_1d([0, 1, 3], [0, .5, 2, 3], [2., 8.], kind="extensive")
    np.testing.assert_allclose(transfer.values, [1, 5, 4])
    assert np.sum(transfer.values) == pytest.approx(10)


def test_vector_signed_remap_handles_cancellation_and_zeros():
    transfer = conservative_remap_1d([0, 1, 2], [0, 2], [[-1., 0., 2.], [1., 0., 4.]])
    np.testing.assert_array_equal(transfer.values, [[0., 0., 3.]])
    np.testing.assert_array_equal(transfer.source_total, [0., 0., 6.])


@pytest.mark.parametrize("target", [[-.1, 1, 3], [0, 1, 3.1], [.1, 1, 3], [0, 1, 2.9]])
def test_extrapolation_and_cropping_are_rejected(target):
    with pytest.raises(ValueError, match="complete domains"):
        conservative_remap_1d([0, 1, 3], target, [2., 4.])


@pytest.mark.parametrize("edges", [[0, 0, 1], [1, 0], [0, np.nan, 1], [0, np.inf], [0]])
def test_invalid_mesh_edges_rejected(edges):
    with pytest.raises(ValueError):
        conservative_remap_1d(edges, [0, 1], [1.])


def test_refined_nonuniform_grid_conserves_nonnegative_amount():
    # Exact cell averages of q=1+3x^2; arbitrary target partitions conserve its integral.
    source = np.linspace(0, 1, 65) ** 1.3
    target = np.linspace(0, 1, 98) ** .7
    density = np.diff(source + source ** 3) / np.diff(source)
    transfer = conservative_remap_1d(source, target, density)
    assert transfer.target_total == pytest.approx(2., rel=5e-14)
    assert np.min(transfer.values) >= 0
    assert transfer.input_sha256 == conservative_remap_1d(source, target, density).input_sha256
    with pytest.raises(ValueError):
        transfer.values[0] = 1


def test_uniform_expansion_geometry_mass_and_surface_scaling():
    stretch = 1.01
    reference, expanded = foil(), foil(stretch * np.eye(3))
    assert expanded.jacobian == pytest.approx(stretch ** 3)
    assert expanded.density_kg_m3 == pytest.approx(7900 / stretch ** 3)
    assert expanded.volume_m3 == pytest.approx(reference.volume_m3 * stretch ** 3)
    assert expanded.density_kg_m3 * expanded.volume_m3 == pytest.approx(reference.mass_kg)
    np.testing.assert_allclose(expanded.face_areas_m2, reference.face_areas_m2 * stretch ** 2)
    np.testing.assert_allclose(expanded.face_normals, reference.face_normals)
    np.testing.assert_allclose(expanded.vertices_m, reference.vertices_m * stretch)


def test_translation_has_no_mass_surface_or_normals_effect():
    reference = foil()
    translated = foil(translation_m=(.1, .2, -.3))
    np.testing.assert_allclose(translated.vertices_m, reference.vertices_m + [.1, .2, -.3])
    np.testing.assert_array_equal(translated.face_normals, reference.face_normals)
    np.testing.assert_array_equal(translated.face_areas_m2, reference.face_areas_m2)
    assert translated.mass_kg == reference.mass_kg


def test_sheared_faces_satisfy_nanson_and_closed_surface_identity():
    F = np.array([[1.02, .12, .03], [0., .99, .02], [0., 0., 1.01]])
    geometry = foil(F)
    np.testing.assert_allclose(np.sum(geometry.face_normals * geometry.face_areas_m2[:, None], axis=0), 0, atol=1e-20)
    # Divergence theorem independently checks surface areas, centers, normals and volume.
    surface_volume = np.sum(geometry.face_areas_m2 * np.sum(geometry.face_centers_m * geometry.face_normals, axis=1)) / 3
    assert surface_volume == pytest.approx(geometry.volume_m3, rel=1e-14)
    assert np.max(np.abs(np.linalg.norm(geometry.face_normals, axis=1) - 1)) < 1e-14


def test_two_centered_free_foils_have_analytic_gap_change():
    stretch = 1.01
    lower = foil(stretch * np.eye(3), center_reference_m=(0., -.0026, 0.))
    upper = foil(stretch * np.eye(3), center_reference_m=(0., .0026, 0.))
    assert projected_gap_y_m(lower, upper) == pytest.approx(.005 - .0002 * (stretch - 1))


def test_overlapping_pair_rejected():
    with pytest.raises(ValueError, match="positive separated"):
        projected_gap_y_m(foil(), foil())


@pytest.mark.parametrize("F", [np.diag([1, 1, -1]), np.diag([1, 1, 0]),
                               np.diag([1e-7, 1, 1]), np.diag([1e7, 1, 1]), np.ones((2, 2))])
def test_invalid_or_ill_conditioned_geometry_rejected(F):
    with pytest.raises(ValueError):
        foil(F)


def test_covector_preserves_line_integral_and_piola_preserves_flux():
    F = np.array([[1.1, .2, 0.], [0., .9, .1], [0., 0., 1.03]])
    E0, dX, B0, dA0 = np.array([1., 2., 3.]), np.array([4., -1., 2.]), np.array([.2, .3, .5]), np.array([1., 2., -3.])
    E = transform_field(E0, F, rule="covector")
    dx = transform_field(dX, F, rule="vector")
    B = transform_field(B0, F, rule="hdiv")
    dA = np.linalg.det(F) * np.linalg.solve(F.T, dA0)
    assert np.dot(E, dx) == pytest.approx(np.dot(E0, dX), rel=1e-14)
    assert np.dot(B, dA) == pytest.approx(np.dot(B0, dA0), rel=1e-14)


def handoff(**kwargs):
    arguments = dict(
        geometries=[foil(center_reference_m=(0., -.0026, 0.)), foil(center_reference_m=(0., .0026, 0.))],
        temperatures_k=[300., 300.], admissible_temperature_range_k=[293.15, 800.],
        transfers=[conservative_remap_1d([0, 1, 2], [0, .5, 2], [1., 2.])],
        mechanics_scope="uniform_affine_free_benchmark", support_description="free benchmark; no physical mount specified",
    )
    arguments.update(kwargs)
    return qualify_field_handoff(**arguments)


def test_handoff_is_benchmark_only_and_does_not_authorize_production():
    receipt = handoff()
    assert receipt["status"] == "qualified_affine_benchmark_interface_only"
    assert receipt["production_field_solve_authorized"] is False
    assert receipt["y_projected_gap_m"] == pytest.approx(.005)


@pytest.mark.parametrize("temperatures", [[290.], [801.], [np.inf], [300, 301]])
def test_handoff_blocks_invalid_or_nonuniform_temperature(temperatures):
    with pytest.raises(ValueError):
        handoff(temperatures_k=temperatures)


@pytest.mark.parametrize("override", [dict(mechanics_scope="nonuniform_transient"), dict(support_description=""),
                                      dict(transfers=[]), dict(transfers=[{"conservation_passed": True}])])
def test_handoff_requires_real_computed_controls(override):
    with pytest.raises(ValueError):
        handoff(**override)


def test_tampered_transfer_result_and_geometry_mass_rejected():
    transfer = conservative_remap_1d([0, 1, 2], [0, .5, 2], [1., 2.])
    with pytest.raises(ValueError, match="recomputed conservative"):
        handoff(transfers=[replace(transfer, values=np.zeros(2))])
    with pytest.raises(ValueError, match="recomputed map"):
        handoff(geometries=[replace(foil(), mass_kg=42.)])
