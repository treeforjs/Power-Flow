"""Independent Q1 volume/surface, midplane gap and conservative-transfer limits."""
import numpy as np
import pytest

from mhdlab.solid_foil_geometry import (
    deform_hex_mesh, extrude_extensive_xy, midplane_gap,
    parallel_plate_shape_sensitivity, remap_extensive_xy,
)


def mesh(nx=2, ny=2, nz=2, center_y=0.):
    x, y, z = np.linspace(-.0025, .0025, nx + 1), center_y + np.linspace(-.0001, .0001, ny + 1), np.linspace(-.0125, .0125, nz + 1)
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    nodes = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    index = np.arange(len(nodes)).reshape((nz + 1, ny + 1, nx + 1))
    cells = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                cells.append([index[k, j, i], index[k, j, i + 1], index[k, j + 1, i + 1], index[k, j + 1, i],
                              index[k + 1, j, i], index[k + 1, j, i + 1], index[k + 1, j + 1, i + 1], index[k + 1, j + 1, i]])
    return nodes, np.array(cells)


def geometry(nodes, cells, displacement=None, **kwargs):
    arguments = dict(reference_density_kg_m3=7900., cell_temperature_k=np.full(len(cells), 300.),
                     admissible_temperature_range_k=[293.15, 773.15])
    arguments.update(kwargs)
    return deform_hex_mesh(nodes, cells, np.zeros_like(nodes) if displacement is None else displacement, **arguments)


def test_identity_volume_mass_boundary_area_and_face_ids():
    nodes, cells = mesh()
    result = geometry(nodes, cells)
    np.testing.assert_allclose(result.cell_volumes_m3, .005 * .0002 * .025 / 8, rtol=1e-14)
    np.testing.assert_allclose(result.gauss_jacobian_ratios, 1.)
    assert np.sum(result.cell_mass_kg) == pytest.approx(7900 * .005 * .0002 * .025)
    assert np.sum(result.triangle_areas_m2) == pytest.approx(2 * (.005 * .025 + .005 * .0002 + .0002 * .025))
    assert set(result.material_face_ids) == {"x_minus", "x_plus", "y_minus", "y_plus", "z_minus", "z_plus"}
    assert len(result.boundary_quads) == 24
    assert result.receipt()["production_field_handoff_authorized"] is False


def test_affine_volume_density_surface_and_translation_analytic_limit():
    nodes, cells = mesh(3, 2, 4)
    F = np.array([[1.01, .03, .001], [0., 1.02, 0.], [0., 0., .99]])
    translation = np.array([.01, .02, .03])
    result = geometry(nodes, cells, nodes @ (F - np.eye(3)).T + translation)
    determinant = np.linalg.det(F)
    np.testing.assert_allclose(result.cell_volumes_m3 / result.reference_cell_volumes_m3, determinant, rtol=2e-13)
    np.testing.assert_allclose(result.cell_density_kg_m3, 7900 / determinant, rtol=2e-13)
    np.testing.assert_allclose(result.cell_density_kg_m3 * result.cell_volumes_m3, result.cell_mass_kg, rtol=1e-15)
    assert result.controls["q1_surface_volume_relative_residual"] < 5e-12


def test_warped_q1_faces_preserve_boundary_vector_closure_and_volume():
    nodes, cells = mesh(3, 2, 4)
    x, y, z = nodes.T
    displacement = np.column_stack((4 * y * z, 2 * x * z, x * y))
    result = geometry(nodes, cells, displacement)
    np.testing.assert_allclose(np.sum(result.triangle_normals * result.triangle_areas_m2[:, :, None], axis=(0, 1)), 0., atol=1e-18)
    surface_volume = np.sum(result.quad_gauss_points_m * result.quad_gauss_area_vectors_m2) / 3
    assert surface_volume == pytest.approx(np.sum(result.cell_volumes_m3), rel=1e-12)
    assert np.min(result.gauss_jacobian_ratios) > 0


@pytest.mark.parametrize("factor", [0., -1., 1e-10])
def test_inverted_or_collapsed_geometry_rejected(factor):
    nodes, cells = mesh(1, 1, 1)
    displacement = np.zeros_like(nodes)
    displacement[:, 1] = (factor - 1) * nodes[:, 1]
    with pytest.raises(ValueError, match="Jacobian"):
        geometry(nodes, cells, displacement)


@pytest.mark.parametrize("temperature", [200., 774., np.nan])
def test_invalid_solid_temperature_is_never_clipped(temperature):
    nodes, cells = mesh(1, 1, 1)
    with pytest.raises(ValueError):
        geometry(nodes, cells, cell_temperature_k=[temperature])


def test_duplicate_cell_and_bad_connectivity_rejected():
    nodes, cells = mesh(1, 1, 1)
    with pytest.raises(ValueError, match="duplicate material"):
        geometry(nodes, np.vstack((cells, cells)))
    repeated = cells.copy()
    repeated[0, 1] = repeated[0, 0]
    with pytest.raises(ValueError, match="repeated"):
        geometry(nodes, repeated)


def test_uniform_own_center_expansion_has_analytic_midplane_gap():
    stretch = 1.01
    results = []
    for center in (-.0026, .0026):
        nodes, cells = mesh(3, 2, 3, center_y=center)
        displacement = (stretch - 1) * (nodes - [0., center, 0.])
        results.append(geometry(nodes, cells, displacement))
    gap = midplane_gap(*results)
    np.testing.assert_allclose(gap["gap_m"], .005 - .0002 * (stretch - 1), rtol=1e-14)
    assert gap["observation_z_m"] == 0.


def test_midplane_node_roundoff_signs_do_not_delete_intersection_segments():
    results = []
    for center in (-.0026, .0026):
        nodes, cells = mesh(8, 4, 4, center_y=center)
        displacement = np.zeros_like(nodes)
        midpoint = nodes[:, 2] == 0.
        # Opposite residual signs on nominally on-plane edges previously added
        # a false third intersection point and silently discarded triangles.
        displacement[midpoint, 2] = np.resize([-1e-16, 1e-16, 1e-16], midpoint.sum())
        results.append(geometry(nodes, cells, displacement))
    gap = midplane_gap(*results)
    np.testing.assert_allclose(gap["gap_m"], .005, rtol=1e-14)
    assert gap["x_m"][0] == pytest.approx(-.0025)
    assert gap["x_m"][-1] == pytest.approx(.0025)


def test_profile_does_not_extrapolate_across_a_real_hole_or_accept_fold():
    from mhdlab.solid_foil_geometry import _profile
    segments = np.array([[[0., 0.], [1., 1.]], [[1. + 1e-8, 1.], [2., 2.]]])
    with pytest.raises(ValueError, match="hole"):
        _profile(segments, [1. + 5e-9])
    overlapping = np.array([[[0., 0.], [1., 1.]], [[0., 1.], [1., 2.]]])
    with pytest.raises(ValueError, match="multiple y"):
        _profile(overlapping, [.5])


def test_nonuniform_gap_and_shape_sensitivity_are_conditional():
    lower_nodes, cells = mesh(4, 2, 4, center_y=-.0026)
    upper_nodes, _ = mesh(4, 2, 4, center_y=.0026)
    displacement = np.zeros_like(upper_nodes)
    displacement[:, 1] = -1e-6 * (1 - (upper_nodes[:, 0] / .0025) ** 2)
    gap = midplane_gap(geometry(lower_nodes, cells), geometry(upper_nodes, cells, displacement))
    assert gap["minimum_gap_m"] == pytest.approx(.005 - 1e-6)
    assert gap["maximum_gap_m"] == pytest.approx(.005)
    estimate = parallel_plate_shape_sensitivity(-110e3, gap["gap_m"])
    np.testing.assert_allclose(estimate["electric_field_magnitude_V_m"], 110e3 / gap["gap_m"])
    assert estimate["field_solve_performed"] is False


def test_contact_and_uncovered_observation_plane_rejected():
    nodes, cells = mesh(2, 2, 2)
    result = geometry(nodes, cells)
    with pytest.raises(ValueError, match="nonpositive"):
        midplane_gap(result, result)
    with pytest.raises(ValueError, match="observation plane"):
        midplane_gap(result, result, z_m=1.)


def test_tensor_product_energy_remap_and_axial_distribution_conserve_totals():
    # Source cell totals correspond to a constant volumetric/reference-area amount 6.
    sx, sy, tx, ty = [0., 1., 3.], [0., 2., 3.], [0., .5, 2., 3.], [0., 1., 3.]
    source = 6 * np.outer(np.diff(sy), np.diff(sx))
    remapped = remap_extensive_xy(sx, sy, tx, ty, source)
    np.testing.assert_allclose(remapped["values"], 6 * np.outer(np.diff(ty), np.diff(tx)))
    assert remapped["source_total"] == remapped["target_total"] == 54.
    extruded = extrude_extensive_xy(remapped["values"], [-.0125, -.01, 0., .0125])
    assert extruded.shape == (3, 2, 3)
    np.testing.assert_allclose(extruded.sum(axis=0), remapped["values"], rtol=1e-15)
    assert extruded.sum() == pytest.approx(54.)


def test_tensor_transfer_rejects_extrapolation_and_temperature_style_wrong_shape():
    with pytest.raises(ValueError, match="complete domains"):
        remap_extensive_xy([0, 1], [0, 1], [0, 1.1], [0, 1], [[1.]])
    with pytest.raises(ValueError, match="shape"):
        remap_extensive_xy([0, 1], [0, 1], [0, 1], [0, 1], [1.])


@pytest.mark.parametrize("gap", [0., -1., np.nan])
def test_parallel_plate_formula_never_continues_through_contact(gap):
    with pytest.raises(ValueError):
        parallel_plate_shape_sensitivity(1., gap)
