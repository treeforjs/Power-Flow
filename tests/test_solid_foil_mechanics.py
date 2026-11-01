"""Analytical patch, thermal, conservation and time/mesh limits for Q1 solids."""
import hashlib
import numpy as np
import pytest
from scipy.linalg import eigh

from mhdlab.thermal_expansion import SolidExpansionMaterial, load_ss304_material
from mhdlab.solid_foil_mechanics import (
    LOCAL_VERTEX_BITS, structured_hex_mesh, isotropic_constitutive_matrix,
    assemble_foil_operators, evaluate_foil_state, transient_foil,
)


def material(variable=False):
    return SolidExpansionMaterial(
        np.array([300., 800.]), np.array([0., .01]),
        np.array([300., 800.]), np.array([200e9, 160e9 if variable else 200e9]),
        8000., 300., "manufactured constant alpha20e-6/K, declared E")


def foil(nx=2, ny=1, nz=2):
    return structured_hex_mesh(np.linspace(-.0025, .0025, nx + 1),
                               np.linspace(.0025, .0027, ny + 1),
                               np.linspace(-.0125, .0125, nz + 1))


def kwargs(variable=False):
    return dict(material=material(variable), stress_free_temperature_K=300., poisson_ratio=.3)


def test_mesh_order_partition_and_immutable_geometry():
    mesh = foil(2, 2, 2)
    assert mesh.cell_shape == (2, 2, 2)
    assert mesh.nodes_m.shape == (27, 3)
    np.testing.assert_array_equal(mesh.cells[0], [0, 1, 4, 3, 9, 10, 13, 12])
    np.testing.assert_array_equal(LOCAL_VERTEX_BITS, [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                                                  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])
    assert mesh.reference_volumes_m3.sum() == pytest.approx(.005 * .0002 * .025)
    np.testing.assert_allclose(mesh.shape_gradients_m_inv.sum(axis=2), 0., atol=2e-12)
    with pytest.raises(ValueError):
        mesh.nodes_m.setflags(write=True)


def test_constitutive_tensor_and_engineering_shear():
    E, nu = 200e9, .3
    D = isotropic_constitutive_matrix(E, nu)
    strain = np.array([.001, -.0002, .0003, .0004, -.0005, .0006])
    tensor = np.array([[strain[0], strain[3] / 2, strain[5] / 2],
                       [strain[3] / 2, strain[1], strain[4] / 2],
                       [strain[5] / 2, strain[4] / 2, strain[2]]])
    stress = E / (1 + nu) * tensor + E * nu / ((1 + nu) * (1 - 2 * nu)) * np.trace(tensor) * np.eye(3)
    expected = stress[[0, 1, 2, 0, 1, 0], [0, 1, 2, 1, 2, 2]]
    np.testing.assert_allclose(D @ strain, expected, rtol=2e-15)
    assert (D @ strain) @ strain == pytest.approx(np.sum(stress * tensor), rel=2e-15)


def test_affine_patch_gauss_strain_stress_and_interior_equilibrium():
    mesh = foil(2, 2, 2)
    gradient = np.array([[.001, .0003, -.0002], [.0001, -.0004, .0002], [.0002, -.0001, .0006]])
    u = mesh.nodes_m @ gradient.T + [1e-6, -2e-6, 3e-6]
    T = np.full(len(mesh.cells), 300.)
    state = evaluate_foil_state(mesh, u, T, **kwargs())
    expected = [gradient[0, 0], gradient[1, 1], gradient[2, 2],
                gradient[0, 1] + gradient[1, 0], gradient[1, 2] + gradient[2, 1],
                gradient[0, 2] + gradient[2, 0]]
    np.testing.assert_allclose(state["gauss_strain"], np.broadcast_to(expected, (8, 8, 6)), atol=3e-17)
    ops = assemble_foil_operators(mesh, T, **kwargs())
    force = (ops["stiffness"] @ u.ravel()).reshape(-1, 3)
    interior = np.flatnonzero(np.all((mesh.nodes_m > mesh.nodes_m.min(axis=0)) &
                                    (mesh.nodes_m < mesh.nodes_m.max(axis=0)), axis=1))
    assert len(interior) == 1
    assert np.max(abs(force[interior])) / np.max(abs(force)) < 2e-13
    assert state["current_mass_kg"] == pytest.approx(state["reference_mass_kg"], rel=2e-15)


@pytest.mark.parametrize("motion", ["translation", "infinitesimal_rotation"])
def test_rigid_motion_has_no_linearized_strain_or_stress(motion):
    mesh = foil()
    if motion == "translation":
        u = np.broadcast_to([1e-6, 2e-6, -3e-6], mesh.nodes_m.shape)
    else:
        omega = np.array([1e-4, -2e-4, 3e-4])
        u = np.cross(omega, mesh.nodes_m - mesh.nodes_m.mean(axis=0))
    state = evaluate_foil_state(mesh, u, np.full(len(mesh.cells), 300.), **kwargs())
    np.testing.assert_allclose(state["gauss_strain"], 0., atol=3e-17)
    np.testing.assert_allclose(state["gauss_stress_Pa"], 0., atol=1e-5)


def test_uniform_free_thermal_affine_state_is_equilibrium_and_conserves_mass():
    mesh = foil(2, 2, 2)
    T = np.full(len(mesh.cells), 400.)
    eps = .002
    u = eps * (mesh.nodes_m - mesh.nodes_m.mean(axis=0))
    state = evaluate_foil_state(mesh, u, T, **kwargs())
    np.testing.assert_allclose(state["gauss_stress_Pa"], 0., atol=2e-5)
    ops = assemble_foil_operators(mesh, T, **kwargs())
    residual = ops["stiffness"] @ u.ravel() - ops["thermal_force_N"].ravel()
    assert np.max(abs(residual)) / np.max(abs(ops["thermal_force_N"])) < 1e-12
    np.testing.assert_allclose(state["gauss_jacobian"], (1 + eps)**3, rtol=2e-15)
    np.testing.assert_allclose(state["gauss_density_kg_m3"], 8000 / (1 + eps)**3, rtol=2e-15)
    assert ops["nodal_reference_mass_kg"].sum() == pytest.approx(8000 * .005 * .0002 * .025)


def test_fully_fixed_uniform_heating_has_hydrostatic_trial_stress_and_reactions():
    mesh = foil(2, 1, 1)  # Every node lies on an end plane: all DOFs fixed.
    temperatures = np.array([np.full(len(mesh.cells), 300.), np.full(len(mesh.cells), 400.)])
    result = transient_foil(mesh, [0., 1e-9], temperatures, boundary="both_end_clamped", **kwargs())
    np.testing.assert_array_equal(result["displacement_m"], 0.)
    expected = np.broadcast_to([-1e9, -1e9, -1e9, 0., 0., 0.], result["gauss_stress_Pa"][-1].shape)
    np.testing.assert_allclose(result["gauss_stress_Pa"][-1], expected, rtol=2e-15, atol=1e-6)
    np.testing.assert_allclose(result["reaction_N"].sum(axis=1), 0., atol=1e-9)
    assert np.max(abs(result["reaction_N"][-1])) > 0
    np.testing.assert_array_equal(result["stress_work_increment_J"], 0.)
    assert result["elastic_energy_J"][-1] > 0
    assert result["heat_debit_J"] is None


@pytest.mark.parametrize("boundary", ["free", "short_plane_clamped", "both_end_clamped"])
def test_zero_temperature_change_is_exact_rest(boundary):
    mesh = foil()
    result = transient_foil(mesh, [0., 1e-9, 2e-9], np.full((3, len(mesh.cells)), 300.), boundary=boundary, **kwargs())
    for key in ("displacement_m", "velocity_m_s", "acceleration_m_s2", "gauss_strain", "gauss_stress_Pa", "reaction_N"):
        np.testing.assert_array_equal(result[key], 0.)
    np.testing.assert_allclose(result["current_mass_kg"], result["reference_mass_kg"], rtol=2e-15)
    if boundary == "short_plane_clamped":
        np.testing.assert_array_equal(mesh.nodes_m[result["constrained_node_indices"], 2], -.0125)


def test_3d_rigid_translation_is_preserved_by_free_transient_solver():
    mesh = foil()
    times = np.linspace(0., 2e-9, 11)
    velocity = np.broadcast_to([1., -2., 3.], mesh.nodes_m.shape)
    result = transient_foil(mesh, times, np.full((len(times), len(mesh.cells)), 300.),
                            boundary="free", initial_velocity_m_s=velocity, **kwargs())
    np.testing.assert_allclose(result["displacement_m"], times[:, None, None] * velocity, rtol=2e-11, atol=2e-20)
    np.testing.assert_allclose(result["velocity_m_s"], np.broadcast_to(velocity, result["velocity_m_s"].shape), rtol=2e-11)


def test_nonuniform_heating_variable_E_closes_discrete_stress_work_not_heat():
    mesh = foil(2, 2, 2)
    times = np.linspace(0., 1e-8, 21)
    spatial = np.linspace(.2, 1., len(mesh.cells))
    temperatures = 300 + 80 * np.linspace(0., 1., len(times))[:, None] * spatial
    result = transient_foil(mesh, times, temperatures, boundary="free", **kwargs(variable=True))
    scale = max(float(np.max(result["kinetic_energy_J"])), float(np.max(abs(result["stress_work_increment_J"]))))
    assert np.max(abs(result["mechanical_balance_residual_J"])) / scale < 2e-10
    momentum = np.einsum("tni,n->ti", result["velocity_m_s"], result["nodal_reference_mass_kg"])
    momentum_scale = np.sum(abs(result["velocity_m_s"]) * result["nodal_reference_mass_kg"][None, :, None])
    assert np.max(abs(momentum)) / momentum_scale < 2e-12
    assert np.max(abs(result["gauss_stress_Pa"])) > 0
    assert np.max(abs(result["displacement_m"])) > 0
    assert result["physical_coupled_ready"] is False
    # Changing constitutive parameters adds elastic energy without an invented heat debit.
    assert np.max(abs(np.diff(result["elastic_energy_J"]) - result["stress_work_increment_J"])) > 0


def test_unforced_3d_newmark_energy_and_time_refinement():
    # A cube separates the time integrator check from conditioning of a single
    # extreme-aspect-ratio foil element. Target geometry has separate patch,
    # source-heating, rigid-mode and continuum-wave checks above/below.
    mesh = structured_hex_mesh([0., .001], [0., .001], [0., .001])
    T = np.full(len(mesh.cells), 300.)
    ops = assemble_foil_operators(mesh, T, **kwargs())
    mass = np.repeat(ops["nodal_reference_mass_kg"], 3)
    eigenvalues, eigenvectors = eigh(ops["stiffness"].toarray(), np.diag(mass))
    index = np.flatnonzero(eigenvalues > eigenvalues.max() * 1e-8)[0]
    omega = np.sqrt(eigenvalues[index])
    mode = eigenvectors[:, index]
    mode = 1e-8 * mode / np.max(abs(mode))
    errors = []
    for count in (11, 21, 41):
        times = np.linspace(0., .7 / omega, count)
        result = transient_foil(mesh, times, np.full((count, len(mesh.cells)), 300.),
                                boundary="free", initial_displacement_m=mode.reshape(-1, 3), **kwargs())
        exact = np.cos(omega * times[-1]) * mode
        errors.append(np.max(abs(result["displacement_m"][-1].ravel() - exact)))
        energy = result["kinetic_energy_J"] + result["elastic_energy_J"]
        np.testing.assert_allclose(energy, energy[0], rtol=2e-10)
    assert min(np.array(errors[:-1]) / errors[1:]) > 3.95


def test_3d_mesh_refinement_against_free_axial_continuum_wave_with_nu_zero():
    length = .005
    speed = np.sqrt(200e9 / 8000)
    times = np.linspace(0., .37 * length / speed, 201)
    errors = []
    for n in (2, 4, 8):
        mesh = foil(n, 1, 1)
        u0 = np.zeros_like(mesh.nodes_m)
        u0[:, 0] = 1e-8 * np.cos(np.pi * (mesh.nodes_m[:, 0] + length / 2) / length)
        result = transient_foil(mesh, times, np.full((len(times), len(mesh.cells)), 300.),
                                material=material(), stress_free_temperature_K=300., poisson_ratio=0.,
                                boundary="free", initial_displacement_m=u0)
        exact = np.cos(np.pi * speed * times[-1] / length) * u0
        errors.append(np.max(abs(result["displacement_m"][-1] - exact)))
    assert min(np.array(errors[:-1]) / errors[1:]) > 3.7


def test_temperature_shape_allows_native_zyx_order_and_rejects_hot_cell():
    mesh = foil(2, 2, 2)
    temperatures = np.full((2,) + mesh.cell_shape, 300.)
    result = transient_foil(mesh, [0., 1e-9], temperatures, boundary="free", **kwargs())
    assert result["cell_temperature_K"].shape == (2, 8)
    temperatures[1, 0, 0, 0] = 1000.
    with pytest.raises(ValueError, match="no clipping"):
        transient_foil(mesh, [0., 1e-9], temperatures, boundary="free", material=load_ss304_material(),
                       stress_free_temperature_K=300., poisson_ratio=.3)


def test_transient_result_binds_geometry_not_just_cell_count():
    first = foil(1, 1, 1)
    second = structured_hex_mesh([-.002, .002], first.y_edges_m, first.z_edges_m)
    hashes = []
    for mesh in (first, second):
        result = transient_foil(mesh, [0., 1e-9], np.full((2, 1), 300.), boundary="free", **kwargs())
        expected = hashlib.sha256(np.asarray(mesh.nodes_m, dtype="<f8").tobytes() +
                                  np.asarray(mesh.cells, dtype="<i8").tobytes()).hexdigest()
        assert result["mesh_geometry_sha256"] == expected
        hashes.append(result["mesh_geometry_sha256"])
    assert len(first.cells) == len(second.cells)
    assert hashes[0] != hashes[1]


@pytest.mark.parametrize("nu", [-1., .5, np.nan, np.inf])
def test_invalid_poisson_ratio_rejected(nu):
    with pytest.raises(ValueError):
        isotropic_constitutive_matrix(200e9, nu)


def test_large_initial_strain_and_mismatched_clamp_state_rejected():
    mesh = foil()
    T = np.full((2, len(mesh.cells)), 300.)
    with pytest.raises(ValueError, match="2 percent"):
        transient_foil(mesh, [0., 1e-9], T, boundary="free", initial_displacement_m=.03 * mesh.nodes_m, **kwargs())
    with pytest.raises(ValueError, match="homogeneous constraints"):
        transient_foil(mesh, [0., 1e-9], T, boundary="short_plane_clamped", initial_velocity_m_s=np.ones_like(mesh.nodes_m), **kwargs())


def test_invalid_time_and_nonfinite_states_rejected():
    mesh = foil()
    T = np.full((2, len(mesh.cells)), 300.)
    with pytest.raises(ValueError):
        transient_foil(mesh, [0., 0.], T, boundary="free", **kwargs())
    with pytest.raises(ValueError):
        transient_foil(mesh, [0., 1e-9], T, boundary="free", initial_displacement_m=np.full_like(mesh.nodes_m, np.nan), **kwargs())


@pytest.mark.parametrize("dt", [1e-200, 1e-160, 1e200])
def test_extreme_timestep_arithmetic_rejected(dt):
    mesh = foil()
    with pytest.raises(ValueError, match="Newmark arithmetic range"):
        transient_foil(mesh, [0., dt], np.full((2, len(mesh.cells)), 300.), boundary="free", **kwargs())
