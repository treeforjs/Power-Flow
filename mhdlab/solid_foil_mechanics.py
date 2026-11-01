"""Three-dimensional, small-strain Q1 solid mechanics for declared benchmarks.

Each solve owns one disconnected rectangular foil. The feed/short assembly is
not invented: clamped end planes are explicit brackets, not measured supports.
Temperature-dependent E and integrated isotropic eigenstrain come from the
bounded material record. Poisson ratio is always an explicit caller assumption.
No magnetic force, plasticity, melting, heat debit or caloric EOS is included.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from itertools import product
import warnings

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import MatrixRankWarning, spsolve


VOIGT_COMPONENTS = ("xx", "yy", "zz", "xy", "yz", "xz")
LOCAL_VERTEX_BITS = np.array([
    [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
    [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
], dtype=np.int64)
MAX_PRINCIPAL_SMALL_STRAIN = .02


def _finite(value, name):
    result = np.asarray(value, dtype=float)
    if not result.size or not np.isfinite(result).all():
        raise ValueError(f"{name} must be nonempty and finite")
    return result


def _readonly(value):
    result = np.ascontiguousarray(value)
    return np.frombuffer(result.tobytes(), dtype=result.dtype).reshape(result.shape)


def _edges(value, name):
    result = _finite(value, name)
    if result.ndim != 1 or len(result) < 2 or np.any(np.diff(result) <= 0):
        raise ValueError(f"{name} must be strictly increasing edges")
    return result


def isotropic_constitutive_matrix(young_modulus_Pa, poisson_ratio):
    """Isotropic D for engineering strains [xx,yy,zz,2exy,2eyz,2exz].

    Stress entries are physical shear stresses, so ``stress @ strain`` is the
    tensor double contraction without extra factors. E may be an array; D then
    has E.shape + (6,6). No grade-qualified Poisson data are silently supplied.
    """
    modulus = _finite(young_modulus_Pa, "Young modulus")
    nu = float(poisson_ratio)
    if np.any(modulus <= 0) or not np.isfinite(nu) or not -1 < nu < .5:
        raise ValueError("positive E and explicit Poisson ratio in (-1,.5) required")
    mu = 1 / (2 * (1 + nu))
    lam = nu / ((1 + nu) * (1 - 2 * nu))
    unit = np.zeros((6, 6))
    unit[:3, :3] = lam
    unit[np.arange(3), np.arange(3)] += 2 * mu
    unit[3:, 3:] = np.eye(3) * mu
    result = modulus[..., None, None] * unit
    return _finite(result, "constitutive matrix")


@dataclass(frozen=True)
class StructuredHexMesh:
    """Reference mesh: global nodes/cells in (z,y,x) order, x fastest.

    Local vertices are 000,100,110,010,001,101,111,011 in xyz bits.
    All Gauss quantities use the same eight quadrature-point ordering.
    """
    x_edges_m: np.ndarray
    y_edges_m: np.ndarray
    z_edges_m: np.ndarray
    nodes_m: np.ndarray
    cells: np.ndarray
    cell_shape: tuple[int, int, int]
    dof_indices: np.ndarray
    shape_gradients_m_inv: np.ndarray  # [cell,gauss,node,xyz]
    B_m_inv: np.ndarray  # [cell,gauss,Voigt,element_dof]
    gauss_weights_m3: np.ndarray
    reference_volumes_m3: np.ndarray


def mesh_geometry_sha256(mesh):
    """Bind reference xyz coordinates and canonical connectivity byte-for-byte."""
    digest = hashlib.sha256()
    digest.update(np.asarray(mesh.nodes_m, dtype="<f8").tobytes(order="C"))
    digest.update(np.asarray(mesh.cells, dtype="<i8").tobytes(order="C"))
    return digest.hexdigest()


def structured_hex_mesh(x_edges_m, y_edges_m, z_edges_m):
    """Construct a nonuniform Cartesian Q1 mesh and cache exact geometry."""
    x, y, z = (_edges(value, name) for value, name in (
        (x_edges_m, "x edges"), (y_edges_m, "y edges"), (z_edges_m, "z edges")))
    nx, ny, nz = len(x) - 1, len(y) - 1, len(z) - 1
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    nodes = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    cells, dimensions = [], []
    for k, j, i in product(range(nz), range(ny), range(nx)):
        bits = LOCAL_VERTEX_BITS + [i, j, k]
        cells.append((bits[:, 2] * (ny + 1) + bits[:, 1]) * (nx + 1) + bits[:, 0])
        dimensions.append([x[i + 1] - x[i], y[j + 1] - y[j], z[k + 1] - z[k]])
    cells = np.array(cells, dtype=np.int64)
    dimensions = np.array(dimensions)
    gauss = np.array([[gx, gy, gz] for gz, gy, gx in product((-1 / np.sqrt(3), 1 / np.sqrt(3)), repeat=3)])
    signs = 2 * LOCAL_VERTEX_BITS - 1
    natural_gradient = np.empty((8, 8, 3))
    for d in range(3):
        other = [p for p in range(3) if p != d]
        natural_gradient[:, :, d] = signs[None, :, d] * np.prod(
            1 + gauss[:, None, other] * signs[None, :, other], axis=-1) / 8
    gradients = natural_gradient[None, ...] * (2 / dimensions)[:, None, None, :]
    volumes = np.prod(dimensions, axis=1)
    weights = np.broadcast_to(volumes[:, None] / 8, (len(cells), 8)).copy()
    _finite(gradients, "shape gradients")
    if not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError("reference cell volumes must be finite and positive")
    B = np.zeros((len(cells), 8, 6, 24))
    for a in range(8):
        dx, dy, dz = (gradients[:, :, a, d] for d in range(3))
        B[:, :, 0, 3 * a] = dx
        B[:, :, 1, 3 * a + 1] = dy
        B[:, :, 2, 3 * a + 2] = dz
        B[:, :, 3, 3 * a] = dy
        B[:, :, 3, 3 * a + 1] = dx
        B[:, :, 4, 3 * a + 1] = dz
        B[:, :, 4, 3 * a + 2] = dy
        B[:, :, 5, 3 * a] = dz
        B[:, :, 5, 3 * a + 2] = dx
    dofs = (cells[..., None] * 3 + np.arange(3)).reshape(-1, 24)
    return StructuredHexMesh(*[_readonly(a) for a in (x, y, z, nodes, cells)],
                             (nz, ny, nx), _readonly(dofs), _readonly(gradients),
                             _readonly(B), _readonly(weights), _readonly(volumes))


def _cell_values(mesh, value, name):
    result = _finite(value, name)
    if result.shape == mesh.cell_shape:
        result = result.reshape(-1)
    if result.shape != (len(mesh.cells),):
        raise ValueError(f"{name} must have shape (ncells,) or (nz,ny,nx)")
    return result


def _material_values(mesh, temperature, material, stress_free_temperature_K):
    temperature = _cell_values(mesh, temperature, "cell temperatures")
    material.check_temperature(temperature)
    eps = material.strain(temperature, stress_free_temperature_K=stress_free_temperature_K)
    modulus = material.modulus(temperature)
    rho = material.density_at_stress_free_temperature(stress_free_temperature_K)
    return temperature, eps, modulus, rho


class _AssemblyCache:
    def __init__(self, mesh, poisson_ratio):
        self.mesh = mesh
        self.D_unit = isotropic_constitutive_matrix(1., poisson_ratio)
        self.ke_unit = np.einsum("eqai,ab,eqbj,eq->eij", mesh.B_m_inv, self.D_unit,
                                 mesh.B_m_inv, mesh.gauss_weights_m3, optimize=True)
        self.fe_unit = np.einsum("eqai,a,eq->ei", mesh.B_m_inv,
                                 self.D_unit @ [1., 1., 1., 0., 0., 0.],
                                 mesh.gauss_weights_m3, optimize=True)
        self.rows = np.broadcast_to(mesh.dof_indices[:, :, None], self.ke_unit.shape).ravel()
        self.cols = np.broadcast_to(mesh.dof_indices[:, None, :], self.ke_unit.shape).ravel()

    def operators(self, modulus, eps):
        ndof = self.mesh.nodes_m.size
        values = self.ke_unit * modulus[:, None, None]
        K = coo_matrix((values.ravel(), (self.rows, self.cols)), shape=(ndof, ndof)).tocsr()
        force = np.zeros(ndof)
        np.add.at(force, self.mesh.dof_indices.ravel(),
                  (self.fe_unit * (modulus * eps)[:, None]).ravel())
        _finite(K.data, "stiffness")
        _finite(force, "thermal force")
        return K, force


def assemble_foil_operators(mesh, cell_temperature_K, *, material,
                            stress_free_temperature_K, poisson_ratio):
    """Assemble 2x2x2 Gauss K/f_th and row-sum lumped reference mass.

    Public independent evaluation; the transient solver caches geometric
    matrices internally and only reweights temperature-dependent properties.
    Matrices use interleaved xyz DOFs; nodal forces have shape (nnodes,3).
    """
    _, eps, modulus, rho = _material_values(mesh, cell_temperature_K, material, stress_free_temperature_K)
    K, force = _AssemblyCache(mesh, poisson_ratio).operators(modulus, eps)
    mass = np.zeros(len(mesh.nodes_m))
    np.add.at(mass, mesh.cells.ravel(), np.repeat(rho * mesh.reference_volumes_m3 / 8, 8))
    return {"stiffness": K, "thermal_force_N": force.reshape(-1, 3),
            "nodal_reference_mass_kg": mass, "cell_modulus_Pa": modulus,
            "cell_thermal_strain": eps, "reference_density_kg_m3": rho}


def _evaluate(mesh, displacement, eps, modulus, rho, D_unit):
    u = _finite(displacement, "displacement")
    if u.shape != mesh.nodes_m.shape:
        raise ValueError("displacement must have shape (nnodes,3)")
    local = u[mesh.cells]
    strain = np.einsum("eqai,ei->eqa", mesh.B_m_inv, local.reshape(-1, 24), optimize=True)
    gradient = np.einsum("eai,eqaj->eqij", local, mesh.shape_gradients_m_inv, optimize=True)
    symmetric = .5 * (gradient + np.swapaxes(gradient, -1, -2))
    _finite(strain, "Gauss strain")
    if np.max(np.abs(np.linalg.eigvalsh(symmetric))) > MAX_PRINCIPAL_SMALL_STRAIN:
        raise ValueError("2 percent principal small-strain benchmark guard exceeded")
    F = gradient + np.eye(3)
    J = np.linalg.det(F)
    if not np.isfinite(J).all() or np.any(J <= 0):
        raise ValueError("invalid deformed Gauss Jacobian")
    elastic = strain - eps[:, None, None] * np.array([1., 1., 1., 0., 0., 0.])
    stress = np.einsum("ab,eqb,e->eqa", D_unit, elastic, modulus, optimize=True)
    gauss_energy = .5 * np.sum(stress * elastic, axis=-1) * mesh.gauss_weights_m3
    volume = np.sum(J * mesh.gauss_weights_m3, axis=1)
    density = rho / J
    _finite(stress, "Gauss stress")
    _finite(gauss_energy, "elastic energy")
    _finite(density, "deformed density")
    return {"gauss_strain": strain, "gauss_stress_Pa": stress,
            "gauss_elastic_strain": elastic, "gauss_elastic_energy_J": gauss_energy,
            "elastic_energy_J": float(gauss_energy.sum()), "gauss_deformation_gradient": F,
            "gauss_jacobian": J, "gauss_density_kg_m3": density,
            "cell_volume_m3": volume, "cell_density_kg_m3": rho * mesh.reference_volumes_m3 / volume,
            "reference_mass_kg": float(rho * mesh.reference_volumes_m3.sum()),
            "current_mass_kg": float(np.sum(density * J * mesh.gauss_weights_m3))}


def evaluate_foil_state(mesh, displacement_m, cell_temperature_K, *, material,
                        stress_free_temperature_K, poisson_ratio):
    """Evaluate linearized strains/stresses and total geometry mass at Gauss points.

    Mass uses the total geometric Jacobian once. Constitutive energy uses the
    reference integration volume; this is a small-strain elastic trial, not a
    finite-strain objective material or a phase/yield qualification.
    """
    _, eps, modulus, rho = _material_values(mesh, cell_temperature_K, material, stress_free_temperature_K)
    return _evaluate(mesh, displacement_m, eps, modulus, rho,
                     isotropic_constitutive_matrix(1., poisson_ratio))


def _constraint_nodes(mesh, boundary):
    if boundary == "free":
        return np.array([], dtype=np.int64)
    z = mesh.nodes_m[:, 2]
    if boundary == "short_plane_clamped":
        return np.flatnonzero(z == mesh.z_edges_m[0])
    if boundary == "both_end_clamped":
        return np.flatnonzero((z == mesh.z_edges_m[0]) | (z == mesh.z_edges_m[-1]))
    raise ValueError("boundary must explicitly be free, short_plane_clamped or both_end_clamped")


def transient_foil(mesh, times_s, cell_temperature_K, *, material,
                   stress_free_temperature_K, poisson_ratio, boundary,
                   initial_displacement_m=None, initial_velocity_m_s=None):
    """Implicit average-acceleration Newmark for one Q1 solid foil.

    T is [time,cell] or [time,z,y,x]; values are uniform within each element.
    Homogeneous displacement constraints are benchmark end-plane brackets.
    Per-step stress work uses trapezoidal stress contracted with engineering
    strain increment and reference Gauss weights. Its negative equals kinetic
    energy change for this unforced, immobile-support discrete system. Changing
    E/eigenstrain changes elastic energy separately; it is NOT a heat debit.
    """
    times = _finite(times_s, "times")
    if times.ndim != 1 or len(times) < 2 or times[0] < 0 or np.any(np.diff(times) <= 0):
        raise ValueError("finite increasing nonnegative times required")
    temp = _finite(cell_temperature_K, "temperature history")
    if temp.shape == (len(times),) + mesh.cell_shape:
        temp = temp.reshape(len(times), -1)
    if temp.shape != (len(times), len(mesh.cells)):
        raise ValueError("temperature history must have shape (nt,ncells) or (nt,nz,ny,nx)")
    material.check_temperature(temp)
    eps = material.strain(temp, stress_free_temperature_K=stress_free_temperature_K)
    modulus = material.modulus(temp)
    rho = material.density_at_stress_free_temperature(stress_free_temperature_K)
    cache = _AssemblyCache(mesh, poisson_ratio)
    mass_nodes = np.zeros(len(mesh.nodes_m))
    np.add.at(mass_nodes, mesh.cells.ravel(), np.repeat(rho * mesh.reference_volumes_m3 / 8, 8))
    if not np.isfinite(mass_nodes).all() or np.any(mass_nodes <= 0):
        raise ValueError("positive finite reference nodal masses required")
    mass = np.repeat(mass_nodes, 3)
    fixed_nodes = _constraint_nodes(mesh, boundary)
    fixed = (fixed_nodes[:, None] * 3 + np.arange(3)).ravel()
    active = np.setdiff1d(np.arange(len(mass)), fixed)
    u = np.zeros((len(times), len(mass)))
    v = np.zeros_like(u)
    a = np.zeros_like(u)
    for target, values, name in ((u, initial_displacement_m, "initial displacement"),
                                 (v, initial_velocity_m_s, "initial velocity")):
        if values is not None:
            values = _finite(values, name)
            if values.shape != mesh.nodes_m.shape or np.any(values.ravel()[fixed] != 0):
                raise ValueError("initial state must have nodal xyz shape and obey homogeneous constraints")
            target[0] = values.ravel()
    K, force = cache.operators(modulus[0], eps[0])
    a[0, active] = ((force - K @ u[0]) / mass)[active]
    _finite(a[0], "initial acceleration")
    reactions = np.zeros_like(u)
    reactions[0, fixed] = (K @ u[0] - force)[fixed]
    residual_max = np.zeros(len(times))
    residual_max[0] = np.max(np.abs((mass * a[0] + K @ u[0] - force)[active]), initial=0)
    states = [_evaluate(mesh, u[0].reshape(-1, 3), eps[0], modulus[0], rho, cache.D_unit)]
    for step, dt in enumerate(np.diff(times), 1):
        dt = float(dt)
        dt_squared = dt * dt
        if (not np.isfinite(dt_squared) or dt_squared <= 0 or
                dt_squared < 4 / np.finfo(float).max):
            raise ValueError("time interval outside finite Newmark arithmetic range")
        coefficient = 4 / dt_squared
        K, force = cache.operators(modulus[step], eps[step])
        predictor = u[step - 1] + dt * v[step - 1] + .25 * dt_squared * a[step - 1]
        rhs = force + coefficient * mass * predictor
        A = K + diags(coefficient * mass)
        _finite(A.data, "effective stiffness")
        _finite(rhs, "Newmark right hand side")
        if len(active):
            with warnings.catch_warnings():
                warnings.simplefilter("error", MatrixRankWarning)
                u[step, active] = spsolve(A[active][:, active].tocsc(), rhs[active], permc_spec="MMD_AT_PLUS_A")
        a[step] = coefficient * (u[step] - predictor)
        v[step] = v[step - 1] + .5 * dt * (a[step - 1] + a[step])
        a[step, fixed] = 0
        v[step, fixed] = 0
        _finite(u[step], "Newmark displacement")
        _finite(v[step], "Newmark velocity")
        _finite(a[step], "Newmark acceleration")
        residual = mass * a[step] + K @ u[step] - force
        reactions[step, fixed] = residual[fixed]
        residual_max[step] = np.max(np.abs(residual[active]), initial=0)
        states.append(_evaluate(mesh, u[step].reshape(-1, 3), eps[step], modulus[step], rho, cache.D_unit))
    result = {key: np.array([state[key] for state in states]) for key in states[0]}
    kinetic = .5 * np.sum(mass[None, :] * v * v, axis=1)
    strain, stress = result["gauss_strain"], result["gauss_stress_Pa"]
    stress_work = np.sum(.5 * (stress[1:] + stress[:-1]) * np.diff(strain, axis=0), axis=-1) * mesh.gauss_weights_m3
    _finite(kinetic, "kinetic energy")
    _finite(stress_work, "stress work")
    result.update({
        "times_s": times.copy(), "displacement_m": u.reshape(len(times), -1, 3),
        "mesh_geometry_sha256": mesh_geometry_sha256(mesh),
        "velocity_m_s": v.reshape(len(times), -1, 3), "acceleration_m_s2": a.reshape(len(times), -1, 3),
        "reaction_N": reactions.reshape(len(times), -1, 3), "nodal_reference_mass_kg": mass_nodes,
        "reference_mass_kg": float(rho * mesh.reference_volumes_m3.sum()),
        "reference_density_kg_m3": float(rho),
        "cell_reference_mass_kg": rho * mesh.reference_volumes_m3,
        "cell_reference_volumes_m3": mesh.reference_volumes_m3,
        "cell_temperature_K": temp.copy(), "cell_modulus_Pa": modulus, "cell_thermal_strain": eps,
        "poisson_ratio": float(poisson_ratio), "poisson_ratio_basis": "explicit benchmark assumption",
        "material_provenance": material.provenance,
        "stress_free_temperature_K": float(stress_free_temperature_K),
        "gauss_weights_m3": mesh.gauss_weights_m3, "voigt_components": VOIGT_COMPONENTS,
        "strain_shear_convention": "engineering shear; stress uses physical shear components",
        "kinetic_energy_J": kinetic, "gauss_stress_work_increment_J": stress_work,
        "stress_work_increment_J": stress_work.sum(axis=(1, 2)),
        "mechanical_balance_residual_J": np.diff(kinetic) + stress_work.sum(axis=(1, 2)),
        "free_dof_equilibrium_residual_max_N": residual_max,
        "constrained_node_indices": fixed_nodes, "boundary": boundary,
        "scope": "one_disconnected_3d_Q1_small_strain_elastic_trial",
        "support_basis": "declared benchmark bracket; actual feed/short supports unspecified",
        "heat_debit_J": None, "physical_coupled_ready": False,
    })
    return result
