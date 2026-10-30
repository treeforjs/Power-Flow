"""Auditable benchmark geometry and conservative material/field transfers.

Coordinates are SI and use x = foil width, y = gap direction, z = axial
direction.  Affine geometry does not determine supports or solve mechanics.
``transform_field`` implements geometric push-forwards only: it does not solve
Maxwell equations on a changed domain.  Nonuniform transient geometry remains
ineligible for field handoff until its mechanics and mesh transfer are qualified.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Literal

import numpy as np


CONSERVATION_RTOL = 5e-13
MIN_JACOBIAN = 1e-8
MAX_CONDITION_NUMBER = 1e6
REFERENCE_DIMENSIONS_M = (0.005, 0.0002, 0.025)
REFERENCE_GAP_M = 0.005
FACE_NAMES = ("x_minus", "x_plus", "y_minus", "y_plus", "z_minus", "z_plus")


def _array(value, name: str) -> np.ndarray:
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} must not contain masked values")
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf" or not raw.size:
        raise ValueError(f"{name} must be a nonempty real array")
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _readonly(array: np.ndarray) -> np.ndarray:
    shape = np.asarray(array).shape
    array = np.ascontiguousarray(array, dtype=float)
    return np.frombuffer(array.tobytes(), dtype=float).reshape(shape)


def _positive_scalar(value, name: str) -> float:
    array = _array(value, name)
    if array.ndim != 0 or array <= 0:
        raise ValueError(f"{name} must be a positive finite scalar")
    return float(array)


def _edges(value, name: str) -> np.ndarray:
    edges = _array(value, name)
    if edges.ndim != 1 or edges.size < 2 or np.any(np.diff(edges) <= 0):
        raise ValueError(f"{name} must be strictly increasing cell edges")
    return edges


def _hash(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        digest.update(str(array.shape).encode("ascii"))
        digest.update(np.asarray(array, dtype="<f8").tobytes())
    return digest.hexdigest()


@dataclass(frozen=True)
class ConservativeTransfer:
    """Computed transfer and immutable evidence used by the handoff guard."""

    values: np.ndarray
    source_edges: np.ndarray
    target_edges: np.ndarray
    source_values: np.ndarray
    kind: str
    relative_conservation_error: float
    source_total: np.ndarray
    target_total: np.ndarray
    input_sha256: str

    def receipt(self) -> dict:
        return {
            "representation": "piecewise_constant_overlap_integrals",
            "kind": self.kind,
            "source_total": self.source_total.tolist(),
            "target_total": self.target_total.tolist(),
            "relative_conservation_error": self.relative_conservation_error,
            "required_relative_tolerance": CONSERVATION_RTOL,
            "coverage": "identical_complete_domain_no_extrapolation",
            "input_sha256": self.input_sha256,
        }


def conservative_remap_1d(source_edges, target_edges, values, *,
                          kind: Literal["density", "extensive"] = "density") -> ConservativeTransfer:
    """Conservatively remap piecewise constant cells on the same complete domain.

    ``density`` represents amount per unit coordinate length; its width integral
    is conserved. ``extensive`` represents the amount already integrated over
    each source cell; its sum is conserved. Trailing value dimensions are allowed.
    Endpoints must be exactly equal: cropping and extrapolation are not inferred.
    This is a spatial remap, not a temporal temperature interpolator or a
    thermodynamic energy inversion. Temperature alone is not conserved energy.
    """
    source = _edges(source_edges, "source_edges")
    target = _edges(target_edges, "target_edges")
    data = _array(values, "values")
    if kind not in ("density", "extensive"):
        raise ValueError("kind must be density or extensive")
    if data.ndim < 1 or data.shape[0] != source.size - 1:
        raise ValueError("values first dimension must equal the source cell count")
    if source[0] != target[0] or source[-1] != target[-1]:
        raise ValueError("complete domains must match exactly; uncovered cells/extrapolation prohibited")
    source_width = np.diff(source)
    target_width = np.diff(target)
    overlap = np.maximum(0.0, np.minimum(target[1:, None], source[None, 1:])
                         - np.maximum(target[:-1, None], source[None, :-1]))
    if not np.allclose(overlap.sum(axis=0), source_width, rtol=CONSERVATION_RTOL, atol=0):
        raise ValueError("source coverage failed")
    if not np.allclose(overlap.sum(axis=1), target_width, rtol=CONSERVATION_RTOL, atol=0):
        raise ValueError("target coverage failed")
    trailing = (1,) * (data.ndim - 1)
    if kind == "density":
        density = data
        source_amount = data * source_width.reshape((-1,) + trailing)
    else:
        source_amount = data
        density = data / source_width.reshape((-1,) + trailing)
    target_amount = np.tensordot(overlap, density, axes=(1, 0))
    mapped = target_amount / target_width.reshape((-1,) + trailing) if kind == "density" else target_amount
    source_total, target_total = source_amount.sum(axis=0), target_amount.sum(axis=0)
    scale = np.abs(source_amount).sum(axis=0)
    difference = np.abs(target_total - source_total)
    relative = np.divide(difference, scale, out=np.zeros_like(difference), where=scale > 0)
    if np.any((scale == 0) & (difference != 0)):
        raise ValueError("zero-total transfer produced an amount")
    error = float(np.max(relative))
    if not np.all(np.isfinite(mapped)) or not np.isfinite(error) or error > CONSERVATION_RTOL:
        raise ValueError("conservative transfer failed numerical conservation controls")
    if np.all(data >= 0) and np.any(mapped < 0):
        raise ValueError("transfer violated positivity")
    return ConservativeTransfer(
        _readonly(mapped), _readonly(source), _readonly(target), _readonly(data), kind,
        error, _readonly(np.asarray(source_total)), _readonly(np.asarray(target_total)),
        _hash(source, target, data),
    )


def checked_deformation_gradient(deformation_gradient) -> tuple[np.ndarray, float]:
    """Require invertible, orientation-preserving, numerically usable geometry."""
    gradient = _array(deformation_gradient, "deformation_gradient")
    if gradient.shape != (3, 3):
        raise ValueError("deformation_gradient must have shape (3, 3)")
    jacobian = float(np.linalg.det(gradient))
    if not np.isfinite(jacobian) or jacobian <= MIN_JACOBIAN:
        raise ValueError("nonpositive or collapsed geometry Jacobian")
    if np.linalg.cond(gradient) > MAX_CONDITION_NUMBER:
        raise ValueError("deformation_gradient is too ill-conditioned")
    return gradient, jacobian


@dataclass(frozen=True)
class AffineFoilGeometry:
    """A parallelepiped geometry interface, with material mass conservation."""

    deformation_gradient: np.ndarray
    center_reference_m: np.ndarray
    translation_m: np.ndarray
    reference_dimensions_m: np.ndarray
    vertices_m: np.ndarray
    face_centers_m: np.ndarray
    face_normals: np.ndarray
    face_areas_m2: np.ndarray
    reference_density_kg_m3: float
    density_kg_m3: float
    jacobian: float
    volume_m3: float
    mass_kg: float

    def receipt(self) -> dict:
        return {
            "geometry_type": "affine_rectangular_foil_benchmark",
            "deformation_gradient": self.deformation_gradient.tolist(),
            "center_reference_m": self.center_reference_m.tolist(),
            "translation_m": self.translation_m.tolist(),
            "reference_dimensions_m": self.reference_dimensions_m.tolist(),
            "vertices_m": self.vertices_m.tolist(),
            "faces": {
                name: {"center_m": center.tolist(), "outward_normal": normal.tolist(), "area_m2": float(area)}
                for name, center, normal, area in zip(FACE_NAMES, self.face_centers_m, self.face_normals, self.face_areas_m2)
            },
            "jacobian": self.jacobian,
            "density_kg_m3": self.density_kg_m3,
            "volume_m3": self.volume_m3,
            "mass_kg": self.mass_kg,
            "density_law": "rho_reference / det(F)",
            "support_condition": "not determined by geometry",
        }


def affine_foil_geometry(deformation_gradient, *, center_reference_m=(0., 0., 0.),
                         translation_m=(0., 0., 0.),
                         reference_dimensions_m=REFERENCE_DIMENSIONS_M,
                         reference_density_kg_m3: float) -> AffineFoilGeometry:
    """Construct x = center_reference + translation + F @ (X - center_reference).

    A zero translation with F = stretch*I expands around each foil's own center.
    That is an explicitly chosen free benchmark origin, not an experimental
    restraint. Surface vectors use Nanson's formula J F^-T n dA. Density changes
    with total volume Jacobian, preserving the reference material mass exactly.
    """
    gradient, jacobian = checked_deformation_gradient(deformation_gradient)
    center = _array(center_reference_m, "center_reference_m")
    translation = _array(translation_m, "translation_m")
    dimensions = _array(reference_dimensions_m, "reference_dimensions_m")
    if center.shape != (3,) or translation.shape != (3,) or dimensions.shape != (3,) or np.any(dimensions <= 0):
        raise ValueError("center, translation and positive dimensions must each have shape (3,)")
    rho0 = _positive_scalar(reference_density_kg_m3, "reference_density_kg_m3")
    signs = np.array([[x, y, z] for x in (-1., 1.) for y in (-1., 1.) for z in (-1., 1.)])
    local_vertices = signs * dimensions / 2
    vertices = local_vertices @ gradient.T + center + translation
    normals_reference = np.repeat(np.eye(3), 2, axis=0) * np.tile([-1., 1.], 3)[:, None]
    centers = (normals_reference * dimensions / 2) @ gradient.T + center + translation
    areas_reference = np.repeat(np.prod(dimensions) / dimensions, 2)
    area_vectors = normals_reference @ (jacobian * np.linalg.inv(gradient)) * areas_reference[:, None]
    areas = np.linalg.norm(area_vectors, axis=1)
    normals = area_vectors / areas[:, None]
    volume0 = float(np.prod(dimensions))
    volume = volume0 * jacobian
    density = rho0 / jacobian
    if not np.isclose(density * volume, rho0 * volume0, rtol=CONSERVATION_RTOL, atol=0):
        raise ValueError("geometry density mapping failed mass conservation")
    return AffineFoilGeometry(
        _readonly(gradient), _readonly(center), _readonly(translation), _readonly(dimensions),
        _readonly(vertices), _readonly(centers), _readonly(normals), _readonly(areas),
        rho0, density, jacobian, volume, rho0 * volume0,
    )


def projected_gap_y_m(lower: AffineFoilGeometry, upper: AffineFoilGeometry) -> float:
    """Return conservative separation of the two complete y-projected solids.

    This equals the physical normal gap for aligned, unsheared foils. With shear
    or rotation it is a conservative y-clearance, not the closest-point distance.
    A nonpositive value rejects this simple interface, including possible overlap.
    """
    gap = float(np.min(upper.vertices_m[:, 1]) - np.max(lower.vertices_m[:, 1]))
    if gap <= 0:
        raise ValueError("foils do not have a positive separated y-projected gap")
    return gap


def transform_field(values, deformation_gradient, *, rule: Literal["vector", "covector", "hdiv"]) -> np.ndarray:
    """Geometric push-forward of explicit representations; never an EM solve.

    ``vector``: v = F v0 (material line/displacement representation).
    ``covector``: e = F^-T e0 (H(curl)/line-integral representation, e.g. E).
    ``hdiv``: b = F b0 / J (Piola flux representation, e.g. B or current flux).
    These preserve their respective geometric integrals. Physical fields after
    deformation still require a compatible field solve and boundary/circuit data.
    """
    gradient, jacobian = checked_deformation_gradient(deformation_gradient)
    data = _array(values, "values")
    if data.ndim < 1 or data.shape[-1] != 3:
        raise ValueError("field values must end with a three-component dimension")
    if rule == "vector":
        result = data @ gradient.T
    elif rule == "covector":
        result = data @ np.linalg.inv(gradient)
    elif rule == "hdiv":
        result = data @ gradient.T / jacobian
    else:
        raise ValueError("rule must explicitly be vector, covector or hdiv")
    return _readonly(result)


def qualify_field_handoff(*, geometries, temperatures_k, admissible_temperature_range_k,
                           transfers, mechanics_scope: str, support_description: str) -> dict:
    """Gate a *benchmark* geometry handoff; never authorize production integration.

    Transfers are re-evaluated from their source values and grids, avoiding trust
    in caller-supplied conservation claims. This interface deliberately rejects
    actual nonuniform transient mechanical deformation. Its FE equilibrium,
    dynamic response, mesh Jacobians, interface conservation and EM boundary
    mapping require separate qualification before a field solve may consume it.
    """
    if mechanics_scope not in ("uniform_affine_free_benchmark", "manufactured_affine_benchmark"):
        raise ValueError("mechanics scope is not qualified for this affine benchmark interface")
    if not isinstance(support_description, str) or not support_description.strip():
        raise ValueError("explicit support description is required")
    temperature = _array(temperatures_k, "temperatures_k")
    bounds = _array(admissible_temperature_range_k, "admissible_temperature_range_k")
    if bounds.shape != (2,) or bounds[0] <= 0 or bounds[0] >= bounds[1]:
        raise ValueError("admissible temperature range must be ordered positive kelvin bounds")
    if np.any(temperature < bounds[0]) or np.any(temperature > bounds[1]):
        raise ValueError("temperature exceeds qualified solid material range; no clipping or elastic extrapolation")
    if mechanics_scope == "uniform_affine_free_benchmark" and np.ptp(temperature) != 0:
        raise ValueError("uniform free benchmark requires uniform temperature")
    geometries, transfers = tuple(geometries), tuple(transfers)
    if not geometries or not transfers:
        raise ValueError("geometry and computed conservation transfers are required")
    geometry_receipts = []
    for geometry in geometries:
        if not isinstance(geometry, AffineFoilGeometry):
            raise ValueError("only computed affine benchmark geometry is supported")
        rebuilt = affine_foil_geometry(
            geometry.deformation_gradient, center_reference_m=geometry.center_reference_m,
            translation_m=geometry.translation_m, reference_dimensions_m=geometry.reference_dimensions_m,
            reference_density_kg_m3=geometry.reference_density_kg_m3,
        )
        for name in ("vertices_m", "face_centers_m", "face_normals", "face_areas_m2"):
            if not np.array_equal(getattr(geometry, name), getattr(rebuilt, name)):
                raise ValueError("geometry receipt differs from its recomputed map")
        for name in ("density_kg_m3", "jacobian", "volume_m3", "mass_kg"):
            if getattr(geometry, name) != getattr(rebuilt, name):
                raise ValueError("geometry mass/Jacobian receipt differs from recomputed map")
        geometry_receipts.append(rebuilt.receipt())
    if len(geometries) == 2:
        gap = projected_gap_y_m(geometries[0], geometries[1])
    else:
        gap = None
    transfer_receipts = []
    for transfer in transfers:
        if not isinstance(transfer, ConservativeTransfer):
            raise ValueError("computed ConservativeTransfer objects are required")
        rebuilt = conservative_remap_1d(transfer.source_edges, transfer.target_edges,
                                       transfer.source_values, kind=transfer.kind)
        if not np.array_equal(transfer.values, rebuilt.values) or transfer.input_sha256 != rebuilt.input_sha256:
            raise ValueError("transfer output differs from its recomputed conservative map")
        transfer_receipts.append(rebuilt.receipt())
    return {
        "status": "qualified_affine_benchmark_interface_only",
        "production_field_solve_authorized": False,
        "mechanics_scope": mechanics_scope,
        "support_description": support_description,
        "temperature_range_k": [float(np.min(temperature)), float(np.max(temperature))],
        "admissible_temperature_range_k": bounds.tolist(),
        "geometries": geometry_receipts,
        "y_projected_gap_m": gap,
        "transfers": transfer_receipts,
        "remaining_field_work": "compatible EM/thermal mesh and boundary/circuit solve; mapped fields are not a solved state",
    }
