"""Guarded Q1 foil geometry and reference-cell conservative transfer.

Nodes/cells follow the mechanics convention: structured (z,y,x), x fastest;
each hex is [000,100,110,010,001,101,111,011] in local x,y,z bits. SI throughout.
Quadrature controls describe benchmark geometry, never authorize a field solve.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np

from .thermal_expansion_mapping import CONSERVATION_RTOL, conservative_remap_1d


_BITS = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)
_SIGNS = 2 * _BITS - 1
_GAUSS = np.array(list(product((-1 / np.sqrt(3), 1 / np.sqrt(3)), repeat=3)))
_CHECK = np.vstack((_GAUSS, _SIGNS, np.zeros((1, 3))))
_FACES = {"x_minus": [0, 4, 7, 3], "x_plus": [1, 2, 6, 5],
          "y_minus": [0, 1, 5, 4], "y_plus": [3, 7, 6, 2],
          "z_minus": [0, 3, 2, 1], "z_plus": [4, 5, 6, 7]}
_QPOINTS = np.array(list(product((-1 / np.sqrt(3), 1 / np.sqrt(3)), repeat=2)))
_QSIGNS = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])


def _finite(value, name):
    original = np.asarray(value)
    if np.ma.isMaskedArray(value) or original.dtype.kind not in "fiu" or original.size == 0:
        raise ValueError(f"{name} must be a nonempty finite real array")
    result = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return result


def _frozen(value):
    array = np.asarray(value)
    return np.frombuffer(np.ascontiguousarray(array).tobytes(), dtype=array.dtype).reshape(array.shape)


def _shape_derivatives(points):
    result = np.empty((len(points), 8, 3))
    for direction in range(3):
        others = [axis for axis in range(3) if axis != direction]
        result[:, :, direction] = (_SIGNS[:, direction][None, :] / 8
            * np.prod(1 + points[:, None, others] * _SIGNS[None, :, others], axis=-1))
    return result


def _determinants(vertices, points):
    derivatives = _shape_derivatives(points)
    jacobian = np.einsum("eia,qib->eqab", vertices, derivatives)
    return np.linalg.det(jacobian)


@dataclass(frozen=True)
class DeformedHexGeometry:
    reference_nodes_m: np.ndarray
    nodes_m: np.ndarray
    cells: np.ndarray
    reference_cell_volumes_m3: np.ndarray
    cell_volumes_m3: np.ndarray
    cell_density_kg_m3: np.ndarray
    cell_mass_kg: np.ndarray
    gauss_jacobian_ratios: np.ndarray
    boundary_quads: np.ndarray
    boundary_cell_ids: np.ndarray
    material_face_ids: tuple[str, ...]
    boundary_triangles_m: np.ndarray
    triangle_normals: np.ndarray
    triangle_areas_m2: np.ndarray
    quad_gauss_points_m: np.ndarray
    quad_gauss_area_vectors_m2: np.ndarray
    controls: dict

    def receipt(self):
        return {
            "schema": "q1-solid-foil-geometry-v1", "units": "SI",
            "reference_volume_m3": float(np.sum(self.reference_cell_volumes_m3)),
            "current_volume_m3": float(np.sum(self.cell_volumes_m3)),
            "material_mass_kg": float(np.sum(self.cell_mass_kg)),
            "minimum_gauss_jacobian_ratio": float(np.min(self.gauss_jacobian_ratios)),
            "triangulated_boundary_area_m2": float(np.sum(self.triangle_areas_m2)),
            "material_face_ids": list(self.material_face_ids),
            "controls": dict(self.controls), "production_field_handoff_authorized": False,
            "surface_representation": "oriented two-triangle quads; Q1 Gauss surface vectors supplied separately",
            "limits": "Jacobian sampling is not a global injectivity proof; no support or thermal-work qualification inferred",
        }


def deform_hex_mesh(nodes_m, cells, displacement_m, *, reference_density_kg_m3,
                    cell_temperature_k, admissible_temperature_range_k) -> DeformedHexGeometry:
    """Map Q1 reference material geometry and enforce phase, mass and Jacobian guards.

    Volume uses 2x2x2 exact Q1 determinant quadrature. Positivity is checked at
    all Gauss points, corners and centers. This does not certify the Jacobian
    everywhere or exclude distant self-contact. Boundary scalar triangle areas
    are exact for their declared triangulation; warped Q1 faces are not planar.
    Separate Q1 surface quadrature preserves the volume-divergence identity.
    """
    reference = _finite(nodes_m, "nodes_m")
    displacement = _finite(displacement_m, "displacement_m")
    connectivity = np.asarray(cells)
    if reference.ndim != 2 or reference.shape[1] != 3 or displacement.shape != reference.shape:
        raise ValueError("nodes and displacements must have equal (node,3) shape")
    if connectivity.dtype.kind not in "iu" or connectivity.ndim != 2 or connectivity.shape[1] != 8 or not len(connectivity):
        raise ValueError("cells must contain eight integer node indices")
    if np.any(connectivity < 0) or np.any(connectivity >= len(reference)) or any(len(set(row)) != 8 for row in connectivity):
        raise ValueError("invalid or repeated cell node index")
    connectivity = connectivity.astype(np.int64)
    if len({tuple(sorted(row)) for row in connectivity}) != len(connectivity):
        raise ValueError("duplicate material cell")
    density = _finite(reference_density_kg_m3, "reference_density_kg_m3")
    if density.ndim != 0 or density <= 0:
        raise ValueError("positive scalar reference density required")
    temperature = _finite(cell_temperature_k, "cell_temperature_k")
    bounds = _finite(admissible_temperature_range_k, "admissible_temperature_range_k")
    if temperature.shape != (len(connectivity),) or bounds.shape != (2,) or not 0 < bounds[0] < bounds[1]:
        raise ValueError("one temperature per material cell and positive ordered range required")
    if np.any(temperature < bounds[0]) or np.any(temperature > bounds[1]):
        raise ValueError("cell temperature outside qualified solid range; no clipping or extrapolation")
    current = reference + displacement
    determinant0 = _determinants(reference[connectivity], _CHECK)
    determinant = _determinants(current[connectivity], _CHECK)
    if np.any(determinant0 <= 0) or np.any(~np.isfinite(determinant0)):
        raise ValueError("reference Q1 cell has nonpositive or invalid sampled Jacobian")
    ratios = determinant / determinant0
    if np.any(ratios <= 1e-8) or np.any(~np.isfinite(ratios)):
        raise ValueError("deformed Q1 cell has nonpositive or collapsed sampled Jacobian")
    volume0, volume = determinant0[:, :8].sum(axis=1), determinant[:, :8].sum(axis=1)
    mass = float(density) * volume0
    current_density = mass / volume
    mass_error = float(np.max(np.abs(current_density * volume - mass) / mass))
    if mass_error > CONSERVATION_RTOL:
        raise ValueError("cell mass conservation failed")

    faces = {}
    for cell_id, cell in enumerate(connectivity):
        for face_id, local in _FACES.items():
            face = cell[local]
            faces.setdefault(tuple(sorted(face)), []).append((cell_id, face_id, face))
    if any(len(owners) > 2 for owners in faces.values()):
        raise ValueError("nonmanifold hex boundary")
    boundary = [owners[0] for owners in faces.values() if len(owners) == 1]
    if not boundary:
        raise ValueError("empty solid boundary")
    cell_ids, face_ids, quads = zip(*boundary)
    quads = np.array(quads)
    vertices = current[quads]
    triangles = vertices[:, [[0, 1, 2], [0, 2, 3]], :]
    vectors = np.cross(triangles[:, :, 1] - triangles[:, :, 0], triangles[:, :, 2] - triangles[:, :, 0]) / 2
    areas = np.linalg.norm(vectors, axis=-1)
    if np.any(areas <= 0):
        raise ValueError("collapsed boundary triangle")
    normals = vectors / areas[:, :, None]
    qshape = np.prod(1 + _QPOINTS[:, None, :] * _QSIGNS[None, :, :], axis=-1) / 4
    derivatives = np.stack((
        _QSIGNS[None, :, 0] * (1 + _QPOINTS[:, None, 1] * _QSIGNS[None, :, 1]) / 4,
        _QSIGNS[None, :, 1] * (1 + _QPOINTS[:, None, 0] * _QSIGNS[None, :, 0]) / 4), axis=-1)
    qcoords = np.einsum("qi,fia->fqa", qshape, vertices)
    tangents = np.einsum("qib,fia->fqab", derivatives, vertices)
    qvectors = np.cross(tangents[:, :, :, 0], tangents[:, :, :, 1])
    total_area = float(areas.sum())
    closure = float(np.linalg.norm(vectors.sum(axis=(0, 1))) / total_area)
    qvolume = float(np.sum(qcoords * qvectors) / 3)
    divergence_error = abs(qvolume - float(volume.sum())) / float(volume.sum())
    if closure > CONSERVATION_RTOL or divergence_error > 5e-12:
        raise ValueError("oriented boundary vector-area/volume conservation failed")
    triangulated_volume = float(np.sum(np.mean(triangles, axis=2) * vectors) / 3)
    controls = {
        "cell_mass_relative_error": mass_error, "boundary_vector_area_relative_residual": closure,
        "q1_surface_volume_relative_residual": divergence_error,
        "triangulated_surface_volume_relative_difference": (triangulated_volume - float(volume.sum())) / float(volume.sum()),
        "jacobian_sample_locations": "8 Gauss + 8 corners + center",
        "boundary_self_contact_checked": False, "physical_coupled_ready": False,
    }
    arrays = [reference, current, connectivity, volume0, volume, current_density, mass, ratios[:, :8], quads,
              np.array(cell_ids), triangles, normals, areas, qcoords, qvectors]
    frozen = [_frozen(value) for value in arrays]
    return DeformedHexGeometry(*frozen[:10], tuple(face_ids), *frozen[10:], controls)


def _section_segments(geometry, face_id, z_m):
    selected = [index for index, name in enumerate(geometry.material_face_ids) if name == face_id]
    segments = []
    tolerance = 1e-13
    for triangle in geometry.boundary_triangles_m[selected].reshape(-1, 3, 3):
        relative = (triangle[:, 2] - z_m).copy()
        # Classify first, then intersect using the same classification. Mixing
        # on-plane tolerance with unsnapped signs creates a false third point
        # when a nominal on-plane edge has opposite roundoff-sized z residuals.
        relative[np.abs(relative) <= tolerance] = 0.
        if np.all(relative == 0):
            raise ValueError("inner surface is coplanar with observation plane; gap profile ambiguous")
        points = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            if relative[i] == 0:
                points.append(triangle[i, :2])
            if relative[i] * relative[j] < 0:
                fraction = -relative[i] / (relative[j] - relative[i])
                points.append(triangle[i, :2] + fraction * (triangle[j, :2] - triangle[i, :2]))
        unique = []
        for point in points:
            if not any(np.linalg.norm(point - other) < tolerance for other in unique):
                unique.append(point)
        if len(unique) > 2:
            raise ValueError("triangle-plane intersection has inconsistent topology")
        if len(unique) == 2:
            if abs(unique[0][0] - unique[1][0]) <= tolerance:
                raise ValueError("inner-face profile is vertical/multivalued at observation plane")
            segments.append(np.array(sorted(unique, key=lambda p: p[0])))
    if not segments:
        raise ValueError("observation plane does not cut a resolvable inner-face profile")
    return np.array(segments)


def _profile(segments, x):
    values = []
    for query in x:
        candidates = []
        for a, b in segments:
            if a[0] - 1e-13 <= query <= b[0] + 1e-13:
                # Only roundoff-close endpoints are snapped; never extend a
                # segment's slope beyond its physical support to fill a hole.
                sample = min(max(query, a[0]), b[0])
                candidates.append(a[1] + (sample - a[0]) * (b[1] - a[1]) / (b[0] - a[0]))
        if not candidates:
            raise ValueError(f"inner-face profile has a hole at x={query:.17g} m")
        if max(candidates) - min(candidates) > 1e-12:
            raise ValueError(f"inner-face profile has multiple y values at x={query:.17g} m")
        values.append(float(np.mean(candidates)))
    return np.array(values)


def midplane_gap(lower, upper, *, z_m=0.):
    """Piecewise-linear physical y gap at z=0, on common projected x support.

    Uses the explicitly triangulated inner faces. No extrapolation beyond either
    face, no closest-distance claim, and no continuation through contact.
    """
    z = float(_finite(z_m, "z_m"))
    a, b = _section_segments(lower, "y_plus", z), _section_segments(upper, "y_minus", z)
    left, right = max(a[:, :, 0].min(), b[:, :, 0].min()), min(a[:, :, 0].max(), b[:, :, 0].max())
    if right <= left:
        raise ValueError("no common inner-face x support")
    knots = np.unique(np.r_[left, right, a[:, :, 0].ravel(), b[:, :, 0].ravel()])
    knots = knots[(knots >= left) & (knots <= right)]
    # Midpoints certify segment coverage between every breakpoint.
    x = np.unique(np.r_[knots, (knots[:-1] + knots[1:]) / 2])
    lower_y, upper_y = _profile(a, x), _profile(b, x)
    gap = upper_y - lower_y
    if np.any(gap <= 0):
        raise ValueError("nonpositive midplane gap/contact")
    return {"x_m": _frozen(x), "lower_y_m": _frozen(lower_y), "upper_y_m": _frozen(upper_y),
            "gap_m": _frozen(gap), "minimum_gap_m": float(gap.min()), "maximum_gap_m": float(gap.max()),
            "observation_z_m": z, "representation": "triangulated inner-face y gap on common physical x support",
            "plane_classification_and_endpoint_tolerance_m": 1e-13,
            "production_field_handoff_authorized": False}


def parallel_plate_shape_sensitivity(voltage_V, gap_m, *, reference_gap_m=.005):
    """Conditional |E|=|V|/gap and d|E|/dg; excludes edge/fringing/inductive fields."""
    voltage, gap, initial_gap = (_finite(voltage_V, "voltage_V"), _finite(gap_m, "gap_m"),
                                 _finite(reference_gap_m, "reference_gap_m"))
    if voltage.ndim != 0 or initial_gap.ndim != 0 or initial_gap <= 0 or np.any(gap <= 0):
        raise ValueError("scalar voltage and positive gaps required")
    return {"electric_field_magnitude_V_m": _frozen(abs(voltage) / gap),
            "derivative_V_m2": _frozen(-abs(voltage) / gap ** 2),
            "field_ratio_to_reference": _frozen(initial_gap / gap),
            "assumption": "parallel plates at supplied voltage; conditional shape sensitivity only",
            "field_solve_performed": False, "production_field_handoff_authorized": False}


def remap_extensive_xy(source_x_edges_m, source_y_edges_m, target_x_edges_m, target_y_edges_m, values_J_or_kg):
    """Tensor-product conservative remap of material-cell totals, shape (ny,nx).

    Apply separately to energy in joules and mass in kilograms. Temperature is
    not extensive and cannot replace an enthalpy/EOS energy map. All grids are
    reference coordinates and must span exactly the same complete domain.
    """
    values = _finite(values_J_or_kg, "values_J_or_kg")
    if values.ndim != 2:
        raise ValueError("extensive xy values must have shape (ny,nx)")
    xmap = conservative_remap_1d(source_x_edges_m, target_x_edges_m, values.T, kind="extensive")
    ymap = conservative_remap_1d(source_y_edges_m, target_y_edges_m, xmap.values.T, kind="extensive")
    before, after = float(values.sum()), float(ymap.values.sum())
    scale = float(np.abs(values).sum())
    residual = abs(after - before) / scale if scale else 0.
    if residual > CONSERVATION_RTOL:
        raise ValueError("tensor-product conservation failure")
    return {"values": ymap.values, "source_total": before, "target_total": after,
            "relative_conservation_error": residual, "x_transfer": xmap.receipt(), "y_transfer": ymap.receipt(),
            "coordinate_system": "reference material coordinates", "depth_applied": False}


def extrude_extensive_xy(values_J_or_kg, z_edges_m):
    """Distribute a full-height xy total over axial cells; depth is used once."""
    values, z = _finite(values_J_or_kg, "values_J_or_kg"), _finite(z_edges_m, "z_edges_m")
    if values.ndim != 2 or z.ndim != 1 or len(z) < 2 or np.any(np.diff(z) <= 0):
        raise ValueError("xy totals and strictly increasing axial edges required")
    output = np.diff(z)[:, None, None] / (z[-1] - z[0]) * values[None, :, :]
    if not np.allclose(output.sum(axis=0), values, rtol=CONSERVATION_RTOL, atol=0.):
        raise ValueError("axial extensive distribution failed conservation")
    return _frozen(output)
