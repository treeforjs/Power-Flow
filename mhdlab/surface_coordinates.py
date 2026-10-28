"""Explicit face-local to physical coordinates for stationary MFEM faces.

The local sampler always puts its surface tangent in column 0, signed outward
normal in column 1, and current direction in column 2. Side faces exchange the
first two physical axes. This orthogonal component permutation changes neither
lengths nor kinetic energy and supplies no artificial drift or speed scaling.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from .flash_surface import COORDINATE_SYSTEM

MFEM_SIDE_COORDINATE_SYSTEM = "u=foil_thickness;v=foil_width;z=current"


def validate_mfem_face_axes(face: Mapping) -> tuple[int, int, int]:
    """Return physical xyz's indices into local arrays, after strict checking.

    Legacy broad faces have local coordinates equal to physical xyz. Side
    faces require both physical axis declarations and use local (u,v,z) equal
    to physical (y,x,z). No general rotation or moving material is inferred.
    """
    if not isinstance(face, Mapping):
        raise ValueError("MFEM face axes require a metadata mapping")
    frame = face.get("coordinate_system")
    if frame == COORDINATE_SYSTEM:
        if (face.get("physical_in_plane_axis", "x") != "x"
                or face.get("physical_normal_axis", "y") != "y"):
            raise ValueError("MFEM broad-face physical axes conflict with its coordinate system")
        return (0, 1, 2)
    if frame == MFEM_SIDE_COORDINATE_SYSTEM:
        if (face.get("physical_in_plane_axis") != "y"
                or face.get("physical_normal_axis") != "x"):
            raise ValueError("MFEM side face requires explicit physical in-plane y and normal x axes")
        return (1, 0, 2)
    raise ValueError("unsupported MFEM face coordinate system")


def to_physical_coordinates(values, face: Mapping) -> np.ndarray:
    """Copy finite local position or velocity rows into physical xyz columns."""
    permutation = validate_mfem_face_axes(face)
    if np.ma.isMaskedArray(values):
        raise ValueError("face coordinate values cannot be masked")
    array = np.asarray(values)
    if array.ndim != 2 or array.shape[1] != 3 or array.dtype.kind not in "iuf":
        raise ValueError("face coordinate values require real numeric shape(n,3)")
    result = np.asarray(array, dtype=float)[:, permutation].copy()
    if not np.isfinite(result).all():
        raise ValueError("face coordinate values must remain finite")
    return result
