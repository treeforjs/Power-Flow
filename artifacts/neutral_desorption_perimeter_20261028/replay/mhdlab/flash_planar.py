"""Configured FLASH export retaining resolved width/current surface patches.

This offline reference rescans snapshots for each configured current interval.
It preserves local patch averages rather than reducing the whole current extent
to one temperature before applying a nonlinear release law.
"""

from __future__ import annotations

from dataclasses import replace
from numbers import Real
from typing import Sequence

import numpy as np

from .flash_extract import FlashLayout, FlashSnapshot, PlanarFace, extract_surface_history
from .planar_release import PlanarSurfaceHistory


EXTRACTION_SCHEMA = "mhdlab-flash-planar-extraction-v1"


def extract_planar_surface_history(
    snapshots: Sequence[FlashSnapshot],
    layout: FlashLayout,
    face: PlanarFace,
    z_edges_m,
    *,
    source_kind: str,
    run_id: str,
    source_reference: str,
    surface_region: str,
    flash_version: str = "",
) -> tuple[PlanarSurfaceHistory, dict]:
    """Return local (time,current,width) temperatures and their full provenance.

    The supplied current partition must cover the face's configured interval
    exactly. Each patch is a finite-volume spatial approximation that requires
    x/z refinement; its temperatures are not values at the material boundary.
    Two-dimensional raw runs retain their explicit invariant-current assumption.
    """
    if not isinstance(layout, FlashLayout) or not isinstance(face, PlanarFace):
        raise ValueError("layout and face must be validated configuration records")
    raw = np.asarray(z_edges_m, dtype=object)
    if (raw.ndim != 1 or raw.size < 2 or any(isinstance(v, (bool, np.bool_))
            or not isinstance(v, Real) for v in raw)):
        raise ValueError("z_edges_m must be a strictly increasing finite real partition")
    z = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(z)) or np.any(np.diff(z) <= 0):
        raise ValueError("z_edges_m must be strictly increasing and finite")
    if (z[0], z[-1]) != face.current_interval_m:
        raise ValueError("the current partition must exactly cover the configured face interval")
    inputs = tuple(snapshots)
    histories, extraction = [], []
    source_identity = None
    for lo, hi in zip(z[:-1], z[1:]):
        local_face = replace(face, current_interval_m=(float(lo), float(hi)))
        history, manifest = extract_surface_history(
            inputs, layout, local_face, source_kind=source_kind, run_id=run_id,
            source_reference=source_reference, surface_region=surface_region, flash_version=flash_version,
        )
        identity = [(item["path"], item["sha256"]) for item in manifest["snapshots"]]
        if source_identity is None:
            source_identity = identity
        elif identity != source_identity:
            raise ValueError("raw source files changed between current-patch extractions")
        histories.append(history)
        extraction.append(manifest)
    # This container now retains the current partition. Per-patch reductions and
    # their original intervals remain in the nested extraction records, instead
    # of claiming that a whole-current averaged trace has been resolved again.
    attributes = dict(histories[0].attributes)
    attributes["current_averaging_interval_m_json"] = "null"
    attributes["resolved_patch_extraction"] = {
        "schema": EXTRACTION_SCHEMA,
        "current_axis_in_raw_data": "current" in layout.physical_axes.values(),
        "projection": "local x/z patch area averages; no whole-current temperature reduction",
        "patch_extraction_manifests": extraction,
    }
    result = PlanarSurfaceHistory(
        histories[0].time_s, histories[0].x_edges_m, z,
        np.stack([history.temperature_k for history in histories], axis=1), attributes,
    )
    manifest = {
        "schema": EXTRACTION_SCHEMA, "source_attributes": result.attributes,
        "z_edges_m": z.tolist(), "current_patches": len(z) - 1,
        "snapshot_identity": [{"path": path, "sha256": sha} for path, sha in source_identity],
        "limitations": ["stationary plane sampled inside material, not an extrapolated boundary value",
                        "patch-temperature averages require width/current/time refinement",
                        "reference exporter rescans raw files for each current patch"],
    }
    return result, manifest


def as_planar_boundary_temperature_proxy(surface: PlanarSurfaceHistory, *, reference: str) -> PlanarSurfaceHistory:
    """Retain all patch/depth provenance while recording a boundary approximation."""
    if not isinstance(surface, PlanarSurfaceHistory) or surface.source_attributes["temperature_location"] != "material_sample_plane":
        raise ValueError("only a material sample plane can be designated a boundary-temperature proxy")
    if not isinstance(reference, str) or not reference.strip():
        raise ValueError("a boundary-temperature proxy requires an explicit evidence reference")
    attributes = surface.attributes
    attributes["temperature_location"] = "boundary_temperature_proxy"
    attributes["boundary_temperature_reference"] = reference
    return PlanarSurfaceHistory(surface.time_s, surface.x_edges_m, surface.z_edges_m,
                                surface.temperature_k, attributes)
