"""Explicit SI sample-table handoff from an MFEM Joule surface evaluator.

This reads neither native MFEM grid functions nor arbitrary volume datasets.
The native producer must evaluate the material-side face on its mesh, resolve
width/current patches, and declare its units and temperature-zero convention.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
from numbers import Real

import numpy as np

from .flash_surface import COORDINATE_SYSTEM, SURFACE_MODEL
from .planar_release import MFEM_SURFACE_SCHEMA, PlanarSurfaceHistory
from .planar_release_io import write_planar_surface_history

CONFIG_SCHEMA = "mhdlab-mfem-surface-samples-v1"
COLUMNS = ("time_s", "z_index", "x_index", "temperature_value")


def _json(raw):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result
    def invalid(value):
        raise ValueError(f"nonfinite JSON constant: {value}")
    return json.loads(raw, object_pairs_hook=pairs, parse_constant=invalid)


def _text(value, name, *, empty=False):
    if not isinstance(value, str) or (not empty and not value.strip()):
        raise ValueError(f"{name} requires explicit text")
    return value


def _evidence(value, name):
    value = _text(value, name)
    if value.strip().upper().startswith("REPLACE_") or value.strip().lower() in {
            "unavailable", "pending", "unknown", "none", "n/a", "tbd"}:
        raise ValueError(f"{name} cannot retain placeholder evidence")
    return value


def _number(value, name, *, positive=False, nonnegative=False):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} requires a real number")
    try:
        value = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} requires finite numeric data") from exc
    if not math.isfinite(value) or positive and value <= 0 or nonnegative and value < 0:
        raise ValueError(f"{name} has invalid finite range")
    return value


def _partition(value, name, *, nonnegative=False):
    if not isinstance(value, list) or len(value) < 2:
        raise ValueError(f"{name} requires an explicit ordered partition")
    values = np.array([_number(x, name, nonnegative=nonnegative) for x in value])
    if np.any(np.diff(values) <= 0) or not np.isfinite(values[-1] - values[0]):
        raise ValueError(f"{name} must be strictly increasing with finite span")
    return values


def _hash(raw):
    return hashlib.sha256(raw).hexdigest()


def load_mfem_surface_samples(path: str | Path) -> tuple[PlanarSurfaceHistory, dict]:
    """Read a fully populated patch table with a hash-bound, explicit SI mapping.

    Each table value is already a local surface-patch average from the native
    producer. No source area, optical depth, exposure, current length or fitted
    geometric factor multiplies the temperature. No spatial averaging occurs.
    """
    path = Path(path).resolve()
    raw = path.read_bytes()
    config = _json(raw)
    if config.get("schema") != CONFIG_SCHEMA:
        raise ValueError("unsupported MFEM sample configuration schema")
    status = config.get("configuration_status")
    if status == "template":
        raise ValueError("template configuration cannot execute")
    if status not in {"configured", "manufactured"}:
        raise ValueError("configuration_status must be configured or manufactured")
    source = config["source"]
    kind = source["source_kind"]
    if (status, kind) not in {("configured", "mfem"), ("manufactured", "synthetic")}:
        raise ValueError("configured MFEM and manufactured synthetic provenance must remain separate")
    if source["heating_backend"] != "mfem_joule":
        raise ValueError("source must explicitly identify the MFEM Joule backend")
    for name in ("run_id", "source_reference", "unit_mapping_reference", "source_clock_reference"):
        _evidence(source[name], name)
    for name in ("mfem_version", "mfem_source_revision"):
        if kind == "mfem":
            _evidence(source[name], name)
        else:
            _text(source[name], name, empty=True)
    if kind == "mfem" and re.fullmatch(r"[0-9a-f]{40}", source["mfem_source_revision"]) is None:
        raise ValueError("mfem_source_revision must identify a full source commit")
    if source["native_field"] != "T":
        raise ValueError("the configured native temperature field must be T")
    if source["si_unit_mapping"] != {"length": "m", "time": "s", "temperature_value": "K"}:
        raise ValueError("sample tables require explicit SI units; no implicit unit conversion")
    limits = source["solver_scope"]
    if not isinstance(limits, dict) or limits.get("stationary_material") is not True:
        raise ValueError("this surface contract requires stationary material")
    for name in ("conductivity_temperature_feedback", "hall_induction", "material_motion"):
        if type(limits.get(name)) is not bool:
            raise ValueError("solver_scope must explicitly declare feedback, Hall and material-motion coverage")
    if limits["material_motion"]:
        raise ValueError("moving material requires a different surface contract")
    _evidence(limits.get("qualification_reference"), "solver qualification reference")
    temperature = config["temperature"]
    basis = temperature["value_basis"]
    if basis not in {"absolute_kelvin", "temperature_increment_kelvin"}:
        raise ValueError("temperature value_basis must declare absolute or increment kelvin")
    offset = _number(temperature["initial_temperature_k"], "initial_temperature_k", nonnegative=True)
    if basis == "absolute_kelvin" and offset != 0:
        raise ValueError("absolute kelvin data cannot receive an initial-temperature offset")
    if basis == "temperature_increment_kelvin" and offset <= 0:
        raise ValueError("increment data require an explicit positive initial temperature")
    _evidence(temperature["zero_convention_reference"], "temperature zero convention reference")
    location = temperature["location"]
    evaluation = temperature["evaluation"]
    if location == "surface_boundary":
        if evaluation != "one_sided_material_boundary_trace_patch_average":
            raise ValueError("DG boundary temperature requires a one-sided material trace, not a volume mean")
    elif location in {"material_sample_plane", "boundary_temperature_proxy"}:
        if evaluation != "resolved_material_sample_plane_patch_average":
            raise ValueError("material-plane temperature requires its sampled-plane association")
    else:
        raise ValueError("unsupported MFEM temperature location")
    face = config["face"]
    if face["coordinate_system"] != COORDINATE_SYSTEM or face["surface_model"] != SURFACE_MODEL:
        raise ValueError("this adapter requires one stationary planar face in the declared axes")
    _evidence(face["surface_region"], "surface_region")
    _number(face["normal_position_m"], "normal_position_m")
    if type(face["outward_normal_sign"]) is not int or face["outward_normal_sign"] not in {-1, 1}:
        raise ValueError("outward_normal_sign must be -1 or +1 along the gap normal")
    if type(face["material_boundary_attribute"]) is not int or face["material_boundary_attribute"] <= 0:
        raise ValueError("face requires an explicit positive material boundary attribute")
    for key in ("material_selection_reference", "full_width_coverage_reference"):
        _evidence(face[key], key)
    if face["current_direction_representation"] != "resolved_patches":
        raise ValueError("current-averaged traces cannot be inflated into resolved current patches")
    time = _partition(config["time_s"], "time_s", nonnegative=True)
    x = _partition(config["x_edges_m"], "x_edges_m")
    z = _partition(config["z_edges_m"], "z_edges_m")
    table = config["table"]
    if re.fullmatch(r"[0-9a-f]{64}", table["sha256"]) is None:
        raise ValueError("sample table requires its exact SHA256")
    table_path = Path(_text(table["path"], "table path"))
    table_path = (path.parent / table_path).resolve() if not table_path.is_absolute() else table_path.resolve()
    table_raw = table_path.read_bytes()
    if _hash(table_raw) != table["sha256"]:
        raise ValueError("sample table SHA256 mismatch")
    rows = csv.DictReader(io.StringIO(table_raw.decode("utf-8"), newline=""))
    if tuple(rows.fieldnames or ()) != COLUMNS:
        raise ValueError(f"sample table columns must be exactly {COLUMNS}")
    result = np.full((len(time), len(z) - 1, len(x) - 1), np.nan)
    time_indices = {float(value): index for index, value in enumerate(time)}
    count = 0
    for row in rows:
        if None in row or any(value is None for value in row.values()):
            raise ValueError("malformed sample table row")
        try:
            stamp = float(row["time_s"])
            if not math.isfinite(stamp) or stamp not in time_indices:
                raise ValueError("sample timestamp is not an explicitly declared snapshot time")
            if any(re.fullmatch(r"0|[1-9][0-9]*", row[key]) is None for key in ("x_index", "z_index")):
                raise ValueError("patch indices require nonnegative integers")
            ix, iz = int(row["x_index"]), int(row["z_index"])
            if ix >= len(x) - 1 or iz >= len(z) - 1:
                raise ValueError("sample index lies outside the declared face partition")
            it = time_indices[stamp]
            if not np.isnan(result[it, iz, ix]):
                raise ValueError("duplicate time/current/width sample")
            value = float(row["temperature_value"])
            converted = value + offset
            if not math.isfinite(converted) or converted <= 0:
                raise ValueError("converted absolute temperature must be finite and positive")
            result[it, iz, ix] = converted
            count += 1
        except (OverflowError, KeyError) as exc:
            raise ValueError("invalid numeric sample row") from exc
    if np.isnan(result).any():
        raise ValueError("sample table must cover every declared time/current/width patch exactly once")
    attributes = {
        "schema": MFEM_SURFACE_SCHEMA, "source_kind": kind, "heating_backend": "mfem_joule",
        **{key: source[key] for key in ("mfem_version", "mfem_source_revision", "run_id", "source_reference",
                                      "unit_mapping_reference", "source_clock_reference")},
        "coordinate_system": COORDINATE_SYSTEM, "surface_model": SURFACE_MODEL,
        "surface_region": face["surface_region"], "temperature_location": location,
        "material_sample_depth_m": _number(temperature["material_sample_depth_m"], "material_sample_depth_m", nonnegative=True),
        "material_sample_depth_units": "m", "boundary_temperature_reference": temperature["boundary_temperature_reference"],
        "current_averaging_interval_m_json": "null", "current_uniformity_reference": "",
        "temperature_value_basis": basis, "initial_temperature_offset_k": offset,
        "temperature_zero_convention_reference": temperature["zero_convention_reference"],
        "temperature_evaluation": evaluation, "mfem_face": face, "mfem_solver_scope": limits,
        "sample_table_sha256": table["sha256"], "sample_configuration_sha256": _hash(raw),
    }
    history = PlanarSurfaceHistory(time, x, z, result, attributes)
    manifest = {"schema": "mhdlab-mfem-surface-extraction-v1", "configuration_status": status,
        "configuration": {"path": str(path), "sha256": _hash(raw)},
        "sample_table": {"path": str(table_path), "sha256": table["sha256"], "rows": count},
        "source_attributes": history.attributes, "temperature_offset_applied_once_k": offset,
        "native_mfem_gridfunction_parsed": False, "new_heating_simulation_performed": False,
        "limitations": ["Native producer owns material selection, boundary evaluation and SI mapping.",
                        "No physical validity is established by matching metadata or hashes.",
                        "Local patch averages require x/z/time refinement before nonlinear release.",
                        "Temperatures are not an optical depth, neutral temperature or emitting extent."]}
    return history, manifest


def export_mfem_surface(path: str | Path, output_dir: str | Path) -> dict:
    """Write a canonical resolved planar history to a fresh output directory."""
    output = Path(output_dir).resolve()
    if output.exists():
        raise ValueError("use a fresh output directory")
    history, manifest = load_mfem_surface_samples(path)
    output.mkdir(parents=True)
    target = output / "surface.h5"
    write_planar_surface_history(target, history)
    manifest["surface_history"] = {"path": str(target), "sha256": _hash(target.read_bytes())}
    (output / "extraction_manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return manifest
