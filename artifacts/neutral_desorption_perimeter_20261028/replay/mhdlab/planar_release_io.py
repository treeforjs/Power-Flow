"""Validated SI persistence for resolved planar heating and delivered sources."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, fields
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from .hdf_integrity import local_dataset
from .planar_release import (
    PLANAR_RELEASE_SCHEMA, PLANAR_SURFACE_SCHEMA, PlanarReleaseHistory, PlanarSurfaceHistory,
)
from .thermal_release import BindingReservoir


TIME_MODEL = "piecewise_frozen_midpoint_temperature_exact_first_order_local_depletion"
WIDTH_UNIFORM_MODEL = "width_redistribution_preserving_each_interval_reservoir_current_patch_number_and_energy"
FULL_FACE_UNIFORM_MODEL = "full_face_redistribution_preserving_each_interval_reservoir_number_and_energy"
SURFACE_UNITS = {"time_s": "s", "x_edges_m": "m", "z_edges_m": "m", "temperature_k": "K"}
RELEASE_UNITS = {
    "time_edges_s": "s", "x_edges_m": "m", "z_edges_m": "m", "temperature_k": "K",
    "coverage_m2": "m^-2", "released_m2": "m^-2",
    "translational_energy_j_m2": "J m^-2", "internal_energy_j_m2": "J m^-2",
    "binding_energy_j_m2": "J m^-2", "total_wall_energy_debit_j_m2": "J m^-2",
    "inventory_balance_residual_m2": "m^-2",
}


def plain_metadata(value: Any) -> Any:
    """Copy nested immutable metadata into JSON-compatible containers."""
    if isinstance(value, Mapping):
        return {key: plain_metadata(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain_metadata(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def _text(value: Any, name: str) -> str:
    if isinstance(value, bytes):
        try:
            value = value.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError(f"{name} must be UTF-8 text") from exc
    if not isinstance(value, str):
        raise ValueError(f"{name} must be a string")
    return value


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON metadata key: {key}")
        result[key] = value
    return result


def _nonfinite_constant(value):
    raise ValueError(f"nonfinite JSON metadata constant: {value}")


def _json_attribute(handle, name):
    if name not in handle.attrs:
        raise ValueError(f"missing required metadata: {name}")
    return json.loads(_text(handle.attrs[name], name), object_pairs_hook=_unique_object,
                      parse_constant=_nonfinite_constant)


def _attribute_equals(handle, name, expected):
    if name not in handle.attrs or _text(handle.attrs[name], name) != expected:
        raise ValueError(f"unsupported or missing {name}")


def _dataset(handle, path, units):
    item = local_dataset(handle, path)
    if item.dtype.kind not in "iuf":
        raise ValueError(f"{path} must be a real numeric dataset")
    if "units" not in item.attrs or _text(item.attrs["units"], f"{path}.units") != units:
        raise ValueError(f"{path} must declare units {units!r}")
    return item[...]


def _write_dataset(handle, name, value, units):
    dataset = handle.create_dataset(name, data=value)
    dataset.attrs["units"] = units


def _check_derived(handle, name, expected, units):
    value = _dataset(handle, name, units)
    if (value.shape != expected.shape or not np.all(np.isfinite(value))
            or not np.allclose(value, expected, rtol=32 * np.finfo(float).eps, atol=0.0)):
        raise ValueError(f"derived dataset disagrees with resolved source: {name}")


def write_planar_surface_history(path: str | Path, history: PlanarSurfaceHistory) -> None:
    """Write local (time,z,x) temperatures with all sampling/provenance metadata."""
    if not isinstance(history, PlanarSurfaceHistory):
        raise ValueError("history must be a PlanarSurfaceHistory")
    validated = PlanarSurfaceHistory(**{name: getattr(history, name) for name in SURFACE_UNITS},
                                     source_attributes=history.attributes)
    try:
        with h5py.File(path, "w") as handle:
            handle.attrs["schema"] = PLANAR_SURFACE_SCHEMA
            handle.attrs["source_attributes_json"] = json.dumps(validated.attributes, sort_keys=True, allow_nan=False)
            for name, units in SURFACE_UNITS.items():
                _write_dataset(handle, name, getattr(validated, name), units)
            _write_dataset(handle, "patch_areas_m2", validated.patch_areas_m2, "m^2")
    except (OSError, TypeError) as exc:
        raise ValueError(f"could not write planar surface {path}: {exc}") from exc


def load_planar_surface_history(path: str | Path) -> PlanarSurfaceHistory:
    """Load and revalidate geometry, local temperature association, and patch areas."""
    try:
        with h5py.File(path, "r") as handle:
            _attribute_equals(handle, "schema", PLANAR_SURFACE_SCHEMA)
            history = PlanarSurfaceHistory(
                **{name: _dataset(handle, name, units) for name, units in SURFACE_UNITS.items()},
                source_attributes=_json_attribute(handle, "source_attributes_json"),
            )
            _check_derived(handle, "patch_areas_m2", history.patch_areas_m2, "m^2")
            return history
    except (OSError, KeyError, TypeError, UnicodeError) as exc:
        raise ValueError(f"could not read planar surface {path}: {exc}") from exc


def write_planar_release_history(
    path: str | Path,
    history: PlanarReleaseHistory,
    *,
    include_full_face_uniform: bool = False,
) -> None:
    """Write local source densities and explicitly distinguished uniform controls."""
    if not isinstance(history, PlanarReleaseHistory):
        raise ValueError("history must be a PlanarReleaseHistory")
    if not isinstance(include_full_face_uniform, (bool, np.bool_)):
        raise ValueError("include_full_face_uniform must be a boolean")
    validated = PlanarReleaseHistory(
        **{name: getattr(history, name) for name in RELEASE_UNITS},
        reservoirs=history.reservoirs, surface_attributes=history.surface_attributes,
        parameter_provenance=history.parameter_provenance,
    )
    try:
        with h5py.File(path, "w") as handle:
            handle.attrs["schema"] = PLANAR_RELEASE_SCHEMA
            handle.attrs["time_model"] = TIME_MODEL
            handle.attrs["width_uniform_model"] = WIDTH_UNIFORM_MODEL
            handle.attrs["full_face_uniform_included"] = bool(include_full_face_uniform)
            handle.attrs["parameter_provenance"] = validated.parameter_provenance
            handle.attrs["surface_attributes_json"] = json.dumps(plain_metadata(validated.surface_attributes), sort_keys=True, allow_nan=False)
            handle.attrs["reservoirs_json"] = json.dumps([asdict(pool) for pool in validated.reservoirs], sort_keys=True, allow_nan=False)
            for name, units in RELEASE_UNITS.items():
                _write_dataset(handle, name, getattr(validated, name), units)
            _write_dataset(handle, "patch_areas_m2", validated.patch_areas_m2, "m^2")
            for name, value in validated.matched_uniform_width().items():
                _write_dataset(handle, f"matched_uniform_width/{name}", value, RELEASE_UNITS[name])
            if include_full_face_uniform:
                handle.attrs["full_face_uniform_model"] = FULL_FACE_UNIFORM_MODEL
                for name, value in validated.matched_uniform_full_face().items():
                    _write_dataset(handle, f"matched_uniform_full_face/{name}", value, RELEASE_UNITS[name])
    except (OSError, TypeError) as exc:
        raise ValueError(f"could not write planar release {path}: {exc}") from exc


def load_planar_release_history(path: str | Path) -> PlanarReleaseHistory:
    """Revalidate stored inventory/energy and regenerate every derived control."""
    try:
        with h5py.File(path, "r") as handle:
            _attribute_equals(handle, "schema", PLANAR_RELEASE_SCHEMA)
            _attribute_equals(handle, "time_model", TIME_MODEL)
            _attribute_equals(handle, "width_uniform_model", WIDTH_UNIFORM_MODEL)
            included = handle.attrs.get("full_face_uniform_included")
            if not isinstance(included, (bool, np.bool_)):
                raise ValueError("full_face_uniform_included must be an explicit boolean")
            records = _json_attribute(handle, "reservoirs_json")
            expected_keys = {item.name for item in fields(BindingReservoir)}
            if (not isinstance(records, list) or not records
                    or any(not isinstance(record, dict) or set(record) != expected_keys for record in records)):
                raise ValueError("reservoir records must include every declared binding parameter")
            if "parameter_provenance" not in handle.attrs:
                raise ValueError("missing parameter_provenance")
            history = PlanarReleaseHistory(
                **{name: _dataset(handle, name, units) for name, units in RELEASE_UNITS.items()},
                reservoirs=tuple(BindingReservoir(**record) for record in records),
                surface_attributes=_json_attribute(handle, "surface_attributes_json"),
                parameter_provenance=_text(handle.attrs["parameter_provenance"], "parameter_provenance"),
            )
            _check_derived(handle, "patch_areas_m2", history.patch_areas_m2, "m^2")
            for name, expected in history.matched_uniform_width().items():
                _check_derived(handle, f"matched_uniform_width/{name}", expected, RELEASE_UNITS[name])
            if included:
                _attribute_equals(handle, "full_face_uniform_model", FULL_FACE_UNIFORM_MODEL)
                for name, expected in history.matched_uniform_full_face().items():
                    _check_derived(handle, f"matched_uniform_full_face/{name}", expected, RELEASE_UNITS[name])
            elif "matched_uniform_full_face" in handle or "full_face_uniform_model" in handle.attrs:
                raise ValueError("undeclared full-face uniform control is present")
            return history
    except (OSError, KeyError, TypeError, UnicodeError) as exc:
        raise ValueError(f"could not read planar release {path}: {exc}") from exc
