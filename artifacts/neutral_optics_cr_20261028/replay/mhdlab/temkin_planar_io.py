"""Self-contained SI interchange for local coverage-dependent planar sources."""

from __future__ import annotations

from dataclasses import asdict, fields
import json
from pathlib import Path

import h5py
import numpy as np

from .hdf_integrity import local_dataset
from .planar_release_io import (
    _attribute_equals, _check_derived, _json_attribute, _text,
    _write_dataset, plain_metadata,
)
from .temkin_history import REGIME_DEFINITIONS, RESIDUAL_INTERPRETATION, TIME_MODEL
from .temkin_io import NUMERICAL_QUALIFICATIONS
from .temkin_planar import (
    COVERAGE_MODEL, FULL_FACE_UNIFORM_MODEL, MIN_POSITIVE_UNIFORM_MOMENT_DENSITY, SCHEMA, WIDTH_UNIFORM_MODEL,
    TemkinPlanarReleaseHistory, _memory_limit,
)
from .temkin_release import LAW_IDENTITY, TemkinReservoir


RELEASE_UNITS = {
    "time_edges_s": "s", "x_edges_m": "m", "z_edges_m": "m", "temperature_k": "K",
    "coverage_m2": "m^-2", "released_m2": "m^-2",
    "translational_energy_j_m2": "J m^-2", "internal_energy_j_m2": "J m^-2",
    "binding_energy_j_m2": "J m^-2", "total_wall_energy_debit_j_m2": "J m^-2",
    "inventory_balance_residual_m2": "m^-2", "inversion_iterations": "1",
    "inversion_log_time_residual": "1", "inversion_regime": "1",
}
UNITS = RELEASE_UNITS
IDENTITIES = {
    "schema": SCHEMA, "law_identity": LAW_IDENTITY, "time_model": TIME_MODEL,
    "coverage_model": COVERAGE_MODEL, "width_uniform_model": WIDTH_UNIFORM_MODEL,
    "inversion_log_time_residual_interpretation": RESIDUAL_INTERPRETATION,
}
PLANAR_NUMERICAL_QUALIFICATIONS = {
    **NUMERICAL_QUALIFICATIONS,
    "minimum_positive_uniform_moment_density_float64": MIN_POSITIVE_UNIFORM_MOMENT_DENSITY,
}


def _stored_dataset(handle, name, units):
    """Validate local storage and units without materializing its array."""
    item = local_dataset(handle, name)
    if item.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a real numeric dataset")
    if "units" not in item.attrs or _text(item.attrs["units"], f"{name}.units") != units:
        raise ValueError(f"{name} must declare units {units!r}")
    return item


def _preflight(handle, nr, include_full_face):
    """Reject oversized or inconsistent HDF shapes before allocating arrays."""
    datasets = {name: _stored_dataset(handle, name, units) for name, units in RELEASE_UNITS.items()}
    for name in ("time_edges_s", "x_edges_m", "z_edges_m"):
        if len(datasets[name].shape) != 1 or datasets[name].shape[0] < 2:
            raise ValueError(f"{name} must be a vector containing at least two edges")
    nt = datasets["time_edges_s"].shape[0] - 1
    nx = datasets["x_edges_m"].shape[0] - 1
    nz = datasets["z_edges_m"].shape[0] - 1
    _memory_limit(nt, nr, nz, nx)
    for name, dataset in datasets.items():
        if name.endswith("edges_s") or name in ("x_edges_m", "z_edges_m"):
            continue
        expected = ((nt, nz, nx) if name == "temperature_k" else
                    ((nt + 1, nr, nz, nx) if name == "coverage_m2" else (nt, nr, nz, nx)))
        if dataset.shape != expected:
            raise ValueError(f"{name} shape disagrees with the declared time, reservoir and patch dimensions")
    if _stored_dataset(handle, "patch_areas_m2", "m^2").shape != (nz, nx):
        raise ValueError("patch_areas_m2 shape disagrees with the planar geometry")
    for group in ("matched_uniform_width", "matched_uniform_full_face") if include_full_face else ("matched_uniform_width",):
        for name in ("released_m2", "translational_energy_j_m2", "internal_energy_j_m2",
                     "binding_energy_j_m2", "total_wall_energy_debit_j_m2"):
            if _stored_dataset(handle, f"{group}/{name}", RELEASE_UNITS[name]).shape != (nt, nr, nz, nx):
                raise ValueError(f"{group}/{name} shape disagrees with the resolved source")
    return datasets


def write_temkin_planar_release_history(
    path: str | Path,
    history: TemkinPlanarReleaseHistory,
    *,
    include_full_face_uniform: bool = False,
) -> None:
    """Persist the identified local law, numerical qualifications and controls."""
    if not isinstance(history, TemkinPlanarReleaseHistory):
        raise ValueError("history must be a distinct TemkinPlanarReleaseHistory")
    if not isinstance(include_full_face_uniform, (bool, np.bool_)):
        raise ValueError("include_full_face_uniform must be a boolean")
    validated = TemkinPlanarReleaseHistory(**{field.name: getattr(history, field.name)
                                             for field in fields(TemkinPlanarReleaseHistory)})
    try:
        with h5py.File(path, "w") as handle:
            for name, value in IDENTITIES.items():
                handle.attrs[name] = value
            for name, value in PLANAR_NUMERICAL_QUALIFICATIONS.items():
                handle.attrs[name] = value
            handle.attrs["inversion_regime_definitions_json"] = json.dumps(REGIME_DEFINITIONS, sort_keys=True, allow_nan=False)
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
    except (OSError, TypeError) as error:
        raise ValueError(f"could not write Temkin planar release {path}: {error}") from error


def load_temkin_planar_release_history(path: str | Path) -> TemkinPlanarReleaseHistory:
    """Replay local kinetics/energy and recheck all derived spatial controls.

    Every dataset must live wholly inside this file. Neither a constant-hazard
    history nor missing qualification metadata is accepted by inference.
    """
    try:
        with h5py.File(path, "r") as handle:
            for name, value in IDENTITIES.items():
                _attribute_equals(handle, name, value)
            if _json_attribute(handle, "inversion_regime_definitions_json") != REGIME_DEFINITIONS:
                raise ValueError("missing or unsupported Temkin planar inversion regime definitions")
            for name, value in PLANAR_NUMERICAL_QUALIFICATIONS.items():
                if name not in handle.attrs or handle.attrs[name] != value:
                    raise ValueError(f"missing or unsupported Temkin planar numerical qualification: {name}")
            included = handle.attrs.get("full_face_uniform_included")
            if not isinstance(included, (bool, np.bool_)):
                raise ValueError("full_face_uniform_included must be an explicit boolean")
            records = _json_attribute(handle, "reservoirs_json")
            expected_keys = {field.name for field in fields(TemkinReservoir)}
            if (not isinstance(records, list) or not records
                    or any(not isinstance(record, dict) or set(record) != expected_keys for record in records)):
                raise ValueError("Temkin reservoir records must include every declared kinetic and energy parameter")
            if "parameter_provenance" not in handle.attrs:
                raise ValueError("missing parameter_provenance")
            if included:
                _attribute_equals(handle, "full_face_uniform_model", FULL_FACE_UNIFORM_MODEL)
            elif "matched_uniform_full_face" in handle or "full_face_uniform_model" in handle.attrs:
                raise ValueError("undeclared full-face uniform control is present")
            datasets = _preflight(handle, len(records), included)
            history = TemkinPlanarReleaseHistory(
                **{name: dataset[...] for name, dataset in datasets.items()},
                reservoirs=tuple(TemkinReservoir(**record) for record in records),
                surface_attributes=_json_attribute(handle, "surface_attributes_json"),
                parameter_provenance=_text(handle.attrs["parameter_provenance"], "parameter_provenance"),
                law_identity=_text(handle.attrs["law_identity"], "law_identity"),
                time_model=_text(handle.attrs["time_model"], "time_model"),
                coverage_model=_text(handle.attrs["coverage_model"], "coverage_model"),
                schema=_text(handle.attrs["schema"], "schema"),
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
    except (OSError, KeyError, TypeError, UnicodeError, OverflowError) as error:
        raise ValueError(f"could not read Temkin planar release {path}: {error}") from error
