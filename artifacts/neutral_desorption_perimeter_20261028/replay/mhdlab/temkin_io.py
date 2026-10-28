"""Self-contained, law-specific HDF5 interchange for Temkin source histories."""

from __future__ import annotations

from dataclasses import asdict, fields
import json
from pathlib import Path

import h5py
import numpy as np

from .hdf_integrity import local_dataset
from .temkin_history import (
    SCHEMA, TIME_MODEL, UNIFORM_MODEL, REGIME_DEFINITIONS, RESIDUAL_INTERPRETATION,
    REPLAY_RELATIVE_TOLERANCE, TemkinReleaseHistory, _plain,
)
from .temkin_release import (
    LAW_IDENTITY, INVERSION_LOG_RESIDUAL_TOLERANCE, MAX_DIMENSIONLESS_COVERAGE_EFFECT, TemkinReservoir,
)


UNITS = {
    "time_edges_s": "s", "x_edges_m": "m", "temperature_k": "K", "coverage_m2": "m^-2",
    "released_m2": "m^-2", "translational_energy_j_m2": "J m^-2", "internal_energy_j_m2": "J m^-2",
    "binding_energy_j_m2": "J m^-2", "total_wall_energy_debit_j_m2": "J m^-2",
    "inventory_balance_residual_m2": "m^-2", "inversion_iterations": "1",
    "inversion_log_time_residual": "1", "inversion_regime": "1",
}
NUMERICAL_QUALIFICATIONS = {
    "maximum_dimensionless_coverage_effect": MAX_DIMENSIONLESS_COVERAGE_EFFECT,
    "inversion_log_residual_tolerance": INVERSION_LOG_RESIDUAL_TOLERANCE,
    "history_replay_relative_tolerance": REPLAY_RELATIVE_TOLERANCE,
}


def _object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON metadata key: {key}")
        result[key] = value
    return result


def _nonfinite(value):
    raise ValueError(f"nonfinite JSON metadata is unsupported: {value}")


def _json(value):
    return json.loads(value, object_pairs_hook=_object, parse_constant=_nonfinite)


def write_temkin_history(path: str | Path, history: TemkinReleaseHistory) -> None:
    """Write a qualified Temkin record; first-order histories are not accepted."""
    if not isinstance(history, TemkinReleaseHistory):
        raise ValueError("history must be a distinct TemkinReleaseHistory")
    validated = TemkinReleaseHistory(**{field.name: getattr(history, field.name) for field in fields(TemkinReleaseHistory)})
    try:
        with h5py.File(path, "w") as handle:
            for key, value in (("schema", SCHEMA), ("law_identity", LAW_IDENTITY), ("time_model", TIME_MODEL),
                               ("uniform_model", UNIFORM_MODEL), ("parameter_provenance", validated.parameter_provenance),
                               ("inversion_log_time_residual_interpretation", RESIDUAL_INTERPRETATION)):
                handle.attrs[key] = value
            handle.attrs["inversion_regime_definitions_json"] = json.dumps(REGIME_DEFINITIONS, sort_keys=True, allow_nan=False)
            handle.attrs["surface_attributes_json"] = json.dumps(_plain(validated.surface_attributes), sort_keys=True, allow_nan=False)
            handle.attrs["reservoirs_json"] = json.dumps([asdict(pool) for pool in validated.reservoirs], sort_keys=True, allow_nan=False)
            for name, value in NUMERICAL_QUALIFICATIONS.items():
                handle.attrs[name] = value
            for name, units in UNITS.items():
                handle.create_dataset(name, data=getattr(validated, name)).attrs["units"] = units
            uniform = handle.create_group("matched_uniform")
            for name, values in validated.matched_uniform().items():
                uniform.create_dataset(name, data=values).attrs["units"] = UNITS[name]
    except (OSError, TypeError) as error:
        raise ValueError(f"could not write Temkin history {path}: {error}") from error


def load_temkin_history(path: str | Path) -> TemkinReleaseHistory:
    """Validate metadata, arrays, local storage, kinetics and inversion records."""
    try:
        with h5py.File(path, "r") as handle:
            for key, expected in (("schema", SCHEMA), ("law_identity", LAW_IDENTITY), ("time_model", TIME_MODEL),
                                  ("uniform_model", UNIFORM_MODEL),
                                  ("inversion_log_time_residual_interpretation", RESIDUAL_INTERPRETATION)):
                if handle.attrs.get(key) != expected:
                    raise ValueError(f"missing or unsupported Temkin {key}")
            if _json(handle.attrs["inversion_regime_definitions_json"]) != REGIME_DEFINITIONS:
                raise ValueError("missing or unsupported Temkin inversion regime definitions")
            for name, value in NUMERICAL_QUALIFICATIONS.items():
                if handle.attrs.get(name) != value:
                    raise ValueError(f"missing or unsupported Temkin numerical qualification: {name}")
            arrays = {}
            for name, units in UNITS.items():
                dataset = local_dataset(handle, name)
                if dataset.attrs.get("units") != units:
                    raise ValueError(f"{name} must declare units {units!r}")
                arrays[name] = dataset[...]
            supplied = _json(handle.attrs["reservoirs_json"])
            if not isinstance(supplied, list) or any(not isinstance(pool, dict) for pool in supplied):
                raise ValueError("reservoirs_json must contain a list of Temkin parameter records")
            history = TemkinReleaseHistory(
                **arrays, reservoirs=tuple(TemkinReservoir(**pool) for pool in supplied),
                surface_attributes=_json(handle.attrs["surface_attributes_json"]),
                parameter_provenance=handle.attrs["parameter_provenance"],
                law_identity=handle.attrs["law_identity"], time_model=handle.attrs["time_model"],
            )
            for name, expected in history.matched_uniform().items():
                dataset = local_dataset(handle, f"matched_uniform/{name}")
                values = dataset[...]
                if (dataset.attrs.get("units") != UNITS[name] or values.shape != expected.shape
                        or values.dtype.kind not in "iuf" or not np.all(np.isfinite(values))
                        or not np.allclose(values, expected, rtol=16 * np.finfo(float).eps, atol=0)):
                    raise ValueError(f"matched Temkin source disagrees with resolved moments: {name}")
            return history
    except (OSError, KeyError, TypeError, OverflowError) as error:
        raise ValueError(f"could not load Temkin history {path}: {error}") from error
