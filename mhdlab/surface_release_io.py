"""Versioned boundary-source output for a separately owned gap solver."""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

import h5py

from .surface_release_history import SurfaceReleaseHistory
from .hdf_integrity import local_dataset
from .thermal_release import BindingReservoir


SCHEMA = "mhdlab-surface-release-v1"
TIME_MODEL = "piecewise_frozen_midpoint_temperature_exact_first_order_depletion"
UNIFORM_MODEL = "spatial_redistribution_preserving_interval_reservoir_number_and_energy_moments"
UNITS = {
    "time_edges_s": "s", "x_edges_m": "m", "temperature_k": "K",
    "coverage_m2": "m^-2", "released_m2": "m^-2",
    "translational_energy_j_m2": "J m^-2", "internal_energy_j_m2": "J m^-2",
    "binding_energy_j_m2": "J m^-2", "total_wall_energy_debit_j_m2": "J m^-2",
    "inventory_balance_residual_m2": "m^-2",
}


def write_release_history(path: str | Path, history: SurfaceReleaseHistory) -> None:
    """Store interval-integrated source quantities, without sampling particles."""
    if not isinstance(history, SurfaceReleaseHistory):
        raise ValueError("history must be a SurfaceReleaseHistory")
    # Revalidate storage even if a caller has changed NumPy's read-only flags.
    validated = SurfaceReleaseHistory(
        **{name: getattr(history, name) for name in UNITS},
        reservoirs=history.reservoirs, surface_attributes=history.surface_attributes,
        parameter_provenance=history.parameter_provenance,
    )
    with h5py.File(path, "w") as handle:
        handle.attrs["schema"] = SCHEMA
        handle.attrs["time_model"] = TIME_MODEL
        handle.attrs["uniform_model"] = UNIFORM_MODEL
        handle.attrs["parameter_provenance"] = validated.parameter_provenance
        handle.attrs["surface_attributes_json"] = json.dumps(dict(validated.surface_attributes), sort_keys=True)
        handle.attrs["reservoirs_json"] = json.dumps([asdict(pool) for pool in validated.reservoirs], sort_keys=True)
        for name, units in UNITS.items():
            dataset = handle.create_dataset(name, data=getattr(validated, name))
            dataset.attrs["units"] = units
        uniform = handle.create_group("matched_uniform")
        for name, values in validated.matched_uniform().items():
            dataset = uniform.create_dataset(name, data=values)
            dataset.attrs["units"] = UNITS[name]


def load_release_history(path: str | Path) -> SurfaceReleaseHistory:
    """Read the source contract and verify its inventory and energy conventions.

    Uniform moments are regenerated from the resolved source. Their stored
    values are convenience products and are checked against that regeneration.
    """
    import numpy as np

    try:
        with h5py.File(path, "r") as handle:
            for key, expected in (("schema", SCHEMA), ("time_model", TIME_MODEL), ("uniform_model", UNIFORM_MODEL)):
                if handle.attrs.get(key) != expected:
                    raise ValueError(f"unsupported or missing {key}")
            arrays = {}
            for name, units in UNITS.items():
                dataset = local_dataset(handle, name)
                if dataset.attrs.get("units") != units:
                    raise ValueError(f"{name} must declare units {units!r}")
                arrays[name] = dataset[...]
            history = SurfaceReleaseHistory(
                **arrays,
                reservoirs=tuple(BindingReservoir(**pool) for pool in json.loads(handle.attrs["reservoirs_json"])),
                surface_attributes=json.loads(handle.attrs["surface_attributes_json"]),
                parameter_provenance=handle.attrs["parameter_provenance"],
            )
            for name, expected in history.matched_uniform().items():
                key = f"matched_uniform/{name}"
                dataset = local_dataset(handle, key)
                values = dataset[...]
                if (dataset.attrs.get("units") != UNITS[name]
                        or values.shape != expected.shape
                        or values.dtype.kind not in "iuf" or not np.all(np.isfinite(values))
                        or not np.allclose(values, expected, rtol=16 * np.finfo(float).eps, atol=0)):
                    raise ValueError(f"matched source dataset disagrees with resolved source: {name}")
            return history
    except (OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError(f"could not read release history {path}: {exc}") from exc
