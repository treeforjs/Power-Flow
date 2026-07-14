"""Density-temperature electrical-conductivity tables."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class ConductivityTable:
    density_kg_m3: np.ndarray
    temperature_k: np.ndarray
    conductivity_s_m: np.ndarray
    source: str | None = None
    metadata: dict[str, Any] | None = None

    @classmethod
    def from_file(cls, path: str | Path) -> "ConductivityTable":
        path = Path(path)
        suffix = path.suffix.lower()
        if suffix == ".csv":
            return cls.from_csv(path)
        if suffix == ".npz":
            return cls.from_npz(path)
        if suffix in {".h5", ".hdf5"}:
            return cls.from_hdf5(path)
        raise ValueError(f"unsupported conductivity table format: {path}")

    @classmethod
    def from_csv(cls, path: str | Path) -> "ConductivityTable":
        path = Path(path)
        rows = []
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            required = {"density_kg_m3", "temperature_k", "conductivity_s_m"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")
            for row in reader:
                rows.append(
                    (
                        int(row.get("ensemble_index") or row.get("member") or 0),
                        float(row["density_kg_m3"]),
                        float(row["temperature_k"]),
                        float(row["conductivity_s_m"]),
                    )
                )
        if not rows:
            raise ValueError(f"conductivity table is empty: {path}")
        return cls._from_long_rows(rows, source=str(path))

    @classmethod
    def from_npz(cls, path: str | Path) -> "ConductivityTable":
        path = Path(path)
        with np.load(path) as data:
            return cls(
                density_kg_m3=np.asarray(data["density_kg_m3"], dtype=float),
                temperature_k=np.asarray(data["temperature_k"], dtype=float),
                conductivity_s_m=np.asarray(data["conductivity_s_m"], dtype=float),
                source=str(path),
            ).validated()

    @classmethod
    def from_hdf5(cls, path: str | Path) -> "ConductivityTable":
        path = Path(path)
        try:
            import h5py
        except ImportError as exc:
            raise RuntimeError("h5py is required to read HDF5 conductivity tables") from exc
        with h5py.File(path, "r") as h5:
            metadata = dict(h5.attrs.items())
            return cls(
                density_kg_m3=np.asarray(h5["density_kg_m3"][...], dtype=float),
                temperature_k=np.asarray(h5["temperature_k"][...], dtype=float),
                conductivity_s_m=np.asarray(h5["conductivity_s_m"][...], dtype=float),
                source=str(path),
                metadata=metadata,
            ).validated()

    @classmethod
    def _from_long_rows(cls, rows: list[tuple[int, float, float, float]], source: str | None = None) -> "ConductivityTable":
        members = np.asarray(sorted({r[0] for r in rows}), dtype=int)
        density = np.asarray(sorted({r[1] for r in rows}), dtype=float)
        temperature = np.asarray(sorted({r[2] for r in rows}), dtype=float)
        member_index = {value: i for i, value in enumerate(members)}
        density_index = {value: i for i, value in enumerate(density)}
        temperature_index = {value: i for i, value in enumerate(temperature)}
        values = np.full((members.size, density.size, temperature.size), np.nan, dtype=float)
        for member, rho, temp, sigma in rows:
            values[member_index[member], density_index[rho], temperature_index[temp]] = sigma
        if np.isnan(values).any():
            raise ValueError(f"conductivity table is not a complete rectangular grid: {source}")
        if members.size == 1:
            values = values[0]
        return cls(density, temperature, values, source=source, metadata={"members": members.tolist()}).validated()

    def validated(self) -> "ConductivityTable":
        density = np.asarray(self.density_kg_m3, dtype=float)
        temperature = np.asarray(self.temperature_k, dtype=float)
        conductivity = np.asarray(self.conductivity_s_m, dtype=float)
        if density.ndim != 1 or temperature.ndim != 1:
            raise ValueError("density and temperature axes must be one-dimensional")
        if density.size < 2 or temperature.size < 2:
            raise ValueError("conductivity table axes must each contain at least two points")
        if np.any(density <= 0.0) or np.any(temperature <= 0.0):
            raise ValueError("conductivity table axes must be positive")
        if np.any(np.diff(density) <= 0.0) or np.any(np.diff(temperature) <= 0.0):
            raise ValueError("conductivity table axes must be strictly increasing")
        expected = (density.size, temperature.size)
        if conductivity.ndim == 3:
            expected = (conductivity.shape[0], *expected)
        if conductivity.shape != expected:
            raise ValueError(f"conductivity table shape {conductivity.shape} does not match axes {expected}")
        if np.any(~np.isfinite(conductivity)) or np.any(conductivity <= 0.0):
            raise ValueError("conductivity table values must be finite and positive")
        return self

    def select_member(self, ensemble_index: int | None = None, statistic: str | None = None) -> np.ndarray:
        values = np.asarray(self.conductivity_s_m, dtype=float)
        if values.ndim == 2:
            return values
        if statistic:
            stat = statistic.lower()
            if stat in {"median", "p50"}:
                return np.median(values, axis=0)
            if stat == "mean":
                return np.mean(values, axis=0)
            if stat in {"min", "lower"}:
                return np.min(values, axis=0)
            if stat in {"max", "upper"}:
                return np.max(values, axis=0)
            raise ValueError(f"unsupported conductivity ensemble statistic: {statistic}")
        index = int(ensemble_index or 0)
        if index < 0 or index >= values.shape[0]:
            raise IndexError(f"conductivity ensemble index {index} outside 0..{values.shape[0] - 1}")
        return values[index]

    def interpolate(
        self,
        density_kg_m3: np.ndarray,
        temperature_k: np.ndarray,
        ensemble_index: int | None = None,
        statistic: str | None = None,
        log_interpolation: bool = True,
    ) -> np.ndarray:
        values = self.select_member(ensemble_index=ensemble_index, statistic=statistic)
        rho_q = np.asarray(density_kg_m3, dtype=float)
        temp_q = np.asarray(temperature_k, dtype=float)
        if log_interpolation:
            return np.exp(
                _interp2d(
                    np.log(self.density_kg_m3),
                    np.log(self.temperature_k),
                    np.log(values),
                    np.log(np.maximum(rho_q, self.density_kg_m3[0])),
                    np.log(np.maximum(temp_q, self.temperature_k[0])),
                )
            )
        return _interp2d(self.density_kg_m3, self.temperature_k, values, rho_q, temp_q)


def _interp2d(x_axis: np.ndarray, y_axis: np.ndarray, values: np.ndarray, xq: np.ndarray, yq: np.ndarray) -> np.ndarray:
    x = np.clip(xq, x_axis[0], x_axis[-1])
    y = np.clip(yq, y_axis[0], y_axis[-1])
    ix = np.clip(np.searchsorted(x_axis, x, side="right") - 1, 0, x_axis.size - 2)
    iy = np.clip(np.searchsorted(y_axis, y, side="right") - 1, 0, y_axis.size - 2)
    x0 = x_axis[ix]
    x1 = x_axis[ix + 1]
    y0 = y_axis[iy]
    y1 = y_axis[iy + 1]
    tx = np.divide(x - x0, x1 - x0, out=np.zeros_like(x, dtype=float), where=x1 != x0)
    ty = np.divide(y - y0, y1 - y0, out=np.zeros_like(y, dtype=float), where=y1 != y0)
    v00 = values[ix, iy]
    v10 = values[ix + 1, iy]
    v01 = values[ix, iy + 1]
    v11 = values[ix + 1, iy + 1]
    return (
        (1.0 - tx) * (1.0 - ty) * v00
        + tx * (1.0 - ty) * v10
        + (1.0 - tx) * ty * v01
        + tx * ty * v11
    )
