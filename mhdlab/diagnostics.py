"""Streaming diagnostic writers and readers."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


FIELD_KEYS = (
    "temperature_k",
    "specific_enthalpy_j_kg",
    "pressure_pa",
    "bx_t",
    "by_t",
    "jz_a_m2",
    "conductivity_s_m",
    "joule_heating_w_m3",
    "surface_displacement_m",
    "electron_density_m3",
    "total_neutral_density_m3",
    "en_td",
)


class HDF5DiagnosticWriter:
    """Append sampled fields to an HDF5 file as each diagnostic slice is made."""

    def __init__(
        self,
        path: str | Path,
        shape: tuple[int, int],
        species: list[str],
        dtype: np.dtype | str = np.float32,
        compression: str | None = "gzip",
        compression_level: int | None = 4,
        metadata: dict[str, Any] | None = None,
    ):
        try:
            import h5py
        except ImportError as exc:
            raise RuntimeError("h5py is required for HDF5 diagnostics; install requirements.txt") from exc

        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.shape = tuple(int(v) for v in shape)
        self.dtype = np.dtype(dtype)
        self.species = list(species)
        self._sample_count = 0
        self._h5 = h5py.File(self.path, "w")
        self._compression_kwargs = _compression_kwargs(compression, compression_level)

        self._h5.attrs["schema"] = "mhdlab-diagnostics-v1"
        self._h5.attrs["created_utc"] = datetime.now(timezone.utc).isoformat()
        self._h5.attrs["completed_samples"] = 0
        self._h5.attrs["dtype"] = str(self.dtype)
        if metadata:
            meta_group = self._h5.create_group("metadata")
            for key, value in metadata.items():
                _write_attr(meta_group.attrs, key, value)

        self._time = self._h5.create_dataset("time_s", shape=(0,), maxshape=(None,), dtype="f8", chunks=(1,))
        self._fields = self._h5.create_group("fields")
        self._field_datasets = {
            key: self._create_time_dataset(self._fields, key)
            for key in FIELD_KEYS
        }
        self._species = self._h5.create_group("species")
        self._species_datasets = {}
        for name in self.species:
            group = self._species.create_group(name)
            self._species_datasets[name] = group.create_dataset(
                "density_m3",
                shape=(0, *self.shape),
                maxshape=(None, *self.shape),
                chunks=(1, *self.shape),
                dtype=self.dtype,
                **self._compression_kwargs,
            )
        self._h5.flush()

    @property
    def sample_count(self) -> int:
        return self._sample_count

    def append(self, sample: dict[str, Any]) -> int:
        index = self._sample_count
        self._time.resize((index + 1,))
        self._time[index] = float(sample["time_s"])
        for key, dataset in self._field_datasets.items():
            _append_array(dataset, index, sample[key], self.shape)
        for name, dataset in self._species_datasets.items():
            density = sample["species_density_m3"][name]
            _append_array(dataset, index, density, self.shape)
        self._sample_count += 1
        self._h5.attrs["completed_samples"] = self._sample_count
        self._h5.flush()
        return index

    def close(self) -> None:
        if self._h5:
            self._h5.flush()
            self._h5.close()

    def __enter__(self) -> "HDF5DiagnosticWriter":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def _create_time_dataset(self, group, name: str):
        return group.create_dataset(
            name,
            shape=(0, *self.shape),
            maxshape=(None, *self.shape),
            chunks=(1, *self.shape),
            dtype=self.dtype,
            **self._compression_kwargs,
        )


def load_hdf5_arrays(path: str | Path) -> dict[str, np.ndarray]:
    """Load HDF5 diagnostics into the legacy flat array dictionary shape."""
    try:
        import h5py
    except ImportError as exc:
        raise RuntimeError("h5py is required to read HDF5 diagnostics; install requirements.txt") from exc

    arrays: dict[str, np.ndarray] = {}
    with h5py.File(path, "r") as h5:
        arrays["time_s"] = h5["time_s"][...]
        for key, dataset in h5["fields"].items():
            arrays[key] = dataset[...]
        if "species" in h5:
            for species, group in h5["species"].items():
                arrays[f"species_{species}_density_m3"] = group["density_m3"][...]
    return arrays


def _compression_kwargs(compression: str | None, compression_level: int | None) -> dict[str, Any]:
    if compression is None or str(compression).lower() in {"", "none", "false"}:
        return {}
    kwargs: dict[str, Any] = {"compression": str(compression)}
    if str(compression).lower() == "gzip" and compression_level is not None:
        kwargs["compression_opts"] = int(compression_level)
    return kwargs


def _append_array(dataset, index: int, values: Any, shape: tuple[int, int]) -> None:
    arr = np.asarray(values)
    if arr.shape != shape:
        raise ValueError(f"diagnostic field {dataset.name} has shape {arr.shape}, expected {shape}")
    dataset.resize((index + 1, *shape))
    dataset[index, :, :] = arr


def _write_attr(attrs, key: str, value: Any) -> None:
    if isinstance(value, (str, int, float, bool, np.number)):
        attrs[key] = value
    elif value is None:
        attrs[key] = ""
    else:
        attrs[key] = str(value)
