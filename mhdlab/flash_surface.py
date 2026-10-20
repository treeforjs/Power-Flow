"""Validated SI surface histories exported from an external FLASH calculation.

This interchange format describes one stationary planar foil face per file. It
is deliberately not a raw FLASH AMR reader: extraction of a material surface,
FLASH unit conversion, and provenance must happen in a separate exporter. It
does not invoke the reduced solver, desorption, neutral transport, or CR model.
Width means are area means per unit length in the current (z) direction.
"""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from .hdf_integrity import local_dataset


SCHEMA = "mhdlab-flash-surface-v2"
COORDINATE_SYSTEM = "x=foil_width;y=gap_normal;z=current"
SURFACE_MODEL = "stationary_planar_single_face"
DATASET_UNITS = {"time_s": "s", "x_edges_m": "m", "temperature_k": "K"}


def _string(value: Any, name: str, *, allow_empty: bool = False) -> str:
    if isinstance(value, bytes):
        try:
            value = value.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError(f"{name} must be a UTF-8 string") from exc
    if not isinstance(value, str) or (not allow_empty and not value.strip()):
        raise ValueError(f"{name} must be a nonempty string")
    return value


def _array(value: Any, name: str, rank: int) -> np.ndarray:
    def contains_bool(item):
        if isinstance(item, (bool, np.bool_)):
            return True
        if isinstance(item, (list, tuple)):
            return any(contains_bool(part) for part in item)
        if isinstance(item, np.ndarray):
            return item.dtype.kind == "b" or (item.dtype.kind == "O" and any(contains_bool(part) for part in item.flat))
        return False
    try:
        raw = np.asarray(value)
        if contains_bool(value) or np.ma.isMaskedArray(value) or raw.dtype.kind not in "iuf":
            raise ValueError(f"{name} must contain unmasked real numbers without bool/complex/string coercion")
        result = np.array(raw, dtype=np.float64, copy=True)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must contain finite numeric values") from exc
    if result.ndim != rank or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must have rank {rank} and contain only finite values")
    return np.frombuffer(result.tobytes(), dtype=float).reshape(result.shape)


@dataclass(frozen=True)
class FlashSurfaceHistory:
    """SI temperatures, cell edges, and provenance for one named foil face.

    Arrays have shapes time_s=(nt,), x_edges_m=(nx+1,), and
    temperature_k=(nt,nx). Temperatures are surface-bin averages, treated as
    piecewise constant within each bin for width integration. The x axis
    spans the full foil width; no moving surface or additional face is implied.
    Synthetic inputs are supported only when explicitly labeled as such.
    """

    time_s: np.ndarray
    x_edges_m: np.ndarray
    temperature_k: np.ndarray
    source_kind: str
    run_id: str
    source_reference: str
    surface_region: str
    flash_version: str = ""
    temperature_location: str = "surface_boundary"
    material_sample_depth_m: float = 0.0
    boundary_temperature_reference: str = ""
    current_averaging_interval_m: tuple[float, float] | None = None
    current_uniformity_reference: str = ""

    def __post_init__(self) -> None:
        time = _array(self.time_s, "time_s", 1)
        edges = _array(self.x_edges_m, "x_edges_m", 1)
        temperature = _array(self.temperature_k, "temperature_k", 2)
        if time.size < 2 or np.any(time < 0.0) or np.any(time[1:] <= time[:-1]):
            raise ValueError("time_s must contain at least two nonnegative, strictly increasing times")
        if edges.size < 2 or np.any(edges[1:] <= edges[:-1]):
            raise ValueError("x_edges_m must contain at least two strictly increasing edges")
        with np.errstate(over="ignore"):
            width = edges[-1] - edges[0]
        if not np.isfinite(width):
            raise ValueError("x_edges_m must span a finite foil width")
        if temperature.shape != (time.size, edges.size - 1):
            raise ValueError("temperature_k shape must be (len(time_s), len(x_edges_m)-1)")
        if np.any(temperature <= 0.0):
            raise ValueError("temperature_k must be strictly positive in kelvin")
        for name, value in (("time_s", time), ("x_edges_m", edges), ("temperature_k", temperature)):
            object.__setattr__(self, name, value)
        for name in ("source_kind", "run_id", "source_reference", "surface_region", "flash_version",
                     "temperature_location", "boundary_temperature_reference", "current_uniformity_reference"):
            object.__setattr__(self, name, _string(getattr(self, name), name,
                                                  allow_empty=name in {"flash_version", "boundary_temperature_reference", "current_uniformity_reference"}))
        if self.source_kind not in {"flash", "synthetic"}:
            raise ValueError("source_kind must be 'flash' or 'synthetic'")
        if self.source_kind == "flash" and not self.flash_version.strip():
            raise ValueError("flash_version must be nonempty for a FLASH source")
        if self.temperature_location not in {"surface_boundary", "material_sample_plane", "boundary_temperature_proxy"}:
            raise ValueError("temperature_location must identify a boundary, material plane, or justified boundary proxy")
        depth = self.material_sample_depth_m
        if isinstance(depth, (bool, np.bool_)) or not isinstance(depth, (int, float, np.integer, np.floating)):
            raise ValueError("material_sample_depth_m must be a finite nonnegative real number")
        depth = float(depth)
        if not np.isfinite(depth) or depth < 0:
            raise ValueError("material_sample_depth_m must be finite and nonnegative")
        object.__setattr__(self, "material_sample_depth_m", depth)
        if self.temperature_location == "surface_boundary" and depth != 0:
            raise ValueError("surface_boundary requires zero material_sample_depth_m")
        if self.temperature_location != "surface_boundary" and depth == 0:
            raise ValueError("material plane and proxy temperatures require a positive sampling depth")
        if self.temperature_location == "material_sample_plane" and self.boundary_temperature_reference:
            raise ValueError("a material sample plane cannot carry a boundary conversion reference")
        if self.temperature_location == "boundary_temperature_proxy" and not self.boundary_temperature_reference.strip():
            raise ValueError("a boundary-temperature proxy requires an explicit justification reference")
        if self.current_averaging_interval_m is not None:
            try:
                interval = tuple(self.current_averaging_interval_m)
            except TypeError as exc:
                raise ValueError("current_averaging_interval_m must be a finite increasing pair or None") from exc
            if (len(interval) != 2 or any(isinstance(v, (bool, np.bool_))
                    or not isinstance(v, (int, float, np.integer, np.floating)) for v in interval)):
                raise ValueError("current_averaging_interval_m must be a finite increasing real pair")
            interval = tuple(float(v) for v in interval)
            if not np.all(np.isfinite(interval)) or interval[1] <= interval[0] or not np.isfinite(interval[1] - interval[0]):
                raise ValueError("current_averaging_interval_m must span a finite positive length")
            object.__setattr__(self, "current_averaging_interval_m", interval)

    @property
    def attributes(self) -> dict[str, str | float]:
        return {
            "schema": SCHEMA,
            "source_kind": self.source_kind,
            "flash_version": self.flash_version,
            "run_id": self.run_id,
            "source_reference": self.source_reference,
            "coordinate_system": COORDINATE_SYSTEM,
            "surface_model": SURFACE_MODEL,
            "surface_region": self.surface_region,
            "temperature_location": self.temperature_location,
            "material_sample_depth_m": self.material_sample_depth_m,
            "material_sample_depth_units": "m",
            "boundary_temperature_reference": self.boundary_temperature_reference,
            "current_averaging_interval_m_json": json.dumps(self.current_averaging_interval_m),
            "current_uniformity_reference": self.current_uniformity_reference,
        }

    def temperature_at(self, time_s: float | np.ndarray) -> np.ndarray:
        """Interpolate each bin linearly in time; never extrapolate."""
        try:
            query = np.asarray(time_s, dtype=np.float64)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("query time_s must be finite numeric values") from exc
        if query.ndim > 1 or not np.all(np.isfinite(query)) or query.size == 0:
            raise ValueError("query time_s must be a finite scalar or nonempty one-dimensional array")
        if np.any(query < self.time_s[0]) or np.any(query > self.time_s[-1]):
            raise ValueError("query time_s lies outside the surface history; extrapolation is not allowed")
        return np.stack(
            [np.interp(query, self.time_s, column) for column in self.temperature_k.T],
            axis=-1,
        )

    def thermal_contrast(
        self,
        *,
        edge_fraction: float = 0.1,
        center_fraction: float = 0.2,
        time_s: float | np.ndarray | None = None,
    ) -> dict[str, np.ndarray]:
        """Return exact overlap-weighted temperatures and thermal contrasts.

        Each edge occupies edge_fraction of the full width. The centered
        interval occupies center_fraction. Intervals may touch but not overlap.
        The combined edge mean weights the two equal-width edge intervals.
        These metrics alone say nothing about desorption or spectral emission.
        """
        try:
            edge = float(edge_fraction)
            center = float(center_fraction)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("region fractions must be finite scalars") from exc
        if not np.isfinite(edge) or not np.isfinite(center) or edge <= 0.0 or center <= 0.0:
            raise ValueError("region fractions must be finite and strictly positive")
        if 2.0 * edge + center > 1.0:
            raise ValueError("edge and center intervals must not overlap or exceed the foil width")
        start, stop = self.x_edges_m[[0, -1]]
        width = stop - start
        midpoint = start + 0.5 * width
        temperature = self.temperature_k if time_s is None else self.temperature_at(time_s)

        def mean(lo: float, hi: float) -> np.ndarray:
            overlap = np.maximum(
                0.0,
                np.minimum(self.x_edges_m[1:], hi) - np.maximum(self.x_edges_m[:-1], lo),
            )
            covered_width = np.sum(overlap)
            if covered_width <= 0.0:
                raise ValueError("region width is too small to resolve numerically on this x axis")
            return np.sum(temperature * (overlap / covered_width), axis=-1)

        left = mean(start, start + edge * width)
        right = mean(stop - edge * width, stop)
        middle = mean(midpoint - 0.5 * center * width, midpoint + 0.5 * center * width)
        combined = 0.5 * (left + right)
        return {
            "time_s": self.time_s.copy() if time_s is None else np.asarray(time_s, dtype=float),
            "left_edge_mean_temperature_k": left,
            "right_edge_mean_temperature_k": right,
            "combined_edge_mean_temperature_k": combined,
            "center_mean_temperature_k": middle,
            "edge_minus_center_k": combined - middle,
            "edge_to_center_ratio": combined / middle,
        }


def load_surface_history(path: str | Path) -> FlashSurfaceHistory:
    """Read the canonical SI interchange file and reject incomplete metadata."""
    try:
        with h5py.File(path, "r") as handle:
            attributes = {}
            for name in (
                "schema", "source_kind", "flash_version", "run_id", "source_reference",
                "coordinate_system", "surface_model", "surface_region",
                "temperature_location", "material_sample_depth_units", "boundary_temperature_reference",
                "current_averaging_interval_m_json", "current_uniformity_reference",
            ):
                if name not in handle.attrs:
                    raise ValueError(f"missing required file attribute: {name}")
                attributes[name] = _string(handle.attrs[name], name,
                                            allow_empty=name in {"flash_version", "boundary_temperature_reference", "current_uniformity_reference"})
            if "material_sample_depth_m" not in handle.attrs:
                raise ValueError("missing required file attribute: material_sample_depth_m")
            attributes["material_sample_depth_m"] = handle.attrs["material_sample_depth_m"]
            for name, expected in (
                ("schema", SCHEMA), ("coordinate_system", COORDINATE_SYSTEM), ("surface_model", SURFACE_MODEL),
                ("material_sample_depth_units", "m"),
            ):
                if attributes[name] != expected:
                    raise ValueError(f"{name} must equal {expected!r}")
            arrays = {}
            for name, units in DATASET_UNITS.items():
                dataset = local_dataset(handle, name)
                if "units" not in dataset.attrs or _string(dataset.attrs["units"], f"{name}.units") != units:
                    raise ValueError(f"{name}.units must equal {units!r}; no unit inference is performed")
                arrays[name] = dataset[...]
            return FlashSurfaceHistory(
                **arrays,
                **{key: attributes[key] for key in (
                    "source_kind", "flash_version", "run_id", "source_reference", "surface_region",
                    "temperature_location", "material_sample_depth_m", "boundary_temperature_reference",
                    "current_uniformity_reference",
                )},
                current_averaging_interval_m=json.loads(attributes["current_averaging_interval_m_json"]),
            )
    except (OSError, TypeError, KeyError) as exc:
        raise ValueError(f"could not read surface history {path}: {exc}") from exc


def write_surface_history(path: str | Path, history: FlashSurfaceHistory) -> None:
    """Write one validated face in canonical SI units; never infer CGS units."""
    if not isinstance(history, FlashSurfaceHistory):
        raise ValueError("history must be a FlashSurfaceHistory")
    # Revalidate even if a caller has overridden NumPy's read-only array flag.
    validated = FlashSurfaceHistory(
        history.time_s, history.x_edges_m, history.temperature_k,
        history.source_kind, history.run_id, history.source_reference,
        history.surface_region, history.flash_version,
        history.temperature_location, history.material_sample_depth_m, history.boundary_temperature_reference,
        history.current_averaging_interval_m, history.current_uniformity_reference,
    )
    try:
        with h5py.File(path, "w") as handle:
            for key, value in validated.attributes.items():
                handle.attrs[key] = value
            for name, units in DATASET_UNITS.items():
                dataset = handle.create_dataset(name, data=getattr(validated, name))
                dataset.attrs["units"] = units
    except (OSError, TypeError) as exc:
        raise ValueError(f"could not write surface history {path}: {exc}") from exc


def thermal_width_report(
    history: FlashSurfaceHistory,
    *,
    edge_fraction: float = 0.1,
    center_fraction: float = 0.2,
    time_s: float | None = None,
) -> dict[str, Any]:
    """Build a JSON-compatible, explicitly labeled thermal-only report."""
    metrics = history.thermal_contrast(
        edge_fraction=edge_fraction, center_fraction=center_fraction, time_s=time_s,
    )
    return {
        **history.attributes,
        "interpretation": "temperature contrast only, at the declared location; no desorption or spectrum prediction",
        "foil_width_m": float(history.x_edges_m[-1] - history.x_edges_m[0]),
        "edge_fraction_per_side": float(edge_fraction),
        "center_fraction": float(center_fraction),
        "metrics": {name: np.asarray(value).tolist() for name, value in metrics.items()},
    }


def as_boundary_temperature_proxy(history: FlashSurfaceHistory, *, reference: str) -> FlashSurfaceHistory:
    """Explicitly designate a material-plane temperature as a boundary proxy.

    Values and sampling depth are preserved. ``reference`` must identify the
    gradient, convergence, analytic fixture or other evidence supporting this
    approximation. The operation records a modeling decision; it does not
    extrapolate the temperature to the material boundary or validate the claim.
    """
    if not isinstance(history, FlashSurfaceHistory) or history.temperature_location != "material_sample_plane":
        raise ValueError("only a material_sample_plane history can be designated a boundary proxy")
    justification = _string(reference, "boundary temperature proxy reference")
    return FlashSurfaceHistory(
        history.time_s, history.x_edges_m, history.temperature_k, history.source_kind,
        history.run_id, history.source_reference, history.surface_region, history.flash_version,
        "boundary_temperature_proxy", history.material_sample_depth_m, justification,
        history.current_averaging_interval_m, history.current_uniformity_reference,
    )


def as_current_uniform_temperature_proxy(history: FlashSurfaceHistory, *, reference: str) -> FlashSurfaceHistory:
    """Record evidence for applying a current-averaged T to a nonlinear source.

    Prefer separately exported thin current-direction slices when temperature
    varies there. This designation preserves the averaging interval and values;
    it does not recover unresolved variations or establish uniformity itself.
    """
    if not isinstance(history, FlashSurfaceHistory) or history.current_averaging_interval_m is None:
        raise ValueError("only a current-averaged history requires this proxy designation")
    justification = _string(reference, "current uniformity reference")
    return FlashSurfaceHistory(
        history.time_s, history.x_edges_m, history.temperature_k, history.source_kind,
        history.run_id, history.source_reference, history.surface_region, history.flash_version,
        history.temperature_location, history.material_sample_depth_m, history.boundary_temperature_reference,
        history.current_averaging_interval_m, justification,
    )


def validate_surface_attributes(attributes: Mapping) -> dict[str, str | float]:
    """Validate every semantic key using the temperature-history record rules.

    Unit-sized arrays exercise only the record's metadata validation here; no
    physical temperature history is constructed or returned by this interface.
    This keeps file consumers on the same depth, timing and provenance rules.
    """
    if not isinstance(attributes, Mapping):
        raise ValueError("surface attributes must be a mapping")
    fields = ("source_kind", "run_id", "source_reference", "surface_region", "flash_version",
              "temperature_location", "material_sample_depth_m", "boundary_temperature_reference",
              "current_uniformity_reference")
    required = {*fields, "schema", "coordinate_system", "surface_model", "material_sample_depth_units",
                "current_averaging_interval_m_json"}
    missing = required - attributes.keys()
    if missing:
        raise ValueError(f"missing required surface attributes: {', '.join(sorted(missing))}")
    for key, expected in (("schema", SCHEMA), ("coordinate_system", COORDINATE_SYSTEM),
                          ("surface_model", SURFACE_MODEL), ("material_sample_depth_units", "m")):
        if attributes[key] != expected:
            raise ValueError(f"surface attribute {key} must equal {expected!r}")
    raw_interval = attributes["current_averaging_interval_m_json"]
    if not isinstance(raw_interval, str):
        raise ValueError("current_averaging_interval_m_json must be an explicit JSON string")
    try:
        interval = json.loads(raw_interval)
    except (TypeError, ValueError) as exc:
        raise ValueError("invalid current-averaging interval metadata") from exc
    record = FlashSurfaceHistory(
        [0.0, 1.0], [0.0, 1.0], [[1.0], [1.0]],
        **{key: attributes[key] for key in fields}, current_averaging_interval_m=interval,
    )
    return record.attributes
