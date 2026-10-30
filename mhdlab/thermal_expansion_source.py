"""Hash-bound, read-only thermal input for solid-expansion benchmarks.

The saved native volume CSV contains upper-foil cell averages in absolute K.
The eight face tables contain temperature increments.  Lower volume was not
exported: its optional reconstruction is explicit and checked against all four
independently exported lower face histories.  It is never called native data.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import io
import json
from pathlib import Path

import numpy as np


MATCHED_256_SOURCE_SHA256 = "df79133957684abbe5678b10b0c4b7625cbb0c8b7e17d31712410decd5efca05"
_BULK_COLUMNS = ("time_s,x_left_m,x_right_m,depth_left_m,depth_right_m,"
                 "temperature_k,current_density_a_m2,joule_power_w_m3")
_FACE_COLUMNS = "time_s,z_index,x_index,temperature_value"


def _readonly(value):
    value = np.asarray(value, dtype=float)
    value.setflags(write=False)
    return value


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _partition(values, name):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all() or np.any(np.diff(values) <= 0):
        raise ValueError(f"{name} must be finite and strictly increasing")
    return values


@dataclass(frozen=True)
class MatchedThermalHistory:
    time_s: np.ndarray
    x_edges_m: np.ndarray
    depth_edges_m: np.ndarray
    upper_y_edges_m: np.ndarray
    lower_y_edges_m: np.ndarray
    upper_temperature_k: np.ndarray  # [time, ascending y, ascending x]
    lower_temperature_k: np.ndarray | None
    native_fields: tuple[dict, ...]
    face_time_s: np.ndarray
    face_temperature_k: dict[str, np.ndarray]
    manifest: dict

    def admissibility(self, minimum_k: float, maximum_k: float) -> dict:
        """Report invalid cells; do not clip, discard or interpolate through them."""
        if not np.isfinite([minimum_k, maximum_k]).all() or minimum_k <= 0 or maximum_k < minimum_k:
            raise ValueError("invalid material temperature interval")
        bad = ((self.upper_temperature_k < minimum_k) |
               (self.upper_temperature_k > maximum_k))
        first = np.flatnonzero(bad.any(axis=(1, 2)))
        face_bad = np.sum([
            ((values < minimum_k) | (values > maximum_k)).sum(axis=1)
            for values in self.face_temperature_k.values()], axis=0)
        first_face = np.flatnonzero(face_bad)
        stop = min(float(self.time_s[first[0]]) if len(first) else np.inf,
                   float(self.face_time_s[first_face[0]]) if len(first_face) else np.inf)
        common_safe = self.time_s[self.time_s < stop]
        return {
            "admissible_temperature_k": [minimum_k, maximum_k],
            "upper_cell_temperature_min_k": float(self.upper_temperature_k.min()),
            "upper_cell_temperature_max_k": float(self.upper_temperature_k.max()),
            "invalid_upper_cell_counts": bad.sum(axis=(1, 2)).tolist(),
            "first_invalid_bulk_time_s": None if not len(first) else float(self.time_s[first[0]]),
            "all_bulk_frames_admissible": not bool(bad.any()),
            "exported_face_temperature_max_k": max(float(v.max()) for v in self.face_temperature_k.values()),
            "invalid_exported_face_patch_counts": face_bad.tolist(),
            "first_invalid_surface_time_s": None if not len(first_face) else float(self.face_time_s[first_face[0]]),
            "last_common_safe_bulk_time_s": None if not len(common_safe) else float(common_safe[-1]),
            "all_exported_temperatures_admissible": not bool(bad.any() or face_bad.any()),
            "lower_volume_basis": self.manifest["lower_volume_basis"],
            "qualification": "cell averages do not bound unresolved nodal or subcell maxima",
        }


def load_matched_history(directory: str | Path, *,
                         expected_source_sha256: str = MATCHED_256_SOURCE_SHA256,
                         lower_volume: str = "unavailable",
                         surface_symmetry_tolerance_k: float = 1e-4) -> MatchedThermalHistory:
    """Load exact-time temperature/Az/rate inputs from one native run.

    No EM solve or material validity inference occurs here. Every consumed file
    is verified against the externally pinned source_run.json manifest, and
    every native field against its own field-history record. No extrapolation,
    temporal interpolation, spatial projection or added bulk 300 K is allowed.
    ``mirror_with_surface_validation`` requests a symmetry benchmark, not a
    recovered lower volume solution.  Default behavior retains missing volume.
    """
    if lower_volume not in {"unavailable", "mirror_with_surface_validation"}:
        raise ValueError("lower_volume must explicitly select unavailable or mirror_with_surface_validation")
    if not np.isfinite(surface_symmetry_tolerance_k) or surface_symmetry_tolerance_k <= 0:
        raise ValueError("surface symmetry tolerance must be finite and positive")
    root = Path(directory).resolve()
    source_path = root / "source_run.json"
    source_hash = _sha256(source_path)
    if source_hash != expected_source_sha256:
        raise ValueError("source_run.json SHA256 differs from the requested matched run")
    run = json.loads(source_path.read_text())
    if run.get("schema") != "mhdlab-mfem-edge-current-run-v1" or run.get("passed") is not True:
        raise ValueError("input is not an accepted native edge-current run")
    hashes = run["outputs_sha256"]
    consumed = {"source_run.json": source_hash}

    def checked(name, additional_hash=None):
        path = (root / name).resolve()
        if not path.is_relative_to(root):
            raise ValueError("input reference escapes matched run directory")
        digest = _sha256(path)
        if digest != hashes.get(name) or (additional_hash is not None and digest != additional_hash):
            raise ValueError(f"SHA256 mismatch for {name}")
        consumed[name] = digest
        return path

    def metadata(name):
        return json.loads(checked(name).read_text())

    config = metadata("configuration.json")
    if config != run["configuration"]:
        raise ValueError("configuration differs from source-run record")
    expected_geometry = {"foil_width_m": 0.005, "foil_thickness_m": 0.0002,
                         "gap_m": 0.005, "current_length_m": 0.025}
    if any(config.get(name) != value for name, value in expected_geometry.items()):
        raise ValueError("input does not match the requested 5 mm / 0.2 mm / 5 mm / 25 mm geometry")
    if config.get("initial_temperature_k") != 300.0:
        raise ValueError("native producer hard-codes 300 K bulk reference; configuration must agree")
    nx, ny = config["width_cells"], config["thickness_cells"]
    if type(nx) is not int or type(ny) is not int or min(nx, ny) < 1:
        raise ValueError("native cell counts must be positive integers")
    width, thickness, gap = config["foil_width_m"], config["foil_thickness_m"], config["gap_m"]
    bulk_meta = metadata("foil_snapshots.json")
    if (bulk_meta.get("schema") != "mhdlab-mfem-foil-snapshots-v1" or
            bulk_meta.get("sampling") != "cell_average" or
            bulk_meta.get("units", {}).get("temperature_k") != "K" or
            bulk_meta.get("geometry", {}).get("foil") != "upper" or
            bulk_meta["geometry"].get("depth_reference") != "inward_from_upper_gap_face"):
        raise ValueError("unsupported bulk sampling, units or electrode basis")
    bulk_path = checked("foil_snapshots.csv", bulk_meta["csv_sha256"])
    with bulk_path.open() as stream:
        if stream.readline().strip() != _BULK_COLUMNS:
            raise ValueError("unexpected native bulk columns")
        cells = np.loadtxt(stream, delimiter=",", ndmin=2)
    if cells.shape[1] != 8 or not np.isfinite(cells).all() or np.any(cells[:, 5] <= 0):
        raise ValueError("nonfinite or invalid native bulk values")
    times = _partition(bulk_meta["field_snapshot_times_s"], "bulk times")
    x_edges = np.linspace(-width / 2, width / 2, nx + 1)
    d_edges = np.linspace(0, thickness, ny + 1)
    # Producer emits MFEM element order, which is not row-major grid order.
    it = np.searchsorted(times, cells[:, 0])
    ix = np.rint((cells[:, 1] + width / 2) / (width / nx)).astype(int)
    iy = np.rint(cells[:, 3] / (thickness / ny)).astype(int)
    if np.any(it >= len(times)) or np.any(ix < 0) or np.any(ix >= nx) or np.any(iy < 0) or np.any(iy >= ny):
        raise ValueError("bulk row outside declared grid/time support")
    if not np.allclose(times[it], cells[:, 0], rtol=0, atol=1e-20):
        raise ValueError("bulk row time differs from declared native field time")
    coordinates = np.column_stack((x_edges[ix], x_edges[ix + 1], d_edges[iy], d_edges[iy + 1]))
    if not np.allclose(coordinates, cells[:, 1:5], rtol=0, atol=1e-14):
        raise ValueError("bulk cells do not tile the configured native grid")
    flat = (it * ny + iy) * nx + ix
    if len(flat) != len(times) * ny * nx or len(np.unique(flat)) != len(flat):
        raise ValueError("duplicate or missing native bulk cells")
    upper = np.empty((len(times), ny, nx))
    upper[it, iy, ix] = cells[:, 5]  # Already absolute; DO NOT add 300 K.
    fields = metadata("em_field_snapshots.json")
    if (fields.get("schema") != "mhdlab-mfem-Az-field-history-v1" or
            fields.get("potential_units") != "T m" or fields.get("rate_units") != "V/m"):
        raise ValueError("unsupported native electromagnetic units or schema")
    for key in ("mesh", "axes"):
        checked(fields[f"{key}_file"], fields[f"{key}_sha256"])
    field_times = np.array([entry["time_s"] for entry in fields["snapshots"]])
    if field_times.shape != times.shape or not np.allclose(field_times, times, rtol=0, atol=1e-20):
        raise ValueError("temperature and native fields have different times")
    native_fields = []
    for entry in fields["snapshots"]:
        result = dict(entry)
        for key in ("potential", "rate"):
            result[f"{key}_path"] = str(checked(entry[f"{key}_file"], entry[f"{key}_sha256"]))
        native_fields.append(result)

    faces = {}
    face_times = None
    for prefix in ("", "anode_"):
        for face in ("top", "bottom", "left", "right"):
            name = prefix + face
            meta = metadata(name + "_surface.json")
            if meta.get("schema") != "mhdlab-mfem-surface-samples-v1":
                raise ValueError("unsupported surface schema")
            basis = meta["temperature"]
            if (basis.get("value_basis") != "temperature_increment_kelvin" or
                    basis.get("initial_temperature_k") != 300.0):
                raise ValueError("native surface increments require explicit 300 K offset")
            current_times = _partition(meta["time_s"], "surface times")
            if face_times is None:
                face_times = current_times
            elif not np.array_equal(face_times, current_times):
                raise ValueError("electrode surface times differ")
            path = checked(meta["table"]["path"], meta["table"]["sha256"])
            with path.open() as stream:
                if stream.readline().strip() != _FACE_COLUMNS:
                    raise ValueError("unexpected surface columns")
                rows = np.loadtxt(stream, delimiter=",", ndmin=2)
            count = nx if face in {"top", "bottom"} else ny
            if rows.shape != (len(face_times) * count, 4) or not np.isfinite(rows).all():
                raise ValueError("incomplete or invalid surface table")
            fi = np.searchsorted(face_times, rows[:, 0])
            indices = rows[:, 2].astype(int)
            if (np.any(fi >= len(face_times)) or np.any(rows[:, 1] != 0) or
                    np.any(indices != rows[:, 2]) or np.any(indices < 0) or np.any(indices >= count)):
                raise ValueError("invalid surface indices")
            if not np.array_equal(face_times[fi], rows[:, 0]):
                raise ValueError("surface table timestamp differs from metadata")
            if len(np.unique(fi * count + indices)) != len(rows):
                raise ValueError("duplicate surface samples")
            temperatures = np.empty((len(face_times), count))
            temperatures[fi, indices] = rows[:, 3] + 300.0
            if np.any(temperatures <= 0):
                raise ValueError("invalid absolute surface temperature")
            faces[name] = _readonly(temperatures)
    if not np.isin(times, face_times).all():
        raise ValueError("native bulk times absent from electrode faces")
    symmetry = {
        "upper_top_lower_bottom_k": float(np.max(np.abs(faces["top"] - faces["anode_bottom"]))),
        "upper_bottom_lower_top_k": float(np.max(np.abs(faces["bottom"] - faces["anode_top"]))),
        "left_reflected_depth_k": float(np.max(np.abs(faces["left"] - faces["anode_left"][:, ::-1]))),
        "right_reflected_depth_k": float(np.max(np.abs(faces["right"] - faces["anode_right"][:, ::-1]))),
    }
    lower = None
    if lower_volume == "mirror_with_surface_validation":
        if max(symmetry.values()) > surface_symmetry_tolerance_k:
            raise ValueError("independent lower surface histories fail requested mirror tolerance")
        lower = _readonly(upper[:, ::-1, :].copy())
    manifest = {
        "schema": "mhdlab-matched-thermal-mechanics-input-v1",
        "directory": str(root), "source_run_sha256": source_hash,
        "files_sha256": consumed, "geometry_m": expected_geometry,
        "time_count": len(times), "upper_shape": list(upper.shape),
        "upper_volume_basis": "native upper cell averages, absolute kelvin",
        "lower_volume_basis": ("unavailable; native lower volume was not exported" if lower is None else
                               "explicit reflection of native upper cell averages; lower surface symmetry checked; lower bulk unverified"),
        "surface_basis": "all eight directly exported faces, increment plus explicit 300 K",
        "surface_symmetry_max_abs_k": symmetry,
        "surface_symmetry_tolerance_k": surface_symmetry_tolerance_k,
        "native_fields_verified": True, "interpolation": "none; exact native matched times",
        "mechanical_supports": "unspecified; electrical feed/short and thermal insulation do not prescribe displacement or traction",
        "axial_representation": "thermal/EM cross-section invariant along 25 mm; no axial temperature resolution",
        "short_z_m": -0.0125, "feed_z_m": 0.0125, "observation_z_m": 0.0,
        "phase_change_in_source": False, "physical_material_prediction": False,
    }
    return MatchedThermalHistory(
        _readonly(times), _readonly(x_edges), _readonly(d_edges),
        _readonly(gap / 2 + d_edges), _readonly(-gap / 2 - d_edges[::-1]),
        _readonly(upper), lower, tuple(native_fields), _readonly(face_times), faces, manifest)
