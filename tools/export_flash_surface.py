#!/usr/bin/env python3
"""Export an explicitly configured raw scalar plane to the SI interchange.

The manufactured demonstration uses the same JSON loader and extraction path
as external files. It exercises the interface without representing a FLASH run
or an inferred temperature at the physical material boundary.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mhdlab.flash_extract import FlashLayout, FlashSnapshot, PlanarFace, extract_surface_history
from mhdlab.flash_surface import write_surface_history


CONFIG_SCHEMA = "mhdlab-flash-extract-config-v1"
_REQUIRED = {"schema", "configuration_status", "layout", "face", "snapshots", "source_metadata"}
_OPTIONAL = {"description"}
_SOURCE_REQUIRED = {"source_kind", "run_id", "source_reference", "surface_region"}
_SOURCE_OPTIONAL = {"flash_version"}


@dataclass(frozen=True)
class ExtractionConfiguration:
    path: Path
    sha256: str
    status: str
    layout: FlashLayout
    face: PlanarFace
    snapshots: tuple[FlashSnapshot, ...]
    source_metadata: dict[str, str]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _reject_constant(value):
    raise ValueError(f"nonfinite JSON constant is not permitted: {value}")


def _object(value, name):
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a JSON object")
    return value


def _keys(value, required, optional, name):
    _object(value, name)
    missing, unknown = required - value.keys(), value.keys() - required - optional
    if missing:
        raise ValueError(f"{name} is missing required keys: {', '.join(sorted(missing))}")
    if unknown:
        raise ValueError(f"{name} has unknown keys: {', '.join(sorted(unknown))}")


def _construct(record, values, name):
    try:
        return record(**_object(values, name))
    except TypeError as error:
        raise ValueError(f"invalid {name}: {error}") from error


def load_configuration(path: str | Path) -> ExtractionConfiguration:
    """Resolve snapshot paths against the configuration's own directory."""
    path = Path(path).resolve()
    raw = path.read_bytes()
    values = json.loads(raw.decode("utf-8"), object_pairs_hook=_unique_object,
                        parse_constant=_reject_constant)
    _keys(values, _REQUIRED, _OPTIONAL, "configuration")
    if values["schema"] != CONFIG_SCHEMA:
        raise ValueError(f"schema must equal {CONFIG_SCHEMA!r}")
    status = values["configuration_status"]
    if status == "template":
        raise ValueError("template configuration cannot execute; replace placeholders and explicitly set configuration_status")
    if status not in {"configured", "manufactured"}:
        raise ValueError("configuration_status must be 'configured', 'manufactured', or non-executable 'template'")
    if "description" in values and not isinstance(values["description"], str):
        raise ValueError("description must be a string")
    metadata = _object(values["source_metadata"], "source_metadata")
    _keys(metadata, _SOURCE_REQUIRED, _SOURCE_OPTIONAL, "source_metadata")
    if any(not isinstance(value, str) for value in metadata.values()):
        raise ValueError("source_metadata values must be strings")
    if status == "manufactured" and metadata["source_kind"] != "synthetic":
        raise ValueError("manufactured configurations must declare source_kind='synthetic'")
    layout = _construct(FlashLayout, values["layout"], "layout")
    face = _construct(PlanarFace, values["face"], "face")
    entries = values["snapshots"]
    if not isinstance(entries, list) or len(entries) < 2:
        raise ValueError("snapshots must be a JSON list with at least two explicit records")
    snapshots = []
    for index, entry in enumerate(entries):
        entry = dict(_object(entry, f"snapshot {index}"))
        if not isinstance(entry.get("path"), str) or not entry["path"].strip():
            raise ValueError(f"snapshot {index} path must be an explicit nonempty string")
        candidate = Path(entry["path"])
        entry["path"] = (path.parent / candidate).resolve() if not candidate.is_absolute() else candidate.resolve()
        snapshots.append(_construct(FlashSnapshot, entry, f"snapshot {index}"))
    return ExtractionConfiguration(path, hashlib.sha256(raw).hexdigest(), status, layout, face,
                                   tuple(snapshots), dict(metadata))


def export_configuration(path: str | Path, output_dir: str | Path) -> dict:
    """Write the same validated export for both external and demo inputs."""
    configuration = load_configuration(path)
    output = Path(output_dir).resolve()
    surface_path = output / "surface.h5"
    manifest_path = output / "extraction_manifest.json"
    inputs = {configuration.path, *(snapshot.path for snapshot in configuration.snapshots)}
    if inputs & {surface_path, manifest_path}:
        raise ValueError("export destinations must not overwrite the configuration or a raw snapshot")
    history, manifest = extract_surface_history(configuration.snapshots, configuration.layout,
                                                configuration.face, **configuration.source_metadata)
    output.mkdir(parents=True, exist_ok=True)
    temporary_paths = []
    try:
        with tempfile.NamedTemporaryFile(dir=output, prefix=".surface-", suffix=".h5", delete=False) as stream:
            temporary_surface = Path(stream.name)
        temporary_paths.append(temporary_surface)
        write_surface_history(temporary_surface, history)
        with h5py.File(temporary_surface, "a") as handle:
            handle.attrs["extraction_config_path"] = str(configuration.path)
            handle.attrs["extraction_config_sha256"] = configuration.sha256
            handle.attrs["extraction_manifest_path"] = str(manifest_path)
        manifest["configuration"] = {"path": str(configuration.path), "sha256": configuration.sha256,
                                     "schema": CONFIG_SCHEMA, "configuration_status": configuration.status}
        manifest["surface_history"] = {"path": str(surface_path), "sha256": _sha256(temporary_surface),
                                       "file_size_bytes": temporary_surface.stat().st_size}
        manifest["extraction_manifest_path"] = str(manifest_path)
        with tempfile.NamedTemporaryFile(mode="w", dir=output, prefix=".manifest-", suffix=".json",
                                         encoding="utf-8", delete=False) as stream:
            temporary_manifest = Path(stream.name)
            temporary_paths.append(temporary_manifest)
            json.dump(manifest, stream, indent=2, allow_nan=False)
            stream.write("\n")
        os.replace(temporary_surface, surface_path)
        os.replace(temporary_manifest, manifest_path)
    finally:
        for temporary in temporary_paths:
            temporary.unlink(missing_ok=True)
    return manifest


def write_demonstration(output_dir: str | Path) -> Path:
    """Generate a small manufactured Cartesian block layout, never a run claim."""
    output = Path(output_dir).resolve()
    raw_dir = output / "manufactured_snapshots"
    raw_dir.mkdir(parents=True, exist_ok=True)
    edges = np.linspace(-.0025, .0025, 7)
    centers = .5 * (edges[:-1] + edges[1:])
    boxes = np.asarray([[[-.0025, 0], [-20e-6, 0]], [[0, .0025], [-20e-6, 0]]])
    snapshots = []
    for index, (time_ns, factor) in enumerate(((0, 0), (60, .5), (120, 1))):
        filename = f"manufactured_snapshot_{index:04d}.h5"
        profile = 300 + factor * (100 + 900 * (centers / .0025) ** 4)
        temperature = np.stack([profile.reshape(2, 3) - 10 * factor,
                                profile.reshape(2, 3)], axis=1)
        with h5py.File(raw_dir / filename, "w") as handle:
            handle.attrs["source_kind"] = "synthetic"
            handle.attrs["description"] = "Manufactured scalar blocks for exporter verification; no simulation run"
            handle.create_dataset("bounding box", data=boxes).attrs["units"] = "m"
            handle.create_dataset("node type", data=np.ones(2, dtype=np.int32))
            handle.create_dataset("manufactured_temperature", data=temperature).attrs["units"] = "K"
            handle.create_dataset("manufactured_material_fraction", data=np.ones_like(temperature)).attrs["units"] = "1"
        snapshots.append({"path": f"manufactured_snapshots/{filename}", "time_value": time_ns,
                          "time_unit": "ns", "time_scale_to_s": 1e-9, "time_offset_s": 0})
    values = {
        "schema": CONFIG_SCHEMA, "configuration_status": "manufactured",
        "description": "Manufactured scalar extraction fixture; temperatures are sampled 5 micrometres inside a declared face, not extrapolated to its boundary.",
        "layout": {"bounding_box_dataset": "bounding box", "node_type_dataset": "node type",
                   "temperature_dataset": "manufactured_temperature",
                   "material_fraction_dataset": "manufactured_material_fraction",
                   "bounding_box_axes": ["x", "y"], "scalar_axes": ["block", "y", "x"],
                   "physical_axes": {"x": "width", "y": "normal"}, "leaf_node_value": 1,
                   "nonleaf_node_values": [2], "coordinate_unit": "m", "coordinate_scale_to_m": 1,
                   "temperature_unit": "K", "temperature_scale_to_k": 1, "temperature_offset_k": 0,
                   "min_material_fraction": 1 - 1e-8, "geometry": "cartesian",
                   "data_location": "cell_centered_no_guard_cells",
                   "axis_signs": {"x": 1, "y": 1}, "axis_offsets_m": {"x": 0, "y": 0}},
        "face": {"x_edges_m": edges.tolist(), "face_normal_coordinate_m": 0,
                 "outward_sign": 1, "material_sample_depth_m": 5e-6,
                 "current_interval_m": [0, .025], "geometry_atol_m": 1e-15},
        "snapshots": snapshots,
        "source_metadata": {"source_kind": "synthetic", "run_id": "manufactured-scalar-extraction",
                            "source_reference": "Manufactured input generated alongside this configuration; no external run",
                            "surface_region": "declared manufactured planar face", "flash_version": ""},
    }
    configuration_path = output / "demo_config.json"
    configuration_path.write_text(json.dumps(values, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return configuration_path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--config", type=Path, help="explicit raw extraction JSON; relative snapshots resolve beside it")
    source.add_argument("--demo", action="store_true", help="write manufactured inputs and pass their JSON through the same exporter")
    parser.add_argument("--output-dir", type=Path, default=Path("runs/flash_extract_reference"))
    args = parser.parse_args(argv)
    try:
        configuration_path = write_demonstration(args.output_dir) if args.demo else args.config
        manifest = export_configuration(configuration_path, args.output_dir)
    except (OSError, ValueError, TypeError, OverflowError, UnicodeError) as error:
        parser.error(str(error))
    print(json.dumps({"source_kind": manifest["source_kind"],
                      "temperature_location": manifest["temperature_location"],
                      "material_sample_depth_m": manifest["material_sample_depth_m"],
                      "surface_history": manifest["surface_history"],
                      "extraction_manifest_path": manifest["extraction_manifest_path"]}, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
