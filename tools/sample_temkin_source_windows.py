#!/usr/bin/env python3
"""Export paired launch packets conditioned on original source-emission windows."""
from __future__ import annotations

import argparse
from dataclasses import fields
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mhdlab.temkin_planar_io import load_temkin_planar_release_history
from mhdlab.temkin_source_window import STATUSES
from mhdlab.temkin_windowed_source import SourceWindowSelection, condition_temkin_planar_source
from mhdlab.temkin_windowed_boundary import sample_temkin_windowed_emission_packets

CODE_PATHS = (
    "tools/sample_temkin_source_windows.py", "mhdlab/temkin_windowed_source.py",
    "mhdlab/temkin_windowed_boundary.py", "mhdlab/temkin_source_window.py",
    "mhdlab/temkin_emission_window.py", "mhdlab/temkin_planar.py",
    "mhdlab/temkin_planar_io.py", "mhdlab/temkin_history.py", "mhdlab/temkin_release.py",
    "mhdlab/temkin_io.py", "mhdlab/thermal_release.py", "mhdlab/planar_release.py",
    "mhdlab/planar_release_io.py", "mhdlab/planar_boundary.py", "mhdlab/thermal_boundary.py",
    "mhdlab/surface_release_history.py", "mhdlab/flash_surface.py",
    "mhdlab/surface_sources.py", "mhdlab/hdf_integrity.py", "mhdlab/constants.py",
)
CONFIG_KEYS = {
    "schema", "parameter_status", "conditioning_provenance", "mass_status", "mass_provenance",
    "state_masses_kg", "selections", "log_mass_tolerance", "face_normal_coordinate_m",
    "outward_normal_sign", "packet_count", "seed",
}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def nonfinite(value):
        raise ValueError(f"nonfinite JSON constant: {value}")

    return json.loads(Path(path).read_bytes(), object_pairs_hook=pairs, parse_constant=nonfinite)


def configuration(path):
    raw = path.read_bytes()
    config = read_json(path)
    if (not isinstance(config, dict) or set(config) != CONFIG_KEYS
            or config["schema"] != "mhdlab-temkin-windowed-packet-config-v1"):
        raise ValueError("configuration requires exactly the mhdlab-temkin-windowed-packet-config-v1 fields")
    for name in ("parameter_status", "mass_status"):
        if not isinstance(config[name], str) or config[name] not in STATUSES:
            raise ValueError(f"{name} must declare a supported input label")
    for name in ("conditioning_provenance", "mass_provenance"):
        if not isinstance(config[name], str) or not config[name].strip():
            raise ValueError(f"{name} must be nonempty")
    for name, lower, upper in (("packet_count", 1, 1_000_000), ("seed", 0, 2**64 - 1)):
        value = config[name]
        if isinstance(value, bool) or not isinstance(value, int) or not lower <= value <= upper:
            raise ValueError(f"{name} must be an integer in [{lower}, {upper}]")
    if not isinstance(config["state_masses_kg"], dict) or not config["state_masses_kg"]:
        raise ValueError("state_masses_kg must map exact state IDs to positive masses")
    for state, mass in config["state_masses_kg"].items():
        if (not isinstance(state, str) or not state.strip() or isinstance(mass, bool)
                or not isinstance(mass, (int, float)) or not math.isfinite(mass) or mass <= 0):
            raise ValueError("each exact-state mass must be finite and positive")
    if not isinstance(config["selections"], list) or not 1 <= len(config["selections"]) <= 10_000:
        raise ValueError("selections must contain 1 to 10000 explicit original-interval windows")
    selections = []
    for row in config["selections"]:
        if not isinstance(row, dict) or set(row) != {"interval_index", "start_s", "stop_s"}:
            raise ValueError("each selection requires interval_index, start_s and stop_s")
        selections.append(SourceWindowSelection(**row))
    return raw, config, tuple(selections)


def source_summary(path, history, source_hash, status):
    if path is None:
        if status != "unvalidated":
            raise ValueError("without a hash-associated source summary, external parameters require status unvalidated")
        return None
    summary = read_json(path)
    if (not isinstance(summary, dict) or summary.get("schema") not in {
            "mhdlab-temkin-planar-run-v1", "mhdlab-joint-source-run-v1", "mhdlab-coverage-source-run-v1"}
            or summary.get("parameter_status") != status
            or summary.get("parameter_provenance") != history.parameter_provenance
            or not isinstance(summary.get("output_sha256"), dict)
            or source_hash not in summary["output_sha256"].values()):
        raise ValueError("source summary must hash-associate the exact source, parameter status and provenance")
    return {"sha256": sha(path), "parameter_status": status,
            "parameter_provenance": history.parameter_provenance, "schema": summary["schema"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--source-run-summary", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    try:
        code_hashes = {name: sha(ROOT / name) for name in CODE_PATHS}
        raw, config, selections = configuration(args.config)
        config_hash = hashlib.sha256(raw).hexdigest()
        source_hash = sha(args.source)
        summary_hash = sha(args.source_run_summary) if args.source_run_summary else None
        history = load_temkin_planar_release_history(args.source)
        association = source_summary(args.source_run_summary, history, source_hash, config["parameter_status"])
        source = condition_temkin_planar_source(
            history, selections, log_mass_tolerance=config["log_mass_tolerance"],
            status=config["parameter_status"], provenance=config["conditioning_provenance"],
        )
        generator = np.random.default_rng(config["seed"])
        draws = generator.random((config["packet_count"], 7))
        if np.any(draws <= 0) or np.any(draws >= 1):
            raise ValueError("generated variates are not strictly open; choose a different explicit seed")
        packets = sample_temkin_windowed_emission_packets(
            source, config["state_masses_kg"], draws,
            face_normal_coordinate_m=config["face_normal_coordinate_m"],
            outward_normal_sign=config["outward_normal_sign"],
        )
        manifest = packets.manifest()
        arrays = {field.name: getattr(packets, field.name) for field in fields(packets)
                  if isinstance(getattr(packets, field.name), np.ndarray)}
        if set(arrays) - set(manifest["array_units"]):
            raise ValueError("packet output is missing array units")
        payload = dict(arrays)
        payload.update(state_ids=np.asarray(packets.state_ids, dtype=str),
                       reservoir_names=np.asarray(packets.reservoir_names, dtype=str),
                       metadata_json=np.asarray(json.dumps(manifest, sort_keys=True, allow_nan=False)))
        report = {
            "schema": "mhdlab-temkin-windowed-packet-export-v1",
            "parameter_status": config["parameter_status"],
            "source_parameter_provenance": history.parameter_provenance,
            "conditioning_provenance": config["conditioning_provenance"],
            "mass_status": config["mass_status"], "mass_provenance": config["mass_provenance"],
            "input_source_sha256": source_hash, "input_config_sha256": config_hash,
            "input_summary_sha256": summary_hash, "source_run_association": association,
            "source_sha256": code_hashes, "conditioned_source": source.manifest(), "packets": manifest,
            "sampling": {"seed": config["seed"], "packet_count": config["packet_count"],
                         "bit_generator": type(generator.bit_generator).__name__, "open_uniform_columns": 7},
            "runtime": {"python": sys.version,
                        **{name: importlib.metadata.version(name) for name in ("numpy", "scipy", "h5py")}},
            "host_budget_accepted": False,
            "limits": [
                "Source-emission selections condition launch times. Earlier launches can still radiate during later photon recording gates.",
                "Interval origin, window start, local delay and log delay retain their separate meanings; rounded absolute time has its own precision flags.",
                "Source and mass status declarations preserve provenance and do not validate physical water parameters, calibration or state branching.",
                "Unequal packet weights and signed number/energy residuals require an explicit downstream acceptance and ledger policy.",
                "Sampled translation and expected 2*kB*T source energy are separate; no velocity adjustment or wall cooling is performed.",
                "Staged output assumes a coordinated single writer; ordinary existing paths are refused.",
            ],
        }
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".temkin-windowed-", dir=args.output_dir.parent) as stage:
            output = Path(stage) / "result"
            output.mkdir()
            (output / "input_config.json").write_bytes(raw)
            shutil.copyfile(args.source, output / "source_input.h5")
            if args.source_run_summary:
                shutil.copyfile(args.source_run_summary, output / "source_run_summary.json")
            np.savez_compressed(output / "boundary_packets.npz", **payload)
            np.save(output / "uniforms.npy", draws, allow_pickle=False)
            if (sha(args.config) != config_hash or sha(args.source) != source_hash
                    or sha(output / "source_input.h5") != source_hash
                    or (args.source_run_summary and (sha(args.source_run_summary) != summary_hash
                        or sha(output / "source_run_summary.json") != summary_hash))
                    or any(sha(ROOT / name) != digest for name, digest in code_hashes.items())):
                raise ValueError("input or implementation changed during export")
            report["output_sha256"] = {p.name: sha(p) for p in output.iterdir() if p.is_file()}
            (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during export; refusing replacement")
            os.rename(output, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, OverflowError, RuntimeError) as error:
        parser.error(str(error))
    print(json.dumps({"output_dir": str(args.output_dir), "parameter_status": config["parameter_status"],
                      "packet_count": config["packet_count"], "components": len(source.components),
                      "host_budget_accepted": False}))


if __name__ == "__main__":
    main()
