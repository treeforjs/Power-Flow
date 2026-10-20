#!/usr/bin/env python3
"""Inspect deterministic source moments associated with a saved packet export."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mhdlab.temkin_planar_io import load_temkin_planar_release_history
from mhdlab.temkin_windowed_boundary import (
    ALLOCATION_MODEL, SCHEMA as PACKET_SCHEMA, VELOCITY_MODEL, _ARRAY_UNITS,
)
from mhdlab.temkin_release import LAW_IDENTITY
from mhdlab.temkin_windowed_source import condition_temkin_planar_source
from mhdlab.windowed_joint_source_moments import windowed_joint_source_moments
from tools.sample_temkin_source_windows import (
    CODE_PATHS as PRODUCER_CODE_PATHS, configuration, read_json, sha, source_summary,
)

CODE_PATHS = (
    "tools/inspect_windowed_source_moments.py", "mhdlab/window_time_moments.py",
    "mhdlab/windowed_joint_source_moments.py", "mhdlab/emission_time_moments.py",
    "mhdlab/temkin_emission_time.py", "tools/sample_temkin_source_windows.py",
    "mhdlab/temkin_windowed_source.py", "mhdlab/temkin_windowed_boundary.py",
    "mhdlab/temkin_source_window.py", "mhdlab/temkin_emission_window.py",
    "mhdlab/temkin_planar.py", "mhdlab/temkin_planar_io.py", "mhdlab/temkin_history.py",
    "mhdlab/temkin_release.py", "mhdlab/temkin_io.py", "mhdlab/thermal_release.py",
    "mhdlab/planar_release.py", "mhdlab/planar_release_io.py", "mhdlab/planar_boundary.py",
    "mhdlab/thermal_boundary.py", "mhdlab/surface_release_history.py", "mhdlab/flash_surface.py",
    "mhdlab/surface_sources.py", "mhdlab/hdf_integrity.py", "mhdlab/constants.py",
)
REQUIRED_FILES = {"input_config.json", "source_input.h5", "boundary_packets.npz", "uniforms.npy"}


def decode_packet_metadata(value):
    def pairs(items):
        result = {}
        for key, item in items:
            if key in result:
                raise ValueError(f"duplicate packet metadata JSON key: {key}")
            result[key] = item
        return result

    def nonfinite(value):
        raise ValueError(f"nonfinite packet metadata JSON constant: {value}")

    if value.shape != () or value.dtype.kind != "U":
        raise ValueError("packet metadata must be a scalar Unicode JSON string")
    return json.loads(str(value), object_pairs_hook=pairs, parse_constant=nonfinite)


def validate_packet_structure(path, metadata, source, packet_count):
    identities = {"schema": PACKET_SCHEMA, "law_identity": LAW_IDENTITY,
                  "allocation_model": ALLOCATION_MODEL, "velocity_model": VELOCITY_MODEL}
    if (any(metadata.get(name) != value for name, value in identities.items())
            or metadata.get("array_units") != _ARRAY_UNITS
            or type(metadata.get("packets")) is not int or metadata["packets"] != packet_count
            or type(metadata.get("required_minimum_n")) is not int
            or metadata["required_minimum_n"] != len(source.components)):
        raise ValueError("packet metadata units, identities or dimensions differ from the pinned packet contract")
    with np.load(path, allow_pickle=False) as packets:
        if set(packets.files) != set(_ARRAY_UNITS) | {"metadata_json", "state_ids", "reservoir_names"}:
            raise ValueError("packet arrays differ from the pinned packet contract")
        if decode_packet_metadata(packets["metadata_json"]) != metadata:
            raise ValueError("packet metadata disagrees with its report")
        for name in _ARRAY_UNITS:
            values = packets[name]
            if name in {"x_edges_m", "z_edges_m", "time_edges_s"}:
                expected = getattr(source, name)
                if not np.array_equal(values, expected):
                    raise ValueError(f"packet source edges disagree: {name}")
                shape = expected.shape
            elif (name.startswith("component_") and name != "component_index") or name == "packet_counts":
                shape = (len(source.components),)
            elif name in {"resolved_position_m", "uniform_width_position_m", "velocity_m_s", "momentum_kg_m_s"}:
                shape = (packet_count, 3)
            else:
                shape = (packet_count,)
            if name.endswith("_index") or name == "packet_counts":
                valid_kind = values.dtype.kind in "iu"
            elif name in {"global_time_precision_qualified", "component_probability_underflowed", "component_probability_representation_qualified"}:
                valid_kind = values.dtype.kind == "b"
            else:
                valid_kind = values.dtype.kind == "f"
            if values.shape != shape or not valid_kind or not np.all(np.isfinite(values)):
                raise ValueError(f"packet array shape/type/value disagrees with the pinned contract: {name}")
        for name in ("state_ids", "reservoir_names"):
            if packets[name].shape != (packet_count,) or packets[name].dtype.kind not in "US":
                raise ValueError(f"packet state/name array shape/type disagrees: {name}")


def associated_source(directory):
    report_path = directory / "report.json"
    report_hash = sha(report_path)
    report = read_json(report_path)
    if (not isinstance(report, dict) or report.get("schema") != "mhdlab-temkin-windowed-packet-export-v1"
            or report.get("host_budget_accepted") is not False):
        raise ValueError("input must be a source-window packet export with no accepted host transaction")
    hashes = report.get("output_sha256")
    if (not isinstance(hashes, dict) or set(hashes) not in (REQUIRED_FILES, REQUIRED_FILES | {"source_run_summary.json"})
            or any(not isinstance(value, str) or len(value) != 64 for value in hashes.values())):
        raise ValueError("packet report must bind exactly the expected local export files")
    for name, digest in hashes.items():
        path = directory / name
        if path.is_symlink() or not path.is_file() or sha(path) != digest:
            raise ValueError(f"packet export file is missing, indirect or changed: {name}")
    expected_producer = {name: sha(ROOT / name) for name in PRODUCER_CODE_PATHS}
    if report.get("source_sha256") != expected_producer:
        raise ValueError("packet producer implementation differs from the pinned current contract")
    raw, config, selections = configuration(directory / "input_config.json")
    source_hash = hashes["source_input.h5"]
    if (report.get("input_source_sha256") != source_hash
            or report.get("input_config_sha256") != hashlib.sha256(raw).hexdigest()
            or report.get("parameter_status") != config["parameter_status"]
            or report.get("conditioning_provenance") != config["conditioning_provenance"]
            or report.get("mass_status") != config["mass_status"] or report.get("mass_provenance") != config["mass_provenance"]):
        raise ValueError("packet report disagrees with its exact input files or status/provenance")
    summary_path = directory / "source_run_summary.json" if "source_run_summary.json" in hashes else None
    if report.get("input_summary_sha256") != (hashes["source_run_summary.json"] if summary_path else None):
        raise ValueError("packet report disagrees with its source-summary association")
    history = load_temkin_planar_release_history(directory / "source_input.h5")
    association = source_summary(summary_path, history, source_hash, config["parameter_status"])
    if report.get("source_run_association") != association or report.get("source_parameter_provenance") != history.parameter_provenance:
        raise ValueError("source parameter provenance or summary association differs")
    source = condition_temkin_planar_source(
        history, selections, log_mass_tolerance=config["log_mass_tolerance"],
        status=config["parameter_status"], provenance=config["conditioning_provenance"],
    )
    manifest = source.manifest()
    packet_metadata = report.get("packets")
    if (not isinstance(packet_metadata, dict) or report.get("conditioned_source") != manifest
            or packet_metadata.get("source") != manifest or packet_metadata.get("host_budget_accepted") is not False):
        raise ValueError("saved conditioned source differs from the original canonical history selections")
    validate_packet_structure(directory / "boundary_packets.npz", packet_metadata, source, config["packet_count"])
    return source, report, {**hashes, "report.json": report_hash}


def refinement_changes(coarse, refined):
    changes = {}
    for hypothesis in ("resolved", "matched_uniform_width"):
        changes[hypothesis] = {}
        for name in coarse["units"]:
            a, b = coarse[hypothesis].get(name), refined[hypothesis].get(name)
            if a is None or b is None:
                continue
            if isinstance(a, bool) or isinstance(b, bool) or not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
                continue
            difference = math.fsum((float(b), -float(a)))
            relative = difference / abs(b) if b else None
            reason = "refined_value_is_zero" if b == 0 else "represented"
            if relative is not None and not math.isfinite(relative):
                relative, reason = None, "relative_change_overflow"
            changes[hypothesis][name] = {"signed_change": difference, "absolute_change": abs(difference),
                                         "relative_to_refined": relative, "relative_representation": reason}
    return {"coarse_order": coarse["quadrature_order"], "refined_order": refined["quadrature_order"],
            "changes": changes, "interpretation": "Changes between finite quadrature estimates; not certified error bounds or experimental uncertainty."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet-export", type=Path, required=True)
    parser.add_argument("--orders", type=int, nargs="+", default=[8, 16, 32])
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    if (not 2 <= len(args.orders) <= 5 or any(not 2 <= value <= 4096 for value in args.orders)
            or any(a >= b for a, b in zip(args.orders, args.orders[1:]))):
        parser.error("--orders requires 2 to 5 strictly increasing integers in [2,4096]")
    try:
        code_hashes = {name: sha(ROOT / name) for name in CODE_PATHS}
        source, packet_report, input_hashes = associated_source(args.packet_export)
        estimates = [windowed_joint_source_moments(source, quadrature_order=order, include_component_records=True)
                     for order in args.orders]
        refinements = [refinement_changes(a, b) for a, b in zip(estimates, estimates[1:])]
        report = {
            "schema": "mhdlab-windowed-source-moment-inspection-v1",
            "parameter_status": source.status, "parameter_provenance": source.parameter_provenance,
            "conditioning_provenance": source.provenance,
            "mass_status": packet_report["mass_status"], "mass_provenance": packet_report["mass_provenance"],
            "input_export_sha256": input_hashes, "source_sha256": code_hashes, "orders": args.orders,
            "refinements": refinements, "host_budget_accepted": False,
            "runtime": {"python": sys.version, **{name: importlib.metadata.version(name) for name in ("numpy", "scipy", "h5py")}},
            "limits": [
                "Reuses the exact original packet export and validates its canonical source selections; no source history evolution or packet sampling is run.",
                "Moments condition original SOURCE emission windows, separate from photon recording gates or diagnostic timing origins.",
                "Probability quadrature refinements are changes in finite estimates; numerical source-law qualification and sampling uncertainty are separate.",
                "Conditional-mean translational energy variation excludes thermal scatter; full kinetic-energy variance includes within-component flux-Maxwellian scatter.",
                "Source and mass status declarations retain provenance without validating physical water parameters, coverage or calibration.",
                "No host transaction, wall feedback, gap, collision, radiation or charged-particle evolution is performed.",
                "Staged output assumes a coordinated single writer; ordinary existing paths are refused.",
            ],
        }
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".windowed-moments-", dir=args.output_dir.parent) as stage:
            output = Path(stage) / "result"
            output.mkdir()
            for name in input_hashes:
                shutil.copyfile(args.packet_export / name, output / ("packet_export_" + name))
            for order, estimate in zip(args.orders, estimates):
                (output / f"moments_order_{order:04d}.json").write_text(json.dumps(estimate, indent=2, allow_nan=False) + "\n")
            if (any(sha(args.packet_export / name) != digest or sha(output / ("packet_export_" + name)) != digest
                    for name, digest in input_hashes.items())
                    or any(sha(ROOT / name) != digest for name, digest in code_hashes.items())):
                raise ValueError("input export or implementation changed during inspection")
            report["output_sha256"] = {path.name: sha(path) for path in output.iterdir() if path.is_file()}
            (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during inspection; refusing replacement")
            os.rename(output, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, OverflowError, RuntimeError, zipfile.BadZipFile) as error:
        parser.error(str(error))
    print(json.dumps({"output_dir": str(args.output_dir), "parameter_status": source.status,
                      "orders": args.orders, "components": len(source.components), "host_budget_accepted": False}))


if __name__ == "__main__":
    main()
