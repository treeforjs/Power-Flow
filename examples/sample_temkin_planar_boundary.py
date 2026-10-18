#!/usr/bin/env python3
"""Export paired 3D Temkin launches from an explicitly resolved planar source.

This runnable boundary example creates no trajectories or collisions. Face
coordinate and orientation are explicit; state masses must agree with source IDs.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mhdlab.temkin_planar_io import load_temkin_planar_release_history
from mhdlab.temkin_planar_boundary import sample_temkin_planar_emission_packets
from tools.sample_surface_boundary import load_mass_config
from tools.sample_temkin_boundary import _sha, write_packets


CODE_PATHS = (
    "mhdlab/temkin_planar.py", "mhdlab/temkin_planar_io.py", "mhdlab/temkin_planar_boundary.py",
    "mhdlab/planar_release.py", "mhdlab/planar_release_io.py", "mhdlab/planar_boundary.py",
    "examples/sample_temkin_planar_boundary.py",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--source-run-summary", type=Path)
    parser.add_argument("--masses", type=Path, required=True)
    parser.add_argument("--interval", type=int, required=True)
    parser.add_argument("--packets", type=int, required=True)
    parser.add_argument("--face-normal-coordinate-m", type=float, required=True)
    parser.add_argument("--normal-sign", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.seed < 0 or args.packets < 1:
        parser.error("--seed must be nonnegative and --packets positive")
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    try:
        source_hash, mass_hash = _sha(args.source), _sha(args.masses)
        config = load_mass_config(args.masses)
        history = load_temkin_planar_release_history(args.source)
        if config["outward_normal_sign"] != args.normal_sign:
            raise ValueError("normal-sign conflicts with the explicit mass/geometry configuration")
        actual_length = float(history.z_edges_m[-1] - history.z_edges_m[0])
        if not np.isclose(config["current_direction_length_m"], actual_length, rtol=8 * np.finfo(float).eps, atol=0):
            raise ValueError("mass configuration current length conflicts with the actual planar extent")
        source_run_record = None
        if args.source_run_summary is not None:
            digest = _sha(args.source_run_summary)
            record = json.loads(args.source_run_summary.read_bytes())
            if (not isinstance(record, dict) or not isinstance(record.get("output_sha256"), dict)
                    or record.get("schema") != "mhdlab-temkin-planar-run-v1"
                    or record["output_sha256"].get("temkin_planar_release.h5") != source_hash
                    or record.get("parameter_status") not in {"manufactured", "sourced", "unvalidated"}
                    or record.get("parameter_provenance") != history.parameter_provenance):
                raise ValueError("source-run summary must match the planar source hash and provenance")
            source_run_record = {"path": str(args.source_run_summary.resolve()), "sha256": digest,
                                 "parameter_status": record["parameter_status"],
                                 "parameter_provenance": record["parameter_provenance"]}
        rng = np.random.default_rng(args.seed)
        draws = rng.random((args.packets, 7))
        redraws = 0
        while np.any(draws == 0):
            selected = draws == 0
            redraws += int(np.count_nonzero(selected))
            draws[selected] = rng.random(int(np.count_nonzero(selected)))
        packets = sample_temkin_planar_emission_packets(history, config["state_masses_kg"], args.interval, draws,
                                                        face_normal_coordinate_m=args.face_normal_coordinate_m,
                                                        outward_normal_sign=args.normal_sign)
        if _sha(args.source) != source_hash or _sha(args.masses) != mass_hash:
            raise ValueError("a source or mass file changed during sampling")
        if source_run_record is not None and _sha(args.source_run_summary) != source_run_record["sha256"]:
            raise ValueError("source-run summary changed during sampling")
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".temkin-planar-packets-", dir=args.output_dir.parent) as stage:
            output = Path(stage) / "result"
            output.mkdir()
            summary = write_packets(output, packets, source_path=args.source, mass_path=args.masses,
                                    source_hash=source_hash, mass_hash=mass_hash, config=config, seed=args.seed,
                                    redraws=redraws, generator=type(rng.bit_generator).__name__,
                                    source_run_record=source_run_record, additional_code_paths=CODE_PATHS)
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during sampling; refusing replacement")
            os.rename(output, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, OverflowError, RuntimeError) as exc:
        parser.error(str(exc))
    print(json.dumps({"output_dir": str(args.output_dir), "packets": summary["packets"],
                      "components": packets.required_minimum_n,
                      "total_packet_particles": summary["total_packet_particles"]}, allow_nan=False))


if __name__ == "__main__":
    main()
