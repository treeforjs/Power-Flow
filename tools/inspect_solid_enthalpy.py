#!/usr/bin/env python3
"""Inspect an explicitly qualified source-evaluated solid enthalpy branch."""
from __future__ import annotations

import argparse
import csv
from dataclasses import fields
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mhdlab.solid_enthalpy import load_mills_2004_solid_enthalpy

CODE_PATHS = ("tools/inspect_solid_enthalpy.py", "mhdlab/solid_enthalpy.py")
UNITS = {"temperature_k": "K", "cp_j_kg_k": "J kg^-1 K^-1", "printed_enthalpy_j_kg": "J kg^-1",
         "delta_enthalpy_from_reference_j_kg": "J kg^-1", "enthalpy_j_kg": "J kg^-1", "roundtrip_temperature_k": "K"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def load_configuration(path):
    def nonfinite(value):
        raise ValueError(f"nonfinite JSON constant: {value}")
    config = json.loads(path.read_bytes(), object_pairs_hook=_pairs, parse_constant=nonfinite)
    if (not isinstance(config, dict) or set(config) != {"schema", "model", "temperatures_k"}
            or config["schema"] != "mhdlab-solid-enthalpy-inspection-config-v1"):
        raise ValueError("configuration requires exactly the mhdlab-solid-enthalpy-inspection-config-v1 fields")
    required = {"temperature_min_k", "temperature_max_k", "reference_temperature_k", "reference_convention",
                "expected_material_family", "expected_phase", "material_qualification", "phase_qualification"}
    if not isinstance(config["model"], dict) or set(config["model"]) != required:
        raise ValueError("model must explicitly supply every source, material, phase, domain and reference qualification")
    temperatures = config["temperatures_k"]
    if not isinstance(temperatures, list) or not temperatures or len(temperatures) > 100000:
        raise ValueError("temperatures_k must be a nonempty one-dimensional list within the inspection budget")
    # The qualified model checks numerical types, finiteness and supported range.
    return config


def plot_values(values, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    temperature = values["temperature_k"]
    axes[0].plot(temperature, values["cp_j_kg_k"], "o-")
    axes[0].set(xlabel="Temperature [K]", ylabel="Mass-specific Cp [J/(kg K)]")
    axes[1].plot(temperature, np.asarray(values["printed_enthalpy_j_kg"]) / 1000, "o-", label="Printed Eq. 11")
    axes[1].plot(temperature, np.asarray(values["enthalpy_j_kg"]) / 1000, "s--", label="Declared reference")
    axes[1].set(xlabel="Temperature [K]", ylabel="Mass-specific enthalpy [kJ/kg]")
    axes[1].legend()
    for axis in axes:
        axis.grid(alpha=.25)
    figure.suptitle("Mills evaluated solid correlation; illustrative sampling, no EOS selection")
    figure.savefig(path, dpi=160)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/mills_enthalpy_example.json")
    parser.add_argument("--manifest", type=Path, default=ROOT / "data/materials/mills_2004_source_manifest.json")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    try:
        config_hash, manifest_hash = sha(args.config), sha(args.manifest)
        code_hashes = {name: sha(ROOT / name) for name in CODE_PATHS}
        config = load_configuration(args.config)
        model = load_mills_2004_solid_enthalpy(args.manifest, **config["model"])
        result = model.evaluate(config["temperatures_k"])
        values = {field.name: np.asarray(getattr(result, field.name)).tolist() for field in fields(result)}
        if np.asarray(result.temperature_k).ndim != 1:
            raise ValueError("temperatures_k must be one-dimensional")
        values["roundtrip_temperature_k"] = np.asarray(model.temperature_from_enthalpy_k(result.enthalpy_j_kg)).tolist()
        report = {"schema": "mhdlab-solid-enthalpy-inspection-v1", "model_metadata": model.metadata,
                  "sample_selection_status": "illustrative temperatures and explicitly selected supported branch/reference",
                  "input_config_sha256": config_hash, "input_manifest_sha256": manifest_hash,
                  "source_sha256": code_hashes, "units": UNITS, "values": values,
                  "maximum_roundtrip_temperature_difference_k": float(np.max(np.abs(np.asarray(values["roundtrip_temperature_k"]) - result.temperature_k))),
                  "limits": ["Source-evaluated mass-specific Cp and enthalpy; no actual foil composition, condition or phase calibration.",
                             "Cp is not Cv; this is not EOS internal energy, volumetric heat capacity or a temperature-evolution solve.",
                             "Printed reference offset is retained separately from an explicitly requested re-zeroed reference.",
                             "No density, liquid/fusion branch, extrapolation, transport-property stitching or runtime host selection.",
                             "Staged output assumes a coordinated single writer; ordinary pre-existing paths are refused."]}
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".solid-enthalpy-", dir=args.output_dir.parent) as stage:
            output = Path(stage) / "result"
            output.mkdir()
            (output / "input_config.json").write_bytes(args.config.read_bytes())
            with (output / "enthalpy.csv").open("w", newline="") as stream:
                writer = csv.writer(stream, lineterminator="\n")
                writer.writerow(values)
                writer.writerows(zip(*values.values()))
            if args.plot:
                plot_values(values, output / "enthalpy.png")
            if (sha(args.config) != config_hash or sha(args.manifest) != manifest_hash
                    or any(sha(ROOT / name) != value for name, value in code_hashes.items())):
                raise ValueError("input or implementation changed during inspection")
            report["output_sha256"] = {p.name: sha(p) for p in output.iterdir() if p.is_file()}
            (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during inspection; refusing replacement")
            os.rename(output, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, OverflowError, RuntimeError) as exc:
        parser.error(str(exc))
    print(json.dumps({"output_dir": str(args.output_dir), "status": "evaluated_source_correlation",
                      "temperature_samples": len(values["temperature_k"])}))


if __name__ == "__main__":
    main()
