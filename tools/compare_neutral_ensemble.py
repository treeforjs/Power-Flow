#!/usr/bin/env python3
"""Compare complete neutral forward cases with explicitly bounded figure readouts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "mhdlab/neutral_spectral_inference.py"
# A separately installed editable package can point at another checkout. This
# command loads its own source explicitly, without changing global installation.
SPEC = importlib.util.spec_from_file_location("_neutral_spectral_inference_local", MODULE)
INFERENCE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = INFERENCE
SPEC.loader.exec_module(INFERENCE)

CONFIG_SCHEMA = "mhdlab-neutral-ensemble-config-v1"
RUN_SCHEMA = "mhdlab-neutral-ensemble-run-v1"
CODE_PATHS = ("mhdlab/neutral_spectral_inference.py", "tools/compare_neutral_ensemble.py")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def load_config(path):
    def invalid(value):
        raise ValueError(f"nonfinite JSON constant {value!r}")
    config = json.loads(Path(path).read_bytes(), object_pairs_hook=_object, parse_constant=invalid)
    required = {"schema", "observation", "cases"}
    if not isinstance(config, dict) or not required <= set(config) or set(config)-required-{"calibration_scales", "manufactured_truth"}:
        raise ValueError("config requires schema, observation, cases; only calibration_scales and manufactured_truth are optional")
    if config["schema"] != CONFIG_SCHEMA:
        raise ValueError(f"configuration schema must be {CONFIG_SCHEMA}")
    if not isinstance(config["observation"], dict) or not isinstance(config["cases"], list):
        raise ValueError("observation must be an object and cases must be a list")
    observed = dict(config["observation"])
    observed["curves"] = tuple(INFERENCE.ReadoutEnvelope(**item) for item in observed["curves"])
    observation = INFERENCE.EnsembleObservation(**observed)
    cases = []
    for values in config["cases"]:
        case = dict(values)
        case["curves"] = tuple(INFERENCE.ProjectedCurve(**item) for item in case["curves"])
        cases.append(INFERENCE.NeutralForwardCase(**case))
    if "manufactured_truth" in config and (observation.status != "manufactured" or any(c.status != "manufactured" for c in cases)):
        raise ValueError("manufactured_truth is allowed only for manufactured observations and cases")
    return config, observation, cases


def write_outputs(output, config, result, input_hash, code_hashes, *, plot):
    (output / "input_config.json").write_text(json.dumps(config, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    with (output / "case_scale_comparisons.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = ["case_id", "hypothesis", "C_abs", "max_normalized_departure", "unbounded_departure", "readout_compatible"]
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(result["case_scale_comparisons"])
    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        cases = [case["case_id"] for case in result["case_metadata"]]
        scales = result["calibration"]["scales"]
        matrix = np.array([[np.nan if row["max_normalized_departure"] is None else row["max_normalized_departure"]
                            for row in result["case_scale_comparisons"] if row["case_id"] == case] for case in cases])
        figure, axes = plt.subplots(1, 2, figsize=(11, max(4.2, .3*len(cases))), constrained_layout=True)
        display = axes[0].imshow(np.log10(1+matrix), aspect="auto", interpolation="none", origin="lower")
        axes[0].set(yticks=range(len(cases)), yticklabels=cases, xlabel="Calibration scenario index", title="log10(1 + worst readout departure)")
        figure.colorbar(display, ax=axes[0], label="Compatibility score transform")
        for hypothesis in sorted({r["hypothesis"] for r in result["hypothesis_scale_summaries"]}):
            rows = [r for r in result["hypothesis_scale_summaries"] if r["hypothesis"] == hypothesis]
            axes[1].plot(scales, [np.nan if r["minimum_normalized_departure"] is None else r["minimum_normalized_departure"] for r in rows],
                         ".-", label=hypothesis)
        axes[1].axhline(1, color="gray", linestyle="--", label="Inside every envelope")
        axes[1].set(xscale="log", yscale="symlog", xlabel="C_abs: multiplier on nominal radiance readout", ylabel="Minimum worst readout departure",
                    title="Supplied finite cases; no confidence interval")
        axes[1].set_ylim(bottom=0)
        axes[1].legend(fontsize=7)
        figure.suptitle(f"Neutral ensemble: {result['observation_status']} / {result['mapping_status']}")
        figure.savefig(output / "neutral_ensemble_comparison.png", dpi=160)
        plt.close(figure)
    report = {"schema": RUN_SCHEMA, "input_config_sha256": input_hash, "source_sha256": code_hashes,
              "comparison": result,
              "output_sha256": {path.name: sha(path) for path in sorted(output.iterdir()) if path.is_file()}}
    if "manufactured_truth" in config:
        report["manufactured_truth_not_used_by_comparison"] = config["manufactured_truth"]
    (output / "comparison_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    try:
        input_hash = sha(args.config)
        code_hashes = {name: sha(ROOT/name) for name in CODE_PATHS}
        config, observation, cases = load_config(args.config)
        result = INFERENCE.compare_neutral_ensemble(observation, cases, config.get("calibration_scales"))
        if sha(args.config) != input_hash or any(sha(ROOT/name) != value for name, value in code_hashes.items()):
            raise ValueError("input or implementation changed during comparison")
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".neutral-ensemble-", dir=args.output_dir.parent) as directory:
            output = Path(directory)/"result"
            output.mkdir()
            report = write_outputs(output, config, result, input_hash, code_hashes, plot=args.plot)
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during comparison")
            os.rename(output, args.output_dir)
        print(json.dumps({"output_dir": str(args.output_dir), "diagnostics": report["comparison"]["diagnostics"]}, indent=2, allow_nan=False))
    except (OSError, ValueError, TypeError, KeyError, OverflowError) as error:
        parser.exit(2, f"neutral ensemble comparison: {error}\n")


if __name__ == "__main__":
    main()
