#!/usr/bin/env python3
"""Inspect a reconstructed published current with explicit time-coordinate choice.

The default input is the recovered shot 15210 Fig. C.11 polyline and its version-1
sidecar. Select --preserve-source-axis or --time-offset-s explicitly. An offset
is a chosen coordinate translation, not evidence of a recovered physical trigger.
This reader checks the CSV hash, declared SI columns, source-unit conversions,
display bounds and legend flags. It does not alter the waveform, infer a circuit,
resample, smooth, differentiate, extrapolate or run a solver.

JSON is printed to stdout. Files are written only through explicit --json-out
and --plot paths. Optional --integrate START_S STOP_S uses selected SI times
and reports exact integrals of the piecewise-linear interpolation, not measured
charge or dissipated heat.

Example:
  python tools/inspect_current_drive.py --preserve-source-axis \
      --json-out runs/current_inspection/summary.json \
      --plot runs/current_inspection/current.png
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mhdlab.drive_trace import DriveTrace, load_current_csv


def _mapping(value, name: str) -> dict:
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be an object")
    return value


def _text(value, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def _number(value, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    return float(value)


def _integer(value, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _pair(value, name: str, *, allow_equal: bool = False) -> list[float]:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{name} must contain two numeric bounds")
    result = [_number(v, name) for v in value]
    if result[1] < result[0] or (not allow_equal and result[1] == result[0]):
        raise ValueError(f"{name} bounds must increase")
    return result


def _near(actual: float, expected: float, name: str, *, absolute: float = 1e-8) -> None:
    # The recovery CSV rounds native coordinate conversions to 12 significant
    # digits. This tolerance checks representation, not physical uncertainty.
    if not math.isclose(actual, expected, rel_tol=2e-11, abs_tol=absolute):
        raise ValueError(f"{name} disagrees with the declared recovery provenance")


def _csv_number(row: dict, key: str) -> float:
    try:
        result = float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"CSV column {key} must contain finite numbers") from exc
    if not math.isfinite(result):
        raise ValueError(f"CSV column {key} must contain finite numbers")
    return result


def _reject_json_constant(value: str):
    raise ValueError(f"Nonfinite JSON constant {value} is not supported")


def inspect_current(
    csv_path: Path, provenance_path: Path, *, mode: str, offset_s: float,
    integrate: tuple[float, float] | None = None,
) -> tuple[dict, DriveTrace, dict]:
    """Read the known recovery schema and return a summary plus plotting data."""
    if mode not in {"preserve_source_axis", "explicit_offset"}:
        raise ValueError("An explicit time-coordinate mode is required")
    offset_s = _number(offset_s, "time offset")
    if mode == "preserve_source_axis" and offset_s != 0:
        raise ValueError("Preserving the source axis requires a zero offset")
    provenance_bytes = provenance_path.read_bytes()
    metadata = _mapping(json.loads(provenance_bytes, parse_constant=_reject_json_constant), "provenance")
    if type(metadata.get("schema_version")) is not int or metadata["schema_version"] != 1:
        raise ValueError("Only recovery provenance schema_version 1 is supported")
    if metadata.get("status") != "reconstructed_published_vector_polyline":
        raise ValueError("Expected reconstructed_published_vector_polyline provenance")
    source = _mapping(metadata.get("source"), "source")
    experiment = _mapping(metadata.get("experiment"), "experiment")
    observable = _mapping(metadata.get("observable"), "observable")
    output = _mapping(metadata.get("output"), "output")
    vector = _mapping(metadata.get("vector_path"), "vector_path")
    calibration = _mapping(metadata.get("calibration"), "calibration")
    time_axis = _mapping(calibration.get("time_axis"), "time_axis")
    current_axis = _mapping(calibration.get("current_axis"), "current_axis")
    display = _mapping(calibration.get("digitization_bounds"), "digitization_bounds")
    gate = _mapping(metadata.get("VUV_gate"), "VUV_gate")
    if (observable.get("quantity") != "published_average_current"
            or observable.get("original_units") != "MA" or observable.get("output_units") != "A"
            or time_axis.get("label") != "Time (ns)" or current_axis.get("legend") != "Average Current [MA]"):
        raise ValueError("Unsupported or inconsistent source units; expected ns/MA source and SI time_s/current_A columns")
    if _number(current_axis.get("A_per_MA"), "declared A_per_MA scale") != 1e6:
        raise ValueError("Declared current scale must be 1e6 A per MA")
    if observable.get("raw_oscilloscope_data") is not False:
        raise ValueError("Reconstructed current provenance cannot claim raw oscilloscope data")
    required = {
        "source_vertex_index", "page_x_pt", "page_y_pt", "time_ns", "time_s", "current_MA", "current_A",
        "time_digitization_lower_ns", "time_digitization_upper_ns",
        "current_digitization_lower_A", "current_digitization_upper_A", "occluded_by_legend",
    }
    declared_columns = output.get("columns")
    if (not isinstance(declared_columns, list) or not all(isinstance(v, str) and v for v in declared_columns)
            or len(set(declared_columns)) != len(declared_columns) or not required.issubset(declared_columns)):
        raise ValueError("Provenance must declare the known SI, display-bound and legend-flag columns")
    expected_hash = _text(output.get("csv_sha256"), "csv_sha256")
    if not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
        raise ValueError("csv_sha256 must be a lowercase SHA-256 digest")
    csv_bytes = csv_path.read_bytes()
    digest = hashlib.sha256(csv_bytes).hexdigest()
    if digest != expected_hash:
        raise ValueError("CSV hash mismatch with the recovery provenance")
    reader = csv.DictReader(io.StringIO(csv_bytes.decode("utf-8-sig")), strict=True)
    rows = list(reader)
    if reader.fieldnames != declared_columns:
        raise ValueError("CSV columns disagree with the recovery provenance")
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError("Malformed CSV row")
    samples = _integer(output.get("samples"), "samples")
    if len(rows) != samples or samples < 2 or _integer(vector.get("retained_vertices"), "retained_vertices") != samples:
        raise ValueError("CSV sample count disagrees with the retained vertices")
    time_bound_ns = _number(display.get("time_plus_minus_ns"), "time display bound")
    current_bound_a = _number(display.get("current_plus_minus_A"), "current display bound")
    if time_bound_ns <= 0 or current_bound_a <= 0:
        raise ValueError("Display bounds must be positive")
    legend = vector.get("later_legend_bounding_rectangle_pt")
    if not isinstance(legend, list) or len(legend) != 4:
        raise ValueError("Expected the source legend rectangle")
    legend = [_number(v, "legend rectangle") for v in legend]
    if not legend[0] < legend[2] or not legend[1] < legend[3]:
        raise ValueError("Legend rectangle is invalid")
    time_s, current_a, flags, indices, low_a, high_a = [], [], [], [], [], []
    for row in rows:
        t_s, t_ns = _csv_number(row, "time_s"), _csv_number(row, "time_ns")
        current, current_ma = _csv_number(row, "current_A"), _csv_number(row, "current_MA")
        _near(t_s * 1e9, t_ns, "CSV SI/source time units", absolute=1e-9)
        _near(current_ma * 1e6, current, "CSV SI/source current units")
        t_low, t_high = _csv_number(row, "time_digitization_lower_ns"), _csv_number(row, "time_digitization_upper_ns")
        i_low, i_high = _csv_number(row, "current_digitization_lower_A"), _csv_number(row, "current_digitization_upper_A")
        if not t_low < t_ns < t_high or not i_low < current < i_high:
            raise ValueError("CSV display bounds do not enclose the source vertex")
        for actual, expected, label in [(t_low, t_ns - time_bound_ns, "lower time bound"),
                                        (t_high, t_ns + time_bound_ns, "upper time bound"),
                                        (i_low, current - current_bound_a, "lower current bound"),
                                        (i_high, current + current_bound_a, "upper current bound")]:
            _near(actual, expected, label)
        if row["occluded_by_legend"] not in {"0", "1"}:
            raise ValueError("CSV legend flags must be exactly 0 or 1")
        flag = int(row["occluded_by_legend"])
        x, y = _csv_number(row, "page_x_pt"), _csv_number(row, "page_y_pt")
        if flag != int(legend[0] <= x <= legend[2] and legend[1] <= y <= legend[3]):
            raise ValueError("CSV legend flag disagrees with the source rectangle")
        try:
            index = int(row["source_vertex_index"])
        except ValueError as exc:
            raise ValueError("Source vertex indices must be nonnegative integers") from exc
        if index < 0 or (indices and index <= indices[-1]):
            raise ValueError("Source vertex indices must strictly increase")
        time_s.append(t_s); current_a.append(current); flags.append(flag); indices.append(index)
        low_a.append(i_low); high_a.append(i_high)
    if sum(flags) != _integer(vector.get("vertices_behind_legend"), "vertices_behind_legend"):
        raise ValueError("CSV legend count disagrees with the provenance")
    for actual, expected in zip([time_s[0] * 1e9, time_s[-1] * 1e9], _pair(output.get("time_extent_ns"), "time extent")):
        _near(actual, expected, "source time extent")
    for actual, expected in zip([min(current_a), max(current_a)], _pair(output.get("current_extent_A"), "current extent", allow_equal=True)):
        _near(actual, expected, "source current extent")
    figure = _text(source.get("figure"), "source figure")
    facility = _text(experiment.get("facility"), "facility")
    shot = experiment.get("shot")
    if isinstance(shot, bool) or not isinstance(shot, (str, int)) or not str(shot).strip():
        raise ValueError("Shot must have an explicit identity")
    title = _text(source.get("title"), "source title")
    trace = load_current_csv(csv_path, time_column="time_s", current_column="current_A",
                             time_scale_to_s=1.0, current_scale_to_a=1.0,
                             diagnostic_time_offset_s=offset_s, source_kind="reconstructed",
                             shot_id=f"{facility} {shot}", source_reference=title,
                             figure_reference=f"Figure {figure}, {current_axis['legend']}")
    if trace.csv_import_metadata["sha256"] != digest:
        raise ValueError("CSV changed during inspection")
    if gate.get("kind") != "published_gate_rectangle_not_breakdown":
        raise ValueError("Expected a separately identified published VUV gate rectangle")
    gate_source = np.array(_pair(gate.get("plotted_time_bounds_ns"), "VUV gate time bounds")) * 1e-9
    gate_selected = gate_source + offset_s
    if not np.all(np.isfinite(gate_selected)) or not gate_selected[0] < gate_selected[1]:
        raise ValueError("Selected offset cannot represent the gate interval")
    peak = int(np.argmax(trace.current_a))
    negative = trace.current_a < 0
    breakdown_ns = gate.get("actual_breakdown_time_ns")
    breakdown_s = None if breakdown_ns is None else _number(breakdown_ns, "source breakdown time") * 1e-9
    limits = metadata.get("limits", [])
    if not isinstance(limits, list) or not all(isinstance(value, str) for value in limits):
        raise ValueError("Provenance limits must be a list of strings")
    report = {
        "schema_version": 1, "scope": "inspection_of_reconstructed_published_current",
        "source_kind": "reconstructed",
        "input": {"csv": str(csv_path.resolve()), "csv_sha256": digest,
                  "provenance": str(provenance_path.resolve()), "provenance_sha256": hashlib.sha256(provenance_bytes).hexdigest(),
                  "source": source, "experiment": experiment,
                  "source_current_legend": current_axis["legend"], "source_time_axis_label": time_axis["label"],
                  "selected_columns": {"time": "time_s", "current": "current_A"},
                  "CSV_units": {"time": "s", "current": "A"}},
        "alignment": {"mode": mode, "applied_time_offset_s": offset_s,
                      "rule": "selected_time_s = source_CSV_time_s + chosen_offset_s",
                      "source_time_origin_reference": observable.get("time_origin_reference"),
                      "physical_trigger_alignment_established": False,
                      "offset_interpretation": "chosen coordinate translation; no recovered physical trigger is asserted"},
        "waveform": {"samples": samples, "source_time_support_s": [time_s[0], time_s[-1]],
                     "selected_time_support_s": [float(trace.time_s[0]), float(trace.time_s[-1])],
                     "minimum_current_A": float(min(trace.current_a)), "maximum_current_A": float(max(trace.current_a)),
                     "peak": {"interpretation": "largest_retained_vertex_of_piecewise_linear_display_curve",
                              "current_A": float(trace.current_a[peak]), "source_time_s": time_s[peak],
                              "selected_time_s": float(trace.time_s[peak]), "source_vertex_index": indices[peak],
                              "occluded_by_legend": bool(flags[peak])},
                     "negative_values": {"count": int(np.count_nonzero(negative)),
                                         "minimum_current_A": float(min(trace.current_a)) if np.any(negative) else None},
                     "legend_covered_vertices": sum(flags), "signed_values_preserved": True,
                     "resampling_or_smoothing_applied": False, "extrapolation_supplied": False},
        "display_bounds": {"time_plus_minus_s": time_bound_ns * 1e-9,
                           "current_plus_minus_A": current_bound_a, "physical_uncertainty": False,
                           "interpretation": _text(display.get("interpretation"), "display-bound interpretation")},
        "VUV_gate": {"kind": gate["kind"], "source_time_bounds_s": gate_source.tolist(),
                     "selected_time_bounds_s": gate_selected.tolist(),
                     "source_page_rectangle_pt": gate.get("page_rectangle_pt"),
                     "actual_breakdown_time_s": breakdown_s,
                     "selected_breakdown_time_s": None if breakdown_s is None else breakdown_s + offset_s,
                     "exact_trigger_reference": gate.get("exact_trigger_reference"),
                     "table_4_1_source_time_ns": gate.get("table_4_1_published_acquisition_time_ns"),
                     "table_4_1_relative_to_breakdown": gate.get("table_4_1_relative_to_breakdown"),
                     "gate_edge_meaning": gate.get("gate_edge_meaning")},
        "limits": limits + [
            "This inspection checks internal file/provenance consistency; it does not independently authenticate the published source or diagnostic calibration.",
            "The explicit coordinate choice does not recover a physical trigger or identify breakdown from the gate geometry.",
            "No circuit fit, derivative, resampled waveform, FLASH boundary, or solver execution is supplied.",
        ],
    }
    if integrate is not None:
        start, stop = (_number(v, "integration bound") for v in integrate)
        report["interpolation_integrals"] = {
            "selected_interval_s": [start, stop], "source_interval_s": [start - offset_s, stop - offset_s],
            "integral_I_A_s": trace.integrate_current_a_s(start, stop),
            "integral_I_squared_A2_s": trace.integrate_current_squared_a2_s(start, stop),
            "interpretation": "piecewise_linear_interpolation_diagnostics", "is_heat": False,
        }
    plot_data = {"legend_flags": np.asarray(flags, dtype=bool), "lower_current_A": np.asarray(low_a),
                 "upper_current_A": np.asarray(high_a)}
    return report, trace, plot_data


def make_plot(report: dict, trace: DriveTrace, data: dict, destination: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(9, 5), layout="constrained")
    time_ns = trace.time_s * 1e9
    axis.fill_between(time_ns, data["lower_current_A"], data["upper_current_A"], color="firebrick", alpha=.2,
                      label="Published display envelope; not physical uncertainty")
    axis.plot(time_ns, trace.current_a, ".-", color="firebrick", lw=.8, ms=3,
              label="Encoded current vertices and linear interpolation")
    hidden = data["legend_flags"]
    axis.scatter(time_ns[hidden], trace.current_a[hidden], marker="s", s=18, facecolors="none", edgecolors="black",
                 label="Encoded behind source legend")
    gate = np.asarray(report["VUV_gate"]["selected_time_bounds_s"]) * 1e9
    axis.axvspan(*gate, alpha=.25, color="blueviolet", label="VUV gate geometry; not breakdown")
    axis.set(xlabel="Selected time coordinate [ns]", ylabel="Reconstructed current [A]",
             title=f"{trace.shot_id}; source figure {report['input']['source']['figure']}")
    axis.grid(alpha=.2)
    axis.legend(fontsize=8)
    figure.suptitle(f"Chosen time offset {trace.diagnostic_time_offset_s:g} s; physical trigger alignment unknown", fontsize=10)
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=150)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", type=Path, default=ROOT / "data/dissertation/shot_15210_current.csv")
    parser.add_argument("--provenance", type=Path, default=ROOT / "data/dissertation/shot_15210_current.json")
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--preserve-source-axis", action="store_true")
    choice.add_argument("--time-offset-s", type=float)
    parser.add_argument("--integrate", type=float, nargs=2, metavar=("START_S", "STOP_S"))
    parser.add_argument("--json-out", type=Path, help="explicit summary destination; stdout JSON is always emitted")
    parser.add_argument("--plot", type=Path, help="explicit plot destination; no resampled trace is written")
    args = parser.parse_args()
    try:
        inputs = {args.csv.resolve(), args.provenance.resolve()}
        outputs = [p.resolve() for p in (args.json_out, args.plot) if p is not None]
        if any(path in inputs for path in outputs):
            raise ValueError("Output paths must not overwrite either input file")
        if len(set(outputs)) != len(outputs):
            raise ValueError("JSON and plot output paths must be distinct")
        report, trace, plot_data = inspect_current(
            args.csv, args.provenance,
            mode="preserve_source_axis" if args.preserve_source_axis else "explicit_offset",
            offset_s=0.0 if args.preserve_source_axis else args.time_offset_s,
            integrate=None if args.integrate is None else tuple(args.integrate),
        )
        rendered = json.dumps(report, indent=2, allow_nan=False) + "\n"
        if args.plot is not None:
            make_plot(report, trace, plot_data, args.plot)
        if args.json_out is not None:
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            args.json_out.write_text(rendered, encoding="utf-8")
    except (ValueError, OSError, csv.Error, ImportError) as exc:
        parser.error(str(exc))
    print(rendered, end="")


if __name__ == "__main__":
    main()
