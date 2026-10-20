"""CLI checks for provenance, explicit time alignment and published-curve diagnostics."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools/inspect_current_drive.py"
CSV_SOURCE = ROOT / "data/dissertation/shot_15210_current.csv"
PROVENANCE_SOURCE = ROOT / "data/dissertation/shot_15210_current.json"


def invoke(*arguments, cwd=None):
    return subprocess.run([sys.executable, str(TOOL), *map(str, arguments)],
                          cwd=cwd or ROOT, capture_output=True, text=True, timeout=30)


def report_for(*arguments, cwd=None):
    result = invoke(*arguments, cwd=cwd)
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def copied_inputs(tmp_path):
    csv_path = tmp_path / CSV_SOURCE.name
    provenance_path = tmp_path / PROVENANCE_SOURCE.name
    csv_path.write_bytes(CSV_SOURCE.read_bytes())
    provenance_path.write_bytes(PROVENANCE_SOURCE.read_bytes())
    return csv_path, provenance_path


def source_rows(path=CSV_SOURCE):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def rewrite_csv_and_hash(csv_path, provenance_path, rows):
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    provenance = json.loads(provenance_path.read_text())
    provenance["output"]["csv_sha256"] = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    provenance_path.write_text(json.dumps(provenance), encoding="utf-8")


def assert_rejected(*arguments):
    result = invoke(*arguments)
    assert result.returncode != 0
    assert result.stderr.strip(), "A rejected request should explain its failure"


def test_default_inputs_require_explicit_source_axis_and_are_read_only(tmp_path):
    before = {p: p.read_bytes() for p in (CSV_SOURCE, PROVENANCE_SOURCE)}
    assert list(tmp_path.iterdir()) == []
    report = report_for("--preserve-source-axis", cwd=tmp_path)
    assert {p: p.read_bytes() for p in before} == before
    assert list(tmp_path.iterdir()) == [], "Default inspection should only write JSON to stdout"
    assert report["source_kind"] == "reconstructed"
    assert report["input"]["csv_sha256"] == hashlib.sha256(before[CSV_SOURCE]).hexdigest()
    assert report["input"]["provenance_sha256"] == hashlib.sha256(before[PROVENANCE_SOURCE]).hexdigest()
    alignment = report["alignment"]
    assert alignment["mode"] == "preserve_source_axis"
    assert alignment["applied_time_offset_s"] == 0
    assert alignment["source_time_origin_reference"] is None
    assert alignment["physical_trigger_alignment_established"] is False
    waveform = report["waveform"]
    assert waveform["samples"] == 118
    assert waveform["legend_covered_vertices"] == 20
    assert waveform["source_time_support_s"] == pytest.approx([2.249774141568e-9, 198.2496266277e-9])
    assert waveform["selected_time_support_s"] == waveform["source_time_support_s"]
    assert waveform["minimum_current_A"] == pytest.approx(-258.644689053)
    assert waveform["negative_values"]["count"] > 0
    assert waveform["negative_values"]["minimum_current_A"] == waveform["minimum_current_A"]
    peak = waveform["peak"]
    assert peak["current_A"] == pytest.approx(847978.195498)
    assert peak["source_time_s"] == pytest.approx(136.749994024e-9)
    assert peak["selected_time_s"] == peak["source_time_s"]
    assert peak["source_vertex_index"] == 82
    assert peak["occluded_by_legend"] in (False, 0)
    bounds = report["display_bounds"]
    assert bounds["physical_uncertainty"] is False
    assert bounds["time_plus_minus_s"] == pytest.approx(.327196320613e-9)
    assert bounds["current_plus_minus_A"] == pytest.approx(2000.10274860)
    gate = report["VUV_gate"]
    assert gate["kind"] == "published_gate_rectangle_not_breakdown"
    assert gate["source_time_bounds_s"] == pytest.approx([52.0545737354e-9, 55.0545612339e-9])
    assert gate["selected_time_bounds_s"] == gate["source_time_bounds_s"]
    assert gate["actual_breakdown_time_s"] is None
    assert gate["exact_trigger_reference"] is None
    assert "interpolation_integrals" not in report


def test_explicit_offset_moves_current_and_gate_but_does_not_establish_trigger():
    offset = -50.e-9
    rows = source_rows()
    # A single original segment is a linear function. Its signed area and
    # squared-current area have independent analytic integrals; no physics
    # model or numerical quadrature is needed for this oracle.
    a, b = rows[0], rows[1]
    t0, t1 = float(a["time_s"]), float(b["time_s"])
    i0, i1 = float(a["current_A"]), float(b["current_A"])
    report = report_for("--time-offset-s=-5e-8", "--integrate",
                        format(t0 + offset, ".24f"), format(t1 + offset, ".24f"))
    alignment = report["alignment"]
    assert alignment["mode"] == "explicit_offset"
    assert alignment["applied_time_offset_s"] == offset
    assert alignment["physical_trigger_alignment_established"] is False
    waveform = report["waveform"]
    assert waveform["selected_time_support_s"] == pytest.approx(
        [v + offset for v in waveform["source_time_support_s"]], rel=1.e-12, abs=1.e-20)
    assert waveform["peak"]["selected_time_s"] == pytest.approx(waveform["peak"]["source_time_s"] + offset)
    gate = report["VUV_gate"]
    assert gate["selected_time_bounds_s"] == pytest.approx([v + offset for v in gate["source_time_bounds_s"]])
    assert gate["actual_breakdown_time_s"] is None
    assert gate["exact_trigger_reference"] is None
    integral = report["interpolation_integrals"]
    assert integral["interpretation"] == "piecewise_linear_interpolation_diagnostics"
    assert integral["is_heat"] is False
    assert integral["selected_interval_s"] == pytest.approx([t0 + offset, t1 + offset])
    assert integral["source_interval_s"] == pytest.approx([t0, t1])
    assert integral["integral_I_A_s"] == pytest.approx((t1 - t0) * (i0 + i1) / 2, rel=1.e-9)
    assert integral["integral_I_A_s"] < 0, "Signed early current must not be floored or absolutized"
    assert integral["integral_I_squared_A2_s"] == pytest.approx(
        (t1 - t0) * (i0 * i0 + i0 * i1 + i1 * i1) / 3, rel=1.e-9)


@pytest.mark.parametrize("arguments", [(), ("--preserve-source-axis", "--time-offset-s=0")])
def test_alignment_choice_is_required_and_mutually_exclusive(arguments):
    assert_rejected(*arguments)


@pytest.mark.parametrize("offset", ["nan", "inf", "-inf"])
def test_nonfinite_time_offset_rejected(offset):
    assert_rejected(f"--time-offset-s={offset}")


def test_csv_hash_mismatch_rejected(tmp_path):
    csv_path, provenance_path = copied_inputs(tmp_path)
    csv_path.write_bytes(csv_path.read_bytes() + b"\n")
    assert_rejected("--csv", csv_path, "--provenance", provenance_path, "--preserve-source-axis")


@pytest.mark.parametrize("alteration", ["current_units", "current_scale", "status"])
def test_changed_provenance_semantics_rejected(tmp_path, alteration):
    csv_path, provenance_path = copied_inputs(tmp_path)
    metadata = json.loads(provenance_path.read_text())
    if alteration == "current_units":
        metadata["observable"]["original_units"] = "kA"
    elif alteration == "current_scale":
        metadata["calibration"]["current_axis"]["A_per_MA"] = 1.e3
    else:
        metadata["status"] = "original_oscilloscope_waveform"
    provenance_path.write_text(json.dumps(metadata), encoding="utf-8")
    assert_rejected("--csv", csv_path, "--provenance", provenance_path, "--preserve-source-axis")


@pytest.mark.parametrize("field", ["current_A", "time_s"])
def test_si_column_inconsistency_rejected_even_with_updated_csv_hash(tmp_path, field):
    csv_path, provenance_path = copied_inputs(tmp_path)
    rows = source_rows(csv_path)
    rows[30][field] = str(float(rows[30][field]) * 1.0001)
    rewrite_csv_and_hash(csv_path, provenance_path, rows)
    assert_rejected("--csv", csv_path, "--provenance", provenance_path, "--preserve-source-axis")


@pytest.mark.parametrize("field", ["time_digitization_lower_ns", "current_digitization_lower_A"])
def test_display_envelope_inconsistency_rejected_after_csv_rehash(tmp_path, field):
    csv_path, provenance_path = copied_inputs(tmp_path)
    rows = source_rows(csv_path)
    rows[30][field] = str(float(rows[30][field]) * .95)
    rewrite_csv_and_hash(csv_path, provenance_path, rows)
    assert_rejected("--csv", csv_path, "--provenance", provenance_path, "--preserve-source-axis")


@pytest.mark.parametrize("alteration", ["nonbinary", "swapped"])
def test_invalid_legend_flags_rejected_even_with_updated_hash_and_unchanged_count(tmp_path, alteration):
    csv_path, provenance_path = copied_inputs(tmp_path)
    rows = source_rows(csv_path)
    visible = next(row for row in rows if row["occluded_by_legend"] == "0")
    covered = [row for row in rows if row["occluded_by_legend"] == "1"]
    if alteration == "swapped":
        # Keep twenty binary flags while mislabelling their geometric locations.
        visible["occluded_by_legend"] = "1"
        covered[0]["occluded_by_legend"] = "0"
    else:
        # Keep the sum at twenty; nonbinary values cannot stand in for flags.
        visible["occluded_by_legend"] = "2"
        covered[0]["occluded_by_legend"] = "0"
        covered[1]["occluded_by_legend"] = "0"
    assert sum(int(row["occluded_by_legend"]) for row in rows) == 20
    rewrite_csv_and_hash(csv_path, provenance_path, rows)
    assert_rejected("--csv", csv_path, "--provenance", provenance_path, "--preserve-source-axis")


@pytest.mark.parametrize("flag", ["--json-out", "--plot"])
@pytest.mark.parametrize("input_target", ["csv", "provenance"])
def test_output_cannot_overwrite_either_input(tmp_path, flag, input_target):
    csv_path, provenance_path = copied_inputs(tmp_path)
    before = {path: path.read_bytes() for path in (csv_path, provenance_path)}
    target = csv_path if input_target == "csv" else provenance_path
    assert_rejected("--csv", csv_path, "--provenance", provenance_path,
                    "--preserve-source-axis", flag, target)
    assert {path: path.read_bytes() for path in before} == before


@pytest.mark.parametrize("interval", [(60.e-9, 50.e-9), (0., 50.e-9), (50.e-9, 200.e-9)])
def test_reversed_or_outside_support_integrals_rejected(interval):
    assert_rejected("--preserve-source-axis", "--integrate", *map(repr, interval))


def test_explicit_report_and_plot_outputs(tmp_path):
    output = tmp_path / "inspection.json"
    plot = tmp_path / "inspection.png"
    report = report_for("--preserve-source-axis", "--json-out", output, "--plot", plot)
    assert json.loads(output.read_text()) == report
    assert plot.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert plot.stat().st_size > 1000
