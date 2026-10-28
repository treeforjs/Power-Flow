"""Fresh-output provenance and independently checked manufactured readout scan."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools/compare_neutral_ensemble.py"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def config_file(directory):
    config = json.loads((ROOT / "configs/neutral_ensemble_manufactured.json").read_bytes())
    path = directory / "request.json"
    path.write_text(json.dumps(config))
    return path, config


def execute(path, output, *extra):
    return subprocess.run([sys.executable, str(TOOL), "--config", str(path),
                           "--output-dir", str(output), *extra], capture_output=True, text=True)


@pytest.fixture(scope="module")
def manufactured_scan(tmp_path_factory):
    directory = tmp_path_factory.mktemp("neutral-ensemble-cli")
    path, config = config_file(directory)
    output = directory / "comparison"
    process = execute(path, output, "--plot")
    assert process.returncode == 0, process.stderr
    report = json.loads((output / "comparison_report.json").read_bytes())
    return path, config, output, report, json.loads(process.stdout)


def test_fresh_output_records_current_source_input_and_saved_hashes(manufactured_scan):
    path, config, output, report, printed = manufactured_scan
    assert report["schema"] == "mhdlab-neutral-ensemble-run-v1"
    assert report["input_config_sha256"] == sha(path)
    assert json.loads((output / "input_config.json").read_bytes()) == config
    assert printed["diagnostics"] == report["comparison"]["diagnostics"]
    assert (output / "neutral_ensemble_comparison.png").is_file()
    from PIL import Image
    with Image.open(output / "neutral_ensemble_comparison.png") as figure:
        assert figure.width > 1000 and figure.height > 500
    for name, digest in report["source_sha256"].items():
        assert sha(ROOT / name) == digest
    for name, digest in report["output_sha256"].items():
        assert sha(output / name) == digest
    assert not list(output.parent.glob(".neutral-ensemble-*"))


def test_every_saved_score_matches_independent_asymmetric_envelope_calculation(manufactured_scan):
    _, config, output, report, _ = manufactured_scan
    curve = config["observation"]["curves"][0]
    source_cases = {case["case_id"]: case for case in config["cases"]}
    with (output / "case_scale_comparisons.csv").open() as stream:
        csv_rows = list(csv.DictReader(stream))
    saved_rows = report["comparison"]["case_scale_comparisons"]
    assert len(csv_rows) == len(saved_rows) == 41*len(source_cases)
    for saved, csv_row in zip(saved_rows, csv_rows):
        scale = saved["C_abs"]
        values = np.asarray(source_cases[saved["case_id"]]["curves"][0]["values"])
        center = scale*np.asarray(curve["center"])
        low, high = scale*np.asarray(curve["lower"]), scale*np.asarray(curve["upper"])
        delta = values-center
        independent_score = np.max(np.abs(delta)/np.where(delta >= 0, high-center, center-low))
        assert saved["max_normalized_departure"] == pytest.approx(independent_score, rel=2e-14)
        assert saved["readout_compatible"] == (independent_score <= 1)
        assert float(csv_row["max_normalized_departure"]) == saved["max_normalized_departure"]
    truth = config["manufactured_truth"]
    target = next(row for row in saved_rows if row["case_id"] == truth["case_id"] and row["C_abs"] == truth["C_abs"])
    assert target["max_normalized_departure"] == 0


def test_recorded_truth_is_never_used_by_comparison(tmp_path, manufactured_scan):
    _, _, _, original, _ = manufactured_scan
    path, config = config_file(tmp_path)
    config["manufactured_truth"] = {"case_id": "does_not_exist", "C_abs": 99999}
    path.write_text(json.dumps(config))
    output = tmp_path / "comparison"
    process = execute(path, output)
    assert process.returncode == 0, process.stderr
    report = json.loads((output / "comparison_report.json").read_bytes())
    assert report["comparison"] == original["comparison"]
    assert report["manufactured_truth_not_used_by_comparison"] == config["manufactured_truth"]


def test_existing_result_is_never_overwritten(manufactured_scan):
    path, _, output, _, _ = manufactured_scan
    before = sha(output / "comparison_report.json")
    process = execute(path, output)
    assert process.returncode == 2
    assert "must be a new path" in process.stderr
    assert sha(output / "comparison_report.json") == before


def test_unavailable_physical_mapping_rejects_without_output(tmp_path):
    path, config = config_file(tmp_path)
    config["observation"]["mapping_status"] = "unavailable"
    path.write_text(json.dumps(config))
    output = tmp_path / "comparison"
    process = execute(path, output)
    assert process.returncode == 2
    assert "withhold inference" in process.stderr
    assert not output.exists()


def test_duplicate_json_keys_and_nonfinite_json_reject(tmp_path):
    path = tmp_path / "request.json"
    output = tmp_path / "comparison"
    path.write_text('{"schema":"a","schema":"b"}')
    process = execute(path, output)
    assert process.returncode == 2 and "duplicate JSON key" in process.stderr
    path.write_text('{"invalid": NaN}')
    process = execute(path, output)
    assert process.returncode == 2 and "nonfinite JSON constant" in process.stderr
    assert not output.exists()
