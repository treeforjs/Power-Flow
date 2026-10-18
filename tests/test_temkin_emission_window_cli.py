"""Independent source association, saved-budget and local-time export checks."""

import csv
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.source_convergence import manufactured_planar_surface
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_planar_io import write_temkin_planar_release_history
from mhdlab.temkin_release import TemkinReservoir


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools/inspect_temkin_emission_window.py"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def configuration(directory):
    config = json.loads((ROOT / "configs/temkin_window_manufactured.json").read_bytes())
    config["probabilities"] = [.1, .5, .9]
    path = directory / "parameters.json"
    path.write_text(json.dumps(config) + "\n")
    return path, config


def execute(output, *arguments):
    return subprocess.run([sys.executable, str(TOOL), *map(str, arguments), "--output-dir", str(output)],
                          capture_output=True, text=True)


@pytest.fixture(scope="module")
def demo_export(tmp_path_factory):
    directory = tmp_path_factory.mktemp("temkin-window-cli")
    path, config = configuration(directory)
    output = directory / "demo"
    process = execute(output, "--demo", "--config", path, "--plot", "--require-represented-budget")
    assert process.returncode == 0, process.stderr
    return path, config, output, json.loads((output / "report.json").read_bytes())


def test_manufactured_export_hashes_labels_qualified_budgets_and_local_time_columns(demo_export):
    path, config, output, report = demo_export
    assert report["schema"] == "mhdlab-temkin-window-inspection-v1"
    assert report["parameter_status"] == "manufactured"
    assert report["parameter_provenance"] == config["parameter_provenance"]
    assert report["input_config_sha256"] == sha(path) == sha(output / "input_config.json")
    assert report["input_source_sha256"] is None
    assert report["input_summary_sha256"] is None
    assert report["window"]["reservoir"] == config["component"]["reservoir"]
    assert report["window"]["status"] == "manufactured"
    assert report["window"]["all_budget_values_qualified"] is True
    assert report["window"]["host_budget_accepted"] is False
    budget = report["represented_budget"]
    assert budget["host_budget_accepted"] is False
    assert budget["values"] == {key: value["value"] for key, value in report["window"]["quantities"].items()}
    assert budget["quantity_units"] == report["window"]["quantity_units"]
    assert budget["quantity_units"]["released_number"] == "parent"
    assert budget["quantity_units"]["translational_energy_j"] == "J"
    for suffix in ("j_m2", "j"):
        values = [budget["values"][name + "_" + suffix] for name in
                  ("total_wall_energy_debit", "translational_energy", "internal_energy", "binding_energy")]
        residual = math.fsum((values[0], -values[1], -values[2], -values[3]))
        assert report["window"]["energy_balance_residual_" + suffix] == residual
        assert budget["energy_balance_residual_" + suffix] == residual
    assert all(report["runtime"][key] for key in ("python", "numpy", "scipy", "h5py"))
    assert (output / "source_window.png").is_file()
    assert set(report["output_sha256"]) == {entry.name for entry in output.iterdir() if entry.name != "report.json"}
    for name, digest in report["output_sha256"].items():
        assert sha(output / name) == digest
    for name, digest in report["source_sha256"].items():
        assert sha(ROOT / name) == digest
    with (output / "quantiles.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert set(rows[0]) == {"probability", "log_local_delay_s", "local_delay_s", "cdf_at_local_delay", "representation"}
    assert [float(row["probability"]) for row in rows] == config["probabilities"]
    previous = 0.
    for row in rows:
        delay = float(row["local_delay_s"])
        assert row["representation"] == "represented_local_delay"
        assert previous < delay < config["window"]["stop_s"] - config["window"]["start_s"]
        assert float(row["log_local_delay_s"]) == pytest.approx(math.log(delay), abs=1e-12)
        assert float(row["cdf_at_local_delay"]) == pytest.approx(float(row["probability"]), rel=1e-7)
        previous = delay
    assert not list(output.parent.glob(".temkin-window-*"))


def test_closed_hazard_quantiles_and_partial_number_from_saved_cli_inputs(tmp_path):
    path, config = configuration(tmp_path)
    config["component"]["reservoir"].update(attempt_frequency_s=2.5, activation_energy_ev=0.)
    config["component"].update(parent_duration_s=2., temperature_k=700., patch_area_m2=2e-6)
    config["window"].update(start_s=.3, stop_s=.9)
    path.write_text(json.dumps(config))
    output = tmp_path / "constant"
    process = execute(output, "--demo", "--config", path)
    assert process.returncode == 0, process.stderr
    report = json.loads((output / "report.json").read_bytes())
    left, right, rate = .3, .9, 2.5
    fraction = -math.expm1(-rate * (right - left))
    number = config["component"]["initial_coverage_m2"] * 2e-6 * math.exp(-rate * left) * fraction
    assert report["window"]["quantities"]["released_number"]["value"] == pytest.approx(number, rel=2e-11)
    with (output / "quantiles.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        probability = float(row["probability"])
        expected = -math.log1p(-probability * fraction) / rate
        assert float(row["local_delay_s"]) == pytest.approx(expected, rel=2e-11)


def test_unrepresented_energy_is_reported_and_required_budget_prevents_output(tmp_path):
    path, config = configuration(tmp_path)
    config["component"]["reservoir"].update(attempt_frequency_s=1., activation_energy_ev=0.,
                                               desorption_energy_ev=0., internal_energy_ev=0.)
    config["component"].update(initial_coverage_m2=1e-300, temperature_k=1e-100,
                                parent_duration_s=1., patch_area_m2=1.)
    config["window"].update(start_s=0., stop_s=1.)
    path.write_text(json.dumps(config))
    output = tmp_path / "diagnostic"
    process = execute(output, "--demo", "--config", path)
    assert process.returncode == 0, process.stderr
    report = json.loads((output / "report.json").read_bytes())
    assert report["represented_budget"] is None
    assert report["window"]["quantities"]["released_number"]["qualified"] is True
    energy = report["window"]["quantities"]["translational_energy_j"]
    assert energy["qualified"] is False
    assert energy["exact_zero"] is False
    assert energy["value"] == 0.
    assert math.isfinite(energy["log_value"])
    required = tmp_path / "required"
    process = execute(required, "--demo", "--config", path, "--require-represented-budget")
    assert process.returncode != 0
    assert "represented number/energy budget" in process.stderr
    assert not required.exists()


def test_lost_linear_delays_keep_log_seconds_and_explicit_flags(tmp_path):
    path, config = configuration(tmp_path)
    duration = float(np.nextafter(0., 1.)) * 4
    config["component"]["reservoir"].update(attempt_frequency_s=1., activation_energy_ev=0.)
    config["component"].update(parent_duration_s=duration, patch_area_m2=1.)
    config["window"].update(start_s=0., stop_s=duration)
    path.write_text(json.dumps(config))
    output = tmp_path / "log-times"
    process = execute(output, "--demo", "--config", path)
    assert process.returncode == 0, process.stderr
    with (output / "quantiles.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert any(row["representation"] == "linear_local_delay_unqualified_use_log" for row in rows)
    for row in rows:
        assert math.isfinite(float(row["log_local_delay_s"]))
        if row["representation"] == "linear_local_delay_unqualified_use_log":
            assert row["local_delay_s"] == ""
            assert row["cdf_at_local_delay"] == ""


@pytest.fixture(scope="module")
def canonical_source(tmp_path_factory):
    directory = tmp_path_factory.mktemp("window-source")
    base = manufactured_planar_surface(width_bins=2, current_bins=1, snapshot_intervals=2,
                                      duration_s=2., width_m=.003, current_length_m=.004)
    surface = PlanarSurfaceHistory(np.array([10., 11., 12.]), np.array([-.001, 0., .002]),
                                   np.array([0., .004]), np.tile([[[700., 900.]]], (3, 1, 1)), base.attributes)
    pool = TemkinReservoir(name="manufactured_component", state_id="M:manufactured-component",
                          attempt_frequency_s=.5, activation_energy_ev=0., desorption_energy_ev=.6,
                          internal_energy_ev=.2, coverage_coefficient=.35, reference_coverage_m2=1e19,
                          max_coverage_m2=5e18, kinetic_provenance="Manufactured constant-hazard source",
                          energy_provenance="Manufactured independent energy references")
    history = evolve_temkin_planar_release(surface, np.array([[[1e18, 2e18]]]), [pool],
                                          parameter_provenance="Manufactured canonical source for independent window CLI checks")
    path = directory / "source.h5"
    write_temkin_planar_release_history(path, history)
    return path, history


def external_config(directory, *, status="unvalidated"):
    path, config = configuration(directory)
    config["parameter_status"] = status
    config["component"] = {"interval_index": 1, "reservoir_index": 0, "current_index": 0, "width_index": 1}
    config["window"].update(start_s=.2, stop_s=.8)
    path.write_text(json.dumps(config))
    return path, config


def test_external_component_uses_depleted_density_actual_area_and_separate_origin(tmp_path, canonical_source):
    source, history = canonical_source
    path, config = external_config(tmp_path)
    output = tmp_path / "selected"
    process = execute(output, "--source", source, "--config", path)
    assert process.returncode == 0, process.stderr
    report = json.loads((output / "report.json").read_bytes())
    assert report["parameter_status"] == "unvalidated"
    assert report["input_source_sha256"] == sha(source) == sha(output / "source_input.h5")
    assert report["input_summary_sha256"] is None
    context, physical = report["source_component"], report["window"]
    assert context["indices"] == config["component"]
    assert context["original_interval_start_s"] == 11.
    assert context["original_interval_stop_s"] == 12.
    assert context["x_bounds_m"] == [0., .002]
    assert context["z_bounds_m"] == [0., .004]
    assert context["surface_attributes"]["source_kind"] == "synthetic"
    assert context["source_parameter_provenance"] == history.parameter_provenance
    assert physical["initial_coverage_m2"] == history.coverage_m2[1, 0, 0, 1]
    assert physical["initial_coverage_m2"] == pytest.approx(2e18 * math.exp(-.5), rel=2e-13)
    assert physical["temperature_k"] == 900.
    assert physical["parent_duration_s"] == 1.
    assert physical["patch_area_m2"] == history.patch_areas_m2[0, 1]
    assert physical["reservoir"] == asdict(history.reservoirs[0])
    assert physical["start_s"] == .2 and physical["stop_s"] == .8
    assert context["whole_interval_released_m2"] == history.released_m2[1, 0, 0, 1]
    count = physical["initial_coverage_m2"] * physical["patch_area_m2"] * math.exp(-.5 * .2) * -math.expm1(-.5 * .6)
    assert physical["quantities"]["released_number"]["value"] == pytest.approx(count, rel=2e-11)
    assert physical["quantities"]["released_m2"]["value"] < context["whole_interval_released_m2"]


@pytest.mark.parametrize("fault", [None, "hash", "provenance", "status", "schema"])
def test_external_parameter_label_needs_an_exact_source_summary_association(tmp_path, canonical_source, fault):
    source, history = canonical_source
    path, _ = external_config(tmp_path, status="manufactured")
    summary = {"schema": "mhdlab-temkin-planar-run-v1", "parameter_status": "manufactured",
               "parameter_provenance": history.parameter_provenance,
               "output_sha256": {"renamed_source.h5": sha(source)}}
    if fault == "hash":
        summary["output_sha256"]["renamed_source.h5"] = "0" * 64
    elif fault == "provenance":
        summary["parameter_provenance"] = "Different source parameters"
    elif fault == "status":
        summary["parameter_status"] = "sourced"
    elif fault == "schema":
        summary["schema"] = "unsupported-run"
    summary_path = tmp_path / "summary.json"
    summary_path.write_text(json.dumps(summary))
    output = tmp_path / "selected"
    process = execute(output, "--source", source, "--config", path, "--source-run-summary", summary_path)
    if fault is None:
        assert process.returncode == 0, process.stderr
        report = json.loads((output / "report.json").read_bytes())
        assert report["parameter_status"] == "manufactured"
        assert report["input_summary_sha256"] == sha(summary_path) == sha(output / "source_run_summary.json")
    else:
        assert process.returncode != 0
        assert not output.exists()


@pytest.mark.parametrize("fault", ["missing_summary", "boolean_index", "out_of_range", "wrong_component_keys"])
def test_external_selection_guards_leave_no_partial_output(tmp_path, canonical_source, fault):
    source, _ = canonical_source
    path, config = external_config(tmp_path)
    if fault == "missing_summary":
        config["parameter_status"] = "manufactured"
    elif fault == "boolean_index":
        config["component"]["interval_index"] = True
    elif fault == "out_of_range":
        config["component"]["width_index"] = 99
    else:
        config["component"]["patch_area_m2"] = 1.
    path.write_text(json.dumps(config))
    output = tmp_path / "selected"
    process = execute(output, "--source", source, "--config", path)
    assert process.returncode != 0
    assert not output.exists()


@pytest.mark.parametrize("existing", ["directory", "file", "symlink"])
def test_existing_output_is_preserved(tmp_path, existing):
    output = tmp_path / "window"
    sentinel = b"Preserve completed source-window work"
    if existing == "directory":
        output.mkdir()
        (output / "sentinel").write_bytes(sentinel)
    elif existing == "file":
        output.write_bytes(sentinel)
    else:
        output.symlink_to(tmp_path / "missing-target", target_is_directory=True)
    process = execute(output, "--demo")
    assert process.returncode != 0
    assert "new path" in process.stderr
    if existing == "directory":
        assert list(output.iterdir()) == [output / "sentinel"]
        assert (output / "sentinel").read_bytes() == sentinel
    elif existing == "file":
        assert output.read_bytes() == sentinel
    else:
        assert output.is_symlink()
        assert not (tmp_path / "missing-target").exists()


@pytest.mark.parametrize("kind", ["schema", "duplicate", "nonfinite", "unknown_top", "unknown_window",
                                  "boolean_area", "zero_width", "outside_parent", "closed_probability",
                                  "boolean_probability", "unsorted_probabilities", "parameter_label"])
def test_invalid_demo_inputs_leave_no_result(tmp_path, kind):
    path, config = configuration(tmp_path)
    if kind == "schema":
        config["schema"] = "unsupported-window"
    elif kind == "nonfinite":
        config["window"]["stop_s"] = float("nan")
    elif kind == "unknown_top":
        config["timestamp_s"] = 0.
    elif kind == "unknown_window":
        config["window"]["photon_gate"] = True
    elif kind == "boolean_area":
        config["component"]["patch_area_m2"] = True
    elif kind == "zero_width":
        config["window"]["stop_s"] = config["window"]["start_s"]
    elif kind == "outside_parent":
        config["window"]["stop_s"] = 2 * config["component"]["parent_duration_s"]
    elif kind == "closed_probability":
        config["probabilities"] = [0., .5]
    elif kind == "boolean_probability":
        config["probabilities"] = [False, .5]
    elif kind == "unsorted_probabilities":
        config["probabilities"] = [.9, .1]
    elif kind == "parameter_label":
        config["parameter_status"] = "sourced"
    text = json.dumps(config)
    if kind == "duplicate":
        text = text.replace('"schema":', f'"schema": {json.dumps(config["schema"])}, "schema":', 1)
    path.write_text(text)
    output = tmp_path / "window"
    process = execute(output, "--demo", "--config", path)
    assert process.returncode != 0
    assert not output.exists()
    assert not list(tmp_path.glob(".temkin-window-*"))


@pytest.mark.parametrize("arguments", [[], ["--demo", "--source", "missing.h5"],
                                         ["--source", "missing.h5"], ["--demo", "--source-run-summary", "missing.json"]])
def test_source_selection_cannot_be_ambiguous_or_implicitly_labeled(tmp_path, arguments):
    output = tmp_path / "window"
    process = execute(output, *arguments)
    assert process.returncode != 0
    assert not output.exists()
