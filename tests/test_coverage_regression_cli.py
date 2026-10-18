"""Independent fresh-output, source-label and saved-objective regression checks."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools/fit_surface_coverage.py"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def config_file(directory):
    config = json.loads((ROOT / "configs/coverage_regression_manufactured.json").read_bytes())
    path = directory / "parameters.json"
    path.write_text(json.dumps(config))
    return path, config


def execute(path, output, *extra):
    return subprocess.run([sys.executable, str(TOOL), "--config", str(path),
                           "--output-dir", str(output), *extra], capture_output=True, text=True)


@pytest.fixture(scope="module")
def manufactured_fit(tmp_path_factory):
    directory = tmp_path_factory.mktemp("coverage-cli")
    path, config = config_file(directory)
    output = directory / "fit"
    process = execute(path, output, "--plot")
    assert process.returncode == 0, process.stderr
    report = json.loads((output / "fit_report.json").read_bytes())
    return path, config, output, report, json.loads(process.stdout)


def test_manufactured_off_grid_fit_records_definition_status_hashes_and_plot(manufactured_fit):
    path, config, output, report, printed = manufactured_fit
    assert report["schema"] == "mhdlab-surface-coverage-fit-run-v1"
    assert report["input_config_sha256"] == sha(path)
    assert json.loads((output / "input_config.json").read_bytes()) == config
    assert (output / "coverage_fit.png").is_file()
    assert set(report["output_sha256"]) == {item.name for item in output.iterdir()
                                           if item.name != "fit_report.json"}
    for name, digest in report["output_sha256"].items():
        assert sha(output / name) == digest
    for name, digest in report["source_sha256"].items():
        assert sha(ROOT / name) == digest
    assert printed["best_fit"] == report["fit"]["best_fit"]
    assert printed["diagnostics"] == report["fit"]["diagnostics"]
    assert report["input_labels"]["prediction_status"] == "manufactured"
    assert report["input_labels"]["observation_status"] == "manufactured"
    assert report["input_labels"]["monolayer_provenance"] == config["grid"]["monolayer_provenance"]
    fit, truth = report["fit"], config["manufactured_truth"]
    assert fit["monolayer"]["areal_density_m2"] == config["grid"]["monolayer_areal_density_m2"]
    assert fit["monolayer"]["definition"] == config["grid"]["monolayer_definition"]
    assert fit["monolayer"]["provenance"] == config["grid"]["monolayer_provenance"]
    best = fit["best_fit"]
    assert best["coverage_ml"] == pytest.approx(truth["coverage_ml"], rel=2e-11)
    assert best["coverage_ml"] not in config["grid"]["coverage_ml"]
    assert best["calibration_scale"] == pytest.approx(truth["calibration_scale"], rel=2e-12)
    assert best["areal_inventory_m2"] == pytest.approx(truth["areal_inventory_m2"], rel=2e-11)
    assert best["coverage_boundary"] == "interior"
    assert best["chi2"] < 1e-22
    assert fit["diagnostics"]["extrapolation_performed"] is False
    assert fit["bin_support"]["signal_units"] == config["observation"]["signal_units"]
    assert fit["bin_support"]["bin_axis_units"] == config["observation"]["bin_axis_units"]
    assert report["manufactured_truth_not_used_by_fit"] == truth
    assert not list(output.parent.glob(".coverage-fit-*"))


def test_best_fit_csv_means_scaled_signal_with_unscaled_known_background(manufactured_fit):
    _, config, output, report, _ = manufactured_fit
    best = report["fit"]["best_fit"]
    knots = np.asarray(config["grid"]["coverage_ml"])
    templates = np.asarray(config["grid"]["predictions"])
    raw = np.array([np.interp(best["coverage_ml"], knots, column) for column in templates.T])
    scaled = best["calibration_scale"] * raw
    background = np.asarray(config["observation"]["known_background"])
    np.testing.assert_allclose(best["prediction_signal"], scaled, rtol=3e-14)
    np.testing.assert_allclose(best["prediction_with_background"], scaled + background, rtol=3e-14)
    with (output / "best_fit_bins.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == len(background)
    for index, row in enumerate(rows):
        assert int(row["bin_index"]) == index
        assert float(row["bin_left"]) == config["grid"]["bin_edges"][index]
        assert float(row["bin_right"]) == config["grid"]["bin_edges"][index + 1]
        assert float(row["observed_signal"]) == config["observation"]["values"][index]
        assert float(row["known_background"]) == background[index]
        assert float(row["scaled_prediction_signal"]) == best["prediction_signal"][index]
        assert float(row["prediction_with_background"]) == best["prediction_with_background"][index]
        assert float(row["observed_minus_prediction"]) == (config["observation"]["values"][index]
                                                         - best["prediction_with_background"][index])


def test_saved_profile_and_objective_retain_full_covariance(manufactured_fit):
    _, config, output, report, _ = manufactured_fit
    profile = report["fit"]["grid_profile"]
    covariance = np.asarray(config["observation"]["covariance"])
    target = np.asarray(config["observation"]["values"]) - config["observation"]["known_background"]
    prior_mean, prior_sigma = (config["calibration"][name] for name in ("prior_mean", "prior_sigma"))
    with (output / "coverage_profile.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == len(config["grid"]["coverage_ml"])
    for index, row in enumerate(rows):
        model = np.asarray(config["grid"]["predictions"][index])
        scale = profile["calibration_scale"][index]
        residual = target - scale * model
        data_cost = residual @ np.linalg.solve(covariance, residual)
        prior_cost = ((scale - prior_mean) / prior_sigma)**2
        assert profile["data_chi2"][index] == pytest.approx(data_cost, rel=8e-14, abs=1e-26)
        assert profile["prior_chi2"][index] == pytest.approx(prior_cost, rel=8e-14, abs=1e-26)
        assert profile["chi2"][index] == pytest.approx(data_cost + prior_cost, rel=8e-14, abs=1e-26)
        for key, values in profile.items():
            if isinstance(values[index], bool):
                assert row[key] == str(values[index])
            else:
                assert float(row[key]) == values[index]


@pytest.mark.parametrize("truth_case", ["removed", "altered"])
def test_fit_is_independent_of_recorded_manufactured_truth(tmp_path, manufactured_fit, truth_case):
    _, _, _, original_report, _ = manufactured_fit
    path, config = config_file(tmp_path)
    if truth_case == "removed":
        del config["manufactured_truth"]
    else:
        config["manufactured_truth"] = {"coverage_ml": 1e8, "calibration_scale": 0.,
                                        "areal_inventory_m2": 0., "usage": "Deliberately incorrect validation label"}
    path.write_text(json.dumps(config))
    output = tmp_path / "fit"
    process = execute(path, output)
    assert process.returncode == 0, process.stderr
    report = json.loads((output / "fit_report.json").read_bytes())
    assert report["fit"] == original_report["fit"]
    if truth_case == "removed":
        assert "manufactured_truth_not_used_by_fit" not in report
    else:
        assert report["manufactured_truth_not_used_by_fit"] == config["manufactured_truth"]


def test_full_covariance_reflection_changes_fit_as_independent_gls_predicts(tmp_path):
    path, config = config_file(tmp_path)
    del config["manufactured_truth"]
    config["grid"]["coverage_ml"] = [0., .5, 1.]
    config["grid"]["predictions"] = [[0., 0.], [.5, 1.], [1., 2.]]
    for record in (config["grid"], config["observation"]):
        record["bin_edges"] = [110., 120., 130.]
        record["bin_support"] = "Two identical supplied manufactured detector bins"
    config["observation"]["values"] = [1., 1.]
    config["observation"]["known_background"] = [0., 0.]
    config["calibration"] = {"mode": "fixed", "fixed_scale": 1., "provenance": "Explicit manufactured unit scale"}
    coverage = []
    for label, covariance in (("full", [[1., .8], [.8, 1.]]), ("diagonal", [[1., 0.], [0., 1.]])):
        config["observation"]["covariance"] = covariance
        path.write_text(json.dumps(config))
        output = tmp_path / label
        process = execute(path, output)
        assert process.returncode == 0, process.stderr
        report = json.loads((output / "fit_report.json").read_bytes())
        best = report["fit"]["best_fit"]
        direction, target = np.array([1., 2.]), np.array([1., 1.])
        inverse_direction = np.linalg.solve(covariance, direction)
        expected = target @ inverse_direction / (direction @ inverse_direction)
        assert best["coverage_ml"] == pytest.approx(expected, rel=8e-14)
        residual = target - best["coverage_ml"] * direction
        assert best["chi2"] == pytest.approx(residual @ np.linalg.solve(covariance, residual), rel=8e-14)
        np.testing.assert_allclose(best["whitened_residual"], np.linalg.solve(np.linalg.cholesky(covariance), residual), rtol=8e-14)
        assert report["fit"]["diagnostics"]["covariance_condition_number"] == pytest.approx(np.linalg.cond(covariance), rel=3e-14)
        coverage.append(best["coverage_ml"])
    np.testing.assert_allclose(coverage, [1 / 3, .6], rtol=8e-14)


def test_observed_labels_remain_explicit_after_manufactured_truth_is_removed(tmp_path):
    path, config = config_file(tmp_path)
    del config["manufactured_truth"]
    config["grid"]["status"] = "owner_supplied"
    config["observation"]["status"] = "reconstructed"
    path.write_text(json.dumps(config))
    output = tmp_path / "fit"
    process = execute(path, output)
    assert process.returncode == 0, process.stderr
    report = json.loads((output / "fit_report.json").read_bytes())
    assert report["input_labels"]["prediction_status"] == "owner_supplied"
    assert report["input_labels"]["observation_status"] == "reconstructed"
    assert report["fit"]["inputs"]["prediction_status"] == "owner_supplied"
    assert report["fit"]["inputs"]["observation_status"] == "reconstructed"
    assert "manufactured_truth_not_used_by_fit" not in report


@pytest.mark.parametrize("existing", ["directory", "file", "symlink"])
def test_existing_destination_is_preserved_without_replacement(tmp_path, existing):
    path, _ = config_file(tmp_path)
    output = tmp_path / "fit"
    sentinel = b"Preserve prior recovered evidence"
    if existing == "directory":
        output.mkdir()
        (output / "sentinel").write_bytes(sentinel)
    elif existing == "file":
        output.write_bytes(sentinel)
    else:
        output.symlink_to(tmp_path / "missing-target", target_is_directory=True)
    process = execute(path, output)
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
    assert not list(tmp_path.glob(".coverage-fit-*"))


@pytest.mark.parametrize("kind", ["schema", "duplicate", "nonfinite", "unknown_top", "unknown_nested",
                                  "support", "units", "edges", "axis", "asymmetric_covariance",
                                  "singular_covariance", "bad_prior", "missing_layer_definition",
                                  "boolean_inventory", "observed_truth", "owner_truth"])
def test_invalid_or_ambiguous_inputs_leave_no_result_or_staging_directory(tmp_path, kind):
    path, config = config_file(tmp_path)
    if kind == "schema":
        config["schema"] = "unsupported-fit"
    elif kind == "nonfinite":
        config["observation"]["values"][0] = float("nan")
    elif kind == "unknown_top":
        config["temperature_k"] = 600.
    elif kind == "unknown_nested":
        config["grid"]["activation_energy_ev"] = 1.
    elif kind == "support":
        config["observation"]["bin_support"] = "Different aperture or gate"
    elif kind == "units":
        config["observation"]["signal_units"] = "photons per second"
    elif kind == "edges":
        config["observation"]["bin_edges"][1] += 1.
    elif kind == "axis":
        config["observation"]["bin_axis_units"] = "m"
    elif kind == "asymmetric_covariance":
        config["observation"]["covariance"][0][1] = .4
    elif kind == "singular_covariance":
        config["observation"]["covariance"] = np.ones((4, 4)).tolist()
    elif kind == "bad_prior":
        config["calibration"]["prior_sigma"] = 0.
    elif kind == "missing_layer_definition":
        del config["grid"]["monolayer_definition"]
    elif kind == "boolean_inventory":
        config["grid"]["monolayer_areal_density_m2"] = True
    elif kind == "observed_truth":
        config["observation"]["status"] = "observed"
    elif kind == "owner_truth":
        config["grid"]["status"] = "owner_supplied"
    text = json.dumps(config)
    if kind == "duplicate":
        text = text.replace('"schema":', f'"schema": {json.dumps(config["schema"])}, "schema":', 1)
    path.write_text(text)
    output = tmp_path / "fit"
    process = execute(path, output)
    assert process.returncode != 0
    assert not output.exists()
    assert not list(tmp_path.glob(".coverage-fit-*"))
