"""Independent export checks for conditional upstream source observables."""

import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import h5py
import numpy as np
import pytest

from mhdlab.temkin_planar_io import load_temkin_planar_release_history


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools/inspect_joint_source_moments.py"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def execute(output, *arguments):
    return subprocess.run(
        [sys.executable, str(TOOL), *map(str, arguments), "--output-dir", str(output)],
        capture_output=True, text=True,
    )


def configuration(tmp_path, **changes):
    config = json.loads((ROOT / "configs/joint_source_manufactured.json").read_text())
    config["geometry"].update(width_bins=3, current_bins=1, snapshot_intervals=2)
    config["subdivisions_per_interval"] = 1
    config.update(changes)
    path = tmp_path / "joint_parameters.json"
    path.write_text(json.dumps(config) + "\n")
    return path, config


@pytest.fixture(scope="module")
def demo_export(tmp_path_factory):
    work = tmp_path_factory.mktemp("joint-cli")
    path, config = configuration(work)
    output = work / "demo"
    result = execute(output, "--demo", "--config", path, "--quadrature-order", "4", "--plot")
    assert result.returncode == 0, result.stderr
    return path, config, output, json.loads((output / "report.json").read_text())


def test_manufactured_export_hashes_and_self_contained_source(demo_export):
    path, config, output, report = demo_export
    assert report["schema"] == "mhdlab-joint-source-run-v1"
    assert report["parameter_status"] == "manufactured"
    assert report["parameter_provenance"] == config["parameter_provenance"]
    assert report["input_config_sha256"] == sha(path) == sha(output / "input_config.json")
    assert report["input_source_sha256"] is None
    assert report["input_summary_sha256"] is None
    assert (output / "joint_moments.png").is_file()
    assert set(report["output_sha256"]) == {item.name for item in output.iterdir() if item.name != "report.json"}
    for name, digest in report["output_sha256"].items():
        assert sha(output / name) == digest
    for name, digest in report["source_sha256"].items():
        assert sha(ROOT / name) == digest
    history = load_temkin_planar_release_history(output / "source_input.h5")
    assert history.temperature_k.shape == (2, 1, 3)
    assert history.parameter_provenance == config["parameter_provenance"]
    assert np.all(history.released_m2 > 0)
    assert report["joint_moments"]


def test_export_retains_thermal_mixture_and_joint_width_control(demo_export):
    _, _, output, report = demo_export
    history = load_temkin_planar_release_history(output / "source_input.h5")
    joint = report["joint_moments"]
    assert joint["schema"] == "mhdlab-joint-source-moments-v1"
    assert joint["units"]["cov_time_width_centered_square_s_m2"] == "s m^2"
    assert joint["units"]["mean_expected_translational_energy_j"] == "J"
    width = history.x_edges_m[-1] - history.x_edges_m[0]
    center = (history.x_edges_m[-1] + history.x_edges_m[0]) / 2
    boltzmann = 1.380649e-23
    for index, pool in enumerate(joint["reservoirs"]):
        weights = history.released_m2[:, index] * history.patch_areas_m2[None, :, :]
        expected_number = weights.sum()
        expected_energy = np.sum(weights * (2 * boltzmann * history.temperature_k)) / expected_number
        resolved, uniform = pool["resolved"], pool["uniform_width"]
        assert pool["state_id"] == history.reservoirs[index].state_id
        np.testing.assert_allclose(pool["released_number"], expected_number, rtol=2e-14)
        for stats in (resolved, uniform):
            np.testing.assert_allclose(stats["released_number"], expected_number, rtol=2e-14)
            np.testing.assert_allclose(stats["mean_expected_translational_energy_j"], expected_energy, rtol=2e-14)
        for name in ("mean_time_s", "variance_time_s2", "mean_z_m", "variance_z_m2",
                     "mean_expected_translational_energy_j", "variance_expected_translational_energy_j2",
                     "cov_time_z_s_m", "cov_time_expected_translational_energy_s_j"):
            assert resolved[name] == uniform[name]
        np.testing.assert_allclose(uniform["mean_x_m"], center, atol=1e-18)
        np.testing.assert_allclose(uniform["variance_x_m2"], width * width / 12, rtol=2e-14)
        for name in ("cov_time_x_s_m", "cov_time_width_centered_square_s_m2",
                     "cov_expected_translational_energy_width_centered_square_j_m2"):
            assert uniform[name] == 0
        assert resolved["cov_expected_translational_energy_width_centered_square_j_m2"] > 0
        assert resolved["cov_time_width_centered_square_s_m2"] != 0
        for group_name in ("conditional_width", "conditional_interval"):
            groups = pool[group_name]
            for hypothesis in ("resolved", "uniform_width"):
                counts = np.array([group[hypothesis]["released_number"] for group in groups])
                np.testing.assert_allclose(counts.sum(), expected_number, rtol=2e-14)
                for mean in ("mean_time_s", "mean_x_m", "mean_z_m", "mean_expected_translational_energy_j"):
                    mixed = math.fsum(number * group[hypothesis][mean] for number, group in zip(counts, groups)
                                      if number > 0) / expected_number
                    scale = {"mean_time_s": history.time_edges_s[-1] - history.time_edges_s[0],
                             "mean_x_m": width,
                             "mean_z_m": history.z_edges_m[-1] - history.z_edges_m[0],
                             "mean_expected_translational_energy_j": expected_energy}[mean]
                    np.testing.assert_allclose(mixed, pool[hypothesis][mean], rtol=5e-13,
                                               atol=16 * np.finfo(float).eps * scale)


def test_conditional_csv_reproduces_report_without_losing_units_or_nulls(demo_export):
    _, _, output, report = demo_export
    for file_name, group_name, lower, upper in (
        ("conditional_width.csv", "conditional_width", "x_left_m", "x_right_m"),
        ("conditional_interval.csv", "conditional_interval", "time_start_s", "time_stop_s"),
    ):
        with (output / file_name).open() as stream:
            rows = list(csv.DictReader(stream))
        expected_rows = sum(2 * len(pool[group_name]) for pool in report["joint_moments"]["reservoirs"])
        assert len(rows) == expected_rows
        for pool in report["joint_moments"]["reservoirs"]:
            for index, group in enumerate(pool[group_name]):
                for hypothesis in ("resolved", "uniform_width"):
                    selected = [row for row in rows if row["state_id"] == pool["state_id"]
                                and int(row["bin_index"]) == index and row["hypothesis"] == hypothesis]
                    assert len(selected) == 1
                    row = selected[0]
                    bounds = json.loads(row["bounds_json"])
                    assert bounds == {key: value for key, value in group.items()
                                      if key not in {"resolved", "uniform_width"}}
                    assert bounds[lower] == group[lower]
                    assert bounds[upper] == group[upper]
                    for name, expected in group[hypothesis].items():
                        if expected is None:
                            assert row[name] in ("", "null")
                        else:
                            assert float(row[name]) == expected


def test_source_parameters_remain_unreported_without_associated_run_summary(tmp_path, demo_export):
    _, _, demo, _ = demo_export
    source = tmp_path / "renamed_source.h5"
    shutil.copyfile(demo / "source_input.h5", source)
    output = tmp_path / "inspected"
    result = execute(output, "--source", source, "--quadrature-order", "4")
    assert result.returncode == 0, result.stderr
    report = json.loads((output / "report.json").read_text())
    assert report["parameter_status"] == "unreported"
    assert report["input_source_sha256"] == sha(source) == sha(output / "source_input.h5")
    assert report["input_config_sha256"] is None
    assert report["input_summary_sha256"] is None
    assert not (output / "surface_input.h5").exists()
    assert not (output / "input_config.json").exists()


@pytest.mark.parametrize("fault", [None, "hash", "provenance", "reservoirs"])
def test_parameter_status_requires_hash_and_provenance_association(tmp_path, demo_export, fault):
    _, config, demo, _ = demo_export
    source = demo / "source_input.h5"
    summary = {
        "schema": "mhdlab-temkin-planar-run-v1", "parameter_status": "unvalidated",
        "parameter_provenance": config["parameter_provenance"],
        "output_sha256": {"temkin_planar_release.h5": sha(source)},
    }
    if fault == "hash":
        summary["output_sha256"]["temkin_planar_release.h5"] = "0" * 64
    elif fault == "provenance":
        summary["parameter_provenance"] = "other experiment"
    elif fault == "reservoirs":
        summary["reservoirs"] = []
    summary_path = tmp_path / "source_summary.json"
    summary_path.write_text(json.dumps(summary))
    output = tmp_path / "inspected"
    result = execute(output, "--source", source, "--source-run-summary", summary_path, "--quadrature-order", "4")
    if fault is None:
        assert result.returncode == 0, result.stderr
        report = json.loads((output / "report.json").read_text())
        assert report["parameter_status"] == "unvalidated"
        assert report["input_summary_sha256"] == sha(summary_path)
    else:
        assert result.returncode != 0
        assert not output.exists()


def test_external_linked_source_data_is_rejected_before_generation(tmp_path, demo_export):
    _, _, demo, _ = demo_export
    source = tmp_path / "linked.h5"
    shutil.copyfile(demo / "source_input.h5", source)
    external = tmp_path / "external.h5"
    with h5py.File(source, "r+") as stream, h5py.File(external, "w") as dependency:
        stream.copy("temperature_k", dependency)
        del stream["temperature_k"]
        stream["temperature_k"] = h5py.ExternalLink(str(external), "temperature_k")
    output = tmp_path / "inspected"
    result = execute(output, "--source", source, "--quadrature-order", "4")
    assert result.returncode != 0
    assert not output.exists()


@pytest.mark.parametrize("existing", ["directory", "file", "symlink"])
def test_existing_destination_and_contents_are_preserved(tmp_path, existing):
    output = tmp_path / "joint"
    sentinel = b"retain previous experimental work"
    if existing == "directory":
        output.mkdir()
        (output / "sentinel").write_bytes(sentinel)
    elif existing == "file":
        output.write_bytes(sentinel)
    else:
        output.symlink_to(tmp_path / "missing-target", target_is_directory=True)
    result = execute(output, "--demo")
    assert result.returncode != 0
    assert "new path" in result.stderr
    if existing == "directory":
        assert list(output.iterdir()) == [output / "sentinel"]
        assert (output / "sentinel").read_bytes() == sentinel
    elif existing == "file":
        assert output.read_bytes() == sentinel
    else:
        assert output.is_symlink()
        assert not (tmp_path / "missing-target").exists()


@pytest.mark.parametrize("arguments", [[], ["--demo", "--source", "missing.h5"],
                                         ["--demo", "--quadrature-order", "1"]])
def test_ambiguous_source_or_unsupported_quadrature_creates_no_result(tmp_path, arguments):
    output = tmp_path / "joint"
    result = execute(output, *arguments)
    assert result.returncode != 0
    assert not output.exists()


@pytest.mark.parametrize("kind", ["duplicate", "nonfinite"])
def test_duplicate_or_nonfinite_json_creates_no_result(tmp_path, kind):
    path, config = configuration(tmp_path)
    if kind == "duplicate":
        text = path.read_text().replace('"schema":', f'"schema": {json.dumps(config["schema"])}, "schema":', 1)
    else:
        config["geometry"]["duration_s"] = float("nan")
        text = json.dumps(config)
    path.write_text(text)
    output = tmp_path / "joint"
    result = execute(output, "--demo", "--config", path)
    assert result.returncode != 0
    assert not output.exists()


@pytest.mark.parametrize("kind", ["status", "boolean_geometry", "coverage_shape", "unknown_field"])
def test_demo_cannot_mislabel_or_silently_reinterpret_input(tmp_path, kind):
    path, config = configuration(tmp_path)
    if kind == "status":
        config["parameter_status"] = "sourced"
    elif kind == "boolean_geometry":
        config["geometry"]["duration_s"] = True
    elif kind == "coverage_shape":
        config["reservoirs"][0]["initial_coverage_m2"] = [1e18]
    else:
        config["density_kg_m3"] = 7900
    path.write_text(json.dumps(config))
    output = tmp_path / "joint"
    result = execute(output, "--demo", "--config", path)
    assert result.returncode != 0
    assert not output.exists()


def test_invalid_source_schema_is_rejected_without_output(tmp_path):
    path = tmp_path / "unqualified.h5"
    with h5py.File(path, "w") as stream:
        stream.attrs["schema"] = "unqualified"
    output = tmp_path / "joint"
    result = execute(output, "--source", path)
    assert result.returncode != 0
    assert not output.exists()
