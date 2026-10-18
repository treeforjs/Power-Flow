"""Bounded executable workflow checks using manufactured heating and laws."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.planar_release_io import load_planar_release_history, load_planar_surface_history, write_planar_surface_history


SCRIPT = Path(__file__).resolve().parents[1] / "tools" / "run_planar_release.py"


def config_file(tmp_path, *, coverage=1e18, length=.025):
    path = tmp_path / "parameters.json"
    value = {
        "schema": "mhdlab-surface-release-config-v1", "parameter_status": "manufactured",
        "parameter_provenance": "manufactured planar CLI law; no material fit",
        "current_direction_length_m": length, "subdivisions_per_interval": 1,
        "reservoirs": [{"name": "verification-pool", "state_id": "X:specified", "attempt_frequency_s": 1e3,
                        "activation_energy_ev": .2, "desorption_energy_ev": .3, "internal_energy_ev": .1,
                        "initial_coverage_m2": coverage}],
    }
    path.write_text(json.dumps(value))
    return path


def surface_file(tmp_path, *, plane=False):
    options = {"temperature_location": "material_sample_plane", "material_sample_depth_m": 10e-6} if plane else {}
    base = FlashSurfaceHistory([0., 1e-7], [0., .005], [[300.], [900.]], "synthetic", "external-fixture",
                               "manufactured external planar input", "test face", **options)
    surface = PlanarSurfaceHistory([0., 1e-7], [0., .001, .005], [0., .01, .025],
                                   np.array([[[300., 400.], [500., 600.]], [[800., 1000.], [1200., 1400.]]]),
                                   base.attributes)
    path = tmp_path / "external.h5"
    write_planar_surface_history(path, surface)
    return path


def invoke(*arguments):
    return subprocess.run([sys.executable, str(SCRIPT), *(str(value) for value in arguments)],
                          text=True, capture_output=True)


def test_demo_writes_labeled_source_controls_hashes_csv_and_scientific_plot(tmp_path):
    config = config_file(tmp_path)
    output = tmp_path / "demo"
    result = invoke("--demo", "--config", config, "--output-dir", output,
                    "--include-full-face-uniform", "--plot")
    assert result.returncode == 0, result.stderr
    stdout = json.loads(result.stdout)
    assert stdout["source_kind"] == "synthetic"
    assert stdout["parameter_status"] == "manufactured"
    summary = json.loads((output / "summary.json").read_text())
    assert summary["uniform_controls"]["default"] == "width_only"
    assert summary["uniform_controls"]["full_face_included"] is True
    assert summary["input_surface"] is None
    assert "mhdlab/hdf_integrity.py" in summary["source_sha256"]
    for name, digest in summary["artifacts_sha256"].items():
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == digest
    surface = load_planar_surface_history(output / "surface_input.h5")
    assert np.ptp(surface.temperature_k[-1, :, 0]) > 0
    assert np.ptp(surface.temperature_k[-1, 0, :]) > 0
    history = load_planar_release_history(output / "planar_release.h5")
    np.testing.assert_allclose(summary["released_particles"], history.extensive_totals().sum(axis=0), rtol=2e-15, atol=0.0)
    with (output / "integrated_source.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == summary["intervals"]
    assert float(rows[0]["released_particles"]) >= 0.0
    assert (output / "overview.png").stat().st_size > 1000


def test_external_surface_uses_exact_per_reservoir_z_x_inventory(tmp_path):
    initial = [[1e18, 2e18], [3e18, 4e18]]
    config = config_file(tmp_path, coverage=initial)
    source = surface_file(tmp_path)
    output = tmp_path / "external_run"
    result = invoke("--surface", source, "--config", config, "--output-dir", output)
    assert result.returncode == 0, result.stderr
    history = load_planar_release_history(output / "planar_release.h5")
    np.testing.assert_array_equal(history.coverage_m2[0, 0], initial)
    summary = json.loads((output / "summary.json").read_text())
    assert summary["input_surface"]["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert summary["uniform_controls"]["full_face_included"] is False


def test_material_plane_requires_explicit_proxy_reference_and_retains_depth(tmp_path):
    source = surface_file(tmp_path, plane=True)
    config = config_file(tmp_path)
    output = tmp_path / "proxy_run"
    rejected = invoke("--surface", source, "--config", config, "--output-dir", output)
    assert rejected.returncode != 0
    assert "not a boundary temperature" in rejected.stderr
    assert not (output / "planar_release.h5").exists()
    accepted = invoke("--surface", source, "--config", config, "--output-dir", output,
                      "--boundary-proxy-reference", "specified gradient/convergence evidence")
    assert accepted.returncode == 0, accepted.stderr
    history = load_planar_release_history(output / "planar_release.h5")
    assert history.surface_attributes["material_sample_depth_m"] == 10e-6
    assert history.surface_attributes["boundary_temperature_reference"] == "specified gradient/convergence evidence"
    assert history.surface_attributes["temperature_location"] == "boundary_temperature_proxy"


def test_external_input_requires_config_and_rejects_implicit_inventory_broadcast(tmp_path):
    source = surface_file(tmp_path)
    missing = invoke("--surface", source, "--output-dir", tmp_path / "missing_config")
    assert missing.returncode != 0
    assert "requires an explicit --config" in missing.stderr
    config = config_file(tmp_path, coverage=[1e18, 2e18])
    invalid = invoke("--surface", source, "--config", config, "--output-dir", tmp_path / "wrong_shape")
    assert invalid.returncode != 0
    assert "exact (nz,nx)" in invalid.stderr


def test_config_length_cannot_rescale_surface_geometry_or_overwrite_inputs(tmp_path):
    source = surface_file(tmp_path)
    config = config_file(tmp_path, length=.1)
    result = invoke("--surface", source, "--config", config, "--output-dir", tmp_path / "bad_length")
    assert result.returncode != 0
    assert "conflicts with the actual planar z extent" in result.stderr
    config = config_file(tmp_path)
    source.rename(tmp_path / "surface_input.h5")
    source = tmp_path / "surface_input.h5"
    checksum = hashlib.sha256(source.read_bytes()).hexdigest()
    result = invoke("--surface", source, "--config", config, "--output-dir", tmp_path)
    assert result.returncode != 0
    assert "must not overwrite" in result.stderr
    assert hashlib.sha256(source.read_bytes()).hexdigest() == checksum
