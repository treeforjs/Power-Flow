"""End-to-end checks for the distinct coverage-dependent source workflow."""

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from mhdlab.flash_surface import FlashSurfaceHistory, write_surface_history
from mhdlab.surface_release_io import load_release_history
from mhdlab.temkin_io import load_temkin_history
from mhdlab.temkin_release import LAW_IDENTITY


ROOT = Path(__file__).resolve().parents[1]


def run(*args):
    return subprocess.run([sys.executable, str(ROOT / "tools/run_temkin_release.py"), *map(str, args)],
                          cwd=ROOT, text=True, capture_output=True, check=False)


def config(tmp_path, **changes):
    values = json.loads((ROOT / "configs/temkin_manufactured.json").read_text())
    values.update(changes)
    path = tmp_path / "parameters.json"
    path.write_text(json.dumps(values))
    return path


def small_surface(tmp_path, **changes):
    surface = FlashSurfaceHistory([0, 1e-7], [-0.002, 0, 0.002],
                                  [[600, 800], [700, 1000]], "synthetic", "manufactured",
                                  "Numerical interface check", "stationary face")
    path = tmp_path / "surface.h5"
    write_surface_history(path, replace(surface, **changes))
    return path


def test_manufactured_controls_conserve_budgets_and_pin_output_bytes(tmp_path):
    output = tmp_path / "result"
    p = run("--demo", "--constant-barrier-control", "--output-dir", output)
    assert p.returncode == 0, p.stderr
    summary = json.loads((output / "summary.json").read_text())
    history = load_temkin_history(output / "temkin_release.h5")
    baseline = load_release_history(output / "constant_barrier_control.h5")
    assert summary["parameter_status"] == "manufactured"
    assert summary["law_identity"] == LAW_IDENTITY == history.law_identity
    length = summary["current_direction_length_m"]
    counts = history.width_integral().sum(axis=0) * length
    np.testing.assert_allclose(summary["resolved_source"]["released_particles"], counts, rtol=1e-14)
    np.testing.assert_allclose(summary["resolved_source"]["wall_energy_debit_j"],
                               history.width_integral("total_wall_energy_debit_j_m2").sum(axis=0) * length,
                               rtol=1e-14)
    # At the same frozen temperature, the nonnegative coverage correction can
    # only accelerate cumulative release relative to the zero-correction law.
    assert np.all(history.coverage_m2[-1] <= baseline.coverage_m2[-1] * (1 + 1e-12))
    assert np.any(history.coverage_m2[-1] < baseline.coverage_m2[-1])
    for name, expected in summary["output_sha256"].items():
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == expected
    for name, expected in summary["source_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected
    # Distinct model identity prevents silently applying the ordinary
    # first-order conditional packet-time CDF to this file.
    with pytest.raises(ValueError, match="schema"):
        load_release_history(output / "temkin_release.h5")


def test_material_plane_requires_named_boundary_proxy(tmp_path):
    surface = small_surface(tmp_path, temperature_location="material_sample_plane", material_sample_depth_m=1e-6)
    parameters = config(tmp_path)
    output = tmp_path / "result"
    denied = run("--surface", surface, "--config", parameters, "--output-dir", output)
    assert denied.returncode != 0
    assert not output.exists()
    allowed = run("--surface", surface, "--config", parameters, "--output-dir", output,
                  "--boundary-proxy-reference", "Manufactured interface qualification only")
    assert allowed.returncode == 0, allowed.stderr
    history = load_temkin_history(output / "temkin_release.h5")
    assert history.surface_attributes["temperature_location"] == "boundary_temperature_proxy"
    assert history.surface_attributes["material_sample_depth_m"] == 1e-6
    summary = json.loads((output / "summary.json").read_text())
    assert summary["surface_input_sha256"] == hashlib.sha256(surface.read_bytes()).hexdigest()


def test_current_average_cannot_silently_drive_local_nonlinear_release(tmp_path):
    surface = small_surface(tmp_path, current_averaging_interval_m=(0.0, 0.025))
    parameters = config(tmp_path)
    output = tmp_path / "result"
    denied = run("--surface", surface, "--config", parameters, "--output-dir", output)
    assert denied.returncode != 0
    assert not output.exists()
    allowed = run("--surface", surface, "--config", parameters, "--output-dir", output,
                  "--current-uniformity-reference", "Manufactured current-uniform surface")
    assert allowed.returncode == 0, allowed.stderr
    assert load_temkin_history(output / "temkin_release.h5").surface_attributes["current_uniformity_reference"]


@pytest.mark.parametrize("change", ["inventory", "schema", "boolean_subdivisions"])
def test_unqualified_inputs_do_not_leave_success_artifacts(tmp_path, change):
    path = config(tmp_path)
    values = json.loads(path.read_text())
    if change == "inventory":
        values["reservoirs"][0]["initial_coverage_m2"] = 3e18
    elif change == "schema":
        values["schema"] = "mhdlab-surface-release-config-v1"
    else:
        values["subdivisions_per_interval"] = True
    path.write_text(json.dumps(values))
    output = tmp_path / "result"
    p = run("--demo", "--config", path, "--output-dir", output)
    assert p.returncode != 0
    assert not output.exists()


def test_existing_outputs_remain_byte_identical(tmp_path):
    output = tmp_path / "result"
    output.mkdir()
    sentinel = output / "sentinel"
    sentinel.write_bytes(b"Retain this result exactly\x00\xff")
    before = sentinel.read_bytes()
    p = run("--demo", "--output-dir", output)
    assert p.returncode != 0
    assert "existing" in p.stderr
    assert list(output.iterdir()) == [sentinel]
    assert sentinel.read_bytes() == before
