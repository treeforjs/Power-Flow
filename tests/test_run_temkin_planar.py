"""Planar workflow verifies extensive units and safe source qualification."""

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from mhdlab.planar_release_io import write_planar_surface_history
from mhdlab.temkin_planar_io import load_temkin_planar_release_history
from mhdlab.temkin_release import LAW_IDENTITY
from tools.run_planar_release import manufactured_surface


ROOT = Path(__file__).resolve().parents[1]


def test_manufactured_planar_workflow_uses_actual_patch_areas_and_pins_outputs(tmp_path):
    output = tmp_path / "result"
    completed = subprocess.run([sys.executable, str(ROOT / "tools/run_temkin_planar.py"),
                                "--demo", "--include-full-face-uniform", "--output-dir", str(output)],
                               cwd=ROOT, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    summary = json.loads((output / "summary.json").read_text())
    history = load_temkin_planar_release_history(output / "temkin_planar_release.h5")
    assert summary["law_identity"] == LAW_IDENTITY == history.law_identity
    assert summary["parameter_status"] == "manufactured"
    np.testing.assert_allclose(summary["face_area_m2"], .005 * .025, rtol=2e-15)
    counts = history.extensive_totals().sum(axis=0)
    np.testing.assert_allclose(summary["released_particles"], counts, rtol=3e-15)
    for name, digest in summary["source_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest
    for name, digest in summary["output_sha256"].items():
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == digest
    width_control = history.matched_uniform_width()["released_m2"]
    widths = np.diff(history.x_edges_m)
    np.testing.assert_allclose(np.sum(width_control * widths, axis=-1),
                               np.sum(history.released_m2 * widths, axis=-1), rtol=4e-15)
    assert not np.allclose(history.temperature_k[-1, 0], history.temperature_k[-1, -1])


def test_inconsistent_configured_current_length_cannot_rescale_resolved_source(tmp_path):
    original = manufactured_surface()
    surface = replace(original, z_edges_m=original.z_edges_m * 2)
    path = tmp_path / "surface.h5"
    write_planar_surface_history(path, surface)
    output = tmp_path / "result"
    completed = subprocess.run([sys.executable, str(ROOT / "tools/run_temkin_planar.py"),
                                "--surface", str(path), "--config", str(ROOT / "configs/temkin_manufactured.json"),
                                "--output-dir", str(output)], cwd=ROOT, capture_output=True, text=True)
    assert completed.returncode != 0
    assert "extent" in completed.stderr
    assert not output.exists()
