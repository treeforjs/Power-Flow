import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from mhdlab.flash_surface import FlashSurfaceHistory, write_surface_history
from mhdlab.surface_release_io import load_release_history


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools/run_surface_release.py"
CONFIG = ROOT / "configs/surface_release_manufactured.json"


def test_manufactured_pipeline_writes_conservative_labeled_source(tmp_path):
    result = subprocess.run([sys.executable, str(TOOL), "--demo", "--subdivisions", "1",
                             "--output-dir", str(tmp_path)], check=True, capture_output=True, text=True)
    assert json.loads(result.stdout)["parameter_status"] == "manufactured"
    summary = json.loads((tmp_path / "summary.json").read_text())
    history = load_release_history(tmp_path / "surface_release.h5")
    assert summary["surface_attributes"]["source_kind"] == "synthetic"
    assert summary["surface_attributes"]["temperature_location"] == "surface_boundary"
    assert np.all(history.coverage_m2[-1] <= history.coverage_m2[0])
    residual = np.abs(summary["particle_balance_residual"]) / np.array(summary["initial_particles"])
    assert np.max(residual) < 1e-13
    assert (tmp_path / "integrated_source.csv").is_file()


def test_plane_and_current_averages_require_separate_explicit_proxy_references(tmp_path):
    source = tmp_path / "plane.h5"
    surface = FlashSurfaceHistory([0, 1e-8], [0, 0.005], [[500], [900]], "synthetic",
                                 "test", "manufactured constant-normal/current fixture", "face", "",
                                 "material_sample_plane", 1e-5, "", (0, 0.025), "")
    write_surface_history(source, surface)
    base = [sys.executable, str(TOOL), "--surface", str(source), "--config", str(CONFIG),
            "--output-dir", str(tmp_path / "release")]
    missing = subprocess.run(base, capture_output=True, text=True)
    assert missing.returncode != 0
    assert "material sample plane" in missing.stderr
    only_boundary = subprocess.run(base + ["--boundary-proxy-reference", "manufactured field is constant in normal"],
                                   capture_output=True, text=True)
    assert only_boundary.returncode != 0
    assert "current-averaged" in only_boundary.stderr
    subprocess.run(base + ["--boundary-proxy-reference", "manufactured field is constant in normal",
                           "--current-uniformity-reference", "manufactured field is constant along current"],
                   check=True, capture_output=True, text=True)
    restored = load_release_history(tmp_path / "release/surface_release.h5")
    assert restored.surface_attributes["temperature_location"] == "boundary_temperature_proxy"
    assert restored.surface_attributes["material_sample_depth_m"] == 1e-5
    assert json.loads(restored.surface_attributes["current_averaging_interval_m_json"]) == [0, 0.025]
    assert restored.surface_attributes["current_uniformity_reference"]
