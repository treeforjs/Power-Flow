"""End-to-end packet artifact checks with manufactured source/mass inputs."""

import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np

from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.surface_release_history import evolve_surface_release
from mhdlab.surface_release_io import write_release_history
from mhdlab.thermal_release import BindingReservoir, SurfaceInventory


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools/sample_surface_boundary.py"
CONFIG = ROOT / "configs/manufactured_boundary_masses.json"


def source_file(path):
    face = FlashSurfaceHistory([0, .25], [0, .001, .004], [[500, 900], [500, 900]],
                               "synthetic", "manufactured-cli-check", "analytic input", "one gap-facing face")
    pool = BindingReservoir("manufactured", "M:manufactured-ground", 2, 0, .3, .1)
    history = evolve_surface_release(face, SurfaceInventory([[1e18, 2e18]]), [pool],
                                     parameter_provenance="manufactured CLI verification parameters")
    write_release_history(path, history)
    return history


def test_cli_writes_safe_paired_packet_contract_and_conserves_source_components(tmp_path):
    source = tmp_path / "surface_release.h5"
    history = source_file(source)
    output = tmp_path / "packets"
    result = subprocess.run([sys.executable, str(TOOL), "--source", str(source), "--masses", str(CONFIG),
                             "--interval", "0", "--packets", "32", "--seed", "7", "--normal-sign", "-1",
                             "--output-dir", str(output)], check=True, capture_output=True, text=True)
    assert json.loads(result.stdout)["mass_status"] == "manufactured"
    summary = json.loads((output / "summary.json").read_text())
    assert summary["surface_attributes"]["source_kind"] == "synthetic"
    assert summary["mass_config"]["mass_status"] == "manufactured"
    assert summary["sampling"]["outward_normal_sign"] == -1
    assert summary["input_files"]["source_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    archive = output / "boundary_packets.npz"
    assert summary["packet_file_sha256"] == hashlib.sha256(archive.read_bytes()).hexdigest()
    with np.load(archive, allow_pickle=False) as packets:
        assert packets["state_ids"].dtype.kind == "U"
        assert packets["state_ids"].tolist() == ["M:manufactured-ground"] * 32
        assert packets["schema"].item() == summary["schema"]
        assert np.all(packets["velocity_m_s"][:, 1] < 0)
        np.testing.assert_array_equal(packets["momentum_kg_m_s"], packets["velocity_m_s"] * packets["mass_kg"][:, None])
        assert np.any(packets["resolved_x_m"] != packets["uniform_x_m"])
        np.testing.assert_allclose(packets["time_s"], packets["time_start_s"] + packets["time_offset_s"], rtol=0, atol=0)
        for width in (0, 1):
            selected = packets["source_width_index"] == width
            expected = history.released_m2[0, 0, width] * np.diff(history.x_edges_m)[width] * .025
            assert math.isclose(math.fsum(packets["particle_weight"][selected]), expected, rel_tol=2e-15)
        assert math.isclose(math.fsum(packets["particle_weight"]), summary["total_released_particles"], rel_tol=2e-15)
        np.testing.assert_array_equal(packets["component_translational_energy_sampling_residual_j"],
                                       packets["component_sampled_translational_energy_j"]
                                       - packets["component_expected_translational_energy_j"])
    assert "sampling error" in summary["thermal_energy"]["interpretation"]
    assert "without resampling" in summary["comparison"]
    assert all(abs(value) < 1e-13 * summary["total_released_particles"]
               for value in summary["component_particle_balance_residual"])


def test_cli_seed_reproduces_arrays_and_unknown_mass_state_fails_before_writing(tmp_path):
    source = tmp_path / "source.h5"
    source_file(source)
    outputs = [tmp_path / "first", tmp_path / "second"]
    for output in outputs:
        subprocess.run([sys.executable, str(TOOL), "--source", str(source), "--masses", str(CONFIG),
                        "--interval", "0", "--packets", "16", "--seed", "42", "--output-dir", str(output)],
                       check=True, capture_output=True, text=True)
    with np.load(outputs[0] / "boundary_packets.npz", allow_pickle=False) as first:
        with np.load(outputs[1] / "boundary_packets.npz", allow_pickle=False) as second:
            assert first.files == second.files
            for key in first.files:
                np.testing.assert_array_equal(first[key], second[key])
    config = json.loads(CONFIG.read_text())
    config["state_masses_kg"] = {"M:unknown-state": 3e-26}
    masses = tmp_path / "unknown.json"
    masses.write_text(json.dumps(config))
    rejected_output = tmp_path / "rejected"
    rejected = subprocess.run([sys.executable, str(TOOL), "--source", str(source), "--masses", str(masses),
                               "--interval", "0", "--packets", "16", "--output-dir", str(rejected_output)],
                              capture_output=True, text=True)
    assert rejected.returncode != 0
    assert "missing exact-state mass" in rejected.stderr
    assert not rejected_output.exists()
