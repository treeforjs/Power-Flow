"""Safe, law-labelled boundary export without downstream evolution."""

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.temkin_history import evolve_temkin_release
from mhdlab.temkin_io import write_temkin_history
from mhdlab.temkin_release import TemkinReservoir, LAW_IDENTITY
from mhdlab.temkin_emission_time import TIME_IDENTITY
from mhdlab.thermal_release import SurfaceInventory


ROOT = Path(__file__).resolve().parents[1]


def make_source(tmp_path):
    surface = FlashSurfaceHistory([0., 1e-8], [-.0025, 0., .0025], [[700., 900.], [700., 900.]],
                                  "synthetic", "export-test", "manufactured source", "stationary face")
    pool = TemkinReservoir("manufactured", "M:manufactured-ground", 1e12, .8, .6, .05,
                           .35, 2e18, 2e18, "manufactured rate", "manufactured independent cost")
    history = evolve_temkin_release(surface, SurfaceInventory([[1e18, 5e17]]), [pool],
                                    parameter_provenance="manufactured export integration")
    path = tmp_path / "source.h5"
    write_temkin_history(path, history)
    return path, history


def test_weighted_packet_cli_exports_only_safe_arrays_and_preserves_law_and_moments(tmp_path):
    source, history = make_source(tmp_path)
    output = tmp_path / "result"
    mass_path = ROOT / "configs/manufactured_boundary_masses.json"
    run_record = tmp_path / "source_summary.json"
    run_record.write_text(json.dumps({"schema": "mhdlab-temkin-run-v1", "parameter_status": "manufactured",
                                      "parameter_provenance": history.parameter_provenance,
                                      "output_sha256": {"temkin_release.h5": hashlib.sha256(source.read_bytes()).hexdigest()}}))
    command = [sys.executable, str(ROOT / "tools/sample_temkin_boundary.py"),
               "--source", str(source), "--masses", str(mass_path), "--interval", "0",
               "--source-run-summary", str(run_record), "--packets", "2048", "--seed", "13", "--output-dir", str(output)]
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    summary = json.loads((output / "summary.json").read_text())
    assert summary["law_identity"] == LAW_IDENTITY
    assert summary["time_model"] == TIME_IDENTITY
    assert summary["mass_configuration"]["mass_status"] == "manufactured"
    assert summary["source_parameter_status"] == "manufactured"
    assert summary["source_run_record"]["sha256"] == hashlib.sha256(run_record.read_bytes()).hexdigest()
    assert summary["input_source"]["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert summary["input_masses"]["sha256"] == hashlib.sha256(mass_path.read_bytes()).hexdigest()
    assert summary["packet_file_sha256"] == hashlib.sha256((output / "boundary_packets.npz").read_bytes()).hexdigest()
    for name, digest in summary["source_code_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest
    with np.load(output / "boundary_packets.npz", allow_pickle=False) as archive:
        for name in archive.files:
            assert archive[name].dtype.kind != "O"
        metadata = json.loads(str(archive["metadata_json"]))
        assert metadata["law_identity"] == LAW_IDENTITY
        assert metadata["time_model"] == TIME_IDENTITY
        assert np.all(archive["velocity_m_s"][:, 1] > 0)
        assert len(archive["particle_weight"]) == 2048
        length = summary["mass_configuration"]["current_direction_length_m"]
        for width in range(2):
            selected = archive["source_width_index"] == width
            expected = history.released_m2[0, 0, width] * np.diff(history.x_edges_m)[width] * length
            np.testing.assert_allclose(archive["particle_weight"][selected].sum(), expected, rtol=4e-15)
        physical_energy = np.sum(.5 * archive["mass_kg"] * archive["particle_weight"]
                                 * np.sum(archive["velocity_m_s"] ** 2, axis=1))
        np.testing.assert_allclose(physical_energy, summary["thermal_energy"]["sampled_translational_energy_j"], rtol=4e-15)
    before = {path.name: path.read_bytes() for path in output.iterdir()}
    denied = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    assert denied.returncode != 0
    assert before == {path.name: path.read_bytes() for path in output.iterdir()}


def test_mass_state_mismatch_leaves_no_apparently_successful_export(tmp_path):
    source, _ = make_source(tmp_path)
    mass_config = json.loads((ROOT / "configs/manufactured_boundary_masses.json").read_text())
    mass_config["state_masses_kg"] = {"M:unrelated-state": 3e-26}
    masses = tmp_path / "mass.json"
    masses.write_text(json.dumps(mass_config))
    output = tmp_path / "result"
    completed = subprocess.run([sys.executable, str(ROOT / "tools/sample_temkin_boundary.py"),
                                "--source", str(source), "--masses", str(masses), "--interval", "0",
                                "--packets", "20", "--output-dir", str(output)],
                               cwd=ROOT, capture_output=True, text=True)
    assert completed.returncode != 0
    assert "mass" in completed.stderr.lower()
    assert not output.exists()


def test_unassociated_source_status_cannot_be_presented_as_bound_provenance(tmp_path):
    source, history = make_source(tmp_path)
    record = tmp_path / "source_summary.json"
    record.write_text(json.dumps({"schema": "mhdlab-temkin-run-v1", "parameter_status": "sourced",
                                 "parameter_provenance": history.parameter_provenance,
                                 "output_sha256": {"temkin_release.h5": "0" * 64}}))
    output = tmp_path / "result"
    command = [sys.executable, str(ROOT / "tools/sample_temkin_boundary.py"), "--source", str(source),
               "--masses", str(ROOT / "configs/manufactured_boundary_masses.json"), "--interval", "0",
               "--packets", "20", "--source-run-summary", str(record), "--output-dir", str(output)]
    denied = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    assert denied.returncode != 0
    assert "source hash" in denied.stderr
    assert not output.exists()
    # Without associated evidence, this output states that the parameter label
    # is unreported; it never infers "sourced" from a mass or heating label.
    position = command.index("--source-run-summary")
    del command[position:position + 2]
    allowed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    assert allowed.returncode == 0, allowed.stderr
    summary = json.loads((output / "summary.json").read_text())
    assert summary["source_parameter_status"] is None
    assert summary["source_run_record"] is None
