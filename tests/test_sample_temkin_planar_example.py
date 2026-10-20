"""End-to-end 3D launch export preserves local patch inventories and units."""

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_planar_io import SCHEMA, write_temkin_planar_release_history
from mhdlab.temkin_release import TemkinReservoir


ROOT = Path(__file__).resolve().parents[1]


def test_planar_example_conserves_each_component_and_serializes_geometry(tmp_path):
    time, x, z = [0., 1e-8], [-.0025, 0., .0025], [0., .01, .025]
    temperature = np.array([[[650., 800.], [900., 1000.]]] * 2)
    base = FlashSurfaceHistory(time, x, temperature[:, 0], "synthetic", "planar-export", "manufactured", "face")
    surface = PlanarSurfaceHistory(time, x, z, temperature, base.attributes)
    pool = TemkinReservoir("manufactured", "M:manufactured-ground", 1e12, .8, .6, .05,
                           .35, 2e18, 2e18, "manufactured rate", "manufactured cost")
    history = evolve_temkin_planar_release(surface, np.full((1, 2, 2), 1e18), [pool],
                                           parameter_provenance="manufactured 3D export check")
    source = tmp_path / "source.h5"
    write_temkin_planar_release_history(source, history)
    record = tmp_path / "source_summary.json"
    record.write_text(json.dumps({"schema": "mhdlab-temkin-planar-run-v1", "parameter_status": "manufactured",
                                 "parameter_provenance": history.parameter_provenance,
                                 "output_sha256": {"temkin_planar_release.h5": hashlib.sha256(source.read_bytes()).hexdigest()}}))
    output = tmp_path / "packets"
    completed = subprocess.run([sys.executable, str(ROOT / "examples/sample_temkin_planar_boundary.py"),
                                "--source", str(source), "--source-run-summary", str(record),
                                "--masses", str(ROOT / "configs/manufactured_boundary_masses.json"),
                                "--interval", "0", "--packets", "512", "--face-normal-coordinate-m", "0",
                                "--normal-sign", "1", "--output-dir", str(output)],
                               cwd=ROOT, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    summary = json.loads((output / "summary.json").read_text())
    assert summary["source_schema"] == SCHEMA
    assert summary["source_parameter_status"] == "manufactured"
    assert summary["array_units"]["resolved_position_m"] == "m"
    assert summary["array_units"]["source_current_index"] == "1"
    assert summary["packet_file_sha256"] == hashlib.sha256((output / "boundary_packets.npz").read_bytes()).hexdigest()
    for name, digest in summary["source_code_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest
    with np.load(output / "boundary_packets.npz", allow_pickle=False) as arrays:
        for name in arrays.files:
            assert arrays[name].dtype.kind != "O"
        metadata = json.loads(str(arrays["metadata_json"]))
        assert metadata["source_history_schema"] == SCHEMA
        assert len(metadata["row_surface_attributes"]) == 2
        np.testing.assert_array_equal(arrays["resolved_position_m"][:, 1:], arrays["uniform_width_position_m"][:, 1:])
        np.testing.assert_array_equal(arrays["resolved_position_m"][:, 1], np.zeros(512))
        for current in range(2):
            for width in range(2):
                selected = ((arrays["source_current_index"] == current)
                            & (arrays["source_width_index"] == width))
                assert np.any(selected)
                expected = history.patch_extensive()[0, 0, current, width]
                np.testing.assert_allclose(arrays["particle_weight"][selected].sum(), expected, rtol=4e-15)
        np.testing.assert_allclose(summary["total_packet_particles"], history.extensive_totals()[0, 0], rtol=4e-15)
