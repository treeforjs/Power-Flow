"""Integration checks for inspecting saved source-window packet associations."""

import hashlib
import importlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest

from mhdlab.constants import KB
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_planar_io import write_temkin_planar_release_history
from mhdlab.temkin_release import TemkinReservoir


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools/inspect_windowed_source_moments.py"
PRODUCER = ROOT / "tools/sample_temkin_source_windows.py"
PARAMETER_PROVENANCE = "Manufactured inspection source: constant hazard and two local temperatures"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    return path


def argv(export, output, orders=(4, 8, 16)):
    return [str(SCRIPT), "--packet-export", str(export), "--orders", *map(str, orders),
            "--output-dir", str(output)]


def run_cli(export, output, orders=(4, 8, 16)):
    return subprocess.run([sys.executable, "-W", "error", *argv(export, output, orders)],
                          cwd=ROOT, text=True, capture_output=True, timeout=45)


def assert_no_result(process, output):
    assert process.returncode != 0, process.stdout
    assert not output.exists() and not output.is_symlink()
    assert not list(output.parent.glob(".windowed-moments-*"))


@pytest.fixture(scope="module")
def packet_export(tmp_path_factory):
    directory = tmp_path_factory.mktemp("moment-inspection-source")
    base = FlashSurfaceHistory([0., 1e-7], [0., 1.], [[300.], [300.]], "synthetic",
                               "moment-inspection-test", "Manufactured prescribed temperature",
                               "Manufactured stationary boundary face")
    surface = PlanarSurfaceHistory([0., 1e-7], [-.003, -.001, .002], [.001, .004],
                                  [[[400., 900.]], [[400., 900.]]], base.attributes)
    pool = TemkinReservoir("inspection-pool", "M:inspection", 2e7, 0., .07, .04,
                           0., 1e18, 1e18, "Manufactured constant hazard",
                           "Manufactured supplied internal and desorption costs")
    history = evolve_temkin_planar_release(surface, [[[1.2e17, 2.4e17]]], [pool],
                                           parameter_provenance=PARAMETER_PROVENANCE)
    source = directory / "source.h5"
    write_temkin_planar_release_history(source, history)
    summary = write_json(directory / "source_summary.json", {
        "schema": "mhdlab-temkin-planar-run-v1", "parameter_status": "manufactured",
        "parameter_provenance": PARAMETER_PROVENANCE, "output_sha256": {source.name: sha(source)},
    })
    config = write_json(directory / "config.json", {
        "schema": "mhdlab-temkin-windowed-packet-config-v1", "parameter_status": "manufactured",
        "conditioning_provenance": "Manufactured disjoint original-interval emission windows",
        "mass_status": "manufactured", "mass_provenance": "Explicit manufactured state mass",
        "state_masses_kg": {"M:inspection": 3.4e-26},
        "selections": [{"interval_index": 0, "start_s": 1e-8, "stop_s": 3e-8},
                       {"interval_index": 0, "start_s": 6e-8, "stop_s": 9e-8}],
        "log_mass_tolerance": 1e-7, "face_normal_coordinate_m": 0.,
        "outward_normal_sign": 1, "packet_count": 64, "seed": 1019,
    })
    output = directory / "packet-export"
    process = subprocess.run([sys.executable, "-W", "error", str(PRODUCER), "--source", str(source),
                              "--source-run-summary", str(summary), "--config", str(config),
                              "--output-dir", str(output)], cwd=ROOT, text=True,
                             capture_output=True, timeout=45)
    assert process.returncode == 0, process.stderr
    return output


@pytest.fixture(scope="module")
def inspected(packet_export, tmp_path_factory):
    output = tmp_path_factory.mktemp("moment-inspection-result") / "output"
    process = run_cli(packet_export, output)
    assert process.returncode == 0, process.stderr
    report = json.loads((output / "report.json").read_text())
    estimates = {order: json.loads((output / f"moments_order_{order:04d}.json").read_text())
                 for order in report["orders"]}
    return output, report, estimates


def test_inspection_preserves_exact_association_and_hashes(packet_export, inspected):
    output, report, estimates = inspected
    old = json.loads((packet_export / "report.json").read_text())
    assert report["schema"] == "mhdlab-windowed-source-moment-inspection-v1"
    assert report["parameter_status"] == report["mass_status"] == "manufactured"
    assert report["parameter_provenance"] == PARAMETER_PROVENANCE
    assert report["conditioning_provenance"] == old["conditioning_provenance"]
    assert report["mass_provenance"] == old["mass_provenance"]
    assert report["orders"] == [4, 8, 16]
    assert report["host_budget_accepted"] is False
    assert set(report["runtime"]) == {"python", "numpy", "scipy", "h5py"}
    assert set(report["input_export_sha256"]) == set(old["output_sha256"]) | {"report.json"}
    for name, digest in report["input_export_sha256"].items():
        assert sha(packet_export / name) == digest
        assert sha(output / ("packet_export_" + name)) == digest
    assert set(report["output_sha256"]) == {path.name for path in output.iterdir()} - {"report.json"}
    for name, digest in report["output_sha256"].items():
        assert sha(output / name) == digest
    for name, digest in report["source_sha256"].items():
        assert sha(ROOT / name) == digest
    for order, estimate in estimates.items():
        assert estimate["schema"] == "mhdlab-windowed-joint-source-moments-v1"
        assert estimate["quadrature_order"] == order
        assert estimate["source"] == old["conditioned_source"]
        assert estimate["status"] == "manufactured" and estimate["host_budget_accepted"] is False
        assert len(estimate["components"]) == 4
        assert len(estimate["per_selection"]) == 2 and len(estimate["per_reservoir"]) == 1
        for component in estimate["components"]:
            assert component["interval_origin_s"] == 0.
            selected = old["conditioned_source"]["selections"][component["selection_index"]]
            assert component["window_start_offset_s"] == selected["start_s"]
            assert component["window_stop_offset_s"] == selected["stop_s"]


def test_analytic_mixture_and_finite_quadrature_changes(inspected):
    _, report, estimates = inspected
    rate = 2e7
    windows = [(1e-8, 3e-8), (6e-8, 9e-8)]
    masses, means, variances = [], [], []
    for start, stop in windows:
        length = stop-start
        denominator = math.expm1(rate*length)
        mean_delay = 1/rate-length/denominator
        second_delay = 2/rate**2-(length**2+2*length/rate)/denominator
        masses.append(math.exp(-rate*start)-math.exp(-rate*stop))
        means.append(start+mean_delay)
        variances.append(second_delay-mean_delay**2)
    total_fraction = math.fsum(masses)
    probabilities = [mass/total_fraction for mass in masses]
    exact_time = math.fsum(p*mean for p, mean in zip(probabilities, means))
    exact_variance = math.fsum(p*(var+(mean-exact_time)**2)
                               for p, mean, var in zip(probabilities, means, variances))
    finest = estimates[16]
    resolved = finest["resolved"]
    uniform = finest["matched_uniform_width"]
    assert resolved["released_number"] == pytest.approx(2.88e12*total_fraction, rel=3e-10, abs=0.)
    assert resolved["mean_time_s"] == pytest.approx(exact_time, rel=3e-10, abs=0.)
    assert resolved["variance_time_s2"] == pytest.approx(exact_variance, rel=3e-9, abs=0.)
    assert abs(resolved["mean_time_s"]-exact_time) < abs(estimates[4]["resolved"]["mean_time_s"]-exact_time)
    assert resolved["mean_x_m"] == pytest.approx(-.000125, rel=3e-14, abs=0.)
    assert resolved["variance_x_m2"] == pytest.approx(1.8177083333333333e-6, rel=3e-14, abs=0.)
    assert uniform["mean_x_m"] == pytest.approx(-.0005, rel=3e-14, abs=0.)
    assert uniform["variance_x_m2"] == pytest.approx(.005**2/12, rel=3e-14, abs=0.)
    assert resolved["mean_z_m"] == pytest.approx(.0025, rel=3e-14, abs=0.)
    assert resolved["variance_z_m2"] == pytest.approx(.003**2/12, rel=3e-14, abs=0.)
    assert resolved["mean_expected_translational_energy_j"] == pytest.approx(1550*KB, rel=3e-14, abs=0.)
    assert resolved["variance_expected_translational_energy_j2"] == pytest.approx(187500*KB**2, rel=3e-14, abs=0.)
    assert resolved["mean_conditional_kinetic_energy_variance_j2"] == pytest.approx(1295000*KB**2, rel=3e-14, abs=0.)
    assert resolved["variance_kinetic_energy_j2"] == pytest.approx(1482500*KB**2, rel=3e-14, abs=0.)
    assert resolved["second_raw_kinetic_energy_j2"] == pytest.approx(3885000*KB**2, rel=3e-14, abs=0.)
    assert resolved["cov_expected_translational_energy_x_j_m"] == pytest.approx(.46875*KB, rel=3e-14, abs=0.)
    assert uniform["cov_expected_translational_energy_x_j_m"] == 0.
    for name in ("released_number", "mean_time_s", "variance_time_s2", "mean_z_m", "variance_z_m2",
                 "mean_expected_translational_energy_j", "variance_expected_translational_energy_j2",
                 "variance_kinetic_energy_j2", "second_raw_kinetic_energy_j2"):
        assert uniform[name] == resolved[name]
    for refinement, (coarse, fine) in zip(report["refinements"], [(4, 8), (8, 16)]):
        assert (refinement["coarse_order"], refinement["refined_order"]) == (coarse, fine)
        assert "not certified" in refinement["interpretation"]
        for hypothesis in ("resolved", "matched_uniform_width"):
            for name, change in refinement["changes"][hypothesis].items():
                a, b = estimates[coarse][hypothesis][name], estimates[fine][hypothesis][name]
                assert change["signed_change"] == b-a
                assert change["absolute_change"] == abs(b-a)
                if b:
                    assert change["relative_to_refined"] == (b-a)/abs(b)
                else:
                    assert change["relative_to_refined"] is None


def test_inspector_never_samples_or_advances_source(packet_export, tmp_path, monkeypatch):
    module = importlib.import_module("tools.inspect_windowed_source_moments")
    producer = importlib.import_module("tools.sample_temkin_source_windows")
    boundary = importlib.import_module("mhdlab.temkin_windowed_boundary")
    evolution = importlib.import_module("mhdlab.temkin_planar")
    before = {path.name: sha(path) for path in packet_export.iterdir()}

    def forbidden(*args, **kwargs):
        raise AssertionError("inspection invoked sampling or source evolution")

    monkeypatch.setattr(producer, "sample_temkin_windowed_emission_packets", forbidden)
    monkeypatch.setattr(boundary, "sample_temkin_windowed_emission_packets", forbidden)
    monkeypatch.setattr(evolution, "evolve_temkin_planar_release", forbidden)
    output = tmp_path / "moments"
    monkeypatch.setattr(sys, "argv", argv(packet_export, output, (4, 8)))
    module.main()
    assert (output / "report.json").is_file()
    assert before == {path.name: sha(path) for path in packet_export.iterdir()}


def rewrite_npz(directory, report, transform):
    target = directory / "boundary_packets.npz"
    with np.load(target, allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files}
    transform(arrays)
    np.savez_compressed(target, **arrays)
    report["output_sha256"][target.name] = sha(target)


TAMPERS = ["input_hash", "external_path", "symlink", "producer_code", "report_status",
           "source_provenance", "conditioned_manifest", "npz_duplicate_metadata",
           "mass_provenance", "summary_hash", "source_run_association", "array_units",
           "missing_array", "wrong_shape"]


@pytest.mark.parametrize("tamper", TAMPERS)
def test_altered_export_is_rejected(packet_export, tmp_path, tamper):
    export = tmp_path / "export"
    shutil.copytree(packet_export, export)
    report = json.loads((export / "report.json").read_text())
    if tamper == "input_hash":
        with (export / "uniforms.npy").open("ab") as stream:
            stream.write(b"changed")
    elif tamper == "external_path":
        report["output_sha256"]["../outside.npy"] = report["output_sha256"].pop("uniforms.npy")
    elif tamper == "symlink":
        outside = tmp_path / "outside.npy"
        shutil.copyfile(export / "uniforms.npy", outside)
        (export / "uniforms.npy").unlink()
        (export / "uniforms.npy").symlink_to(outside)
    elif tamper == "producer_code":
        first = next(iter(report["source_sha256"]))
        report["source_sha256"][first] = "0"*64
    elif tamper == "report_status":
        report["parameter_status"] = "sourced"
    elif tamper == "source_provenance":
        report["source_parameter_provenance"] = "Detached provenance"
    elif tamper == "conditioned_manifest":
        report["conditioned_source"]["parameter_provenance"] = "Detached source association"
        report["packets"]["source"] = report["conditioned_source"]
        rewrite_npz(export, report, lambda arrays: arrays.update(metadata_json=np.asarray(json.dumps(report["packets"]))))
    elif tamper == "npz_duplicate_metadata":
        def duplicate(arrays):
            text = str(arrays["metadata_json"])
            schema = json.dumps(report["packets"]["schema"])
            arrays["metadata_json"] = np.asarray('{"schema":' + schema + ',' + text[1:])
        rewrite_npz(export, report, duplicate)
    elif tamper == "mass_provenance":
        report["mass_provenance"] = "Detached mass provenance"
    elif tamper == "summary_hash":
        report["input_summary_sha256"] = "0"*64
    elif tamper == "source_run_association":
        report["source_run_association"]["parameter_status"] = "unvalidated"
    elif tamper == "array_units":
        report["packets"]["array_units"]["temperature_k"] = "eV"
        rewrite_npz(export, report, lambda arrays: arrays.update(metadata_json=np.asarray(json.dumps(report["packets"]))))
    elif tamper == "missing_array":
        rewrite_npz(export, report, lambda arrays: arrays.pop("velocity_m_s"))
    elif tamper == "wrong_shape":
        rewrite_npz(export, report, lambda arrays: arrays.update(velocity_m_s=arrays["velocity_m_s"][:, :2]))
    write_json(export / "report.json", report)
    output = tmp_path / "result"
    process = run_cli(export, output, (4, 8))
    assert_no_result(process, output)
    assert "error:" in process.stderr


@pytest.mark.parametrize("orders", [(4,), (4, 4), (8, 4), (1, 4), (2, 4, 8, 16, 32, 64), ("4.5", 8)])
def test_invalid_order_requests_are_rejected(packet_export, tmp_path, orders):
    output = tmp_path / "result"
    process = run_cli(packet_export, output, orders)
    assert_no_result(process, output)
    assert "--orders" in process.stderr


@pytest.mark.parametrize("kind", ["directory", "file", "broken_symlink"])
def test_existing_output_paths_are_preserved(packet_export, tmp_path, kind):
    output = tmp_path / "result"
    if kind == "directory":
        output.mkdir()
        (output / "sentinel").write_text("retain")
    elif kind == "file":
        output.write_text("retain")
    else:
        output.symlink_to(tmp_path / "absent")
    process = run_cli(packet_export, output, (4, 8))
    assert process.returncode != 0 and "new path" in process.stderr
    if kind == "directory":
        assert (output / "sentinel").read_text() == "retain"
    elif kind == "file":
        assert output.read_text() == "retain"
    else:
        assert output.is_symlink() and not output.exists()
    assert not list(tmp_path.glob(".windowed-moments-*"))


def test_partial_write_failure_rolls_back_stage(packet_export, tmp_path, monkeypatch):
    module = importlib.import_module("tools.inspect_windowed_source_moments")
    original = Path.write_text
    output = tmp_path / "result"

    def fail_moment(path, *args, **kwargs):
        if path.name.startswith("moments_order_"):
            raise OSError("manufactured moment write failure")
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, "write_text", fail_moment)
    monkeypatch.setattr(sys, "argv", argv(packet_export, output, (4, 8)))
    with pytest.raises(SystemExit) as failure:
        module.main()
    assert failure.value.code == 2
    assert not output.exists() and not list(tmp_path.glob(".windowed-moments-*"))


@pytest.mark.parametrize("change", ["input", "implementation"])
def test_changes_during_staging_are_rejected(packet_export, tmp_path, monkeypatch, change):
    module = importlib.import_module("tools.inspect_windowed_source_moments")
    export = tmp_path / "export"
    shutil.copytree(packet_export, export)
    output = tmp_path / "result"
    original_write = Path.write_text
    original_sha = module.sha
    staged = False

    def mark_stage(path, *args, **kwargs):
        nonlocal staged
        result = original_write(path, *args, **kwargs)
        if path.name.startswith("moments_order_") and not staged:
            staged = True
            if change == "input":
                with (export / "uniforms.npy").open("ab") as stream:
                    stream.write(b"changed after copy")
        return result

    def changed_hash(path):
        if change == "implementation" and staged and Path(path) == SCRIPT:
            return "0"*64
        return original_sha(path)

    monkeypatch.setattr(Path, "write_text", mark_stage)
    monkeypatch.setattr(module, "sha", changed_hash)
    monkeypatch.setattr(sys, "argv", argv(export, output, (4, 8)))
    with pytest.raises(SystemExit) as failure:
        module.main()
    assert failure.value.code == 2 and staged
    assert not output.exists() and not list(tmp_path.glob(".windowed-moments-*"))


def test_output_appearing_during_staging_is_preserved(packet_export, tmp_path, monkeypatch):
    module = importlib.import_module("tools.inspect_windowed_source_moments")
    original = Path.write_text
    output = tmp_path / "result"

    def competing_writer(path, *args, **kwargs):
        result = original(path, *args, **kwargs)
        if path.name.startswith("moments_order_") and not output.exists():
            output.mkdir()
            original(output / "sentinel", "retain")
        return result

    monkeypatch.setattr(Path, "write_text", competing_writer)
    monkeypatch.setattr(sys, "argv", argv(packet_export, output, (4, 8)))
    with pytest.raises(SystemExit) as failure:
        module.main()
    assert failure.value.code == 2
    assert (output / "sentinel").read_text() == "retain"
    assert not list(tmp_path.glob(".windowed-moments-*"))
