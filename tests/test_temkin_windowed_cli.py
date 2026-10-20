"""Portable integration checks for paired window-conditioned source exports."""

import hashlib
import importlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import h5py
import numpy as np
import pytest

from mhdlab.constants import KB, QE
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.temkin_emission_window import ReducedTemkinEmissionWindow
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_planar_io import write_temkin_planar_release_history
from mhdlab.temkin_release import TemkinReservoir


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools/sample_temkin_source_windows.py"
PARAMETER_PROVENANCE = "Manufactured CLI source with explicit independent pools"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def configuration():
    return {
        "schema": "mhdlab-temkin-windowed-packet-config-v1", "parameter_status": "manufactured",
        "conditioning_provenance": "Manufactured original-interval source-emission selections",
        "mass_status": "manufactured", "mass_provenance": "Explicit manufactured exact-state masses",
        "state_masses_kg": {"M:a": 3e-26, "M:b": 4e-26},
        "selections": [{"interval_index": 0, "start_s": 3e-9, "stop_s": 15e-9},
                       {"interval_index": 1, "start_s": 2e-9, "stop_s": 18e-9}],
        "log_mass_tolerance": 1e-7, "face_normal_coordinate_m": .001,
        "outward_normal_sign": -1, "packet_count": 96, "seed": 2025,
    }


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")
    return path


@pytest.fixture(scope="module")
def source_fixture(tmp_path_factory):
    directory = tmp_path_factory.mktemp("windowed-cli-source")
    base = FlashSurfaceHistory([0., 4e-8], [0., 1.], [[300.], [300.]], "synthetic",
                               "windowed-cli-test", "Manufactured fixed-patch temperature history",
                               "Manufactured stationary single face")
    temperatures = np.array([[[500., 800.], [650., 950.]], [[600., 900.], [750., 1050.]],
                             [[700., 1000.], [850., 1150.]]])
    surface = PlanarSurfaceHistory([0., 2e-8, 4e-8], [-.002, -.0005, .002], [0., .004, .01],
                                   temperatures, base.attributes)
    pools = [TemkinReservoir("a", "M:a", 2e8, .12, .08, 0., .3, 1e18, 1e18,
                             "Manufactured independent pool a", "Manufactured supplied wall cost"),
             TemkinReservoir("b", "M:b", 1e8, .18, .11, .03, .2, 1e18, 1e18,
                             "Manufactured independent pool b", "Manufactured supplied internal/wall costs")]
    initial = np.array([[[2e17, 1.1e17], [3e17, 2.2e17]], [[1.5e17, 2.5e17], [1e17, 1.8e17]]])
    history = evolve_temkin_planar_release(surface, initial, pools, parameter_provenance=PARAMETER_PROVENANCE)
    source = directory / "source.h5"
    write_temkin_planar_release_history(source, history)
    summary = write_json(directory / "source_summary.json", {
        "schema": "mhdlab-temkin-planar-run-v1", "parameter_status": "manufactured",
        "parameter_provenance": PARAMETER_PROVENANCE, "output_sha256": {"source.h5": sha(source)},
    })
    return {"source": source, "summary": summary, "history": history, "directory": directory}


def command(source, config, output, summary):
    arguments = [str(SCRIPT), "--source", str(source), "--config", str(config), "--output-dir", str(output)]
    if summary is not None:
        arguments += ["--source-run-summary", str(summary)]
    return arguments


def run_cli(source_fixture, tmp_path, *, config=None, summary="default", output=None, raw=None):
    config_path = tmp_path / "config.json"
    if raw is None:
        write_json(config_path, configuration() if config is None else config)
    else:
        config_path.write_text(raw)
    summary = source_fixture["summary"] if summary == "default" else summary
    output = tmp_path / "output" if output is None else output
    process = subprocess.run([sys.executable, "-W", "error", *command(source_fixture["source"], config_path, output, summary)],
                              cwd=ROOT, text=True, capture_output=True, timeout=30)
    return process, output, config_path


def assert_failed_without_stage(process, output):
    assert process.returncode != 0, process.stdout
    assert not output.exists() and not output.is_symlink()
    assert not list(output.parent.glob(".temkin-windowed-*"))


def read_export(output):
    report = json.loads((output / "report.json").read_text())
    with np.load(output / "boundary_packets.npz", allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files}
    return report, arrays, np.load(output / "uniforms.npy", allow_pickle=False)


@pytest.fixture(scope="module")
def exported(source_fixture, tmp_path_factory):
    directory = tmp_path_factory.mktemp("windowed-cli-export")
    process, output, config = run_cli(source_fixture, directory)
    assert process.returncode == 0, process.stderr
    return {"output": output, "config": config, "data": read_export(output)}


def test_export_hash_association_units_and_saved_open_variates(source_fixture, exported):
    output = exported["output"]
    report, arrays, uniforms = exported["data"]
    assert report["schema"] == "mhdlab-temkin-windowed-packet-export-v1"
    assert report["parameter_status"] == report["mass_status"] == "manufactured"
    assert report["source_parameter_provenance"] == PARAMETER_PROVENANCE
    assert report["input_source_sha256"] == sha(source_fixture["source"]) == sha(output / "source_input.h5")
    assert report["input_config_sha256"] == sha(exported["config"]) == sha(output / "input_config.json")
    assert report["input_summary_sha256"] == sha(source_fixture["summary"]) == sha(output / "source_run_summary.json")
    assert report["source_run_association"]["parameter_status"] == "manufactured"
    assert report["host_budget_accepted"] is False
    assert set(report["output_sha256"]) == {"input_config.json", "source_input.h5", "source_run_summary.json",
                                            "boundary_packets.npz", "uniforms.npy"}
    for name, digest in report["output_sha256"].items():
        assert sha(output / name) == digest
    for name, digest in report["source_sha256"].items():
        assert sha(ROOT / name) == digest
    assert uniforms.shape == (96, 7) and np.all((uniforms > 0) & (uniforms < 1))
    np.testing.assert_array_equal(uniforms, np.random.default_rng(2025).random((96, 7)))
    metadata = json.loads(str(arrays["metadata_json"]))
    assert metadata == report["packets"]
    assert set(arrays) - {"metadata_json", "state_ids", "reservoir_names"} <= set(metadata["array_units"])


def test_saved_packets_preserve_original_timing_mixture_pairing_and_distinct_energy_ledgers(source_fixture, exported):
    report, arrays, uniforms = exported["data"]
    history, config = source_fixture["history"], configuration()
    component = arrays["component_index"]
    interval, pool, zi, xi = (arrays[name] for name in
                             ("interval_index", "reservoir_index", "source_current_index", "source_width_index"))
    temperatures = history.temperature_k[interval, zi, xi]
    np.testing.assert_array_equal(arrays["temperature_k"], temperatures)
    np.testing.assert_array_equal(arrays["interval_origin_s"], history.time_edges_s[interval])
    selections = config["selections"]
    starts = np.array([selections[int(index)]["start_s"] for index in arrays["selection_index"]])
    stops = np.array([selections[int(index)]["stop_s"] for index in arrays["selection_index"]])
    np.testing.assert_array_equal(arrays["window_start_offset_s"], starts)
    np.testing.assert_array_equal(arrays["window_stop_offset_s"], stops)
    delay = arrays["local_time_delay_s"]
    assert np.all((delay > 0) & (delay < stops-starts))
    np.testing.assert_allclose(arrays["log_local_time_delay_s"], np.log(delay), rtol=2e-14, atol=0.)
    np.testing.assert_allclose(arrays["best_effort_time_s"], history.time_edges_s[interval]+starts+delay,
                               rtol=4e-15, atol=0.)
    assert np.all(arrays["global_time_precision_qualified"])
    expected_states = np.array([history.reservoirs[int(index)].state_id for index in pool])
    np.testing.assert_array_equal(arrays["state_ids"], expected_states)
    np.testing.assert_array_equal(arrays["mass_kg"], [config["state_masses_kg"][state] for state in expected_states])
    resolved, uniform = arrays["resolved_position_m"], arrays["uniform_width_position_m"]
    np.testing.assert_array_equal(resolved[:, 1:], uniform[:, 1:])
    assert np.all(resolved[:, 1] == config["face_normal_coordinate_m"])
    assert np.any(resolved[:, 0] != uniform[:, 0])
    width = np.diff(history.x_edges_m)
    np.testing.assert_allclose(resolved[:, 0], history.x_edges_m[xi]+width[xi]*uniforms[:, 2], rtol=2e-13, atol=2e-18)
    np.testing.assert_allclose(uniform[:, 0], history.x_edges_m[0]+np.ptp(history.x_edges_m)*uniforms[:, 2],
                               rtol=2e-13, atol=2e-18)
    np.testing.assert_allclose(resolved[:, 2], history.z_edges_m[zi]+np.diff(history.z_edges_m)[zi]*uniforms[:, 6],
                               rtol=2e-13, atol=2e-18)
    velocity = arrays["velocity_m_s"].astype(np.longdouble)
    assert np.all(velocity[:, 1] < 0)
    np.testing.assert_array_equal(arrays["momentum_kg_m_s"], arrays["mass_kg"][:, None]*arrays["velocity_m_s"])
    energy = arrays["mass_kg"].astype(np.longdouble) * np.sum(velocity**2, axis=1, dtype=np.longdouble) / 2
    from_draws = np.longdouble(KB)*temperatures.astype(np.longdouble)*(
        -np.log1p(-uniforms[:, 3].astype(np.longdouble))-np.log(uniforms[:, 4].astype(np.longdouble)))
    np.testing.assert_allclose(energy, from_draws, rtol=4e-14, atol=0.)
    source_components = report["conditioned_source"]["components"]
    expected_counts = np.array([item["source_window"]["quantities"]["released_number"]["value"]
                                for item in source_components])
    np.testing.assert_array_equal(arrays["component_released_particles"], expected_counts)
    assert np.all(arrays["packet_counts"] > 0)
    for index, count in enumerate(expected_counts):
        selected = component == index
        law_record = source_components[index]["source_window"]["time_window"]
        law = ReducedTemkinEmissionWindow(**{name: law_record[name] for name in (
            "coupling", "log_initial_rate_s", "parent_duration_s", "start_s", "stop_s", "log_mass_tolerance")})
        np.testing.assert_allclose(delay[selected], law.quantile(uniforms[selected, 1]), rtol=3e-12, atol=0.)
        weights = arrays["particle_weight"][selected].astype(np.longdouble)
        assert np.count_nonzero(selected) == arrays["packet_counts"][index]
        np.testing.assert_allclose(np.sum(weights), count, rtol=4e-15, atol=0.)
        expected_energy = count * 2 * KB * float(temperatures[selected][0])
        sampled_energy = np.sum(weights*energy[selected], dtype=np.longdouble)
        np.testing.assert_allclose(arrays["component_expected_translational_energy_j"][index], expected_energy,
                                   rtol=4e-14, atol=0.)
        np.testing.assert_allclose(arrays["component_sampled_translational_energy_j"][index], sampled_energy,
                                   rtol=4e-14, atol=0.)
        r = history.reservoirs[int(pool[selected][0])]
        np.testing.assert_allclose(arrays["component_expected_binding_energy_j"][index],
                                   count*r.desorption_energy_ev*QE, rtol=4e-14, atol=0.)
        np.testing.assert_allclose(arrays["component_expected_internal_energy_j"][index],
                                   count*r.internal_energy_ev*QE, rtol=4e-14, atol=0.)
    expected = arrays["component_expected_translational_energy_j"]
    sampled = arrays["component_sampled_translational_energy_j"]
    residual = arrays["component_translational_energy_sampling_residual_j"]
    assert np.any(sampled != expected)
    assert np.all(np.abs(residual-(sampled-expected)) <= 1e-13*np.maximum(expected,sampled))


def test_explicit_seed_reproduces_arrays_without_assuming_archive_timestamps(source_fixture, exported, tmp_path):
    process, output, _ = run_cli(source_fixture, tmp_path)
    assert process.returncode == 0, process.stderr
    _, repeated, uniforms = read_export(output)
    _, original, original_uniforms = exported["data"]
    np.testing.assert_array_equal(uniforms, original_uniforms)
    assert set(repeated) == set(original)
    for name in repeated:
        np.testing.assert_array_equal(repeated[name], original[name])


def test_missing_source_summary_requires_and_retains_unvalidated_parameter_status(source_fixture, tmp_path):
    config = configuration()
    config["parameter_status"] = "unvalidated"
    process, output, _ = run_cli(source_fixture, tmp_path, config=config, summary=None)
    assert process.returncode == 0, process.stderr
    report, _, _ = read_export(output)
    assert report["parameter_status"] == "unvalidated"
    assert report["mass_status"] == "manufactured"
    assert report["source_run_association"] is None and report["input_summary_sha256"] is None
    assert report["conditioned_source"]["status"] == "unvalidated"
    assert not (output / "source_run_summary.json").exists()


def test_manufactured_parameters_cannot_be_asserted_without_the_source_summary(source_fixture, tmp_path):
    process, output, _ = run_cli(source_fixture, tmp_path, summary=None)
    assert_failed_without_stage(process, output)
    assert "unvalidated" in process.stderr


@pytest.mark.parametrize("kind", ["hash", "status", "provenance"])
def test_source_summary_must_bind_the_exact_hash_status_and_provenance(source_fixture, tmp_path, kind):
    summary = json.loads(source_fixture["summary"].read_text())
    if kind == "hash":
        summary["output_sha256"]["source.h5"] = "0" * 64
    elif kind == "status":
        summary["parameter_status"] = "sourced"
    else:
        summary["parameter_provenance"] = "Detached from the actual original source"
    path = write_json(tmp_path / "summary.json", summary)
    process, output, _ = run_cli(source_fixture, tmp_path, summary=path)
    assert_failed_without_stage(process, output)
    assert "hash-associate" in process.stderr


@pytest.mark.parametrize("kind", ["extra", "bool_count", "unknown_mass_status", "nonfinite_tolerance",
                                 "bool_start", "selection_extra", "bool_mass", "missing_exact_state",
                                 "species_fallback", "zero_normal", "too_few_packets"])
def test_invalid_configuration_or_allocation_is_rejected_before_export(source_fixture, tmp_path, kind):
    config = configuration()
    if kind == "extra": config["extra"] = 1
    elif kind == "bool_count": config["packet_count"] = True
    elif kind == "unknown_mass_status": config["mass_status"] = "unknown"
    elif kind == "nonfinite_tolerance": config["log_mass_tolerance"] = float("inf")
    elif kind == "bool_start": config["selections"][0]["start_s"] = True
    elif kind == "selection_extra": config["selections"][0]["extra"] = 1
    elif kind == "bool_mass": config["state_masses_kg"]["M:a"] = True
    elif kind == "missing_exact_state": del config["state_masses_kg"]["M:b"]
    elif kind == "species_fallback": config["state_masses_kg"] = {"M": 3e-26}
    elif kind == "zero_normal": config["outward_normal_sign"] = 0
    else: config["packet_count"] = 1
    process, output, _ = run_cli(source_fixture, tmp_path, config=config)
    assert_failed_without_stage(process, output)


@pytest.mark.parametrize("nested", [False, True])
def test_duplicate_json_keys_are_rejected_at_every_config_level(source_fixture, tmp_path, nested):
    text = json.dumps(configuration())
    if nested:
        text = text.replace('"state_masses_kg": {', '"state_masses_kg": {"M:a": 9e-26,', 1)
    else:
        text = text[:-1] + ', "seed": 19}'
    process, output, _ = run_cli(source_fixture, tmp_path, raw=text)
    assert_failed_without_stage(process, output)
    assert "duplicate JSON key" in process.stderr


@pytest.mark.parametrize("kind", ["directory", "file", "broken_symlink"])
def test_existing_output_paths_are_preserved_including_broken_symlinks(source_fixture, tmp_path, kind):
    target = tmp_path / "occupied"
    if kind == "directory":
        target.mkdir()
        (target / "sentinel").write_text("keep")
    elif kind == "file": target.write_text("keep")
    else: target.symlink_to(tmp_path / "missing-target")
    process, _, _ = run_cli(source_fixture, tmp_path, output=target)
    assert process.returncode != 0 and "new path" in process.stderr
    if kind == "directory": assert (target / "sentinel").read_text() == "keep"
    elif kind == "file": assert target.read_text() == "keep"
    else: assert target.is_symlink() and not target.exists()
    assert not list(tmp_path.glob(".temkin-windowed-*"))


def inprocess_inputs(source_fixture, tmp_path):
    source, summary = tmp_path / "source.h5", tmp_path / "summary.json"
    shutil.copyfile(source_fixture["source"], source)
    shutil.copyfile(source_fixture["summary"], summary)
    config = write_json(tmp_path / "config.json", configuration())
    return source, config, summary, tmp_path / "output"


def test_staged_array_write_failure_leaves_no_partial_export(source_fixture, tmp_path, monkeypatch, capsys):
    cli = importlib.import_module("tools.sample_temkin_source_windows")
    source, config, summary, output = inprocess_inputs(source_fixture, tmp_path)
    original = cli.np.save

    def fail_uniforms(path, *args, **kwargs):
        if Path(path).name == "uniforms.npy":
            raise OSError("manufactured staged write failure")
        return original(path, *args, **kwargs)

    monkeypatch.setattr(cli.np, "save", fail_uniforms)
    monkeypatch.setattr(sys, "argv", command(source, config, output, summary))
    with pytest.raises(SystemExit) as error:
        cli.main()
    assert error.value.code == 2
    assert "manufactured staged write failure" in capsys.readouterr().err
    assert not output.exists() and not list(tmp_path.glob(".temkin-windowed-*"))


@pytest.mark.parametrize("kind", ["source", "config", "summary", "code"])
def test_input_and_implementation_changes_during_staging_are_rejected(source_fixture, tmp_path, monkeypatch, capsys, kind):
    cli = importlib.import_module("tools.sample_temkin_source_windows")
    source, config, summary, output = inprocess_inputs(source_fixture, tmp_path)
    original_save, original_sha = cli.np.savez_compressed, cli.sha
    staged = {"written": False}

    def save_then_change(*args, **kwargs):
        original_save(*args, **kwargs)
        staged["written"] = True
        if kind == "source":
            with h5py.File(source, "a") as handle:
                handle.attrs["manufactured_race_marker"] = "changed after source copy"
        elif kind == "config": config.write_bytes(config.read_bytes() + b"\n")
        elif kind == "summary": summary.write_bytes(summary.read_bytes() + b"\n")

    def guarded_sha(path):
        if kind == "code" and staged["written"] and Path(path).resolve() == SCRIPT:
            return "0" * 64
        return original_sha(path)

    monkeypatch.setattr(cli.np, "savez_compressed", save_then_change)
    monkeypatch.setattr(cli, "sha", guarded_sha)
    monkeypatch.setattr(sys, "argv", command(source, config, output, summary))
    with pytest.raises(SystemExit) as error:
        cli.main()
    assert error.value.code == 2 and staged["written"]
    assert "changed during export" in capsys.readouterr().err
    assert not output.exists() and not list(tmp_path.glob(".temkin-windowed-*"))


def test_output_created_during_staging_is_not_replaced(source_fixture, tmp_path, monkeypatch, capsys):
    cli = importlib.import_module("tools.sample_temkin_source_windows")
    source, config, summary, output = inprocess_inputs(source_fixture, tmp_path)
    original_save = cli.np.save

    def save_then_competing_writer(path, *args, **kwargs):
        result = original_save(path, *args, **kwargs)
        if Path(path).name == "uniforms.npy":
            output.mkdir()
            (output / "other-writer").write_text("keep")
        return result

    monkeypatch.setattr(cli.np, "save", save_then_competing_writer)
    monkeypatch.setattr(sys, "argv", command(source, config, output, summary))
    with pytest.raises(SystemExit) as error:
        cli.main()
    assert error.value.code == 2
    assert "output path appeared" in capsys.readouterr().err
    assert (output / "other-writer").read_text() == "keep"
    assert not list(tmp_path.glob(".temkin-windowed-*"))
