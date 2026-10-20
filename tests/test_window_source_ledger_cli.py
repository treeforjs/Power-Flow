"""Fresh saved-study export checks with independent source/energy controls."""
from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
from fractions import Fraction
import copy
import csv
import hashlib
import importlib
import io
import json
import math
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from mhdlab.constants import KB, QE
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_planar_io import write_temkin_planar_release_history
from mhdlab.temkin_release import TemkinReservoir


ROOT = Path(__file__).resolve().parents[1]
PROVENANCE = "Fresh manufactured unequal-area two-pool ledger export fixture"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def read_json(path):
    return json.loads(path.read_text())


def invoke(command, study, output):
    stdout, stderr = io.StringIO(), io.StringIO()
    with patch.object(sys, "argv", [str(command.__file__), "--study", str(study), "--output-dir", str(output)]):
        with redirect_stdout(stdout), redirect_stderr(stderr):
            try:
                command.main()
                code = 0
            except SystemExit as error:
                code = error.code
    return SimpleNamespace(returncode=code, stdout=stdout.getvalue(), stderr=stderr.getvalue())


@pytest.fixture(scope="module")
def fresh_saved_study(tmp_path_factory):
    """Create new source/moments once; export tests never rerun those paths."""
    directory = tmp_path_factory.mktemp("fresh-window-ledger-study")
    times = [0., 1e-7]
    x = [-.003, -.001, .002]
    z = [.001, .002, .006]
    temperatures = np.array([[350., 700.], [500., 900.]])
    initial = np.array([[[1e17, 2e17], [3e17, 4e17]],
                        [[2e17, 4e17], [6e17, 8e17]]])
    pools = (
        TemkinReservoir("sites-a", "M:state0", 1e7, 0., .11, .02, 0., 1e18, 1e18,
                        "Manufactured first-order pool A", "Manufactured independent reference A"),
        TemkinReservoir("sites-b", "N:state1", 2e7, 0., .23, .07, 0., 1e18, 1e18,
                        "Manufactured first-order pool B", "Manufactured independent reference B"),
    )
    definition = directory/"definition.json"
    write_json(definition, {"parameter_status": "manufactured", "x_edges_m": x, "z_edges_m": z,
                            "initial_m2": initial.tolist(), "temperature_k": temperatures.tolist(),
                            "rates_s": [pool.attempt_frequency_s for pool in pools]})
    base = FlashSurfaceHistory(times, [0., 1.], [[400.], [400.]], "synthetic", "fresh-ledger-export-test",
                               "Manufactured stationary patch temperatures", "Manufactured boundary")
    surface = PlanarSurfaceHistory(times, x, z, np.stack((temperatures, temperatures)), base.attributes)
    history = evolve_temkin_planar_release(surface, initial, pools, parameter_provenance=PROVENANCE)
    source = directory/"fresh.h5"
    write_temkin_planar_release_history(source, history)
    resolution = {"width_bins": 2, "current_bins": 2, "snapshot_intervals": 1, "subdivisions_per_interval": 1}
    evidence = directory/"association.json"
    write_json(evidence, {"schema": "mhdlab-preserved-source-association-v1",
                          "parameter_status": "manufactured", "parameter_provenance": PROVENANCE,
                          "source_sha256": sha(source), "resolution": resolution,
                          "preserved_inputs_sha256": {definition.name: sha(definition)},
                          "source_evolution_repeated": False})
    config = directory/"config.json"
    windows = [{"start_s": 1e-8, "stop_s": 4e-8}, {"start_s": 6e-8, "stop_s": 9e-8}]
    write_json(config, {
        "schema": "mhdlab-windowed-source-resolution-config-v1", "parameter_status": "manufactured",
        "conditioning_provenance": "Fresh manufactured source-emission restriction; no photon gate",
        "comparison_provenance": "Fresh saved single case, with no convergence comparison",
        "relative_time_tolerance": 1e-12, "log_mass_tolerance": 1e-7, "windows": windows,
        "orders": [4, 8], "cases": [{"case_id": "fresh", "source_file": str(source),
            "source_sha256": sha(source), "evidence_file": str(evidence), "evidence_sha256": sha(evidence),
            "resolution": resolution}], "comparisons": [],
    })
    comparison = importlib.import_module("tools.compare_windowed_source_resolution")
    output = directory/"saved-study"
    with patch.object(sys, "argv", [str(comparison.__file__), "--config", str(config), "--output-dir", str(output)]):
        comparison.main()
    return SimpleNamespace(path=output, initial=initial, temperatures=temperatures, x=np.array(x), z=np.array(z),
                           pools=pools, windows=windows)


@pytest.fixture(autouse=True)
def forbid_old_execution(fresh_saved_study, monkeypatch):
    """Canonical constructor validation remains allowed; drivers do not."""
    def forbidden(*args, **kwargs):
        pytest.fail("saved-study export invoked an old evolution, moments, or sampling path")
    targets = {
        "mhdlab.temkin_planar": ["evolve_temkin_planar_release"],
        "mhdlab.temkin_history": ["evolve_temkin_release"],
        "mhdlab.windowed_joint_source_moments": ["windowed_joint_source_moments"],
        "tools.compare_windowed_source_resolution": ["inspect_case", "windowed_joint_source_moments"],
        "mhdlab.thermal_boundary": ["sample_emission_packets"],
        "mhdlab.planar_boundary": ["sample_planar_emission_packets"],
        "mhdlab.temkin_boundary": ["sample_temkin_emission_packets"],
        "mhdlab.temkin_planar_boundary": ["sample_temkin_planar_emission_packets"],
        "mhdlab.temkin_windowed_boundary": ["sample_temkin_windowed_emission_packets"],
    }
    for module, names in targets.items():
        loaded = importlib.import_module(module)
        for name in names:
            monkeypatch.setattr(loaded, name, forbidden)
    monkeypatch.setattr(np.random, "default_rng", forbidden)


@pytest.fixture
def command():
    return importlib.import_module("tools.export_window_source_ledger")


@pytest.fixture
def study(fresh_saved_study, tmp_path):
    destination = tmp_path/"study"
    shutil.copytree(fresh_saved_study.path, destination)
    return destination


def assert_rejected(process, output):
    assert process.returncode != 0, process.stdout
    assert not output.exists() and not output.is_symlink()
    assert not [p for p in output.parent.iterdir() if p.name.startswith(".window-source-ledger-")]


def rebound_json(study, filename, value):
    """Rebind outer bytes, so semantic validation rather than SHA alone is tested."""
    write_json(study/filename, value)
    report = read_json(study/"report.json")
    report["output_sha256"][filename] = sha(study/filename)
    if filename == "input_config.json":
        report["input_config_sha256"] = sha(study/filename)
    write_json(study/"report.json", report)


def quantity(group, name):
    record = group["quantities"][name]
    assert record["available"] is True
    assert record["projection"]["finite"] is True
    assert record["projection"]["underflowed"] is False
    assert record["projection"]["overflowed"] is False
    return record["value"]


def analytic_window(fixture, window):
    duration = window["stop_s"]-window["start_s"]
    areas = np.diff(fixture.z)[:, None]*np.diff(fixture.x)[None, :]
    counts = np.array([
        initial*areas*math.exp(-pool.attempt_frequency_s*window["start_s"])
        * (-math.expm1(-pool.attempt_frequency_s*duration))
        for initial, pool in zip(fixture.initial, fixture.pools)
    ])
    return {
        "released_number": counts,
        "translational_energy_j": counts*2*KB*fixture.temperatures[None],
        "internal_energy_j": counts*np.array([pool.internal_energy_ev for pool in fixture.pools])[:, None, None]*QE,
        "binding_energy_j": counts*np.array([pool.desorption_energy_ev for pool in fixture.pools])[:, None, None]*QE,
    }


def test_new_ledgers_match_independent_counts_and_all_energy_paths(command, study, tmp_path, fresh_saved_study):
    output = tmp_path/"ledger"
    process = invoke(command, study, output)
    assert process.returncode == 0, process.stderr
    result = read_json(output/"case_fresh.json")
    assert result["case_id"] == "fresh" and result["parameter_status"] == "manufactured"
    assert result["parameter_provenance"] == PROVENANCE
    assert len(result["windows"]) == 2
    for row, window in zip(result["windows"], fresh_saved_study.windows):
        assert row["window"] == window
        ledger = row["ledger"]
        assert ledger["schema"] == "mhdlab-window-source-ledger-v1"
        expected = analytic_window(fresh_saved_study, window)
        expected["total_wall_energy_debit_j"] = (expected["translational_energy_j"]
            + expected["internal_energy_j"] + expected["binding_energy_j"])
        for hypothesis in ("resolved", "matched_uniform_width"):
            total = ledger[hypothesis]["total"]
            for name, values in expected.items():
                assert quantity(total, name) == pytest.approx(float(np.sum(values)), rel=1e-12, abs=0.)
            assert quantity(total, "total_wall_energy_debit_j") == pytest.approx(
                math.fsum(quantity(total, name) for name in (
                    "translational_energy_j", "internal_energy_j", "binding_energy_j")), rel=1e-14, abs=0.)
            for group in ledger[hypothesis]["by_state"]:
                index = [pool.state_id for pool in fresh_saved_study.pools].index(group["key"]["state_id"])
                for name, values in expected.items():
                    assert quantity(group, name) == pytest.approx(float(np.sum(values[index])), rel=1e-12, abs=0.)
        for name in expected:
            assert ledger["resolved"]["total"]["quantities"][name]["exact_sum"] == (
                ledger["matched_uniform_width"]["total"]["quantities"][name]["exact_sum"])


def test_unequal_area_width_control_preserves_full_face_mixture(command, study, tmp_path, fresh_saved_study):
    output = tmp_path/"ledger"
    process = invoke(command, study, output)
    assert process.returncode == 0, process.stderr
    result = read_json(output/"case_fresh.json")
    fractions = np.diff(fresh_saved_study.x)/(fresh_saved_study.x[-1]-fresh_saved_study.x[0])
    face_area = (fresh_saved_study.x[-1]-fresh_saved_study.x[0])*(fresh_saved_study.z[-1]-fresh_saved_study.z[0])
    for entry, window in zip(result["windows"], fresh_saved_study.windows):
        ledger = entry["ledger"]
        expected = analytic_window(fresh_saved_study, window)
        resolved = {group["key"]["width_index"]: group for group in ledger["resolved"]["by_width"]}
        uniform = {group["key"]["width_index"]: group for group in ledger["matched_uniform_width"]["by_width"]}
        for index in range(2):
            for name, values in expected.items():
                assert quantity(resolved[index], name) == pytest.approx(float(np.sum(values[:, :, index])), rel=1e-12, abs=0.)
                assert quantity(uniform[index], name) == pytest.approx(float(np.sum(values))*fractions[index], rel=1e-12, abs=0.)
        assert quantity(resolved[0], "released_number") != pytest.approx(quantity(uniform[0], "released_number"), rel=1e-2, abs=0.)
        total = ledger["resolved"]["total"]
        assert total["area"]["value_m2"] == pytest.approx(face_area, rel=1e-14, abs=0.)
        for name in expected:
            per_area = total["quantities"][name]["per_area"]
            assert per_area["available"] is True
            assert per_area["value"] == pytest.approx(quantity(total, name)/face_area, rel=1e-14, abs=0.)


def test_copies_old_study_bytes_and_binds_new_files(command, study, tmp_path):
    output = tmp_path/"ledger"
    original = {path.name: path.read_bytes() for path in study.iterdir() if path.is_file()}
    process = invoke(command, study, output)
    assert process.returncode == 0, process.stderr
    report = read_json(output/"report.json")
    assert report["schema"] == "mhdlab-window-source-ledger-export-v1"
    assert report["source_evolution_repeated"] is False
    assert report["host_budget_accepted"] is False
    for name, content in original.items():
        assert (output/"study"/name).read_bytes() == content
        assert (study/name).read_bytes() == content
    for name, expected in report["output_sha256"].items():
        assert sha(output/name) == expected
    assert (output/"ledger.csv").is_file()


def test_flat_csv_rows_match_ledger_groups_and_use_numeric_flags(command, study, tmp_path):
    output = tmp_path/"ledger"
    process = invoke(command, study, output)
    assert process.returncode == 0, process.stderr
    result = read_json(output/"case_fresh.json")
    with (output/"ledger.csv").open(newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == list(command.CSV_FIELDS)
        rows = list(reader)
    expected = {}
    for window in result["windows"]:
        for hypothesis in ("resolved", "matched_uniform_width"):
            for group_type, collection in window["ledger"][hypothesis].items():
                groups = [collection] if group_type == "total" else collection
                for index, group in enumerate(groups):
                    for name, amount in group["quantities"].items():
                        expected[window["window_index"], hypothesis, group_type, index, name] = (group, amount)
    assert len(rows) == len(expected)
    for row in rows:
        key = (int(row["window_index"]), row["hypothesis"], row["group_type"], int(row["group_index"]), row["quantity"])
        group, amount = expected.pop(key)
        assert row["case_id"] == "fresh"
        assert row["available"] == "1"
        assert float(row["value"]) == amount["value"]
        assert float(row["per_area_value"]) == amount["per_area"]["value"]
        assert float(row["area_m2"]) == group["area"]["value_m2"]
        for prefix in ("", "area_", "qualified_subtotal_", "per_area_", "per_area_qualified_subtotal_"):
            for flag in command.PROJECTION_FLAGS:
                assert row[prefix+flag] in {"0", "1"}
        for field, value in row.items():
            assert not value.startswith(("{", "[")), f"structured record leaked into CSV field{field}"
    assert not expected


def test_csv_availability_is_separate_from_positive_underflow_and_overflow(command, tmp_path):
    def projected(value, *, positive=False, zero=False, underflow=False, overflow=False):
        return {"value": value, "exact_sum": {"numerator": "1", "denominator": "1"},
                "projection": {"exact_zero": zero, "positive": positive, "negative": False,
                    "finite": not overflow, "underflowed": underflow, "overflowed": overflow,
                    "rounded": underflow or overflow, "subnormal": False,
                    "exact_rounding_error": None}}
    def amount(number, *, available=True):
        return {"available": available,
                **(copy.deepcopy(number) if available else {"value": None, "exact_sum": None, "projection": None}),
                "qualified_subtotal": projected(0., zero=True),
                "per_area": {"available": available,
                    **(copy.deepcopy(number) if available else {"value": None, "exact_sum": None, "projection": None}),
                    "qualified_subtotal": projected(0., zero=True)},
                "unqualified": {"count": 0 if available else 1, "catalog_indices": [] if available else [0]}}
    quantities = {
        "released_number": amount(projected(0., positive=True, underflow=True)),
        "translational_energy_j": amount(projected(None, positive=True, overflow=True)),
        "internal_energy_j": amount(projected(None), available=False),
        "binding_energy_j": amount(projected(0., zero=True)),
        "total_wall_energy_debit_j": amount(projected(1e-18, positive=True)),
    }
    group = {"key": {}, "area": {"value_m2": 1., "projection": projected(1., positive=True)["projection"]},
             "quantities": quantities}
    ledger = {"resolved": {"total": group}, "matched_uniform_width": {"total": copy.deepcopy(group)}}
    output = tmp_path/"serialization.csv"
    command.write_ledger_csv(output, [{"case_id": "serialization-control", "window_index": 0,
                                      "window": None, "ledger": ledger}])
    with output.open(newline="") as handle:
        rows = {row["quantity"]: row for row in csv.DictReader(handle) if row["hypothesis"] == "resolved"}
    tiny = rows["released_number"]
    assert (tiny["available"], tiny["value"], tiny["positive"], tiny["underflowed"], tiny["finite"]) == (
        "1", "0.0", "1", "1", "1")
    assert tiny["projection_reason"] == "si_projection_underflow"
    huge = rows["translational_energy_j"]
    assert (huge["available"], huge["value"], huge["positive"], huge["overflowed"], huge["finite"]) == (
        "1", "", "1", "1", "0")
    assert huge["projection_reason"] == "si_projection_overflow"
    missing = rows["internal_energy_j"]
    assert missing["available"] == "0" and missing["value"] == "" and missing["positive"] == ""
    assert missing["unqualified_count"] == "1" and missing["qualified_subtotal_value"] == "0.0"
    assert missing["availability_reason"] == "unqualified_contributors"
    zero = rows["binding_energy_j"]
    assert (zero["available"], zero["value"], zero["exact_zero"], zero["positive"], zero["underflowed"]) == (
        "1", "0.0", "1", "0", "0")


@pytest.mark.parametrize("mutation", ["source_hash", "code_pin", "status", "conditioned_status", "state_id",
                                      "missing_positive_component", "plan", "provenance"])
def test_saved_associations_are_validated_after_outer_rehash(command, study, tmp_path, mutation):
    report = read_json(study/"report.json")
    if mutation == "source_hash":
        with (study/"source_fresh.h5").open("ab") as handle:
            handle.write(b"changed-after-recording")
    elif mutation == "code_pin":
        name = next(iter(report["source_sha256"]))
        report["source_sha256"][name] = "0"*64
        write_json(study/"report.json", report)
    elif mutation == "status":
        config = read_json(study/"input_config.json")
        config["parameter_status"] = "sourced"
        report["parameter_status"] = "sourced"
        write_json(study/"report.json", report)
        rebound_json(study, "input_config.json", config)
    elif mutation in {"conditioned_status", "state_id", "missing_positive_component", "plan"}:
        case = read_json(study/"case_fresh.json")
        if mutation == "conditioned_status":
            case["windows"][0]["conditioned_source"]["status"] = "sourced"
        elif mutation == "state_id":
            case["windows"][0]["conditioned_source"]["source_reservoirs"][0]["state_id"] = "different:state"
        elif mutation == "missing_positive_component":
            case["windows"][0]["conditioned_source"]["components"].pop()
        else:
            case["plan"]["window_indices"][0] = 9
        rebound_json(study, "case_fresh.json", case)
    else:
        evidence = read_json(study/"evidence_fresh.json")
        evidence["parameter_provenance"] = "Different physical source association"
        write_json(study/"evidence_fresh.json", evidence)
        config = read_json(study/"input_config.json")
        config["cases"][0]["evidence_sha256"] = sha(study/"evidence_fresh.json")
        report["cases"][0]["evidence_sha256"] = config["cases"][0]["evidence_sha256"]
        report["output_sha256"]["evidence_fresh.json"] = sha(study/"evidence_fresh.json")
        write_json(study/"report.json", report)
        rebound_json(study, "input_config.json", config)
    output = tmp_path/"ledger"
    assert_rejected(invoke(command, study, output), output)


@pytest.mark.parametrize("mutation", ["plan_boolean_index", "pool_boolean_coefficient", "evidence_boolean_count"])
def test_boolean_numeric_association_aliases_reject(command, study, tmp_path, mutation):
    """JSON booleans must not alias equal-valued numeric canonical metadata."""
    if mutation in {"plan_boolean_index", "pool_boolean_coefficient"}:
        case = read_json(study/"case_fresh.json")
        if mutation == "plan_boolean_index":
            assert type(case["plan"]["window_indices"][0]) is int
            assert case["plan"]["window_indices"][0] == 0
            case["plan"]["window_indices"][0] = False
        else:
            pool = case["windows"][0]["conditioned_source"]["source_reservoirs"][0]
            assert type(pool["coverage_coefficient"]) is float and pool["coverage_coefficient"] == 0.
            pool["coverage_coefficient"] = False
        rebound_json(study, "case_fresh.json", case)
    else:
        evidence = read_json(study/"evidence_fresh.json")
        assert type(evidence["resolution"]["snapshot_intervals"]) is int
        assert evidence["resolution"]["snapshot_intervals"] == 1
        evidence["resolution"]["snapshot_intervals"] = True
        write_json(study/"evidence_fresh.json", evidence)
        config = read_json(study/"input_config.json")
        config["cases"][0]["evidence_sha256"] = sha(study/"evidence_fresh.json")
        report = read_json(study/"report.json")
        report["cases"][0]["evidence_sha256"] = config["cases"][0]["evidence_sha256"]
        report["output_sha256"]["evidence_fresh.json"] = sha(study/"evidence_fresh.json")
        write_json(study/"report.json", report)
        rebound_json(study, "input_config.json", config)
    output = tmp_path/"ledger"
    assert_rejected(invoke(command, study, output), output)


@pytest.mark.parametrize("mutation", ["duplicate_json", "nonfinite_json", "unsafe_output_key", "unsafe_case_id", "input_symlink"])
def test_unsafe_or_ambiguous_saved_inputs_reject(command, study, tmp_path, mutation):
    if mutation == "duplicate_json":
        path = study/"report.json"
        raw = path.read_text()
        path.write_text('{"schema":"duplicate",'+raw[1:])
    elif mutation == "nonfinite_json":
        path = study/"report.json"
        raw = path.read_text()
        path.write_text('{"unvalidated_float":NaN,'+raw[1:])
    elif mutation == "unsafe_output_key":
        report = read_json(study/"report.json")
        report["output_sha256"]["../escape.json"] = "1"*64
        write_json(study/"report.json", report)
    elif mutation == "unsafe_case_id":
        config = read_json(study/"input_config.json")
        config["cases"][0]["case_id"] = "../escape"
        rebound_json(study, "input_config.json", config)
    else:
        source = study/"source_fresh.h5"
        replacement = tmp_path/"same-source.h5"
        shutil.move(source, replacement)
        source.symlink_to(replacement)
    output = tmp_path/"ledger"
    assert_rejected(invoke(command, study, output), output)


def test_existing_output_is_preserved(command, study, tmp_path):
    output = tmp_path/"existing"
    output.mkdir()
    marker = output/"keep.txt"
    marker.write_text("preserve these bytes")
    process = invoke(command, study, output)
    assert process.returncode != 0
    assert marker.read_text() == "preserve these bytes"
    assert set(output.iterdir()) == {marker}


def test_output_symlink_is_refused(command, study, tmp_path):
    target = tmp_path/"target"
    target.mkdir()
    output = tmp_path/"ledger-link"
    output.symlink_to(target, target_is_directory=True)
    assert invoke(command, study, output).returncode != 0
    assert output.is_symlink() and not list(target.iterdir())


def test_late_input_mutation_rolls_back_new_destination(command, study, tmp_path, monkeypatch):
    original_copy = command.shutil.copyfile
    changed = []
    def change_after_copy(source, destination, *args, **kwargs):
        result = original_copy(source, destination, *args, **kwargs)
        source = Path(source)
        if source.parent == study and not changed:
            with source.open("ab") as handle:
                handle.write(b" ")
            changed.append(source)
        return result
    monkeypatch.setattr(command.shutil, "copyfile", change_after_copy)
    output = tmp_path/"ledger"
    process = invoke(command, study, output)
    assert changed, "late mutation hook was not exercised"
    assert_rejected(process, output)
