"""Externally pinned, purely algebraic admission of fresh source-ledger exports."""
from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import copy
import csv
import hashlib
import importlib
import io
import json
import os
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_planar_io import write_temkin_planar_release_history
from mhdlab.temkin_release import TemkinReservoir


ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def driver(module_name, arguments):
    module = importlib.import_module(module_name)
    with patch.object(sys, "argv", [str(module.__file__), *map(str, arguments)]):
        module.main()


@pytest.fixture(scope="module")
def fresh_exports(tmp_path_factory):
    """Generate new examples once; admission itself must consume only JSON/CSV."""
    preserved_fixture = os.environ.get("MITL_ADMISSION_CLI_FIXTURE_ROOT")
    if preserved_fixture is not None:
        directory = Path(preserved_fixture)
        if directory.exists():
            index = read_json(directory/"fixture_index.json")
            assert index["schema"] == "mhdlab-new-admission-cli-fixture-v1"
            result = {name: directory/relative for name, relative in index["exports"].items()}
            assert set(result) == {"study", "packet"}
            assert all(sha(path/"report.json") == index["report_sha256"][name] for name, path in result.items())
            return result
        directory.mkdir(parents=True)
    else:
        directory = tmp_path_factory.mktemp("fresh-ledger-admission-inputs")
    times = [0., 1e-7]
    base = FlashSurfaceHistory(times, [0., 1.], [[400.], [400.]], "synthetic", "fresh-admission-fixture",
                               "Manufactured constant patch temperatures", "Manufactured stationary face")
    surface = PlanarSurfaceHistory(times, [-.003, -.001, .002], [0., .001, .004],
                                  [[[400., 700.], [500., 900.]]]*2, base.attributes)
    pools = [
        TemkinReservoir("pool-a", "M:a", 1e7, 0., .12, .03, 0., 1e18, 1e18,
                        "Manufactured first-order source A", "Manufactured supplied energy A"),
        TemkinReservoir("pool-b", "M:b", 2e7, 0., .24, .07, 0., 1e18, 1e18,
                        "Manufactured first-order source B", "Manufactured supplied energy B"),
    ]
    history = evolve_temkin_planar_release(surface,
        np.array([[[1e17, 2e17], [2e17, 4e17]], [[0., 2e17], [3e17, 4e17]]]), pools,
        parameter_provenance="Fresh manufactured source for pure export admission")
    source = directory/"source.h5"
    write_temkin_planar_release_history(source, history)
    definition = directory/"definition.json"
    write_json(definition, {"parameter_status": "manufactured", "provenance": history.parameter_provenance})
    summary = directory/"summary.json"
    write_json(summary, {"schema": "mhdlab-temkin-planar-run-v1", "parameter_status": "manufactured",
                        "parameter_provenance": history.parameter_provenance,
                        "output_sha256": {source.name: sha(source)}})
    resolution = {"width_bins": 2, "current_bins": 2, "snapshot_intervals": 1, "subdivisions_per_interval": 1}
    evidence = directory/"evidence.json"
    write_json(evidence, {"schema": "mhdlab-preserved-source-association-v1", "parameter_status": "manufactured",
                         "parameter_provenance": history.parameter_provenance, "source_sha256": sha(source),
                         "resolution": resolution, "preserved_inputs_sha256": {definition.name: sha(definition)},
                         "source_evolution_repeated": False})
    windows = [{"start_s": 1e-8, "stop_s": 3e-8}, {"start_s": 6e-8, "stop_s": 9e-8}]
    study_config = directory/"study_config.json"
    write_json(study_config, {
        "schema": "mhdlab-windowed-source-resolution-config-v1", "parameter_status": "manufactured",
        "conditioning_provenance": "Manufactured source windows for pure admission",
        "comparison_provenance": "Manufactured single preserved source without convergence inference",
        "relative_time_tolerance": 1e-12, "log_mass_tolerance": 1e-7, "windows": windows, "orders": [4, 8],
        "cases": [{"case_id": "fresh", "source_file": str(source), "source_sha256": sha(source),
                   "evidence_file": str(evidence), "evidence_sha256": sha(evidence), "resolution": resolution}],
        "comparisons": [],
    })
    study = directory/"study"
    driver("tools.compare_windowed_source_resolution", ["--config", study_config, "--output-dir", study])
    study_export = directory/"study-ledger"
    driver("tools.export_window_source_ledger", ["--study", study, "--output-dir", study_export])
    packet_config = directory/"packet_config.json"
    write_json(packet_config, {
        "schema": "mhdlab-temkin-windowed-packet-config-v1", "parameter_status": "manufactured",
        "conditioning_provenance": "Manufactured composite source selection for pure admission",
        "mass_status": "manufactured", "mass_provenance": "Manufactured two-parent mass declarations",
        "state_masses_kg": {"M:a": 3e-26, "M:b": 6e-26},
        "selections": [{"interval_index": 0, "start_s": window["start_s"], "stop_s": window["stop_s"]}
                       for window in windows],
        "log_mass_tolerance": 1e-7, "face_normal_coordinate_m": 0., "outward_normal_sign": 1,
        "packet_count": 32, "seed": 19341,
    })
    packet = directory/"packet"
    driver("tools.sample_temkin_source_windows", ["--source", source, "--source-run-summary", summary,
        "--config", packet_config, "--output-dir", packet])
    packet_export = directory/"packet-ledger"
    driver("tools.export_packet_source_ledger", ["--packet-export", packet, "--output-dir", packet_export])
    result = {"study": study_export, "packet": packet_export}
    if preserved_fixture is not None:
        write_json(directory/"fixture_index.json", {
            "schema": "mhdlab-new-admission-cli-fixture-v1",
            "exports": {name: path.relative_to(directory).as_posix() for name, path in result.items()},
            "report_sha256": {name: sha(path/"report.json") for name, path in result.items()},
            "interpretation": "Fresh fixture generation completed once; later focused invocations read these exact exports only",
        })
    return result


@pytest.fixture(autouse=True)
def pure_admission_only(fresh_exports, monkeypatch):
    """A new admission must not revalidate inputs through physical machinery."""
    def forbidden(*args, **kwargs):
        pytest.fail("admission invoked a physical constructor, HDF/array loader, law, source, moment, sampler, or ledger producer")
    targets = {
        "mhdlab.temkin_planar": ["evolve_temkin_planar_release"],
        "mhdlab.temkin_history": ["evolve_temkin_release"],
        "mhdlab.temkin_release": ["release_step", "energy_budgets"],
        "mhdlab.temkin_planar_io": ["load_temkin_planar_release_history"],
        "mhdlab.temkin_windowed_source": ["condition_temkin_planar_source"],
        "mhdlab.source_window_plan": ["plan_source_windows"],
        "mhdlab.window_source_ledger": ["window_source_ledger"],
        "mhdlab.windowed_joint_source_moments": ["windowed_joint_source_moments"],
        "mhdlab.temkin_windowed_boundary": ["sample_temkin_windowed_emission_packets"],
        "tools.compare_windowed_source_resolution": ["configuration", "load_case", "inspect_case", "windowed_joint_source_moments"],
        "tools.inspect_windowed_source_moments": ["associated_source", "windowed_joint_source_moments"],
        "tools.export_window_source_ledger": ["validate_study", "inspect_case", "window_source_ledger"],
        "tools.export_packet_source_ledger": ["associated_source", "window_source_ledger"],
    }
    for module_name, names in targets.items():
        module = importlib.import_module(module_name)
        for name in names:
            monkeypatch.setattr(module, name, forbidden)
    classes = {
        "mhdlab.temkin_release": ["TemkinReservoir"],
        "mhdlab.temkin_planar": ["TemkinPlanarReleaseHistory"],
        "mhdlab.temkin_history": ["TemkinReleaseHistory"],
        "mhdlab.temkin_source_window": ["FrozenTemkinSourceWindow"],
        "mhdlab.temkin_emission_window": ["ReducedTemkinEmissionWindow"],
        "mhdlab.temkin_windowed_source": ["SourceWindowSelection", "TemkinPlanarSourceWindows"],
        "mhdlab.source_window_plan": ["SourceTimeWindow", "SourceWindowPlan"],
    }
    for module_name, names in classes.items():
        module = importlib.import_module(module_name)
        for name in names:
            monkeypatch.setattr(getattr(module, name), "__post_init__", forbidden)
    monkeypatch.setattr(h5py, "File", forbidden)
    monkeypatch.setattr(np, "load", forbidden)
    monkeypatch.setattr(np.random, "default_rng", forbidden)


@pytest.fixture
def command():
    return importlib.import_module("tools.admit_source_ledger_export")


@pytest.fixture(params=("study", "packet"))
def exported(request, fresh_exports, tmp_path):
    path = tmp_path/"ledger-export"
    shutil.copytree(fresh_exports[request.param], path)
    return request.param, path


def invoke(command, export, output, *, expected=None, omit_pin=False):
    args = [str(command.__file__), "--export", str(export), "--output-dir", str(output)]
    if not omit_pin:
        args.extend(["--expected-report-sha256", sha(export/"report.json") if expected is None else expected])
    stdout, stderr = io.StringIO(), io.StringIO()
    with patch.object(sys, "argv", args), redirect_stdout(stdout), redirect_stderr(stderr):
        try:
            command.main()
            code = 0
        except SystemExit as error:
            code = error.code
    return SimpleNamespace(returncode=code, stdout=stdout.getvalue(), stderr=stderr.getvalue())


def assert_rejected(process, output):
    assert process.returncode != 0, process.stdout
    assert not output.exists() and not output.is_symlink()
    assert not list(output.parent.glob(".source-ledger-admission-*"))


def ledger_path(kind):
    return "case_fresh.json" if kind == "study" else "ledger.json"


def get_ledger(kind, document, index=0):
    return document["windows"][index]["ledger"] if kind == "study" else document


def mutate_ledger(export, kind, change):
    name = ledger_path(kind)
    document = read_json(export/name)
    change(get_ledger(kind, document))
    rebound_file(export, name, document)


def rebound_file(export, name, document):
    write_json(export/name, document)
    report = read_json(export/"report.json")
    report["output_sha256"][name] = sha(export/name)
    write_json(export/"report.json", report)


def refresh_csv_projection(export, kind):
    """Match a tampered ledger's scalar CSV without regenerating its physics.

    This pure test-setup projection prevents a CSV mismatch from being the only
    reason that an incorrect source identity or group association is refused.
    """
    document = read_json(export/ledger_path(kind))
    records = ([{"case_id": document["case_id"], **row} for row in document["windows"]]
               if kind == "study" else [{"case_id": "packet_source", "window_index": None,
                                          "window": None, "ledger": document}])
    importlib.import_module("tools.export_window_source_ledger").write_ledger_csv(export/"ledger.csv", records)
    report = read_json(export/"report.json")
    report["output_sha256"]["ledger.csv"] = sha(export/"ledger.csv")
    write_json(export/"report.json", report)


def pointer(document, reference):
    assert reference == "" or reference.startswith("/")
    for part in reference.split("/")[1:]:
        key = part.replace("~1", "/").replace("~0", "~")
        document = document[int(key)] if isinstance(document, list) else document[key]
    return document


def test_admits_both_export_schemas_without_physical_calls(command, exported, tmp_path):
    kind, export = exported
    output = tmp_path/"admitted"
    pin = sha(export/"report.json")
    process = invoke(command, export, output, expected=pin)
    assert process.returncode == 0, process.stderr
    admission = read_json(output/"admission.json")
    handoff = read_json(output/"owner_handoff.json")
    assert admission["schema"] == "mhdlab-source-ledger-export-admission-v1"
    assert handoff["schema"] == "mhdlab-source-ledger-owner-handoff-v1"
    assert admission["expected_report_sha256"] == pin
    assert admission["input_export_schema"] == read_json(export/"report.json")["schema"]
    for result in (admission, handoff):
        assert result["host_budget_accepted"] is False
        assert result["owner_implementation_verified"] is False
    sources = handoff["sources"]
    assert len(sources) == (2 if kind == "study" else 1)
    for source in sources:
        ledger_reference = source["ledger"]
        source_reference = source["original_conditioned_source"]
        assert sha(output/ledger_reference["path"]) == ledger_reference["sha256"]
        assert sha(output/source_reference["path"]) == source_reference["sha256"]
        ledger_document = read_json(output/ledger_reference["path"])
        source_document = read_json(output/source_reference["path"])
        ledger = pointer(ledger_document, ledger_reference["json_pointer"])
        original = pointer(source_document, source_reference["json_pointer"])
        assert ledger["components_included"] is True
        joint = source["component_joint_source"]
        assert joint["component_count"] == len(ledger["components"]) == len(original["components"])
        assert pointer(ledger_document, joint["ledger_components_pointer"]) == ledger["components"]
        assert pointer(source_document, joint["source_components_pointer"]) == original["components"]
        raw_source = source["original_source"]
        assert sha(output/raw_source["path"]) == raw_source["sha256"]
        assert {row["key"]["state_id"] for row in ledger["resolved"]["by_state"]} == {"M:a", "M:b"}
        assert source["source_association"] == ledger["source_association"]
        assert source["geometry"] == ledger["geometry"]
        assert source["resolved_total"] == ledger["resolved"]["total"]
        assert source["matched_uniform_width_total"] == ledger["matched_uniform_width"]["total"]
        assert source["validation_summary"]["algebraically_consistent"] is True
        assert source["validation_summary"]["source_association_bound"] is True


def test_preserves_exact_export_bytes_and_admission_hashes(command, exported, tmp_path):
    _, export = exported
    output = tmp_path/"admitted"
    original = {path.relative_to(export).as_posix(): path.read_bytes() for path in export.rglob("*") if path.is_file()}
    process = invoke(command, export, output)
    assert process.returncode == 0, process.stderr
    report = read_json(output/"admission.json")
    assert set(report["input_export_sha256"]) == set(original)
    for name, content in original.items():
        assert (output/"export"/name).read_bytes() == content
        assert (export/name).read_bytes() == content
        assert sha(export/name) == report["input_export_sha256"][name]
    for name, digest in report["output_sha256"].items():
        assert sha(output/name) == digest


@pytest.mark.parametrize("expected", ["0"*64, "A"*64, "f"*63, "g"*64])
def test_required_external_report_pin_rejects_wrong_or_invalid_values(command, exported, tmp_path, expected):
    _, export = exported
    output = tmp_path/"admitted"
    assert_rejected(invoke(command, export, output, expected=expected), output)


def test_missing_external_report_pin_rejects(command, exported, tmp_path):
    _, export = exported
    output = tmp_path/"admitted"
    assert_rejected(invoke(command, export, output, omit_pin=True), output)


def test_rebound_report_cannot_use_the_previous_trusted_pin(command, exported, tmp_path):
    kind, export = exported
    pin = sha(export/"report.json")
    mutate_ledger(export, kind, lambda ledger: ledger["resolved"]["total"]["quantities"]["released_number"].update(value=1.))
    assert sha(export/"report.json") != pin
    output = tmp_path/"admitted"
    assert_rejected(invoke(command, export, output, expected=pin), output)


@pytest.mark.parametrize("mutation", ["count", "projection", "boolean", "catalog", "area", "current", "state",
                                      "negative_zero", "summary_only", "missing_component", "source_provenance"])
def test_rebound_report_pin_does_not_bypass_original_source_or_exact_algebra(command, exported, tmp_path, mutation):
    kind, export = exported
    def change(ledger):
        total = ledger["resolved"]["total"]
        if mutation == "count":
            total["quantities"]["released_number"]["exact_sum"]["numerator"] = "1"
        elif mutation == "projection":
            total["quantities"]["released_number"]["projection"]["underflowed"] = True
        elif mutation == "boolean":
            total["quantities"]["released_number"]["available"] = 1
        elif mutation == "catalog":
            ledger["unqualified_contributors"].append({"catalog_index": 0, "quantity_name": "released_number"})
        elif mutation == "area":
            total["area"]["exact_m2"]["numerator"] = "1"
        elif mutation == "current":
            ledger["resolved"]["patch_state"][0]["key"]["current_index"] = 99
        elif mutation == "state":
            ledger["resolved"]["patch_state"][0]["key"]["state_id"] = "undeclared:state"
        elif mutation == "negative_zero":
            assert ledger["geometry"]["z_edges_m"][0] == 0.
            ledger["geometry"]["z_edges_m"][0] = -0.
        elif mutation == "summary_only":
            ledger["components_included"] = False
            ledger.pop("components")
        elif mutation == "missing_component":
            ledger["components"].pop()
        else:
            ledger["source_association"]["parameter_provenance"] = "Different source association"
    mutate_ledger(export, kind, change)
    if mutation in {"projection", "current", "state"}:
        refresh_csv_projection(export, kind)
    output = tmp_path/"admitted"
    assert_rejected(invoke(command, export, output), output)


def test_csv_must_match_exact_ledger_projections_after_rebinding(command, exported, tmp_path):
    _, export = exported
    path = export/"ledger.csv"
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        fields, rows = reader.fieldnames, list(reader)
    rows[0]["value"] = "12345"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    report = read_json(export/"report.json")
    report["output_sha256"]["ledger.csv"] = sha(path)
    write_json(export/"report.json", report)
    output = tmp_path/"admitted"
    assert_rejected(invoke(command, export, output), output)


@pytest.mark.parametrize("mutation", ["case_id", "window_index", "window_label", "crosswindow_ledger"])
def test_case_and_window_identity_cannot_be_rebound_to_another_source(command, fresh_exports, tmp_path, mutation):
    export = tmp_path/"ledger-export"
    shutil.copytree(fresh_exports["study"], export)
    document = read_json(export/"case_fresh.json")
    if mutation == "case_id":
        document["case_id"] = "unassociated_case"
    elif mutation == "window_index":
        document["windows"][0]["window_index"] = 1
    elif mutation == "window_label":
        document["windows"][0]["window"] = copy.deepcopy(document["windows"][1]["window"])
    else:
        document["windows"][0]["ledger"] = copy.deepcopy(document["windows"][1]["ledger"])
    rebound_file(export, "case_fresh.json", document)
    refresh_csv_projection(export, "study")
    output = tmp_path/"admitted"
    assert_rejected(invoke(command, export, output), output)


@pytest.mark.parametrize("mutation", ["duplicate_json", "nonfinite_json", "duplicate_report", "nonfinite_report",
                                      "unsafe_path", "input_symlink"])
def test_unsafe_inputs_reject_even_with_current_report_pin(command, exported, tmp_path, mutation):
    kind, export = exported
    if mutation in {"duplicate_json", "nonfinite_json", "duplicate_report", "nonfinite_report"}:
        path = export/("report.json" if mutation.endswith("report") else ledger_path(kind))
        original = path.read_text()
        prefix = '"schema":"duplicate",' if mutation.startswith("duplicate") else '"invalid":NaN,'
        path.write_text("{"+prefix+original[1:])
        if path.name != "report.json":
            report = read_json(export/"report.json")
            report["output_sha256"][path.name] = sha(path)
            write_json(export/"report.json", report)
    elif mutation == "unsafe_path":
        report = read_json(export/"report.json")
        report["output_sha256"]["../escape.json"] = "1"*64
        write_json(export/"report.json", report)
    else:
        original = export/ledger_path(kind)
        copied = tmp_path/"external-ledger.json"
        shutil.move(original, copied)
        original.symlink_to(copied)
    output = tmp_path/"admitted"
    assert_rejected(invoke(command, export, output), output)


def test_oversize_json_rejects_with_a_bounded_read(command, exported, tmp_path, monkeypatch):
    _, export = exported
    limit = 128
    assert (export/"report.json").stat().st_size > limit
    monkeypatch.setattr(command, "MAX_JSON_BYTES", limit)
    output = tmp_path/"admitted"
    process = invoke(command, export, output)
    assert_rejected(process, output)
    assert "Traceback" not in process.stderr


def test_deep_json_rejects_cleanly(command, exported, tmp_path):
    _, export = exported
    (export/"report.json").write_text('{"deep":'+"["*2000+"0"+"]"*2000+"}")
    output = tmp_path/"admitted"
    process = invoke(command, export, output)
    assert_rejected(process, output)
    assert "Traceback" not in process.stderr


@pytest.mark.parametrize("location", ["report", "ledger"])
def test_overflowing_json_number_literal_rejects_cleanly(command, exported, tmp_path, location):
    kind, export = exported
    path = export/("report.json" if location == "report" else ledger_path(kind))
    original = path.read_text()
    path.write_text('{"unused_overflow":1e999,'+original[1:])
    if location == "ledger":
        report = read_json(export/"report.json")
        report["output_sha256"][path.name] = sha(path)
        write_json(export/"report.json", report)
    output = tmp_path/"admitted"
    process = invoke(command, export, output)
    assert_rejected(process, output)
    assert "Traceback" not in process.stderr


@pytest.mark.parametrize("kind", ["directory", "file", "broken_symlink"])
def test_existing_output_is_preserved(command, exported, tmp_path, kind):
    _, export = exported
    output = tmp_path/"admitted"
    if kind == "directory":
        output.mkdir()
        (output/"keep").write_bytes(b"preserve")
    elif kind == "file":
        output.write_bytes(b"preserve")
    else:
        output.symlink_to(tmp_path/"missing")
    assert invoke(command, export, output).returncode != 0
    if kind == "directory":
        assert (output/"keep").read_bytes() == b"preserve"
    elif kind == "file":
        assert output.read_bytes() == b"preserve"
    else:
        assert output.is_symlink()


@pytest.mark.parametrize("location", ["original", "copy"])
def test_late_export_mutation_rolls_back(command, exported, tmp_path, monkeypatch, location):
    _, export = exported
    original_copy = command.shutil.copyfile
    changed = []
    def mutate_after_copy(source, destination, *args, **kwargs):
        result = original_copy(source, destination, *args, **kwargs)
        source = Path(source)
        if export in source.parents and not changed:
            target = source if location == "original" else Path(destination)
            with target.open("ab") as stream:
                stream.write(b" ")
            changed.append(target)
        return result
    monkeypatch.setattr(command.shutil, "copyfile", mutate_after_copy)
    output = tmp_path/"admitted"
    process = invoke(command, export, output)
    assert changed, "late mutation hook was not exercised"
    assert_rejected(process, output)
