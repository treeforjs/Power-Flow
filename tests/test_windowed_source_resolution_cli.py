"""Independent integration checks for preserved source-resolution comparisons."""

import copy
from collections.abc import Mapping
from dataclasses import replace
import hashlib
import importlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from mhdlab.constants import KB, QE
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_planar_io import write_temkin_planar_release_history
from mhdlab.temkin_release import TemkinReservoir


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools/compare_windowed_source_resolution.py"
PROVENANCE = "Manufactured preserved resolution fixture with two constant-temperature regions"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    return path


@pytest.fixture(scope="module")
def preserved_sources(tmp_path_factory):
    directory = tmp_path_factory.mktemp("source-resolution-inputs")
    definition = write_json(directory / "fixture_definition.json", {
        "status": "manufactured", "initial_density_m2": 1.5e17,
        "temperature_regions_k": [400., 800.], "rate_regions_s": [1e7, 2e7],
        "interpretation": "Independent piecewise-constant region subdivision equality control",
    })
    base = FlashSurfaceHistory([0., 1e-7], [0., 1.], [[300.], [300.]], "synthetic",
                               "preserved-source-resolution-test", "Manufactured prescribed boundary",
                               "Manufactured stationary single face")
    pool = TemkinReservoir("resolution-pool", "M:resolution", 4e7, (KB/QE)*400*math.log(4),
                           .09, .02, 0., 1e18, 1e18, "Manufactured two known release rates",
                           "Manufactured supplied internal and desorption costs")
    records = {}
    for case_id, edges, temperatures, initial in (
        ("coarse", [-.003, -.001, .002], [400., 800.], [1.5e17, 1.5e17]),
        ("fine", [-.003, -.002, -.001, .0005, .002], [400., 400., 800., 800.], [1.5e17]*4),
        ("nonuniform", [-.003, -.001, .002], [400., 800.], [1e17, 2e17]),
        ("shifted", [-.002, 0., .003], [400., 800.], [1.5e17, 1.5e17]),
    ):
        surface = PlanarSurfaceHistory([0., 1e-7], edges, [.001, .004],
                                      [[temperatures], [temperatures]], base.attributes)
        history = evolve_temkin_planar_release(surface, [[initial]], [pool], parameter_provenance=PROVENANCE)
        source = directory / f"{case_id}.h5"
        write_temkin_planar_release_history(source, history)
        resolution = {"width_bins": len(edges)-1, "current_bins": 1,
                      "snapshot_intervals": 1, "subdivisions_per_interval": 1}
        evidence = write_json(directory / f"{case_id}.json", {
            "schema": "mhdlab-preserved-source-association-v1", "parameter_status": "manufactured",
            "parameter_provenance": PROVENANCE, "source_sha256": sha(source), "resolution": resolution,
            "preserved_inputs_sha256": {definition.name: sha(definition)}, "source_evolution_repeated": False,
            "audit_note": "New small manufactured fixture; exact release rates are independent controls",
        })
        records[case_id] = {"case_id": case_id, "source_file": str(source), "source_sha256": sha(source),
                            "evidence_file": str(evidence), "evidence_sha256": sha(evidence),
                            "resolution": resolution}
    return records


def configuration(records, *, comparisons=True):
    return {
        "schema": "mhdlab-windowed-source-resolution-config-v1", "parameter_status": "manufactured",
        "conditioning_provenance": "Manufactured absolute SOURCE emission windows, independent of photon gates",
        "comparison_provenance": "Manufactured exact patch subdivision control; no physical convergence claim",
        "relative_time_tolerance": 1e-12, "log_mass_tolerance": 1e-7,
        "windows": [{"start_s": 1e-8, "stop_s": 3e-8}, {"start_s": 6e-8, "stop_s": 9e-8}],
        "orders": [4, 8, 16],
        "cases": [copy.deepcopy(records["coarse"]), copy.deepcopy(records["fine"])],
        "comparisons": ([{"estimate_case_id": "coarse", "reference_case_id": "fine", "family": "width"}]
                        if comparisons else []),
    }


def argv(config, output):
    return [str(SCRIPT), "--config", str(config), "--output-dir", str(output)]


def run_cli(tmp_path, config, *, output=None, raw=None):
    path = tmp_path / "config.json"
    if raw is None:
        write_json(path, config)
    else:
        path.write_text(raw)
    output = tmp_path / "result" if output is None else output
    process = subprocess.run([sys.executable, "-W", "error", *argv(path, output)],
                             cwd=ROOT, text=True, capture_output=True, timeout=60)
    return process, output, path


def assert_no_output(process, output):
    assert process.returncode != 0, process.stdout
    assert not output.exists() and not output.is_symlink()
    assert not list(output.parent.glob(".windowed-source-resolution-*"))


@pytest.fixture(scope="module")
def compared(preserved_sources, tmp_path_factory):
    directory = tmp_path_factory.mktemp("source-resolution-result")
    process, output, config_path = run_cli(directory, configuration(preserved_sources))
    assert process.returncode == 0, process.stderr
    return output, json.loads((output / "report.json").read_text()), config_path


def test_preserved_inputs_and_fixed_order_differences(preserved_sources, compared):
    output, report, config = compared
    assert report["schema"] == "mhdlab-windowed-source-resolution-report-v1"
    assert report["parameter_status"] == "manufactured"
    assert report["source_evolution_repeated"] is False and report["host_budget_accepted"] is False
    assert sha(config) == report["input_config_sha256"] == sha(output / "input_config.json")
    for case in report["cases"]:
        name = case["case_id"]
        assert case == preserved_sources[name]
        assert sha(output / f"source_{name}.h5") == case["source_sha256"]
        assert sha(output / f"evidence_{name}.json") == case["evidence_sha256"]
    assert set(report["output_sha256"]) == {p.name for p in output.iterdir()} - {"report.json"}
    for name, digest in report["output_sha256"].items():
        assert sha(output / name) == digest
    for name, digest in report["source_sha256"].items():
        assert sha(ROOT / name) == digest
    cases = {name: json.loads((output / f"case_{name}.json").read_text()) for name in ("coarse", "fine")}
    assert len(report["source_comparisons"]) == 1
    differences = report["source_comparisons"][0]["differences"]
    assert {(d["window_index"], d["quadrature_order"]) for d in differences} == {
        (window, order) for window in range(2) for order in (4, 8, 16)}
    selected_metrics = ("released_number", "mean_time_s", "variance_time_s2", "mean_x_m", "variance_x_m2",
                        "mean_expected_translational_energy_j", "variance_kinetic_energy_j2",
                        "cov_time_expected_translational_energy_s_j", "cov_expected_translational_energy_x_j_m")
    for difference in differences:
        window = difference["window_index"]
        index = [4, 8, 16].index(difference["quadrature_order"])
        a = cases["coarse"]["windows"][window]["estimates"][index]
        b = cases["fine"]["windows"][window]["estimates"][index]
        for hypothesis in ("resolved", "matched_uniform_width"):
            for metric in selected_metrics:
                assert a[hypothesis][metric] == pytest.approx(b[hypothesis][metric], rel=2e-12, abs=0.)
            number = difference["estimate_minus_reference"][hypothesis]["released_number"]
            assert number["available"] is True
            assert number["signed_change"]["value"] == a[hypothesis]["released_number"]-b[hypothesis]["released_number"]
        assert "mean_time_offset_s" not in difference["estimate_minus_reference"]["resolved"]


def test_analytic_absolute_windows_and_paired_width_null(compared):
    output, _, _ = compared
    case = json.loads((output / "case_coarse.json").read_text())
    assert case["plan"]["window_indices"] == [0, 1]
    for entry, (start, stop) in zip(case["windows"], [(1e-8, 3e-8), (6e-8, 9e-8)]):
        length = stop-start
        numbers, means, variances = [], [], []
        for initial_number, rate in [(9e11, 1e7), (1.35e12, 2e7)]:
            exponential = math.expm1(rate*length)
            numbers.append(initial_number*(math.exp(-rate*start)-math.exp(-rate*stop)))
            means.append(start+1/rate-length/exponential)
            variances.append(1/rate**2-length**2*math.exp(rate*length)/exponential**2)
        count = math.fsum(numbers)
        p = numbers[1]/count
        mean_time = (1-p)*means[0]+p*means[1]
        variance_time = (1-p)*variances[0]+p*variances[1]+p*(1-p)*(means[1]-means[0])**2
        estimate = entry["estimates"][-1]
        resolved, uniform = estimate["resolved"], estimate["matched_uniform_width"]
        assert resolved["released_number"] == pytest.approx(count, rel=3e-10, abs=0.)
        assert resolved["mean_time_s"] == pytest.approx(mean_time, rel=3e-10, abs=0.)
        assert resolved["variance_time_s2"] == pytest.approx(variance_time, rel=3e-9, abs=0.)
        assert resolved["mean_x_m"] == pytest.approx(-.002+.0025*p, rel=3e-12, abs=0.)
        assert uniform["mean_x_m"] == pytest.approx(-.0005, rel=3e-14, abs=0.)
        assert uniform["variance_x_m2"] == pytest.approx(.005**2/12, rel=3e-14, abs=0.)
        assert resolved["mean_expected_translational_energy_j"] == pytest.approx((800+800*p)*KB, rel=3e-12, abs=0.)
        assert resolved["variance_kinetic_energy_j2"] == pytest.approx(
            (640000*p*(1-p)+320000+960000*p)*KB**2, rel=3e-12, abs=0.)
        assert resolved["second_raw_kinetic_energy_j2"] == pytest.approx((960000+2880000*p)*KB**2, rel=3e-12, abs=0.)
        assert resolved["cov_expected_translational_energy_x_j_m"] == pytest.approx(2*p*(1-p)*KB, rel=3e-12, abs=0.)
        assert resolved["cov_time_expected_translational_energy_s_j"] == pytest.approx(
            p*(1-p)*(means[1]-means[0])*800*KB, rel=3e-9, abs=0.)
        assert uniform["cov_expected_translational_energy_x_j_m"] == 0.
        for metric in ("released_number", "mean_time_s", "variance_time_s2", "mean_z_m", "variance_z_m2",
                       "mean_expected_translational_energy_j", "variance_kinetic_energy_j2",
                       "cov_time_expected_translational_energy_s_j"):
            assert uniform[metric] == resolved[metric]
        assert len(entry["conditioned_source"]["components"]) == 2
        for change, (coarse, fine) in zip(entry["quadrature_changes"], [(4, 8), (8, 16)]):
            assert (change["coarse_order"], change["refined_order"]) == (coarse, fine)
            a, b = [e["resolved"]["variance_time_s2"] for e in entry["estimates"]
                    if e["quadrature_order"] in (coarse, fine)]
            assert change["refined_minus_coarse"]["resolved"]["variance_time_s2"]["signed_change"]["value"] == b-a


def test_single_case_and_relative_input_paths(preserved_sources, tmp_path):
    config = configuration(preserved_sources, comparisons=False)
    config["cases"] = [config["cases"][0]]
    for key in ("source_file", "evidence_file"):
        original = Path(config["cases"][0][key])
        shutil.copyfile(original, tmp_path / original.name)
        config["cases"][0][key] = original.name
    process, output, _ = run_cli(tmp_path, config)
    assert process.returncode == 0, process.stderr
    report = json.loads((output / "report.json").read_text())
    assert report["source_comparisons"] == [] and len(report["cases"]) == 1


@pytest.mark.parametrize("mutation", ["duplicate", "unknown", "boolean_tolerance", "nonfinite", "outside_window",
                                      "duplicate_orders", "unsafe_id", "two_knobs", "wrong_family", "window_work_limit"])
def test_invalid_configuration_is_rejected(preserved_sources, tmp_path, mutation):
    config = configuration(preserved_sources)
    raw = None
    if mutation == "duplicate":
        raw = '{"schema":' + json.dumps(config["schema"]) + ',' + json.dumps(config)[1:]
    elif mutation == "unknown":
        config["auto_evolve"] = True
    elif mutation == "boolean_tolerance":
        config["relative_time_tolerance"] = True
    elif mutation == "nonfinite":
        raw = json.dumps(config).replace('"log_mass_tolerance": 1e-07', '"log_mass_tolerance": NaN')
        assert "NaN" in raw
    elif mutation == "outside_window":
        config["windows"] = [{"start_s": 0., "stop_s": 2e-7}]
    elif mutation == "duplicate_orders":
        config["orders"] = [8, 8]
    elif mutation == "unsafe_id":
        config["cases"][0]["case_id"] = "../escape"
    elif mutation == "two_knobs":
        config["cases"][1]["resolution"]["current_bins"] = 2
    elif mutation == "wrong_family":
        config["comparisons"][0]["family"] = "current"
    elif mutation == "window_work_limit":
        config["windows"] *= 9
    process, output, _ = run_cli(tmp_path, config, raw=raw)
    assert_no_output(process, output)
    assert "error:" in process.stderr


@pytest.mark.parametrize("mutation", ["source_hash", "evidence_hash", "status", "provenance", "resolution",
                                      "repeated", "empty_input_hashes", "source_symlink", "evidence_symlink"])
def test_detached_source_evidence_is_rejected(preserved_sources, tmp_path, mutation):
    config = configuration(preserved_sources, comparisons=False)
    config["cases"] = [config["cases"][0]]
    record = config["cases"][0]
    if mutation == "source_hash":
        record["source_sha256"] = "0"*64
    elif mutation == "evidence_hash":
        record["evidence_sha256"] = "0"*64
    elif mutation.endswith("symlink"):
        key = "source_file" if mutation == "source_symlink" else "evidence_file"
        link = tmp_path / Path(record[key]).name
        link.symlink_to(record[key])
        record[key] = str(link)
    else:
        evidence = json.loads(Path(record["evidence_file"]).read_text())
        if mutation == "status":
            evidence["parameter_status"] = "sourced"
        elif mutation == "provenance":
            evidence["parameter_provenance"] = "Detached parameter provenance"
        elif mutation == "resolution":
            evidence["resolution"]["width_bins"] = 99
        elif mutation == "repeated":
            evidence["source_evolution_repeated"] = True
        elif mutation == "empty_input_hashes":
            evidence["preserved_inputs_sha256"] = {}
        path = write_json(tmp_path / "altered_evidence.json", evidence)
        record["evidence_file"], record["evidence_sha256"] = str(path), sha(path)
    process, output, _ = run_cli(tmp_path, config)
    assert_no_output(process, output)


@pytest.mark.parametrize("source_case", ["nonuniform", "shifted"])
def test_mechanically_incompatible_sources_are_rejected(preserved_sources, tmp_path, source_case):
    config = configuration(preserved_sources)
    config["cases"][0] = copy.deepcopy(preserved_sources[source_case])
    config["comparisons"][0]["estimate_case_id"] = source_case
    process, output, _ = run_cli(tmp_path, config)
    assert_no_output(process, output)
    assert "uniform initial" in process.stderr if source_case == "nonuniform" else "domains" in process.stderr


def test_cli_does_not_advance_sources_or_sample_packets(preserved_sources, tmp_path, monkeypatch):
    module = importlib.import_module("tools.compare_windowed_source_resolution")
    evolution = importlib.import_module("mhdlab.temkin_planar")
    sampling = importlib.import_module("mhdlab.temkin_windowed_boundary")

    def forbidden(*args, **kwargs):
        raise AssertionError("preserved comparison invoked evolution or sampling")

    monkeypatch.setattr(evolution, "evolve_temkin_planar_release", forbidden)
    monkeypatch.setattr(sampling, "sample_temkin_windowed_emission_packets", forbidden)
    path = write_json(tmp_path / "config.json", configuration(preserved_sources))
    output = tmp_path / "result"
    monkeypatch.setattr(sys, "argv", argv(path, output))
    module.main()
    assert (output / "report.json").is_file()


@pytest.mark.parametrize("kind", ["directory", "file", "broken_symlink"])
def test_existing_output_is_preserved(preserved_sources, tmp_path, kind):
    output = tmp_path / "result"
    if kind == "directory":
        output.mkdir()
        (output / "sentinel").write_text("retain")
    elif kind == "file":
        output.write_text("retain")
    else:
        output.symlink_to(tmp_path / "absent")
    process, _, _ = run_cli(tmp_path, configuration(preserved_sources), output=output)
    assert process.returncode != 0 and "new path" in process.stderr
    if kind == "directory":
        assert (output / "sentinel").read_text() == "retain"
    elif kind == "file":
        assert output.read_text() == "retain"
    else:
        assert output.is_symlink() and not output.exists()


@pytest.mark.parametrize("changed_file", ["source_file", "evidence_file"])
def test_first_input_changed_after_its_check_is_rejected(preserved_sources, tmp_path, monkeypatch, changed_file):
    module = importlib.import_module("tools.compare_windowed_source_resolution")
    config = configuration(preserved_sources)
    original_path = Path(config["cases"][0][changed_file])
    local = tmp_path / original_path.name
    shutil.copyfile(original_path, local)
    config["cases"][0][changed_file] = str(local)
    path = write_json(tmp_path / "config.json", config)
    original_write = Path.write_text
    output = tmp_path / "result"

    def late_change(destination, *args, **kwargs):
        result = original_write(destination, *args, **kwargs)
        if destination.name == "case_fine.json":
            with local.open("ab") as stream:
                stream.write(b"changed after earlier case check")
        return result

    monkeypatch.setattr(Path, "write_text", late_change)
    monkeypatch.setattr(sys, "argv", argv(path, output))
    with pytest.raises(SystemExit) as error:
        module.main()
    assert error.value.code == 2
    assert not output.exists() and not list(tmp_path.glob(".windowed-source-resolution-*"))


def test_partial_output_write_rolls_back(preserved_sources, tmp_path, monkeypatch):
    module = importlib.import_module("tools.compare_windowed_source_resolution")
    path = write_json(tmp_path / "config.json", configuration(preserved_sources))
    output = tmp_path / "result"
    original_write = Path.write_text

    def failed_write(destination, *args, **kwargs):
        if destination.name == "case_fine.json":
            raise OSError("manufactured case write failure")
        return original_write(destination, *args, **kwargs)

    monkeypatch.setattr(Path, "write_text", failed_write)
    monkeypatch.setattr(sys, "argv", argv(path, output))
    with pytest.raises(SystemExit) as error:
        module.main()
    assert error.value.code == 2
    assert not output.exists() and not list(tmp_path.glob(".windowed-source-resolution-*"))


def test_common_interval_snapshot_config_accepts_compensating_counts(preserved_sources, tmp_path):
    module = importlib.import_module("tools.compare_windowed_source_resolution")
    config = configuration(preserved_sources)
    for case, snapshots, subdivisions in zip(config["cases"], (8, 16), (64, 32)):
        case["resolution"] = {"width_bins": 2, "current_bins": 1,
                              "snapshot_intervals": snapshots, "subdivisions_per_interval": subdivisions}
    config["comparisons"][0]["family"] = "snapshot_common_intervals"
    path = write_json(tmp_path / "config.json", config)
    _, parsed, _ = module.configuration(path)
    assert parsed["comparisons"][0]["family"] == "snapshot_common_intervals"


def test_common_interval_snapshot_sources_require_identical_time_edges():
    module = importlib.import_module("tools.compare_windowed_source_resolution")
    common = {"reservoirs": (), "parameter_provenance": "Manufactured shared fixture", "surface_attributes": {},
              "x_edges_m": np.array([0., 1.]), "z_edges_m": np.array([0., 1.]),
              "coverage_m2": np.ones((3, 1, 1, 1))}
    a = SimpleNamespace(**common, time_edges_s=np.array([0., .5, 1.]))
    b = SimpleNamespace(**common, time_edges_s=np.array([0., .4, 1.]))
    with pytest.raises(ValueError, match="identical|same|interval"):
        module.shared_inputs(a, b, "snapshot_common_intervals")


def test_changes_preserve_local_time_difference_and_relative_underflow():
    module = importlib.import_module("tools.compare_windowed_source_resolution")

    def estimate(local_mean):
        anchor = {"interval_origin_s": 1e300, "window_start_offset_s": 0.,
                  "mean_time_offset_s": local_mean}
        return {"units": {"mean_time_s": "s", "mean_time_offset_s": "s"},
                "resolved": {"mean_time_s": 1e300, "mean_time_offset_s": local_mean},
                "matched_uniform_width": {"mean_time_s": 1e300, "mean_time_offset_s": local_mean},
                "mean_time_anchor": anchor, "matched_uniform_width_mean_time_anchor": dict(anchor)}

    comparison = module.changes(estimate(2e-300), estimate(1e-300))
    for hypothesis in ("resolved", "matched_uniform_width"):
        assert "mean_time_offset_s" not in comparison[hypothesis]
        difference = comparison[hypothesis]["mean_time_s"]
        assert difference["signed_change"]["value"] == 1e-300
        assert difference["signed_change"]["underflowed"] is False
        assert difference["relative_to_reference"]["value"] == 0.
        assert difference["relative_to_reference"]["underflowed"] is True
        assert difference["relative_to_reference"]["exact_represented_input_fraction"]["numerator"] != "0"


def plain_metadata(value):
    if isinstance(value, Mapping):
        return {name: plain_metadata(item) for name, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain_metadata(item) for item in value]
    return value


@pytest.fixture(scope="module")
def common_grid_sources(tmp_path_factory):
    directory = tmp_path_factory.mktemp("source-resolution-common-grid")
    base = FlashSurfaceHistory([0., 1e-7], [0., 1.], [[300.], [300.]], "synthetic",
                               "common-grid-resolution-test", "Manufactured common-grid source",
                               "Manufactured stationary single face")
    surface = PlanarSurfaceHistory([0., .5e-7, 1e-7], [-.003, -.001, .002], [.001, .004],
                                  [[[400., 800.]], [[400., 800.]], [[400., 800.]]], base.attributes)
    pool = TemkinReservoir("resolution-pool", "M:resolution", 4e7, (KB/QE)*400*math.log(4),
                           .09, .02, 0., 1e18, 1e18, "Manufactured common-grid rates",
                           "Manufactured supplied energy costs")
    history = evolve_temkin_planar_release(surface, [[[1.5e17, 1.5e17]]], [pool], parameter_provenance=PROVENANCE)
    definition = write_json(directory / "common_grid_definition.json", {
        "status": "manufactured", "time_s": [0., .5e-7, 1e-7], "x_edges_m": [-.003, -.001, .002],
        "z_edges_m": [.001, .004], "temperature_k": [400., 800.], "initial_coverage_m2": 1.5e17,
        "attempt_frequency_s": 4e7, "activation_energy_ev": (KB/QE)*400*math.log(4),
        "coverage_coefficient": 0., "desorption_energy_ev": .09, "internal_energy_ev": .02,
    })
    records, histories = {}, {}
    for name, snapshots, subdivisions in [("coarse", 1, 2), ("fine", 2, 1)]:
        attributes = plain_metadata(history.surface_attributes)
        original_attributes = plain_metadata(history.surface_attributes)
        attributes["followup_common_grid_resampling"] = {
            "original_snapshot_intervals": snapshots,
            "equivalent_subdivisions_per_original_interval": subdivisions,
            "actual_integrator_snapshot_intervals": 2, "actual_integrator_subdivisions": 1,
            "additional_recovered_information": False, "temporal_smoothing": False,
            "original_snapshot_attributes": original_attributes,
        }
        derived = replace(history, surface_attributes=attributes)
        source = directory / f"{name}.h5"
        write_temkin_planar_release_history(source, derived)
        counts = {"width_bins": 2, "current_bins": 1, "snapshot_intervals": snapshots,
                  "subdivisions_per_interval": subdivisions}
        evidence = write_json(directory / f"{name}.json", {
            "schema": "mhdlab-preserved-source-association-v1", "parameter_status": "manufactured",
            "parameter_provenance": PROVENANCE, "source_sha256": sha(source), "resolution": counts,
            "preserved_inputs_sha256": {definition.name: sha(definition)},
            "source_evolution_repeated": False,
        })
        records[name] = {"case_id": name, "source_file": str(source), "source_sha256": sha(source),
                         "evidence_file": str(evidence), "evidence_sha256": sha(evidence), "resolution": counts}
        histories[name] = derived
    return records, histories


@pytest.mark.parametrize("mutation", ["none", "forged_original_count", "changed_original_attributes"])
def test_common_grid_metadata_qualifies_only_declared_count_changes(common_grid_sources, tmp_path, mutation):
    records, histories = common_grid_sources
    config = configuration(records)
    config["comparisons"][0]["family"] = "snapshot_common_intervals"
    if mutation != "none":
        attributes = plain_metadata(histories["fine"].surface_attributes)
        remap = attributes["followup_common_grid_resampling"]
        if mutation == "forged_original_count":
            remap["original_snapshot_intervals"] = 99
        else:
            remap["original_snapshot_attributes"]["run_id"] = "Different original physical history"
        changed = replace(histories["fine"], surface_attributes=attributes)
        source = tmp_path / "changed_source.h5"
        write_temkin_planar_release_history(source, changed)
        record = config["cases"][1]
        record["source_file"], record["source_sha256"] = str(source), sha(source)
        evidence = json.loads(Path(record["evidence_file"]).read_text())
        evidence["source_sha256"] = sha(source)
        evidence_path = write_json(tmp_path / "changed_evidence.json", evidence)
        record["evidence_file"], record["evidence_sha256"] = str(evidence_path), sha(evidence_path)
    process, output, _ = run_cli(tmp_path, config)
    if mutation == "none":
        assert process.returncode == 0, process.stderr
        report = json.loads((output / "report.json").read_text())
        assert report["source_comparisons"][0]["family"] == "snapshot_common_intervals"
        for case in report["cases"]:
            assert case["resolution"]["snapshot_intervals"]*case["resolution"]["subdivisions_per_interval"] == 2
    else:
        assert_no_output(process, output)
        expected = "resampling metadata" if mutation == "forged_original_count" else "surface metadata"
        assert expected in process.stderr


@pytest.mark.parametrize("control", ["current_refinement", "width_refinement", "changed_physical_metadata"])
def test_current_partition_metadata_changes_only_for_current_family(control):
    module = importlib.import_module("tools.compare_windowed_source_resolution")

    def source(z_edges):
        return SimpleNamespace(
            reservoirs=("same manufactured pool",), parameter_provenance="Same manufactured physics",
            surface_attributes={"physical_drive_reference": "Same manufactured drive",
                                "current_patch_intervals_m": [list(pair) for pair in zip(z_edges, z_edges[1:])]},
            time_edges_s=np.array([0., 1.]), x_edges_m=np.array([0., 1.]), z_edges_m=np.asarray(z_edges),
            coverage_m2=np.ones((2, 1, len(z_edges)-1, 1)),
        )

    a, b = source([0., 1.]), source([0., .4, 1.])
    if control == "changed_physical_metadata":
        b.surface_attributes["physical_drive_reference"] = "Different physical drive"
    if control == "current_refinement":
        module.shared_inputs(a, b, "current")
    else:
        with pytest.raises(ValueError, match="metadata|coordinates"):
            module.shared_inputs(a, b, "width" if control == "width_refinement" else "current")

