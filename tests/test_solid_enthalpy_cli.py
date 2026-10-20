"""Independent Mills-expression, bounded-reference and fresh-output checks."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools/inspect_solid_enthalpy.py"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def configuration(tmp_path, **changes):
    config = json.loads((ROOT / "configs/mills_enthalpy_example.json").read_text())
    config.update(changes)
    path = tmp_path / "enthalpy_parameters.json"
    path.write_text(json.dumps(config) + "\n")
    return path, config


def execute(path, output, *extra):
    return subprocess.run(
        [sys.executable, str(TOOL), "--config", str(path), "--output-dir", str(output), *extra],
        capture_output=True, text=True,
    )


@pytest.mark.parametrize("convention", ["printed_eq11", "rezeroed_eq11"])
def test_evaluated_source_export_preserves_printed_and_declared_references(tmp_path, convention):
    path, config = configuration(tmp_path)
    config["model"]["reference_convention"] = convention
    config["temperatures_k"] = [300., float(np.nextafter(300., np.inf)), 550., 1000., 1500.]
    path.write_text(json.dumps(config) + "\n")
    output = tmp_path / "enthalpy"
    result = execute(path, output, "--plot")
    assert result.returncode == 0, result.stderr
    report = json.loads((output / "report.json").read_text())
    assert report["schema"] == "mhdlab-solid-enthalpy-inspection-v1"
    assert report["input_config_sha256"] == sha(path) == sha(output / "input_config.json")
    assert report["input_manifest_sha256"] == sha(ROOT / "data/materials/mills_2004_source_manifest.json")
    assert (output / "enthalpy.png").is_file()
    assert set(report["output_sha256"]) == {item.name for item in output.iterdir() if item.name != "report.json"}
    for name, digest in report["output_sha256"].items():
        assert sha(output / name) == digest
    for name, digest in report["source_sha256"].items():
        assert sha(ROOT / name) == digest
    metadata = report["model_metadata"]
    assert metadata["parameter_status"] == "evaluated_source_correlation"
    assert metadata["cp_is_cv"] is False
    assert metadata["enthalpy_is_internal_energy"] is False
    assert metadata["actual_dissertation_foil_validated"] is False
    assert metadata["density_or_fusion_supplied"] is False
    assert metadata["extrapolation_or_clamping"] is False
    assert metadata["adopted_as_solver_default"] is False
    assert metadata["material_family"] == "austenitic stainless steel"
    assert metadata["phase"] == "solid"
    assert metadata["reference_convention"] == convention
    assert metadata["units"]["cp_j_kg_k"] == "J kg^-1 K^-1"
    assert metadata["units"]["enthalpy_j_kg"] == "J kg^-1"
    assert metadata["supported_temperature_range_k"] == [300., 1500.]
    values = report["values"]
    temperature = np.array(config["temperatures_k"])
    cp = 472 + .136 * temperature - 2.82e6 / temperature**2
    printed = 472 * temperature + .068 * temperature**2 + 2.82e6 / temperature - 156000
    reference = config["model"]["reference_temperature_k"]
    delta = ((temperature - reference)
             * (472 + .068 * (temperature + reference) - 2.82e6 / (temperature * reference)))
    np.testing.assert_array_equal(values["temperature_k"], temperature)
    np.testing.assert_allclose(values["cp_j_kg_k"], cp, rtol=2e-15)
    np.testing.assert_allclose(values["printed_enthalpy_j_kg"], printed, rtol=2e-15)
    np.testing.assert_allclose(values["delta_enthalpy_from_reference_j_kg"], delta, rtol=2e-15)
    np.testing.assert_allclose(values["enthalpy_j_kg"], printed if convention == "printed_eq11" else delta, rtol=2e-15)
    np.testing.assert_allclose(values["roundtrip_temperature_k"], temperature, rtol=8e-16)
    assert values["roundtrip_temperature_k"][0] == temperature[0]
    assert values["roundtrip_temperature_k"][-1] == temperature[-1]
    if convention == "rezeroed_eq11":
        assert values["enthalpy_j_kg"][0] == 0
        assert values["enthalpy_j_kg"][1] > 0
    else:
        assert values["enthalpy_j_kg"][0] != 0
    assert metadata["printed_enthalpy_at_reference_j_kg"] != 0
    with (output / "enthalpy.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == len(temperature)
    assert set(rows[0]) == set(values)
    for index, row in enumerate(rows):
        for key, expected in values.items():
            assert float(row[key]) == expected[index]


@pytest.mark.parametrize("kind", ["query_below", "query_above", "boolean_query", "unsupported_phase",
                                  "unqualified_family", "reference_outside", "source_domain", "unknown_field"])
def test_inputs_cannot_extrapolate_substitute_phase_or_drop_qualifications(tmp_path, kind):
    path, config = configuration(tmp_path)
    if kind == "query_below":
        config["temperatures_k"] = [float(np.nextafter(config["model"]["temperature_min_k"], -np.inf))]
    elif kind == "query_above":
        config["temperatures_k"] = [float(np.nextafter(config["model"]["temperature_max_k"], np.inf))]
    elif kind == "boolean_query":
        config["temperatures_k"] = [True]
    elif kind == "unsupported_phase":
        config["model"]["expected_phase"] = "liquid"
    elif kind == "unqualified_family":
        config["model"]["expected_material_family"] = "304L"
    elif kind == "reference_outside":
        config["model"]["reference_temperature_k"] = 298.
    elif kind == "source_domain":
        config["model"]["temperature_max_k"] = 1701.
    else:
        config["model"]["density_kg_m3"] = 7900.
    path.write_text(json.dumps(config))
    output = tmp_path / "enthalpy"
    result = execute(path, output)
    assert result.returncode != 0
    assert not output.exists()


def test_changed_source_manifest_requires_review_instead_of_automatic_adoption(tmp_path):
    path, _ = configuration(tmp_path)
    manifest = tmp_path / "changed_manifest.json"
    source = (ROOT / "data/materials/mills_2004_source_manifest.json").read_bytes()
    manifest.write_bytes(source + b"\n")
    output = tmp_path / "enthalpy"
    result = execute(path, output, "--manifest", str(manifest))
    assert result.returncode != 0
    assert "SHA256" in result.stderr
    assert not output.exists()


@pytest.mark.parametrize("existing", ["directory", "file", "symlink"])
def test_existing_destination_and_contents_are_preserved(tmp_path, existing):
    path, _ = configuration(tmp_path)
    output = tmp_path / "enthalpy"
    sentinel = b"retain prior material evidence"
    if existing == "directory":
        output.mkdir()
        (output / "sentinel").write_bytes(sentinel)
    elif existing == "file":
        output.write_bytes(sentinel)
    else:
        output.symlink_to(tmp_path / "missing-target", target_is_directory=True)
    result = execute(path, output)
    assert result.returncode != 0
    assert "new path" in result.stderr
    if existing == "directory":
        assert list(output.iterdir()) == [output / "sentinel"]
        assert (output / "sentinel").read_bytes() == sentinel
    elif existing == "file":
        assert output.read_bytes() == sentinel
    else:
        assert output.is_symlink()
        assert not (tmp_path / "missing-target").exists()


@pytest.mark.parametrize("kind", ["duplicate", "nonfinite"])
def test_duplicate_or_nonfinite_json_creates_no_result(tmp_path, kind):
    path, config = configuration(tmp_path)
    if kind == "duplicate":
        text = path.read_text().replace('"schema":', f'"schema": {json.dumps(config["schema"])}, "schema":', 1)
    else:
        config["temperatures_k"] = [float("nan")]
        text = json.dumps(config)
    path.write_text(text)
    output = tmp_path / "enthalpy"
    result = execute(path, output)
    assert result.returncode != 0
    assert not output.exists()
