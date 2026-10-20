"""Bounded SI transport properties with explicit material and source context."""

from dataclasses import FrozenInstanceError
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from mhdlab.material_transport import BoundedTransportTable, load_flynn_304_transport


ROOT = Path(__file__).resolve().parents[1]
SOURCE_CSV = ROOT / "data/materials/flynn_1968_304_transport.csv"
SOURCE_JSON = ROOT / "data/materials/flynn_1968_304_transport.json"
VALUE_FIELDS = ("temperature_k", "resistivity_ohm_m", "thermal_conductivity_w_m_k",
                "electrical_conductivity_s_m")


def table(**overrides):
    parameters = dict(temperature_k=[300., 500., 900.],
                      resistivity_ohm_m=[2e-6, 6e-6, 10e-6],
                      thermal_conductivity_w_m_k=[5., 15., 35.],
                      material_grade="manufactured material",
                      sample_condition="explicit mathematical fixture; no physical sample",
                      phase="solid", parameter_status="manufactured",
                      source_reference="manufactured linear SI verification values",
                      source_data_sha256=hashlib.sha256(b"manufactured fixture").hexdigest())
    parameters.update(overrides)
    return BoundedTransportTable(**parameters)


def load(csv_path=SOURCE_CSV, metadata_path=SOURCE_JSON, **overrides):
    parameters = dict(expected_material_grade="AISI Type 304")
    parameters.update(overrides)
    return load_flynn_304_transport(csv_path, metadata_path, **parameters)


def copy_source(tmp_path):
    csv_path, metadata_path = tmp_path / SOURCE_CSV.name, tmp_path / SOURCE_JSON.name
    csv_path.write_bytes(SOURCE_CSV.read_bytes())
    metadata_path.write_bytes(SOURCE_JSON.read_bytes())
    return csv_path, metadata_path


def read_metadata(path):
    return json.loads(path.read_text())


def write_metadata(path, metadata):
    path.write_text(json.dumps(metadata, allow_nan=False))


def set_nested(record, keys, value):
    for key in keys[:-1]:
        record = record[key]
    record[keys[-1]] = value


def refresh_digest(csv_path, metadata_path):
    metadata = read_metadata(metadata_path)
    metadata["output"]["csv_sha256"] = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    write_metadata(metadata_path, metadata)


def write_rows(path, rows):
    with path.open("w", newline="") as handle:
        csv.writer(handle).writerows(rows)


def test_manufactured_piecewise_linear_rho_and_heat_conduction_then_reciprocal_sigma():
    model = table()
    values = model.evaluate([300., 400., 500., 700., 900.])
    np.testing.assert_array_equal(values.temperature_k, [300., 400., 500., 700., 900.])
    np.testing.assert_allclose(values.resistivity_ohm_m, np.array([2., 4., 6., 8., 10.]) * 1e-6,
                               rtol=3e-15, atol=0.)
    np.testing.assert_allclose(values.thermal_conductivity_w_m_k, [5., 10., 15., 25., 35.],
                               rtol=3e-15, atol=0.)
    np.testing.assert_allclose(values.electrical_conductivity_s_m,
                               1. / (np.array([2., 4., 6., 8., 10.]) * 1e-6), rtol=3e-15, atol=0.)
    assert values.electrical_conductivity_s_m[1] != pytest.approx((1 / 2e-6 + 1 / 6e-6) / 2)
    assert model.representation == "piecewise_linear_in_temperature"
    assert model.parameter_status == "manufactured"


def test_scalar_and_unsorted_repeated_vector_queries_preserve_shape_and_values():
    model = table()
    scalar = model.evaluate(400.)
    for name in VALUE_FIELDS:
        assert type(getattr(scalar, name)) is float
    assert scalar.resistivity_ohm_m == pytest.approx(4e-6, rel=3e-15, abs=0.)
    assert scalar.thermal_conductivity_w_m_k == pytest.approx(10., rel=3e-15, abs=0.)
    vector = model.evaluate([700., 300., 700.])
    np.testing.assert_array_equal(vector.temperature_k, [700., 300., 700.])
    np.testing.assert_allclose(vector.resistivity_ohm_m, [8e-6, 2e-6, 8e-6], rtol=3e-15, atol=0.)


def test_table_and_evaluation_arrays_are_copied_and_immutable():
    temperatures = np.array([300., 500., 900.])
    resistivity = np.array([2e-6, 6e-6, 10e-6])
    model = table(temperature_k=temperatures, resistivity_ohm_m=resistivity)
    temperatures[0], resistivity[0] = 1., 1.
    np.testing.assert_array_equal(model.temperature_k, [300., 500., 900.])
    np.testing.assert_array_equal(model.resistivity_ohm_m, [2e-6, 6e-6, 10e-6])
    for name in ("temperature_k", "resistivity_ohm_m", "thermal_conductivity_w_m_k"):
        with pytest.raises(ValueError):
            getattr(model, name).setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        model.material_grade = "different material"
    query = np.array([400., 700.])
    result = model.evaluate(query)
    query[0] = 500.
    np.testing.assert_array_equal(result.temperature_k, [400., 700.])
    for name in VALUE_FIELDS:
        with pytest.raises(ValueError):
            getattr(result, name).setflags(write=True)


def test_nested_source_metadata_is_preserved_and_deeply_immutable():
    source = {"sample": {"composition": None, "history": ["unknown", {"measured": False}]}}
    model = table(source_metadata=source)
    source["sample"]["history"][1]["measured"] = True
    assert model.source_metadata["sample"]["composition"] is None
    assert model.source_metadata["sample"]["history"][1]["measured"] is False
    with pytest.raises(TypeError):
        model.source_metadata["sample"]["composition"] = "assumed alloy"
    with pytest.raises(TypeError):
        model.source_metadata["sample"]["history"][1]["measured"] = True
    with pytest.raises((TypeError, AttributeError)):
        model.source_metadata["sample"]["history"].append("extra")


def test_metadata_export_is_independent_and_json_serializable():
    model = table(source_metadata={"sample": {"composition": None, "labels": ["unknown"]}})
    metadata = model.metadata
    json.dumps(metadata, allow_nan=False)
    metadata["units"]["resistivity_ohm_m"] = "ohm cm"
    metadata["source_metadata"]["sample"]["labels"].append("assumed")
    assert model.metadata["units"]["resistivity_ohm_m"] == "ohm m"
    assert model.source_metadata["sample"]["labels"] == ("unknown",)


def test_derived_electrical_conductivity_cannot_be_supplied_independently():
    result = table().evaluate(400.)
    with pytest.raises(TypeError):
        type(result)(temperature_k=400., resistivity_ohm_m=4e-6,
                     thermal_conductivity_w_m_k=10., electrical_conductivity_s_m=1.)
    with pytest.raises(FrozenInstanceError):
        result.electrical_conductivity_s_m = 1.


def test_tiny_temperature_bin_uses_convex_interpolation_without_slope_overflow():
    smallest = float(np.nextafter(0., 1.))
    model = table(temperature_k=[smallest, 3 * smallest],
                  resistivity_ohm_m=[1., 2.], thermal_conductivity_w_m_k=[1., 2.])
    values = model.evaluate(2 * smallest)
    assert values.resistivity_ohm_m == pytest.approx(1.5, rel=3e-15, abs=0.)
    assert values.thermal_conductivity_w_m_k == pytest.approx(1.5, rel=3e-15, abs=0.)
    assert values.electrical_conductivity_s_m == pytest.approx(2 / 3, rel=3e-15, abs=0.)


def test_unrepresentable_reciprocal_raises_on_evaluation_instead_of_returning_infinity():
    model = table(resistivity_ohm_m=[1e-320, 2e-320, 3e-320])
    with pytest.raises(ValueError):
        model.evaluate(400.)


def test_wide_temperature_bin_preserves_underflowed_weights_with_finite_property_contribution():
    # The manufactured property is exactly y(T)=T on this interval.  Its tiny
    # interpolation weight underflows, but weight*upper_property is finite.
    model = table(temperature_k=[1e-300, 1e300], resistivity_ohm_m=[1e-300, 1e300],
                  thermal_conductivity_w_m_k=[1e-300, 1e300])
    result = model.evaluate(2e-300)
    assert result.resistivity_ohm_m == pytest.approx(2e-300, rel=6e-13, abs=0.)
    assert result.thermal_conductivity_w_m_k == pytest.approx(2e-300, rel=6e-13, abs=0.)
    assert result.electrical_conductivity_s_m == pytest.approx(5e299, rel=6e-13, abs=0.)


def test_constant_maximum_finite_properties_remain_finite_when_interpolated():
    maximum = np.finfo(float).max
    model = table(temperature_k=[1., 2.], resistivity_ohm_m=[maximum, maximum],
                  thermal_conductivity_w_m_k=[maximum, maximum])
    result = model.evaluate(1.5)
    assert result.resistivity_ohm_m == maximum
    assert result.thermal_conductivity_w_m_k == maximum
    assert np.isfinite(result.electrical_conductivity_s_m)
    assert result.electrical_conductivity_s_m > 0.


@pytest.mark.parametrize("overrides", [
    {"temperature_k": []}, {"temperature_k": [300.]}, {"temperature_k": [300., 300., 900.]},
    {"temperature_k": [300., 200., 900.]}, {"temperature_k": [0., 500., 900.]},
    {"temperature_k": [-1., 500., 900.]}, {"temperature_k": [300., np.inf, 900.]},
    {"temperature_k": [300., np.nan, 900.]}, {"temperature_k": [[300., 500., 900.]]},
    {"temperature_k": [300., True, 900.]}, {"temperature_k": [300., 500.+0j, 900.]},
    {"temperature_k": ["300", "500", "900"]}, {"resistivity_ohm_m": [1e-6, 2e-6]},
    {"resistivity_ohm_m": [0., 2e-6, 3e-6]}, {"resistivity_ohm_m": [-1., 2e-6, 3e-6]},
    {"resistivity_ohm_m": [np.nan, 2e-6, 3e-6]}, {"resistivity_ohm_m": [True, 2e-6, 3e-6]},
    {"resistivity_ohm_m": [1j, 2e-6, 3e-6]}, {"resistivity_ohm_m": ["1e-6", "2e-6", "3e-6"]},
    {"thermal_conductivity_w_m_k": [1., 2.]}, {"thermal_conductivity_w_m_k": [1., 0., 3.]},
    {"thermal_conductivity_w_m_k": [1., -2., 3.]}, {"thermal_conductivity_w_m_k": [1., np.inf, 3.]},
    {"thermal_conductivity_w_m_k": [[1., 2., 3.]]}, {"phase": "liquid"}, {"phase": "Solid"},
    {"phase": np.array(["solid"])}, {"parameter_status": np.array(["manufactured"])},
    {"parameter_status": "validated"}, {"parameter_status": "measured"},
    {"material_grade": ""}, {"sample_condition": " "}, {"source_reference": ""},
    {"material_grade": None}, {"sample_condition": True}, {"source_reference": 1968},
    {"source_data_sha256": ""}, {"source_data_sha256": "not-a-digest"},
    {"source_data_sha256": "z" * 64}, {"source_data_sha256": "A" * 64}, {"source_data_sha256": None},
    {"source_metadata": {"invalid": np.nan}}, {"source_metadata": {"invalid": 1j}},
    {"source_metadata": []}, {"source_metadata": {"invalid": {1, 2}}},
    {"temperature_k": np.ma.array([300., 500., 900.], mask=[False, True, False])},
    {"resistivity_ohm_m": np.ma.array([1., 2., 3.], mask=[True, False, False])},
])
def test_table_rejects_invalid_shapes_values_and_missing_context(overrides):
    with pytest.raises(ValueError):
        table(**overrides)


@pytest.mark.parametrize("query", [
    299.99, 900.01, np.nextafter(300., -np.inf), np.nextafter(900., np.inf),
    [], [[400.]], [300., 901.], np.nan, np.inf, True, 400.+0j, "400",
    [300., False], [300., 500.+0j], ["300", "500"],
    np.ma.array([300., 500.], mask=[False, True]),
])
def test_evaluation_never_extrapolates_clamps_or_coerces_invalid_types(query):
    with pytest.raises(ValueError):
        table().evaluate(query)


def test_flynn_source_anchors_and_endpoints_are_preserved_with_explicit_labels():
    model = load()
    assert model.material_grade == "AISI Type 304"
    assert model.phase == "solid"
    assert model.parameter_status == "recovered_published_smoothed_table"
    assert "10.6028/NBS.RPT.9869" in model.source_reference
    assert model.source_data_sha256 == hashlib.sha256(SOURCE_CSV.read_bytes()).hexdigest()
    assert model.sample_condition.strip()
    np.testing.assert_allclose(model.temperature_k, [253.15, 273.15, 293.15, 313.15, 333.15,
                                                    353.15, 373.15, 393.15, 413.15, 423.15], rtol=0., atol=0.)
    np.testing.assert_allclose(model.resistivity_ohm_m[[0, 2, -1]], [6.82e-7, 7.17e-7, 8.19e-7], rtol=0., atol=0.)
    np.testing.assert_array_equal(model.thermal_conductivity_w_m_k[[0, 2, -1]], [13.9, 14.8, 17.9])
    ends = model.evaluate([253.15, 423.15])
    np.testing.assert_array_equal(ends.resistivity_ohm_m, [6.82e-7, 8.19e-7])
    np.testing.assert_array_equal(ends.thermal_conductivity_w_m_k, [13.9, 17.9])


def test_flynn_midbin_properties_are_interpolated_before_reciprocal_conversion():
    values = load().evaluate(283.15)
    expected_rho = (7.00e-7 + 7.17e-7) / 2
    assert values.resistivity_ohm_m == pytest.approx(expected_rho, rel=3e-15, abs=0.)
    assert values.thermal_conductivity_w_m_k == pytest.approx(14.55, rel=3e-15, abs=0.)
    assert values.electrical_conductivity_s_m == pytest.approx(1 / expected_rho, rel=3e-15, abs=0.)
    assert values.electrical_conductivity_s_m != pytest.approx((1 / 7e-7 + 1 / 7.17e-7) / 2, rel=1e-6)


def test_full_recovery_metadata_retains_unknown_sample_and_uncertainty_meaning():
    metadata = load().source_metadata
    assert metadata["observable"]["sample_composition"] is None
    assert metadata["observable"]["sample_thermal_history"] is None
    assert metadata["source"]["author"] == "D. R. Flynn"
    assert metadata["source_QA"]["render_recipe"]["used_PyMuPDF_version"] == "1.28.2"
    assert metadata["source_error_statements"]["confidence_level"] is None
    assert metadata["source_error_statements"]["standard_deviation_interpretation"] is False
    assert metadata["source_error_statements"]["between_heat_or_foil_variability_quantified"] is False
    assert "unknown" in load().sample_condition.lower()
    with pytest.raises(TypeError):
        metadata["observable"]["sample_composition"] = "assumed matching foil"


@pytest.mark.parametrize("grade", ["304", "AISI Type 304L", "AISI Type 316", "SS304", "", None])
def test_loader_requires_exact_explicit_material_grade(grade):
    with pytest.raises(ValueError):
        load(expected_material_grade=grade)


def test_loader_does_not_supply_a_default_expected_grade():
    with pytest.raises(TypeError):
        load_flynn_304_transport(SOURCE_CSV, SOURCE_JSON)


def test_loader_detects_exact_csv_content_change(tmp_path):
    csv_path, metadata_path = copy_source(tmp_path)
    csv_path.write_text(csv_path.read_text().replace("13.9", "13.8"))
    with pytest.raises(ValueError):
        load(csv_path, metadata_path)


def test_unchanged_source_pair_can_be_loaded_from_an_explicit_local_copy(tmp_path):
    csv_path, metadata_path = copy_source(tmp_path)
    model = load(csv_path, metadata_path)
    np.testing.assert_array_equal(model.temperature_k, load().temperature_k)
    assert model.source_data_sha256 == hashlib.sha256(csv_path.read_bytes()).hexdigest()


@pytest.mark.parametrize("keys,value", [
    (("schema",), "generic-table"), (("status",), "manufactured"),
    (("source", "doi"), "10.1234/wrong"), (("source", "report"), "Unknown report"),
    (("source", "pdf_sha256"), "0" * 64), (("source", "table_text_sha256"), "0" * 64),
    (("observable", "material"), "AISI Type 316, one supplied solid plate sample"),
    (("observable", "material"), "AISI Type 304, liquid"),
    (("observable", "sample_composition"), "assumed 304 composition"),
    (("observable", "sample_thermal_history"), "assumed annealed"),
    (("observable", "temperature_source_unit"), "K"), (("observable", "temperature_si_unit"), "degC"),
    (("observable", "resistivity_unit"), "ohm cm"), (("observable", "conductivity_unit"), "W cm^-1 K^-1"),
    (("observable", "supported_temperature_range_degC"), [-20, 160]),
    (("observable", "supported_temperature_range_K"), [253.15, 433.15]),
    (("conversion", "resistivity_scale_ohm_m"), "1e-6"),
    (("recovery", "rows"), 9), (("output", "rows"), 9),
    (("output", "adopted_as_runtime_defaults"), True),
    (("output", "columns"), ["temperature_K", "resistivity_ohm_m", "thermal_conductivity_W_m_K"]),
    (("source_error_statements", "resistivity", "plus_minus_percent"), 0.),
    (("source_error_statements", "conductivity", "percent_range"), [0, 0]),
    (("source_error_statements", "evidence", "printed_page"), 13),
    (("source_error_statements", "evidence", "pdf_page_1based"), 17),
    (("source_error_statements", "resistivity", "interpretation"), "one standard deviation"),
    (("source_error_statements", "conductivity", "interpretation"), "one standard deviation"),
])
def test_loader_rejects_incompatible_units_scales_ranges_and_source_context(tmp_path, keys, value):
    csv_path, metadata_path = copy_source(tmp_path)
    metadata = read_metadata(metadata_path)
    set_nested(metadata, keys, value)
    write_metadata(metadata_path, metadata)
    with pytest.raises(ValueError):
        load(csv_path, metadata_path)


@pytest.mark.parametrize("row,column,value", [
    (1, 2, "253.16"), (1, 4, "6.83e-7"), (2, 0, "1"),
    (2, 1, "1"), (2, 2, "274.15"), (3, 3, "7.18"),
    (4, 4, "nan"), (4, 5, "-1"), (5, 5, "inf"),
])
def test_loader_checks_decimal_source_to_si_conversions_after_digest_refresh(tmp_path, row, column, value):
    csv_path, metadata_path = copy_source(tmp_path)
    rows = list(csv.reader(csv_path.read_text().splitlines()))
    rows[row][column] = value
    write_rows(csv_path, rows)
    refresh_digest(csv_path, metadata_path)
    with pytest.raises(ValueError):
        load(csv_path, metadata_path)


def test_named_source_adapter_rejects_coherent_table_changes_even_with_updated_hash(tmp_path):
    csv_path, metadata_path = copy_source(tmp_path)
    rows = list(csv.reader(csv_path.read_text().splitlines()))
    rows[3][3], rows[3][4] = "7.18", "7.18e-7"
    write_rows(csv_path, rows)
    refresh_digest(csv_path, metadata_path)
    with pytest.raises(ValueError):
        load(csv_path, metadata_path)


@pytest.mark.parametrize("change", ["omit_row", "extra_column", "wrong_header", "duplicate_header", "missing_cell"])
def test_loader_checks_shape_and_headers_after_digest_refresh(tmp_path, change):
    csv_path, metadata_path = copy_source(tmp_path)
    rows = list(csv.reader(csv_path.read_text().splitlines()))
    if change == "omit_row":
        rows.pop(3)
    elif change == "extra_column":
        rows[0].append("unrelated")
        for row in rows[1:]:
            row.append("0")
    elif change == "wrong_header":
        rows[0][4] = "resistivity_ohm_cm"
    elif change == "duplicate_header":
        rows[0][4] = rows[0][3]
    else:
        rows[2].pop()
    write_rows(csv_path, rows)
    refresh_digest(csv_path, metadata_path)
    with pytest.raises(ValueError):
        load(csv_path, metadata_path)


def test_loader_rejects_duplicate_metadata_keys_and_nonfinite_json_constants(tmp_path):
    csv_path, metadata_path = copy_source(tmp_path)
    original = metadata_path.read_text()
    duplicate = original.rstrip()[:-1] + ',"status":"recovered_published_smoothed_table"}'
    assert json.loads(duplicate)["status"] == "recovered_published_smoothed_table"
    metadata_path.write_text(duplicate)
    with pytest.raises(ValueError):
        load(csv_path, metadata_path)
    metadata = json.loads(original)
    metadata["source_error_statements"]["resistivity"]["plus_minus_percent"] = float("nan")
    metadata_path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError):
        load(csv_path, metadata_path)
