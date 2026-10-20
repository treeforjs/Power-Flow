"""Published-table guards, source anchors and exact unit/rounding recovery."""

import csv
from decimal import Decimal
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools/recover_flynn_304_transport.py"
CSV = ROOT / "data/materials/flynn_1968_304_transport.csv"
METADATA = ROOT / "data/materials/flynn_1968_304_transport.json"
spec = importlib.util.spec_from_file_location("recover_flynn_304_transport", TOOL)
recovery = importlib.util.module_from_spec(spec)
spec.loader.exec_module(recovery)

# The PDF's actual table text layer, including its visually checked OCR defects.
# No external PDF or optional rendering package is needed for these unit tests.
TABLE_TEXT = (
    'Temperature\n°C\nElectrical Resistivity\nQm\nThermal Conductivity\nW m-1\n°C" 1\n'
    '-20\n6.82 x 10* 7\n13.9\n0\n7.00\n14.3\n20\n7.17\n14.8\n'
    '40\n7.33\n15.2\n60\n7.50\n15.7\n80\n7.66\n16.2\n'
    '100\n7.81\n16.7\n120\n7.97\n17.1\n140\n8.12\n17.7\n150\n8.19\n17-9\n'
)


def published_rows():
    with CSV.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def test_actual_text_fixture_matches_pinned_source_and_parses_ten_rows():
    assert hashlib.sha256(TABLE_TEXT.encode()).hexdigest() == recovery.TABLE_TEXT_SHA256
    cells = recovery.parse_table_text(TABLE_TEXT)
    assert len(cells) == 10
    assert cells[0] == ("-20", "6.82", "13.9")
    assert cells[1] == ("0", "7.00", "14.3")
    assert cells[8] == ("140", "8.12", "17.7")
    assert cells[-1] == ("150", "8.19", "17.9")


@pytest.mark.parametrize("text", [
    TABLE_TEXT.replace("Electrical Resistivity", "Electrical Conductivity"),
    TABLE_TEXT.replace("°C\n", "K\n", 1),
    TABLE_TEXT.replace("0\n7.00\n14.3\n", ""),
    TABLE_TEXT + "160\n8.26\n18.2\n",
    TABLE_TEXT.replace("7.33", "7.330"),
    TABLE_TEXT.replace("7.33", "7.34"),
    TABLE_TEXT.replace("7.81", "nan"),
    TABLE_TEXT.replace("10* 7", "10* 8"),
    TABLE_TEXT.replace("17-9", "17.9"),
])
def test_changed_table_units_count_precision_or_cells_are_rejected(text):
    with pytest.raises(ValueError):
        recovery.parse_table_text(text)


def test_source_identity_guard_also_rejects_text_only_whitespace_change():
    # The parser tolerates harmless line whitespace; source identity remains
    # stricter, so a new extraction/rendering cannot pass without review.
    changed = TABLE_TEXT + "\n"
    assert recovery.parse_table_text(changed) == recovery.parse_table_text(TABLE_TEXT)
    with pytest.raises(ValueError, match="Table text differs"):
        recovery.guard_digest("Table text", changed.encode(), recovery.TABLE_TEXT_SHA256)


def test_an_unreviewed_renderer_cannot_silently_change_source_crop_guards():
    with pytest.raises(ValueError, match="has not been visually inspected"):
        recovery.render_reference("unreviewed")


def test_si_values_and_published_decimal_precision():
    rows = recovery.convert_cells(recovery.parse_table_text(TABLE_TEXT))
    assert rows[0] == {"source_row": "1", "temperature_source_degC": "-20",
        "temperature_K": "253.15", "resistivity_source_mantissa": "6.82",
        "resistivity_ohm_m": "6.82e-7", "thermal_conductivity_W_m_K": "13.9"}
    assert rows[1]["temperature_K"] == "273.15"
    assert rows[1]["resistivity_source_mantissa"] == "7.00"
    assert rows[1]["resistivity_ohm_m"] == "7.00e-7"
    assert rows[2]["temperature_K"] == "293.15"
    assert rows[2]["resistivity_ohm_m"] == "7.17e-7"
    assert rows[-1]["temperature_K"] == "423.15"
    assert rows[-1]["resistivity_ohm_m"] == "8.19e-7"
    assert rows[-1]["thermal_conductivity_W_m_K"] == "17.9"
    assert recovery.csv_bytes(rows) == CSV.read_bytes()


def test_conversion_cannot_be_used_to_add_unpublished_temperatures():
    cells = recovery.parse_table_text(TABLE_TEXT)
    cells.append(("200", "8.50", "19.0"))
    with pytest.raises(ValueError, match="restricted to the ten"):
        recovery.convert_cells(cells)


def test_published_artifact_anchors_and_source_error_semantics():
    metadata = json.loads(METADATA.read_text())
    rows = published_rows()
    assert len(rows) == metadata["output"]["rows"] == 10
    assert metadata["output"]["csv_sha256"] == hashlib.sha256(CSV.read_bytes()).hexdigest()
    assert metadata["status"] == "recovered_published_smoothed_table"
    assert metadata["source"]["pdf_sha256"] == "ff6afb078d512aba2a07b40228a3997857f8bb5448eba465dca9c3d06c74edfa"
    assert metadata["source"]["printed_page"] == 12
    assert metadata["source"]["pdf_page_1based"] == 16
    assert [r["temperature_source_degC"] for r in rows] == ["-20", "0", "20", "40", "60", "80", "100", "120", "140", "150"]
    assert min(Decimal(r["temperature_K"]) for r in rows) == Decimal("253.15")
    assert max(Decimal(r["temperature_K"]) for r in rows) == Decimal("423.15")
    observable = metadata["observable"]
    assert observable["sample_composition"] is None
    assert observable["sample_thermal_history"] is None
    assert observable["resistivity_unit"] == "ohm m"
    assert observable["conductivity_unit"] == "W m^-1 K^-1"
    assert metadata["conversion"]["resistivity_scale_ohm_m"] == "1e-7"
    errors = metadata["source_error_statements"]
    assert errors["resistivity"]["plus_minus_percent"] == .5
    assert errors["conductivity"]["percent_range"] == [2, 3]
    assert errors["coverage_factor"] is None and errors["confidence_level"] is None
    assert errors["standard_deviation_interpretation"] is False
    assert errors["printed_rounding_is_separate_from_source_error"] is True
    assert metadata["output"]["adopted_as_runtime_defaults"] is False
    assert metadata["recovery"]["extrapolation"] is False
    assert metadata["recovery"]["interpolation"] is False


def test_changed_pdf_is_rejected_before_optional_dependencies_or_output(tmp_path):
    pdf = tmp_path / "wrong.pdf"
    pdf.write_bytes(b"%PDF-1.4\nNot the inspected report")
    output = tmp_path / "table.csv"
    metadata = tmp_path / "table.json"
    result = subprocess.run([sys.executable, str(TOOL), "--pdf", str(pdf),
        "--csv", str(output), "--metadata", str(metadata)], capture_output=True, text=True, timeout=30)
    assert result.returncode == 2
    assert "PDF differs from the inspected source" in result.stderr
    assert not output.exists() and not metadata.exists()


def test_source_destination_collision_is_rejected_without_overwrite(tmp_path):
    pdf = tmp_path / "source.pdf"
    raw = b"untouched source"
    pdf.write_bytes(raw)
    result = subprocess.run([sys.executable, str(TOOL), "--pdf", str(pdf),
        "--csv", str(pdf), "--metadata", str(tmp_path / "table.json")],
        capture_output=True, text=True, timeout=30)
    assert result.returncode == 2
    assert "destinations must be different" in result.stderr
    assert pdf.read_bytes() == raw


def test_optional_source_pdf_integration(tmp_path):
    pdf = Path(os.environ.get("FLYNN_304_SOURCE_PDF",
        "/workspace/research_sources/ss304_properties/flynn_1968_nbs9869.pdf"))
    if not pdf.is_file():
        pytest.skip("Validated external source PDF is not available")
    pytest.importorskip("pymupdf", reason="Optional recovery dependency is not installed")
    rows, details = recovery.recover(pdf, tmp_path / "qa")
    assert recovery.csv_bytes(rows) == CSV.read_bytes()
    assert [r["xref"] for r in details["native_rasters"]] == [73, 74, 75]
    reference = recovery.render_reference(details["render_recipe"]["used_PyMuPDF_version"])
    assert details["source_crops"][0]["decoded_RGB_sha256"] == reference[0]
    assert details["source_crops"][1]["decoded_RGB_sha256"] == reference[1]
    for crop in details["source_crops"]:
        path = Path(crop["written_path"])
        assert path.is_file()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == crop["encoded_PNG_sha256"]
