#!/usr/bin/env python3
"""Recover Flynn's ten-row, published smoothed AISI 304 transport table.

Source: NBS Report 9869 (1968), printed p. 12 / PDF p. 16. The inspected
scan, page content, native rasters, table text and rendered crop are pinned.
The text layer has known OCR errors: the shared 10^-7 resistivity scale,
the ohm symbol, and the last conductivity cell require visual confirmation.
Only those documented corrections are accepted. No interpolation, fit,
extrapolation or runtime material model is supplied.

Recovery dependency: PyMuPDF (inspected renderings: 1.26.6 and 1.28.2).
Parsing and conversion helpers use only the Python standard library.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal
import hashlib
import io
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
PDF_SHA256 = "ff6afb078d512aba2a07b40228a3997857f8bb5448eba465dca9c3d06c74edfa"
PAGE_INDEX = 15
PAGE_CONTENT_SHA256 = "a34325dc2bc302f1a841b75b46861cba035d9d05e84750ea861d23a712f4931c"
TABLE_CLIP_PT = (90, 275, 545, 450)
CONTEXT_CLIP_PT = (65, 95, 545, 455)
TABLE_TEXT_SHA256 = "249300bc3e3ccd2509b2f3c55983f080c3a061aa0bf0f03013952e1e11544167"
TABLE_RGB_SHA256 = "b75599a7aca7e4a541c0af0653eccd7d5451a68a9d71c20040da5674675f5c4e"
CONTEXT_RGB_SHA256 = "9d01220aa958479ecec989a81ec21ab7b3a492a308953ba882f3f06d5d957118"
INSPECTED_RENDERERS = {
    "1.26.6": (TABLE_RGB_SHA256, CONTEXT_RGB_SHA256),
    "1.28.2": ("e7dac76a1f17f89152b43f7df069590c488ac092f51b24c133c8069a8075d18f",
               "c01c85ca34088a264e686b11aa7707a45b8b363ba282364ec12cb2928df43b2d"),
}
NATIVE_RASTERS = (
    (73, 1044, 1487, 3, "6e9bd8d0e84f053713d3258d864741cf3cfa44ec874964bce69b99e0fcd92eec"),
    (74, 3132, 4461, 3, "603cec50d22dd8257e5896738f33c675985e1b009af41ee236a6cc7e925e871b"),
    (75, 3132, 4461, 1, "a3a169a602e3a76daa45aa0cd93464d1ef9025ff8d7536b409c7f0b7afe135f2"),
)
HEADER = ("Temperature", "°C", "Electrical Resistivity", "Qm",
          "Thermal Conductivity", "W m-1", '°C" 1')
# Every cell below was checked against the actual page image. Resistivity's
# 10^-7 ohm m scale is printed once in the first row and applies to the column.
INSPECTED_CELLS = (
    ("-20", "6.82", "13.9"), ("0", "7.00", "14.3"),
    ("20", "7.17", "14.8"), ("40", "7.33", "15.2"),
    ("60", "7.50", "15.7"), ("80", "7.66", "16.2"),
    ("100", "7.81", "16.7"), ("120", "7.97", "17.1"),
    ("140", "8.12", "17.7"), ("150", "8.19", "17.9"),
)
FIELDS = ("source_row", "temperature_source_degC", "temperature_K",
          "resistivity_source_mantissa", "resistivity_ohm_m",
          "thermal_conductivity_W_m_K")


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def guard_digest(label: str, data: bytes, expected: str) -> None:
    if sha256(data) != expected:
        raise ValueError(f"{label} differs from the inspected source; inspect before changing guards")


def render_reference(version: str) -> tuple[str, str]:
    try:
        return INSPECTED_RENDERERS[version]
    except KeyError as error:
        raise ValueError(f"PyMuPDF renderer {version} has not been visually inspected") from error


def parse_table_text(text: str) -> list[tuple[str, str, str]]:
    """Parse the inspected OCR layout, accepting only documented visual repairs."""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if tuple(lines[:len(HEADER)]) != HEADER:
        raise ValueError("Table headings or OCR unit labels differ from the inspected page")
    cells = lines[len(HEADER):]
    if len(cells) != 30:
        raise ValueError("Expected exactly ten source rows with three cells each")
    if cells[1] != "6.82 x 10* 7" or cells[-1] != "17-9":
        raise ValueError("Expected the documented first-row scale and final-cell OCR forms")
    cells[1] = "6.82"
    cells[-1] = "17.9"
    parsed = []
    for offset in range(0, len(cells), 3):
        temperature, resistivity, conductivity = cells[offset:offset + 3]
        if not (re.fullmatch(r"-?\d+", temperature)
                and re.fullmatch(r"\d+\.\d{2}", resistivity)
                and re.fullmatch(r"\d+\.\d", conductivity)):
            raise ValueError("Unexpected numeric cell or changed printed precision")
        parsed.append((temperature, resistivity, conductivity))
    if tuple(parsed) != INSPECTED_CELLS:
        raise ValueError("A table cell differs from the visually verified transcription")
    return parsed


def convert_cells(cells: list[tuple[str, str, str]]) -> list[dict[str, str]]:
    """Convert exactly the verified cells; retain all original decimal precision."""
    if tuple(cells) != INSPECTED_CELLS:
        raise ValueError("Conversion is restricted to the ten inspected published rows")
    rows = []
    for number, (temperature, resistivity, conductivity) in enumerate(cells, 1):
        rows.append({
            "source_row": str(number),
            "temperature_source_degC": temperature,
            "temperature_K": format(Decimal(temperature) + Decimal("273.15"), ".2f"),
            "resistivity_source_mantissa": resistivity,
            "resistivity_ohm_m": format(Decimal(resistivity) * Decimal("1e-7"), ".2e"),
            "thermal_conductivity_W_m_K": conductivity,
        })
    return rows


def csv_bytes(rows: list[dict[str, str]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def recover(pdf: Path, qa_dir: Path | None = None) -> tuple[list[dict[str, str]], dict]:
    guard_digest("PDF", pdf.read_bytes(), PDF_SHA256)
    try:
        import pymupdf
    except ImportError as error:
        raise RuntimeError("PDF recovery needs PyMuPDF; install the recovery dependency") from error

    with pymupdf.open(pdf) as document:
        if len(document) != 20:
            raise ValueError("Expected the inspected 20-page NBS report")
        page = document[PAGE_INDEX]
        guard_digest("Page content", page.read_contents(), PAGE_CONTENT_SHA256)
        if tuple(page.rect) != (0.0, 0.0, 563.0, 802.0):
            raise ValueError("Source page dimensions changed")
        images = page.get_images(full=True)
        if [(image[0], image[1]) for image in images] != [(73, 0), (74, 75)]:
            raise ValueError("Expected the two inspected page rasters and their mask")
        text = page.get_text("text", clip=pymupdf.Rect(TABLE_CLIP_PT), sort=False)
        guard_digest("Table text", text.encode("utf-8"), TABLE_TEXT_SHA256)
        rows = convert_cells(parse_table_text(text))
        qa_records = []
        table_digest, context_digest = render_reference(pymupdf.VersionBind)
        for name, clip, shape, digest in (
            ("table_crop", TABLE_CLIP_PT, (910, 350, 3), table_digest),
            ("uncertainty_and_table", CONTEXT_CLIP_PT, (960, 720, 3), context_digest),
        ):
            pixmap = page.get_pixmap(matrix=pymupdf.Matrix(2, 2), clip=pymupdf.Rect(clip),
                                     colorspace=pymupdf.csRGB, alpha=False)
            if (pixmap.width, pixmap.height, pixmap.n) != shape:
                raise ValueError("Rendered source crop dimensions changed")
            guard_digest(name, pixmap.samples, digest)
            record = {"name": name, "page_clip_pt": list(clip), "pixel_shape_width_height_channels": list(shape),
                      "decoded_RGB_sha256": digest, "written_path": None, "encoded_PNG_sha256": None}
            if qa_dir is not None:
                qa_dir.mkdir(parents=True, exist_ok=True)
                destination = qa_dir / f"flynn_1968_304_{name}.png"
                raw = pixmap.tobytes("png")
                destination.write_bytes(raw)
                record.update(written_path=str(destination.resolve()), encoded_PNG_sha256=sha256(raw))
            qa_records.append(record)

        # Render before separately decoding native Pixmaps. In MuPDF 1.28.2,
        # decoding the standalone soft mask first changes later compositing
        # through the document image cache. Both fresh renderings were inspected;
        # the native image and mask identity checks still all run before return.
        native_records = []
        for xref, width, height, channels, digest in NATIVE_RASTERS:
            pixmap = pymupdf.Pixmap(document, xref)
            if (pixmap.width, pixmap.height, pixmap.n) != (width, height, channels):
                raise ValueError("Native raster dimensions or channels changed")
            guard_digest(f"Native raster {xref}", pixmap.samples, digest)
            native_records.append({"xref": xref, "width_px": width, "height_px": height,
                                   "channels": channels, "decoded_samples_sha256": digest})

    details = {
        "native_rasters": native_records,
        "render_recipe": {"matrix_scale": 2, "colorspace": "RGB", "alpha": False,
                          "inspected_PyMuPDF_versions": list(INSPECTED_RENDERERS),
                          "used_PyMuPDF_version": pymupdf.VersionBind,
                          "operation_order": "Render source crops before separately decoding native image and mask Pixmaps; avoid mask-cache-dependent compositing",
                          "version_effect": "These inspected renderers differ in composited RGB crop pixels while source PDF, native decoded rasters and extracted table text match; each rendering has a separate visually checked digest"},
        "source_crops": qa_records,
    }
    return rows, details


def metadata_for(rows: list[dict[str, str]], csv_path: Path, raw_csv: bytes, details: dict) -> dict:
    return {
        "schema": "mhdlab-published-smoothed-transport-table-v1",
        "status": "recovered_published_smoothed_table",
        "source": {"author": "D. R. Flynn", "report_date": "1968-06-25", "report": "NBS Report 9869",
            "title": "Thermal conductivity and electrical resistivity of a sample of AISI type 304 stainless steel",
            "doi": "10.6028/NBS.RPT.9869",
            "url": "https://nvlpubs.nist.gov/nistpubs/Legacy/RPT/nbsreport9869.pdf",
            "pdf_sha256": PDF_SHA256, "pdf_pages": 20, "printed_page": 12, "pdf_page_1based": 16,
            "page_content_sha256": PAGE_CONTENT_SHA256,
            "table_text_sha256": TABLE_TEXT_SHA256, "table_text_clip_pt": list(TABLE_CLIP_PT),
            "table_text_method": "PyMuPDF get_text('text', clip=Rect(...), sort=False), UTF-8 exact bytes",
            "permission_context": "Official public archive; PDF p. 3 bears public-release approval dated October 9, 2015; no new license assigned."},
        "recovery": {"method": "Guarded source-text parsing with visually verified cell transcription and an explicit OCR correction whitelist",
            "visual_inspection_date": "2026-10-18", "rows": len(rows),
            "known_OCR_repairs": [
                {"source_text": "Qm", "visually_verified_print": "ohm m", "role": "Resistivity unit label"},
                {"source_text": "6.82 x 10* 7", "visually_verified_print": "6.82 x 10^-7", "role": "First row establishes shared column scale; following mantissas retain it"},
                {"source_text": "17-9", "visually_verified_print": "17.9", "role": "150 degC conductivity cell"},
                {"source_text": 'W m-1 / °C" 1', "visually_verified_print": "W m^-1 degC^-1", "role": "Conductivity temperature-gradient unit"}],
            "missing_rows_invented": False, "graph_digitization": False, "interpolation": False,
            "extrapolation": False, "fit": False},
        "observable": {"material": "AISI Type 304, one supplied solid plate sample",
            "value_kind": "Published smoothed results from measurements, not original acquisition samples",
            "sample_composition": None, "sample_thermal_history": None,
            "sample_context_evidence": {"printed_page": 1, "pdf_page_1based": 5},
            "actual_dissertation_foil_match": "Unknown; source sample is not the dissertation foil",
            "temperature_source_unit": "degC", "temperature_si_unit": "K",
            "supported_temperature_range_degC": [-20, 150],
            "supported_temperature_range_K": [253.15, 423.15],
            "resistivity_unit": "ohm m", "conductivity_unit": "W m^-1 K^-1"},
        "conversion": {"temperature": "T_K = T_degC + 273.15, decimal arithmetic",
            "resistivity": "Printed mantissa x 10^-7 ohm m; shared scale visually established in the first row",
            "resistivity_scale_ohm_m": "1e-7",
            "conductivity": "One degree Celsius temperature difference equals one kelvin; numerical conductivity unchanged",
            "rounding": "Integer source temperatures, two-decimal resistivity mantissas and one-decimal conductivities preserved; converted K has exact .15 offset"},
        "source_error_statements": {
            "evidence": {"printed_page": 12, "pdf_page_1based": 16},
            "resistivity": {"plus_minus_percent": 0.5, "interpretation": "Author believes error no greater than this value"},
            "conductivity": {"percent_range": [2, 3], "interpretation": "Author estimates uncertainty no greater than two to three percent; chiefly uncorrected insulation heat loss"},
            "coverage_factor": None, "confidence_level": None,
            "standard_deviation_interpretation": False,
            "between_heat_or_foil_variability_quantified": False,
            "printed_rounding_is_separate_from_source_error": True},
        "validation_limits": ["No values outside 253.15-423.15 K are supported by this table",
            "Unknown chemical analysis and thermal history; no 304L/316 substitution",
            "No high-temperature phase, enthalpy, heat-capacity or pulsed-heating calibration",
            "Not a runtime property model or a material default; no interpolation function supplied"],
        "source_QA": details,
        "output": {"csv_filename": csv_path.name, "csv_sha256": sha256(raw_csv),
            "rows": len(rows), "columns": list(FIELDS), "adopted_as_runtime_defaults": False},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pdf", type=Path, required=True, help="Validated NBS Report 9869 PDF")
    parser.add_argument("--csv", type=Path, default=ROOT / "data/materials/flynn_1968_304_transport.csv")
    parser.add_argument("--metadata", type=Path, default=ROOT / "data/materials/flynn_1968_304_transport.json")
    parser.add_argument("--qa-dir", type=Path, help="Optional local directory for inspected source crops")
    args = parser.parse_args()
    try:
        if args.csv.resolve() == args.metadata.resolve() or args.pdf.resolve() in (
                args.csv.resolve(), args.metadata.resolve()):
            raise ValueError("Source PDF, CSV and metadata destinations must be different")
        rows, details = recover(args.pdf, args.qa_dir)
        raw_csv = csv_bytes(rows)
        metadata = metadata_for(rows, args.csv, raw_csv, details)
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        args.metadata.parent.mkdir(parents=True, exist_ok=True)
        args.csv.write_bytes(raw_csv)
        args.metadata.write_text(json.dumps(metadata, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(2, f"Recovery failed: {error}\n")
    print(json.dumps({"rows": len(rows), "csv": str(args.csv), "csv_sha256": sha256(raw_csv),
                      "metadata": str(args.metadata)}, sort_keys=True))


if __name__ == "__main__":
    main()
