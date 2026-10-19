#!/usr/bin/env python3
"""Recover the finite plotted drive insets of dissertation Figures 4.28/4.29.

These are reconstructed published legacy Allegra simulation drives. Their display
clock and literal current coordinates are not an experimental clock registration
or a recovered original simulation input array.
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path
from typing import Any

PDF_SHA256 = "1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7"
FIGURES = (
    {"figure": "4.28", "page_index": 118, "printed_page": 94, "xref": 2847, "snapshot_ns": 40,
     "rgb_sha256": "458b3ced8db4501db3fcbcc90939ce5f336a7ec83eb97ee259c66f7112bf07c8",
     "reference_png_sha256": "1f9c939c8b72947da2629ccad92c42dda04503719c3cdb89fe0e921bc5ce5488",
     "page_content_sha256": "934979de84040dc7970be53747773d2f8c78fedad9ac2b7116757a8052449103",
     "expected_retained_columns": 178},
    {"figure": "4.29", "page_index": 120, "printed_page": 96, "xref": 2870, "snapshot_ns": 60,
     "rgb_sha256": "ac25b035329213a8d0a354b8496b1674c2be24ceaeff65d95a48ef2e6376b34f",
     "reference_png_sha256": "5ae90775c30a8fad657f74d5813c150b2a4284d1a0fd8075ebc4fcd21b88e89a",
     "page_content_sha256": "17080e390125995d993f2fa580e2fcef6cf4cb166a267bf08a9e1a9798a8d43e",
     "expected_retained_columns": 262},
)
SOURCE_RGB = (178, 34, 34)
TIME_TICK_PIXELS = (100.5, 184.5, 268.5, 352.5, 436.5)
TIME_TICK_VALUES = (0.0, 20.0, 40.0, 60.0, 80.0)
CURRENT_TICK_PIXELS = (379.5, 313.5, 247.5, 182.5, 116.5, 50.5)
CURRENT_TICK_VALUES = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
PANEL_X_MIN, PANEL_X_MAX = 90, 477
PANEL_Y_MIN, PANEL_Y_MAX = 19, 411
READOUT_PIXELS = 2.0
FIELDS = (
    "source_identifier", "figure", "snapshot_time_s", "pixel_x", "pixel_y_median",
    "run_first_pixel_y", "run_last_pixel_y", "matched_color_pixels",
    "time_ns", "time_readout_lower_ns", "time_readout_upper_ns",
    "current_MA", "current_readout_lower_MA", "current_readout_upper_MA",
    "time_s", "time_readout_lower_s", "time_readout_upper_s",
    "current_A", "current_readout_lower_A", "current_readout_upper_A",
    "left_display_clip_column", "thin_run_endpoint_candidate",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def split_runs(indices: list[int]) -> list[list[int]]:
    runs: list[list[int]] = []
    for value in indices:
        if not runs or value != runs[-1][-1]+1:
            runs.append([value])
        else:
            runs[-1].append(value)
    return [[run[0], run[-1]] for run in runs]


def select_column(pixels: Any, x: int) -> tuple[dict[str, Any] | None, str, list[list[int]]]:
    """Recover one exact-colour run without endpoint padding or gap filling."""
    matches = [y for y in range(PANEL_Y_MIN, PANEL_Y_MAX+1)
               if tuple(int(c) for c in pixels[y, x]) == SOURCE_RGB]
    runs = split_runs(matches)
    if len(runs) != 1:
        return None, "no_exact_colour_run" if not runs else "multiple_exact_colour_runs", runs
    first, last = runs[0]
    if not 2 <= len(matches) <= 8:
        return None, "run_width_outside_2_to_8_pixels", runs
    if first <= PANEL_Y_MIN or last >= PANEL_Y_MAX:
        return None, "run_touches_vertical_panel_ROI_boundary", runs
    return {"pixel_x": x, "pixel_y_median": (first+last)/2,
            "run_first_pixel_y": first, "run_last_pixel_y": last,
            "matched_color_pixels": len(matches),
            "left_display_clip_column": int(x == PANEL_X_MIN),
            "thin_run_endpoint_candidate": int(len(matches) == 2)}, "retained", runs


def linear_calibration(pixels: tuple[float, ...], values: tuple[float, ...]) -> dict[str, Any]:
    x = [Fraction(str(value)) for value in pixels]
    y = [Fraction(str(value)) for value in values]
    n = len(x)
    slope = (n*sum(a*b for a, b in zip(x, y))-sum(x)*sum(y))/(n*sum(a*a for a in x)-sum(x)**2)
    intercept = (sum(y)-slope*sum(x))/n
    residuals = [b-(slope*a+intercept) for a, b in zip(x, y)]
    return {"tick_pixel_centres": list(pixels), "tick_values": list(values),
            "slope": float(slope), "intercept": float(intercept),
            "slope_exact": [str(slope.numerator), str(slope.denominator)],
            "intercept_exact": [str(intercept.numerator), str(intercept.denominator)],
            "tick_residuals": [float(r) for r in residuals],
            "max_abs_tick_residual": float(max(abs(r) for r in residuals)),
            "fit": "ordinary_least_squares_exact_rational_tick_coordinates"}


def validate_ticks(pixels: Any) -> None:
    time = split_runs([x for x in range(90, 478) if max(int(c) for c in pixels[417, x]) < 128])
    current = split_runs([y for y in range(19, 412) if max(int(c) for c in pixels[y, 84]) < 128])
    if [(a+b)/2 for a, b in time] != list(TIME_TICK_PIXELS):
        raise ValueError("native time tick anchors do not match configured calibration")
    if [(a+b)/2 for a, b in reversed(current)] != list(CURRENT_TICK_PIXELS):
        raise ValueError("native current tick anchors do not match configured calibration")


def transform(calibration: dict[str, Any], pixel: float) -> float:
    return calibration["slope"]*pixel+calibration["intercept"]


def make_rows(pixels: Any, figure: dict[str, Any], calibration: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows, retained, excluded = [], [], []
    time, current = calibration["time"], calibration["current"]
    identifier = f"figure_{figure['figure'].replace('.', '_')}_simulation_drive"
    for x in range(PANEL_X_MIN, PANEL_X_MAX+1):
        selected, reason, runs = select_column(pixels, x)
        raw = {"pixel_x": x, "exact_colour_runs_inclusive_y": runs, "reason": reason}
        if selected is None:
            excluded.append(raw)
            continue
        t = transform(time, x)
        tlo = transform(time, x-READOUT_PIXELS)-time["max_abs_tick_residual"]
        thi = transform(time, x+READOUT_PIXELS)+time["max_abs_tick_residual"]
        value = transform(current, selected["pixel_y_median"])
        vlo = transform(current, selected["run_last_pixel_y"]+READOUT_PIXELS)-current["max_abs_tick_residual"]
        vhi = transform(current, selected["run_first_pixel_y"]-READOUT_PIXELS)+current["max_abs_tick_residual"]
        rows.append({"source_identifier": identifier, "figure": figure["figure"],
                     "snapshot_time_s": figure["snapshot_ns"]*1e-9, **selected,
                     "time_ns": t, "time_readout_lower_ns": tlo, "time_readout_upper_ns": thi,
                     "current_MA": value, "current_readout_lower_MA": vlo, "current_readout_upper_MA": vhi,
                     "time_s": t*1e-9, "time_readout_lower_s": tlo*1e-9, "time_readout_upper_s": thi*1e-9,
                     "current_A": value*1e6, "current_readout_lower_A": vlo*1e6,
                     "current_readout_upper_A": vhi*1e6})
        retained.append(raw)
    return rows, {"retained_columns": retained, "excluded_columns": excluded,
                  "retained_count": len(retained), "excluded_count": len(excluded),
                  "catalogued_column_count": len(retained)+len(excluded)}


def recover(pdf: Path) -> list[dict[str, Any]]:
    """Read only the pinned native source; no output writes or other producers."""
    if sha256(pdf.read_bytes()) != PDF_SHA256:
        raise ValueError("source PDF SHA256 mismatch")
    import numpy as np
    import pymupdf
    from PIL import Image

    calibration = {"time": linear_calibration(TIME_TICK_PIXELS, TIME_TICK_VALUES),
                   "current": linear_calibration(CURRENT_TICK_PIXELS, CURRENT_TICK_VALUES)}
    results = []
    with pymupdf.open(pdf) as document:
        if len(document) != 203:
            raise ValueError("source PDF page count mismatch")
        for figure in FIGURES:
            page = document[figure["page_index"]]
            if sha256(page.read_contents()) != figure["page_content_sha256"]:
                raise ValueError("source figure page content SHA256 mismatch")
            if figure["xref"] not in {entry[0] for entry in page.get_images(full=True)}:
                raise ValueError("source image xref is not associated with the figure page")
            encoded = document.extract_image(figure["xref"])["image"]
            image = Image.open(io.BytesIO(encoded))
            if image.mode != "RGB" or image.size != (975, 975):
                raise ValueError("native source image must be 975x975 RGB")
            if sha256(image.tobytes()) != figure["rgb_sha256"]:
                raise ValueError("native source RGB SHA256 mismatch")
            pixels = np.asarray(image)
            validate_ticks(pixels)
            rows, catalog = make_rows(pixels, figure, calibration)
            if len(rows) != figure["expected_retained_columns"]:
                raise ValueError("retained native drive support changed")
            results.append({"figure": figure, "rows": rows, "catalog": catalog, "image": image.copy(),
                            "calibration": calibration, "encoded_png_sha256": sha256(encoded),
                            "parser_version": pymupdf.VersionBind, "page_text": page.get_text()})
    return results


def make_qa(result: dict[str, Any], path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)
    rows = result["rows"]
    axes[0].imshow(result["image"])
    axes[0].set_xlim(20, 495)
    axes[0].set_ylim(480, 0)
    axes[0].plot([r["pixel_x"] for r in rows], [r["pixel_y_median"] for r in rows], ".", color="navy", markersize=1.4)
    axes[0].set_xlabel("Native source pixel x")
    axes[0].set_ylabel("Native source pixel y")
    axes[0].set_title("Native current-drive inset and retained exact-colour runs")
    times = [r["time_ns"] for r in rows]
    axes[1].plot(times, [r["current_MA"] for r in rows], color="firebrick", label="recovered plotted drive")
    axes[1].fill_between(times, [r["current_readout_lower_MA"] for r in rows],
                         [r["current_readout_upper_MA"] for r in rows], color="firebrick", alpha=.2,
                         label="full-run + ±2 px + tick-fit current bounds")
    axes[1].errorbar(times[::10], [r["current_MA"] for r in rows][::10],
                     xerr=READOUT_PIXELS*abs(result["calibration"]["time"]["slope"]),
                     fmt="none", color="navy", alpha=.4, capsize=1, label="±2 px time readout (every tenth sample)")
    axes[1].set_xlabel("Published simulation display Time [ns]; clock unregistered")
    axes[1].set_ylabel("Literal plotted Current [MA]")
    axes[1].grid(alpha=.2)
    axes[1].legend(fontsize=8)
    figure.suptitle(f"Dissertation Figure {result['figure']['figure']}: legacy Allegra simulation "
                   f"drive inset at {result['figure']['snapshot_ns']} ns snapshot\nFinite raster support; no endpoint padding or time registration")
    figure.savefig(path.with_suffix(".png"), dpi=160)
    figure.savefig(path.with_suffix(".pdf"))
    plt.close(figure)


def metadata(result: dict[str, Any], csv_path: Path, qa_paths: list[Path], geometry: Path) -> dict[str, Any]:
    figure, rows = result["figure"], result["rows"]
    identifier = rows[0]["source_identifier"]
    return {
        "schema": "mhdlab-dissertation-simulation-drive-recovery-v1", "source_identifier": identifier,
        "source": {"title": "Characterizing Desorbed Surface Contaminants Leading to Plasma Formation in Power Feeds of Pulsed Power Accelerators",
                   "item_url": "https://deepblue.lib.umich.edu/items/9ccfa886-0b19-4bd8-b3b2-bb685b9ddefc",
                   "pdf_sha256": PDF_SHA256, "figure": figure["figure"], "printed_page": figure["printed_page"],
                   "pdf_page": figure["page_index"]+1, "pdf_page_index": figure["page_index"],
                   "page_content_sha256": figure["page_content_sha256"], "page_text": result["page_text"],
                   "embedded_image_xref": figure["xref"], "native_image_size_pixels": [975, 975], "native_image_mode": "RGB",
                   "reference_extracted_PNG_sha256": figure["reference_png_sha256"],
                   "extracted_PNG_sha256": result["encoded_png_sha256"], "decoded_RGB_pixels_sha256": figure["rgb_sha256"],
                   "hash_validation": "PDF, page content and decoded RGB mandatory; lossless PNG encoding may vary",
                   "parser": "PyMuPDF", "parser_version": result["parser_version"]},
        "simulation_context": {"code_label_in_dissertation": "Allegra", "drive_shot_caption_label": 17133,
                               "caption_snapshot_time_s": figure["snapshot_ns"]*1e-9,
                               "observation_type": "reconstructed published legacy simulation drive inset",
                               "experimental_measurement": False, "FLASH_result": False,
                               "original_simulation_input_array_recovered": False,
                               "source_display_clock": "literal Time [ns] axis in the published simulation inset",
                               "experimental_clock_alignment_established": False,
                               "time_shift_or_scale_fit_applied": False, "current_repair_or_normalization_applied": False,
                               "VUV_gate_recovered": False,
                               "snapshot_is_observed_waveform_end": False},
        "unit_interpretation": {"horizontal_axis_label": "Time [ns]", "vertical_axis_label": "Current [MA]",
                                "time_scale_to_s": 1e-9, "current_scale_to_A": 1e6,
                                "value_role": "literal plotted axis coordinates; not a validated original simulation input",
                                "original_simulation_input_normalization_verified": False,
                                "negative_display_times_preserved": True,
                                "readout_bounds_clipped_to_nonnegative_current": False},
        "calibration": {**result["calibration"], "time_unit": "ns", "current_unit": "MA",
                        "tick_detection": {"dark_channel_threshold_exclusive": 128, "time_tick_tail_row": 417,
                                           "current_tick_tail_column": 84,
                                           "pixel_coordinate_convention": "native zero-based integer pixel centres"},
                        "time_fit_input_sha256": sha256(json.dumps([TIME_TICK_PIXELS, TIME_TICK_VALUES]).encode()),
                        "current_fit_input_sha256": sha256(json.dumps([CURRENT_TICK_PIXELS, CURRENT_TICK_VALUES]).encode())},
        "selection": {"panel_ROI_inclusive": {"x": [PANEL_X_MIN, PANEL_X_MAX], "y": [PANEL_Y_MIN, PANEL_Y_MAX]},
                      "native_frame_strokes": {"left_x": [88, 89], "right_x": [478, 479], "bottom_y": [412, 413], "top_y": [17, 18]},
                      "source_RGB": list(SOURCE_RGB), "column_stride": 1,
                      "rule": "one contiguous exact RGB run of 2–8 pixels within the vertical ROI; single-pixel caps and absent/ambiguous columns excluded",
                      **result["catalog"], "left_display_support_clipped": True,
                      "left_clip_qualification": "curve reaches the first native interior column; its beginning outside the displayed frame is unknown",
                      "right_support_qualification": "finite visible snapshot trace; final retained column is not a recovered pulse tail or an inserted snapshot endpoint",
                      "interpolation": False, "smoothing": False, "extrapolation": False,
                      "endpoint_padding": False, "zero_time_or_current_inserted": False},
        "readout_bounds": {"correlated_display_bounds_only": True, "independent_noise_or_confidence_interval": False,
                           "additional_pixel_readout_half_width": READOUT_PIXELS,
                           "time_recipe": "native column +/-2 pixels, mapped by fitted time axis and expanded by maximum absolute time tick-fit residual",
                           "current_recipe": "entire matched run [first_y,last_y], extended by +/-2 pixels, mapped by fitted current axis and expanded by maximum absolute current tick-fit residual",
                           "bounds_are_observed_source_support": False,
                           "physical_or_simulation_uncertainty_quantified": False},
        "output": {"csv_filename": csv_path.name, "csv_sha256": sha256(csv_path.read_bytes()), "columns": list(FIELDS),
                   "samples": len(rows), "line_ending": "LF", "time_support_ns": [rows[0]["time_ns"], rows[-1]["time_ns"]],
                   "time_support_s": [rows[0]["time_s"], rows[-1]["time_s"]],
                   "current_extent_MA": [min(r["current_MA"] for r in rows), max(r["current_MA"] for r in rows)],
                   "generator": "tools/recover_dissertation_simulation_drive.py", "generator_sha256": sha256(Path(__file__).read_bytes()),
                   "source_geometry_filename": geometry.name, "source_geometry_sha256": sha256(geometry.read_bytes()),
                   "QA_files_sha256": {p.name: sha256(p.read_bytes()) for p in qa_paths}},
        "limits": ["Published legacy Allegra simulation raster; distinct from the measured Figure C.19 shot-current trace.",
                   "The displayed time origin, trigger relation and experimental clock registration are unknown; negative plotted times are retained.",
                   "Caption snapshot times label the figures and do not determine an observed waveform tail or justify endpoint padding.",
                   "The left trace is display-clipped; no missing beginning, zero crossing, pulse onset or continuation is reconstructed.",
                   "Literal ns/MA conversions preserve plotted coordinates; original simulation input arrays and normalization remain unavailable.",
                   "Samples and full-run/readout envelopes are correlated raster observations, not independent experimental noise.",
                   "No current renormalization, time fit, VUV gate, magnetic/current closure, physical heating input or FLASH result is inferred."],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--qa-dir", required=True, type=Path)
    args = parser.parse_args()
    targets = [args.output_dir/f"figure_{f['figure'].replace('.', '_')}_simulation_drive.{suffix}"
               for f in FIGURES for suffix in ("csv", "json")]
    if any(path.exists() for path in targets):
        raise ValueError("refusing to overwrite existing recovered datasets")
    results = recover(args.pdf)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.qa_dir.mkdir(parents=True, exist_ok=True)
    summary = []
    for result in results:
        stem = f"figure_{result['figure']['figure'].replace('.', '_')}_simulation_drive"
        csv_path = args.output_dir/f"{stem}.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
            writer.writeheader()
            writer.writerows({key: format(value, ".12g") if isinstance(value, float) else value for key, value in row.items()} for row in result["rows"])
        qa = args.qa_dir/stem
        make_qa(result, qa)
        geometry = args.qa_dir/f"{stem}_source_geometry.json"
        geometry.write_text(json.dumps({"schema": "mhdlab-simulation-drive-source-geometry-v1", "figure": result["figure"],
                                        "source_RGB": list(SOURCE_RGB), "calibration": result["calibration"],
                                        "selection": result["catalog"]}, indent=2, allow_nan=False)+"\n")
        document = metadata(result, csv_path, [qa.with_suffix(".png"), qa.with_suffix(".pdf")], geometry)
        metadata_path = args.output_dir/f"{stem}.json"
        metadata_path.write_text(json.dumps(document, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
        summary.append({"figure": result["figure"]["figure"], "samples": len(result["rows"]),
                        "csv": str(csv_path), "csv_sha256": sha256(csv_path.read_bytes()),
                        "metadata": str(metadata_path), "metadata_sha256": sha256(metadata_path.read_bytes()),
                        "time_support_ns": document["output"]["time_support_ns"],
                        "current_extent_MA": document["output"]["current_extent_MA"],
                        "QA_files_sha256": document["output"]["QA_files_sha256"]})
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
