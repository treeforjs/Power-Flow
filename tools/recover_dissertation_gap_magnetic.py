#!/usr/bin/env python3
"""Recover conservative signed B_x central-gap raster profiles, Figs. 4.28/4.29.

These are published legacy Allegra simulation displays. The primary quantity is
the literal signed B_x [100 Tesla] axis value. Conditional tesla columns do not
verify original electromagnetic normalization or supply physical field inputs.
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
    {"figure": "4.28", "page_index": 118, "printed_page": 94, "xref": 2847, "time_s": 40e-9,
     "rgb_sha256": "458b3ced8db4501db3fcbcc90939ce5f336a7ec83eb97ee259c66f7112bf07c8",
     "reference_png_sha256": "1f9c939c8b72947da2629ccad92c42dda04503719c3cdb89fe0e921bc5ce5488",
     "page_content_sha256": "934979de84040dc7970be53747773d2f8c78fedad9ac2b7116757a8052449103"},
    {"figure": "4.29", "page_index": 120, "printed_page": 96, "xref": 2870, "time_s": 60e-9,
     "rgb_sha256": "ac25b035329213a8d0a354b8496b1674c2be24ceaeff65d95a48ef2e6376b34f",
     "reference_png_sha256": "5ae90775c30a8fad657f74d5813c150b2a4284d1a0fd8075ebc4fcd21b88e89a",
     "page_content_sha256": "17080e390125995d993f2fa580e2fcef6cf4cb166a267bf08a9e1a9798a8d43e"},
)
PROFILES = (("foil_edge", (255, 105, 180)), ("foil_middle", (186, 85, 211)))
MAGNETIC_TICK_PIXELS = (606.5, 706.5, 806.5, 906.5)
MAGNETIC_TICK_VALUES = (-1.0, -0.5, 0.0, 0.5)
GAP_TICK_PIXELS = (56.5, 135.5, 214.5, 293.5, 372.5)
GAP_TICK_VALUES = (4.0, 2.0, 0.0, -2.0, -4.0)
PANEL_X_MIN, PANEL_X_MAX = 570, 957
PANEL_ROW_MIN, PANEL_ROW_MAX = 19, 411
INTERIOR_ROW_MIN, INTERIOR_ROW_MAX = 150, 280
READOUT_PIXELS = 2
CONDITIONAL_TESLA_SCALE = 100
FIELDS = (
    "figure", "profile", "time_s", "pixel_y", "pixel_x_median", "run_first_pixel_x",
    "run_last_pixel_x", "matched_color_pixels", "gap_coordinate_mm",
    "gap_coordinate_readout_lower_mm", "gap_coordinate_readout_upper_mm",
    "Bx_literal_axis_value", "Bx_literal_readout_lower_axis_value", "Bx_literal_readout_upper_axis_value",
    "Bx_literal_conditional_T", "Bx_literal_conditional_readout_lower_T", "Bx_literal_conditional_readout_upper_T",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def exact_pair(value: Fraction) -> list[str]:
    return [str(value.numerator), str(value.denominator)]


def split_runs(indices: list[int]) -> list[list[int]]:
    runs: list[list[int]] = []
    for value in indices:
        if not runs or value != runs[-1][-1]+1:
            runs.append([value])
        else:
            runs[-1].append(value)
    return [[run[0], run[-1]] for run in runs]


def row_context(y: int) -> str:
    if y <= 23:
        return "upper_frame_margin_excluded"
    if y <= 85:
        return "upper_legend_margin_excluded"
    if y < INTERIOR_ROW_MIN:
        return "upper_outer_crossing_jump_margin_excluded"
    if y <= INTERIOR_ROW_MAX:
        return "central_gap_support"
    if y <= 329:
        return "lower_crossing_jump_margin_excluded"
    return "lower_outer_frame_support_excluded"


def select_row(pixels: Any, y: int, rgb: tuple[int, int, int]) -> tuple[dict[str, Any] | None, str, list[list[int]]]:
    """Select a single exact-colour central run; retain every other raw run too."""
    matches = [x for x in range(PANEL_X_MIN, PANEL_X_MAX+1)
               if tuple(int(c) for c in pixels[y, x]) == rgb]
    runs = split_runs(matches)
    if not INTERIOR_ROW_MIN <= y <= INTERIOR_ROW_MAX:
        return None, "outside_conservative_central_gap_support", runs
    if len(runs) != 1:
        return None, "no_exact_colour_run" if not runs else "multiple_exact_colour_runs", runs
    first, last = runs[0]
    if not 2 <= len(matches) <= 8:
        return None, "run_width_outside_2_to_8_pixels", runs
    if first <= PANEL_X_MIN or last >= PANEL_X_MAX:
        return None, "run_touches_panel_ROI_boundary", runs
    return {"pixel_y": y, "pixel_x_median": (first+last)/2,
            "run_first_pixel_x": first, "run_last_pixel_x": last,
            "matched_color_pixels": len(matches)}, "retained", runs


def linear_calibration(pixels: tuple[float, ...], values: tuple[float, ...]) -> dict[str, Any]:
    x = [Fraction(str(value)) for value in pixels]
    y = [Fraction(str(value)) for value in values]
    n = len(x)
    slope = (n*sum(a*b for a, b in zip(x, y))-sum(x)*sum(y))/(n*sum(a*a for a in x)-sum(x)**2)
    intercept = (sum(y)-slope*sum(x))/n
    residuals = [b-slope*a-intercept for a, b in zip(x, y)]
    maximum = max(abs(value) for value in residuals)
    return {"tick_pixel_centres": list(pixels), "tick_values": list(values),
            "slope": float(slope), "intercept": float(intercept),
            "slope_exact": exact_pair(slope), "intercept_exact": exact_pair(intercept),
            "tick_residuals": [float(value) for value in residuals],
            "tick_residuals_exact": [exact_pair(value) for value in residuals],
            "max_abs_tick_residual": float(maximum), "max_abs_tick_residual_exact": exact_pair(maximum),
            "fit": "ordinary_least_squares_exact_rational_tick_coordinates"}


def map_exact(calibration: dict[str, Any], pixel: int | float) -> Fraction:
    return (Fraction(*map(int, calibration["slope_exact"]))*Fraction(str(pixel))
            + Fraction(*map(int, calibration["intercept_exact"])))


def validate_ticks(pixels: Any) -> None:
    magnetic = split_runs([x for x in range(PANEL_X_MIN, PANEL_X_MAX+1)
                           if max(int(c) for c in pixels[417, x]) < 128])
    gap = split_runs([y for y in range(PANEL_ROW_MIN, PANEL_ROW_MAX+1)
                      if max(int(c) for c in pixels[y, 561]) < 128])
    if [(a+b)/2 for a, b in magnetic] != list(MAGNETIC_TICK_PIXELS):
        raise ValueError("magnetic native tick anchors do not match configured calibration")
    if [(a+b)/2 for a, b in gap] != list(GAP_TICK_PIXELS):
        raise ValueError("gap native tick anchors do not match configured calibration")


def make_rows(pixels: Any, figure: dict[str, Any], calibration: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows, catalog = [], {}
    bfit, gfit = calibration["magnetic"], calibration["gap"]
    bres = Fraction(*map(int, bfit["max_abs_tick_residual_exact"]))
    gres = Fraction(*map(int, gfit["max_abs_tick_residual_exact"]))
    for profile, rgb in PROFILES:
        retained, excluded = [], []
        for y in range(PANEL_ROW_MAX, PANEL_ROW_MIN-1, -1):
            selected, reason, runs = select_row(pixels, y, rgb)
            raw = {"pixel_y": y, "exact_colour_runs_inclusive_x": runs,
                   "reason": reason, "display_context": row_context(y)}
            if selected is None:
                excluded.append(raw)
                continue
            b = map_exact(bfit, selected["pixel_x_median"])
            blo = map_exact(bfit, selected["run_first_pixel_x"]-READOUT_PIXELS)-bres
            bhi = map_exact(bfit, selected["run_last_pixel_x"]+READOUT_PIXELS)+bres
            gap = map_exact(gfit, y)
            glo = map_exact(gfit, y+READOUT_PIXELS)-gres
            ghi = map_exact(gfit, y-READOUT_PIXELS)+gres
            rows.append({"figure": figure["figure"], "profile": profile, "time_s": figure["time_s"], **selected,
                         "gap_coordinate_mm": float(gap), "gap_coordinate_readout_lower_mm": float(glo),
                         "gap_coordinate_readout_upper_mm": float(ghi), "Bx_literal_axis_value": float(b),
                         "Bx_literal_readout_lower_axis_value": float(blo), "Bx_literal_readout_upper_axis_value": float(bhi),
                         "Bx_literal_conditional_T": float(CONDITIONAL_TESLA_SCALE*b),
                         "Bx_literal_conditional_readout_lower_T": float(CONDITIONAL_TESLA_SCALE*blo),
                         "Bx_literal_conditional_readout_upper_T": float(CONDITIONAL_TESLA_SCALE*bhi)})
            retained.append(raw)
        catalog[profile] = {"RGB": list(rgb), "retained_rows": retained, "excluded_rows": excluded,
                            "retained_count": len(retained), "excluded_count": len(excluded),
                            "catalogued_row_count": len(retained)+len(excluded)}
    return rows, catalog


def recover(pdf: Path) -> list[dict[str, Any]]:
    """Recover only the new magnetic profiles from the pinned native source."""
    if sha256(pdf.read_bytes()) != PDF_SHA256:
        raise ValueError("source PDF SHA256 mismatch")
    import numpy as np
    import pymupdf
    from PIL import Image

    calibration = {"magnetic": linear_calibration(MAGNETIC_TICK_PIXELS, MAGNETIC_TICK_VALUES),
                   "gap": linear_calibration(GAP_TICK_PIXELS, GAP_TICK_VALUES)}
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
            if any(catalog[name]["retained_count"] != 131 for name, _ in PROFILES):
                raise ValueError("selected central-gap native support changed")
            results.append({"figure": figure, "rows": rows, "catalog": catalog, "calibration": calibration,
                            "image": image.copy(), "encoded_png_sha256": sha256(encoded),
                            "page_text": page.get_text(), "parser_version": pymupdf.VersionBind})
    return results


def make_qa(result: dict[str, Any], path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(12, 7), constrained_layout=True)
    axes[0].imshow(result["image"])
    axes[0].set_xlim(495, 972)
    axes[0].set_ylim(480, 0)
    axes[0].axhspan(INTERIOR_ROW_MIN, INTERIOR_ROW_MAX, color="royalblue", alpha=.10)
    for profile, rgb in PROFILES:
        selected = [row for row in result["rows"] if row["profile"] == profile]
        axes[0].plot([row["pixel_x_median"] for row in selected], [row["pixel_y"] for row in selected],
                     ".", color="navy" if profile == "foil_edge" else "darkred", markersize=1.5)
        colour = tuple(c/255 for c in rgb)
        gap = [row["gap_coordinate_mm"] for row in selected]
        axes[1].plot([row["Bx_literal_axis_value"] for row in selected], gap, color=colour, label=profile)
        axes[1].fill_betweenx(gap, [row["Bx_literal_readout_lower_axis_value"] for row in selected],
                             [row["Bx_literal_readout_upper_axis_value"] for row in selected], color=colour, alpha=.2)
    axes[0].set_title("Native magnetic panel; retained central rows shaded")
    axes[0].set_xlabel("Native source pixel x")
    axes[0].set_ylabel("Native source pixel y")
    axes[1].set_xlabel("Signed literal B_x [100 Tesla]; normalization unverified")
    axes[1].set_ylabel("Distance from Middle of AK Gap [mm]")
    axes[1].legend()
    axes[1].grid(alpha=.2)
    figure.suptitle(f"Dissertation Figure {result['figure']['figure']}: legacy Allegra simulation, "
                   f"{result['figure']['time_s']*1e9:g} ns nominal snapshot\nCorrelated full-run + ±2 px + tick-fit readout bounds; coarse central gap only")
    figure.savefig(path.with_suffix(".png"), dpi=160)
    figure.savefig(path.with_suffix(".pdf"))
    plt.close(figure)


def metadata(result: dict[str, Any], csv_path: Path, qa_paths: list[Path], geometry: Path) -> dict[str, Any]:
    source, rows = result["figure"], result["rows"]
    return {
        "schema": "mhdlab-dissertation-magnetic-gap-recovery-v1",
        "source": {"title": "Characterizing Desorbed Surface Contaminants Leading to Plasma Formation in Power Feeds of Pulsed Power Accelerators",
                   "item_url": "https://deepblue.lib.umich.edu/items/9ccfa886-0b19-4bd8-b3b2-bb685b9ddefc",
                   "pdf_sha256": PDF_SHA256, "figure": source["figure"], "printed_page": source["printed_page"],
                   "pdf_page": source["page_index"]+1, "pdf_page_index": source["page_index"],
                   "page_content_sha256": source["page_content_sha256"], "page_text": result["page_text"],
                   "embedded_image_xref": source["xref"], "native_image_size_pixels": [975, 975], "native_image_mode": "RGB",
                   "reference_extracted_PNG_sha256": source["reference_png_sha256"], "extracted_PNG_sha256": result["encoded_png_sha256"],
                   "decoded_RGB_pixels_sha256": source["rgb_sha256"],
                   "hash_validation": "PDF, page content and decoded RGB mandatory; lossless PNG encoding may vary",
                   "parser": "PyMuPDF", "parser_version": result["parser_version"]},
        "simulation_context": {"code_label_in_dissertation": "Allegra", "drive_shot_caption_label": 17133,
                               "caption_snapshot_time_s": source["time_s"],
                               "snapshot_time_representation": "Direct decimal-to-float seconds for nominal 40 ns/60 ns caption labels; no experimental clock registration.",
                               "observation_type": "published legacy simulation coarse central-gap B_x lineout",
                               "experimental_measurement": False, "FLASH_result": False,
                               "initial_foil_thickness_m": 200e-6, "initial_AK_gap_m": 5e-3,
                               "context_source": "section 4.5 and Figures 4.27–4.29, printed pages 93–96",
                               "actual_material_boundaries_recovered": False, "edge_width_offset_m": None,
                               "full_vector_field_recovered": False, "experimental_clock_alignment_established": False},
        "unit_interpretation": {"source_magnetic_axis_label": "B_x [100 Tesla]",
                                "source_gap_axis_label": "Distance from Middle of AK Gap [mm]",
                                "primary_recovered_quantity": "signed literal B_x plotted axis value",
                                "literal_axis_unit": "100 Tesla", "literal_conditional_scale_to_T": CONDITIONAL_TESLA_SCALE,
                                "original_electromagnetic_normalization_verified": False,
                                "original_simulation_output_units_confirmed": False,
                                "normalization_or_legend_error_resolved": False,
                                "unit_repair_or_empirical_rescaling_applied": False,
                                "component_signs_preserved": True, "coordinate_frame_rotation_applied": False,
                                "SI_columns_qualification": "conditional literal printed factor 100 only; not a validated physical field input"},
        "calibration": {**result["calibration"], "magnetic_unit": "literal B_x [100 Tesla]", "gap_unit": "mm",
                        "tick_detection": {"dark_channel_threshold_exclusive": 128, "magnetic_tick_tail_row": 417,
                                           "gap_tick_tail_column": 561, "pixel_coordinate_convention": "native zero-based integer pixel centres"},
                        "magnetic_fit_input_sha256": sha256(json.dumps([MAGNETIC_TICK_PIXELS, MAGNETIC_TICK_VALUES]).encode()),
                        "gap_fit_input_sha256": sha256(json.dumps([GAP_TICK_PIXELS, GAP_TICK_VALUES]).encode())},
        "selection": {"panel_ROI_inclusive": {"x": [PANEL_X_MIN, PANEL_X_MAX], "y": [PANEL_ROW_MIN, PANEL_ROW_MAX]},
                      "retained_central_y_inclusive": [INTERIOR_ROW_MIN, INTERIOR_ROW_MAX], "row_stride": 1,
                      "mask_basis": "conservative central-gap display support separated from visible crossings, sheets and legend; no material-boundary reconstruction",
                      "rule": "one contiguous exact RGB run of 2–8 pixels away from ROI boundaries",
                      "profiles": result["catalog"],
                      "display_context_qualification": "Exclusion context labels describe conservative display margins, not a claim that every omitted row is a crossing, sheet or occluded pixel.",
                      "interpolation": False, "smoothing": False, "extrapolation": False,
                      "surface_or_endpoint_padding": False, "symmetry_imposed": False},
        "readout_bounds": {"correlated_display_bounds_only": True, "independent_noise_or_confidence_interval": False,
                           "additional_pixel_readout_half_width": READOUT_PIXELS,
                           "magnetic_recipe": "entire matched [first_x,last_x] run extended by +/-2 pixels, exact fitted-axis map, then max absolute magnetic tick-fit residual",
                           "gap_recipe": "native row +/-2 pixels, exact fitted-axis map, then max absolute gap tick-fit residual",
                           "physical_or_simulation_uncertainty_quantified": False},
        "output": {"csv_filename": csv_path.name, "csv_sha256": sha256(csv_path.read_bytes()), "columns": list(FIELDS),
                   "line_ending": "LF", "numeric_projection": "Exact rational coordinates, bounds and conditional scale converted once to float, then CSV formatted to 12 significant digits.",
                   "samples": len(rows), "samples_by_profile": {name: sum(row["profile"] == name for row in rows) for name, _ in PROFILES},
                   "gap_coordinate_extent_mm": [min(row["gap_coordinate_mm"] for row in rows), max(row["gap_coordinate_mm"] for row in rows)],
                   "signed_Bx_literal_extent_by_profile": {name: [min(row["Bx_literal_axis_value"] for row in rows if row["profile"] == name), max(row["Bx_literal_axis_value"] for row in rows if row["profile"] == name)] for name, _ in PROFILES},
                   "generator": "tools/recover_dissertation_gap_magnetic.py", "generator_sha256": sha256(Path(__file__).read_bytes()),
                   "source_geometry_filename": geometry.name, "source_geometry_sha256": sha256(geometry.read_bytes()),
                   "QA_files_sha256": {path.name: sha256(path.read_bytes()) for path in qa_paths}},
        "limits": ["Published legacy Allegra simulation raster, not measured electromagnetic data or an actual FLASH result.",
                   "The central B_x component is signed negative on both coordinate halves; no sign flip from gap-coordinate sign is applied.",
                   "Only finite conservative central-gap support is recovered. Exteriors, crossings, sharp sheets, legend regions and exact surface fields remain excluded.",
                   "Nominal 200 µm foil spans about 7.9 native pixels on this gap axis; no skin-depth profile, derivative, Ampere/current or field-energy closure is supplied.",
                   "Original arrays, mesh, material fractions, moving boundaries, full vector fields and electromagnetic normalization are unavailable.",
                   "Conditional tesla columns do not establish physical magnitude, Hall transport, species magnetization, gap/host field acceptance or a physical coordinate-frame mapping.",
                   "Two nominal snapshots are not a registered field time history; no trigger or current-to-field alignment is inferred.",
                   "Readout bounds and adjacent raster samples are correlated display quantities, not physical uncertainty or independent observations."],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--qa-dir", required=True, type=Path)
    args = parser.parse_args()
    targets = [args.output_dir/f"figure_{source['figure'].replace('.', '_')}_magnetic_gap.{suffix}"
               for source in FIGURES for suffix in ("csv", "json")]
    if any(path.exists() for path in targets):
        raise ValueError("refusing to overwrite existing recovered datasets")
    results = recover(args.pdf)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.qa_dir.mkdir(parents=True, exist_ok=True)
    summary = []
    for result in results:
        stem = f"figure_{result['figure']['figure'].replace('.', '_')}_magnetic_gap"
        csv_path = args.output_dir/f"{stem}.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
            writer.writeheader()
            writer.writerows({key: format(value, ".12g") if isinstance(value, float) else value for key, value in row.items()} for row in result["rows"])
        qa = args.qa_dir/stem
        make_qa(result, qa)
        geometry = args.qa_dir/f"{stem}_source_geometry.json"
        geometry.write_text(json.dumps({"schema": "mhdlab-magnetic-gap-source-geometry-v1", "figure": result["figure"],
                                        "calibration": result["calibration"], "profiles": result["catalog"]}, indent=2, allow_nan=False)+"\n")
        document = metadata(result, csv_path, [qa.with_suffix(".png"), qa.with_suffix(".pdf")], geometry)
        metadata_path = args.output_dir/f"{stem}.json"
        metadata_path.write_text(json.dumps(document, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
        summary.append({"figure": result["figure"]["figure"], "samples": len(result["rows"]),
                        "csv": str(csv_path), "csv_sha256": sha256(csv_path.read_bytes()),
                        "metadata": str(metadata_path), "metadata_sha256": sha256(metadata_path.read_bytes()),
                        "gap_coordinate_extent_mm": document["output"]["gap_coordinate_extent_mm"],
                        "signed_Bx_literal_extent_by_profile": document["output"]["signed_Bx_literal_extent_by_profile"],
                        "QA_files_sha256": document["output"]["QA_files_sha256"]})
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
