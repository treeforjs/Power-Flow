#!/usr/bin/env python3
"""Recover signed interior J_z raster lineouts from dissertation Figs. 4.28/4.29.

The source is a published legacy Allegra simulation, not experimental data.
Literal unit conversions remain conditional on the unrecovered normalization.
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
    {"figure": "4.28", "page_index": 118, "printed_page": 94, "xref": 2847, "time_ns": 40,
     "rgb_sha256": "458b3ced8db4501db3fcbcc90939ce5f336a7ec83eb97ee259c66f7112bf07c8",
     "reference_png_sha256": "1f9c939c8b72947da2629ccad92c42dda04503719c3cdb89fe0e921bc5ce5488",
     "page_content_sha256": "934979de84040dc7970be53747773d2f8c78fedad9ac2b7116757a8052449103"},
    {"figure": "4.29", "page_index": 120, "printed_page": 96, "xref": 2870, "time_ns": 60,
     "rgb_sha256": "ac25b035329213a8d0a354b8496b1674c2be24ceaeff65d95a48ef2e6376b34f",
     "reference_png_sha256": "5ae90775c30a8fad657f74d5813c150b2a4284d1a0fd8075ebc4fcd21b88e89a",
     "page_content_sha256": "17080e390125995d993f2fa580e2fcef6cf4cb166a267bf08a9e1a9798a8d43e"},
)
PROFILES = (("foil_edge", (154, 205, 50)), ("foil_middle", (0, 128, 0)))
CURRENT_TICK_PIXELS = (110.5, 196.5, 283.5, 369.5, 456.5)
CURRENT_TICK_VALUES = (-0.8, -0.6, -0.4, -0.2, 0.0)
GAP_TICK_PIXELS = (900.5, 851.5, 802.5, 752.5, 703.5, 653.5, 604.5, 555.5, 505.5)
GAP_TICK_VALUES = (2.40, 2.45, 2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80)
PANEL_ROW_MIN, PANEL_ROW_MAX = 505, 900
INTERIOR_ROW_MIN, INTERIOR_ROW_MAX = 640, 780
PANEL_X_MIN, PANEL_X_MAX = 90, 477
READOUT_PIXELS = 2.0
INITIAL_FACE_MM = 2.5
CONDITIONAL_SI_SCALE = 1e10
FIELDS = (
    "figure", "profile", "time_s", "pixel_y", "pixel_x_median", "run_first_pixel_x",
    "run_last_pixel_x", "matched_color_pixels", "gap_coordinate_mm",
    "gap_coordinate_readout_lower_mm", "gap_coordinate_readout_upper_mm",
    "nominal_initial_depth_um", "nominal_initial_depth_readout_lower_um",
    "nominal_initial_depth_readout_upper_um", "Jz_literal_MA_per_cm2",
    "Jz_literal_readout_lower_MA_per_cm2", "Jz_literal_readout_upper_MA_per_cm2",
    "Jz_literal_conditional_A_per_m2", "Jz_literal_conditional_readout_lower_A_per_m2",
    "Jz_literal_conditional_readout_upper_A_per_m2",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def split_runs(indices: list[int]) -> list[list[int]]:
    runs: list[list[int]] = []
    for value in indices:
        if not runs or value != runs[-1][-1] + 1:
            runs.append([value])
        else:
            runs[-1].append(value)
    return [[run[0], run[-1]] for run in runs]


def select_row(pixels: Any, y: int, rgb: tuple[int, int, int]) -> tuple[dict[str, Any] | None, str, list[list[int]]]:
    """Select one bounded, unambiguous exact-colour interior run; never fill gaps."""
    matches = [x for x in range(PANEL_X_MIN, PANEL_X_MAX + 1)
               if tuple(int(c) for c in pixels[y, x]) == rgb]
    runs = split_runs(matches)
    if not INTERIOR_ROW_MIN <= y <= INTERIOR_ROW_MAX:
        return None, "outside_conservative_interior_display_mask", runs
    if len(runs) != 1:
        return None, "no_exact_colour_run" if not runs else "multiple_exact_colour_runs", runs
    first, last = runs[0]
    if not 2 <= len(matches) <= 8:
        return None, "run_width_outside_2_to_8_pixels", runs
    if first <= PANEL_X_MIN or last >= PANEL_X_MAX:
        return None, "run_touches_panel_ROI_boundary", runs
    return {"pixel_y": y, "pixel_x_median": (first + last) / 2,
            "run_first_pixel_x": first, "run_last_pixel_x": last,
            "matched_color_pixels": len(matches)}, "retained", runs


def linear_calibration(pixels: tuple[float, ...], values: tuple[float, ...]) -> dict[str, Any]:
    """Ordinary least squares in exact rational arithmetic on declared tick centres."""
    x = [Fraction(str(v)) for v in pixels]
    z = [Fraction(str(v)) for v in values]
    n = len(x)
    slope = (n * sum(a*b for a, b in zip(x, z)) - sum(x)*sum(z)) / (n*sum(a*a for a in x)-sum(x)**2)
    intercept = (sum(z)-slope*sum(x))/n
    residuals = [b-(slope*a+intercept) for a, b in zip(x, z)]
    return {"tick_pixel_centres": list(pixels), "tick_values": list(values),
            "slope": float(slope), "intercept": float(intercept),
            "slope_exact": [str(slope.numerator), str(slope.denominator)],
            "intercept_exact": [str(intercept.numerator), str(intercept.denominator)],
            "tick_residuals": [float(r) for r in residuals],
            "max_abs_tick_residual": float(max(abs(r) for r in residuals)),
            "fit": "ordinary_least_squares_exact_rational_tick_coordinates"}


def validate_ticks(pixels: Any) -> None:
    """Bind the configured anchors to separately visible dark tick tails."""
    current = split_runs([x for x in range(90, 478) if max(int(c) for c in pixels[904, x]) < 128])
    gap = split_runs([y for y in range(505, 902) if max(int(c) for c in pixels[y, 84]) < 128])
    if [(a+b)/2 for a, b in current] != list(CURRENT_TICK_PIXELS):
        raise ValueError("J_z native tick anchors do not match configured calibration")
    if [(a+b)/2 for a, b in reversed(gap)] != list(GAP_TICK_PIXELS):
        raise ValueError("gap native tick anchors do not match configured calibration")


def transform(calibration: dict[str, Any], pixel: float) -> float:
    return calibration["slope"] * pixel + calibration["intercept"]


def make_rows(pixels: Any, figure: dict[str, Any], calibration: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    catalog: dict[str, Any] = {}
    jfit, gfit = calibration["current"], calibration["gap"]
    for profile, rgb in PROFILES:
        retained, excluded = [], []
        for y in range(PANEL_ROW_MAX, PANEL_ROW_MIN-1, -1):
            selected, reason, runs = select_row(pixels, y, rgb)
            raw = {"pixel_y": y, "exact_colour_runs_inclusive_x": runs, "reason": reason}
            if selected is None:
                excluded.append(raw)
                continue
            # Bounds include the complete matched run, an additional +/-2px readout,
            # and the largest tick-fit residual. They are correlated display bounds.
            x = selected["pixel_x_median"]
            jlo = transform(jfit, selected["run_first_pixel_x"]-READOUT_PIXELS)-jfit["max_abs_tick_residual"]
            jhi = transform(jfit, selected["run_last_pixel_x"]+READOUT_PIXELS)+jfit["max_abs_tick_residual"]
            gap = transform(gfit, y)
            glo = transform(gfit, y+READOUT_PIXELS)-gfit["max_abs_tick_residual"]
            ghi = transform(gfit, y-READOUT_PIXELS)+gfit["max_abs_tick_residual"]
            jz = transform(jfit, x)
            rows.append({"figure": figure["figure"], "profile": profile,
                         "time_s": figure["time_ns"]*1e-9, **selected,
                         "gap_coordinate_mm": gap,
                         "gap_coordinate_readout_lower_mm": glo,
                         "gap_coordinate_readout_upper_mm": ghi,
                         "nominal_initial_depth_um": (gap-INITIAL_FACE_MM)*1000,
                         "nominal_initial_depth_readout_lower_um": (glo-INITIAL_FACE_MM)*1000,
                         "nominal_initial_depth_readout_upper_um": (ghi-INITIAL_FACE_MM)*1000,
                         "Jz_literal_MA_per_cm2": jz,
                         "Jz_literal_readout_lower_MA_per_cm2": jlo,
                         "Jz_literal_readout_upper_MA_per_cm2": jhi,
                         "Jz_literal_conditional_A_per_m2": jz*CONDITIONAL_SI_SCALE,
                         "Jz_literal_conditional_readout_lower_A_per_m2": jlo*CONDITIONAL_SI_SCALE,
                         "Jz_literal_conditional_readout_upper_A_per_m2": jhi*CONDITIONAL_SI_SCALE})
            retained.append(raw)
        catalog[profile] = {"RGB": list(rgb), "retained_rows": retained, "excluded_rows": excluded,
                            "retained_count": len(retained), "excluded_count": len(excluded),
                            "catalogued_row_count": len(retained)+len(excluded)}
    return rows, catalog


def recover(pdf: Path) -> list[dict[str, Any]]:
    """Read the pinned source and recover both new figures; does not write files."""
    if sha256(pdf.read_bytes()) != PDF_SHA256:
        raise ValueError("source PDF SHA256 mismatch")
    import numpy as np
    import pymupdf
    from PIL import Image

    calibration = {"current": linear_calibration(CURRENT_TICK_PIXELS, CURRENT_TICK_VALUES),
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
            if any(catalog[p]["retained_count"] != 141 for p, _ in PROFILES):
                raise ValueError("selected native interior row support changed")
            if not all(row["Jz_literal_readout_upper_MA_per_cm2"] < 0 for row in rows):
                raise ValueError("selected current readout envelopes must remain signed negative")
            results.append({"figure": figure, "image": image.copy(), "rows": rows, "catalog": catalog,
                            "calibration": calibration, "encoded_png_sha256": sha256(encoded),
                            "page_text": page.get_text(), "parser_version": pymupdf.VersionBind})
    return results


def make_qa(result: dict[str, Any], path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(12, 7), constrained_layout=True)
    axes[0].imshow(result["image"])
    axes[0].set_xlim(45, 492)
    axes[0].set_ylim(930, 485)
    axes[0].axhspan(INTERIOR_ROW_MIN, INTERIOR_ROW_MAX, color="royalblue", alpha=.10)
    for profile, rgb in PROFILES:
        selected = [r for r in result["rows"] if r["profile"] == profile]
        colour = tuple(c/255 for c in rgb)
        axes[0].plot([r["pixel_x_median"] for r in selected], [r["pixel_y"] for r in selected],
                     ".", color="crimson" if profile == "foil_edge" else "navy", markersize=1.5)
        depth = [r["nominal_initial_depth_um"] for r in selected]
        axes[1].plot([r["Jz_literal_MA_per_cm2"] for r in selected], depth, color=colour, label=profile)
        axes[1].fill_betweenx(depth, [r["Jz_literal_readout_lower_MA_per_cm2"] for r in selected],
                             [r["Jz_literal_readout_upper_MA_per_cm2"] for r in selected], color=colour, alpha=.2)
    axes[0].set_title("Native J panel; conservative retained rows shaded")
    axes[0].set_xlabel("Native source pixel x")
    axes[0].set_ylabel("Native source pixel y")
    axes[1].set_xlabel("Signed literal J_z [MA/cm²]; normalization unverified")
    axes[1].set_ylabel("Nominal initial depth [µm] from gap coordinate 2.5 mm")
    axes[1].legend()
    axes[1].grid(alpha=.2)
    figure.suptitle(f"Dissertation Figure {result['figure']['figure']}: legacy Allegra simulation, "
                   f"{result['figure']['time_ns']} ns\nCorrelated full-run + ±2 px + tick-fit readout envelopes")
    figure.savefig(path.with_suffix(".png"), dpi=160)
    figure.savefig(path.with_suffix(".pdf"))
    plt.close(figure)


def metadata(result: dict[str, Any], csv_path: Path, qa_paths: list[Path], geometry: Path) -> dict[str, Any]:
    figure = result["figure"]
    rows = result["rows"]
    return {
        "schema": "mhdlab-dissertation-current-depth-recovery-v1",
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
        "simulation_context": {"code_label_in_dissertation": "Allegra", "drive_shot": 17133,
                               "time_s": figure["time_ns"]*1e-9, "observation_type": "published legacy simulation depth lineout",
                               "experimental_measurement": False, "FLASH_result": False,
                               "initial_foil_thickness_m": 200e-6, "initial_AK_gap_m": 5e-3,
                               "dielectric_support_included": False,
                               "context_source": "section 4.5 and Figures 4.27–4.29, printed pages 93–96",
                               "edge_width_offset_m": None, "actual_material_boundaries_recovered": False,
                               "coordinate_frame_rotation_applied": False, "component_signs_preserved": True},
        "unit_interpretation": {"source_current_axis_label": "J_z [MA/cm²]",
                                "source_gap_axis_label": "Distance from Middle of AK Gap [mm]",
                                "primary_recovered_quantity": "signed literal plotted J_z ordinate",
                                "literal_conditional_scale_to_A_per_m2": CONDITIONAL_SI_SCALE,
                                "literal_conversion_derivation": "1e6 A / (1e-2 m)^2 = 1e10 A/m²",
                                "original_simulation_output_units_confirmed": False,
                                "original_normalization_verified": False,
                                "normalization_or_legend_error_resolved": False,
                                "unit_repair_or_empirical_rescaling_applied": False,
                                "SI_columns_qualification": "conditional literal printed-unit conversion only; not a validated dimensional benchmark"},
        "calibration": {**result["calibration"],
                        "current_unit": "literal MA/cm²", "gap_unit": "mm",
                        "tick_detection": {"dark_channel_threshold_exclusive": 128, "current_tick_tail_row": 904,
                                           "gap_tick_tail_column": 84, "pixel_coordinate_convention": "native zero-based integer pixel centres"},
                        "current_fit_input_sha256": sha256(json.dumps([CURRENT_TICK_PIXELS, CURRENT_TICK_VALUES]).encode()),
                        "gap_fit_input_sha256": sha256(json.dumps([GAP_TICK_PIXELS, GAP_TICK_VALUES]).encode())},
        "selection": {"panel_ROI_inclusive": {"x": [PANEL_X_MIN, PANEL_X_MAX], "y": [PANEL_ROW_MIN, PANEL_ROW_MAX]},
                      "retained_interior_y_inclusive": [INTERIOR_ROW_MIN, INTERIOR_ROW_MAX], "row_stride": 1,
                      "mask_basis": "conservative display interior with margins from visible corners/ramps; not material-boundary reconstruction",
                      "run_rule": "one contiguous exact RGB run of 2–8 pixels, separated from ROI boundaries",
                      "profiles": result["catalog"], "interpolation": False, "smoothing": False,
                      "extrapolation": False, "surface_or_endpoint_padding": False},
        "readout_bounds": {"correlated_display_bounds_only": True, "independent_noise_or_confidence_interval": False,
                           "additional_pixel_readout_half_width": READOUT_PIXELS,
                           "current_recipe": "entire matched run [first_x,last_x], extended by +/-2 pixels, mapped by fitted axis and expanded by max absolute current tick-fit residual",
                           "gap_recipe": "native row +/-2 pixels, mapped by fitted axis and expanded by max absolute gap tick-fit residual",
                           "physical_or_simulation_uncertainty_quantified": False},
        "nominal_initial_depth": {"definition": "1000 * (gap_coordinate_mm - 2.5) µm",
                                  "initial_positive_gap_face_mm": INITIAL_FACE_MM,
                                  "basis": "initial 5 mm symmetric gap and 200 µm foil in section 4.5/Figure 4.27",
                                  "actual_moving_face_coordinate": False, "initial_geometry_reference_only": True},
        "output": {"csv_filename": csv_path.name, "csv_sha256": sha256(csv_path.read_bytes()), "columns": list(FIELDS),
                   "samples": len(rows), "samples_by_profile": {p: sum(r["profile"] == p for r in rows) for p, _ in PROFILES},
                   "gap_coordinate_extent_mm": [min(r["gap_coordinate_mm"] for r in rows), max(r["gap_coordinate_mm"] for r in rows)],
                   "nominal_initial_depth_extent_um": [min(r["nominal_initial_depth_um"] for r in rows), max(r["nominal_initial_depth_um"] for r in rows)],
                   "signed_Jz_literal_extent_by_profile": {p: [min(r["Jz_literal_MA_per_cm2"] for r in rows if r["profile"] == p), max(r["Jz_literal_MA_per_cm2"] for r in rows if r["profile"] == p)] for p, _ in PROFILES},
                   "generator": "tools/recover_dissertation_depth_current.py", "generator_sha256": sha256(Path(__file__).read_bytes()),
                   "source_geometry_filename": geometry.name, "source_geometry_sha256": sha256(geometry.read_bytes()),
                   "QA_files_sha256": {p.name: sha256(p.read_bytes()) for p in qa_paths}},
        "limits": ["Original simulation arrays, mesh, material fractions, component-frame mapping and normalization are unavailable.",
                   "The signed negative J_z values are preserved; no absolute-value replacement or sign rotation is applied.",
                   "Only a conservative smooth interior is retained. Legends, exterior/zero branches, corners, boundary ramps and every other unselected native row are catalogued and excluded.",
                   "Raster samples and their readout bounds are correlated, not independent experimental observations.",
                   "No moving-face coordinate, surface value, boundary peak, whole-depth integral or full cross-section is recovered.",
                   "No current-integral, Ampere-law, J-squared heating, conductivity, energy closure or drive-clock alignment is inferred.",
                   "The conditional SI columns do not verify the original printed normalization or authorize physical heating inputs.",
                   "Magnetic gap panels remain pending: their coarse foil support is not recovered here."],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--qa-dir", required=True, type=Path)
    args = parser.parse_args()
    targets = [args.output_dir/f"figure_{f['figure'].replace('.', '_')}_current_depth.{suffix}"
               for f in FIGURES for suffix in ("csv", "json")]
    if any(p.exists() for p in targets):
        raise ValueError("refusing to overwrite existing recovered datasets")
    results = recover(args.pdf)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.qa_dir.mkdir(parents=True, exist_ok=True)
    summary = []
    for result in results:
        stem = f"figure_{result['figure']['figure'].replace('.', '_')}_current_depth"
        csv_path = args.output_dir/f"{stem}.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows({k: format(v, ".12g") if isinstance(v, float) else v for k, v in row.items()} for row in result["rows"])
        qa = args.qa_dir/stem
        make_qa(result, qa)
        geometry = args.qa_dir/f"{stem}_source_geometry.json"
        geometry.write_text(json.dumps({"schema": "mhdlab-current-depth-source-geometry-v1", "figure": result["figure"],
                                        "calibration": result["calibration"], "profiles": result["catalog"]}, indent=2, allow_nan=False)+"\n")
        qa_paths = [qa.with_suffix(".png"), qa.with_suffix(".pdf")]
        document = metadata(result, csv_path, qa_paths, geometry)
        metadata_path = args.output_dir/f"{stem}.json"
        metadata_path.write_text(json.dumps(document, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
        summary.append({"figure": result["figure"]["figure"], "samples": len(result["rows"]),
                        "csv": str(csv_path), "csv_sha256": sha256(csv_path.read_bytes()),
                        "metadata": str(metadata_path), "metadata_sha256": sha256(metadata_path.read_bytes()),
                        "QA_files_sha256": document["output"]["QA_files_sha256"]})
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
