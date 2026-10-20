#!/usr/bin/env python3
"""Recover signed electromagnetic width profiles from dissertation Fig. 4.30.

The source is the native 975 x 975 RGB raster on printed page 98 / PDF page 123.
Its exact legends are B_x [100 Tesla] (purple) and J_z [MA/cm²] (green). The
shared ordinate says Signal Height. Preserve that ordinate as the primary
recovery. Separate SI columns apply only the literal legend interpretation:
one B ordinate unit is 100 T; one J ordinate unit is 1e6 A / 1e-4 m² = 1e10 A/m².
Original output units/normalization are unavailable. No unit repair, sign flip,
field rotation or normalization to a drive current is performed.

Inspected fixed calibration, 2026-10-18; native pixel origin upper left:
* x=122.5,258.5,395.5,531.5,668.5,804.5,941.5 -> -3,-2,-1,0,1,2,3 mm.
* y=98.5,254.5,410.5,566.5,721.5,877.5 -> 0.2,0,-0.2,-0.4,-0.6,-0.8.
* B_x RGB=(186,85,211), y ROI 340:700; J_z RGB=(154,205,50), y ROI 430:835.
* Sample columns 199:864 at stride 5, plus 531. Exclude crossing neighbourhoods
  x=215:225 and 838:848, and legend rectangle x=392:672,y=744:867. Exact sharp
  edge jumps and the outermost approximately 63 um on each side are excluded.
* Accept one contiguous 2--24 exact-colour pixel run per column and use its
  median. Wider vertical spans on steep smooth curves are retained with their
  full visible-run envelope, rather than treated as precise point estimates.

Position envelopes combine estimated +/-2 native pixels and tick-fit residual.
Ordinate envelopes include the exact-colour run extent, +/-2 pixels and tick-fit
residual. These quantify display readout only, not physical uncertainty.
The curve is a 40 ns, 12.5 um subsurface legacy simulation lineout driven by shot
17133. A single-depth width profile does not determine an integrated current,
current continuity, an Ampere-law check, or full magnetic/current vector fields.

Dependencies: numpy, pillow, pymupdf. Optional --qa-plot needs matplotlib.
Example:
  python tools/recover_dissertation_em_profiles.py --pdf /path/to/smtrevor_1.pdf \
      --qa-plot runs/recovery/figure_4_30_electromagnetic_qa.png
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from PIL import Image
import pymupdf


ROOT = Path(__file__).resolve().parents[1]
PDF_SHA256 = "1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7"
REFERENCE_PNG_SHA256 = "eafdd82cc2204c891aeb6cfe350e628561d163b44f3814d01c57fd78a2727627"
RGB_PIXELS_SHA256 = "05cf38abf7718721bd5a22f9f4b80f005810941fb0d75d765ad715a4e5fd2a90"
X_TICKS_PX = np.array([122.5, 258.5, 395.5, 531.5, 668.5, 804.5, 941.5])
X_TICKS_MM = np.arange(-3.0, 4.0)
Y_TICKS_PX = np.array([98.5, 254.5, 410.5, 566.5, 721.5, 877.5])
SIGNAL_TICKS = np.array([0.2, 0.0, -0.2, -0.4, -0.6, -0.8])
PROFILES = (
    {"quantity": "B_x", "legend": "B_x [100 Tesla]", "rgb": (186, 85, 211),
     "y_roi_inclusive": (340, 700), "literal_scale_to_SI": 100.0, "SI_unit": "T"},
    {"quantity": "J_z", "legend": "J_z [MA/cm²]", "rgb": (154, 205, 50),
     "y_roi_inclusive": (430, 835), "literal_scale_to_SI": 1e10, "SI_unit": "A/m2"},
)
CROSSING_EXCLUSIONS = ((215, 225), (838, 848))
LEGEND_RECTANGLE = (392, 744, 672, 867)
READOUT_PX = 2.0
FIELDS = (
    "quantity", "source_legend", "pixel_x", "pixel_y_median", "exact_color_run_first_y",
    "exact_color_run_last_y", "matched_exact_color_pixels", "width_offset_mm",
    "width_offset_digitization_lower_mm", "width_offset_digitization_upper_mm",
    "plotted_signal_height", "signal_height_digitization_lower", "signal_height_digitization_upper",
    "value_from_literal_legend_SI", "literal_legend_value_digitization_lower_SI",
    "literal_legend_value_digitization_upper_SI", "literal_legend_SI_unit",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_image(pdf: Path) -> tuple[np.ndarray, str]:
    if sha256(pdf.read_bytes()) != PDF_SHA256:
        raise ValueError("PDF differs from the inspected dissertation")
    with pymupdf.open(pdf) as document:
        if len(document) != 203:
            raise ValueError("Expected 203 PDF pages")
        page = document[122]
        if "Figure 4.30" not in page.get_text() or 2890 not in {entry[0] for entry in page.get_images(full=True)}:
            raise ValueError("The inspected Figure 4.30 image is absent from PDF page 123")
        encoded = document.extract_image(2890)["image"]
    image = Image.open(io.BytesIO(encoded))
    if image.mode != "RGB" or image.size != (975, 975) or sha256(image.tobytes()) != RGB_PIXELS_SHA256:
        raise ValueError("Native RGB pixels changed; inspect and recalibrate before recovery")
    return np.asarray(image), sha256(encoded)


def recover(image: np.ndarray) -> tuple[list[dict], dict]:
    x_fit = np.polyfit(X_TICKS_PX, X_TICKS_MM, 1)
    signal_fit = np.polyfit(Y_TICKS_PX, SIGNAL_TICKS, 1)
    x_residual = X_TICKS_MM - np.polyval(x_fit, X_TICKS_PX)
    signal_residual = SIGNAL_TICKS - np.polyval(signal_fit, Y_TICKS_PX)
    x_bound = abs(float(x_fit[0])) * READOUT_PX + float(max(abs(x_residual)))
    signal_calibration_bound = float(max(abs(signal_residual)))
    columns = sorted(set(range(199, 865, 5)) | {531})
    rows, excluded = [], []
    for profile in PROFILES:
        lower, upper = profile["y_roi_inclusive"]
        for pixel_x in columns:
            if any(left <= pixel_x <= right for left, right in CROSSING_EXCLUSIONS):
                excluded.append({"quantity": profile["quantity"], "pixel_x": pixel_x,
                                 "reason": "inspected_curve_crossing_neighbourhood"})
                continue
            matches = np.flatnonzero(np.all(image[lower:upper + 1, pixel_x] == profile["rgb"], axis=1)) + lower
            x0, y0, x1, y1 = LEGEND_RECTANGLE
            if x0 <= pixel_x <= x1:
                matches = matches[(matches < y0) | (matches > y1)]
            reason = None
            if len(matches) < 2:
                reason = "fewer_than_two_unoccluded_exact_color_pixels"
            elif len(matches) > 24:
                reason = "overly_broad_or_near_vertical_run"
            elif np.any(np.diff(matches) != 1):
                reason = "multiple_separated_exact_color_runs"
            elif matches[0] == lower or matches[-1] == upper:
                reason = "run_touches_ROI_boundary"
            if reason:
                excluded.append({"quantity": profile["quantity"], "pixel_x": pixel_x,
                                 "reason": reason, "matching_pixels_y": matches.tolist()})
                continue
            pixel_y = float(np.median(matches))
            x_mm = float(np.polyval(x_fit, pixel_x))
            signal = float(np.polyval(signal_fit, pixel_y))
            signal_lower = float(np.polyval(signal_fit, matches[-1] + READOUT_PX)) - signal_calibration_bound
            signal_upper = float(np.polyval(signal_fit, matches[0] - READOUT_PX)) + signal_calibration_bound
            if not signal_lower < signal < signal_upper < 0:
                raise ValueError("Expected a signed negative EM curve with ordered display envelopes")
            scale = profile["literal_scale_to_SI"]
            rows.append({
                "quantity": profile["quantity"], "source_legend": profile["legend"],
                "pixel_x": pixel_x, "pixel_y_median": pixel_y,
                "exact_color_run_first_y": int(matches[0]), "exact_color_run_last_y": int(matches[-1]),
                "matched_exact_color_pixels": len(matches), "width_offset_mm": x_mm,
                "width_offset_digitization_lower_mm": x_mm - x_bound,
                "width_offset_digitization_upper_mm": x_mm + x_bound,
                "plotted_signal_height": signal,
                "signal_height_digitization_lower": signal_lower, "signal_height_digitization_upper": signal_upper,
                "value_from_literal_legend_SI": signal * scale,
                "literal_legend_value_digitization_lower_SI": signal_lower * scale,
                "literal_legend_value_digitization_upper_SI": signal_upper * scale,
                "literal_legend_SI_unit": profile["SI_unit"],
            })
    if not rows:
        raise ValueError("No unambiguous EM profile recovered")
    calibration = {
        "inspected_date": "2026-10-18", "coordinate_origin": "upper_left_native_raster",
        "method": "least_squares_linear_fit_to_inspected_tick_centres",
        "x_axis": {"label": "Distance from Foil-Middle Location [mm]", "tick_pixels_x": X_TICKS_PX.tolist(),
                   "tick_values_mm": X_TICKS_MM.tolist(), "mm_per_pixel": float(x_fit[0]),
                   "intercept_mm": float(x_fit[1]), "tick_residuals_mm": x_residual.tolist()},
        "y_axis": {"label": "Signal Height", "tick_pixels_y": Y_TICKS_PX.tolist(),
                   "tick_values": SIGNAL_TICKS.tolist(), "signal_height_per_pixel": float(signal_fit[0]),
                   "intercept": float(signal_fit[1]), "tick_residuals": signal_residual.tolist()},
        "selection": {"profiles": [{**p, "rgb": list(p["rgb"]), "y_roi_inclusive": list(p["y_roi_inclusive"])} for p in PROFILES],
                      "sample_columns_before_exclusions": columns,
                      "crossing_excluded_x_intervals_inclusive": [list(interval) for interval in CROSSING_EXCLUSIONS],
                      "legend_excluded_rectangle_xyxy_inclusive": list(LEGEND_RECTANGLE),
                      "accepted_contiguous_run_pixels": [2, 24],
                      "line_position": "median_of_single_contiguous_unoccluded_exact_color_run",
                      "interpolation": False, "extrapolation": False, "symmetry_enforced": False},
        "digitization_envelopes": {"estimated_coordinate_readout_pixels_plus_minus": READOUT_PX,
                                   "width_offset_mm_plus_minus": x_bound,
                                   "ordinate_rule": "convert full exact-colour run extent expanded by +/-2 pixels; add maximum absolute tick-fit residual",
                                   "ordinate_calibration_residual_bound": signal_calibration_bound,
                                   "meaning": "display readout only; not physical uncertainty or confidence intervals"},
        "excluded_samples": excluded,
    }
    return rows, calibration


def make_qa(image: np.ndarray, rows: list[dict], path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, (source_ax, profile_ax) = plt.subplots(1, 2, figsize=(12, 5.7), layout="constrained")
    source_ax.imshow(image, origin="upper")
    source_ax.set(xlim=(0, 975), ylim=(975, 0), title="Native Fig. 4.30; black recovered medians")
    source_ax.axis("off")
    for profile in PROFILES:
        selected = [r for r in rows if r["quantity"] == profile["quantity"]]
        source_ax.scatter([r["pixel_x"] for r in selected], [r["pixel_y_median"] for r in selected], s=5, color="black")
        color = np.array(profile["rgb"]) / 255.0
        x = [r["width_offset_mm"] for r in selected]
        signal = [r["plotted_signal_height"] for r in selected]
        lower = [r["signal_height_digitization_lower"] for r in selected]
        upper = [r["signal_height_digitization_upper"] for r in selected]
        # Plot points only: no line or envelope is drawn across excluded crossings.
        profile_ax.errorbar(x, signal, yerr=[np.subtract(signal, lower), np.subtract(upper, signal)],
                            fmt=".", ms=3, lw=0.5, color=color, label=profile["legend"])
    profile_ax.set(xlabel="Distance from foil middle [mm]", ylabel="Published Signal Height",
                   title="Signed source ordinates; literal units remain unverified")
    profile_ax.grid(alpha=0.2)
    profile_ax.legend(fontsize=8)
    figure.suptitle("Legacy 40 ns / 12.5 µm lineout; crossings and boundary jumps excluded", fontsize=10)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=150)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--csv", type=Path, default=ROOT / "data/dissertation/figure_4_30_electromagnetic.csv")
    parser.add_argument("--metadata", type=Path, default=ROOT / "data/dissertation/figure_4_30_electromagnetic.json")
    parser.add_argument("--qa-plot", type=Path)
    args = parser.parse_args()
    image, encoded_hash = load_image(args.pdf)
    rows, calibration = recover(image)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format(value, ".12g") if isinstance(value, float) else value for key, value in row.items()})
    metadata = {
        "schema_version": 1, "status": "partially_digitized_published_simulation_with_unverified_legend_units",
        "source": {"author": "Trevor Johannes Smith", "year": 2025,
                   "title": "Characterizing Desorbed Surface Contaminants Leading to Plasma Formation in Power Feeds of Pulsed Power Accelerators",
                   "item_url": "https://deepblue.lib.umich.edu/items/9ccfa886-0b19-4bd8-b3b2-bb685b9ddefc",
                   "pdf_sha256": PDF_SHA256, "figure": "4.30", "printed_page": 98, "pdf_page": 123,
                   "embedded_image_xref": 2890, "native_image_size_pixels": [975, 975], "native_image_mode": "RGB",
                   "reference_extracted_PNG_sha256": REFERENCE_PNG_SHA256, "extracted_PNG_sha256": encoded_hash,
                   "decoded_RGB_pixels_sha256": RGB_PIXELS_SHA256,
                   "hash_validation": "PDF and decoded native RGB pixels mandatory; lossless PNG encoding may vary"},
        "simulation_context": {"code_label_in_dissertation": "Allegra", "drive_shot": 17133,
                               "time_s": 40e-9, "depth_from_AK_gap_m": 12.5e-6,
                               "initial_foil_thickness_m": 200e-6, "initial_foil_width_m": 5e-3,
                               "initial_AK_gap_m": 5e-3, "dielectric_support_included": False,
                               "context_source": "section 4.5 and Figures 4.27–4.30, printed pages 93–98",
                               "observation_type": "published legacy simulation width lineout",
                               "experimental_measurement": False, "FLASH_result": False, "surface_lineout": False,
                               "drive_15210_used": False},
        "unit_interpretation": {
            "primary_recovered_quantity": "signed plotted Signal Height; original legend strings retained",
            "normalization_statement_search_scope": "section 4.5, printed pages 93–98, and Figure 4.30 legend/caption",
            "explicit_normalization_recipe_found": False,
            "literal_legend_conversions": [
                {"quantity": "B_x", "source_legend": "B_x [100 Tesla]", "assumption": "one ordinate unit denotes 100 tesla",
                 "scale_to_SI": 100.0, "SI_unit": "T"},
                {"quantity": "J_z", "source_legend": "J_z [MA/cm²]", "assumption": "one ordinate unit denotes one megaampere per square centimetre",
                 "scale_to_SI": 1e10, "SI_unit": "A/m2", "derivation": "1e6 A / (1e-2 m)^2"},
            ],
            "original_simulation_output_units_confirmed": False,
            "normalization_or_legend_error_resolved": False,
            "unit_repair_or_empirical_rescaling_applied": False,
            "source_component_signs_preserved": True,
            "vector_coordinate_rotation_applied": False,
            "SI_columns_qualification": "conditional literal-legend interpretation only; not a validated dimensional benchmark",
        },
        "calibration": calibration,
        "output": {"csv_filename": args.csv.name, "csv_sha256": sha256(args.csv.read_bytes()), "columns": list(FIELDS),
                   "samples": len(rows), "samples_by_quantity": {p["quantity"]: sum(r["quantity"] == p["quantity"] for r in rows) for p in PROFILES},
                   "width_offset_extent_mm": [min(r["width_offset_mm"] for r in rows), max(r["width_offset_mm"] for r in rows)],
                   "signal_height_extent_by_quantity": {p["quantity"]: [min(r["plotted_signal_height"] for r in rows if r["quantity"] == p["quantity"]),
                                                                          max(r["plotted_signal_height"] for r in rows if r["quantity"] == p["quantity"])] for p in PROFILES},
                   "generator": "tools/recover_dissertation_em_profiles.py", "generator_sha256": sha256(Path(__file__).read_bytes())},
        "dimensional_consistency_limits": [
            "The shared ordinate is labeled Signal Height; a separate source normalization recipe and original simulation-unit export are unavailable.",
            "The B_x [100 Tesla] convention is interpreted literally as 100 T per ordinate unit; this interpretation is explicit and remains unverified.",
            "The J_z conversion uses the printed cm² denominator exactly. No alternative denominator or factor-of-100 adjustment is inferred or applied.",
            "A one-depth width lineout does not supply the full current-density cross-section; multiplying it by foil thickness would introduce an unsupported depth-uniformity assumption.",
            "No total-current integral, Ampere-law agreement, magnetic divergence check, Joule-power closure or normalization to the separate drive history is claimed.",
            "Only named components are recovered. Other components, their gradients and a mapping of the original vector frame to a host solver are not supplied.",
        ],
        "limits": [
            "Original simulation arrays, mesh, material fractions, unit conventions and current drive are not recovered.",
            "The exact-color medians are correlated published-raster samples, not independent experimental observations.",
            "Near-vertical foil-edge jumps, outermost edge regions, exterior curves, inspected crossings and legend glyphs are excluded; gaps are not interpolated.",
            "Current-density and field peaks at excluded boundaries are not recovered; extrema in this CSV describe only retained samples.",
            "Display envelopes include line-run extent and calibration/readout estimates, not physical simulation or measurement uncertainty.",
            "The temperature tool and data are separate frozen artifacts and are not modified by this recovery.",
        ],
    }
    args.metadata.parent.mkdir(parents=True, exist_ok=True)
    args.metadata.write_text(json.dumps(metadata, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    if args.qa_plot:
        make_qa(image, rows, args.qa_plot)
    print(f"Recovered {len(rows)} signed source-ordinate samples: {args.csv}")
    print(f"Literal-legend SI conversion remains unverified; metadata: {args.metadata}")


if __name__ == "__main__":
    main()
