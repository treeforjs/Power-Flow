#!/usr/bin/env python3
"""Recover the two central profiles in dissertation Figure 4.20, reproducibly.

This is a narrowly calibrated raster digitizer, not a general figure recognizer.
The source is Trevor Smith's 2025 dissertation, printed page 84 / PDF page 109
(zero-based index 108). Its embedded RGB image, xref 2707, is 1040 x 692 pixels.
All coordinates below refer to that original embedded image, with origin at its
upper-left corner. Do not calibrate using a screenshot or rescaled rendering.

Manual calibration, inspected on 2026-10-18:
* x pixels 154, 269, 383 correspond to log10(density / cm^-3) = 10, 11, 12;
  a least-squares straight-line calibration handles the unequal tick spacing.
* y pixel 356 is 0 mm from the cathode; y pixel 82 is 5 mm.
* H Lyman-alpha: exact RGB (255, 127, 14), inclusive central-line x ROI 210:261.
* H2 Lyman band at 161 nm: RGB (0, 191, 255), inclusive x ROI 305:361.
* Sample rows 87 through 351 inclusive, at stride 5 (the last sampled row is 347).
  Use the median of a single contiguous run of 3--8 exact-color pixels; skip
  thinner, clipped, separated, or overly broad matches as ambiguous.

The +/-3 horizontal pixels and +/-1 vertical pixel are conservative digitization
estimates, not physical error bars or recovered published uncertainty. These
curves show plotted transition densities, not ground-state neutral densities.
Samples on the same raster curve are correlated and are not raw observations.

Core dependencies: pymupdf, pillow, numpy. Optional --qa-plot needs matplotlib.
Example (run from repository root):
  python tools/recover_dissertation_profiles.py --pdf /path/to/smtrevor_1.pdf \
      --qa-plot runs/recovery/figure_4_20_qa.png
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pymupdf
from PIL import Image


PDF_SHA256 = "1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7"
IMAGE_SHA256 = "55d3d42a310319769d6d231fafed42f925bbb8c22620f965a5bd015239c35d25"
RGB_PIXELS_SHA256 = "1abbc44c3138785c31690a0732378e8fb7c8e14114cf0d85ea7ca4e3c0337f04"
PDF_PAGE_INDEX = 108
IMAGE_XREF = 2707
IMAGE_SIZE = (1040, 692)
X_TICKS = np.array([154.0, 269.0, 383.0])
LOG10_DENSITY_TICKS = np.array([10.0, 11.0, 12.0])
Y_TICKS = np.array([356.0, 82.0])
GAP_TICKS_MM = np.array([0.0, 5.0])
HORIZONTAL_DIGITIZATION_PX = 3.0
VERTICAL_DIGITIZATION_PX = 1.0
PROFILES = (
    {"species": "H", "transition": "Lyman-alpha", "wavelength_nm": 121.6,
     "rgb": (255, 127, 14), "x_roi_inclusive": (210, 261)},
    {"species": "H2", "transition": "Lyman_band", "wavelength_nm": 161.0,
     "rgb": (0, 191, 255), "x_roi_inclusive": (305, 361)},
)
FIELDS = (
    "species", "transition", "wavelength_nm", "pixel_x_median", "pixel_y",
    "matched_exact_color_pixels", "gap_mm", "gap_digitization_lower_mm",
    "gap_digitization_upper_mm", "log10_transition_density_cm3",
    "nominal_transition_density_cm3", "transition_density_digitization_lower_cm3",
    "transition_density_digitization_upper_cm3",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_source_image(pdf: Path) -> tuple[np.ndarray, str]:
    """Reject a different PDF or changed image before using the fixed calibration."""
    digest = sha256(pdf.read_bytes())
    if digest != PDF_SHA256:
        raise ValueError(f"PDF SHA256 mismatch: expected {PDF_SHA256}, got {digest}")
    with pymupdf.open(pdf) as document:
        if len(document) != 203:
            raise ValueError("The calibrated dissertation must contain 203 PDF pages")
        page = document[PDF_PAGE_INDEX]
        if "Figure 4.20" not in page.get_text():
            raise ValueError("PDF page 109 does not identify Figure 4.20")
        if IMAGE_XREF not in {item[0] for item in page.get_images(full=True)}:
            raise ValueError("The calibrated image xref is absent from PDF page 109")
        embedded = document.extract_image(IMAGE_XREF)
        image_bytes = embedded["image"]
    image = Image.open(io.BytesIO(image_bytes))
    if image.size != IMAGE_SIZE or image.mode != "RGB":
        raise ValueError(f"Expected RGB image {IMAGE_SIZE}, got {image.mode} {image.size}")
    # PNG encoding can vary between PyMuPDF builds; the calibrated RGB pixels
    # must match regardless of lossless encoder details.
    if sha256(image.tobytes()) != RGB_PIXELS_SHA256:
        raise ValueError("Embedded image pixels changed; inspect calibration before extraction")
    return np.asarray(image), digest


def recover_profiles(image: np.ndarray) -> tuple[list[dict], dict]:
    density_slope, density_intercept = np.polyfit(X_TICKS, LOG10_DENSITY_TICKS, 1)
    gap_slope, gap_intercept = np.polyfit(Y_TICKS, GAP_TICKS_MM, 1)
    rows: list[dict] = []
    skipped: list[dict] = []
    for profile in PROFILES:
        lower, upper = profile["x_roi_inclusive"]
        for pixel_y in range(87, 352, 5):
            matches = np.flatnonzero(
                np.all(image[pixel_y, lower:upper + 1] == profile["rgb"], axis=1)
            ) + lower
            reason = None
            if len(matches) < 3:
                reason = "fewer_than_three_exact_color_pixels"
            elif len(matches) > 8:
                reason = "overly_broad_exact_color_run"
            elif np.any(np.diff(matches) != 1):
                reason = "multiple_separated_exact_color_runs"
            elif matches[0] == lower or matches[-1] == upper:
                reason = "run_touches_roi_boundary"
            if reason:
                skipped.append({"species": profile["species"], "pixel_y": pixel_y,
                                "reason": reason, "matching_pixels": matches.tolist()})
                continue
            pixel_x = float(np.median(matches))
            log_density = float(density_slope * pixel_x + density_intercept)
            log_estimate = float(abs(density_slope) * HORIZONTAL_DIGITIZATION_PX)
            gap_mm = float(gap_slope * pixel_y + gap_intercept)
            gap_estimate = float(abs(gap_slope) * VERTICAL_DIGITIZATION_PX)
            rows.append({
                "species": profile["species"], "transition": profile["transition"],
                "wavelength_nm": profile["wavelength_nm"], "pixel_x_median": pixel_x,
                "pixel_y": pixel_y, "matched_exact_color_pixels": len(matches),
                "gap_mm": gap_mm, "gap_digitization_lower_mm": gap_mm - gap_estimate,
                "gap_digitization_upper_mm": gap_mm + gap_estimate,
                "log10_transition_density_cm3": log_density,
                "nominal_transition_density_cm3": 10.0 ** log_density,
                "transition_density_digitization_lower_cm3": 10.0 ** (log_density - log_estimate),
                "transition_density_digitization_upper_cm3": 10.0 ** (log_density + log_estimate),
            })
    if not rows:
        raise ValueError("No unambiguous central-line pixels recovered")
    return rows, {
        "density_axis": {"type": "log10", "tick_pixels_x": X_TICKS.tolist(),
                         "tick_log10_density_cm3": LOG10_DENSITY_TICKS.tolist(),
                         "method": "least_squares_linear_fit",
                         "log10_density_per_pixel": float(density_slope),
                         "log10_density_intercept": float(density_intercept),
                         "tick_fit_residual_log10": (
                             LOG10_DENSITY_TICKS - (density_slope * X_TICKS + density_intercept)
                         ).tolist()},
        "gap_axis": {"type": "linear", "tick_pixels_y": Y_TICKS.tolist(),
                     "tick_gap_mm": GAP_TICKS_MM.tolist(),
                     "gap_mm_per_pixel": float(gap_slope),
                     "gap_mm_intercept": float(gap_intercept)},
        "profiles": [{**p, "rgb": list(p["rgb"]),
                      "x_roi_inclusive": list(p["x_roi_inclusive"])} for p in PROFILES],
        "row_sampling": {"first_pixel_y": 87, "last_allowed_pixel_y": 351,
                         "stride_pixels": 5, "last_sampled_pixel_y": 347},
        "accepted_exact_color_run_length_pixels": [3, 8],
        "accepted_run_requirement": "single contiguous run strictly inside specified ROI",
        "digitization_estimates": {
            "horizontal_pixels_plus_minus": HORIZONTAL_DIGITIZATION_PX,
            "vertical_pixels_plus_minus": VERTICAL_DIGITIZATION_PX,
            "log10_density_plus_minus": float(abs(density_slope) * HORIZONTAL_DIGITIZATION_PX),
            "gap_mm_plus_minus": float(abs(gap_slope) * VERTICAL_DIGITIZATION_PX),
            "meaning": "conservative digitization estimates; not physical uncertainty or confidence intervals",
        },
        "skipped_rows": skipped,
    }


def make_qa_plot(image: np.ndarray, rows: list[dict], destination: Path) -> None:
    """Optional visual check; source pixels and recovered values share calibration."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (source_ax, value_ax) = plt.subplots(1, 2, figsize=(10, 5), layout="constrained")
    source_ax.imshow(image, origin="upper")
    source_ax.set(xlim=(70, 435), ylim=(418, 55), xlabel="Embedded image x [pixel]",
                  ylabel="Embedded image y [pixel]", title="Source central curves + recovered pixels")
    for profile in PROFILES:
        selected = [r for r in rows if r["species"] == profile["species"]]
        color = np.asarray(profile["rgb"]) / 255.0
        source_ax.scatter([r["pixel_x_median"] for r in selected],
                          [r["pixel_y"] for r in selected], s=12, facecolors="none",
                          edgecolors="black", linewidths=0.5)
        density = np.asarray([r["nominal_transition_density_cm3"] for r in selected])
        gap = np.asarray([r["gap_mm"] for r in selected])
        lower = np.asarray([r["transition_density_digitization_lower_cm3"] for r in selected])
        upper = np.asarray([r["transition_density_digitization_upper_cm3"] for r in selected])
        value_ax.plot(density, gap, color=color, label=f'{profile["species"]}: {profile["transition"]}')
        value_ax.errorbar(density[::5], gap[::5],
                          xerr=np.vstack((density - lower, upper - density))[:, ::5],
                          yerr=[r["gap_digitization_upper_mm"] - r["gap_mm"]
                                for r in selected][::5],
                          fmt=".", color=color, capsize=2, linewidth=0.8)
    value_ax.set(xscale="log", xlim=(1e10, 2e12), ylim=(0, 5),
                 xlabel=r"Plotted transition density [emitters cm$^{-3}$]",
                 ylabel="Distance from cathode [mm]", title="Recovered profiles")
    value_ax.grid(alpha=0.2, which="both")
    value_ax.legend(loc="lower right")
    fig.suptitle("Figure 4.20, shot 15750: digitized published transition densities")
    fig.supxlabel("Bars indicate digitization estimates only; points on each curve are correlated.", fontsize=9)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pdf", type=Path, required=True, help="Exact dissertation PDF, retained outside version control")
    parser.add_argument("--output-dir", type=Path, default=Path("data/dissertation"))
    parser.add_argument("--qa-plot", type=Path, help="Optional QA image, e.g. runs/recovery/figure_4_20_qa.png")
    args = parser.parse_args()
    image, pdf_digest = load_source_image(args.pdf)
    rows, calibration = recover_profiles(image)
    metadata = {
        "schema_version": 1,
        "artifact": "figure_4_20_profiles.csv",
        "generated_by": "tools/recover_dissertation_profiles.py",
        "source": {
            "author": "Trevor Smith", "publication_year": 2025,
            "title": "Characterizing Desorbed Surface Contaminants Leading to Plasma Formation in Power Feeds of Pulsed Power Accelerators",
            "repository_item_url": "https://deepblue.lib.umich.edu/items/9ccfa886-0b19-4bd8-b3b2-bb685b9ddefc",
            "pdf_url": "https://backend.production.deepblue-documents.lib.umich.edu/server/api/core/bitstreams/24e682a8-896a-410d-85bb-ca4654a9db02/content",
            "pdf_filename": "smtrevor_1.pdf", "pdf_sha256": pdf_digest,
            "pdf_page_count": 203, "printed_page": 84, "pdf_page_1_based": 109,
            "pdf_page_index_0_based": PDF_PAGE_INDEX, "figure": "4.20",
            "embedded_image_xref": IMAGE_XREF,
            "reference_extraction_png_sha256": IMAGE_SHA256,
            "embedded_rgb_pixels_sha256": RGB_PIXELS_SHA256,
            "embedded_image_size_pixels": list(IMAGE_SIZE), "embedded_image_mode": "RGB",
        },
        "experiment": {"facility": "Mykonos", "shot": 15750, "foil_material": "SS304",
                       "foil_thickness_um": 200, "relative_breakdown_timing": "before"},
        "observable": {
            "quantity": "plotted_transition_density",
            "units": "emitters/cm^3", "spatial_coordinate": "distance_from_cathode_mm",
            "physical_interpretation": "Published H Lyman-alpha and H2 161 nm Lyman-band central transition-density profiles across the anode-cathode gap",
            "is_total_ground_state_neutral_density": False,
            "is_raw_detector_data": False,
            "points_are_independent_observations": False,
            "physical_uncertainty_bands_recovered": False,
            "limitations": [
                "Samples from the same raster curve are correlated digitized constraints, not independent raw observations.",
                "The spatial coordinate spans the electrode gap, not the foil width.",
                "Digitization bounds cover estimated pixel reading precision only; no statistical confidence interpretation.",
                "Published shaded uncertainty bands and dashed density bounds are deliberately not extracted.",
                "Conversion to ground-state species density requires a separately validated population/emission model and experimental diagnostic normalization.",
            ],
        },
        "calibration": calibration,
        "row_count": len(rows),
        "row_count_per_species": {p["species"]: sum(r["species"] == p["species"] for r in rows) for p in PROFILES},
        "csv_columns": {
            "pixel_x_median": "Median exact-color x coordinate in the native embedded image.",
            "pixel_y": "Sampled y coordinate, increasing downwards in the native embedded image.",
            "matched_exact_color_pixels": "Number of pixels in the accepted contiguous central-line run.",
            "gap_mm": "Distance from cathode, from the linear image y-axis calibration.",
            "log10_transition_density_cm3": "Base-10 logarithm of the numerical transition density expressed in emitters/cm^3.",
            "nominal_transition_density_cm3": "Central plotted transition density obtained from the logarithmic image x axis.",
            "digitization_lower_upper_columns": "Bounds from +/-3 x pixels or +/-1 y pixel; estimates only, not physical errors.",
            "wavelength_nm": "Approximate displayed transition wavelength; 121.6 nm for H Lyman-alpha, 161 nm for H2 Lyman band.",
        },
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / metadata["artifact"]
    with csv_path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    json_path = args.output_dir / "figure_4_20_profiles.json"
    json_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if args.qa_plot:
        make_qa_plot(image, rows, args.qa_plot)
    print(json.dumps({"csv": str(csv_path), "metadata": str(json_path),
                      "row_count": len(rows), "skipped_rows": len(calibration["skipped_rows"]),
                      "qa_plot": str(args.qa_plot) if args.qa_plot else None}))


if __name__ == "__main__":
    main()
