#!/usr/bin/env python3
"""Recover the temperature curve in dissertation Figure 4.30 from fixed pixels.

The source is the embedded 975 x 975 RGB raster on printed page 98 / PDF page
123 of Trevor Johannes Smith's 2025 dissertation. This is a 40 ns legacy
simulation profile at 12.5 micrometres depth, not a measurement or a surface
temperature. The legend explicitly labels temperature in eV. Convert k_B T in
eV to kelvin using the exact SI elementary charge and Boltzmann constant.

Inspected calibration, 2026-10-18 (native pixel origin at upper left):
* x tick centres 122.5, 258.5, 395.5, 531.5, 668.5, 804.5, 941.5 map to
  -3, -2, -1, 0, 1, 2, 3 mm from the foil middle.
* y tick centres 98.5, 254.5, 410.5, 566.5, 721.5, 877.5 map to shared signal
  heights 0.2, 0, -0.2, -0.4, -0.6, -0.8. Temperature uses this axis in eV.
* Temperature is exact RGB (255, 140, 0). Sample columns 194 through 869 at
  stride 5, plus column 531 nearest the centre. Within y=90:260 inclusive,
  accept one contiguous run of 2--8 exact-colour pixels and use its median.
* Exclude the near-vertical foil-edge jumps and every exterior zero segment.
  Do not extrapolate a boundary peak, interpolate missing columns or impose
  left/right symmetry. The steepest actual edge values remain unrecovered.

Coordinate +/-2 pixel estimates describe digitization, not experimental or
simulation uncertainty. Nearby raster samples are correlated.

Dependencies: numpy, pillow, pymupdf. Optional --qa-plot needs matplotlib.
Example:
  python tools/recover_dissertation_heating.py --pdf /path/to/smtrevor_1.pdf \
      --qa-plot runs/recovery/figure_4_30_heating_qa.png
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
PDF_PAGE_INDEX = 122
IMAGE_XREF = 2890
IMAGE_SIZE = (975, 975)
TEMPERATURE_RGB = (255, 140, 0)
X_TICKS_PX = np.array([122.5, 258.5, 395.5, 531.5, 668.5, 804.5, 941.5])
X_TICKS_MM = np.arange(-3.0, 4.0)
Y_TICKS_PX = np.array([98.5, 254.5, 410.5, 566.5, 721.5, 877.5])
Y_TICKS_EV = np.array([0.2, 0.0, -0.2, -0.4, -0.6, -0.8])
ELEMENTARY_CHARGE_J_PER_EV = 1.602176634e-19
BOLTZMANN_J_PER_K = 1.380649e-23
K_PER_EV = ELEMENTARY_CHARGE_J_PER_EV / BOLTZMANN_J_PER_K
DIGITIZATION_PX = 2.0
FIELDS = (
    "pixel_x", "pixel_y_median", "exact_color_run_first_y", "exact_color_run_last_y",
    "matched_exact_color_pixels", "width_offset_mm", "width_offset_digitization_lower_mm",
    "width_offset_digitization_upper_mm", "temperature_eV", "temperature_digitization_lower_eV",
    "temperature_digitization_upper_eV", "temperature_K", "temperature_digitization_lower_K",
    "temperature_digitization_upper_K",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_image(pdf: Path) -> tuple[np.ndarray, str]:
    """Check the source identity before applying this one-raster calibration."""
    digest = sha256(pdf.read_bytes())
    if digest != PDF_SHA256:
        raise ValueError(f"PDF SHA256 mismatch: expected {PDF_SHA256}, got {digest}")
    with pymupdf.open(pdf) as document:
        if len(document) != 203:
            raise ValueError("The calibrated source must contain 203 pages")
        page = document[PDF_PAGE_INDEX]
        if "Figure 4.30" not in page.get_text():
            raise ValueError("PDF page 123 does not identify Figure 4.30")
        if IMAGE_XREF not in {entry[0] for entry in page.get_images(full=True)}:
            raise ValueError("The calibrated image xref is absent from PDF page 123")
        encoded = document.extract_image(IMAGE_XREF)["image"]
    image = Image.open(io.BytesIO(encoded))
    if image.mode != "RGB" or image.size != IMAGE_SIZE:
        raise ValueError(f"Expected RGB {IMAGE_SIZE}, got {image.mode} {image.size}")
    if sha256(image.tobytes()) != RGB_PIXELS_SHA256:
        raise ValueError("Native pixels changed; inspect and recalibrate before extraction")
    # The decoded pixel hash is authoritative: lossless PNG encoder bytes can
    # differ across library versions even when every calibrated pixel is equal.
    return np.asarray(image), sha256(encoded)


def recover(image: np.ndarray) -> tuple[list[dict], dict]:
    x_slope, x_intercept = np.polyfit(X_TICKS_PX, X_TICKS_MM, 1)
    t_slope, t_intercept = np.polyfit(Y_TICKS_PX, Y_TICKS_EV, 1)
    dx = abs(float(x_slope)) * DIGITIZATION_PX
    dt = abs(float(t_slope)) * DIGITIZATION_PX
    columns = sorted(set(range(194, 870, 5)) | {531})
    rows, skipped = [], []
    for pixel_x in columns:
        matches = np.flatnonzero(
            np.all(image[90:261, pixel_x] == TEMPERATURE_RGB, axis=1)
        ) + 90
        reason = None
        if len(matches) < 2:
            reason = "fewer_than_two_exact_color_pixels"
        elif len(matches) > 8:
            reason = "overly_broad_or_near_vertical_run"
        elif np.any(np.diff(matches) != 1):
            reason = "multiple_separated_exact_color_runs"
        elif matches[0] == 90 or matches[-1] == 260:
            reason = "run_touches_roi_boundary"
        if reason:
            skipped.append({"pixel_x": pixel_x, "reason": reason,
                            "matching_pixels_y": matches.tolist()})
            continue
        pixel_y = float(np.median(matches))
        x_mm = float(x_slope * pixel_x + x_intercept)
        t_ev = float(t_slope * pixel_y + t_intercept)
        if t_ev <= 0:
            raise ValueError("Accepted temperature must be inside the nonzero foil profile")
        rows.append({
            "pixel_x": pixel_x, "pixel_y_median": pixel_y,
            "exact_color_run_first_y": int(matches[0]),
            "exact_color_run_last_y": int(matches[-1]),
            "matched_exact_color_pixels": len(matches),
            "width_offset_mm": x_mm,
            "width_offset_digitization_lower_mm": x_mm - dx,
            "width_offset_digitization_upper_mm": x_mm + dx,
            "temperature_eV": t_ev,
            "temperature_digitization_lower_eV": t_ev - dt,
            "temperature_digitization_upper_eV": t_ev + dt,
            "temperature_K": t_ev * K_PER_EV,
            "temperature_digitization_lower_K": (t_ev - dt) * K_PER_EV,
            "temperature_digitization_upper_K": (t_ev + dt) * K_PER_EV,
        })
    if not rows:
        raise ValueError("No unambiguous temperature line was recovered")
    calibration = {
        "inspected_date": "2026-10-18", "coordinate_origin": "upper_left_native_raster",
        "method": "least_squares_linear_fit_to_inspected_tick_centres",
        "x_axis": {"label": "Distance from Foil-Middle Location [mm]",
                   "tick_pixels_x": X_TICKS_PX.tolist(), "tick_values_mm": X_TICKS_MM.tolist(),
                   "mm_per_pixel": float(x_slope), "intercept_mm": float(x_intercept),
                   "tick_residuals_mm": (X_TICKS_MM - np.polyval([x_slope, x_intercept], X_TICKS_PX)).tolist()},
        "y_axis": {"label": "Signal Height", "temperature_legend": "Temperature [eV]",
                   "tick_pixels_y": Y_TICKS_PX.tolist(), "tick_values_eV": Y_TICKS_EV.tolist(),
                   "eV_per_pixel": float(t_slope), "intercept_eV": float(t_intercept),
                   "tick_residuals_eV": (Y_TICKS_EV - np.polyval([t_slope, t_intercept], Y_TICKS_PX)).tolist()},
        "selection": {"exact_RGB": list(TEMPERATURE_RGB), "y_roi_inclusive": [90, 260],
                      "x_first": 194, "x_last": 869, "x_stride": 5,
                      "additional_x_columns": [531], "accepted_contiguous_run_pixels": [2, 8],
                      "line_position": "median_of_single_contiguous_exact_color_run",
                      "interpolation": False, "extrapolation": False, "symmetry_enforced": False},
        "digitization_estimate": {"coordinate_pixels_plus_minus": DIGITIZATION_PX,
                                  "width_offset_mm_plus_minus": dx, "temperature_eV_plus_minus": dt,
                                  "temperature_K_plus_minus": dt * K_PER_EV,
                                  "meaning": "conservative coordinate readout estimate, not physical uncertainty or a confidence interval"},
        "skipped_columns": skipped,
    }
    return rows, calibration


def make_qa(image: np.ndarray, rows: list[dict], path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, (source_ax, profile_ax) = plt.subplots(1, 2, figsize=(12, 5.6), layout="constrained")
    source_ax.imshow(image, origin="upper")
    source_ax.scatter([r["pixel_x"] for r in rows], [r["pixel_y_median"] for r in rows],
                      s=6, color="black", label="Recovered temperature centre")
    source_ax.set_title("Native Fig. 4.30 pixels and recovered samples")
    source_ax.set_xlim(0, 975)
    source_ax.set_ylim(975, 0)
    source_ax.axis("off")
    source_ax.legend(fontsize=7)
    x = np.array([r["width_offset_mm"] for r in rows])
    t = np.array([r["temperature_K"] for r in rows])
    low = np.array([r["temperature_digitization_lower_K"] for r in rows])
    high = np.array([r["temperature_digitization_upper_K"] for r in rows])
    profile_ax.fill_between(x, low, high, color="orange", alpha=0.25,
                            label="±2 pixel vertical readout estimate")
    profile_ax.plot(x, t, ".", ms=3, color="darkorange", label="Recovered plotted temperature")
    profile_ax.set(xlabel="Distance from foil middle [mm]", ylabel="Temperature [K]",
                   title="Legacy simulation, 40 ns, 12.5 µm depth")
    profile_ax.grid(alpha=0.2)
    profile_ax.legend(fontsize=8)
    figure.suptitle("Published simulation profile; foil-edge jumps excluded; no surface-temperature inference", fontsize=10)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=150)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--csv", type=Path, default=ROOT / "data/dissertation/figure_4_30_heating.csv")
    parser.add_argument("--metadata", type=Path, default=ROOT / "data/dissertation/figure_4_30_heating.json")
    parser.add_argument("--qa-plot", type=Path)
    args = parser.parse_args()
    image, encoded_hash = load_image(args.pdf)
    rows, calibration = recover(image)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format(value, ".12g") if isinstance(value, float) else value
                             for key, value in row.items()})
    nearest_middle = min(rows, key=lambda row: abs(row["width_offset_mm"]))
    metadata = {
        "schema_version": 1, "status": "digitized_published_legacy_simulation",
        "source": {"title": "Characterizing Desorbed Surface Contaminants Leading to Plasma Formation in Power Feeds of Pulsed Power Accelerators",
                   "author": "Trevor Johannes Smith", "year": 2025,
                   "item_url": "https://deepblue.lib.umich.edu/items/9ccfa886-0b19-4bd8-b3b2-bb685b9ddefc",
                   "pdf_url": "https://backend.production.deepblue-documents.lib.umich.edu/server/api/core/bitstreams/24e682a8-896a-410d-85bb-ca4654a9db02/content",
                   "pdf_sha256": PDF_SHA256, "figure": "4.30", "printed_page": 98, "pdf_page": 123,
                   "pdf_page_index_zero_based": PDF_PAGE_INDEX, "embedded_image_xref": IMAGE_XREF,
                   "native_image_size_pixels": list(IMAGE_SIZE), "native_image_mode": "RGB",
                   "extracted_png_sha256": encoded_hash, "reference_extracted_png_sha256": REFERENCE_PNG_SHA256,
                   "decoded_RGB_pixels_sha256": RGB_PIXELS_SHA256,
                   "hash_validation": "PDF and decoded RGB pixels are mandatory; lossless PNG encoding may vary"},
        "simulation_context": {"code_label_in_dissertation": "Allegra", "drive_shot": 17133,
                               "foil_material": "SS304", "foil_thickness_m": 200e-6,
                               "foil_thickness_source": "section 4.5, printed page 93",
                               "time_s": 40e-9, "depth_from_AK_gap_m": 12.5e-6,
                               "foil_width_m": 5e-3, "AK_gap_m": 5e-3,
                               "dielectric_support_included": False,
                               "observation_type": "published simulation lineout, not experimental measurement",
                               "surface_temperature": False, "FLASH_result": False},
        "temperature_conversion": {"interpretation": "legend Temperature [eV] is k_B*T expressed in electronvolts",
                                   "equation": "T_K = temperature_eV * elementary_charge_J_per_eV / Boltzmann_J_per_K",
                                   "elementary_charge_J_per_eV": ELEMENTARY_CHARGE_J_PER_EV,
                                   "Boltzmann_J_per_K": BOLTZMANN_J_PER_K, "K_per_eV": K_PER_EV,
                                   "Celsius_offset_applied": False},
        "calibration": calibration,
        "output": {"csv_filename": args.csv.name, "csv_sha256": sha256(args.csv.read_bytes()),
                   "columns": list(FIELDS), "accepted_samples": len(rows),
                   "width_offset_extent_mm": [rows[0]["width_offset_mm"], rows[-1]["width_offset_mm"]],
                   "temperature_extent_K": [min(r["temperature_K"] for r in rows), max(r["temperature_K"] for r in rows)],
                   "nearest_middle_sample": nearest_middle,
                   "generator": "tools/recover_dissertation_heating.py",
                   "generator_sha256": sha256(Path(__file__).read_bytes())},
        "limits": [
            "Raster digitization of a legacy simulation curve, not recovered original simulation arrays or measured data.",
            "The lineout is 12.5 micrometres inside the foil at 40 ns, not the gap-facing surface or a temperature time history.",
            "Foil-edge near-vertical jumps and exterior zero segments are excluded; exact edge maxima and boundary values are not recovered.",
            "Coordinate readout bounds do not cover material, simulation, timing, geometry or experimental uncertainties; neighbouring samples are correlated.",
            "Section 4.5 p93 describes a 700 C surface condition at 40 ns. That text does not calibrate this subsurface eV-labelled curve.",
            "The p93/p95 wording conflict about the 60 ns surface/melt temperature remains unresolved; this 40 ns subsurface extraction does not settle it.",
            "No smoothing, symmetry constraint, interpolation or extrapolation is applied; decimal precision records the pixel calibration, not physical accuracy.",
            "Magnetic-field and current-density curves sharing the signal-height axis are not recovered by this tool.",
        ],
    }
    args.metadata.parent.mkdir(parents=True, exist_ok=True)
    args.metadata.write_text(json.dumps(metadata, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    if args.qa_plot:
        make_qa(image, rows, args.qa_plot)
    print(f"Recovered {len(rows)} temperature samples: {args.csv}")
    print(f"Metadata: {args.metadata}")
    print(f"Nearest-middle temperature: {nearest_middle['temperature_K']:.1f} K; vertical readout estimate ±{calibration['digitization_estimate']['temperature_K_plus_minus']:.1f} K")


if __name__ == "__main__":
    main()
