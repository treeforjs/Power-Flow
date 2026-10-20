#!/usr/bin/env python3
"""Recover interior temperature profiles from dissertation Figs. 4.28 and 4.29.

The lower-right panels contain foil-edge (dark orange) and foil-middle (peru)
temperature profiles at 40 and 60 ns. Both temperature axes are explicitly K;
the vertical coordinate is distance from the middle of the AK gap in mm.
This tool does not use the eV axis of Fig. 4.30.

Inspected native-raster calibration, 2026-10-18; pixel origin upper left:
* Native RGB images are 975 x 975 pixels, printed pages 94/96, PDF pages 119/121.
* Temperature tick centres x=597.5,700.5,802.5,905.5 map to 0,1000,2000,3000 K.
* Gap-centre tick centres y=900.5,851.5,802.5,752.5,703.5,653.5,604.5,555.5,505.5
  map to 2.40,2.45,2.50,2.55,2.60,2.65,2.70,2.75,2.80 mm.
* Exact RGB (255,140,0) is foil-edge; (205,133,63) is foil-middle. Inspect only
  x=610:950 and rows 795,790,...,610, below the legend and away from zero-valued
  exterior segments and boundary ramps. Accept one contiguous 2--8 pixel run.
* The nominal INITIAL cathode gap face is +2.5 mm from the gap centre, with
  nominal back face +2.7 mm, from the reported 5 mm gap and 200 um foil.
  Actual moving material boundaries, mesh and material fractions are unavailable.
  Export the absolute plotted coordinate and a clearly labeled nominal depth;
  do not identify either recovered endpoint as a physical surface.

No smoothing, interpolation, extrapolation, surface value or temperature time
history is constructed. These are correlated raster samples from legacy
simulation figures, not measurements or FLASH output. Readout bounds combine
an estimated +/-2 native pixels with maximum tick-calibration residual; they
are not physical uncertainty or confidence intervals.

Dependencies: numpy, pillow, pymupdf. Optional --qa-dir needs matplotlib.
Example:
  python tools/recover_dissertation_depth_heating.py --pdf /path/to/smtrevor_1.pdf \
      --qa-dir runs/recovery
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
X_TICKS_PX = np.array([597.5, 700.5, 802.5, 905.5])
TEMPERATURE_TICKS_K = np.array([0.0, 1000.0, 2000.0, 3000.0])
Y_TICKS_PX = np.array([900.5, 851.5, 802.5, 752.5, 703.5, 653.5, 604.5, 555.5, 505.5])
GAP_COORDINATE_TICKS_MM = np.array([2.40, 2.45, 2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80])
INITIAL_GAP_FACE_MM = 2.5
INITIAL_BACK_FACE_MM = 2.7
COORDINATE_READOUT_PX = 2.0
PROFILES = (
    {"id": "foil_edge", "legend": "Foil-edge location (cf. Fig. 4.27)", "rgb": (255, 140, 0)},
    {"id": "foil_middle", "legend": "Foil-middle location (cf. Fig. 4.27)", "rgb": (205, 133, 63)},
)
FIGURES = (
    {"figure": "4.28", "stem": "figure_4_28_temperature", "printed_page": 94,
     "pdf_page": 119, "image_xref": 2847, "time_s": 40e-9,
     "reference_png_sha256": "1f9c939c8b72947da2629ccad92c42dda04503719c3cdb89fe0e921bc5ce5488",
     "rgb_sha256": "458b3ced8db4501db3fcbcc90939ce5f336a7ec83eb97ee259c66f7112bf07c8"},
    {"figure": "4.29", "stem": "figure_4_29_temperature", "printed_page": 96,
     "pdf_page": 121, "image_xref": 2870, "time_s": 60e-9,
     "reference_png_sha256": "5ae90775c30a8fad657f74d5813c150b2a4284d1a0fd8075ebc4fcd21b88e89a",
     "rgb_sha256": "ac25b035329213a8d0a354b8496b1674c2be24ceaeff65d95a48ef2e6376b34f"},
)
FIELDS = (
    "profile", "time_s", "pixel_y", "pixel_x_median", "exact_color_run_first_x",
    "exact_color_run_last_x", "matched_exact_color_pixels", "gap_center_distance_mm",
    "gap_center_distance_digitization_lower_mm", "gap_center_distance_digitization_upper_mm",
    "nominal_depth_from_initial_gap_face_um", "nominal_depth_digitization_lower_um",
    "nominal_depth_digitization_upper_um", "temperature_K",
    "temperature_digitization_lower_K", "temperature_digitization_upper_K",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_image(document: pymupdf.Document, specification: dict) -> tuple[np.ndarray, str]:
    page = document[specification["pdf_page"] - 1]
    if f"Figure {specification['figure']}" not in page.get_text():
        raise ValueError("The source page does not identify the calibrated figure")
    xref = specification["image_xref"]
    if xref not in {entry[0] for entry in page.get_images(full=True)}:
        raise ValueError("Calibrated embedded image is absent from the source page")
    encoded = document.extract_image(xref)["image"]
    image = Image.open(io.BytesIO(encoded))
    if image.mode != "RGB" or image.size != (975, 975):
        raise ValueError(f"Expected native RGB 975 x 975 image, got {image.mode} {image.size}")
    if sha256(image.tobytes()) != specification["rgb_sha256"]:
        raise ValueError("Native image pixels changed; inspect and recalibrate before extraction")
    return np.asarray(image), sha256(encoded)


def recover(image: np.ndarray, specification: dict) -> tuple[list[dict], dict]:
    temperature_fit = np.polyfit(X_TICKS_PX, TEMPERATURE_TICKS_K, 1)
    coordinate_fit = np.polyfit(Y_TICKS_PX, GAP_COORDINATE_TICKS_MM, 1)
    temperature_residual = TEMPERATURE_TICKS_K - np.polyval(temperature_fit, X_TICKS_PX)
    coordinate_residual = GAP_COORDINATE_TICKS_MM - np.polyval(coordinate_fit, Y_TICKS_PX)
    temperature_bound = abs(float(temperature_fit[0])) * COORDINATE_READOUT_PX + float(max(abs(temperature_residual)))
    coordinate_bound = abs(float(coordinate_fit[0])) * COORDINATE_READOUT_PX + float(max(abs(coordinate_residual)))
    rows, skipped = [], []
    for profile in PROFILES:
        for pixel_y in range(795, 609, -5):
            matches = np.flatnonzero(np.all(image[pixel_y, 610:951] == profile["rgb"], axis=1)) + 610
            reason = None
            if len(matches) < 2:
                reason = "fewer_than_two_exact_color_pixels"
            elif len(matches) > 8:
                reason = "overly_broad_or_boundary_run"
            elif np.any(np.diff(matches) != 1):
                reason = "multiple_separated_exact_color_runs"
            elif matches[0] == 610 or matches[-1] == 950:
                reason = "run_touches_ROI_boundary"
            if reason:
                skipped.append({"profile": profile["id"], "pixel_y": pixel_y,
                                "reason": reason, "matching_pixels_x": matches.tolist()})
                continue
            pixel_x = float(np.median(matches))
            temperature = float(np.polyval(temperature_fit, pixel_x))
            coordinate = float(np.polyval(coordinate_fit, pixel_y))
            nominal_depth = (coordinate - INITIAL_GAP_FACE_MM) * 1000.0
            if temperature <= 0 or not INITIAL_GAP_FACE_MM < coordinate < INITIAL_BACK_FACE_MM:
                raise ValueError("Accepted point is outside the positive-temperature nominal foil interior")
            rows.append({
                "profile": profile["id"], "time_s": specification["time_s"],
                "pixel_y": pixel_y, "pixel_x_median": pixel_x,
                "exact_color_run_first_x": int(matches[0]), "exact_color_run_last_x": int(matches[-1]),
                "matched_exact_color_pixels": len(matches), "gap_center_distance_mm": coordinate,
                "gap_center_distance_digitization_lower_mm": coordinate - coordinate_bound,
                "gap_center_distance_digitization_upper_mm": coordinate + coordinate_bound,
                "nominal_depth_from_initial_gap_face_um": nominal_depth,
                "nominal_depth_digitization_lower_um": nominal_depth - coordinate_bound * 1000.0,
                "nominal_depth_digitization_upper_um": nominal_depth + coordinate_bound * 1000.0,
                "temperature_K": temperature,
                "temperature_digitization_lower_K": temperature - temperature_bound,
                "temperature_digitization_upper_K": temperature + temperature_bound,
            })
    if not rows:
        raise ValueError("No unambiguous interior temperature profile recovered")
    calibration = {
        "inspected_date": "2026-10-18", "coordinate_origin": "upper_left_native_raster",
        "panel": "lower_right", "method": "least_squares_linear_fit_to_inspected_tick_centres",
        "temperature_axis": {"label": "Temperature [K]", "tick_pixels_x": X_TICKS_PX.tolist(),
                             "tick_values_K": TEMPERATURE_TICKS_K.tolist(),
                             "K_per_pixel": float(temperature_fit[0]), "intercept_K": float(temperature_fit[1]),
                             "tick_residuals_K": temperature_residual.tolist(), "eV_conversion_applied": False},
        "gap_coordinate_axis": {"label": "Distance from Middle of AK Gap [mm]",
                                "tick_pixels_y": Y_TICKS_PX.tolist(),
                                "tick_values_mm": GAP_COORDINATE_TICKS_MM.tolist(),
                                "mm_per_pixel": float(coordinate_fit[0]), "intercept_mm": float(coordinate_fit[1]),
                                "tick_residuals_mm": coordinate_residual.tolist()},
        "selection": {"profiles": [{**p, "rgb": list(p["rgb"])} for p in PROFILES],
                      "x_ROI_inclusive": [610, 950], "sampled_rows": list(range(795, 609, -5)),
                      "accepted_contiguous_run_pixels": [2, 8],
                      "line_position": "median_of_single_contiguous_exact_color_run",
                      "interpolation": False, "extrapolation": False, "smoothing": False,
                      "boundary_or_vacuum_values_recovered": False},
        "digitization_bounds": {"coordinate_readout_pixels_plus_minus": COORDINATE_READOUT_PX,
                                "combination": "absolute_calibration_slope_times_pixel_estimate_plus_maximum_absolute_tick_fit_residual",
                                "temperature_K_plus_minus": temperature_bound,
                                "gap_center_distance_mm_plus_minus": coordinate_bound,
                                "nominal_depth_um_plus_minus": coordinate_bound * 1000.0,
                                "meaning": "coordinate readout envelope; not physical uncertainty or confidence intervals"},
        "displayed_melt_reference": {"legend": "SS304 Melt", "inspected_median_pixel_x": 767.0,
                                     "line_RGB": [128, 128, 128],
                                     "temperature_K": float(np.polyval(temperature_fit, 767.0)),
                                     "temperature_digitization_plus_minus_K": temperature_bound,
                                     "role": "figure annotation only; not an independently sourced material EOS or adopted melting criterion"},
        "skipped_rows": skipped,
    }
    return rows, calibration


def make_qa(image: np.ndarray, rows: list[dict], calibration: dict, specification: dict, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, (source_ax, profile_ax) = plt.subplots(1, 2, figsize=(11, 5.4), layout="constrained")
    source_ax.imshow(image, origin="upper")
    source_ax.set(xlim=(490, 975), ylim=(970, 490), title=f"Native Fig. {specification['figure']} thermal panel")
    source_ax.axis("off")
    for profile, marker in zip(PROFILES, ("o", "s")):
        selected = [r for r in rows if r["profile"] == profile["id"]]
        source_ax.scatter([r["pixel_x_median"] for r in selected], [r["pixel_y"] for r in selected],
                          s=7, color="black", marker=marker)
        temperature = np.array([r["temperature_K"] for r in selected])
        depth = np.array([r["nominal_depth_from_initial_gap_face_um"] for r in selected])
        low = np.array([r["temperature_digitization_lower_K"] for r in selected])
        high = np.array([r["temperature_digitization_upper_K"] for r in selected])
        color = np.array(profile["rgb"]) / 255.0
        profile_ax.fill_betweenx(depth, low, high, alpha=0.2, color=color)
        profile_ax.plot(temperature, depth, marker, ms=3, color=color, label=profile["id"].replace("_", " "))
    profile_ax.axvline(calibration["displayed_melt_reference"]["temperature_K"], color="gray", ls="--",
                       label="Displayed SS304 melt annotation")
    profile_ax.set(xlabel="Temperature [K]", ylabel="Nominal depth from initial gap face [µm]",
                   title=f"{specification['time_s'] * 1e9:g} ns; interior samples only")
    profile_ax.grid(alpha=0.2)
    profile_ax.legend(fontsize=8)
    figure.suptitle("Legacy simulation; actual material boundaries and physical uncertainty unrecovered", fontsize=10)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=150)
    plt.close(figure)


def save(rows: list[dict], calibration: dict, specification: dict, encoded_hash: str, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / (specification["stem"] + ".csv")
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format(value, ".12g") if isinstance(value, float) else value for key, value in row.items()})
    metadata = {
        "schema_version": 1, "status": "partially_digitized_published_legacy_simulation",
        "source": {"author": "Trevor Johannes Smith", "year": 2025,
                   "title": "Characterizing Desorbed Surface Contaminants Leading to Plasma Formation in Power Feeds of Pulsed Power Accelerators",
                   "item_url": "https://deepblue.lib.umich.edu/items/9ccfa886-0b19-4bd8-b3b2-bb685b9ddefc",
                   "pdf_sha256": PDF_SHA256, "figure": specification["figure"],
                   "printed_page": specification["printed_page"], "pdf_page": specification["pdf_page"],
                   "pdf_page_index_zero_based": specification["pdf_page"] - 1,
                   "embedded_image_xref": specification["image_xref"],
                   "native_image_size_pixels": [975, 975], "native_image_mode": "RGB",
                   "reference_extracted_PNG_sha256": specification["reference_png_sha256"],
                   "extracted_PNG_sha256": encoded_hash, "decoded_RGB_pixels_sha256": specification["rgb_sha256"],
                   "hash_validation": "PDF and decoded native RGB pixels mandatory; lossless PNG encoding may vary"},
        "simulation_context": {"code_label_in_dissertation": "Allegra", "drive_shot": 17133,
                               "time_s": specification["time_s"], "foil_material": "SS304",
                               "initial_foil_thickness_m": 200e-6, "initial_gap_m": 5e-3,
                               "dielectric_support_included": False,
                               "context_source": "section 4.5 and Figure 4.27, printed page 93 / PDF page 118",
                               "observation_type": "published simulation lineout, not experimental data",
                               "FLASH_result": False, "temperature_time_history": False,
                               "drive_15210_used": False},
        "spatial_coordinates": {"absolute_coordinate": "positive cathode side, distance from middle of AK gap in mm",
                                "nominal_initial_gap_face_mm": INITIAL_GAP_FACE_MM,
                                "nominal_initial_back_face_mm": INITIAL_BACK_FACE_MM,
                                "nominal_depth_equation": "1000 * (gap_center_distance_mm - 2.5)",
                                "nominal_depth_units": "micrometres",
                                "boundary_evidence": "5 mm initial AK gap and 200 micrometre foil in section 4.5; geometry in Figure 4.27",
                                "actual_time_dependent_material_boundary_recovered": False,
                                "boundary_motion_or_material_fraction_data_available": False,
                                "profile_lateral_locations": "foil-edge and foil-middle labels as in Figure 4.27; precise transverse sampling offsets not numerically specified"},
        "calibration": calibration,
        "output": {"csv_filename": csv_path.name, "csv_sha256": sha256(csv_path.read_bytes()),
                   "columns": list(FIELDS), "samples": len(rows),
                   "samples_by_profile": {p["id"]: sum(r["profile"] == p["id"] for r in rows) for p in PROFILES},
                   "nominal_depth_extent_um": [min(r["nominal_depth_from_initial_gap_face_um"] for r in rows),
                                               max(r["nominal_depth_from_initial_gap_face_um"] for r in rows)],
                   "temperature_extent_K_by_profile": {p["id"]: [min(r["temperature_K"] for r in rows if r["profile"] == p["id"]),
                                                                  max(r["temperature_K"] for r in rows if r["profile"] == p["id"])] for p in PROFILES},
                   "generator": "tools/recover_dissertation_depth_heating.py",
                   "generator_sha256": sha256(Path(__file__).read_bytes())},
        "source_wording_limits": [
            {"source": "section 4.5 printed page 93", "statement": "40 ns foil-middle gap-facing surface at 700 C; 60 ns foil-middle gap-facing surface at 1600 C"},
            {"source": "section 4.5 printed page 95", "statement": "60 ns foil-edge exceeds melt while foil-middle remains well below melt"},
            {"source": "Figures 4.28 and 4.29 thermal panels", "statement": "Axes explicitly label Temperature [K]; these labels determine the recovered temperature units, with no eV or Celsius conversion"},
        ],
        "limits": [
            "Only the two reliable interior temperature curves are digitized; current, magnetic field and current-density panels are not recovered here.",
            "No original simulation arrays, grid, EOS, material fractions or material-boundary trajectories were recovered.",
            "Nominal depth is measured from the initial 2.5 mm gap-face coordinate; it is not depth from an inferred evolving surface.",
            "Vacuum/exterior zero segments, boundary ramps, sharp surface peaks and the legend region are excluded; no surface temperatures or boundary values are assigned.",
            "The foil-edge profile label does not specify its precise transverse distance from the geometrical edge; it need not equal the maximum in Figure 4.30.",
            "Pixel readout bounds and tick-fit residuals do not include uncertainty in source materials, simulation, geometry, current drive or timing.",
            "Samples on each raster curve are correlated; numerical decimal precision records calibration rather than physical accuracy.",
            "The two figures are snapshots from a simulation driven by shot 17133, not a time-resolved experimental temperature series or a thermal prediction for shot 15210.",
            "The plotted SS304 melt line is retained only as a figure annotation, without adopting it as a material property.",
            "Narrative temperature/melt wording is not used to recalibrate or override the Kelvin axes; the surface wording conflict remains unresolved.",
        ],
    }
    json_path = output_dir / (specification["stem"] + ".json")
    json_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Recovered {len(rows)} interior temperature samples: {csv_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/dissertation")
    parser.add_argument("--qa-dir", type=Path)
    args = parser.parse_args()
    if sha256(args.pdf.read_bytes()) != PDF_SHA256:
        raise ValueError("PDF differs from the inspected dissertation")
    with pymupdf.open(args.pdf) as document:
        if len(document) != 203:
            raise ValueError("Expected 203 PDF pages")
        # Extract and validate both native images before creating any output.
        sources = [(specification, *load_image(document, specification)) for specification in FIGURES]
    for specification, image, encoded_hash in sources:
        rows, calibration = recover(image, specification)
        save(rows, calibration, specification, encoded_hash, args.output_dir)
        if args.qa_dir:
            make_qa(image, rows, calibration, specification,
                    args.qa_dir / (specification["stem"] + "_qa.png"))


if __name__ == "__main__":
    main()
