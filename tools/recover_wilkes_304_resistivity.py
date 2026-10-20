#!/usr/bin/env python3
"""Recover ten unambiguous specimen-2 markers from Wilkes 1997 Figure 1.

Only the 0.203 mm 304 foil has an absolute-resistivity figure. Filled diamonds
identify heating; enclosed white diamonds identify cooling. Three overlapping
marker clusters are excluded, and Figure 2 percentage deviations are not
converted into absolute resistivity for the 0.305 mm foil.

Fixed native-pixel ROIs were visually inspected on 2026-10-18. Original PDF,
page content and decoded monochrome rasters are guarded. A 2D affine tick fit
accounts for scan skew. No regression/reference line, interpolated point or
temperature nominalization is recovered. Pixel envelopes are readout estimates,
not physical measurement errors. Dependencies: numpy, scipy; PyMuPDF
for PDF extraction, matplotlib only for --qa-overlay.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from scipy import ndimage


ROOT = Path(__file__).resolve().parents[1]
PDF_SHA256 = "5738cdc4407bf6d09da8a8f2379105ff7a32c6d77c6fbb495d5391be67d154e7"
PAGE_RECORDS = (
    (4, 22, (2528, 3327), "c93f29f295f1f0c5961c8fb2dab4439776a90e8f34543a3c152cb83b81e61be3", "733f674c1fa08847c0c603a90bee36172b9cbbaca3e6b8a7e89543b7e1dc5335"),
    (6, 32, (2528, 3322), "143435ad898153d555c16cdb731cd1516f671ef32971797f628aa5c722e3ff1b", "81c0100d553ada5509f26a78d3fa4a9231afb430ff7fd2d39558b0fb7f2a3bf5"),
    (7, 37, (2528, 3328), "c93f29f295f1f0c5961c8fb2dab4439776a90e8f34543a3c152cb83b81e61be3", "417942d97b92488f967543a42f38599ec7dc411b1455c1deb26dfb1c36bc01cf"),
)
X_TICK_ROIS = ((680, 715), (870, 900), (1060, 1100), (1255, 1285),
               (1445, 1485), (1640, 1680), (1830, 1870))
Y_TICK_ROIS = ((2635, 2660), (2505, 2530), (2375, 2400), (2240, 2270),
               (2110, 2140), (1980, 2005), (1850, 1875))
INSPECTED_ANCHORS = np.array([
    [692, 2647, 0, 60], [885.5, 2646, 100, 60], [1078, 2644.5, 200, 60],
    [1269.5, 2643, 300, 60], [1462, 2642, 400, 60], [1655.5, 2641, 500, 60],
    [1849, 2640, 600, 60], [692, 2647, 0, 60], [692, 2517, 0, 70],
    [691, 2384, 0, 80], [691, 2252.5, 0, 90], [690, 2122.5, 0, 100],
    [690, 1992.5, 0, 110], [689, 1860.5, 0, 120],
], dtype=float)
MARKERS = (
    ("C1", "cooling", 742, 2474), ("C2", "cooling", 764, 2460),
    ("C3", "cooling", 791, 2445), ("H1", "heating", 835, 2420),
    ("C4", "cooling", 887, 2390), ("H2", "heating", 986, 2336),
    ("H3", "heating", 1187, 2244), ("H4", "heating", 1384, 2159),
    ("H5", "heating", 1587, 2083), ("H6", "heating", 1665, 2057),
)
EXCLUDED_OVERLAP_ROIS = (
    {"name": "near_200_degC", "xyxy": [1058, 2258, 1118, 2315]},
    {"name": "near_300_degC", "xyxy": [1256, 2168, 1320, 2225]},
    {"name": "near_400_degC", "xyxy": [1455, 2088, 1518, 2145]},
)
READOUT_PX = 4.0
FIELDS = ("point_id", "source_figure", "specimen_id", "material_grade",
    "nominal_foil_thickness_mm", "branch", "pixel_x", "pixel_y",
    "temperature_source_degC", "temperature_K", "temperature_readout_lower_degC",
    "temperature_readout_upper_degC", "resistivity_source_microohm_cm",
    "resistivity_ohm_m", "resistivity_readout_lower_ohm_m",
    "resistivity_readout_upper_ohm_m", "marker_white_hole_pixels",
    "marker_black_depth_pixels", "marker_centre_plateau_pixels")


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def guard_digest(label: str, raw: bytes, expected: str) -> None:
    if sha256(raw) != expected:
        raise ValueError(f"{label} differs from the inspected source")


def load_pages(pdf: Path) -> tuple[dict[int, np.ndarray], list[dict]]:
    raw = pdf.read_bytes()
    guard_digest("PDF", raw, PDF_SHA256)
    try:
        import pymupdf
    except ImportError as error:
        raise RuntimeError("PDF extraction needs the optional PyMuPDF recovery dependency") from error
    pages, evidence = {}, []
    with pymupdf.open(stream=raw, filetype="pdf") as document:
        if len(document) != 15:
            raise ValueError("Expected the inspected 15-page Wilkes report")
        for number, xref, dimensions, content_hash, pixel_hash in PAGE_RECORDS:
            page = document[number - 1]
            guard_digest(f"Page {number} content", page.read_contents(), content_hash)
            images = page.get_images(full=True)
            if len(images) != 1 or images[0][0] != xref or images[0][1] != 0 or images[0][4] != 1:
                raise ValueError("Expected the inspected one-bit monochrome native scan")
            pix = pymupdf.Pixmap(document, xref)
            if (pix.width, pix.height, pix.n) != (*dimensions, 1):
                raise ValueError("Native scan dimensions or channels changed")
            guard_digest(f"Page {number} native pixels", pix.samples, pixel_hash)
            image = np.frombuffer(pix.samples, dtype=np.uint8).reshape(dimensions[1], dimensions[0]).copy()
            if not np.all((image == 0) | (image == 255)):
                raise ValueError("Expected black/white decoded scan pixels")
            pages[number] = image
            evidence.append({"pdf_page_1based": number, "native_image_xref": xref,
                "width_height_px": list(dimensions), "native_bits_per_component": 1,
                "decoded_channels": 1, "page_content_sha256": content_hash,
                "decoded_gray_pixels_sha256": pixel_hash})
    return pages, evidence


def _centre_of_run(indices: np.ndarray) -> float:
    if not 3 <= len(indices) <= 7 or np.any(np.diff(indices) != 1):
        raise ValueError("Inspected axis stroke is not one narrow contiguous run")
    return float(np.mean(indices))


def tick_anchors(image: np.ndarray) -> np.ndarray:
    anchors = []
    for value, (lo, hi) in zip(range(0, 601, 100), X_TICK_ROIS):
        scores = (image[2625:2667, lo:hi] < 128).sum(axis=0)
        x = _centre_of_run(np.flatnonzero(scores > 25) + lo)
        sample_x = np.arange(round(x) - 18, round(x) + 19)
        sample_x = sample_x[abs(sample_x - x) > 7]
        scores = (image[2630:2660, sample_x] < 128).sum(axis=1)
        y = _centre_of_run(np.flatnonzero(scores >= .7 * max(scores)) + 2630)
        anchors.append((x, y, value, 60))
    for value, (lo, hi) in zip(range(60, 121, 10), Y_TICK_ROIS):
        scores = (image[lo:hi, 673:711] < 128).sum(axis=1)
        y = _centre_of_run(np.flatnonzero(scores > 25) + lo)
        sample_y = np.arange(round(y) - 18, round(y) + 19)
        sample_y = sample_y[abs(sample_y - y) > 7]
        scores = (image[sample_y, 680:705] < 128).sum(axis=0)
        x = _centre_of_run(np.flatnonzero(scores >= .7 * max(scores)) + 680)
        anchors.append((x, y, 0, value))
    result = np.asarray(anchors, dtype=float)
    if not np.array_equal(result, INSPECTED_ANCHORS):
        raise ValueError("Tick coordinates differ from the inspected calibration")
    return result


def affine_fit(anchors: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    anchors = np.asarray(anchors, dtype=float)
    if anchors.ndim != 2 or anchors.shape[1] != 4 or len(anchors) < 3 or not np.all(np.isfinite(anchors)):
        raise ValueError("Calibration requires finite [x,y,temperature,resistivity] anchors")
    design = np.column_stack((anchors[:, :2], np.ones(len(anchors))))
    coefficients, _, rank, _ = np.linalg.lstsq(design, anchors[:, 2:], rcond=None)
    if rank != 3:
        raise ValueError("Calibration anchors do not span two native-pixel dimensions")
    return coefficients, anchors[:, 2:] - design @ coefficients


def marker_centre(image: np.ndarray, branch: str, seed_x: int, seed_y: int) -> tuple[float, float, dict]:
    if branch not in ("heating", "cooling"):
        raise ValueError("A branch requires the inspected filled/open diamond legend")
    x0, y0 = seed_x - 16, seed_y - 16
    region = image[y0:seed_y + 17, x0:seed_x + 17]
    if region.shape != (33, 33):
        raise ValueError("Marker ROI lies outside native scan")
    white = region > 128
    labels, count = ndimage.label(white, structure=ndimage.generate_binary_structure(2, 1))
    holes = []
    for label in range(1, count + 1):
        ys, xs = np.where(labels == label)
        if len(xs) >= 4 and not (np.any(xs == 0) or np.any(ys == 0) or np.any(xs == 32) or np.any(ys == 32)):
            holes.append((ys, xs))
    distance = ndimage.distance_transform_cdt(~white, metric="taxicab")
    depth = int(np.max(distance[10:23, 10:23]))
    plateau = 0
    if branch == "cooling":
        if len(holes) != 1 or not 40 <= len(holes[0][0]) <= 65 or depth > 6:
            raise ValueError("Cooling ROI lacks one clearly enclosed white diamond centre")
        ys, xs = holes[0]
        x, y = float(xs.mean() + x0), float(ys.mean() + y0)
        spread_x = spread_y = 0.0
        hole_pixels = len(xs)
    else:
        if holes or not 8 <= depth <= 14:
            raise ValueError("Heating ROI lacks a distinct filled diamond core")
        local_y, local_x = np.where(distance[10:23, 10:23] == depth)
        if not 1 <= len(local_x) <= 9 or np.any((local_x == 0) | (local_x == 12) | (local_y == 0) | (local_y == 12)):
            raise ValueError("Filled-marker core is broad or at the centre-search boundary")
        x, y = float(local_x.mean() + x0 + 10), float(local_y.mean() + y0 + 10)
        spread_x = float(np.ptp(local_x) / 2)
        spread_y = float(np.ptp(local_y) / 2)
        plateau, hole_pixels = len(local_x), 0
    if abs(x - seed_x) > 6 or abs(y - seed_y) > 6:
        raise ValueError("Marker centre exceeds the inspected seed neighbourhood")
    quality = {"roi_xyxy_exclusive": [x0, y0, seed_x + 17, seed_y + 17],
        "roi_gray_pixels_sha256": sha256(region.tobytes()), "white_hole_pixels": hole_pixels,
        "black_taxicab_depth_px": depth, "max_depth_plateau_pixels": plateau,
        "centre_spread_x_px": spread_x, "centre_spread_y_px": spread_y,
        "centre_method": "centroid_of_enclosed_white_face" if branch == "cooling" else "centroid_of_maximum_taxicab_depth_plateau",
        "black_threshold": 128, "white_component_connectivity": 4}
    return x, y, quality


def recover_markers(image: np.ndarray) -> tuple[list[dict], dict]:
    anchors = tick_anchors(image)
    coefficients, residuals = affine_fit(anchors)
    maximum_residual = np.max(abs(residuals), axis=0)
    rows, picks = [], []
    for point_id, branch, seed_x, seed_y in MARKERS:
        x, y, quality = marker_centre(image, branch, seed_x, seed_y)
        if any(a["xyxy"][0] <= x <= a["xyxy"][2] and a["xyxy"][1] <= y <= a["xyxy"][3] for a in EXCLUDED_OVERLAP_ROIS):
            raise ValueError("A selected marker entered an explicitly excluded overlap region")
        temperature, resistivity = np.array([x, y, 1]) @ coefficients
        pixel_bounds = READOUT_PX + np.array([quality["centre_spread_x_px"], quality["centre_spread_y_px"]])
        bound_t, bound_rho = pixel_bounds @ abs(coefficients[:2]) + maximum_residual
        rows.append({"point_id": point_id, "source_figure": "1", "specimen_id": 2,
            "material_grade": "AISI Type 304", "nominal_foil_thickness_mm": .203, "branch": branch,
            "pixel_x": x, "pixel_y": y, "temperature_source_degC": float(temperature),
            "temperature_K": float(temperature + 273.15),
            "temperature_readout_lower_degC": float(temperature - bound_t),
            "temperature_readout_upper_degC": float(temperature + bound_t),
            "resistivity_source_microohm_cm": float(resistivity), "resistivity_ohm_m": float(resistivity * 1e-8),
            "resistivity_readout_lower_ohm_m": float((resistivity - bound_rho) * 1e-8),
            "resistivity_readout_upper_ohm_m": float((resistivity + bound_rho) * 1e-8),
            "marker_white_hole_pixels": quality["white_hole_pixels"],
            "marker_black_depth_pixels": quality["black_taxicab_depth_px"],
            "marker_centre_plateau_pixels": quality["max_depth_plateau_pixels"]})
        picks.append({"point_id": point_id, "branch": branch, "seed_xy": [seed_x, seed_y],
                      "selected_centre_xy": [x, y], **quality})
    if len(rows) != 10 or sum(r["branch"] == "heating" for r in rows) != 6:
        raise ValueError("Expected six heating and four cooling markers")
    details = {"calibration": {"method": "least_squares_2D_affine_fit_to_native_scan_tick_centres",
        "coordinate_origin": "upper_left_native_scan", "anchor_columns": ["pixel_x", "pixel_y", "temperature_degC", "resistivity_microohm_cm"],
        "anchors": anchors.tolist(), "x_axis_label": "Temperature, °C",
        "y_axis_label": "Electrical Resistivity, micro ohm-cm", "coefficients_rows_x_y_constant_columns_temperature_resistivity": coefficients.tolist(),
        "residuals_temperature_degC_resistivity_microohm_cm": residuals.tolist(),
        "maximum_absolute_residuals": maximum_residual.tolist()},
        "marker_picks": picks, "readout": {"estimated_centre_readout_pixels_plus_minus": READOUT_PX,
        "envelope_method": "Propagate rectangular native-pixel centre bounds plus filled-core plateau half-span through absolute affine coefficients, then add maximum axis-fit residual",
        "physical_measurement_uncertainty": False, "confidence_interval": False,
        "marker_symbol_size_is_not_physical_scatter": True}}
    return rows, details


def csv_bytes(rows: list[dict]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({k: format(v, ".12g") if isinstance(v, float) else v for k, v in row.items()})
    return stream.getvalue().encode("utf-8")


def make_qa(image: np.ndarray, rows: list[dict], details: dict, destination: Path) -> dict:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    figure, (left, right) = plt.subplots(1, 2, figsize=(15, 6.5), layout="constrained")
    left.imshow(image, cmap="gray", vmin=0, vmax=255)
    for exclusion in EXCLUDED_OVERLAP_ROIS:
        x0, y0, x1, y1 = exclusion["xyxy"]
        left.add_patch(Rectangle((x0, y0), x1-x0, y1-y0, facecolor="orange", edgecolor="darkorange", alpha=.35))
    for branch, color in [("heating", "red"), ("cooling", "dodgerblue")]:
        selected = [r for r in rows if r["branch"] == branch]
        left.scatter([r["pixel_x"] for r in selected], [r["pixel_y"] for r in selected],
                     marker="o", s=90, facecolors="none", edgecolors=color, linewidths=1.2, label=branch)
        right.errorbar([r["temperature_source_degC"] for r in selected],
            [r["resistivity_source_microohm_cm"] for r in selected],
            xerr=np.array([[r["temperature_source_degC"]-r["temperature_readout_lower_degC"], r["temperature_readout_upper_degC"]-r["temperature_source_degC"]] for r in selected]).T,
            yerr=np.array([[(r["resistivity_ohm_m"]-r["resistivity_readout_lower_ohm_m"])*1e8, (r["resistivity_readout_upper_ohm_m"]-r["resistivity_ohm_m"])*1e8] for r in selected]).T,
            fmt="o", linestyle="none", color=color, capsize=2, ms=4, label=branch)
    for row in rows:
        left.annotate(row["point_id"], (row["pixel_x"]+9, row["pixel_y"]-15), color="blue", fontsize=8)
    left.set(xlim=(430, 1990), ylim=(2865, 1770), xlabel="Native x [pixel]", ylabel="Native y [pixel]",
             title="Figure 1: selected marker centres; orange overlaps excluded")
    left.legend(loc="lower right", fontsize=8)
    right.set(xlabel="Source plotted temperature [°C]", ylabel="Resistivity [micro-ohm cm]",
              title="Specimen 2, 0.203 mm: ten recovered markers")
    right.grid(alpha=.25); right.legend()
    figure.suptitle("Published measurement markers; envelopes describe pixel readout, not physical errors")
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return {"path": str(destination.resolve()), "png_sha256": sha256(destination.read_bytes()),
            "public_source_image_redistributed_in_repository": False}


def metadata_for(rows: list[dict], details: dict, evidence: list[dict], csv_path: Path, raw_csv: bytes, qa: dict | None) -> dict:
    return {"schema": "mhdlab-wilkes-304-marker-recovery-v1", "status": "partially_digitized_published_measurement_markers",
        "source": {"title": "Thermophysical properties of stainless steel foils", "authors": ["Kenneth E. Wilkes", "Joseph P. Strizak", "Fred J. Weaver", "John E. Besser", "David L. Smith"],
            "year": 1997, "report": "ORNL/CP-96813", "doi": "10.2172/658177", "url": "https://www.osti.gov/servlets/purl/658177",
            "pdf_sha256": PDF_SHA256, "pdf_pages": 15, "figure": "1", "pdf_page_1based": 6,
            "printed_page": None, "page_number_note": "Unnumbered manuscript; PDF page is one-based",
            "native_evidence": evidence, "permission_context": "Public official OSTI scan; distribution unlimited and retained government license printed; no new license assigned"},
        "observable": {"quantity": "electrical_resistivity", "not_electrical_conductivity": True,
            "original_units": "micro-ohm cm", "SI_units": "ohm m", "SI_per_original_unit": 1e-8,
            "temperature_original_unit": "degC", "temperature_SI_unit": "K", "K_offset_from_degC": 273.15,
            "value_kind": "Digitized centres of published measurement markers, not raw acquisition samples or the plotted regression line",
            "primary_reported_temperature_domain": "Room temperature to about 500 degC and cooling back to room temperature",
            "retained_temperature_coordinate_bounds_degC": [min(r["temperature_source_degC"] for r in rows), max(r["temperature_source_degC"] for r in rows)],
            "coordinate_bounds_are_continuous_curve_support": False,
            "temperature_points_snapped_to_nominal_setpoints": False},
        "material": {"reported_grade": "AISI Type 304", "specimen_id": 2, "nominal_thickness_mm": .203,
            "specimen_table_evidence": {"table": "I", "pdf_page_1based": 4}, "chemical_analysis": None,
            "heat_treatment": "Authors say rolled foils most likely cold worked; no independently established treatment state",
            "actual_dissertation_foil_match": "Unknown; distinct source specimens",
            "other_304_specimen": {"specimen_id": 1, "nominal_thickness_mm": .305,
                "absolute_resistivity_recovery_status": "not_recovered", "reason": "Figure 2 identifies this specimen but gives deviation percent from CINDAS, not absolute resistivity; no reference reconstruction performed"}},
        "selection": {"inspected_date": "2026-10-18", "retained_markers": 10, "heating_markers": 6, "cooling_markers": 4,
            "legend_mapping": {"filled_diamond": "Heating", "open_diamond": "Cooling"},
            "ROI_half_width_pixels": 16, "excluded_overlap_groups": EXCLUDED_OVERLAP_ROIS,
            "excluded_overlap_hidden_point_count": None,
            "excluded_features": ["Three mixed/merged marker clusters", "Regression line", "CINDAS Type 304/321/201 reference lines", "All Figure 2 deviation markers"],
            "original_marker_count_recovered": False, "smoothing": False, "resampling": False, "interpolation": False, "extrapolation": False},
        "figure_2_status": {"pdf_page_1based": 7, "ordinate": "Deviation, percent", "absolute_resistivity_reconstructed": False,
            "heating_cooling_legend_available": False, "specimen_legend": {"open_diamond": "Specimen 1, Type 304", "open_square": "Specimen 2, Type 304"},
            "reason_not_used": "A percentage-deviation marker is not an absolute-resistivity measurement; no assumed CINDAS baseline or branch assignment"},
        "source_error_context": {"pointwise_physical_error_bars_recovered": False,
            "source_thickness_resolution_mm": .0013, "resolution_evidence_pdf_page": 7,
            "reported_resolution_effect_range_percent": [.4, 2.5], "range_scope": "Thickest and thinnest specimens across the report, not a pointwise uncertainty for specimen 2",
            "comparison_to_CINDAS_percent": [1.5, 2], "comparison_is_measurement_uncertainty": False,
            "between_heat_variability_quantified": False, "confidence_level": None,
            "Flynn_error_statement_transferred": False},
        "recovery": details, "QA_overlay": qa,
        "limitations": ["Only specimen 2 has recovered absolute values; specimen 1 absolute curve unresolved",
            "Higher-temperature cooling and mixed clusters are excluded; no recovered full hysteresis loop",
            "Figure bounds and sparse markers do not authorize a continuous material model",
            "Published marker coordinates and pixel readout are not original acquisition precision or physical error bars",
            "No phase, high-field, nanosecond or actual-foil calibration; no conductivity curve recovered"],
        "output": {"csv_filename": csv_path.name, "csv_sha256": sha256(raw_csv), "rows": len(rows), "columns": list(FIELDS), "adopted_as_runtime_defaults": False}}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--csv", type=Path, default=ROOT / "data/materials/wilkes_1997_304_resistivity.csv")
    parser.add_argument("--metadata", type=Path, default=ROOT / "data/materials/wilkes_1997_304_resistivity.json")
    parser.add_argument("--qa-overlay", type=Path)
    args = parser.parse_args()
    try:
        destinations = [args.csv.resolve(), args.metadata.resolve()] + ([args.qa_overlay.resolve()] if args.qa_overlay else [])
        if len(set(destinations)) != len(destinations) or args.pdf.resolve() in destinations:
            raise ValueError("Source PDF and output destinations must all be different")
        pages, evidence = load_pages(args.pdf)
        rows, details = recover_markers(pages[6])
        qa = make_qa(pages[6], rows, details, args.qa_overlay) if args.qa_overlay else None
        raw_csv = csv_bytes(rows)
        metadata = metadata_for(rows, details, evidence, args.csv, raw_csv, qa)
        args.csv.parent.mkdir(parents=True, exist_ok=True); args.csv.write_bytes(raw_csv)
        args.metadata.parent.mkdir(parents=True, exist_ok=True)
        args.metadata.write_text(json.dumps(metadata, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(2, f"Recovery failed: {error}\n")
    print(json.dumps({"rows": len(rows), "heating": 6, "cooling": 4,
                      "csv_sha256": sha256(raw_csv), "metadata": str(args.metadata)}, sort_keys=True))


if __name__ == "__main__":
    main()
