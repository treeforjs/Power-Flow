#!/usr/bin/env python3
"""Recover partial, literal transition-density lineouts from Figures 4.20/4.21.

This is a guarded raster readout, not an inversion to counts or radiance. The
source PDF stays outside version control. Legend-covered and ambiguous columns
remain missing. Display envelopes describe glyph extent, not measurement noise.
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from PIL import Image


PDF_SHA256 = "1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7"
RGB_SHA256 = "1abbc44c3138785c31690a0732378e8fb7c8e14114cf0d85ea7ca4e3c0337f04"
PAGE_SHA256 = "f0e6e0336f8165ccf7f8f3bd71159cbb5f84e54d9b8251be5867645b13a865fb"
FIGURES = {
    "4.20": {"page_index": 108, "xref": 2707, "page_sha256": PAGE_SHA256,
             "rgb_sha256": RGB_SHA256, "shot": 15750, "thickness_um": 200,
             "x_ticks": [555.5, 632.5, 710.5, 788.5, 866.5],
             "y_ticks": [607.5, 549.5, 491.5, 433.5], "density_ticks": [0, 2e11, 4e11, 6e11]},
    "4.21": {"page_index": 109, "xref": 2721,
             "page_sha256": "4519fa7e1b81b7a5c284d2a9fc922ae6784d5e19d52677fd91a654ec8b829f37",
             "rgb_sha256": "206c076b53ba73e6ccc25fcee714dca5eaf66997f9588ba9e96e4cb7c21b1e92",
             "shot": 15765, "thickness_um": 400,
             "x_ticks": [527.5, 608.5, 689.5, 770.5, 851.5],
             "y_ticks": [605.5, 550.5, 494.5, 438.5], "density_ticks": [0, 1e10, 2e10, 3e10]},
}
SERIES = {"top": (128, 128, 0), "middle": (178, 34, 34), "bottom": (106, 90, 205)}
FIELDS = (
    "series", "support_segment_id", "pixel_x", "pixel_y_median", "matched_exact_color_pixels",
    "wavelength_nm", "wavelength_lower_nm", "wavelength_upper_nm",
    "nominal_transition_density_cm3", "transition_density_lower_cm3",
    "transition_density_upper_cm3",
)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def linear_map(pixels: list[float], values: list[float]) -> tuple[float, float, float]:
    """Fit exact rational native tick centres; round coefficients only once."""
    x, y = [Fraction(v) for v in pixels], [Fraction(v) for v in values]
    n = len(x)
    slope = (n * sum(a * b for a, b in zip(x, y)) - sum(x) * sum(y)) / (
        n * sum(a * a for a in x) - sum(x) ** 2
    )
    intercept = (sum(y) - slope * sum(x)) / n
    residual = max(abs(b - (slope * a + intercept)) for a, b in zip(x, y))
    return float(slope), float(intercept), float(residual)


def load_image(pdf: Path, figure: str = "4.20") -> np.ndarray:
    if digest(pdf.read_bytes()) != PDF_SHA256:
        raise ValueError("Dissertation PDF hash mismatch; fixed pixel calibration is unsafe")
    import pymupdf
    config = FIGURES[figure]
    with pymupdf.open(pdf) as document:
        if len(document) != 203:
            raise ValueError("Expected 203-page dissertation")
        page = document[config["page_index"]]
        if digest(page.read_contents()) != config["page_sha256"] or f"Figure {figure}" not in page.get_text():
            raise ValueError(f"Figure {figure} source-page identity mismatch")
        if config["xref"] not in {entry[0] for entry in page.get_images(full=True)}:
            raise ValueError(f"Expected embedded Figure {figure} image is absent")
        image = Image.open(io.BytesIO(document.extract_image(config["xref"])["image"]))
        if image.mode != "RGB" or image.size != (1040, 692) or digest(image.tobytes()) != config["rgb_sha256"]:
            raise ValueError(f"Figure {figure} RGB identity mismatch")
        return np.asarray(image)


def recover(image: np.ndarray, figure: str = "4.20") -> tuple[list[dict], dict]:
    if image.shape != (692, 1040, 3):
        raise ValueError("Use the native RGB image without scaling")
    config = FIGURES[figure]
    x_ticks = config["x_ticks"]
    wavelengths = [120, 130, 140, 150, 160]
    y_ticks = config["y_ticks"]
    density_ticks = config["density_ticks"]
    xs, xi, xr = linear_map(x_ticks, wavelengths)
    ys, yi, yr = linear_map(y_ticks, density_ticks)
    rows, omitted = [], []
    for series, rgb in SERIES.items():
        segment_id, previous_x = -1, None
        for x in range(514, 905, 4):
            # The complete legend rectangle is excluded, including faint curves
            # visible through it. Never undo its alpha compositing.
            matches = np.flatnonzero(np.all(image[414:614, x] == rgb, axis=1)) + 414
            matches = matches[~((517 <= x) & (x <= 795) & (matches >= 420) & (matches <= 525))]
            reason = None
            if not len(matches):
                reason = "no_unobscured_exact_color"
            elif len(matches) > 8:
                reason = "overly_broad_stroke"
            elif np.any(np.diff(matches) != 1):
                reason = "separated_runs"
            elif matches[0] == 414 or matches[-1] == 613:
                reason = "selection_boundary_contact"
            elif 517 <= x <= 795 and (matches[0] == 526 or matches[-1] == 419):
                reason = "legend_boundary_contact"
            if reason:
                omitted.append({"series": series, "pixel_x": x, "reason": reason})
                continue
            y = float(np.median(matches))
            if previous_x is None or x != previous_x + 4:
                segment_id += 1
            previous_x = x
            rows.append({
                "series": series, "support_segment_id": segment_id, "pixel_x": x, "pixel_y_median": y,
                "matched_exact_color_pixels": len(matches),
                "wavelength_nm": xs * x + xi,
                "wavelength_lower_nm": xs * (x - 1) + xi - xr,
                "wavelength_upper_nm": xs * (x + 1) + xi + xr,
                "nominal_transition_density_cm3": ys * y + yi,
                "transition_density_lower_cm3": ys * (matches[-1] + 1) + yi - yr,
                "transition_density_upper_cm3": ys * (matches[0] - 1) + yi + yr,
            })
    if not rows:
        raise ValueError("No unambiguous lineout pixels recovered")
    return rows, {
        "axis_calibration": {
            "method": "exact_rational_least_squares_of_inspected_tick_centres",
            "wavelength": {"tick_pixels": x_ticks, "tick_values_nm": wavelengths,
                           "slope_nm_per_pixel": xs, "intercept_nm": xi, "max_residual_nm": xr},
            "density": {"tick_pixels": y_ticks, "tick_values_emitters_cm3": density_ticks,
                        "slope_per_pixel": ys, "intercept": yi, "max_residual": yr},
        },
        "series_rgb": {key: list(value) for key, value in SERIES.items()},
        "selection": {"pixel_x_first": 514, "pixel_x_last": 902, "stride_pixels": 4,
                      "pixel_y_inclusive": [414, 613], "legend_exclusion_xyxy": [517, 420, 795, 525],
                      "accepted_contiguous_run_pixels": [1, 8],
                      "antialiased_colors_used": False, "interpolated_missing_columns": False},
        "omitted_columns": omitted,
        "readout_envelopes": "full exact-color stroke plus one native pixel and maximum tick-fit residual; correlated glyph extent, not confidence or noise",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--figure", choices=tuple(FIGURES), default="4.20")
    parser.add_argument("--output-dir", type=Path, default=Path("data/dissertation"))
    parser.add_argument("--qa-plot", type=Path, help="Local source-pixel overlay; retain outside public artifacts")
    args = parser.parse_args()
    config = FIGURES[args.figure]
    image = load_image(args.pdf, args.figure)
    rows, recovery = recover(image, args.figure)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = "figure_" + args.figure.replace(".", "_") + "_wavelength_lineouts"
    csv_path = args.output_dir / (stem + ".csv")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "schema": "mhdlab-published-vuv-transition-density-lineouts-v1",
        "source": {"pdf_sha256": PDF_SHA256, "page_content_sha256": config["page_sha256"],
                   "native_rgb_sha256": config["rgb_sha256"], "figure": args.figure, "printed_page": config["page_index"] - 24,
                   "pdf_page_1based": config["page_index"] + 1, "native_image_xref": config["xref"], "native_size_px": [1040, 692]},
        "experiment": {"shot": config["shot"], "foil_material": "SS304", "foil_thickness_um": config["thickness_um"],
                       "diagnostic": "VUV gated spectroscopy", "gate_duration_ns": 3,
                       "signed_gate_start_s": None, "signed_gate_end_s": None,
                       "time_registration": "before-breakdown context only; gate endpoints not inferred from Table 4.1"},
        "observable": {"quantity": "plotted_transition_density", "axis_label_literal": "Transition Density [Emitters/cm^3]",
                       "units": "emitters/cm^3", "coordinate": "wavelength_nm_as_displayed",
                       "series_literal": {"top": "Top Lineout", "middle": "Middle Lineout", "bottom": "Bottom Lineout"},
                       "spatial_support": "top/middle/bottom dashed bands in figure; quantitative collection kernel unresolved",
                       "raw_counts": False, "absolute_radiance": False, "ground_neutral_density": False,
                       "independent_measurements": False, "noise_covariance_recovered": False,
                       "integrals_supported": False, "missing_segments_retained": True,
                       "yacora_curve_recovered": False},
        "recovery": recovery,
        "row_count": len(rows),
        "row_count_by_series": {series: sum(row["series"] == series for row in rows) for series in SERIES},
        "csv_sha256": digest(csv_path.read_bytes()),
        "limitations": [
            "The ordinate is a derived published quantity. No unique count/radiance inverse mapping is established.",
            "Legend coverage, curve overlap and steep/wide strokes omit columns. No interpolation or complete line integral is supplied.",
            "Sampling stride is neither a detector-bin width nor the instrumental line-spread function.",
            "Upper-state assignment and per-bin transition reduction remain unresolved; a scalar calibration factor cannot repair them.",
        ],
    }
    if args.qa_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(12, 7))
        ax.imshow(image)
        for series in SERIES:
            selected = [row for row in rows if row["series"] == series]
            ax.scatter([r["pixel_x"] for r in selected], [r["pixel_y_median"] for r in selected],
                       s=12, facecolors="none", edgecolors="black", linewidths=0.6)
        ax.set(xlim=(490, 925), ylim=(640, 390), title=f"Figure {args.figure}: retained literal display readouts")
        args.qa_plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.qa_plot, dpi=180)
        plt.close(fig)
        metadata["local_source_overlay_sha256"] = digest(args.qa_plot.read_bytes())
    sidecar = args.output_dir / (stem + ".json")
    sidecar.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"csv": str(csv_path), "sidecar": str(sidecar), "rows": len(rows),
                      "by_series": metadata["row_count_by_series"]}))


if __name__ == "__main__":
    main()
