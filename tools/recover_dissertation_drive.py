#!/usr/bin/env python3
"""Recover shot 15210's published current polyline from dissertation Fig. C.11.

Figure A.10 (printed p111 / PDF p136) was inspected first: its current ramp
is clipped above approximately 0.24 MA. The same-shot Fig. C.11 (printed p157 /
PDF p182) contains a full vector curve. This tool retains its encoded line
vertices, including vertices behind the later legend, and labels that overlap.
These are published display commands, not original oscilloscope samples.

Inspected calibration, 2026-10-18: native PDF page coordinates in points, origin
at upper left; eleven time ticks 0:20:200 ns and six current ticks 0:.2:1 MA.
The exact red stroke is drawing 87, sequence 174, 119 connected line segments.
Only original vertices inside the displayed axes are retained. No endpoint is
inserted, no extrapolation is performed, and no pre-current zero is assumed.
The source time origin/trigger is unknown and is not shifted. Display-envelope
bounds use half the published stroke width plus maximum axis-fit residual;
they are neither measurement uncertainty nor uncertainty in the physical drive.

Dependencies: numpy, pymupdf, pillow. Optional --qa-plot needs matplotlib.
Example:
  python tools/recover_dissertation_drive.py --pdf /path/to/smtrevor_1.pdf \
      --qa-plot runs/recovery/shot_15210_current_qa.png
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


ROOT = Path(__file__).resolve().parents[1]
PDF_SHA256 = "1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7"
PAGE_INDEX = 181
PAGE_CONTENT_SHA256 = "1d12ec2ed170923d2a27a3082a0317a50ce048d8f25a7601b21bd65d6c36e50f"
PATH_SHA256 = "26428b0374930703af2f9d64f5f175de32605125b6d3b955ddf8bdac321bd549"
A10_IMAGE_XREF = 3001
A10_PNG_SHA256 = "57020d7de286661007e01d109334877d9c1d6344544db7933a699e7e3e9ebac9"
A10_RGB_SHA256 = "ca96300fd9a56b4f676f2b18e9d25a75b4eef5b0395861538b1967ff2492cf5f"
DRAWING_INDEX = 87
X_TICK_DRAWINGS = [5, 7, 10, 13, 16, 19, 23, 27, 31, 35, 39]
Y_TICK_DRAWINGS = [51, 55, 59, 63, 67, 71]
X_TICKS_PT = np.array([110.60025024414062, 150.42379760742188, 190.24766540527344,
    230.0712127685547, 269.89508056640625, 309.7186279296875, 349.5425109863281,
    389.3660888671875, 429.19122314453125, 469.01348876953125, 508.8357238769531])
Y_TICKS_PT = np.array([517.0767211914062, 452.0224609375, 386.96820068359375,
    321.9139099121094, 256.86029052734375, 191.8050537109375])
TIME_TICKS_NS = np.arange(0., 201., 20.)
CURRENT_TICKS_MA = np.arange(6) * .2
FIELDS = ("source_vertex_index", "page_x_pt", "page_y_pt", "time_ns", "time_s",
    "time_digitization_lower_ns", "time_digitization_upper_ns", "current_MA", "current_A",
    "current_digitization_lower_A", "current_digitization_upper_A", "occluded_by_legend")


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(record: dict) -> bytes:
    return json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def recover(pdf: Path) -> tuple[list[dict], dict]:
    if sha256(pdf.read_bytes()) != PDF_SHA256:
        raise ValueError("PDF differs from the inspected dissertation")
    with pymupdf.open(pdf) as document:
        if len(document) != 203:
            raise ValueError("Expected 203 PDF pages")
        page = document[PAGE_INDEX]
        text = page.get_text()
        if "Figure C.11" not in text or "Average Current [MA]" not in text or "Time (ns)" not in text:
            raise ValueError("Expected shot 15210 current figure and inspected axis units")
        if page.get_images(full=True):
            raise ValueError("Expected vector-only C.11, without an embedded native raster")
        if sha256(page.read_contents()) != PAGE_CONTENT_SHA256:
            raise ValueError("Figure C.11 content stream changed")
        drawings = page.get_drawings()
        path = drawings[DRAWING_INDEX]
        if path["type"] != "s" or path["seqno"] != 174 or len(path["items"]) != 119:
            raise ValueError("Expected the inspected current stroke")
        points = [tuple(path["items"][0][1])]
        for item in path["items"]:
            if item[0] != "l" or tuple(item[1]) != points[-1]:
                raise ValueError("Current must be one continuous straight-segment polyline")
            points.append(tuple(item[2]))
        geometry = {"color": list(path["color"]), "width": path["width"], "points": points}
        if sha256(canonical(geometry)) != PATH_SHA256:
            raise ValueError("Current vector geometry changed; inspect before recalibrating")
        x_ticks = np.array([drawings[i]["items"][0][1].x for i in X_TICK_DRAWINGS])
        y_ticks = np.array([drawings[i]["items"][0][1].y for i in Y_TICK_DRAWINGS])
        if not np.array_equal(x_ticks, X_TICKS_PT) or not np.array_equal(y_ticks, Y_TICKS_PT):
            raise ValueError("Inspected tick-mark coordinates changed")
        axes = tuple(drawings[1]["rect"])
        legend = tuple(drawings[158]["rect"])
        gate = tuple(drawings[2]["rect"])
        # Validate the inspected, clipped A.10 raster only as a paired reference.
        if A10_IMAGE_XREF not in {r[0] for r in document[135].get_images(full=True)}:
            raise ValueError("Paired Fig. A.10 image is absent")
        a10_encoded = document.extract_image(A10_IMAGE_XREF)["image"]
        a10_image = Image.open(io.BytesIO(a10_encoded))
        if a10_image.mode != "RGB" or a10_image.size != (1024, 512) or sha256(a10_image.tobytes()) != A10_RGB_SHA256:
            raise ValueError("Paired native A.10 pixels changed")

    time_fit = np.polyfit(X_TICKS_PT, TIME_TICKS_NS, 1)
    current_fit = np.polyfit(Y_TICKS_PT, CURRENT_TICKS_MA, 1)
    time_residual = TIME_TICKS_NS - np.polyval(time_fit, X_TICKS_PT)
    current_residual = CURRENT_TICKS_MA - np.polyval(current_fit, Y_TICKS_PT)
    half_width = path["width"] / 2.
    time_bound = abs(float(time_fit[0])) * half_width + float(max(abs(time_residual)))
    current_bound = abs(float(current_fit[0])) * half_width + float(max(abs(current_residual)))
    rows, excluded = [], []
    for index, (x, y) in enumerate(points):
        if not (axes[0] <= x <= axes[2] and axes[1] <= y <= axes[3]):
            excluded.append({"source_vertex_index": index, "page_point": [x, y], "reason": "outside_displayed_axes"})
            continue
        time = float(np.polyval(time_fit, x))
        current = float(np.polyval(current_fit, y))
        rows.append({"source_vertex_index": index, "page_x_pt": x, "page_y_pt": y,
            "time_ns": time, "time_s": time * 1.e-9,
            "time_digitization_lower_ns": time - time_bound,
            "time_digitization_upper_ns": time + time_bound,
            "current_MA": current, "current_A": current * 1.e6,
            "current_digitization_lower_A": (current - current_bound) * 1.e6,
            "current_digitization_upper_A": (current + current_bound) * 1.e6,
            "occluded_by_legend": int(legend[0] <= x <= legend[2] and legend[1] <= y <= legend[3])})
    if len(rows) != 118 or any(b["time_ns"] <= a["time_ns"] for a, b in zip(rows, rows[1:])):
        raise ValueError("Expected 118 original vertices with strictly increasing time")
    details = {
        "vector_path": {"drawing_index_zero_based": DRAWING_INDEX, "drawing_sequence": path["seqno"],
            "RGB": list(path["color"]), "stroke_width_pt": path["width"], "segments": 119,
            "encoded_vertices": len(points), "retained_vertices": len(rows),
            "canonical_geometry_sha256": PATH_SHA256,
            "canonical_geometry_hash_format": "sorted compact JSON with color, width, points keys; native float coordinates",
            "excluded_vertices": excluded, "axes_rectangle_pt": list(axes),
            "later_legend_bounding_rectangle_pt": list(legend),
            "occlusion_test": "vertex inside later legend bounding rectangle; conservative for rounded corners",
            "vertices_behind_legend": sum(r["occluded_by_legend"] for r in rows)},
        "calibration": {"inspected_date": "2026-10-18", "coordinates": "native_PDF_points_origin_upper_left",
            "method": "least_squares_linear_fit_to_fixed_inspected_tick_mark_centres",
            "time_axis": {"label": "Time (ns)", "tick_page_x_pt": X_TICKS_PT.tolist(),
                "tick_values_ns": TIME_TICKS_NS.tolist(), "ns_per_pt": float(time_fit[0]),
                "intercept_ns": float(time_fit[1]), "fit_residuals_ns": time_residual.tolist()},
            "current_axis": {"label": "Signal Height", "legend": "Average Current [MA]",
                "tick_page_y_pt": Y_TICKS_PT.tolist(), "tick_values_MA": CURRENT_TICKS_MA.tolist(),
                "MA_per_pt": float(current_fit[0]), "intercept_MA": float(current_fit[1]),
                "fit_residuals_MA": current_residual.tolist(), "A_per_MA": 1.e6},
            "digitization_bounds": {"interpretation": "conservative displayed stroke envelope plus maximum tick-fit residual; not physical uncertainty or confidence intervals",
                "coordinate_half_stroke_pt": half_width, "time_plus_minus_ns": time_bound,
                "current_plus_minus_A": current_bound * 1.e6,
                "native_pixel_readout_uncertainty_applicable": False}},
        "VUV_gate": {"kind": "published_gate_rectangle_not_breakdown", "drawing_index_zero_based": 2,
            "page_rectangle_pt": list(gate), "plotted_time_bounds_ns": np.polyval(time_fit, [gate[0], gate[2]]).tolist(),
            "display_coordinate_bound_ns": time_bound,
            "table_4_1_published_acquisition_time_ns": 55, "table_4_1_relative_to_breakdown": "After",
            "actual_breakdown_time_ns": None, "exact_trigger_reference": None,
            "gate_edge_meaning": "geometrical extent of published rectangle; electronics edge definition not recovered"},
        "paired_A10": {"figure": "A.10", "printed_page": 111, "pdf_page": 136,
            "embedded_image_xref": A10_IMAGE_XREF, "image_size_pixels": [1024, 512],
            "reference_encoded_PNG_sha256": A10_PNG_SHA256,
            "extracted_encoded_PNG_sha256": sha256(a10_encoded), "decoded_RGB_pixels_sha256": A10_RGB_SHA256,
            "role": "inspected same-shot reference; current ramp clips at upper axis near 0.24 MA, so C.11 supplies the recovered current"},
    }
    return rows, details


def make_qa(pdf: Path, rows: list[dict], details: dict, destination: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    clip = pymupdf.Rect(80, 175, 525, 568)
    with pymupdf.open(pdf) as document:
        pix = document[PAGE_INDEX].get_pixmap(matrix=pymupdf.Matrix(2, 2), clip=clip, alpha=False)
        pixels = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    figure, (left, right) = plt.subplots(1, 2, figsize=(13, 5.5), layout="constrained")
    left.imshow(pixels, extent=(clip.x0, clip.x1, clip.y1, clip.y0))
    for occluded, marker in ((0, "."), (1, "s")):
        selected = [r for r in rows if r["occluded_by_legend"] == occluded]
        left.scatter([r["page_x_pt"] for r in selected], [r["page_y_pt"] for r in selected],
            marker=marker, s=12 if occluded else 8, facecolors="none" if occluded else "black",
            edgecolors="black", linewidths=.6, label="Encoded behind legend" if occluded else "Recovered vertices")
    left.set(xlabel="Native PDF x [pt]", ylabel="Native PDF y [pt]", title="C.11 source and encoded current vertices")
    left.legend(loc="lower right", fontsize=8)
    times = np.array([r["time_ns"] for r in rows])
    current = np.array([r["current_MA"] for r in rows])
    bound = details["calibration"]["digitization_bounds"]["current_plus_minus_A"] * 1.e-6
    right.fill_between(times, current - bound, current + bound, color="firebrick", alpha=.2,
        label="Published-stroke envelope")
    right.plot(times, current, ".-", ms=3, lw=.8, color="firebrick", label="Original polyline vertices")
    gate = details["VUV_gate"]["plotted_time_bounds_ns"]
    right.axvspan(*gate, color="blueviolet", alpha=.3, label="VUV gate, not breakdown")
    hidden = [r for r in rows if r["occluded_by_legend"]]
    right.scatter([r["time_ns"] for r in hidden], [r["current_MA"] for r in hidden],
        marker="s", s=18, facecolors="none", edgecolors="black", label="Encoded behind legend")
    right.set(xlabel="Source plotted time [ns]; trigger origin unknown", ylabel="Average current [MA]",
        title="Shot 15210: published display curve")
    right.grid(alpha=.2)
    right.legend(fontsize=8)
    figure.suptitle("Reconstructed published vector polyline; original waveform and physical error bars unavailable", fontsize=11)
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=150)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--csv", type=Path, default=ROOT / "data/dissertation/shot_15210_current.csv")
    parser.add_argument("--metadata", type=Path, default=ROOT / "data/dissertation/shot_15210_current.json")
    parser.add_argument("--qa-plot", type=Path)
    args = parser.parse_args()
    rows, details = recover(args.pdf)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format(value, ".12g") if isinstance(value, float) else value for key, value in row.items()})
    metadata = {"schema_version": 1, "status": "reconstructed_published_vector_polyline",
        "source": {"author": "Trevor Johannes Smith", "year": 2025,
            "title": "Characterizing Desorbed Surface Contaminants Leading to Plasma Formation in Power Feeds of Pulsed Power Accelerators",
            "item_url": "https://deepblue.lib.umich.edu/items/9ccfa886-0b19-4bd8-b3b2-bb685b9ddefc",
            "pdf_filename": "smtrevor_1.pdf", "pdf_sha256": PDF_SHA256,
            "figure": "C.11", "printed_page": 157, "pdf_page": 182, "pdf_page_index_zero_based": PAGE_INDEX,
            "page_content_stream_xrefs": [4366], "decoded_page_content_stream_sha256": PAGE_CONTENT_SHA256,
            "native_image_pixel_hash": None, "native_image_pixel_hash_applicability": "not applicable: C.11 contains vector paths and no embedded raster",
            "parser_reference_version": "PyMuPDF 1.26.6", "identity_guards": "PDF, content stream, vector geometry and tick coordinates are checked"},
        "experiment": {"facility": "Mykonos", "shot": 15210, "foil_material": "SS304", "foil_thickness_um": 200,
            "thickness_evidence": ["Table 4.1", "Figure 4.23 caption", "Figure A.10 caption", "Figure C.11 in-plot title"],
            "unresolved_caption_conflicts": ["Figure B.9 caption labels shot 15210 as 500 micrometres", "Figure C.11 caption labels shot 15210 as 500 micrometres, while its own plot title says 200 micrometres"]},
        "observable": {"quantity": "published_average_current", "original_units": "MA", "output_units": "A",
            "time_coordinate": "unshifted source plotted time in ns", "time_origin_reference": None,
            "trigger_or_current_start_offset_recovered": False, "raw_oscilloscope_data": False,
            "independent_samples": False, "source_plot_simplification_or_filtering": "unknown",
            "negative_early_current_vertices_preserved": True},
        **details,
        "output": {"csv_filename": args.csv.name, "csv_sha256": sha256(args.csv.read_bytes()),
            "columns": list(FIELDS), "samples": len(rows),
            "time_extent_ns": [rows[0]["time_ns"], rows[-1]["time_ns"]],
            "current_extent_A": [min(r["current_A"] for r in rows), max(r["current_A"] for r in rows)],
            "generator": "tools/recover_dissertation_drive.py", "generator_sha256": sha256(Path(__file__).read_bytes())},
        "limits": [
            "Encoded display vertices are recovered from the published PDF, not original acquisition samples or a fitted drive waveform.",
            "Vertices behind the later white legend are explicitly flagged; they are encoded source geometry and are not inferred from neighbouring pixels.",
            "The two encoded endpoints outside the displayed axes are excluded; no points at time zero or 200 ns are inserted.",
            "No pre-pulse zero, physical time origin, tail, symmetry, resampling, smoothing or extrapolation is imposed.",
            "Small negative early display values are retained without flooring; a physical drive consumer must choose and document its treatment.",
            "Display-envelope bounds quantify plotted line width and axis-fit residual only, not diagnostic calibration, shot variability or model uncertainty.",
            "The VUV gate rectangle is independent metadata; no breakdown event time is inferred from it or from the voltage curve.",
            "Decimal precision preserves coordinate conversion and does not establish physical measurement accuracy.",
        ]}
    args.metadata.parent.mkdir(parents=True, exist_ok=True)
    args.metadata.write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if args.qa_plot:
        make_qa(args.pdf, rows, details, args.qa_plot)
    print(f"Recovered {len(rows)} published current vertices: {args.csv}")


if __name__ == "__main__":
    main()
