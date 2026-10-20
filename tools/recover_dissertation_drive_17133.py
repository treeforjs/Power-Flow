#!/usr/bin/env python3
"""Recover the published shot 17133 current vertices from dissertation Figure C.19.

The source is a vector display curve, rather than original acquisition samples.
Calibration was inspected on 2026-10-19 from eleven time ticks 0:20:200 ns and six
shared ordinate ticks 0:.2:1. The current legend supplies MA; the shared ordinate
label remains literally Signal Height. The plotted time origin is unknown.

Only original encoded vertices inside the displayed axes are retained. Their
signs and ordering remain unchanged; no boundary intersections, zero endpoints,
smoothing, resampling, extrapolation or trigger alignment are supplied. A later
white legend is checked for occlusion even if the resulting count is zero.
Bounds describe the displayed stroke and tick fit, not physical uncertainty.

Dependencies: numpy and pymupdf; optional QA output also needs matplotlib.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pymupdf


ROOT = Path(__file__).resolve().parents[1]
PDF_SHA256 = "1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7"
PAGE_INDEX = 189
PAGE_CONTENT_XREFS = [4991]
PAGE_CONTENT_SHA256 = "0a6619c5541d64d19d1534e738268c6d7868eca658152967a2f8c8dfe9c57edb"
PATH_SHA256 = "cebaa12de7cf1187cf0a8f2766f773ba165408d4d548e96ab8ebe3709d85117a"
DRAWING_INDEX = 20
DRAWING_SEQUENCE = 55
AXES_DRAWING = 1
GATE_DRAWING = 2
LEGEND_OCCLUDER_DRAWING = 27
LEGEND_SAMPLE_DRAWING = 28
X_TICK_DRAWINGS = tuple(range(3, 14))
Y_TICK_DRAWINGS = tuple(range(14, 20))
X_TICKS_PT = np.array([110.60025024414062, 150.42379760742188, 190.24766540527344,
    230.0712127685547, 269.89508056640625, 309.7186279296875, 349.5425109863281,
    389.3660888671875, 429.19122314453125, 469.01348876953125, 508.8357238769531])
Y_TICKS_PT = np.array([477.6004333496094, 430.1172790527344, 382.6341552734375,
    335.15130615234375, 287.66815185546875, 240.18499755859375])
TIME_TICKS_NS = np.array([0., 20., 40., 60., 80., 100., 120., 140., 160., 180., 200.])
CURRENT_TICKS_MA = np.array([0., .2, .4, .6, .8, 1.])
AXES_PT = (110.60025024414062, 240.18499755859375, 508.8389892578125, 501.3420104980469)
GATE_PT = (162.44944763183594, 240.18499755859375, 168.42311096191406, 501.3420104980469)
LEGEND_PT = (411.8297424316406, 243.43499755859375, 505.5889892578125, 286.68536376953125)
FIELDS = ("source_vertex_index", "page_x_pt", "page_y_pt", "time_ns", "time_s",
    "time_digitization_lower_ns", "time_digitization_upper_ns", "current_MA", "current_A",
    "current_digitization_lower_A", "current_digitization_upper_A", "occluded_by_legend")


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(record: dict) -> bytes:
    return json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def recover(pdf: Path) -> tuple[list[dict], dict]:
    """Validate the inspected source and return only its retained native vertices."""
    if sha256(pdf.read_bytes()) != PDF_SHA256:
        raise ValueError("PDF differs from the inspected dissertation")
    with pymupdf.open(pdf) as document:
        if len(document) != 203:
            raise ValueError("Expected 203 dissertation pages")
        page = document[PAGE_INDEX]
        text = page.get_text()
        required_text = ("Figure C.19", "Shot 17133", "Average Current [MA", "Time (ns)", "Signal Height",
                         "SS304, 200-\u03bcm Foil", "200-\u00b5m\u2013SS304")
        if any(value not in text for value in required_text):
            raise ValueError("Expected the inspected shot 17133 figure, material, thickness and units")
        if page.get_images(full=True):
            raise ValueError("Expected vector-only Figure C.19 without a native raster")
        if page.get_contents() != PAGE_CONTENT_XREFS or sha256(page.read_contents()) != PAGE_CONTENT_SHA256:
            raise ValueError("Figure C.19 content stream changed")
        drawings = page.get_drawings()
        if len(drawings) != 33:
            raise ValueError("Expected 33 inspected source drawings")
        path = drawings[DRAWING_INDEX]
        if path["type"] != "s" or path["seqno"] != DRAWING_SEQUENCE or len(path["items"]) != 133:
            raise ValueError("Expected the inspected 133-segment red current stroke")
        points = [tuple(path["items"][0][1])]
        for item in path["items"]:
            if item[0] != "l" or tuple(item[1]) != points[-1]:
                raise ValueError("Current must be one continuous straight-segment polyline")
            points.append(tuple(item[2]))
        geometry = {"color": list(path["color"]), "width": path["width"], "points": points}
        if sha256(canonical(geometry)) != PATH_SHA256:
            raise ValueError("Current vector geometry changed; inspect before recalibrating")
        x_ticks = np.array([drawings[index]["items"][0][1].x for index in X_TICK_DRAWINGS])
        y_ticks = np.array([drawings[index]["items"][0][1].y for index in Y_TICK_DRAWINGS])
        if not np.array_equal(x_ticks, X_TICKS_PT) or not np.array_equal(y_ticks, Y_TICKS_PT):
            raise ValueError("Inspected tick coordinates changed")
        axes = tuple(drawings[AXES_DRAWING]["rect"])
        gate = tuple(drawings[GATE_DRAWING]["rect"])
        legend = tuple(drawings[LEGEND_OCCLUDER_DRAWING]["rect"])
        if axes != AXES_PT or gate != GATE_PT or legend != LEGEND_PT:
            raise ValueError("Inspected axes, gate or later legend geometry changed")
        occluder = drawings[LEGEND_OCCLUDER_DRAWING]
        sample = drawings[LEGEND_SAMPLE_DRAWING]
        if (occluder["seqno"] != 63 or occluder["fill"] != (1., 1., 1.)
                or occluder["seqno"] <= path["seqno"] or sample["seqno"] != 65
                or sample["color"] != path["color"] or len(sample["items"]) != 1):
            raise ValueError("The later white occluder and red legend sample must remain distinct")
        literal_legend = next(line for line in text.splitlines() if line.startswith("Average Current [MA"))

    time_fit = np.polyfit(X_TICKS_PT, TIME_TICKS_NS, 1)
    current_fit = np.polyfit(Y_TICKS_PT, CURRENT_TICKS_MA, 1)
    time_residual = TIME_TICKS_NS-np.polyval(time_fit, X_TICKS_PT)
    current_residual = CURRENT_TICKS_MA-np.polyval(current_fit, Y_TICKS_PT)
    half_width = float(path["width"])/2.
    time_bound = abs(float(time_fit[0]))*half_width+float(max(abs(time_residual)))
    current_bound = abs(float(current_fit[0]))*half_width+float(max(abs(current_residual)))
    rows, excluded = [], []
    for index, (x, y) in enumerate(points):
        if not (axes[0] <= x <= axes[2] and axes[1] <= y <= axes[3]):
            excluded.append({"source_vertex_index": index, "page_point": [x, y], "reason": "outside_displayed_axes"})
            continue
        time = float(np.polyval(time_fit, x))
        current = float(np.polyval(current_fit, y))
        rows.append({"source_vertex_index": index, "page_x_pt": x, "page_y_pt": y,
            "time_ns": time, "time_s": time*1e-9,
            "time_digitization_lower_ns": time-time_bound, "time_digitization_upper_ns": time+time_bound,
            "current_MA": current, "current_A": current*1e6,
            "current_digitization_lower_A": (current-current_bound)*1e6,
            "current_digitization_upper_A": (current+current_bound)*1e6,
            "occluded_by_legend": int(legend[0] <= x <= legend[2] and legend[1] <= y <= legend[3])})
    if (len(rows) != 132 or [row["source_vertex_index"] for row in rows] != list(range(1, 133))
            or any(b["time_ns"] <= a["time_ns"] for a, b in zip(rows, rows[1:]))):
        raise ValueError("Expected 132 original vertices with strictly increasing time")
    details = {
        "vector_path": {"drawing_index_zero_based": DRAWING_INDEX, "drawing_sequence": path["seqno"],
            "RGB": list(path["color"]), "stroke_width_pt": path["width"], "segments": 133,
            "encoded_vertices": len(points), "retained_vertices": len(rows),
            "canonical_geometry_sha256": PATH_SHA256,
            "canonical_geometry_hash_format": "sorted compact JSON with color, width, points keys; native float coordinates",
            "excluded_vertices": excluded, "axes_rectangle_pt": list(axes),
            "later_legend_drawing_index_zero_based": LEGEND_OCCLUDER_DRAWING,
            "later_legend_drawing_sequence": occluder["seqno"],
            "later_legend_bounding_rectangle_pt": list(legend),
            "legend_current_sample_drawing_index_zero_based": LEGEND_SAMPLE_DRAWING,
            "legend_current_sample_drawing_sequence": sample["seqno"],
            "occlusion_test": "vertex inside later white legend bounding rectangle; conservative for rounded corners",
            "vertices_behind_legend": sum(row["occluded_by_legend"] for row in rows)},
        "calibration": {"inspected_date": "2026-10-19", "coordinates": "native_PDF_points_origin_upper_left",
            "method": "least_squares_linear_fit_to_fixed_inspected_tick_mark_centres",
            "time_axis": {"label": "Time (ns)", "tick_drawing_indices_zero_based": list(X_TICK_DRAWINGS),
                "tick_page_x_pt": X_TICKS_PT.tolist(), "tick_values_ns": TIME_TICKS_NS.tolist(),
                "ns_per_pt": float(time_fit[0]), "intercept_ns": float(time_fit[1]),
                "fit_residuals_ns": time_residual.tolist(), "seconds_per_ns": 1e-9},
            "current_axis": {"label": "Signal Height", "legend": "Average Current [MA]",
                "literal_extracted_legend_text": literal_legend,
                "legend_identifier_interpretation": "MA is directly published; the normalized reader identifier restores the closing bracket omitted by PDF text extraction",
                "tick_drawing_indices_zero_based": list(Y_TICK_DRAWINGS),
                "tick_page_y_pt": Y_TICKS_PT.tolist(), "tick_values_MA": CURRENT_TICKS_MA.tolist(),
                "MA_per_pt": float(current_fit[0]), "intercept_MA": float(current_fit[1]),
                "fit_residuals_MA": current_residual.tolist(), "A_per_MA": 1e6,
                "shared_ordinate_interpretation": "Signal Height is the literal shared label; only the red Average Current legend supplies MA",
                "zero_current_tick_page_y_pt": float(Y_TICKS_PT[0]),
                "frame_bottom_is_zero_current": False},
            "digitization_bounds": {"interpretation": "conservative displayed stroke envelope plus maximum tick-fit residual; not physical uncertainty or confidence intervals",
                "coordinate_half_stroke_pt": half_width, "time_plus_minus_ns": time_bound,
                "current_plus_minus_A": current_bound*1e6, "native_pixel_readout_uncertainty_applicable": False}},
        "VUV_gate": {"kind": "published_gate_rectangle_not_breakdown", "drawing_index_zero_based": GATE_DRAWING,
            "page_rectangle_pt": list(gate), "plotted_time_bounds_ns": np.polyval(time_fit, [gate[0], gate[2]]).tolist(),
            "display_coordinate_bound_ns": time_bound, "table_4_1_published_acquisition_time_ns": None,
            "table_4_1_relative_to_breakdown": None, "actual_breakdown_time_ns": None, "exact_trigger_reference": None,
            "gate_edge_meaning": "geometrical extent of published rectangle; electronics edge definition not recovered"},
    }
    return rows, details


def make_qa(pdf: Path, rows: list[dict], details: dict, png: Path | None, pdf_out: Path | None) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    clip = pymupdf.Rect(80., 222., 530., 525.)
    with pymupdf.open(pdf) as document:
        pix = document[PAGE_INDEX].get_pixmap(matrix=pymupdf.Matrix(2., 2.), clip=clip, alpha=False)
        pixels = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    figure, (native, recovered) = plt.subplots(1, 2, figsize=(13, 5.5), layout="constrained")
    native.imshow(pixels, extent=(clip.x0, clip.x1, clip.y1, clip.y0))
    native.scatter([row["page_x_pt"] for row in rows], [row["page_y_pt"] for row in rows],
                   s=7, color="black", label="Original current vertices")
    legend = details["vector_path"]["later_legend_bounding_rectangle_pt"]
    native.add_patch(plt.Rectangle((legend[0], legend[1]), legend[2]-legend[0], legend[3]-legend[1],
                                   fill=False, edgecolor="navy", linewidth=.8, label="Later white legend bounds"))
    native.set(xlabel="Native PDF x [pt]", ylabel="Native PDF y [pt]", title="C.19 source and encoded current vertices")
    native.legend(loc="lower right", fontsize=8)
    times = np.array([row["time_ns"] for row in rows])
    current = np.array([row["current_MA"] for row in rows])
    bound = details["calibration"]["digitization_bounds"]["current_plus_minus_A"]*1e-6
    recovered.fill_between(times, current-bound, current+bound, color="firebrick", alpha=.2,
                           label="Display envelope; not physical uncertainty")
    recovered.plot(times, current, ".-", ms=3, lw=.8, color="firebrick", label="Original polyline vertices")
    recovered.axvspan(*details["VUV_gate"]["plotted_time_bounds_ns"], color="blueviolet", alpha=.3,
                      label="VUV gate geometry; not breakdown")
    recovered.set(xlabel="Source plotted time [ns]; trigger origin unknown", ylabel="Average current [MA]",
                  title="Shot 17133: published display curve")
    recovered.grid(alpha=.2)
    recovered.legend(fontsize=8)
    figure.suptitle("Reconstructed source vector; raw waveform and physical error bars unavailable", fontsize=11)
    for destination in (png, pdf_out):
        if destination is not None:
            destination.parent.mkdir(parents=True, exist_ok=True)
            figure.savefig(destination, dpi=150)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--csv", type=Path, default=ROOT/"data/dissertation/shot_17133_current.csv")
    parser.add_argument("--metadata", type=Path, default=ROOT/"data/dissertation/shot_17133_current.json")
    parser.add_argument("--qa-plot", type=Path)
    parser.add_argument("--qa-pdf", type=Path)
    args = parser.parse_args()
    try:
        outputs = [path.resolve() for path in (args.csv, args.metadata, args.qa_plot, args.qa_pdf) if path is not None]
        if len(set(outputs)) != len(outputs) or args.pdf.resolve() in outputs:
            raise ValueError("PDF input and output paths must be distinct")
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
                "pdf_filename": "smtrevor_1.pdf", "pdf_sha256": PDF_SHA256, "figure": "C.19",
                "printed_page": 165, "pdf_page": 190, "pdf_page_index_zero_based": PAGE_INDEX,
                "page_content_stream_xrefs": PAGE_CONTENT_XREFS, "decoded_page_content_stream_sha256": PAGE_CONTENT_SHA256,
                "native_image_pixel_hash": None,
                "native_image_pixel_hash_applicability": "not applicable: C.19 contains vector paths and no embedded raster",
                "parser_reference_version": "PyMuPDF 1.28.2", "parser_used_version": pymupdf.VersionBind,
                "identity_guards": "Exact PDF, page content, connected vector geometry, native tick coordinates, axes, gate and later legend are checked"},
            "experiment": {"facility": "Mykonos", "shot": 17133, "foil_material": "SS304", "foil_thickness_um": 200,
                "thickness_evidence": ["Figure C.19 in-plot title", "Figure C.19 caption"], "unresolved_caption_conflicts": []},
            "observable": {"quantity": "published_average_current", "original_units": "MA", "output_units": "A",
                "time_coordinate": "unshifted source plotted time in ns", "time_origin_reference": None,
                "trigger_or_current_start_offset_recovered": False, "raw_oscilloscope_data": False,
                "independent_samples": False, "source_plot_simplification_or_filtering": "unknown",
                "signed_current_vertices_preserved": True, "negative_early_current_vertices_preserved": True,
                "retained_negative_vertex_count": sum(row["current_A"] < 0 for row in rows),
                "negative_current_handling": "no sign floor is applied; the inspected retained C.19 current vertices are all positive"},
            **details,
            "output": {"csv_filename": args.csv.name, "csv_sha256": sha256(args.csv.read_bytes()), "columns": list(FIELDS),
                "samples": len(rows), "time_extent_ns": [rows[0]["time_ns"], rows[-1]["time_ns"]],
                "current_extent_A": [min(row["current_A"] for row in rows), max(row["current_A"] for row in rows)],
                "generator": "tools/recover_dissertation_drive_17133.py", "generator_sha256": sha256(Path(__file__).read_bytes())},
            "limits": [
                "Encoded display vertices are recovered from the published PDF, not original acquisition samples or a fitted drive waveform.",
                "Later white legend overlap is explicitly tested and flagged, including when no current vertices lie behind it.",
                "The two original vertices outside the displayed axes are excluded; no boundary intersections or points at zero or 200 ns are inserted.",
                "No physical time origin, pre-pulse zero, tail, symmetry, smoothing, resampling or extrapolation is imposed.",
                "Signed values are preserved without flooring; zero negative retained vertices is a property of this source trace, not a clipping rule.",
                "The plotted vertices are correlated display geometry; original filtering, simplification, calibration uncertainty and shot variability are unknown.",
                "Display bounds quantify stroke width and tick-fit residuals only, not physical uncertainty or confidence intervals.",
                "The VUV gate is an independent published rectangle. Breakdown, trigger reference and Table 4.1 context are not inferred from this page.",
                "Decimal precision preserves coordinate conversion and does not establish physical measurement accuracy.",
            ]}
        args.metadata.parent.mkdir(parents=True, exist_ok=True)
        args.metadata.write_text(json.dumps(metadata, indent=2, allow_nan=False)+"\n", encoding="utf-8")
        if args.qa_plot is not None or args.qa_pdf is not None:
            make_qa(args.pdf, rows, details, args.qa_plot, args.qa_pdf)
    except (OSError, ValueError, ImportError) as error:
        parser.error(str(error))
    peak = max(rows, key=lambda row: row["current_A"])
    print(json.dumps({"shot": 17133, "samples": len(rows), "csv": str(args.csv),
        "csv_sha256": metadata["output"]["csv_sha256"], "metadata_sha256": sha256(args.metadata.read_bytes()),
        "generator_sha256": metadata["output"]["generator_sha256"], "peak": peak,
        "time_extent_ns": metadata["output"]["time_extent_ns"], "current_extent_A": metadata["output"]["current_extent_A"],
        "legend_covered_vertices": details["vector_path"]["vertices_behind_legend"],
        "VUV_gate_time_bounds_ns": details["VUV_gate"]["plotted_time_bounds_ns"]}, indent=2))


if __name__ == "__main__":
    main()
