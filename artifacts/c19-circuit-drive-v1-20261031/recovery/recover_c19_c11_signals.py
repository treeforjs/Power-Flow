"""Recover the published resistive-voltage and light-signal vector curves of dissertation Figures C.19 and C.11.

Source: T. J. Smith, dissertation (2025), PDF sha256 1d7b1f9c...ed7 (the same PDF that the repository current
recoveries tools/recover_dissertation_drive_17133.py and tools/recover_dissertation_drive.py check).

  Figure C.19 (PDF page index 189), Mykonos shot 17133:
    drawing 21  green   "V_bank - L_0 dI/dt [MV]"   resistive voltage, MV on the shared Signal Height ordinate
    drawing 22  orange  "Normalized APD 2"           dashed, dimensionless
  Figure C.11 (PDF page index 181), Mykonos shot 15210:
    drawing 23  green   "V_bank - L_0 dI/dt [MV]"   (the extracted legend text has lost its opening bracket)
    drawing 24  blue    "2x Normalized SVS Intensity"

Every curve is one continuous straight-segment stroke. Only original encoded vertices inside the displayed axes
are kept; signs and order are preserved; nothing is smoothed, resampled, shifted or extrapolated. Vertices inside the
later white legend rectangle are flagged. The time and ordinate calibrations use the same fixed inspected tick
centres as the repository current recoveries (refit here from the native tick drawings and compared with them).
Bounds = half stroke width mapped through the axis fit + maximum tick-fit residual: a display envelope only, not
physical uncertainty. The plotted time origin is the unshifted published coordinate; it is not a trigger reference.

Usage: python -I recover_c19_c11_signals.py PDF OUTDIR
"""
from __future__ import annotations

import csv, hashlib, json, sys
from pathlib import Path

import numpy as np
import pymupdf

PDF_SHA256 = "1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7"
TIME_TICKS_NS = np.arange(0., 201., 20.)
ORD_TICKS = np.array([0., .2, .4, .6, .8, 1.])

FIGURES = {
    "C.19": {"page_index": 189, "shot": 17133, "content_sha256": "0a6619c5541d64d19d1534e738268c6d7868eca658152967a2f8c8dfe9c57edb",
             "n_drawings": 33, "x_ticks": tuple(range(3, 14)), "y_ticks": tuple(range(14, 20)), "axes": 1, "legend": 27,
             "current_drawing": 20,
             "curves": {21: ("resistive_voltage", "V_bank - L_0 dI/dt [MV]", "MV", 1e6, "V"),
                        22: ("apd2_normalized", "Normalized APD 2", "normalized", 1.0, "1")}},
    "C.11": {"page_index": 181, "shot": 15210, "content_sha256": "1d12ec2ed170923d2a27a3082a0317a50ce048d8f25a7601b21bd65d6c36e50f",
             "n_drawings": 37, "x_ticks": tuple(range(5, 16)), "y_ticks": tuple(range(16, 22)), "axes": 1, "legend": 29,
             "current_drawing": 22,
             "curves": {23: ("resistive_voltage", "V_bank - L_0 dI/dt [MV]", "MV", 1e6, "V"),
                        24: ("svs_2x_normalized", "2x Normalized SVS Intensity", "normalized", 1.0, "1")}},
}


def sha256(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def polyline(d):
    pts = [tuple(d["items"][0][1])]
    for it in d["items"]:
        if it[0] != "l" or tuple(it[1]) != pts[-1]:
            raise ValueError("expected one continuous straight-segment polyline")
        pts.append(tuple(it[2]))
    return pts


def main(pdf: Path, out: Path) -> None:
    if sha256(pdf.read_bytes()) != PDF_SHA256:
        raise ValueError("PDF differs from the inspected dissertation")
    out.mkdir(parents=True, exist_ok=True)
    manifest = {"schema": "c19-circuit-drive-v1-signal-recovery-v1", "source_pdf_sha256": PDF_SHA256,
                "generator": Path(__file__).name, "generator_sha256": sha256(Path(__file__).read_bytes()),
                "status": "vector_recovered_published_display_not_raw_scope_data", "figures": {}}
    with pymupdf.open(pdf) as doc:
        for fig, spec in FIGURES.items():
            page = doc[spec["page_index"]]
            if sha256(page.read_contents()) != spec["content_sha256"]:
                raise ValueError(f"{fig} content stream changed")
            dr = page.get_drawings()
            if len(dr) != spec["n_drawings"]:
                raise ValueError(f"{fig} drawing count changed")
            xt = np.array([dr[i]["items"][0][1].x for i in spec["x_ticks"]])
            yt = np.array([dr[i]["items"][0][1].y for i in spec["y_ticks"]])
            tfit = np.polyfit(xt, TIME_TICKS_NS, 1); yfit = np.polyfit(yt, ORD_TICKS, 1)
            tres = TIME_TICKS_NS - np.polyval(tfit, xt); yres = ORD_TICKS - np.polyval(yfit, yt)
            axes = tuple(dr[spec["axes"]]["rect"]); legend = tuple(dr[spec["legend"]]["rect"])
            fig_rec = {"page_index_zero_based": spec["page_index"], "shot": spec["shot"],
                       "page_content_sha256": spec["content_sha256"],
                       "time_axis": {"ns_per_pt": float(tfit[0]), "intercept_ns": float(tfit[1]),
                                     "max_abs_tick_residual_ns": float(np.max(np.abs(tres)))},
                       "ordinate_axis": {"label": "Signal Height", "units_per_pt": float(yfit[0]), "intercept": float(yfit[1]),
                                         "max_abs_tick_residual": float(np.max(np.abs(yres)))},
                       "axes_rectangle_pt": list(axes), "legend_rectangle_pt": list(legend), "curves": {}}
            for idx, (name, legend_text, unit, si_scale, si_unit) in spec["curves"].items():
                d = dr[idx]
                pts = polyline(d)
                geometry = {"color": list(d["color"]), "width": d["width"], "points": pts}
                gsha = sha256(json.dumps(geometry, sort_keys=True, separators=(",", ":")).encode())
                hw = float(d["width"]) / 2.
                tb = abs(float(tfit[0])) * hw + float(np.max(np.abs(tres)))
                yb = abs(float(yfit[0])) * hw + float(np.max(np.abs(yres)))
                rows, excluded = [], []
                for k, (x, y) in enumerate(pts):
                    if not (axes[0] <= x <= axes[2] and axes[1] <= y <= axes[3]):
                        excluded.append({"source_vertex_index": k, "page_point": [x, y], "reason": "outside_displayed_axes"})
                        continue
                    t = float(np.polyval(tfit, x)); v = float(np.polyval(yfit, y))
                    rows.append({"source_vertex_index": k, "page_x_pt": x, "page_y_pt": y, "time_ns": t, "time_s": t * 1e-9,
                                 "time_bound_ns": tb, "value_plotted": v, "value_SI": v * si_scale,
                                 "value_bound_SI": yb * si_scale,
                                 "occluded_by_legend": int(legend[0] <= x <= legend[2] and legend[1] <= y <= legend[3])})
                times = [r["time_ns"] for r in rows]
                monotone = all(b > a for a, b in zip(times, times[1:]))
                stem = f"fig_{fig.replace('.', '_')}_shot_{spec['shot']}_{name}"
                with (out / (stem + ".csv")).open("w", newline="") as h:
                    w = csv.DictWriter(h, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
                vals = np.array([r["value_SI"] for r in rows])
                fig_rec["curves"][name] = {
                    "drawing_index_zero_based": idx, "drawing_sequence": d["seqno"], "RGB": list(d["color"]),
                    "stroke_width_pt": d["width"], "dashes": d.get("dashes"), "legend_text": legend_text,
                    "plotted_unit": unit, "SI_unit": si_unit, "SI_per_plotted_unit": si_scale,
                    "encoded_vertices": len(pts), "retained_vertices": len(rows), "excluded_vertices": excluded,
                    "vertices_behind_legend": sum(r["occluded_by_legend"] for r in rows),
                    "time_strictly_increasing": monotone, "canonical_geometry_sha256": gsha,
                    "time_extent_ns": [min(times), max(times)], "value_extent_SI": [float(vals.min()), float(vals.max())],
                    "time_bound_ns": tb, "value_bound_SI": yb * si_scale,
                    "csv": stem + ".csv", "csv_sha256": sha256((out / (stem + ".csv")).read_bytes())}
            manifest["figures"][fig] = fig_rec
    manifest["limits"] = [
        "Published display vertices, not oscilloscope samples; original filtering, plot simplification and scope calibration are unknown.",
        "V_bank - L_0 dI/dt is the experimenters' resistive-voltage construction (dissertation Eq. 4.1-4.2, pp. 59-60): V_bank is measured at the "
        "water-line output, L_0 is the total inductance from there through the target (value not published), dI/dt from calibrated B-dot probes; "
        "signals are propagation-time corrected by the experimenters. Its L_0 value and timing corrections are not recovered.",
        "The plotted time coordinate is unshifted; it is the same coordinate as the repository current recoveries of the same figures.",
        "APD and SVS curves are normalized light signals; only their timing (onset) is a comparison quantity.",
        "Bounds are display envelopes (half stroke + tick residual), not physical uncertainty."]
    (out / "signal_recovery_manifest.json").write_text(json.dumps(manifest, indent=1))
    for fig, rec in manifest["figures"].items():
        for name, c in rec["curves"].items():
            print(fig, name, c["encoded_vertices"], c["retained_vertices"], c["time_extent_ns"], c["value_extent_SI"],
                  "monotone", c["time_strictly_increasing"], "legend", c["vertices_behind_legend"])


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
