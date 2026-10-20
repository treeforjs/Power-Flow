#!/usr/bin/env python3
"""Recover published calibration display geometry from dissertation Figs.4.18/19.

The exact inspected PDF, page streams, native drawing geometry and tick marks
are guarded. Original encoded vertices and measured-marker centres are retained;
no spectral resampling, fit, radiance conversion or absolute correction occurs.
The literal ordinate unit is ambiguous: microW/(sr nm mm), not silently mm^2.
Optional source-pixel QA overlays belong outside the repository/delivery archive.
PyMuPDF is required only for PDF recovery; --help and isolated tests need none.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PDF_SHA256 = "1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7"
LITERAL_UNIT = "µW/(sr·nm·mm)"
SCHEMA = "mhdlab-dissertation-calibration-vector-v1"
SOURCE_URL = "https://backend.production.deepblue-documents.lib.umich.edu/server/api/core/bitstreams/24e682a8-896a-410d-85bb-ca4654a9db02/content"
PAGE_STREAM_SHA256 = {
    105: "5111aee2cd4e9be7e1c386e375055d5796ba459d2919ca452c893402f9bf0b22",
    106: "2009bea2193d3a2a4997398ce00ad4c00408a7cac1b3530feb3eae4b3a724a38",
}
NATIVE_DRAWINGS_SHA256 = {
    105: "43a3be1acb11d1073c10a64bc84e90a8f7b88d4f76abedf69f0dd0d57dd39611",
    106: "8c46adff20b1ab5e93a5568d00444ba9a0611c528521c147568043ea9391fcd7",
}
SPECS = {
    "4.18": dict(page=105, curve=17, x_ticks=list(range(2, 9)), y_ticks=list(range(9, 17)),
                 x_values=list(range(110, 171, 10)), y_values=list(range(0, 35001, 5000)),
                 legend=22, logarithmic_y=False),
    "4.19": dict(page=106, curve=99, band=2, markers=list(range(3, 39)),
                 x_ticks=list(range(39, 45)), y_ticks=list(range(45, 51)),
                 x_values=list(range(0, 101, 20)), y_values=list(range(1, 7)),
                 legend=104, logarithmic_y=True),
}
LAMP_FIELDS = ("source_vertex_index", "page_x_pt", "page_y_pt", "wavelength_nm", "plotted_intensity",
               "wavelength_display_lower_nm", "wavelength_display_upper_nm",
               "intensity_display_lower", "intensity_display_upper", "occluded_by_legend")
GAIN_FIELDS = ("series", "record_kind", "source_path_index", "source_element_index", "page_x_pt",
               "page_y_pt", "gain", "log10_plotted_intensity", "plotted_intensity",
               "gain_display_lower", "gain_display_upper", "intensity_display_lower",
               "intensity_display_upper", "gain_region", "occluded_by_legend")


def sha256(payload):
    return hashlib.sha256(payload).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _geometry_value(value):
    if isinstance(value, (tuple, list)):
        return [_geometry_value(item) for item in value]
    if hasattr(value, "x") and hasattr(value, "y"):
        return [float(value.x), float(value.y)]
    if hasattr(value, "x0"):
        return [float(value.x0), float(value.y0), float(value.x1), float(value.y1)]
    return value


def native_geometry(path):
    """Stable explicit subset of native drawing fields, without renderer seqno."""
    return {key: _geometry_value(path[key]) for key in
            ("type", "color", "fill", "width", "closePath", "items")}


def drawings_digest(drawings):
    return sha256(canonical([native_geometry(path) for path in drawings]))


def polyline_vertices(path):
    if not path["items"] or path["items"][0][0] != "l":
        raise ValueError("expected nonempty native line-segment geometry")
    points = [tuple(path["items"][0][1])]
    for item in path["items"]:
        if item[0] != "l" or tuple(item[1]) != points[-1]:
            raise ValueError("expected one continuous original polyline; no joining or interpolation allowed")
        points.append(tuple(item[2]))
    return points


def marker_center(path):
    if len(path["items"]) != 8 or any(item[0] != "c" for item in path["items"]):
        raise ValueError("expected the inspected eight-cubic-segment circular marker")
    rect = tuple(path["rect"])
    if abs((rect[2]-rect[0])-6) > 5e-5 or abs((rect[3]-rect[1])-6) > 5e-5:
        raise ValueError("expected six-point marker diameter")
    return (rect[0]+rect[2])/2, (rect[1]+rect[3])/2


def axis_fit(coordinates, values):
    coordinates, values = np.asarray(coordinates, float), np.asarray(values, float)
    if (coordinates.ndim != 1 or values.shape != coordinates.shape or len(values) < 2
            or not np.all(np.isfinite(coordinates)) or not np.all(np.isfinite(values))
            or len(np.unique(coordinates)) != len(coordinates)):
        raise ValueError("finite distinct tick coordinates and paired values required")
    # Exact least squares on the binary tick inputs avoids LAPACK/NumPy-version
    # differences in final CSV digits. This fits coordinates, not physical data.
    xs = [Fraction.from_float(float(value)) for value in coordinates]
    ys = [Fraction.from_float(float(value)) for value in values]
    x_mean, y_mean = sum(xs)/len(xs), sum(ys)/len(ys)
    slope = sum((x-x_mean)*(y-y_mean) for x,y in zip(xs,ys))/sum((x-x_mean)**2 for x in xs)
    intercept = y_mean-slope*x_mean
    fit = np.array([float(slope),float(intercept)])
    residual = np.array([float(y-(slope*x+intercept)) for x,y in zip(xs,ys)])
    return fit, residual


def _inside(point, rectangle):
    return rectangle[0] <= point[0] <= rectangle[2] and rectangle[1] <= point[1] <= rectangle[3]


def calibrated_points(points, x_fit, y_fit, x_residual, y_residual, half_extent, logarithmic_y):
    dx = abs(float(x_fit[0]))*half_extent + float(max(abs(x_residual)))
    dy = abs(float(y_fit[0]))*half_extent + float(max(abs(y_residual)))
    for x, y in points:
        x_value = float(x_fit[0])*float(x)+float(x_fit[1])
        ordinate = float(y_fit[0])*float(y)+float(y_fit[1])
        if logarithmic_y:
            yield x_value, ordinate, 10**ordinate, x_value-dx, x_value+dx, 10**(ordinate-dy), 10**(ordinate+dy)
        else:
            yield x_value, ordinate, ordinate, x_value-dx, x_value+dx, ordinate-dy, ordinate+dy


def _calibration(drawings, spec):
    x_points = [float(drawings[i]["items"][0][1].x) for i in spec["x_ticks"]]
    y_points = [float(drawings[i]["items"][0][1].y) for i in spec["y_ticks"]]
    x_fit, x_residual = axis_fit(x_points, spec["x_values"])
    y_fit, y_residual = axis_fit(y_points, spec["y_values"])
    details = {
        "coordinate_system": "native PDF points, upper-left origin; no raster calibration",
        "method": "exact-rational least-squares linear map from fixed native major-tick coordinates, coefficients rounded once to float64",
        "x": {"tick_drawing_indices": spec["x_ticks"], "tick_coordinates_pt": x_points,
              "tick_values": spec["x_values"], "slope": float(x_fit[0]), "intercept": float(x_fit[1]),
              "residuals": x_residual.tolist()},
        "y": {"tick_drawing_indices": spec["y_ticks"], "tick_coordinates_pt": y_points,
              "tick_values": spec["y_values"], "mapped_quantity": "log10 plotted intensity" if spec["logarithmic_y"] else "plotted intensity",
              "slope": float(y_fit[0]), "intercept": float(y_fit[1]), "residuals": y_residual.tolist()},
        "display_envelope": "Half stroke width or marker radius plus maximum tick-fit residual; visible glyph extent only, not marker-centre error, physical calibration uncertainty, confidence interval or independent observations",
    }
    return (x_fit, y_fit, x_residual, y_residual), details


def recover(pdf):
    """Return both (rows, metadata) records after validating exact source bytes."""
    payload = Path(pdf).read_bytes()
    if sha256(payload) != PDF_SHA256:
        raise ValueError("PDF SHA256 differs from the inspected dissertation")
    import pymupdf
    recovered = {}
    with pymupdf.open(stream=payload, filetype="pdf") as document:
        if len(document) != 203:
            raise ValueError("expected 203 dissertation PDF pages")
        for figure, spec in SPECS.items():
            page = document[spec["page"]-1]
            if (sha256(page.read_contents()) != PAGE_STREAM_SHA256[spec["page"]]
                    or page.get_images(full=True) or f"Figure {figure}" not in page.get_text()
                    or "mm)]" not in page.get_text()):
                raise ValueError("calibration page content, vector-only format or literal unit changed")
            drawings = page.get_drawings()
            if drawings_digest(drawings) != NATIVE_DRAWINGS_SHA256[spec["page"]]:
                raise ValueError("native drawing geometry changed; inspect before recalibration")
            calibration, axis_metadata = _calibration(drawings, spec)
            legend = tuple(drawings[spec["legend"]]["rect"])
            curve = drawings[spec["curve"]]
            vertices = polyline_vertices(curve)
            path_indices = [spec["curve"]]
            rows = []
            if figure == "4.18":
                if len(vertices) != 243:
                    raise ValueError("expected 243 original lamp-curve vertices")
                values = calibrated_points(vertices, *calibration, curve["width"]/2, False)
                for index, ((x,y), (wavelength, _, intensity, xlo,xhi,ylo,yhi)) in enumerate(zip(vertices,values)):
                    rows.append(dict(zip(LAMP_FIELDS,(index,x,y,wavelength,intensity,xlo,xhi,ylo,yhi,int(_inside((x,y),legend))))))
            else:
                if len(vertices) != 101:
                    raise ValueError("expected 101 original fit vertices")
                groups = [("exponential_fit", "polyline_vertex", spec["curve"], vertices, curve["width"]/2)]
                for index in spec["markers"]:
                    groups.append(("measured_marker", "marker_center", index, [marker_center(drawings[index])], 3.))
                boundary = polyline_vertices(drawings[spec["band"]])
                if len(boundary) != 203 or boundary[0] != boundary[-1]:
                    raise ValueError("expected 203 encoded vertices including repeated closure of fit-error polygon")
                groups.append(("fit_error_boundary", "polygon_boundary_vertex", spec["band"], boundary, .5))
                path_indices.extend(spec["markers"]+[spec["band"]])
                for series, kind, path_index, points, extent in groups:
                    values = calibrated_points(points, *calibration, extent, True)
                    for index, ((x,y),(gain,log_intensity,intensity,xlo,xhi,ylo,yhi)) in enumerate(zip(points,values)):
                        region = "below_measured_gain_range" if xhi<20 else "above_measured_gain_range" if xlo>55 else "within_or_at_reported_measured_range"
                        values_row = (series,kind,path_index,-1 if kind=="marker_center" else index,x,y,gain,
                                      log_intensity,intensity,xlo,xhi,ylo,yhi,region,int(_inside((x,y),legend)))
                        rows.append(dict(zip(GAIN_FIELDS,values_row)))
                if Counter(row["series"] for row in rows) != {"exponential_fit":101,"measured_marker":36,"fit_error_boundary":203}:
                    raise ValueError("unexpected gain series counts")
            metadata = {
                "schema": SCHEMA, "status": "recovered_published_vector_display_data",
                "source": {"author":"Trevor Johannes Smith","year":2025,"figure":figure,
                           "printed_page":spec["page"]-25,"pdf_page_1based":spec["page"],
                           "pdf_url":SOURCE_URL,"pdf_sha256":PDF_SHA256,
                           "page_content_sha256":PAGE_STREAM_SHA256[spec["page"]],
                           "all_native_drawings_sha256":NATIVE_DRAWINGS_SHA256[spec["page"]]},
                "native_hash_format": "sorted compact UTF-8 JSON; each path has type,color,fill,width,closePath,items; coordinates native floats; no renderer sequence numbers",
                "paths": [{"drawing_index":index,"canonical_path_sha256":sha256(canonical(native_geometry(drawings[index]))),
                           "items":len(drawings[index]["items"]),"stroke_width_pt":drawings[index]["width"]} for index in path_indices],
                "calibration": axis_metadata,
                "observable": {"x_label":"Wavelength [nm]" if figure=="4.18" else "Gain",
                               "y_label":"Intensity [µW/(sr·nm·mm)]","literal_ordinate_unit":LITERAL_UNIT,
                               "quantity_kind":"encoded lamp display curve" if figure=="4.18" else "encoded gain display fit, marker centres and uninterpreted band boundary",
                               "is_raw_detector_count":False,"is_calibration_certificate":False,
                               "absolute_factor_resolved":False,"SI_radiance_conversion_applied":False},
                "recovery": {"original_vertices_only":True,"resampling":False,"smoothing":False,"new_fit":False,
                             "legend_drawing_index":spec["legend"],"legend_rectangle_pt":list(legend),
                             "legend_symbols_excluded":True,"occluded_rows":sum(row["occluded_by_legend"] for row in rows),
                             "rows":len(rows),"rows_are_independent_measurements":False},
                "unresolved_absolute_calibration": [
                    "Literal single-mm ordinate conflicts with an ordinary area-based spectral radiance interpretation and quoted mm² sr étendue; no mm² replacement adopted",
                    "Original lamp certificate/traceability and detector-count-to-photon calibration files not recovered",
                    "Gain-marker exposure assignments, normalization and numerical original fit coefficients/covariance not supplied in inspected section",
                    "Mykonos numerical étendue, effective per-spatial-bin geometry and wavelength throughput factors not supplied in inspected section",
                    "No numerical absolute correction inferred from these displayed curves",
                ],
                "source_context": {"lamp":"McPherson 632 D2","lamp_slit_width_um":50,
                                   "reported_lamp_etendue_mm2_sr":.037,"reported_gain_measurement_range":[20,55],
                                   "reported_lamp_gate_ns":[3,5,10,100],"shot_gate_ns":3,
                                   "reported_spectral_sampling_nm_per_pixel":.03,
                                   "per_marker_gate_assignment":None,"fit_error_statistical_definition":None,
                                   "gain_extrapolation_is_validation":False},
            }
            if figure == "4.18":
                metadata["source_context"].update(caption_wavelength_range_nm=[100,200], plotted_tick_range_nm=[110,170])
                backwards = [{"from_vertex":index,"to_vertex":index+1,
                              "wavelength_from_nm":first["wavelength_nm"],"wavelength_to_nm":second["wavelength_nm"]}
                             for index,(first,second) in enumerate(zip(rows,rows[1:]))
                             if second["wavelength_nm"] <= first["wavelength_nm"]]
                metadata["recovery"].update(wavelength_strictly_increasing=not backwards,
                                          nonincreasing_wavelength_segments=backwards,
                                          interpolation_ready=False,
                                          source_vertex_order_preserved=True)
            else:
                metadata["recovery"].update(series_counts=dict(Counter(row["series"] for row in rows)),
                                         fit_band_interpretation="Literal Fit Error polygon boundary; no confidence level or probability distribution inferred",
                                         excluded_legend_marker_drawing=107,
                                         marker_radius_pt=3., marker_centre_method="midpoint of native circular-glyph bounding rectangle")
            recovered[figure] = rows, metadata
    return recovered


def csv_bytes(rows, fields):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: format(value,".17g") if isinstance(value,float) else value for key,value in row.items()})
    return stream.getvalue().encode("utf-8")


def make_qa(pdf, figure, rows, destination):
    import pymupdf
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    clip = pymupdf.Rect(70,90,510,430)
    payload = Path(pdf).read_bytes()
    if sha256(payload) != PDF_SHA256:
        raise ValueError('QA source PDF changed after recovery')
    with pymupdf.open(stream=payload,filetype='pdf') as document:
        pix = document[SPECS[figure]["page"]-1].get_pixmap(matrix=pymupdf.Matrix(2,2),clip=clip,alpha=False)
        pixels = np.frombuffer(pix.samples,np.uint8).reshape(pix.height,pix.width,pix.n)
    fig,(left,right)=plt.subplots(1,2,figsize=(13,5),layout="constrained")
    left.imshow(pixels,extent=(clip.x0,clip.x1,clip.y1,clip.y0))
    if figure=="4.18":
        left.scatter([r['page_x_pt'] for r in rows],[r['page_y_pt'] for r in rows],s=4,color='black')
        right.plot([r['wavelength_nm'] for r in rows],[r['plotted_intensity'] for r in rows],'.-',ms=2,lw=.6)
        right.set_xlabel('Published wavelength [nm]')
    else:
        for series,color in (("exponential_fit","blue"),("measured_marker","red"),("fit_error_boundary","gray")):
            selected=[r for r in rows if r['series']==series]
            left.scatter([r['page_x_pt'] for r in selected],[r['page_y_pt'] for r in selected],s=5,facecolors='none',edgecolors=color,linewidths=.5)
            right.plot([r['gain'] for r in selected],[r['plotted_intensity'] for r in selected],'.' if series=='measured_marker' else '-',ms=3,lw=.7,label=series)
        right.set_yscale('log'); right.set_xlabel('Published gain setting'); right.legend(fontsize=7)
    left.set(title='Native source pixels and recovered geometry',xlabel='PDF x [pt]',ylabel='PDF y [pt]')
    right.set(title='Published display data; absolute factor unresolved',ylabel='Literal ordinate [µW/(sr·nm·mm)]')
    fig.suptitle(f'Figure {figure}; no original detector counts or absolute correction')
    fig.savefig(destination,dpi=150); plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pdf',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,default=ROOT/'data/dissertation')
    parser.add_argument('--qa-dir',type=Path)
    args=parser.parse_args()
    results=recover(args.pdf)
    if args.qa_dir:
        args.qa_dir.mkdir(parents=True,exist_ok=True)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    outputs=[]
    for figure,(rows,metadata) in results.items():
        stem='figure_4_18_calibration' if figure=='4.18' else 'figure_4_19_gain'
        output=args.output_dir/(stem+'.csv')
        sidecar=args.output_dir/(stem+'.json')
        if args.pdf.resolve() in (output.resolve(),sidecar.resolve()):
            raise ValueError('source PDF cannot be an output destination')
        payload=csv_bytes(rows,LAMP_FIELDS if figure=='4.18' else GAIN_FIELDS)
        metadata['output']={'csv_filename':output.name,'csv_sha256':sha256(payload),'columns':list(LAMP_FIELDS if figure=='4.18' else GAIN_FIELDS),'rows':len(rows)}
        if args.qa_dir:
            qa=args.qa_dir/(stem+'_qa.png')
            make_qa(args.pdf,figure,rows,qa)
            metadata['QA_overlay']={'path':str(qa),'sha256':sha256(qa.read_bytes()),'contains_source_pixels':True,'redistribute_in_public_archive':False}
        output.write_bytes(payload)
        sidecar.write_text(json.dumps(metadata,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
        outputs.append({'figure':figure,'rows':len(rows),'csv':str(output),'sha256':sha256(payload)})
    print(json.dumps(outputs,indent=2))


if __name__=='__main__':
    main()
