#!/usr/bin/env python3
"""Recover a partial observed visible SVS lineout from dissertation Fig. 4.15.

This is deterministic selection of literal green pixels in the inspected native
image, not recovery of raw counts, spectral bins, radiance or independent data.
The source PDF stays external; source-pixel QA must stay outside delivery files.
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


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "mhdlab-published-observed-visible-spectrum-v1"
PDF_SHA256 = "1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7"
RGB_SHA256 = "27a9beba2c5aa9edb0de720ff0f50866157eb845b12c1318864279f1b16d5276"
PAGE_CONTENT_SHA256 = "064dfb491db122b454c86a3ee2a09e74fe40b3bdeb1c4928247969c180ae5442"
SOURCE_URL = "https://backend.production.deepblue-documents.lib.umich.edu/server/api/core/bitstreams/24e682a8-896a-410d-85bb-ca4654a9db02/content"
IMAGE_SIZE = (1040, 1320)
GREEN_RGB = (60, 179, 113)
X_TICK_ROIS = ((162,795,169,803),(254,795,261,803),(347,795,354,803),
               (439,795,446,803),(532,795,539,803),(625,795,632,803),
               (718,795,725,803),(811,795,818,803),(904,795,911,803))
Y_TICK_ROIS = ((938,542,946,551),(938,624,946,633),(938,706,946,715))
X_TICK_COORDS = (163.5,256.5,349.5,442.5,535.5,628.5,721.5,814.5,907.5)
Y_TICK_COORDS = (545.5,627.5,710.5)
X_TICK_VALUES = tuple(range(350,751,50))
Y_TICK_VALUES = (150,100,50)
LEGEND_RECTANGLE = (626,510,926,611)
Y_ROI = (516,788)
SAMPLE_COLUMNS = tuple(range(220,841,4))
FIELDS = ("pixel_x", "pixel_y", "matched_green_pixels", "wavelength_nm",
          "observed_intensity_arb", "wavelength_readout_lower_nm", "wavelength_readout_upper_nm",
          "intensity_readout_lower_arb", "intensity_readout_upper_arb")


def sha256(payload):
    return hashlib.sha256(payload).hexdigest()


def axis_fit(pixel_coords, tick_values):
    """Fit tick geometry exactly, then round coefficients once to float64."""
    xs, ys = np.asarray(pixel_coords, float), np.asarray(tick_values, float)
    if (xs.ndim != 1 or ys.shape != xs.shape or xs.size < 2
            or not np.all(np.isfinite(xs)) or not np.all(np.isfinite(ys))
            or np.unique(xs).size != xs.size):
        raise ValueError("finite distinct tick coordinates and paired values required")
    xf = [Fraction.from_float(float(x)) for x in xs]
    yf = [Fraction.from_float(float(y)) for y in ys]
    xm, ym = sum(xf)/len(xf), sum(yf)/len(yf)
    a = sum((x-xm)*(y-ym) for x,y in zip(xf,yf))/sum((x-xm)**2 for x in xf)
    b = ym-a*xm
    return np.array([float(a),float(b)]), np.array([float(y-a*x-b) for x,y in zip(xf,yf)])


def load_source_image(pdf):
    """Decode only the byte-identical inspected dissertation/native image."""
    payload = Path(pdf).read_bytes()
    if sha256(payload) != PDF_SHA256:
        raise ValueError("PDF SHA256 differs from inspected official dissertation")
    try:
        import pymupdf
    except ImportError:
        import fitz as pymupdf
    with pymupdf.open(stream=payload, filetype="pdf") as document:
        if len(document) != 203:
            raise ValueError("unexpected PDF page count")
        page = document[99]
        if sha256(page.read_contents()) != PAGE_CONTENT_SHA256:
            raise ValueError("PDF page content differs from inspected Figure 4.15")
        if 2554 not in [item[0] for item in page.get_images(full=True)]:
            raise ValueError("expected native Figure 4.15 image absent")
        pixmap = pymupdf.Pixmap(document,2554)
        if (pixmap.width,pixmap.height) != IMAGE_SIZE or pixmap.n != 3 or pixmap.alpha:
            raise ValueError("expected native 1040 by 1320 RGB image without alpha")
        image = np.frombuffer(pixmap.samples,dtype=np.uint8).reshape(1320,1040,3).copy()
    _guard_image(image)
    return image, {"url":SOURCE_URL, "pdf_sha256":PDF_SHA256,
                   "page_content_sha256":PAGE_CONTENT_SHA256, "figure":"4.15",
                   "printed_page":75, "pdf_page_1based":100, "native_image_xref":2554,
                   "native_image_size_px":list(IMAGE_SIZE), "native_rgb_sha256":RGB_SHA256}


def _guard_image(image):
    if (not isinstance(image,np.ndarray) or image.dtype != np.uint8
            or image.shape != (IMAGE_SIZE[1],IMAGE_SIZE[0],3)
            or sha256(image.tobytes(order="C")) != RGB_SHA256):
        raise ValueError("native RGB image shape, type or SHA256 differs from inspected source")


def green_mask(image, exclusions=()):
    """Return literal series-color pixels; rectangle endpoints are inclusive."""
    if not isinstance(image,np.ndarray) or image.dtype != np.uint8 or image.ndim != 3 or image.shape[2] != 3:
        raise ValueError("RGB uint8 image required")
    result = np.all(image == GREEN_RGB,axis=2)
    for rectangle in exclusions:
        if len(rectangle) != 4 or any(isinstance(v,bool) or not isinstance(v,(int,np.integer)) for v in rectangle):
            raise ValueError("integer inclusive exclusion rectangle required")
        x0,y0,x1,y1 = rectangle
        if not (0 <= x0 <= x1 < image.shape[1] and 0 <= y0 <= y1 < image.shape[0]):
            raise ValueError("exclusion rectangle outside image")
        result[y0:y1+1,x0:x1+1] = False
    return result


def column_center(mask, pixel_x, y_roi, min_pixels=3, max_pixels=25):
    """Select one bounded contiguous color run, rejecting ambiguous columns."""
    if not isinstance(mask,np.ndarray) or mask.ndim != 2 or mask.dtype != np.bool_:
        raise ValueError("two-dimensional boolean mask required")
    if (isinstance(pixel_x,bool) or not isinstance(pixel_x,(int,np.integer))
            or not 0 <= pixel_x < mask.shape[1] or len(y_roi) != 2
            or any(isinstance(v,bool) or not isinstance(v,(int,np.integer)) for v in y_roi)
            or not 0 <= y_roi[0] < y_roi[1] < mask.shape[0]
            or isinstance(min_pixels,bool) or isinstance(max_pixels,bool)
            or not isinstance(min_pixels,int) or not isinstance(max_pixels,int)
            or not 1 <= min_pixels <= max_pixels):
        raise ValueError("valid integer column, inclusive ROI and run-size bounds required")
    y0,y1 = y_roi
    rows = np.flatnonzero(mask[y0:y1+1,pixel_x])+y0
    if rows.size == 0:
        raise ValueError("no_green_pixels")
    if rows.size < min_pixels:
        raise ValueError("too_few_green_pixels")
    if rows.size > max_pixels:
        raise ValueError("broad_projected_green_run")
    if np.any(np.diff(rows) != 1):
        raise ValueError("separated_green_runs")
    if rows[0] == y0 or rows[-1] == y1:
        raise ValueError("green_run_touches_roi_boundary")
    return {"pixel_y":float(np.median(rows)),"matched_pixels":int(rows.size),
            "y_first":int(rows[0]),"y_last":int(rows[-1])}


def _tick_centres(image, rois, coordinate):
    centres=[]
    for x0,y0,x1,y1 in rois:
        yy,xx = np.where(np.all(image[y0:y1,x0:x1] < 45,axis=2))
        if xx.size == 0:
            raise ValueError("inspected black tick absent")
        centres.append(float(np.median(xx+x0 if coordinate == "x" else yy+y0)))
    return centres


def recover_spectrum(image):
    """Return retained native-pixel samples and all explicit omitted columns."""
    _guard_image(image)
    xt = _tick_centres(image,X_TICK_ROIS,"x")
    yt = _tick_centres(image,Y_TICK_ROIS,"y")
    if tuple(xt) != X_TICK_COORDS or tuple(yt) != Y_TICK_COORDS:
        raise ValueError("native tick centres differ from inspected calibration")
    xfit,xres = axis_fit(xt,X_TICK_VALUES)
    yfit,yres = axis_fit(yt,Y_TICK_VALUES)
    mask = green_mask(image,(LEGEND_RECTANGLE,))
    dx = abs(float(xfit[0]))+float(max(abs(xres)))
    dyres = float(max(abs(yres)))
    rows,skipped = [],[]
    for x in SAMPLE_COLUMNS:
        try:
            point = column_center(mask,x,Y_ROI)
        except ValueError as error:
            skipped.append({"pixel_x":x,"reason":str(error),
                            "matched_green_pixels":int(np.count_nonzero(mask[Y_ROI[0]:Y_ROI[1]+1,x]))})
            continue
        wavelength = float(xfit[0])*x+float(xfit[1])
        intensity = float(yfit[0])*point["pixel_y"]+float(yfit[1])
        # Decreasing y map: full visible run plus one pixel gives an envelope.
        ilo = float(yfit[0])*(point["y_last"]+1)+float(yfit[1])-dyres
        ihi = float(yfit[0])*(point["y_first"]-1)+float(yfit[1])+dyres
        rows.append(dict(zip(FIELDS,(x,point["pixel_y"],point["matched_pixels"],wavelength,intensity,
                                     wavelength-dx,wavelength+dx,ilo,ihi))))
    details = {
        "coordinate_system":"native RGB integer pixel centres, upper-left origin; no rescaling",
        "axis_calibration":{
            "method":"exact-rational least-squares on inspected black major-tick centres; coefficients rounded once to float64",
            "roi_convention":"tick ROI x0,y0,x1,y1 with exclusive upper bounds",
            "x":{"tick_rois_px":X_TICK_ROIS,"tick_centres_px":xt,"tick_values_nm":X_TICK_VALUES,
                 "slope_nm_per_pixel":float(xfit[0]),"intercept_nm":float(xfit[1]),"residuals_nm":xres.tolist()},
            "y":{"tick_rois_px":Y_TICK_ROIS,"tick_centres_px":yt,"tick_values_arb":Y_TICK_VALUES,
                 "slope_arb_per_pixel":float(yfit[0]),"intercept_arb":float(yfit[1]),"residuals_arb":yres.tolist()}},
        "sampling":{"green_rgb":GREEN_RGB,"column_start_px":220,"column_stop_inclusive_px":840,
                    "column_stride_px":4,"candidate_columns":len(SAMPLE_COLUMNS),
                    "y_roi_inclusive_px":Y_ROI,"legend_exclusion_inclusive_px":LEGEND_RECTANGLE,
                    "min_green_pixels":3,"max_green_pixels":25,"retained_samples":len(rows),
                    "skipped_columns":skipped,"interpolation":False,"smoothing":False,
                    "rule":"literal color, one contiguous run per column; exclude broad, separate, absent, small or boundary-touching runs"},
        "readout_envelopes":{
            "wavelength":"plus/minus one native pixel and maximum wavelength tick-fit residual",
            "intensity":"entire selected literal-green pixel run padded one pixel at each end plus maximum intensity tick-fit residual",
            "interpretation":"raster glyph/readout extents only; not physical error bars, independent noise, confidence intervals or covariance"},
    }
    return rows,details


def serialize_csv(rows):
    buffer=io.StringIO(newline="")
    writer=csv.DictWriter(buffer,fieldnames=FIELDS,lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key:format(value,".17g") if isinstance(value,float) else value for key,value in row.items()})
    return buffer.getvalue().encode("utf-8")


def recover(pdf):
    image,source=load_source_image(pdf)
    rows,details=recover_spectrum(image)
    metadata={
        "schema":SCHEMA,"parameter_status":"reconstructed_published_observed_visible_trace",
        "source":source,
        "observable":{"diagnostic":"SVS visible streak spectroscopy","panel":"middle linear lineout panel",
                      "series_literal":"SVS Lineout 3: 330 ns","axis_label_literal":"Intensity [arb.]",
                      "wavelength_unit_literal":"nm","intensity_unit":"arbitrary",
                      "raw_detector_counts":False,"absolute_radiance":False,
                      "transition_density":False,"inferred_covariance":False},
        "experimental_context":{"facility":"Mykonos","shot":15765,"foil_material_literal":"SS304",
                                "foil_thickness_um_caption":400,"visible_lineout_time_ns_literal":330,
                                "visible_lineout_time_origin":"not recovered as a measured trigger or breakdown reference",
                                "visible_exposure_or_averaging_width_ns":None,
                                "same_shot_vuv_context_only":{"time_before_breakdown_ns_table_4_1":29,"diagnostic_gate_ns":3,
                                                              "is_this_visible_acquisition":False},
                                "collection_region":"SVS side-on gap view; distinct from VUV collection",
                                "collection_geometry_reported":{
                                    "source_printed_page":71,"source_pdf_page_1based":96,
                                    "setup_figure":"4.13","setup_figure_printed_page":72,
                                    "setup_figure_shot_range":[15750,15802],
                                    "focal_point_distance_from_cathode_um":450,
                                    "focal_spot_radius_um":175,
                                    "nearest_reported_collection_distance_um":275,
                                    "primary_lens_focal_length_mm":70.3,
                                    "object_distance_from_primary_optic_mm":70.3,
                                    "reported_second_focusing_optic_distance_mm":40.2,
                                    "reported_demagnification":0.57,
                                    "spatial_kernel_or_etendue_inferred":False},
                                "source_emission_gate_inferred":False},
        "recovery":details,
        "limitations":["Partial correlated display samples from a published raster, not original independent detector samples.",
                       "Broad green strokes at sharp features are omitted; peak maxima, missing segments and integrated line strengths are not recovered.",
                       "Spectral bin widths, instrumental resolution, raw counts, exposure averaging, background, gain and covariance are unavailable here.",
                       "The 330-ns visible label and the same-shot 29-ns-before-breakdown VUV record are separate diagnostic contexts.",
                       "No absolute amplitude correction, source-emission gate, temperature, species identification or surface coverage is inferred."]}
    return rows,metadata


def qa_plot(image,rows,metadata,path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    fig,(left,right)=plt.subplots(1,2,figsize=(15,5),layout="constrained")
    left.imshow(image)
    left.scatter([r["pixel_x"] for r in rows],[r["pixel_y"] for r in rows],s=8,c="magenta",label="retained centres")
    x0,y0,x1,y1=LEGEND_RECTANGLE
    left.add_patch(Rectangle((x0,y0),x1-x0+1,y1-y0+1,fill=False,color="black",ls="--"))
    for item in metadata["recovery"]["sampling"]["skipped_columns"]:
        left.axvline(item["pixel_x"],color="red",alpha=.2,lw=.5)
    left.set(xlim=(130,945),ylim=(805,495),title="Native source pixels + selections (source-local QA)")
    left.legend(fontsize="small")
    xx=np.array([r["wavelength_nm"] for r in rows]); yy=np.array([r["observed_intensity_arb"] for r in rows])
    right.vlines(xx,[r["intensity_readout_lower_arb"] for r in rows],
                 [r["intensity_readout_upper_arb"] for r in rows],color="lightgray",lw=1)
    right.scatter(xx,yy,s=8,color=np.array(GREEN_RGB)/255)
    right.set(xlabel="Visible wavelength [nm]",ylabel="Intensity [arb.]",
              title="Partial 330-ns observed lineout; glyph envelopes")
    right.grid(alpha=.2)
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=160)
    plt.close(fig)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf",type=Path,required=True,help="external exact inspected dissertation PDF")
    parser.add_argument("--csv",type=Path,default=ROOT/"data/dissertation/observed_spectrum_4_15.csv")
    parser.add_argument("--metadata",type=Path,default=ROOT/"data/dissertation/observed_spectrum_4_15.json")
    parser.add_argument("--qa-plot",type=Path,help="source-pixel QA outside repository/delivery archive")
    args=parser.parse_args(argv)
    paths=[args.pdf,args.csv,args.metadata]+([args.qa_plot] if args.qa_plot else [])
    if len({p.resolve() for p in paths}) != len(paths):
        parser.error("source and output paths must be distinct")
    if args.qa_plot and args.qa_plot.resolve().is_relative_to(ROOT):
        parser.error("source-pixel QA must be outside repository")
    rows,metadata=recover(args.pdf)
    payload=serialize_csv(rows)
    metadata["artifact"]={"csv_file":args.csv.name,"csv_sha256":sha256(payload),
                          "rows":len(rows),"columns":list(FIELDS),"numeric_format":".17g float64; integer pixels/counts"}
    if args.qa_plot:
        image,_=load_source_image(args.pdf)
        qa_plot(image,rows,metadata,args.qa_plot)
        metadata["source_local_qa"]={"file":args.qa_plot.name,"sha256":sha256(args.qa_plot.read_bytes()),
                                     "contains_source_pixels":True,"redistribution_in_delivery":False}
    args.csv.parent.mkdir(parents=True,exist_ok=True)
    args.metadata.parent.mkdir(parents=True,exist_ok=True)
    args.csv.write_bytes(payload)
    args.metadata.write_text(json.dumps(metadata,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"retained_samples":len(rows),"skipped_columns":len(metadata["recovery"]["sampling"]["skipped_columns"]),
                      "csv_sha256":sha256(payload),"metadata_sha256":sha256(args.metadata.read_bytes())}))


if __name__ == "__main__":
    main()
