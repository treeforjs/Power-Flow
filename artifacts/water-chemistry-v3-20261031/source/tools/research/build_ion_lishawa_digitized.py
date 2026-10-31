"""Digitized cross sections from Lishawa et al., J. Chem. Phys. 93, 3196 (1990), Figs. 5 and 7 (DTIC reprint
ADA232204, native 300 dpi 1-bit page scans extracted with pdfimages; pages 7 and 8 of the PDF).

The source tabulates no cross sections; these numbers are DIGITIZED FROM FIGURES. Marker centres (pixel
coordinates relative to the crops given below) were found with fig_markers.py (fixed-radius Hough transform for
open circles, erosion cores for filled markers) and fig_digitize.py (connected-component bounding boxes), and
every accepted marker was checked visually on an overlay. Overlapping markers whose individual centres could
not be separated unambiguously were NOT digitized; their approximate locations are listed in the notes.
Axis calibration: least-squares line through the tick-mark pixel centres listed below.

Usage: python3 -I build_ion_lishawa_digitized.py <extracted_dir>
"""
import json
import os
import sys

SRC = "ion-lishawa-1990-h2oplus-h2o"
CIT = "C. R. Lishawa, R. A. Dressler, J. A. Gardner, R. H. Salter, E. Murad, J. Chem. Phys. 93, 3196-3206 (1990)"


def fit(pairs):
    n = len(pairs)
    sx = sum(v for v, p in pairs)
    sp = sum(p for v, p in pairs)
    sxx = sum(v * v for v, p in pairs)
    sxp = sum(v * p for v, p in pairs)
    b = (n * sxp - sx * sp) / (n * sxx - sx * sx)  # px per unit
    a = (sp - b * sx) / n                          # px at value 0
    resid = max(abs(p - (a + b * v)) for v, p in pairs)
    return a, b, resid


# Fig. 5: crop (380,2080)-(1160,2960) of PDF page 7 scan
F5_XT = [(1, 51.5), (2, 71.5), (3, 92.5), (4, 112.5), (5, 133.0), (7, 174.0), (8, 195.0), (9, 216.0), (10, 236.5),
         (11, 256.5), (12, 277.0), (13, 298.0), (14, 318.5), (16, 358.5), (17, 379.5), (19, 421.0), (20, 441.5),
         (21, 461.5), (23, 502.5), (24, 523.5), (26, 564.0), (29, 626.0), (31, 667.0), (32, 687.5), (33, 708.5),
         (34, 729.0)]
F5_YT = [(2, 685.5), (5, 613.0), (6, 589.5), (7, 565.5), (8, 541.5), (9, 517.5), (10, 493.5), (11, 470.0),
         (12, 446.0), (13, 422.0), (14, 398.5), (15, 374.5), (16, 351.0), (17, 327.5), (18, 303.5), (19, 279.5),
         (21, 232.0), (22, 208.5), (23, 184.5), (24, 160.5), (25, 137.0), (26, 113.0), (28, 65.5), (29, 41.5)]
F5_CIRCLES = [(56, 185), (56, 262), (56, 275), (82, 529), (87, 578), (111, 565), (139, 587), (197, 639), (225, 670),
              (231, 647), (256, 672), (257, 642), (314, 696), (317, 669), (348, 685), (371, 698), (426, 712),
              (435, 702), (546, 709), (595, 717), (652, 721)]
F5_TRIANGLES = [(59.0, 44.0, ""), (70.0, 296.0, ""), (78.0, 418.5, ""), (134.0, 618.5, ""), (188.0, 684.5, ""),
                (238.5, 682.5, "overlaps open circles; centre from erosion core shifted 1.5 px up"),
                (294.0, 711.5, ""),
                (339.5, 718.5, "sits on the 15 eV axis tick; centre from erosion core shifted 1.5 px up"),
                (404.0, 721.0, "")]
# Fig. 7: crop (150,2050)-(1150,2960) of PDF page 8 scan
F7_XT = [(1, 160.5), (2, 181.0), (3, 201.0), (4, 222.0), (5, 242.5), (6, 263.0), (7, 283.0), (8, 303.5), (9, 325.0),
         (10, 345.0), (11, 365.5), (12, 385.5), (13, 406.0), (14, 426.5), (15, 447.0), (16, 467.0), (17, 488.0),
         (18, 509.0), (19, 529.0), (20, 549.5), (21, 569.5), (22, 590.5), (23, 611.0), (24, 631.5), (25, 652.0),
         (26, 672.5), (27, 693.0), (28, 714.0), (29, 734.5), (30, 754.5), (31, 775.5), (32, 796.0), (33, 816.5),
         (34, 837.0)]
F7_YT = [(2, 748.5), (4, 724.5), (8, 677.0), (10, 653.0), (12, 629.5), (14, 605.0), (16, 581.0), (18, 557.5),
         (20, 533.5), (22, 510.0), (24, 486.5), (26, 462.5), (28, 438.5), (30, 414.5), (32, 391.5), (34, 367.5),
         (36, 344.5), (38, 320.5), (40, 296.5), (44, 249.0), (46, 225.5), (48, 201.5), (50, 178.0), (52, 154.0),
         (54, 130.5), (56, 106.5), (58, 82.5)]
F7_CIRCLES = [(164, 600), (165, 582), (191, 611), (197, 668), (220, 647), (249, 648), (334, 674), (340, 659),
              (365, 670), (366, 653), (395, 670), (426, 666), (456, 658), (533, 706), (541, 687), (600, 675),
              (699, 722), (760, 706)]
F7_DIAMONDS = [(156.0, 278.5, "overlaps an open circle near 42e-16 cm2; centre from erosion core"),
               (164.5, 506.5, ""), (174.0, 618.5, ""), (232.0, 656.0, ""), (321.5, 689.0, ""), (503.5, 694.0, ""),
               (550.0, 703.0, "")]

F5_UNRESOLVED = ("Not digitized (overlapping markers): open circle partly hidden behind the triangle near E~4.2 eV, "
                 "sigma~6.2; two fused open circles near E~6.6 eV, sigma~3.9-4.2; overlapping open circle/triangle "
                 "groups near E~22 eV (sigma~0.9-1.0), ~25 eV (second circle, sigma~0.9) and ~27.5 eV (sigma~0.6-0.9); "
                 "the triangle overlapped by circles near E~3.6 eV, sigma~6.3.")
F7_UNRESOLVED = ("Not digitized (overlapping markers): the open circle near E~0.8 eV, sigma~42 (behind a diamond); "
                 "stacked diamonds near E~2.1 eV (sigma~9.3-10.7); circle/diamond groups near E~4 eV (sigma~9.5), "
                 "E~6.7-7.6 eV (sigma~5.5-9.5), E~11 eV (diamonds sigma~5-6), E~13.5 eV (sigma~6-8), "
                 "E~15.6-16.7 eV (diamonds sigma~6), E~21.6-22 eV (sigma~4.5-7), E~24.7 eV (sigma~5.5-8) and "
                 "E~27 eV (sigma~6-7).")


def convert(markers, xt, yt):
    ax, bx, rx = fit(xt)
    ay, by, ry = fit(yt)
    pts = []
    for m in markers:
        px, py = m[0], m[1]
        flag = m[2] if len(m) > 2 else ""
        e = (px - ax) / bx
        s = (py - ay) / by
        pts.append((round(e, 2), round(s, 2), flag, px, py))
    pts.sort()
    return pts, (ax, bx, rx, ay, by, ry)


def table(tid, fig, page, process, pts, cal, unit_px_x, unit_px_y, unresolved, extra):
    ax, bx, rx, ay, by, ry = cal
    dx = 1.5 / bx
    dy = 1.5 / abs(by)
    return {
        "table_id": tid, "source_id": SRC, "citation": CIT, "doi": "10.1063/1.458852",
        "locator": "J. Chem. Phys. 93, %s (DTIC reprint ADA232204, PDF page %d)" % (fig, page),
        "process": process, "quantity": "cross_section_m2", "x_name": "energy_eV_cm",
        "x": [p[0] for p in pts],
        "y": [round(p[1] * 1e-20, 24) for p in pts],
        "y_uncertainty": None,
        "digitization_uncertainty": {"x_eV": round(dx, 3), "y_m2": round(dy * 1e-20, 24),
                                     "basis": "+/-1.5 px marker-centre reading at 300 dpi; calibration max tick residual x %.2f px, y %.2f px" % (rx, ry)},
        "point_flags": [p[2] for p in pts],
        "marker_pixel_centres_crop": [[p[3], p[4]] for p in pts],
        "axis_calibration": {"x_px_at_0": round(ax, 3), "x_px_per_eV": round(bx, 4),
                             "y_px_at_0": round(ay, 3), "y_px_per_1e-16cm2": round(by, 4)},
        "original_units": "1e-16 cm2 (figure axis); x is centre-of-mass collision energy in eV",
        "conversion_to_SI": "sigma[m2] = sigma[1e-16 cm2] * 1e-20",
        "tag": "measured",
        "extraction_method": ("DIGITIZED FROM FIGURE (no tabulated cross sections exist in the source): automated marker "
                              "detection on the native 300 dpi 1-bit scan (fig_markers.py / fig_digitize.py) with visual "
                              "verification of every accepted marker; overlapping markers that could not be separated were "
                              "omitted."),
        "notes": (extra + " The source gives no numerical uncertainty for these cross sections and shows no error bars in "
                  "this figure. Each point is an individual measurement as plotted (not averaged or smoothed here). "
                  "Primary ions were made by 20-30 eV electron impact; possible vibrational excitation not excluded. "
                  + unresolved),
    }


def main():
    out = sys.argv[1]
    p5c, cal5 = convert(F5_CIRCLES, F5_XT, F5_YT)
    p5t, _ = convert(F5_TRIANGLES, F5_XT, F5_YT)
    p7c, cal7 = convert(F7_CIRCLES, F7_XT, F7_YT)
    p7d, _ = convert(F7_DIAMONDS, F7_XT, F7_YT)
    print("Fig5 calibration", cal5)
    print("Fig7 calibration", cal7)
    tabs = [
        table("ion-lishawa1990-fig5-proton-transfer-3a-sigma-digitized", "p. 3200, Fig. 5 (open circles)", 7,
              "H2O+ + D2(18)O -> OH + HD2(18)O+ (proton transfer, reaction 3a)", p5c, cal5, 0, 0, F5_UNRESOLVED,
              "Cross sections corrected by the source for a 20.8% collection efficiency (isotropic product ions). The "
              "source states proton transfer is >90% of oxonium production and that, above ~4 eV, reaction 3a exceeds 3b "
              "(isotope effect), reversed below ~4 eV. Target is D2-18O (mass 22), not H2O."),
        table("ion-lishawa1990-fig5-deuteron-transfer-3b-sigma-digitized", "p. 3200, Fig. 5 (filled triangles)", 7,
              "D2O+ + D2(18)O -> OD + D3(18)O+ (deuteron transfer, reaction 3b)", p5t, cal5, 0, 0, F5_UNRESOLVED,
              "Cross sections corrected by the source for a 20.8% collection efficiency. Isotopic analogue; shows the "
              "isotope effect relative to 3a above ~4 eV."),
        table("ion-lishawa1990-fig7-charge-exchange-4a-sigma-digitized", "p. 3201, Fig. 7 (open circles)", 8,
              "H2O+ + D2(18)O -> H2O + D2(18)O+ (symmetric charge exchange, reaction 4a)", p7c, cal7, 0, 0, F7_UNRESOLVED,
              "Cross sections corrected by the source for a 20.8% collection efficiency. The source finds no isotope "
              "effect between 4a and 4b and states that above 1 eV charge exchange is the dominant channel with a "
              "cross section between 6 and 10e-16 cm2 (conclusion, p. 3205)."),
        table("ion-lishawa1990-fig7-charge-exchange-4b-sigma-digitized", "p. 3201, Fig. 7 (filled diamonds)", 8,
              "D2(18)O+ + H2O -> D2(18)O + H2O+ (symmetric charge exchange, reaction 4b)", p7d, cal7, 0, 0, F7_UNRESOLVED,
              "Cross sections corrected by the source for a 20.8% collection efficiency. Only isolated diamonds were "
              "digitized; most diamonds overlap other markers."),
    ]
    for t in tabs:
        with open(os.path.join(out, t["table_id"] + ".json"), "w") as f:
            json.dump(t, f, indent=1)
        print("wrote", t["table_id"], list(zip(t["x"], [round(v * 1e20, 2) for v in t["y"]])))


if __name__ == "__main__":
    main()
