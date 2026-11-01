#!/usr/bin/env python3
"""Digitization of the measured dissociative-recombination cross sections in two retained figures.

usage: dr_digitize.py rosen|jensen IMAGE.png OUT_DIR
  rosen   Rosen et al., Faraday Discuss. 115, 295 (2000), Fig. 4: H2O+ (solid squares, CRYRING), log-log,
          x 1e-4..1e2 eV, y (<v sigma>/v) 1e-18..1e-11 cm^2; native 300 ppi image extracted with pdfimages.
  jensen  Jensen et al., ApJ 543, 764 (2000), Fig. 3: H3O+ "this work (x10)" (open circles), log-log,
          x 1e-4..~50 eV, y 1e-20..~1e-11 cm^2; native 800 ppi stencil extracted with pdfimages.

Method (both figures):
 1. Ink mask = minority pixel value. Plot frame = the longest horizontal and vertical ink lines.
 2. Axis calibration from the tick marks attached to the frame: tick positions (column/row centres of the
    ink runs leaving the frame line), majors = the longest class; decade spacing from consecutive majors;
    minor ticks assigned to log10(2..9) of their decade; least-squares line pixel = a + b log10(value).
    Two independent calibrations: (A) bottom and left frame ticks, (B) top and right frame ticks.
    The calibration residual RMS and the A/B difference enter the reading uncertainty.
 3. Markers. rosen: filled squares = components of the ink eroded by a (s-2) x (s-2) square (s = marker
    size measured on an isolated marker); open squares and triangles have no filled core and drop out.
    jensen: open circles = peaks of the normalized cross-correlation with a ring template of the measured
    radius (non-maximum suppression; peaks shared with the inverted-triangle and square series are
    removed by a shape test: ring coverage >= 0.8 and straight-edge test). Every detection is drawn on an
    overlay image for visual review; the review result is recorded in the output.
 4. Values: E = 10^((x - a_x)/b_x), sigma = 10^((y - a_y)/b_y) cm^2 (jensen: divided by 10, the stated
    x10 of the H3O+ series), converted to m^2. Reading uncertainty per point: marker-centre uncertainty
    (+-1 px) combined with the calibration residual and the A/B difference, in decades.
Writes OUT_DIR/<name>-points.json (raw pixel centres, both calibrations, values, uncertainties) and
OUT_DIR/<name>-overlay.png.
"""
import json
import math
import sys

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage, signal


def ink_mask(path):
    a = np.array(Image.open(path).convert('L')) < 128
    if a.mean() > 0.5:
        a = ~a
    return a


def skew_angle(ink):
    """Rotation of the plot frame (radians, positive = bottom line descends to the right), from the centre row of
    the thick bottom frame line measured per column and fitted by least squares."""
    h, w = ink.shape
    rows = ink.sum(axis=1)
    rb = int(np.argmax(rows[h // 2:])) + h // 2
    xs, ys = [], []
    for x in range(int(0.2 * w), int(0.95 * w), 7):
        col = ink[rb - 40:rb + 40, x]
        idx = np.nonzero(col)[0]
        if idx.size == 0:
            continue
        runs = np.split(idx, np.nonzero(np.diff(idx) > 1)[0] + 1)
        r = max(runs, key=len)
        if 3 <= len(r) <= 12:
            xs.append(x)
            ys.append(rb - 40 + r.mean())
    xs, ys = np.array(xs), np.array(ys)
    k, c = np.polyfit(xs, ys, 1)
    resid = ys - (k * xs + c)
    keep = np.abs(resid) < 2.5
    k, c = np.polyfit(xs[keep], ys[keep], 1)
    return math.atan(k), int(keep.sum())


def frame(ink):
    h, w = ink.shape
    rows = [(r, longest_run(ink[r, :])) for r in range(h)]
    cols = [(c, longest_run(ink[:, c])) for c in range(w)]
    rmax = max(v for _, v in rows)
    cmax = max(v for _, v in cols)
    hr = [r for r, v in rows if v > 0.8 * rmax]
    vc = [c for c, v in cols if v > 0.8 * cmax]
    top, bottom = group_edges(hr)
    left, right = group_edges(vc)
    return top, bottom, left, right


def longest_run(v):
    best = cur = 0
    for x in v:
        cur = cur + 1 if x else 0
        best = max(best, cur)
    return best


def group_edges(idx):
    """first and last group of consecutive indices: return (inner edge of first group, inner edge of last group)."""
    groups = []
    for i in idx:
        if groups and i - groups[-1][-1] <= 1:
            groups[-1].append(i)
        else:
            groups.append([i])
    return groups[0][-1], groups[-1][0]


def ticks_along(ink, line, lo, hi, direction, horizontal, maxlen):
    """Tick runs leaving a frame line. horizontal=True: frame is a row 'line', ticks are vertical runs
    in columns lo..hi going in 'direction' (-1 up, +1 down). Positions are classified by run length first
    (major = long class, minor = short class) and grouped per class; minor groups wider than 3 px (merged
    neighbours near the top of a decade) are discarded. Returns [(centre, length)]."""
    lengths = []
    for k in range(lo, hi + 1):
        n = 0
        gap = 0
        for d in range(1, maxlen):
            p = line + direction * d
            v = ink[p, k] if horizontal else ink[k, p]
            if v:
                n = d
                gap = 0
            else:
                gap += 1
                if gap > 2:
                    break
        if n >= maxlen - 3:
            n = 0   # reaches the search window: a data feature (error bar, frame), not a tick
        lengths.append(n)
    Lmax = max(lengths)
    out = []
    for cls in ('major', 'minor'):
        cur = []
        groups = []
        for i, n in enumerate(lengths):
            ok = (n >= 0.75 * Lmax) if cls == 'major' else (3 <= n < 0.75 * Lmax)
            if ok:
                cur.append((lo + i, n))
            elif cur:
                groups.append(cur)
                cur = []
        if cur:
            groups.append(cur)
        wmax_minor = max(3, int(round(0.0016 * max(ink.shape))))   # 3 px at 976 px, 6 px at 3667 px
        for g in groups:
            if cls == 'minor' and len(g) > wmax_minor:
                continue
            if cls == 'major' and len(g) > 3 * wmax_minor:
                continue
            out.append((sum(p for p, _ in g) / len(g), max(n for _, n in g)))
    return sorted(out)


TICK_CLASSES = {   # (major length range, minor length range) in px, measured on each extracted image
    'rosen': ((12, 16), (5, 9)),
    'jensen': ((24, 32), (13, 22)),
}


def calibrate(ticks, first_major_value, reverse, edge=None, classes=None):
    """ticks [(pixel, length)]; returns (a, b, rms_px, details). value increases with pixel unless reverse.
    Majors/minors are classified relative to the median tick length (data features such as error bars that
    touch the frame have other lengths and drop out); majors may be every one or two decades; minors are
    assigned to log10(k), k = 1..9, of their decade; outliers (> 3 px) are rejected and the line refitted.
    edge = (pixel, value) of a labelled frame edge anchors the decade of the first major."""
    ticks = sorted(ticks)
    (M0, M1), (m0_, m1_) = classes
    Lm = (m0_ + m1_) / 2
    majors = [p for p, L in ticks if M0 <= L <= M1]
    minors = [p for p, L in ticks if m0_ <= L <= m1_]
    # majors are equally spaced: drop the member of any too-close pair that breaks the regular spacing
    while len(majors) > 2:
        diffs = [b - a for a, b in zip(majors, majors[1:])]
        med = float(np.median(diffs))
        bad = [i for i, d in enumerate(diffs) if d < 0.7 * med]
        if not bad:
            break
        i = bad[0]
        cand = []
        for drop in (i, i + 1):
            rest = majors[:drop] + majors[drop + 1:]
            dd = [b - a for a, b in zip(rest, rest[1:])]
            cand.append((float(np.std(dd)), drop))
        majors.pop(min(cand)[1])
    if reverse:
        majors = majors[::-1]
    gaps = [abs(b - a) for a, b in zip(majors, majors[1:])]
    gap = float(np.median(gaps))
    m0, m1 = majors[0], majors[1]
    between = [p for p in minors if min(m0, m1) < p < max(m0, m1)]
    dec_per_major = 1 if len(between) <= 10 else 2
    ppd = gap / dec_per_major
    sign = -1 if reverse else 1
    d0 = math.log10(first_major_value)
    if edge is not None:
        d0 = math.log10(edge[1]) + round(abs(majors[0] - edge[0]) / ppd)
    pts = [(p, d0 + i * dec_per_major) for i, p in enumerate(majors)]
    for p in minors:
        best = None
        for i, mp in enumerate(majors):
            f = sign * (p - mp) / ppd
            if f <= 0 or f >= dec_per_major:
                continue
            idec = math.floor(f + 0.02)
            rem = f - idec
            k = max(1, round(10 ** rem))
            if k <= 9 and abs(math.log10(k) + idec - f) < 0.035:
                best = d0 + i * dec_per_major + idec + math.log10(k)
            break
        if best is None and majors:
            f = sign * (majors[0] - p) / ppd   # before the first major (towards smaller values)
            if 0 < f < 1:
                k = max(1, round(10 ** (1 - f)))
                if k <= 9 and abs(math.log10(k) - (1 - f)) < 0.035:
                    best = d0 - 1 + math.log10(k)
        if best is not None:
            pts.append((p, best))
    P = np.array([p for p, _ in pts], dtype=float)
    X = np.array([x for _, x in pts], dtype=float)
    keep = np.ones(len(P), dtype=bool)
    for _ in range(3):
        A = np.vstack([np.ones(keep.sum()), X[keep]]).T
        (a, b), *_ = np.linalg.lstsq(A, P[keep], rcond=None)
        res = P - (a + b * X)
        new = np.abs(res) <= 3.0
        if (new == keep).all():
            break
        keep = new
    rms = float(np.sqrt(np.mean(res[keep] ** 2)))
    return float(a), float(b), rms, {'majors': len(majors), 'minors_used': int(keep.sum()) - len(majors), 'rejected': int((~keep).sum()),
                                    'pixels_per_decade': ppd, 'decades_per_major': dec_per_major, 'median_minor_length_px': Lm}


def calibrate_snapped(ticks, ref, classes):
    """Independent fit of one axis side whose majors are partly hidden: every tick of either length class is
    assigned the nearest valid tick value (decade + log10 k, k = 1..9) predicted by the reference line of the
    opposite side (accepted within 0.02 decades), then fitted by least squares with outlier rejection."""
    a_ref, b_ref = ref
    (M0, M1), (m0_, m1_) = classes
    P, X = [], []
    for p, L in ticks:
        if not (M0 <= L <= M1 or m0_ <= L <= m1_):
            continue
        x = (p - a_ref) / b_ref
        dec = math.floor(x)
        best = min((dec + math.log10(k) for k in range(1, 10)), key=lambda v: abs(v - x))
        best = min(best, dec + 1.0, key=lambda v: abs(v - x))
        if abs(best - x) < 0.02:
            P.append(p)
            X.append(best)
    P, X = np.array(P, dtype=float), np.array(X, dtype=float)
    keep = np.ones(len(P), dtype=bool)
    for _ in range(3):
        A = np.vstack([np.ones(keep.sum()), X[keep]]).T
        (a, b), *_ = np.linalg.lstsq(A, P[keep], rcond=None)
        res = P - (a + b * X)
        new = np.abs(res) <= 3.0
        if (new == keep).all():
            break
        keep = new
    rms = float(np.sqrt(np.mean(res[keep] ** 2)))
    return float(a), float(b), rms, {'ticks_used': int(keep.sum()), 'rejected': int((~keep).sum()), 'assignment': 'snapped to the opposite side'}


def main(which, path, out):
    ink = ink_mask(path)
    theta, nfit = skew_angle(ink)
    # deskew only when the frame is visibly tilted: the longest horizontal ink run (the bottom frame
    # line) is broken into pieces shorter than 90 % of the longest vertical-frame-based width estimate
    rows_run = max(longest_run(ink[r, :]) for r in range(ink.shape[0]))
    cols_span = ink.shape[1]
    deskew = rows_run < 0.75 * cols_span
    if deskew:
        # rotate the ink image about its centre so that the frame is axis-aligned (nearest neighbour)
        img = Image.fromarray((ink * 255).astype(np.uint8))
        img = img.rotate(math.degrees(theta), resample=Image.NEAREST, fillcolor=0)
        ink = np.array(img) > 127
    h, w = ink.shape
    top, bottom, left, right = frame(ink)
    maxlen = max(h, w) // 40
    tb = ticks_along(ink, bottom, left + 3, right - 3, -1, True, maxlen)
    tt = ticks_along(ink, top, left + 3, right - 3, +1, True, maxlen)
    tl = ticks_along(ink, left, top + 3, bottom - 3, +1, False, maxlen)
    tr = ticks_along(ink, right, top + 3, bottom - 3, -1, False, maxlen)
    # axis anchors read from the tick labels: both figures have 0.0001 eV at the left frame edge;
    # Rosen: the lowest y major is 1e-18 (bottom frame below it is not a decade); Jensen: 1e-20 at the bottom edge
    xedge = (left, 1e-4)
    if which == 'rosen':
        ykw = {'first_major_value': 1e-18}
    else:
        ykw = {'first_major_value': 1e-20, 'edge': (bottom, 1e-20)}
    cal = {}
    for name, (tx, ty) in {'A_bottom_left': (tb, tl), 'B_top_right': (tt, tr)}.items():
        try:
            ax, bx, rx, dx = calibrate(tx, 1e-4, False, edge=xedge, classes=TICK_CLASSES[which])
        except (IndexError, ValueError):
            ref = cal['A_bottom_left']['x']
            ax, bx, rx, dx = calibrate_snapped(tx, (ref[0], ref[1]), TICK_CLASSES[which])
        try:
            ay, by, ry, dy = calibrate(ty, ykw['first_major_value'], True, edge=ykw.get('edge'), classes=TICK_CLASSES[which])
        except (IndexError, ValueError):
            ref = cal['A_bottom_left']['y']
            ay, by, ry, dy = calibrate_snapped(ty, (ref[0], ref[1]), TICK_CLASSES[which])
        cal[name] = {'x': [ax, bx, rx, dx], 'y': [ay, by, ry, dy]}
    # markers
    det = []
    if which == 'rosen':
        lab, n = ndimage.label(ink)
        # marker size from components that are filled squares (fill ratio ~1, size 8..20)
        sizes = []
        for sl in ndimage.find_objects(lab):
            hh, ww = sl[0].stop - sl[0].start, sl[1].stop - sl[1].start
            sub = ink[sl]
            if 6 <= hh <= 24 and abs(hh - ww) <= 1 and sub.mean() > 0.93:
                sizes.append(hh)
        s = int(np.median(sizes))
        # filled markers only: opening with an (s-2)^2 square removes open squares, triangles, error bars
        # and text strokes; overlapping filled squares stay as merged bands
        opened = ndimage.binary_opening(ink, structure=np.ones((s - 2, s - 2)))
        opened[:, :left + s] = False
        opened[:, right - s:] = False
        opened[:top + s, :] = False
        opened[bottom - s:, :] = False
        # band trace: per column the centres of the vertical runs, binned at one marker width
        cols = {}
        for x in range(w):
            ys = np.nonzero(opened[:, x])[0]
            if ys.size == 0:
                continue
            runs = np.split(ys, np.nonzero(np.diff(ys) > 1)[0] + 1)
            cols[x] = [float(r.mean()) for r in runs]
        xs = sorted(cols)
        bins = []
        cur = []
        for x in xs:
            if cur and x - cur[0] >= s:
                bins.append(cur)
                cur = []
            cur.append(x)
        if cur:
            bins.append(cur)
        for bcols in bins:
            ys = [y for x in bcols for y in cols[x]]
            ys.sort()
            ymed = ys[len(ys) // 2]
            spread = (ys[-1] - ys[0]) / 2.0
            det.append((float(np.mean(bcols)), float(ymed), float(spread), len(bcols)))
        meta = {'marker_size_px': s, 'opening_kernel': s - 2, 'method': 'filled-marker band trace, bins of one marker width'}
    else:
        # empirical templates cut from isolated markers of each series in this figure (bounding boxes
        # measured on the extracted image): H3O+ this work (circle with its error bar), D3O+ (square),
        # previous ASTRID H3O+ (inverted triangle), HD2O+ (triangle), Mul et al. (filled circle)
        boxes = {'circle': (98, 688, 30, 30), 'circle_b': (221, 1081, 30, 30), 'square': (392, 921, 29, 29),
                 'nabla': (139, 822, 27, 27), 'triangle': (637, 985, 27, 27), 'filled': (457, 1846, 36, 35)}
        W = 41
        f = ink.astype(float)
        box = np.ones((W, W))
        energy = signal.fftconvolve(f, box, mode='same')
        scores = {}
        for name, (r0, c0, hh, ww) in boxes.items():
            t = np.zeros((W, W))
            sub = ink[r0:r0 + hh, c0:c0 + ww].astype(float)
            t[(W - hh) // 2:(W - hh) // 2 + hh, (W - ww) // 2:(W - ww) // 2 + ww] = sub
            num = signal.fftconvolve(f, t[::-1, ::-1], mode='same')
            scores[name] = num / np.sqrt(np.maximum(energy, 1e-9) * (t ** 2).sum())
        circ = np.maximum(scores['circle'], scores['circle_b'])
        other = np.maximum.reduce([ndimage.maximum_filter(scores[k], size=7) for k in ('nabla', 'filled')])
        mx = ndimage.maximum_filter(circ, size=21)
        peaks = np.argwhere((circ == mx) & (circ > 0.60) & (circ > other))
        # the H3O+ (x10) series band: polyline through the series read on the figure (deskewed pixels),
        # +-BAND px; the D3O+ and HD2O+ series lie outside it for E >= 0.1 eV
        path_disp = [(1000, 250), (1100, 320), (1175, 340), (1250, 370), (1300, 430), (1345, 465), (1370, 530), (1385, 590),
                     (1410, 630), (1440, 650), (1465, 725), (1490, 760), (1520, 785), (1545, 800), (1570, 850), (1600, 850),
                     (1620, 800), (1640, 785), (1680, 720), (1710, 700), (1740, 680), (1780, 660), (1820, 665), (1860, 700),
                     (1880, 750)]
        bpath = [(x * 1.8335, y * 1.8335) for x, y in path_disp]
        BAND = 110
        def band_dist(px, py):
            best = 1e9
            for (x0, y0), (x1, y1) in zip(bpath, bpath[1:]):
                vx, vy = x1 - x0, y1 - y0
                t = max(0.0, min(1.0, ((px - x0) * vx + (py - y0) * vy) / (vx * vx + vy * vy)))
                best = min(best, math.hypot(px - x0 - t * vx, py - y0 - t * vy))
            return best
        for py, px in peaks:
            if px < bpath[0][0] or px > bpath[-1][0]:
                continue
            if band_dist(px, py) <= BAND:
                det.append((float(px), float(py), 0.0, 1))
        meta = {'method': 'empirical-template cosine similarity inside the H3O+ series band (+-%d px), circle score > 0.60 and above the inverted-triangle and filled-circle scores, NMS 21 px; reviewed on overlays, edits listed in manual_edits' % BAND,
                'templates': boxes, 'band_path_display_px': path_disp, 'band_scale': 1.8335}
    # values
    pts = []
    for px, py, spread, ncol in sorted(det):
        rec = {'px': px, 'py': py, 'band_half_spread_px': spread, 'columns': ncol}
        vals = {}
        for name, c in cal.items():
            ax, bx = c['x'][0], c['x'][1]
            ay, by = c['y'][0], c['y'][1]
            lx = (px - ax) / bx
            ly = (py - ay) / by
            vals[name] = (lx, ly)
        lxA, lyA = vals['A_bottom_left']
        lxB, lyB = vals['B_top_right']
        ux = math.sqrt((1 / abs(cal['A_bottom_left']['x'][1])) ** 2 * (1 + cal['A_bottom_left']['x'][2] ** 2) + (lxA - lxB) ** 2)
        uy = math.sqrt((1 / abs(cal['A_bottom_left']['y'][1])) ** 2 * (1 + cal['A_bottom_left']['y'][2] ** 2 + spread ** 2) + (lyA - lyB) ** 2)
        scale = 0.1 if which == 'jensen' else 1.0
        rec.update({'log10_E_A': lxA, 'log10_E_B': lxB, 'log10_sigma_cm2_A': lyA, 'log10_sigma_cm2_B': lyB,
                    'E_eV': 10 ** lxA, 'sigma_m2': 10 ** lyA * 1e-4 * scale, 'u_log10_E': ux, 'u_log10_sigma': uy})
        pts.append(rec)
    json.dump({'figure': which, 'image': path, 'image_size': [w, h], 'skew_rad': theta, 'skew_fit_columns': nfit, 'deskewed': deskew,
               'frame': {'top': top, 'bottom': bottom, 'left': left, 'right': right},
               'calibration': cal, 'detection': meta, 'points': pts}, open('%s/%s-points.json' % (out, which), 'w'), indent=1)
    im = Image.fromarray(((~ink) * 255).astype(np.uint8)).convert('RGB')   # ink dark on white, deskewed frame
    dr = ImageDraw.Draw(im)
    for k, p in enumerate(pts):
        x, y = p['px'], p['py']
        dr.ellipse([x - 4, y - 4, x + 4, y + 4], outline=(255, 0, 0), width=2)
        dr.text((x + 6, y - 14), str(k), fill=(0, 0, 255))
    im.save('%s/%s-overlay.png' % (out, which))
    print(which, 'skew %.5f rad (deskewed %s)' % (theta, deskew), 'frame', top, bottom, left, right, 'markers', len(pts))
    for name, c in cal.items():
        print('  calibration', name, 'x: b=%.3f px/decade rms %.2f px %s' % (c['x'][1], c['x'][2], c['x'][3]),
              'y: b=%.3f rms %.2f %s' % (c['y'][1], c['y'][2], c['y'][3]))


if __name__ == '__main__':
    main(*sys.argv[1:4])
