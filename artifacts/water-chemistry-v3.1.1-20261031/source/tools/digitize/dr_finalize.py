#!/usr/bin/env python3
"""Final digitized DR cross sections (Rosen 2000 Fig. 4 H2O+, Jensen 2000 Fig. 3 H3O+) with reading
uncertainties; writes data/extracted/rec-dr-*-sigma-digitized.json of the v3.1 component.

usage: dr_finalize.py DIGITIZATION_DIR
Inputs: <dir>/rosen-points.json and <dir>/jensen-points.json (dr_digitize.py), the extracted figure images,
and the manual review edits below (every edit is listed in the output).

Repeat digitization (second, independent pass):
  Rosen: the band trace is repeated with an (s-3)^2 opening kernel and bins offset by half a marker width;
         the pass-1/pass-2 difference interpolated at the pass-1 energies is the repeatability.
  Jensen: each marker centre is re-measured as the centroid of the marker's own ink within its 30x30 box
          (pass 2) versus the template-correlation peak or manual reading (pass 1); the band trace above 6 eV
          is repeated with bins offset by half a bin.
Per-point uncertainty (decades, 1 sigma, combined in quadrature): calibration (fit RMS and the difference
between the bottom/left and top/right calibrations), +-1 px marker centre, band half-spread (band traces),
and the pass-1/pass-2 difference.
"""
import json
import math
import os
import sys

import numpy as np
from PIL import Image
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dr_digitize as dd  # noqa: E402

COMP = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

JENSEN_EDITS = {
    'remove_px': [[2376, 778, 'inverted triangle (previous ASTRID)'], [2705, 1134, 'inverted triangle (previous ASTRID)'],
                  [2747, 1414, 'duplicate of the circle at (2748, 1394)'], [2748, 1395, 'duplicate of the circle at (2748, 1394)'],
                  [3009, 1238, 'inverted triangle (previous ASTRID)'], [3023, 1299, 'error-bar cap of the circle at (3024, 1316)'],
                  [3041, 1404, 'error-bar cap of the circle at (3042, 1362)'], [3085, 1251, 'inverted triangle (previous ASTRID)'],
                  [2850, 1590, 'error-bar cap']],
    'add_px': [[2325, 718, 2], [2595, 1104, 2], [2781, 1384, 2], [2920, 1514, 2], [2937, 1492, 2],
               [3062, 1313, 4], [3083, 1318, 4], [3104, 1318, 4]],   # [x, y, reading uncertainty px]
    'markers_from_eV': 0.15, 'band_from_px': 3150, 'band_to_px': 3440, 'band_rows': [1150, 1420],
    'band_note': 'above 6 eV the H3O+ circles (this work) and the previous ASTRID inverted triangles overlap; '
                 'the source states agreement with the previous ASTRID data outside the 0.6-5 eV valley, so the '
                 'upper band (both series, both x10) is traced together there',
}


def load_ink(path, deskew):
    ink = dd.ink_mask(path)
    th = 0.0
    if deskew:
        th, _ = dd.skew_angle(ink)
        img = Image.fromarray((ink * 255).astype(np.uint8)).rotate(math.degrees(th), resample=Image.NEAREST, fillcolor=0)
        ink = np.array(img) > 127
    return ink, th


def to_values(px, py, cal, scale):
    out = {}
    for name, c in cal.items():
        ax, bx = c['x'][0], c['x'][1]
        ay, by = c['y'][0], c['y'][1]
        out[name] = ((px - ax) / bx, (py - ay) / by)
    lxA, lyA = out['A_bottom_left']
    lxB, lyB = out['B_top_right']
    return lxA, lyA + math.log10(scale), lxA - lxB, lyA - lyB


def band_trace(opened, x0, x1, binw, offset, rows=None):
    pts = []
    xs = [x for x in range(x0, x1) if opened[:, x].any()] if rows is None else list(range(x0, x1))
    start = x0 + offset
    edges = list(range(start, x1 + binw, binw))
    for a, b in zip(edges, edges[1:]):
        ys = []
        cols = []
        for x in range(max(a, x0), min(b, x1)):
            col = opened[:, x] if rows is None else opened[rows[0]:rows[1], x]
            idx = np.nonzero(col)[0]
            if idx.size == 0:
                continue
            if rows is not None:
                idx = idx + rows[0]
            ys.extend(idx.tolist())
            cols.append(x)
        if len(cols) < max(2, binw // 3):
            continue
        ys.sort()
        pts.append((float(np.mean(cols)), float(ys[len(ys) // 2]), (ys[int(0.84 * (len(ys) - 1))] - ys[int(0.16 * (len(ys) - 1))]) / 2.0))
    return pts


def interp_log(xs, ys, x):
    for (a, ya), (b, yb) in zip(zip(xs, ys), zip(xs[1:], ys[1:])):
        if a <= x <= b:
            return ya + (yb - ya) * (x - a) / (b - a)
    return None


def finalize_rosen(d):
    R = json.load(open(os.path.join(d, 'rosen-points.json')))
    cal = R['calibration']
    ink, _ = load_ink(R['image'], False)
    s = R['detection']['marker_size_px']
    fr = R['frame']
    pts1 = []
    for p in R['points']:
        lx, ly, dxab, dyab = to_values(p['px'], p['py'], cal, 1.0)
        pts1.append({'px': p['px'], 'py': p['py'], 'lx': lx, 'ly': ly, 'spread_px': p['band_half_spread_px'], 'dxab': dxab, 'dyab': dyab})
    # pass 2: smaller kernel, bins offset by half a marker width
    op2 = ndimage.binary_opening(ink, structure=np.ones((s - 3, s - 3)))
    op2[:, :fr['left'] + s] = False
    op2[:, fr['right'] - s:] = False
    op2[:fr['top'] + s, :] = False
    op2[fr['bottom'] - s:, :] = False
    x0 = int(min(p['px'] for p in R['points'])) - s
    x1 = int(max(p['px'] for p in R['points'])) + s
    t2 = band_trace(op2, x0, x1, s, s // 2)
    l2 = [to_values(px, py, cal, 1.0) for px, py, _ in t2]
    xs2 = [v[0] for v in l2]
    ys2 = [v[1] for v in l2]
    by = abs(cal['A_bottom_left']['y'][1])
    bx = abs(cal['A_bottom_left']['x'][1])
    ry = cal['A_bottom_left']['y'][2]
    rx = cal['A_bottom_left']['x'][2]
    out = []
    own = 0.34 * s   # 16-84 % half-range of the rows of a single marker of size s
    for p in pts1:
        rep = interp_log(xs2, ys2, p['lx'])
        drep = (p['ly'] - rep) if rep is not None else None
        scatter = math.sqrt(max(0.0, p['spread_px'] ** 2 - own ** 2))
        uy = math.sqrt((1 + ry ** 2 + scatter ** 2) / by ** 2 + p['dyab'] ** 2 + (drep or 0.0) ** 2)
        ux = math.sqrt((1 + rx ** 2) / bx ** 2 + p['dxab'] ** 2)
        out.append(dict(p, repeat_diff_log10=drep, u_log10_sigma=uy, u_log10_E=ux))
    return out, {'pass2_points': len(t2), 'kernel_pass1': s - 2, 'kernel_pass2': s - 3, 'bin_px': s, 'pass2_bin_offset_px': s // 2}


def finalize_jensen(d):
    J = json.load(open(os.path.join(d, 'jensen-points.json')))
    cal = J['calibration']
    ink, th = load_ink(J['image'], True)
    E = JENSEN_EDITS
    rem = E['remove_px']
    pts = []
    for p in J['points']:
        if any(abs(p['px'] - x) <= 3 and abs(p['py'] - y) <= 3 for x, y, _ in rem):
            continue
        pts.append((p['px'], p['py'], 1.0, 'template'))
    for x, y, u in E['add_px']:
        pts.append((float(x), float(y), float(u), 'manual'))
    by = abs(cal['A_bottom_left']['y'][1])
    bx = abs(cal['A_bottom_left']['x'][1])
    ry, rx = cal['A_bottom_left']['y'][2], cal['A_bottom_left']['x'][2]
    out = []
    for px, py, u, how in sorted(pts):
        if px >= E['band_from_px']:
            continue
        lx, ly, dxab, dyab = to_values(px, py, cal, 0.1)
        if 10 ** lx < E['markers_from_eV']:
            continue
        # pass 2: centroid of the marker ink in its 30x30 box
        y0, x0 = int(py) - 15, int(px) - 15
        sub = ink[y0:y0 + 31, x0:x0 + 31]
        cy, cx = ndimage.center_of_mass(sub)
        cx, cy = cx + x0, cy + y0
        lx2, ly2, _, _ = to_values(cx, cy, cal, 0.1)
        drep = ly - ly2
        uy = math.sqrt((u ** 2 + ry ** 2) / by ** 2 + dyab ** 2 + drep ** 2)
        ux = math.sqrt((u ** 2 + rx ** 2) / bx ** 2 + dxab ** 2 + (lx - lx2) ** 2)
        out.append({'px': px, 'py': py, 'how': how, 'lx': lx, 'ly': ly, 'pass2_px': cx, 'pass2_py': cy, 'repeat_diff_log10': drep,
                    'u_log10_sigma': uy, 'u_log10_E': ux})
    # band trace above 6 eV (upper band only)
    rows = E['band_rows']
    binw = 18
    t1 = band_trace(ink, E['band_from_px'], E['band_to_px'], binw, 0, rows=rows)
    t2 = band_trace(ink, E['band_from_px'], E['band_to_px'], binw, binw // 2, rows=rows)
    l2 = [to_values(px, py, cal, 0.1) for px, py, _ in t2]
    xs2, ys2 = [v[0] for v in l2], [v[1] for v in l2]
    for px, py, spread in t1:
        lx, ly, dxab, dyab = to_values(px, py, cal, 0.1)
        rep = interp_log(xs2, ys2, lx)
        drep = (ly - rep) if rep is not None else 0.0
        scatter = math.sqrt(max(0.0, spread ** 2 - (0.34 * 30) ** 2))   # minus the extent of one 30 px marker
        uy = math.sqrt((1 + ry ** 2 + scatter ** 2) / by ** 2 + dyab ** 2 + drep ** 2)
        out.append({'px': px, 'py': py, 'how': 'band', 'lx': lx, 'ly': ly, 'band_half_spread_px': spread, 'repeat_diff_log10': drep,
                    'u_log10_sigma': uy, 'u_log10_E': math.sqrt((1 + rx ** 2) / bx ** 2 + dxab ** 2)})
    out.sort(key=lambda r: r['lx'])
    return out, {'skew_rad': th, 'band_bin_px': binw}


def write(table_id, source_id, locator, process, pts, meta, notes, cutoff_lo):
    sel = [p for p in pts if 10 ** p['lx'] >= cutoff_lo]
    # the tabulated range ends before the first bin whose reading uncertainty exceeds a factor of 2
    dropped = []
    for i, p in enumerate(sel):
        if p['u_log10_sigma'] > math.log10(2.0):
            dropped = sel[i:]
            sel = sel[:i]
            break
    meta = dict(meta, dropped_tail=[{'E_eV': 10 ** p['lx'], 'sigma_m2': 10 ** p['ly'] * 1e-4, 'u_log10_sigma': p['u_log10_sigma']} for p in dropped])
    # strictly increasing energies (merge points closer than 0.002 decades)
    merged = []
    for p in sel:
        if merged and p['lx'] - merged[-1]['lx'] < 0.002:
            q = merged[-1]
            q['ly'] = 0.5 * (q['ly'] + p['ly'])
            q['u_log10_sigma'] = max(q['u_log10_sigma'], p['u_log10_sigma'])
            continue
        merged.append(dict(p))
    rec = {
        'table_id': table_id, 'source_id': source_id, 'locator': locator, 'process': process, 'quantity': 'cross_section_m2',
        'x_name': 'energy_eV_cm', 'x': [10 ** p['lx'] for p in merged], 'y': [10 ** p['ly'] * 1e-4 for p in merged],
        'y_uncertainty_relative_1sigma': [10 ** p['u_log10_sigma'] - 1 for p in merged],
        'x_uncertainty_log10': [p['u_log10_E'] for p in merged],
        'original_units': 'cm^2 (figure ordinate)', 'conversion_to_SI': 'sigma[m2] = sigma[cm2] * 1e-4',
        'tag': 'measured', 'extraction_method': 'figure digitization (measured, figure-digitized): ' + meta['method'],
        'digitization': meta, 'points_detail': merged, 'notes': notes,
    }
    path = os.path.join(COMP, 'data/extracted', table_id + '.json')
    json.dump(rec, open(path, 'w'), indent=1)
    u = [r for r in rec['y_uncertainty_relative_1sigma']]
    print(table_id, len(rec['x']), 'points %.4g-%.4g eV' % (rec['x'][0], rec['x'][-1]), 'reading uncertainty median %.3f max %.3f' % (sorted(u)[len(u) // 2], max(u)))


def main(d):
    rp, rmeta = finalize_rosen(d)
    R = json.load(open(os.path.join(d, 'rosen-points.json')))
    write('rec-dr-rosen2000-h2oplus-sigma-digitized', 'rec-dr-rosen-2000-h2oplus-cryring',
          'Rosen et al., Faraday Discuss. 115, 295 (2000), Fig. 4 (printed p. 300), solid squares (CRYRING), (<v sigma>/v) in cm^2 vs centre-of-mass energy; '
          'native 300 ppi figure image extracted from the retained PDF (sha256 fb9893fd1e828c452a7ae672b3d360687ffbb7e65a6e5cbc3e2ff16191db0d7d)',
          'H2O+ + e -> neutral products (total DR)', rp,
          {'method': 'filled-marker band trace (opening with an (s-2)^2 square; open squares, triangles and error bars removed), median per marker-width bin; '
                     'log-log axes calibrated from the frame ticks (bottom/left and top/right independently)',
           'calibration': R['calibration'], 'frame': R['frame'], **rmeta},
          'Ordinate is <v sigma>/v as plotted by the source (merged-beam rate divided by the relative velocity, not fully deconvolved; the source states '
          '~20 % total uncertainty at 1 sigma and up to 30 % statistical above 20 eV). Values below 2 eV are not used by the tables (the stated '
          'power laws are kept there); they verify the join. Reading uncertainty excludes the source uncertainty.', 0.5)
    jp, jmeta = finalize_jensen(d)
    J = json.load(open(os.path.join(d, 'jensen-points.json')))
    write('rec-dr-jensen2000-h3oplus-sigma-digitized', 'rec-dr-jensen-2000-h3oplus-astrid',
          'Jensen et al., ApJ 543, 764 (2000), Fig. 3, H3O+ "this work (x10)" open circles (divided by 10 here), sigma in cm^2 vs E; '
          'native 800 ppi stencil extracted from the retained PDF (sha256 a32df2e0d3ac3aa327b7e7f0bc0170169ad217b7fea7d526e0028625a92c1ba0), deskewed',
          'H3O+ + e -> neutral products (total DR)', jp,
          {'method': 'empirical-template detection of the H3O+ circles inside the series band, reviewed on overlays with the listed manual edits; '
                     'above 6 eV upper-band trace (circles and previous ASTRID points, see edits note); log-log axes calibrated from the frame ticks',
           'calibration': J['calibration'], 'frame': J['frame'], 'manual_edits': JENSEN_EDITS, **jmeta},
          'Figure values are the H3O+ series multiplied by 10 by the source; divided by 10 here. The source states sigma ~ E^-1.15 at low energy, '
          'structure 0.2-1 eV, a valley 0.6-5 eV where the earlier ASTRID data are higher (background not subtracted), and a strong peak near '
          '10 eV; cross sections were measured from 0 to 20-35 eV; ~ the H3O+ series ends near 20 eV in this figure.', 0.15)


if __name__ == '__main__':
    main(sys.argv[1])
