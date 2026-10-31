#!/usr/bin/python -I
"""Plots of a Hall v3.2 pre-run time series (series.jsonl of test_v32 prerun) as white-background PNG files.
Standard library only (zlib PNG writer, 5x7 bitmap font); no plotting package is needed.

  plot_series.py SERIES_JSONL OUT_DIR [TITLE]
"""
import json, math, struct, sys, zlib
from pathlib import Path

FONT = {  # 5x7 glyphs, rows top to bottom, 5 bits each (MSB = left)
    'A': (0x0E, 0x11, 0x11, 0x1F, 0x11, 0x11, 0x11), 'B': (0x1E, 0x11, 0x11, 0x1E, 0x11, 0x11, 0x1E), 'C': (0x0E, 0x11, 0x10, 0x10, 0x10, 0x11, 0x0E),
    'D': (0x1E, 0x11, 0x11, 0x11, 0x11, 0x11, 0x1E), 'E': (0x1F, 0x10, 0x10, 0x1E, 0x10, 0x10, 0x1F), 'F': (0x1F, 0x10, 0x10, 0x1E, 0x10, 0x10, 0x10),
    'G': (0x0E, 0x11, 0x10, 0x17, 0x11, 0x11, 0x0F), 'H': (0x11, 0x11, 0x11, 0x1F, 0x11, 0x11, 0x11), 'I': (0x0E, 0x04, 0x04, 0x04, 0x04, 0x04, 0x0E),
    'J': (0x07, 0x02, 0x02, 0x02, 0x02, 0x12, 0x0C), 'K': (0x11, 0x12, 0x14, 0x18, 0x14, 0x12, 0x11), 'L': (0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x1F),
    'M': (0x11, 0x1B, 0x15, 0x15, 0x11, 0x11, 0x11), 'N': (0x11, 0x11, 0x19, 0x15, 0x13, 0x11, 0x11), 'O': (0x0E, 0x11, 0x11, 0x11, 0x11, 0x11, 0x0E),
    'P': (0x1E, 0x11, 0x11, 0x1E, 0x10, 0x10, 0x10), 'Q': (0x0E, 0x11, 0x11, 0x11, 0x15, 0x12, 0x0D), 'R': (0x1E, 0x11, 0x11, 0x1E, 0x14, 0x12, 0x11),
    'S': (0x0F, 0x10, 0x10, 0x0E, 0x01, 0x01, 0x1E), 'T': (0x1F, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04), 'U': (0x11, 0x11, 0x11, 0x11, 0x11, 0x11, 0x0E),
    'V': (0x11, 0x11, 0x11, 0x11, 0x11, 0x0A, 0x04), 'W': (0x11, 0x11, 0x11, 0x15, 0x15, 0x15, 0x0A), 'X': (0x11, 0x11, 0x0A, 0x04, 0x0A, 0x11, 0x11),
    'Y': (0x11, 0x11, 0x0A, 0x04, 0x04, 0x04, 0x04), 'Z': (0x1F, 0x01, 0x02, 0x04, 0x08, 0x10, 0x1F),
    '0': (0x0E, 0x11, 0x13, 0x15, 0x19, 0x11, 0x0E), '1': (0x04, 0x0C, 0x04, 0x04, 0x04, 0x04, 0x0E), '2': (0x0E, 0x11, 0x01, 0x02, 0x04, 0x08, 0x1F),
    '3': (0x1F, 0x02, 0x04, 0x02, 0x01, 0x11, 0x0E), '4': (0x02, 0x06, 0x0A, 0x12, 0x1F, 0x02, 0x02), '5': (0x1F, 0x10, 0x1E, 0x01, 0x01, 0x11, 0x0E),
    '6': (0x06, 0x08, 0x10, 0x1E, 0x11, 0x11, 0x0E), '7': (0x1F, 0x01, 0x02, 0x04, 0x08, 0x08, 0x08), '8': (0x0E, 0x11, 0x11, 0x0E, 0x11, 0x11, 0x0E),
    '9': (0x0E, 0x11, 0x11, 0x0F, 0x01, 0x02, 0x0C), '.': (0, 0, 0, 0, 0, 0x0C, 0x0C), '-': (0, 0, 0, 0x1F, 0, 0, 0), '+': (0, 0x04, 0x04, 0x1F, 0x04, 0x04, 0),
    '(': (0x02, 0x04, 0x08, 0x08, 0x08, 0x04, 0x02), ')': (0x08, 0x04, 0x02, 0x02, 0x02, 0x04, 0x08), '/': (0x01, 0x01, 0x02, 0x04, 0x08, 0x10, 0x10),
    '%': (0x18, 0x19, 0x02, 0x04, 0x08, 0x13, 0x03), '_': (0, 0, 0, 0, 0, 0, 0x1F), ',': (0, 0, 0, 0, 0x0C, 0x04, 0x08), ':': (0, 0x0C, 0x0C, 0, 0x0C, 0x0C, 0),
    '|': (0x04,) * 7, '=': (0, 0, 0x1F, 0, 0x1F, 0, 0), '^': (0x04, 0x0A, 0x11, 0, 0, 0, 0), ' ': (0,) * 7,
}
COLORS = [(31, 119, 180), (214, 39, 40), (44, 160, 44), (255, 127, 14), (148, 103, 189), (140, 86, 75)]


class Canvas:
    def __init__(self, w, h):
        self.w, self.h = w, h
        self.px = bytearray([255] * (w * h * 3))

    def set(self, x, y, c):
        if 0 <= x < self.w and 0 <= y < self.h:
            i = 3 * (y * self.w + x)
            self.px[i:i + 3] = bytes(c)

    def line(self, x0, y0, x1, y1, c, width=1):
        x0, y0, x1, y1 = int(round(x0)), int(round(y0)), int(round(x1)), int(round(y1))
        dx, dy = abs(x1 - x0), -abs(y1 - y0)
        sx, sy = (1 if x0 < x1 else -1), (1 if y0 < y1 else -1)
        err = dx + dy
        while True:
            for ox in range(-(width // 2), width - width // 2):
                for oy in range(-(width // 2), width - width // 2):
                    self.set(x0 + ox, y0 + oy, c)
            if x0 == x1 and y0 == y1:
                break
            e2 = 2 * err
            if e2 >= dy:
                err += dy; x0 += sx
            if e2 <= dx:
                err += dx; y0 += sy

    def text(self, x, y, s, c=(0, 0, 0), scale=2, vertical=False):
        for n, ch in enumerate(s):
            g = FONT.get(ch.upper(), FONT[' '])
            for r, bits in enumerate(g):
                for col in range(5):
                    if bits & (1 << (4 - col)):
                        for a in range(scale):
                            for b in range(scale):
                                if vertical:
                                    self.set(x + r * scale + a, y - (n * 6 + col) * scale - b, c)
                                else:
                                    self.set(x + (n * 6 + col) * scale + a, y + r * scale + b, c)

    def save(self, path):
        raw = b''.join(b'\x00' + bytes(self.px[3 * y * self.w:3 * (y + 1) * self.w]) for y in range(self.h))
        def chunk(t, d):
            return struct.pack('>I', len(d)) + t + d + struct.pack('>I', zlib.crc32(t + d) & 0xffffffff)
        Path(path).write_bytes(b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', self.w, self.h, 8, 2, 0, 0, 0)) +
                               chunk(b'IDAT', zlib.compress(raw, 9)) + chunk(b'IEND', b''))


def nice_ticks(lo, hi, n=6):
    if hi <= lo:
        hi = lo + 1
    step = (hi - lo) / n
    mag = 10 ** math.floor(math.log10(step))
    step = min((m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= step), default=10 * mag)
    t = math.ceil(lo / step) * step
    out = []
    while t <= hi + 1e-12 * step:
        out.append(t); t += step
    return out


def fmt(v):
    if v == 0:
        return '0'
    a = abs(v)
    if 1e-3 <= a < 1e5:
        s = ('%.4f' % v).rstrip('0').rstrip('.')
        return s
    return '%.2e' % v


def plot(path, title, xlabel, ylabel, series, logy=False, W=1100, H=650):
    c = Canvas(W, H)
    L, R, T, B = 130, 30, 60, 90
    xs = [x for _, pts, _ in series for x, _ in pts]
    ys = [y for _, pts, _ in series for _, y in pts if (y > 0 or not logy)]
    if not xs or not ys:
        return
    x0, x1 = min(xs), max(xs)
    if logy:
        y0, y1 = math.floor(math.log10(min(ys))), math.ceil(math.log10(max(ys)))
        if y1 == y0:
            y1 += 1
        yt = list(range(int(y0), int(y1) + 1))
        fy = lambda y: math.log10(y)
    else:
        y0, y1 = min(ys), max(ys)
        pad = 0.05 * (y1 - y0) if y1 > y0 else (abs(y0) * 0.05 + 1e-30)
        y0, y1 = y0 - pad, y1 + pad
        yt = nice_ticks(y0, y1)
        fy = lambda y: y
    X = lambda x: L + (x - x0) / ((x1 - x0) or 1) * (W - L - R)
    Y = lambda y: H - B - (fy(y) - y0) / ((y1 - y0) or 1) * (H - T - B)
    grid, axis = (225, 225, 225), (0, 0, 0)
    for t in nice_ticks(x0, x1, 8):
        px = X(t); c.line(px, T, px, H - B, grid); c.line(px, H - B, px, H - B + 6, axis)
        s = fmt(t); c.text(int(px - len(s) * 6), H - B + 12, s)
    for t in yt:
        py = Y(10 ** t) if logy else Y(t)
        c.line(L, py, W - R, py, grid); c.line(L - 6, py, L, py, axis)
        s = ('1E%d' % t) if logy else fmt(t); c.text(L - 10 - len(s) * 12, int(py - 7), s)
    c.line(L, T, L, H - B, axis); c.line(L, H - B, W - R, H - B, axis); c.line(L, T, W - R, T, axis); c.line(W - R, T, W - R, H - B, axis)
    for k, (name, pts, color) in enumerate(series):
        prev = None
        for x, y in pts:
            if logy and y <= 0:
                prev = None; continue
            p = (X(x), Y(y))
            if prev:
                c.line(prev[0], prev[1], p[0], p[1], color, 2)
            prev = p
        c.line(L + 20, T + 18 + 22 * k, L + 60, T + 18 + 22 * k, color, 3)
        c.text(L + 70, T + 11 + 22 * k, name, color)
    c.text(L, 18, title, scale=3)
    c.text(int((L + W - R) / 2 - len(xlabel) * 6), H - 40, xlabel)
    c.text(22, int((T + H - B) / 2 + len(ylabel) * 6), ylabel, vertical=True)
    c.save(path)


def main(series_path, out_dir, title='Hall v3.2 pre-run'):
    rows = [json.loads(l) for l in open(series_path)]
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    t = [r['time_s'] * 1e9 for r in rows]
    S = lambda f: [(x, f(r)) for x, r in zip(t, rows)]
    files = []
    def P(name, *a, **k):
        plot(out / name, *a, **k); files.append(str(out / name))
    P('terminal-voltage.png', title + ': terminal voltage', 'time (ns)', 'voltage (kV)', [('terminal V', S(lambda r: r['terminal_V'] / 1e3), COLORS[0])])
    P('drive-current.png', title + ': drive current (C19)', 'time (ns)', 'current (kA)', [('I(t)', S(lambda r: r['drive_current_at_t_A'] / 1e3), COLORS[0])])
    P('field-energy.png', title + ': field energy', 'time (ns)', 'energy (J)', [('W_E', S(lambda r: r['W_E_J']), COLORS[0]), ('W_B', S(lambda r: r['W_B_J']), COLORS[1])])
    faces = ['lower_inner', 'lower_outer', 'upper_inner', 'upper_outer']
    P('surface-field.png', title + ': max surface normal field', 'time (ns)', 'field (MV/m)',
      [(f, S(lambda r, f=f: r['surface_field'][f]['max_abs_E_normal_V_m'] / 1e6), COLORS[i]) for i, f in enumerate(faces)])
    P('foil-temperature.png', title + ': foil temperatures', 'time (ns)', 'temperature (K)',
      [('lower max', S(lambda r: r['lower_foil']['max_T_K']), COLORS[0]), ('lower mean', S(lambda r: r['lower_foil']['mean_T_K']), COLORS[2]),
       ('upper max', S(lambda r: r['upper_foil']['max_T_K']), COLORS[1]), ('upper mean', S(lambda r: r['upper_foil']['mean_T_K']), COLORS[3])])
    P('gauss-divb-ratios.png', title + ': guard ratios', 'time (ns)', 'ratio', [
        ('gauss |D-L|/bound', S(lambda r: r['gauss']['guard_ratio_max_abs_D_minus_L_over_bound']), COLORS[0]),
        ('gauss v3.1 form |D|/bound', S(lambda r: r['gauss']['v31_form_ratio_max_abs_D_over_bound']), COLORS[1]),
        ('gauss |D|/(bound+|L|)', S(lambda r: r['gauss']['implied_ratio_max_abs_D_over_bound_plus_abs_L']), COLORS[2]),
        ('divB/bound', S(lambda r: r['divB_over_bound_max_max']), COLORS[3])], logy=True)
    P('gauss-ledger.png', title + ': admitted gauss ledger', 'time (ns)', 'charge (C)', [
        ('max |L|', S(lambda r: r['gauss']['max_abs_L_C']), COLORS[0]), ('max |D|', S(lambda r: r['gauss']['max_abs_D_C']), COLORS[1])], logy=True)
    print(json.dumps({'plots': files}))


if __name__ == '__main__':
    if len(sys.argv) not in (3, 4):
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
