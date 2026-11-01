"""Reconstruct a rotated (landscape) table from pdftotext -bbox word boxes.

Usage: python3 -I song_rotated_table_words.py <bbox_html> [xmax_limit]
Keeps only words whose box is taller than wide (text rotated by 90 degrees),
clusters them into table columns by y-centre and into table rows by x-centre,
and prints a grid. Output is for transcription checking; it does not alter values.
"""
import html
import re
import sys

path = sys.argv[1]
xlim = float(sys.argv[2]) if len(sys.argv) > 2 else 1e9
words = []
for m in re.finditer(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">([^<]*)</word>', open(path, encoding="utf-8").read()):
    x0, y0, x1, y1 = map(float, m.group(1, 2, 3, 4))
    t = html.unescape(m.group(5))
    if (y1 - y0) > (x1 - x0) * 1.05 or len(t) <= 2:
        if x1 <= xlim:
            words.append(((x0 + x1) / 2, (y0 + y1) / 2, t, x0, y0, x1, y1))


def cluster(vals, tol):
    vals = sorted(set(round(v, 1) for v in vals))
    groups = []
    for v in vals:
        if groups and v - groups[-1][-1] <= tol:
            groups[-1].append(v)
        else:
            groups.append([v])
    return [sum(g) / len(g) for g in groups]


cols = cluster([w[1] for w in words], 3.0)
rows = cluster([w[0] for w in words], 3.0)


def nearest(v, cs):
    return min(range(len(cs)), key=lambda i: abs(cs[i] - v))


grid = {}
for w in words:
    r = nearest(w[0], rows)
    c = nearest(w[1], cols)
    grid.setdefault((r, c), []).append(w[2])
print("columns (y centres, page order bottom-to-top = table left-to-right):")
for ci in sorted(range(len(cols)), key=lambda i: -cols[i]):
    cells = [" ".join(grid.get((ri, ci), [])) for ri in range(len(rows))]
    if "--aligned" in sys.argv:
        print(f"y={cols[ci]:7.1f}: " + " | ".join((c if c else "_") for c in cells))
    else:
        print(f"y={cols[ci]:7.1f}: " + " | ".join(c for c in cells if c))
print("row x centres:", " ".join(f"{r:.1f}" for r in rows))
