"""Locate axis frame and dark connected components in a 1-bit scanned figure crop (figure digitization aid).

Usage: python3 -I fig_components.py <image.png> <x0> <y0> <x1> <y1>
Prints frame line positions (crop coordinates) and the bounding box, centroid, pixel count and fill ratio of
each connected dark component inside the crop. Writes nothing.
"""
import sys
from collections import deque

from PIL import Image


def main():
    fn = sys.argv[1]
    x0, y0, x1, y1 = (int(v) for v in sys.argv[2:6])
    im = Image.open(fn).convert("L").crop((x0, y0, x1, y1))
    w, h = im.size
    px = im.load()
    dark = [[px[x, y] < 128 for x in range(w)] for y in range(h)]
    colsum = [sum(dark[y][x] for y in range(h)) for x in range(w)]
    rowsum = [sum(dark[y][x] for x in range(w)) for y in range(h)]
    vcols = [x for x in range(w) if colsum[x] > 0.5 * h]
    hrows = [y for y in range(h) if rowsum[y] > 0.5 * w]
    print("vertical frame columns:", vcols)
    print("horizontal frame rows:", hrows)
    seen = [[False] * w for _ in range(h)]
    comps = []
    for y in range(h):
        for x in range(w):
            if dark[y][x] and not seen[y][x]:
                q = deque([(x, y)])
                seen[y][x] = True
                pts = []
                while q:
                    cx, cy = q.popleft()
                    pts.append((cx, cy))
                    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1), (1, -1), (-1, 1)):
                        nx, ny = cx + dx, cy + dy
                        if 0 <= nx < w and 0 <= ny < h and dark[ny][nx] and not seen[ny][nx]:
                            seen[ny][nx] = True
                            q.append((nx, ny))
                xs = [p[0] for p in pts]
                ys = [p[1] for p in pts]
                bw, bh = max(xs) - min(xs) + 1, max(ys) - min(ys) + 1
                comps.append((min(xs), min(ys), bw, bh, sum(xs) / len(xs), sum(ys) / len(ys), len(pts)))
    for c in sorted(comps, key=lambda c: c[4]):
        bx, by, bw, bh, cx, cy, n = c
        if n < 15 or bw > 0.5 * w or bh > 0.5 * h:
            continue
        print("bbox=(%d,%d,%d,%d) centre_bbox=(%.1f,%.1f) centroid=(%.1f,%.1f) n=%d fill=%.2f" % (
            bx, by, bw, bh, bx + (bw - 1) / 2, by + (bh - 1) / 2, cx, cy, n, n / float(bw * bh)))


if __name__ == "__main__":
    main()
