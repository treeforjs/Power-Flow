"""Marker detection for 1-bit scanned scatter plots with open circles and filled markers (digitization aid).

Filled markers: binary erosion with a square of half-width k removes thin strokes (rings, axes, ticks); remaining
connected components give filled-marker cores (centre = core bbox centre).
Open circles: fixed-radius circular Hough transform on dark pixels outside dilated filled-marker regions.

Usage: python3 -I fig_markers.py <image.png> <crop x0 y0 x1 y1> <plot box px_left py_top px_right py_bottom>
       <exclude boxes "x0,y0,x1,y1;..."> <ring_radius_px> <erode_k> <hough_threshold_fraction> <overlay_out.png>
Coordinates (except crop) are crop-relative. Prints detections; writes only the overlay PNG given.
"""
import math
import sys
from collections import deque

from PIL import Image, ImageDraw


def components(mask, w, h, minn):
    seen = [[False] * w for _ in range(h)]
    res = []
    for y in range(h):
        for x in range(w):
            if mask[y][x] and not seen[y][x]:
                q = deque([(x, y)])
                seen[y][x] = True
                pts = []
                while q:
                    a, b = q.popleft()
                    pts.append((a, b))
                    for dx in (-1, 0, 1):
                        for dy in (-1, 0, 1):
                            nx, ny = a + dx, b + dy
                            if 0 <= nx < w and 0 <= ny < h and mask[ny][nx] and not seen[ny][nx]:
                                seen[ny][nx] = True
                                q.append((nx, ny))
                if len(pts) >= minn:
                    xs = [p[0] for p in pts]
                    ys = [p[1] for p in pts]
                    res.append((min(xs), min(ys), max(xs), max(ys), len(pts)))
    return res


def main():
    a = sys.argv
    im = Image.open(a[1]).convert("L").crop(tuple(int(v) for v in a[2:6]))
    pl, pt, pr, pb = (int(v) for v in a[6:10])
    excl = [tuple(int(v) for v in b.split(",")) for b in a[10].split(";")] if a[10] else []
    R = float(a[11])
    k = int(a[12])
    thr = float(a[13])
    outpng = a[14]
    w, h = im.size
    px = im.load()

    def ok(x, y):
        if not (pl + 4 < x < pr - 4 and pt + 4 < y < pb - 4):
            return False
        return not any(e[0] <= x <= e[2] and e[1] <= y <= e[3] for e in excl)

    dark = [[px[x, y] < 128 and ok(x, y) for x in range(w)] for y in range(h)]
    # erosion
    er = [[False] * w for _ in range(h)]
    for y in range(k, h - k):
        for x in range(k, w - k):
            if dark[y][x]:
                good = True
                for yy in range(y - k, y + k + 1):
                    row = dark[yy]
                    for xx in range(x - k, x + k + 1):
                        if not row[xx]:
                            good = False
                            break
                    if not good:
                        break
                er[y][x] = good
    filled = components(er, w, h, 4)
    print("# filled-marker cores (crop px centre)")
    fcent = []
    for (x0, y0, x1, y1, n) in filled:
        cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0
        fcent.append((cx, cy, x1 - x0 + 1, y1 - y0 + 1, n))
        print("F %.1f %.1f core=%dx%d n=%d" % (cx, cy, x1 - x0 + 1, y1 - y0 + 1, n))
    # mask filled regions (dilate cores by k+2)
    dk = k + 3
    vote_ok = [row[:] for row in dark]
    for (x0, y0, x1, y1, n) in filled:
        for y in range(max(0, y0 - dk), min(h, y1 + dk + 1)):
            for x in range(max(0, x0 - dk), min(w, x1 + dk + 1)):
                vote_ok[y][x] = False
    acc = {}
    offs = []
    for t in range(0, 360, 4):
        ang = math.radians(t)
        offs.append((int(round(R * math.cos(ang))), int(round(R * math.sin(ang)))))
    offs = list(set(offs))
    for y in range(h):
        for x in range(w):
            if vote_ok[y][x]:
                for (dx, dy) in offs:
                    key = (x - dx, y - dy)
                    acc[key] = acc.get(key, 0) + 1
    # score = fraction of ring sample points that are dark
    def ring_score(cx, cy):
        hits = 0
        for (dx, dy) in offs:
            xx, yy = cx + dx, cy + dy
            if 0 <= xx < w and 0 <= yy < h and px[xx, yy] < 128:
                hits += 1
        return hits / float(len(offs))

    cands = sorted(((v, c) for c, v in acc.items()), reverse=True)
    circles = []
    for v, (cx, cy) in cands:
        if v < 0.5 * len(offs):
            break
        if any((cx - ox) ** 2 + (cy - oy) ** 2 < (0.9 * R) ** 2 for (ox, oy, s) in circles):
            continue
        s = ring_score(cx, cy)
        if s >= thr and (px[cx, cy] >= 128 or thr < 0):  # open centre
            circles.append((cx, cy, s))
    print("# open circles (crop px centre, ring dark fraction)")
    for (cx, cy, s) in sorted(circles):
        print("C %d %d score=%.2f" % (cx, cy, s))
    ov = im.convert("RGB")
    d = ImageDraw.Draw(ov)
    for (cx, cy, s) in circles:
        d.ellipse([cx - R, cy - R, cx + R, cy + R], outline=(255, 0, 0))
        d.point((cx, cy), fill=(255, 0, 0))
    for (cx, cy, bw, bh, n) in fcent:
        d.rectangle([cx - 2, cy - 2, cx + 2, cy + 2], outline=(0, 0, 255))
    ov = ov.resize((w * 2, h * 2), Image.NEAREST)
    ov.save(outpng)


if __name__ == "__main__":
    main()
