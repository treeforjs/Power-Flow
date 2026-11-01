"""Digitization aid for 1-bit scanned scatter plots: list dark connected components inside a plot frame and
convert their bounding-box centres to data coordinates with a linear frame calibration.

Usage: python3 -I fig_digitize.py <image.png> <crop x0 y0 x1 y1> <frame px_left px_right data_left data_right>
       <frame py_bottom py_top data_bottom data_top> <exclude boxes as x0,y0,x1,y1;...> [margin]
Coordinates are crop-relative. Prints one line per component; writes nothing.
"""
import sys
from collections import deque

from PIL import Image


def main():
    a = sys.argv
    fn = a[1]
    cx0, cy0, cx1, cy1 = (int(v) for v in a[2:6])
    pxl, pxr, dxl, dxr = (float(v) for v in a[6:10])
    pyb, pyt, dyb, dyt = (float(v) for v in a[10:14])
    excl = []
    if len(a) > 14 and a[14]:
        for box in a[14].split(";"):
            excl.append(tuple(int(v) for v in box.split(",")))
    margin = int(a[15]) if len(a) > 15 else 4
    im = Image.open(fn).convert("L").crop((cx0, cy0, cx1, cy1))
    w, h = im.size
    px = im.load()

    def inside(x, y):
        if not (min(pxl, pxr) + margin < x < max(pxl, pxr) - margin and min(pyt, pyb) + margin < y < max(pyt, pyb) - margin):
            return False
        for (ex0, ey0, ex1, ey1) in excl:
            if ex0 <= x <= ex1 and ey0 <= y <= ey1:
                return False
        return True

    dark = [[px[x, y] < 128 and inside(x, y) for x in range(w)] for y in range(h)]
    seen = [[False] * w for _ in range(h)]
    out = []
    for y in range(h):
        for x in range(w):
            if dark[y][x] and not seen[y][x]:
                q = deque([(x, y)])
                seen[y][x] = True
                pts = []
                while q:
                    qx, qy = q.popleft()
                    pts.append((qx, qy))
                    for dx in (-1, 0, 1):
                        for dy in (-1, 0, 1):
                            nx, ny = qx + dx, qy + dy
                            if 0 <= nx < w and 0 <= ny < h and dark[ny][nx] and not seen[ny][nx]:
                                seen[ny][nx] = True
                                q.append((nx, ny))
                if len(pts) < 12:
                    continue
                xs = [p[0] for p in pts]
                ys = [p[1] for p in pts]
                bx, by = min(xs), min(ys)
                bw, bh = max(xs) - bx + 1, max(ys) - by + 1
                mx, my = bx + (bw - 1) / 2.0, by + (bh - 1) / 2.0
                X = dxl + (mx - pxl) * (dxr - dxl) / (pxr - pxl)
                Y = dyb + (my - pyb) * (dyt - dyb) / (pyt - pyb)
                out.append((X, Y, bx, by, bw, bh, len(pts), len(pts) / float(bw * bh)))
    for X, Y, bx, by, bw, bh, n, fill in sorted(out):
        print("x=%7.3f y=%7.3f  bbox=(%d,%d,%d,%d) n=%d fill=%.2f" % (X, Y, bx, by, bw, bh, n, fill))


if __name__ == "__main__":
    main()
