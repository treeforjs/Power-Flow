#!/usr/bin/env python3
"""Tile page images into a labelled contact sheet for visual page finding.

Usage:
    metal_contact_sheet.py <out.png> <cols> <thumb_width> <img1> [<img2> ...]
"""
import os
import sys

from PIL import Image, ImageDraw


def main() -> int:
    if len(sys.argv) < 5:
        print(__doc__)
        return 2
    out, cols, tw = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    paths = sys.argv[4:]
    thumbs = []
    for p in paths:
        im = Image.open(p).convert("L")
        h = int(im.height * tw / im.width)
        thumbs.append((os.path.basename(p), im.resize((tw, h))))
    th = max(t.height for _, t in thumbs) + 14
    rows = (len(thumbs) + cols - 1) // cols
    sheet = Image.new("L", (cols * tw, rows * th), 255)
    draw = ImageDraw.Draw(sheet)
    for i, (name, t) in enumerate(thumbs):
        x, y = (i % cols) * tw, (i // cols) * th
        sheet.paste(t, (x, y + 14))
        draw.text((x + 2, y), name, fill=0)
    sheet.save(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
