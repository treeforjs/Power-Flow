#!/usr/bin/env python3
"""Gridded, enlarged crops of the deskewed figure for manual marker reading.
usage: grid_crops.py IMAGE.png OUT_DIR PREFIX OFFSET x0,y0,x1,y1 [...]
The image is deskewed exactly as dr_digitize.py does (same skew estimate), so read pixel coordinates are in
the calibrated frame. Grid lines every 10 px (labels every 50 px), shifted by OFFSET px for the second pass."""
import math
import sys

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import dr_digitize as dd  # noqa: E402

path, out, prefix, off = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
ink = dd.ink_mask(path)
th, _ = dd.skew_angle(ink)
img = Image.fromarray((ink * 255).astype(np.uint8)).rotate(math.degrees(th), resample=Image.NEAREST, fillcolor=0)
ink = np.array(img) > 127
base = Image.fromarray(((~ink) * 255).astype(np.uint8)).convert('RGB')
Z = 2  # overridden per crop: 6 for windows narrower than 100 px
for k, box in enumerate(sys.argv[5:]):
    x0, y0, x1, y1 = map(int, box.split(','))
    Z = 6 if x1 - x0 < 100 else 2
    c = base.crop((x0, y0, x1, y1)).resize(((x1 - x0) * Z, (y1 - y0) * Z), Image.NEAREST)
    d = ImageDraw.Draw(c)
    for x in range(x0 - (x0 - off) % 10 + 10, x1, 10):
        major = (x - off) % 50 == 0
        d.line([((x - x0) * Z, 0), ((x - x0) * Z, (y1 - y0) * Z)], fill=(0, 150, 0) if major else (190, 225, 190), width=1)
        if major:
            d.text(((x - x0) * Z + 2, 2), str(x), fill=(0, 110, 0))
    for y in range(y0 - (y0 - off) % 10 + 10, y1, 10):
        major = (y - off) % 50 == 0
        d.line([(0, (y - y0) * Z), ((x1 - x0) * Z, (y - y0) * Z)], fill=(0, 150, 0) if major else (190, 225, 190), width=1)
        if major:
            d.text((2, (y - y0) * Z + 2), str(y), fill=(0, 110, 0))
    c.save('%s/%s-%d.png' % (out, prefix, k))
print('skew', th)
