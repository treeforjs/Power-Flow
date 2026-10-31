"""Check extracted LXCat cross-section tables against their raw text download.

Usage: python3 -I verify_lxcat_cs.py <extracted_dir> <sources_dir>
For each lxcat-*-h2o-NN-*.json table: re-reads the raw block by its recorded
line range, re-tokenizes the numeric table, and requires exact token equality,
float equality, and that the total number of data lines in each raw file equals
the number of extracted points from that file.
"""
import glob
import hashlib
import json
import os
import re
import sys

ext, src = sys.argv[1], sys.argv[2]
DASH = re.compile(r"^-{5,}\s*$")
per_file = {}
bad = 0
for p in sorted(glob.glob(os.path.join(ext, "lxcat-*-h2o-[0-9][0-9]-*.json"))):
    t = json.load(open(p))
    raw = os.path.join(src, t["source_id"], "Cross section.txt")
    data = open(raw, "rb").read()
    assert hashlib.sha256(data).hexdigest() == t["lxcat"]["raw_file_sha256"], p
    lines = data.decode().replace("\r\n", "\n").split("\n")
    m = re.search(r"lines (\d+)-(\d+)", t["locator"])
    a, b = int(m.group(1)), int(m.group(2))
    blk = lines[a - 1:b]
    d = [i for i, l in enumerate(blk) if DASH.match(l)]
    rows = [l.split() for l in blk[d[0] + 1:d[1]]]
    xs = [r[0] for r in rows]
    ys = [r[1] for r in rows]
    ok = (xs == t["x_text"] and ys == t["y_text"]
          and [float(v) for v in xs] == t["x"] and [float(v) for v in ys] == t["y"]
          and len(t["x"]) == len(t["y"]))
    ok = ok and ("PROCESS: " + t["lxcat"]["process_label"]) in blk
    if not ok:
        bad += 1
        print("MISMATCH", p)
    per_file.setdefault(raw, 0)
    per_file[raw] += len(xs)
for raw, n in per_file.items():
    lines = open(raw).read().split("\n")
    # count numeric data lines between dash pairs
    cnt, inside = 0, False
    for l in lines:
        if DASH.match(l):
            inside = not inside
            continue
        if inside and l.strip():
            cnt += 1
    status = "OK" if cnt == n else "COUNT MISMATCH"
    print(status, os.path.basename(os.path.dirname(raw)), "raw data lines", cnt, "extracted points", n)
print("tables with mismatches:", bad)
