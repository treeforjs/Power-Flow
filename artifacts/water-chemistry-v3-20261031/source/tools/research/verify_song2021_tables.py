"""Second check of the Song 2021 transcriptions against the PDF text layer.

Usage: python3 -I verify_song2021_tables.py <extracted_dir> <table1_layout.txt> <table11_layout.txt> <table8_layout.txt> <p9_bbox.html>

table1/table11/table8 files: `pdftotext -layout` output cropped to the table region.
p9_bbox.html: `pdftotext -bbox` output of PDF page 9 (rotated Table 2).
Every transcribed token must equal the corresponding text-layer token.
"""
import html
import json
import os
import re
import sys

ext, t1, t11, t8, p9 = sys.argv[1:6]
NUM = re.compile(r"^\d+(\.\d+)?$")
fails = 0


def check(name, got, want):
    global fails
    ok = got == want
    print(("OK  " if ok else "FAIL"), name, len(got), "tokens")
    if not ok:
        fails += 1
        for i, (a, b) in enumerate(zip(got, want)):
            if a != b:
                print("   first difference at", i, a, b)
                break


def two_pair_rows(path):
    left_x, left_y, right_x, right_y = [], [], [], []
    for line in open(path, encoding="utf-8"):
        toks = line.split()
        if toks and all(NUM.match(t) for t in toks) and len(toks) in (2, 4):
            left_x.append(toks[0]); left_y.append(toks[1])
            if len(toks) == 4:
                right_x.append(toks[2]); right_y.append(toks[3])
    return left_x + right_x, left_y + right_y


def load(tid):
    return json.load(open(os.path.join(ext, tid + ".json"), encoding="utf-8"))


x, y = two_pair_rows(t1)
r = load("song2021-table1-total-scattering")
check("Table 1 energies", r["x_text"], x)
check("Table 1 TCS", r["y_text"], y)
check("Table 1 floats", [float(v) for v in r["x_text"]], r["x"])

x, y = two_pair_rows(t11)
r = load("song2021-table11-O1S0-production")
check("Table 11 energies", r["x_text"], x)
check("Table 11 values", r["y_text"], y)

vals = [t for t in open(t8, encoding="utf-8").read().split() if re.match(r"^\d\.\d\d$", t)]
r = load("song2021-table8-electronic-state-energies")
check("Table 8 energies", r["y_text"], vals)

# Table 2 from word boxes: rotated words; columns are y bands, rows are x positions
words = []
for m in re.finditer(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">([^<]*)</word>',
                     open(p9, encoding="utf-8").read()):
    x0, y0, x1, y1 = map(float, m.group(1, 2, 3, 4))
    t = html.unescape(m.group(5))
    if (y1 - y0) > (x1 - x0) and re.match(r"^\d\.\d{3}$", t) and x1 < 320:
        words.append(((y0 + y1) / 2, (x0 + x1) / 2, t))
bands = {}
for yc, xc, t in words:
    key = round(yc)
    k = next((k for k in bands if abs(k - key) <= 3), key)
    bands.setdefault(k, []).append((xc, t))
cols = [[t for _, t in sorted(v)] for _, v in sorted(bands.items(), key=lambda kv: -kv[0])]
energies = ["2", "4", "6", "8", "10", "15", "20", "30", "50", "100"]
if len(cols) != 20:
    print("FAIL Table 2: expected 20 value columns, found", len(cols))
    fails += 1
else:
    for i, E in enumerate(energies):
        r = load(f"song2021-table2-elastic-dcs-{E}eV")
        check(f"Table 2 {E} eV DCS", r["y_text"], cols[2 * i])
        check(f"Table 2 {E} eV delta", r["y_uncertainty_text"], cols[2 * i + 1])
        assert len(r["x"]) == len(r["y"]) == len(r["y_uncertainty"])
print("failures:", fails)
