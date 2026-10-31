"""List the vector drawings of dissertation Figures C.19 (PDF page index 189) and C.11 (index 181)."""
import sys, pymupdf
pdf = sys.argv[1]
with pymupdf.open(pdf) as doc:
    for idx in (189, 181):
        page = doc[idx]
        print("=== page index", idx)
        for line in page.get_text().splitlines():
            print("   text:", line)
        for i, d in enumerate(page.get_drawings()):
            items = d["items"]
            print(i, d["type"], "seq", d["seqno"], "color", d.get("color"), "fill", d.get("fill"),
                  "width", d.get("width"), "dashes", d.get("dashes"), "n_items", len(items),
                  "kinds", "".join(sorted(set(it[0] for it in items))), "rect", tuple(round(v, 2) for v in d["rect"]))
