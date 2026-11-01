"""Validate extracted table JSON files against the agreed schema.

Usage: python3 -I validate_tables.py <extracted_dir> [prefix ...]
Only files whose names start with one of the prefixes are checked
(default prefixes: ion-, neutral-, rec-).
"""
import json
import math
import os
import sys

REQUIRED = ["table_id", "source_id", "locator", "process", "quantity", "x_name", "x", "y",
            "y_uncertainty", "original_units", "conversion_to_SI", "tag", "extraction_method", "notes"]
TAGS = {"measured", "recommended", "calculated", "model_derived"}


def check(path):
    errs = []
    with open(path) as f:
        d = json.load(f)
    for k in REQUIRED:
        if k not in d:
            errs.append("missing key " + k)
    if errs:
        return errs
    base = os.path.basename(path)[:-5]
    if d["table_id"] != base:
        errs.append("table_id %r != file name %r" % (d["table_id"], base))
    if d["tag"] not in TAGS:
        errs.append("bad tag " + str(d["tag"]))
    x, y = d["x"], d["y"]
    if not isinstance(x, list) or not isinstance(y, list):
        errs.append("x/y must be lists")
    elif len(x) != len(y):
        errs.append("len(x)=%d != len(y)=%d" % (len(x), len(y)))
    u = d["y_uncertainty"]
    if isinstance(u, list) and len(u) != len(y):
        errs.append("len(y_uncertainty) != len(y)")
    for v in (y if isinstance(y, list) else []):
        if isinstance(v, (int, float)) and (math.isnan(v) or math.isinf(v)):
            errs.append("non-finite y")
    return errs


def main():
    root = sys.argv[1]
    prefixes = tuple(sys.argv[2:]) or ("ion-", "neutral-", "rec-")
    bad = 0
    for fn in sorted(os.listdir(root)):
        if fn.endswith(".json") and fn.startswith(prefixes):
            e = check(os.path.join(root, fn))
            print(("OK   " if not e else "FAIL ") + fn + ("" if not e else ": " + "; ".join(e)))
            bad += bool(e)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
