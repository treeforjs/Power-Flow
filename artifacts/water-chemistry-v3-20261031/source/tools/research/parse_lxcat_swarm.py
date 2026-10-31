"""Parse a legacy LXCat swarm/transport text download into extracted-table JSON.

Usage:
  python3 -I parse_lxcat_swarm.py <lxcat_txt> <out_json> <spec_json>

<spec_json> is a small JSON file holding the per-table metadata that cannot be
read from the LXCat text (table_id, source_id, quantity, tag, notes, ...).
The numeric columns are copied exactly as printed: each value is stored both
as the original token string ("x_text"/"y_text") and as the float parsed from
that token. Nothing is interpolated, rounded, smoothed or rescaled.

The file must contain exactly one SPECIES/PROCESS block; otherwise the script
stops with an error so that blocks are never merged silently.
"""
import hashlib
import json
import re
import sys


def parse_blocks(text):
    lines = text.splitlines()
    blocks = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("SPECIES:"):
            hdr = {}
            j = i
            while j < len(lines) and not re.match(r"^-{5,}\s*$", lines[j]):
                m = re.match(r"^([A-Z.]+):\s*(.*)$", lines[j])
                if m:
                    key = m.group(1).rstrip(".")
                    hdr.setdefault(key, []).append(m.group(2).rstrip())
                j += 1
            if j >= len(lines):
                raise SystemExit("unterminated header block")
            j += 1
            rows = []
            while j < len(lines) and not re.match(r"^-{5,}\s*$", lines[j]):
                toks = lines[j].split()
                if len(toks) != 2:
                    raise SystemExit(f"unexpected data line: {lines[j]!r}")
                rows.append(toks)
                j += 1
            if j >= len(lines):
                raise SystemExit("unterminated data table")
            blocks.append({"header": hdr, "rows": rows})
            i = j + 1
        else:
            i += 1
    return blocks


def main():
    src, out, spec_path = sys.argv[1], sys.argv[2], sys.argv[3]
    raw = open(src, "rb").read()
    text = raw.decode("utf-8")
    with open(spec_path, encoding="utf-8") as fh:
        spec = json.load(fh)
    blocks = parse_blocks(text)
    if len(blocks) != 1:
        raise SystemExit(f"expected exactly one block, found {len(blocks)}")
    b = blocks[0]
    hdr = b["header"]
    x_text = [r[0] for r in b["rows"]]
    y_text = [r[1] for r in b["rows"]]
    x = [float(v) for v in x_text]
    y = [float(v) for v in y_text]
    if any(x[k + 1] <= x[k] for k in range(len(x) - 1)):
        monotonic = False
    else:
        monotonic = True
    table = {
        "table_id": spec["table_id"],
        "source_id": spec["source_id"],
        "locator": spec["locator"],
        "process": spec.get("process", "e in H2O swarm"),
        "quantity": spec["quantity"],
        "x_name": "E_over_N_Td",
        "x": x,
        "y": y,
        "y_uncertainty": None,
        "threshold_eV": None,
        "original_units": spec["original_units"],
        "conversion_to_SI": spec["conversion_to_SI"],
        "tag": spec["tag"],
        "extraction_method": "verbatim text download",
        "notes": spec["notes"],
        "x_text": x_text,
        "y_text": y_text,
        "y_original": y,
        "n_points": len(x),
        "x_strictly_increasing": monotonic,
        "lxcat_header": {k: v for k, v in hdr.items()},
        "raw_file": spec["raw_file"],
        "raw_file_sha256": hashlib.sha256(raw).hexdigest(),
    }
    for k in ("citation_as_given_by_lxcat", "primary_citation", "primary_doi",
              "crossref_record", "temperature_K", "pressure", "derived_quantities_note"):
        if k in spec:
            table[k] = spec[k]
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(table, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    print(out, len(x), x[0], x[-1])


if __name__ == "__main__":
    main()
