"""Assemble data/sources/ION-NEUTRAL-SOURCES.json from per-source registry records, the raw source directories
(sha256 of every file) and the extracted tables (prefixes ion-, neutral-, rec-).

Usage: python3 -I build_ion_neutral_index.py <registry_dir> <sources_dir> <extracted_dir> <output_json>
"""
import datetime
import hashlib
import json
import os
import sys

PREFIXES = ("ion-", "neutral-", "rec-")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    reg, src, ext, out = sys.argv[1:5]
    tables = {}
    for fn in sorted(os.listdir(ext)):
        if fn.endswith(".json") and fn.startswith(PREFIXES):
            with open(os.path.join(ext, fn)) as f:
                t = json.load(f)
            xs = [v for v in t.get("x", []) if isinstance(v, (int, float))]
            tables[t["table_id"]] = {
                "table_id": t["table_id"], "file": "data/extracted/" + fn, "sha256": sha256(os.path.join(ext, fn)),
                "source_id": t["source_id"], "process": t["process"], "quantity": t["quantity"],
                "x_name": t["x_name"], "n_points": len(t.get("y", [])),
                "x_range": [min(xs), max(xs)] if xs else None, "tag": t["tag"],
                "extraction_method": t["extraction_method"], "locator": t["locator"],
            }
    sources, missing, warnings = [], [], []
    seen_src = set()
    for fn in sorted(os.listdir(reg)):
        if not fn.endswith(".json"):
            continue
        with open(os.path.join(reg, fn)) as f:
            r = json.load(f)
        if fn.startswith("missing-"):
            r["id"] = fn[:-5]
            missing.append(r)
            continue
        sid = r["source_id"]
        seen_src.add(sid)
        d = os.path.join(src, sid)
        files = []
        if os.path.isdir(d):
            for root, _, names in os.walk(d):
                for n in sorted(names):
                    p = os.path.join(root, n)
                    files.append({"file": os.path.relpath(p, src), "bytes": os.path.getsize(p), "sha256": sha256(p)})
        else:
            warnings.append("source directory missing for " + sid)
        cf = r.get("crossref_file")
        for one in ([c.strip() for c in cf.split(";")] if isinstance(cf, str) else (cf or [])):
            if one and not os.path.exists(os.path.join(d, one)):
                warnings.append("crossref file not found for %s: %s" % (sid, one))
        tids = sorted(set(r.get("tables", [])) | {k for k, v in tables.items() if v["source_id"] == sid})
        for t in tids:
            if t not in tables:
                warnings.append("registry of %s lists unknown table %s" % (sid, t))
        entry = dict(r)
        entry["source_directory"] = "data/sources/" + sid
        entry["raw_files"] = files
        entry["tables"] = [tables[t] for t in tids if t in tables]
        sources.append(entry)
    for name in sorted(os.listdir(src)):
        if name.startswith(PREFIXES) and os.path.isdir(os.path.join(src, name)) and name not in seen_src:
            warnings.append("source directory without registry record: " + name)
    for k, v in tables.items():
        if v["source_id"] not in seen_src:
            warnings.append("table %s references unregistered source %s" % (k, v["source_id"]))
    missing.sort(key=lambda m: (str(m.get("channel")), m["id"]))
    doc = {
        "title": "Heavy-particle reaction data sources for the water-vapour PIC/MCC chemistry component (ions, neutral fragments, recombination)",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "scope": "Directories and tables with prefixes ion-, neutral-, rec- only (electron-impact data are indexed separately).",
        "conventions": {
            "tags": "measured | recommended | calculated | model_derived (per table; source-level tag is the dominant character)",
            "digitized_tables": "Any table whose extraction_method begins 'DIGITIZED FROM FIGURE' was read from a figure because no tabulation exists.",
            "expression_tables": "Tables with empty x/y carry the source's fit/rate expression parameters in an 'expression' object; no values were computed from expressions.",
            "sha256": "of every file in each source directory as stored at index time",
        },
        "channels": {
            "1": "H2O+ + H2O proton transfer and symmetric charge transfer",
            "2": "ion mobility/momentum transfer in H2O and fragment-ion reactions with H2O",
            "3": "negative-ion reactions, collisional and associative detachment",
            "4": "electron-ion dissociative recombination",
            "5": "ion-ion recombination",
            "6": "neutral fragments with H2O (transport and fast neutral reactions)",
        },
        "sources": sources,
        "n_sources": len(sources),
        "n_tables": sum(len(s["tables"]) for s in sources),
        "missing": missing,
        "index_warnings": warnings,
    }
    with open(out, "w") as f:
        json.dump(doc, f, indent=1, ensure_ascii=False)
    print("sources", len(sources), "tables", doc["n_tables"], "missing", len(missing))
    for w in warnings:
        print("WARN", w)


if __name__ == "__main__":
    main()
