"""Parse a legacy LXCat cross-section text download into one JSON table per process block.

Usage:
  python3 -I parse_lxcat_cs.py <cross_section_txt> <source_id> <db_short> <tag> <out_dir> [--dry-run]

The numeric columns are copied verbatim: every value is kept both as the
original token (x_text / y_text) and as the float parsed from that token.
No interpolation, smoothing, unit change or point removal is performed.
"""
import hashlib
import json
import os
import re
import sys

KEYWORDS = {"ELASTIC", "EFFECTIVE", "EXCITATION", "IONIZATION", "ATTACHMENT", "ROTATION", "ROTATIONAL"}
DASHES = re.compile(r"^-{5,}\s*$")
META = ("SPECIES:", "PROCESS:", "PARAM.:", "COMMENT:", "UPDATED:", "COLUMNS:")


def slug(text, maxlen=48):
    s = text.lower()
    s = s.replace("+", "p").replace("->", "to").replace("^", "")
    s = re.sub(r"[^a-z0-9]+", "-", s).strip("-")
    return s[:maxlen].strip("-")


def parse(path):
    with open(path, encoding="utf-8", errors="strict") as fh:
        lines = fh.read().split("\n")
    header = {}
    group_comment = []
    blocks = []
    i = 0
    n = len(lines)
    in_db_header = False
    while i < n:
        line = lines[i]
        if line.startswith("DATABASE:"):
            in_db_header = True
        if in_db_header:
            m = re.match(r"^([A-Z][A-Z .]+):\s*(.*)$", line)
            if m and m.group(1) in ("DATABASE", "PERMLINK", "DESCRIPTION", "CONTACT", "HOW TO REFERENCE"):
                key = m.group(1)
                header[key] = m.group(2).strip()
                last_key = key
            elif line.startswith(" ") and line.strip() and header:
                header[last_key] = (header[last_key] + " " + line.strip()).strip()
            elif line.startswith("xxxxx") and header:
                in_db_header = False
            i += 1
            continue
        if line.startswith("COMMENT:") and not blocks and not any(lines[j].strip() in KEYWORDS for j in range(max(0, i - 3), i)):
            # data-group comment printed before the species banner
            group_comment.append(line[len("COMMENT:"):].strip())
            j = i + 1
            while j < n and lines[j].startswith(" ") and lines[j].strip():
                group_comment.append(lines[j].strip())
                j += 1
            i = j
            continue
        if line.strip() in KEYWORDS and line == line.rstrip() and line.strip() == line:
            start = i
            kind = line.strip()
            line2 = lines[i + 1]
            j = i + 2
            param_tokens = None
            if j < n and re.match(r"^\s*[-+0-9.eE]+(\s+[-+0-9.eE]+)*\s*$", lines[j]):
                param_tokens = lines[j].split()
                j += 1
            meta = {"SPECIES": None, "PROCESS": None, "PARAM.": None, "UPDATED": None, "COLUMNS": None}
            comments = []
            while j < n and not DASHES.match(lines[j]):
                l = lines[j]
                key = next((k for k in META if l.startswith(k)), None)
                if key == "COMMENT:":
                    comments.append(l[len("COMMENT:"):].strip())
                elif key:
                    meta[key[:-1]] = l[len(key):].strip()
                elif l.strip():
                    if comments:
                        comments[-1] = comments[-1] + " " + l.strip()
                    else:
                        comments.append(l.strip())
                j += 1
            if j >= n:
                raise SystemExit(f"unterminated block at line {start + 1}")
            j += 1
            xs, ys = [], []
            data_start = j
            while j < n and not DASHES.match(lines[j]):
                toks = lines[j].split()
                if len(toks) != 2:
                    raise SystemExit(f"bad data line {j + 1}: {lines[j]!r}")
                xs.append(toks[0])
                ys.append(toks[1])
                j += 1
            if j >= n:
                raise SystemExit(f"unterminated data table at line {data_start + 1}")
            blocks.append({
                "kind": kind,
                "line2": line2,
                "param_tokens": param_tokens,
                "meta": meta,
                "comments": comments,
                "x_text": xs,
                "y_text": ys,
                "line_first": start + 1,
                "line_last": j + 1,
            })
            i = j + 1
            continue
        i += 1
    return header, group_comment, blocks


def reaction_from_label(label):
    if label is None:
        return None
    body = re.sub(r",\s*(Elastic|Effective|Excitation|Ionization|Attachment|Rotation|Rotational)\s*$", "", label)
    body = re.sub(r"(?<![A-Za-z0-9])E(?![A-Za-z0-9(])", "e", body)
    return body


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    dry = "--dry-run" in sys.argv
    path, source_id, db_short, tag, out_dir = args
    header, group_comment, blocks = parse(path)
    with open(path, "rb") as fh:
        file_sha = hashlib.sha256(fh.read()).hexdigest()
    written = []
    used = set()
    for k, b in enumerate(blocks, start=1):
        label = b["meta"]["PROCESS"]
        x = [float(t) for t in b["x_text"]]
        y = [float(t) for t in b["y_text"]]
        threshold = None
        mass_ratio = None
        stat_weight_ratio = None
        pt = b["param_tokens"]
        if pt:
            if b["kind"] in ("ELASTIC", "EFFECTIVE"):
                mass_ratio = float(pt[0])
            else:
                threshold = float(pt[0])
            if len(pt) > 1:
                stat_weight_ratio = float(pt[1])
        base = f"lxcat-{db_short}-h2o-{k:02d}-{b['kind'].lower()}-{slug(b['line2'].split('->')[-1] if '->' in b['line2'] else b['line2'])}"
        if (b["comments"] and b["kind"] in ("ELASTIC",)) or base in used:
            base = base + "-" + slug(b["comments"][0] if b["comments"] else str(k), 24)
        comment_text = " ".join(b["comments"]).lower()
        process = reaction_from_label(label)
        if process is not None and "->" not in process and b["kind"] in ("ELASTIC", "EFFECTIVE"):
            process = process + " -> e + H2O"
        is_total = ("total ionization" in comment_text or "total ionisation" in comment_text
                    or "beb ionisation cross section" in comment_text)
        if b["kind"] == "IONIZATION" and is_total:
            base = base + "-total"
            process = "e + H2O -> 2e + (all cations), total ionization; LXCat label: " + str(label)
        if b["kind"] == "ATTACHMENT" and "total attachment" in comment_text:
            base = base + "-total"
            process = "e + H2O -> (all anions), total attachment; LXCat label: " + str(label)
        used.add(base)
        notes = []
        if b["comments"]:
            notes.append("LXCat COMMENT: " + " | ".join(b["comments"]))
        if group_comment:
            notes.append("LXCat data-group comment: " + " ".join(group_comment))
        notes.append("LXCat UPDATED: " + str(b["meta"]["UPDATED"]))
        notes.append("LXCat PARAM.: " + str(b["meta"]["PARAM."]))
        if b["kind"] == "ATTACHMENT" and pt is None:
            notes.append("Attachment block: LXCat format has no threshold line.")
        if b["kind"] in ("EXCITATION", "IONIZATION") and threshold == 0.0:
            notes.append("LXCat third line gives energy loss 0 eV; for emission/production blocks this is a placeholder, not a physical threshold.")
        rec = {
            "table_id": base,
            "source_id": source_id,
            "locator": f"LXCat {header.get('DATABASE', db_short)}; PROCESS: {label}; block {k} of {len(blocks)}; lines {b['line_first']}-{b['line_last']} of 'Cross section.txt'",
            "process": process,
            "quantity": "cross_section_m2",
            "x_name": "energy_eV",
            "x": x,
            "y": y,
            "y_uncertainty": None,
            "threshold_eV": threshold,
            "original_units": "energy in eV; cross section in m2 (LXCat COLUMNS: " + str(b["meta"]["COLUMNS"]) + ")",
            "conversion_to_SI": "none; values already in eV and m2",
            "tag": tag,
            "extraction_method": "verbatim text download",
            "notes": " ".join(notes),
            "lxcat": {
                "database": header.get("DATABASE"),
                "permlink": header.get("PERMLINK"),
                "collision_type": b["kind"],
                "line2": b["line2"],
                "third_line_tokens": pt,
                "mass_ratio_m_over_M": mass_ratio,
                "statistical_weight_ratio": stat_weight_ratio,
                "species": b["meta"]["SPECIES"],
                "process_label": label,
                "param": b["meta"]["PARAM."],
                "updated": b["meta"]["UPDATED"],
                "comments": b["comments"],
                "data_group_comment": group_comment,
                "raw_file": "Cross section.txt",
                "raw_file_sha256": file_sha,
            },
            "x_text": b["x_text"],
            "y_text": b["y_text"],
        }
        out = os.path.join(out_dir, base + ".json")
        if not dry:
            if os.path.exists(out):
                raise SystemExit(f"refusing to overwrite {out}")
            with open(out, "w", encoding="utf-8") as fh:
                json.dump(rec, fh, indent=1)
                fh.write("\n")
        written.append((base, b["kind"], label, threshold, len(x), x[0] if x else None, x[-1] if x else None))
    for w in written:
        print(*w, sep=" | ")
    print("database header:", json.dumps(header))
    print("group comment:", group_comment)


if __name__ == "__main__":
    main()
