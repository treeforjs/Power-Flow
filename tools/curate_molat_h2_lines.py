#!/usr/bin/env python3
"""Verify retained MOLAT raw files and reproduce the H2 B/C line-list CSV.

Only radiative bound-bound transitions are handled; this is not a CR solver.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path


HEADER = ["upper_state", "upper_v", "upper_j", "lower_state", "lower_v",
          "lower_j", "einstein_a_per_s", "transition_wavenumber_per_cm",
          "vacuum_wavelength_nm", "source_file"]
STATES = {"BCan94": ("B", {-1, 1}), "CmCan94": ("C_minus", {0}),
          "CpCan94": ("C_plus", {-1, 1})}


def read_source(path, expected_hash):
    content = path.read_bytes()
    if hashlib.sha256(content).hexdigest() != expected_hash:
        raise ValueError(f"Source checksum differs from inspected MOLAT file: {path}")
    lines = content.decode("ascii").splitlines()
    if lines[1].split() != ["VU", "JU", "VL", "JL", "A", "(s-1)", "TR.E.(cm-1)"]:
        raise ValueError(f"Unexpected source columns: {path}")
    declared = int(lines[2])
    rows = [line.split() for line in lines[3:] if line.strip()]
    if len(rows) != declared or any(len(row) != 6 for row in rows):
        raise ValueError(f"Row count/shape mismatch: {path}")
    state, allowed_dj = STATES[path.name]
    seen = set()
    result = []
    for row in rows:
        vu, ju, vl, jl = map(int, row[:4])
        key = (vu, ju, vl, jl)
        a, wavenumber = map(Decimal, row[4:])
        if key in seen or min(key) < 0 or ju - jl not in allowed_dj:
            raise ValueError(f"Invalid/duplicate level labels: {path}: {row}")
        if not a.is_finite() or not wavenumber.is_finite() or a < 0 or wavenumber <= 0:
            raise ValueError(f"Invalid radiative data: {path}: {row}")
        seen.add(key)
        with localcontext() as context:
            context.prec = 15
            wavelength = Decimal("1e7") / wavenumber
        result.append([state, vu, ju, "X", vl, jl, row[4], row[5], str(wavelength), path.name])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path,
                        default=Path(__file__).resolve().parents[1] / "data/cr_model/molat_h2")
    parser.add_argument("--check", action="store_true", help="Compare instead of overwriting CSV")
    args = parser.parse_args()
    metadata = json.loads((args.data_dir / "metadata.json").read_text())
    result = []
    for entry in metadata["sources"]:
        rows = read_source(args.data_dir / "raw" / entry["file"], entry["sha256"])
        if len(rows) != entry["declared_line_count"]:
            raise ValueError(f"Metadata count differs for {entry['file']}")
        result.extend(rows)
    target = args.data_dir / "transitions.csv"
    if args.check:
        with target.open(newline="") as stream:
            current = list(csv.reader(stream))
        expected = [HEADER] + [[str(value) for value in row] for row in result]
        if current != expected:
            raise SystemExit("Normalized CSV differs from the retained source data.")
    else:
        with target.open("w", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(HEADER)
            writer.writerows(result)
    print(f"PASS: {len(result)} H2 radiative transitions match source counts and checksums.")


if __name__ == "__main__":
    main()
