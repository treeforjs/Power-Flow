#!/usr/bin/env python3
"""Compare curated Song Tables 12/13 against the original PDF's numeric cells.

Requires PyMuPDF. The PDF is an external input, not distributed with this tool.
This uses geometric word extraction independently of the manual transcription.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re

import pymupdf


def numeric_cell(words, left, right, row_y):
    selected = sorted(
        (word for word in words
         if left < word[0] < right and abs(word[1] - row_y) < 3),
        key=lambda word: word[0],
    )
    text = " ".join(word[4] for word in selected)
    # The PDF's multiplication-sign glyph extracts as the character '3'.
    text = re.sub(r"\s+3\s+10", "e", text).replace("−", "-").replace(" ", "")
    return None if not text or set(text) == {"."} else Decimal(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument(
        "--data-dir", type=Path,
        default=Path(__file__).resolve().parents[1] / "data/cross_sections/song2021",
    )
    args = parser.parse_args()
    metadata = json.loads((args.data_dir / "metadata.json").read_text())
    digest = hashlib.sha256(args.pdf.read_bytes()).hexdigest()
    if digest != metadata["source_pdf_sha256"]:
        raise SystemExit("PDF SHA-256 differs from the inspected source; inspect its layout first.")
    pdf = pymupdf.open(args.pdf)
    positive = ["H2O+", "OH+", "O+", "O{2+}", "H2+", "H+", None]
    negative = ["H-", "O-", "OH-"]
    count = 0
    for entry in metadata["tables"]:
        if entry["source_table"] == 12:
            index = positive.index(entry["observed_ion"])
            page = pdf[20]
            energy_bounds = (35, 70)
            value_bounds = [(100, 160), (175, 222), (240, 280), (295, 350),
                            (365, 410), (426, 475), (489, 540)][index]
            low, high, factor = 145, 460, Decimal("1e-20")
        else:
            index = negative.index(entry["observed_ion"])
            page = pdf[21]
            energy_bounds = [(115, 155), (250, 285), (370, 410)][index]
            value_bounds = [(193, 214), (312, 334), (432, 454)][index]
            low, high = 160, 470
            factor = Decimal(["1e-22", "1e-23", "1e-24"][index])
        words = page.get_text("words")
        anchors = [word for word in words
                   if energy_bounds[0] < word[0] < energy_bounds[1]
                   and low < word[1] < high
                   and re.fullmatch(r"\d+(?:\.\d+)?", word[4])]
        expected = {}
        for word in anchors:
            value = numeric_cell(words, *value_bounds, word[1])
            if value is not None:
                expected[Decimal(word[4])] = value * factor
        with (args.data_dir / entry["file"]).open(newline="") as stream:
            rows = list(csv.DictReader(stream))
            actual = {Decimal(row["energy_ev"]): Decimal(row["cross_section_m2"])
                      for row in rows}
        if len(actual) != len(rows):
            raise SystemExit(f"Duplicate energy samples: {entry['file']}")
        if actual != expected:
            raise SystemExit(f"PDF/table mismatch: {entry['file']}")
        if len(actual) != entry["count"]:
            raise SystemExit(f"Metadata count mismatch: {entry['file']}")
        count += len(actual)
    print(f"PASS: {len(metadata['tables'])} tables, {count} numeric cross sections match PDF.")


if __name__ == "__main__":
    main()
