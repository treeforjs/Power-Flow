#!/usr/bin/env python3
"""Reconstruct the abstract text from saved OpenAlex work JSON files.

Usage:
    metal_openalex_abstract.py <openalex_json> [<openalex_json> ...]

OpenAlex stores abstracts as an inverted index {word: [positions]}; this
script inverts it back to running text. Input files are treated as data only.
"""
import json
import sys


def abstract(path: str) -> str:
    with open(path, encoding="utf-8") as f:
        w = json.load(f)
    inv = w.get("abstract_inverted_index") or {}
    if not inv:
        return ""
    pos = {}
    for word, idxs in inv.items():
        for i in idxs:
            pos[i] = word
    return " ".join(pos[i] for i in sorted(pos))


def main() -> int:
    for p in sys.argv[1:]:
        print(f"== {p}")
        text = abstract(p)
        print(text if text else "(no abstract in record)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
