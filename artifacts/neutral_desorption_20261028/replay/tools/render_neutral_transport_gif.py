#!/usr/bin/env python3
"""Render saved H/H2/H2O/OH number-density maps as a four-panel GIF."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.meta_path[:] = [finder for finder in sys.meta_path
                    if type(finder).__module__ != "_editable_skbc_mhdlab"]
from mhdlab.neutral_transport_render import render_transport


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="saved transport NPZ")
    parser.add_argument("--metadata", required=True, type=Path, help="model/provenance JSON sidecar")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--fps", type=float, default=10.)
    parser.add_argument("--dpi", type=int, default=110)
    parser.add_argument("--color-span-decades", type=float, help="limit each fixed color scale to this span below its maximum; lower positive values retain the lowest color")
    args = parser.parse_args()
    try:
        report = render_transport(args.input, args.metadata, args.output_dir, fps=args.fps, dpi=args.dpi,
                                  color_span_decades=args.color_span_decades)
    except (OSError, ValueError, TypeError, KeyError) as exc:
        parser.exit(2, f"neutral transport GIF: {exc}\n")
    print(json.dumps({"gif": str(args.output_dir / "neutral_transport.gif"),
                      "frames": report["frame_count"], "status": report["status"]}))


if __name__ == "__main__":
    main()
