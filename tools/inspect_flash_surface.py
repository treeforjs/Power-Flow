"""Inspect a canonical SI surface history, not raw FLASH plot or checkpoint files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mhdlab.flash_surface import load_surface_history, thermal_width_report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("surface_hdf5", type=Path)
    parser.add_argument("--edge-fraction", type=float, default=0.1, help="fraction of foil width at EACH edge")
    parser.add_argument("--center-fraction", type=float, default=0.2)
    parser.add_argument("--time-s", type=float, help="optional interpolated time in seconds; no extrapolation")
    args = parser.parse_args()
    try:
        history = load_surface_history(args.surface_hdf5)
        report = thermal_width_report(
            history,
            edge_fraction=args.edge_fraction,
            center_fraction=args.center_fraction,
            time_s=args.time_s,
        )
    except ValueError as exc:
        parser.error(str(exc))
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
