"""Command-line driver for the Power Flow MHD/CR prototype."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from mhdlab import load_config, run_from_config


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the reduced MHD / neutral / CR prototype.")
    parser.add_argument("--config", default="configs/mykonos.yaml", help="YAML/JSON run configuration")
    parser.add_argument(
        "--mhd-backend",
        choices=["auto", "numpy", "cuda", "cupy", "gpu"],
        default=None,
        help="Override mhd.backend from the config",
    )
    parser.add_argument(
        "--neutral-backend",
        choices=["auto", "numpy", "cuda", "cupy", "gpu"],
        default=None,
        help="Override neutrals.backend from the config",
    )
    parser.add_argument("--print-run-dir", action="store_true", help="Print only the output run directory")
    parser.add_argument("--quiet", action="store_true", help="Suppress progress messages")
    args = parser.parse_args(argv)

    config = load_config(args.config)
    config.setdefault("project_root", str(Path(__file__).resolve().parent))
    if args.mhd_backend:
        config.setdefault("mhd", {})["backend"] = args.mhd_backend
    if args.neutral_backend:
        config.setdefault("neutrals", {})["backend"] = args.neutral_backend
    progress = None if args.quiet or args.print_run_dir else _print_progress
    result = run_from_config(config, progress=progress)
    if args.print_run_dir:
        print(result["run_dir"])
    else:
        print(
            json.dumps(
                {
                    "run_dir": result["run_dir"],
                    "fields_hdf5": result["fields_hdf5"],
                    "sample_count": result["sample_count"],
                    "rl_fit": result["rl_fit"],
                    "spectra": result["spectra"],
                },
                indent=2,
            )
        )
    return 0


def _print_progress(message: str) -> None:
    print(f"[power-flow] {message}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
