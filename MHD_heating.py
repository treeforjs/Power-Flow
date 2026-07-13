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
    parser.add_argument("--print-run-dir", action="store_true", help="Print only the output run directory")
    args = parser.parse_args(argv)

    config = load_config(args.config)
    config.setdefault("project_root", str(Path(__file__).resolve().parent))
    result = run_from_config(config)
    if args.print_run_dir:
        print(result["run_dir"])
    else:
        print(json.dumps({"run_dir": result["run_dir"], "rl_fit": result["rl_fit"], "spectra": result["spectra"]}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
