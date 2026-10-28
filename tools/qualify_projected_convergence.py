"""Inspect independent projected-feature refinement and uncertainty evidence."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.meta_path[:] = [finder for finder in sys.meta_path
                    if type(finder).__module__ != "_editable_skbc_mhdlab"]
from mhdlab.projected_convergence import qualify_projected_convergence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("use a fresh output file")
    raw = args.config.read_bytes()
    result = qualify_projected_convergence(json.loads(raw))
    result["input_sha256"] = hashlib.sha256(raw).hexdigest()
    result["source_sha256"] = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        for name in ("tools/qualify_projected_convergence.py", "mhdlab/projected_convergence.py")}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(args.output), "resolution_criteria_met": result["resolution_criteria_met"],
                      "contrast_criteria_met": result["contrast_criteria_met"]}))


if __name__ == "__main__":
    main()
