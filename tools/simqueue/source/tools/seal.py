#!/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen
"""Write <runtime>/sealed-files.json: sha256 of every runtime file.

switchover.py refuses to install when a sealed file differs.
Usage: seal.py [--check]
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

RUNTIME = Path(__file__).resolve().parents[2]
SKIP_NAMES = {"sealed-files.json"}


def files() -> dict[str, str]:
    result = {}
    for path in sorted(RUNTIME.rglob("*")):
        if not path.is_file() or "__pycache__" in path.parts or path.name in SKIP_NAMES:
            continue
        result[str(path.relative_to(RUNTIME))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def main() -> int:
    current = files()
    target = RUNTIME / "sealed-files.json"
    if "--check" in sys.argv:
        recorded = json.loads(target.read_text())["files"]
        changed = sorted(k for k in set(recorded) | set(current) if recorded.get(k) != current.get(k))
        print("sealed files match" if not changed else "changed: " + ", ".join(changed))
        return 0 if not changed else 1
    target.write_text(json.dumps({
        "schema": "simqueue-sealed-files-v1",
        "runtime": str(RUNTIME),
        "sealed_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "interpreter": "/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen",
        "files": current,
    }, indent=2) + "\n")
    print(f"sealed {len(current)} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
