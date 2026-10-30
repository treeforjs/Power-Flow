#!/usr/bin/env python3
"""Copy one atomically published native checkpoint; never mutate the live output."""
import hashlib
import json
from pathlib import Path
import sys
import time

source, destination = map(Path, sys.argv[1:3])
destination.mkdir(parents=True, exist_ok=True)
for attempt in range(5):
    try:
        metadata = (source / "restart.meta").read_bytes()
        name = metadata.split(None, 1)[0].decode("ascii")
        if Path(name).name != name or not name.startswith("checkpoint_") or not name.endswith(".bin"):
            raise ValueError("Invalid published checkpoint name")
        with (source / name).open("rb") as stream:
            data = stream.read()
        if (source / "restart.meta").read_bytes() != metadata:
            continue
        target = destination / name
        target.write_bytes(data)
        (destination / "restart.meta").write_bytes(metadata)
        record = {"source_output_directory": str(source), "copied_checkpoint": str(target),
                  "checkpoint_bytes": len(data), "checkpoint_sha256": hashlib.sha256(data).hexdigest(),
                  "restart_meta_sha256": hashlib.sha256(metadata).hexdigest(),
                  "copy_unix_time_s": time.time(), "stable_metadata_match": True,
                  "read_only_source": True, "attempt": attempt + 1}
        (destination / "copy-provenance.json").write_text(json.dumps(record, indent=2) + "\n")
        print(json.dumps(record))
        break
    except FileNotFoundError:
        continue
else:
    raise SystemExit("Could not capture a stable checkpoint within five attempts")
