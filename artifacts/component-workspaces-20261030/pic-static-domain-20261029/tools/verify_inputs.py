#!/usr/bin/env python3
"""Read-only hash guard. Does not launch, stop or signal any process."""
import hashlib
import json
from pathlib import Path
import sys

manifest=Path(__file__).resolve().parents[1]/'input-bindings.json'
data=json.loads(manifest.read_text())
failures=[]
for name,expected in data['sha256'].items():
    h=hashlib.sha256()
    try:
        with Path(name).open('rb') as stream:
            for block in iter(lambda:stream.read(1048576),b''):
                h.update(block)
        if h.hexdigest()!=expected:
            failures.append(name)
    except OSError:
        failures.append(name)
print(json.dumps({'passed':not failures,'files_checked':len(data['sha256']),
                  'failed_paths':failures,'manifest':str(manifest)}))
sys.exit(bool(failures))
