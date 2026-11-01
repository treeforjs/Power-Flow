#!/usr/bin/env python3
"""Read-only check that the proposed adoption patch still matches its target."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

if len(sys.argv) != 2:
    raise SystemExit('usage: verify_base.py TARGET_SOURCE_DIRECTORY (read-only)')
a = Path(__file__).resolve().parent
target = Path(sys.argv[1]).resolve()
m = json.loads((a / 'base-sha256.json').read_text())
for name, expected in m['base_files_sha256'].items():
    actual = hashlib.sha256((target / name).read_bytes()).hexdigest()
    if actual != expected:
        raise SystemExit('BASE CHANGED: ' + name + '; rebase/review; do not force the patch')
subprocess.run(['git', 'apply', '--check', str(a / 'v3.2-audit-fixes.patch')], cwd=target, check=True)
print('PASS exact source hashes and git apply --check; target not modified')
