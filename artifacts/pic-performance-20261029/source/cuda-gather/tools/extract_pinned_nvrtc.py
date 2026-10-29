#!/usr/bin/env python3
"""Extract the user-authorized pinned NVRTC wheel only into this investigation."""
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import stat
import time
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
RESULTS = Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/cuda-gather-investigation')
proposal = json.loads((RESULTS / 'toolchain-proposal.json').read_text())
package = proposal['nvrtc_wheel']
assert os.sched_getaffinity(0) == {6} and os.getpriority(os.PRIO_PROCESS, 0) >= 15
cache = ROOT / 'tools/packages'
cache.mkdir(parents=True, exist_ok=True)
wheel = cache / package['filename']
if not wheel.exists():
    partial = wheel.with_suffix('.partial')
    with urllib.request.urlopen(package['url'], timeout=60) as source, partial.open('wb') as out:
        shutil.copyfileobj(source, out, length=1024 * 1024)
    partial.replace(wheel)
digest = hashlib.sha256(wheel.read_bytes()).hexdigest()
if wheel.stat().st_size != package['size'] or digest != package['digests']['sha256']:
    raise RuntimeError('Pinned wheel hash/size mismatch; not extracting')
destination = ROOT / 'tools/nvrtc'
with zipfile.ZipFile(wheel) as archive:
    members = archive.infolist()
    total = sum(x.file_size for x in members)
    if total > proposal['disk_budget']['extracted_budget_bytes']:
        raise RuntimeError(f'Unpacked cap exceeded: {total} bytes')
    for member in members:
        path = PurePosixPath(member.filename)
        if path.is_absolute() or '..' in path.parts or stat.S_ISLNK(member.external_attr >> 16):
            raise RuntimeError('Unsafe wheel path/type')
    archive.extractall(destination)
files = {str(x.relative_to(destination)): {'size': x.stat().st_size, 'sha256': hashlib.sha256(x.read_bytes()).hexdigest()}
         for x in destination.rglob('*') if x.is_file()}
receipt = {'status': 'extracted_inside_assigned_prefix', 'user_authorized': True,
           'gpu_initialized_or_executed': False, 'unix_s': time.time(),
           'wheel': str(wheel), 'wheel_bytes': wheel.stat().st_size, 'wheel_sha256': digest,
           'extracted_prefix': str(destination), 'extracted_bytes': total,
           'peak_archive_plus_extracted_bytes': wheel.stat().st_size + total,
           'logical_cpu': sorted(os.sched_getaffinity(0)), 'nice': os.getpriority(os.PRIO_PROCESS, 0),
           'files': files}
(RESULTS / 'nvrtc-extraction.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({k: v for k, v in receipt.items() if k != 'files'}, indent=2))
