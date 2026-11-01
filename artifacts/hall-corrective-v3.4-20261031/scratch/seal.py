#!/usr/bin/python -I
"""Seal Hall corrective v3.4: write SOURCE-SHA256.json into the runtime directory and remove write permission from the
sealed files and their directories. Usage: seal.py BUILD_LABEL"""
import hashlib, json, os, stat, sys
from pathlib import Path

R = Path('/home/treefor/projects/runtimes/hall-corrective-v3.4-20261031')
O = Path('/home/treefor/projects/.workspace/outputs/hall-corrective-v3.4-20261031')
label = sys.argv[1]
SKIP_DIRS = {R / 'build', R / 'coordinator-v7.1-copy/build', R / 'scratch'}   # scratch: operational helpers, not sealed
SKIP_FILES = {R / 'HANDOFF.txt', R / 'SOURCE-SHA256.json'}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


files, dirs = {}, set()
for root, ds, fs in os.walk(R):
    root = Path(root)
    ds[:] = sorted(d for d in ds if root / d not in SKIP_DIRS and d != '__pycache__')
    for f in sorted(fs):
        p = root / f
        if p in SKIP_FILES:
            continue
        files[str(p.relative_to(R))] = sha(p)
        if root != R and root != R / 'coordinator-v7.1-copy':
            dirs.add(root)
br = json.loads((O / label / 'build/build-receipt.json').read_text())
assert br['status'] == 'passed'
for p, h in br['own_sources_sha256'].items():
    assert files[str(Path(p).relative_to(R))] == h, p
manifest = {'schema': 'hall-v3.4-sealed-source-v1', 'runtime': str(R), 'build_label': label,
            'build_receipt': str(O / label / 'build/build-receipt.json'), 'binaries_sha256': br['binaries_sha256'],
            'isa_flags': br.get('isa_flags', [])}
iso = json.loads((O / 'isa-hook-rebuild.json').read_text())
manifest['binaries_identical_to_b13'] = iso['b14_vs_b13']['binaries_identical']
manifest['coordinator_builds'] = {
    'c4 (coordinator v5 T1-T7 parity; = c3)': iso['c4_vs_c3']['binary_sha256'],
    'r8 (coordinator v7.1 copy, Hall v3.4 + native v7.1; = r6)': iso['r8_vs_r6']['binary_sha256'],
    'r9 (coordinator v7.1 copy, Hall v3.4 + native v7; = r7)': iso['r9_vs_r7']['binary_sha256']}
manifest['files'] = dict(sorted(files.items()))
out = R / 'SOURCE-SHA256.json'
out.write_text(json.dumps(manifest, indent=1) + '\n')
W = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH
for rel in files:
    p = R / rel
    p.chmod(p.stat().st_mode & ~W)
out.chmod(out.stat().st_mode & ~W)
for d in sorted(dirs, key=lambda x: -len(x.parts)):
    d.chmod(d.stat().st_mode & ~W)
bad = [rel for rel, h in files.items() if sha(R / rel) != h or os.access(R / rel, os.W_OK)]
print(json.dumps({'manifest': str(out), 'manifest_sha256': sha(out), 'files': len(files), 'read_only_dirs': len(dirs), 'problems': bad}))
