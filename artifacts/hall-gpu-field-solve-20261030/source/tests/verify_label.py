"""Recompute every sha256 recorded in one label's receipts and compare (read-only)."""
import hashlib, json, sys
from pathlib import Path
root = Path('/home/treefor/sim-workspace/powerflow/outputs/hall-gpu-resident-v1-20261030/field') / sys.argv[1]
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''): h.update(c)
    return h.hexdigest()
checked, bad = 0, []
def check(path, expected):
    global checked
    checked += 1
    if sha(path) != expected: bad.append(str(path))
b = json.loads((root / 'build.json').read_text())
for n, h in b['binaries'].items(): check(root / 'build' / n, h)
check(b['input_sha256_manifest'], b['input_sha256_manifest_sha256'])
for p, h in json.loads(Path(b['input_sha256_manifest']).read_text()).items():
    if p.startswith(str(root)): check(p, h)
for ph in sorted(root.glob('meshes/*/*-phase.json')):
    r = json.loads(ph.read_text())
    for c in r.get('commands', []): check(c['log'], c['log_sha256'])
    if 'result' in r: check(r['result'], r['result_sha256'])
    for f in r.get('fixtures', {}).values(): check(f['file'], f['sha256'])
    for k, v in r.get('bit_identity', {}).items():
        for run, e in v['processes'].items():
            check(ph.parent / run / k.split('/')[0] / ('x_gpu_%s.bin' % k.split('/')[1]), e['x_file_sha256'])
print(json.dumps({'label': sys.argv[1], 'checked': checked, 'mismatched': bad}, indent=1))
sys.exit(1 if bad else 0)
