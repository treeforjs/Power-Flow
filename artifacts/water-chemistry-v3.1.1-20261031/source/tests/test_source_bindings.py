"""Verify immutable source data and read-in-place dependencies against v3.1's seal."""
import hashlib
import json
from pathlib import Path

r = Path(__file__).resolve().parents[1]
old = Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1-20261031')
seal = json.loads((old / 'SEAL.json').read_text())
def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

checked = {}
for name, expected in seal['pinned_dependencies'].items():
    actual = sha(Path(name))
    assert actual == expected, ('dependency changed', name, expected, actual)
    checked[name] = actual
for name, expected in seal['files'].items():
    if name.startswith('data/sources/') or name.startswith('data/extracted/'):
        actual = sha(r / name)
        assert actual == expected, ('source data changed', name, expected, actual)
        checked[name] = actual
(r / 'evidence/immutable-inputs.json').write_text(json.dumps(checked, indent=1) + '\n')
print('PASS unchanged sealed source-data/dependency bindings:', len(checked))
