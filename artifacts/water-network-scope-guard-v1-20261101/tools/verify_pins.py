import hashlib
import json
from pathlib import Path

r=Path(__file__).resolve().parents[1]
p=json.loads((r/'INPUTS.json').read_text())
actual={str(f.relative_to(r)) for d in ('include','src','tests','tools') for f in (r/d).rglob('*')
        if f.is_file() and '__pycache__' not in f.parts}
assert actual==set(p['source_files']), 'source set changed'
for section in ('source_files','dependencies'):
    for path,expected in p[section].items():
        assert hashlib.sha256((r/path).read_bytes()).hexdigest()==expected, 'pinned input changed: '+path
print('PASS source/table/frozen dependency identities')
