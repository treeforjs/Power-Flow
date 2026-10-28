"""Hash unchanged shared implementations and owned physical reference data."""
from pathlib import Path
import hashlib
import json

ROOT=Path('/workspace/Power-Flow')
OUT=Path(__file__).parent
excluded={'mhdlab/runner.py','mhdlab/__init__.py','mhdlab/mhd.py','mhdlab/materials.py','mhdlab/conductivity_tables.py'}
paths=set(ROOT.glob('mhdlab/**/*.py'))
for subtree in ['src','data/cr_model','data/cross_sections']:
    paths.update(path for path in (ROOT/subtree).rglob('*') if path.is_file())
paths.update([ROOT/'CMakeLists.txt',ROOT/'configs/water_fragment_mechanism.yaml'])
records={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
         for path in sorted(paths) if str(path.relative_to(ROOT)) not in excluded and path.exists()}
record={'schema':'mitl-mhd-retirement-independent-protected-pins-v1','files':records,
        'source_or_test_execution':False}
target=OUT/'protected-before.json'
assert not target.exists()
target.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
print(json.dumps({'file_count':len(records),'record_sha256':hashlib.sha256(target.read_bytes()).hexdigest()},sort_keys=True))
