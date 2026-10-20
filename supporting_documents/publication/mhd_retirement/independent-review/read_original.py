"""Read-only baseline inspection for the MHD retirement boundary."""
import ast
import hashlib
import json
from pathlib import Path

ROOT=Path('/workspace/Power-Flow')
OUT=Path(__file__).parent
paths=['mhdlab/runner.py','mhdlab/__init__.py','mhdlab/config.py','configs/mykonos.yaml',
       'tests/test_cr_data.py','tests/test_geometry.py','tests/test_neutrals.py']
record={'schema':'mitl-mhd-retirement-original-read-v1','files':{},'source_or_test_execution':False}
for name in paths:
    raw=(ROOT/name).read_bytes()
    record['files'][name]={'sha256':hashlib.sha256(raw).hexdigest()}
    if name=='mhdlab/runner.py':
        tree=ast.parse(raw)
        record['helpers']={node.name:{'ast':ast.dump(node,include_attributes=False),
            'source':ast.get_source_segment(raw.decode(),node)} for node in tree.body
            if isinstance(node,ast.FunctionDef) and node.name in ('_grid_shape_from_config','_thermionic_electron_energy_ev')}
    if name in ('mhdlab/__init__.py','mhdlab/config.py','tests/test_cr_data.py'):
        record['files'][name]['source']=raw.decode()
    if name=='configs/mykonos.yaml':
        lines=raw.decode().splitlines(True)
        start=next(i for i,line in enumerate(lines) if line=='cr_model:\n')
        stop=next((i for i in range(start+1,len(lines)) if lines[i].strip() and not lines[i][0].isspace()),len(lines))
        record['cr_model_text']=''.join(lines[start:stop]).rstrip()+'\n'
(OUT/'original-read.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
print(json.dumps({'schema':record['schema'],'paths':list(record['files']),'helpers':list(record['helpers']),
                  'record_sha256':hashlib.sha256((OUT/'original-read.json').read_bytes()).hexdigest()},sort_keys=True))
