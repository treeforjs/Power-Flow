"""Preserve and repair a pre-execution new-test symbol typo."""
from pathlib import Path
import ast
import hashlib
import json

R=Path('/workspace/Power-Flow')
W=Path('/workspace/research_work/mitl-github-push-20261020T0601/mhd-retirement')
P=R/'tests/test_legacy_mhd_retirement.py'
initial_record=(W/'source-freeze.json').read_bytes()
record=json.loads(initial_record)
current=P.read_text()
initial=current.replace('"ConstantConductivityHalfSpace",','"HalfSpaceHeating",')
assert hashlib.sha256(initial.encode()).hexdigest()==record['owned_file_hashes']['tests/test_legacy_mhd_retirement.py']['after_sha256']
assert current.count('"ConstantConductivityHalfSpace",')==1
compile(current,str(P),'exec')
tree=ast.parse((R/'mhdlab/heating_benchmarks.py').read_text())
assert any(isinstance(node,ast.ClassDef) and node.name=='ConstantConductivityHalfSpace' for node in tree.body)
(W/'source-freeze.initial.json').write_bytes(initial_record)
(W/'test_legacy_mhd_retirement.initial.py').write_text(initial)
record['owned_file_hashes']['tests/test_legacy_mhd_retirement.py']['after_sha256']=hashlib.sha256(P.read_bytes()).hexdigest()
record['preexecution_test_symbol_correction']='Text-only class name correction HalfSpaceHeating -> ConstantConductivityHalfSpace; no target execution occurred before this correction.'
(W/'source-freeze.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
print(json.dumps({'source_freeze_sha256':hashlib.sha256((W/'source-freeze.json').read_bytes()).hexdigest(),
                  'new_test_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),
                  'initial_test_sha256':hashlib.sha256(initial.encode()).hexdigest(),
                  'target_executions':0},sort_keys=True))
