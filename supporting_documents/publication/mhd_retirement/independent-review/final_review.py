"""Independent file/AST/metadata review; no target imports or tests."""
import ast
import hashlib
import json
from pathlib import Path

ROOT=Path('/workspace/Power-Flow')
WORK=Path('/workspace/research_work/mitl-github-push-20261020T0601/mhd-retirement')
OUT=Path(__file__).parent
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
original=json.loads((OUT/'original-read.json').read_text())
freeze=json.loads((WORK/'source-freeze.json').read_text())
assert sha(WORK/'source-freeze.json')=='fc3ba952b90456e2e27c4836baffe0327b00aa0d5fd2382bbaabcd6a5281cdcb'
for path, record in freeze['owned_file_hashes'].items():
    actual=sha(ROOT/path) if (ROOT/path).exists() else None
    assert actual==record['after_sha256'],path
    if record['before_sha256'] is not None:
        assert sha(WORK/'before-files'/path)==record['before_sha256'],path

runner_raw=(ROOT/'mhdlab/runner.py').read_text()
runner=ast.parse(runner_raw)
functions={node.name:node for node in runner.body if isinstance(node,ast.FunctionDef)}
assert set(functions)==set(original['helpers'])
for name,node in functions.items():
    assert ast.dump(node,include_attributes=False)==original['helpers'][name]['ast'],name
    assert ast.get_source_segment(runner_raw,node)==original['helpers'][name]['source'],name
imports=[]
for node in runner.body:
    if isinstance(node,ast.Import):
        imports.extend(alias.name for alias in node.names)
    elif isinstance(node,ast.ImportFrom):
        imports.append('.'*node.level+(node.module or ''))
assert set(imports)=={'__future__','numpy','.constants','.geometry'}
assert all(isinstance(node,(ast.Import,ast.ImportFrom,ast.FunctionDef)) or
           isinstance(node,ast.Expr) and isinstance(node.value,ast.Constant) and isinstance(node.value.value,str)
           for node in runner.body)

init=ast.parse((ROOT/'mhdlab/__init__.py').read_text())
init_imports=[node for node in init.body if isinstance(node,(ast.Import,ast.ImportFrom))]
assert len(init_imports)==1
imp=init_imports[0]
assert isinstance(imp,ast.ImportFrom) and imp.level==1 and imp.module=='config'
assert [(alias.name,alias.asname) for alias in imp.names]==[('load_config',None)]
exports=[node for node in init.body if isinstance(node,ast.Assign) and
         any(isinstance(target,ast.Name) and target.id=='__all__' for target in node.targets)]
assert len(exports)==1 and ast.literal_eval(exports[0].value)==['load_config']

config_before=(WORK/'before-files/configs/mykonos.yaml').read_bytes()
lines=config_before.splitlines(keepends=True)
start=lines.index(b'cr_model:\n')
stop=next(i for i in range(start+1,len(lines)) if lines[i].strip() and not lines[i][:1].isspace())
subsection=b''.join(lines[start:stop])
assert (ROOT/'configs/cr_reference_library.yaml').read_bytes()==subsection
new_config=(ROOT/'configs/cr_reference_library.yaml').read_text()
for field in ('references','reaction_overview','rovibrational_manifest'):
    line=next(line for line in new_config.splitlines() if line.strip().startswith(field+':'))
    relative=line.split(':',1)[1].strip()
    assert (ROOT/'configs'/relative).resolve().is_file()
test_before=original['files']['tests/test_cr_data.py']['source']
assert (ROOT/'tests/test_cr_data.py').read_text()==test_before.replace('configs/mykonos.yaml','configs/cr_reference_library.yaml')

independent=json.loads((OUT/'protected-before.json').read_text())['files']
assert all(sha(ROOT/path)==pin for path,pin in independent.items())
preserved=json.loads((WORK/'preserved-before.json').read_text())
assert sha(WORK/'preserved-before.json')==freeze['preserved_inventory_sha256']
assert all(sha(ROOT/path)==pin for path,pin in preserved.items())
assert len(preserved)==308

test_source=(ROOT/'tests/test_legacy_mhd_retirement.py').read_text()
test_tree=ast.parse(test_source)
test_names=[node.name for node in test_tree.body if isinstance(node,ast.FunctionDef) and node.name.startswith('test_')]
assert len(test_names)==4
assert '"-I"' in test_source and 'MetaPathFinder' in test_source
assert 'importlib.util.find_spec' in test_source and 'run_from_config' in test_source
assert 'source' not in functions and 'run_from_config' not in functions
for module,symbols in {'flash_surface':['FlashSurfaceHistory','load_surface_history'],
        'heating_benchmarks':['ConstantConductivityHalfSpace'],
        'hall_induction_benchmark':['PeriodicHallInduction'],
        'thermal_release':['BindingReservoir','release_step']}.items():
    target_tree=ast.parse((ROOT/'mhdlab'/f'{module}.py').read_bytes())
    names={node.name for node in target_tree.body if isinstance(node,(ast.FunctionDef,ast.ClassDef))}
    assert set(symbols)<=names,(module,symbols)
    assert f'mhdlab.{module}' in test_source
initial_test=(WORK/'test_legacy_mhd_retirement.initial.py').read_text()
assert test_source==initial_test.replace('HalfSpaceHeating','ConstantConductivityHalfSpace')

record={'schema':'mitl-legacy-mhd-retirement-independent-review-v1',
        'status':'accepted_read_only_review_no_substantive_issue_identified',
        'source_freeze_sha256':sha(WORK/'source-freeze.json'),
        'owned_after_sha256':{path:value['after_sha256'] for path,value in freeze['owned_file_hashes'].items()},
        'checks':{'twelve_removed_paths_absent':len(freeze['removed_paths'])==12,
                  'only_two_runner_helpers_and_minimal_imports':True,
                  'independent_original_helper_AST_and_source_unchanged':True,
                  'package_only_exports_load_config_without_runner_import':True,
                  'cr_subsection_exact_original_bytes_and_same_directory':True,
                  'cr_test_only_filename_changed':True,
                  'independent_protected_103_files_unchanged':True,
                  'author_protected_308_files_unchanged':True,
                  'focused_tests_cover_clean_imports_absence_and_reference_interfaces':True,
                  'preexecution_test_symbol_repair_is_only_class_name':True},
        'focused_test_names':test_names,'findings':[],
        'limits':['Static file/AST/hash review only; no target modules imported and no tests or models executed.',
                  'Runtime verification remains the parent-owned focused test command; no science suite replay is requested.'],
        'target_imports':0,'tests_executed':0,'model_or_source_evolution_calls':0,
        'public_files_edited':0,'commits_or_status_mutations':0}
(OUT/'completion.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':record['status'],'completion_sha256':sha(OUT/'completion.json'),
                  'protected_independent':len(independent),'protected_author':len(preserved),
                  'focused_test_functions':len(test_names)},sort_keys=True))
