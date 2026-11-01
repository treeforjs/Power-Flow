from pathlib import Path
import json,hashlib,shutil,difflib,sys
R=Path(__file__).resolve().parent;O=Path('/home/treefor/projects/.workspace/outputs/gas-target-owner-v1.4-binding-20261101')
cpu=sys.argv[1:]==['--cpu-milestone'];assert cpu or not sys.argv[1:]
A=Path('/home/treefor/projects/Power-Flow/artifacts')/('gas-candidate-binding-cpu-milestone-20261101' if cpu else 'gas-candidate-binding-v1.4-20261101')
assert not A.exists()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
base=json.loads((R/'BASELINE-v1.3.json').read_text());B=Path(base['base'])
assert all(sha(B/n)==v for n,v in base['source_files'].items())
files={str(p.relative_to(R)):sha(p) for top in ['include','src','device','tests','tools'] for p in sorted((R/top).rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
identity=hashlib.sha256(json.dumps(files,sort_keys=True,separators=(',',':')).encode()).hexdigest()
(R/'SOURCE.json').write_text(json.dumps({'identity':identity,'files':files,'base_file_pins':base['source_files']},indent=2)+'\n')
diff=[]
for n in sorted(set(base['source_files'])|set(files)):
    if not n.startswith(('include/','src/','device/','tests/','tools/')):continue
    a=(B/n).read_text().splitlines(True) if n in base['source_files'] else []
    b=(R/n).read_text().splitlines(True) if n in files else []
    if a!=b:diff.extend(difflib.unified_diff(a,b,fromfile='a/'+n if a else '/dev/null',tofile='b/'+n if b else '/dev/null'))
(R/'candidate-binding.diff').write_text(''.join(diff))
receipts={}
for label in ['h2','s2']+([] if cpu else ['g1']):
    q=json.loads((O/label/'receipt.json').read_text())
    assert q['passed'] and q['inputs_unchanged'] and q['binary_unchanged'] and q['frozen_root_verified']
    assert all(sha(Path(n))==v for n,v in q['input_hashes'].items())
    assert sha(O/label/'test_guards')==q['binary_sha256']
    assert all(c['rc']==0 for c in q['runs'])
    receipts[label]=q
    assert 'PASS 317 checks' in (O/label/'host.log').read_text()
    if label=='g1':
        assert q['state_parity']
        for lane in ['device','verified']:assert 'PASS 172 checks' in (O/label/(lane+'.log')).read_text()
for key in ['state ','motion_state ','publication_state ']:
    a=next(s for s in (O/'h2/host.log').read_text().splitlines() if s.startswith(key))
    b=next(s for s in (O.parent/'gas-target-owner-v1.3-publication-20261101/g1/host.log').read_text().splitlines() if s.startswith(key))
    assert a==b
A.mkdir(parents=True)
def copy(p,q):q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
for top in ['include','src','device','tests','tools','baseline','reference-v1.1']:
    for p in sorted((R/top).rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:copy(p,A/p.relative_to(R))
for n in ['BASELINE-v1.3.json','SOURCE.json','candidate-binding.diff','HANDOFF.txt','BINDING-CONTRACT.txt','PROPOSED-TARGET-FEEDBACK.txt','preserve-binding.py']:copy(R/n,A/n)
evidence={}
for label in ['h1','s1','h2','s2']+([] if cpu else ['g1']):
    for p in sorted((O/label).rglob('*')):
        if not p.is_file():continue
        n=str(p.relative_to(O));evidence[n]={'bytes':p.stat().st_size,'sha256':sha(p)}
        if p.suffix in ('.json','.log','.txt') and p.stat().st_size<=262144 and p.name not in ['processes-before.txt','processes-after.txt']:copy(p,A/'evidence'/n)
(A/'evidence-manifest.json').write_text(json.dumps(evidence,indent=2)+'\n')
q={'schema':'gas-candidate-binding-v1.4-component','source_identity':identity,'runtime':str(R),'evidence_root':str(O),
 'status':'CPU-and-sanitizer-milestone-GPU-pending' if cpu else 'qualified-read-only-binding-component',
 'cpu_job':390,'sanitizer_job':391,'gpu_job':392,'gpu_pending':cpu,'cpu_checks':317,'sanitizer_checks':317,
 'frozen_userspace_verified':True,'published_numerical_identity_unchanged_from_v1.3':True,
 'scope':'Read-only transient owner lifetime/generation/mutation-revision/step. No publication authority; content digest binding still required.',
 'error_guards':['failed device begin initialization requires abandon','failed partial key setup requires abandon','revision overflow refuses before mutation'],
 'preliminary_evidence':'381/382 (h1/s1) passed313 controls in host environment; outer launcher omitted. They are retained only as preliminary evidence, superseded by frozen390/391 final fixture317.',
 'device_fault_scope':'Test-only link wrapper injects exceptions before/after actual copy; no actual driver OOM recovery claim. '+('Pending GPU until392 passes.' if cpu else 'Shared-GPU392 passed before/after-copy refusal controls.'),
 'gpu_checks':None if cpu else {'host':317,'device':172,'verified':172},
 'untested_limit':'2^64 mutation-revision overflow is source guarded but not runtime exercised.',
 'production_qualified':False,'full_system_energy_closure':False,
 'local_only':'Test binaries and large numerical payloads remain at hashed runtime paths; full-machine process listings are excluded.'}
(A/'QUALIFICATION.json').write_text(json.dumps(q,indent=2)+'\n')
assert all(sha(R/n)==v and sha(A/n)==v for n,v in files.items())
sealed={str(p.relative_to(A)):sha(p) for p in sorted(A.rglob('*')) if p.is_file()}
(A/'SEAL.json').write_text(json.dumps({'schema':'gas-candidate-binding-v1.4-seal','source_identity':identity,'files':sealed},indent=2)+'\n')
print(A);print('seal',sha(A/'SEAL.json'))
