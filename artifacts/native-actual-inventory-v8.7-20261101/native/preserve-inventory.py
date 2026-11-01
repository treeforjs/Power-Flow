from pathlib import Path
import json,hashlib,shutil,difflib,importlib.util
N=Path(__file__).resolve().parent;B=N.parent/'native-full-time-v8.5-host-target-feedback-20261101';V6=N.parent/'native-full-time-v8.6-candidate-population-20261101'
O=Path('/home/treefor/projects/.workspace/outputs');A=Path('/home/treefor/projects/Power-Flow/artifacts/native-actual-inventory-v8.7-20261101')
assert not A.exists()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
m=json.loads((N/'SOURCE-SHA256.json').read_text());bm=json.loads((B/'SOURCE-SHA256.json').read_text());assert m['source_identity']=='d795695aca20ddef0efadedaebea8baa6c994e819514c0258bbac290fb66103d'
assert all(sha(N/n)==v for n,v in m['files'].items()) and all(sha(B/n)==v for n,v in bm['files'].items())
spec=importlib.util.spec_from_file_location('nt',N/'tooling.py');nt=importlib.util.module_from_spec(spec);spec.loader.exec_module(nt);deps=nt.sealed_inputs()
receipts={}
for mode,label in [('normal','cpu-c2'),('sanitized','asan-c1')]:
 p=O/'native-v8.7-reservoir-stock-20261101'/label/'receipt.json';r=json.loads(p.read_text());assert r['passed'] and all(q['rc']==0 for q in r['runs']) and r['dependencies']==deps
 assert all(sha(N/f)==h for f,h in r['source_inputs'].items())
 tag='sanitized-' if mode=='sanitized' else ''
 assert all(sha(N/'build'/(tag+f))==h for f,h in r['binary_sha256'].items())
 receipts[mode]={'path':str(p),'sha256':sha(p),'job':r['job'],'stock_checks':r['stock_controls']['checks'],'population_checks':r['population_controls']['checks']}
A.mkdir()
def copy(p,q):q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
for top in ['source','inputs','forms','tests']:
 for p in (N/top).rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts:copy(p,A/'native'/p.relative_to(N))
for f in ['SOURCE-SHA256.json','POPULATION-BASE.json','STOCK-BASE.json','tooling.py','build-components.py','qualify-candidate-population.py','qualify-reservoir-stocks.py','PRESERVATION-NOTES.txt','preserve-inventory.py']:copy(N/f,A/'native'/f)
diff=[]
for name in sorted(set(bm['files'])|set(m['files'])):
 old=(B/name).read_text().splitlines(True) if name in bm['files'] else []
 new=(N/name).read_text().splitlines(True) if name in m['files'] else []
 if old!=new:diff.extend(difflib.unified_diff(old,new,fromfile='a/'+name if old else '/dev/null',tofile='b/'+name if new else '/dev/null'))
(A/'adoption-from-v8.5.diff').write_text(''.join(diff));copy(B/'SOURCE-SHA256.json',A/'base-v8.5-source.json');copy(V6/'SOURCE-SHA256.json',A/'precursor-v8.6-source.json')
evidence={}
for root,labels in [(O/'native-v8.6-candidate-population-20261101',['cpu-c1','cpu-c2','cpu-c3','build-normal-job450-failed-fixture','build-normal-job452-fixture-flight']), (O/'native-v8.7-reservoir-stock-20261101',['cpu-c1','cpu-c2','asan-c1','build-normal','build-sanitized','build-normal-job457-namespace-typo'])]:
 for label in labels:
  for p in (root/label).rglob('*'):
   if not p.is_file():continue
   key=root.name+'/'+str(p.relative_to(root));evidence[key]={'bytes':p.stat().st_size,'sha256':sha(p)}
   if p.suffix in ('.json','.log','.txt') and p.stat().st_size<=262144 and p.name not in ('processes-before.txt','processes-after.txt'):copy(p,A/'evidence'/key)
(A/'evidence-manifest.json').write_text(json.dumps(evidence,indent=2)+'\n')
q={'schema':'native-v8.7-actual-inventory-cpu-reference','source':m['source_identity'],'normal_archive':sha(N/'build/libpic_native_v8.a'),'sanitized_archive':sha(N/'build/sanitized-libpic_native_v8.a'),'receipts':receipts,'production_qualified':False,'strict_energy_owner_gate_open':False,'full_system_closure_qualified':False,'device_view_qualified':False,'physical_model_qualified':False,'cumulative_nonzero_internal_sink_fixture':False,'packet_U_fixture':'manufactured authored native birth, not coupled source','local_only':'binaries, large checkpoint payloads and full machine process listings; originals remain at evidence paths with hashes; no performance claim'}
(A/'QUALIFICATION.json').write_text(json.dumps(q,indent=2)+'\n')
assert all(sha(N/f)==h and sha(A/'native'/f)==h for f,h in m['files'].items())
files={str(p.relative_to(A)):sha(p) for p in sorted(A.rglob('*')) if p.is_file()};(A/'SEAL.json').write_text(json.dumps({'schema':'native-v8.7-actual-inventory-seal','files':files,'source':m['source_identity']},indent=2)+'\n')
print(A);print(sha(A/'SEAL.json'))
