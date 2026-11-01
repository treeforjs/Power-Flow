from pathlib import Path
import hashlib,json,shutil,difflib,importlib.util
R=Path(__file__).resolve().parent;RT=R.parent
N=RT/'pic-reaction-transactions-20261029/native-full-time-v8.5-host-target-feedback-20261101'
B=RT/'pic-reaction-transactions-20261029/native-full-time-v8.4-water-scope-20261101'
O=Path('/home/treefor/projects/.workspace/outputs');A=Path('/home/treefor/projects/Power-Flow/artifacts/native-host-target-feedback-v8.5-20261101')
assert not A.exists()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
m=json.loads((N/'SOURCE-SHA256.json').read_text());base=json.loads((N/'SCOPE-BASE.json').read_text())
assert all(sha(N/n)==v for n,v in m['files'].items()) and all(sha(B/n)==v for n,v in base['files'].items())
assert m['source_identity']=='d996c286997cdfa4ffdf485da77433083dff1c3f0aacfeeb1058585b32971b89'
spec=importlib.util.spec_from_file_location("native_tooling",N/"tooling.py");tooling=importlib.util.module_from_spec(spec);spec.loader.exec_module(tooling)
current_deps=tooling.sealed_inputs()
receipts={}
for kind,path in [('cpu',O/'native-gas-target-feedback-v1-20261101/h5/receipt.json'),('native_sanitized',O/'native-v8.5-host-target-feedback-20261101/s2/receipt.json'),('all_sanitized',O/'native-gas-target-feedback-v1-20261101/s1/receipt.json')]:
 q=json.loads(path.read_text());assert q['passed'];assert all(x['rc']==0 for x in q['runs'])
 for name,digest in q['source_inputs'].items():
  p=Path(name) if Path(name).is_absolute() else N/name
  assert sha(p)==digest,(kind,name)
 assert q['dependencies']==current_deps
 if kind=='native_sanitized':
  for name,digest in q['binary_sha256'].items():assert sha(N/'build'/('sanitized-'+name))==digest
 else:assert sha(path.parent/'test_targets')==q['binary_sha256']
 receipts[kind]={'path':str(path),'sha256':sha(path),'job':q['job']}
A.mkdir()
def copy(p,q):q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
for top in ['source','inputs','forms','tests']:
 for p in (N/top).rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts:copy(p,A/'native'/p.relative_to(N))
for name in ['SOURCE-SHA256.json','SCOPE-BASE.json','tooling.py','build-components.py','qualify-native-feedback.py','IMPLEMENTATION-NOTES.txt']:
 copy(N/name,A/'native'/name)
for top in ['include','tests','tools']:
 for p in (R/top).rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts:copy(p,A/'provider'/p.relative_to(R))
for name in ['README.txt','COMMON-STAGE-PLAN.txt','DEPENDENCIES.json','preserve.py']:copy(R/name,A/'provider'/name)
dep=json.loads((R/'DEPENDENCIES.json').read_text());copy(Path(dep['binding_header']),A/'provider/pinned-dependency/native_gas_stage.hpp')
diff=[]
for name in sorted(set(base['files'])|set(m['files'])):
 a=(B/name).read_text().splitlines(True) if name in base['files'] else []
 b=(N/name).read_text().splitlines(True) if name in m['files'] else []
 if a!=b:diff.extend(difflib.unified_diff(a,b,fromfile='a/'+name if a else '/dev/null',tofile='b/'+name if b else '/dev/null'))
(A/'native-feedback-adoption.diff').write_text(''.join(diff))
copy(N/'pre-callback-owner-guard/PREIMAGE.json',A/'development-callback-guard-preimage.json')
evidence={}
for root,labels in [(O/'native-gas-target-feedback-v1-20261101',['h1','h2','h3','h5','s1']), (O/'native-v8.5-host-target-feedback-20261101',['build-normal','build-normal-pre-callback-owner-guard','build-sanitized','build-sanitized-pre-callback-owner-guard','s1','s2'])]:
 for label in labels:
  for p in (root/label).rglob('*'):
   if not p.is_file():continue
   name=root.name+'/'+str(p.relative_to(root));evidence[name]={'bytes':p.stat().st_size,'sha256':sha(p)}
   if p.suffix in ['.json','.log','.txt'] and p.stat().st_size<=262144 and p.name not in ['processes-before.txt','processes-after.txt']:
    copy(p,A/'evidence'/name)
(A/'evidence-manifest.json').write_text(json.dumps(evidence,indent=2)+'\n')
q={'schema':'native-host-target-feedback-v8.5-cpu-reference','native_source':m['source_identity'],'normal_archive_sha256':sha(N/'build/libpic_native_v8.a'),'sanitized_archive_sha256':sha(N/'build/sanitized-libpic_native_v8.a'),'receipts':receipts,'mode':'CPU diagnostic component','production_qualified':False,'full_system_energy_closure':False,'device_feedback_qualified':False,'physical_model_validated':False,'maxwellian_moment_approximation':True,'strict_target_energy_owner_guard':'closed','configuration_checkpoint_schema':'v8.5-host-target-feedback-v1; old v8.4 refused','claims':['actual immutable native inputs and final target-byte binding','actual gas lifetime/revision/step binding','zero-correction physical records and ledgers equal','controlled nonzero moment injection to stage2','before/at-anchor restored feedback roots and tamper refusal','exact catalogue energy convention registration guard'],'excluded_evidence':'Full-machine process listings, binaries and large checkpoint numerical payloads remain local at hashed evidence paths. No performance claim.'}
(A/'QUALIFICATION.json').write_text(json.dumps(q,indent=2)+'\n')
assert all(sha(N/n)==v and sha(A/'native'/n)==v for n,v in m['files'].items())
files={str(p.relative_to(A)):sha(p) for p in sorted(A.rglob('*')) if p.is_file()}
(A/'SEAL.json').write_text(json.dumps({'schema':'native-host-target-feedback-v8.5-seal','files':files,'native_source':m['source_identity']},indent=2)+'\n')
print(A);print(sha(A/'SEAL.json'))
