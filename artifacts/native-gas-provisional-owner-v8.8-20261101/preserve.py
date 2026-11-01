from pathlib import Path
import json,hashlib,shutil,difflib
RT=Path('/home/treefor/projects/runtimes');N=RT/'pic-reaction-transactions-20261029/native-full-time-v8.8-provisional-target-owner-20261101';B=N.parent/'native-full-time-v8.7-reservoir-stock-20261101'
C=RT/'native-gas-common-stage-v1.2-provisional-owner-20261101';CB=RT/'native-gas-common-stage-v1.1-reference-accounting-20261101'
O=Path('/home/treefor/projects/.workspace/outputs/native-gas-common-stage-v1.2-provisional-owner-20261101');A=Path('/home/treefor/projects/Power-Flow/artifacts/native-gas-provisional-owner-v8.8-20261101')
BACK=Path('/home/treefor/Documents/Codex/reconstructed-source-20261029/native-gas-provisional-owner-v8.8-20261101')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def copy(p,q):q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
m=json.loads((N/'SOURCE-SHA256.json').read_text());bm=json.loads((B/'SOURCE-SHA256.json').read_text());assert m['source_identity']=='9538ae6d7f1cffb4c0f446eba7df0449234a34eeb14169a1127d18e52853d6a7'
assert all(sha(N/f)==h for f,h in m['files'].items()) and all(sha(B/f)==h for f,h in bm['files'].items())
receipts={}
for label in ['cpu-c12','asan-c3']:
 p=O/label/'receipt.json';r=json.loads(p.read_text());assert r['passed'] and r['native_source']==m['source_identity'] and all(x['rc']==0 for x in r['runs'])
 assert all(sha(p)==h for p,h in r['source_inputs'].items()) and all(sha(p)==h for p,h in r['compiler_headers'].items())
 assert sha(O/label/'test_common')==r['binary_sha256']
 tag='sanitized-' if r['mode']=='sanitized' else '';assert sha(N/'build'/(tag+'libpic_native_v8.a'))==r['native_archive_sha256']
 receipts[label]={'path':str(p),'sha256':sha(p),'job':r['job'],'controls':len(r['results']),'archive':r['native_archive_sha256']}
p=O/'allocation-sweep-c3/receipt.json';r=json.loads(p.read_text());assert r['passed'] and len(r['checks'])==r['measured_allocations']==1243
assert [v['offset']for v in r['checks']]==list(range(1243)) and all(sha(p)==h for p,h in r['source_pins'].items())
receipts['sweep']={'path':str(p),'sha256':sha(p),'job':r['job'],'persistent_failure_positions':1243,'retry_after_each_failure':True}
assert not A.exists() and not BACK.exists();A.mkdir()
for top in ['source','forms','inputs','tests']:
 for p in (N/top).rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts:copy(p,A/'native'/p.relative_to(N))
for name in ['SOURCE-SHA256.json','FORK-BASE.json','tooling.py','build-components.py']:
 if (N/name).exists():copy(N/name,A/'native'/name)
for top in ['include','tests','tools']:
 for p in (C/top).rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts:copy(p,A/'bridge'/p.relative_to(C))
for name in ['DEPENDENCIES.json','NATIVE-BUILD.json','PROVISIONAL-OWNER-INTERFACE.txt']:
 copy(C/name,A/'bridge'/name)
def diff(old,new,files):
 out=[]
 for f in sorted(files):
  a=(old/f).read_text().splitlines(True) if (old/f).exists() else [];b=(new/f).read_text().splitlines(True)if(new/f).exists()else[]
  if a!=b:out.extend(difflib.unified_diff(a,b,fromfile='a/'+f if a else '/dev/null',tofile='b/'+f if b else '/dev/null'))
 return ''.join(out)
(A/'adoption-native-from-v8.7.diff').write_text(diff(B,N,set(bm['files'])|set(m['files'])))
bridge_files={str(p.relative_to(x))for x in[C,CB]for t in['include','tests','tools']for p in(x/t).rglob('*')if p.is_file()and'__pycache__'not in p.parts}
(A/'adoption-bridge-from-v1.1.diff').write_text(diff(CB,C,bridge_files));copy(B/'SOURCE-SHA256.json',A/'base-native-v8.7.json')
manifest={}
for root in[O,Path('/home/treefor/projects/.workspace/outputs/native-v8.8-provisional-target-owner-20261101')]:
 for p in root.rglob('*'):
  if not p.is_file() or p.name=='preserve-provisional-draft.py':continue
  key=root.name+'/'+str(p.relative_to(root));manifest[key]={'bytes':p.stat().st_size,'sha256':sha(p)}
  is_sweep='allocation-sweep-c3' in p.parts
  selected=not is_sweep or p.name=='receipt.json' or any(v in p.parts for v in['allocation-seal-0','allocation-seal-173','allocation-seal-256','allocation-seal-1242'])
  if selected and p.suffix in['.json','.log','.txt','.gdb']and p.stat().st_size<=524288 and p.name not in['processes-before.txt','processes-after.txt']:copy(p,A/'evidence'/key)
for job in [532,535,541,544,551,553,554,555,557,562,567,568,570,573,579,586,587,589]:
 p=Path('/home/treefor/.local/state/simqueue/logs')/(str(job).zfill(6)+'.log')
 if p.exists():copy(p,A/'evidence/scheduler'/p.name)
(A/'evidence-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
q={'schema':'native-gas-provisional-owner-cpu-v8.8','native_source':m['source_identity'],'receipts':receipts,'production_qualified':False,'physical_model_qualified':False,'device_coupled_path_qualified':False,'common_checkpoint_restore_qualified':False,'coupled_bookkeeping_scope':'actual populated host-only two-stage native/gas inventories; no births/sinks/cell-crossing motion','raw_physical_exact_balance_claim':False,'local_reaction_EM_continuity_qualified':False,'accounting':'actual stock closure after explicit guarded numerical event defects and species-reference/process-Q operand conversion','policy':'existing unknown-rate and minimum-network guards retained; component eventful tests use explicit diagnostic network scope','capability':'private genuine stage/final gas witnesses, complete two-stage transcript including empty stages, accepted-pair continuation; ordinary public and diagnostic preparation cannot open coupled gate','continuation':'live same-owner accepted pair only; failure permits retry; successful publication invalidates old pair; unattested restored pairs refuse','resources':'frozen userspace, one scheduled host CPU per job, explicit memory; no GPU work and no performance claim','local_only':'binaries and complete large per-offset fault-sweep outputs remain at hashed evidence paths; no full-machine process listings published'}
(A/'QUALIFICATION.json').write_text(json.dumps(q,indent=2)+'\n')
(A/'HANDOFF.txt').write_text('''This release qualifies a bounded CPU native/gas ownership transaction, not a full water-breakdown simulation.

The concrete bridge quarantines eventful native preparations until independent actual carrier, pending neutral-parcel, packet-internal, gas correction and reference accounting closes. Both stages, including empty ones, are required. The final accepted native and gas states publish only after all allocations/hashing and genuine capability preflight. Live accepted-pair continuation supports retry and subsequent steps; restored/unattested pairs are refused.

Physical raw sums are reported separately. Exact accounted closure includes bounded event numerical defects and explicit species-reference/process-Q conversion operands; it does not assert exact raw physical conservation or establish missing reaction rates. Nucleus accounting is H/O only. Global charge closure does not prove nodal charge continuity: charged targets can be removed at their own positions while charged products are created at the projectile position. Actual Hall census/charge admission must refuse unsupported relocation until a genuine reaction current and field-work owner are qualified. Target distributions remain the declared moment-based model. The strict ordinary no-target-owner path, unknown-rate/scope/resource guards and CPU-only refusal remain.

Qualification includes actual 1,686 stage-1 and 2,759 stage-2 events, 1,049 neutral reservoir releases, an actual second accepted common step, empty stage transcript, stale/foreign/destroyed/moved owner controls, unchanged accepted RNG/state on refusal, and no allocations in accepted publication. Persistent failure is injected at every one of 1,243 measured final-seal allocation positions; every failed attempt rolls back and then succeeds on a retry from the same accepted pair, producing the identical final common receipt as the uninterrupted control. ASan/UBSan repeats functional and selected failure controls. This is a final-seal fault sweep, not a claim that every earlier preparation allocation path is exhaustively tested.

Two inherited frozen-JSON hazards were reproduced: container destruction allocates a flattening stack; null-to-object indexing changes the tag before allocating storage. The native final seal and its JSON copy/build path now construct object storage first and drain guarded partial trees during unwinding. No vendored dependency was changed. An intermediate proof-output indexing typo and a test assertion for a nonexistent JSON key were corrected before final qualification; failed evidence is retained.

Next adoption: a separate authenticated prepared-motion snapshot and gas conservative-census transport must add births, ordinary sinks and cell crossings; current bridge refuses them. Genuine source/Hall/material owners, common checkpoint restore, complete chemistry network and physical convergence remain separate gates. The new RNG2/kind-32 dynamic-CDF release is not included here: this transaction release remains on its pinned original chemistry table/layout. No old production run was resumed.
''')
copy(Path(__file__),A/'preserve.py')
assert all(sha(A/'native'/f)==h for f,h in m['files'].items())
files={str(p.relative_to(A)):sha(p)for p in sorted(A.rglob('*'))if p.is_file()};(A/'SEAL.json').write_text(json.dumps({'schema':'native-gas-provisional-owner-v8.8-seal','source':m['source_identity'],'files':files},indent=2)+'\n')
shutil.copytree(A,BACK/'artifact')
for p in [N/'build/libpic_native_v8.a',N/'build/sanitized-libpic_native_v8.a',N/'build/kinetic-boundary.o',N/'build/sanitized-kinetic-boundary.o',O/'cpu-c12/test_common',O/'asan-c3/test_common']:
 copy(p,BACK/'local-binaries'/p.parent.name/p.name)
assert all(sha(p)==sha(BACK/'artifact'/p.relative_to(A))for p in A.rglob('*')if p.is_file())
(BACK/'BACKUP-SHA256.json').write_text(json.dumps({str(p.relative_to(BACK)):sha(p)for p in sorted(BACK.rglob('*'))if p.is_file()},indent=2)+'\n')
print(A);print('seal',sha(A/'SEAL.json'));print('backup',BACK)
