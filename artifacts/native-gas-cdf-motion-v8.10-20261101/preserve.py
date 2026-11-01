from pathlib import Path
import json,hashlib,shutil,difflib
RT=Path('/home/treefor/projects/runtimes');N=RT/'pic-reaction-transactions-20261029/native-full-time-v8.10-cdf-motion-20261101';B=N.parent/'native-full-time-v8.9-motion-inventory-20261101'
C=RT/'native-gas-common-stage-v1.4-cdf-motion-20261101';CB=RT/'native-gas-common-stage-v1.3.1-motion-unwind-20261101';IM=RT/'native-gas-adapter-v1.3-cdf-20261101'
O=Path('/home/treefor/projects/.workspace/outputs/native-gas-common-stage-v1.4-cdf-motion-20261101');NO=O.parent/'native-v8.10-cdf-motion-20261101'
A=Path('/home/treefor/projects/Power-Flow/artifacts/native-gas-cdf-motion-v8.10-20261101');BACK=Path('/home/treefor/Documents/Codex/reconstructed-source-20261029/native-gas-cdf-motion-v8.10-20261101')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  while b:=f.read(1<<20):h.update(b)
 return h.hexdigest()
def copy(p,q):q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
m=json.loads((N/'SOURCE-SHA256.json').read_text());bm=json.loads((B/'SOURCE-SHA256.json').read_text())
assert m['source_identity']=='7888c6beb64f40a98cefeba9712267be43e15983a4da7298e7bfafb59b36362e'
assert all(sha(N/f)==h for f,h in m['files'].items()) and all(sha(B/f)==h for f,h in bm['files'].items())
receipts={};cpu=None
for label in ['cpu-c2','sanitized-final']:
 p=O/label/'receipt.json';r=json.loads(p.read_text());assert r['passed'] and r['native_source']==m['source_identity'] and all(x['rc']==0 for x in r['runs'])
 assert all(sha(p)==h for p,h in r['source_inputs'].items()) and all(sha(p)==h for p,h in r['compiler_headers'].items())
 assert sha(O/label/'test_common')==r['binary_sha256'];tag='sanitized-'if r['mode']=='sanitized'else''
 assert sha(N/'build'/(tag+'libpic_native_v8.a'))==r['native_archive_sha256']
 cdf=json.loads((O/label/'actual-cdf-common/result.json').read_text());assert cdf['passed'] and cdf['events']==19 and cdf['releases']==38 and len(cdf['refusals'])==2
 receipts[label]={'path':str(p),'sha256':sha(p),'job':r['job'],'controls':len(r['results']),'archive':r['native_archive_sha256'],'cdf':{k:cdf[k]for k in ['checks','events','releases','refusals']}}
 if r['mode']=='normal':cpu=r
retry={}
for label,phase,count,control in [('seal-sweep-final','seal',1472,'normal'),('motion-sweep-final','motion',1263,'motion-empty-flows')]:
 p=O/label/'receipt.json';r=json.loads(p.read_text());assert r['passed'] and len(r['checks'])==r['measured_allocations']==count
 assert [x['offset']for x in r['checks']]==list(range(count)) and all(sha(p)==h for p,h in r['source_pins'].items())
 baseline=cpu['results'][control]['record']['sha256'];seen=set()
 for x in r['checks']:
  f=O/label/('allocation-'+phase+'-'+str(x['offset']))/'result.json';assert sha(f)==x['result_sha256'];j=json.loads(f.read_text());assert j['retry_accepted'];seen.add(j['retry_receipt_sha256'])
 assert seen=={baseline},(phase,seen,baseline)
 receipts[label]={'path':str(p),'sha256':sha(p),'job':r['job'],'persistent_failure_positions':count,'retry_equals_uninterrupted':True};retry[phase]={'positions':count,'retry_common_receipt':baseline}
assert not A.exists() and not BACK.exists();A.mkdir()
for top in ['source','forms','inputs','tests']:
 for p in (N/top).rglob('*'):
  if p.is_file()and'__pycache__'not in p.parts:copy(p,A/'native'/p.relative_to(N))
for name in ['SOURCE-SHA256.json','COMBINED-FORK-BASE.json','COMBINED-IDENTITY-NOTES.txt','adoption-cdf.patch','tooling.py','build-components.py']:copy(N/name,A/'native'/name)
for top in ['include','tests','tools']:
 for p in(C/top).rglob('*'):
  if p.is_file()and'__pycache__'not in p.parts:copy(p,A/'bridge'/p.relative_to(C))
for name in ['DEPENDENCIES.json','NATIVE-BUILD.json','NATIVE-BUILD-SANITIZED.json','COMBINED-FORK-BASE.json']:copy(C/name,A/'bridge'/name)
for name in ['include/native_gas_stage.hpp','tests/test_kind32.cpp']:copy(IM/name,A/'immutable-adapter'/name)
def diff(old,new,files):
 out=[]
 for f in sorted(files):
  a=(old/f).read_text().splitlines(True)if(old/f).exists()else[];b=(new/f).read_text().splitlines(True)if(new/f).exists()else[]
  if a!=b:out.extend(difflib.unified_diff(a,b,fromfile='a/'+f if a else'/dev/null',tofile='b/'+f if b else'/dev/null'))
 return ''.join(out)
(A/'adoption-native-from-v8.9.diff').write_text(diff(B,N,set(bm['files'])|set(m['files'])))
files={str(p.relative_to(x))for x in[C,CB]for t in['include','tests','tools']for p in(x/t).rglob('*')if p.is_file()and'__pycache__'not in p.parts}
(A/'adoption-bridge-from-v1.3.1.diff').write_text(diff(CB,C,files));copy(B/'SOURCE-SHA256.json',A/'base-native-v8.9.json')
manifest={}
for root in[O,NO]:
 for p in root.rglob('*'):
  if not p.is_file()or p==Path(__file__)or p.name.startswith('core')or p.name in ['processes-before.txt','processes-after.txt']:continue
  key=root.name+'/'+str(p.relative_to(root));manifest[key]={'bytes':p.stat().st_size,'sha256':sha(p)}
  sweep=next((q for q in ['seal-sweep-final','motion-sweep-final']if q in p.parts),None)
  select=not sweep or p.name=='receipt.json' or any(q in p.parts for q in ['allocation-seal-0','allocation-seal-312','allocation-seal-1471','allocation-motion-0','allocation-motion-64','allocation-motion-1262'])
  if select and p.suffix in['.json','.log','.txt']and p.stat().st_size<=524288:copy(p,A/'evidence'/key)
for job in[662,668,674,675,676]:
 p=Path('/home/treefor/.local/state/simqueue/logs')/(str(job).zfill(6)+'.log');copy(p,A/'evidence/scheduler'/p.name)
(A/'evidence-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');(A/'RETRY-IDENTITY.json').write_text(json.dumps(retry,indent=2)+'\n')
normal=cpu['results']['normal'];(A/'TIMING-BASELINE.json').write_text(json.dumps({'scope':'single scheduled CPU under concurrent machine load; component timing, no exclusive speedup claim','carriers':normal['census_carriers'],'events':[x['events']for x in normal['record']['stages']],'phases':normal['phase_seconds'],'source_sha256':sha(O/'cpu-c2/normal/result.json'),'note':'motion_* are nested wall-only subphases of admission_accept'},indent=2)+'\n')
q={'schema':'native-gas-cdf-motion-cpu-v8.10','native_source':m['source_identity'],'receipts':receipts,'production_qualified':False,'device_coupled_path_qualified':False,'common_checkpoint_restore_qualified':False,'full_chemistry_network_qualified':False,'local_reaction_EM_continuity_qualified':False,'thermal_routing_of_gas_corrections_qualified':False,'physical_raw_exact_balance_claim':False,'scope':'RNG2/blob2 kind32 CDF adoption with the authentic motion/strict target-owner transaction; immutable retained stage stocks, populations, admitted exports, actual motion and exact per-exit gas transport receipts','diagnostic_CDF':'actual 19 events/38 H fragments with explicit diagnostic missing-channel admission; minimum-water scope and selected unsupported CDF energy independently refuse','accounting':'actual independent stored carrier/reservoir reference and kinetic operands; explicit numerical-defect/conversion ledgers; held conductor/exterior corrections','limitations':['source-conditioned H2 central-quantile approximation remains explicit; fixed-seed records change deliberately with RNG2','all missing process groups/channels and outside-support guards remain closed in strict mode','ordinary bridge refuses merge, returned-subcell and wall-response birth mappings','root must validate actual final native nodal charge against Hall endpoints and every nonzero conductor correction before publication','component fixtures do not qualify measured field/heating/source trajectories or source-held layers'],'preservation':'compiled binaries and complete per-offset outputs remain local with hashes; no unrelated process listings; 662 missing-MFEM-link harness failure retained'}
(A/'QUALIFICATION.json').write_text(json.dumps(q,indent=2)+'\n')
(A/'HANDOFF.txt').write_text('''This CPU release combines the qualified RNG2/dynamic-CDF kernel with the authentic native/gas motion transaction. It is not a full physical water-breakdown run. Numerical source identity and archive hashes are in QUALIFICATION.json; the source-conditioned H2 descriptor remains separately model-qualified only within its declared support. Fixed-seed trajectories intentionally differ from RNG1.

The actual CDF fixture admits 19 kind-32 reactions and releases 38 H fragments, with the strict target-energy ownership transaction in place. It explicitly admits otherwise missing encounter coverage for this diagnostic. Separate controls preserve minimum-water incompleteness refusal and low-energy selected CDF_SUPPORT whole-step rollback. No missing physics is silently assigned zero.

CoupledTrial keeps the same construction/call ordering as v1.3.1: stage1, actual native motion, gas birth/census/sink flows, stage2, complete all-owner preparation, final native/gas seal, all-owner preflight, no-fail publication. Immutable historical getters stage_stocks(1/2), stage_export(1/2), stage_population(1/2), accepted_motion_snapshot() and motion_transport_receipt() now retain actual bounded objects created by their owners. They confer no publication authority. Root must retain/preflight that exact CoupledTrial and apply the final Hall-charge and per-exit conductor-correction guards. Signed exact local exit shares may cancel in totals; use the retained individual receipt. Gas correction exits remain held energy/momentum, not implicitly Hall heat.

CPU 668 and ASan/UBSan 674 pass the populated CDF and inherited scope/target-energy/reference/stored-stock controls. Persistent failures cover all 1,263 measured nonempty-motion allocations and all 1,472 measured final-seal allocations (675/676), with identical uninterrupted receipt on every retry and unchanged accepted native/RNG/gas state on refusal. This is not a global all-code allocation-failure guarantee. Timings are concurrent-load component diagnostics, not an isolated speedup claim. Underlying combined device execution is not qualified by these CPU common-owner tests.

Build 662 successfully built the native/archive and common/CDF programs, but its small immutable-adapter test link omitted MFEM libraries. Only that harness link was fixed; 668 reused the same archive. The SANITIZED build-cache marker deliberately refuses reuse; actual sanitized archive pins are in its final qualification receipt. Frozen dependencies and all model/resource guards are unchanged.

A versioned paired-checkpoint successor is separate. This release still refuses unattested restored pairs. Actual root Hall/source/arrival ownership restoration and source-held electron-layer lifecycle are additional integration gates. No previous production run was resumed.
''')
copy(Path(__file__),A/'preserve.py');assert all(sha(A/'native'/f)==h for f,h in m['files'].items());files={str(p.relative_to(A)):sha(p)for p in sorted(A.rglob('*'))if p.is_file()};(A/'SEAL.json').write_text(json.dumps({'schema':'native-gas-cdf-motion-v8.10-seal','source':m['source_identity'],'files':files},indent=2)+'\n')
shutil.copytree(A,BACK/'artifact')
for p in[N/'build/libpic_native_v8.a',N/'build/sanitized-libpic_native_v8.a',N/'build/kinetic-boundary.o',N/'build/sanitized-kinetic-boundary.o',O/'cpu-c2/test_common',O/'cpu-c2/test_common_cdf',O/'sanitized-final/test_common',O/'sanitized-final/test_common_cdf']:copy(p,BACK/'local-binaries'/p.parent.name/p.name)
assert all(sha(p)==sha(BACK/'artifact'/p.relative_to(A))for p in A.rglob('*')if p.is_file());(BACK/'BACKUP-SHA256.json').write_text(json.dumps({str(p.relative_to(BACK)):sha(p)for p in sorted(BACK.rglob('*'))if p.is_file()},indent=2)+'\n');print(A);print('seal',sha(A/'SEAL.json'));print('backup',BACK)
