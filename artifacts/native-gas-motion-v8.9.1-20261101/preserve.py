from pathlib import Path
import json,hashlib,shutil,difflib
RT=Path('/home/treefor/projects/runtimes');N=RT/'pic-reaction-transactions-20261029/native-full-time-v8.9-motion-inventory-20261101';B=N.parent/'native-full-time-v8.8-provisional-target-owner-20261101'
C=RT/'native-gas-common-stage-v1.3.1-motion-unwind-20261101';CB=RT/'native-gas-common-stage-v1.2-provisional-owner-20261101'
O=Path('/home/treefor/projects/.workspace/outputs/native-gas-common-stage-v1.3.1-motion-unwind-20261101');OLD=O.parent/'native-gas-common-stage-v1.3-motion-20261101';NO=O.parent/'native-v8.9-motion-inventory-20261101'
A=Path('/home/treefor/projects/Power-Flow/artifacts/native-gas-motion-v8.9.1-20261101');BACK=Path('/home/treefor/Documents/Codex/reconstructed-source-20261029/native-gas-motion-v8.9.1-20261101')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  while b:=f.read(1<<20):h.update(b)
 return h.hexdigest()
def copy(p,q):q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
m=json.loads((N/'SOURCE-SHA256.json').read_text());bm=json.loads((B/'SOURCE-SHA256.json').read_text())
assert m['source_identity']=='8e8d9553ccac604e7e69ecc3354f124e71643a17e8fbc85b9883dd4f798897db'
assert all(sha(N/f)==h for f,h in m['files'].items()) and all(sha(B/f)==h for f,h in bm['files'].items())
receipts={};cpu=None
for label in ['cpu-final','sanitized-final']:
 p=O/label/'receipt.json';r=json.loads(p.read_text());assert r['passed'] and r['native_source']==m['source_identity'] and all(x['rc']==0 for x in r['runs'])
 assert all(sha(p)==h for p,h in r['source_inputs'].items()) and all(sha(p)==h for p,h in r['compiler_headers'].items())
 assert sha(O/label/'test_common')==r['binary_sha256'];tag='sanitized-'if r['mode']=='sanitized'else''
 assert sha(N/'build'/(tag+'libpic_native_v8.a'))==r['native_archive_sha256']
 receipts[label]={'path':str(p),'sha256':sha(p),'job':r['job'],'controls':len(r['results']),'archive':r['native_archive_sha256']}
 if r['mode']=='normal':cpu=r
retry={}
for label,phase,count,control in [('seal-sweep-final','seal',1463,'normal'),('motion-sweep-final','motion',1255,'motion-empty-flows')]:
 p=O/label/'receipt.json';r=json.loads(p.read_text());assert r['passed'] and len(r['checks'])==r['measured_allocations']==count
 assert [x['offset']for x in r['checks']]==list(range(count)) and all(sha(p)==h for p,h in r['source_pins'].items())
 baseline=cpu['results'][control]['record']['sha256'];seen=set()
 for x in r['checks']:
  f=O/label/('allocation-'+phase+'-'+str(x['offset']))/'result.json';assert sha(f)==x['result_sha256'];j=json.loads(f.read_text());assert j['retry_accepted'];seen.add(j['retry_receipt_sha256'])
 assert seen=={baseline},(phase,seen,baseline)
 receipts[label]={'path':str(p),'sha256':sha(p),'job':r['job'],'persistent_failure_positions':count,'retry_equals_uninterrupted':True}
 retry[phase]={'positions':count,'retry_common_receipt':baseline}
assert not A.exists() and not BACK.exists();A.mkdir()
for top in ['source','forms','inputs','tests']:
 for p in (N/top).rglob('*'):
  if p.is_file()and'__pycache__'not in p.parts:copy(p,A/'native'/p.relative_to(N))
for name in ['SOURCE-SHA256.json','MOTION-FORK-BASE.json','tooling.py','build-components.py']:copy(N/name,A/'native'/name)
for top in ['include','tests','tools']:
 for p in(C/top).rglob('*'):
  if p.is_file()and'__pycache__'not in p.parts:copy(p,A/'bridge'/p.relative_to(C))
for name in ['DEPENDENCIES.json','NATIVE-BUILD.json','NATIVE-BUILD-SANITIZED.json','MOTION-FORK-BASE.json','UNWIND-FORK-BASE.json']:copy(C/name,A/'bridge'/name)
def diff(old,new,files):
 out=[]
 for f in sorted(files):
  a=(old/f).read_text().splitlines(True)if(old/f).exists()else[];b=(new/f).read_text().splitlines(True)if(new/f).exists()else[]
  if a!=b:out.extend(difflib.unified_diff(a,b,fromfile='a/'+f if a else'/dev/null',tofile='b/'+f if b else'/dev/null'))
 return ''.join(out)
(A/'adoption-native-from-v8.8.diff').write_text(diff(B,N,set(bm['files'])|set(m['files'])))
files={str(p.relative_to(x))for x in[C,CB]for t in['include','tests','tools']for p in(x/t).rglob('*')if p.is_file()and'__pycache__'not in p.parts}
(A/'adoption-bridge-from-v1.2.diff').write_text(diff(CB,C,files));copy(B/'SOURCE-SHA256.json',A/'base-native-v8.8.json')
manifest={}
for root in[O,NO,OLD]:
 for p in root.rglob('*'):
  if not p.is_file()or p==Path(__file__) or p.name.startswith('core') or p.name in ['processes-before.txt','processes-after.txt']:continue
  key=root.name+'/'+str(p.relative_to(root));manifest[key]={'bytes':p.stat().st_size,'sha256':sha(p)}
  sweep=next((q for q in ['seal-sweep-final','motion-sweep-final']if q in p.parts),None)
  select=not sweep or p.name=='receipt.json' or any(q in p.parts for q in ['allocation-seal-0','allocation-seal-312','allocation-seal-1462','allocation-motion-0','allocation-motion-64','allocation-motion-1254'])
  if select and p.suffix in['.json','.log','.txt']and p.stat().st_size<=524288 and p.name not in['processes-before.txt','processes-after.txt']:copy(p,A/'evidence'/key)
for job in[622,626,629,630,635,636,637,638,650,653,655,656,657]:
 p=Path('/home/treefor/.local/state/simqueue/logs')/(str(job).zfill(6)+'.log');copy(p,A/'evidence/scheduler'/p.name)
(A/'evidence-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');(A/'RETRY-IDENTITY.json').write_text(json.dumps(retry,indent=2)+'\n')
normal=cpu['results']['normal'];timing={'source_result':str(O/'cpu-final/normal/result.json'),'source_sha256':sha(O/'cpu-final/normal/result.json'),'scope':'single scheduled core under concurrent machine load; component diagnostic, not production or exclusive speedup','carriers':normal['census_carriers'],'events':[x['events']for x in normal['record']['stages']],'phases':normal['phase_seconds'],'note':'motion_* are subphases of admission_accept; their CPU time is unmeasured/null; do not double-count'}
(A/'TIMING-BASELINE.json').write_text(json.dumps(timing,indent=2)+'\n')
q={'schema':'native-gas-authentic-motion-cpu-v8.9.1','native_source':m['source_identity'],'receipts':receipts,'production_qualified':False,'device_coupled_path_qualified':False,'common_checkpoint_restore_qualified':False,'full_chemistry_network_qualified':False,'local_reaction_EM_continuity_qualified':False,'thermal_routing_of_gas_corrections_qualified':False,'physical_raw_exact_balance_claim':False,'scope':'actual immutable prepared native before/start/after carrier and packet-energy snapshots, native births/first-hit sinks, gas census-conservative transport, complete stage1-motion-stage2 transcript and live final census','accounting':'H/O inventories; guarded numerical defects and explicit reference operand conversion; gas correction entry totals plus separate held conductor/exterior transfer totals','limitations':['same corrected old chemistry/RNG1 binding as v8.8; RNG2/CDF merge is a separate successor','ordinary bridge still refuses merge, returned-subcell and wall-response birth mappings','per-exit correction receipt is not retained by this bridge; root must refuse conductor sinks until per-exit routing/guard successor is adopted','actual final candidate must pass root Hall nodal charge comparison before full-owner publication','authored source inputs and reference motion fixtures are component controls, not measured physical heating trajectories'],'resource_policy':'frozen userspace, one scheduled CPU,3GiB normal/8GiB sanitizer; no GPU or uncontrolled timing','local_only':'compiled binaries and complete per-offset fault outputs are hash-listed but excluded from Git; no full-machine process command listings'}
(A/'QUALIFICATION.json').write_text(json.dumps(q,indent=2)+'\n')
(A/'HANDOFF.txt').write_text('''The v1.3 bridge failed a persistent allocation at motion offset 1163: a long-key JSON lookup allocated after constructing a populated temporary. v1.3.1 binds the array/element references before creating that object. The actual successful common receipt is unchanged, and native v8.9 source/archive is unchanged. The failed job 638 and trace 650 are retained; only the corrected bridge has passed the exhaustive motion sweep.

This release extends the bounded native/gas ownership transaction to authentic host motion inventories, ordinary source births, native first-hit conductor/exterior sinks, and conservative cell transport. It is not a full physical water-breakdown run.

Construction: build gas Configuration with Transport::census_conservative; set native coupled_target_contract_sha256=coupled_contract(base,options,gas), then enable host_target_feedback and compute its target contract. Construct CoupledTrial on the actual PrivateFullTimeStep and same live gas Owner/previous Continuation. Run stage1; prepare actual source births and native transport trial; admit genuine total-field motion; call coupled.accept_transport(move(motion)); run stage2; prepare all source/Hall future roots; seal common native/gas; obtain candidate_boundary() and final_candidate_population(); compare actual final census against Hall motion endpoints under the unchanged nodal bound; preflight every owner; perform only no-fail accepted publication. Root owns the field/source/material transaction.

motion_snapshot(token,byte_cap,carrier_cap) returns one immutable native-built snapshot per exact token. A copied/equal-attempt/moved-from/consumed snapshot is not authority. The optional zero carrier_cap means the common bound. Capacity checks are identical before and after cache hits, including zero-count stream and sink receipt cardinalities. The snapshot contains actual before-birth, after-birth start and motion-end carriers, packet-owned internal energy, authentic birth/sink receipts and native motion/state bindings. final_candidate_population() is retained from the actual stage2 preparation and verified against final native row identity, not reconstructed from event records.

Gas flow records use every actual positive binary carrier weight; no packet count is projected to one double. Its census transport handles exact evacuated-cell residue without erasing it. Nonzero correction exits remain held transfers, with conductor and exterior totals separate; neither is silently Hall heat and historical wall totals must not be counted twice. v1.3.1 does not retain individual correction exits: signed shares may cancel in aggregate, so the root integration currently must refuse conductor sinks. The separate v1.4 successor retains the genuine per-exit receipt for an exact local guard. Source birth/internal stock, incident sink kinetic energy and held gas corrections have distinct owners.

Global H/O/charge/energy accounting does not prove local Maxwell continuity. Charged targets can be removed at their position while products appear at the projectile position. Root must apply the actual final census-versus-Hall endpoint charge guard after stage2, before any publication; the existing Hall preflight alone does not inspect stage2 particles. Unsupported reaction relocation/current/work remains a refusal.

CPU and ASan/UBSan controls include real four-birth/two-sink motion, carried packet energy, crossing, separate nonzero conductor/exterior correction exits from an evacuated cell, authentic final population, cached cap refusal, equal-token/copy/move/consumed controls, same-digest gas mutation, rollback and two-step continuation. Persistent faults cover every measured final-seal allocation (1,463) and every allocation in a nonempty four-birth/two-sink motion (1,255). Every failed attempt preserves accepted native/RNG/gas and produces the exact uninterrupted common receipt on retry. This is not an exhaustive all-code OOM claim.

Component timings are in TIMING-BASELINE.json, with CPU/wall scope and subphase nesting explicit. Frozen source/RNG1 dependencies are unchanged; a corrected RNG2/dynamic-CDF combined successor is separate. Unknown-rate/minimum-network, target-owner, resource and CPU-only guards are preserved. Full field/source/material integration, network completeness, sheath/grid/material convergence, thermal routing and common checkpoint restore remain separate gates.

No unattested restored pair can acquire Continuation. A future trusted common loader must validate atomic multiowner manifests, actual component bytes/pins, clocks/roots and actual combined stock/transcript invariants before privately minting a fresh live accepted-pair capability; serialized hashes alone cannot authorize it. No previous production run was resumed.
''')
copy(Path(__file__),A/'preserve.py');assert all(sha(A/'native'/f)==h for f,h in m['files'].items())
files={str(p.relative_to(A)):sha(p)for p in sorted(A.rglob('*'))if p.is_file()};(A/'SEAL.json').write_text(json.dumps({'schema':'native-gas-motion-v8.9.1-seal','source':m['source_identity'],'files':files},indent=2)+'\n')
shutil.copytree(A,BACK/'artifact')
for p in[N/'build/libpic_native_v8.a',N/'build/sanitized-libpic_native_v8.a',N/'build/kinetic-boundary.o',N/'build/sanitized-kinetic-boundary.o',O/'cpu-final/test_common',O/'sanitized-final/test_common']:copy(p,BACK/'local-binaries'/p.parent.name/p.name)
assert all(sha(p)==sha(BACK/'artifact'/p.relative_to(A))for p in A.rglob('*')if p.is_file())
(BACK/'BACKUP-SHA256.json').write_text(json.dumps({str(p.relative_to(BACK)):sha(p)for p in sorted(BACK.rglob('*'))if p.is_file()},indent=2)+'\n')
print(A);print('seal',sha(A/'SEAL.json'));print('backup',BACK)
