from pathlib import Path
import hashlib,json,shutil
R=Path('/home/treefor/projects/runtimes/gas-target-owner-v1.6-reservoir-release-20261101')
OUT=Path('/home/treefor/projects/.workspace/outputs')/R.name
REPO=Path('/home/treefor/projects/Power-Flow')
ART=REPO/'artifacts/gas-reservoir-release-v1.6-20261101'
BACK=Path('/home/treefor/Documents/Codex/reconstructed-source-20261029')/ART.name
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2)+'\n')
assert not ART.exists() and not BACK.exists()
receipts={}
for label,job in [('c2',478),('s1',479),('g1',480)]:
 rec=json.loads((OUT/label/'receipt.json').read_text());assert int(rec['job'])==job and rec['passed'] and rec['frozen_root_verified'] and rec['inputs_unchanged'] and rec['binary_unchanged']
 for p,h in rec['input_hashes'].items():assert sha(Path(p))==h,p
 receipts[label]=rec
assert receipts['g1']['state_parity']
assert 'configuration identity' in (OUT/'c2/legacy-release-checkpoint-refusal.log').read_text()
ART.mkdir()
for name in ['include','src','device','tests','tools','baseline','reference-v1.1','integration']:
 shutil.copytree(R/name,ART/name,ignore=shutil.ignore_patterns('__pycache__'))
for name in ['BASELINE-v1.5.json','RELEASE-CONTRACT.txt']:shutil.copy2(R/name,ART/name)
shutil.copy2(Path(__file__),ART/'preserve-release.py')
for label in ['c1','c2','s1','g1']:
 d=ART/'evidence'/label;d.mkdir(parents=True)
 for p in (OUT/label).iterdir():
  if p.suffix in ('.json','.log','.checkpoint'):shutil.copy2(p,d/p.name)
for job in [472,478,479,480]:shutil.copy2(Path('/home/treefor/.local/state/simqueue/logs')/f'{job:06}.log',ART/'evidence'/f'{job:06}.log')
checks={label:{p.stem:p.read_text().splitlines()[-1] for p in (OUT/label).glob('*.log') if p.stem in ('host','device','verified')} for label in receipts}
probe=Path('/home/treefor/projects/.workspace/outputs/native-gas-common-stage-v1.1-reference-accounting-20261101/cpu-c1/controls/result.json')
actual=json.loads(probe.read_text());assert actual['passed'] and all(s['energy_exact_zero'] and s['momentum_exact_zero'] and s['reference_operands_independently_verified'] for s in actual['stages'])
shutil.copy2(probe,ART/'evidence/actual-native-two-stage-result.json')
write(ART/'RESULTS.json',{'jobs':{'cpu':478,'sanitized':479,'gpu':480},'checks':checks,'actual_v15_release_mirror_preimage_reproduced':True,'old_checkpoint_refused_by_configuration_identity':True,'device_state_parity':True,'native_two_stage_probe':{'job':476,'result_path':str(probe),'sha256':sha(probe),'source':actual['native_source'],'common_accepted_publication':False},'semantic_change':'pending reservoir kinetic/momentum mirror subtracts signed releases','internal_energy_ownership_changed':False,'legacy_checkpoint_migration':False,'physical_breakdown_qualified':False,'speedup_qualified':False})
(ART/'README.txt').write_text('Gas pending-reservoir release correction\n\nActual native stage-2 releases exposed the inherited v1.5 defect: G_RES_E/P kept products that had already become carriers. The shared host/device core now subtracts actual signed release kinetic energy and momentum. Internal energy remains separately owned.\n\nQualification covers partial and split release, empty/pending stock, full drain, rollback, inverse signed bookkeeping, non-release isolation and checkpoint continuation. Actual v1.5 executables reproduce the stale mirror. New executables refuse actual old checkpoints by the changed configuration identity; no migration is implemented. CPU, ASan/UBSan and scheduled GPU/device-verified paths pass; state digits agree between CPU and GPU within this version. Old/new state digests intentionally differ.\n\nA separate real native/gas two-stage diagnostic probe closes actual carrier+stored-parcel charge, H/O inventory, energy and momentum exactly with independently derived chemical-reference conversion operands. Its result is preserved here; it does not establish joint accepted publication, a complete reaction network, a physical breakdown run or a speedup. Frozen input and allocation receipts are included.\n')
sources={str(p.relative_to(ART)):sha(p) for name in ['include','src','device','tests','tools','baseline','reference-v1.1','integration'] for p in sorted((ART/name).rglob('*')) if p.is_file()}
write(ART/'SOURCE.json',{'identity':hashlib.sha256((json.dumps(sources,sort_keys=True,indent=2)+'\n').encode()).hexdigest(),'files':sources})
files={str(p.relative_to(ART)):sha(p) for p in sorted(ART.rglob('*')) if p.is_file()}
write(ART/'SEAL.json',{'identity':hashlib.sha256((json.dumps(files,sort_keys=True,indent=2)+'\n').encode()).hexdigest(),'files':files})
shutil.copytree(ART,BACK)
for p in ART.rglob('*'):
 if p.is_file():assert p.read_bytes()==(BACK/p.relative_to(ART)).read_bytes()
print(json.dumps({'artifact':str(ART),'backup':str(BACK),'seal_sha256':sha(ART/'SEAL.json'),'files':len(files)+1,'checks':checks},indent=2))
