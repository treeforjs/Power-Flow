from pathlib import Path
import hashlib,json,shutil
R=Path('/home/treefor/projects/runtimes/pic-water-common-commit-v8.0.3-publication-20261101');O=Path('/home/treefor/projects/.workspace/outputs')/R.name
REPO=Path('/home/treefor/projects/Power-Flow');A=REPO/'artifacts/common-prepared-publication-v8.0.3-20261101';B=Path('/home/treefor/Documents/Codex/reconstructed-source-20261029')/A.name
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,sort_keys=True,indent=2)+'\n')
assert not A.exists() and not B.exists()
for label in ['c6','s1']:
 br=json.loads((O/('build-'+label)/'build-receipt.json').read_text());assert br['status']=='passed' and sha(Path(br['binary']))==br['binary_sha256']
 for path,h in br['own_sources_sha256'].items():assert sha(Path(path))==h,path
 for path,h in br['dependency_sha256'].items():
  p=Path(path)
  if p.is_relative_to('/usr') or p.is_relative_to('/etc'):p=R.parent/'arch-base-20261028/root'/p.relative_to('/')
  assert sha(p)==h,path
 for label2 in ([label+'-checks',label+'-publication-guard',label+'-replacement',label+'-lifetime'] if label=='c6' else [label+'-checks',label+'-publication-guard',label+'-unchanged-source-replacement',label+'-retained-source-lifetime']):
  rec=json.loads((O/label2/'run-receipt.json').read_text());assert rec['status']=='passed' and all(rec[k] for k in ['sealed_unchanged_after_run','own_sources_unchanged_after_run','binary_unchanged_after_run'])
assert json.loads((O/'s1-sanitized-result.json').read_text())['passed']
parity=json.loads((O/'c6-baseline-parity-p2/comparison.json').read_text());assert parity['passed'] and all(c['bit_identical'] for c in parity['comparisons'])
continuation=json.loads((O/'c6-continuation-results.json').read_text());assert continuation['passed']==continuation['total']==12
A.mkdir()
for name in ['source','tests','tools','vendor']:shutil.copytree(R/name,A/name,ignore=shutil.ignore_patterns('__pycache__'))
for name in ['run_v8.py','BASELINE-v8.0.2.json','run_v8.py.c1-preimage.txt','driver_v8_checks.c2-preimage.txt']:shutil.copy2(R/name,A/name)
shutil.copy2(Path(__file__),A/'preserve-publication.py')
shutil.copytree(O/'replacement-preimage-c4-source',A/'preimages/unchanged-source-c4')
old=R.parent/'pic-water-common-commit-v8.0.2-publication-preimage-20261101';oldout=O.parent/old.name
for name in ['source','tests','vendor']:shutil.copytree(old/name,A/'preimages/allocation-v8.0.2'/name,ignore=shutil.ignore_patterns('__pycache__'))
for name in ['run_v8.py','BASELINE.json']:
 if (old/name).is_file():shutil.copy2(old/name,A/'preimages/allocation-v8.0.2'/name)
pre=json.loads((oldout/'observed-p1-trace/preimage-receipt.json').read_text());assert pre['passed'] and pre['mixed_accepted_state_observed']
for name in ['observed-p1','observed-p1-trace','build-p1']:
 for p in (oldout/name).glob('*'):
  if p.suffix in ('.json','.log'):
   q=A/'evidence/allocation-preimage'/name/p.name;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
for p in O.rglob('*'):
 if not p.is_file() or 'replacement-preimage-c4-source' in p.parts or 'scratch' in p.parts or 'frames' in p.parts or 'checkpoints' in p.parts:continue
 allowed=p.suffix in ('.json','.log','.txt','.jsonl') or p.name in ('ledger.bin','deterministic.bin','step-digests.txt')
 if allowed and p.stat().st_size<=8*1024**2:
  q=A/'evidence/current'/p.relative_to(O);q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
for job in [460,464,466,483,486,490,496,501,502,507,508,510]:
 q=A/'evidence'/f'{job:06}.log';q.parent.mkdir(exist_ok=True);shutil.copy2(Path('/home/treefor/.local/state/simqueue/logs')/q.name,q)
checks={label:json.loads((O/(label+'-checks/checks-result.json')).read_text()) for label in ['c6','s1']}
(A/'README.txt').write_text('Prepared all-owner publication and source identity\n\nThe actual v8.0.2 coordinator could terminate after an allocation failure with Hall advanced and arrival unchanged (466). Arrival v1.2 and source v7.2 prepare all future states, hashes and storage before publication. The complete publication window now uses no allocations and verifies every resulting source/arrival hash.\n\nA second actual preimage (496) replaces an inactive source with an equal-state owner bound to a different thermal view. The old preflight accepted it, producing an incorrect recorded source hash. This successor retains exact source objects through the trial and checks every source identity before any source/thermal access, including inactive slices; retained lifetime prevents address reuse. Trial destruction after the accepted owner set is destroyed is tested.\n\nNormal502 and rebuilt coordinator/source/Hall ASan+UBSan507 pass all53 inherited controls plus whole-window allocation, replacement and detached-lifetime tests. The pinned native archive and external libraries were not rebuilt in this sanitizer run. Twelve continuation checks and six wrapper integrity controls pass. Actual qualified v8.0.2/c2 versus successor populated20/adaptive trajectories and ledgers are bit-identical (510); no timing claim is made.\n\nThe source v7.2 thermal forwarding wrapper is accepted only with its exact reviewed bytes and pinned original hash; exactly one defining thermal header is still required. Initial binding/fixture/provenance failures and their source preimages are retained. Initial comparison508 correctly refused outdated baseline c1 runner metadata;510 uses the matching qualified c2 receipt of the same binary.\n\nThis is the existing parent-chemistry coordinator with corrected transaction behavior, not adoption of the new full native/gas breakdown network or B3 emission bridge. It does not qualify a physical full run, complete chemistry, material response, sheath convergence or a speedup. Source and evidence are backed up outside projects. Large checkpoint payloads and build binaries remain in the local receipted outputs; receipts and bounded raw comparison ledgers are preserved here.\n')
sources={str(p.relative_to(A)):sha(p) for name in ['source','tests','tools','vendor'] for p in sorted((A/name).rglob('*')) if p.is_file()};sources['run_v8.py']=sha(A/'run_v8.py')
save(A/'SOURCE.json',{'files':sources,'identity':hashlib.sha256(json.dumps(sources,sort_keys=True).encode()).hexdigest()})
save(A/'RESULTS.json',{'normal_job':502,'sanitizer_job':507,'parity_job':510,'allocation_preimage_job':466,'source_replacement_preimage_job':496,'normal_checks':checks['c6'].get('passed'),'sanitizer_checks':checks['s1'].get('passed'),'continuation_checks':12,'publication_allocation_count':0,'detached_source_lifetime_qualified':True,'parity_fixtures':[c.get('fixture') for c in parity['comparisons']],'physical_breakdown_qualified':False,'speedup_qualified':False})
files={str(p.relative_to(A)):sha(p) for p in sorted(A.rglob('*')) if p.is_file()};save(A/'SEAL.json',{'files':files,'identity':hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest()})
shutil.copytree(A,B)
for p in A.rglob('*'):
 if p.is_file():assert p.read_bytes()==(B/p.relative_to(A)).read_bytes()
print(json.dumps({'artifact':str(A),'backup':str(B),'files':len(files)+1,'seal_sha256':sha(A/'SEAL.json')},indent=2))
