from pathlib import Path
import hashlib,json,tarfile,io,shutil
W=Path('/home/treefor/sim-workspace/powerflow');R=Path('/home/treefor/projects/Power-Flow');O=R/'artifacts/cathode-field-and-water-review-20261029';items={}
def add(p,key):
 if p.is_file() and p.suffix not in {'.pdf','.o','.a','.so'} and p.read_bytes()[:4]!=b'\x7fELF':items[key]=p
A=W/'outputs/cathode-emission-audit-20261029'
for p in A.rglob('*'):
 if 'build' not in p.relative_to(A).parts and '__pycache__' not in p.parts:add(p,'diagnostics/'+p.relative_to(A).as_posix())
F=W/'outputs/pic-profile/cuda-gather-investigation/actual-midpoint-fixture'
for p in F.iterdir():add(p,'fixture/'+p.name)
S=Path('/home/treefor/projects/runtimes/pic-static-domain-20261029')
for name in ['CMakeLists.txt','README.txt','input-bindings.json']:
 add(S/name,'static-source/'+name)
for sub in ['source','tools']:
 for p in (S/sub).rglob('*'):
  if '__pycache__' not in p.parts:add(p,'static-source/'+p.relative_to(S).as_posix())
G=Path('/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029')
for name in ['batch_snapshot.hpp','batch_snapshot.cpp','gather_push_core.hpp']:
 add(G/'source'/name,'snapshot-support/'+name)
for name in ['host_guard.py','run-frozen-with-host-proc']:add(G/'tools'/name,'snapshot-support/'+name)
P=W/'outputs/supplied-chicago-papers-20261029'
for p in P.rglob('*'):
 if 'tmp' not in p.relative_to(P).parts and '__pycache__' not in p.parts and p.suffix in {'.json','.txt','.py'} and p.name!='paper.txt':add(p,'literature-review/'+p.relative_to(P).as_posix())
P=W/'outputs/pic-profile/fragmentation-raikou-plan'
for p in P.rglob('*'):
 if not any(x in p.relative_to(P).parts for x in ['source-pages','source-documents','__pycache__']) and p.suffix in {'.json','.txt','.py','.csv','.sh'}:add(p,'water-data-handoff/'+p.relative_to(P).as_posix())
for p in (W/'outputs/voltage-pulse-review-20261029').iterdir():add(p,'voltage/'+p.name)
for name in ['gap-pulse-midplane110kV.csv','gap-pulse-midplane110kV.json','gap-pulse-220kV-anchored.csv','gap-pulse-220kV-anchored.json','electron-mapping-resume.json','pic-profiling-resume.json']:
 add(W/'outputs'/name,'state-and-inputs/'+name)
archive=O/'source-and-diagnostic-evidence.tar.gz'
manifest_path=O/'members-sha256.json'
if archive.exists():
 manifest=json.loads(manifest_path.read_text());items={k:Path(v['source']) for k,v in manifest.items()}
else:
 manifest={k:{'source':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in sorted(items.items())};manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
 with tarfile.open(archive,'w:gz',compresslevel=6) as t:
  for k,p in sorted(items.items()):
   data=p.read_bytes();assert hashlib.sha256(data).hexdigest()==manifest[k]['sha256'],str(p)
   info=tarfile.TarInfo(k);info.size=len(data);info.mode=0o644;t.addfile(info,io.BytesIO(data))
with tarfile.open(archive,'r:gz') as t:
 members=t.getmembers();assert len(members)==len(manifest)
 for member in members:
  assert member.isfile();data=t.extractfile(member).read();assert hashlib.sha256(data).hexdigest()==manifest[member.name]['sha256']
archive_bytes=archive.stat().st_size;archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest();parts=[]
with archive.open('rb') as f:
 n=0
 while data:=f.read(48*1024**2):
  n+=1;q=O/f'source-and-diagnostic-evidence.tar.gz.part{n:03d}';q.write_bytes(data);parts.append({'name':q.name,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
joined=hashlib.sha256()
for part in parts:
 data=(O/part['name']).read_bytes();assert hashlib.sha256(data).hexdigest()==part['sha256'];joined.update(data)
assert joined.hexdigest()==archive_sha256
archive.unlink()
q=O/'surface_births_and_later_cloud.png';shutil.copy2(A/'preview/surface_births_and_later_cloud.png',q)
shutil.copy2(W/'outputs/voltage-pulse-review-20261029/voltage_waveform_used.png',O/'voltage_waveform_used.png')
shutil.copy2(W/'outputs/supplied-chicago-papers-20261029/H2O-readiness.txt',O/'H2O-readiness.txt')
receipt={'files':len(manifest),'uncompressed_bytes':sum(x['bytes'] for x in manifest.values()),'archive_bytes':archive_bytes,'archive_sha256':archive_sha256,'ordered_parts':parts,'joined_parts_sha256_verified':True,'all_members_reopened_and_hash_verified':True,'full_third_party_PDF_or_text_or_renders_included':False,'third_party_binary_libraries_included':False,'long_PIC_resumed':False}
(O/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
for n in ['README.txt']:
 (O/n).write_text('Surface source/cloud audit, fixed-state electric exterior sensitivity, water-data readiness and supplied CHICAGO paper scope.\n\nThe original9mm strict-capture comparison remains FAILED and preserved. The separate CPU causal audit explains its fresh-solve nodal differences, without changing solver/source/physical guards. Larger static domains ran under a separately documented continuation; no particles advanced.9-to25mm fixed-state vector E differences are14.03% inside the finite gap and39.90% overall. Field-domain/sheath/material/full3D convergence and full H2O kinetics are not established.\n\nThe archive retains own simulation fixtures, source, receipts, ledgers, plots, source-data hashes and authored literature/data summaries. Full external PDFs/text/page renders and compiled libraries are excluded. Join the numbered gzip-tar parts in order to restore the archive. Every archive member was reopened and hash-verified; the joined part stream reproduces the original archive hash. The run remains explicitly paused, dependencies remain frozen, and changed source/field inputs require a new supported-start history.\n')
B=Path('/home/treefor/sim-workspace/source-backups/reconstructed-source-20261029/cathode-field-water-review-20261029');B.mkdir(exist_ok=False)
for p in O.iterdir():
 if p.is_file():shutil.copy2(p,B/p.name)
print(json.dumps(receipt,indent=2))
