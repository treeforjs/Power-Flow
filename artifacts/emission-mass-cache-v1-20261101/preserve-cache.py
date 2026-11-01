from pathlib import Path
import json,hashlib,shutil
R=Path('/home/treefor/projects/runtimes/emission-mass-cache-v1-20261101');O=Path('/home/treefor/projects/.workspace/outputs')/R.name
A=Path('/home/treefor/projects/Power-Flow/artifacts/emission-mass-cache-v1-20261101');B=Path('/home/treefor/Documents/Codex/reconstructed-source-20261029')/A.name
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,sort_keys=True,indent=2)+'\n')
assert not A.exists() and not B.exists();final={}
for label,job in [('c4',530),('s2',531)]:
 d=json.loads((O/label/'receipt.json').read_text());assert d['passed'] and d['inputs_unchanged'] and int(d['job'])==job
 for path,h in d['input_hashes'].items():assert sha(Path(path))==h,path
 final[label]=d
assert (O/'c4/legacy-new.stdout').read_bytes()==(O/'c4/legacy-old.stdout').read_bytes()
A.mkdir()
for name in ['fast','baseline-fast','reference-kernel','tests','tools']:shutil.copytree(R/name,A/name)
shutil.copy2(R/'BASELINE.json',A/'BASELINE.json');shutil.copy2(Path(__file__),A/'preserve-cache.py')
for label in ['c1','c2','c3','s1','c4','s2']:
 D=A/'evidence'/label;D.mkdir(parents=True)
 for p in (O/label).iterdir():
  if p.suffix in ('.json','.stdout','.stderr') and (label=='c4' or not p.name.startswith('legacy-')):shutil.copy2(p,D/p.name)
for job in [516,519,521,522,530,531]:shutil.copy2(Path('/home/treefor/.local/state/simqueue/logs')/f'{job:06}.log',A/'evidence'/f'{job:06}.log')
(A/'README.txt').write_text('Exact emission cache for declared electron mass\n\nThe qualified source kernel supports declared electron mass and provenance, but the previous fast wrapper sent every nonlegacy mass directly to the scalar implementation. This successor permits bounded retarding-spectrum reuse and the existing guarded extracting fast path for those configurations. Cache keys include exact mass bits and policy identity. Barrier curvature, action and underflow bounds, Richardson factor and decay width use the declared mass; returned-distribution comparison checks both mass and identity. The scalar per-configuration self-check and fallback remain active.\n\nCPU530 and ASan/UBSan531 pass independent serialized field/node/refusal comparisons, cache-key isolation, invalid-configuration checks, explicit fallback and fast-path admission. Tests include default quadrature with the actual mass/provenance from the qualified registry fixture, adjacent mass bits, same-mass different policy and a diagnostic doubled mass. Old and new legacy spectra are byte-identical, including default quadrature. Raw legacy traces and hashes are retained. Initial516 rejected an uninitialized work-function value in the new test fixture at compile time; no numerical guard was changed.\n\nNo reference kernel, frozen userspace or active simulation was changed. This is an isolated component ready for explicit source-owner adoption, not a coupled performance qualification. No speedup is claimed from shared-core wall times; functional fast-path use and byte parity are established.\n')
source={str(p.relative_to(A)):sha(p) for name in ['fast','baseline-fast','reference-kernel','tests','tools'] for p in sorted((A/name).rglob('*')) if p.is_file()};save(A/'SOURCE.json',{'files':source,'identity':hashlib.sha256(json.dumps(source,sort_keys=True).encode()).hexdigest()})
save(A/'RESULTS.json',{'cpu_job':530,'sanitizer_job':531,'checks':{x:(O/x/'tests.stdout').read_text().strip() for x in ['c4','s2']},'legacy_trace_sha256':final['c4']['legacy_trace_sha256'],'actual_registry_mass_fixture':json.loads((A/'tests/qualified-registry-mass.json').read_text()),'fast_path_used_for_actual_registry_default_quadrature':True,'source_owner_adopted':False,'speedup_qualified':False})
files={str(p.relative_to(A)):sha(p) for p in sorted(A.rglob('*')) if p.is_file()};save(A/'SEAL.json',{'files':files,'identity':hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest()})
shutil.copytree(A,B)
for p in A.rglob('*'):
 if p.is_file():assert p.read_bytes()==(B/p.relative_to(A)).read_bytes()
print(json.dumps({'artifact':str(A),'backup':str(B),'files':len(files)+1,'source_identity':json.loads((A/'SOURCE.json').read_text())['identity'],'seal_sha256':sha(A/'SEAL.json')},indent=2))
