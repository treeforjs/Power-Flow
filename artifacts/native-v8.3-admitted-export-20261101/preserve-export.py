"""Preserve the admitted-export component and bounded qualification evidence."""
from tooling import *
import shutil,sys
cpu_only=sys.argv[1:]==['--cpu-milestone']
need(not sys.argv[1:] or cpu_only,'usage: preserve-export.py [--cpu-milestone]')
A=Path('/home/treefor/projects/Power-Flow/artifacts')/('native-v8.3-admitted-export-cpu-milestone-20261101' if cpu_only else 'native-v8.3-admitted-export-20261101')
need(not A.exists(),'artifact exists; never overwrite a sealed component')
manifest=json.loads((P/'SOURCE-SHA256.json').read_text())
identity=manifest['source_identity'];archive=sha256(P/'build/libpic_native_v8.a');inputs=sealed_inputs()
need(all(sha256(P/n)==v for n,v in manifest['files'].items()),'source changed')
fixtures=tree_hashes(P/'tests')
audits={}
for label,n in ([('sealed-cpu',14)] if cpu_only else [('sealed-cpu',14),('final-gpu',4)]):
    q=json.loads((OUT/label/'audit-result.json').read_text())
    need(q['passed'] and q['source_identity']==identity and q['library_sha256']==archive,'audit pins')
    need(len(q['commands'])==n and q['fixture_pins']==fixtures,'incomplete or stale audit fixtures')
    need(q['dependencies']==inputs,'audit dependencies')
    for name,digest in q['binary_pins'].items():need(sha256(P/name)==digest,'audit binary changed')
    for c in q['commands']:need(c['returncode']==0 and sha256(OUT/label/c['log'])==c['log_sha256'],'audit command evidence')
    audits[label]=q
q=json.loads((OUT/'sealed-export/result.json').read_text())
need(q['passed'] and q['source_identity']==identity and q['archive_sha256']==archive,'export pins')
need(q['fixture_pins']==fixtures and q['dependencies']==inputs and q['binary_sha256']==sha256(P/'build/test_v83_export'),'export fixtures/binary')
controls=json.loads((OUT/'sealed-export/controls/export-controls.json').read_text())
need(controls['passed'] and controls['checks']==12241,'export controls incomplete')
san=json.loads((OUT/'sanitized-export/result.json').read_text())
need(san['passed'] and san['source_identity']==identity and san['fixture_pins']==fixtures and san['dependencies']==inputs,'sanitizer pins')
need(len(san['commands'])==5 and san['archive_sha256']==sha256(P/'build/sanitized-libpic_native_v8.a'),'sanitizer archive/completeness')
for c in san['commands']:need(c['returncode']==0 and sha256(OUT/'sanitized-export'/c['log'])==c['log_sha256'],'sanitizer evidence')
A.mkdir(parents=True)
def copy(src,dst):
    dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
for name in manifest['files']:copy(P/name,A/name)
for f in sorted((P/'tests').rglob('*')):
    if f.is_file() and '__pycache__' not in f.parts:copy(f,A/f.relative_to(P))
for name in ['tooling.py','build-components.py','audit-tests.py','check-export.py','qualify-export.py',
             'qualify-export-next.py','qualify-export-binding.py','qualify-export-sealed.py','qualify-export-sanitized.py',
             'preserve-export.py','SOURCE-SHA256.json','BIRTH-MEMO-BASE.json','ADOPTION.json','admitted-export.diff',
             'API-ADMITTED-EXPORT.txt','HANDOFF.txt','EXPORT-RESUME.json']:
    copy(P/name,A/name)
labels=['build-core-initial','build-normal','build-sanitized','build-sanitized-limit3G','first-export','first-cpu','final-export','final-cpu',
        'binding-export','binding-cpu','sealed-export','sealed-cpu','sanitized-export','final-gpu']
local={};exclusions=[]
for label in labels:
    for f in sorted((OUT/label).rglob('*')):
        if not f.is_file():continue
        rel=str(f.relative_to(OUT));local[rel]={'bytes':f.stat().st_size,'sha256':sha256(f)}
        if f.name in ['processes-before.txt','processes-after.txt']:
            exclusions.append({'path':rel,'reason':'unrelated full process command lines remain local only'});continue
        if f.suffix in ('.json','.log','.txt') and f.stat().st_size<=262144:copy(f,A/'evidence'/rel)
(A/'evidence-manifest.json').write_text(json.dumps(local,indent=2)+'\n')
(A/'LOCAL-ONLY-EVIDENCE.json').write_text(json.dumps({'excluded':exclusions,'policy':'Large checkpoints, binaries and raw payloads remain at hashed runtime evidence paths; full machine process listings are never published.'},indent=2)+'\n')
qualification={'schema':'native-v8.3-admitted-export-component-v1','source_identity':identity,'archive_sha256':archive,
 'runtime':str(P),'evidence_root':str(OUT),'status':'CPU-and-sanitizer-milestone-GPU-pending' if cpu_only else 'qualified-diagnostic-component-not-coupled-closure',
 'dependencies':inputs,'cpu_job':366,'gpu_job':363,'sanitizer_job':385,'cpu_commands':14,'gpu_commands':0 if cpu_only else 4,'pending_gpu':cpu_only,
 'export_assertions':controls['checks'],'sanitizer_commands':5,
 'controls':['raw/admitted prefix separation; inventory suffix refusal; GAS nondebit',
 'real whole and partial target debits; exact nonbinary multiple-carrier remainder',
 'tampered bindings/counts/events refused; empty stage export',
 'fresh immutable retry object; moved/foreign/published native token preflight',
 'abandoned outer leaves accepted state/RNG unchanged',
 'export digests survive before/at-anchor restart; changed retained commitment and old schema refused',
 'GPU export parity and device Stage B controls pending' if cpu_only else 'actual host/device_verified/device and partition export digest parity; device Stage B controls'],
 'limits':['No gas feedback, gas authority, common-owner publication, target energy closure, or physical breakdown run.',
 'Actual projectile/target overlap and carried internal consumption are guarded; only manufactured flag/refusal controls qualify these export branches.',
 'No peak-memory or performance qualification; compact target gather grows from 5 to 11 words per inspected target.',
 'No backwards checkpoint migration; source/ABI identity changes require complete dependent rebuild.'],
 'prior_fixture_failure':'Job 361 failed compilation of new test due to sha256 versus state_sha256 member typo; corrected test only; numerical archive unchanged. Job367 compiler exceeded inherited3GiB address-space limit; job385 passed with7GiB compiler AS inside8GiB scheduler memory. No physics/resource execution guard changed.',
 'production_qualified':False,'strict_target_owner_gate_open':False,'full_system_energy_closure':False}
(A/'QUALIFICATION.json').write_text(json.dumps(qualification,indent=2)+'\n')
need(inputs==sealed_inputs(),'dependencies changed during preservation')
need(all(sha256(P/n)==v and sha256(A/n)==v for n,v in manifest['files'].items()),'source changed during copy')
need(sha256(P/'build/libpic_native_v8.a')==archive,'archive changed')
files={str(f.relative_to(A)):sha256(f) for f in sorted(A.rglob('*')) if f.is_file()}
(A/'SEAL.json').write_text(json.dumps({'schema':'native-v8.3-admitted-export-seal-v1','source_identity':identity,'archive_sha256':archive,'files':files},indent=2)+'\n')
print(A)
print('seal',sha256(A/'SEAL.json'))
