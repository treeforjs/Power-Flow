"""Preserve a qualified component without copying large numerical payloads or binaries.
All actual runtime evidence remains in OUT; this copies source and reviewed receipts.
Never overwrites an artifact and refuses incomplete or mismatched qualification.
"""
from tooling import *
import shutil
import sys

cpu_milestone = sys.argv[1:] == ['--cpu-milestone']
need(not sys.argv[1:] or cpu_milestone, 'usage: preserve-successor.py [--cpu-milestone]')
A = Path('/home/treefor/projects/Power-Flow/artifacts') / (
    'native-v8.2-birth-memo-cpu-milestone-20261101' if cpu_milestone else 'native-v8.2-birth-memo-20261101')
need(not A.exists(), 'artifact already exists; do not overwrite a preserved release')
manifest = json.loads((P/'SOURCE-SHA256.json').read_text())
identity = manifest['source_identity']
archive = sha256(P/'build/libpic_native_v8.a')
inputs = sealed_inputs()
for name, digest in manifest['files'].items():
    need(sha256(P/name) == digest, 'source changed: '+name)

audits = {}
audit_cases = [('final-cpu', 14)]
if not cpu_milestone:
    audit_cases += [('final-gpu-core', 2), ('final-gpu-stageb', 2)]
for label, commands in audit_cases:
    q = json.loads((OUT/label/'audit-result.json').read_text())
    need(q['passed'] and q['source_identity'] == identity and q['library_sha256'] == archive,
         'audit source/archive mismatch: '+label)
    need(len(q['commands']) == commands, 'incomplete control suite: '+label)
    for c in q['commands']:
        need(c['returncode'] == 0 and sha256(OUT/label/c['log']) == c['log_sha256'], 'audit log mismatch')
    audits[label] = q

diffs = {}
diff_cases = [('functional-differential-c3', 8), ('functional-anchor-differential-c3', 17)]
if not cpu_milestone:
    diff_cases += [('owner-timing-c3', 17)]
for label, steps in diff_cases:
    q = json.loads((OUT/label/'differential-result.json').read_text())
    need(q['passed'] and not q['differences'], 'differential did not pass: '+label)
    need(q['build']['candidate_identity'] == identity and q['build']['candidate_library'] == archive,
         'differential candidate identity mismatch')
    need(q['build']['fixture_sha256'] == sha256(P/'tests/test_birth_memo_diff.cpp'), 'differential fixture changed')
    for run in q['runs']:
        r = json.loads((OUT/label/run['run']/'report.json').read_text())
        need(r['steps'] == steps and not r['failures'], 'differential run incomplete')
        need(sha256(OUT/label/run['run']/'physical.bin') == run['physical_sha256'], 'physical payload changed')
    diffs[label] = q

A.mkdir(parents=True)
def copy(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
for name in manifest['files']:
    copy(P/name, A/name)
for directory in ['tests', 'tools', 'gas-api-inventory']:
    for f in sorted((P/directory).rglob('*')):
        if f.is_file() and '__pycache__' not in f.parts:
            copy(f, A/f.relative_to(P))
for name in ['SOURCE-SHA256.json','SUCCESSOR-BASE.json','birth-memo.diff','EQUIVALENCE.txt',
             'HANDOFF.txt','BIRTH-MEMO-RESUME.json','tooling.py','build-components.py',
             'audit-tests.py','verify-final.py','build-differential.py','run-differential.py',
             'qualify-birth-memo.py','qualify-differential.py','qualify-differential-anchors.py','preserve-successor.py']:
    copy(P/name, A/name)
labels = ['build-normal','build-differential','build-differential-prior-fixture','build-differential-before-anchor-fields',
          *audits, 'functional-differential', *diffs]
evidence_manifest = {}
for label in labels:
    for f in sorted((OUT/label).rglob('*')):
        if not f.is_file():
            continue
        rel = str(f.relative_to(OUT))
        evidence_manifest[rel] = {'bytes':f.stat().st_size, 'sha256':sha256(f)}
        # Keep top-level command receipts/logs and bounded summary JSON. Full
        # checkpoints and raw physical bytes stay at their hashed runtime paths.
        if f.suffix in ('.json','.log','.txt') and f.stat().st_size <= 262144:
            copy(f, A/'evidence'/rel)
(A/'evidence-manifest.json').write_text(json.dumps(evidence_manifest, indent=2)+'\n')
qualification = {
    'schema':'powerflow-native-birth-memo-qualification-v1',
    'status':'cpu-milestone-gpu-qualification-pending' if cpu_milestone else 'qualified-component-not-full-coupled-run',
    'source_identity':identity,'archive_sha256':archive,'runtime':str(P),
    'evidence_root':str(OUT),'dependencies':inputs,
    'cpu_commands':14,'gpu_commands':0 if cpu_milestone else 4,
    'pending_gpu_jobs':[321,322] if cpu_milestone else [],
    'differential_cases':{k:{'passed':v['passed'], 'runs':v['runs'],
        'source_identity_bound_seal_differences':v['source_identity_bound_seal_differences'],
        'timing_uncontaminated':v['timing_uncontaminated']} for k,v in diffs.items()},
    'timing_scope':'Manufactured zero-field ordinary source births/sinks. Only owner-timing can qualify timing, and only if its measured SMT sibling check passed. No production forecast.',
    'coverage_limit':'Cap test proves early refusal, not memo branch instrumentation. Stage B separately covers wall/returned/bridge/internal-energy/merge variants.',
    'unchanged_guards':'All inherited audit corrections, frozen dependencies and strict eventful target-energy-owner refusal remain. No gas feedback or physical model changes.',
    'fixture_failures':json.loads((P/'BIRTH-MEMO-RESUME.json').read_text())['fixture_failures'],
    'production_qualified':False,'full_system_energy_closure':False}
(A/'QUALIFICATION.json').write_text(json.dumps(qualification,indent=2)+'\n')
need(sealed_inputs() == inputs, 'dependencies changed during preservation')
for name, digest in manifest['files'].items():
    need(sha256(P/name) == digest and sha256(A/name) == digest, 'source copy changed: '+name)
need(sha256(P/'build/libpic_native_v8.a') == archive, 'archive changed during preservation')
files = {str(f.relative_to(A)):sha256(f) for f in sorted(A.rglob('*')) if f.is_file()}
(A/'SEAL.json').write_text(json.dumps({'schema':'powerflow-native-birth-memo-seal-v1',
    'source_identity':identity,'archive_sha256':archive,
    'qualification_sha256':sha256(A/'QUALIFICATION.json'),'files':files},indent=2)+'\n')
print(A)
