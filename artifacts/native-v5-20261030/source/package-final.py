"""Packages native v5 after verify-final.py and long-run.py have passed:
immutable snapshot + source archive, filled COORDINATOR-V4-LINK.txt, handoff."""
import shutil, tarfile
from tooling import *

resource_guard()
F = OUT / 'final'
v = json.loads((F / 'verification-result.json').read_text())
manifest = json.loads((P / 'SOURCE-SHA256.json').read_text())
identity = manifest['source_identity']
need(v['status'] == 'passed-bounded-components' and v['source_identity'] == identity, 'verification identity')
for name, digest in manifest['files'].items():
    need(sha256(P / name) == digest, 'source changed after verification: ' + name)
summary = json.loads((OUT / 'long-run/long-run-summary.json').read_text())
need(summary['pins']['v5_source_identity'] == identity, 'long run used another source identity')
lib, adapter = P / 'build/libpic_native_v5.a', P / 'build/libpic_hall_motion_adapter_v5.a'
genesis = json.loads((F / 'ported/genesis.json').read_text())['payload']
text = (P / 'COORDINATOR-V4-LINK.template.txt').read_text()
for key, value in {'@V5_IDENTITY@': identity, '@V5_LIB_SHA@': sha256(lib), '@V5_ADAPTER_SHA@': sha256(adapter),
                   '@V5_GENESIS_PREFIX@': genesis['boundary']['accepted_arena_sha256'],
                   '@V5_GENESIS_INITIAL@': genesis['initial_coupling_prefix'],
                   '@V5_GENESIS_CHAIN@': genesis['metadata']['chain']}.items():
    text = text.replace(key, value)
need('@' not in text.replace('`', ''), 'unfilled template field')
(P / 'COORDINATOR-V4-LINK.txt').write_text(text)
S = P / 'snapshots' / ('candidate005-' + identity[:8])
S.mkdir(parents=True, exist_ok=False)
for name in ['source', 'inputs', 'tests']:
    shutil.copytree(P / name, S / name)
for name in ['SOURCE-SHA256.json', 'API-005.txt', 'COORDINATOR-V4-LINK.txt', 'tooling.py', 'build-components.py',
             'build-baseline-c003.py', 'verify-final.py', 'long-run.py', 'v4-link-check.py', 'package-final.py']:
    shutil.copy2(P / name, S / name)
A = OUT / ('native-full-time-candidate005-source.tar.gz')
with tarfile.open(A, 'w:gz') as archive:
    for path in sorted(S.rglob('*')):
        archive.add(path, arcname=str(path.relative_to(S)), recursive=False)
files = {str(p.relative_to(S)): sha256(p) for p in sorted(S.rglob('*')) if p.is_file()}
(OUT / 'source-archive-manifest.json').write_text(json.dumps(files, indent=2) + '\n')
shutil.copy2(P / 'SOURCE-SHA256.json', OUT / 'SOURCE-SHA256.json')
shutil.copy2(P / 'API-005.txt', OUT / 'API-005.txt')
shutil.copy2(P / 'COORDINATOR-V4-LINK.txt', OUT / 'COORDINATOR-V4-LINK.txt')
artifacts = {str(p): sha256(p) for p in [lib, adapter, P / 'build/test_legacy', P / 'build/test_ported',
                                          P / 'build/test_history', P / 'build/test_long_run']}
artifacts[str(EXACT / 'build/native.o')] = sha256(EXACT / 'build/native.o')
handoff = {'status': 'ready_for_read_only_component_review', 'source_identity': identity,
           'candidate_prefix': str(P), 'immutable_snapshot': str(S), 'source_archive': str(A),
           'archive_sha256': sha256(A), 'base': {'candidate003_identity': v['baseline']['candidate003_identity'],
                                                 'rebuilt_from_archive_sha256': v['baseline']['rebuilt_from_archive_sha256']},
           'libraries': {'core': str(lib), 'hall_adapter': str(adapter)}, 'artifact_sha256': artifacts,
           'verification': str(F / 'verification-result.json'), 'long_run': str(OUT / 'long-run/long-run-summary.json'),
           'v4_link_check': str(OUT / 'v4-link-check/v4-link-check.json'),
           'remaining_unqualified': ['Nonempty full-time disk restore (anchor chain is serializable; restore deferred)',
                                     'Births, sinks/collection, active layers and competing channel sets (H2-H4)',
                                     'Identity tombstone sets grow with events (bounded by explicit cap only)',
                                     'Hall v2 motion record is JSON-sealed and capped at 1e5 carriers',
                                     'GPU block-hash computation (design and host verification only)'],
           'utc': utc()}
(OUT / 'handoff.json').write_text(json.dumps(handoff, indent=2) + '\n')
(OUT / 'READY').write_text('Native v5 ready for read-only component review. Source ' + identity + '.\n')
outputs = {str(p.relative_to(OUT)): sha256(p) for p in sorted(OUT.rglob('*'))
           if p.is_file() and p.name != 'output-manifest.json' and 'superseded-dev' not in p.parts}
(OUT / 'output-manifest.json').write_text(json.dumps(outputs, indent=2) + '\n')
print(json.dumps({'source_identity': identity, 'archive_sha256': handoff['archive_sha256'], 'snapshot': str(S)}, indent=1))
