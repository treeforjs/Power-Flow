"""Native v5 final verification. Run under one Power-Flow core lease:
  core-lease-20261030/any-core.sh --wait env CUDA_VISIBLE_DEVICES= python3 verify-final.py
Rebuilds everything (normal and ASan/UBSan), runs the legacy (368), ported
candidate003 and new v5 suites, checks parity and seals, writes
OUT/final/verification-result.json. Never modifies candidate003 or its evidence."""
import sys, tarfile
from tooling import *

cpus = resource_guard(build=True)
F = OUT / 'final'
F.mkdir(parents=True, exist_ok=True)
r = Runner(F, 'final')
C003_OUT = OUT.parent / 'native-full-time-candidate003'
EXACT_OUT = OUT.parent / 'native-exact-count'

# 1. The v5 baseline is the sealed candidate003 tarball, byte for byte.
c003_manifest = json.loads((C003 / 'SOURCE-SHA256.json').read_text())
for name, digest in c003_manifest['files'].items():
    need(sha256(C003 / name) == digest, 'sealed candidate003 changed: ' + name)
    need(sha256(P / 'reference/candidate003-c045' / name) == digest, 'reference baseline differs: ' + name)
handoff = json.loads((C003_OUT / 'handoff.json').read_text())
need(sha256(C003_OUT / 'native-full-time-candidate003-source.tar.gz') == handoff['archive_sha256'],
     'candidate003 evidence archive changed')
with tarfile.open(C003_OUT / 'native-full-time-candidate003-source.tar.gz') as archive:
    for member in archive.getmembers():
        if member.isfile():
            data = archive.extractfile(member).read()
            need(hashlib.sha256(data).hexdigest() == sha256(P / 'reference/candidate003-c045' / member.name),
                 'reference differs from archive: ' + member.name)
for name in ['inputs/' + p for p in ['collision_scheduler.hpp', 'electrode_heat_feedback.hpp', 'h2o-gas-20261029.json',
                                     'heat-v2-count.hpp', 'water_parent_bridge.hpp']] + \
            ['inputs/hall-motion-v2/' + p for p in ['kinetic_boundary.cpp', 'kinetic_boundary.hpp', 'motion_receipt.hpp',
                                                   'native_arena_motion.cpp']] + \
            ['source/include/raikou/exact_particle_count.h', 'source/include/raikou/mfem_gap_pic.h',
             'source/src/mfem_gap_pic.cpp', 'source/src/mfem_pic_hall_motion_adapter.cpp']:
    need(sha256(P / name) == c003_manifest['files'][name], 'v5 must keep candidate003 bytes: ' + name)
for name in ['kinetic_boundary.cpp', 'kinetic_boundary.hpp', 'motion_receipt.hpp', 'native_arena_motion.cpp']:
    need(sha256(HALL_DRIVER / name) == sha256(P / 'inputs/hall-motion-v2' / name), 'Hall provider differs: ' + name)

# 2. Build normal and sanitized.
# build-components.py runs every compiler through the frozen launcher itself.
with (F / 'build.log').open('w') as log:
    built = subprocess.run([sys.executable, str(P / 'build-components.py'), 'all'], stdout=log,
                           stderr=subprocess.STDOUT, timeout=7200, cwd=str(P))
need(built.returncode == 0, 'build failed; see final/build.log')
manifest = json.loads((P / 'SOURCE-SHA256.json').read_text())
identity = manifest['source_identity']

# 3. Normal suites.
B = P / 'build'
dirs = {k: F / k for k in ['legacy', 'ported', 'history', 'sanitized-legacy', 'sanitized-ported', 'sanitized-history']}
for d in dirs.values():
    d.mkdir(exist_ok=True)
r.run([B / 'test_legacy', P / 'inputs/h2o-gas-20261029.json', dirs['legacy']], 'legacy-run.log')
r.run([B / 'test_ported', dirs['ported']], 'ported-run.log')
r.run([B / 'test_history', dirs['history']], 'history-run.log', timeout=7200)
r.run([B / 'test_digest_bench', F / 'digest-scaling.json'], 'digest-bench.log', timeout=3600)
# 4. Sanitized suites (sanitizer RSS limit instead of RLIMIT_AS).
r.run([B / 'sanitized-test_legacy', P / 'inputs/h2o-gas-20261029.json', dirs['sanitized-legacy']],
      'sanitized-legacy-run.log', sanitized=True, timeout=3600)
r.run([B / 'sanitized-test_ported', dirs['sanitized-ported']], 'sanitized-ported-run.log', sanitized=True, timeout=3600)
r.run([B / 'sanitized-test_history', dirs['sanitized-history'], 'reduced'], 'sanitized-history-run.log',
      sanitized=True, timeout=7200)

# 5. Parity.
load = lambda p: json.loads(Path(p).read_text())
legacy = load(dirs['legacy'] / 'arena-tests.json')
need(legacy == load(EXACT_OUT / 'tests/arena-tests.json') == load(C003_OUT / 'tests/arena-tests.json') ==
     load(dirs['sanitized-legacy'] / 'arena-tests.json'), 'legacy results/proofs differ')
ported = load(dirs['ported'] / 'candidate005-ported-tests.json')
need(ported == load(dirs['sanitized-ported'] / 'candidate005-ported-tests.json'), 'ported normal/sanitized differ')
c003 = load(C003_OUT / 'component-tests/candidate003-tests.json')
sweep = {'allocation failure before private publish leaves token/base reusable',
         'retry after private publication allocation failure succeeds'}
named = lambda j: [x['name'] for x in j['checks'] if x['name'] not in sweep]
need(named(ported) == named(c003) and all(x['passed'] for x in ported['checks']) and
     ported['forbidden_native_solver_calls'] == 0 and ported['forbidden_native_advance_calls'] == 0,
     'ported candidate003 controls differ')
history = load(dirs['history'] / 'v5-history-tests.json')
san_history = load(dirs['sanitized-history'] / 'v5-history-tests.json')
normal_only = {'1e6 resident carriers: first exposure publishes', '1e6 resident carriers: complete native census',
               'Hall v2 motion refuses a 1e6 census (provider cap, not native)'}
h_names = [x['name'] for x in history['checks']]
# The last population check appears once per population; drop the 1e6 occurrence.
trimmed, seen_validate = [], 0
for name in h_names:
    if name in normal_only:
        continue
    if name == 'resident carrier population validates at the configured cap':
        seen_validate += 1
        if seen_validate == 2:
            continue
    trimmed.append(name)
need(all(x['passed'] for x in history['checks'] + san_history['checks']) and
     trimmed == [x['name'] for x in san_history['checks']], 'v5 history normal/sanitized check lists differ')
need(json.loads((P / 'SOURCE-SHA256.json').read_text()) == manifest, 'source changed during verification')
for name, digest in c003_manifest['files'].items():
    need(sha256(C003 / name) == digest, 'sealed candidate003 changed during verification: ' + name)
artifacts = {str(p.relative_to(P)): sha256(p) for p in sorted(B.glob('*')) if p.is_file()}
(F / 'build-artifact-hashes.json').write_text(json.dumps(artifacts, indent=2) + '\n')
result = {
    'status': 'passed-bounded-components',
    'source_identity': identity,
    'baseline': {'candidate003_identity': c003_manifest['source_identity'],
                 'rebuilt_from_archive_sha256': handoff['archive_sha256']},
    'legacy_checks': legacy['check_count'],
    'legacy_identical_to_exact_count_and_candidate003_and_sanitized': True,
    'ported_candidate003_controls': {'check_count': ported['check_count'], 'named_controls': len(named(ported)),
                                     'allocation_failure_sweep': ported['allocation_failure_sweep'],
                                     'normal_sanitized_identical': True},
    'v5_checks': {'normal': history['check_count'], 'sanitized_reduced': san_history['check_count'],
                  'normal_only_controls': sorted(normal_only) + ['resident carrier population validates at the configured cap (1e6)']},
    'v5_measurements': history['measurements'],
    'digest_scaling': load(F / 'digest-scaling.json')['digest_scaling'],
    'forbidden_native_solver_or_advance_calls': history['forbidden_native_solver_or_advance_calls'],
    'sanitizers': 'address/undefined/leak on v5 arena, canonical, adapter, Hall motion objects and tests; '
                  'sealed exact-count native.o and MFEM uninstrumented; hard_rss_limit_mb=3072',
    'resources': {'cpus': cpus, 'nice': 15, 'threads': 1, 'jobs': 1, 'GPU': False,
                  'rlimit_as_bytes_normal': RLIMIT_AS, 'launcher': FROZEN},
    'genuine_joint_step': False,
    'hall_motion': 'genuine Hall v2 provider at identically zero total field; synthetic field/heat digests',
    'nonempty_disk_restore': 'deferred (typed refusal after validating exported anchors/window)',
    'utc': utc()}
(F / 'verification-result.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: result[k] for k in ['status', 'source_identity', 'legacy_checks', 'ported_candidate003_controls',
                                          'v5_checks']}, indent=1))
