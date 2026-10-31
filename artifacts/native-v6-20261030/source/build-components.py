"""Native v6 build: core archive (sealed exact-count native.o + v6 arena and
canonical objects), the v2 JSON Hall adapter archive, the binary motion v3
Hall adapter archive (Hall corrective v2 receipts), genuine Hall motion
objects (sealed 8dddec0 provider and the Hall corrective v2 provider, each
compiled with Hall's flags, no -Werror), the sealed source v5 support objects
for the source/arrival sweep, and every test/benchmark executable, normal and
ASan/UBSan. Two independent compiles run at a time (both leased cores).
Run under core-lease-20261030/both-cores.sh with CUDA_VISIBLE_DEVICES empty.
Usage: build-components.py [normal|sanitized|all] [target ...]"""
import sys
from concurrent.futures import ThreadPoolExecutor
from tooling import *

resource_guard(build=True)
mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
only = set(sys.argv[2:])
B = P / 'build'
B.mkdir(exist_ok=True)
inputs = sealed_inputs()


def source_manifest():
    files = {}
    for top in ['source', 'inputs']:
        for f in sorted((P / top).rglob('*')):
            if f.is_file():
                files[str(f.relative_to(P))] = sha256(f)
    identity = hashlib.sha256(json.dumps(files, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return {'scope': 'source/ and inputs/; tests/evidence separately hashed',
            'formula': "sha256(json.dumps(files, sort_keys=True, separators=(',', ':')))",
            'source_identity': identity, 'files': files}


manifest = source_manifest()
identity = manifest['source_identity']
(P / 'SOURCE-SHA256.json').write_text(json.dumps(manifest, indent=2) + '\n')
NATIVE_O = EXACT / 'build/native.o'
INC = ['-I' + str(P / 'source/include'), '-isystem', str(P / 'inputs'), '-I' + str(VENDOR)] + SYSTEM
INC5 = ['-I' + str(V5 / 'source/include'), '-isystem', str(V5 / 'inputs'), '-I' + str(VENDOR)] + SYSTEM
STRICT = ['-Wall', '-Wextra', '-Wpedantic', '-Werror']
FP = ['-ffp-contract=off', '-fno-fast-math', '-ffunction-sections', '-fdata-sections']
DEF = ['-DPIC_REACTION_NATIVE_SOURCE_SHA256="' + identity + '"']
TESTINC = ['-isystem', str(HALL), '-I' + str(P / 'tests')]
V6INC = ['-I' + str(HALL_V2 / 'source'), '-I' + str(ARRIVAL / 'include'), '-I' + str(P / 'source/adapters'),
         '-I' + str(V3_COMMIT / 'source'), '-DPF_ARRIVAL_EXACT_COUNT_HEADER="raikou/exact_particle_count.h"']
SOURCE_INC = ['-I' + str(SV3 / 'hall-owner/include'), '-I' + str(SV3 / 'hall-owner/src'), '-I' + str(S5_MOD),
              '-I' + str(SV3), '-I' + str(SV2), '-I' + str(SV2 / 'source-kernel/include'), *SYSTEM]


def opt(level, sanitized):
    return (['-O1', '-g'] + ASAN) if sanitized else [level, '-g0']


def run_all(r, jobs):
    """jobs: list of (command, log). Two at a time; any failure raises."""
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(r.run, c, log, 3600) for c, log in jobs]
        for f in futures:
            f.result()


def build(sanitized):
    tag = 'sanitized-' if sanitized else ''
    out = OUT / ('build-sanitized' if sanitized else 'build-normal')
    r = Runner(out, 'build')
    o = lambda name: B / (tag + name + '.o')
    want = lambda name: not only or name in only
    cxx = ['/usr/bin/mpic++', '-std=c++20']
    hall_cc = ['/usr/bin/c++', '-std=c++20'] + (['-O1', '-g'] + ASAN if sanitized else ['-O2', '-DNDEBUG']) + FP
    jobs = []
    if want('core'):
        for src, name in [('source/src/mfem_pic_arena.cpp', 'arena'), ('source/src/mfem_pic_canonical.cpp', 'canonical')]:
            jobs.append((cxx + opt('-O2', sanitized) + FP + INC + DEF + STRICT + ['-c', P / src, '-o', o(name)], tag + name + '-compile.log'))
        jobs.append((cxx + opt('-O2', sanitized) + FP + INC + DEF + STRICT +
                     ['-c', P / 'source/src/mfem_pic_hall_motion_adapter.cpp', '-o', o('hall-adapter')], tag + 'hall-adapter-compile.log'))
        jobs.append((cxx + opt('-O2', sanitized) + FP + INC + DEF + STRICT + ['-I' + str(HALL_V2 / 'source'), '-isystem', str(HALL),
                     '-c', P / 'source/src/mfem_pic_hall_motion_v3_adapter.cpp', '-o', o('hall-v3-adapter')], tag + 'hall-v3-adapter-compile.log'))
        if not sanitized:
            jobs.append((cxx + ['-O0', '-g0'] + FP + INC + DEF + STRICT + ['-c', P / 'tests/test_hall_adapter_compile.cpp', '-o', o('hall-header-traits')],
                         'hall-header-traits-compile.log'))
    if want('hall'):
        jobs.append((hall_cc + INC + ['-c', HALL_DRIVER / 'kinetic_boundary.cpp', '-o', o('kinetic-boundary')], tag + 'kinetic-boundary-compile.log'))
        jobs.append((hall_cc + INC + ['-c', HALL_DRIVER / 'native_arena_motion.cpp', '-o', o('hall-motion')], tag + 'hall-motion-compile.log'))
        jobs.append((hall_cc + ['-I' + str(HALL_V2 / 'source')] + INC + ['-isystem', str(HALL), '-c', HALL_V2 / 'source/motion_provider.cpp',
                     '-o', o('hall-v2-motion')], tag + 'hall-v2-motion-compile.log'))
        if not sanitized:
            # Reference builds against the sealed native v5 headers (v5 parity/benchmarks).
            jobs.append((hall_cc + INC5 + ['-c', HALL_DRIVER / 'native_arena_motion.cpp', '-o', o('v5-hall-motion')], 'v5-hall-motion-compile.log'))
            jobs.append((hall_cc + ['-I' + str(HALL_V2 / 'source')] + INC5 + ['-isystem', str(HALL), '-c', HALL_V2 / 'source/motion_provider.cpp',
                         '-o', o('v5-hall-v2-motion')], 'v5-hall-v2-motion-compile.log'))
    if want('source-support'):
        # Sealed source-kernel sources compiled as the source v5 CMake targets do.
        for name in ['surface_emission', 'micro_layer']:
            jobs.append((['/usr/bin/c++', '-std=gnu++20', '-O2', '-ffp-contract=off', '-fno-fast-math', '-I' + str(SV2 / 'source-kernel/include'),
                          '-c', SV2 / ('source-kernel/src/%s.cpp' % name), '-o', B / (name + '.o')], name + '-compile.log'))
        jobs.append((['/usr/bin/c++', '-std=gnu++20'] + opt('-O2', sanitized) + ['-Wall', '-Wextra', '-Wno-misleading-indentation',
                      '-ffp-contract=off', '-fno-fast-math', '-DHEAT_SOURCE_EXACT_COUNT_HEADER="%s"' % COUNT,
                      '-DPF_ARRIVAL_EXACT_COUNT_HEADER="%s"' % COUNT, '-I' + str(ARRIVAL / 'include'), '-I' + str(P / 'tests')] + SOURCE_INC +
                     ['-c', P / 'tests/source_v5_harness.cpp', '-o', o('source-v5-harness')], tag + 'source-v5-harness-compile.log'))
    tests = [  # name, source, level, wraps, include set, normal only
        ('legacy', 'tests/test_arena.cpp', '-O0', [], 'plain', False),
        ('ported', 'tests/test_candidate005_ported.cpp', '-O0', TRAPS, 'plain', False),
        ('history', 'tests/test_v5_history.cpp', '-O0', TRAPS[:2], 'sealed', False),
        ('flows', 'tests/test_v6_flows.cpp', '-O1', TRAPS[:2], 'v6', False),
        ('source-arrival', 'tests/test_v6_source_arrival.cpp', '-O1', TRAPS[:2], 'v6', False),
        ('parity', 'tests/parity_trace.cpp', '-O2', TRAPS[:2], 'sealed', True),
        ('bench', 'tests/bench_motion_v3.cpp', '-O2', [], 'bench', True),
        ('long-run', 'tests/long_run.cpp', '-O2', TRAPS[:2], 'sealed', True),
        ('digest-bench', 'tests/digest_bench.cpp', '-O2', [], 'sealed', True)]
    for name, src, level, wraps, kind, normal_only in tests:
        if not want(name) or (sanitized and normal_only):
            continue
        inc = {'plain': ['-I' + str(P / 'tests')], 'sealed': TESTINC + ['-DNATIVE_V5'], 'v6': TESTINC + V6INC,
               'bench': ['-I' + str(HALL_V2 / 'source'), '-isystem', str(HALL)]}[kind]
        jobs.append((cxx + opt(level, sanitized) + FP + INC + DEF + STRICT + inc + ['-c', P / src, '-o', o('test-' + name)], tag + name + '-compile.log'))
        if name == 'parity' and not sanitized:
            jobs.append((cxx + ['-O2', '-g0'] + FP + INC5 + ['-DPIC_REACTION_NATIVE_SOURCE_SHA256="' + V5_IDENTITY + '"'] + STRICT +
                         TESTINC + ['-DNATIVE_V5', '-c', P / src, '-o', o('v5-test-parity')], 'v5-parity-compile.log'))
        if name == 'bench' and not sanitized:
            jobs.append((cxx + ['-O2', '-g0'] + FP + INC5 + ['-DPIC_REACTION_NATIVE_SOURCE_SHA256="' + V5_IDENTITY + '"'] + STRICT +
                         inc + ['-DBENCH_V5', '-c', P / src, '-o', o('v5-test-bench')], 'v5-bench-compile.log'))
    run_all(r, jobs)
    # Archives and links.
    core = [o('arena'), o('canonical'), NATIVE_O]
    if want('core'):
        if not sanitized:
            r.run(['/usr/bin/nm', '-u', o('arena')], 'core-undefined-symbols.txt')
            r.run(['/usr/bin/nm', '-u', o('hall-adapter')], 'adapter-undefined-symbols.txt')
            r.run(['/usr/bin/nm', '-u', o('hall-v3-adapter')], 'v3-adapter-undefined-symbols.txt')
            need('hall_driver' not in (out / 'core-undefined-symbols.txt').read_text(), 'core acquired Hall dependency')
            need('hall_driver' in (out / 'adapter-undefined-symbols.txt').read_text(), 'adapter lost genuine Hall binding')
        for lib, objs in [('libpic_native_v6.a', core), ('libpic_hall_motion_adapter_v6.a', [o('hall-adapter')]),
                          ('libpic_hall_motion_v3_adapter_v6.a', [o('hall-v3-adapter')])]:
            path = B / (tag + lib)
            path.unlink(missing_ok=True)
            r.run(['/usr/bin/ar', 'rcs', path, *objs], tag + lib + '-archive.log')
    adapters = [o('hall-adapter'), o('hall-v3-adapter')]
    sealed_hall = [o('hall-motion'), o('kinetic-boundary')]
    v2_hall = [o('hall-v2-motion'), o('kinetic-boundary')]
    links = []
    for name, src, level, wraps, kind, normal_only in tests:
        if not want(name) or (sanitized and normal_only):
            continue
        extra = {'legacy': [], 'ported': [], 'history': sealed_hall + [o('hall-adapter')],
                 'flows': v2_hall + adapters,
                 'source-arrival': [o('source-v5-harness'), B / 'surface_emission.o', B / 'micro_layer.o', o('kinetic-boundary')],
                 'parity': sealed_hall + [o('hall-adapter')], 'bench': v2_hall + adapters,
                 'long-run': sealed_hall + [o('hall-adapter')], 'digest-bench': []}[name]
        exe = B / (tag + 'test_' + name.replace('-', '_'))
        links.append((['/usr/bin/mpic++'] + (ASAN if sanitized else []) + ['-Wl,--wrap=' + w for w in wraps] +
                      ['-Wl,--gc-sections', o('test-' + name)] + extra + core + LIBS + ['-o', exe], tag + name + '-link.log'))
        if name == 'parity' and not sanitized:
            links.append((['/usr/bin/mpic++'] + ['-Wl,--wrap=' + w for w in wraps] + ['-Wl,--gc-sections', o('v5-test-parity'),
                          o('v5-hall-motion'), o('kinetic-boundary'), V5 / 'build/libpic_hall_motion_adapter_v5.a',
                          V5 / 'build/libpic_native_v5.a'] + LIBS + ['-o', B / 'v5_test_parity'], 'v5-parity-link.log'))
        if name == 'bench' and not sanitized:
            links.append((['/usr/bin/mpic++', '-Wl,--gc-sections', o('v5-test-bench'), o('v5-hall-v2-motion'), o('kinetic-boundary'),
                           V5 / 'build/libpic_hall_motion_adapter_v5.a', V5 / 'build/libpic_native_v5.a'] + LIBS +
                          ['-o', B / 'v5_test_bench'], 'v5-bench-link.log'))
    run_all(r, links)
    (out / 'build-result.json').write_text(json.dumps({'source_identity': identity, 'mode': 'sanitized' if sanitized else 'normal',
                                                       'sealed_inputs': inputs, 'utc': utc(),
                                                       'cpus': sorted(os.sched_getaffinity(0))}, indent=2) + '\n')


if mode in ['normal', 'all']:
    build(False)
if mode in ['sanitized', 'all']:
    build(True)
print('source identity', identity)
