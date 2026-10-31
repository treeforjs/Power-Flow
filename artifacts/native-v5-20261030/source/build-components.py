"""Native v5 build: core archive (sealed exact-count native.o + v5 arena and
canonical objects), separate Hall adapter archive, genuine Hall v2 motion
objects for tests, and every test/benchmark executable, normal and ASan/UBSan.
Run under core-lease-20261030/both-cores.sh with CUDA_VISIBLE_DEVICES empty.
Usage: build-components.py [normal|sanitized|all] [target ...]"""
import sys
from tooling import *

resource_guard(build=True)
mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
only = set(sys.argv[2:])
B = P / 'build'
B.mkdir(exist_ok=True)


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
STRICT = ['-Wall', '-Wextra', '-Wpedantic', '-Werror']
FP = ['-ffp-contract=off', '-fno-fast-math', '-ffunction-sections', '-fdata-sections']
DEF = ['-DPIC_REACTION_NATIVE_SOURCE_SHA256="' + identity + '"']
TESTINC = ['-isystem', str(HALL), '-I' + str(P / 'tests')]
HALL_FLAGS = ['/usr/bin/c++', '-std=c++20', '-O2', '-DNDEBUG'] + FP + INC


def opt(level, sanitized):
    return (['-O1', '-g'] + ASAN) if sanitized else [level, '-g0']


def build(sanitized):
    tag = 'sanitized-' if sanitized else ''
    out = OUT / ('build-sanitized' if sanitized else 'build-normal')
    r = Runner(out, 'build')
    o = lambda name: B / (tag + name + '.o')
    want = lambda name: not only or name in only
    cxx = ['/usr/bin/mpic++', '-std=c++20']
    if want('core'):
        for src, name in [('source/src/mfem_pic_arena.cpp', 'arena'), ('source/src/mfem_pic_canonical.cpp', 'canonical')]:
            r.run(cxx + opt('-O2', sanitized) + FP + INC + DEF + STRICT + ['-c', P / src, '-o', o(name)], tag + name + '-compile.log')
        if not sanitized:
            r.run(cxx + ['-O2', '-g0'] + FP + INC + DEF + STRICT + ['-c', P / 'source/src/mfem_pic_hall_motion_adapter.cpp', '-o', o('hall-adapter')], 'hall-adapter-compile.log')
            r.run(cxx + ['-O0', '-g0'] + FP + INC + DEF + STRICT + ['-c', P / 'tests/test_hall_adapter_compile.cpp', '-o', o('hall-header-traits')], 'hall-header-traits-compile.log')
            r.run(['/usr/bin/nm', '-u', o('arena')], 'core-undefined-symbols.txt')
            r.run(['/usr/bin/nm', '-u', o('hall-adapter')], 'adapter-undefined-symbols.txt')
            need('hall_driver' not in (out / 'core-undefined-symbols.txt').read_text(), 'core acquired Hall dependency')
            need('hall_driver' in (out / 'adapter-undefined-symbols.txt').read_text(), 'adapter lost genuine Hall binding')
            for lib in [B / 'libpic_native_v5.a', B / 'libpic_hall_motion_adapter_v5.a']:
                lib.unlink(missing_ok=True)
            r.run(['/usr/bin/ar', 'rcs', B / 'libpic_native_v5.a', NATIVE_O, o('arena'), o('canonical')], 'core-archive.log')
            r.run(['/usr/bin/ar', 'rcs', B / 'libpic_hall_motion_adapter_v5.a', o('hall-adapter')], 'adapter-archive.log')
    if want('hall'):
        # Genuine Hall v2 motion provider, compiled from the Hall checkout with
        # Hall's own flags (no -Werror); used only by tests and the long run.
        for src, name in [(HALL_DRIVER / 'kinetic_boundary.cpp', 'kinetic-boundary'), (HALL_DRIVER / 'native_arena_motion.cpp', 'hall-motion')]:
            flags = HALL_FLAGS if not sanitized else (['/usr/bin/c++', '-std=c++20', '-O1', '-g'] + ASAN + FP + INC)
            r.run(flags + ['-c', src, '-o', o(name)], tag + name + '-compile.log')
    core = [o('arena'), o('canonical'), NATIVE_O]
    adapter = [B / 'hall-adapter.o'] if not sanitized else []
    if sanitized:
        r.run(cxx + ['-O1', '-g'] + ASAN + FP + INC + DEF + STRICT + ['-c', P / 'source/src/mfem_pic_hall_motion_adapter.cpp', '-o', o('hall-adapter')], tag + 'hall-adapter-compile.log') if want('core') else None
        adapter = [o('hall-adapter')]
    hall = [o('hall-motion'), o('kinetic-boundary')]
    tests = [('legacy', 'tests/test_arena.cpp', '-O0', [], False, False),
             ('ported', 'tests/test_candidate005_ported.cpp', '-O0', TRAPS, False, False),
             ('history', 'tests/test_v5_history.cpp', '-O0', TRAPS[:2], True, True),
             ('long-run', 'tests/long_run.cpp', '-O2', TRAPS[:2], True, True),
             ('digest-bench', 'tests/digest_bench.cpp', '-O2', [], False, True)]
    for name, src, level, wraps, with_hall, v5flag in tests:
        if not want(name) or (sanitized and name in ['long-run', 'digest-bench']):
            continue
        flags = cxx + opt(level, sanitized) + FP + INC + DEF + STRICT + (TESTINC if with_hall or v5flag else ['-I' + str(P / 'tests')]) + (['-DNATIVE_V5'] if v5flag else [])
        r.run(flags + ['-c', P / src, '-o', o('test-' + name)], tag + name + '-compile.log')
        link = ['/usr/bin/mpic++'] + (ASAN if sanitized else []) + ['-Wl,--wrap=' + w for w in wraps] + ['-Wl,--gc-sections', o('test-' + name)] + \
               (hall + adapter if with_hall else []) + core + LIBS + ['-o', B / (tag + 'test_' + name.replace('-', '_'))]
        r.run(link, tag + name + '-link.log')
    (out / 'build-result.json').write_text(json.dumps({'source_identity': identity, 'mode': 'sanitized' if sanitized else 'normal',
                                                       'utc': utc(), 'cpus': sorted(os.sched_getaffinity(0))}, indent=2) + '\n')


if mode in ['normal', 'all']:
    build(False)
if mode in ['sanitized', 'all']:
    build(True)
print('source identity', identity)
