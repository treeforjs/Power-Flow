"""Builds the long-run driver against the SEALED candidate003 archives (read only)
plus the genuine Hall v2 motion provider, for the same-fixture comparison."""
import sys
from tooling import *

resource_guard(build=True)
out = OUT / 'baseline-candidate003'
r = Runner(out, 'build')
B = P / 'build' / 'c003'
B.mkdir(parents=True, exist_ok=True)
c003_identity = json.loads((C003 / 'SOURCE-SHA256.json').read_text())['source_identity']
pins = {str(p): sha256(p) for p in [C003 / 'build/libpic_native_candidate003.a', C003 / 'build/libpic_hall_motion_adapter003.a',
                                     HALL_DRIVER / 'native_arena_motion.cpp', HALL_DRIVER / 'kinetic_boundary.cpp',
                                     HALL_DRIVER / 'motion_receipt.hpp', HALL_DRIVER / 'kinetic_boundary.hpp']}
(out / 'input-pins.json').write_text(json.dumps(pins, indent=2) + '\n')
hall_flags = ['/usr/bin/c++', '-std=c++20', '-O2', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math',
              '-ffunction-sections', '-fdata-sections', '-I' + str(C003 / 'source/include'),
              '-isystem', str(C003 / 'inputs'), '-I' + str(VENDOR)] + SYSTEM
for src, name in [(HALL_DRIVER / 'kinetic_boundary.cpp', 'kinetic-boundary'), (HALL_DRIVER / 'native_arena_motion.cpp', 'hall-motion')]:
    r.run(hall_flags + ['-c', src, '-o', B / (name + '.o')], 'compile-' + name + '.log')
strict = ['/usr/bin/mpic++', '-std=c++20', '-O2', '-g0', '-ffp-contract=off', '-fno-fast-math', '-ffunction-sections',
          '-fdata-sections', '-I' + str(C003 / 'source/include'), '-isystem', str(C003 / 'inputs'), '-I' + str(VENDOR)] + SYSTEM + \
         ['-isystem', str(HALL), '-I' + str(P / 'tests'), '-Wall', '-Wextra', '-Wpedantic', '-Werror']
r.run(strict + ['-c', P / 'tests/long_run.cpp', '-o', B / 'long-run.o'], 'compile-long-run.log')
link = ['/usr/bin/mpic++', '-Wl,--wrap=' + TRAPS[0], '-Wl,--wrap=' + TRAPS[1], '-Wl,--gc-sections', B / 'long-run.o',
        B / 'hall-motion.o', B / 'kinetic-boundary.o', C003 / 'build/libpic_hall_motion_adapter003.a',
        C003 / 'build/libpic_native_candidate003.a'] + LIBS + ['-o', B / 'long_run_c003']
r.run(link, 'link-long-run.log')
(out / 'build-result.json').write_text(json.dumps({'source_identity': c003_identity, 'binary': str(B / 'long_run_c003'),
                                                   'binary_sha256': sha256(B / 'long_run_c003'), 'utc': utc()}, indent=2) + '\n')
print('built', B / 'long_run_c003')
