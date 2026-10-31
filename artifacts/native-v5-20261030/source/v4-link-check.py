"""Read-only link check: compile the unmodified coordinator v4 Stage 1 sources
against the native v5 headers and link them with the v5 archives, using the v4
runner's own include order and options with the native path, identity and
library switched to v5 (the exact change list in COORDINATOR-V4-LINK.txt).
Nothing in the v4 directory is written; objects go to build/v4-link-check."""
from tooling import *

resource_guard(build=True)
V4 = RT.parent / 'pic-water-common-commit-v4-20261030'
V3 = RT.parent / 'pic-water-common-commit-v3-20261030'
I = RT.parent / 'pic-water-full-time-subflow-20261029/inputs'
H = HALL.parent
HC = HALL / 'hall-field-correction-20261030'
out = OUT / 'v4-link-check'
r = Runner(out, 'v4-link-check')
B = P / 'build/v4-link-check'
B.mkdir(parents=True, exist_ok=True)
identity = json.loads((P / 'SOURCE-SHA256.json').read_text())['source_identity']
v4_files = {str(p.relative_to(V4)): sha256(p) for p in sorted((V4 / 'source').glob('*')) + sorted((V4 / 'tests').glob('*'))}
common = ['-std=c++20', '-O2', '-ffp-contract=off', '-fno-fast-math', '-ffunction-sections', '-fdata-sections']
isystem = []
for d in (P / 'source/include', I, I / 'include', P / 'inputs', HC, H / 'development',
          MFEM / 'install-mpi/include', MFEM / 'tpl/hypre/include', MFEM / 'tpl/metis/include'):
    isystem += ['-isystem', str(d)]
ours = ['/usr/bin/mpic++', *common, '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-I' + str(V4 / 'source'),
        '-I' + str(V3 / 'source'), '-DPIC_REACTION_NATIVE_SOURCE_SHA256="' + identity + '"', *isystem]
hall = ['/usr/bin/c++', *common, '-DNDEBUG', *isystem]
objects = {'native_arena_motion': (hall, HALL_DRIVER / 'native_arena_motion.cpp'),
           'kinetic_boundary': (hall, HALL_DRIVER / 'kinetic_boundary.cpp'),
           'hall_motion_adapter': (ours, P / 'source/src/mfem_pic_hall_motion_adapter.cpp'),
           'alloc_hook': (ours, V4 / 'tests/alloc_hook.cpp'),
           'test_stage1_joint': (ours, V4 / 'tests/test_stage1_joint.cpp')}
status = {}
for name, (flags, src) in objects.items():
    code = r.run([*flags, '-c', src, '-o', B / (name + '.o')], 'compile-' + name + '.log', check=False)
    status[name] = code
link_code = None
if all(c == 0 for c in status.values()):
    link = ['/usr/bin/mpic++', *['-Wl,--wrap=' + w for w in TRAPS], '-Wl,--gc-sections',
            *[B / (k + '.o') for k in objects], P / 'build/libpic_native_v5.a', *LIBS, '-o', B / 'test_stage1_joint_v5']
    link_code = r.run(link, 'link.log', check=False)
after = {str(p.relative_to(V4)): sha256(p) for p in sorted((V4 / 'source').glob('*')) + sorted((V4 / 'tests').glob('*'))}
need(after == v4_files, 'v4 sources changed during the read-only check')
(out / 'v4-link-check.json').write_text(json.dumps({'v5_identity': identity, 'compile_status': status, 'link_status': link_code,
                                                    'v4_sources_sha256': v4_files, 'utc': utc()}, indent=2) + '\n')
print(json.dumps({'compile': status, 'link': link_code}))
