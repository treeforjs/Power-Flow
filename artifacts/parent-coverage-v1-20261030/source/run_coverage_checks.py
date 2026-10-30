#!/usr/bin/env python3
"""Pinned qualification of parent-channel coverage scope v1.

Run under the core lease wrapper, which holds the lease of the core it binds:
  core-lease-20261030/any-core.sh --wait python3 run_coverage_checks.py
Builds and runs, inside the frozen launcher, one thread each:
  1. test_parent_coverage (normal -O2 and UBSan),
  2. the generated coverage-v1 port of the sealed v3 parent regression
     (normal and UBSan), compared with the sealed baseline exactly as
     run_parent_regression.py compares it.
"""
import hashlib, json, os, resource, subprocess, sys, time
from pathlib import Path

P = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
B = Path('/home/treefor/sim-workspace/powerflow/outputs')
O = B / 'pic-water-parent-coverage-v1-20261030'
V3 = R / 'pic-water-common-commit-v3-20261030'
N = R / 'pic-reaction-transactions-20261029/native-full-time-candidate003'
I = R / 'pic-water-full-time-subflow-20261029/inputs'
M = R / 'mfem-4.10-pic-20261029'
F = R / 'arch-base-20261028/run-frozen'
identity = 'c045c42d2f57082e3fb5e12172c440d2470e681a38df9ac406bd7e063609a08d'

if not __debug__ or os.getenv('PYTHONOPTIMIZE'):
    raise SystemExit('Optimized guards refused')
core = os.getenv('PF_CORE')
if core not in ('6', '7') or os.sched_getaffinity(0) != {int(core)} or os.getpriority(os.PRIO_PROCESS, 0) != 15:
    raise SystemExit('run under any-core.sh: one leased Power-Flow core (6 or 7) at nice 15 required')
resource.setrlimit(resource.RLIMIT_AS, (3 * 1024**3, 3 * 1024**3))
env = {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
       'BLIS_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1', 'OMP_THREAD_LIMIT': '1',
       'CUDA_VISIBLE_DEVICES': '', 'PYTHONNOUSERSITE': '1', 'UBSAN_OPTIONS': 'print_stacktrace=1:halt_on_error=1'}
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

lib = N / 'build/libpic_native_candidate003.a'
if sha(lib) != 'd2568af9fac13f8aa056ae4a9a50bb2a07eee05d7218d1c4b79aef21317d8a20':
    raise RuntimeError('Unqualified native library')
manifest = json.loads((N / 'SOURCE-SHA256.json').read_text())
if manifest['source_identity'] != identity:
    raise RuntimeError('Native source identity differs')
for name, h in manifest['files'].items():
    if sha(N / name) != h:
        raise RuntimeError('Native source file changed: ' + name)
baseline = B / 'pic-water-exact-count-subflow-20261029'
data = R / 'pic-water-exact-count-subflow-20261029/inputs/h2o-gas.json'
binding = baseline / 'run-binding.json'
if sha(data) != '840663cdbcf4c39f1ab6f5407f9ecacb232a07c3c48002a4a182c8187792e4bf':
    raise RuntimeError('Water dataset differs from h2o-gas-20261029-v1')
libs = [M / x for x in ('install-mpi/lib/libmfem.a', 'tpl/hypre/lib/libHYPRE.a', 'tpl/metis/lib/libmetis.a', 'tpl/GKlib/lib/libGKlib.a')]
trees = [P / 'source', P / 'tests', V3 / 'source', V3 / 'tests', N / 'source/include', N / 'inputs', I,
         M / 'install-mpi/include', M / 'tpl/hypre/include', M / 'tpl/metis/include']


def pins():
    files = {p for tree in trees for p in tree.rglob('*') if p.is_file()}
    files.update([Path(__file__).resolve(), P / 'make_v3_port.py', P / 'V3-TEST-SHA256', F, lib, *libs, data, binding,
                  baseline / 'result.json', N / 'SOURCE-SHA256.json', V3 / 'run_parent_regression.py',
                  R / 'arch-base-20261028/root/usr/bin/mpic++', R / 'arch-base-20261028/root/usr/bin/g++'])
    return {str(p): sha(p) for p in sorted(files)}


before = pins()
O.mkdir(exist_ok=True)
(P / 'build').mkdir(exist_ok=True)
A = O / f'attempt-{len(list(O.glob("attempt-*"))) + 1:03d}'
A.mkdir()
(A / 'inputs.json').write_text(json.dumps(before, indent=2) + '\n')
flags = ['-std=c++20', '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-ffp-contract=off', '-fno-fast-math',
         '-ffunction-sections', '-fdata-sections', '-I' + str(P / 'source'), '-I' + str(V3 / 'source'),
         '-DPIC_REACTION_NATIVE_SOURCE_SHA256="' + identity + '"']
for d in (N / 'source/include', I, I / 'include', N / 'inputs', M / 'install-mpi/include', M / 'tpl/hypre/include', M / 'tpl/metis/include'):
    flags += ['-isystem', str(d)]
variants = {'normal': ['-O2'], 'ubsan': ['-O1', '-g', '-fsanitize=undefined', '-fno-sanitize-recover=all']}
tests = {'coverage': P / 'tests/test_parent_coverage.cpp', 'v3_port': P / 'tests/test_v3_parent_regression_port.cpp'}
commands, runs = [], {}
for t, src in tests.items():
    for v, opt in variants.items():
        exe = P / f'build/{src.stem}-{v}'
        out = A / f'{t}-{v}'
        out.mkdir()
        commands.append([str(F), '/usr/bin/mpic++', *flags, *opt, str(src), str(lib), *map(str, libs),
                         '-Wl,--gc-sections', '-lcrypto', '-o', str(exe)])
        run = [str(F), str(exe), str(data), str(binding), str(out)] + ([identity] if t == 'v3_port' else [])
        commands.append(run)
        runs[(t, v)] = out
(A / 'commands.json').write_text(json.dumps({'commands': commands, 'environment': env, 'core': int(core)}, indent=2) + '\n')

start = time.monotonic()
status, error, comparison = 'failed', None, {}
try:
    port_check = subprocess.run([sys.executable, '-I', str(P / 'make_v3_port.py'), '--check'], env=env,
                                capture_output=True, text=True, timeout=30)
    if port_check.returncode != 0:
        raise RuntimeError('Port differs from its generator: ' + port_check.stdout + port_check.stderr)
    for index, command in enumerate(commands):
        with (A / f'command-{index}.log').open('w') as log:
            subprocess.run(command, env=env, cwd=P, stdout=log, stderr=subprocess.STDOUT,
                           timeout=600, check=True)
    old = json.loads((baseline / 'result.json').read_text())

    def events(j):
        return [[{k: v for k, v in e.items() if k != 'native'} for e in r['events']] for r in j['receipts']]

    for v in variants:
        cov = json.loads((runs[('coverage', v)] / 'result.json').read_text())
        if cov['status'] != 'passed' or cov['physical_trajectory_advanced'] or cov['full_coupled_run']:
            raise RuntimeError('Coverage checks failed: ' + v)
        comparison[f'coverage_{v}'] = {'checks': cov['check_count'], 'all_passed': all(c['passed'] for c in cov['checks']),
                                       'coverage_scope_sha256': cov['coverage_scope_sha256']}
        port = json.loads((runs[('v3_port', v)] / 'result.json').read_text())
        c = {'checks': port['check_count'], 'first_events': port['first_events'],
             'all_check_results_equal': port['checks'] == old['checks'],
             'recorded_event_physics_equal': events(port) == events(old),
             'recorded_event_counts': [len(x['events']) for x in port['receipts']],
             'physical_trajectory_advanced': port['physical_trajectory_advanced'], 'full_coupled_run': port['full_coupled_run']}
        if (port['status'] != 'passed' or c['checks'] != 442 or c['first_events'] != old['first_events'] or
                not c['all_check_results_equal'] or not c['recorded_event_physics_equal'] or
                c['physical_trajectory_advanced'] or c['full_coupled_run']):
            raise RuntimeError('Coverage-v1 port of the v3 parent regression differs: ' + v)
        comparison[f'v3_regression_port_{v}'] = c
    a = json.loads((runs[('coverage', 'normal')] / 'result.json').read_text())
    u = json.loads((runs[('coverage', 'ubsan')] / 'result.json').read_text())
    comparison['coverage_normal_ubsan_results_identical'] = a == u
    if a != u:
        raise RuntimeError('Normal and UBSan coverage results differ')
    (A / 'coverage-scope.json').write_text(json.dumps(a['coverage_scope'], indent=2) + '\n')
    if pins() != before:
        raise RuntimeError('Dependencies changed during the checks')
    status = 'passed'
except Exception as exc:
    error = repr(exc)
receipt = {'status': status, 'error': error, 'comparison': comparison, 'elapsed_s': time.monotonic() - start,
           'native_source_identity': identity, 'cpu': int(core), 'nice': 15, 'threads': 1, 'lease': 'any-core.sh',
           'GPU_used': False, 'full_time_joint_executed': False, 'physical_time_advanced_s': 0,
           'scope': 'parent-channel coverage scope v1 component checks and the v3 parent regression through the '
                    'coverage-v1 subflow; frozen-stage operator, no Hall coupling',
           'input_manifest_sha256': sha(A / 'inputs.json')}
(A / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({'attempt': str(A), **receipt}))
if status != 'passed':
    raise SystemExit(1)
