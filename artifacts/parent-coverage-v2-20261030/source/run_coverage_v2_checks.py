#!/usr/bin/env python3
"""Pinned qualification of parent-channel coverage scope v2.

Run under a core lease wrapper, which holds the lease of the core(s) it binds:
  core-lease-20261030/both-cores.sh python3 run_coverage_v2_checks.py   (builds two at a time)
  core-lease-20261030/any-core.sh --wait python3 run_coverage_v2_checks.py
Builds and runs, inside the frozen launcher, one thread per process:
  1. test_parent_coverage_v2 (normal -O2 and UBSan);
  2. the v2 port of the sealed v3 parent regression without high-energy bands
     (must equal the sealed baseline: 442 checks, event physics, 73/128);
  3. the v2 port of the sealed v3 parent regression with all v2 bands (the 8 keV
     entry removed: 440 checks equal to the baseline's, event physics, 73/128);
  4. the v2 port of the coverage-v1 component test (171 checks equal to the
     sealed v1 evidence, ledgers equal apart from the added categories);
  5. the dataset-reference outer3 program (dataset include tree), compared bit
     for bit with the v2 outer3 evaluation;
  6. the 50-digit Decimal BEB/RBEB oracle, compared with the C++ values;
  7. the short exposure benchmark (v1 vs v2).
Receipts: outputs/pic-water-parent-coverage-v2-20261030/attempt-NNN/.
"""
import hashlib, json, os, resource, subprocess, sys, time
from decimal import Decimal
from pathlib import Path

P = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
B = Path('/home/treefor/projects/.workspace/outputs')
O = B / 'pic-water-parent-coverage-v2-20261030'
V1 = R / 'pic-water-parent-coverage-v1-20261030'
V1A = Path('/home/treefor/projects/Power-Flow/artifacts/parent-coverage-v1-20261030')
V3 = R / 'pic-water-common-commit-v3-20261030'
N = R / 'pic-reaction-transactions-20261029/native-full-time-candidate003'
I = R / 'pic-water-full-time-subflow-20261029/inputs'
M = R / 'mfem-4.10-pic-20261029'
DS = R / 'h2o-collision-data-20261029'
F = R / 'arch-base-20261028/run-frozen'
identity = 'c045c42d2f57082e3fb5e12172c440d2470e681a38df9ac406bd7e063609a08d'

if not __debug__ or os.getenv('PYTHONOPTIMIZE'):
    raise SystemExit('Optimized guards refused')
core = os.getenv('PF_CORE')
allowed = {'6': {6}, '7': {7}, '6,7': {6, 7}}
if core not in allowed or os.sched_getaffinity(0) != allowed[core] or os.getpriority(os.PRIO_PROCESS, 0) != 15:
    raise SystemExit('run under any-core.sh or both-cores.sh: leased Power-Flow core(s) 6/7 at nice 15 required')
parallel = 2 if core == '6,7' else 1
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
outer3 = DS / 'data/collisions/h2o-gas-20261029-outer3.json'
binding = baseline / 'run-binding.json'
kernel_binding = P / 'kernel-binding-v2.json'
if sha(data) != '840663cdbcf4c39f1ab6f5407f9ecacb232a07c3c48002a4a182c8187792e4bf':
    raise RuntimeError('Water dataset differs from h2o-gas-20261029-v1')
if sha(outer3) != 'c436020d3f2d07902fb059cfc0a269b2e914a3cbb2789300093d3631df896551':
    raise RuntimeError('Outer3 dataset differs from h2o-gas-20261029-outer3-v5')
# The v2 build includes the v1 runtime sources (exact_sum.hpp and, in tests, the
# v1 subflow); they must equal the sealed v1 artifact copy.
for name in ('exact_sum.hpp', 'parent_coverage_scope.hpp', 'parent_coverage_subflow.hpp'):
    if sha(V1 / 'source' / name) != sha(V1A / 'source/source' / name):
        raise RuntimeError('v1 runtime source differs from the sealed artifact: ' + name)
v1_evidence = V1A / 'evidence/attempt-001/coverage-normal/result.json'
libs = [M / x for x in ('install-mpi/lib/libmfem.a', 'tpl/hypre/lib/libHYPRE.a', 'tpl/metis/lib/libmetis.a', 'tpl/GKlib/lib/libGKlib.a')]
trees = [P / 'source', P / 'tests', P / 'tools', V1 / 'source', V3 / 'source', V3 / 'tests', N / 'source/include', N / 'inputs', I,
         DS / 'include', M / 'install-mpi/include', M / 'tpl/hypre/include', M / 'tpl/metis/include']


def pins():
    files = {p for tree in trees for p in tree.rglob('*') if p.is_file()}
    files.update([Path(__file__).resolve(), P / 'V3-TEST-SHA256', P / 'V1-TEST-SHA256', kernel_binding, F, lib, *libs, data, outer3,
                  binding, baseline / 'result.json', v1_evidence, V1A / 'source/tests/test_parent_coverage.cpp',
                  R / 'pic-water-exact-count-subflow-20261029/kernel-binding.json', N / 'SOURCE-SHA256.json',
                  R / 'arch-base-20261028/root/usr/bin/mpic++', R / 'arch-base-20261028/root/usr/bin/g++'])
    return {str(p): sha(p) for p in sorted(files)}


before = pins()
O.mkdir(exist_ok=True)
(P / 'build').mkdir(exist_ok=True)
A = O / f'attempt-{len(list(O.glob("attempt-*"))) + 1:03d}'
A.mkdir()
(A / 'inputs.json').write_text(json.dumps(before, indent=2) + '\n')
common = ['-std=c++20', '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-ffp-contract=off', '-fno-fast-math',
          '-ffunction-sections', '-fdata-sections']
flags = common + ['-I' + str(P / 'source'), '-I' + str(V1 / 'source'), '-I' + str(V3 / 'source'),
                  '-DPIC_REACTION_NATIVE_SOURCE_SHA256="' + identity + '"']
for d in (N / 'source/include', I, I / 'include', N / 'inputs', M / 'install-mpi/include', M / 'tpl/hypre/include', M / 'tpl/metis/include'):
    flags += ['-isystem', str(d)]
link = [str(lib), *map(str, libs), '-Wl,--gc-sections', '-lcrypto']
variants = {'normal': ['-O2'], 'ubsan': ['-O1', '-g', '-fsanitize=undefined', '-fno-sanitize-recover=all']}
tests = {'coverage_v2': P / 'tests/test_parent_coverage_v2.cpp', 'v3_port_a': P / 'tests/test_v3_parent_regression_port_v2a.cpp',
         'v3_port_b': P / 'tests/test_v3_parent_regression_port_v2b.cpp', 'v1_port': P / 'tests/test_v1_coverage_port_v2.cpp'}
builds, runs, outdirs = [], [], {}
for t, src in tests.items():
    for v, opt in variants.items():
        exe = P / f'build/{src.stem}-{v}'
        out = A / f'{t}-{v}'
        out.mkdir()
        builds.append([str(F), '/usr/bin/mpic++', *flags, *opt, str(src), *link, '-o', str(exe)])
        args = {'coverage_v2': [str(data), str(binding), str(outer3), str(kernel_binding), str(out)],
                'v3_port_a': [str(data), str(binding), str(out), identity], 'v3_port_b': [str(data), str(binding), str(out), identity],
                'v1_port': [str(data), str(binding), str(out)]}[t]
        runs.append([str(F), str(exe), *args])
        outdirs[(t, v)] = out
ref_exe = P / 'build/dataset_reference_outer3-normal'
builds.append([str(F), '/usr/bin/mpic++', *common, '-isystem', str(DS / 'include'), '-O2', str(P / 'tests/dataset_reference_outer3.cpp'),
               '-o', str(ref_exe)])
runs.append([str(F), str(ref_exe), str(outer3), str(A / 'dataset-reference-outer3.json')])
bench_exe = P / 'build/bench_exposure_v2-normal'
builds.append([str(F), '/usr/bin/mpic++', *flags, '-O2', str(P / 'tests/bench_exposure_v2.cpp'), *link, '-o', str(bench_exe)])
runs.append([str(F), str(bench_exe), str(data), str(binding), str(outer3), str(kernel_binding), str(A)])
runs.append([str(F), '/usr/bin/python3', '-I', str(P / 'tools/oracle_beb_rbeb.py'), str(outer3), str(A / 'oracle-beb-rbeb.json')])
(A / 'commands.json').write_text(json.dumps({'builds': builds, 'runs': runs, 'environment': env, 'core': core,
                                             'parallel_builds': parallel}, indent=2) + '\n')


start = time.monotonic()
status, error, comparison, timing = 'failed', None, {}, {}
try:
    for script in ('tools/make_kernel_binding.py', 'tools/make_ports.py'):
        r = subprocess.run([sys.executable, '-I', str(P / script), '--check'], env=env, capture_output=True, text=True, timeout=60)
        if r.returncode != 0:
            raise RuntimeError(script + ' check failed: ' + r.stdout + r.stderr)
    t0 = time.monotonic()
    pending = list(enumerate(builds))
    while pending:
        batch, pending = pending[:parallel], pending[parallel:]
        procs = []
        for index, command in batch:
            log = (A / f'build-{index}.log').open('w')
            procs.append((index, subprocess.Popen(command, env=env, cwd=P, stdout=log, stderr=subprocess.STDOUT), log))
        for index, proc, log in procs:
            rc = proc.wait(timeout=1800);log.close()
            if rc != 0:
                raise RuntimeError(f'build {index} failed; see build-{index}.log')
    timing['builds_s'] = time.monotonic() - t0
    for index, command in enumerate(runs):
        t1 = time.monotonic()
        with (A / f'run-{index}.log').open('w') as log:
            subprocess.run(command, env=env, cwd=P, stdout=log, stderr=subprocess.STDOUT, timeout=1800, check=True)
        timing[f'run_{index}_s'] = time.monotonic() - t1

    old = json.loads((baseline / 'result.json').read_text())
    events = lambda j: [[{k: v for k, v in e.items() if k != 'native'} for e in r['events']] for r in j['receipts']]
    res = lambda t, v: json.loads((outdirs[(t, v)] / 'result.json').read_text())
    # 1. coverage v2 component checks
    for v in variants:
        c = res('coverage_v2', v)
        if c['status'] != 'passed' or c['physical_trajectory_advanced'] or c['full_coupled_run'] or not all(x['passed'] for x in c['checks']):
            raise RuntimeError('Coverage v2 checks failed: ' + v)
        comparison[f'coverage_v2_{v}'] = {'checks': c['check_count'], 'all_passed': True, 'coverage_scope_sha256': c['coverage_scope_sha256'],
                                          'arena_channel_binding': c['arena_channel_binding']}
    a, u = res('coverage_v2', 'normal'), res('coverage_v2', 'ubsan')
    comparison['coverage_v2_normal_ubsan_results_identical'] = a == u
    if a != u:
        raise RuntimeError('Normal and UBSan coverage v2 results differ')
    # 2./3. v3 regression ports
    names = [c['name'] for c in old['checks']]
    pair = ['unsupported incident energy cannot be silentlyignored', 'support error leaves all owners unchanged']
    at = [i for i in range(len(names) - 1) if names[i:i + 2] == pair]
    if len(at) != 2:
        raise RuntimeError('baseline unsupported-energy loop not found')
    trimmed = old['checks'][:at[1]] + old['checks'][at[1] + 2:]
    for t, expected, count in (('v3_port_a', old['checks'], 442), ('v3_port_b', trimmed, 440)):
        for v in variants:
            port = res(t, v)
            c = {'checks': port['check_count'], 'first_events': port['first_events'],
                 'all_check_results_equal': port['checks'] == expected, 'recorded_event_physics_equal': events(port) == events(old),
                 'recorded_event_counts': [len(x['events']) for x in port['receipts']],
                 'physical_trajectory_advanced': port['physical_trajectory_advanced'], 'full_coupled_run': port['full_coupled_run']}
            if (port['status'] != 'passed' or c['checks'] != count or c['first_events'] != old['first_events'] or
                    not c['all_check_results_equal'] or not c['recorded_event_physics_equal'] or
                    c['physical_trajectory_advanced'] or c['full_coupled_run']):
                raise RuntimeError(f'{t} differs from the sealed baseline: {v}')
            comparison[f'{t}_{v}'] = c
    comparison['v3_port_b_removed_baseline_entries'] = {'indices': [at[1], at[1] + 1], 'names': pair,
                                                        'reason': '8000 eV is covered in v2 (dataset outer3 BEB band)'}
    # 4. v1 component test port
    v1 = json.loads(v1_evidence.read_text())
    strip = lambda x: ({k: strip(w) for k, w in x.items() if k != 'categories'} if isinstance(x, dict) else
                       [strip(w) for w in x] if isinstance(x, list) else x)
    for v in variants:
        port = res('v1_port', v)
        c = {'checks': port['check_count'], 'all_check_results_equal': port['checks'] == v1['checks'],
             'ledger_receipts_equal_without_categories': strip(port['receipts']) == strip(v1['receipts']),
             'coverage_bins_equal': port['coverage_bins'] == v1['coverage_bins']}
        if port['status'] != 'passed' or c['checks'] != 171 or not all(c[k] for k in c if k != 'checks'):
            raise RuntimeError('v1 component port differs: ' + v)
        comparison[f'v1_port_{v}'] = c
    # 5. dataset reference, bit for bit
    ref = json.loads((A / 'dataset-reference-outer3.json').read_text())
    mine = json.loads((outdirs[('coverage_v2', 'normal')] / 'outer3-grid.json').read_text())
    same = len(ref) == len(mine) and all(x['energy'] == y['energy'] and x['parent'] == y['parent'] and x['valence'] == y['valence']
                                         for x, y in zip(ref, mine))
    comparison['outer3_dataset_reference_bitwise'] = {'rows': len(ref), 'equal': same}
    if not same:
        raise RuntimeError('outer3 rate differs from the dataset implementation')
    # 6. Decimal oracle
    orc = json.loads((A / 'oracle-beb-rbeb.json').read_text())['rows']
    cpp = a['report']['oracle_values']
    worst = {'rbeb_orbital': 0.0, 'beb_extended_orbital': 0.0, 'beb_dataset_arithmetic_orbital': 0.0, 'sums': 0.0}
    rel = lambda x, d: 0.0 if Decimal(d) == 0 and x == 0 else float(abs(Decimal(repr(x)) / Decimal(d) - 1)) if Decimal(d) != 0 else float('inf')
    if len(orc) != len(cpp):
        raise RuntimeError('oracle grid differs')
    for o, c in zip(orc, cpp):
        if o['energy_eV'] != c['energy_eV']:
            raise RuntimeError('oracle energies differ')
        for po, pcpp in zip(o['orbitals'], c['orbitals']):
            worst['rbeb_orbital'] = max(worst['rbeb_orbital'], rel(pcpp['rbeb_m2'], po['rbeb_m2']))
            worst['beb_extended_orbital'] = max(worst['beb_extended_orbital'], rel(pcpp['beb_extended_m2'], po['beb_m2']))
            worst['beb_dataset_arithmetic_orbital'] = max(worst['beb_dataset_arithmetic_orbital'], rel(pcpp['beb_dataset_arithmetic_m2'], po['beb_m2']))
        for k in ('beb_parent_m2', 'rbeb_parent_m2', 'beb_valence_m2', 'rbeb_valence_m2'):
            worst['sums'] = max(worst['sums'], rel(c[k], o[k]))
    comparison['decimal_oracle_worst_relative'] = worst
    if worst['rbeb_orbital'] > 2e-15 or worst['beb_extended_orbital'] > 2e-15 or worst['sums'] > 2e-15 or worst['beb_dataset_arithmetic_orbital'] > 1e-14:
        raise RuntimeError('C++ BEB/RBEB values differ from the Decimal oracle')
    # 7. benchmark
    comparison['benchmark'] = json.loads((A / 'benchmark.json').read_text())
    (A / 'coverage-scope-v2.json').write_text(json.dumps(a['coverage_scope'], indent=2) + '\n')
    (A / 'physics-report.json').write_text(json.dumps(a['report'], indent=2) + '\n')
    if pins() != before:
        raise RuntimeError('Dependencies changed during the checks')
    status = 'passed'
except Exception as exc:
    error = repr(exc)
receipt = {'status': status, 'error': error, 'comparison': comparison, 'timing_s': timing, 'elapsed_s': time.monotonic() - start,
           'native_source_identity': identity, 'cpu': core, 'nice': 15, 'threads_per_process': 1, 'parallel_builds': parallel,
           'lease': 'both-cores.sh' if core == '6,7' else 'any-core.sh', 'GPU_used': False, 'full_time_joint_executed': False,
           'physical_time_advanced_s': 0,
           'scope': 'parent-channel coverage scope v2 component checks, v3 parent regression and v1 component ports through the '
                    'coverage-v2 subflow, dataset-reference and Decimal-oracle comparisons, exposure benchmark; frozen-stage operator, '
                    'no Hall coupling',
           'input_manifest_sha256': sha(A / 'inputs.json')}
(A / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({'attempt': str(A), 'status': status, 'error': error}))
if status != 'passed':
    raise SystemExit(1)
