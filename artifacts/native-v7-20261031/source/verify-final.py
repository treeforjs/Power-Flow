"""Native v7 final verification with receipts. Start restricted to cores 6-7:
  taskset -c 6,7 nice -n 15 env CUDA_VISIBLE_DEVICES= OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \\
      MKL_NUM_THREADS=1 python3 verify-final.py LABEL
Builds run under core-lease-20261030/both-cores.sh (build-components.py all); every CPU
test/benchmark process runs under core-lease-20261030/any-core.sh --wait (one core,
nice 15, one thread, no GPU) inside the frozen launcher, at most two at a time; the GPU
forms test runs under the two GPU locks (testing-gpu.lock, testing-gpu-<UUID>.lock).
Checks parity with the sealed native v6 evidence; writes OUT/LABEL/verification-result.json.
Never modifies native v6, its evidence or any sealed input."""
import sys
from concurrent.futures import ThreadPoolExecutor
from tooling import *

LEASE = Path('/home/treefor/projects/runtimes/core-lease-20261030')
label = sys.argv[1] if len(sys.argv) > 1 else 'final'
cpus = os.sched_getaffinity(0)
need(cpus <= {6, 7} and os.getpriority(os.PRIO_PROCESS, 0) >= 15, 'start under taskset -c 6,7 nice 15')
need(os.environ.get('CUDA_VISIBLE_DEVICES', None) == '', 'CUDA_VISIBLE_DEVICES must be empty')
F = OUT / label
F.mkdir(parents=True, exist_ok=False)
records = []
v6_before = tree_hashes(V6)
v6_out_before = tree_hashes(V6_OUT / 'final')
inputs_before = sealed_inputs()


def execute(name, argv, timeout=7200):
    log = F / (name + '.log')
    start = time.monotonic()
    with log.open('w') as h:
        try:
            rc = subprocess.run([str(a) for a in argv], stdout=h, stderr=subprocess.STDOUT, timeout=timeout,
                                cwd=str(P)).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    rec = {'name': name, 'argv': [str(a) for a in argv], 'returncode': rc, 'elapsed_s': round(time.monotonic() - start, 3),
           'log': str(log), 'log_sha256': sha256(log), 'utc': utc()}
    records.append(rec)
    (F / 'receipt-commands.json').write_text(json.dumps(records, indent=1) + '\n')
    print(name, rc, rec['elapsed_s'], flush=True)
    return rc


def cpu_run(name, argv, sanitized=False, rlimit=RLIMIT_AS, timeout=7200):
    env = ['env', 'CUDA_VISIBLE_DEVICES=']
    if sanitized:
        env += ['%s=%s' % kv for kv in ASAN_ENV.items()]
        prefix = [FROZEN]
    else:
        prefix = [FROZEN, 'prlimit', '--as=%d' % rlimit, '--']
    return execute(name, [LEASE / 'any-core.sh', '--wait', *env, *prefix, *argv], timeout)


def gpu_run(name, argv, timeout=1800):
    return execute(name, ['flock', LEASES / 'testing-gpu.lock', 'flock', LEASES / ('testing-gpu-%s.lock' % GPU_UUID),
                          LEASE / 'any-core.sh', '--wait', 'env', 'CUDA_VISIBLE_DEVICES=' + GPU_UUID, FROZEN, *argv], timeout)


# 1. Build everything (normal + ASan/UBSan + GPU forms test) under both leased cores.
rc = execute('build', [LEASE / 'both-cores.sh', 'env', 'CUDA_VISIBLE_DEVICES=', sys.executable, P / 'build-components.py', 'all'],
             timeout=7200)
need(rc == 0, 'build failed')
manifest = json.loads((P / 'SOURCE-SHA256.json').read_text())
identity = manifest['source_identity']
B = P / 'build'
names_ = ['legacy', 'ported', 'history', 'flows', 'source', 'forms-cpu', 'forms-native', 'work-gross', 'flows-long',
          'parity', 'bench', 'gpu', 'san-legacy', 'san-ported', 'san-history', 'san-flows', 'san-source', 'san-forms-cpu',
          'san-forms-native', 'san-work-gross']
d = {k: F / k for k in names_}
for x in d.values():
    x.mkdir()
G = 3 * 1024 ** 3
H2O = P / 'inputs/h2o-gas-20261029.json'
runs = [('legacy', [B / 'test_legacy', H2O, d['legacy']], False, G),
        ('ported', [B / 'test_ported', d['ported']], False, G),
        ('history', [B / 'test_history', d['history']], False, G),
        ('flows', [B / 'test_flows', d['flows']], False, G),
        ('source-arrival', [B / 'test_source_arrival', d['source']], False, G),
        ('forms-cpu', [B / 'test_forms_cpu', d['forms-cpu'], 1000000], False, G),
        ('forms-native', [B / 'test_forms_native', d['forms-native']], False, G),
        ('work-gross', [B / 'test_work_gross', d['work-gross']], False, G),
        ('flows-long-200000', [B / 'test_flows', d['flows-long'], 'long=200000'], False, G),
        ('parity-v7', [B / 'test_parity', d['parity'] / 'v7.json', 60], False, G),
        ('parity-v6', [B / 'v6_test_parity', d['parity'] / 'v6.json', 60], False, G),
        ('san-legacy', [B / 'sanitized-test_legacy', H2O, d['san-legacy']], True, None),
        ('san-ported', [B / 'sanitized-test_ported', d['san-ported']], True, None),
        ('san-history', [B / 'sanitized-test_history', d['san-history'], 'reduced'], True, None),
        ('san-flows', [B / 'sanitized-test_flows', d['san-flows'], 'reduced'], True, None),
        ('san-source-arrival', [B / 'sanitized-test_source_arrival', d['san-source'], 'reduced'], True, None),
        ('san-forms-cpu', [B / 'sanitized-test_forms_cpu', d['san-forms-cpu'], 100000], True, None),
        ('san-forms-native', [B / 'sanitized-test_forms_native', d['san-forms-native'], 'reduced'], True, None),
        ('san-work-gross', [B / 'sanitized-test_work_gross', d['san-work-gross']], True, None)]
with ThreadPoolExecutor(max_workers=2) as pool:
    futures = {name: pool.submit(cpu_run, name, argv, san, lim or G) for name, argv, san, lim in runs}
    codes = {name: f.result() for name, f in futures.items()}
need(all(c == 0 for c in codes.values()), 'test failures: %s' % {k: v for k, v in codes.items() if v != 0})
# Independent rational verification of the exact sums.
rc = cpu_run('verify-exact', ['/usr/bin/python', '-I', P / 'forms/tests/verify_exact.py', d['forms-cpu'] / 'adversarial.txt'])
need(rc == 0, 'rational verification failed')
# GPU: CPU vs GPU bit identity, two fresh processes.
for k in [1, 2]:
    rc = gpu_run('gpu-run%d' % k, [B / 'test_forms_gpu', NVRTC, GPU_UUID, 1000000])
    need(rc == 0, 'GPU forms test failed')
# Benchmarks one at a time (timing; the 1e6 runs need more address space).
for n in [10000, 100000, 1000000]:
    for tag, binary in [('v7', B / 'test_bench'), ('v6', B / 'v6_test_bench')]:
        rc = cpu_run('bench-%s-%d' % (tag, n), [binary, d['bench'] / ('%s-%d.json' % (tag, n)), n, 5 if n < 1000000 else 3],
                     rlimit=12 * 1024 ** 3)
        need(rc == 0, 'benchmark failed')

# 2. Parity and results.
load = lambda p: json.loads(Path(p).read_text())
V6F = V6_OUT / 'final'
legacy = load(d['legacy'] / 'arena-tests.json')
need(legacy == load(V6F / 'legacy/arena-tests.json') == load(d['san-legacy'] / 'arena-tests.json'),
     'legacy results/proofs differ from native v6')


def flat(x, p=''):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from flat(v, p + '/' + k)
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from flat(v, p + '/' + str(i))
    else:
        yield p, x


legacy_file_diffs = {}
for f in sorted((V6F / 'legacy').glob('*.json')):
    a, b = dict(flat(load(d['legacy'] / f.name))), dict(flat(load(f)))
    legacy_file_diffs[f.name] = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
need(all(set(v) <= {'/payload/source_sha256', '/sha256'} for v in legacy_file_diffs.values()),
     'legacy outputs differ beyond the source identity: %s' % legacy_file_diffs)
ported = load(d['ported'] / 'candidate005-ported-tests.json')
v6_ported = load(V6F / 'ported/candidate005-ported-tests.json')
names = lambda j: [x['name'] for x in j['checks']]
need(names(ported) == names(v6_ported) and ported == load(d['san-ported'] / 'candidate005-ported-tests.json') and
     all(x['passed'] for x in ported['checks']) and ported['allocation_failure_sweep'] == v6_ported['allocation_failure_sweep'],
     'ported controls differ from native v6')
g7 = load(d['ported'] / 'genesis.json')['payload']
g6 = load(V6F / 'ported/genesis.json')['payload']
genesis_diff = sorted(k for k in set(dict(flat(g7))) | set(dict(flat(g6))) if dict(flat(g7)).get(k) != dict(flat(g6)).get(k))
need(set(genesis_diff) <= {'/genesis/history/anchor_root_sha256', '/metadata/history/anchor_root_sha256', '/source_sha256'},
     'genesis checkpoint differs beyond source identity/anchor root: %s' % genesis_diff)
history = load(d['history'] / 'v5-history-tests.json')
v6_history = load(V6F / 'history/v5-history-tests.json')
new_history = [n for n in names(history) if n not in names(v6_history)]
removed_history = [n for n in names(v6_history) if n not in names(history)]
need(all(x['passed'] for x in history['checks']) and
     removed_history == ['Hall canonical census digest of the native census',
                         'complete motion record remains available to the caller for the step'],
     'history controls differ from v6 beyond the ported motion checks')
v6_flows = load(V6F / 'flows/v6-flows-tests.json')
flows_ported = {'removed_v6_names': [n for n in names(v6_flows) if n not in names(load(d['flows'] / 'v6-flows-tests.json'))],
                'new_v7_names': [n for n in names(load(d['flows'] / 'v6-flows-tests.json')) if n not in names(v6_flows)]}
san_history = load(d['san-history'] / 'v5-history-tests.json')
flows = load(d['flows'] / 'v6-flows-tests.json')
san_flows = load(d['san-flows'] / 'v6-flows-tests.json')
normal_only = {'long run: resident memory bounded after warm-up (< 8 MiB growth)'}
need(flows['status'] == 'passed' and san_flows['status'] == 'passed' and
     [n for n in names(flows) if n not in normal_only] == names(san_flows), 'flows normal/sanitized check lists differ')
source = load(d['source'] / 'v6-source-arrival-tests.json')
san_source = load(d['san-source'] / 'v6-source-arrival-tests.json')
v6_source = load(V6F / 'source/v6-source-arrival-tests.json')
need(source['status'] == san_source['status'] == 'passed' and source['births'] == v6_source['births'] and
     source['exact_stream_equalities'] == v6_source['exact_stream_equalities'], 'source/arrival sweep differs from v6')
forms_native = load(d['forms-native'] / 'v7-forms-native.json')
work = load(d['work-gross'] / 'v7-work-gross.json')
san_work = load(d['san-work-gross'] / 'v7-work-gross.json')
need(work['status'] == san_work['status'] == 'passed', 'work tolerance controls failed')
long_run = load(d['flows-long'] / 'v6-flows-long-run.json')
p7, p6 = load(d['parity'] / 'v7.json'), load(d['parity'] / 'v6.json')
step_diff = sorted({k for x, y in zip(p7['steps'], p6['steps']) for k in x if x[k] != y[k]})
need(p7['genesis'] == p6['genesis'] and len(p7['steps']) == len(p6['steps']) == 60 and
     all(x['physical_state_sha256'] == y['physical_state_sha256'] for x, y in zip(p7['steps'], p6['steps'])) and
     set(step_diff) <= {'accepted_arena_sha256', 'anchor_head_sha256', 'cumulative_bytes'},
     'physical-state parity trace differs: %s' % step_diff)


def results_of(log):
    return [line.strip() for line in log.read_text().splitlines() if line.startswith('RESULT ')]


gpu1, gpu2 = results_of(F / 'gpu-run1.log'), results_of(F / 'gpu-run2.log')
strip_timing = lambda lines: [' '.join(t for t in l.split() if 'ms' not in t and 'wall' not in t and 'us_per' not in t) for l in lines]
need(strip_timing(gpu1) == strip_timing(gpu2) and any('gpu_checks=' in l and 'failures=0' in l for l in gpu1),
     'GPU runs differ or failed')
bench = {}
for n in [10000, 100000, 1000000]:
    for tag in ['v7', 'v6']:
        j = load(d['bench'] / ('%s-%d.json' % (tag, n)))
        keys = sorted(k for k in j['runs'][0] if k.endswith('_s'))
        med = {k: sorted(r[k] for r in j['runs'])[len(j['runs']) // 2] for k in keys}
        bench['%s-%d' % (tag, n)] = {'median_s': med, 'VmHWM_KiB': j['VmHWM_KiB'], 'repeats': len(j['runs'])}
need(tree_hashes(V6) == v6_before and tree_hashes(V6_OUT / 'final') == v6_out_before, 'native v6 or its evidence changed')
need(sealed_inputs() == inputs_before, 'sealed inputs changed during verification')
need(json.loads((P / 'SOURCE-SHA256.json').read_text()) == manifest, 'v7 source changed during verification')
artifacts = {str(p.relative_to(P)): sha256(p) for p in sorted(B.glob('*')) if p.is_file()}
(F / 'build-artifact-hashes.json').write_text(json.dumps(artifacts, indent=2) + '\n')
result = {
    'status': 'passed',
    'source_identity': identity,
    'base_v6_identity': V6_IDENTITY,
    'legacy_checks': legacy['check_count'],
    'legacy_arena_tests_identical_to_native_v6_and_sanitized': True,
    'legacy_output_differences_vs_v6': legacy_file_diffs,
    'ported_candidate003_controls': {'check_count': ported['check_count'], 'identical_names_to_v6': True,
                                     'allocation_failure_sweep': ported['allocation_failure_sweep']},
    'v5_history_controls': {'normal': history['check_count'], 'sanitized_reduced': san_history['check_count'],
                            'new_v7_checks': new_history, 'replaced_v6_checks': removed_history},
    'v6_flows_name_changes': flows_ported,
    'genesis_checkpoint_differences_vs_v6': genesis_diff,
    'physical_state_parity_trace': {'steps': 60, 'anchor_interval': 4, 'physical_state_identical': True,
                                    'differing_identity_fields': step_diff},
    'v6_flows_ported': {'normal': flows['check_count'], 'sanitized_reduced': san_flows['check_count'],
                        'normal_only': sorted(normal_only), 'measurements': flows['measurements']},
    'v6_source_arrival': {k: source[k] for k in ['check_count', 'cases', 'admitted_steps', 'births', 'exact_stream_equalities',
                                                  'source_refused', 'arrival_refused', 'exposure_refused', 'surface_births',
                                                  'interior_births', 'active_layer_births', 'multi_step_committed_cases']},
    'v6_source_arrival_sanitized_reduced': {k: san_source[k] for k in ['check_count', 'exact_stream_equalities', 'births']},
    'forms_cpu': results_of(F / 'forms-cpu.log'),
    'forms_cpu_sanitized': results_of(F / 'san-forms-cpu.log')[-1:],
    'exact_rational_verification': results_of(F / 'verify-exact.log'),
    'forms_native': {'check_count': forms_native['check_count'], 'measurements': forms_native['measurements']},
    'work_gross_v1': {'check_count': work['check_count'], 'sanitized': san_work['check_count'], 'cases': work['cases']},
    'gpu': gpu1,
    'long_run_200000': {k: v for k, v in long_run['measurements']['long_run'].items()},
    'benchmark': bench,
    'sealed_inputs': inputs_before,
    'resources': {'orchestrator_cpus': sorted(cpus), 'build': 'both-cores.sh (cores 6,7), two compiles at a time',
                  'runs': 'any-core.sh --wait (one core each, at most two concurrent), nice 15, one thread, no GPU',
                  'gpu': 'testing-gpu.lock + testing-gpu-%s.lock, any-core.sh --wait' % GPU_UUID,
                  'rlimit_as_tests': G, 'rlimit_as_benchmarks': 12 * 1024 ** 3, 'launcher': FROZEN},
    'utc': utc()}
(F / 'verification-result.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: result[k] for k in ['status', 'source_identity', 'legacy_checks', 'ported_candidate003_controls',
                                          'v5_history_controls', 'physical_state_parity_trace']}, indent=1))
