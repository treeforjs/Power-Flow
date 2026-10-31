"""Native v6 final verification with receipts. Start restricted to cores 6-7:
  taskset -c 6,7 nice -n 15 env CUDA_VISIBLE_DEVICES= OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
      MKL_NUM_THREADS=1 python3 verify-final.py LABEL
Builds run under core-lease-20261030/both-cores.sh (build-components.py all); every
test/benchmark process runs under core-lease-20261030/any-core.sh --wait (one core,
nice 15, one thread, no GPU) inside the frozen launcher, at most two at a time.
Checks parity with sealed native v5 evidence; writes OUT/LABEL/verification-result.json.
Never modifies native v5, its evidence or any sealed input."""
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
v5_before = tree_hashes(V5)
v5_out_before = tree_hashes(V5_OUT / 'final')
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


# 1. Build everything (normal + ASan/UBSan) under both leased cores.
rc = execute('build', [LEASE / 'both-cores.sh', 'env', 'CUDA_VISIBLE_DEVICES=', sys.executable, P / 'build-components.py', 'all'],
             timeout=7200)
need(rc == 0, 'build failed')
manifest = json.loads((P / 'SOURCE-SHA256.json').read_text())
identity = manifest['source_identity']
B = P / 'build'
d = {k: F / k for k in ['legacy', 'ported', 'history', 'flows', 'source', 'flows-long', 'parity', 'bench',
                         'san-legacy', 'san-ported', 'san-history', 'san-flows', 'san-source']}
for x in d.values():
    x.mkdir()
G = 3 * 1024 ** 3
runs = [('legacy', [B / 'test_legacy', P / 'inputs/h2o-gas-20261029.json', d['legacy']], False, G),
        ('ported', [B / 'test_ported', d['ported']], False, G),
        ('history', [B / 'test_history', d['history']], False, G),
        ('flows', [B / 'test_flows', d['flows']], False, G),
        ('source-arrival', [B / 'test_source_arrival', d['source']], False, G),
        ('flows-long-200000', [B / 'test_flows', d['flows-long'], 'long=200000'], False, G),
        ('parity-v6', [B / 'test_parity', d['parity'] / 'v6.json', 60], False, G),
        ('parity-v5', [B / 'v5_test_parity', d['parity'] / 'v5.json', 60], False, G),
        ('san-legacy', [B / 'sanitized-test_legacy', P / 'inputs/h2o-gas-20261029.json', d['san-legacy']], True, None),
        ('san-ported', [B / 'sanitized-test_ported', d['san-ported']], True, None),
        ('san-history', [B / 'sanitized-test_history', d['san-history'], 'reduced'], True, None),
        ('san-flows', [B / 'sanitized-test_flows', d['san-flows'], 'reduced'], True, None),
        ('san-source-arrival', [B / 'sanitized-test_source_arrival', d['san-source'], 'reduced'], True, None)]
with ThreadPoolExecutor(max_workers=2) as pool:
    futures = {name: pool.submit(cpu_run, name, argv, san, lim or G) for name, argv, san, lim in runs}
    codes = {name: f.result() for name, f in futures.items()}
need(all(c == 0 for c in codes.values()), 'test failures: %s' % {k: v for k, v in codes.items() if v != 0})
# Benchmarks one at a time (timing; the 1e6 runs need more address space).
for n in [10000, 100000, 1000000]:
    for tag, binary in [('v6', B / 'test_bench'), ('v5', B / 'v5_test_bench')]:
        rc = cpu_run('bench-%s-%d' % (tag, n), [binary, d['bench'] / ('%s-%d.json' % (tag, n)), n, 3 if n < 1000000 else 2],
                     rlimit=12 * 1024 ** 3)
        need(rc == 0, 'benchmark failed')

# 2. Parity and results.
load = lambda p: json.loads(Path(p).read_text())
legacy = load(d['legacy'] / 'arena-tests.json')
need(legacy == load(V5_OUT / 'final/legacy/arena-tests.json') == load(d['san-legacy'] / 'arena-tests.json'),
     'legacy results/proofs differ from native v5')
ported = load(d['ported'] / 'candidate005-ported-tests.json')
v5_ported = load(V5_OUT / 'final/ported/candidate005-ported-tests.json')
names = lambda j: [x['name'] for x in j['checks']]
need(names(ported) == names(v5_ported) and ported == load(d['san-ported'] / 'candidate005-ported-tests.json') and
     all(x['passed'] for x in ported['checks']) and ported['allocation_failure_sweep'] == v5_ported['allocation_failure_sweep'],
     'ported controls differ from native v5')
g6 = load(d['ported'] / 'genesis.json')['payload']
g5 = load(V5_OUT / 'final/ported/genesis.json')['payload']
genesis_diff = []


def diff(a, b, path=''):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if a.get(k) != b.get(k):
                diff(a.get(k), b.get(k), path + '/' + k)
    elif a != b:
        genesis_diff.append(path)


diff(g6, g5)
need(sorted(genesis_diff) == ['/genesis/history/anchor_root_sha256', '/metadata/history/anchor_root_sha256', '/source_sha256'],
     'genesis checkpoint differs beyond source identity/anchor root: %s' % genesis_diff)
history = load(d['history'] / 'v5-history-tests.json')
v5_history = load(V5_OUT / 'final/history/v5-history-tests.json')
need(names(history) == names(v5_history) and all(x['passed'] for x in history['checks']), 'v5 history controls differ')
san_history = load(d['san-history'] / 'v5-history-tests.json')
need(names(san_history) == names(load(V5_OUT / 'final/sanitized-history/v5-history-tests.json')),
     'sanitized v5 history controls differ')
flows = load(d['flows'] / 'v6-flows-tests.json')
san_flows = load(d['san-flows'] / 'v6-flows-tests.json')
normal_only = {'long run: resident memory bounded after warm-up (< 8 MiB growth)'}
need(flows['status'] == 'passed' and san_flows['status'] == 'passed' and
     [n for n in names(flows) if n not in normal_only] == names(san_flows), 'flows normal/sanitized check lists differ')
source = load(d['source'] / 'v6-source-arrival-tests.json')
san_source = load(d['san-source'] / 'v6-source-arrival-tests.json')
need(source['status'] == san_source['status'] == 'passed', 'source/arrival sweep failed')
long_run = load(d['flows-long'] / 'v6-flows-long-run.json')
p6, p5 = load(d['parity'] / 'v6.json'), load(d['parity'] / 'v5.json')
step_diff = sorted({k for x, y in zip(p6['steps'], p5['steps']) for k in x if x[k] != y[k]})
need(p6['genesis'] == p5['genesis'] and len(p6['steps']) == len(p5['steps']) == 60 and step_diff == ['anchor_head_sha256'],
     'identity parity trace differs: %s' % step_diff)
bench = {}
for n in [10000, 100000, 1000000]:
    for tag in ['v6', 'v5']:
        j = load(d['bench'] / ('%s-%d.json' % (tag, n)))
        keys = sorted(j['runs'][0])
        med = {k: sorted(r[k] for r in j['runs'])[len(j['runs']) // 2] for k in keys if k != 'VmRSS_KiB'}
        bench['%s-%d' % (tag, n)] = {'median_s': med, 'VmHWM_KiB': j['VmHWM_KiB'], 'repeats': len(j['runs'])}
v5_after = tree_hashes(V5)
need(v5_after == v5_before and tree_hashes(V5_OUT / 'final') == v5_out_before, 'native v5 or its evidence changed')
need(sealed_inputs() == inputs_before, 'sealed inputs changed during verification')
need(json.loads((P / 'SOURCE-SHA256.json').read_text()) == manifest, 'v6 source changed during verification')
artifacts = {str(p.relative_to(P)): sha256(p) for p in sorted(B.glob('*')) if p.is_file()}
(F / 'build-artifact-hashes.json').write_text(json.dumps(artifacts, indent=2) + '\n')
result = {
    'status': 'passed',
    'source_identity': identity,
    'base_v5_identity': V5_IDENTITY,
    'legacy_checks': legacy['check_count'],
    'legacy_identical_to_native_v5_and_sanitized': True,
    'ported_candidate003_controls': {'check_count': ported['check_count'], 'identical_names_to_v5': True,
                                     'allocation_failure_sweep': ported['allocation_failure_sweep']},
    'v5_history_controls': {'normal': history['check_count'], 'sanitized_reduced': san_history['check_count'],
                            'identical_names_to_v5': True},
    'genesis_checkpoint_differences_vs_v5': sorted(genesis_diff),
    'identity_parity_trace': {'steps': 60, 'anchor_interval': 4, 'differing_fields': step_diff,
                              'v6_anchor_root': p6['anchor_root_sha256'], 'v5_anchor_root': p5['anchor_root_sha256']},
    'v6_flows': {'normal': flows['check_count'], 'sanitized_reduced': san_flows['check_count'],
                 'normal_only': sorted(normal_only), 'measurements': flows['measurements']},
    'v6_source_arrival': {k: source[k] for k in ['check_count', 'cases', 'admitted_steps', 'births', 'exact_stream_equalities',
                                                  'source_refused', 'arrival_refused', 'exposure_refused', 'surface_births',
                                                  'interior_births', 'active_layer_births', 'multi_step_committed_cases']},
    'v6_source_arrival_sanitized_reduced': {k: san_source[k] for k in ['check_count', 'exact_stream_equalities', 'births']},
    'long_run_200000': {k: v for k, v in long_run['measurements']['long_run'].items()},
    'admission_benchmark': bench,
    'sealed_inputs': inputs_before,
    'resources': {'orchestrator_cpus': sorted(cpus), 'build': 'both-cores.sh (cores 6,7), two compiles at a time',
                  'runs': 'any-core.sh --wait (one core each, at most two concurrent), nice 15, one thread, no GPU',
                  'rlimit_as_tests': G, 'rlimit_as_benchmarks': 12 * 1024 ** 3, 'launcher': FROZEN},
    'utc': utc()}
(F / 'verification-result.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: result[k] for k in ['status', 'source_identity', 'legacy_checks', 'ported_candidate003_controls',
                                          'v5_history_controls', 'identity_parity_trace']}, indent=1))
