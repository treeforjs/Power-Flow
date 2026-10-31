"""Native v7.1 final verification with receipts. Start restricted to cores 6-7:
  taskset -c 6,7 nice -n 15 env CUDA_VISIBLE_DEVICES= OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \\
      MKL_NUM_THREADS=1 python3 verify-final.py LABEL
Builds run under core-lease-20261030/both-cores.sh (build-components.py all); every CPU
test process runs under scripts/lease.sh (one free core of 6, 7 and the lent cores 0, 1;
nice 15, one thread, no GPU) inside the frozen launcher, at most three at a time (the GPU sequence runs beside them); GPU processes take the
two GPU locks first, then a core lease. Host mode: every suite, compared with the native v7
evidence; ASan/UBSan (host mode). Device modes (device_verified, device): the same suites
(history reduced) compared with host mode. Timing: native per-step host cost, native v7
(reference build of the same driver) vs v7.1 host vs v7.1 device, 1e4/1e5/1e6 carriers,
200 um and 50 um meshes. Never modifies native v7, its evidence or any sealed input."""
import sys
import threading
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
records_lock = threading.Lock()
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
    with records_lock:
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
    return execute(name, [P / 'scripts/lease.sh', *env, *prefix, *argv], timeout)


def gpu_run(name, argv, mode=None, timeout=3600):
    env = ['env', 'CUDA_VISIBLE_DEVICES=' + GPU_UUID]
    if mode:
        env += ['PF_V71_EXECUTION=' + mode, 'PF_V71_UUID=' + GPU_UUID, 'PF_V71_NVRTC=' + str(NVRTC)]
    return execute(name, ['flock', LEASES / 'testing-gpu.lock', 'flock', LEASES / ('testing-gpu-%s.lock' % GPU_UUID),
                          P / 'scripts/lease.sh', *env, FROZEN, *argv], timeout)


v7_before = tree_hashes(V7)
v7_out_before = tree_hashes(V7_OUT / 'final')
rc = execute('build', [LEASE / 'both-cores.sh', 'env', 'CUDA_VISIBLE_DEVICES=', sys.executable, P / 'build-components.py', 'all'],
             timeout=7200)
need(rc == 0, 'build failed')
manifest = json.loads((P / 'SOURCE-SHA256.json').read_text())
identity = manifest['source_identity']
B = P / 'build'
G = 3 * 1024 ** 3
H2O = P / 'inputs/h2o-gas-20261029.json'
SUITES = [('ported', 'test_ported', []), ('work', 'test_work_gross', []), ('forms-native', 'test_forms_native', []),
          ('flows', 'test_flows', []), ('source', 'test_source_arrival', []), ('history-reduced', 'test_history', ['reduced'])]
d = {}


def out_dir(name):
    d[name] = F / name
    d[name].mkdir(exist_ok=True)
    return d[name]


runs = [('legacy', [B / 'test_legacy', H2O, out_dir('legacy')], False, G),
        ('history', [B / 'test_history', out_dir('history')], False, G),
        ('forms-cpu', [B / 'test_forms_cpu', out_dir('forms-cpu'), 1000000], False, G),
        ('flows-long-200000', [B / 'test_flows', out_dir('flows-long'), 'long=200000'], False, G),
        ('parity-host', [B / 'test_parity', out_dir('parity') / 'host.json', 60], False, G),
        ('parity-v7', [B / 'v7_test_parity', d['parity'] / 'v7.json', 60], False, G),
        ('san-legacy', [B / 'sanitized-test_legacy', H2O, out_dir('san-legacy')], True, None),
        ('san-ported', [B / 'sanitized-test_ported', out_dir('san-ported')], True, None),
        ('san-history', [B / 'sanitized-test_history', out_dir('san-history'), 'reduced'], True, None),
        ('san-flows', [B / 'sanitized-test_flows', out_dir('san-flows'), 'reduced'], True, None),
        ('san-source-arrival', [B / 'sanitized-test_source_arrival', out_dir('san-source'), 'reduced'], True, None),
        ('san-forms-cpu', [B / 'sanitized-test_forms_cpu', out_dir('san-forms-cpu'), 100000], True, None),
        ('san-forms-native', [B / 'sanitized-test_forms_native', out_dir('san-forms-native'), 'reduced'], True, None),
        ('san-work-gross', [B / 'sanitized-test_work_gross', out_dir('san-work-gross')], True, None)]
for name, exe, extra in SUITES:
    runs.append(('host-' + name, [B / exe, out_dir('host-' + name), *extra], False, G))
for mode in ['device_verified', 'device']:
    for name, exe, extra in SUITES:
        out_dir('%s-%s' % (mode, name))


def gpu_phase():
    """GPU forms parity (two fresh processes) and the device-mode suites; runs beside the CPU tests."""
    codes = {}
    for k in [1, 2]:
        codes['gpu-forms-run%d' % k] = gpu_run('gpu-forms-run%d' % k, [B / 'test_forms_gpu', NVRTC, GPU_UUID, 1000000])
    for mode in ['device_verified', 'device']:
        for name, exe, extra in SUITES:
            tag = '%s-%s' % (mode, name)
            codes[tag] = gpu_run(tag, [B / exe, d[tag], *extra], mode)
        codes[mode + '-parity'] = gpu_run(mode + '-parity', [B / 'test_parity', d['parity'] / ('%s.json' % mode), 60], mode)
    return codes


with ThreadPoolExecutor(max_workers=4) as pool:
    gpu_future = pool.submit(gpu_phase)
    futures = {name: pool.submit(cpu_run, name, argv, san, lim or G) for name, argv, san, lim in runs}
    codes = {name: f.result() for name, f in futures.items()}
    codes.update(gpu_future.result())
need(all(c == 0 for c in codes.values()), 'test failures: %s' % {k: v for k, v in codes.items() if v != 0})
rc = cpu_run('verify-exact', ['/usr/bin/python', '-I', P / 'forms/tests/verify_exact.py', d['forms-cpu'] / 'adversarial.txt'])
need(rc == 0, 'rational verification failed')
timing = out_dir('timing')
for mesh in [45, 180]:
    for n, steps in [(10000, 8), (100000, 5), (1000000, 4)]:
        tag = '%d-%d' % (mesh, n)
        need(cpu_run('timing-v7-host-' + tag, [B / 'v7_test_step_timing', timing / ('v7-host-%s.json' % tag), 'host', n, steps,
                                               mesh, NVRTC, GPU_UUID], rlimit=12 * 1024 ** 3) == 0, 'timing failed')
        need(cpu_run('timing-v71-host-' + tag, [B / 'test_step_timing', timing / ('v71-host-%s.json' % tag), 'host', n, steps,
                                                mesh, NVRTC, GPU_UUID], rlimit=12 * 1024 ** 3) == 0, 'timing failed')
        need(gpu_run('timing-v71-device-' + tag, [B / 'test_step_timing', timing / ('v71-device-%s.json' % tag), 'device', n,
                                                  steps, mesh, NVRTC, GPU_UUID]) == 0, 'timing failed')
    # Device lane with the v7 streaming field/current digests (hashed on the host): the cost the block form removes.
    need(gpu_run('timing-v71-device_stream-%d-10000' % mesh, [B / 'test_step_timing', timing / ('v71-device_stream-%d-10000.json' % mesh),
                                                            'device_stream', 10000, 8, mesh, NVRTC, GPU_UUID]) == 0, 'timing failed')

# 2. Parity and results.
load = lambda p: json.loads(Path(p).read_text())
V7F = V7_OUT / 'final'


def flat(x, p=''):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from flat(v, p + '/' + k)
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from flat(v, p + '/' + str(i))
    else:
        yield p, x


def measured(k):
    """Run-dependent measurements (memory, elapsed seconds) recorded by the suites."""
    leaf = k.rsplit('/', 1)[-1]
    return k == '/seconds' or ('/measurements/' in k and ('KiB' in leaf or leaf.endswith('_s') or leaf == 'seconds'))


def differences(a, b, ignore=(), exact=(), contains=()):
    """Paths whose values differ, other than measurements and the given paths."""
    fa, fb = dict(flat(a)), dict(flat(b))
    return sorted(k for k in set(fa) | set(fb) if fa.get(k) != fb.get(k) and not measured(k)
                  and not (ignore and k.endswith(tuple(ignore))) and k not in exact and not any(c in k for c in contains))


# Checkpoint fields bound to the source identity: the file digest, source_sha256 and the history anchors.
IDENTITY_BOUND = {'exact': ('/sha256', '/payload/source_sha256'), 'contains': ('anchor_root_sha256', '/history/anchors/')}


OUTPUT = {'ported': 'candidate005-ported-tests.json', 'work': 'v7-work-gross.json', 'forms-native': 'v7-forms-native.json',
          'flows': 'v6-flows-tests.json', 'source': 'v6-source-arrival-tests.json', 'history-reduced': 'v5-history-tests.json'}
V7_NAME = {'ported': 'ported', 'work': 'work-gross', 'forms-native': 'forms-native', 'flows': 'flows', 'source': 'source'}
host_vs_v7 = {}
for name, v7dir in V7_NAME.items():
    host_vs_v7[name] = differences(load(V7F / v7dir / OUTPUT[name]), load(d['host-' + name] / OUTPUT[name]))
host_vs_v7['history'] = differences(load(V7F / 'history' / OUTPUT['history-reduced']), load(d['history'] / OUTPUT['history-reduced']))
for name, v7dir, ours in [('ported', 'ported', 'host-ported'), ('forms-native', 'forms-native', 'host-forms-native'),
                          ('flows', 'flows', 'host-flows'), ('history', 'history', 'history'), ('legacy', 'legacy', 'legacy')]:
    # Checkpoints: equal except the fields bound to the source identity (source_sha256, anchors, file digest).
    for f in sorted((V7F / v7dir).glob('*.json')):
        if f.name != OUTPUT.get(name, '') and f.name != 'arena-tests.json':
            host_vs_v7['%s/%s' % (name, f.name)] = differences(load(f), load(d[ours] / f.name), **IDENTITY_BOUND)
need(all(not v for v in host_vs_v7.values()), 'host-mode outputs differ from native v7: %s' % host_vs_v7)
need(load(V7F / 'legacy/arena-tests.json') == load(d['legacy'] / 'arena-tests.json') == load(d['san-legacy'] / 'arena-tests.json'),
     'legacy results differ from native v7')
device_vs_host = {}
for mode in ['device_verified', 'device']:
    for name in OUTPUT:
        device_vs_host['%s-%s' % (mode, name)] = differences(load(d['host-' + name] / OUTPUT[name]),
                                                            load(d['%s-%s' % (mode, name)] / OUTPUT[name]))
        # Every other file the suite wrote (checkpoints): identical bytes.
        for f in sorted(d['host-' + name].glob('*.json')):
            if f.name not in (OUTPUT[name], 'execution.json'):
                same = f.read_bytes() == (d['%s-%s' % (mode, name)] / f.name).read_bytes()
                device_vs_host['%s-%s/%s' % (mode, name, f.name)] = [] if same else ['bytes differ']
need(all(not v for v in device_vs_host.values()), 'device-mode outputs differ from host mode: %s' % device_vs_host)
executions = {k: load(d[k] / 'execution.json') for k in d if k.startswith('device')}
for mode in ['device_verified', 'device']:
    executions[mode + '-parity'] = load(d['parity'] / ('%s.json.execution.json' % mode))
p_host, p_v7 = load(d['parity'] / 'host.json'), load(d['parity'] / 'v7.json')
parity = {'v7_vs_host': differences(p_v7, p_host, ignore=('/seconds', 'anchor_head_sha256', 'anchor_root_sha256', 'source_identity'))}
for mode in ['device_verified', 'device']:
    parity[mode + '_vs_host'] = differences(p_host, load(d['parity'] / ('%s.json' % mode)))
need(not parity['v7_vs_host'] and not parity['device_verified_vs_host'] and not parity['device_vs_host'],
     'parity traces differ: %s' % parity)


def results_of(log):
    return [line.strip() for line in log.read_text().splitlines() if line.startswith('RESULT ')]


g1, g2 = results_of(F / 'gpu-forms-run1.log'), results_of(F / 'gpu-forms-run2.log')
strip = lambda lines: [' '.join(t for t in l.split() if 'ms' not in t and 'wall' not in t and 'us_per' not in t) for l in lines]
need(strip(g1) == strip(g2) and any('gpu_checks=' in l and 'failures=0' in l for l in g1), 'GPU forms runs differ or failed')
timings = {p.stem: load(p) for p in sorted(timing.glob('*.json'))}
need(tree_hashes(V7) == v7_before and tree_hashes(V7_OUT / 'final') == v7_out_before, 'native v7 or its evidence changed')
need(sealed_inputs() == inputs_before, 'sealed inputs changed during verification')
need(json.loads((P / 'SOURCE-SHA256.json').read_text()) == manifest, 'v7.1 source changed during verification')
artifacts = {str(p.relative_to(P)): sha256(p) for p in sorted(B.glob('*')) if p.is_file()}
(F / 'build-artifact-hashes.json').write_text(json.dumps(artifacts, indent=2) + '\n')
result = {'status': 'passed', 'source_identity': identity, 'base_v7_identity': V7_IDENTITY,
          'host_vs_v7_differences': host_vs_v7, 'device_vs_host_differences': device_vs_host, 'executions': executions,
          'parity': parity, 'forms_cpu': results_of(F / 'forms-cpu.log'),
          'exact_rational_verification': results_of(F / 'verify-exact.log'), 'gpu_forms': g1, 'timings': timings,
          'sanitized': {k: (F / (k + '.log')).read_text().strip().splitlines()[-1] for k in
                        ['san-legacy', 'san-ported', 'san-history', 'san-flows', 'san-source-arrival', 'san-forms-cpu',
                         'san-forms-native', 'san-work-gross']},
          'long_run_200000': load(d['flows-long'] / 'v6-flows-long-run.json')['measurements']['long_run'],
          'suite_summaries': {rec['name']: Path(rec['log']).read_text().strip().splitlines()[-1] for rec in records
                              if rec['name'] != 'build' and Path(rec['log']).read_text().strip()},
          'sealed_inputs': inputs_before, 'utc': utc()}
(F / 'verification-result.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: result[k] for k in ['status', 'source_identity', 'parity']}, indent=1))
