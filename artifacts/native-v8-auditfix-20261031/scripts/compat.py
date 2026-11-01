"""Native v8 compatible mode: the v7.1 suites (host mode; device modes with --device) run against the v8
build and compared with the native v7.1 evidence (outputs; checkpoints except the fields bound to the source
identity). usage: python3 scripts/compat.py LABEL [--device] [--device-only] [--queued] [--compare-only] [suite ...]
Every CPU process runs under scripts/lease.sh inside the frozen launcher (one core of 6, 7, 0, 1; nice 15);
GPU processes take the two GPU locks first. --device-only skips the host runs; --queued runs the device processes
directly (inside a SimQueue job, which holds the core and the GPU share). Device outputs are compared with the
native v7.1 device evidence of the same mode. --compare-only runs nothing and compares the outputs already in
compat-LABEL (CPU only)."""
import sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent))
from tooling import *

args = [a for a in sys.argv[1:] if not a.startswith('--')]
label, only = args[0], set(args[1:])
device = '--device' in sys.argv or '--device-only' in sys.argv
device_only = '--device-only' in sys.argv
queued = '--queued' in sys.argv
compare_only = '--compare-only' in sys.argv
F = OUT / ('compat-' + label)
F.mkdir(parents=True, exist_ok=True)
B = P / 'build'
G = 3 * 1024 ** 3
V71F = V71_OUT / 'final'
H2O = P / 'inputs/h2o-gas-20261029.json'
records = []


def execute(name, argv, timeout=7200):
    log = F / (name + '.log')
    start = time.monotonic()
    with log.open('w') as h:
        rc = subprocess.run([str(a) for a in argv], stdout=h, stderr=subprocess.STDOUT, timeout=timeout, cwd=str(P)).returncode
    records.append({'name': name, 'argv': [str(a) for a in argv], 'returncode': rc,
                    'elapsed_s': round(time.monotonic() - start, 3), 'log_sha256': sha256(log), 'utc': utc()})
    print(name, rc, round(time.monotonic() - start, 1), flush=True)
    return rc


def cpu(name, argv):
    return execute(name, [P / 'scripts/lease.sh', 'env', 'CUDA_VISIBLE_DEVICES=', FROZEN, 'prlimit', '--as=%d' % G, '--', *argv])


def gpu(name, argv, mode):
    env = ['env', 'CUDA_VISIBLE_DEVICES=' + GPU_UUID, 'PF_V71_EXECUTION=' + mode, 'PF_V71_UUID=' + GPU_UUID, 'PF_V71_NVRTC=' + str(NVRTC)]
    if queued:
        return execute(name, ['nice', '-n', '15', *env, 'OMP_NUM_THREADS=1', FROZEN, *argv])
    return execute(name, ['flock', LEASES / 'testing-gpu.lock', 'flock', LEASES / ('testing-gpu-%s.lock' % GPU_UUID),
                          P / 'scripts/lease.sh', *env, FROZEN, *argv])


def d(name):
    p = F / name
    p.mkdir(exist_ok=True)
    return p


SUITES = [('ported', 'test_ported', []), ('work', 'test_work_gross', []), ('forms-native', 'test_forms_native', []),
          ('flows', 'test_flows', []), ('source', 'test_source_arrival', []), ('history-reduced', 'test_history', ['reduced'])]
runs = [('legacy', [B / 'test_legacy', H2O, d('legacy')]), ('history', [B / 'test_history', d('history')]),
        ('parity-host', [B / 'test_parity', d('parity') / 'host.json', 60])]
runs += [('host-' + n, [B / e, d('host-' + n), *x]) for n, e, x in SUITES]
runs = [r for r in runs if (not only or r[0] in only) and not device_only and not compare_only]
with ThreadPoolExecutor(max_workers=3) as pool:
    codes = dict(zip([r[0] for r in runs], pool.map(lambda r: cpu(*r), runs)))
if device and not compare_only:
    for mode in ['device_verified', 'device']:
        for n, e, x in SUITES:
            if not only or ('host-' + n) in only:
                codes['%s-%s' % (mode, n)] = gpu('%s-%s' % (mode, n), [B / e, d('%s-%s' % (mode, n)), *x], mode)
load = lambda p: json.loads(Path(p).read_text())


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
    leaf = k.rsplit('/', 1)[-1]
    return k == '/seconds' or ('/measurements/' in k and ('KiB' in leaf or leaf.endswith('_s') or leaf == 'seconds'))


def differences(a, b, ignore=(), exact=(), contains=()):
    fa, fb = dict(flat(a)), dict(flat(b))
    return sorted(k for k in set(fa) | set(fb) if fa.get(k) != fb.get(k) and not measured(k)
                  and not (ignore and k.endswith(tuple(ignore))) and k not in exact and not any(c in k for c in contains))


IDENTITY_BOUND = {'exact': ('/sha256', '/payload/source_sha256'), 'contains': ('anchor_root_sha256', '/history/anchors/')}
OUTPUT = {'ported': 'candidate005-ported-tests.json', 'work': 'v7-work-gross.json', 'forms-native': 'v7-forms-native.json',
          'flows': 'v6-flows-tests.json', 'source': 'v6-source-arrival-tests.json', 'history-reduced': 'v5-history-tests.json'}
diff = {}
for name in OUTPUT:
    if (F / ('host-' + name) / OUTPUT[name]).exists():
        diff[name] = differences(load(V71F / ('host-' + name) / OUTPUT[name]), load(F / ('host-' + name) / OUTPUT[name]))
        for f in sorted((V71F / ('host-' + name)).glob('*.json')):
            if f.name not in (OUTPUT[name], 'execution.json'):
                diff['host-%s/%s' % (name, f.name)] = differences(load(f), load(F / ('host-' + name) / f.name), **IDENTITY_BOUND)
for mode in ['device_verified', 'device']:
    for name in OUTPUT:
        here, there = F / ('%s-%s' % (mode, name)), V71F / ('%s-%s' % (mode, name))
        if (here / OUTPUT[name]).exists() and (there / OUTPUT[name]).exists():
            diff['%s-%s' % (mode, name)] = differences(load(there / OUTPUT[name]), load(here / OUTPUT[name]))
            for f in sorted(there.glob('*.json')):
                if f.name not in (OUTPUT[name], 'execution.json'):
                    diff['%s-%s/%s' % (mode, name, f.name)] = differences(load(f), load(here / f.name), **IDENTITY_BOUND) \
                        if (here / f.name).exists() else ['missing']
for name in ['history', 'legacy']:
    if (F / name).exists() and any((F / name).iterdir()):
        for f in sorted((V71F / name).glob('*.json')):
            if not (F / name / f.name).exists():
                diff['%s/%s' % (name, f.name)] = ['missing']
            else:
                diff['%s/%s' % (name, f.name)] = differences(load(f), load(F / name / f.name), **IDENTITY_BOUND) \
                    if f.name != 'arena-tests.json' else ([] if load(f) == load(F / name / f.name) else ['differs'])
if (F / 'parity/host.json').exists():
    diff['parity'] = differences(load(V71F / 'parity/host.json'), load(F / 'parity/host.json'),
                                 ignore=('/seconds', 'anchor_head_sha256', 'anchor_root_sha256', 'source_identity'))
summary = {'codes': codes, 'differences': {k: v for k, v in diff.items() if v}, 'compared': sorted(diff),
           'records': records, 'utc': utc()}
(F / 'compat-result.json').write_text(json.dumps(summary, indent=1) + '\n')
print(json.dumps({'codes': codes, 'differing': {k: v[:6] for k, v in diff.items() if v}}, indent=1))
