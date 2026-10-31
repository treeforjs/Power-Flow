#!/usr/bin/python -I
"""Build and run the resolution-series preparation programs against the sealed Hall corrective v3.1 sources.

  build:  core-lease-20261030/any-core.sh --wait /usr/bin/python -I build_run.py build LABEL [--sanitize]
  run:    <lease wrapper> /usr/bin/python -I build_run.py run LABEL NAME BINARY [ARGS ...]

Leases: CPU runs on core 6 or 7 (any-core.sh) or on one lent core 0, 1 or 3 (flock testing-core<K>.lock, taskset, nice 15,
refused when pic-profile/RETURN-CORES-0-3 exists). GPU runs: SimQueue (simq submit -p powerflow -c 1 --gpus 1 --gpu-mem X ...), which
orders cores and the GPU itself, or flock testing-gpu.lock flock testing-gpu-<UUID>.lock any-core.sh --wait ... (GPU locks first, then
the core lease; refused unless both are held).
ARGS may use {TRACE} {NVRTC} {UUID} {OUT}. Every compile and run executes inside arch-base-20261028/run-frozen.
Hall corrective v3.1 (sources, tests), Hall 8dddec0, hall-gpu-resident-v1 common/push, native candidate003 and the C19
trace are hashed before and after every stage and are only read.
Outputs: /home/treefor/projects/.workspace/outputs/resolution-series-prep-20261031/LABEL/.
"""
from pathlib import Path
import fcntl, hashlib, json, os, platform, subprocess, sys, time

P = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
OUT = Path('/home/treefor/projects/.workspace/outputs/resolution-series-prep-20261031')
V31 = R / 'hall-corrective-v3.1-20261031'
N = R / 'pic-reaction-transactions-20261029/native-full-time-candidate003'
H = R / 'hall-coupled-transport-20261029/source'
HD = H / 'development/hall-canonical-driver-20261029'
I = R / 'pic-water-full-time-subflow-20261029/inputs'
M = R / 'mfem-4.10-pic-20261029'
G = R / 'hall-gpu-resident-v1-20261030'
F = R / 'arch-base-20261028/run-frozen'
LEASE_DIR = Path('/home/treefor/projects/.workspace/outputs/pic-profile')
LEASES = {6: LEASE_DIR / 'testing-core.lock', 7: LEASE_DIR / 'testing-core7.lock', 0: LEASE_DIR / 'testing-core0.lock',
          1: LEASE_DIR / 'testing-core1.lock', 3: LEASE_DIR / 'testing-core3.lock'}
LENT = {0, 1, 3}
RETURN = LEASE_DIR / 'RETURN-CORES-0-3'
GPU_UUID = 'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9'
GPU_LOCKS = [LEASE_DIR / 'testing-gpu.lock', LEASE_DIR / ('testing-gpu-%s.lock' % GPU_UUID)]
NVRTC = R / 'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/lib/libnvrtc.so.12'
NVRTC_INC = R / 'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/include'
CUDA_INC = R / 'mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include'
HALL_COMMIT = '8dddec0610982778bf02bdfe697b70a37129c012'
IDENTITY = 'c045c42d2f57082e3fb5e12172c440d2470e681a38df9ac406bd7e063609a08d'
NATIVE_LIB = N / 'build/libpic_native_candidate003.a'
NATIVE_LIB_SHA = 'd2568af9fac13f8aa056ae4a9a50bb2a07eee05d7218d1c4b79aef21317d8a20'
TRACE = H / 'development/hall-finite-foil-em-20261029/data/shot_17133_current.csv'
TRACE_SHA = '2c54d2dbaec6b6c8fcf2e45b5abd837aa72999f5fd14217f530277ffd5949076'
LIBS = [M / 'install-mpi/lib/libmfem.a', M / 'tpl/hypre/lib/libHYPRE.a', M / 'tpl/metis/lib/libmetis.a', M / 'tpl/GKlib/lib/libGKlib.a']
ENV = {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
       'PYTHONNOUSERSITE': '1'}
if os.getenv('CUDA_VISIBLE_DEVICES') is not None:   # SimQueue GPU allocation (the programs select the device by UUID)
    ENV['CUDA_VISIBLE_DEVICES'] = os.environ['CUDA_VISIBLE_DEVICES']
# Hall v3.1 sealed record of its own sources (README.txt: sealed); hashed and compared before/after every stage.
V31_FILES = sorted(p for d in ('source', 'tests') for p in (V31 / d).rglob('*') if p.is_file()) + [V31 / 'README.txt', V31 / 'run_hall_v31.py']


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fail(m):
    raise RuntimeError(m)


def held(lock):
    with open(lock, 'a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(handle, fcntl.LOCK_UN)
        return False


def simq_context():
    """SimQueue job context (SimQueue orders cores, GPU and memory itself)."""
    if not os.getenv('SIMQ_JOB_ID'):
        return None
    cpus = os.getenv('SIMQ_CPU_LIST', '')
    want = set()
    for part in cpus.split(','):
        if '-' in part:
            a, b = part.split('-'); want |= set(range(int(a), int(b) + 1))
        elif part:
            want.add(int(part))
    return {'job': os.getenv('SIMQ_JOB_ID'), 'partition': os.getenv('SIMQ_PARTITION'), 'cpu_list': sorted(want),
            'gpus': os.getenv('SIMQ_GPUS'), 'gpu_mode': os.getenv('SIMQ_GPU_MODE'), 'gpu_mem': os.getenv('SIMQ_GPU_MEM'),
            'borrow': os.getenv('SIMQ_BORROW'), 'mem': os.getenv('SIMQ_MEM')}


def admission(stage):
    cores = os.sched_getaffinity(0)
    if any(os.getenv(k) != '1' for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']) and not os.getenv('SIMQ_JOB_ID'):
        fail('one thread per numerical library required')
    if Path('/proc/sys/kernel/osrelease').read_text().strip() != (R / 'arch-base-20261028/kernel.release').read_text().strip():
        fail('kernel release differs from the frozen runtime')
    sq = simq_context()
    base = {'stage': stage, 'cpu_affinity': sorted(cores), 'nice': os.getpriority(os.PRIO_PROCESS, 0), 'PF_CORE': os.getenv('PF_CORE'),
            'python': sys.executable, 'python_version': platform.python_version(), 'frozen_launcher_sha256': sha(F), 'time_unix': time.time()}
    if sq:
        # SimQueue: affinity must be the allocated cores; a GPU stage needs an allocated GPU; never cores 2, 4, 5 of other
        # partitions unless borrowed through SimQueue
        if set(cores) != set(sq['cpu_list']) or len(cores) != 1 and stage != 'build':
            fail('SimQueue affinity %s differs from the allocation %s' % (sorted(cores), sq['cpu_list']))
        if stage == 'gpu' and not sq['gpus']:
            fail('GPU stage without a SimQueue GPU allocation')
        base['simqueue'] = sq
        return base
    if len(cores) != 1 or os.getpriority(os.PRIO_PROCESS, 0) < 15:
        fail('exactly one leased core at nice >= 15 is admitted; affinity %s' % sorted(cores))
    c = next(iter(cores))
    if c not in LEASES:
        fail('core %d is not a Power-Flow or lent core' % c)
    if c in LENT and stage == 'gpu':
        fail('GPU stages run on cores 6/7 only')
    if stage == 'gpu' and not all(held(lk) for lk in GPU_LOCKS):
        fail('both GPU lock files must be held by the launching wrapper')
    if c in LENT and RETURN.exists():
        fail('lent cores 0-3 have been recalled')
    if not held(LEASES[c]):
        fail('core %d lease is not held by the launching wrapper' % c)
    base['core_lease_held'] = True
    return base


def sealed_pins():
    rec = {}
    manifest = json.loads((N / 'SOURCE-SHA256.json').read_text())
    if manifest['source_identity'] != IDENTITY:
        fail('native source identity differs')
    for name, digest in manifest['files'].items():
        if sha(N / name) != digest:
            fail('native sealed file changed: ' + name)
    if sha(NATIVE_LIB) != NATIVE_LIB_SHA:
        fail('native candidate003 library differs')
    rec['native'] = {'source_identity': IDENTITY, 'library_sha256': NATIVE_LIB_SHA, 'files_verified': len(manifest['files'])}
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=H, text=True).strip()
    dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=H, text=True)
    if head != HALL_COMMIT or dirty:
        fail('Hall repository is not the clean corrective commit')
    rec['hall'] = {'commit': head, 'git_clean': True}
    rec['hall_v31_sources'] = {str(p): sha(p) for p in V31_FILES}
    rec['gpu_sealed_sources'] = {str(p): sha(p) for p in sorted((G / 'common').rglob('*')) + sorted((G / 'push').rglob('*')) if p.is_file()}
    if sha(TRACE) != TRACE_SHA:
        fail('C19 trace differs')
    rec['trace_sha256'] = TRACE_SHA
    rec['nvrtc'] = {'library': str(NVRTC), 'sha256': sha(NVRTC)}
    return rec


def own_sources():
    # compiled sources and this launcher (analysis scripts are not part of a build)
    files = sorted(p for p in P.rglob('*') if p.is_file() and p.suffix in ('.hpp', '.cpp', '.cu', '.cuh', '.h')
                   and 'build' not in p.relative_to(P).parts and 'merge-split' not in p.relative_to(P).parts) + [Path(__file__).resolve()]
    return {str(p): sha(p) for p in files}


def execute(name, cmd, log_dir, timeout, cwd=None):
    log = log_dir / (name + '.log')
    start = time.monotonic()
    with log.open('w') as h:
        try:
            rc = subprocess.run([str(F), *map(str, cmd)], cwd=cwd or P, env=ENV, stdout=h, stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    r = {'name': name, 'argv': [str(F), *map(str, cmd)], 'returncode': rc, 'elapsed_s': time.monotonic() - start,
         'log': str(log), 'log_sha256': sha(log)}
    print(name, rc, round(r['elapsed_s'], 1), flush=True)
    return r


def isystem():
    out = []
    for d in [N / 'source/include', I, I / 'include', N / 'inputs', H / 'development', M / 'install-mpi/include',
              M / 'tpl/hypre/include', M / 'tpl/metis/include']:
        out += ['-isystem', str(d)]
    return out


COMMON = ['-std=c++20', '-O2', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math', '-ffunction-sections', '-fdata-sections']
GPU_INC = ['-I' + str(G), '-I' + str(G / 'push'), '-isystem', str(CUDA_INC), '-isystem', str(NVRTC_INC), '-Wno-c++26-extensions']
SRC = V31 / 'source'


def flags(gpu=False):
    f = ['/usr/bin/c++', *COMMON, '-I' + str(SRC), '-I' + str(V31 / 'tests'), '-I' + str(G), *isystem()]
    return f + (GPU_INC if gpu else [])


# The Hall v3.1 objects are compiled from the sealed sources with the flags of run_hall_v31.py (target test_resident).
OBJECTS = {
    'motion_provider_v31': (flags(), SRC / 'motion_provider.cpp'),
    'kinetic_boundary': (['/usr/bin/c++', *COMMON, *isystem()], HD / 'kinetic_boundary.cpp'),
    'gpu_backends_v31': (flags(True), SRC / 'gpu/gpu_backends.cpp'),
    'hall_gpu_push_v31': (['/usr/bin/c++', *COMMON, '-I' + str(SRC), *GPU_INC], SRC / 'gpu/hall_gpu_push_v3.cpp'),
    'hall_gpu_common': (['/usr/bin/c++', *COMMON, *GPU_INC], G / 'common/hall_gpu_common.cpp'),
    'hall_gpu_primitives': (['/usr/bin/c++', *COMMON, *GPU_INC], G / 'common/hall_gpu_primitives.cpp'),
    'resident_field': (flags(True), SRC / 'gpu/resident_field.cpp'),
    'main_field_scaling': (flags(True), P / 'hall-scaling/field_scaling.cpp'),
    'main_test_gmg': (flags(True) + ['-I' + str(P / 'precond')], P / 'precond/test_gmg.cpp'),
}
HALL_OBJECTS = [k for k in OBJECTS if not k.startswith('main_')]
TARGETS = {'field_scaling': HALL_OBJECTS + ['main_field_scaling'], 'gmg_solve': HALL_OBJECTS + ['main_test_gmg']}


SAN = ['-fsanitize=address,undefined', '-fno-sanitize-recover=undefined', '-fno-omit-frame-pointer', '-g']


def build(label, sanitize=False):
    adm = admission('build')
    out = OUT / label / 'build'
    out.mkdir(parents=True, exist_ok=False)
    b = P / 'build' / label
    b.mkdir(parents=True, exist_ok=False)
    pins = sealed_pins()
    before = own_sources()
    receipt = {'schema': 'resolution-series-build-receipt-v1', 'label': label, 'admission': adm, 'sealed': pins, 'commands': [],
               'sanitizers': SAN if sanitize else []}
    status = 'passed'
    try:
        for k, (fl, src) in OBJECTS.items():
            r = execute('compile-' + k, [*fl, *(SAN if sanitize else []), '-MD', '-MF', b / (k + '.d'), '-c', src, '-o', b / (k + '.o')], out, 3600)
            receipt['commands'].append(r)
            if r['returncode'] != 0:
                fail(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-8000:])
        for name, parts in TARGETS.items():
            r = execute('link-' + name, ['/usr/bin/mpic++', *(SAN if sanitize else []), '-Wl,--gc-sections', *[b / (p + '.o') for p in parts], NATIVE_LIB, *LIBS, '-lcrypto', '-ldl',
                                         '-o', b / name], out, 1800)
            receipt['commands'].append(r)
            if r['returncode'] != 0:
                fail(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-8000:])
        deps = set()
        for k in OBJECTS:
            text = (b / (k + '.d')).read_text().replace('\\\n', ' ')
            for line in text.splitlines():
                if ':' in line:
                    deps |= {str(Path(x).resolve()) for x in line.split(':', 1)[1].split()}
        receipt['dependency_sha256'] = {d: sha(d) for d in sorted(deps)}
        receipt['binaries_sha256'] = {n: sha(b / n) for n in TARGETS}
        receipt['build_dir'] = str(b)
        receipt['own_sources_sha256'] = before
        if own_sources() != before or sealed_pins() != pins:
            fail('sources or sealed inputs changed during the build')
    except Exception as e:
        status = 'failed'
        receipt['error'] = str(e)
    receipt['status'] = status
    (out / 'build-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    print(json.dumps({'build': str(out), 'status': status, 'error': receipt.get('error', '')[:6000]}))
    return 0 if status == 'passed' else 1


def load_build(label):
    br = json.loads((OUT / label / 'build/build-receipt.json').read_text())
    if br['status'] != 'passed':
        fail('build did not pass')
    for n, d in br['binaries_sha256'].items():
        if sha(Path(br['build_dir']) / n) != d:
            fail('binary differs from its build receipt: ' + n)
    return br


def run(label, name, binary, args):
    gpu = (binary == 'gmg_solve' and '{NVRTC}' in args) or (len(args) > 0 and args[0].startswith('gpu-'))
    adm = admission('gpu' if gpu else 'cpu')
    br = load_build(label)
    pins = sealed_pins()
    out = OUT / label / ('%s-%s' % (name, time.strftime('%H%M%S')))
    out.mkdir(parents=True)
    sub = {'{TRACE}': str(TRACE), '{NVRTC}': str(NVRTC), '{UUID}': GPU_UUID, '{OUT}': str(out)}

    def fill(a):
        for k, v in sub.items():
            a = a.replace(k, v)
        return a
    argv = [Path(br['build_dir']) / binary] + [fill(a) for a in args]
    r = execute(name, argv, out, 12 * 3600)
    receipt = {'schema': 'resolution-series-run-receipt-v1', 'label': label, 'name': name, 'admission': adm, 'sealed': pins,
               'build_receipt': str(OUT / label / 'build/build-receipt.json'), 'binary_sha256': sha(argv[0]),
               'own_sources_unchanged_since_build': own_sources() == br['own_sources_sha256'], 'command': r,
               'gpu_locks_held_by_wrapper': [held(lk) for lk in GPU_LOCKS] if gpu and not os.getenv('SIMQ_JOB_ID') else None}
    if sealed_pins() != pins:
        fail('sealed inputs changed during the run')
    (out / 'run-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'returncode': r['returncode']}))
    return 0 if r['returncode'] == 0 else 1


if __name__ == '__main__':
    a = sys.argv
    if len(a) < 3 or a[1] not in ('build', 'run') or '/' in a[2]:
        raise SystemExit(__doc__)
    if a[1] == 'build':
        raise SystemExit(build(a[2], '--sanitize' in a[3:]))
    raise SystemExit(run(a[2], a[3], a[4], a[5:]))
