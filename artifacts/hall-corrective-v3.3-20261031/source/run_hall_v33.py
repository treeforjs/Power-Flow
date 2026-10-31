#!/usr/bin/python -I
"""Build and run Hall corrective v3.3 (v3.2 + double-word magnetic update and seed) and its qualification.

  build:  core-lease-20261030/both-cores.sh /usr/bin/python -I run_hall_v33.py build LABEL [TARGET ...]
          or, on the loaned CPU-only cores: flock testing-core<K>.lock ... taskset -c <K,...> nice -n 15 \
              /usr/bin/python -I run_hall_v33.py build LABEL [TARGET ...]
  run:    core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_hall_v33.py run LABEL NAME BINARY [ARGS ...]
          (or one loaned core 0, 1 or 3 with its lease held)
  gpu:    flock <testing-gpu.lock> flock <testing-gpu-UUID.lock> core-lease-20261030/any-core.sh --wait \
              /usr/bin/python -I run_hall_v33.py gpu LABEL NAME BINARY [ARGS ...]   (or inside a SimQueue powerflow job)

ARGS may use the placeholders {TRACE} {NVRTC} {UUID} {OUT} (the run directory). Every compile and run executes
inside arch-base-20261028/run-frozen. Outputs: /home/treefor/projects/.workspace/outputs/hall-corrective-v3.3-20261031/
LABEL/. Sealed inputs (Hall 8dddec0, Hall corrective v3 and v3.1, hall-gpu-resident-v1 common/push, native candidate003)
are hashed before and after every stage and are only read. Loaned cores 0, 1, 3 take CPU-only stages only and no new
stage starts on them once pic-profile/RETURN-CORES-0-3 exists.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import fcntl, hashlib, json, os, platform, subprocess, sys, time

P = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
OUT = Path('/home/treefor/projects/.workspace/outputs/hall-corrective-v3.3-20261031')
N = R / 'pic-reaction-transactions-20261029/native-full-time-candidate003'
H = R / 'hall-coupled-transport-20261029/source'
HD = H / 'development/hall-canonical-driver-20261029'
I = R / 'pic-water-full-time-subflow-20261029/inputs'
M = R / 'mfem-4.10-pic-20261029'
G = R / 'hall-gpu-resident-v1-20261030'
V3 = R / 'hall-corrective-v3-20261030'
V31 = R / 'hall-corrective-v3.1-20261031'
V32 = R / 'hall-corrective-v3.2-20261031'
N7 = R / 'pic-reaction-transactions-20261029/native-full-time-candidate007'   # forms/ read and compiled in place
FORMS = N7 / 'forms'
F = R / 'arch-base-20261028/run-frozen'
LEASE_DIR = Path('/home/treefor/projects/.workspace/outputs/pic-profile')
LEASES = {6: LEASE_DIR / 'testing-core.lock', 7: LEASE_DIR / 'testing-core7.lock', 0: LEASE_DIR / 'testing-core0.lock',
          1: LEASE_DIR / 'testing-core1.lock', 3: LEASE_DIR / 'testing-core3.lock'}   # 0, 1, 3: temporary loan, CPU-only
LOANED = {0, 1}   # core 3 returned to the Spheral session
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
GPU_SEALED = [G / 'common/hall_gpu_common.cpp', G / 'common/hall_gpu_common.hpp', G / 'common/hall_gpu_primitives.cpp',
              G / 'common/hall_gpu_primitives.hpp', G / 'common/device/deterministic_primitives.cuh', G / 'push/hall_gpu_push.cpp',
              G / 'push/hall_gpu_push.hpp', G / 'push/push_core.hpp', G / 'push/push_abi.hpp', G / 'push/push_bodies.hpp',
              G / 'push/push_kernels.cu', G / 'push/push_fold.hpp', G / 'push/hall_adapter.hpp']
ENV = {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
       'PYTHONNOUSERSITE': '1'}


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


def admission(stage):
    cores = os.sched_getaffinity(0)
    if not (cores <= {6, 7} or cores <= LOANED) or not cores or os.getpriority(os.PRIO_PROCESS, 0) < 15:
        fail('only cores 6/7 (or loaned cores 0, 1, 3) at nice >= 15 are admitted; affinity %s' % sorted(cores))
    if stage == 'gpu' and cores & LOANED:
        fail('GPU stages run on cores 6/7 only')
    if cores & LOANED and RETURN.exists():
        fail('loaned cores 0-3 are being returned (RETURN-CORES-0-3 exists)')
    if stage != 'build' and len(cores) != 1:
        fail('run stages need exactly one leased core')
    for c in cores:
        if not held(LEASES[c]):
            fail('core %d lease is not held by the launching wrapper' % c)
    if stage == 'gpu' and not all(held(lk) for lk in GPU_LOCKS):
        fail('both GPU lock files must be held')
    if any(os.getenv(k) != '1' for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']):
        fail('one thread per numerical library required')
    if Path('/proc/sys/kernel/osrelease').read_text().strip() != (R / 'arch-base-20261028/kernel.release').read_text().strip():
        fail('kernel release differs from the frozen runtime')
    return {'stage': stage, 'cpu_affinity': sorted(cores), 'nice': os.getpriority(os.PRIO_PROCESS, 0),
            'core_leases_held': {str(c): held(LEASES[c]) for c in LEASES}, 'gpu_locks_held': [held(lk) for lk in GPU_LOCKS],
            'PF_CORE': os.getenv('PF_CORE'), 'python': sys.executable, 'python_version': platform.python_version(),
            'frozen_launcher_sha256': sha(F), 'time_unix': time.time()}


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
    rec['gpu_sealed_sources'] = {str(p): sha(p) for p in GPU_SEALED}
    rec['hall_v3_sources'] = {str(p): sha(p) for d in ('source', 'tests') for p in sorted((V3 / d).rglob('*')) if p.is_file()}
    rec['hall_v3_1_sources'] = {str(p): sha(p) for d in ('source', 'tests') for p in sorted((V31 / d).rglob('*')) if p.is_file()}
    rec['hall_v3_2_sources'] = {str(p): sha(p) for d in ('source', 'tests') for p in sorted((V32 / d).rglob('*')) if p.is_file()}
    if sha(TRACE) != TRACE_SHA:
        fail('C19 trace differs')
    rec['native_v7_forms_sources'] = {str(p): sha(p) for p in sorted(FORMS.rglob('*')) if p.is_file()}
    for doc in ('API-007.txt', 'HALL-LINK.txt'):
        rec['native_v7_forms_sources'][str(N7 / doc)] = sha(N7 / doc)
    rec['trace_sha256'] = TRACE_SHA
    rec['nvrtc'] = {'library': str(NVRTC), 'sha256': sha(NVRTC)}
    return rec


def own_sources():
    files = sorted(p for d in (P / 'source', P / 'tests') for p in d.rglob('*') if p.is_file() and p.suffix in ('.hpp', '.cpp', '.cu', '.cuh', '.py'))
    return {str(p): sha(p) for p in files + [Path(__file__).resolve()]}


def execute(name, cmd, log_dir, timeout, cwd=None, env=None):
    log = log_dir / (name + '.log')
    start = time.monotonic()
    with log.open('w') as h:
        try:
            rc = subprocess.run([str(F), *map(str, cmd)], cwd=cwd or P, env=env or ENV, stdout=h, stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    r = {'name': name, 'argv': [str(F), *map(str, cmd)], 'returncode': rc, 'elapsed_s': time.monotonic() - start,
         'log': str(log), 'log_sha256': sha(log)}
    print(name, rc, round(r['elapsed_s'], 1), flush=True)
    return r


def isystem():
    out = []
    for d in [N / 'source/include', I, I / 'include', N / 'inputs', H / 'development', FORMS / 'include', M / 'install-mpi/include',
              M / 'tpl/hypre/include', M / 'tpl/metis/include']:
        out += ['-isystem', str(d)]
    return out


COMMON = ['-std=c++20', '-O2', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math', '-ffunction-sections', '-fdata-sections']
GPU_INC = ['-I' + str(G), '-I' + str(G / 'push'), '-isystem', str(CUDA_INC), '-isystem', str(NVRTC_INC), '-isystem', str(FORMS / 'include'),
           '-Wno-c++26-extensions']


def flags(tree, gpu=False, extra=()):
    src = (V3 if tree == 'v3' else P) / 'source'
    f = ['/usr/bin/c++', *COMMON, '-I' + str(src), '-I' + str(P / 'tests'), '-I' + str(G), *isystem()]
    if gpu:
        f += GPU_INC
    return f + list(extra)


def objects_for(tree, gpu):
    src = (V3 if tree == 'v3' else P) / 'source'
    o = {'motion_provider_' + tree: (flags(tree), src / 'motion_provider.cpp'),
         'kinetic_boundary': (['/usr/bin/c++', *COMMON, *isystem()], HD / 'kinetic_boundary.cpp'),
         'forms_cpu_v7': (['/usr/bin/c++', *COMMON, '-isystem', str(FORMS / 'include')], FORMS / 'src/forms_cpu.cpp')}
    if gpu:
        o.update({'gpu_backends_' + tree: (flags(tree, True), src / 'gpu/gpu_backends.cpp'),
                  'hall_gpu_push_' + tree: (['/usr/bin/c++', *COMMON, '-I' + str(src), *GPU_INC], src / 'gpu/hall_gpu_push_v3.cpp'),
                  'hall_gpu_common': (['/usr/bin/c++', *COMMON, *GPU_INC], G / 'common/hall_gpu_common.cpp'),
                  'hall_gpu_primitives': (['/usr/bin/c++', *COMMON, *GPU_INC], G / 'common/hall_gpu_primitives.cpp')})
        if tree == 'v32':
            o['resident_field'] = (flags(tree, True), src / 'gpu/resident_field.cpp')
    return o


# target -> (tree, gpu, main source, extra flags)
TARGETS = {
    'test_hall_v3': ('v32', False, 'tests/test_hall_v3.cpp', []),
    'test_resident': ('v32', True, 'tests/test_resident.cpp', []),
    'test_v32': ('v32', True, 'tests/test_v32.cpp', []),
}


def build(label, targets):
    adm = admission('build')
    out = OUT / label / 'build'
    out.mkdir(parents=True, exist_ok=False)
    b = P / 'build' / label
    b.mkdir(parents=True, exist_ok=False)
    pins = sealed_pins()
    before = own_sources()
    objs, links = {}, []
    for t in targets:
        tree, gpu, main, extra = TARGETS[t]
        o = objects_for(tree, gpu)
        objs.update(o)
        objs['main_' + t] = (flags(tree, gpu, extra), P / main)
        links.append((t, ['main_' + t, *o.keys()], ['-ldl'] if gpu else []))
    receipt = {'schema': 'hall-corrective-v3.3-build-receipt-v1', 'label': label, 'targets': targets, 'admission': adm, 'sealed': pins, 'commands': []}
    status = 'passed'
    try:
        jobs = [('compile-' + k, [*fl, '-MD', '-MF', b / (k + '.d'), '-c', src, '-o', b / (k + '.o')]) for k, (fl, src) in objs.items()]
        workers = 2 if len(adm['cpu_affinity']) == 2 else 1
        with ThreadPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(lambda j: execute(j[0], j[1], out, 3600), jobs))
        receipt['commands'] += results
        for r in results:
            if r['returncode'] != 0:
                fail(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-8000:])
        link = ['/usr/bin/mpic++', '-Wl,--gc-sections']
        tail = [NATIVE_LIB, *LIBS, '-lcrypto']
        for name, parts, extra in links:
            r = execute('link-' + name, [*link, *[b / (p + '.o') for p in parts], *tail, *extra, '-o', b / name], out, 1800)
            receipt['commands'].append(r)
            if r['returncode'] != 0:
                fail(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-8000:])
        deps = set()
        for k in objs:
            text = (b / (k + '.d')).read_text().replace('\\\n', ' ')
            for line in text.splitlines():
                if ':' in line:
                    deps |= {str(Path(x).resolve()) for x in line.split(':', 1)[1].split()}
        receipt['dependency_sha256'] = {d: sha(d) for d in sorted(deps)}
        receipt['binaries_sha256'] = {n: sha(b / n) for n, _, _ in links}
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


def run(stage, label, name, binary, args):
    adm = admission(stage)
    br = load_build(label)
    pins = sealed_pins()
    out = OUT / label / ('%s-%s-%s' % (stage, name, time.strftime('%H%M%S')))
    out.mkdir(parents=True)
    sub = {'{TRACE}': str(TRACE), '{NVRTC}': str(NVRTC), '{UUID}': GPU_UUID, '{OUT}': str(out)}
    def fill(a):
        for k, v in sub.items():
            a = a.replace(k, v)
        return a
    argv = [Path(br['build_dir']) / binary] + [fill(a) for a in args]
    r = execute(name, argv, out, 12 * 3600)
    receipt = {'schema': 'hall-corrective-v3.3-run-receipt-v1', 'stage': stage, 'label': label, 'name': name, 'admission': adm,
               'sealed': pins, 'build_receipt': str(OUT / label / 'build/build-receipt.json'), 'binary_sha256': sha(argv[0]),
               'own_sources_unchanged_since_build': own_sources() == br['own_sources_sha256'], 'command': r}
    if sealed_pins() != pins:
        fail('sealed inputs changed during the run')
    (out / 'run-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'returncode': r['returncode']}))
    return 0 if r['returncode'] == 0 else 1


if __name__ == '__main__':
    a = sys.argv
    if len(a) < 3 or a[1] not in ('build', 'run', 'gpu') or '/' in a[2]:
        raise SystemExit(__doc__)
    if a[1] == 'build':
        raise SystemExit(build(a[2], a[3:] or list(TARGETS)))
    raise SystemExit(run('gpu' if a[1] == 'gpu' else 'cpu', a[2], a[3], a[4], a[5:]))
