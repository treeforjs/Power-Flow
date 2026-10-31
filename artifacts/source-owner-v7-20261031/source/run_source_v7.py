#!/usr/bin/python -I
"""Build and run the heat/source owner v7 qualification (receipted).

  build: core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_source_v7.py build LABEL [normal|ubsan|asan]
  run:   core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_source_v7.py run LABEL [SANITIZER_LABEL ...]
  gpu:   flock <testing-gpu.lock> flock <testing-gpu-UUID.lock> core-lease-20261030/any-core.sh --wait
             /usr/bin/python -I run_source_v7.py gpu LABEL
Cores 0, 1 and 3 (lent for CPU-only work) are admitted too when the caller holds testing-core<K>.lock, pins the
process with taskset at nice 15 and RETURN-CORES-0-3 does not exist; never cores 2, 4, 5.

Binaries per build:
  test_source_v7          v6.1 -> v7 port (A02 probes, digest identity at benchmark settings), trial interval bound,
                          owner-level surface physics (classes, phase map, births, wall responses, ESD release, ledgers)
  test_source_v6_on_v7    the sealed source v6 qualification (pic-electrode-heat-source-v6-20261030/tests/test_source_v6.cpp,
                          unchanged) compiled against owner v7 through compat-v6/ (every v7 option at its v6.1 value)
  test_surface_v7         surface physics without the owner: deterministic math, kernel oracles (Murphy-Good/Forbes FN,
                          Richardson-Schottky, MG temperature factor), Furman-Pivi yields and spectra against the cited data,
                          1e6-event exact ledger closure, RNG determinism over 1/2/5/8 partitions, parameter registry
  test_wall_gpu_v7        (normal builds only) NVRTC wall kernel and fixed-order reduction vs the CPU reference, bitwise
Every compile and run step executes inside arch-base-20261028/run-frozen. Outputs:
/home/treefor/projects/.workspace/outputs/source-owner-v7-20261031/LABEL/. Sealed inputs are hashed and only read.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import fcntl, hashlib, json, os, platform, subprocess, sys, time

P = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
OUT = Path('/home/treefor/projects/.workspace/outputs/source-owner-v7-20261031')
S61 = R / 'pic-electrode-heat-source-v6.1-20261031'
S6 = R / 'pic-electrode-heat-source-v6-20261030'
N = R / 'pic-reaction-transactions-20261029/native-full-time-candidate003'
H = R / 'hall-coupled-transport-20261029/source'
HD = H / 'development/hall-canonical-driver-20261029'
HB = HD / 'bindings/heat-v2'
I = R / 'pic-water-full-time-subflow-20261029/inputs'
M = R / 'mfem-4.10-pic-20261029'
G = R / 'hall-gpu-resident-v1-20261030'
HV3 = R / 'hall-corrective-v3-20261030'
S5 = R / 'pic-electrode-heat-source-v5-20261030/native/mfem_electrode_heat_source_v5'
V3K = R / 'pic-electrode-heat-feedback-20261029/native/mfem_electrode_heat_source_v3'
A1 = R / 'pic-arrival-owner-v1-20261030'
A11 = R / 'pic-arrival-owner-v1.1-20261030'
F = R / 'arch-base-20261028/run-frozen'
V7N = P / 'native/mfem_electrode_heat_source_v7'
V61N = S61 / 'native/mfem_electrode_heat_source_v6_1'
LEASE_DIR = Path('/home/treefor/projects/.workspace/outputs/pic-profile')
LEASES = {0: LEASE_DIR / 'testing-core0.lock', 1: LEASE_DIR / 'testing-core1.lock', 3: LEASE_DIR / 'testing-core3.lock',
          6: LEASE_DIR / 'testing-core.lock', 7: LEASE_DIR / 'testing-core7.lock'}
RETURN_FLAG = LEASE_DIR / 'RETURN-CORES-0-3'
GPU_LOCKS = [LEASE_DIR / 'testing-gpu.lock', LEASE_DIR / 'testing-gpu-GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9.lock']
GPU_UUID = 'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9'
NVRTC_DIR = R / 'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc'
NVRTC = NVRTC_DIR / 'lib/libnvrtc.so.12'
CUDA_INC = R / 'mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include'
NATIVE_LIB = N / 'build/libpic_native_candidate003.a'
LIBS = [M / 'install-mpi/lib/libmfem.a', M / 'tpl/hypre/lib/libHYPRE.a', M / 'tpl/metis/lib/libmetis.a', M / 'tpl/GKlib/lib/libGKlib.a']
TRACE = H / 'development/hall-finite-foil-em-20261029/data/shot_17133_current.csv'
ENV = {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'PYTHONNOUSERSITE': '1'}
VARIANTS = {
    'normal': ['-std=c++20', '-O2', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math'],
    'ubsan': ['-std=c++20', '-O1', '-g', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math', '-fsanitize=undefined', '-fno-sanitize-recover=all', '-fno-omit-frame-pointer'],
    'asan': ['-std=c++20', '-O1', '-g', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math', '-fsanitize=address', '-fno-omit-frame-pointer'],
}
LINK_SAN = {'normal': [], 'ubsan': ['-fsanitize=undefined'], 'asan': ['-fsanitize=address']}


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


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
    if not cores <= set(LEASES) or os.getpriority(os.PRIO_PROCESS, 0) < 15:
        raise RuntimeError('only cores 0, 1, 3, 6, 7 at nice >= 15 are admitted')
    if cores & {0, 1, 3} and RETURN_FLAG.exists():
        raise RuntimeError('RETURN-CORES-0-3 exists: no new work on cores 0-3')
    for c in cores:
        if not held(LEASES[c]):
            raise RuntimeError('core %d lease is not held' % c)
    rec = {'stage': stage, 'cpu_affinity': sorted(cores), 'nice': os.getpriority(os.PRIO_PROCESS, 0), 'python': platform.python_version(),
           'frozen_launcher_sha256': sha(F), 'time_unix': time.time()}
    if stage == 'gpu':
        for lock in GPU_LOCKS:
            if not held(lock):
                raise RuntimeError('GPU lease %s is not held' % lock.name)
        rec['gpu_locks'] = [str(x) for x in GPU_LOCKS]
    return rec


def sealed():
    dirs = (S5, V3K, HB, A1 / 'include', A11 / 'include', S6 / 'native', S6 / 'binding', S6 / 'tests', S61 / 'native', S61 / 'binding', G / 'common')
    files = [p for d in dirs for p in sorted(d.rglob('*')) if p.is_file() and p.suffix in ('.hpp', '.h', '.cpp', '.cuh', '.json', '.txt')]
    return {str(p): sha(p) for p in files + [NATIVE_LIB]}


def own():
    files = sorted(p for p in P.rglob('*') if p.is_file() and p.suffix in ('.hpp', '.h', '.cpp', '.cu', '.py', '.txt') and 'build' not in p.parts)
    return {str(p): sha(p) for p in files}


def execute(name, cmd, log_dir, timeout, env=None):
    log = log_dir / (name + '.log')
    t = time.monotonic()
    with log.open('w') as h:
        try:
            rc = subprocess.run([str(F), *map(str, cmd)], cwd=P, env=env or ENV, stdout=h, stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    print(name, rc, round(time.monotonic() - t, 1), flush=True)
    return {'name': name, 'argv': [str(F), *map(str, cmd)], 'returncode': rc, 'elapsed_s': time.monotonic() - t, 'log': str(log), 'log_sha256': sha(log)}


def owner_flags(base, compat=False):
    inc = (['-I' + str(P / 'compat-v6')] if compat else []) + ['-I' + str(R), '-I' + str(P / 'binding'), '-I' + str(V7N), '-I' + str(V61N),
           '-I' + str(HB / 'source-kernel/include'), '-I' + str(HV3 / 'source'), '-I' + str(G), '-I' + str(HB), '-I' + str(V3K),
           '-I' + str(A11 / 'include'), '-I' + str(A1 / 'include'),
           '-DHEAT_SOURCE_EXACT_COUNT_HEADER="raikou/exact_particle_count.h"', '-DPF_ARRIVAL_EXACT_COUNT_HEADER="raikou/exact_particle_count.h"']
    for d in [N / 'source/include', I, I / 'include', N / 'inputs', H / 'development', M / 'install-mpi/include', M / 'tpl/hypre/include', M / 'tpl/metis/include']:
        inc += ['-isystem', str(d)]
    return ['/usr/bin/c++', *base, *inc]


def light_flags(base):
    return ['/usr/bin/c++', *base, '-I' + str(V7N), '-I' + str(V61N), '-I' + str(HB / 'source-kernel/include'), '-I' + str(A1 / 'include'),
            '-isystem', str(N / 'source/include'), '-isystem', str(N / 'inputs')]


def gpu_flags(base):
    return ['/usr/bin/c++', *base, '-Wno-c++26-extensions', '-I' + str(G), '-I' + str(V7N), '-I' + str(A1 / 'include'), '-isystem', str(N / 'source/include'),
            '-isystem', str(N / 'inputs'),
            '-isystem', str(CUDA_INC), '-isystem', str(NVRTC_DIR / 'include')]


EXES = {
    'test_source_v7': ('owner', ['test_source_v7']),
    'test_source_v6_on_v7': ('owner', ['test_source_v6_on_v7']),
    'test_surface_v7': ('light', ['test_surface_v7']),
    'test_wall_gpu_v7': ('gpu', ['test_wall_gpu_v7', 'hall_gpu_common']),
}


def build(label, variant):
    adm = admission('build')
    if variant not in VARIANTS:
        raise SystemExit('variant must be one of ' + ', '.join(VARIANTS))
    out = OUT / label / 'build'
    out.mkdir(parents=True, exist_ok=False)
    b = P / 'build' / label
    b.mkdir(parents=True, exist_ok=False)
    pins, before = sealed(), own()
    base = VARIANTS[variant]
    se_flags = ['/usr/bin/c++', *base, '-I' + str(HB / 'source-kernel/include')]
    objs = {'test_source_v7': (owner_flags(base), P / 'tests/test_source_v7.cpp'),
            'test_source_v6_on_v7': (owner_flags(base, True), S6 / 'tests/test_source_v6.cpp'),
            'test_surface_v7': (light_flags(base), P / 'tests/test_surface_v7.cpp'),
            'fast_emission': ([*se_flags, '-I' + str(V61N)], V61N / 'fast_emission.cpp'),
            'emission_hot_v7': (se_flags, V7N / 'emission_hot_v7.cpp'),
            'motion_provider_v3': (owner_flags(base), HV3 / 'source/motion_provider.cpp'),
            'kinetic_boundary': (owner_flags(base), HD / 'kinetic_boundary.cpp'),
            'surface_emission': (se_flags, HB / 'source-kernel/src/surface_emission.cpp'),
            'micro_layer': (se_flags, HB / 'source-kernel/src/micro_layer.cpp')}
    exes = ['test_source_v7', 'test_source_v6_on_v7', 'test_surface_v7']
    if variant == 'normal':
        objs['test_wall_gpu_v7'] = (gpu_flags(base), P / 'tests/test_wall_gpu_v7.cpp')
        objs['hall_gpu_common'] = (gpu_flags(base), G / 'common/hall_gpu_common.cpp')
        exes.append('test_wall_gpu_v7')
    receipt = {'schema': 'source-v7-build-receipt-v1', 'label': label, 'variant': variant, 'admission': adm, 'sealed_sha256': pins, 'commands': []}
    status = 'passed'
    try:
        jobs = [('compile-' + k, [*f, '-MD', '-MF', b / (k + '.d'), '-c', src, '-o', b / (k + '.o')]) for k, (f, src) in objs.items()]
        with ThreadPoolExecutor(max_workers=max(1, len(adm['cpu_affinity']))) as pool:
            results = list(pool.map(lambda j: execute(j[0], j[1], out, 3600), jobs))
        receipt['commands'] += results
        for r in results:
            if r['returncode'] != 0:
                raise RuntimeError(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-8000:])
        sealed_objs = [b / (k + '.o') for k in ('fast_emission', 'emission_hot_v7', 'motion_provider_v3', 'kinetic_boundary', 'surface_emission', 'micro_layer')]
        links = {'test_source_v7': ['/usr/bin/mpic++', '-Wl,--gc-sections', *LINK_SAN[variant], b / 'test_source_v7.o', *sealed_objs, NATIVE_LIB, *LIBS, '-lcrypto'],
                 'test_source_v6_on_v7': ['/usr/bin/mpic++', '-Wl,--gc-sections', *LINK_SAN[variant], b / 'test_source_v6_on_v7.o', *sealed_objs, NATIVE_LIB, *LIBS, '-lcrypto'],
                 'test_surface_v7': ['/usr/bin/c++', *LINK_SAN[variant], b / 'test_surface_v7.o', b / 'fast_emission.o', b / 'emission_hot_v7.o', b / 'surface_emission.o', '-lcrypto'],
                 'test_wall_gpu_v7': ['/usr/bin/c++', b / 'test_wall_gpu_v7.o', b / 'hall_gpu_common.o', '-ldl', '-lcrypto']}
        for exe in exes:
            r = execute('link-' + exe, [*links[exe], '-o', b / exe], out, 1800)
            receipt['commands'].append(r)
            if r['returncode'] != 0:
                raise RuntimeError('link failed:\n' + Path(r['log']).read_text()[-8000:])
        deps = set()
        for k in objs:
            for line in (b / (k + '.d')).read_text().replace('\\\n', ' ').splitlines():
                if ':' in line:
                    deps |= {str(Path(x).resolve()) for x in line.split(':', 1)[1].split()}
        receipt['dependency_sha256'] = {d: sha(d) for d in sorted(deps)}
        receipt['sealed_v6_owner_headers_in_binaries'] = sorted(d for d in deps if '/pic-electrode-heat-source-v6-20261030/native/' in d)
        receipt['binary_sha256'] = {exe: sha(b / exe) for exe in exes}
        receipt['build_dir'] = str(b)
        receipt['own_sha256'] = before
        if sealed() != pins or own() != before:
            raise RuntimeError('sources changed during the build')
    except Exception as e:
        status = 'failed'
        receipt['error'] = str(e)
    receipt['status'] = status
    (out / 'build-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    print(json.dumps({'build': str(out), 'status': status, 'error': receipt.get('error', '')[:3000]}))
    return 0 if status == 'passed' else 1


def usable(lab):
    br = json.loads((OUT / lab / 'build/build-receipt.json').read_text())
    for exe, digest in br['binary_sha256'].items():
        if br['status'] != 'passed' or sha(Path(br['build_dir']) / exe) != digest:
            raise RuntimeError('build %s not usable' % lab)
    return br


def summarize(f):
    j = json.loads(f.read_text())
    return {'failed': j['failed'], 'passed': sum(1 for c in j['cases'] if c['passed']), 'cases': len(j['cases']), 'error': j.get('error'),
            'names': [{'group': c['group'], 'name': c['name'], 'passed': c['passed']} for c in j['cases']]}


def run(label, extra):
    adm = admission('run')
    out = OUT / label / ('run-' + time.strftime('%H%M%S'))
    out.mkdir(parents=True)
    summary = {'schema': 'source-v7-run-receipt-v1', 'label': label, 'admission': adm, 'sealed_sha256': sealed(), 'own_sha256': own(), 'commands': []}
    plan = [(label, 'normal')] + [(x, 'sanitizer-' + json.loads((OUT / x / 'build/build-receipt.json').read_text())['variant']) for x in extra]
    for lab, sub in plan:
        br = usable(lab)
        env = dict(ENV)
        if 'asan' in sub:
            env['ASAN_OPTIONS'] = 'detect_leaks=0:abort_on_error=1'   # MPI/MFEM static state is not freed at exit
        for exe, args, result in (('test_surface_v7', [], 'surface-v7.json'), ('test_source_v7', [], 'source-v7.json'),
                                  ('test_source_v6_on_v7', [TRACE], 'source-v6.json')):
            d = out / (sub + '-' + exe)
            d.mkdir()
            r = execute(sub + '-' + exe, [Path(br['build_dir']) / exe, *args, d], out, 7200, env)
            summary['commands'].append(r)
            if (d / result).exists():
                summary[sub + '/' + exe] = summarize(d / result)
    summary['all_passed'] = all(c['returncode'] == 0 for c in summary['commands'])
    (out / 'run-receipt.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'all_passed': summary['all_passed']}))
    return 0 if summary['all_passed'] else 1


def gpu(label):
    adm = admission('gpu')
    br = usable(label)
    out = OUT / label / ('gpu-' + time.strftime('%H%M%S'))
    out.mkdir(parents=True)
    env = dict(ENV, CUDA_VISIBLE_DEVICES=GPU_UUID)
    summary = {'schema': 'source-v7-gpu-receipt-v1', 'label': label, 'admission': adm, 'sealed_sha256': sealed(), 'own_sha256': own(),
               'nvrtc_library': str(NVRTC), 'nvrtc_library_sha256': sha(NVRTC), 'device_uuid_required': GPU_UUID, 'commands': []}
    d = out / 'test_wall_gpu_v7'
    d.mkdir()
    r = execute('gpu-test_wall_gpu_v7', [Path(br['build_dir']) / 'test_wall_gpu_v7', NVRTC, GPU_UUID, V7N, A1 / 'include', d], out, 3600, env)
    summary['commands'].append(r)
    if (d / 'wall-gpu-v7.json').exists():
        summary['test_wall_gpu_v7'] = summarize(d / 'wall-gpu-v7.json')
    summary['all_passed'] = r['returncode'] == 0
    (out / 'gpu-receipt.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({'gpu': str(out), 'all_passed': summary['all_passed']}))
    return 0 if summary['all_passed'] else 1


if __name__ == '__main__':
    a = sys.argv
    if len(a) < 3 or a[1] not in ('build', 'run', 'gpu') or '/' in a[2]:
        raise SystemExit(__doc__)
    if a[1] == 'build':
        raise SystemExit(build(a[2], a[3] if len(a) > 3 else 'normal'))
    if a[1] == 'run':
        raise SystemExit(run(a[2], a[3:]))
    raise SystemExit(gpu(a[2]))
