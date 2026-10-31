#!/usr/bin/python -I
"""Build and run the SS304 phase thermal owner v1.1 and its qualification (receipted, frozen launcher).

  build:  <lease wrapper> /usr/bin/python -I run_owner_v1_1.py build LABEL [--sanitize]
  run:    <lease wrapper> /usr/bin/python -I run_owner_v1_1.py run LABEL NAME BINARY [ARGS ...]
  gpu:    flock <testing-gpu.lock> flock <testing-gpu-UUID.lock> <lease wrapper> \
              /usr/bin/python -I run_owner_v1_1.py gpu LABEL NAME BINARY [ARGS ...]

Lease wrappers: core-lease-20261030/any-core.sh --wait (cores 6/7), or for the lent cores 0/1/3 (CPU stages only)
  flock <pic-profile/testing-core<K>.lock> taskset -c <K> nice -n 15 env OMP_NUM_THREADS=1 ...
ARGS may use {OUT} (run directory), {DATA} (runtime data/), {NVRTC}, {UUID}. Every compile and run executes inside
arch-base-20261028/run-frozen. Outputs: /home/treefor/projects/.workspace/outputs/ss304-phase-thermal-owner-v1.1-20261031/
LABEL/. Sealed inputs (heat-v2 material headers, hall-gpu-resident-v1 common/, NVRTC) are hashed before and after.
"""
from pathlib import Path
import fcntl, hashlib, json, os, platform, subprocess, sys, time

P = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
OUT = Path('/home/treefor/projects/.workspace/outputs/ss304-phase-thermal-owner-v1.1-20261031')
F = R / 'arch-base-20261028/run-frozen'
G = R / 'hall-gpu-resident-v1-20261030'
SEALED_MATERIAL = R / 'hall-coupled-transport-20261029/source/development/hall-canonical-driver-20261029/bindings/heat-v2/material-v2'
NVRTC = R / 'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/lib/libnvrtc.so.12'
NVRTC_INC = R / 'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/include'
CUDA_INC = R / 'mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include'
LEASE_DIR = Path('/home/treefor/projects/.workspace/outputs/pic-profile')
LEASES = {6: LEASE_DIR / 'testing-core.lock', 7: LEASE_DIR / 'testing-core7.lock', 0: LEASE_DIR / 'testing-core0.lock',
          1: LEASE_DIR / 'testing-core1.lock', 3: LEASE_DIR / 'testing-core3.lock'}
GPU_UUID = 'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9'
GPU_LOCKS = [LEASE_DIR / 'testing-gpu.lock', LEASE_DIR / ('testing-gpu-%s.lock' % GPU_UUID)]
V1 = R / 'ss304-phase-thermal-owner-v1-20261031'
SEALED = [SEALED_MATERIAL / 'sealed1800.hpp', SEALED_MATERIAL / 'ss304_thermal_v2.hpp', G / 'common/hall_gpu_common.cpp',
          G / 'common/hall_gpu_common.hpp', NVRTC] + sorted((V1 / 'source').glob('*.h*'))
ENV = {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
       'PYTHONNOUSERSITE': '1'}
COMMON = ['/usr/bin/c++', '-std=c++20', '-O2', '-ffp-contract=off', '-fno-fast-math', '-Wall', '-Wextra']
SAN = ['-fsanitize=address,undefined', '-fno-sanitize-recover=undefined', '-fno-omit-frame-pointer', '-g']
GPU_INC = ['-I' + str(G), '-isystem', str(CUDA_INC), '-isystem', str(NVRTC_INC), '-Wno-c++26-extensions']


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def held(lock):
    with open(lock, 'a') as h:
        try:
            fcntl.flock(h, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(h, fcntl.LOCK_UN)
        return False


def admission(stage):
    cores = os.sched_getaffinity(0)
    if len(cores) != 1 or not cores <= set(LEASES) or os.getpriority(os.PRIO_PROCESS, 0) < 15:
        raise SystemExit('one leased core (6, 7, or lent 0/1/3) at nice >= 15 required; affinity %s' % sorted(cores))
    c = next(iter(cores))
    if c in (0, 1, 3) and (LEASE_DIR / 'RETURN-CORES-0-3').exists():
        raise SystemExit('cores 0-3 are being returned; no new work there')
    if stage == 'gpu' and c not in (6, 7):
        raise SystemExit('GPU stages run on cores 6/7 only')
    if not held(LEASES[c]):
        raise SystemExit('core %d lease is not held by the launching wrapper' % c)
    if stage == 'gpu' and not all(held(lk) for lk in GPU_LOCKS):
        raise SystemExit('both GPU lock files must be held')
    return {'stage': stage, 'core': c, 'nice': os.getpriority(os.PRIO_PROCESS, 0), 'gpu_locks_held': [held(lk) for lk in GPU_LOCKS],
            'python': platform.python_version(), 'frozen_launcher_sha256': sha(F), 'time_unix': time.time()}


def own_sources():
    files = sorted(p for d in ('source', 'tests', 'data') for p in (P / d).rglob('*') if p.is_file())
    return {str(p.relative_to(P)): sha(p) for p in files + [Path(__file__).resolve()]}


def sealed():
    return {str(p): sha(p) for p in SEALED}


def execute(name, cmd, log_dir, timeout, env=None):
    log = log_dir / (name + '.log')
    t = time.monotonic()
    with log.open('w') as h:
        try:
            rc = subprocess.run([str(F), *map(str, cmd)], cwd=P, env=env or ENV, stdout=h, stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    r = {'name': name, 'argv': [str(F), *map(str, cmd)], 'returncode': rc, 'elapsed_s': time.monotonic() - t, 'log': str(log), 'log_sha256': sha(log)}
    print(name, rc, round(r['elapsed_s'], 1), flush=True)
    return r


def build(label, sanitize):
    adm = admission('build')
    out = OUT / label / 'build'
    out.mkdir(parents=True, exist_ok=False)
    b = P / 'build' / label
    b.mkdir(parents=True, exist_ok=False)
    pins, before = sealed(), own_sources()
    flags = COMMON + (SAN if sanitize else []) + ['-Isource', '-Itests']
    jobs = [('test_law', [*flags, '-isystem', str(SEALED_MATERIAL), 'tests/test_law.cpp', '-lcrypto', '-o', b / 'test_law']),
            ('test_owner', [*flags, 'tests/test_owner.cpp', '-lcrypto', '-o', b / 'test_owner']),
            ('test_v11', [*flags, 'tests/test_v11.cpp', '-lcrypto', '-o', b / 'test_v11']),
            ('test_regression', [*flags, 'tests/test_regression.cpp', '-lcrypto', '-o', b / 'test_regression']),
            ('nvrtc_archs', [*flags, '-isystem', str(NVRTC_INC), '-Wno-c++26-extensions', 'tests/nvrtc_archs.cpp', '-ldl', '-lcrypto', '-o', b / 'nvrtc_archs']),
            ('hall_gpu_common', [*COMMON[:-2], *(SAN if sanitize else []), *GPU_INC, '-c', str(G / 'common/hall_gpu_common.cpp'), '-o', b / 'hall_gpu_common.o']),
            ('gpu_engine', [*flags, *GPU_INC, '-c', 'source/gpu_engine.cpp', '-o', b / 'gpu_engine.o']),
            ('test_gpu_o', [*flags, *GPU_INC, '-c', 'tests/test_gpu.cpp', '-o', b / 'test_gpu.o']),
            ('test_gpu', [*COMMON, *(SAN if sanitize else []), b / 'test_gpu.o', b / 'gpu_engine.o', b / 'hall_gpu_common.o', '-ldl', '-lcrypto', '-o', b / 'test_gpu'])]
    rec = {'schema': 'ss304-phase-owner-v1.1-build-receipt-v1', 'label': label, 'sanitize': sanitize, 'admission': adm, 'sealed_sha256': pins,
           'own_sources_sha256': before, 'commands': []}
    status = 'passed'
    for name, cmd in jobs:
        r = execute(name, cmd, out, 3600)
        rec['commands'].append(r)
        if r['returncode'] != 0:
            status = 'failed'
            rec['error'] = Path(r['log']).read_text()[-6000:]
            break
    if status == 'passed':
        rec['binaries_sha256'] = {n: sha(b / n) for n in ('test_law', 'test_owner', 'test_v11', 'test_regression', 'nvrtc_archs', 'test_gpu')}
        rec['compiler'] = subprocess.run([str(F), '/usr/bin/c++', '--version'], capture_output=True, text=True, env=ENV).stdout.splitlines()[0]
        if own_sources() != before or sealed() != pins:
            status = 'failed'
            rec['error'] = 'sources or sealed inputs changed during the build'
    rec['status'] = status
    rec['build_dir'] = str(b)
    (out / 'build-receipt.json').write_text(json.dumps(rec, indent=1) + '\n')
    print(json.dumps({'build': str(out), 'status': status}))
    return 0 if status == 'passed' else 1


def run(stage, label, name, binary, args):
    adm = admission(stage)
    br = json.loads((OUT / label / 'build/build-receipt.json').read_text())
    if br['status'] != 'passed':
        raise SystemExit('build did not pass')
    exe = Path(br['build_dir']) / binary
    if sha(exe) != br['binaries_sha256'][binary]:
        raise SystemExit('binary differs from its build receipt')
    pins = sealed()
    out = OUT / label / ('%s-%s-%s' % (stage, name, time.strftime('%H%M%S')))
    out.mkdir(parents=True)
    sub = {'{OUT}': str(out), '{DATA}': str(P / 'data'), '{NVRTC}': str(NVRTC), '{UUID}': GPU_UUID}
    def fill(a):
        for k, v in sub.items():
            a = a.replace(k, v)
        return a
    env = dict(ENV)
    if br.get('sanitize'):
        # the CUDA driver maps memory in the ASan shadow gap and keeps allocations until process exit
        env['ASAN_OPTIONS'] = 'protect_shadow_gap=0:detect_leaks=0' if stage == 'gpu' else 'detect_leaks=1'
        env['UBSAN_OPTIONS'] = 'print_stacktrace=1:halt_on_error=1'
    prefix = []
    if br.get('sanitize') and stage == 'gpu':
        # under ASan the dlopen interceptor hides libnvrtc's $ORIGIN, so libnvrtc-builtins is found by path
        prefix = ['/usr/bin/env', 'LD_LIBRARY_PATH=' + str(NVRTC.parent)]
    r = execute(name, [*prefix, exe, *[fill(a) for a in args]], out, 12 * 3600, env)
    rec = {'schema': 'ss304-phase-owner-v1.1-run-receipt-v1', 'stage': stage, 'label': label, 'name': name, 'admission': adm,
           'binary_sha256': sha(exe), 'build_receipt': str(OUT / label / 'build/build-receipt.json'),
           'own_sources_unchanged_since_build': own_sources() == br['own_sources_sha256'], 'sealed_unchanged': sealed() == pins,
           'environment': env, 'command': r}
    (out / 'run-receipt.json').write_text(json.dumps(rec, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'returncode': r['returncode']}))
    return 0 if r['returncode'] == 0 else 1


if __name__ == '__main__':
    a = sys.argv
    if len(a) < 3 or a[1] not in ('build', 'run', 'gpu') or '/' in a[2]:
        raise SystemExit(__doc__)
    if a[1] == 'build':
        raise SystemExit(build(a[2], '--sanitize' in a[3:]))
    raise SystemExit(run('gpu' if a[1] == 'gpu' else 'cpu', a[2], a[3], a[4], a[5:]))
