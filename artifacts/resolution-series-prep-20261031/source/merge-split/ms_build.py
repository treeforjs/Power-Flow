#!/usr/bin/python -I
"""Build and run the merge/split component and its tests (receipted).

  build:  python -I ms_build.py build LABEL            (run under a core lease, e.g. lent core 3)
  run:    python -I ms_build.py run LABEL NAME BINARY [ARGS ...]

Every compile and every test executes inside arch-base-20261028/run-frozen with cwd in this
directory. ARGS may use the placeholders {OUT} (run directory), {NVRTC} (pinned NVRTC 12.8
library) and {UUID} (GPU device UUID). Binaries: build/LABEL/. Logs and receipts:
/home/treefor/projects/.workspace/outputs/resolution-series-prep-20261031/merge-split/LABEL/.
Read-only inputs (philox core, native ExactParticleCount header, hall-gpu-resident-v1 common,
NVRTC) are hashed before and after each stage. Core leases and GPU locks are taken by the
caller (flock/taskset), not by this script.
"""
from pathlib import Path
import hashlib, json, os, subprocess, sys, time

P = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
OUT = Path('/home/treefor/projects/.workspace/outputs/resolution-series-prep-20261031/merge-split')
FROZEN = R / 'arch-base-20261028/run-frozen'
PHILOX = R / 'pic-arrival-owner-v1-20261030/include/pf_arrival'
NATIVE_INC = R / 'pic-reaction-transactions-20261029/native-full-time-candidate006/source/include'
COMMON = R / 'hall-gpu-resident-v1-20261030/common'
NVRTC = R / 'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/lib/libnvrtc.so.12'
NVRTC_INC = R / 'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/include'
CUDA_INC = R / 'mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include'
UUID = 'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9'
SEALED = [PHILOX / 'philox_core.h', NATIVE_INC / 'raikou/exact_particle_count.h', COMMON / 'hall_gpu_common.cpp',
          COMMON / 'hall_gpu_common.hpp', COMMON / 'hall_gpu_primitives.cpp', COMMON / 'hall_gpu_primitives.hpp',
          COMMON / 'device/deterministic_primitives.cuh', NVRTC]
ENV = {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'}

CXX = ['/usr/bin/c++', '-std=c++20']
FP = ['-ffp-contract=off', '-fno-fast-math']
WARN = ['-Wall', '-Wextra']
INC = ['-I' + str(P / 'include'), '-I' + str(P / 'tests'), '-I' + str(P / 'gpu'), '-isystem', str(PHILOX), '-isystem', str(NATIVE_INC)]
GPU_INC = ['-I' + str(COMMON), '-isystem', str(CUDA_INC), '-isystem', str(NVRTC_INC), '-Wno-c++26-extensions', '--embed-dir=' + str(PHILOX)]
OPT = ['-O2']
SAN = ['-O1', '-g', '-fsanitize=address,undefined', '-fno-omit-frame-pointer', '-fno-sanitize-recover=all']

# object name -> (flags, source)
OBJECTS = {
    'merge_split': (CXX + OPT + FP + WARN + INC, P / 'src/merge_split.cpp'),
    'relocation': (CXX + OPT + FP + WARN + INC, P / 'src/relocation.cpp'),
    'test_cpu': (CXX + OPT + FP + WARN + INC, P / 'tests/test_cpu.cpp'),
    'merge_split_san': (CXX + SAN + FP + WARN + INC, P / 'src/merge_split.cpp'),
    'relocation_san': (CXX + SAN + FP + WARN + INC, P / 'src/relocation.cpp'),
    'test_cpu_san': (CXX + SAN + FP + WARN + INC, P / 'tests/test_cpu.cpp'),
    'ms_gpu': (CXX + OPT + FP + WARN + INC + GPU_INC, P / 'gpu/ms_gpu.cpp'),
    'hall_gpu_common': (CXX + OPT + FP + GPU_INC, COMMON / 'hall_gpu_common.cpp'),
    'hall_gpu_primitives': (CXX + OPT + FP + GPU_INC, COMMON / 'hall_gpu_primitives.cpp'),
    'test_gpu': (CXX + OPT + FP + WARN + INC + GPU_INC, P / 'tests/test_gpu.cpp'),
}
TARGETS = {
    'test_cpu': (['test_cpu', 'merge_split', 'relocation'], ['-lcrypto'], []),
    'test_cpu_asan': (['test_cpu_san', 'merge_split_san', 'relocation_san'], ['-lcrypto'], ['-fsanitize=address,undefined']),
    'test_gpu': (['test_gpu', 'ms_gpu', 'merge_split', 'relocation', 'hall_gpu_common', 'hall_gpu_primitives'], ['-ldl', '-lcrypto'], []),
}


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def own_sources():
    files = sorted(x for d in ('include', 'src', 'gpu', 'tests') for x in (P / d).rglob('*') if x.is_file())
    return {str(x): sha(x) for x in files + [Path(__file__).resolve()]}


def sealed():
    return {str(x): sha(x) for x in SEALED}


def admission():
    return {'cpu_affinity': sorted(os.sched_getaffinity(0)), 'nice': os.getpriority(os.PRIO_PROCESS, 0), 'pid': os.getpid(),
            'time_unix': time.time(), 'frozen_launcher_sha256': sha(FROZEN),
            'return_cores_flag': Path('/home/treefor/projects/.workspace/outputs/pic-profile/RETURN-CORES-0-3').exists()}


def execute(name, cmd, log_dir, timeout):
    log = log_dir / (name + '.log')
    t = time.monotonic()
    with log.open('w') as h:
        try:
            rc = subprocess.run([str(FROZEN), *map(str, cmd)], cwd=P, env=ENV, stdout=h, stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    r = {'name': name, 'argv': [str(FROZEN), *map(str, cmd)], 'returncode': rc, 'elapsed_s': time.monotonic() - t, 'log': str(log), 'log_sha256': sha(log)}
    print(name, rc, round(r['elapsed_s'], 1), flush=True)
    return r


def build(label, targets):
    out = OUT / label / 'build'
    out.mkdir(parents=True, exist_ok=False)
    b = P / 'build' / label
    b.mkdir(parents=True, exist_ok=False)
    pins, before = sealed(), own_sources()
    receipt = {'schema': 'merge-split-build-receipt-v1', 'label': label, 'admission': admission(), 'sealed_sha256': pins, 'commands': []}
    status = 'passed'
    try:
        need = sorted({o for t in targets for o in TARGETS[t][0]})
        for o in need:
            fl, src = OBJECTS[o]
            r = execute('compile-' + o, [*fl, '-MD', '-MF', b / (o + '.d'), '-c', src, '-o', b / (o + '.o')], out, 3600)
            receipt['commands'].append(r)
            if r['returncode'] != 0:
                raise RuntimeError(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-6000:])
        for t in targets:
            objs, libs, extra = TARGETS[t]
            r = execute('link-' + t, ['/usr/bin/c++', *extra, *[b / (o + '.o') for o in objs], *libs, '-o', b / t], out, 1800)
            receipt['commands'].append(r)
            if r['returncode'] != 0:
                raise RuntimeError(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-6000:])
        deps = set()
        for o in need:
            text = (b / (o + '.d')).read_text().replace('\\\n', ' ')
            for line in text.splitlines():
                if ':' in line:
                    deps |= {str(Path(x).resolve()) for x in line.split(':', 1)[1].split()}
        receipt['dependency_sha256'] = {d: sha(d) for d in sorted(deps)}
        receipt['binaries_sha256'] = {t: sha(b / t) for t in targets}
        receipt['build_dir'] = str(b)
        receipt['own_sources_sha256'] = before
        if own_sources() != before or sealed() != pins:
            raise RuntimeError('sources or read-only inputs changed during the build')
    except Exception as e:
        status = 'failed'
        receipt['error'] = str(e)
    receipt['status'] = status
    (out / 'build-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    print(json.dumps({'build': str(out), 'status': status, 'error': receipt.get('error', '')[:6000]}))
    return 0 if status == 'passed' else 1


def run(label, name, binary, args):
    br = json.loads((OUT / label / 'build/build-receipt.json').read_text())
    if br['status'] != 'passed':
        raise SystemExit('build did not pass')
    exe = Path(br['build_dir']) / binary
    if sha(exe) != br['binaries_sha256'][binary]:
        raise SystemExit('binary differs from its build receipt')
    out = OUT / label / ('%s-%s' % (name, time.strftime('%H%M%S')))
    out.mkdir(parents=True)
    sub = {'{OUT}': str(out), '{NVRTC}': str(NVRTC), '{UUID}': UUID}
    def fill(a):
        for k, v in sub.items():
            a = a.replace(k, v)
        return a
    pins = sealed()
    adm = admission()
    r = execute(name, [exe, *[fill(a) for a in args]], out, 4 * 3600)
    receipt = {'schema': 'merge-split-run-receipt-v1', 'label': label, 'name': name, 'admission': adm, 'binary_sha256': sha(exe),
               'own_sources_unchanged_since_build': own_sources() == br['own_sources_sha256'], 'sealed_unchanged': sealed() == pins,
               'command': r}
    (out / 'run-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'returncode': r['returncode']}))
    return 0 if r['returncode'] == 0 else 1


if __name__ == '__main__':
    a = sys.argv
    if len(a) >= 3 and a[1] == 'build' and '/' not in a[2]:
        raise SystemExit(build(a[2], a[3:] or list(TARGETS)))
    if len(a) >= 5 and a[1] == 'run' and '/' not in a[2]:
        raise SystemExit(run(a[2], a[3], a[4], a[5:]))
    raise SystemExit(__doc__)
