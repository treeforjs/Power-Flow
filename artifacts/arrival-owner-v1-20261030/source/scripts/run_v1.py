#!/usr/bin/env python3
"""Build and run the arrival owner v1 qualification; write receipts.

Builds run under  core-lease-20261030/both-cores.sh run-frozen <argv>   (cores 6,7; nice 15).
CPU runs run under core-lease-20261030/any-core.sh --wait env CUDA_VISIBLE_DEVICES= run-frozen <argv>
                   (core 6 or 7, nice 15, one thread, no GPU).
GPU runs run under flock testing-gpu.lock flock testing-gpu-<UUID>.lock any-core.sh --wait
                   env CUDA_VISIBLE_DEVICES=<UUID> run-frozen <argv>.
Every step uses a working directory under /home/treefor/projects (run-frozen keeps $PWD).
Start this script restricted to cores 6-7:  taskset -c 6,7 python3 -I run_v1.py LABEL [stage ...]
Stages: build, core, source, native, gpu, oracle, probe (default: all but probe, in that order).
Each stage appends to OUT/LABEL/receipt.json (logs and outputs in OUT/LABEL/).
"""
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

R = Path('/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030')
BUILD = R / 'build'
OUT_ROOT = Path('/home/treefor/projects/.workspace/outputs/pic-arrival-owner-v1-20261030')
LEASE = Path('/home/treefor/projects/runtimes/core-lease-20261030')
LOCKS = Path('/home/treefor/projects/.workspace/outputs/pic-profile')
GPU_UUID = 'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9'
GPU_LOCKS = [LOCKS / 'testing-gpu.lock', LOCKS / ('testing-gpu-%s.lock' % GPU_UUID)]
FROZEN = '/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
CXX = '/usr/bin/c++'
MPICXX = '/usr/bin/mpic++'

# Sealed prerequisites (read only; hashes verified before every stage).
COUNT = Path('/home/treefor/projects/runtimes/pic-water-exact-count-subflow-20261029/native/source/include/raikou/exact_particle_count.h')
COUNT_SHA = '46e01fe22e67cc9deb0da9320edef4041f4b3ef433de8ce1eca484463b405501'
S5 = Path('/home/treefor/projects/runtimes/pic-electrode-heat-source-v5-20261030')
S5_MOD = S5 / 'native/mfem_electrode_heat_source_v5'
S5_SEAL = Path('/home/treefor/projects/Power-Flow/artifacts/source-owner-v5-20261030/source-seal-v5.json')
V3 = Path('/home/treefor/projects/runtimes/pic-electrode-heat-feedback-20261029/native/mfem_electrode_heat_source_v3')
V2 = Path('/home/treefor/projects/runtimes/pic-electrode-heat-feedback-20261029/native/mfem_electrode_heat_feedback_v2')
MFEM = Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
MFEM_LIBS = [MFEM / 'install-mpi/lib/libmfem.a', MFEM / 'tpl/hypre/lib/libHYPRE.a', MFEM / 'tpl/metis/lib/libmetis.a', MFEM / 'tpl/GKlib/lib/libGKlib.a']
MFEM_INC = ['-isystem', str(MFEM / 'install-mpi/include'), '-isystem', str(MFEM / 'tpl/hypre/include'), '-isystem', str(MFEM / 'tpl/metis/include')]
N5 = Path('/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-full-time-candidate005')
N5_LIB = N5 / 'build/libpic_native_v5.a'
N5_LIB_SHA = '715729114e5a7d3975bd7ee7db92c570e61f78c907f3e999812b34e7a74159f2'
N5_IDENTITY = 'cca72822c9b61900aae7572f29c5d57314f45aef17c4b01c9fd2346c6fa3c883'
N5_VENDOR = Path('/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/source/vendor')
GPU = Path('/home/treefor/projects/runtimes/hall-gpu-resident-v1-20261030')
NVRTC_DIR = Path('/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc')
NVRTC = NVRTC_DIR / 'lib/libnvrtc.so.12'
CUDA_INC = Path('/home/treefor/projects/runtimes/mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include')

FP = ['-ffp-contract=off', '-fno-fast-math']
STRICT = ['-Wall', '-Wextra', '-Wpedantic', '-Werror']
SAN = ['-O1', '-g', '-fno-omit-frame-pointer', '-fsanitize=address,undefined', '-fno-sanitize-recover=all']
SAN_ENV = ['ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:abort_on_error=1', 'UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1']
COUNT_DEF = '-DPF_ARRIVAL_EXACT_COUNT_HEADER="%s"' % COUNT


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def prerequisites():
    """Hashes of every sealed input; refuses if a pinned one changed."""
    if sha(COUNT) != COUNT_SHA:
        raise RuntimeError('ExactParticleCount header changed')
    if sha(N5_LIB) != N5_LIB_SHA:
        raise RuntimeError('native v5 archive changed')
    seal = json.loads(S5_SEAL.read_text())
    changed = [f['path'] for f in seal['files'] if sha(S5 / f['path']) != f['sha256']]
    if changed:
        raise RuntimeError('source v5 sealed files changed: %s' % changed[:4])
    gpu = {str(p): sha(p) for p in sorted((GPU / 'common').rglob('*')) if p.is_file()}
    own = {str(p): sha(p) for p in sorted(R.rglob('*')) if p.is_file() and 'build' not in p.parts and '__pycache__' not in p.parts}
    return {'exact_count_header': {'path': str(COUNT), 'sha256': COUNT_SHA},
            'source_v5_seal': {'path': str(S5_SEAL), 'sha256': sha(S5_SEAL), 'source_identity': seal['source_identity'], 'files_verified': len(seal['files'])},
            'native_v5_archive': {'path': str(N5_LIB), 'sha256': N5_LIB_SHA, 'source_identity': N5_IDENTITY},
            'gpu_common': gpu, 'nvrtc': {'path': str(NVRTC), 'sha256': sha(NVRTC)}, 'component_files': own}


class Runner:
    def __init__(self, label):
        self.out = OUT_ROOT / label
        self.out.mkdir(parents=True, exist_ok=True)
        (self.out / 'logs').mkdir(exist_ok=True)
        self.receipt_path = self.out / 'receipt.json'
        self.receipt = json.loads(self.receipt_path.read_text()) if self.receipt_path.exists() else {'label': label, 'stages': []}

    def step(self, name, argv, kind, env=(), timeout=3600, cwd=None):
        if kind == 'build':
            cmd = [str(LEASE / 'both-cores.sh'), FROZEN, *argv]
        elif kind == 'cpu':
            cmd = [str(LEASE / 'any-core.sh'), '--wait', 'env', 'CUDA_VISIBLE_DEVICES=', *env, FROZEN, *argv]
        elif kind == 'gpu':
            cmd = ['flock', str(GPU_LOCKS[0]), 'flock', str(GPU_LOCKS[1]), str(LEASE / 'any-core.sh'), '--wait', 'env',
                   'CUDA_VISIBLE_DEVICES=' + GPU_UUID, *env, FROZEN, *argv]
        else:
            raise ValueError(kind)
        log = self.out / 'logs' / (name + '.log')
        t0 = time.time()
        with log.open('w') as f:
            try:
                rc = subprocess.run(cmd, cwd=str(cwd or R), stdout=f, stderr=subprocess.STDOUT, timeout=timeout).returncode
            except subprocess.TimeoutExpired:
                rc = 'timeout'
        entry = {'name': name, 'argv': cmd, 'returncode': rc, 'elapsed_s': round(time.time() - t0, 3), 'log': str(log), 'log_sha256': sha(log)}
        print(name, rc, entry['elapsed_s'], flush=True)
        return entry

    def stage(self, name, entries, extra=None):
        record = {'stage': name, 'utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'host': platform.node(),
                  'kernel': platform.release(), 'prerequisites': prerequisites(), 'steps': entries}
        if extra:
            record.update(extra)
        self.receipt['stages'].append(record)
        self.receipt_path.write_text(json.dumps(self.receipt, indent=1) + '\n')
        failed = [e['name'] for e in entries if e['returncode'] != 0]
        if failed:
            raise SystemExit('stage %s failed: %s' % (name, failed))


def core_cmd(sanitized):
    opt = SAN if sanitized else ['-O2']
    out = BUILD / ('test_core_sanitized' if sanitized else 'test_core')
    return [CXX, '-std=c++20', *opt, *STRICT, *FP, '-I' + str(R / 'include'), COUNT_DEF, str(R / 'tests/test_core.cpp'), '-lcrypto', '-o', str(out)]


def source_cmds(sanitized):
    tag = 'sanitized' if sanitized else 'normal'
    opt = SAN if sanitized else ['-O2']
    obj = BUILD / ('source-' + tag)
    inc = ['-I' + str(V3 / 'hall-owner/include'), '-I' + str(V3 / 'hall-owner/src'), '-I' + str(S5_MOD), '-I' + str(V3), '-I' + str(V2),
           '-I' + str(V2 / 'source-kernel/include'), *MFEM_INC]
    count = '-DHEAT_SOURCE_EXACT_COUNT_HEADER="%s"' % COUNT
    # Sealed supporting sources compiled exactly as the source v5 CMake targets do (no sanitizer, no -Werror).
    cmds = [[CXX, '-std=gnu++20', '-O2', *FP, '-I' + str(V2 / 'source-kernel/include'), '-c', str(V2 / 'source-kernel/src/surface_emission.cpp'), '-o', str(obj / 'surface_emission.o')],
            [CXX, '-std=gnu++20', '-O2', *FP, '-I' + str(V2 / 'source-kernel/include'), '-c', str(V2 / 'source-kernel/src/micro_layer.cpp'), '-o', str(obj / 'micro_layer.o')],
            [CXX, '-std=gnu++20', '-O2', *FP, *inc, '-c', str(V3 / 'hall-owner/src/mfem_gap_pic.cpp'), '-o', str(obj / 'mfem_gap_pic.o')],
            [CXX, '-std=gnu++20', *opt, '-Wall', '-Wextra', '-Wno-misleading-indentation', *FP, count, COUNT_DEF, '-I' + str(R / 'include'), *inc, '-c', str(R / 'tests/test_source_v5.cpp'), '-o', str(obj / 'test_source_v5.o')]]
    link = [CXX, *(['-fsanitize=address,undefined'] if sanitized else []), str(obj / 'test_source_v5.o'), str(obj / 'mfem_gap_pic.o'), str(obj / 'surface_emission.o'),
            str(obj / 'micro_layer.o'), *map(str, MFEM_LIBS), '/usr/lib/libmpi.so', '-lcrypto', '-o', str(BUILD / ('test_source_v5_' + tag))]
    return obj, cmds + [link]


def native_cmds(sanitized):
    tag = 'sanitized' if sanitized else 'normal'
    opt = SAN if sanitized else ['-O2']
    inc = ['-I' + str(N5 / 'source/include'), '-isystem', str(N5 / 'inputs'), '-I' + str(N5_VENDOR), *MFEM_INC]
    # The native include path provides the same ExactParticleCount bytes (hash pinned).
    obj = BUILD / ('native-' + tag)
    comp = [MPICXX, '-std=c++20', *opt, '-Wall', '-Wextra', *FP, '-DPF_ARRIVAL_EXACT_COUNT_HEADER="raikou/exact_particle_count.h"', '-I' + str(R / 'include'), *inc,
            '-c', str(R / 'tests/test_native_census.cpp'), '-o', str(obj / 'test_native_census.o')]
    link = [MPICXX, *(['-fsanitize=address,undefined'] if sanitized else []), '-Wl,--gc-sections', str(obj / 'test_native_census.o'), str(N5_LIB), *map(str, MFEM_LIBS), '-lcrypto',
            '-o', str(BUILD / ('test_native_census_' + tag))]
    return obj, [comp, link]


def gpu_cmd():
    return [CXX, '-std=c++20', '-O2', '-Wall', '-Wextra', '-Wno-c++26-extensions', *FP, '-I' + str(R / 'include'), '-I' + str(GPU), '-isystem', str(CUDA_INC),
            '-isystem', str(NVRTC_DIR / 'include'), str(R / 'tests/gpu_parity.cpp'), str(GPU / 'common/hall_gpu_common.cpp'), '-ldl', '-lcrypto', '-o', str(BUILD / 'gpu_parity')]


def stage_build(run, which):
    BUILD.mkdir(exist_ok=True)
    entries = []
    if 'core' in which:
        entries.append(run.step('build-core-normal', core_cmd(False), 'build'))
        entries.append(run.step('build-core-sanitized', core_cmd(True), 'build'))
    for name, fn in [('source', source_cmds), ('native', native_cmds)]:
        if name not in which:
            continue
        for sanitized in (False, True):
            obj, cmds = fn(sanitized)
            obj.mkdir(parents=True, exist_ok=True)
            # One lease for the whole target: a shell script runs the commands in order.
            script = obj / 'build.sh'
            script.write_text('#!/bin/sh\nset -eux\n' + ''.join(' '.join("'%s'" % a for a in c) + '\n' for c in cmds))
            script.chmod(0o755)
            entries.append(run.step('build-%s-%s' % (name, 'sanitized' if sanitized else 'normal'), ['/bin/sh', str(script)], 'build', timeout=5400))
    if 'gpu' in which:
        entries.append(run.step('build-gpu-parity', gpu_cmd(), 'build'))
    binaries = {str(p): sha(p) for p in sorted(BUILD.glob('test_*')) + sorted(BUILD.glob('gpu_parity')) if p.is_file()}
    run.stage('build', entries, {'targets': sorted(which), 'binaries': binaries})


def stage_core(run):
    entries = [run.step('core-normal', [str(BUILD / 'test_core'), str(run.out / 'core-normal')], 'cpu'),
               run.step('core-sanitized', [str(BUILD / 'test_core_sanitized'), str(run.out / 'core-sanitized')], 'cpu', env=SAN_ENV)]
    run.stage('core', entries, {'binaries': {p: sha(BUILD / p) for p in ['test_core', 'test_core_sanitized']}})


def stage_source(run):
    entries = [run.step('source-normal', [str(BUILD / 'test_source_v5_normal'), str(run.out / 'source-normal')], 'cpu', timeout=7200),
               run.step('source-sanitized', [str(BUILD / 'test_source_v5_sanitized'), str(run.out / 'source-sanitized')], 'cpu', env=SAN_ENV, timeout=14400)]
    run.stage('source', entries)


def stage_native(run):
    entries = [run.step('native-normal', [str(BUILD / 'test_native_census_normal'), str(run.out / 'native-normal')], 'cpu'),
               run.step('native-sanitized', [str(BUILD / 'test_native_census_sanitized'), str(run.out / 'native-sanitized')], 'cpu', env=SAN_ENV)]
    run.stage('native', entries)


def stage_probe(run):
    """Large neutral stocks against the sealed source v5 (sealed-code behaviour, recorded only)."""
    entries = []
    for n0 in ['1e20', '1e22', '1e25', '1e30', '1e100', '1e200']:
        for rate in ['1e300', '1e14']:
            e = run.step('probe-%s-%s' % (n0, rate), [str(BUILD / 'test_source_v5_normal'), str(run.out / 'probe'), 'probe', n0, rate, '1'], 'cpu')
            e['expected_returncode'] = 'recorded'
            entries.append(e)
    record = {'stage': 'probe', 'utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'steps': entries,
              'note': 'one neutral release per run; nonzero exit = sealed source/heat engine abort, recorded not gated'}
    run.receipt['stages'].append(record)
    run.receipt_path.write_text(json.dumps(run.receipt, indent=1) + '\n')


def stage_gpu(run, records):
    entries = []
    for k in (1, 2):
        entries.append(run.step('gpu-parity-process%d' % k, [str(BUILD / 'gpu_parity'), str(run.out / ('gpu-process%d' % k)), str(NVRTC), GPU_UUID,
                                                            str(records), '3'], 'gpu'))
    run.stage('gpu', entries, {'records': records, 'gpu_uuid': GPU_UUID})


def stage_oracle(run):
    entries = [run.step('oracle', ['/usr/bin/python3', '-I', str(R / 'tests/oracle.py'), str(run.out / 'core-normal/oracle-input.json'),
                                   str(run.out / 'oracle-result.json')], 'cpu')]
    run.stage('oracle', entries)


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    if any(c not in (6, 7) for c in os.sched_getaffinity(0)):
        raise SystemExit('restrict the driver to cores 6-7 (taskset -c 6,7)')
    label = sys.argv[1]
    stages = sys.argv[2:] or ['build', 'core', 'source', 'native', 'gpu', 'oracle']
    run = Runner(label)
    for s in stages:
        if s.startswith('build'):
            which = set(s.split(':')[1].split(',')) if ':' in s else {'core', 'source', 'native', 'gpu'}
            stage_build(run, which)
        elif s == 'core':
            stage_core(run)
        elif s == 'source':
            stage_source(run)
        elif s == 'native':
            stage_native(run)
        elif s.startswith('gpu'):
            stage_gpu(run, int(s.split(':')[1]) if ':' in s else 1 << 20)
        elif s == 'oracle':
            stage_oracle(run)
        elif s == 'probe':
            stage_probe(run)
        else:
            raise SystemExit('unknown stage ' + s)


if __name__ == '__main__':
    main()
