#!/usr/bin/env python3
"""Qualification of chemistry v3: builds, CPU checks (normal + sanitized), CPU/GPU parity, swarm validation.

usage: taskset -c 6,7 python3 -I run_checks.py LABEL [stage ...]
stages: inputs tables build core sanitized gpu swarm report   (default: all, in that order)
Receipts and logs: /home/treefor/projects/.workspace/outputs/pic-water-chemistry-v3-20261031/LABEL/

Builds and CPU runs execute inside the frozen launcher on one leased core (tools/lease_any.sh: cores 6, 7 and the
lent cores 0, 1, 3 with their lease files; nice 15, one thread, no GPU). Swarm jobs use tools/run_jobs.py (cores 6, 7 and the lent cores
0, 1, 3 with their lease files). The GPU stage holds testing-gpu.lock and the device lock for minutes only.
"""
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

R = Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3-20261031')
OUT_ROOT = Path('/home/treefor/projects/.workspace/outputs/pic-water-chemistry-v3-20261031')
LEASE = Path('/home/treefor/projects/runtimes/core-lease-20261030/any-core.sh')
FROZEN = '/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
LOCKS = Path('/home/treefor/projects/.workspace/outputs/pic-profile')
GPU_UUID = 'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9'
NVRTC = Path('/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/lib/libnvrtc.so.12')
PINNED = {
    # read-in-place dependencies (sha256 at qualification time)
    '/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include/pf_arrival/philox_core.h': None,
    '/home/treefor/projects/runtimes/hall-gpu-resident-v1-20261030/common/hall_gpu_common.cpp': None,
    '/home/treefor/projects/runtimes/hall-gpu-resident-v1-20261030/common/hall_gpu_common.hpp': None,
    str(NVRTC): None,
}
EN_LIST = [30, 60, 100, 150, 200, 300, 500, 700, 1000, 1500, 2000, 3000]


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


class Run:
    def __init__(self, label):
        self.out = OUT_ROOT / label
        self.out.mkdir(parents=True, exist_ok=True)
        self.rp = self.out / 'receipt.json'
        self.receipt = json.loads(self.rp.read_text()) if self.rp.exists() else {'label': label, 'host': platform.node(), 'stages': []}

    def step(self, name, argv, kind='cpu', timeout=7200):
        log = self.out / ('%s.log' % name)
        if kind == 'cpu':
            cmd = [str(R / 'tools/lease_any.sh'), 'env', 'CUDA_VISIBLE_DEVICES=', 'OMP_NUM_THREADS=1', FROZEN, *argv]
        elif kind == 'gpu':
            cmd = ['flock', str(LOCKS / 'testing-gpu.lock'), 'flock', str(LOCKS / ('testing-gpu-%s.lock' % GPU_UUID)), str(LEASE), '--wait', 'env',
                   'CUDA_VISIBLE_DEVICES=' + GPU_UUID, FROZEN, *argv]
        elif kind == 'host':
            cmd = list(argv)
        t0 = time.time()
        with open(log, 'w') as f:
            p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=str(R), timeout=timeout)
        e = {'name': name, 'argv': cmd, 'exit': p.returncode, 'seconds': round(time.time() - t0, 2), 'log': str(log)}
        print('%-28s exit %d  %.1f s' % (name, p.returncode, e['seconds']), flush=True)
        return e

    def stage(self, name, entries, extra=None):
        self.receipt['stages'].append({'stage': name, 'utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'entries': entries,
                                       'passed': all(e.get('exit', 0) == 0 for e in entries), **(extra or {})})
        self.rp.write_text(json.dumps(self.receipt, indent=1))


def stage_inputs(run):
    pins = {p: sha(p) for p in PINNED}
    man = json.load(open(R / 'data/tables/chemistry-v3-manifest.json'))
    files = {}
    for k, v in man['input_files_sha256'].items():
        p = k if os.path.isabs(k) else str(R / k)
        files[k] = {'expected': v, 'actual': sha(p)}
    bad = [k for k, v in files.items() if v['expected'] != v['actual']]
    own = {str(p.relative_to(R)): sha(p) for p in sorted(R.rglob('*')) if p.is_file() and 'build' not in p.parts and '__pycache__' not in p.parts
           and 'data' not in p.parts}
    run.stage('inputs', [{'name': 'table input files unchanged', 'exit': 1 if bad else 0, 'changed': bad}],
              {'pinned_dependencies': pins, 'table_inputs': files, 'component_sources': own})


def stage_tables(run):
    tmp = run.out / 'tables'
    tmp.mkdir(exist_ok=True)
    e1 = run.step('spec-regenerate', ['python3', '-I', str(R / 'tools/make_production_spec.py'), str(tmp / 'chemistry-v3-spec.json')], 'host')
    same_spec = (tmp / 'chemistry-v3-spec.json').read_bytes() == (R / 'data/tables/chemistry-v3-spec.json').read_bytes()
    e2 = run.step('tables-rebuild', ['python3', '-I', str(R / 'tools/build_tables.py'), str(R / 'data/tables/chemistry-v3-spec.json'), str(tmp / 'chemistry-v3')], 'host')
    a, b = sha(tmp / 'chemistry-v3.blob'), sha(R / 'data/tables/chemistry-v3.blob')
    run.stage('tables', [e1, e2, {'name': 'spec regenerated byte-identical', 'exit': 0 if same_spec else 1},
                         {'name': 'blob rebuilt byte-identical', 'exit': 0 if a == b else 1, 'rebuilt': a, 'shipped': b}], {'blob_sha256': b})


def stage_build(run):
    bdir = R / 'build/bin'
    bdir.mkdir(parents=True, exist_ok=True)
    es = [run.step('build-test_core', [str(R / 'tools/dev_build.sh'), 'normal', str(R / 'tests/test_core.cpp'), str(bdir / 'test_core')]),
          run.step('build-test_core_san', [str(R / 'tools/dev_build.sh'), 'sanitized', str(R / 'tests/test_core.cpp'), str(bdir / 'test_core_sanitized')]),
          run.step('build-swarm_mc', [str(R / 'tools/dev_build.sh'), 'normal', str(R / 'tests/swarm_mc.cpp'), str(bdir / 'swarm_mc')]),
          run.step('build-gpu_parity', [str(R / 'tools/dev_build_gpu.sh'), str(bdir / 'gpu_parity')])]
    run.stage('build', es, {'binaries': {p.name: sha(p) for p in sorted(bdir.iterdir()) if p.is_file()}})


def stage_core(run):
    blob, chk = R / 'data/tables/chemistry-v3.blob', R / 'data/tables/chemistry-v3-checks.json'
    e = run.step('core-normal', [str(R / 'build/bin/test_core'), str(blob), str(chk), str(run.out / 'core-normal'), '1000000'], timeout=10800)
    run.stage('core', [e], {'blob_sha256': sha(blob)})


def stage_sanitized(run):
    blob, chk = R / 'data/tables/chemistry-v3.blob', R / 'data/tables/chemistry-v3-checks.json'
    e = run.step('core-sanitized', ['env', 'ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:abort_on_error=1', 'UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1',
                                    str(R / 'build/bin/test_core_sanitized'), str(blob), str(chk), str(run.out / 'core-sanitized'), '100000'], timeout=14400)
    run.stage('sanitized', [e])


def stage_gpu(run):
    blob = R / 'data/tables/chemistry-v3.blob'
    es = [run.step('gpu-parity-process%d' % k, [str(R / 'build/bin/gpu_parity'), str(blob), str(run.out / ('gpu-process%d' % k)), str(NVRTC), GPU_UUID, '50000', '2'],
                   'gpu', timeout=3600) for k in (1, 2)]
    a = json.loads((run.out / 'gpu-process1/gpu-parity.json').read_text())
    b = json.loads((run.out / 'gpu-process2/gpu-parity.json').read_text())
    same = [c.get('record_digest_gpu') for c in a['cases']] == [c.get('record_digest_gpu') for c in b['cases']]
    run.stage('gpu', es + [{'name': 'two fresh processes give identical device digests', 'exit': 0 if same else 1}])


def stage_swarm(run):
    sdir = run.out / 'swarm'
    sdir.mkdir(exist_ok=True)
    # isotropic-elastic variant of the production set (angular-treatment comparison)
    spec = json.load(open(R / 'data/tables/chemistry-v3-spec.json'))
    for p in spec['groups'][0]['processes']:
        if p['kind'] == 'ELASTIC':
            p['angular'] = 'ISO'
            p.pop('sigma_is_momentum_transfer', None)
            p['id'] += ' [isotropic variant]'
    spec['id'] += '-isotropic-elastic-variant'
    json.dump(spec, open(sdir / 'v3-iso-spec.json', 'w'), indent=1)
    e0 = run.step('swarm-iso-tables', ['python3', '-I', str(R / 'tools/build_tables.py'), str(sdir / 'v3-iso-spec.json'), str(sdir / 'v3-iso')], 'host')
    sets = {'v3': str(R / 'data/tables/chemistry-v3.blob'), 'v3-iso': str(sdir / 'v3-iso.blob')}
    for s in ('muroranit', 'triniti', 'hayashi', 'phelps', 'morgan'):
        src = next(R.glob('data/sources/lxcat-%s-h2o-cs-*/Cross section.txt' % s))
        subprocess.run(['python3', '-I', str(R / 'tools/make_lxcat_spec.py'), str(src), s, str(sdir / ('lxcat-%s-spec.json' % s))], check=True, cwd=str(R))
        subprocess.run(['python3', '-I', str(R / 'tools/build_tables.py'), str(sdir / ('lxcat-%s-spec.json' % s)), str(sdir / ('lxcat-%s' % s))], check=True, cwd=str(R))
        sets['lxcat-' + s] = str(sdir / ('lxcat-%s.blob' % s))
    chunks = [EN_LIST[0:4], EN_LIST[4:8], EN_LIST[8:12]]
    jobs = []
    for name, blob in sets.items():
        for k, ch in enumerate(chunks):
            jobs.append({'name': '%s-%d' % (name, k), 'cwd': str(R), 'log': str(sdir / ('%s-%d.log' % (name, k))),
                         'argv': [str(R / 'build/bin/swarm_mc'), blob, str(sdir / ('%s-%d.json' % (name, k))), ','.join(map(str, ch)), '3.22e22', '300', '2000', '20261031']})
    json.dump(jobs, open(sdir / 'jobs.json', 'w'), indent=1)
    e1 = run.step('swarm-jobs', ['python3', '-I', str(R / 'tools/run_jobs.py'), str(sdir / 'jobs.json'), str(sdir / 'jobs-receipt.json')], 'host', timeout=21600)
    e2 = run.step('swarm-table', ['python3', '-I', str(R / 'tools/swarm_table.py'), str(sdir), str(run.out / 'SWARM-VALIDATION')], 'host')
    run.stage('swarm', [e0, e1, e2], {'sets': {k: sha(v) for k, v in sets.items()}})


def main():
    label = sys.argv[1]
    stages = sys.argv[2:] or ['inputs', 'tables', 'build', 'core', 'sanitized', 'gpu', 'swarm']
    run = Run(label)
    fns = {'inputs': stage_inputs, 'tables': stage_tables, 'build': stage_build, 'core': stage_core, 'sanitized': stage_sanitized,
           'gpu': stage_gpu, 'swarm': stage_swarm}
    for s in stages:
        fns[s](run)
    ok = all(st['passed'] for st in run.receipt['stages'])
    print('RECEIPT', run.rp, 'passed' if ok else 'FAILED')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
