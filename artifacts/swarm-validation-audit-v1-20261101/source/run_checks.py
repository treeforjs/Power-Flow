#!/usr/bin/env python3
"""Qualification of chemistry v3: builds, CPU checks (normal + sanitized), CPU/GPU parity, swarm validation.

usage: taskset -c 6,7 python3 -I run_checks.py LABEL [stage ...]
stages: inputs tables policy build core sanitized gpu cost bench swarm   (default: all, in that order)
Receipts and logs: /home/treefor/projects/.workspace/outputs/pic-water-chemistry-v3.2-20261031/LABEL/

Builds and CPU runs execute inside the frozen launcher on one leased core (tools/lease_any.sh: cores 6, 7 and the
lent cores 0, 1 with their lease files; nice 15, one thread, no GPU). Swarm jobs use tools/run_jobs.py (cores 6, 7 and the lent cores
0, 1 with their lease files). The GPU stage holds testing-gpu.lock and the device lock for minutes only.
"""
import hashlib
import json
import importlib.util
import os
import platform
import signal
import subprocess
import sys
import time
from pathlib import Path

R = Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.2-20261031')
OUT_ROOT = Path('/home/treefor/projects/.workspace/outputs/pic-water-chemistry-v3.2-20261031')
FROZEN = '/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
LOCKS = Path('/home/treefor/projects/.workspace/outputs/pic-profile')
GPU_UUID = 'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9'
GPU_LEASE = '/home/treefor/projects/runtimes/simqueue-0.6.2-20261031/bin/gpu-lease'
V31 = Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1-20261031')   # read-only baseline (sealed)
V3_SWARM = Path('/home/treefor/projects/.workspace/outputs/pic-water-chemistry-v3.1-20261031/qual-001/SWARM-VALIDATION.json')   # = v3 values
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


_qspec = importlib.util.spec_from_file_location('qualification_integrity', R / 'tools/qualification_integrity.py')
QI = importlib.util.module_from_spec(_qspec)
_qspec.loader.exec_module(QI)


class Run:
    def __init__(self, label):
        self.out = OUT_ROOT / label
        self.out.mkdir(parents=True, exist_ok=True)
        self.rp = self.out / 'receipt.json'
        self.receipt = json.loads(self.rp.read_text()) if self.rp.exists() else {'schema': QI.SCHEMA, 'label': label, 'host': platform.node(), 'stages': []}
        QI.need(self.receipt.get('schema') == QI.SCHEMA, 'legacy receipt requires a new qualification label')
        self.binding = None

    def step(self, name, argv, kind='cpu', timeout=7200):
        log = self.out / ('%s.log' % name)
        if kind == 'cpu':
            cmd = [str(R / 'tools/lease_any.sh'), 'env', 'CUDA_VISIBLE_DEVICES=', 'OMP_NUM_THREADS=1', FROZEN, *argv]
        elif kind == 'gpu':
            # the node's GPU lease helper: GPU lease files first, then a core lease (fair line with SimQueue GPU jobs)
            cmd = [GPU_LEASE, '--core', 'any', '--label', 'pf-chem-v3.2 ' + name, '--', 'env', 'CUDA_VISIBLE_DEVICES=' + GPU_UUID, FROZEN, *argv]
        elif kind == 'host':
            cmd = list(argv)
        t0 = time.time()
        with open(log, 'w') as f:
            p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=str(R), start_new_session=True)
            try:
                p.wait(timeout=timeout)
            except BaseException:
                # Allow run_jobs to reap its own child sessions before a forced stop.
                # This group belongs to this invocation; never signal a guessed PID.
                try: os.killpg(p.pid, signal.SIGTERM)
                except ProcessLookupError: pass
                try:
                    p.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    try: os.killpg(p.pid, signal.SIGKILL)
                    except ProcessLookupError: pass
                    p.wait(timeout=5)
                raise
        e = {'name': name, 'argv': cmd, 'exit': p.returncode, 'seconds': round(time.time() - t0, 2), 'log': str(log)}
        print('%-28s exit %d  %.1f s' % (name, p.returncode, e['seconds']), flush=True)
        return e

    def begin_stage(self, name):
        dep_roots = [V31 / 'include',
            Path('/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include'),
            Path('/home/treefor/projects/runtimes/hall-gpu-resident-v1-20261030/common'),
            Path('/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/include'),
            Path('/home/treefor/projects/runtimes/mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include')]
        deps = sorted(set(PINNED) | {str(p) for d in dep_roots for p in d.rglob('*') if p.is_file()} |
                      {str(V31 / 'build/bin/gpu_parity'), str(V31 / 'data/tables/chemistry-v3.1.blob')})
        self.dependencies = deps
        self.binding = QI.begin_stage(R, deps, self.receipt, name)
        os.environ['PF_QUAL_BINDING_SHA256'] = self.binding['before']['sha256']

    def stage(self, name, entries, extra=None):
        binding = QI.finish_stage(R, self.dependencies, self.binding)
        self.receipt['stages'].append({'stage': name, 'utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'entries': entries,
                                       'passed': bool(entries) and all(type(e.get('exit')) is int and e['exit'] == 0 for e in entries), 'binding': binding, **(extra or {})})
        self.rp.write_text(json.dumps(self.receipt, indent=1))


def stage_inputs(run):
    pins = {p: sha(p) for p in PINNED}
    man = json.load(open(R / 'data/tables/chemistry-v3.2-manifest.json'))
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
    e1 = run.step('spec-regenerate', ['python3', '-I', str(R / 'tools/make_production_spec.py'), str(tmp / 'chemistry-v3.2-spec.json')], 'host')
    same_spec = (tmp / 'chemistry-v3.2-spec.json').read_bytes() == (R / 'data/tables/chemistry-v3.2-spec.json').read_bytes()
    e2 = run.step('tables-rebuild', ['python3', '-I', str(R / 'tools/build_tables.py'), str(R / 'data/tables/chemistry-v3.2-spec.json'), str(tmp / 'chemistry-v3.2')], 'host')
    a, b = sha(tmp / 'chemistry-v3.2.blob'), sha(R / 'data/tables/chemistry-v3.2.blob')
    run.stage('tables', [e1, e2, {'name': 'spec regenerated byte-identical', 'exit': 0 if same_spec else 1},
                         {'name': 'blob rebuilt byte-identical', 'exit': 0 if a == b else 1, 'rebuilt': a, 'shipped': b}], {'blob_sha256': b})


def stage_policy(run):
    tmp = run.out / 'tables'
    tmp.mkdir(exist_ok=True)
    e1 = run.step('processes-report', ['python3', '-I', str(R / 'tools/report_processes.py'), str(R / 'data/tables/chemistry-v3.2-manifest.json'), str(tmp / 'PROCESSES.txt')], 'host')
    e2 = run.step('coverage-policy', ['python3', '-I', str(R / 'tools/make_coverage_policy.py'), str(R / 'data/tables/chemistry-v3.2-manifest.json'), str(tmp / 'policy.json')], 'host')
    same1 = (tmp / 'PROCESSES.txt').read_bytes() == (R / 'data/tables/PROCESSES.txt').read_bytes()
    same2 = (tmp / 'policy.json').read_bytes() == (R / 'data/tables/chemistry-v3.2-coverage-policy.json').read_bytes()
    pol = json.load(open(R / 'data/tables/chemistry-v3.2-coverage-policy.json'))
    run.stage('policy', [e1, e2, {'name': 'PROCESSES.txt regenerated byte-identical', 'exit': 0 if same1 else 1},
                         {'name': 'coverage policy regenerated byte-identical', 'exit': 0 if same2 else 1}],
              {'policy_sha256': sha(R / 'data/tables/chemistry-v3.2-coverage-policy.json'), 'admitted_unbounded': len(pol['admitted_unbounded_regions']),
               'refusing_unbounded': len(pol['refusing_unbounded_regions']), 'bounded': len(pol['bounded_regions'])})


def stage_build(run):
    bdir = R / 'build/bin'
    bdir.mkdir(parents=True, exist_ok=True)
    es = [run.step('build-test_core', [str(R / 'tools/dev_build.sh'), 'normal', str(R / 'tests/test_core.cpp'), str(bdir / 'test_core')]),
          run.step('build-test_core_san', [str(R / 'tools/dev_build.sh'), 'sanitized', str(R / 'tests/test_core.cpp'), str(bdir / 'test_core_sanitized')]),
          run.step('build-swarm_mc', [str(R / 'tools/dev_build.sh'), 'normal', str(R / 'tests/swarm_mc.cpp'), str(bdir / 'swarm_mc')]),
          run.step('build-gpu_parity', [str(R / 'tools/dev_build_gpu.sh'), str(bdir / 'gpu_parity')])]
    run.stage('build', es, {'binaries': {p.name: sha(p) for p in sorted(bdir.iterdir()) if p.is_file()},
                            'PF_BUILD_ISA_FLAGS': os.environ.get('PF_BUILD_ISA_FLAGS', '')})


def stage_core(run):
    blob, chk = R / 'data/tables/chemistry-v3.2.blob', R / 'data/tables/chemistry-v3.2-checks.json'
    e = run.step('core-normal', ['/usr/bin/time', '-v', str(R / 'build/bin/test_core'), str(blob), str(chk), str(run.out / 'core-normal'), '1000000'], timeout=10800)
    run.stage('core', [e], {'blob_sha256': sha(blob)})


def stage_sanitized(run):
    blob, chk = R / 'data/tables/chemistry-v3.2.blob', R / 'data/tables/chemistry-v3.2-checks.json'
    e = run.step('core-sanitized', ['env', 'ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:abort_on_error=1', 'UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1',
                                    str(R / 'build/bin/test_core_sanitized'), str(blob), str(chk), str(run.out / 'core-sanitized'), '100000'], timeout=14400)
    run.stage('sanitized', [e])


def stage_gpu(run):
    blob = R / 'data/tables/chemistry-v3.2.blob'
    es = [run.step('gpu-parity-process%d' % k, [str(R / 'build/bin/gpu_parity'), str(blob), str(run.out / ('gpu-process%d' % k)), str(NVRTC), GPU_UUID, '25000', '2'],
                   'gpu', timeout=3600) for k in (1, 2)]
    a = json.loads((run.out / 'gpu-process1/gpu-parity.json').read_text())
    b = json.loads((run.out / 'gpu-process2/gpu-parity.json').read_text())
    same = [c.get('record_digest_gpu') for c in a['cases']] == [c.get('record_digest_gpu') for c in b['cases']]
    run.stage('gpu', es + [{'name': 'two fresh processes give identical device digests', 'exit': 0 if same else 1}])


def stage_cost(run):
    # sampling cost per collision against v3.1 on the same GPU and CPU: v3.1 code + v3.1 tables (its sealed binary),
    # v3.2 code + v3.1 tables, and v3.2 code + v3.2 tables (the gpu stage); cpu_ms / gpu_kernel_ms per real event
    v31_blob = V31 / 'data/tables/chemistry-v3.1.blob'
    es = [run.step('cost-v31code-v31tables', [str(V31 / 'build/bin/gpu_parity'), str(v31_blob), str(run.out / 'cost-v31code-v31tables'), str(NVRTC), GPU_UUID, '25000', '2'],
                   'gpu', timeout=3600),
          run.step('cost-v32code-v31tables', [str(R / 'build/bin/gpu_parity'), str(v31_blob), str(run.out / 'cost-v32code-v31tables'), str(NVRTC), GPU_UUID, '25000', '2'],
                   'gpu', timeout=3600)]
    summ = {}
    for name, d in (('v3.1 code, v3.1 tables', 'cost-v31code-v31tables'), ('v3.2 code, v3.1 tables', 'cost-v32code-v31tables'),
                    ('v3.2 code, v3.2 tables', 'gpu-process1')):
        f = run.out / d / 'gpu-parity.json'
        if not f.exists():
            continue
        cs = [c for c in json.loads(f.read_text())['cases'] if c.get('case') == 'step' and not c['split_launch'] and c['block'] == 128 and c['repeat'] == 1]
        rows = [{'seed': c['seed'], 'events': c['events'], 'cpu_ms': c['cpu_ms'], 'gpu_kernel_ms': c['gpu_kernel_ms'],
                 'cpu_us_per_event': 1e3 * c['cpu_ms'] / c['events'], 'gpu_ns_per_event': 1e6 * c['gpu_kernel_ms'] / c['events']} for c in cs]
        ev = sum(r['events'] for r in rows)
        summ[name] = {'cases': rows, 'events': ev, 'cpu_us_per_event': 1e3 * sum(r['cpu_ms'] for r in rows) / ev,
                      'gpu_ns_per_event': 1e6 * sum(r['gpu_kernel_ms'] for r in rows) / ev}
    (run.out / 'cost-per-collision.json').write_text(json.dumps(summ, indent=1))
    run.stage('cost', es, {'cost_per_collision': {k: {'cpu_us_per_event': v['cpu_us_per_event'], 'gpu_ns_per_event': v['gpu_ns_per_event'],
                                                      'events': v['events']} for k, v in summ.items()}})


def stage_bench(run):
    # CPU sampling cost per collision on identical inputs: v3.1 code (v3.1 headers, read in place) and v3.2 code, the
    # same three populations as tests/gpu_parity.cpp; minimum over interleaved rounds (tools/bench_cost.sh)
    bdir, out = R / 'build/bin', run.out / 'bench'
    out.mkdir(exist_ok=True)
    arr = '/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include/pf_arrival'
    flags = ['-std=c++20', '-O2', '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-ffp-contract=off', '-fno-fast-math'] + os.environ.get('PF_BUILD_ISA_FLAGS', '').split()
    es = [run.step('build-bench_gen', ['/usr/bin/c++', *flags, '-DBENCH_GEN', '-I', str(R / 'include/pf_chem_v3'), '-I', arr, str(R / 'tests/bench_step.cpp'), '-lcrypto',
                                       '-o', str(bdir / 'bench_gen')]),
          run.step('build-bench_v32', ['/usr/bin/c++', *flags, '-I', str(R / 'include/pf_chem_v3'), '-I', arr, str(R / 'tests/bench_step.cpp'), '-lcrypto',
                                       '-o', str(bdir / 'bench_v32')]),
          run.step('build-bench_v31code', ['/usr/bin/c++', *flags, '-I', str(V31 / 'include/pf_chem_v3'), '-I', arr, str(R / 'tests/bench_step.cpp'), '-lcrypto',
                                           '-o', str(bdir / 'bench_v31code')])]
    for seed, cells in ((1000, 97), (4711, 31), (99, 7)):
        es.append(run.step('bench-population-%d' % seed, [str(bdir / 'bench_gen'), 'gen', str(R / 'data/tables/chemistry-v3.2.blob'), '25000', str(cells), str(seed),
                                                          str(out / ('pop-%d.bin' % seed))]))
    log = out / 'bench-cost.txt'
    t0 = time.time()
    with open(log, 'w') as f:
        p = subprocess.run([str(R / 'tools/lease_any.sh'), 'sh', str(R / 'tools/bench_cost.sh'), str(out), '6'], stdout=f, stderr=subprocess.STDOUT, cwd=str(R))
    es.append({'name': 'bench-cost', 'exit': p.returncode, 'seconds': round(time.time() - t0, 2), 'log': str(log)})
    print('%-28s exit %d  %.1f s' % ('bench-cost', p.returncode, es[-1]['seconds']), flush=True)
    run.stage('bench', es)


def stage_swarm(run):
    sdir = run.out / 'swarm'
    sdir.mkdir(exist_ok=True)
    # isotropic-elastic variant of the production set (angular-treatment comparison)
    spec = json.load(open(R / 'data/tables/chemistry-v3.2-spec.json'))
    for p in spec['groups'][0]['processes']:
        if p['kind'] == 'ELASTIC':
            p['angular'] = 'ISO'
            p.pop('sigma_is_momentum_transfer', None)
            p['id'] += ' [isotropic variant]'
    spec['id'] += '-isotropic-elastic-variant'
    json.dump(spec, open(sdir / 'v3-iso-spec.json', 'w'), indent=1)
    e0 = run.step('swarm-iso-tables', ['python3', '-I', str(R / 'tools/build_tables.py'), str(sdir / 'v3-iso-spec.json'), str(sdir / 'v3-iso')], 'host')
    sets = {'v3': str(R / 'data/tables/chemistry-v3.2.blob'), 'v3-iso': str(sdir / 'v3-iso.blob')}
    for s in ('muroranit', 'triniti', 'hayashi', 'phelps', 'morgan'):
        src = next(R.glob('data/sources/lxcat-%s-h2o-cs-*/Cross section.txt' % s))
        subprocess.run(['python3', '-I', str(R / 'tools/make_lxcat_spec.py'), str(src), s, str(sdir / ('lxcat-%s-spec.json' % s))], check=True, cwd=str(R))
        subprocess.run(['python3', '-I', str(R / 'tools/build_tables.py'), str(sdir / ('lxcat-%s-spec.json' % s)), str(sdir / ('lxcat-%s' % s))], check=True, cwd=str(R))
        sets['lxcat-' + s] = str(sdir / ('lxcat-%s.blob' % s))
    chunks = [EN_LIST[0:4], EN_LIST[4:8], EN_LIST[8:12]]
    jobs = []
    for name, blob in sets.items():
        for k, ch in enumerate(chunks):
            jobs.append({'name': '%s-%d' % (name, k), 'cwd': str(R), 'log': str(sdir / ('%s-%d.log' % (name, k))), 'output': str(sdir / ('%s-%d.json' % (name, k))),
                         'argv': [str(R / 'build/bin/swarm_mc'), blob, str(sdir / ('%s-%d.json' % (name, k))), ','.join(map(str, ch)), '3.22e22', '300', '2000', '20261031']})
    json.dump(jobs, open(sdir / 'jobs.json', 'w'), indent=1)
    e1 = run.step('swarm-jobs', ['python3', '-I', str(R / 'tools/run_jobs.py'), str(sdir / 'jobs.json'), str(sdir / 'jobs-receipt.json')], 'host', timeout=21600)
    e2 = run.step('swarm-table', ['python3', '-I', str(R / 'tools/swarm_table.py'), str(sdir), str(run.out / 'SWARM-VALIDATION')], 'host')
    # v3.2: the electron set changed (rotational J=2/3 tails, ion-pair tails, attachment above the resonances). Re-validation
    # (decision of 2026-10-31): six seeds of v3 (v3.1 tables) and v3.2, the Monte Carlo run directly at every measured E/N;
    # each v3.2 mean within the v3 multi-seed mean +- max(2 SD, 1 percentage point)
    mdir = run.out / 'swarm-multiseed'
    mdir.mkdir(exist_ok=True)
    seeds = ','.join(str(x) for x in range(20261101, 20261107))
    ms = R / 'tools/swarm_multiseed.py'
    subprocess.run(['python3', '-I', str(ms), 'jobs', 'v3', str(V31 / 'data/tables/chemistry-v3.1.blob'), str(mdir), seeds, str(mdir / 'jobs-v3.json')], check=True, cwd=str(R))
    subprocess.run(['python3', '-I', str(ms), 'jobs', 'v32', str(R / 'data/tables/chemistry-v3.2.blob'), str(mdir), seeds, str(mdir / 'jobs-v32.json')], check=True, cwd=str(R))
    allj = json.load(open(mdir / 'jobs-v3.json')) + json.load(open(mdir / 'jobs-v32.json'))
    allj = [x for pair in zip(allj[:6], allj[6:]) for x in pair]   # interleave the sets (partial runs stay balanced)
    json.dump(allj, open(mdir / 'jobs.json', 'w'), indent=1)
    e3 = run.step('swarm-multiseed-jobs', ['python3', '-I', str(R / 'tools/run_jobs.py'), str(mdir / 'jobs.json'), str(mdir / 'jobs-receipt.json')], 'host', timeout=43200)
    e4 = run.step('swarm-multiseed-summary-v3', ['python3', '-I', str(ms), 'summary', str(mdir), 'v3', str(mdir / 'summary-v3.json')], 'host')
    e5 = run.step('swarm-multiseed-summary-v32', ['python3', '-I', str(ms), 'summary', str(mdir), 'v32', str(mdir / 'summary-v32.json')], 'host')
    e6 = run.step('swarm-revalidation', ['python3', '-I', str(ms), 'compare', str(mdir / 'summary-v32.json'), str(mdir / 'summary-v3.json'), str(V3_SWARM),
                                         str(run.out / 'SWARM-REVALIDATION')], 'host')
    proof = json.load(open(run.out / 'SWARM-REVALIDATION.json'))
    run.stage('swarm', [e0, e1, e2, e3, e4, e5, e6], {'sets': {k: sha(v) for k, v in sets.items()},
              'revalidation_sha256': sha(run.out / 'SWARM-REVALIDATION.json'),
              'job_plan_sha256': {d: sha(run.out / d / 'jobs.json') for d in ('swarm','swarm-multiseed')},
              'scientific_criterion_status': proof['scientific_criterion_status']})


def main():
    label = sys.argv[1]
    stages = sys.argv[2:] or ['inputs', 'tables', 'policy', 'build', 'core', 'sanitized', 'gpu', 'cost', 'bench', 'swarm']
    run = Run(label)
    fns = {'inputs': stage_inputs, 'tables': stage_tables, 'policy': stage_policy, 'build': stage_build, 'core': stage_core, 'sanitized': stage_sanitized,
           'gpu': stage_gpu, 'cost': stage_cost, 'bench': stage_bench, 'swarm': stage_swarm}
    for s in stages:
        run.begin_stage(s)
        fns[s](run)
    ok = all(st['passed'] for st in run.receipt['stages'])
    print('RECEIPT', run.rp, 'passed' if ok else 'FAILED')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
