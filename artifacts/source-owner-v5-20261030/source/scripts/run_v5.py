#!/usr/bin/env python3
"""Configure, build and run heat/source owner v5 qualification; write receipts.

Builds run under  core-lease-20261030/both-cores.sh run-frozen <argv>   (cores 6,7; -j2; nice 15).
Runs   run under  core-lease-20261030/any-core.sh --wait env CUDA_VISIBLE_DEVICES= run-frozen <argv>
                  (core 6 or 7, nice 15, one thread, no GPU).
Every step uses a working directory under /home/treefor/projects (run-frozen chdirs to $PWD).
Start this script at nice 0 restricted to cores 6-7:  taskset -c 6,7 python3 run_v5.py [stage ...]
Stages: build run oracle (default: all, in that order). build = configure + build + compile-reject.
Each invocation writes a new OUT/attempt-NNN directory with logs, outputs and receipt.json.
"""
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
import time
from pathlib import Path

V5 = Path('/home/treefor/projects/runtimes/pic-electrode-heat-source-v5-20261030')
MOD = V5 / 'native/mfem_electrode_heat_source_v5'
PORT = V5 / 'audit-port'
BUILD = V5 / 'build'
# Workspace alias; resolved to the physical directory, which is the one bound inside the frozen sandbox.
WORKSPACE = Path('/home/treefor/sim-workspace/powerflow').resolve()
OUT = WORKSPACE / 'outputs/pic-electrode-heat-source-v5-20261030'
LEASE = Path('/home/treefor/projects/runtimes/core-lease-20261030')
FROZEN = '/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
FROZEN_ROOT = '/home/treefor/projects/runtimes/arch-base-20261028/root'
WT = Path('/home/treefor/projects/runtimes/pic-electrode-heat-feedback-20261029')
V4_HEAD = 'c65ac51874f069b4258e454537530eaff7755e66'
COUNT_HEADER = '/home/treefor/projects/runtimes/pic-water-exact-count-subflow-20261029/native/source/include/raikou/exact_particle_count.h'
COUNT_SHA = '46e01fe22e67cc9deb0da9320edef4041f4b3ef433de8ce1eca484463b405501'
MFEM_DIR = '/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/cmake/mfem'
AUDIT_SRC = Path('/home/treefor/projects/Power-Flow/artifacts/pre-integration-audit-20261030/reviews/source')
V4_PORT = Path('/home/treefor/projects/Power-Flow/artifacts/audit-reproducer-ports-20261030/source-v4')
SAN_ENV = ['ASAN_OPTIONS=detect_leaks=1:halt_on_error=1', 'UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1']
TARGETS = ['electrode-heat-source-v5-tests', 'electrode-heat-source-v5-sanitized', 'source_v5_port_normal',
           'source_v5_port_sanitized', 'min_repro_a08_gaps_v5_normal', 'min_repro_a08_gaps_v5_sanitized']
REJECT = {1: r'read-only', 2: r'private within this context', 3: r'no matching function|private within this context',
          4: r'no matching function|cannot convert', 5: r"has no member named .native.", 6: r'private within this context',
          7: r'no matching function|private within this context', 8: r'private within this context',
          9: r'use of deleted function', 10: r'use of deleted function', 11: r'read-only'}


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def host_path(p):
    if p.startswith('/home/treefor/projects') or p.startswith(str(WORKSPACE.parent.parent)):
        return p
    return FROZEN_ROOT + p


def tree_hash(roots):
    files = sorted(str(p) for r in roots for p in Path(r).rglob('*') if p.is_file())
    h = hashlib.sha256()
    for f in files:
        h.update(f'{sha(f)}  {f}\n'.encode())
    return h.hexdigest(), len(files)


class Runner:
    def __init__(self, attempt):
        self.A = attempt
        self.steps = []

    def step(self, name, argv, kind, extra_env=(), timeout=3600, expect=0):
        if kind == 'build':
            cmd = [str(LEASE / 'both-cores.sh'), FROZEN, *argv]
        else:
            cmd = [str(LEASE / 'any-core.sh'), '--wait', 'env', 'CUDA_VISIBLE_DEVICES=', *extra_env, FROZEN, *argv]
        env = dict(os.environ)
        env['PWD'] = str(V5)
        t0 = time.monotonic()
        wall = time.time()
        try:
            r = subprocess.run(cmd, cwd=V5, env=env, capture_output=True, text=True, timeout=timeout)
            code, out, err = r.returncode, r.stdout, r.stderr
        except subprocess.TimeoutExpired as e:
            code = 'timeout'
            out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or '')
            err = e.stderr.decode() if isinstance(e.stderr, bytes) else (e.stderr or '')
        elapsed = time.monotonic() - t0
        log = self.A / 'logs' / f'{name}.log'
        log.write_text(f'$ {shlex.join(cmd)}\n--- stdout ---\n{out}\n--- stderr ---\n{err}')
        cpu = re.search(r'on CPU(\d)|cpu=(\d)', out)
        rec = {'name': name, 'argv': cmd, 'cwd': str(V5), 'exit_code': code, 'expected_exit_code': expect,
               'elapsed_s': round(elapsed, 3), 'started_unix': round(wall, 3), 'log': str(log),
               'lease': 'both-cores (6,7)' if kind == 'build' else 'any-core (6 or 7)',
               'cpu_reported': next((g for g in cpu.groups() if g), None) if cpu else None}
        self.steps.append(rec)
        print(f'{name}: exit {code} ({elapsed:.1f} s)', flush=True)
        return code, out, err


def deps(path):
    text = Path(path).read_text().replace('\\\n', ' ')
    files = []
    for line in text.splitlines():
        if ':' not in line:
            continue
        for p in line.split(':', 1)[1].split():
            if p not in files:
                files.append(p)
    return files


def v4_state():
    head = subprocess.run(['git', '-C', str(WT), 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(['git', '-C', str(WT), 'status', '--porcelain', '--', 'native'], capture_output=True, text=True).stdout
    th, n = tree_hash([WT / 'native', WT / 'build-heat-source-v4'])
    ph, pn = tree_hash([V4_PORT])
    # Same pipeline (and locale) as v4-baseline-before.txt.
    shell = subprocess.run(f'find {WT}/native {WT}/build-heat-source-v4 -type f -print0 | sort -z | xargs -0 sha256sum | sha256sum',
                           shell=True, capture_output=True, text=True).stdout.split()[0]
    return {'head': head, 'baseline_pipeline_tree_sha256': shell, 'native_status_porcelain': dirty, 'native_and_build_v4_tree_sha256': th, 'files': n,
            'source_v4_port_artifacts_tree_sha256': ph, 'source_v4_port_files': pn}


def main(stages):
    assert os.sched_getaffinity(0) <= {6, 7}, 'start under taskset -c 6,7'
    assert os.getpriority(os.PRIO_PROCESS, 0) == 0, 'start at nice 0; each leased step adds nice 15'
    before = v4_state()
    assert before['head'] == V4_HEAD and before['native_status_porcelain'] == '', before
    assert sha(COUNT_HEADER) == COUNT_SHA
    OUT.mkdir(parents=True, exist_ok=True)
    n = 1
    while (OUT / f'attempt-{n:03d}').exists():
        n += 1
    A = OUT / f'attempt-{n:03d}'
    for d in ('logs',):
        (A / d).mkdir(parents=True)
    R = Runner(A)
    res = {}

    if 'build' in stages:
        # One both-cores lease: configure, build all targets, compile-reject cases.
        R.step('configure-build-reject', ['/bin/bash', str(V5 / 'scripts/build.sh'), str(A)], 'build', timeout=10800)
        reject = {}
        for case in range(12):
            logp, exitp = A / 'reject' / f'case-{case:02d}.log', A / 'reject' / f'case-{case:02d}.exit'
            if not (logp.exists() and exitp.exists()):
                continue
            log = logp.read_text()
            code = int(exitp.read_text().strip())
            first = next((l for l in log.splitlines() if ': error:' in l), '')
            ok = code == 0 if case == 0 else (code != 0 and 'compile_reject.cpp' in first.split(': error:')[0] and re.search(REJECT[case], log) is not None)
            reject[case] = {'exit_code': code, 'first_error': first, 'expected_pattern': REJECT.get(case), 'as_expected': ok}
        res['compile_reject'] = reject
        res['compile_reject_all_as_expected'] = len(reject) == 12 and all(v['as_expected'] for v in reject.values())

    if 'run' in stages:
        tests = BUILD / 'native/mfem_electrode_heat_source_v5'
        R.step('tests-normal', [str(tests / 'electrode-heat-source-v5-tests'), str(A / 'checks-normal')], 'run')
        R.step('tests-sanitized', [str(tests / 'electrode-heat-source-v5-sanitized'), str(A / 'checks-sanitized')], 'run', SAN_ENV)
        R.step('port-normal', [str(BUILD / 'audit-port/source_v5_port_normal'), str(A / 'port-normal')], 'run')
        R.step('port-sanitized', [str(BUILD / 'audit-port/source_v5_port_sanitized'), str(A / 'port-sanitized')], 'run', SAN_ENV)
        for tag, binary, env in (('normal', 'min_repro_a08_gaps_v5_normal', ()), ('sanitized', 'min_repro_a08_gaps_v5_sanitized', SAN_ENV)):
            code, out, err = R.step(f'minrepro-{tag}', [str(BUILD / 'audit-port' / binary)], 'run', env)
            (A / f'minrepro-{tag}.txt').write_text(out)
            res[f'minrepro_{tag}'] = {'exit_code': code, 'stdout': out.splitlines()}

    if 'oracle' in stages:
        py = ['/usr/bin/python3', '-I']
        for tag in ('normal', 'sanitized'):
            code, out, err = R.step(f'count-oracle-{tag}', [*py, str(MOD / 'count_oracle.py'), str(A / f'checks-{tag}')], 'run')
            (A / f'count-oracle-{tag}.json').write_text(out)
            R.step(f'precision-oracle-{tag}', [*py, str(MOD / 'precision_oracle.py'), str(A / f'checks-{tag}/precision-dump.json'),
                                                str(A / f'precision-oracle-{tag}.json')], 'run')
            R.step(f'port-fraction-oracle-{tag}', [*py, str(PORT / 'fraction_oracle_v5.py'), str(A / f'port-{tag}/exact-counts.json'),
                                                    str(A / f'port-{tag}/port-results.json'), str(AUDIT_SRC), str(A / f'port-fraction-oracle-{tag}.json')], 'run')

    # ---- parity, inputs, results
    def load(p):
        try:
            return json.loads(Path(p).read_text())
        except Exception as e:
            return {'unreadable': str(e)}
    parity = {}
    for a, b in (('checks-normal', 'checks-sanitized'), ('port-normal', 'port-sanitized')):
        da, db = A / a, A / b
        if da.exists() and db.exists():
            names = sorted({p.name for p in da.iterdir()} | {p.name for p in db.iterdir()})
            parity[a] = {nm: ((da / nm).exists() and (db / nm).exists() and (da / nm).read_bytes() == (db / nm).read_bytes()) for nm in names}
    if (A / 'minrepro-normal.txt').exists() and (A / 'minrepro-sanitized.txt').exists():
        parity['minrepro-stdout'] = (A / 'minrepro-normal.txt').read_bytes() == (A / 'minrepro-sanitized.txt').read_bytes()
    inputs = {}
    for d in sorted(BUILD.rglob('*.o.d')):
        for p in deps(d):
            hp = host_path(p)
            inputs[p] = {'sha256': sha(hp) if os.path.isfile(hp) else None, 'hashed_at': hp}
    for p in [*sorted(MOD.iterdir()), *sorted(PORT.iterdir()), *sorted((V5 / 'scripts').iterdir()), V5 / 'CMakeLists.txt']:
        if p.is_file():
            inputs[str(p)] = {'sha256': sha(p), 'hashed_at': str(p)}
    for p in ['/usr/lib/libmpi.so', '/usr/lib/libcrypto.so', '/usr/sbin/c++', '/usr/bin/python3', '/usr/bin/cmake']:
        hp = os.path.realpath(FROZEN_ROOT + p)
        inputs[p] = {'sha256': sha(hp), 'hashed_at': hp}
    for p in [FROZEN, f'{FROZEN_ROOT}/../kernel.release', str(LEASE / 'any-core.sh'), str(LEASE / 'both-cores.sh'),
              *(str(BUILD / t) for t in ('native/mfem_electrode_heat_source_v5/libsource-v5-native-support.a',
                                         'native/mfem_electrode_heat_source_v5/v3-owner/v2-owner/libheat-boundary-source-kernel.a')),
              *(str(BUILD / 'native/mfem_electrode_heat_source_v5' / t) for t in ('electrode-heat-source-v5-tests', 'electrode-heat-source-v5-sanitized')),
              *(str(BUILD / 'audit-port' / t) for t in TARGETS[2:])]:
        hp = os.path.realpath(p)
        if os.path.isfile(hp):
            inputs[p] = {'sha256': sha(hp), 'hashed_at': hp}
    after = v4_state()
    baseline = (OUT / 'v4-baseline-before.txt').read_text()
    receipt = {
        'schema': 'source-owner-v5-qualification-receipt-v1', 'attempt': str(A), 'stages': stages,
        'kernel': Path('/proc/sys/kernel/osrelease').read_text().strip(), 'frozen_launcher': FROZEN,
        'frozen_launcher_sha256': sha(FROZEN), 'threads': 1, 'nice': 15, 'cores': '6-7 only', 'gpu': 'none (CUDA_VISIBLE_DEVICES empty)',
        'v4_state_before': before, 'v4_state_after': after,
        'v4_unchanged': before == after and after['baseline_pipeline_tree_sha256'] in baseline,
        'steps': R.steps, **res, 'normal_sanitized_byte_parity': parity,
        'tests_normal': load(A / 'checks-normal/results.json'), 'tests_sanitized': load(A / 'checks-sanitized/results.json'),
        'count_oracle_normal': load(A / 'count-oracle-normal.json'),
        'precision_oracle_normal': {k: v for k, v in load(A / 'precision-oracle-normal.json').items() if k != 'checks'},
        'precision_oracle_sanitized_summary': {k: load(A / 'precision-oracle-sanitized.json').get(k) for k in ('passed', 'failed')},
        'port_fraction_oracle_normal': load(A / 'port-fraction-oracle-normal.json'),
        'port_verdicts': [{k: p.get(k) for k in ('id', 'audit_finding', 'v5_verdict')} for p in load(A / 'port-normal/port-results.json').get('probes', [])],
        'inputs_sha256': inputs,
    }
    (A / 'receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    print('receipt:', A / 'receipt.json')
    bad = [s['name'] for s in R.steps if s['exit_code'] != 0]
    print('unexpected step exits:', bad or 'none')
    print('v4 unchanged:', receipt['v4_unchanged'])
    return 0 if not bad and receipt['v4_unchanged'] else 1


if __name__ == '__main__':
    stages = sys.argv[1:] or ['build', 'run', 'oracle']
    sys.exit(main(stages))
