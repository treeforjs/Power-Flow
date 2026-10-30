#!/usr/bin/env python3
"""Versioned CPU7 runner for the SS304 melt-ablation proposal v1 sources.

Successor to ss304-melt-ablation-proposal-v1-20261030/source/build_check.py,
which asserts the pre-reboot stopped production process and frozen cgroup.
This runner leaves the v1 sources untouched, computes the same source identity,
uses explicit checks that cannot be disabled by optimization, and guards the
paused production run by its accepted checkpoint bytes instead of a process.
"""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

V1 = Path('/home/treefor/projects/runtimes/ss304-melt-ablation-proposal-v1-20261030')
S = V1 / 'source'
O = Path('/home/treefor/sim-workspace/powerflow/outputs') / V1.name
FROZEN = '/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
PAUSED = Path('/home/treefor/sim-workspace/powerflow/outputs/full60-gpu-manual-pause-checkpoint/checkpoint_437462.bin')
PAUSED_SHA256 = '0a5a91ed93821e2487167b69e463dc1d56ee33a8df03707c96bb7bfc69c0ced3'
RUNNER_VERSION = 'ss304-melt-ablation-runner-v2'


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def production():
    running = []
    for p in Path('/proc').iterdir():
        if p.name.isdigit():
            try:
                if 'heating-history-pic' in (p / 'cmdline').read_bytes().decode(errors='replace'):
                    running.append(int(p.name))
            except OSError:
                pass
    r = dict(paused_checkpoint=str(PAUSED), paused_checkpoint_sha256=sha(PAUSED),
             production_processes=sorted(running))
    require(r['paused_checkpoint_sha256'] == PAUSED_SHA256, 'paused production checkpoint changed')
    require(not running, 'a production PIC process is running; component checks refused')
    return r


def source_files():
    return {str(p.relative_to(S)): sha(p) for p in sorted(S.rglob('*')) if p.is_file()}


def main():
    require(os.sched_getaffinity(0) == {7}, 'run with taskset -c 7')
    require(os.getpriority(os.PRIO_PROCESS, 0) == 15, 'run with nice -n 15')
    require(len(sys.argv) == 2 and sys.argv[1] in {'normal', 'ubsan'}, 'usage: build_check.py normal|ubsan')
    mode = sys.argv[1]
    lock = (O.parent / 'pic-profile/testing-core7.lock').open('a')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print('CPU7 lease busy; no computation started')
        return 75
    resource.setrlimit(resource.RLIMIT_AS, (3 * 1024**3, 3 * 1024**3))
    before = production()
    deps = json.loads((O / 'prerequisite-bindings.json').read_text())
    for original, pin in deps['inputs'].items():
        require(sha(Path(original)) == pin['sha256'] == sha(S / pin['copy']), f'prerequisite changed: {original}')
    files = source_files()
    identity = hashlib.sha256(json.dumps(files, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    attempt = 1
    while (O / f'{mode}-runner-v2-{attempt:03d}').exists():
        attempt += 1
    A = O / f'{mode}-runner-v2-{attempt:03d}'
    A.mkdir()
    binary = A / 'test_proposal'
    cmd = [FROZEN, '/usr/bin/g++', '-std=c++20', '-O1', '-fno-fast-math',
           '-ffp-contract=off', '-Wall', '-Wextra', '-Wpedantic', '-Werror',
           '-I' + str(S / 'vendor'),
           '-DABLATION_PROPOSAL_SOURCE_SHA256="' + identity + '"']
    if mode == 'ubsan':
        cmd += ['-fsanitize=undefined', '-fno-sanitize-recover=all']
    cmd += [str(S / 'test_proposal.cpp'), '-lcrypto', '-o', str(binary)]
    env = {'PATH': '/usr/bin:/bin', 'LANG': 'C', 'LC_ALL': 'C',
           'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
           'MKL_NUM_THREADS': '1', 'BLIS_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1',
           'UBSAN_OPTIONS': 'halt_on_error=1:print_stacktrace=1'}
    commands = {'compile': cmd, 'run': [FROZEN, str(binary), str(A)], 'environment': env}
    (A / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
    started = time.monotonic()
    receipt = dict(runner=RUNNER_VERSION, runner_sha256=sha(Path(__file__).resolve()),
                   mode=mode, source_identity=identity, source_files=files,
                   cpu=7, nice=15, exclusive_core7_lease=True, frozen=True,
                   threads=1, production_before=before, full_integration=False,
                   new_GPU_work=False, signals_sent=False)
    try:
        with (A / 'compile.log').open('w') as log:
            c = subprocess.run(cmd, cwd=V1, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=48)
        receipt['compile_returncode'] = c.returncode
        if c.returncode:
            print((A / 'compile.log').read_text()[-10000:])
            return c.returncode
        receipt['binary_sha256'] = sha(binary)
        with (A / 'run.log').open('w') as log:
            r = subprocess.run(commands['run'], cwd=V1, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=5)
        receipt['run_returncode'] = r.returncode
        print((A / 'run.log').read_text().strip())
        return r.returncode
    finally:
        receipt['elapsed_s'] = time.monotonic() - started
        receipt['production_after'] = production()
        receipt['sources_unchanged'] = files == source_files()
        (A / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        require(receipt['production_after'] == before, 'production state changed during checks')
        require(receipt['sources_unchanged'], 'v1 sources changed during checks')
        print(json.dumps({'attempt': str(A), 'elapsed_s': receipt['elapsed_s'],
                          'compile': receipt.get('compile_returncode'),
                          'run': receipt.get('run_returncode')}))


if __name__ == '__main__':
    raise SystemExit(main())
