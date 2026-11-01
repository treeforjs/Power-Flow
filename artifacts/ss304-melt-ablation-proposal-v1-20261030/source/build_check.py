#!/usr/bin/env python3
"""Bounded, exclusive CPU7 component checks; no integration or GPU execution."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

P = Path(__file__).resolve().parent.parent
S = P / 'source'
O = Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs') / P.name
FROZEN = '/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def production():
    p = Path('/proc/278902')
    fields = (p / 'stat').read_text().rsplit(')', 1)[1].split()
    events = dict(x.split() for x in Path('/sys/fs/cgroup/powerflow-pic-gpu-reserved/cgroup.events').read_text().splitlines())
    r = dict(pid=278902, state=fields[0], start_ticks=int(fields[19]),
             executable=os.readlink(p / 'exe'), executable_sha256=sha(p / 'exe'),
             cgroup_events=events)
    assert r['state'] == 'T' and r['start_ticks'] == 5288339
    assert events == {'populated': '1', 'frozen': '1'}
    assert r['executable_sha256'] == 'cddd9f1867b50f471b6e8ce49e698c45414b568acb6e42c5a4364d5c5fb37797'
    return r


def main():
    assert __debug__ and not os.getenv('PYTHONOPTIMIZE')
    assert os.sched_getaffinity(0) == {7} and os.getpriority(os.PRIO_PROCESS, 0) == 15
    mode = sys.argv[1]
    assert mode in {'normal', 'ubsan'}
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
        assert sha(Path(original)) == pin['sha256'] == sha(S / pin['copy'])
    files = {str(p.relative_to(S)): sha(p) for p in sorted(S.rglob('*')) if p.is_file()}
    identity = hashlib.sha256(json.dumps(files, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    attempt = 1
    while (O / f'{mode}-{attempt:03d}').exists():
        attempt += 1
    A = O / f'{mode}-{attempt:03d}'
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
    receipt = dict(mode=mode, source_identity=identity, source_files=files,
                   cpu=7, nice=15, exclusive_core7_lease=True, frozen=True,
                   threads=1, production_before=before, full_integration=False,
                   new_GPU_work=False, signals_sent=False)
    result = 1
    try:
        with (A / 'compile.log').open('w') as log:
            c = subprocess.run(cmd, cwd=P, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=48)
        receipt['compile_returncode'] = c.returncode
        if c.returncode:
            print((A / 'compile.log').read_text()[-10000:])
            return c.returncode
        receipt['binary_sha256'] = sha(binary)
        with (A / 'run.log').open('w') as log:
            r = subprocess.run(commands['run'], cwd=P, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=5)
        receipt['run_returncode'] = r.returncode
        print((A / 'run.log').read_text().strip())
        result = r.returncode
        return result
    finally:
        receipt['elapsed_s'] = time.monotonic() - started
        receipt['production_after'] = production()
        assert receipt['production_after'] == before
        assert files == {str(p.relative_to(S)): sha(p) for p in sorted(S.rglob('*')) if p.is_file()}
        (A / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        print(json.dumps({'attempt': str(A), 'elapsed_s': receipt['elapsed_s'],
                          'compile': receipt.get('compile_returncode'),
                          'run': receipt.get('run_returncode')}))


if __name__ == '__main__':
    raise SystemExit(main())
