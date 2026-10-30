#!/usr/bin/env python3
"""Build and run the source-v4 audit reproducer port (normal and ASan/UBSan).

Every compile and run step is wrapped as
  flock <testing-core7.lock> taskset -c 7 nice -n 15 env <one-thread env> run-frozen <argv>
from a working directory under /home/treefor/projects. Builds are serial (one
compiler process at a time). Start this script itself at nice 0 under
`taskset -c 7` so the per-step `nice -n 15` yields exactly 15.

Writes a new numbered attempt directory under OUT and a receipt.json there.
"""
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

SRC = Path('/home/treefor/projects/runtimes/audit-reproducer-ports-20261030/source-v4')
OUT = Path('/home/treefor/sim-workspace/powerflow/outputs/audit-reproducer-ports-20261030/source-v4')
LOCK = '/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core7.lock'
FROZEN = '/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
FROZEN_ROOT = '/home/treefor/projects/runtimes/arch-base-20261028/root'
WT = Path('/home/treefor/projects/runtimes/pic-electrode-heat-feedback-20261029')
V4 = WT / 'native/mfem_electrode_heat_source_v4'
V3 = WT / 'native/mfem_electrode_heat_source_v3'
V2 = WT / 'native/mfem_electrode_heat_feedback_v2'
BUILD = WT / 'build-heat-source-v4'
COUNT_HEADER = '/home/treefor/projects/runtimes/pic-water-exact-count-subflow-20261029/native/source/include/raikou/exact_particle_count.h'
MFEM = Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
AUDIT = Path('/home/treefor/projects/Power-Flow/artifacts/pre-integration-audit-20261030')
AUDIT_SRC = AUDIT / 'reviews/source'
MELT = Path('/home/treefor/projects/runtimes/ss304-melt-ablation-proposal-v1-20261030/source')
MELT_IDENTITY = '908fec5449356342c2f1ba9e6192e530b31ec4e482a2930e3642cc9b2d4c2e73'
V4_HEAD = 'c65ac51874f069b4258e454537530eaff7755e66'

# Exactly the v4 build-tree flags (CMakeFiles/*/flags.make and link.txt).
DEFINE = f'-DHEAT_SOURCE_EXACT_COUNT_HEADER="{COUNT_HEADER}"'
INCLUDES = [f'-I{V3}/hall-owner/include', f'-I{V3}/hall-owner/src', f'-I{V4}', f'-I{V3}', f'-I{V2}',
            f'-I{V2}/source-kernel/include', '-isystem', f'{MFEM}/install-mpi/include',
            '-isystem', f'{MFEM}/tpl/hypre/include', '-isystem', f'{MFEM}/tpl/metis/include']
NORMAL = ['-O3', '-DNDEBUG', '-O2', '-ffp-contract=off', '-fno-fast-math']
SANITIZED = ['-O3', '-DNDEBUG', '-O1', '-g0', '-fno-omit-frame-pointer', '-fsanitize=address,undefined',
             '-fno-sanitize-recover=all', '-ffp-contract=off', '-fno-fast-math']
LIBS = ['-Wl,-rpath', '-Wl,/usr/lib', '-Wl,--enable-new-dtags', str(BUILD / 'libsource-v4-native-support.a'),
        '/usr/lib/libcrypto.so', f'{MFEM}/install-mpi/lib/libmfem.a', f'{MFEM}/tpl/hypre/lib/libHYPRE.a',
        f'{MFEM}/tpl/metis/lib/libmetis.a', f'{MFEM}/tpl/GKlib/lib/libGKlib.a', '-Wl,-rpath', '-Wl,/usr/lib',
        '-Wl,--enable-new-dtags', '/usr/lib/libmpi.so', str(BUILD / 'v3-owner/v2-owner/libheat-boundary-source-kernel.a')]
CXX = '/usr/sbin/c++'
ONE_THREAD = ['OMP_NUM_THREADS=1', 'OPENBLAS_NUM_THREADS=1', 'MKL_NUM_THREADS=1', 'CUDA_VISIBLE_DEVICES=']
SAN_ENV = ['ASAN_OPTIONS=detect_leaks=1:halt_on_error=1', 'UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1']
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
    """Path as seen on the host for a path used inside the frozen sandbox."""
    if p.startswith('/home/treefor/projects') or p.startswith('/home/treefor/sim-workspace'):
        return p
    return FROZEN_ROOT + p


class Runner:
    def __init__(self, attempt):
        self.A = attempt
        self.steps = []

    def step(self, name, argv, extra_env=(), timeout=1800, expect=0):
        cmd = ['flock', LOCK, 'taskset', '-c', '7', 'nice', '-n', '15', 'env', *ONE_THREAD, *extra_env, FROZEN, *argv]
        env = dict(os.environ)
        env['PWD'] = str(SRC)  # run-frozen chdirs to $PWD; keep it under /home/treefor/projects
        started_wall = time.time()
        t0 = time.monotonic()
        try:
            r = subprocess.run(cmd, cwd=SRC, env=env, capture_output=True, text=True, timeout=timeout)
            code, out, err = r.returncode, r.stdout, r.stderr
        except subprocess.TimeoutExpired as e:
            code, out, err = 'timeout', e.stdout or '', e.stderr or ''
            out = out.decode() if isinstance(out, bytes) else out
            err = err.decode() if isinstance(err, bytes) else err
        elapsed = time.monotonic() - t0
        log = self.A / 'logs' / f'{name}.log'
        log.write_text(f'$ {" ".join(cmd)}\n--- stdout ---\n{out}\n--- stderr ---\n{err}')
        rec = {'name': name, 'argv': cmd, 'cwd': str(SRC), 'exit_code': code, 'expected_exit_code': expect,
               'elapsed_s': round(elapsed, 3), 'started_unix': round(started_wall, 3), 'log': str(log)}
        self.steps.append(rec)
        print(f'{name}: exit {code} ({elapsed:.1f} s)', flush=True)
        return code, out, err


def deps(path):
    """Header/source paths from a make-style dependency file."""
    text = Path(path).read_text().replace('\\\n', ' ')
    files = []
    for line in text.splitlines():
        if ':' not in line:
            continue
        for p in line.split(':', 1)[1].split():
            if p not in files:
                files.append(p)
    return files


def main():
    assert os.sched_getaffinity(0) == {7}, 'start under taskset -c 7'
    assert os.getpriority(os.PRIO_PROCESS, 0) == 0, 'start at nice 0; each step adds nice 15'
    head = subprocess.run(['git', '-C', str(WT), 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(['git', '-C', str(WT), 'status', '--porcelain', '--', 'native'], capture_output=True, text=True).stdout
    assert head == V4_HEAD and dirty == '', (head, dirty)
    assert sha(COUNT_HEADER) == '46e01fe22e67cc9deb0da9320edef4041f4b3ef433de8ce1eca484463b405501'
    melt_files = {str(p.relative_to(MELT)): sha(p) for p in sorted(MELT.rglob('*')) if p.is_file()}
    melt_identity = hashlib.sha256(json.dumps(melt_files, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    assert melt_identity == MELT_IDENTITY, melt_identity

    OUT.mkdir(parents=True, exist_ok=True)
    n = 1
    while (OUT / f'attempt-{n:03d}').exists():
        n += 1
    A = OUT / f'attempt-{n:03d}'
    A.mkdir()
    for d in ('logs', 'build', 'normal', 'sanitized', 'melt-normal', 'melt-ubsan'):
        (A / d).mkdir()
    R = Runner(A)
    B = A / 'build'

    # ---- builds (serial)
    port = str(SRC / 'source_v4_port.cpp')
    R.step('build-port-normal', [CXX, *NORMAL, DEFINE, *INCLUDES, '-MD', '-MF', str(B / 'port-normal.d'),
                                 port, '-o', str(B / 'source_v4_port'), *LIBS])
    R.step('build-port-sanitized', [CXX, *SANITIZED, DEFINE, *INCLUDES, '-MD', '-MF', str(B / 'port-sanitized.d'),
                                    port, '-o', str(B / 'source_v4_port_sanitized'), *LIBS])
    mini = str(SRC / 'min_repro_a08_gaps.cpp')
    R.step('build-minrepro-normal', [CXX, *NORMAL, DEFINE, *INCLUDES, '-MD', '-MF', str(B / 'minrepro-normal.d'),
                                     mini, '-o', str(B / 'min_repro_a08_gaps'), *LIBS])
    R.step('build-minrepro-sanitized', [CXX, *SANITIZED, DEFINE, *INCLUDES, mini, '-o', str(B / 'min_repro_a08_gaps_sanitized'), *LIBS])
    reject = {}
    for case in range(0, 12):
        code, out, err = R.step(f'compile-reject-{case:02d}', [CXX, *NORMAL, DEFINE, *INCLUDES, '-fsyntax-only',
                                f'-DREJECT_CASE={case}', '-MD', '-MF', str(B / f'reject-{case:02d}.d'),
                                str(SRC / 'compile_reject.cpp')], expect=0 if case == 0 else 'nonzero')
        first = next((l for l in err.splitlines() if ': error:' in l), '')
        if case == 0:
            ok = code == 0
        else:
            ok = code not in (0, 'timeout') and 'compile_reject.cpp' in first.split(': error:')[0] \
                and re.search(REJECT[case], err) is not None
        reject[case] = {'exit_code': code, 'first_error': first, 'expected_pattern': REJECT.get(case), 'as_expected': ok}
    melt_flags = ['/usr/bin/g++', '-std=c++20', '-O1', '-fno-fast-math', '-ffp-contract=off', '-Wall', '-Wextra',
                  '-Wpedantic', '-Werror', f'-I{MELT}/vendor', f'-I{MELT}', f'-DABLATION_PROPOSAL_SOURCE_SHA256="{MELT_IDENTITY}"']
    R.step('build-melt-normal', [*melt_flags, '-MD', '-MF', str(B / 'melt-normal.d'), str(SRC / 'melt_conflict.cpp'),
                                 '-lcrypto', '-o', str(B / 'melt_conflict')])
    R.step('build-melt-ubsan', [*melt_flags, '-fsanitize=undefined', '-fno-sanitize-recover=all', str(SRC / 'melt_conflict.cpp'),
                                '-lcrypto', '-o', str(B / 'melt_conflict_ubsan')])

    # ---- runs
    if (B / 'source_v4_port').exists():
        R.step('run-port-normal', [str(B / 'source_v4_port'), str(A / 'normal')], timeout=900)
    if (B / 'source_v4_port_sanitized').exists():
        R.step('run-port-sanitized', [str(B / 'source_v4_port_sanitized'), str(A / 'sanitized')], extra_env=SAN_ENV, timeout=1800)
    minirun = {}
    for tag, binary, extra_env in (('normal', 'min_repro_a08_gaps', ()), ('sanitized', 'min_repro_a08_gaps_sanitized', SAN_ENV)):
        if (B / binary).exists():
            code, out, err = R.step(f'run-minrepro-{tag}', [str(B / binary)], extra_env=extra_env, timeout=900)
            minirun[tag] = {'exit_code': code, 'stdout': out.splitlines()}
    counts = A / 'normal' / 'exact-counts.json'
    if counts.exists():
        R.step('run-fraction-oracle', ['/usr/bin/python3', '-I', str(SRC / 'fraction_oracle.py'), str(counts),
                                       str(A / 'normal' / 'port-results.json'), str(AUDIT_SRC), str(A / 'fraction-oracle.json')], timeout=300)
        accepted = json.loads(counts.read_text()).get('c1_accepted', [])
        (A / 'v4-accepted-stock-carriers.txt').write_text(''.join(h + '\n' for h in accepted))
        if accepted and (B / 'melt_conflict').exists():
            R.step('run-melt-normal', [str(B / 'melt_conflict'), str(A / 'v4-accepted-stock-carriers.txt'), str(A / 'melt-normal')], timeout=300)
        if accepted and (B / 'melt_conflict_ubsan').exists():
            R.step('run-melt-ubsan', [str(B / 'melt_conflict_ubsan'), str(A / 'v4-accepted-stock-carriers.txt'), str(A / 'melt-ubsan')],
                   extra_env=['UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1'], timeout=300)

    # ---- inputs (every compiled header from the dependency files, libraries, launcher, audit records)
    inputs = {}
    for d in sorted(B.glob('*.d')):
        for p in deps(d):
            hp = host_path(p)
            if os.path.isfile(hp):
                inputs[p] = {'sha256': sha(hp), 'hashed_at': hp}
            else:
                inputs[p] = {'sha256': None, 'hashed_at': hp, 'note': 'missing'}
    extra = [*(str(p) for p in sorted(SRC.iterdir()) if p.is_file()),
             *(str(p) for p in sorted(V4.rglob('*')) if p.is_file()),
             str(BUILD / 'libsource-v4-native-support.a'), str(BUILD / 'v3-owner/v2-owner/libheat-boundary-source-kernel.a'),
             f'{MFEM}/install-mpi/lib/libmfem.a', f'{MFEM}/tpl/hypre/lib/libHYPRE.a', f'{MFEM}/tpl/metis/lib/libmetis.a',
             f'{MFEM}/tpl/GKlib/lib/libGKlib.a', '/usr/lib/libmpi.so', '/usr/lib/libcrypto.so', CXX, '/usr/bin/g++', '/usr/bin/python3',
             FROZEN, f'{FROZEN_ROOT}/../kernel.release',
             str(AUDIT / 'AUDIT.txt'), *(str(AUDIT_SRC / x) for x in ('source_v3_reproducer.cpp', 'reproducer.log',
             'exact-allocation-evidence.json', 'build_command.json', 'REPRODUCE.txt', 'REPORT.txt', 'findings.json')),
             *(str(MELT / f) for f in melt_files)]
    for p in extra:
        hp = host_path(os.path.realpath(p)) if not p.startswith(FROZEN_ROOT) else p
        if p in ('/usr/lib/libmpi.so', '/usr/lib/libcrypto.so', CXX, '/usr/bin/g++', '/usr/bin/python3'):
            hp = os.path.realpath(FROZEN_ROOT + p)
        inputs[p] = {'sha256': sha(hp) if os.path.isfile(hp) else None, 'hashed_at': hp}

    # ---- results
    def load(p):
        try:
            return json.loads(Path(p).read_text())
        except Exception as e:  # recorded, not hidden
            return {'unreadable': str(e)}
    normal = load(A / 'normal' / 'port-results.json')
    parity = {}
    for f in ('port-results.json', 'exact-counts.json'):
        a, b = A / 'normal' / f, A / 'sanitized' / f
        parity[f] = (a.exists() and b.exists() and a.read_bytes() == b.read_bytes())
    receipt = {
        'schema': 'audit-reproducer-port-receipt-v1', 'component': 'electrode-heat-source-v4', 'v4_head': V4_HEAD,
        'reproducer': str(AUDIT_SRC / 'source_v3_reproducer.cpp'), 'attempt': str(A),
        'cpu': 7, 'nice': 15, 'threads': 1, 'lease': LOCK, 'frozen_launcher': FROZEN, 'build_parallelism': 1,
        'kernel': Path('/proc/sys/kernel/osrelease').read_text().strip(),
        'full_joint_test': False, 'v4_sources_or_build_tree_modified': False,
        'melt_source_identity': melt_identity,
        'steps': R.steps, 'compile_reject': reject,
        'compile_reject_all_as_expected': all(v['as_expected'] for v in reject.values()),
        'normal_sanitized_byte_parity': parity,
        'per_counterexample': [{k: p.get(k) for k in ('id', 'audit_finding', 'audit_input', 'v3_behavior', 'v4_verdict')}
                               for p in normal.get('probes', [])],
        'port_results_normal': normal,
        'port_results_sanitized_status': 'see sanitized/port-results.json',
        'fraction_oracle': load(A / 'fraction-oracle.json'),
        'minimal_reproducer_a08_gaps': minirun,
        'melt_conflict_normal': load(A / 'melt-normal' / 'melt-conflict.json'),
        'melt_conflict_ubsan': load(A / 'melt-ubsan' / 'melt-conflict.json'),
        'inputs_sha256': inputs,
    }
    (A / 'receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    print('receipt:', A / 'receipt.json')
    bad = [s['name'] for s in R.steps if not (s['exit_code'] == 0 or (s['expected_exit_code'] == 'nonzero' and s['exit_code'] not in (0, 'timeout')))]
    print('unexpected step exits:', bad or 'none')
    return 0 if not bad else 1


if __name__ == '__main__':
    sys.exit(main())
