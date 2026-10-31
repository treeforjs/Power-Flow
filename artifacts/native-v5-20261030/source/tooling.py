"""Shared build/run helpers for native v5. Every compiler and test process runs
through the frozen launcher with cwd under /home/treefor/projects, inside the
caller-held Power-Flow core lease (core-lease-20261030/both-cores.sh for builds,
any-core.sh --wait for single-threaded runs), at nice 15 with one thread and no
GPU. Normal processes run under RLIMIT_AS 3 GiB (prlimit inside the launcher);
ASan/UBSan processes cannot reserve their shadow mapping under RLIMIT_AS, so they
use the sanitizer's own hard_rss_limit_mb=3072 instead (recorded per command)."""
from pathlib import Path
import json, os, subprocess, hashlib, time, datetime

P = Path(__file__).resolve().parent
RT = P.parent
C003 = RT / 'native-full-time-candidate003'
EXACT = RT / 'native-exact-count'
OUT = Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/'
           'reaction-transactions-20261029/native-full-time-candidate005')
FROZEN = '/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
MFEM = Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
VENDOR = RT / 'source/vendor'
HALL = Path('/home/treefor/projects/runtimes/hall-coupled-transport-20261029/source/development')
HALL_DRIVER = HALL / 'hall-canonical-driver-20261029'
RLIMIT_AS = 3 * 1024 ** 3
LIBS = [str(MFEM / 'install-mpi/lib/libmfem.a'), str(MFEM / 'tpl/hypre/lib/libHYPRE.a'),
        str(MFEM / 'tpl/metis/lib/libmetis.a'), str(MFEM / 'tpl/GKlib/lib/libGKlib.a'), '-lcrypto']
SYSTEM = []
for d in ['install-mpi/include', 'tpl/hypre/include', 'tpl/metis/include']:
    SYSTEM += ['-isystem', str(MFEM / d)]
TRAPS = ['_ZN6raikou12mfem_backend6GapPic12solve_fieldsEv', '_ZN6raikou12mfem_backend6GapPic7advanceEd',
         '_ZN6raikou12mfem_backend6GapPic13advance_emptyEd',
         '_ZN6raikou12mfem_backend6GapPic16advance_observedEdRKSt8functionIFvRKNS0_11GapParticleEidEE']
ASAN = ['-fsanitize=address,undefined', '-fno-omit-frame-pointer']
ASAN_ENV = {'ASAN_OPTIONS': 'detect_leaks=1:abort_on_error=1:hard_rss_limit_mb=3072',
            'UBSAN_OPTIONS': 'halt_on_error=1:print_stacktrace=1'}


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def resource_guard(build=False):
    cpus = os.sched_getaffinity(0)
    allowed = [{6}, {7}] + ([{6, 7}] if build else [])
    need(cpus in allowed, f'Power-Flow core lease required (affinity {sorted(cpus)})')
    need(os.getpriority(os.PRIO_PROCESS, 0) == 15, 'nice 15 required')
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
        need(os.environ.get(name) == '1', 'one thread required: ' + name)
    need(os.environ.get('CUDA_VISIBLE_DEVICES', None) == '', 'CUDA_VISIBLE_DEVICES must be empty')
    return sorted(cpus)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Runner:
    def __init__(self, out, label):
        self.out = Path(out)
        self.out.mkdir(parents=True, exist_ok=True)
        self.label = label
        self.records = []

    def run(self, command, log, timeout=1800, env=None, sanitized=False, limit=True, check=True):
        prefix = [FROZEN]
        if limit and not sanitized:
            prefix += ['prlimit', f'--as={RLIMIT_AS}', '--']
        full = prefix + [str(x) for x in command]
        environment = dict(os.environ)
        if sanitized:
            environment.update(ASAN_ENV)
        if env:
            environment.update(env)
        start = time.monotonic()
        with (self.out / log).open('w') as handle:
            p = subprocess.run(full, stdout=handle, stderr=subprocess.STDOUT, timeout=timeout,
                               env=environment, cwd=str(P))
        record = {'log': log, 'command': full, 'returncode': p.returncode,
                  'elapsed_s': round(time.monotonic() - start, 3),
                  'cpus': sorted(os.sched_getaffinity(0)), 'nice': os.getpriority(os.PRIO_PROCESS, 0),
                  'rlimit_as': RLIMIT_AS if (limit and not sanitized) else None,
                  'sanitizer_rss_limit_mb': 3072 if sanitized else None}
        self.records.append(record)
        (self.out / (self.label + '-commands.json')).write_text(json.dumps(self.records, indent=2) + '\n')
        if check and p.returncode:
            print((self.out / log).read_text()[-12000:])
            raise RuntimeError(f'failed {log} exit {p.returncode}')
        return p.returncode


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()
