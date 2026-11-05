#!/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen
"""Run the SimQueue test suite while holding the cores' lease files.

Usage: run-tests.py [--cores LIST] [--wait-any [--min-cores N]] [--shared] [--receipt PATH] [PATTERN ...]

Tests run on cores 0-1 (the run refuses to start and stops when
<pic-profile>/RETURN-CORES-0-3 exists) or 6-7 (as any-core.sh); cores 2-5
are never used.  Temporary files go to $SIMQ_TEST_TMPDIR (default:
<outputs>/simq063-test-tmp), never /tmp.  The
directory is used by its resolved path: the outputs directory lies behind a
symbolic link, and the wrapper compares a job's resolved working directory
with its project root, as do the submit fingerprints.  test_switchover.py
needs a directory whose path does not contain "simqueue-0.6." (switchover.py
takes that text in the sandbox unit as "0.6 already installed"), for example
SIMQ_TEST_TMPDIR=/home/treefor/projects/.workspace/outputs/simq062-test-tmp.

The runner takes the flock(2) lease file of every core it uses (the same
testing-core<K>.lock files as any-core.sh and pf-launch, ascending order),
pins itself and every test process to those cores at nice 15, and exports
SIMQ_TEST_CPUS so the tests place their real jobs only there.  Tests that
need lease files of their own use temporary ones.

While the tests run, a guard checks every two seconds for work that needs
the leased cores: an owner job queued or running in the live queue database
(read-only) under a name whose cores include a leased one (0.6.2 partitions:
spheral 0-3, powerflow 6-7; 0.6.3: the pool zenobia, every core, and its
aliases with their old cores), or a process outside the
test run whose CPU affinity lies within the leased cores (work pinned there
without a lease, as Spheral runs and builds are).  It then
stops every test unit (``simqt*``), releases the leases and exits with
status 3.  The same check refuses to start the run.  ``--no-live-guard``
disables it (sandbox use only).

``--shared`` (cores 6-7, shared with Power-Flow work under the lease
protocol): queue jobs of the cores' partition are no reason to refuse or
stop; the queue never starts a job on a core whose lease the run holds and
starts it once the run releases the lease.  Work pinned to the leased cores
without a lease, and the RETURN flag, still stop the run.
``--wait-any --min-cores N`` waits until at least N of the listed cores are
free at the same moment (nothing is held while waiting).

0.7.0 modes:

``--in-job``  inside a live queue job (``simq run -c 2 ... -- run-tests.py
              --in-job``): the job's cores (``SIMQ_CPU_LIST``) are the test
              cores; the job unit already holds their lease files, so none
              is taken and no live guard runs (the queue arbitrates).
``--pure``    no test unit at all (``SIMQ_TEST_CPUS`` unset: unit-level
              tests skip); no lease, no guard.
``--sandbox-cpus LIST``  development only: test units on LIST (for example
              idle SMT siblings) without leases or guard.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import os
import sqlite3
import subprocess
import sys
import threading
import time
import unittest
from contextlib import closing
from pathlib import Path


HERE = Path(__file__).resolve().parent
LEASE_DIR = Path("/home/treefor/projects/.workspace/outputs/pic-profile")
LEASE_OVERRIDES = {6: "testing-core.lock"}
FORBIDDEN = {2, 3, 4, 5}  # Spheral's own cores and the raikou-mfem allocation: never used by these tests
LIVE_DATABASE = Path("/home/treefor/.local/state/simqueue/queue.sqlite3")
GUARDED_PARTITIONS = ("spheral",)
GUARDED_CORES = {0, 1, 2, 3}  # the spheral partition (0.6.2) and alias (0.6.3)
# 0.6.3: zenobia is one pool of every core; the old names are its aliases
CORE_GROUPS = {"spheral": {0, 1, 2, 3}, "powerflow": {6, 7}, "zenobia": set(range(8))}
RETURN_FLAG = LEASE_DIR / "RETURN-CORES-0-3"
DEFAULT_TMPDIR = Path("/home/treefor/projects/.workspace/outputs/simq070-test-tmp")
TEST_UNIT_PATTERN = "simqt*"


def guarded_partitions(cores) -> tuple[str, ...]:
    return tuple(name for name, group in CORE_GROUPS.items() if group & set(cores))


def guarded_cores(cores) -> set[int]:
    """Cores whose pinned work counts against the run: the old partition
    groups of the leased cores (the pool's group would be every core)."""
    found: set[int] = set()
    for name, group in CORE_GROUPS.items():
        if name != "zenobia" and group & set(cores):
            found |= group
    return found or set(GUARDED_CORES)


def return_requested(cores) -> list[dict]:
    if set(cores) & {0, 1, 2, 3} and RETURN_FLAG.exists():
        return [{"flag": str(RETURN_FLAG)}]
    return []


def live_jobs_on(database: Path, partitions=GUARDED_PARTITIONS) -> list[dict]:
    """Queued or active owner jobs of *partitions* in a live queue database (read-only)."""
    if not database.exists() or not partitions:
        return []
    marks = ",".join("?" for _ in partitions)
    try:
        with closing(sqlite3.connect(f"file:{database}?mode=ro", uri=True, timeout=2)) as connection:
            columns = {row[1] for row in connection.execute("PRAGMA table_info(jobs)")}
            partition = "COALESCE(partition, 'spheral')" if "partition" in columns else "'spheral'"
            owner = "AND COALESCE(borrow, 0) = 0" if "borrow" in columns else ""
            rows = connection.execute(
                f"SELECT id, name, state, {partition} FROM jobs "
                f"WHERE state IN ('queued', 'starting', 'running', 'cancelling') AND {partition} IN ({marks}) {owner}",
                tuple(partitions),
            ).fetchall()
    except sqlite3.Error:
        return []
    return [{"id": r[0], "name": r[1], "state": r[2], "partition": r[3]} for r in rows]


def _process_table() -> dict[int, tuple[int, str]]:
    """pid -> (ppid, cgroup) for every readable process."""
    table: dict[int, tuple[int, str]] = {}
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            text = Path(f"/proc/{entry}/stat").read_text()
            ppid = int(text[text.rfind(")") + 2 :].split()[1])
            cgroup = Path(f"/proc/{entry}/cgroup").read_text()
        except (OSError, ValueError, IndexError):
            continue
        table[int(entry)] = (ppid, cgroup)
    return table


def foreign_work_on(cores: list[int], *, own_root: int | None = None) -> list[dict]:
    """Processes of this user, outside this test run, pinned within *cores*,
    or within the spheral cores while touching one of *cores*: work placed
    there without taking a lease.

    A process belongs to the test run when it or an ancestor is the run
    itself or lives in a test unit (``simqt*``), so work that a test job
    moves into a scope of its own (bounded-job.py) is not foreign, or when it
    carries the run's ``SIMQ_TEST_RUN`` token (a sandbox worker that was
    reparented after its starter exited).
    """
    token = f"SIMQ_TEST_RUN={own_root}".encode() if own_root else None
    leased = set(cores)
    table = _process_table()
    cache: dict[int, bool] = {}

    def ours(pid: int) -> bool:
        chain = []
        current = pid
        result = False
        for _ in range(64):
            if current in cache:
                result = cache[current]
                break
            if current == own_root or (current in table and "/simqt" in table[current][1]):
                result = True
                break
            if current not in table or current <= 1:
                break
            chain.append(current)
            current = table[current][0]
        for item in chain:
            cache[item] = result
        return result

    uid = os.getuid()
    found = []
    for pid in table:
        base = Path("/proc") / str(pid)
        try:
            if base.stat().st_uid != uid:
                continue
            command = base.joinpath("cmdline").read_bytes().replace(b"\0", b" ").decode("utf-8", "replace").strip()
            if not command:
                continue  # kernel thread
            allowed = None
            for line in base.joinpath("status").read_text().splitlines():
                if line.startswith("Cpus_allowed_list:"):
                    allowed = set(parse_cores(line.split(":", 1)[1].strip()))
        except (OSError, ValueError):
            continue
        # work pinned within the leased cores, or within the spheral cores and
        # touching a leased one (daemons allowed on 0-7 or everywhere are not)
        if not allowed or not allowed & leased or ours(pid):
            continue
        if not (allowed <= leased or allowed <= guarded_cores(cores)):
            continue
        if token is not None:
            try:
                if token in base.joinpath("environ").read_bytes().split(b"\0"):
                    continue
            except OSError:
                pass
        found.append({"pid": pid, "allowed": sorted(allowed), "command": command[:160]})
    return found


def kill_run_processes(own_root: int) -> None:
    """SIGKILL every other process carrying this run's token (test workers
    started in their own sessions, sandbox workers)."""
    import signal

    token = f"SIMQ_TEST_RUN={own_root}".encode()
    for entry in os.listdir("/proc"):
        if not entry.isdigit() or int(entry) == own_root:
            continue
        try:
            if token in Path(f"/proc/{entry}/environ").read_bytes().split(b"\0"):
                os.kill(int(entry), signal.SIGKILL)
        except (OSError, ProcessLookupError):
            continue


class LiveGuard(threading.Thread):
    """Abort the test run when the live queue wants the leased cores."""

    def __init__(self, database: Path, cores: list[int], interval: float = 2.0, *, shared: bool = False) -> None:
        super().__init__(daemon=True)
        self.database = database
        self.cores = cores
        self.interval = interval
        self.shared = shared

    def run(self) -> None:
        while True:
            jobs = (([] if self.shared else live_jobs_on(self.database, guarded_partitions(self.cores)))
                    + foreign_work_on(self.cores, own_root=os.getpid()) + return_requested(self.cores))
            if jobs:
                print(f"\nLIVE GUARD: other work needs the leased cores: {jobs}; stopping the tests now",
                      file=sys.stderr, flush=True)
                for verb in ("thaw", "stop"):
                    subprocess.run(["systemctl", "--user", verb, TEST_UNIT_PATTERN], stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL, check=False, timeout=60)
                kill_run_processes(os.getpid())
                listed = subprocess.run(["systemctl", "--user", "list-units", "simqt*.slice", "--all", "--no-legend",
                                         "--plain"], text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                                        check=False, timeout=60).stdout
                for line in listed.splitlines():
                    name = line.split()[0] if line.split() else ""
                    if name.startswith("simqt") and name != "simqtest.slice":
                        for verb in ("revert", "stop"):
                            subprocess.run(["systemctl", "--user", verb, name], stdout=subprocess.DEVNULL,
                                           stderr=subprocess.DEVNULL, check=False, timeout=60)
                os._exit(3)  # closes the lease descriptors
            time.sleep(self.interval)


def lease_path(core: int) -> Path:
    return LEASE_DIR / LEASE_OVERRIDES.get(core, f"testing-core{core}.lock")


def parse_cores(text: str) -> list[int]:
    cores: set[int] = set()
    for part in text.split(","):
        low, _, high = part.partition("-")
        cores.update(range(int(low), int(high or low) + 1))
    return sorted(cores)


def free_by_proc_locks(cores: list[int]) -> list[int]:
    """Cores whose lease file nobody holds, read from /proc/locks."""
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    from simqueue.locks import lock_holders

    holders = lock_holders([str(lease_path(core)) for core in cores])
    return [core for core in cores if not holders.get(str(lease_path(core)))]


def take_leases(cores: list[int], *, wait_any: bool, min_cores: int = 1) -> tuple[list[int], list[int]]:
    held: list[int] = []
    fds: list[int] = []
    if wait_any:
        while True:
            if min_cores > 1 and len(free_by_proc_locks(cores)) < min_cores:
                time.sleep(2)  # look without touching the lease files until enough are free
                continue
            for core in cores:
                fd = os.open(lease_path(core), os.O_RDONLY | os.O_CREAT | os.O_APPEND, 0o644)
                try:
                    fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    os.close(fd)
                    continue
                held.append(core)
                fds.append(fd)
            if len(held) >= min_cores:
                return held, fds
            # too few free at once: hold nothing while waiting
            for fd in fds:
                os.close(fd)
            held, fds = [], []
            time.sleep(2)
    for core in cores:
        fd = os.open(lease_path(core), os.O_RDONLY | os.O_CREAT | os.O_APPEND, 0o644)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print(f"waiting for {lease_path(core)} ...", file=sys.stderr, flush=True)
            fcntl.flock(fd, fcntl.LOCK_EX)
        held.append(core)
        fds.append(fd)
    return held, fds


def run_suites(args: argparse.Namespace, held: list[int] | None, lease_files: list[str]) -> int:
    sys.path.insert(0, str(HERE))
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    for pattern in args.patterns:
        suite.addTests(loader.discover(str(HERE / "tests"), pattern=pattern))
    started = time.time()
    result = unittest.TextTestRunner(verbosity=2 if args.verbose else 1).run(suite)
    receipt = {
        "schema": "simqueue-test-receipt-v1",
        "started_unix": started,
        "duration_s": round(time.time() - started, 3),
        "cores_leased": held,
        "lease_files": lease_files,
        "mode": "in-job" if args.in_job else ("pure" if args.pure else ("sandbox" if args.sandbox_cpus else "leases")),
        "queue_job": os.environ.get("SIMQ_JOB_ID"),
        "patterns": args.patterns,
        "tests_run": result.testsRun,
        "failures": len(result.failures),
        "errors": len(result.errors),
        "skipped": len(result.skipped),
        "skipped_detail": [f"{test.id()}: {reason}" for test, reason in result.skipped],
        "failed_detail": [test.id() for test, _ in result.failures + result.errors],
        "successful": result.wasSuccessful(),
        "python": sys.version,
    }
    if args.receipt:
        Path(args.receipt).write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: receipt[k] for k in ("tests_run", "failures", "errors", "skipped", "cores_leased",
                                              "duration_s", "mode")}))
    return 0 if result.wasSuccessful() else 1


def isolate_state(tmpdir: Path) -> None:
    """No test process may reach the live queue: a path that lacks
    SIMQUEUE_STATE_DIR (a child given a reduced environment) falls back to
    XDG_STATE_HOME, which points into the test directory as well."""
    sandbox = tmpdir / f"default-state-{os.getpid()}"
    sandbox.mkdir(parents=True, exist_ok=True)
    os.environ["XDG_STATE_HOME"] = str(sandbox)
    os.environ.setdefault("SIMQUEUE_STATE_DIR", str(sandbox / "simqueue"))


def run_without_leases(args: argparse.Namespace) -> int:
    tmpdir = Path(os.environ.get("SIMQ_TEST_TMPDIR") or DEFAULT_TMPDIR).resolve()
    tmpdir.mkdir(parents=True, exist_ok=True)
    isolate_state(tmpdir)
    os.environ["TMPDIR"] = str(tmpdir)
    import tempfile

    tempfile.tempdir = str(tmpdir)
    os.environ["SIMQ_TEST_RUN"] = str(os.getpid())
    if args.receipts_dir:
        os.environ["SIMQ_TEST_RECEIPTS"] = args.receipts_dir
    held: list[int] | None = None
    if args.in_job:
        cpu_list = os.environ.get("SIMQ_CPU_LIST", "")
        if not cpu_list or not os.environ.get("SIMQ_JOB_ID"):
            print("--in-job needs to run inside a queue job (SIMQ_JOB_ID, SIMQ_CPU_LIST)", file=sys.stderr)
            return 64
        held = parse_cores(cpu_list)
        os.environ["SIMQ_TEST_CPUS"] = ",".join(str(c) for c in held)
    elif args.sandbox_cpus:
        held = parse_cores(args.sandbox_cpus)
        os.environ["SIMQ_TEST_CPUS"] = ",".join(str(c) for c in held)
    else:
        os.environ.pop("SIMQ_TEST_CPUS", None)
    os.environ.pop("SIMQ_TEST_LEASES", None)
    # the suites start their own sandbox workers and units: none of the
    # job's own SIMQ_* placement may leak into them
    keep = ("SIMQ_TEST", "SIMQ_LIVE_DB_COPY", "SIMQ_UPGRADE_UNSEALED")
    for key in [k for k in os.environ if k.startswith("SIMQ_") and not k.startswith(keep)]:
        os.environ.pop(key, None)
    try:
        os.nice(max(0, 10 - os.nice(0)))
    except OSError:
        pass
    return run_suites(args, held, [])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cores", default="0-1")
    parser.add_argument("--wait-any", action="store_true", help="use whichever listed cores are free")
    parser.add_argument("--receipt")
    parser.add_argument("--receipts-dir", help="directory for test observation receipts (SIMQ_TEST_RECEIPTS)")
    parser.add_argument("--no-live-guard", action="store_true")
    parser.add_argument("--shared", action="store_true",
                        help="cores shared under the lease protocol: queue jobs of their partition do not stop the run")
    parser.add_argument("--min-cores", type=int, default=1, help="with --wait-any: wait until this many are free")
    parser.add_argument("--in-job", action="store_true", help="run inside a live queue job on its cores")
    parser.add_argument("--pure", action="store_true", help="no test units, no leases")
    parser.add_argument("--sandbox-cpus", help="development: test units on these CPUs, no leases")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("patterns", nargs="*", default=["test_*.py"])
    args = parser.parse_args()
    if args.in_job or args.pure or args.sandbox_cpus:
        return run_without_leases(args)
    cores = parse_cores(args.cores)
    if FORBIDDEN & set(cores):
        parser.error("cores 2-5 are never used by these tests (Spheral's own cores and raikou-mfem's)")
    tmpdir = Path(os.environ.get("SIMQ_TEST_TMPDIR") or DEFAULT_TMPDIR).resolve()
    tmpdir.mkdir(parents=True, exist_ok=True)
    isolate_state(tmpdir)
    os.environ["TMPDIR"] = str(tmpdir)
    import tempfile

    tempfile.tempdir = str(tmpdir)
    if not args.no_live_guard:
        waiting = return_requested(cores)
        if waiting:
            print(f"cores 0-3 must be returned ({RETURN_FLAG}); not starting the tests", file=sys.stderr)
            return 3
    held, _fds = take_leases(cores, wait_any=args.wait_any, min_cores=min(args.min_cores, len(cores)))
    if not args.no_live_guard:
        waiting = (([] if args.shared else live_jobs_on(LIVE_DATABASE, guarded_partitions(held)))
                   + foreign_work_on(held, own_root=os.getpid())
                   + return_requested(held))
        if waiting:
            print(f"other work needs the cores: {waiting}; not starting the tests", file=sys.stderr)
            return 3
    os.environ["SIMQ_TEST_RUN"] = str(os.getpid())
    if not args.no_live_guard:
        LiveGuard(LIVE_DATABASE, held, shared=args.shared).start()
    if args.receipts_dir:
        os.environ["SIMQ_TEST_RECEIPTS"] = args.receipts_dir
    os.sched_setaffinity(0, held)
    os.nice(max(0, 15 - os.nice(0)))
    os.environ["SIMQ_TEST_CPUS"] = ",".join(str(c) for c in held)
    os.environ["SIMQ_TEST_LEASES"] = ",".join(str(lease_path(c)) for c in held)
    sys.path.insert(0, str(HERE))
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    for pattern in args.patterns:
        suite.addTests(loader.discover(str(HERE / "tests"), pattern=pattern))
    started = time.time()
    result = unittest.TextTestRunner(verbosity=2 if args.verbose else 1).run(suite)
    receipt = {
        "schema": "simqueue-test-receipt-v1",
        "started_unix": started,
        "duration_s": round(time.time() - started, 3),
        "cores_leased": held,
        "lease_files": [str(lease_path(c)) for c in held],
        "patterns": args.patterns,
        "tests_run": result.testsRun,
        "failures": len(result.failures),
        "errors": len(result.errors),
        "skipped": len(result.skipped),
        "skipped_detail": [f"{test.id()}: {reason}" for test, reason in result.skipped],
        "failed_detail": [test.id() for test, _ in result.failures + result.errors],
        "successful": result.wasSuccessful(),
        "python": sys.version,
    }
    if args.receipt:
        Path(args.receipt).write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: receipt[k] for k in ("tests_run", "failures", "errors", "skipped", "cores_leased", "duration_s")}))
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
