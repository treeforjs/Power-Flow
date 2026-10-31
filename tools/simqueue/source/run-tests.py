#!/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen
"""Run the SimQueue test suite while holding the cores' lease files.

Usage: run-tests.py [--cores LIST] [--wait-any] [--receipt PATH] [PATTERN ...]

The runner takes the flock(2) lease file of every core it uses (the same
testing-core<K>.lock files as any-core.sh and pf-launch, ascending order),
pins itself and every test process to those cores at nice 15, and exports
SIMQ_TEST_CPUS so the tests place their real jobs only there.  Tests that
need lease files of their own use temporary ones.

While the tests run, a guard checks every two seconds for work that needs
the leased cores: a job of the live spheral partition (cores 0-3) queued or
running in the live queue database (read-only), or a process outside the
test run whose CPU affinity lies within the leased cores (work pinned there
without a lease, as Spheral runs and builds are).  It then
stops every test unit (``simqt*``), releases the leases and exits with
status 3.  The same check refuses to start the run.  ``--no-live-guard``
disables it (sandbox use only).
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
FORBIDDEN = {4, 5}  # raikou-mfem cores: never used by these tests
LIVE_DATABASE = Path("/home/treefor/.local/state/simqueue/queue.sqlite3")
GUARDED_PARTITIONS = ("spheral",)
GUARDED_CORES = {0, 1, 2, 3}  # the spheral partition
TEST_UNIT_PATTERN = "simqt*"


def live_jobs_on(database: Path, partitions=GUARDED_PARTITIONS) -> list[dict]:
    """Queued or active jobs of *partitions* in a live queue database (read-only)."""
    if not database.exists():
        return []
    marks = ",".join("?" for _ in partitions)
    try:
        with closing(sqlite3.connect(f"file:{database}?mode=ro", uri=True, timeout=2)) as connection:
            columns = {row[1] for row in connection.execute("PRAGMA table_info(jobs)")}
            partition = "COALESCE(partition, 'spheral')" if "partition" in columns else "'spheral'"
            rows = connection.execute(
                f"SELECT id, name, state, {partition} FROM jobs "
                f"WHERE state IN ('queued', 'starting', 'running', 'cancelling') AND {partition} IN ({marks})",
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
        if not (allowed <= leased or allowed <= GUARDED_CORES):
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

    def __init__(self, database: Path, cores: list[int], interval: float = 2.0) -> None:
        super().__init__(daemon=True)
        self.database = database
        self.cores = cores
        self.interval = interval

    def run(self) -> None:
        while True:
            jobs = live_jobs_on(self.database) + foreign_work_on(self.cores, own_root=os.getpid())
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


def take_leases(cores: list[int], *, wait_any: bool) -> tuple[list[int], list[int]]:
    held: list[int] = []
    fds: list[int] = []
    if wait_any:
        while not held:
            for core in cores:
                fd = os.open(lease_path(core), os.O_RDONLY | os.O_CREAT | os.O_APPEND, 0o644)
                try:
                    fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    os.close(fd)
                    continue
                held.append(core)
                fds.append(fd)
            if not held:
                time.sleep(2)
        return held, fds
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cores", default="0-3")
    parser.add_argument("--wait-any", action="store_true", help="use whichever listed cores are free")
    parser.add_argument("--receipt")
    parser.add_argument("--receipts-dir", help="directory for test observation receipts (SIMQ_TEST_RECEIPTS)")
    parser.add_argument("--no-live-guard", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("patterns", nargs="*", default=["test_*.py"])
    args = parser.parse_args()
    cores = parse_cores(args.cores)
    if FORBIDDEN & set(cores):
        parser.error("cores 4-5 belong to the raikou-mfem allocation")
    if not args.no_live_guard:
        waiting = live_jobs_on(LIVE_DATABASE) + foreign_work_on(cores, own_root=os.getpid())
        if waiting:
            print(f"other work needs the cores: {waiting}; not starting the tests", file=sys.stderr)
            return 3
    held, _fds = take_leases(cores, wait_any=args.wait_any)
    os.environ["SIMQ_TEST_RUN"] = str(os.getpid())
    if not args.no_live_guard:
        LiveGuard(LIVE_DATABASE, held).start()
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
