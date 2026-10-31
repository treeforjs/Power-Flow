#!/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen
"""Run the SimQueue test suite while holding the cores' lease files.

Usage: run-tests.py [--cores LIST] [--wait-any] [--receipt PATH] [PATTERN ...]

The runner takes the flock(2) lease file of every core it uses (the same
testing-core<K>.lock files as any-core.sh and pf-launch, ascending order),
pins itself and every test process to those cores at nice 15, and exports
SIMQ_TEST_CPUS so the tests place their real jobs only there.  Tests that
need lease files of their own use temporary ones.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import os
import sys
import time
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
LEASE_DIR = Path("/home/treefor/projects/.workspace/outputs/pic-profile")
LEASE_OVERRIDES = {6: "testing-core.lock"}
FORBIDDEN = {4, 5}  # raikou-mfem cores: never used by these tests


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
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("patterns", nargs="*", default=["test_*.py"])
    args = parser.parse_args()
    cores = parse_cores(args.cores)
    if FORBIDDEN & set(cores):
        parser.error("cores 4-5 belong to the raikou-mfem allocation")
    held, _fds = take_leases(cores, wait_any=args.wait_any)
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
