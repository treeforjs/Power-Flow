#!/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen
"""Hold the lease files of the given cores, pin to them at nice 15, run a command.

Usage: with-leases.py --cores LIST [--wait-any] -- COMMAND [ARGS...]
Used for validation runs; the same lease files as any-core.sh/pf-launch.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import importlib.util

spec = importlib.util.spec_from_file_location("run_tests", Path(__file__).resolve().parents[1] / "run-tests.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def main() -> int:
    args = sys.argv[1:]
    wait_any = "--wait-any" in args
    args = [a for a in args if a != "--wait-any"]
    if len(args) < 4 or args[0] != "--cores" or args[2] != "--":
        print(__doc__, file=sys.stderr)
        return 64
    cores = runner.parse_cores(args[1])
    if runner.FORBIDDEN & set(cores):
        print("cores 4-5 belong to the raikou-mfem allocation", file=sys.stderr)
        return 64
    held, _fds = runner.take_leases(cores, wait_any=wait_any)
    os.sched_setaffinity(0, held)
    os.nice(max(0, 15 - os.nice(0)))
    env = dict(os.environ, SIMQ_TEST_CPUS=",".join(map(str, held)))
    return subprocess.call(args[3:], env=env)


if __name__ == "__main__":
    raise SystemExit(main())
