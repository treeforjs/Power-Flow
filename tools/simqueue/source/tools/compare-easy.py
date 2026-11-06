#!/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen
"""Run the test suite with every planner call checked against the 0.8.0 planner.

Every call of ``simqueue.scheduler.plan`` made in this process -- by the
tests, by ``borrow.plan_node`` and by the GPU line's what-if plans -- runs
twice: 0.8.1 at reservation depth 1 and the 0.8.0 runtime's own ``plan()``.
Their decisions (every start with its allocation, the reservation, every
reason) must be equal; the call returns the depth-1 decision, so the suite
runs as it did with 0.8.0.  A difference is recorded in the receipt and
raised as an AssertionError in the calling test.  Planner calls made by
worker processes that the unit-level tests start are not seen here.

Usage: compare-easy.py [--receipt PATH] [run-tests.py arguments ...]
       (default: --pure and every test module but test_reservations.py,
       whose own tests are about depths above 1, and test_shield_helper.py,
       which cannot be imported without source/privileged)
"""

from __future__ import annotations

import collections
import importlib
import importlib.util
import json
import runpy
import sys
from dataclasses import asdict, fields
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[1]
RUNTIME_080 = Path("/home/treefor/projects/runtimes/simqueue-0.8.0-20261106/source")


def load_080():
    name = "simqueue080"
    spec = importlib.util.spec_from_file_location(name, RUNTIME_080 / "simqueue" / "__init__.py",
                                                  submodule_search_locations=[str(RUNTIME_080 / "simqueue")])
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return importlib.import_module(f"{name}.scheduler")


def summary(decision) -> tuple:
    return ([(r.job_id, a.cpus, a.gpus, a.gpu_claim, a.gpu_exclusive, a.mem, a.borrow, a.placement)
             for r, a in decision.starts], decision.reservation, decision.skipped)


def main() -> int:
    arguments = sys.argv[1:]
    receipt_path = None
    if arguments[:1] == ["--receipt"]:
        receipt_path, arguments = arguments[1], arguments[2:]
    if not any(not a.startswith("-") for a in arguments):
        skip = {"test_reservations.py", "test_shield_helper.py"}  # depths above 1; not importable here
        modules = sorted(p.name for p in (SOURCE / "tests").glob("test_*.py") if p.name not in skip)
        arguments = [*(arguments or ["--pure"]), *modules]
    sys.path.insert(0, str(SOURCE))
    from simqueue import scheduler

    old = load_080()
    ours = scheduler.plan
    counts: collections.Counter = collections.Counter()
    callers: collections.Counter = collections.Counter()
    mismatches: list[dict] = []

    def checked(config, pool, running, queue, now, *, backfill=True, depth=None):
        frame = sys._getframe(1)
        callers[f"{Path(frame.f_code.co_filename).name}:{frame.f_code.co_name}"] += 1
        decision = ours(config, pool, running, queue, now, backfill=backfill, depth=1)
        old_pool = old.Pool(**{f.name: getattr(pool, f.name) for f in fields(old.Pool)}).copy()
        old_running = [old.Running(r.job_id, old.Allocation(**asdict(r.allocation)), r.end_time, r.mem_actual)
                       for r in running]
        old_queue = [old.Request(**asdict(r)) for r in queue]
        reference = old.plan(config, old_pool, old_running, old_queue, now, backfill=backfill)
        counts["calls"] += 1
        counts["queued_jobs"] += len(queue)
        counts["with_reservation"] += decision.reservation is not None
        counts["backfill_refusals"] += sum(1 for r in decision.skipped.values() if r.startswith("backfill would"))
        if summary(decision) != summary(reference):
            mismatches.append({"caller": f"{frame.f_code.co_filename}:{frame.f_lineno}",
                               "ours": repr(summary(decision)), "080": repr(summary(reference))})
            raise AssertionError(f"0.8.1 at depth 1 differs from 0.8.0: {summary(decision)!r} != "
                                 f"{summary(reference)!r}")
        return decision

    scheduler.plan = checked
    sys.argv = [str(SOURCE / "run-tests.py"), *arguments]
    code = 0
    try:
        runpy.run_path(str(SOURCE / "run-tests.py"), run_name="__main__")
    except SystemExit as stop:
        code = int(stop.code or 0)
    receipt = {"schema": "simqueue-easy-comparison-v1", "reference": str(RUNTIME_080), "arguments": arguments,
               "plan_calls_compared": counts["calls"], "queued_jobs_planned": counts["queued_jobs"],
               "calls_with_a_reservation": counts["with_reservation"],
               "backfill_refusals": counts["backfill_refusals"], "mismatches": mismatches,
               "callers": dict(callers.most_common()), "suite_exit_code": code}
    print(json.dumps({k: v for k, v in receipt.items() if k != "callers"}, indent=2))
    if receipt_path:
        Path(receipt_path).write_text(json.dumps(receipt, indent=2) + "\n")
    return 1 if mismatches or code else 0


if __name__ == "__main__":
    raise SystemExit(main())
