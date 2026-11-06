#!/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen
"""Measure the planner's cost per worker pass (0.8.1).

Scenarios on the 16-core node configuration (config/zenobia-16core-node.json):

  live         the owner jobs of a queue database (a read-only copy given
               with --db), queued and running, planned as the worker plans
               them
  synthetic    100 queued jobs of zenobia-like shapes behind a 16-core head,
               14 running jobs ending over the next 10 h, two cores free
  adversarial  92 queued jobs of 92 different shapes that fit the two free
               cores now and end before the head's reservation, so each is
               tested against the seven reservations behind it (the refusal
               cache cannot help)

Each is planned with plan() at depths 1, 8 and 16, by the 0.8.0 planner, and
as one worker pass (plan_node plus the GPU line's two what-if plans).
Prints and writes (--out) median, p95 and max milliseconds.

Usage: bench-plan.py [--db QUEUE.sqlite3] [--repeat N] [--out FILE]
"""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import json
import random
import sqlite3
import statistics
import sys
import time
from contextlib import closing
from dataclasses import asdict, fields
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[1]
RUNTIME = SOURCE.parent
RUNTIME_080 = Path("/home/treefor/projects/runtimes/simqueue-0.8.0-20261106/source")
sys.path.insert(0, str(SOURCE))

from simqueue.borrow import plan_node  # noqa: E402
from simqueue.gpuqueue import ready_except_gpu  # noqa: E402
from simqueue.nodeconfig import load_node_config  # noqa: E402
from simqueue.nodeworker import ACTIVE_STATES, job_request, running_entries  # noqa: E402
from simqueue.scheduler import Allocation, Request, Running, build_pool, plan  # noqa: E402

G = 1024**3
NOW = 1_800_000_000.0


def load_080():
    if not (RUNTIME_080 / "simqueue" / "scheduler.py").exists():
        return None
    name = "simqueue080"
    spec = importlib.util.spec_from_file_location(name, RUNTIME_080 / "simqueue" / "__init__.py",
                                                  submodule_search_locations=[str(RUNTIME_080 / "simqueue")])
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return importlib.import_module(f"{name}.scheduler")


def live(config, database: Path):
    with closing(sqlite3.connect(f"file:{database}?mode=ro", uri=True)) as connection:
        connection.row_factory = sqlite3.Row
        rows = [dict(r) for r in connection.execute(
            "SELECT * FROM jobs WHERE state IN ('queued', 'starting', 'running', 'cancelling') "
            "ORDER BY borrow, priority DESC, queue_order, id")]
    now = time.time()
    active = [j for j in rows if j["state"] in ACTIVE_STATES and not j.get("borrow")]
    running = running_entries(active, config, now)
    queue = [job_request(j, config) for j in rows if j["state"] == "queued" and not j.get("borrow")]
    return running, queue, now


def synthetic(config, rng: random.Random):
    running = [Running(500 + c, Allocation(cpus=(c,), mem=G), NOW + rng.choice([10, 30, 60, 120, 240, 600]) * 60)
               for c in range(12)]
    running += [Running(600 + c, Allocation(cpus=(c,), mem=2 * G), NOW + 8.5 * 3600) for c in (12, 13)]
    queue = [Request(job_id=1, partition="zenobia", ncpus=16, mem=15 * G, time_limit=3 * 3600)]
    for job_id in range(2, 101):
        queue.append(Request(job_id=job_id, partition=rng.choice(["zenobia", "spheral", "powerflow"]),
                             ncpus=rng.choice([1, 1, 2, 3, 4, 4, 8, 16]), mem=rng.choice([1, 2, 3, 4, 6, 8, 12]) * G,
                             time_limit=rng.choice([10, 20, 30, 60, 120, 540, 1170]) * 60,
                             gpus=1 if rng.random() < 0.06 else 0))
    return running, queue, NOW


def adversarial(config, rng: random.Random):
    running = [Running(500 + c, Allocation(cpus=(c,), mem=G), NOW + (c + 1) * 600) for c in range(12)]
    running += [Running(600 + c, Allocation(cpus=(c,), mem=G), NOW + 8.5 * 3600) for c in (12, 13)]
    # the head's shadow is 8.5 h away; seven 4-core jobs reserve the cores
    # freeing up over the next two hours; 92 jobs of 92 shapes fit the two
    # free cores now, end before the head and must be tested against each
    queue = [Request(job_id=1, partition="zenobia", ncpus=16, mem=8 * G, time_limit=3600)]
    queue += [Request(job_id=1 + i, partition="zenobia", ncpus=4, mem=G, time_limit=1800 + i) for i in range(1, 8)]
    for job_id in range(9, 101):
        queue.append(Request(job_id=job_id, partition=rng.choice(["zenobia", "spheral", "powerflow"]),
                             ncpus=rng.choice([1, 2]), mem=rng.choice([1, 2]) * G, time_limit=3600 + 60 * job_id))
    return running, queue, NOW


def measure(function, repeat: int) -> dict:
    timings = []
    for _ in range(repeat):
        began = time.perf_counter()
        function()
        timings.append((time.perf_counter() - began) * 1000.0)
    timings.sort()
    return {"median_ms": round(statistics.median(timings), 3),
            "p95_ms": round(timings[min(len(timings) - 1, int(0.95 * len(timings)))], 3),
            "max_ms": round(timings[-1], 3)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db")
    parser.add_argument("--repeat", type=int, default=200)
    parser.add_argument("--out")
    args = parser.parse_args()
    config = load_node_config(RUNTIME / "config" / "zenobia-16core-node.json")
    old = load_080()
    rng = random.Random(81)
    scenarios = {"synthetic": synthetic(config, rng), "adversarial": adversarial(config, rng)}
    if args.db:
        scenarios = {"live": live(config, Path(args.db)), **scenarios}
    gpu = next(iter(config.gpus), None)
    result: dict = {"config": str(config.path), "repeat": args.repeat, "scenarios": {}}
    for name, (running, queue, now) in scenarios.items():
        pool = build_pool(config, running)
        entry: dict = {"queued": len(queue), "running": len(running), "free_cores": len(pool.free_cpus)}
        for depth in (1, 8, 16):
            decision = plan(config, pool, running, queue, now, depth=depth)
            entry[f"depth_{depth}"] = {
                **measure(lambda: plan(config, pool, running, queue, now, depth=depth), args.repeat),
                "starts": len(decision.starts), "reservations": len(decision.reservations),
                "refused": sum(1 for r in decision.skipped.values() if r.startswith("backfill would"))}
        if old is not None:
            old_pool = old.Pool(**{f.name: getattr(pool, f.name) for f in fields(old.Pool)})
            old_running = [old.Running(r.job_id, old.Allocation(**asdict(r.allocation)), r.end_time, r.mem_actual)
                           for r in running]
            old_queue = [old.Request(**asdict(r)) for r in queue]
            entry["planner_080"] = measure(lambda: old.plan(config, old_pool, old_running, old_queue, now),
                                           args.repeat)

        def worker_pass():
            plan_node(config, pool, running, queue, now, owner_of={})
            if gpu is not None:
                ready_except_gpu(config, gpu, running=running, owners=queue, busy_cpus=(), blocked=(gpu,),
                                 started=set(), now=now)

        entry["worker_pass_depth_8"] = measure(worker_pass, max(20, args.repeat // 4))
        result["scenarios"][name] = entry
    text = json.dumps(result, indent=2)
    print(text)
    if args.out:
        Path(args.out).write_text(text + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
