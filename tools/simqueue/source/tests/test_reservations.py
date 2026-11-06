"""SimQueue 0.8.1: backfill reservations for several blocked jobs (pure; fake clock).

0.8.0 reserved only the first blocked job (EASY backfill).  On zenobia on
2026-11-06 a 16-core priority-10 job held that reservation 8.5 h ahead while
five 4-core priority-5 jobs waited 1.5 h behind it: every core that freed up
went to one of some 30 one-core priority-0 probes (2 h limits), none of which
delayed the distant head, so four cores were never free at once.  These
tests replay that queue, check that depth 1 is exactly the 0.8.0 planner
(the 0.8.0 runtime's own scheduler.py as the reference), that no reserved
head is ever delayed, and the diagnostics, borrowing and GPU-line paths.
"""

from __future__ import annotations

import importlib
import importlib.util
import json
import random
import sys
import time
import unittest
from dataclasses import asdict, fields, replace
from pathlib import Path

from simqueue.nodeconfig import ConfigError, SchedulerConfig, load_node_config, parse_node_config
from simqueue.scheduler import INFINITY, Allocation, Pool, Request, Running, build_pool, plan

RUNTIME = Path(__file__).resolve().parents[2]
RUNTIME_080 = Path("/home/treefor/projects/runtimes/simqueue-0.8.0-20261106/source")
G = 1024**3
MIN = 60
HOUR = 3600
T0 = 1_800_000_000.0  # fake clock


def zenobia16(**extra) -> dict:
    """The 16-core pool as on zenobia (two L3 domains of 8 cores, the GPU,
    the old partition names as aliases), given without sysfs."""
    document = {
        "schema": "simqueue-node-config-v1", "node": "z", "default_partition": "zenobia",
        "memory_schedulable": "28G",
        "topology": {"l3_domains": ["0-7", "8-15"]},
        "gpus": [{"uuid": "GPU-z", "index": 0, "memory": "12G", "budget": "11G"}],
        "partitions": [{"name": "zenobia", "cpus": "0-15", "gpus": ["GPU-z"], "default_mem": "4G",
                        "max_mem": "20G", "default_time": "2-00:00:00", "max_time": "14-00:00:00",
                        "aliases": [{"name": "spheral", "default_cpus": 4, "prefer_cpus": [0, 1, 2, 3]},
                                    {"name": "powerflow", "default_cpus": 2, "prefer_cpus": [6, 7]}]}],
    }
    document.update(extra)
    return document


def three_partitions(**extra) -> dict:
    """Partitions with GPUs shared between two of them (0.6.2 layout shape)."""
    document = {
        "schema": "simqueue-node-config-v1", "node": "t", "default_partition": "a",
        "memory_schedulable": "16G",
        "gpus": [{"uuid": "GPU-x", "index": 0, "memory": "12G", "budget": "10G"}],
        "partitions": [
            {"name": "a", "cpus": "0-3", "default_mem": "1G", "max_mem": "8G",
             "default_time": "1:00:00", "max_time": "10:00:00"},
            {"name": "g", "cpus": "4-5", "gpus": ["GPU-x"], "default_mem": "1G", "max_mem": "8G",
             "default_time": "1:00:00", "max_time": "10:00:00"},
            {"name": "h", "cpus": "6-7", "gpus": ["GPU-x"], "default_mem": "1G", "max_mem": "8G",
             "default_time": "1:00:00", "max_time": "10:00:00"},
        ],
    }
    document.update(extra)
    return document


def parse(document: dict):
    return parse_node_config(document, path=Path("/dev/null"), check_topology=False)


def job(job_id: int, ncpus: int, *, minutes: float | None = 60, mem: int = G, partition: str = "zenobia",
        **more) -> Request:
    return Request(job_id=job_id, partition=partition, ncpus=ncpus, mem=mem,
                   time_limit=None if minutes is None else int(minutes * MIN), **more)


def started_ids(decision) -> list[int]:
    return [request.job_id for request, _ in decision.starts]


def simulate(config, running: list[Running], queue: list[Request], *, depth: int | None, until: float,
             now: float = T0, early: random.Random | None = None):
    """Plan at every job end on a fake clock.  Jobs run to their time limit
    (with *early*, to a random 20-100% of it).  Returns the start time of
    every started job and each pass's (time, decision)."""
    running = list(running)
    ends = {item.job_id: item.end_time for item in running}
    waiting = list(queue)
    starts: dict[int, float] = {}
    passes = []
    while waiting and now <= until:
        decision = plan(config, build_pool(config, running), running, waiting, now, depth=depth)
        passes.append((now, decision))
        for request, allocation in decision.starts:
            limit = float(request.time_limit)
            running.append(Running(request.job_id, allocation, now + limit))
            ends[request.job_id] = now + (limit if early is None else max(MIN, early.uniform(0.2, 1.0) * limit))
            starts[request.job_id] = now
        waiting = [r for r in waiting if r.job_id not in starts]
        upcoming = [ends[item.job_id] for item in running if ends[item.job_id] > now]
        if not upcoming:
            break
        now = min(upcoming)
        running = [item for item in running if ends[item.job_id] > now]
    return starts, passes


class SchedulerConfigTests(unittest.TestCase):
    def test_depth_default_and_bounds(self) -> None:
        self.assertEqual(parse(three_partitions()).scheduler, SchedulerConfig(reservation_depth=8))
        self.assertEqual(parse(three_partitions(scheduler={"reservation_depth": 1, "notes": "x"}))
                         .scheduler.reservation_depth, 1)
        for bad in (0, 65, True, "8", 2.5, None):
            with self.subTest(bad=bad), self.assertRaises(ConfigError):
                parse(three_partitions(scheduler={"reservation_depth": bad}))
        with self.assertRaises(ConfigError):
            parse(three_partitions(scheduler={"depth": 4}))
        with self.assertRaises(ConfigError):
            parse(three_partitions(scheduler=[8]))

    def test_shipped_16_core_config_sets_depth_8(self) -> None:
        document = json.loads((RUNTIME / "config" / "zenobia-16core-node.json").read_text())
        self.assertEqual(document["scheduler"]["reservation_depth"], 8)
        config = load_node_config(RUNTIME / "config" / "zenobia-16core-node.json")
        self.assertEqual(config.scheduler.reservation_depth, 8)

    def test_explicit_depth_overrides_the_configuration(self) -> None:
        config = parse(zenobia16())
        running = [Running(100, Allocation(cpus=tuple(range(16)), mem=G), T0 + HOUR)]
        queue = [job(i, 16) for i in range(1, 12)]
        for depth, count in ((None, 8), (1, 1), (3, 3), (20, 11)):
            decision = plan(config, build_pool(config, running), running, queue, T0, depth=depth)
            self.assertEqual(len(decision.reservations), count, depth)
            self.assertEqual(decision.reservation, decision.reservations[0])


class LiveStarvationScenarioTests(unittest.TestCase):
    """The queue of 2026-11-06 ~02:00Z: head 1402 (priority 10, 16 cores,
    15G, shadow 8.5 h away), five merge-gating jobs (priority 5, 4 cores, 30
    min) and a stream of one-core probes (priority 0, 2 h)."""

    def setUp(self) -> None:
        self.config = parse(zenobia16())
        # four long jobs keep the head's shadow 8.5 h away; twelve probes,
        # started earlier, end every 10 minutes over the next two hours
        self.running = [Running(900 + i, Allocation(cpus=(i,), mem=2 * G), T0 + 8.5 * HOUR) for i in range(4)]
        self.running += [Running(910 + i, Allocation(cpus=(4 + i,), mem=G), T0 + (i + 1) * 10 * MIN)
                         for i in range(12)]
        self.head = Request(job_id=1402, partition="zenobia", ncpus=16, mem=15 * G, time_limit=3 * HOUR)
        self.gating = [job(1500 + i, 4, minutes=30, mem=4 * G) for i in range(5)]
        self.probes = [job(2000 + i, 1, minutes=120, mem=G) for i in range(40)]
        self.queue = [self.head, *self.gating, *self.probes]  # dispatch order: priority 10, 5, 0

    def test_first_pass_reserves_the_head_and_the_gating_jobs(self) -> None:
        pool = build_pool(self.config, self.running)
        decision = plan(self.config, pool, self.running, self.queue, T0)
        self.assertEqual(decision.reservation, (1402, T0 + 8.5 * HOUR))
        reserved = [job_id for job_id, _ in decision.reservations]
        self.assertEqual(reserved, [1402, 1500, 1501, 1502, 1503, 1504, 2000, 2001])
        shadows = dict(decision.reservations)
        # FIFO, each on four cores no reservation ahead of it holds then:
        # 4-7 at 40 min; 4-7 again when job 1500 ends (70 min); 8-11 at 80 min
        # while job 1501 holds 4-7; ...
        self.assertEqual([(shadows[1500 + i] - T0) / MIN for i in range(5)], [40, 70, 80, 100, 110])

    def test_a_probe_on_a_freed_core_is_refused_and_names_the_gating_job(self) -> None:
        running = self.running[:4] + self.running[5:]  # the first probe ended: one free core
        now = T0 + 10 * MIN
        decision = plan(self.config, build_pool(self.config, running), running, self.queue, now)
        self.assertEqual(started_ids(decision), [])
        self.assertEqual(decision.skipped[2000], "backfill would delay job 1500")
        self.assertEqual(decision.delays[2000], 1500)
        self.assertEqual(decision.skipped[1500], "resources busy")
        # 0.8.0 (depth 1) gives the core to the first probe
        easy = plan(self.config, build_pool(self.config, running), running, self.queue, now, depth=1)
        self.assertEqual(started_ids(easy), [2000])
        self.assertEqual(easy.reservations, ((1402, T0 + 8.5 * HOUR),))

    def test_gating_jobs_start_within_one_probe_limit_and_the_head_is_not_delayed(self) -> None:
        starts, passes = simulate(self.config, self.running, self.queue, depth=None, until=T0 + 12 * HOUR)
        shadow = passes[0][1].reservation[1]
        self.assertEqual(shadow, T0 + 8.5 * HOUR)
        gating = [starts[1500 + i] - T0 for i in range(5)]
        self.assertTrue(all(wait <= 2 * HOUR for wait in gating), gating)  # one probe time limit
        self.assertEqual(gating, sorted(gating), "same-priority jobs start in FIFO order")
        self.assertLessEqual(starts[1402], shadow)  # the head is not delayed at all
        # the probes still run in the gaps: the pool is not drained for the head early
        self.assertGreaterEqual(sum(1 for j in starts if j >= 2000), 12)

    def test_depth_1_reproduces_the_starvation(self) -> None:
        starts, passes = simulate(self.config, self.running, self.queue, depth=1, until=T0 + 12 * HOUR)
        self.assertEqual(passes[0][1].reservation, (1402, T0 + 8.5 * HOUR))
        first = min(starts.get(1500 + i, INFINITY) for i in range(5)) - T0
        self.assertGreater(first, 4 * HOUR)  # every freed core went to a probe
        self.assertLessEqual(starts[1402], T0 + 8.5 * HOUR)

    def test_head_never_starts_later_when_jobs_end_early(self) -> None:
        for seed in range(6):
            with self.subTest(seed=seed):
                starts, passes = simulate(self.config, self.running, self.queue, depth=None,
                                          until=T0 + 12 * HOUR, early=random.Random(seed))
                first = passes[0][1].reservation[1]
                head_shadows = [d.reservation[1] for _, d in passes if d.reservation and d.reservation[0] == 1402]
                self.assertEqual(head_shadows, sorted(head_shadows, reverse=True), "head shadow never moves later")
                self.assertLessEqual(starts[1402], first)


class NoStarvationTests(unittest.TestCase):
    """Random closed workloads: the head's start never moves later than its
    first shadow time, every job starts, and every start fits the node."""

    def workload(self, rng: random.Random, config, width: int):
        running = []
        cpus = list(range(width))
        rng.shuffle(cpus)
        for job_id in range(500, 500 + rng.randint(1, 6)):
            take = [cpus.pop() for _ in range(min(len(cpus), rng.randint(1, 4)))]
            if not take:
                break
            running.append(Running(job_id, Allocation(cpus=tuple(sorted(take)), mem=rng.choice([1, 2, 4]) * G),
                                   T0 + rng.choice([10, 30, 60, 120, 300]) * MIN))
        queue = []
        for job_id in range(1, rng.randint(8, 40)):
            queue.append(job(job_id, rng.choice([1, 1, 1, 2, 4, 8, width]), minutes=rng.choice([5, 20, 30, 60, 120, 600]),
                             mem=rng.choice([1, 2, 4, 8]) * G))
        return running, queue

    def test_random_workloads(self) -> None:
        rng = random.Random(81)
        for trial in range(60):
            document = zenobia16() if trial % 2 else three_partitions()
            config = parse(document)
            width = 16 if trial % 2 else 4
            running, queue = self.workload(rng, config, width)
            if trial % 2 == 0:
                queue = [replace(r, partition="a", ncpus=min(r.ncpus, 4)) for r in queue]
            for depth in (1, 2, 8):
                with self.subTest(trial=trial, depth=depth):
                    starts, passes = simulate(config, running, queue, depth=depth, until=T0 + 4000 * HOUR,
                                              early=random.Random(trial))
                    self.assertEqual(sorted(starts), sorted(r.job_id for r in queue), "every job starts")
                    heads: dict[int, float] = {}
                    for now, decision in passes:
                        if decision.reservation is not None:
                            head, shadow = decision.reservation
                            heads.setdefault(head, shadow)
                            self.assertLessEqual(shadow, heads[head])
                        self.assertLessEqual(len(decision.reservations), depth)
                    for head, shadow in heads.items():
                        self.assertLessEqual(starts[head], shadow, f"head {head} started after its shadow")


class FairnessTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = parse(zenobia16())
        # cores 14-15 free; 0-1 free at 1 h; 2-13 at 4 h (the head's shadow)
        self.running = [Running(100, Allocation(cpus=(0, 1), mem=G), T0 + HOUR),
                        Running(101, Allocation(cpus=tuple(range(2, 14)), mem=G), T0 + 4 * HOUR)]

    def test_same_priority_wide_job_is_protected_from_a_narrower_later_one(self) -> None:
        head = job(1, 16, minutes=60)
        wide = job(2, 4, minutes=30)  # 0-1 and 14-15 at 1 h
        narrow = job(3, 1, minutes=120)  # same priority as job 2, submitted after it; ends before the head
        decision = plan(self.config, build_pool(self.config, self.running), self.running, [head, wide, narrow], T0)
        self.assertEqual(decision.reservations[:2], ((1, T0 + 4 * HOUR), (2, T0 + HOUR)))
        self.assertEqual(decision.skipped[3], "backfill would delay job 2")
        easy = plan(self.config, build_pool(self.config, self.running), self.running, [head, wide, narrow], T0,
                    depth=1)
        self.assertEqual(started_ids(easy), [3])  # 0.8.0: it only had to spare the head

    def test_reservations_are_placed_after_the_ones_ahead(self) -> None:
        # job 2 (4 cores, 6 h) cannot run before the head without delaying it:
        # its reservation follows the head's, so short jobs may use the free
        # cores until then (independent shadows would hold them from 1 h on)
        head = job(1, 16, minutes=60)
        long = job(2, 4, minutes=360)
        short = job(3, 2, minutes=90)
        decision = plan(self.config, build_pool(self.config, self.running), self.running, [head, long, short], T0)
        self.assertEqual(decision.reservations[:2], ((1, T0 + 4 * HOUR), (2, T0 + 5 * HOUR)))
        self.assertEqual(started_ids(decision), [3])

    def test_two_reservations_do_not_count_the_same_cores(self) -> None:
        head = job(1, 16, minutes=60)
        first = job(2, 3, minutes=30)
        second = job(3, 3, minutes=30)
        decision = plan(self.config, build_pool(self.config, self.running), self.running, [head, first, second], T0)
        shadows = dict(decision.reservations)
        self.assertEqual(shadows[2], T0 + HOUR)  # three of 0-1, 14-15
        self.assertEqual(shadows[3], T0 + HOUR + 30 * MIN)  # the fourth is not enough: after job 2


class EasyEquivalenceTests(unittest.TestCase):
    """Depth 1 decides exactly as the 0.8.0 planner (its runtime's own
    scheduler.py): starts, allocations, reservation and every reason."""

    @classmethod
    def setUpClass(cls) -> None:
        if not (RUNTIME_080 / "simqueue" / "scheduler.py").exists():
            raise unittest.SkipTest("0.8.0 runtime not present")
        name = "simqueue080"
        if name not in sys.modules:
            spec = importlib.util.spec_from_file_location(
                name, RUNTIME_080 / "simqueue" / "__init__.py",
                submodule_search_locations=[str(RUNTIME_080 / "simqueue")])
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
        cls.old = importlib.import_module(f"{name}.scheduler")
        cls.old_config = importlib.import_module(f"{name}.nodeconfig")

    def convert(self, pool: Pool, running: list[Running], queue: list[Request]):
        old = self.old
        old_pool = old.Pool(**{f.name: getattr(pool, f.name) for f in fields(old.Pool)}).copy()
        old_running = [old.Running(r.job_id, old.Allocation(**asdict(r.allocation)), r.end_time, r.mem_actual)
                       for r in running]
        old_queue = [old.Request(**asdict(r)) for r in queue]
        return old_pool, old_running, old_queue

    @staticmethod
    def summary(decision):
        return ([(r.job_id, a.cpus, a.gpus, a.gpu_claim, a.gpu_exclusive, a.mem) for r, a in decision.starts],
                decision.reservation, decision.skipped)

    def random_case(self, rng: random.Random, config, gpu: str):
        names = list(config.partitions) + list(config.aliases)
        used: set[int] = set()
        running = []
        for job_id in range(300, 300 + rng.randint(0, 10)):
            free = [c for c in config.node_cpus() if c not in used]
            if not free:
                break
            cpus = tuple(sorted(rng.sample(free, rng.randint(1, min(4, len(free))))))
            used |= set(cpus)
            part = config.partition_of(cpus[0])
            gpus = (gpu,) if part is not None and part.gpus and rng.random() < 0.25 else ()
            exclusive = bool(gpus) and rng.random() < 0.5
            claim = config.gpus[gpu].budget_bytes if exclusive else (rng.choice([2, 3, 4]) * G if gpus else 0)
            mem = rng.choice([1, 2, 4, 8]) * G
            end = rng.choice([INFINITY, T0, T0 + 300, T0 + 600, T0 + HOUR, T0 + HOUR, T0 + 3 * HOUR])
            actual = None if rng.random() < 0.5 else rng.randint(0, mem)
            running.append(Running(job_id, Allocation(cpus=cpus, gpus=gpus, gpu_claim=claim, gpu_exclusive=exclusive,
                                                      mem=mem), end, actual))
        external = {c for c in config.node_cpus() if c not in used and rng.random() < 0.1}
        actual_free = None if rng.random() < 0.5 else rng.choice([2, 6, 10, 20]) * G
        pool = build_pool(config, running, external_cpus=external,
                          external_gpus={gpu} if rng.random() < 0.1 else (), actual_free=actual_free)
        pool.soft_cpus = {c for c in pool.free_cpus if rng.random() < 0.15}
        queue = []
        for job_id in range(1, rng.randint(1, 30)):
            name = rng.choice(names)
            part = config.find_partition(name)
            ncpus = rng.choice([1, 1, 2, 3, 4, len(part.cpus) // 2 or 1, len(part.cpus)])
            more: dict = {}
            roll = rng.random()
            if roll < 0.08:
                more["cpu_list"] = tuple(sorted(rng.sample(part.cpus, min(ncpus, len(part.cpus)))))
            elif roll < 0.14 and len(config.topology.domains) > 1:
                more["l3_domains"] = 2
                ncpus = max(2, ncpus)
            elif roll < 0.18 and len(config.topology.domains) > 1:
                more["split_domains"] = (rng.randint(1, 6), rng.randint(1, 6))
                ncpus = sum(more["split_domains"])
            if part.gpus and rng.random() < 0.25:
                more["gpus"] = 1
                more["gpu_mem"] = rng.choice([None, 2 * G, 4 * G])
            queue.append(Request(job_id=job_id, partition=name, ncpus=ncpus, mem=rng.choice([1, 2, 4, 8, 15]) * G,
                                 time_limit=rng.choice([None, 300, 600, 1800, HOUR, 3 * HOUR, 30 * HOUR]),
                                 blocked=rng.random() < 0.03, **more))
        return pool, running, queue

    def test_random_queues_match_080(self) -> None:
        rng = random.Random(801)
        documents = [three_partitions(), zenobia16(),
                     zenobia16(topology={"l3_domains": ["0-7", "8-15"], "l3_placement": False}),
                     zenobia16(topology={"l3_domains": ["0-3", "4-7", "8-11", "12-15"]})]
        cases = 0
        refused = 0  # cases where the backfill test refused a job
        backfilled = 0  # cases where a job started behind the reservation
        for index, document in enumerate(documents):
            config = parse(document)
            old_config = self.old_config.parse_node_config(document, path=Path("/dev/null"), check_topology=False)
            gpu = next(iter(config.gpus))
            for trial in range(1500):
                pool, running, queue = self.random_case(rng, config, gpu)
                backfill = rng.random() > 0.05
                old_pool, old_running, old_queue = self.convert(pool, running, queue)
                new = plan(config, pool, running, queue, T0, backfill=backfill, depth=1)
                old = self.old.plan(old_config, old_pool, old_running, old_queue, T0, backfill=backfill)
                with self.subTest(config=index, trial=trial):
                    self.assertEqual(self.summary(new), self.summary(old))
                    self.assertEqual(new.reservations, () if new.reservation is None else (new.reservation,))
                cases += 1
                refused += any(r.startswith("backfill would delay") for r in old.skipped.values())
                if old.reservation is not None:
                    order = [r.job_id for r in queue]
                    head = order.index(old.reservation[0])
                    backfilled += any(order.index(r.job_id) > head for r, _ in old.starts)
        self.assertEqual(cases, 6000)
        # the backfill test itself was exercised both ways
        self.assertGreater(refused, 300)
        self.assertGreater(backfilled, 300)

    def test_existing_scenarios_with_depth_1_against_080(self) -> None:
        """The planner calls of the other test modules' fixtures, replayed."""
        from test_scheduler import small_config

        config = small_config()
        running = [Running(100, Allocation(cpus=(0, 1), mem=G), T0 + 30 * MIN)]
        for queue in ([job(1, 4, partition="a"), job(2, 2, partition="a", minutes=20)],
                      [job(1, 4, partition="a"), job(2, 1, partition="a", minutes=31)],
                      [job(1, 3, partition="a"), job(2, 1, partition="a", minutes=600),
                       job(3, 1, partition="a", minutes=600)]):
            pool = build_pool(config, running)
            old_pool, old_running, old_queue = self.convert(pool, running, queue)
            self.assertEqual(self.summary(plan(config, pool, running, queue, T0, depth=1)),
                             self.summary(self.old.plan(config, old_pool, old_running, old_queue, T0)))


class NodePlannerTests(unittest.TestCase):
    """borrow.plan_node and the GPU line use the same reservations."""

    def test_plan_node_without_borrow_jobs_is_plan(self) -> None:
        from simqueue.borrow import plan_node

        config = parse(zenobia16())
        rng = random.Random(7)
        for _ in range(100):
            running = [Running(100, Allocation(cpus=tuple(range(rng.randint(4, 15))), mem=G),
                               T0 + rng.choice([10, 60, 600]) * MIN)]
            queue = [job(i, rng.choice([1, 2, 4, 16]), minutes=rng.choice([10, 60, 300])) for i in range(1, 20)]
            a = plan(config, build_pool(config, running), running, queue, T0)
            b = plan_node(config, build_pool(config, running), running, queue, T0, owner_of={})
            self.assertEqual((started_ids(a), a.reservations, a.skipped),
                             (started_ids(b), b.reservations, b.skipped))

    def test_borrowers_use_cores_held_for_a_second_reservation_and_yield_to_it(self) -> None:
        from simqueue.borrow import Borrower, CoreState, plan_node

        document = zenobia16(borrow={"lend_reserved": True})
        config = parse(document)
        running = [Running(1, Allocation(cpus=tuple(range(12)), mem=G), T0 + 4 * HOUR),
                   Running(2, Allocation(cpus=(12,), mem=G), T0 + HOUR)]
        head = job(10, 16, minutes=60)
        second = job(11, 4, minutes=30)  # reservation 2: cores 12-15 at 1 h
        borrow = Request(job_id=20, partition="zenobia", ncpus=3, mem=G, time_limit=5 * HOUR, borrow=True)
        cores = {cpu: CoreState(cpu) for cpu in range(16)}
        owner_of = {c: 1 for c in range(12)} | {12: 2}
        decision = plan_node(config, build_pool(config, running), running, [head, second, borrow], T0,
                             cores=cores, owner_of=owner_of)
        self.assertEqual([j for j, _ in decision.reservations], [10, 11])
        lent = {r.job_id: a for r, a in decision.starts if a.borrow}
        self.assertEqual(lent[20].cpus, (13, 14, 15))  # the free cores held for job 11
        # at 1 h job 2 ends: job 11 starts on 12-15 and preempts the borrower
        now = T0 + HOUR
        running = running[:1]
        borrower = Borrower(job_id=20, home="zenobia", cpus=(13, 14, 15), mem=G, deadline=T0 + 5 * HOUR)
        decision = plan_node(config, build_pool(config, running), running, [head, second], now,
                             borrowers=[borrower], cores=cores, owner_of={c: 1 for c in range(12)})
        self.assertEqual(started_ids(decision), [11])
        self.assertEqual([(p.job_id, p.by_job) for p in decision.preemptions], [(20, 11)])

    def test_gpu_job_that_would_delay_a_later_reservation_is_not_ready(self) -> None:
        from simqueue.gpuqueue import ready_except_gpu

        document = {
            "schema": "simqueue-node-config-v1", "node": "t", "default_partition": "g",
            "memory_schedulable": "16G",
            "gpus": [{"uuid": "GPU-x", "index": 0, "memory": "12G", "budget": "10G"}],
            "partitions": [{"name": "g", "cpus": "0-3", "gpus": ["GPU-x"], "default_mem": "1G", "max_mem": "8G",
                            "default_time": "1:00:00", "max_time": "10:00:00"}],
        }
        running = [Running(1, Allocation(cpus=(0, 1), mem=G), T0 + HOUR)]
        first = job(10, 3, partition="g", minutes=30)    # reservation 1 at 1 h
        second = job(11, 4, partition="g", minutes=30)   # reservation 2 at 1 h 30
        gpu_job = job(12, 1, partition="g", minutes=105, gpus=1)  # still running at 1 h 30
        for depth, ready in ((1, [12]), (8, [])):
            config = replace(parse(document), scheduler=SchedulerConfig(reservation_depth=depth))
            found = ready_except_gpu(config, "GPU-x", running=running, owners=[first, second, gpu_job],
                                     busy_cpus=(), blocked={"GPU-x"}, started=set(), now=T0)
            self.assertEqual(found, ready, depth)


class ParkingTests(unittest.TestCase):
    """Parked jobs (a memory request the node cannot admit) are left out of
    the planner's queue by the worker (test_memory_units.HeadOfLineTests
    runs that with real units); infeasible and held jobs get no reservation
    from the planner itself."""

    def test_infeasible_and_held_jobs_hold_no_reservation(self) -> None:
        config = parse(zenobia16())
        running = [Running(1, Allocation(cpus=tuple(range(16)), mem=G), T0 + HOUR)]
        queue = [job(10, 4, mem=30 * G), job(11, 17), job(12, 16, minutes=60), job(13, 4)]
        decision = plan(config, build_pool(config, running), running, queue, T0)
        self.assertEqual([j for j, _ in decision.reservations], [12, 13])
        self.assertIn("schedulable memory", decision.skipped[10])
        paused = [job(20, 4, blocked=True), job(21, 4)]
        decision = plan(config, build_pool(config, running), running, paused, T0)
        self.assertEqual((decision.reservations, decision.skipped[21]), ((), "held behind a paused ensemble"))


class DiagnosticsTests(unittest.TestCase):
    def worker(self):
        from simqueue.nodeworker import NodeWorker

        worker = NodeWorker.__new__(NodeWorker)
        worker._reasons, worker._reason_written, worker._reserved_since = {}, {}, {}
        worker.gpu = type("G", (), {"reasons": {}})()
        written: dict[int, str] = {}
        worker.store = type("S", (), {"set_message": lambda self, job, text: written.__setitem__(job, text)})()
        return worker, written

    def test_every_reserved_job_says_when_it_is_reserved_to_start(self) -> None:
        worker, written = self.worker()
        skipped = {1: "resources busy", 2: "resources busy", 3: "backfill would delay job 2", 4: "resources busy"}
        worker._record_reasons(skipped, ((1, 4102444800.0), (2, 4102448400.0), (4, INFINITY)), {})
        self.assertEqual(written[1], "pending: resources busy; reserved to start by 2100-01-01T00:00:00Z; "
                                     "head of the queue for 0 min")
        self.assertEqual(written[2], "pending: resources busy; reserved to start by 2100-01-01T01:00:00Z; "
                                     "reservation 2 of the queue for 0 min")
        self.assertEqual(written[3], "pending: backfill would delay job 2")
        self.assertEqual(written[4], "pending: resources busy; reservation 3 of the queue for 0 min")
        # a job that leaves the queue is forgotten; one still waiting keeps its time
        worker._reserved_since[1] -= 600
        worker._record_reasons({1: "resources busy"}, ((1, 4102444800.0),), {})
        self.assertIn("head of the queue for 10 min", written[1])
        self.assertEqual(set(worker._reserved_since), {1})

    def test_single_reservation_form_is_still_accepted(self) -> None:
        worker, written = self.worker()
        worker._record_reasons({5: "resources busy"}, (5, 4102444800.0), {})
        self.assertIn("reserved to start by 2100-01-01T00:00:00Z; head of the queue for", written[5])
        worker._record_reasons({6: "resources busy"}, None, {})
        self.assertEqual(written[6], "pending: resources busy")


class PlanCostTests(unittest.TestCase):
    """plan() stays far below the worker's 1 s poll with ~100 queued and 16
    running jobs, also when many queued jobs fit now and must be tested
    against every reservation (tools/bench-plan.py measures it in detail)."""

    def test_100_queued_16_running(self) -> None:
        config = parse(zenobia16())
        rng = random.Random(3)
        running = [Running(100 + c, Allocation(cpus=(c,), mem=G), T0 + rng.choice([10, 60, 120, 600]) * MIN)
                   for c in range(12)]
        running += [Running(200 + c, Allocation(cpus=(c,), mem=G), T0 + 8 * HOUR) for c in (12, 13)]
        queue = [job(1, 16, mem=15 * G, minutes=180)]
        queue += [job(i, rng.choice([1, 2, 3, 4, 8]), minutes=rng.choice([10, 20, 30, 60, 120, 540]),
                      mem=rng.choice([1, 2, 3, 4, 6, 8]) * G) for i in range(2, 101)]
        pool = build_pool(config, running)
        timings = []
        for _ in range(5):
            began = time.perf_counter()
            decision = plan(config, pool, running, queue, T0)
            timings.append(time.perf_counter() - began)
        self.assertEqual(len(decision.reservations), 8)
        self.assertLess(min(timings), 0.25, timings)


if __name__ == "__main__":
    unittest.main()
