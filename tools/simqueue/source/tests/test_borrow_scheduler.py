"""Borrowing (0.6.1): pure planning, usage classification and configuration.

Every test uses an explicit (fake) clock value and fake utilization figures;
nothing is launched.
"""

from __future__ import annotations

import importlib
import importlib.util
import json
import random
import sys
import unittest
from collections import deque
from pathlib import Path

from simqueue.borrow import Borrower, CoreState, eligible_partitions, in_use_cores, plan_node, protected_cpus
from simqueue.nodeconfig import ConfigError, load_node_config, parse_node_config
from simqueue.scheduler import Allocation, Pool, Request, Running, build_pool, plan
from simqueue.telemetry import (
    HolderGroup,
    ProcInfo,
    UsageTracker,
    classify_external,
    classify_owner,
    descendants,
    parse_stat,
    read_cpu_times,
)


G = 1024**3
RUNTIME = Path(__file__).resolve().parents[2]
ZENOBIA_CONFIG = RUNTIME / "config" / "zenobia-node.json"
RUNTIME_060 = Path("/home/treefor/projects/runtimes/simqueue-0.6.0-20261030/source")


def document(policies=None, *, borrow=None, extra_partitions=(), **extra):
    policies = policies or {}
    doc = {
        "schema": "simqueue-node-config-v1",
        "node": "test",
        "default_partition": "a",
        "memory_schedulable": "16G",
        "allow_smt_siblings": True,
        "gpus": [{"uuid": "GPU-x", "index": 0, "memory": "12G", "budget": "10G"}],
        "smt_siblings": {str(c): [c + 8] for c in range(8)},
        "borrow": {"idle_window_s": 10, "busy_window_s": 2, "freeze_max_s": 1800, **(borrow or {})},
        "partitions": [
            {"name": "a", "cpus": "0-3", "default_mem": "1G", "max_mem": "8G",
             "default_time": "1:00:00", "max_time": "10:00:00", **policies.get("a", {})},
            {"name": "g", "cpus": "4-5", "gpus": ["GPU-x"], "default_mem": "1G", "max_mem": "8G",
             "default_time": "1:00:00", "max_time": "10:00:00", **policies.get("g", {})},
            {"name": "h", "cpus": "6-7", "gpus": ["GPU-x"], "default_mem": "1G", "max_mem": "8G",
             "default_time": "1:00:00", "max_time": "10:00:00", **policies.get("h", {})},
            *extra_partitions,
        ],
    }
    doc.update(extra)
    return doc


def config(policies=None, **kwargs):
    return parse_node_config(document(policies, **kwargs), path=Path("/dev/null"), check_topology=False)


def req(job_id, partition="a", ncpus=1, *, borrow=False, max_cpus=None, borrow_from=None, mem=G, minutes=60,
        exclusive=False, on_preempt="requeue", cpus=None, gpus=0, priority=0):
    return Request(
        job_id=job_id, partition=partition, ncpus=ncpus, cpu_list=cpus, gpus=gpus, mem=mem,
        time_limit=None if minutes is None else minutes * 60, borrow=borrow, max_cpus=max_cpus,
        borrow_from=borrow_from, exclusive_timing=exclusive, on_preempt=on_preempt,
    )


def owner(job_id, cpus, *, mem=G, end=10_000.0):
    return Running(job_id, Allocation(cpus=tuple(cpus), mem=mem), end)


def borrower(job_id, cpus, *, home="h", mode="requeue", classes=None, owners=None, frozen=None, pending=None,
             mem=G, order=None, holder_keys=None, deadline=None, borrow_from=None):
    return Borrower(
        job_id=job_id, home=home, cpus=tuple(cpus), mem=mem, on_preempt=mode,
        placement=dict(classes or {cpu: "free" for cpu in cpus}), owners=dict(owners or {}),
        holder_keys=dict(holder_keys or {}), frozen_since=frozen, pending=pending,
        borrow_from=borrow_from, deadline=deadline, order=job_id if order is None else order,
    )


def run(cfg, queue, *, running=(), borrowers=(), cores=None, external=(), now=0.0, exclusive=(), **kwargs):
    running = list(running)
    pool = build_pool(cfg, running, external_cpus=external)
    owner_of = {cpu: item.job_id for item in running for cpu in item.allocation.cpus}
    return plan_node(cfg, pool, running, list(queue), now, borrowers=list(borrowers), cores=cores or {},
                     owner_of=owner_of, exclusive_owner_cpus=exclusive, **kwargs)


def started(decision):
    return {request.job_id: allocation for request, allocation in decision.starts}


def preempted(decision):
    return {item.job_id: item for item in decision.preemptions}


class BorrowPlacementTests(unittest.TestCase):
    def test_borrower_prefers_home_then_highest_free_cores(self):
        cfg = config()
        decision = run(cfg, [req(10, "a", 2, borrow=True)])
        allocation = started(decision)[10]
        self.assertTrue(allocation.borrow)
        self.assertEqual(allocation.cpus, (2, 3))
        self.assertEqual({cls for _, _, cls in allocation.placement}, {"free"})

    def test_borrower_spills_to_lending_partitions_in_order(self):
        cfg = config()
        decision = run(cfg, [req(10, "h", 3, borrow=True)], running=[owner(1, [6, 7])])
        allocation = started(decision)[10]
        self.assertEqual(allocation.cpus, (1, 2, 3))
        self.assertEqual({part for _, part, _ in allocation.placement}, {"a"})

    def test_max_borrow_cores_is_elastic_between_cores_and_maximum(self):
        cfg = config()
        decision = run(cfg, [req(10, "a", 1, borrow=True, max_cpus=3, borrow_from=("a",))])
        self.assertEqual(started(decision)[10].cpus, (1, 2, 3))
        decision = run(cfg, [req(10, "a", 1, borrow=True, max_cpus=3, borrow_from=("a",))],
                       running=[owner(1, [0, 1, 2])])
        self.assertEqual(started(decision)[10].cpus, (3,))
        decision = run(cfg, [req(10, "a", 2, borrow=True, borrow_from=("a",))], running=[owner(1, [0, 1, 2])])
        self.assertNotIn(10, started(decision))
        self.assertIn("lendable", decision.skipped[10])

    def test_max_borrow_cores_caps_a_job_even_when_more_cores_are_idle(self):
        cfg = config()
        decision = run(cfg, [req(10, "h", 2, borrow=True, max_cpus=2)])
        self.assertEqual(len(started(decision)[10].cpus), 2)

    def test_borrow_from_restricts_lenders(self):
        cfg = config()
        decision = run(cfg, [req(10, "h", 2, borrow=True, borrow_from=("g",))], running=[owner(1, [6, 7])])
        self.assertEqual(started(decision)[10].cpus, (4, 5))
        decision = run(cfg, [req(10, "h", 3, borrow=True, borrow_from=("g",))], running=[owner(1, [6, 7])])
        self.assertNotIn(10, started(decision))

    def test_waiting_owner_job_keeps_its_partition_cores_from_borrowers(self):
        cfg = config()
        queue = [req(1, "a", 4), req(10, "h", 2, borrow=True)]
        decision = run(cfg, queue, running=[owner(2, [0, 1]), owner(3, [6, 7])])
        self.assertNotIn(1, started(decision))
        self.assertEqual(started(decision)[10].cpus, (4, 5))  # not a's free cores 2-3

    def test_owner_held_by_paused_ensemble_does_not_block_lending(self):
        cfg = config()
        paused = Request(job_id=1, partition="a", ncpus=4, mem=G, time_limit=3600, blocked=True)
        decision = run(cfg, [paused, req(10, "h", 2, borrow=True)], running=[owner(3, [6, 7])])
        self.assertEqual(started(decision)[10].cpus, (2, 3))

    def test_lend_false_lends_to_nobody_else_but_home_jobs_may_borrow(self):
        cfg = config({"a": {"lend": False}})
        decision = run(cfg, [req(10, "h", 2, borrow=True)], running=[owner(3, [6, 7])])
        self.assertEqual(started(decision)[10].cpus, (4, 5))
        decision = run(cfg, [req(11, "a", 2, borrow=True)])
        self.assertEqual(started(decision)[11].cpus, (2, 3))

    def test_timing_sensitive_partition_and_its_smt_siblings_are_never_lent(self):
        extra = ({"name": "s", "cpus": "8-11", "default_mem": "1G", "max_mem": "8G",
                  "default_time": "1:00:00", "max_time": "10:00:00"},)
        cfg = config({"a": {"timing_sensitive": True}}, extra_partitions=extra)
        self.assertTrue({0, 1, 2, 3, 8, 9, 10, 11} <= protected_cpus(cfg))
        decision = run(cfg, [req(10, "h", 1, borrow=True, max_cpus=8)], running=[owner(3, [6, 7])])
        self.assertEqual(started(decision)[10].cpus, (4, 5))
        # not even the partition's own borrow jobs
        decision = run(cfg, [req(11, "a", 1, borrow=True, borrow_from=("g",))])
        self.assertEqual(set(started(decision)[11].cpus) & {0, 1, 2, 3, 8, 9, 10, 11}, set())
        self.assertEqual(eligible_partitions(cfg, req(11, "a", 1, borrow=True)), ["g", "h", "s"])

    def test_exclusive_timing_owner_protects_its_cores_and_siblings(self):
        extra = ({"name": "s", "cpus": "14-15", "default_mem": "1G", "max_mem": "8G",
                  "default_time": "1:00:00", "max_time": "10:00:00"},)
        cfg = config(extra_partitions=extra)
        cores = {6: CoreState(6, owner_idle=True), 7: CoreState(7, owner_idle=True)}
        running = [owner(1, [6, 7]), owner(2, [0, 1, 2, 3]), owner(3, [4, 5])]
        queue = [req(10, "a", 1, borrow=True, max_cpus=4)]
        decision = run(cfg, queue, running=running, cores=cores)
        self.assertEqual(started(decision)[10].cpus, (6, 7, 14, 15))
        decision = run(cfg, queue, running=running, cores=cores, exclusive=[6])
        self.assertEqual(started(decision)[10].cpus, (7, 15))  # 6 and its sibling 14 protected
        decision = run(cfg, [req(11, "a", 3, borrow=True, max_cpus=4)], running=running, cores=cores,
                       exclusive=[6])
        self.assertIn("CPU 6 job 1 (exclusive timing)", decision.skipped[11])
        self.assertIn("CPU 14 SMT sibling of an exclusive-timing job", decision.skipped[11])
        # an exclusive-timing owner starting now preempts a borrower on its sibling
        decision = run(cfg, [req(4, "h", 1, exclusive=True)],
                       borrowers=[borrower(20, [14], home="s")])
        self.assertEqual(started(decision)[4].cpus, (6,))
        self.assertIn(20, preempted(decision))

    def test_owner_idle_cores_need_lend_and_are_ranked_after_free_cores(self):
        cfg = config()
        cores = {7: CoreState(7, owner_idle=True)}
        decision = run(cfg, [req(10, "a", 2, borrow=True, borrow_from=("h",))],
                       running=[owner(1, [6, 7]), owner(2, [0, 1, 2])], cores=cores)
        allocation = started(decision)[10]
        self.assertEqual(allocation.placement, ((3, "a", "free"), (7, "h", "owner-idle")))
        cfg = config({"h": {"lend": False}})
        decision = run(cfg, [req(10, "a", 2, borrow=True)], running=[owner(1, [6, 7]), owner(2, [0, 1, 2])],
                       cores=cores)
        self.assertEqual(started(decision)[10].cpus, (3, 5))

    def test_borrowers_always_sort_after_owners(self):
        cfg = config()
        queue = [req(10, "a", 4, borrow=True, priority=100), req(1, "a", 4)]
        decision = run(cfg, queue)
        self.assertEqual(started(decision)[1].cpus, (0, 1, 2, 3))
        self.assertEqual(started(decision)[10].cpus, (4, 5, 6, 7))  # elsewhere, never a's cores
        decision = run(cfg, [req(10, "a", 4, borrow=True, priority=100, borrow_from=("a",)), req(1, "a", 4)])
        self.assertNotIn(10, started(decision))

    def test_borrowers_never_get_a_reservation_or_delay_owners(self):
        cfg = config()
        queue = [req(1, "a", 4, minutes=10), req(10, "a", 1, borrow=True, minutes=None)]
        decision = run(cfg, queue, running=[owner(2, [0, 1], end=600.0)])
        self.assertEqual(decision.reservation, (1, 600.0))
        self.assertEqual(started(decision)[10].cpus, (5,))  # a's free cores are waited for

    def test_gpu_and_explicit_cpu_borrow_requests_are_infeasible(self):
        cfg = config()
        decision = run(cfg, [req(10, "h", 1, borrow=True, gpus=1), req(11, "a", 1, borrow=True, cpus=(1,))])
        self.assertIn("GPU", decision.skipped[10])
        self.assertIn("explicit", decision.skipped[11])

    def test_borrow_memory_must_fit_the_free_memory(self):
        cfg = config()
        decision = run(cfg, [req(10, "a", 1, borrow=True, mem=10 * G)], running=[owner(1, [6], mem=8 * G)])
        self.assertEqual(decision.skipped[10], "borrow: memory busy")

    def test_free_core_busy_with_unleased_work_is_not_lent(self):
        cfg = config()
        cores = {3: CoreState(3, cpu_busy=True), 2: CoreState(2, cpu_busy=True)}
        decision = run(cfg, [req(10, "a", 1, borrow=True, borrow_from=("a",))], cores=cores)
        self.assertEqual(started(decision)[10].cpus, (1,))
        decision = run(cfg, [req(10, "a", 3, borrow=True, borrow_from=("a",))], cores=cores)
        self.assertIn("CPU 2-3 CPU busy outside the queue", decision.skipped[10])
        # owners are not affected, and a borrower already there is not preempted for it
        decision = run(cfg, [req(1, "a", 4)], cores=cores)
        self.assertEqual(started(decision)[1].cpus, (0, 1, 2, 3))
        decision = run(cfg, [], borrowers=[borrower(20, [3])], cores=cores)
        self.assertEqual(decision.preemptions, [])

    def test_paused_queue_starts_and_resumes_nothing_but_still_preempts(self):
        cfg = config()
        cores = {3: CoreState(3, external=True, holder="pid 1 x")}
        decision = run(cfg, [req(10, "a", 1, borrow=True)], borrowers=[borrower(11, [3], frozen=None),
                                                                      borrower(12, [2], mode="freeze", frozen=1.0)],
                       cores=cores, external=[3], starts_allowed=False)
        self.assertEqual(started(decision), {})
        self.assertEqual(decision.resumes, [])
        self.assertIn(11, preempted(decision))

    def test_disabled_borrowing_starts_no_borrow_job(self):
        cfg = config(borrow={"enabled": False})
        decision = run(cfg, [req(10, "a", 1, borrow=True)])
        self.assertIn("disabled", decision.skipped[10])


class PreemptionTests(unittest.TestCase):
    def test_owner_start_preempts_borrower_with_each_mode_and_never_waits(self):
        cfg = config()
        for mode in ("freeze", "requeue", "cancel"):
            with self.subTest(mode=mode):
                decision = run(cfg, [req(1, "a", 4)], borrowers=[borrower(20, [2, 3], home="h", mode=mode)])
                self.assertEqual(started(decision)[1].cpus, (0, 1, 2, 3))
                item = preempted(decision)[20]
                self.assertEqual((item.mode, item.by_job, item.escalated), (mode, 1, False))
                self.assertIn("needed by job 1", item.reason)

    def test_owner_prefers_free_cores_over_borrowed_ones(self):
        cfg = config()
        decision = run(cfg, [req(1, "a", 2)], borrowers=[borrower(20, [0, 1])])
        self.assertEqual(started(decision)[1].cpus, (2, 3))
        self.assertEqual(decision.preemptions, [])

    def test_frozen_borrower_cores_go_to_owners_without_another_action(self):
        cfg = config()
        frozen = borrower(20, [2, 3], mode="freeze", frozen=50.0)
        decision = run(cfg, [req(1, "a", 4)], borrowers=[frozen], now=100.0)
        self.assertEqual(started(decision)[1].cpus, (0, 1, 2, 3))
        self.assertEqual(decision.preemptions, [])
        self.assertEqual(decision.resumes, [])

    def test_frozen_borrower_resumes_when_its_cores_are_lendable_again(self):
        cfg = config()
        frozen = borrower(20, [2, 3], home="h", mode="freeze", frozen=50.0)
        decision = run(cfg, [], borrowers=[frozen], now=100.0)
        self.assertEqual(decision.resumes, [(20, ((2, "a", "free"), (3, "a", "free")))])
        decision = run(cfg, [], borrowers=[frozen], running=[owner(1, [2, 3])], now=100.0)
        self.assertEqual(decision.resumes, [])
        decision = run(cfg, [req(1, "a", 4)], borrowers=[frozen], running=[owner(2, [0])], now=100.0)
        self.assertEqual(decision.resumes, [])  # an owner of a is waiting

    def test_frozen_borrower_cores_stay_reserved_against_other_borrowers(self):
        cfg = config()
        frozen = borrower(20, [2, 3], home="h", mode="freeze", frozen=50.0)
        decision = run(cfg, [req(21, "h", 1, borrow=True, max_cpus=4, borrow_from=("a",))],
                       borrowers=[frozen], running=[owner(1, [6, 7])], now=100.0)
        self.assertEqual(started(decision)[21].cpus, (0, 1))

    def test_stopping_borrower_frees_cores_for_owners_only(self):
        cfg = config()
        stopping = borrower(20, [2, 3], pending="requeue")
        decision = run(cfg, [req(21, "h", 1, borrow=True, max_cpus=4, borrow_from=("a",)), req(1, "a", 4)],
                       borrowers=[stopping], running=[owner(2, [6, 7])])
        self.assertEqual(started(decision)[1].cpus, (0, 1, 2, 3))
        self.assertNotIn(21, started(decision))
        self.assertEqual(decision.preemptions, [])

    def test_memory_needed_by_an_owner_stops_borrowers_and_escalates_freeze(self):
        cfg = config()
        decision = run(cfg, [req(1, "g", 1, mem=8 * G)],
                       borrowers=[borrower(20, [3], mode="freeze", mem=6 * G, order=1),
                                  borrower(21, [2], mode="requeue", mem=1 * G, order=2)],
                       running=[owner(5, [6], mem=4 * G)])
        self.assertIn(1, started(decision))
        items = preempted(decision)
        # newest first: 21 frees 1G (not enough), then 20 (freeze -> requeue)
        self.assertEqual(items[21].mode, "requeue")
        self.assertEqual((items[20].mode, items[20].escalated), ("requeue", True))
        self.assertIn("memory", items[20].reason)

    def test_frozen_too_long_or_near_its_time_limit_is_requeued(self):
        cfg = config()
        decision = run(cfg, [], borrowers=[borrower(20, [2], mode="freeze", frozen=0.0)], running=[owner(1, [2])],
                       now=1801.0)
        self.assertEqual((preempted(decision)[20].mode, preempted(decision)[20].escalated), ("requeue", True))
        decision = run(cfg, [], borrowers=[borrower(20, [2], mode="freeze", frozen=90.0, deadline=150.0)],
                       running=[owner(1, [2])], now=100.0)
        self.assertIn("time limit", preempted(decision)[20].reason)

    def test_lease_taken_on_a_lent_free_core_preempts_immediately(self):
        cfg = config()
        cores = {3: CoreState(3, external=True, holder="pid 9 any-core.sh", holder_key="9:1")}
        decision = run(cfg, [], borrowers=[borrower(20, [3])], cores=cores, external=[3])
        self.assertIn("taken outside the queue", preempted(decision)[20].reason)

    def test_no_borrow_marker_preempts_and_blocks_lending(self):
        cfg = config()
        cores = {3: CoreState(3, no_borrow=True)}
        decision = run(cfg, [req(10, "a", 1, borrow=True, max_cpus=4)], borrowers=[borrower(20, [3])], cores=cores)
        self.assertIn(20, preempted(decision))
        self.assertNotIn(3, started(decision)[10].cpus)


class LeasedIdleTests(unittest.TestCase):
    def cores(self, *, idle=True, busy=False, stale=False, key="77:5"):
        return {cpu: CoreState(cpu, external=True, external_idle=idle, external_busy=busy, stale=stale,
                               holder="pid 77 any-core.sh --wait python t.py", holder_key=key) for cpu in (6, 7)}

    def queue(self):
        return [req(10, "a", 1, borrow=True, borrow_from=("h",))]

    def test_leased_idle_cores_are_lent_only_with_the_policy(self):
        running = [owner(1, [0, 1, 2, 3])]
        decision = run(config(), self.queue(), running=running, cores=self.cores(), external=[6, 7])
        self.assertNotIn(10, started(decision))
        cfg = config({"h": {"lend_leased_idle": True}})
        decision = run(cfg, self.queue(), running=running, cores=self.cores(), external=[6, 7])
        self.assertEqual(started(decision)[10].placement, ((7, "h", "leased-idle"),))
        decision = run(cfg, self.queue(), running=running, cores=self.cores(idle=False), external=[6, 7])
        self.assertNotIn(10, started(decision))
        decision = run(cfg, self.queue(), running=running, cores=self.cores(stale=True), external=[6, 7])
        self.assertNotIn(10, started(decision))

    def test_busy_or_replaced_holder_preempts_the_borrower(self):
        cfg = config({"h": {"lend_leased_idle": True}})
        lent = borrower(20, [7], home="a", classes={7: "leased-idle"}, holder_keys={7: "77:5"})
        decision = run(cfg, [], borrowers=[lent], cores=self.cores(), external=[6, 7])
        self.assertEqual(decision.preemptions, [])
        decision = run(cfg, [], borrowers=[lent], cores=self.cores(busy=True), external=[6, 7])
        self.assertIn("is busy", preempted(decision)[20].reason)
        decision = run(cfg, [], borrowers=[lent], cores=self.cores(key="78:9"), external=[6, 7])
        self.assertIn("changed hands", preempted(decision)[20].reason)
        decision = run(config(), [], borrowers=[lent], cores=self.cores(), external=[6, 7])
        self.assertIn("no longer lends", preempted(decision)[20].reason)

    def test_owner_needing_its_idle_core_back_preempts(self):
        cfg = config()
        lent = borrower(20, [7], home="a", classes={7: "owner-idle"}, owners={7: 1})
        cores = {6: CoreState(6, owner_idle=False), 7: CoreState(7, owner_idle=True)}
        decision = run(cfg, [], borrowers=[lent], running=[owner(1, [6, 7])], cores=cores)
        self.assertEqual(decision.preemptions, [])
        cores = {cpu: CoreState(cpu, owner_busy=True) for cpu in (6, 7)}
        decision = run(cfg, [], borrowers=[lent], running=[owner(1, [6, 7])], cores=cores)
        self.assertIn("needs CPU 7 again", preempted(decision)[20].reason)


class CompatibilityTests(unittest.TestCase):
    """Without borrow jobs, plan_node returns exactly what 0.6.0's plan does."""

    @classmethod
    def setUpClass(cls):
        if not (RUNTIME_060 / "simqueue" / "scheduler.py").exists():
            raise unittest.SkipTest("0.6.0 runtime not present")
        name = "simqueue060"
        spec = importlib.util.spec_from_file_location(
            name, RUNTIME_060 / "simqueue" / "__init__.py", submodule_search_locations=[str(RUNTIME_060 / "simqueue")]
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        cls.old_scheduler = importlib.import_module(f"{name}.scheduler")
        cls.old_nodeconfig = importlib.import_module(f"{name}.nodeconfig")

    def test_random_queues_match_060(self):
        rng = random.Random(61)
        # the 0.6.0 configuration: both parsers accept it
        doc = json.loads((RUNTIME_060.parent / "config" / "zenobia-node.json").read_text())
        new_cfg = parse_node_config(doc, path=Path("/dev/null"), check_topology=False)
        old_cfg = self.old_nodeconfig.parse_node_config(doc, path=Path("/dev/null"), check_topology=False)
        names = list(new_cfg.partitions)
        gpu = list(new_cfg.gpus)[0]
        for trial in range(400):
            running_new, running_old, used = [], [], set()
            for job_id in range(100, 100 + rng.randint(0, 4)):
                part = new_cfg.partitions[rng.choice(names)]
                free = [c for c in part.cpus if c not in used]
                if not free:
                    continue
                cpus = tuple(sorted(rng.sample(free, rng.randint(1, len(free)))))
                used |= set(cpus)
                gpus = (gpu,) if part.gpus and rng.random() < 0.3 else ()
                kwargs = dict(cpus=cpus, gpus=gpus, gpu_claim=new_cfg.gpus[gpu].budget_bytes if gpus else 0,
                              gpu_exclusive=bool(gpus), mem=rng.choice([1, 4, 18]) * G)
                end = rng.choice([float("inf"), 600.0, 3600.0])
                running_new.append(Running(job_id, Allocation(**kwargs), end))
                running_old.append(self.old_scheduler.Running(job_id, self.old_scheduler.Allocation(**kwargs), end))
            external = {c for c in new_cfg.node_cpus() if c not in used and rng.random() < 0.15}
            queue_args = []
            for job_id in range(1, rng.randint(1, 7)):
                part = new_cfg.partitions[rng.choice(names)]
                gpus = 1 if part.gpus and rng.random() < 0.3 else 0
                queue_args.append(dict(
                    job_id=job_id, partition=part.name, ncpus=rng.randint(1, len(part.cpus)), gpus=gpus,
                    gpu_mem=rng.choice([None, 3 * G]) if gpus else None, mem=rng.choice([1, 4, 18]) * G,
                    time_limit=rng.choice([None, 300, 3600]), blocked=rng.random() < 0.05,
                ))
            new_pool = build_pool(new_cfg, running_new, external_cpus=external)
            old_pool = self.old_scheduler.build_pool(old_cfg, running_old, external_cpus=external)
            new = plan_node(new_cfg, new_pool, running_new, [Request(**a) for a in queue_args], 0.0,
                            owner_of={c: r.job_id for r in running_new for c in r.allocation.cpus})
            old = self.old_scheduler.plan(old_cfg, old_pool, running_old,
                                          [self.old_scheduler.Request(**a) for a in queue_args], 0.0)
            with self.subTest(trial=trial):
                self.assertEqual([(r.job_id, a.cpus, a.gpus, a.gpu_claim, a.mem) for r, a in new.starts],
                                 [(r.job_id, a.cpus, a.gpus, a.gpu_claim, a.mem) for r, a in old.starts])
                self.assertEqual(new.reservation, old.reservation)
                self.assertEqual(new.skipped, old.skipped)
                self.assertEqual(new.preemptions, [])

    def test_plan_with_soft_cores_prefers_free_cores_only(self):
        cfg = config()
        pool = build_pool(cfg, [])
        pool.soft_cpus = {0, 1}
        decision = plan(cfg, pool, [], [req(1, "a", 3)], 0.0)
        self.assertEqual(decision.starts[0][1].cpus, (0, 2, 3))  # free 2-3 first, then the lowest soft core


class TelemetryTests(unittest.TestCase):
    borrow = config().borrow

    def test_in_use_cores_rounding(self):
        self.assertEqual(in_use_cores(0.0, 0.1), 0)
        self.assertEqual(in_use_cores(0.09, 0.1), 0)
        self.assertEqual(in_use_cores(0.11, 0.1), 1)
        self.assertEqual(in_use_cores(1.08, 0.1), 1)
        self.assertEqual(in_use_cores(1.2, 0.1), 2)

    def test_external_holder_idle_busy_and_partly_idle(self):
        busy_cpu = {6: 0.9, 7: 0.02}
        group = HolderGroup("1:1", (6, 7), "pid 1 x", usage_idle=1.0, usage_busy=1.0)
        flags = classify_external(group, cpu_busy=busy_cpu, running_borrowed=set(), occupied=set(), borrow=self.borrow)
        self.assertEqual({c: f.external_idle for c, f in flags.items()}, {6: False, 7: True})
        group = HolderGroup("1:1", (6,), "pid 1 x", usage_idle=0.02, usage_busy=0.02)
        flags = classify_external(group, cpu_busy={6: 0.03}, running_borrowed=set(), occupied=set(), borrow=self.borrow)
        self.assertTrue(flags[6].external_idle)
        self.assertFalse(flags[6].external_busy)
        group = HolderGroup("1:1", (6,), "pid 1 x", usage_idle=None, usage_busy=0.0)
        flags = classify_external(group, cpu_busy={6: 0.0}, running_borrowed=set(), occupied=set(), borrow=self.borrow)
        self.assertFalse(flags[6].external_idle)  # not observed for a whole window yet
        # lent to a running borrower: busy above 0.25 cores
        for usage, expected in ((0.2, False), (0.3, True)):
            group = HolderGroup("1:1", (6,), "pid 1 x", usage_idle=0.05, usage_busy=usage)
            flags = classify_external(group, cpu_busy={6: 1.0}, running_borrowed={6}, occupied={6},
                                      borrow=self.borrow)
            self.assertEqual(flags[6].external_busy, expected)
        # a frozen borrower's core is idle again only after the holder's idle window is quiet
        group = HolderGroup("1:1", (6,), "pid 1 x", usage_idle=0.6, usage_busy=0.0)
        flags = classify_external(group, cpu_busy={6: 0.0}, running_borrowed=set(), occupied={6}, borrow=self.borrow)
        self.assertFalse(flags[6].external_idle)
        group = HolderGroup("stale", (6,), "holder gone", usage_idle=None, usage_busy=None, stale=True)
        flags = classify_external(group, cpu_busy={}, running_borrowed=set(), occupied=set(), borrow=self.borrow)
        self.assertTrue(flags[6].stale and flags[6].external and not flags[6].external_idle)

    def test_owner_job_partly_idle_and_needing_cores_back(self):
        flags = classify_owner((4, 5), usage_idle=1.0, usage_busy=1.0, cpu_busy={4: 1.0, 5: 0.0},
                               running_borrowed=set(), occupied=set(), borrow=self.borrow)
        self.assertEqual({c: f.owner_idle for c, f in flags.items()}, {4: False, 5: True})
        flags = classify_owner((4, 5), usage_idle=1.0, usage_busy=1.6, cpu_busy={4: 1.0, 5: 1.0},
                               running_borrowed={5}, occupied={5}, borrow=self.borrow)
        self.assertTrue(flags[5].owner_busy)
        flags = classify_owner((4, 5), usage_idle=1.0, usage_busy=1.1, cpu_busy={4: 1.0, 5: 1.0},
                               running_borrowed={5}, occupied={5}, borrow=self.borrow)
        self.assertFalse(flags[5].owner_busy)
        self.assertTrue(flags[5].owner_idle)

    def test_usage_tracker_tree_deltas_new_members_and_windows(self):
        tracker = UsageTracker(hz=100)

        def table(**ticks):
            return {int(k[1:]): ProcInfo(int(k[1:]), 1, start, value, "R") for k, (start, value) in ticks.items()}

        tracker.observe_tree("t", "id", [10], table(p10=(5, 1000)), 0.0)
        self.assertIsNone(tracker.tree_usage("t", 2.0, 0.0))
        tracker.observe_tree("t", "id", [10], table(p10=(5, 1100)), 1.0)
        tracker.observe_tree("t", "id", [10, 11], table(p10=(5, 1200), p11=(9, 50)), 2.0)
        self.assertAlmostEqual(tracker.tree_usage("t", 2.0, 2.0), (200 + 50) / 100 / 2.0)
        self.assertAlmostEqual(tracker.tree_usage("t", 1.0, 2.0), 1.5)
        # pid 10 reused by another process (new start time): counted as a new member
        tracker.observe_tree("t", "id", [10], table(p10=(77, 30)), 3.0)
        self.assertAlmostEqual(tracker.tree_usage("t", 1.0, 3.0), 0.3)
        self.assertIsNone(tracker.tree_usage("t", 10.0, 3.0))
        self.assertIsNotNone(tracker.tree_usage("t", 10.0, 3.0, full=False))
        # a new identity restarts the series
        tracker.observe_tree("t", "other", [10], table(p10=(77, 40)), 4.0)
        self.assertIsNone(tracker.tree_usage("t", 0.5, 4.0))

    def test_cpu_busy_from_proc_stat(self):
        first = "cpu  1 2 3 4\ncpu0 100 0 100 800 0 0 0 0 0 0\ncpu1 0 0 0 1000 0 0 0 0 0 0\n"
        second = "cpu  1 2 3 4\ncpu0 150 0 150 900 0 0 0 0 0 0\ncpu1 0 0 0 1200 0 0 0 0 0 0\n"
        self.assertEqual(read_cpu_times(first)[0], (200, 1000))
        tracker = UsageTracker(hz=100)
        tracker.observe_cpus(read_cpu_times(first), 0.0)
        tracker.observe_cpus(read_cpu_times(second), 2.0)
        self.assertAlmostEqual(tracker.cpu_busy(0, 2.0, 2.0), 0.5)
        self.assertAlmostEqual(tracker.cpu_busy(1, 2.0, 2.0), 0.0)

    def test_parse_stat_and_descendants(self):
        info = parse_stat(42, "42 (a b) c) S 7 42 42 0 -1 0 0 0 0 0 11 22 3 4 20 0 1 0 999 0 0")
        self.assertEqual((info.ppid, info.ticks, info.start, info.state), (7, 40, 999, "S"))
        table = {pid: ProcInfo(pid, ppid, 1, 0, "S") for pid, ppid in ((1, 0), (2, 1), (3, 2), (4, 1), (5, 99))}
        self.assertEqual(descendants(table, [2]), {2, 3})
        self.assertEqual(descendants(table, [1]), {1, 2, 3, 4})


class ConfigTests(unittest.TestCase):
    def test_defaults_and_zenobia_policy(self):
        cfg = config()
        part = cfg.partitions["a"]
        self.assertEqual((part.lend, part.lend_leased_idle, part.timing_sensitive), (True, False, False))
        self.assertEqual(cfg.borrow.slice, "simqborrow.slice")
        zenobia = load_node_config(ZENOBIA_CONFIG)
        for name in ("spheral", "raikou-mfem", "powerflow"):
            item = zenobia.partitions[name]
            self.assertEqual((item.lend, item.lend_leased_idle, item.timing_sensitive), (True, False, False), name)
        self.assertEqual(zenobia.borrow.slice_cpu_weight, "idle")
        self.assertEqual(zenobia.borrow.unit_cpu_weight, 1)
        self.assertEqual(zenobia.siblings(6), (14,))
        self.assertEqual(zenobia.locks.borrow_lock_paths(6),
                         ("/home/treefor/projects/.workspace/outputs/pic-profile/borrow-core6.lock",))

    def test_invalid_policy_and_borrow_settings_are_rejected(self):
        for policies, borrow in (
            ({"a": {"lend": "yes"}}, None),
            ({"a": {"timing_sensitive": 1}}, None),
            (None, {"slice": "simq-borrow.slice"}),
            (None, {"slice": "simq.service"}),
            (None, {"slice_cpu_weight": "low"}),
            (None, {"default_on_preempt": "kill"}),
            (None, {"freeze_max_s": 0}),
            (None, {"unknown": 1}),
        ):
            with self.subTest(policies=policies, borrow=borrow):
                with self.assertRaises(ConfigError):
                    config(policies, borrow=borrow)


if __name__ == "__main__":
    unittest.main()


class NodeTelemetryTests(unittest.TestCase):
    """The whole sampling pipeline with a fake clock and fake /proc data."""

    def setUp(self):
        import tempfile

        from simqueue.locks import LockHolder

        self.temporary = tempfile.TemporaryDirectory()
        base = Path(self.temporary.name)
        doc = document({"h": {"lend_leased_idle": True}})
        doc["lock_files"] = {"enabled": True, "core_path": str(base / "core{core}.lock")}
        self.cfg = parse_node_config(doc, path=Path("/dev/null"), check_topology=False)
        for cpu in range(8):
            (base / f"core{cpu}.lock").touch()
        self.lease6 = str(base / "core6.lock")
        self.now = 0.0
        self.ticks = {100: 0, 101: 0}  # 100 holds the lease (any-core.sh); 101 is its child
        self.cpu_ticks = {cpu: [0, 0] for cpu in range(16)}
        self.external = {6: [LockHolder(self.lease6, "FLOCK", "WRITE", 99)]}  # locker pid 99 has exited
        from simqueue.telemetry import NodeTelemetry, file_key

        key = file_key(self.lease6)
        self.telemetry = NodeTelemetry(
            self.cfg,
            clock=lambda: self.now,
            read_table=lambda: {pid: ProcInfo(pid, 1 if pid == 100 else 100, 7, t, "S")
                                for pid, t in self.ticks.items()},
            read_cpus=lambda: {cpu: tuple(v) for cpu, v in self.cpu_ticks.items()},
            fd_holders=lambda keys: {k: ({100, 101} if k == key else set()) for k in keys},
            lock_holders=lambda paths: {},
        )

    def tearDown(self):
        self.temporary.cleanup()

    def advance(self, seconds, *, holder_cores=0.0, cpu6_busy=0.0):
        self.now += seconds
        if 101 in self.ticks:
            self.ticks[101] += int(round(holder_cores * seconds * 100))
        for cpu, state in self.cpu_ticks.items():
            busy = cpu6_busy if cpu == 6 else 0.0
            state[0] += int(round(busy * seconds * 100))
            state[1] += int(round(seconds * 100))

    def sample(self, borrowers=()):
        return self.telemetry.sample(owners=[], external=self.external, borrowers=list(borrowers))

    def test_idle_holder_becomes_lendable_then_busy_preempts(self):
        states = self.sample()
        self.assertTrue(states[6].external and not states[6].external_idle)
        self.assertIn("pid 100", states[6].holder)
        for _ in range(12):  # 12 s of an idle holder; idle window is 10 s
            self.advance(1.0)
            states = self.sample()
        self.assertTrue(states[6].external_idle)
        self.assertFalse(states[6].external_busy)
        decision = run(self.cfg, [req(10, "a", 1, borrow=True, borrow_from=("h",))],
                       running=[owner(1, [0, 1, 2, 3]), owner(2, [7])], cores=states, external=[6])
        self.assertEqual(started(decision)[10].placement, ((6, "h", "leased-idle"),))
        lent = borrower(10, [6], home="a", classes={6: "leased-idle"}, holder_keys={6: states[6].holder_key})
        self.advance(1.0, holder_cores=0.0, cpu6_busy=1.0)  # the borrower itself uses CPU 6
        states = self.sample([lent])
        self.assertFalse(states[6].external_busy)
        self.advance(1.0, holder_cores=1.0, cpu6_busy=1.0)  # the holder's child starts working
        self.advance(1.0, holder_cores=1.0, cpu6_busy=1.0)
        states = self.sample([lent])
        self.assertTrue(states[6].external_busy)
        decision = run(self.cfg, [], borrowers=[lent], cores=states, external=[6])
        self.assertIn("is busy", preempted(decision)[10].reason)

    def test_holder_gone_is_stale_and_free_cores_need_an_idle_cpu(self):
        self.ticks = {}
        states = self.sample()
        self.assertTrue(states[6].stale)
        self.assertIn("holder gone", states[6].holder)
        self.assertTrue(states[0].cpu_busy)  # not measured yet
        self.advance(1.0)
        states = self.sample()
        self.assertFalse(states[0].cpu_busy)
        self.cpu_ticks[0][0] += 300
        self.advance(1.0)
        self.assertTrue(self.sample()[0].cpu_busy)
