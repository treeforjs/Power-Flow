"""Pure scheduling decisions: partitions, concurrency, EASY backfill, GPU budget.

All tests use an explicit (fake) clock value; nothing is launched.
"""

from __future__ import annotations

import json
import unittest
from pathlib import Path

from simqueue.nodeconfig import (
    ConfigError,
    load_node_config,
    parse_node_config,
    parse_size,
    parse_time_limit,
)
from simqueue.scheduler import (
    INFINITY,
    Allocation,
    Request,
    Running,
    build_pool,
    earliest_start,
    plan,
)


ZENOBIA_CONFIG = Path(__file__).resolve().parents[2] / "config" / "zenobia-node.json"
GPU = "GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9"
G = 1024**3


def zenobia():
    config = load_node_config(ZENOBIA_CONFIG)
    assert config is not None
    return config


def small_config(**extra):
    document = {
        "schema": "simqueue-node-config-v1",
        "node": "test",
        "default_partition": "a",
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
    return parse_node_config(document, path=Path("/dev/null"), check_topology=False)


def req(job_id, partition="a", ncpus=1, *, minutes=60, gpus=0, gpu_mem=None, mem=G, cpus=None, blocked=False):
    return Request(
        job_id=job_id, partition=partition, ncpus=ncpus, cpu_list=cpus, gpus=gpus,
        gpu_mem=gpu_mem, mem=mem, time_limit=None if minutes is None else minutes * 60, blocked=blocked,
    )


def started_ids(decision):
    return [request.job_id for request, _ in decision.starts]


class ParsingTests(unittest.TestCase):
    def test_slurm_time_formats(self):
        self.assertEqual(parse_time_limit("30"), 1800)
        self.assertEqual(parse_time_limit("1:30"), 90)
        self.assertEqual(parse_time_limit("2:00:00"), 7200)
        self.assertEqual(parse_time_limit("1-12"), 36 * 3600)
        self.assertEqual(parse_time_limit("1-00:30"), 86400 + 1800)
        self.assertEqual(parse_time_limit("2-00:00:00"), 172800)
        self.assertIsNone(parse_time_limit("UNLIMITED"))
        for bad in ("", "0", "1:2:3:4", "x", "1-x"):
            with self.assertRaises(ValueError):
                parse_time_limit(bad)

    def test_slurm_memory_sizes(self):
        self.assertEqual(parse_size("18G"), 18 * G)
        self.assertEqual(parse_size("18GiB"), 18 * G)
        self.assertEqual(parse_size("512M"), 512 * 1024**2)
        self.assertEqual(parse_size("4096"), 4096 * 1024**2)  # bare number = MiB, as sbatch --mem
        self.assertEqual(parse_size("1GB"), 10**9)
        self.assertEqual(parse_size("100B"), 100)
        with self.assertRaises(ValueError):
            parse_size("4X")

    def test_zenobia_config_partitions(self):
        config = zenobia()
        self.assertEqual(config.default_partition, "spheral")
        self.assertEqual(config.partitions["spheral"].cpus, (0, 1, 2, 3))
        self.assertEqual(config.partitions["spheral"].gpus, ())
        self.assertEqual(config.partitions["raikou-mfem"].cpus, (4, 5))
        self.assertEqual(config.partitions["powerflow"].cpus, (6, 7))
        self.assertEqual(config.partitions["raikou-mfem"].gpus, (GPU,))
        self.assertEqual(config.partitions["powerflow"].gpus, (GPU,))
        self.assertEqual(config.partitions["spheral"].default_mem, 18 * G)
        self.assertEqual(config.partitions["spheral"].resolve_default_cpus(), 4)
        self.assertTrue(config.locks.enabled)

    def test_zenobia_lease_paths_match_existing_lock_files(self):
        config = zenobia()
        base = "/home/treefor/projects/.workspace/outputs/pic-profile/"
        self.assertEqual(config.locks.core_lock_paths(6), (base + "testing-core.lock",))
        self.assertEqual(config.locks.core_lock_paths(7), (base + "testing-core7.lock",))
        self.assertEqual(config.locks.core_lock_paths(0), (base + "testing-core0.lock",))
        self.assertEqual(
            config.locks.gpu_lock_paths(config.gpus[GPU]),
            (base + f"testing-gpu-{GPU}.lock", base + "testing-gpu.lock"),
        )

    def test_config_rejects_smt_siblings_and_bad_values(self):
        document = json.loads(ZENOBIA_CONFIG.read_text())
        document["partitions"][2]["cpus"] = "6,14"
        with self.assertRaises(ConfigError):
            parse_node_config(document, path=ZENOBIA_CONFIG)
        document = json.loads(ZENOBIA_CONFIG.read_text())
        document["partitions"][0]["max_mem"] = "64G"
        with self.assertRaises(ConfigError):
            parse_node_config(document, path=ZENOBIA_CONFIG)
        document = json.loads(ZENOBIA_CONFIG.read_text())
        document["default_partition"] = "nope"
        with self.assertRaises(ConfigError):
            parse_node_config(document, path=ZENOBIA_CONFIG)
        document = json.loads(ZENOBIA_CONFIG.read_text())
        document["partitions"][1]["gpus"] = ["GPU-unknown"]
        with self.assertRaises(ConfigError):
            parse_node_config(document, path=ZENOBIA_CONFIG)


class PlanningTests(unittest.TestCase):
    now = 1_000_000.0

    def plan(self, config, queue, running=(), **pool_args):
        running = list(running)
        pool = build_pool(config, running, **pool_args)
        return plan(config, pool, running, queue, self.now)

    def test_partitions_run_concurrently_and_fifo_within_partition(self):
        config = zenobia()
        queue = [
            req(1, "spheral", 4, mem=18 * G),
            req(2, "powerflow", 1, mem=4 * G),
            req(3, "powerflow", 1, mem=4 * G),
            req(4, "powerflow", 1, mem=4 * G),
            req(5, "raikou-mfem", 2, mem=4 * G),
        ]
        decision = self.plan(config, queue)
        # 18 + 4 + 4 GiB are allocated; job 5's 4 GiB would exceed the 28 GiB
        # schedulable memory, so it waits even though cores 4-5 are free.
        self.assertEqual(started_ids(decision), [1, 2, 3])
        self.assertEqual(decision.skipped[5], "resources busy")
        queue[4] = req(5, "raikou-mfem", 2, mem=2 * G)
        decision = self.plan(config, queue)
        self.assertEqual(started_ids(decision), [1, 2, 3, 5])
        allocations = {r.job_id: a for r, a in decision.starts}
        self.assertEqual(allocations[1].cpus, (0, 1, 2, 3))
        self.assertEqual(allocations[2].cpus, (6,))
        self.assertEqual(allocations[3].cpus, (7,))
        self.assertEqual(allocations[5].cpus, (4, 5))
        self.assertEqual(decision.skipped[4], "resources busy")

    def test_running_allocations_and_outside_leases_are_never_reused(self):
        config = zenobia()
        running = [Running(9, Allocation(cpus=(6,), mem=4 * G), self.now + 600)]
        decision = self.plan(config, [req(1, "powerflow", 1), req(2, "powerflow", 1)], running,
                             external_cpus={7})
        self.assertEqual(started_ids(decision), [])
        decision = self.plan(config, [req(1, "powerflow", 1)], running)
        self.assertEqual([a.cpus for _, a in decision.starts], [(7,)])

    def test_priority_orders_before_fifo(self):
        config = small_config()
        # The caller passes the queue in dispatch order (priority DESC, FIFO).
        queue = [req(2, "a", 4), req(1, "a", 4)]
        decision = self.plan(config, queue)
        self.assertEqual(started_ids(decision), [2])

    def test_explicit_cpu_list(self):
        config = small_config()
        decision = self.plan(config, [req(1, "a", 2, cpus=(1, 3))])
        self.assertEqual(decision.starts[0][1].cpus, (1, 3))
        running = [Running(5, Allocation(cpus=(3,)), self.now + 60)]
        decision = self.plan(config, [req(1, "a", 2, cpus=(1, 3)), req(2, "a", 1)], running)
        # Job 1 waits for core 3; job 2 backfills on core 0, which job 1 never needs.
        self.assertEqual(started_ids(decision), [2])
        self.assertEqual(decision.starts[0][1].cpus, (0,))
        self.assertEqual(decision.reservation, (1, self.now + 60))

    def test_infeasible_request_is_skipped_without_blocking(self):
        config = small_config()
        decision = self.plan(config, [req(1, "a", 6), req(2, "a", 1)])
        self.assertEqual(started_ids(decision), [2])
        self.assertIn("4 cores", decision.skipped[1])


class BackfillTests(unittest.TestCase):
    """EASY backfill: a later job may start only if it cannot delay the
    first blocked job, given every job's time limit."""

    now = 2_000_000.0

    def setUp(self):
        self.config = small_config()
        # Two of partition a's four cores are busy until now+30 min.
        self.running = [Running(100, Allocation(cpus=(0, 1), mem=G), self.now + 30 * 60)]

    def decide(self, queue, running=None):
        running = self.running if running is None else running
        pool = build_pool(self.config, running)
        return plan(self.config, pool, list(running), queue, self.now)

    def test_shadow_time_of_blocked_head(self):
        pool = build_pool(self.config, self.running)
        shadow = earliest_start(req(1, "a", 4), self.config.partitions["a"], pool, self.running, self.now)
        self.assertEqual(shadow, self.now + 30 * 60)

    def test_short_job_backfills_before_shadow_time(self):
        decision = self.decide([req(1, "a", 4), req(2, "a", 2, minutes=20)])
        self.assertEqual(started_ids(decision), [2])
        self.assertEqual(decision.reservation, (1, self.now + 30 * 60))

    def test_job_ending_exactly_at_shadow_time_backfills(self):
        decision = self.decide([req(1, "a", 4), req(2, "a", 2, minutes=30)])
        self.assertEqual(started_ids(decision), [2])

    def test_long_job_that_would_delay_the_head_does_not_backfill(self):
        decision = self.decide([req(1, "a", 4), req(2, "a", 1, minutes=31)])
        self.assertEqual(started_ids(decision), [])
        self.assertIn("would delay job 1", decision.skipped[2])

    def test_unlimited_job_never_backfills_into_needed_cores(self):
        decision = self.decide([req(1, "a", 4), req(2, "a", 1, minutes=None)])
        self.assertEqual(started_ids(decision), [])

    def test_long_job_may_use_cores_the_head_does_not_need(self):
        # Head needs 3 of 4 cores at shadow time; one core is spare forever.
        decision = self.decide([req(1, "a", 3), req(2, "a", 1, minutes=600)])
        # Job 1 fits now? 2 free cores < 3 -> blocked; job 2 takes one of the
        # two free cores; at shadow time 3 cores are still free for job 1.
        self.assertEqual(started_ids(decision), [2])
        decision = self.decide([req(1, "a", 3), req(2, "a", 1, minutes=600), req(3, "a", 1, minutes=600)])
        self.assertEqual(started_ids(decision), [2])
        self.assertIn("would delay job 1", decision.skipped[3])

    def test_other_partitions_are_unaffected_by_a_reservation(self):
        decision = self.decide([req(1, "a", 4), req(2, "h", 2, minutes=None)])
        self.assertEqual(started_ids(decision), [2])

    def test_shared_gpu_cannot_backfill_past_a_waiting_exclusive_gpu_job(self):
        running = [
            Running(100, Allocation(cpus=(4,), gpus=("GPU-x",), gpu_claim=4 * G, mem=G), self.now + 600),
        ]
        queue = [
            req(1, "h", 1, gpus=1),  # exclusive: blocked until now+600
            req(2, "h", 1, gpus=1, gpu_mem=2 * G, minutes=5),  # ends before: allowed
            req(3, "g", 1, gpus=1, gpu_mem=2 * G, minutes=60),  # would delay job 1
        ]
        decision = self.decide(queue, running)
        self.assertEqual(started_ids(decision), [2])
        self.assertEqual(decision.reservation, (1, self.now + 600))
        self.assertIn("would delay job 1", decision.skipped[3])

    def test_memory_is_a_counted_resource(self):
        running = [Running(100, Allocation(cpus=(4,), mem=12 * G), self.now + 600)]
        decision = self.decide([req(1, "a", 1, mem=8 * G), req(2, "a", 1, mem=4 * G, minutes=5)], running)
        self.assertEqual(started_ids(decision), [2])
        self.assertEqual(decision.reservation, (1, self.now + 600))

    def test_outside_lease_holder_gives_unknown_shadow_time(self):
        pool = build_pool(self.config, [], external_cpus={0, 1, 2})
        decision = plan(self.config, pool, [], [req(1, "a", 2), req(2, "a", 1, minutes=None)], self.now)
        self.assertEqual(started_ids(decision), [2])
        self.assertEqual(decision.reservation, (1, INFINITY))

    def test_backfill_can_be_disabled(self):
        pool = build_pool(self.config, self.running)
        decision = plan(self.config, pool, list(self.running), [req(1, "a", 4), req(2, "a", 1, minutes=1)],
                        self.now, backfill=False)
        self.assertEqual(started_ids(decision), [])


class GpuBudgetTests(unittest.TestCase):
    now = 3_000_000.0

    def decide(self, queue, running=(), **pool_args):
        config = zenobia()
        running = list(running)
        pool = build_pool(config, running, **pool_args)
        return plan(config, pool, running, queue, self.now)

    def test_default_gpu_use_is_exclusive(self):
        decision = self.decide([req(1, "powerflow", 1, gpus=1), req(2, "raikou-mfem", 1, gpus=1)])
        self.assertEqual(started_ids(decision), [1])
        allocation = decision.starts[0][1]
        self.assertTrue(allocation.gpu_exclusive)
        self.assertEqual(allocation.gpus, (GPU,))
        self.assertEqual(allocation.gpu_claim, 11 * G)

    def test_gpu_memory_claims_share_within_budget(self):
        decision = self.decide([
            req(1, "powerflow", 1, gpus=1, gpu_mem=4 * G),
            req(2, "raikou-mfem", 1, gpus=1, gpu_mem=6 * G),
            req(3, "powerflow", 1, gpus=1, gpu_mem=2 * G),
        ])
        self.assertEqual(started_ids(decision), [1, 2])  # 4 + 6 = 10 <= 11; +2 > 11
        self.assertFalse(decision.starts[0][1].gpu_exclusive)

    def test_exclusive_request_waits_for_shared_users(self):
        running = [Running(7, Allocation(cpus=(6,), gpus=(GPU,), gpu_claim=1 * G, mem=G), self.now + 60)]
        decision = self.decide([req(1, "raikou-mfem", 1, gpus=1)], running)
        self.assertEqual(started_ids(decision), [])

    def test_gpu_lease_held_outside_the_queue_blocks_gpu_jobs_only(self):
        decision = self.decide(
            [req(1, "powerflow", 1, gpus=1, gpu_mem=G), req(2, "powerflow", 1)],
            external_gpus={GPU},
        )
        self.assertEqual(started_ids(decision), [2])

    def test_paused_ensemble_holds_only_its_partition(self):
        decision = self.decide([
            req(1, "powerflow", 1, blocked=True),
            req(2, "powerflow", 1),
            req(3, "spheral", 4, mem=18 * G),
        ])
        self.assertEqual(started_ids(decision), [3])
        self.assertEqual(decision.skipped[2], "held behind a paused ensemble")


if __name__ == "__main__":
    unittest.main()
