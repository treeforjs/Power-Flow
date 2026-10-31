"""Borrowing (0.6.1): database operations and command-line validation.

No unit is started: the store is exercised directly and the CLI with
SIMQUEUE_NO_AUTOSTART=1.
"""

from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

from node_fixture import PYTHON_FROZEN, NodeFixture

from simqueue.store import QueueStore


RUNTIME_060 = Path("/home/treefor/projects/runtimes/simqueue-0.6.0-20261030")


class BorrowStoreTests(NodeFixture, unittest.TestCase):
    needs_units = False

    def submit_job(self, name: str, *, borrow: bool = False, priority: int = 0, partition: str = "pa") -> dict:
        scheduling = {"partition": partition, "priority": priority, "mem_bytes": 1 << 28, "time_limit_s": 600}
        if borrow:
            scheduling.update({"borrow": 1, "on_preempt": "freeze"})
        return self.store.submit(name=name, cwd=self.projects, argv=["/bin/true"], environment={},
                                 env_override_keys=[], requested_cores=1, scheduling=scheduling)

    def claim(self, job_id: int, cpus: str, **kwargs) -> dict | None:
        return self.store.claim_allocation(
            job_id, cpu_list=cpus, gpu_list="", gpu_claim_bytes=0, gpu_exclusive=False, mem_bytes=1 << 28,
            swap_bytes=0, time_limit_s=600, partition="pa", unit_prefix=self.prefix, job_root=self.state / "jobs",
            lock_paths=[], runtime_dir="/x", busy_cpus=set(), **kwargs)

    def test_borrow_jobs_sort_after_owner_jobs_in_every_view(self) -> None:
        borrow = self.submit_job("borrow", borrow=True, priority=50)
        owner = self.submit_job("owner", priority=0)
        low = self.submit_job("owner-low", priority=-5)
        self.assertEqual([j["id"] for j in self.store.queued_jobs()], [owner["id"], low["id"], borrow["id"]])
        self.assertEqual([j["id"] for j in self.store.node_jobs()], [owner["id"], low["id"], borrow["id"]])
        listed = [j["id"] for j in self.store.list_jobs() if j["state"] == "queued"]
        self.assertEqual(listed, [owner["id"], low["id"], borrow["id"]])
        # promote moves a borrow job to the front of the borrow jobs only
        second = self.submit_job("borrow-2", borrow=True)
        result = self.store.promote(second["id"])
        self.assertEqual(result["queue_position"], 3)
        self.assertEqual([j["id"] for j in self.store.queued_jobs()],
                         [owner["id"], low["id"], second["id"], borrow["id"]])
        result = self.store.demote(owner["id"])
        self.assertEqual(result["queue_position"], 2)
        self.assertEqual([j["id"] for j in self.store.queued_jobs()][:2], [low["id"], owner["id"]])

    def test_claim_overlap_rules(self) -> None:
        owner = self.submit_job("owner")
        borrow = self.submit_job("borrow", borrow=True)
        self.assertIsNotNone(self.claim(borrow["id"], "0", borrow=True, placement={"cpus": []}))
        self.store.mark_unit_running(borrow["id"], pid=1, boot_id="b")
        with self.assertRaises(RuntimeError):
            self.claim(owner["id"], "0")  # a running borrower must be preempted first
        self.assertTrue(self.store.mark_preempted(borrow["id"], mode="freeze", reason="test", by_job=owner["id"]))
        self.assertIsNotNone(self.claim(owner["id"], "0"))
        self.store.mark_unit_running(owner["id"], pid=2, boot_id="b")
        other = self.submit_job("borrow-2", borrow=True)
        with self.assertRaises(RuntimeError):
            # a frozen borrower's cores stay reserved for its resume
            self.claim(other["id"], "0", borrow=True, placement={}, owner_idle_cpus={0})
        self.assertTrue(self.store.mark_preempted(borrow["id"], mode="requeue", reason="x"))
        self.assertTrue(self.store.requeue_preempted(borrow["id"], "requeued"))
        with self.assertRaises(RuntimeError):
            self.claim(other["id"], "0", borrow=True, placement={})  # owner's core, not lent as owner-idle
        self.assertIsNotNone(self.claim(other["id"], "0", borrow=True, placement={}, owner_idle_cpus={0}))

    def test_preemption_resume_requeue_and_events(self) -> None:
        owner = self.submit_job("owner")
        borrow = self.submit_job("borrow", borrow=True)
        self.claim(borrow["id"], "1", borrow=True, placement={"cpus": [{"cpu": 1, "partition": "pa", "class": "free"}]})
        self.store.mark_unit_running(borrow["id"], pid=1, boot_id="b")
        self.assertTrue(self.store.mark_preempted(borrow["id"], mode="freeze", reason="cores needed", by_job=owner["id"]))
        self.assertFalse(self.store.mark_preempted(borrow["id"], mode="freeze", reason="again"))
        job = self.store.require(borrow["id"])
        self.assertTrue(job["frozen_at"])
        self.assertEqual(job["preempt_count"], 1)
        self.assertTrue(self.store.mark_resumed(borrow["id"], {"cpus": [{"cpu": 1, "class": "free"}]}))
        self.assertIsNone(self.store.require(borrow["id"])["frozen_at"])
        self.assertTrue(self.store.mark_preempted(borrow["id"], mode="requeue", reason="memory", by_job=owner["id"],
                                                  escalated=True))
        job = self.store.require(borrow["id"])
        self.assertEqual((job["preempt_pending"], job["preempt_count"]), ("requeue", 2))
        token = job["launch_token"]
        self.assertTrue(self.store.requeue_preempted(borrow["id"], "pending: preempted; requeued"))
        job = self.store.require(borrow["id"])
        self.assertEqual((job["state"], job["cpu_list"], job["preempt_pending"], job["unit_name"]),
                         ("queued", None, None, None))
        self.assertNotEqual(job["launch_token"], token)
        self.assertEqual(job["queue_order"], borrow["queue_order"])
        kinds = [event["kind"] for event in self.store.job_events(borrow["id"])]
        self.assertEqual(kinds, ["started", "preempted", "resumed", "preempted", "requeued"])
        escalation = self.store.job_events(borrow["id"])[3]["detail"]
        self.assertTrue(escalation["escalated"])
        reclaimed = [e["kind"] for e in self.store.job_events(owner["id"])]
        self.assertEqual(reclaimed, ["reclaimed", "reclaimed"])
        # owner jobs cannot be preempted
        self.claim(owner["id"], "0")
        self.assertFalse(self.store.mark_preempted(owner["id"], mode="requeue", reason="x"))

    def test_retry_keeps_the_borrow_request(self) -> None:
        borrow = self.submit_job("borrow", borrow=True)
        self.store.cancel(borrow["id"])
        retried = self.store.retry(borrow["id"])
        self.assertEqual((retried["borrow"], retried["on_preempt"]), (1, "freeze"))


class BorrowCliTests(NodeFixture, unittest.TestCase):
    needs_units = False

    def fails(self, *arguments: str) -> str:
        result = self.cli("submit", *arguments)
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        return result.stderr

    def test_borrow_option_validation(self) -> None:
        self.assertIn("cannot use GPUs", self.fails("--borrow", "-p", "pb", "--gpus", "1", "--", "/bin/true"))
        self.assertIn("--cpus cannot be combined", self.fails("--borrow", "--cpus", str(self.part_a[0]), "--",
                                                              "/bin/true"))
        self.assertIn("--on-preempt needs --borrow", self.fails("--on-preempt", "freeze", "--", "/bin/true"))
        self.assertIn("--exclusive-timing is for owner jobs",
                      self.fails("--borrow", "--exclusive-timing", "--", "/bin/true"))
        self.assertIn("unknown partition", self.fails("--borrow-from", "nope", "--", "/bin/true"))
        total = len(self.part_a) + len(self.part_b)
        self.assertIn(f"between 1 and {total}", self.fails("--borrow", "-c", str(total + 1), "--", "/bin/true"))
        self.assertIn("--max-borrow-cores must be between", self.fails("--borrow", "-c", "2", "--max-borrow-cores",
                                                                       "1", "--", "/bin/true"))

    def test_borrow_submission_fields_and_json(self) -> None:
        total = len(self.part_a) + len(self.part_b)
        result = self.cli("submit", "--json", "--opportunistic", "-p", "pa", "--max-borrow-cores", str(total),
                          "--on-preempt", "freeze", "--borrow-from", "pb", "--", "/bin/true")
        self.assertEqual(result.returncode, 0, result.stderr)
        payload = json.loads(result.stdout)
        self.assertEqual((payload["borrow"], payload["borrow_from"], payload["max_borrow_cores"],
                          payload["on_preempt"], payload["requested_cores"]), (True, ["pb"], total, "freeze", 1))
        job = self.store.require(payload["job_id"])
        self.assertEqual((job["borrow"], json.loads(job["borrow_from"]), job["max_borrow_cores"], job["on_preempt"]),
                         (1, ["pb"], total, "freeze"))
        shown = self.cli("show", str(payload["job_id"])).stdout
        self.assertIn("borrow job", shown)
        self.assertIn("On preempt", shown)
        listing = self.cli("list").stdout
        self.assertIn("~pa", listing)
        self.assertIn("~ borrow job", listing)
        owner = self.cli("submit", "--json", "-p", "pb", "--exclusive-timing", "--", "/bin/true")
        self.assertEqual(json.loads(owner.stdout)["exclusive_timing"], True)

    def test_borrow_from_a_partition_that_does_not_lend_is_refused(self) -> None:
        document = self.config_document()
        for part in document["partitions"]:
            if part["name"] == "pb":
                part["lend"] = False
        self.write_config(document)
        self.assertIn("does not lend", self.fails("--borrow-from", "pb", "-p", "pa", "--", "/bin/true"))
        document["partitions"][1]["lend"] = True
        document["partitions"][1]["timing_sensitive"] = True
        self.write_config(document)
        self.assertIn("timing-sensitive", self.fails("--borrow-from", "pb", "-p", "pa", "--", "/bin/true"))

    def test_test_only_shows_borrow_placement(self) -> None:
        result = self.cli("submit", "--test-only", "--borrow", "-p", "pa", "-c", "1", "--", "/bin/true")
        self.assertEqual(result.returncode, 0, result.stderr)
        payload = json.loads(result.stdout)
        self.assertTrue(payload["would_start_now"])
        self.assertEqual(payload["borrow_placement"][0]["class"], "free")
        self.assertTrue(all(lock["mode"] == "shared" and "borrow-core" in lock["path"]
                            for lock in payload["lease_files"]))

    def test_060_submission_replayed_under_061_returns_the_same_job(self) -> None:
        if not (RUNTIME_060 / "source" / "simq").exists():
            self.skipTest("0.6.0 runtime not present")
        state = self.base / "state-060"
        env = dict(self.environment, SIMQUEUE_STATE_DIR=str(state))
        arguments = ("submit", "--json", "--idempotency-key", "fp:1", "-p", "pa", "-c", "1", "--mem", "512M",
                     "-t", "5:00", "--priority", "2", "-e", "X=1", "--", "/bin/echo", "hi")
        old_launcher = self.make_launcher(RUNTIME_060 / "source")
        old = subprocess.run([PYTHON_FROZEN, str(old_launcher), *arguments], cwd=self.projects, env=env, text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        self.assertEqual(old.returncode, 0, old.stderr)
        new = self.cli(*arguments, env=env)  # migrates that database to schema 8, then replays the key
        self.assertEqual(new.returncode, 0, new.stderr)
        first, replay = json.loads(old.stdout), json.loads(new.stdout)
        self.assertEqual((replay["job_id"], replay["created"]), (first["job_id"], False))
        self.assertEqual(replay["request_sha256"], first["request_sha256"])

    def test_legacy_mode_rejects_borrow_options(self) -> None:
        self.config_path.unlink()
        env = dict(self.environment)
        env.pop("SIMQUEUE_CONFIG")
        for options in (("--borrow",), ("--exclusive-timing",), ("--max-borrow-cores", "2")):
            result = self.cli("submit", *options, "--", "/bin/true", env=env)
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertIn("need a node configuration", result.stderr)

    def test_partitions_shows_policy_and_per_core_states(self) -> None:
        holder = self.hold_lease(self.locks / f"core{self.part_b[0]}.lock")
        text = self.cli("partitions", "--sample", "0.2").stdout
        self.assertIn("POLICY", text)
        self.assertIn("lend", text)
        self.assertIn(f"CPU {self.part_b[0]}: external lease idle", text)
        self.assertIn(f"holder pid {holder.pid}", text)
        status = json.loads(self.cli("partitions", "--json", "--sample", "0").stdout)
        cpus = {c["cpu"]: c for p in status["partitions"] for c in p["cpus"]}
        self.assertTrue(cpus[self.part_b[0]]["state"].startswith("lease held: "))  # 0.6.0 key kept
        self.assertEqual(cpus[self.part_a[0]]["state"], "free")
        self.release_lease(holder)


if __name__ == "__main__":
    unittest.main()
