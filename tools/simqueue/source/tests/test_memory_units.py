"""SimQueue 0.7.0 memory tracking with real job units.

Jobs allocate known amounts; the worker samples their cgroups, records the
unit's own memory.peak and CPU time at the end, derives --mem auto from
them, holds starts that do not fit actual memory, and handles memory
pressure.  The admission-guard and pressure tests use a simulated node
(memory.test_node_total): its total is small and only the queue's job units
use it, so a few hundred MiB drive it while the machine itself is never
pressed.
"""

from __future__ import annotations

import json
import os
import signal
import time
import unittest
from pathlib import Path

from node_fixture import PYTHON_FROZEN, NodeFixture

MiB = 1024**2
ALLOC = r'''
import os, sys, time
mib, hold, burn = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
gate = sys.argv[4] if len(sys.argv) > 4 else None
end = time.time() + burn
while time.time() < end:
    pass
block = bytearray(mib << 20)
for i in range(0, len(block), 4096):
    block[i] = 1
print(f"allocated {mib} MiB", flush=True)
deadline = time.time() + hold
while time.time() < deadline:
    if gate and os.path.exists(gate):
        break
    time.sleep(0.05)
del block
print("released", flush=True)
'''


class MemoryUnitFixture(NodeFixture):
    memory: dict = {"sample_seconds": 1, "record_seconds": 2}

    def setUp(self) -> None:  # noqa: N802
        super().setUp()
        (self.projects / "alloc.py").write_text(ALLOC, encoding="utf-8")

    def config_document(self):
        document = super().config_document()
        cpus = ",".join(str(c) for c in (list(self.cpus) or [0, 1]))
        document["partitions"] = [{"name": "pa", "cpus": cpus, "default_mem": "256M", "max_mem": "2G",
                                   "default_time": "10:00", "max_time": "1:00:00", "default_cpus": 1}]
        document["gpus"] = []
        document["lock_files"].pop("gpu_path", None)
        document["lock_files"].pop("gpu_extra", None)
        document["launch"]["idle_seconds"] = 5
        document["memory"] = dict(self.memory)
        document["borrow"] = {"idle_window_s": 1, "busy_window_s": 1, "cpu_idle_threshold": 1.0,
                              "idle_threshold": 1.0, "busy_threshold": 1.0}
        document["github"] = {"enabled": False}
        return document

    def alloc(self, mib: int, hold: float, burn: float = 0.0, gate: Path | None = None) -> list[str]:
        command = ["--", PYTHON_FROZEN, str(self.projects / "alloc.py"), str(mib), str(hold), str(burn)]
        return command + ([str(gate)] if gate is not None else [])

    def wait_state(self, job_id: int, state: str, timeout: float = 60) -> dict:
        return self.wait_until(lambda: self.job(job_id) if self.state_of(job_id) == state else None,
                               timeout=timeout, message=f"job {job_id} {state}")

    def events(self, job_id: int) -> list[dict]:
        return self.store.job_events(job_id)


class MemoryTrackingTests(MemoryUnitFixture, unittest.TestCase):
    min_cpus = 1

    def test_known_allocation_peak_cpu_time_samples_and_auto(self) -> None:
        self.start_worker()
        job_id = self.submit("-p", "pa", "-n", "memtrack", "--mem", "1G", *self.alloc(256, 5, burn=1.0))
        # while it runs: the sample reaches the database and simq show reads the live cgroup
        self.wait_until(lambda: (self.job(job_id).get("mem_current_bytes") or 0) >= 256 * MiB, timeout=30,
                        message="a recorded sample of the allocation")
        show = self.cli("show", str(job_id)).stdout
        self.assertIn("Memory now", show)
        self.assertIn("read from the unit's cgroup", show)
        partitions = json.loads(self.cli("partitions", "--json", "--sample", "0").stdout)
        row = next(r for r in partitions["node_memory"]["jobs"] if r["job_id"] == job_id)
        self.assertGreaterEqual(row["current"], 256 * MiB)
        self.assertEqual(row["request"], 1 << 30)
        job = self.wait_state(job_id, "succeeded")
        peak = int(job["mem_peak_bytes"])
        self.assertGreaterEqual(peak, 256 * MiB)
        self.assertLess(peak, 256 * MiB + 96 * MiB, "peak = the allocation plus interpreter and wrapper")
        self.assertEqual(job["mem_peak_source"], "unit")
        self.assertGreaterEqual(int(job["cpu_usage_usec"]), 800_000, "1 s of busy loop")
        self.assertIsNone(job["mem_current_bytes"])
        series = json.loads(job["mem_series_json"])
        self.assertGreaterEqual(max(b[2] for b in series["b"]), 256)
        show = self.cli("show", str(job_id)).stdout
        self.assertIn("memory.peak recorded by the unit at its end", show)
        self.assertIn("CPU time", show)
        self.assertIn("Memory history", show)
        self.assertIn(f"/1G", self.cli("list").stdout)
        # --mem auto from that peak: 1.25 x, rounded up to 256M
        auto = json.loads(self.cli("submit", "--json", "-p", "pa", "-n", "memtrack", "--mem", "auto",
                                   "--", "/bin/true").stdout)
        wanted = -(-int(peak * 1.25) // (256 * MiB)) * 256 * MiB
        self.assertEqual(auto["mem_bytes"], max(wanted, 256 * MiB))
        self.assertEqual(auto["mem_auto"]["jobs"], [job_id])
        self.wait_state(auto["job_id"], "succeeded")
        self.assertIn("mem 512M", self.worker_log())

    def test_adopted_job_keeps_being_tracked_across_a_worker_restart(self) -> None:
        worker = self.start_worker()
        job_id = self.submit("-p", "pa", "-n", "adopt-mem", "--mem", "512M", *self.alloc(160, 6))
        self.wait_until(lambda: (self.job(job_id).get("mem_current_bytes") or 0) >= 160 * MiB, timeout=30,
                        message="a recorded sample")
        worker.send_signal(signal.SIGTERM)
        worker.wait(timeout=30)
        self.assertEqual(self.state_of(job_id), "running")
        self.start_worker()
        job = self.wait_state(job_id, "succeeded")
        log = self.worker_log()
        self.assertIn(f"adopted job {job_id} ", log)
        self.assertGreaterEqual(int(job["mem_peak_bytes"]), 160 * MiB)
        self.assertEqual(job["mem_peak_source"], "unit")

    def test_out_of_memory_job_records_its_peak_at_the_limit(self) -> None:
        from node_fixture import real_oom_allowed

        if not real_oom_allowed():
            self.skipTest("real cgroup OOM: only with SIMQ_TEST_REAL_OOM=1, never on zenobia (test_oom_safety.py "
                          "covers the OOM path with a mocked unit)")
        self.start_worker()
        job_id = self.submit("-p", "pa", "-n", "oom", "--mem", "128M", *self.alloc(400, 2))
        job = self.wait_state(job_id, "failed")
        self.assertEqual(job["result"], "oom-kill")
        self.assertGreaterEqual(int(job["mem_peak_bytes"]), 120 * MiB)
        auto = json.loads(self.cli("submit", "--json", "-p", "pa", "-n", "oom", "--mem", "auto", "--test-only",
                                   "--", "/bin/true").stdout)
        self.assertEqual(auto["request"]["mem_bytes"], 256 * MiB)  # 1.25 x the 128M limit, at least min_mem


class AdmissionGuardTests(MemoryUnitFixture, unittest.TestCase):
    memory = {"sample_seconds": 1, "record_seconds": 2, "test_node_total": "1G", "admission_margin": "128M",
              "ramp_seconds": 0, "high_water_available": "64M", "high_water_percent": 0,
              "low_water_available": "64M", "low_water_percent": 0}

    def test_start_held_until_actual_memory_is_free(self) -> None:
        self.start_worker()
        gate = self.projects / "release-a"
        a = self.submit("-p", "pa", "-n", "big", "--mem", "800M", *self.alloc(600, 60, gate=gate))
        self.wait_until(lambda: (self.job(a).get("mem_current_bytes") or 0) >= 600 * MiB, timeout=30,
                        message="job a allocated")
        # reservations allow b (800M + 400M of 8G schedulable); actual memory does not
        b = self.submit("-p", "pa", "-n", "small", "--mem", "400M", "--", "/bin/sleep", "1")
        # 0.7.1: the planner itself holds it (actual memory is a planned resource)
        self.wait_until(lambda: "memory busy: needs 400M" in str(self.job(b).get("message") or ""), timeout=20,
                        message="b held for memory")
        self.assertEqual(self.state_of(b), "queued")
        listing = self.cli("list").stdout
        self.assertIn("memory busy: needs 400M", listing)
        gate.write_text("go")
        self.wait_state(a, "succeeded")
        self.wait_state(b, "succeeded", timeout=30)


class HeadOfLineTests(MemoryUnitFixture, unittest.TestCase):
    """0.7.1 regression (seen live on 0.7.0): a head job the memory guard can
    never admit must not hold the queue; small jobs behind it run."""

    memory = {"sample_seconds": 1, "record_seconds": 2, "test_node_total": "1G", "admission_margin": "128M",
              "ramp_seconds": 0, "high_water_available": "64M", "high_water_percent": 0,
              "low_water_available": "64M", "low_water_percent": 0}

    def test_unadmittable_head_job_is_parked_and_small_jobs_run(self) -> None:
        cores = ",".join(str(c) for c in self.cpus)
        result = self.cli("submit", "--json", "-p", "pa", "--cpus", cores, "--mem", "950M", "--priority", "10",
                          "-n", "too-big", "--", "/bin/sleep", "1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("exceeds what testnode can admit now (~896M", result.stderr)
        head = int(json.loads(result.stdout)["job_id"])
        smalls = [self.submit("-p", "pa", "-c", "1", "--mem", "64M", "-t", "5:00", "-n", f"small-{i}",
                              "--", "/bin/sleep", "1") for i in range(3)]
        self.start_worker()
        for job_id in smalls:
            self.wait_state(job_id, "succeeded", timeout=60)
        self.assertEqual(self.state_of(head), "queued")
        self.assertIn("parked: mem 950M exceeds admissible ~896M on testnode", self.job(head)["message"])
        self.assertIn("parked: mem 950M", self.cli("list").stdout)

    def test_head_job_waits_for_memory_running_jobs_free_while_others_backfill(self) -> None:
        self.start_worker()
        gate = self.projects / "release"
        big = self.submit("-p", "pa", "-c", "1", "--mem", "700M", "-n", "holder", *self.alloc(540, 120, gate=gate))
        self.wait_until(lambda: (self.job(big).get("mem_current_bytes") or 0) >= 540 * MiB, timeout=30,
                        message="holder allocated")
        head = self.submit("-p", "pa", "-c", "1", "--mem", "500M", "--priority", "10", "-n", "head",
                           "--", "/bin/sleep", "1")
        self.wait_until(lambda: "memory busy: needs 500M" in str(self.job(head).get("message") or ""),
                        timeout=20, message="head waits for memory")
        small = self.submit("-p", "pa", "-c", "1", "--mem", "64M", "-t", "1:00", "-n", "small",
                            "--", "/bin/sleep", "1")
        self.wait_state(small, "succeeded", timeout=60)  # backfills: it ends before the holder's limit
        gate.write_text("go")
        self.wait_state(big, "succeeded", timeout=60)
        self.wait_state(head, "succeeded", timeout=60)


class MemoryPressureTests(MemoryUnitFixture, unittest.TestCase):
    memory = {"sample_seconds": 1, "record_seconds": 2, "test_node_total": "1G", "admission_margin": "0",
              "ramp_seconds": 0, "high_water_available": "300M", "high_water_percent": 0,
              "low_water_available": "500M", "low_water_percent": 0, "pressure_escalate_seconds": 600}

    def start_borrower(self, policy: str) -> int:
        job_id = self.submit("--borrow", "--on-preempt", policy, "-p", "pa", "-c", "1", "--mem", "64M",
                             "-n", f"borrow-{policy}", "--", "/bin/sleep", "40")
        self.wait_state(job_id, "running", timeout=40)
        return job_id

    def spike(self, gate: Path) -> int:
        # 750M of a 900M limit: far from an OOM (never drive a unit into one), and
        # 1G - 750M - overhead is below the 300M high-water mark
        owner = self.submit("-p", "pa", "-c", "1", "-n", "grower", "--mem", "900M", *self.alloc(750, 60, gate=gate))
        self.wait_until(lambda: "memory pressure: MemAvailable" in self.worker_log(), timeout=30,
                        message="memory pressure")
        return owner

    def test_pressure_stops_dispatch_freezes_a_borrower_and_spares_the_owner(self) -> None:
        if len(self.cpus) < 2:
            self.skipTest("needs two CPUs")
        self.start_worker()
        borrower = self.start_borrower("freeze")
        gate = self.projects / "release-owner"
        owner = self.spike(gate)
        self.wait_until(lambda: self.job(borrower).get("frozen_at"), timeout=20, message="borrower frozen")
        events = self.events(borrower)
        self.assertTrue(any("memory pressure" in json.dumps(e["detail"]) for e in events), events)
        owner_events = [e for e in self.events(owner) if e["kind"] == "memory-pressure"]
        self.assertTrue(owner_events and owner_events[0]["detail"]["culprit"], self.events(owner))
        self.assertIn("owner jobs are not stopped", owner_events[0]["detail"]["action"])
        queued = self.submit("-p", "pa", "-c", "1", "-n", "waits", "--mem", "64M", "--", "/bin/true")
        self.wait_until(lambda: "memory pressure" in str(self.job(queued).get("message") or ""), timeout=20,
                        message="dispatch stopped")
        self.assertEqual(self.state_of(queued), "queued")
        self.assertEqual(self.state_of(owner), "running", "an owner job is never stopped for memory pressure")
        partitions = self.cli("partitions").stdout
        self.assertIn("MEMORY PRESSURE", partitions)
        meta = json.loads(self.store.get_meta("memory_pressure"))
        self.assertEqual(meta["culprit"]["job_id"], owner)
        gate.write_text("go")
        self.wait_state(owner, "succeeded")
        self.wait_until(lambda: "memory pressure over" in self.worker_log(), timeout=20, message="recovery")
        self.wait_state(queued, "succeeded", timeout=30)
        self.assertIsNone(self.store.get_meta("memory_pressure"))
        self.cli("cancel", str(borrower))
        self.wait_until(lambda: self.state_of(borrower) in ("canceled", "succeeded"), timeout=30,
                        message="borrower ended")

    def test_pressure_requeues_a_requeue_policy_borrower(self) -> None:
        if len(self.cpus) < 2:
            self.skipTest("needs two CPUs")
        self.start_worker()
        borrower = self.start_borrower("requeue")
        gate = self.projects / "release-owner"
        owner = self.spike(gate)
        self.wait_until(lambda: int(self.job(borrower).get("preempt_count") or 0) >= 1, timeout=20,
                        message="borrower preempted")
        self.wait_until(lambda: self.state_of(borrower) == "queued", timeout=30, message="borrower requeued")
        self.assertEqual(self.state_of(owner), "running")
        gate.write_text("go")
        self.wait_state(owner, "succeeded")
        self.wait_state(borrower, "running", timeout=40)  # back once the pressure is over
        self.cli("cancel", str(borrower))


class SamplingCostTests(MemoryUnitFixture, unittest.TestCase):
    """CPU cost of the worker while it samples running jobs every second."""

    min_cpus = 1
    memory = {"sample_seconds": 1, "record_seconds": 10}

    @staticmethod
    def cpu_seconds(pid: int) -> float:
        fields = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()
        ticks = os.sysconf("SC_CLK_TCK")
        return sum(int(fields[i]) for i in (11, 12, 13, 14)) / ticks  # utime stime cutime cstime

    def test_sampling_loop_cost(self) -> None:
        worker = self.start_worker()
        jobs = [self.submit("-p", "pa", "-c", "1", "-n", f"sleeper-{i}", "--mem", "64M", "--", "/bin/sleep", "60")
                for i in range(len(self.cpus))]
        for job_id in jobs:
            self.wait_state(job_id, "running")
        time.sleep(3)  # past the first passes
        before, started = self.cpu_seconds(worker.pid), time.monotonic()
        time.sleep(20)
        used, span = self.cpu_seconds(worker.pid) - before, time.monotonic() - started
        rss = int(next(line.split()[1] for line in Path(f"/proc/{worker.pid}/status").read_text().splitlines()
                       if line.startswith("VmRSS:"))) * 1024
        record = {"jobs_sampled": len(jobs), "sample_seconds": 1, "window_s": round(span, 2),
                  "cpu_s": round(used, 4), "cpu_s_per_hour": round(used / span * 3600, 2),
                  "percent_of_one_core": round(100 * used / span, 4), "rss_bytes": rss}
        receipts = os.environ.get("SIMQ_TEST_RECEIPTS")
        if receipts:
            Path(receipts).mkdir(parents=True, exist_ok=True)
            (Path(receipts) / "sampling-loop-cost.json").write_text(json.dumps(record, indent=2) + "\n")
        print(f"\nsampling loop: {record}", flush=True)
        # one wake-up per second with a few file reads: far below 1% of a core
        self.assertLess(used / span, 0.01, record)
        for job_id in jobs:
            self.cli("cancel", str(job_id))


if __name__ == "__main__":
    unittest.main()
