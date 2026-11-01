"""GPU line (0.6.2) with the real worker, claim units and job units.

Real units run only on the CPUs in SIMQ_TEST_CPUS (run-tests.py), with a
private state directory and private lease files; raw ``flock`` users and
``gpu-lease`` helpers are pinned there too.  No real lease file is touched.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import time
import unittest
from pathlib import Path

from node_fixture import FAKE_GPU, PYTHON_FROZEN, NodeFixture, flock_free

from simqueue import gpuline


REPOSITORY = Path(__file__).resolve().parents[1]
GPULINE = REPOSITORY / "simqueue" / "gpuline.py"


class GpuLineFixture(NodeFixture):
    min_cpus = 1
    warn_s = 600

    def config_document(self):
        document = super().config_document()
        document["gpu_line"] = {"poll_s": 0.05, "claim_cpus": list(self.cpus) or [0],
                                "starvation_warn_s": self.warn_s, "claim_hold_max_s": 60}
        return document

    def setUp(self) -> None:
        super().setUp()
        self.l1 = self.locks / "gpu.lock"
        self.l2 = self.locks / f"gpu-{FAKE_GPU}.lock"
        for path in (self.l1, self.l2):
            path.touch()
        self.events = self.projects / "events.txt"
        self.events.write_text("")
        self.pin = ["taskset", "-c", ",".join(str(c) for c in self.cpus)]
        self.cfg = gpuline.parse_line_config(self.config_document())

    # ----- outside users ------------------------------------------------------

    def body(self, tag: str, hold: float) -> list[str]:
        return ["sh", "-c", f'echo "{tag} start $(date +%s.%N)" >> {self.events}; sleep {hold}; '
                            f'echo "{tag} end $(date +%s.%N)" >> {self.events}']

    def spawn(self, argv: list[str]) -> subprocess.Popen:
        # own session: tearDown kills the whole tree (``flock A flock B cmd``
        # forks), so no blocked inner flock outlives a failed test
        proc = subprocess.Popen(argv, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                start_new_session=True)
        self.holders.append(proc)
        return proc

    def tearDown(self) -> None:
        for proc in self.holders:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except OSError:
                pass
        super().tearDown()

    def raw(self, tag: str, hold: float) -> subprocess.Popen:
        return self.spawn([*self.pin, "flock", str(self.l1), "flock", str(self.l2), *self.body(tag, hold)])

    def raw_loop(self, tag: str, count: int, hold: float) -> subprocess.Popen:
        inner = f'echo "{tag} start $(date +%s.%N)" >> {self.events}; sleep {hold}; ' \
                f'echo "{tag} end $(date +%s.%N)" >> {self.events}'
        script = f"for i in $(seq {count}); do flock {self.l1} flock {self.l2} sh -c '{inner}'; done"
        return self.spawn([*self.pin, "bash", "-c", script])

    def helper(self, tag: str, hold: float, *options: str) -> subprocess.Popen:
        return self.spawn([*self.pin, PYTHON_FROZEN, str(GPULINE), "lease", "--config", str(self.config_path), "-q",
                           *options, "--", *self.body(tag, hold)])

    def event_list(self) -> list[tuple[str, str, float]]:
        found = []
        for line in self.events.read_text().splitlines():
            parts = line.split()
            if len(parts) == 3:
                found.append((parts[0], parts[1], float(parts[2])))
        return found

    def first(self, tag: str, kind: str = "start") -> float:
        return min(t for name, k, t in self.event_list() if name == tag and k == kind)

    def l1_waiters(self) -> list[int]:
        return [w.pid for w in gpuline.path_locks([str(self.l1)])[str(self.l1)].waiters]

    def l1_holders(self) -> list[int]:
        return [h.pid for h in gpuline.path_locks([str(self.l1)])[str(self.l1)].holders]

    # ----- queue side -----------------------------------------------------------

    def claims(self) -> list[dict]:
        root = self.state / "claims"
        found = []
        if root.exists():
            for directory in sorted(root.iterdir()):
                status = gpuline.read_json(directory / "status.json")
                if status is not None:
                    found.append(status)
        return found

    def claim_phase(self) -> str | None:
        active = [c for c in self.claims() if c.get("phase") not in ("handed-off", "released", "failed")]
        return active[-1]["phase"] if active else None

    def active_claim_units(self) -> list[str]:
        result = subprocess.run(["systemctl", "--user", "list-units", f"{self.prefix}-gpuclaim-*", "--no-legend",
                                 "--plain", "--state=active,activating,deactivating"], text=True,
                                stdout=subprocess.PIPE, check=False)
        return [line.split()[0] for line in result.stdout.splitlines() if line.strip()]

    def job_start(self, job_id: int) -> float:
        return float(self.read_intervals()[job_id]["start"])

    def stop_worker(self, worker: subprocess.Popen, sig: int = signal.SIGTERM) -> None:
        worker.send_signal(sig)
        worker.wait(timeout=30)


class FairHandoffTests(GpuLineFixture, unittest.TestCase):
    def test_queued_job_starts_at_the_first_release_despite_raw_and_helper_users(self) -> None:
        r0 = self.raw("r0", 4.0)
        self.wait_until(lambda: self.l1_holders(), 10, "raw holder")
        loops = [self.raw_loop("loopA", 3, 0.8), self.raw_loop("loopB", 3, 0.8)]
        self.wait_until(lambda: len(self.l1_waiters()) >= 2, 10, "back-to-back raw waiters")
        job = self.submit("-p", "pb", "--gpus", "1", *self.probe_command(1.0))
        worker = self.start_worker()
        self.wait_until(lambda: self.claim_phase() == "queued", 30, "claim queued in the kernel")
        helpers = [self.helper("helper0", 0.4)]
        time.sleep(0.3)
        helpers.append(self.helper("helper1", 0.4))
        self.wait_until(lambda: self.state_of(job) == "succeeded", 90, "GPU job done")
        for proc in [r0, *helpers, *loops]:
            proc.wait(timeout=90)
        r0_end = self.first("r0", "end")
        job_start = self.job_start(job)
        later = sorted((t, name) for name, kind, t in self.event_list() if kind == "start" and t > r0_end)
        self.assertLess(job_start, later[0][0], f"job started at {job_start}, others: {later[:4]}")
        self.assertLess(job_start - r0_end, 3.0)
        # after the job: the helpers in start order, then the raw users
        self.assertEqual([name for _, name in later[:2]], ["helper0", "helper1"], later)
        log = self.log_text(job)
        self.assertIn("gpu line:  received", log)
        kinds = [e["kind"] for e in self.store.job_events(job)]
        self.assertIn("gpu-handoff", kinds)
        self.assertIn("moved raw waiter", Path(self.cfg.log_path).read_text())
        self.assertEqual(self.active_claim_units(), [])
        self.stop_worker(worker)

    def test_helper_that_started_before_the_job_goes_first(self) -> None:
        r0 = self.raw("r0", 2.5)
        self.wait_until(lambda: self.l1_holders(), 10, "raw holder")
        early = self.helper("early", 0.5)
        self.wait_until(lambda: any(t["kind"] == "helper" and t["state"] == "queued"
                                    for t in gpuline.load_tickets(self.cfg.ticket_dir)), 10, "helper queued")
        time.sleep(0.2)
        job = self.submit("-p", "pb", "--gpus", "1", *self.probe_command(0.5))
        worker = self.start_worker()
        self.wait_until(lambda: self.state_of(job) == "succeeded", 90, "GPU job done")
        early.wait(timeout=30)
        r0.wait(timeout=30)
        self.assertLess(self.first("early", "end"), self.job_start(job) + 0.05)
        self.stop_worker(worker)

    def test_time_limit_starts_when_the_job_starts(self) -> None:
        r0 = self.raw("r0", 8.0)
        self.wait_until(lambda: self.l1_holders(), 10, "raw holder")
        job = self.submit("-p", "pb", "--gpus", "1", "-t", "0:06", *self.probe_command(4.0))
        worker = self.start_worker()
        self.wait_until(lambda: self.state_of(job) in ("succeeded", "failed"), 90, "GPU job done")
        record = self.job(job)
        self.assertEqual(record["state"], "succeeded", record.get("message"))
        self.assertNotEqual(record.get("result"), "timeout")
        r0.wait(timeout=30)
        self.stop_worker(worker)

    def test_shared_job_gets_a_shared_lease_after_the_handoff(self) -> None:
        r0 = self.raw("r0", 2.0)
        self.wait_until(lambda: self.l1_holders(), 10, "raw holder")
        first = self.submit("-p", "pb", "--gpu-mem", "1G", "--", "/bin/sleep", "6")
        worker = self.start_worker()
        self.wait_until(lambda: self.state_of(first) == "running", 60, "shared job running")
        from simqueue.locks import lock_holders

        path = str(self.l1)
        self.wait_until(lambda: lock_holders([path]).get(path), 10, "lease")
        self.wait_until(lambda: all(h.mode == "READ" for h in lock_holders([path])[path]), 10, "shared lease")
        second = self.submit("-p", "pb", "--gpu-mem", "1G", "--", "/bin/sleep", "1")
        self.wait_until(lambda: self.state_of(second) in ("running", "succeeded"), 60, "second shared job")
        self.wait_until(lambda: self.state_of(first) == "succeeded", 60, "first done")
        r0.wait(timeout=30)
        self.stop_worker(worker)

    def test_cancel_of_a_waiting_job_stops_the_claim(self) -> None:
        r0 = self.raw("r0", 30.0)
        self.wait_until(lambda: self.l1_holders(), 10, "raw holder")
        job = self.submit("-p", "pb", "--gpus", "1", *self.probe_command(0.2))
        worker = self.start_worker()
        self.wait_until(lambda: self.claim_phase() == "queued", 30, "claim queued")
        self.assertIn("waiting for the GPU", str(self.job(job)["message"]))
        self.assertEqual(self.cli("cancel", str(job)).returncode, 0)
        self.wait_until(lambda: not self.active_claim_units(), 30, "claim stopped")
        self.wait_until(lambda: not [t for t in gpuline.load_tickets(self.cfg.ticket_dir) if t["kind"] == "queue"],
                        10, "claim ticket removed")
        self.assertEqual(self.claims()[-1]["phase"], "released")
        self.assertEqual(self.l1_waiters(), [])
        r0.kill()
        self.stop_worker(worker)


class RestartTests(GpuLineFixture, unittest.TestCase):
    def test_worker_restart_while_the_job_waits_and_while_it_runs(self) -> None:
        r0 = self.raw("r0", 6.0)
        self.wait_until(lambda: self.l1_holders(), 10, "raw holder")
        loop = self.raw_loop("loop", 3, 0.5)
        job = self.submit("-p", "pb", "--gpus", "1", *self.probe_command(6.0))
        worker = self.start_worker()
        self.wait_until(lambda: self.claim_phase() == "queued", 30, "claim queued")
        # the worker goes away while the job waits; the claim keeps its place
        self.stop_worker(worker)
        self.wait_until(lambda: self.claim_phase() == "held", 30, "claim holds while the worker is down")
        self.assertEqual(self.state_of(job), "queued")
        time.sleep(0.5)
        self.assertNotIn("loop", [name for name, _, _ in self.event_list()])
        worker = self.start_worker()
        self.wait_until(lambda: self.state_of(job) == "running", 30, "job started by the new worker")
        self.assertIn("adopted GPU claim", self.worker_log())
        # the worker goes away while the job runs (SIGKILL); the job keeps the GPU
        self.wait_until(lambda: "gpu line:  received" in self.log_text(job), 10, "hand-off")
        self.kill9(worker)
        time.sleep(1.0)
        self.assertFalse(flock_free(self.l1))
        self.assertEqual(self.state_of(job), "running")
        worker = self.start_worker()
        self.wait_until(lambda: self.state_of(job) == "succeeded", 60, "job recorded by the third worker")
        self.assertIn(f"adopted job {job} ", self.worker_log())
        loop.wait(timeout=60)
        r0.wait(timeout=30)
        job_end = float(self.read_intervals()[job]["end"])
        self.assertTrue(all(t >= job_end - 0.05 for name, kind, t in self.event_list() if name == "loop"))
        self.stop_worker(worker)


class CoreLeaseTests(GpuLineFixture, unittest.TestCase):
    def test_no_deadlock_between_gpu_users_and_core_leases(self) -> None:
        """The queue's claim and job and a ``gpu-lease --core`` helper never
        deadlock with raw users that take the core lease first or the GPU
        files in reversed order.  Raw users of two different orders can
        deadlock with each other whatever the queue does (``core, l1, l2``
        against ``l2, l1`` or ``l1, l2, core``), so each phase has raw users
        of one order only, next to the queue job and the helper."""
        cpu = self.cpus[0]
        core = self.locks / f"core{cpu}.lock"
        core.touch()
        phases = {
            # core lease first, then the GPU files (the order that deadlocked once)
            "corefirst": ["flock", str(core), "flock", str(self.l1), "flock", str(self.l2)],
            # per-GPU file first
            "reversed": ["flock", str(self.l2), "flock", str(self.l1)],
            # documented order: GPU files, then the core lease
            "gpufirst": ["flock", str(self.l1), "flock", str(self.l2), "flock", str(core)],
        }
        worker = self.start_worker()
        for tag, prefix in phases.items():
            r0 = self.spawn([*self.pin, "flock", str(self.l1), *self.body(f"r0-{tag}", 1.5)])
            self.wait_until(lambda: self.l1_holders(), 10, "raw holder")
            raw = self.spawn([*self.pin, *prefix, *self.body(tag, 0.5)])
            helper = self.helper(f"helper-{tag}", 0.5, "--core", str(cpu))
            job = self.submit("-p", "pb", "--gpus", "1", *self.probe_command(0.5))
            deadline = time.monotonic() + 60
            for proc in (r0, raw, helper):
                proc.wait(timeout=max(1, deadline - time.monotonic()))
            self.wait_until(lambda: self.state_of(job) == "succeeded", max(1, deadline - time.monotonic()),
                            f"job of phase {tag}")
        names = sorted({name for name, kind, _ in self.event_list() if kind == "end"})
        self.assertEqual(names, sorted([*phases, *(f"r0-{t}" for t in phases), *(f"helper-{t}" for t in phases)]))
        self.stop_worker(worker)


class DisplayTests(GpuLineFixture, unittest.TestCase):
    warn_s = 2

    def test_partitions_show_and_published_state(self) -> None:
        r0 = self.raw("r0", 30.0)
        self.wait_until(lambda: self.l1_holders(), 10, "raw holder")
        waiter = self.raw("rawwaiter", 0.1)
        self.wait_until(lambda: len(self.l1_waiters()) == 1, 10, "raw waiter")
        job = self.submit("-p", "pb", "--gpus", "1", *self.probe_command(0.2))
        worker = self.start_worker()
        self.wait_until(lambda: self.claim_phase() == "queued", 30, "claim queued")
        time.sleep(2.5)
        state = gpuline.load_queue_state(self.cfg.state_path)
        entries = state["gpus"][FAKE_GPU]["queue"]
        self.assertEqual([e["job_id"] for e in entries], [job])
        self.assertIsNotNone(entries[0]["since_unix"])
        self.assertTrue(state["gpus"][FAKE_GPU]["warnings"])
        output = self.cli("partitions", "--sample", "0")
        self.assertEqual(output.returncode, 0, output.stderr)
        text = output.stdout
        self.assertIn("holder: raw flock holder pid", text)
        self.assertIn(f"queue job {job}", text)
        self.assertIn("WARNING: queue job", text)
        self.assertIn("raw flock waiters", text)
        self.assertIn("queue claim:", text)
        shown = self.cli("show", str(job)).stdout
        self.assertIn("GPU line", shown)
        self.assertIn("waiting for the GPU", shown)
        listed = self.cli("list").stdout
        self.assertIn("waiting for the GPU", listed)
        status = subprocess.run([PYTHON_FROZEN, str(GPULINE), "lease", "--config", str(self.config_path), "--status"],
                                text=True, stdout=subprocess.PIPE, check=False)
        self.assertIn("holder: raw flock holder pid", status.stdout)
        self.assertIn(f"queue job {job}", status.stdout)
        doctor = self.cli("doctor").stdout
        self.assertIn("GPU line", doctor)
        r0.kill()
        self.wait_until(lambda: self.state_of(job) == "succeeded", 60, "job after the holder ended")
        waiter.wait(timeout=30)
        self.stop_worker(worker)
        state = gpuline.load_queue_state(self.cfg.state_path)
        self.assertEqual(state["heartbeat_unix"], 0)  # withdrawn when the worker stopped


if __name__ == "__main__":
    unittest.main()
