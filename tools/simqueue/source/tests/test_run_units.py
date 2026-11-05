"""0.6.3 with real systemd user units: one pool with aliases, and ``simq run``.

Real units run only on the leased test CPUs (``SIMQ_TEST_CPUS``, set by
run-tests.py), with temporary lease files, a per-test unit prefix and a
per-test borrow slice.  The test configuration is one pool partition of the
test CPUs with the aliases ``pa`` and ``pb`` (the two partitions of the
0.6.2-style test configuration, whose cores they prefer).
"""

from __future__ import annotations

import json
import signal
import subprocess
import time
import unittest

from node_fixture import FAKE_GPU, PYTHON_FROZEN

from test_borrow_units import BorrowFixture


PRINT_AND_EXIT = "import sys; print('RUN-OUTPUT', sys.argv[1], flush=True); raise SystemExit(int(sys.argv[1]))"
SLEEPER = "import time; print('RUN-STARTED', flush=True); time.sleep(120)"


class PoolFixture(BorrowFixture):
    def two_partition_document(self):
        return super().config_document()

    def config_document(self):
        document = super().config_document()
        cpus = sorted(set(self.part_a) | set(self.part_b))
        document["default_partition"] = "pool"
        document["partitions"] = [{
            "name": "pool", "cpus": ",".join(map(str, cpus)), "gpus": [FAKE_GPU], "default_cpus": 1,
            "default_mem": "256M", "max_mem": "2G", "default_time": "10:00", "max_time": "1:00:00",
            "aliases": [
                {"name": "pa", "prefer_cpus": list(self.part_a)},
                {"name": "pb", "default_mem": "200M", "prefer_cpus": list(self.part_b)},
            ],
        }]
        return document

    def run_command(self, *arguments: str, timeout: float = 120) -> subprocess.CompletedProcess[str]:
        return self.cli("run", *arguments, timeout=timeout)

    def run_popen(self, *arguments: str) -> subprocess.Popen:
        proc = subprocess.Popen(
            [PYTHON_FROZEN, str(self.launcher), "run", *arguments], cwd=self.projects, env=self.environment,
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True,
        )
        self.holders.append(proc)
        return proc

    def jobs_named(self, name: str) -> list[dict]:
        return [job for job in self.store.list_jobs(include_all=True, limit=1000) if job["name"] == name]


class AliasUnitTests(PoolFixture, unittest.TestCase):
    def test_jobs_recorded_under_old_partition_names_are_adopted_and_scheduled_in_the_pool(self) -> None:
        # 0.6.2-style configuration: pa and pb are partitions of their own
        self.write_config(self.two_partition_document())
        worker = self.start_worker()
        long_b = self.submit("-p", "pb", "-c", str(len(self.part_b)), "-n", "long-b", *self.probe_command(8.0))
        self.wait_state(long_b, "running")
        self.kill9(worker)
        first_a = self.submit("-p", "pa", "-c", "1", "-n", "first-a", *self.probe_command(1.0))
        then_b = self.submit("-p", "pb", "-c", "1", "-n", "then-b", *self.probe_command(1.0))
        self.assertEqual([self.job(j)["partition"] for j in (long_b, first_a, then_b)], ["pb", "pa", "pb"])
        # the pool replaces both partitions; the old names are aliases now
        self.write_config()
        self.start_worker()
        self.wait_until(lambda: f"adopted job {long_b} " in self.worker_log(), message="adoption")
        for job_id in (long_b, first_a, then_b):
            self.wait_state(job_id, "succeeded", timeout=60)
        records = self.read_intervals()
        # then-b was recorded for pb, whose cores long-b held: in the pool it runs
        # on another core while long-b still runs
        self.assertLess(records[then_b]["start"], records[long_b]["end"])
        self.assertTrue(set(records[then_b]["affinity"]) <= set(self.part_a), records[then_b])
        self.assert_no_cpu_overlap()
        self.assertEqual([self.job(j)["partition"] for j in (long_b, first_a, then_b)], ["pb", "pa", "pb"])
        self.assertIn("partition: pb (pool pool)", self.log_text(then_b))
        self.assertIn("aliases pa,pb", self.worker_log())

    def test_borrow_job_under_an_alias_and_an_owner_job_that_needs_the_whole_pool(self) -> None:
        self.start_worker()
        out = self.projects / "borrow-info.json"
        borrow_id = self.submit("--borrow", "-p", "pb", "-c", "1", "--borrow-from", "pa", "--mem", "256M",
                                "-n", "pr-tests", *self.py("info.py", str(out), "6"))
        job = self.wait_state(borrow_id, "running")
        self.wait_until(out.exists, message="borrow job info")
        info = json.loads(out.read_text())
        self.assertEqual(info["env"]["SIMQ_PARTITION"], "pb")
        self.assertEqual(info["env"]["SIMQ_BORROW"], "1")
        self.assertIn(self.borrow_slice, info["cgroup"])
        placement = json.loads(job["placement_json"])
        self.assertEqual([item["partition"] for item in placement["cpus"]], ["pool"])
        self.assertEqual(job["partition"], "pb")
        # an owner job of every pool core starts at once and preempts the borrower
        width = len(set(self.part_a) | set(self.part_b))
        submitted = time.time()
        owner_id = self.submit("-p", "pa", "-c", str(width), "-n", "pilot", *self.probe_command(2))
        owner = self.wait_state(owner_id, "running", "succeeded", timeout=30)
        from datetime import datetime

        started = datetime.fromisoformat(owner["started_at"].replace("Z", "+00:00")).timestamp()
        self.assertLess(started - submitted, 4.0, "the owner job waited for the borrower")
        self.wait_until(lambda: "requeued" in self.kinds(borrow_id), message="requeue")
        self.wait_state(owner_id, "succeeded")
        self.wait_state(borrow_id, "succeeded", timeout=60)


class RunUnitTests(PoolFixture, unittest.TestCase):
    def test_run_streams_the_log_and_returns_the_job_exit_code(self) -> None:
        self.start_worker()
        result = self.run_command("-p", "pa", "-c", "1", "-n", "run-fail", "--", PYTHON_FROZEN, "-c",
                                  PRINT_AND_EXIT, "3")
        self.assertEqual(result.returncode, 3, result.stderr)
        job = self.jobs_named("run-fail")[0]
        self.assertEqual((job["state"], job["exit_code"]), ("failed", 3))
        self.assertIn("RUN-OUTPUT 3", result.stdout)
        self.assertIn(f"===== simq job {job['id']}: run-fail =====", result.stdout)
        self.assertIn(f"===== simq job {job['id']} finished =====", result.stdout)
        self.assertIn(f"simq: job {job['id']} (run-fail) queued in partition pa (pool pool)", result.stderr)
        self.assertIn("exit status 3", result.stderr)
        self.assertNotIn("simq:", result.stdout)  # messages go to stderr only
        result = self.run_command("-p", "pb", "-n", "run-ok", "--", PYTHON_FROZEN, "-c", PRINT_AND_EXIT, "0")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("RUN-OUTPUT 0", result.stdout)
        self.assertEqual(self.jobs_named("run-ok")[0]["state"], "succeeded")
        self.assertEqual(self.log_text(self.jobs_named("run-ok")[0]["id"]).count("RUN-OUTPUT 0"), 1)

    def test_a_replayed_idempotency_key_reattaches_instead_of_submitting_again(self) -> None:
        self.start_worker()
        arguments = ("-p", "pa", "-c", "1", "-n", "ci-step", "--idempotency-key", "ci:run:7", "--",
                     PYTHON_FROZEN, "-c", PRINT_AND_EXIT, "5")
        first = self.run_command(*arguments)
        self.assertEqual(first.returncode, 5, first.stderr)
        again = self.run_command(*arguments)
        self.assertEqual(again.returncode, 5, again.stderr)
        jobs = self.jobs_named("ci-step")
        self.assertEqual(len(jobs), 1)
        self.assertIn(f"attached to existing job {jobs[0]['id']}", again.stderr)
        self.assertIn("RUN-OUTPUT 5", again.stdout)  # the whole log again
        # reattach to a job that is still running, from a second client
        arguments = ("-p", "pa", "-c", "1", "-n", "ci-long", "--idempotency-key", "ci:run:8", "--json", "--",
                     PYTHON_FROZEN, "-c", "import time; time.sleep(3); print('LONG-DONE')")
        one = self.run_popen(*arguments)
        self.wait_until(lambda: self.jobs_named("ci-long") and self.jobs_named("ci-long")[0]["state"] == "running",
                        message="first client's job running")
        two = self.run_command(*arguments)
        out, err = one.communicate(timeout=60)
        self.assertEqual((one.returncode, two.returncode), (0, 0), err + two.stderr)
        self.assertEqual(len(self.jobs_named("ci-long")), 1)
        self.assertIn("LONG-DONE", out)
        self.assertIn("LONG-DONE", two.stdout)
        payload = json.loads(two.stderr.splitlines()[0])
        self.assertEqual((payload["created"], payload["job_id"]), (False, self.jobs_named("ci-long")[0]["id"]))

    def test_sigint_and_sigterm_cancel_the_job(self) -> None:
        self.start_worker()
        for signum, status in ((signal.SIGINT, 130), (signal.SIGTERM, 143)):
            with self.subTest(signal=signum.name):
                name = f"run-{signum.name.lower()}"
                proc = self.run_popen("-p", "pa", "-c", "1", "-n", name, "--", PYTHON_FROZEN, "-c", SLEEPER)
                self.wait_until(lambda: self.jobs_named(name) and "RUN-STARTED" in self.log_text(
                    self.jobs_named(name)[0]["id"]), message="job running")
                job = self.jobs_named(name)[0]
                proc.send_signal(signum)
                out, err = proc.communicate(timeout=60)
                self.assertEqual(proc.returncode, status, err)
                self.assertIn("RUN-STARTED", out)
                self.assertIn(f"{signum.name}: canceling job {job['id']}", err)
                self.assertEqual(self.state_of(int(job["id"])), "canceled")
                self.assertFalse(self.unit_active(str(self.job(int(job["id"]))["unit_name"])))

    def test_signal_while_queued_cancels_before_launch(self) -> None:
        # no worker: the job stays queued
        proc = self.run_popen("-p", "pa", "-c", "1", "-n", "run-queued", "--", PYTHON_FROZEN, "-c", SLEEPER)
        self.wait_until(lambda: self.jobs_named("run-queued"), message="job queued")
        time.sleep(0.5)
        proc.send_signal(signal.SIGTERM)
        out, err = proc.communicate(timeout=60)
        self.assertEqual(proc.returncode, 143, err)
        self.assertIn("worker autostart is disabled", err)
        job = self.jobs_named("run-queued")[0]
        self.assertEqual((job["state"], job["started_at"]), ("canceled", None))


if __name__ == "__main__":
    unittest.main()
