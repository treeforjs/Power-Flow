"""SimQueue 0.7.1 OOM safety (critical fix).

On 2026-11-05 an OOM group kill of a test job unit killed the user manager:
the manager had just spawned (vfork, sharing its mm) a process into the
unit's cgroup, and __oom_kill_process kills every process sharing the
victim's mm.  Job units therefore never make the manager spawn into their
cgroup after the job started (no ExecStopPost, ExecStartPost, ExecStop),
use OOMPolicy=stop (memory.oom.group 0), and are not collected, so the
worker reads a killed job's result from the unit.  No test here drives a
unit into a real OOM: the OOM path is exercised with a mocked unit state.
"""

from __future__ import annotations

import json
import os
import sqlite3
import subprocess
import sys
import time
import unittest
from contextlib import closing
from pathlib import Path

from node_fixture import PYTHON_FROZEN, NodeFixture

from simqueue.nodeconfig import LaunchConfig
from simqueue.units import FORBIDDEN_JOB_PROPERTIES, UnitState, build_run_command

MiB = 1024**2


class CommandTests(unittest.TestCase):
    def test_job_unit_command_has_no_post_start_spawns(self) -> None:
        command = build_run_command(LaunchConfig(), unit="simq-job-1-abc.service", description="d", cwd="/tmp",
                                    cpus=(0, 1), mem_bytes=MiB * 512, swap_bytes=0, runtime_max_s=60, nice=0,
                                    log_path="/tmp/l", spec_path="/tmp/s", setenv={})
        text = " ".join(command)
        for prop in FORBIDDEN_JOB_PROPERTIES:
            self.assertNotIn(f"--property={prop}=", text, prop)
        self.assertIn("--property=OOMPolicy=stop", command)
        self.assertNotIn("OOMPolicy=kill", text)
        self.assertNotIn("--collect", command)
        self.assertIn("--property=KillMode=control-group", command)
        from simqueue.units import JOBWRAP

        self.assertEqual(command[-5:], ["--", LaunchConfig().wrapper_python, str(JOBWRAP), "run", "/tmp/s"])


class UnitPropertyTests(NodeFixture, unittest.TestCase):
    """A real job unit (it never comes near its memory limit)."""

    min_cpus = 1

    def show(self, unit: str, *props: str) -> dict[str, str]:
        result = subprocess.run(["systemctl", "--user", "show", unit, *[f"--property={p}" for p in props]],
                                text=True, stdout=subprocess.PIPE, check=False)
        return dict(line.split("=", 1) for line in result.stdout.splitlines() if "=" in line)

    def test_running_unit_properties_and_failed_unit_is_recorded_then_unloaded(self) -> None:
        self.start_worker()
        job_id = self.submit("-p", "pa", "-c", "1", "--", "/bin/sh", "-c", "sleep 4; exit 3")
        self.wait_until(lambda: self.state_of(job_id) == "running", message="running")
        unit = str(self.job(job_id)["unit_name"])
        props = self.show(unit, "OOMPolicy", "ExecStopPost", "ExecStartPost", "ExecStop", "ControlGroup")
        self.assertEqual(props["OOMPolicy"], "stop")
        for prop in ("ExecStopPost", "ExecStartPost", "ExecStop"):
            self.assertEqual(props.get(prop, ""), "", prop)
        cgroup = Path("/sys/fs/cgroup") / props["ControlGroup"].lstrip("/")
        self.assertEqual((cgroup / "memory.oom.group").read_text().strip(), "0")
        self.wait_until(lambda: self.state_of(job_id) == "failed", timeout=60, message="failed")
        job = self.job(job_id)
        self.assertEqual(job["exit_code"], 3)
        self.assertIsNotNone(job["mem_peak_bytes"])
        # the worker read the failed unit, then unloaded it
        self.wait_until(lambda: self.show(unit, "LoadState").get("LoadState") == "not-found", timeout=20,
                        message="failed unit reset")


class MockedOomTests(NodeFixture, unittest.TestCase):
    """The OOM path from the unit's state to the database, with no real OOM."""

    needs_units = False

    def test_oom_killed_unit_is_recorded_from_its_properties(self) -> None:
        from simqueue import nodeworker
        from simqueue.nodeconfig import load_node_config

        job = self.store.submit(name="oom", cwd=self.projects, argv=["/bin/true"], environment={},
                                env_override_keys=[], requested_cores=1,
                                scheduling={"partition": "pa", "mem_bytes": 128 * MiB})
        job_dir = self.state / "jobs" / "000001-oom"
        job_dir.mkdir(parents=True)
        unit_name = f"{self.prefix}-{job['id']}-{job['launch_token'][:12]}.service"
        (job_dir / "status.json").write_text(json.dumps({"launch_token": job["launch_token"], "phase": "running",
                                                         "boot_id": None}))
        with closing(sqlite3.connect(self.store.path)) as connection:
            connection.execute("UPDATE jobs SET state = 'running', unit_name = ?, job_dir = ?, cpu_list = '0', "
                               "launch_backend = 'systemd-user', started_at = '2026-11-05T00:00:00Z' WHERE id = ?",
                               (unit_name, str(job_dir), job["id"]))
            connection.commit()
        worker = nodeworker.NodeWorker(load_node_config())
        reset_calls = []
        saved = nodeworker.reset_failed
        nodeworker.reset_failed = lambda launch, name: reset_calls.append(name) or True
        try:
            unit = UnitState(name=unit_name, load_state="loaded", active_state="failed", sub_state="failed",
                             result="oom-kill", main_pid=0, exec_main_status=9, exec_main_code=2, invocation_id="x",
                             environment="", memory_peak=128 * MiB, cpu_usage_nsec=2_500_000_000)
            worker._reconcile_unit_job(self.store.require(job["id"]), unit, startup=False)
        finally:
            nodeworker.reset_failed = saved
            worker.store.close()
        row = self.store.require(job["id"])
        self.assertEqual((row["state"], row["result"]), ("failed", "oom-kill"))
        self.assertIn("out of memory", row["message"])
        self.assertEqual((row["mem_peak_bytes"], row["mem_peak_source"], row["cpu_usage_usec"]),
                         (128 * MiB, "unit", 2_500_000))
        self.assertEqual(reset_calls, [unit_name])
        auto = json.loads(self.cli("submit", "--json", "-p", "pa", "-n", "oom", "--mem", "auto", "--test-only",
                                   "--", "/bin/true").stdout)
        self.assertEqual(auto["request"]["mem_bytes"], 256 * MiB)  # 1.25 x the limit it was killed at


class SharedMemoryCleanupTests(NodeFixture, unittest.TestCase):
    needs_units = False

    def test_leftover_mpi_segments_are_removed_and_others_kept(self) -> None:
        from simqueue import nodeworker
        from simqueue.nodeconfig import load_node_config

        shm = self.base / "shm"
        shm.mkdir()
        started = time.time() - 30
        leftover = shm / "sm_segment.zenobia.1000.4b2a0001.0"
        leftover.write_bytes(b"x" * 4096)
        vader = shm / "vader_segment.zenobia.7.1"
        vader.write_bytes(b"y")
        unrelated = shm / "my-app-cache"
        unrelated.write_bytes(b"z")
        older = shm / "sm_segment.zenobia.1000.old.0"
        older.write_bytes(b"o")
        os.utime(older, (started - 3600, started - 3600))
        holder = subprocess.Popen([sys.executable, "-c", f"f=open({str(vader)!r}); import sys; sys.stdin.read()"],
                                  stdin=subprocess.PIPE)
        try:
            time.sleep(0.5)
            worker = nodeworker.NodeWorker(load_node_config())
            worker.shm_dir = shm
            iso = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(started))
            worker._clean_shm({"id": 7, "started_at": iso})
            worker.store.close()
        finally:
            holder.stdin.close()
            holder.wait(timeout=10)
        self.assertFalse(leftover.exists(), "a leftover MPI segment of the job is removed")
        self.assertTrue(vader.exists(), "a segment some process holds stays")
        self.assertTrue(unrelated.exists(), "other files stay")
        self.assertTrue(older.exists(), "segments older than the job stay")
        self.assertIn("removed 1 leftover shared-memory segment(s)", self.worker_log())


class SuiteGuardTests(unittest.TestCase):
    def test_real_oom_is_never_allowed_on_zenobia(self) -> None:
        import socket

        from node_fixture import real_oom_allowed

        saved = os.environ.get("SIMQ_TEST_REAL_OOM")
        os.environ["SIMQ_TEST_REAL_OOM"] = "1"
        try:
            self.assertEqual(real_oom_allowed(), socket.gethostname() != "zenobia")
        finally:
            if saved is None:
                os.environ.pop("SIMQ_TEST_REAL_OOM", None)
            else:
                os.environ["SIMQ_TEST_REAL_OOM"] = saved
        result = subprocess.run([PYTHON_FROZEN, str(Path(__file__).resolve().parents[1] / "run-tests.py"), "--pure",
                                 "test_none_*.py"], env=dict(os.environ, SIMQ_TEST_REAL_OOM="1"),
                                text=True, capture_output=True, check=False)
        if socket.gethostname() == "zenobia":
            self.assertEqual(result.returncode, 64)
            self.assertIn("refused on zenobia", result.stderr)


if __name__ == "__main__":
    unittest.main()
