"""Node scheduler with real systemd user units on the leased test CPUs."""

from __future__ import annotations

import json
import os
import subprocess
import time
import unittest
from pathlib import Path

from node_fixture import FAKE_GPU, PYTHON_FROZEN, NodeFixture, flock_free

from simqueue.config import format_cpu_list, parse_cpu_list


ENV_PROBE = (
    "import json,os\n"
    "keys=sorted(k for k in os.environ if k.startswith(('SIMQ_','CUDA_','SLURM_')))\n"
    "chain=[]; pid=os.getpid()\n"
    "for _ in range(4):\n"
    "    pid=int(open(f'/proc/{pid}/stat').read().rsplit(')',1)[1].split()[1]); chain.append(pid)\n"
    "print('PROBE '+json.dumps({'env':{k:os.environ[k] for k in keys},"
    "'affinity':sorted(os.sched_getaffinity(0)),"
    "'cgroup':open('/proc/self/cgroup').read().strip(),'parents':chain}))\n"
)


def probe_payload(text: str) -> dict:
    for line in text.splitlines():
        if line.startswith("PROBE "):
            return json.loads(line[len("PROBE "):])
    raise AssertionError("probe output missing:\n" + text)


class PlacementTests(NodeFixture, unittest.TestCase):
    def test_placement_environment_unit_and_affinity(self) -> None:
        env = dict(self.environment, SLURM_JOB_ID="4242", SLURM_NTASKS="9")
        job_id = self.submit(
            "-p", "pb", "--gpus", "1", "--gpu-mem", "3G", "--mem", "300M", "-t", "7:00",
            "-e", "SLURM_KEEP=explicit", "--", PYTHON_FROZEN, "-c", ENV_PROBE, env=env,
        )
        self.run_drain()
        job = self.job(job_id)
        self.assertEqual(job["state"], "succeeded", self.log_text(job_id))
        payload = probe_payload(self.log_text(job_id))
        environment = payload["env"]
        cpu = self.part_b[0]
        self.assertEqual(environment["SIMQ_JOB_ID"], str(job_id))
        self.assertEqual(environment["SIMQ_PARTITION"], "pb")
        self.assertEqual(environment["SIMQ_CPU_LIST"], str(cpu))
        self.assertEqual(environment["SIMQ_CPU_COUNT"], "1")
        self.assertEqual(environment["SIMQ_NCPUS"], "1")
        self.assertEqual(environment["SIMQ_GPUS"], FAKE_GPU)
        self.assertEqual(environment["CUDA_VISIBLE_DEVICES"], FAKE_GPU)
        self.assertEqual(environment["SIMQ_GPU_MEM"], str(3 * 1024**3))
        self.assertEqual(environment["SIMQ_GPU_MODE"], "shared")
        self.assertEqual(environment["SIMQ_MEM"], str(300 * 1024**2))
        self.assertEqual(environment["SIMQ_TIME_LIMIT"], "420")
        self.assertEqual(environment["SIMQ_LAUNCH_TOKEN"], job["launch_token"])
        self.assertEqual(environment["SIMQ_UNIT"], job["unit_name"])
        # Inherited SLURM_* would send pf-launch to srun; only -e values pass.
        self.assertNotIn("SLURM_JOB_ID", environment)
        self.assertNotIn("SLURM_NTASKS", environment)
        self.assertEqual(environment["SLURM_KEEP"], "explicit")
        self.assertEqual(payload["affinity"], [cpu])
        self.assertTrue(payload["cgroup"].endswith("/" + job["unit_name"]), payload["cgroup"])
        self.assertEqual(job["partition"], "pb")
        self.assertEqual(job["gpu_list"], FAKE_GPU)
        self.assertEqual(job["result"], "success")

    def test_job_without_gpus_gets_empty_cuda_visible_devices(self) -> None:
        job_id = self.submit("-p", "pa", "-c", "1", "--", PYTHON_FROZEN, "-c", ENV_PROBE)
        self.run_drain()
        environment = probe_payload(self.log_text(job_id))["env"]
        self.assertEqual(environment["CUDA_VISIBLE_DEVICES"], "")
        self.assertEqual(environment["SIMQ_GPUS"], "")
        self.assertEqual(environment["SIMQ_GPU_MODE"], "none")

    def test_job_processes_are_not_worker_descendants(self) -> None:
        job_id = self.submit("-p", "pa", "-c", "1", "--", PYTHON_FROZEN, "-c", ENV_PROBE)
        worker = self.start_worker()
        self.wait_until(lambda: self.state_of(job_id) == "succeeded", 60, "job completion")
        payload = probe_payload(self.log_text(job_id))
        self.assertNotIn(worker.pid, payload["parents"])
        worker.terminate()
        worker.wait(timeout=20)


class ConcurrencyTests(NodeFixture, unittest.TestCase):
    def test_partitions_run_concurrently_fifo_and_no_cpu_overlap(self) -> None:
        width = len(self.part_a)
        first = [self.submit("-p", "pa", "-c", "1", *self.probe_command(2.0)) for _ in range(width)]
        waiting = self.submit("-p", "pa", "-c", "1", *self.probe_command(0.5))
        other = self.submit("-p", "pb", "-c", "1", *self.probe_command(1.0))
        self.run_drain()
        for job_id in [*first, waiting, other]:
            self.assertEqual(self.state_of(job_id), "succeeded", self.log_text(job_id))
        records = self.read_intervals()
        starts = [records[j]["start"] for j in first]
        ends = [records[j]["end"] for j in first]
        if width > 1:
            self.assertLess(max(starts), min(ends), "partition jobs did not run concurrently")
        self.assertGreaterEqual(records[waiting]["start"], min(ends) - 0.05)
        self.assertLess(records[other]["start"], min(ends), "other partition waited")
        used = {tuple(records[j]["affinity"]) for j in first}
        self.assertEqual(len(used), width)
        self.assert_no_cpu_overlap()

    def test_whole_partition_job_waits_for_all_cores(self) -> None:
        if len(self.part_a) < 2:
            self.skipTest("needs a two-core partition")
        short = self.submit("-p", "pa", "-c", "1", *self.probe_command(1.5))
        whole = self.submit("-p", "pa", *self.probe_command(0.3))  # default: all partition cores
        self.run_drain()
        records = self.read_intervals()
        self.assertEqual(records[whole]["affinity"], list(self.part_a))
        self.assertGreaterEqual(records[whole]["start"], records[short]["end"] - 0.05)
        self.assertEqual(self.job(whole)["cpu_list"], format_cpu_list(self.part_a))


class AdoptionTests(NodeFixture, unittest.TestCase):
    def test_kill9_worker_jobs_survive_are_adopted_by_relocated_code(self) -> None:
        long_a = self.submit("-p", "pa", "-c", "1", *self.probe_command(6.0))
        long_b = self.submit("-p", "pb", "-c", "1", *self.probe_command(6.0, 7))
        if len(self.part_a) > 1:
            short = self.submit("-p", "pa", "-c", "1", *self.probe_command(1.0))
            whole = self.submit("-p", "pa", *self.probe_command(0.3))
        else:
            short = None
            whole = self.submit("-p", "pa", *self.probe_command(0.3))
        worker = self.start_worker()
        expected_running = [long_a, long_b] + ([short] if short else [])
        self.wait_until(lambda: all(self.state_of(j) == "running" for j in expected_running), 60,
                        "jobs running")
        units = {j: str(self.job(j)["unit_name"]) for j in expected_running}
        self.kill9(worker)
        killed_at = time.time()
        for unit in units.values():
            self.assertTrue(self.unit_active(unit), f"{unit} died with the worker")
        if short:
            # The short job finishes while no worker is running.
            status = Path(str(self.job(short)["job_dir"])) / "status.json"
            self.wait_until(lambda: json.loads(status.read_text()).get("phase") == "finished"
                            if status.exists() else False, 30, "short job to finish unsupervised")
            self.assertEqual(self.state_of(short), "running")  # database not yet reconciled
        self.assertEqual(self.state_of(whole), "queued")
        relocated = self.copy_runtime("relocated")
        restart_iso = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        restarted = self.start_worker("--drain", launcher=self.make_launcher(relocated))
        restarted.wait(timeout=120)
        self.assertEqual(restarted.returncode, 0)
        a, b = self.job(long_a), self.job(long_b)
        self.assertEqual(a["state"], "succeeded")
        self.assertEqual((b["state"], b["exit_code"]), ("failed", 7))
        self.assertIsNotNone(a["adopted_at"])
        self.assertIsNotNone(b["adopted_at"])
        if short:
            s = self.job(short)
            self.assertEqual(s["state"], "succeeded")
            # The recorded end is when the job really ended (worker down).
            self.assertLessEqual(s["finished_at"], restart_iso)
        self.assertEqual(self.state_of(whole), "succeeded")
        records = self.read_intervals()
        self.assertGreaterEqual(records[whole]["start"], records[long_a]["end"] - 0.05)
        self.assert_no_cpu_overlap()
        log = self.worker_log()
        self.assertIn(f"adopted job {long_a}", log)
        self.assertIn(f"adopted job {long_b}", log)
        from simqueue import __version__

        self.assertIn(f"{__version__}+relocated", log)
        self.assertGreater(records[long_a]["end"], killed_at)

    def test_cancel_while_worker_down_and_cancel_after_adoption(self) -> None:
        first = self.submit("-p", "pa", "-c", "1", "--", "/bin/sleep", "120")
        second = self.submit("-p", "pb", "-c", "1", "--", "/bin/sleep", "120")
        worker = self.start_worker()
        self.wait_until(lambda: self.state_of(first) == "running" and self.state_of(second) == "running",
                        60, "both running")
        self.kill9(worker)
        second_unit = str(self.job(second)["unit_name"])
        result = self.cli("cancel", str(second))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.wait_until(lambda: not self.unit_active(second_unit), 30, "unit stop without worker")
        self.assertEqual(self.state_of(second), "cancelling")
        worker = self.start_worker()
        self.wait_until(lambda: self.state_of(second) == "canceled", 30, "reconciled cancel")
        self.wait_until(lambda: self.job(first)["adopted_at"], 30, "adoption")
        result = self.cli("cancel", str(first))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.wait_until(lambda: self.state_of(first) == "canceled", 30, "adopted cancel")
        self.assertFalse(self.unit_active(str(self.job(first)["unit_name"])))
        worker.terminate()
        worker.wait(timeout=20)

    def test_sigterm_to_worker_leaves_jobs_running(self) -> None:
        job_id = self.submit("-p", "pa", "-c", "1", *self.probe_command(3.0))
        worker = self.start_worker()
        self.wait_until(lambda: self.state_of(job_id) == "running", 60, "running")
        worker.terminate()
        worker.wait(timeout=20)
        unit = str(self.job(job_id)["unit_name"])
        self.assertTrue(self.unit_active(unit))
        self.assertIn("continue independently", self.worker_log())
        self.run_drain()
        self.assertEqual(self.state_of(job_id), "succeeded")


class LeaseTests(NodeFixture, unittest.TestCase):
    def test_outside_lease_holder_is_respected_and_job_units_hold_leases(self) -> None:
        if len(self.part_a) < 2:
            self.skipTest("needs a two-core partition")
        blocked_cpu, free_cpu = self.part_a[0], self.part_a[1]
        holder = self.hold_lease(self.locks / f"core{blocked_cpu}.lock")
        whole = self.submit("-p", "pa", *self.probe_command(0.2))
        single = self.submit("-p", "pa", "-c", "1", *self.probe_command(3.0))
        worker = self.start_worker()
        self.wait_until(lambda: self.state_of(single) == "running", 60, "single-core job start")
        self.assertEqual(self.job(single)["cpu_list"], str(free_cpu))
        self.assertEqual(self.state_of(whole), "queued")
        lease = self.locks / f"core{free_cpu}.lock"
        self.wait_until(lambda: self.lease_holder_pids(lease), 10, "job lease")
        self.assertIn("lease held outside the queue", str(self.job(whole)["message"]))
        self.wait_until(lambda: self.state_of(single) == "succeeded", 60, "single done")
        self.wait_until(lambda: flock_free(lease), 10, "lease release after job")
        time.sleep(1.0)
        self.assertEqual(self.state_of(whole), "queued")  # outside holder still present
        self.release_lease(holder)
        self.wait_until(lambda: self.state_of(whole) == "succeeded", 60, "whole job after release")
        worker.terminate()
        worker.wait(timeout=20)
        self.assert_no_cpu_overlap()

    def test_gpu_lease_outside_queue_blocks_gpu_jobs_only(self) -> None:
        holder = self.hold_lease(self.locks / "gpu.lock")
        gpu_job = self.submit("-p", "pb", "--gpus", "1", *self.probe_command(0.2))
        cpu_job = self.submit("-p", "pb", *self.probe_command(0.2))
        worker = self.start_worker()
        self.wait_until(lambda: self.state_of(cpu_job) == "succeeded", 60, "cpu job")
        self.assertEqual(self.state_of(gpu_job), "queued")
        self.release_lease(holder)
        self.wait_until(lambda: self.state_of(gpu_job) == "succeeded", 60, "gpu job")
        worker.terminate()
        worker.wait(timeout=20)

    def test_exclusive_and_shared_gpu_lease_modes(self) -> None:
        exclusive = self.submit("-p", "pb", "--gpus", "1", "--", "/bin/sleep", "3")
        worker = self.start_worker()
        self.wait_until(lambda: self.state_of(exclusive) == "running", 60, "gpu job")
        from simqueue.locks import lock_holders

        path = str(self.locks / f"gpu-{FAKE_GPU}.lock")
        self.wait_until(lambda: lock_holders([path]).get(path), 10, "gpu lease")
        self.assertEqual(lock_holders([path])[path][0].mode, "WRITE")
        self.wait_until(lambda: self.state_of(exclusive) == "succeeded", 60, "done")
        shared = self.submit("-p", "pb", "--gpu-mem", "1G", "--", "/bin/sleep", "3")
        self.wait_until(lambda: self.state_of(shared) == "running", 60, "shared gpu job")
        self.wait_until(lambda: lock_holders([path]).get(path), 10, "shared gpu lease")
        self.assertEqual(lock_holders([path])[path][0].mode, "READ")
        worker.terminate()
        worker.wait(timeout=20)


class LimitTests(NodeFixture, unittest.TestCase):
    def test_time_limit_is_enforced_by_the_unit(self) -> None:
        job_id = self.submit("-p", "pa", "-c", "1", "-t", "0:03", "--", "/bin/sleep", "60")
        started = time.monotonic()
        self.run_drain()
        job = self.job(job_id)
        self.assertEqual((job["state"], job["result"]), ("failed", "timeout"))
        self.assertIn("time limit", job["message"])
        self.assertLess(time.monotonic() - started, 40)

    def test_out_of_memory_kills_the_whole_job(self) -> None:
        job_id = self.submit(
            "-p", "pa", "-c", "1", "--mem", "96M", "--",
            PYTHON_FROZEN, "-c", "import time; x=bytearray(400*1024*1024); time.sleep(5)",
        )
        self.run_drain()
        job = self.job(job_id)
        self.assertEqual((job["state"], job["result"]), ("failed", "oom-kill"), self.log_text(job_id))
        self.assertIn("out of memory", job["message"])


RESOURCE_RUNTIME = Path("/home/treefor/projects/runtimes/simqueue-resource-20261030")
LIVE_RECORDS = Path("/home/treefor/.local/state/simqueue/resource-scopes")


class BoundedJobCompatibilityTests(NodeFixture, unittest.TestCase):
    """The Spheral bounded-job.py gate checks SIMQ_CPU_LIST / SIMQ_CPU_COUNT
    against the process affinity and needs SIMQ_JOB_ID / SIMQ_LAUNCH_TOKEN."""

    def setUp(self) -> None:
        super().setUp()
        if not (RESOURCE_RUNTIME / "scope-entry.py").exists():
            self.skipTest("simqueue-resource runtime not present")

    def test_scope_entry_gate_accepts_a_06_allocation(self) -> None:
        receipt = self.base / "scope-entry-receipt.json"
        script = (
            'set -eu\n'
            'unit="simqt-gate-${SIMQ_JOB_ID}-$(printf %s "$SIMQ_LAUNCH_TOKEN" | cut -c1-12).scope"\n'
            f'exec systemd-run --user --scope --quiet --unit="$unit" -p MemoryMax=268435456 '
            f'-p MemorySwapMax=0 -p OOMPolicy=kill {PYTHON_FROZEN} {RESOURCE_RUNTIME}/scope-entry.py '
            f'--memory-bytes 268435456 --swap-bytes 0 --unit "$unit" --receipt {receipt} '
            '-- /bin/sh -c \'echo "gate-ok cpus=$SIMQ_CPU_LIST count=$SIMQ_CPU_COUNT"\'\n'
        )
        job_id = self.submit("-p", "pa", "--", "/bin/sh", "-c", script)
        self.run_drain()
        job = self.job(job_id)
        self.assertEqual(job["state"], "succeeded", self.log_text(job_id))
        record = json.loads(receipt.read_text())
        self.assertEqual(record["cpu_affinity"], list(self.part_a))
        self.assertIn(f"gate-ok cpus={format_cpu_list(self.part_a)} count={len(self.part_a)}", self.log_text(job_id))

    def test_unmodified_bounded_job_runs_with_records_redirected(self) -> None:
        """Runs bounded-job.py itself.  Its record directory is hard-coded to
        the live queue state, so the job runs in a private mount namespace
        where that directory is a bind mount of a scratch directory."""
        probe = subprocess.run(
            ["unshare", "--user", "--map-current-user", "--mount", "--keep-caps", "--",
             "/bin/sh", "-c", f"mount --bind {self.base} {LIVE_RECORDS} && echo ok"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False,
        )
        if probe.stdout.strip() != "ok":
            self.skipTest(f"private mount namespace unavailable: {probe.stderr.strip()}")
        before = sorted(p.name for p in LIVE_RECORDS.iterdir())
        records = self.base / "redirected-records"
        records.mkdir()
        script = (
            "set -eu\n"
            f"mount --bind {records} {LIVE_RECORDS}\n"
            f'test "$(stat -c %d:%i {LIVE_RECORDS})" = "$(stat -c %d:%i {records})"\n'
            f"exec {PYTHON_FROZEN} {RESOURCE_RUNTIME}/bounded-job.py --validation "
            "--memory-bytes 268435456 --swap-bytes 0 -- "
            "/bin/sh -c 'echo bounded-ok affinity=$(grep Cpus_allowed_list /proc/self/status | cut -f2)'\n"
        )
        job_id = self.submit(
            "-p", "pa", "--",
            "unshare", "--user", "--map-current-user", "--mount", "--keep-caps", "--",
            "/bin/sh", "-c", script,
        )
        self.run_drain()
        job = self.job(job_id)
        self.assertEqual(job["state"], "succeeded", self.log_text(job_id))
        self.assertIn("bounded-ok", self.log_text(job_id))
        entry = json.loads((records / f"{job_id}-entry.json").read_text())
        final = json.loads((records / f"{job_id}.json").read_text())
        self.assertEqual(entry["cpu_affinity"], list(self.part_a))
        self.assertTrue(final["memory_policy_verified"])
        self.assertTrue(final["scope_empty"])
        self.assertEqual(final["process_returncode"], 0)
        after = sorted(p.name for p in LIVE_RECORDS.iterdir())
        self.assertEqual(before, after, "live resource-scope records were touched")


class EnsembleNodeTests(NodeFixture, unittest.TestCase):
    def manifest(self, name: str, *, budget: str | None, claim: str, members: int, parallel: int) -> Path:
        document = {
            "schema_version": "simqueue-ensemble-v1",
            "name": name,
            "max_parallel": parallel,
            "allow_shared_cwd": True,
            "defaults": {"cwd": str(self.projects), "gpu_memory": claim},
            "members": [
                {"name": f"m{i}", "argv": [PYTHON_FROZEN, "-c", ENV_PROBE]} for i in range(members)
            ],
        }
        if budget is not None:
            document["gpu"] = {"budget": budget}
        path = self.projects / f"{name}.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        return path

    def member_payloads(self, job_id: int) -> list[dict]:
        directory = Path(str(self.store.require_ensemble(job_id)["artifact_dir"])) / "members"
        return [probe_payload(p.read_text()) for p in sorted(directory.glob("*.log"))]

    def test_cpu_ensemble_runs_inside_its_allocation(self) -> None:
        width = len(self.part_a)
        path = self.manifest("cpu-ens", budget=None, claim="0B", members=3, parallel=width)
        result = self.cli("ensemble", "submit", "-p", "pa", str(path), "--json")
        self.assertEqual(result.returncode, 0, result.stderr)
        job_id = int(json.loads(result.stdout)["job_id"])
        self.run_drain()
        job = self.job(job_id)
        self.assertEqual(job["state"], "succeeded", self.log_text(job_id))
        self.assertEqual(job["cpu_list"], format_cpu_list(self.part_a))
        for payload in self.member_payloads(job_id):
            self.assertEqual(payload["env"]["CUDA_VISIBLE_DEVICES"], "")
            self.assertEqual(len(payload["affinity"]), 1)
            self.assertIn(payload["affinity"][0], self.part_a)

    def test_gpu_ensemble_gets_a_shared_gpu_claim_of_its_budget(self) -> None:
        path = self.manifest("gpu-ens", budget="2GiB", claim="1GiB", members=2, parallel=1)
        result = self.cli("ensemble", "submit", "-p", "pb", str(path), "--json")
        self.assertEqual(result.returncode, 0, result.stderr)
        job_id = int(json.loads(result.stdout)["job_id"])
        self.run_drain()
        job = self.job(job_id)
        self.assertEqual(job["state"], "succeeded", self.log_text(job_id))
        self.assertEqual(job["gpu_list"], FAKE_GPU)
        self.assertEqual(job["gpu_claim_bytes"], 2 * 1024**3)
        self.assertEqual(job["gpu_exclusive"], 0)
        for payload in self.member_payloads(job_id):
            self.assertEqual(payload["env"]["CUDA_VISIBLE_DEVICES"], FAKE_GPU)


class EnsembleDrainTests(NodeFixture, unittest.TestCase):
    def test_pause_drains_releases_the_allocation_and_resume_continues_same_id(self) -> None:
        document = {
            "schema_version": "simqueue-ensemble-v1", "name": "drain-ens", "max_parallel": 1,
            "allow_shared_cwd": True, "defaults": {"cwd": str(self.projects), "gpu_memory": "0B"},
            "members": [{"name": f"m{i}", "argv": ["/bin/sleep", "2"]} for i in range(3)],
        }
        path = self.projects / "drain.json"
        path.write_text(json.dumps(document))
        result = self.cli("ensemble", "submit", "-p", "pa", str(path), "--json")
        self.assertEqual(result.returncode, 0, result.stderr)
        job_id = int(json.loads(result.stdout)["job_id"])
        worker = self.start_worker()
        self.wait_until(lambda: self.state_of(job_id) == "running", 60, "ensemble start")
        first_unit = str(self.job(job_id)["unit_name"])
        self.wait_until(lambda: any(m["state"] == "running" for m in self.store.list_ensemble_members(job_id)),
                        30, "first member")
        paused = self.cli("ensemble", "pause", str(job_id))
        self.assertEqual(paused.returncode, 0, paused.stderr)
        self.wait_until(lambda: self.state_of(job_id) == "queued", 60, "drained allocation released")
        job = self.job(job_id)
        self.assertIsNone(job["cpu_list"])
        self.assertIsNone(job["unit_name"])
        self.assertFalse(self.unit_active(first_unit))
        states = [m["state"] for m in self.store.list_ensemble_members(job_id)]
        self.assertEqual(states.count("succeeded"), 1)
        self.assertEqual(states.count("queued"), 2)
        resumed = self.cli("ensemble", "resume", str(job_id))
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        self.wait_until(lambda: self.state_of(job_id) == "succeeded", 60, "resumed ensemble")
        self.assertNotEqual(self.job(job_id)["unit_name"], first_unit)
        self.assertEqual([m["state"] for m in self.store.list_ensemble_members(job_id)], ["succeeded"] * 3)
        worker.terminate()
        worker.wait(timeout=20)


class LegacyJobTests(NodeFixture, unittest.TestCase):
    def test_running_05_job_is_observed_never_killed_and_keeps_its_cores(self) -> None:
        cpu = self.part_a[0]
        legacy = self.store.submit(
            name="legacy-05", cwd=self.projects, argv=["/bin/sleep", "1"],
            environment={}, env_override_keys=[], requested_cores=1,
        )
        token = str(legacy["launch_token"])
        marker = self.base / "legacy-finished"
        process = subprocess.Popen(
            ["taskset", "--cpu-list", str(cpu), "/bin/sh", "-c", f"sleep 4; touch {marker}"],
            env={"SIMQ_LAUNCH_TOKEN": token, "PATH": os.environ["PATH"]},
            start_new_session=True,
        )
        with self.store.connect() as connection:
            connection.execute(
                "UPDATE jobs SET state='running', cpu_list=?, pid=?, pgid=?, started_at=? WHERE id=?",
                (str(cpu), process.pid, process.pid, "2026-10-30T00:00:00Z", legacy["id"]),
            )
        newcomer = self.submit("-p", "pa", "--cpus", str(cpu), *self.probe_command(0.2))
        once = self.cli("worker", "--once")
        self.assertEqual(once.returncode, 0, once.stderr)
        self.assertEqual(self.state_of(int(legacy["id"])), "running")
        self.assertEqual(self.state_of(newcomer), "queued")
        self.assertIn("still running", self.worker_log())
        process.wait(timeout=30)
        self.assertTrue(marker.exists(), "the 0.5 job was killed")
        self.run_drain()
        job = self.job(int(legacy["id"]))
        self.assertEqual(job["state"], "interrupted")
        self.assertIn("exit status unknown", job["message"])
        self.assertEqual(self.state_of(newcomer), "succeeded")


class RefusalTests(NodeFixture, unittest.TestCase):
    needs_units = False

    def test_wrapper_refuses_a_busy_lease_without_running_the_command(self) -> None:
        from simqueue.units import JOBWRAP

        lease = self.locks / "busy.lock"
        holder = self.hold_lease(lease)
        marker = self.base / "must-not-exist"
        spec = {
            "schema": "simqueue-job-spec-v1", "job_id": 1, "name": "x", "launch_token": "ab" * 16,
            "unit": "none.service", "argv": ["/bin/touch", str(marker)], "env": {"PATH": "/usr/bin"},
            "cwd": str(self.projects), "project_root": str(self.projects),
            "log_path": str(self.base / "x.log"), "status_path": str(self.base / "status.json"),
            "stop_post_path": str(self.base / "stop-post.json"),
            "cpu_list": format_cpu_list(sorted(os.sched_getaffinity(0))),
            "locks": [{"path": str(lease), "mode": "exclusive"}], "require_unit": False,
        }
        spec_path = self.base / "spec.json"
        spec_path.write_text(json.dumps(spec))
        result = subprocess.run([PYTHON_FROZEN, str(JOBWRAP), "run", str(spec_path)],
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
        self.assertEqual(result.returncode, 254, result.stdout)
        status = json.loads((self.base / "status.json").read_text())
        self.assertEqual((status["phase"], status["refusal"]), ("refused", "lock-busy"))
        self.assertFalse(marker.exists())
        self.release_lease(holder)
        result = subprocess.run([PYTHON_FROZEN, str(JOBWRAP), "run", str(spec_path)],
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertTrue(marker.exists())

    def test_lock_busy_refusal_requeues_the_same_job(self) -> None:
        from simqueue.nodeconfig import load_node_config
        from simqueue.nodeworker import NodeWorker
        from simqueue.units import UnitState

        job_id = int(self.store.submit(
            name="race", cwd=self.projects, argv=["/bin/true"], environment={},
            env_override_keys=[], requested_cores=1, scheduling={"partition": "pa"},
        )["id"])
        claimed = self.store.claim_allocation(
            job_id, cpu_list=str(self.part_a[0]), gpu_list="", gpu_claim_bytes=0, gpu_exclusive=False,
            mem_bytes=1, swap_bytes=0, time_limit_s=60, partition="pa", unit_prefix=self.prefix,
            job_root=self.state / "jobs", lock_paths=[], runtime_dir="x", busy_cpus=set(),
        )
        job_dir = Path(str(claimed["job_dir"]))
        job_dir.mkdir(parents=True)
        (job_dir / "status.json").write_text(json.dumps({
            "launch_token": claimed["launch_token"], "phase": "refused", "refusal": "lock-busy",
            "reason": "lease file held by another process",
        }))
        worker = NodeWorker(load_node_config())
        unit = UnitState(str(claimed["unit_name"]), "not-found", "inactive", "dead", "", 0, None, None, "", "")
        worker._finalize_unit_job(self.job(job_id), unit, *worker._job_files(self.job(job_id)))
        job = self.job(job_id)
        self.assertEqual(job["state"], "queued")
        self.assertIsNone(job["unit_name"])
        self.assertIn("requeued", job["message"])

    def test_unit_that_never_started_is_requeued_not_duplicated(self) -> None:
        from simqueue.nodeconfig import load_node_config
        from simqueue.nodeworker import NodeWorker
        from simqueue.units import UnitState

        job_id = int(self.store.submit(
            name="crash-before-start", cwd=self.projects, argv=["/bin/true"], environment={},
            env_override_keys=[], requested_cores=1, scheduling={"partition": "pa"},
        )["id"])
        claimed = self.store.claim_allocation(
            job_id, cpu_list=str(self.part_a[0]), gpu_list="", gpu_claim_bytes=0, gpu_exclusive=False,
            mem_bytes=1, swap_bytes=0, time_limit_s=60, partition="pa", unit_prefix=self.prefix,
            job_root=self.state / "jobs", lock_paths=[], runtime_dir="x", busy_cpus=set(),
        )
        old_token = claimed["launch_token"]
        worker = NodeWorker(load_node_config())
        unit = UnitState(str(claimed["unit_name"]), "not-found", "inactive", "dead", "", 0, None, None, "", "")
        worker._finalize_unit_job(self.job(job_id), unit, None, None)
        job = self.job(job_id)
        self.assertEqual(job["state"], "queued")
        self.assertNotEqual(job["launch_token"], old_token)

    def test_claim_refuses_overlapping_allocation(self) -> None:
        ids = [int(self.store.submit(
            name=f"j{i}", cwd=self.projects, argv=["/bin/true"], environment={},
            env_override_keys=[], requested_cores=1, scheduling={"partition": "pa"},
        )["id"]) for i in range(2)]
        common = dict(gpu_list="", gpu_claim_bytes=0, gpu_exclusive=False, mem_bytes=1, swap_bytes=0,
                      time_limit_s=60, partition="pa", unit_prefix=self.prefix, job_root=self.state / "jobs",
                      lock_paths=[], runtime_dir="x", busy_cpus=set())
        self.assertIsNotNone(self.store.claim_allocation(ids[0], cpu_list="0", **common))
        with self.assertRaises(RuntimeError):
            self.store.claim_allocation(ids[1], cpu_list="0", **common)
        self.assertEqual(self.state_of(ids[1]), "queued")


class NodeCliTests(NodeFixture, unittest.TestCase):
    needs_units = False

    def test_submit_without_partition_uses_default_and_05_core_semantics(self) -> None:
        result = self.cli("submit", "--json", "--", "/bin/true")
        self.assertEqual(result.returncode, 0, result.stderr)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["partition"], "pa")
        self.assertEqual(payload["requested_cores"], len(self.part_a))  # 0.5: default = all cores
        result = self.cli("submit", "--json", "-c", "1", "--", "/bin/true")
        self.assertEqual(json.loads(result.stdout)["requested_cores"], 1)
        result = self.cli("submit", "-c", str(len(self.part_a) + 1), "--", "/bin/true")
        self.assertEqual(result.returncode, 2)
        self.assertIn("between 1 and", result.stderr)

    def test_request_validation(self) -> None:
        cases = [
            (("-p", "nope"), "unknown partition"),
            (("-p", "pa", "--mem", "3G"), "--mem must be"),
            (("-p", "pa", "-t", "2:00:00"), "--time exceeds"),
            (("-p", "pa", "--gpus", "1"), "GPU(s)"),
            (("-p", "pb", "--gpu-mem", "11G"), "exceeds the GPU budget"),
            (("-p", "pa", "--cpus", "999"), "not in partition"),
        ]
        for options, message in cases:
            with self.subTest(options=options):
                result = self.cli("submit", *options, "--", "/bin/true")
                self.assertEqual(result.returncode, 2, result.stdout)
                self.assertIn(message, result.stderr)

    def test_partitions_and_priority_commands(self) -> None:
        low = self.submit("-p", "pa", "--", "/bin/true")
        high = self.submit("-p", "pa", "--priority", "5", "--", "/bin/true")
        queued = [int(j["id"]) for j in self.store.queued_jobs()]
        self.assertEqual(queued, [high, low])
        result = self.cli("priority", str(low), "9")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([int(j["id"]) for j in self.store.queued_jobs()], [low, high])
        result = self.cli("demote", str(low))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([int(j["id"]) for j in self.store.queued_jobs()], [high, low])
        result = self.cli("partitions", "--json")
        self.assertEqual(result.returncode, 0, result.stderr)
        status = json.loads(result.stdout)
        self.assertEqual([p["name"] for p in status["partitions"]], ["pa", "pb"])
        self.assertEqual(status["partitions"][0]["queued"], 2)

    def test_idempotent_05_style_replay_keeps_its_fingerprint(self) -> None:
        """A 0.5-syntax keyed submission gets the fingerprint the 0.5 code
        computes, so controllers replaying across the upgrade still match."""
        env = {"PATH": os.environ["PATH"], "A": "1",
               **{k: v for k, v in self.environment.items() if k.startswith("SIMQUEUE_")}}
        result = self.cli("submit", "--json", "-n", "replay", "-c", "1", "--idempotency-key", "k:1",
                          "--", "/bin/true", env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        first = json.loads(result.stdout)
        legacy = subprocess.run(
            [PYTHON_FROZEN, "-c",
             "import json,sys\n"
             "sys.path.insert(0, '/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/source')\n"
             "from pathlib import Path\n"
             "from simqueue.store import _command_request_sha256\n"
             "a=json.loads(sys.argv[1])\n"
             "print(_command_request_sha256(name='replay', cwd=Path(a['cwd']), argv=['/bin/true'],"
             " environment=a['env'], env_override_keys=[], requested_cores=1, retry_of=None, attempt=1))",
             json.dumps({"cwd": str(self.projects),
                         "env": json.loads(self.job(int(first["job_id"]))["env_json"])})],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True,
        )
        self.assertEqual(first["request_sha256"], legacy.stdout.strip())
        again = json.loads(self.cli("submit", "--json", "-n", "replay", "-c", "1", "--idempotency-key", "k:1",
                                    "--", "/bin/true", env=env).stdout)
        self.assertEqual((again["job_id"], again["created"]), (first["job_id"], False))
        changed = self.cli("submit", "--json", "-n", "replay", "-c", "1", "-p", "pb", "--idempotency-key",
                           "k:1", "--", "/bin/true", env=env)
        self.assertEqual(changed.returncode, 2)
        self.assertIn("different submission", changed.stderr)

    def test_legacy_mode_rejects_node_options(self) -> None:
        env = dict(self.environment)
        env.pop("SIMQUEUE_CONFIG")
        env["SIMQUEUE_PCORES"] = format_cpu_list(sorted(os.sched_getaffinity(0))[:2])
        result = self.cli("submit", "-p", "pa", "--", "/bin/true", env=env)
        self.assertEqual(result.returncode, 2)
        self.assertIn("need a node configuration", result.stderr)


if __name__ == "__main__":
    unittest.main()


PF_LAUNCH_11 = Path("/home/treefor/projects/runtimes/pf-launch-v1.1-20261030")
PF_LAUNCH_10 = Path("/home/treefor/projects/Power-Flow/tools/pf-launch")


class PfLaunchInsideAllocationTests(NodeFixture, unittest.TestCase):
    """pf-launch 1.1 inside a real 0.6 job unit uses the allocation and takes
    no leases (its own unit holds them); 1.0 refuses on those same leases."""

    def setUp(self) -> None:
        super().setUp()
        if not (PF_LAUNCH_11 / "pf_launch.py").exists():
            self.skipTest("pf-launch 1.1 not present")
        site = json.loads((PF_LAUNCH_11 / "sites" / "zenobia.json").read_text())
        site["hostname"] = os.uname().nodename
        site["reservations"] = []
        # Only the leased test cores exist for pf-launch here, so neither
        # version can ever place a process outside them.
        site["reserved_cores"] = {"cores": [], "reason": ""}
        site["core_preference"] = list(self.part_a)
        site["leases"] = {"core_path": str(self.locks / "core{core}.lock"),
                          "gpu_path": str(self.locks / "gpu-{uuid}.lock")}
        site["execution"] = {"kind": "host", "prefix": [], "mpi_inside_prefix": False, "wrapper_python": None,
                             "root": None, "manifest": None, "image": None, "image_digest": None, "hash_image": False}
        site["identity"] = {"policy": "record-only", "record": {}}
        site["launcher"]["nvidia_smi"] = "/nonexistent/nvidia-smi"  # no GPU queries in tests
        site["receipt_dir"] = str(self.base / "receipts")
        self.site = self.base / "site.json"
        self.site.write_text(json.dumps(site))

    def pf_job(self, launcher_dir: Path, receipt: Path) -> int:
        # Explicit --cpus = the allocation: pf-launch 1.0 cannot pick other cores.
        return self.submit(
            "-p", "pa", "--", PYTHON_FROZEN, str(launcher_dir / "pf_launch.py"), "--site", str(self.site),
            "--cpus", format_cpu_list(self.part_a),
            "--launcher", "exec", "--gpus", "none", "--receipt", str(receipt), "--",
            "/bin/sh", "-c", "echo PFCHILD affinity=$(grep Cpus_allowed_list /proc/self/status | cut -f2) "
                             "cuda=[$CUDA_VISIBLE_DEVICES] core=$PF_CORE",
        )

    def test_pf_launch_11_uses_the_allocation_without_double_locking(self) -> None:
        receipt = self.base / "pf11.json"
        job_id = self.pf_job(PF_LAUNCH_11, receipt)
        self.run_drain()
        self.assertEqual(self.state_of(job_id), "succeeded", self.log_text(job_id))
        cpus = format_cpu_list(self.part_a)
        self.assertIn(f"PFCHILD affinity={cpus} cuda=[]", self.log_text(job_id))
        data = json.loads(receipt.read_text())
        self.assertEqual(data["selected_cores"]["simqueue"]["SIMQ_JOB_ID"], str(job_id))
        self.assertEqual(data["selected_cores"]["cpus"], list(self.part_a))
        self.assertEqual(data["selected_cores"]["leases"], [])
        self.assertEqual(data["exit_codes"]["pf_launch"], 0)

    def test_pf_launch_10_refuses_on_its_own_jobs_leases(self) -> None:
        if not (PF_LAUNCH_10 / "pf_launch.py").exists():
            self.skipTest("pf-launch 1.0 not present")
        if 'VERSION = "1.0' not in (PF_LAUNCH_10 / "pf_launch.py").read_text():
            self.skipTest(f"{PF_LAUNCH_10} is no longer pf-launch 1.0")
        receipt = self.base / "pf10.json"
        job_id = self.pf_job(PF_LAUNCH_10, receipt)
        self.run_drain()
        job = self.job(job_id)
        self.assertEqual((job["state"], job["exit_code"]), ("failed", 75), self.log_text(job_id))
        self.assertIn("leased", self.log_text(job_id))
