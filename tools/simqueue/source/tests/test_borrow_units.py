"""Borrowing (0.6.1) with real systemd user units.

Real units run only on the leased test CPUs (``SIMQ_TEST_CPUS``, cores 0-3
taken by run-tests.py), with temporary lease files, a per-test unit prefix
and a per-test borrow slice that is reverted afterwards.  The powerflow test
observes the real lease files of cores 6-7 through the scheduler only: its
borrow job sleeps (no CPU work), runs SCHED_IDLE in the idle slice, and is
never placed on a core whose lease holder is busy.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import unittest
from pathlib import Path

from node_fixture import PYTHON_FROZEN, NodeFixture, flock_free

from simqueue.locks import lock_holders


USER_CGROUP = Path(f"/sys/fs/cgroup/user.slice/user-{os.getuid()}.slice/user@{os.getuid()}.service")
REAL_LEASES = {
    6: "/home/treefor/projects/.workspace/outputs/pic-profile/testing-core.lock",
    7: "/home/treefor/projects/.workspace/outputs/pic-profile/testing-core7.lock",
}

INFO = r'''
import json, os, sys, time
out, seconds = sys.argv[1], float(sys.argv[2])
cgroup = next(l[3:] for l in open("/proc/self/cgroup").read().splitlines() if l.startswith("0::"))
record = {"pid": os.getpid(), "policy": os.sched_getscheduler(0), "affinity": sorted(os.sched_getaffinity(0)),
          "cgroup": cgroup, "env": {k: v for k, v in os.environ.items() if k.startswith("SIMQ_")}}
with open(out + ".tmp", "w") as h:
    json.dump(record, h)
os.replace(out + ".tmp", out)
time.sleep(seconds)
'''
BURN = r'''
import sys, time
end = time.time() + float(sys.argv[1])
while time.time() < end:
    pass
'''
TICKER = r'''
import os, sys, time
path, count = sys.argv[1], int(sys.argv[2])
for _ in range(count):
    with open(path, "a") as h:
        h.write(f"{os.getpid()} {time.time()}\n")
    time.sleep(0.1)
'''
# Owner job for the owner-idle test: one busy process pinned to the first
# allocated CPU; once TRIGGER exists, one more on every other allocated CPU.
OWNER_LOAD = r"""
import os, subprocess, sys, time
trigger, seconds = sys.argv[1], float(sys.argv[2])
cpus = sorted(os.sched_getaffinity(0))
burn = "import time\nend=time.time()+%f\nwhile time.time()<end: pass\n"
end = time.time() + seconds
procs = [subprocess.Popen([sys.executable, "-c", burn % seconds])]
os.sched_setaffinity(procs[0].pid, {cpus[0]})
triggered = False
while time.time() < end:
    if not triggered and os.path.exists(trigger):
        triggered = True
        for cpu in cpus[1:]:
            proc = subprocess.Popen([sys.executable, "-c", burn % max(0.5, end - time.time())])
            os.sched_setaffinity(proc.pid, {cpu})
            procs.append(proc)
    time.sleep(0.05)
for proc in procs:
    proc.wait()
"""
# A tool outside the queue: holds core lease files; "burn" starts a busy child.
HOLDER = r'''
import fcntl, subprocess, sys
handles = []
for path in sys.argv[1:]:
    h = open(path, "a"); fcntl.flock(h, fcntl.LOCK_EX); handles.append(h)
print("held", flush=True)
children = []
for line in sys.stdin:
    if line.strip() == "burn":
        children.append(subprocess.Popen([sys.executable, "-c", "while True: pass"]))
        print("burning", flush=True)
for child in children:
    child.kill()
'''


class BorrowFixture(NodeFixture):
    min_cpus = 2
    policies: dict = {}
    borrow_settings: dict = {}

    def setUp(self) -> None:  # noqa: N802
        super().setUp()
        for name, text in (("info.py", INFO), ("burn.py", BURN), ("ticker.py", TICKER),
                           ("owner_load.py", OWNER_LOAD), ("holder.py", HOLDER)):
            (self.projects / name).write_text(text, encoding="utf-8")

    def tearDown(self) -> None:  # noqa: N802
        # systemd refuses to stop a frozen unit: thaw this test's units first.
        subprocess.run(["systemctl", "--user", "thaw", f"{self.prefix}-*"], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=False)
        super().tearDown()
        subprocess.run(["systemctl", "--user", "revert", self.borrow_slice], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=False)
        subprocess.run(["systemctl", "--user", "stop", self.borrow_slice], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=False)

    @property
    def borrow_slice(self) -> str:
        return self.prefix + "b.slice"

    def config_document(self):
        document = super().config_document()
        document["borrow"] = {
            "slice": self.prefix + "b.slice",
            "idle_window_s": 3,
            "busy_window_s": 1,
            "freeze_max_s": None,
            **self.borrow_settings,
        }
        for part in document["partitions"]:
            part.update(self.policies.get(part["name"], {}))
        return document

    # ----- helpers -----------------------------------------------------------

    def py(self, script: str, *arguments: str) -> list[str]:
        return ["--", PYTHON_FROZEN, str(self.projects / script), *arguments]

    def wait_state(self, job_id: int, *states: str, timeout: float = 60) -> dict:
        return self.wait_until(lambda: self.job(job_id) if self.state_of(job_id) in states else None,
                               timeout=timeout, message=f"job {job_id} in {states}")

    def unit_cgroup(self, job_id: int) -> Path:
        job = self.job(job_id)
        return USER_CGROUP / self.borrow_slice / str(job["unit_name"])

    def cpu_usage_usec(self, path: Path) -> int:
        for line in (path / "cpu.stat").read_text().splitlines():
            if line.startswith("usage_usec"):
                return int(line.split()[1])
        raise AssertionError("no usage_usec")

    def events(self, job_id: int) -> list[dict]:
        return self.store.job_events(job_id)

    def kinds(self, job_id: int) -> list[str]:
        return [event["kind"] for event in self.events(job_id)]

    def freezer_state(self, unit: str) -> str:
        return subprocess.run(["systemctl", "--user", "show", "-P", "FreezerState", unit], text=True,
                              stdout=subprocess.PIPE, check=False).stdout.strip()

    def start_holder(self, *paths: Path) -> subprocess.Popen:
        proc = subprocess.Popen([sys.executable, str(self.projects / "holder.py"), *map(str, paths)],
                                stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        self.assertEqual(proc.stdout.readline().strip(), "held")
        self.holders.append(proc)
        return proc


class BorrowUnitTests(BorrowFixture, unittest.TestCase):
    def test_borrow_unit_is_idle_priority_and_holds_only_the_borrow_lease(self) -> None:
        self.start_worker()
        out = self.projects / "info.json"
        job_id = self.submit("--borrow", "-p", "pa", "-c", "1", "-n", "info", *self.py("info.py", str(out), "4"))
        self.wait_until(out.exists, message="borrow probe")
        info = json.loads(out.read_text())
        job = self.job(job_id)
        cpu = int(job["cpu_list"])
        self.assertEqual(cpu, self.part_a[-1])  # home partition, highest free core
        self.assertEqual(info["policy"], os.SCHED_IDLE)
        self.assertEqual(info["affinity"], [cpu])
        self.assertEqual(info["env"]["SIMQ_BORROW"], "1")
        self.assertEqual(info["env"]["SIMQ_ON_PREEMPT"], "requeue")
        self.assertTrue(info["cgroup"].endswith(f"/{self.borrow_slice}/{job['unit_name']}"), info["cgroup"])
        cgroup = USER_CGROUP / self.borrow_slice
        self.assertEqual((cgroup / job["unit_name"] / "cpu.weight").read_text().strip(), "1")
        self.assertEqual((cgroup / "cpu.idle").read_text().strip(), "1")
        ionice = subprocess.run(["ionice", "-p", str(info["pid"])], text=True, stdout=subprocess.PIPE, check=False)
        self.assertEqual(ionice.stdout.strip(), "idle")
        borrow_lease = self.locks / f"borrow-core{cpu}.lock"
        held = lock_holders([str(borrow_lease)]).get(str(borrow_lease), [])
        self.assertTrue(held and all(h.mode == "READ" for h in held), held)
        self.assertTrue(flock_free(self.locks / f"core{cpu}.lock"), "a borrow job must not hold the core lease")
        placement = json.loads(job["placement_json"])
        self.assertEqual(placement["cpus"], [{"class": "free", "cpu": cpu, "partition": "pa"}])
        self.wait_state(job_id, "succeeded")
        self.assertEqual(self.kinds(job_id), ["started", "finished"])

    def test_owner_process_keeps_cpu_priority_over_a_running_borrower(self) -> None:
        """A busy process outside the queue (no lease taken, so the scheduler
        does not preempt) on the borrower's core gets the CPU: the borrow
        slice is idle-weight, the unit CPUWeight=1 and SCHED_IDLE."""
        from simqueue.telemetry import read_cpu_times

        self.start_worker()
        job_id = self.submit("--borrow", "-p", "pa", "-c", "1", "-n", "burner", *self.py("burn.py", "20"))
        job = self.wait_state(job_id, "running")
        cpu = int(job["cpu_list"])
        cgroup = self.unit_cgroup(job_id)
        time.sleep(1.0)
        cpu_before, alone_before = read_cpu_times()[cpu], self.cpu_usage_usec(cgroup)
        time.sleep(1.0)
        cpu_after, alone_after = read_cpu_times()[cpu], self.cpu_usage_usec(cgroup)
        alone = (alone_after - alone_before) / 1e6
        core_busy = (cpu_after[0] - cpu_before[0]) / max(1, cpu_after[1] - cpu_before[1])
        competitor = subprocess.Popen([sys.executable, str(self.projects / "burn.py"), "3.5"])
        os.sched_setaffinity(competitor.pid, {cpu})
        time.sleep(0.5)
        before = self.cpu_usage_usec(cgroup)
        start = time.monotonic()
        time.sleep(2.0)
        borrowed = (self.cpu_usage_usec(cgroup) - before) / 1e6
        span = time.monotonic() - start
        _, status, usage = os.wait4(competitor.pid, 0)
        competitor_cpu = usage.ru_utime + usage.ru_stime
        self.cli("cancel", str(job_id))
        result = {"cpu": cpu, "alone_share": alone, "core_busy_alone": core_busy,
                  "borrower_share_under_contention": borrowed / span, "competitor_cpu_s": competitor_cpu,
                  "competitor_wall_s": 3.5}
        out = os.environ.get("SIMQ_TEST_RECEIPTS")
        if out:
            Path(out).mkdir(parents=True, exist_ok=True)
            (Path(out) / "cpu-priority.json").write_text(json.dumps(result, indent=2))
        self.assertGreater(core_busy, 0.9, "the core was not busy while the borrower ran alone")
        self.assertLess(borrowed / span, 0.05, f"borrower used {borrowed:.3f} s of {span:.2f} s under contention")
        # the competitor shares the leased cores with the test harness itself
        self.assertGreater(competitor_cpu, 0.5 * 3.5, result)
        self.assertGreater(competitor_cpu / 3.5, 20 * borrowed / span, result)


class PreemptionUnitTests(BorrowFixture, unittest.TestCase):
    def test_requeue_preemption_owner_starts_at_once_and_borrower_reruns(self) -> None:
        self.start_worker()
        borrow_id = self.submit("--borrow", "-p", "pa", "-c", str(len(self.part_a)), "--borrow-from", "pa", "-n", "build",
                                *self.probe_command(4))
        self.wait_state(borrow_id, "running")
        time.sleep(0.5)
        submitted = time.time()
        owner_id = self.submit("-p", "pa", "-c", str(len(self.part_a)), "-n", "owner", *self.probe_command(2))
        owner = self.wait_state(owner_id, "running", "succeeded", timeout=30)
        started = owner["started_at"]
        from datetime import datetime

        started_unix = datetime.fromisoformat(started.replace("Z", "+00:00")).timestamp()
        self.assertLess(started_unix - submitted, 4.0, "the owner job waited for the borrower")
        self.wait_until(lambda: "requeued" in self.kinds(borrow_id), message="requeue")
        preempted = next(e for e in self.events(borrow_id) if e["kind"] == "preempted")
        self.assertEqual(preempted["detail"]["mode"], "requeue")
        self.assertEqual(preempted["detail"]["by_job"], owner_id)
        self.wait_state(owner_id, "succeeded")
        job = self.wait_state(borrow_id, "succeeded", timeout=60)
        self.assertEqual(job["preempt_count"], 1)
        self.assertEqual(self.kinds(borrow_id), ["started", "preempted", "requeued", "started", "finished"])
        self.assertIn("reclaimed", [e["kind"] for e in self.store.job_events(owner_id)])
        self.assert_no_cpu_overlap()
        records = self.read_intervals()
        self.assertIn(borrow_id, records)

    def test_freeze_preemption_stops_the_cpu_and_resumes_the_same_process(self) -> None:
        self.start_worker()
        ticks = self.projects / "ticks.txt"
        borrow_id = self.submit("--borrow", "-p", "pa", "-c", str(len(self.part_a)), "--borrow-from", "pa", "--on-preempt",
                                "freeze", "-n", "ticker", *self.py("ticker.py", str(ticks), "60"))
        self.wait_until(lambda: ticks.exists() and len(ticks.read_text().splitlines()) >= 10, message="ticks")
        owner_id = self.submit("-p", "pa", "-c", str(len(self.part_a)), "-n", "owner", *self.probe_command(3))
        self.wait_state(owner_id, "running")
        job = self.wait_until(lambda: self.job(borrow_id) if self.job(borrow_id)["frozen_at"] else None,
                              message="frozen")
        unit = str(job["unit_name"])
        self.assertEqual(self.freezer_state(unit), "frozen")
        self.assertEqual((self.unit_cgroup(borrow_id) / "cgroup.freeze").read_text().strip(), "1")
        count = len(ticks.read_text().splitlines())
        usage = self.cpu_usage_usec(self.unit_cgroup(borrow_id))
        time.sleep(1.5)
        self.assertEqual(len(ticks.read_text().splitlines()), count, "a frozen job made progress")
        self.assertLess(self.cpu_usage_usec(self.unit_cgroup(borrow_id)) - usage, 20_000)
        self.assertEqual(self.state_of(borrow_id), "running")
        listing = self.cli("list").stdout
        self.assertIn("frozen", listing)
        self.wait_state(owner_id, "succeeded")
        job = self.wait_state(borrow_id, "succeeded", timeout=60)
        lines = [line.split() for line in ticks.read_text().splitlines()]
        self.assertEqual(len(lines), 60)
        self.assertEqual(len({pid for pid, _ in lines}), 1, "the frozen job was restarted instead of resumed")
        gaps = [float(b[1]) - float(a[1]) for a, b in zip(lines, lines[1:])]
        self.assertGreater(max(gaps), 2.0)
        self.assertEqual(self.kinds(borrow_id), ["started", "preempted", "resumed", "finished"])
        self.assertEqual(job["preempt_count"], 1)

    def test_cancel_preemption_ends_the_borrower(self) -> None:
        self.start_worker()
        borrow_id = self.submit("--borrow", "-p", "pa", "-c", str(len(self.part_a)), "--borrow-from", "pa", "--on-preempt", "cancel",
                                *self.probe_command(30))
        self.wait_state(borrow_id, "running")
        owner_id = self.submit("-p", "pa", "-c", str(len(self.part_a)), *self.probe_command(1))
        job = self.wait_state(borrow_id, "canceled", timeout=30)
        self.assertEqual(job["result"], "preempted")
        self.assertIn("canceled (--on-preempt cancel)", job["message"])
        self.wait_state(owner_id, "succeeded")
        self.assertEqual(self.kinds(borrow_id), ["started", "preempted", "finished"])

    def test_worker_restart_keeps_a_frozen_borrower_frozen_then_resumes_it(self) -> None:
        worker = self.start_worker()
        ticks = self.projects / "ticks.txt"
        borrow_id = self.submit("--borrow", "-p", "pa", "-c", str(len(self.part_a)), "--borrow-from", "pa", "--on-preempt",
                                "freeze", *self.py("ticker.py", str(ticks), "40"))
        self.wait_until(lambda: ticks.exists() and len(ticks.read_text().splitlines()) >= 5, message="ticks")
        owner_id = self.submit("-p", "pa", "-c", str(len(self.part_a)), *self.probe_command(6))
        self.wait_until(lambda: self.job(borrow_id)["frozen_at"], message="frozen")
        self.kill9(worker)
        unit = str(self.job(borrow_id)["unit_name"])
        self.assertEqual(self.freezer_state(unit), "frozen")
        relocated = self.copy_runtime("restart")
        self.start_worker(launcher=self.make_launcher(relocated))
        self.wait_until(lambda: f"adopted job {borrow_id}" in self.worker_log(), message="adoption")
        self.assertIn("borrow job, frozen", self.worker_log())
        self.assertEqual(self.state_of(owner_id), "running")
        self.assertEqual(self.freezer_state(unit), "frozen")
        self.wait_state(owner_id, "succeeded")
        self.wait_state(borrow_id, "succeeded", timeout=60)
        self.assertEqual(len({line.split()[0] for line in ticks.read_text().splitlines()}), 1)
        self.assertIn("+restart", self.worker_log())

    def test_a_lease_taken_on_a_lent_core_preempts_the_borrower(self) -> None:
        self.start_worker()
        borrow_id = self.submit("--borrow", "-p", "pa", "-c", "1", *self.probe_command(30))
        job = self.wait_state(borrow_id, "running")
        cpu = int(job["cpu_list"])
        taken = time.monotonic()
        holder = self.hold_lease(self.locks / f"core{cpu}.lock")
        self.wait_until(lambda: "preempted" in self.kinds(borrow_id), timeout=15, message="preemption")
        latency = time.monotonic() - taken
        event = next(e for e in self.events(borrow_id) if e["kind"] == "preempted")
        self.assertIn("taken outside the queue", event["detail"]["reason"])
        self.assertLess(latency, 5.0)
        # requeued, then restarted on another core
        self.wait_until(lambda: self.kinds(borrow_id).count("started") == 2, timeout=30, message="restart")
        job = self.wait_state(borrow_id, "running", timeout=30)
        self.assertNotEqual(int(job["cpu_list"]), cpu)
        self.assertIn("requeued", self.kinds(borrow_id))
        self.release_lease(holder)
        self.cli("cancel", str(borrow_id))


class LeasedIdleUnitTests(BorrowFixture, unittest.TestCase):
    policies = {"pb": {"lend_leased_idle": True}, "pa": {"lend": False}}

    def test_idle_lease_holder_lends_its_core_and_busy_holder_preempts(self) -> None:
        self.start_worker()
        filler = self.submit("-p", "pa", "-c", str(len(self.part_a)), *self.probe_command(60))
        self.wait_state(filler, "running")
        leased = self.part_b[-1]
        holder = self.start_holder(self.locks / f"core{leased}.lock")
        borrow_id = self.submit("--borrow", "-p", "pa", "-c", str(len(self.part_b)), "--borrow-from", "pb",
                                *self.probe_command(60))
        job = self.wait_state(borrow_id, "running", timeout=30)
        placement = {item["cpu"]: item for item in json.loads(job["placement_json"])["cpus"]}
        self.assertEqual(placement[leased]["class"], "leased-idle")
        self.assertIn(f"pid {holder.pid}", placement[leased]["holder"])
        status = json.loads(self.cli("partitions", "--json", "--sample", "0.5").stdout)
        cpus = {c["cpu"]: c for p in status["partitions"] for c in p["cpus"]}
        self.assertEqual(cpus[leased]["class"], "external-idle")
        self.assertEqual(cpus[leased]["holders"][0]["pid"], holder.pid)
        self.assertEqual(cpus[leased]["borrowers"][0]["job_id"], borrow_id)
        busy_at = time.monotonic()
        holder.stdin.write("burn\n")
        holder.stdin.flush()
        self.assertEqual(holder.stdout.readline().strip(), "burning")
        self.wait_until(lambda: "preempted" in self.kinds(borrow_id), timeout=20, message="busy preemption")
        latency = time.monotonic() - busy_at
        event = next(e for e in self.events(borrow_id) if e["kind"] == "preempted")
        self.assertIn("is busy", event["detail"]["reason"])
        self.assertLess(latency, 6.0)
        status = json.loads(self.cli("partitions", "--json", "--sample", "0.5").stdout)
        cpus = {c["cpu"]: c for p in status["partitions"] for c in p["cpus"]}
        self.assertEqual(cpus[leased]["class"], "external-busy")
        holder.stdin.close()
        holder.wait(timeout=10)
        self.cli("cancel", str(borrow_id))
        self.cli("cancel", str(filler))
        (self.base / "latency.json").write_text(json.dumps({"busy_preemption_s": latency}))

    def test_lease_classification_with_and_without_sampling(self) -> None:
        holder = self.hold_lease(self.locks / f"core{self.part_b[0]}.lock")
        status = json.loads(self.cli("partitions", "--json", "--sample", "0").stdout)
        cpus = {c["cpu"]: c for p in status["partitions"] for c in p["cpus"]}
        self.assertEqual(cpus[self.part_b[0]]["class"], "external")
        self.assertEqual(cpus[self.part_a[0]]["class"], "free")
        text = self.cli("partitions", "--sample", "0.3").stdout
        self.assertIn(f"holder pid {holder.pid}", text)
        self.assertIn("external lease idle", text)
        self.release_lease(holder)


class OwnerIdleUnitTests(BorrowFixture, unittest.TestCase):
    """One partition with every test CPU: an owner job takes all of them but
    keeps only the first busy; a borrow job of the same partition may use the
    unused ones (owner-idle) until the owner needs them."""

    def config_document(self):
        document = super().config_document()
        document["partitions"] = [
            {"name": "pw", "cpus": ",".join(map(str, self.cpus)), "default_mem": "256M", "max_mem": "2G",
             "default_time": "10:00", "max_time": "1:00:00"},
        ]
        document["default_partition"] = "pw"
        document["gpus"] = []
        document["lock_files"].pop("gpu_extra", None)
        return document

    def run_owner_and_borrower(self, *owner_options: str) -> tuple[int, int, Path]:
        self.start_worker()
        trigger = self.projects / "trigger"
        owner_id = self.submit("-p", "pw", "-c", str(len(self.cpus)), *owner_options,
                               *self.py("owner_load.py", str(trigger), "40"))
        self.wait_state(owner_id, "running")
        borrow_id = self.submit("--borrow", "-p", "pw", "-c", "1", *self.probe_command(60))
        return owner_id, borrow_id, trigger

    def test_idle_core_of_a_running_owner_is_lent_and_reclaimed(self) -> None:
        owner_id, borrow_id, trigger = self.run_owner_and_borrower()
        job = self.wait_state(borrow_id, "running", timeout=30)
        placement = json.loads(job["placement_json"])["cpus"]
        self.assertEqual(len(placement), 1)
        self.assertEqual(placement[0]["class"], "owner-idle")
        self.assertEqual(placement[0]["owner_job"], owner_id)
        self.assertEqual(placement[0]["cpu"], max(self.cpus))
        triggered = time.monotonic()
        trigger.write_text("go")
        self.wait_until(lambda: "preempted" in self.kinds(borrow_id), timeout=20, message="owner reclaim")
        latency = time.monotonic() - triggered
        event = next(e for e in self.events(borrow_id) if e["kind"] == "preempted")
        self.assertIn(f"job {owner_id} needs CPU", event["detail"]["reason"])
        self.assertLess(latency, 6.0)
        out = os.environ.get("SIMQ_TEST_RECEIPTS")
        if out:
            Path(out).mkdir(parents=True, exist_ok=True)
            (Path(out) / "owner-idle-reclaim.json").write_text(json.dumps(
                {"reclaim_latency_s": latency, "placement": placement, "events": self.events(borrow_id)},
                indent=2, default=str))
        for job_id in (borrow_id, owner_id):
            self.cli("cancel", str(job_id))

    def test_exclusive_timing_owner_is_never_lent(self) -> None:
        owner_id, borrow_id, _ = self.run_owner_and_borrower("--exclusive-timing")
        deadline = time.monotonic() + 9.0
        while time.monotonic() < deadline:
            self.assertEqual(self.state_of(borrow_id), "queued")
            time.sleep(0.5)
        self.assertIn("exclusive timing", str(self.job(borrow_id)["message"]))
        shown = self.cli("show", str(owner_id)).stdout
        self.assertIn("Exclusive timing", shown)
        for job_id in (borrow_id, owner_id):
            self.cli("cancel", str(job_id))


class PowerflowLeaseTests(BorrowFixture, unittest.TestCase):
    """Cores 6-7 are observed only through the scheduler and their real lease
    files; the borrow job sleeps and may run there only while the lease
    holder's tree is idle."""

    min_cpus = 2

    def config_document(self):
        document = super().config_document()
        document["lock_files"]["core_overrides"] = {str(cpu): path for cpu, path in REAL_LEASES.items()}
        document["lock_files"]["borrow_path"] = str(self.locks / "borrow-core{core}.lock")
        document["partitions"] = [
            {"name": "pa", "cpus": ",".join(map(str, self.part_a)), "default_mem": "256M", "max_mem": "2G",
             "default_time": "10:00", "max_time": "1:00:00", "lend": False},
            {"name": "pf", "cpus": "6-7", "default_mem": "256M", "max_mem": "2G", "default_time": "10:00",
             "max_time": "1:00:00", "lend_leased_idle": True},
        ]
        document["default_partition"] = "pa"
        return document

    def test_busy_lease_holders_on_cores_6_7_are_never_borrowed(self) -> None:
        before = {cpu: sorted(h.pid for h in lock_holders([path]).get(path, [])) for cpu, path in REAL_LEASES.items()}
        worker = self.start_worker()
        borrow_id = self.submit("--borrow", "-p", "pf", "-c", "1", "--max-borrow-cores", "1", "-n", "pf-sleeper",
                                "-t", "2:00", "--", "/bin/sleep", "3")
        observations = []
        deadline = time.monotonic() + 12.0
        while time.monotonic() < deadline:
            status = json.loads(self.cli("partitions", "--json", "--sample", "1").stdout)
            cpus = {c["cpu"]: c for p in status["partitions"] for c in p["cpus"]}
            job = self.job(borrow_id)
            observations.append({
                "at": time.time(), "state": job["state"], "cpu_list": job["cpu_list"],
                "classes": {cpu: cpus[cpu]["class"] for cpu in (6, 7)},
                "usage": {cpu: cpus[cpu].get("usage_cores") for cpu in (6, 7)},
                "holders": {cpu: cpus[cpu].get("holders", []) for cpu in (6, 7)},
                "message": job["message"],
            })
            for cpu in (6, 7):
                if cpus[cpu]["class"] in ("external-busy", "external-idle"):
                    self.assertTrue(cpus[cpu]["holders"], cpus[cpu])
                    for holder in cpus[cpu]["holders"]:
                        self.assertTrue(holder["command"])
            if job["cpu_list"]:
                placement = json.loads(job["placement_json"])["cpus"]
                for item in placement:
                    self.assertIn(item["cpu"], (6, 7))
                    self.assertIn(item["class"], ("leased-idle", "free"))
            time.sleep(0.5)
        both_busy = all(o["classes"][6] == "external-busy" and o["classes"][7] == "external-busy"
                        for o in observations)
        if both_busy:
            self.assertTrue(all(o["state"] == "queued" for o in observations),
                            "a borrow job started on a core whose lease holder was busy")
        self.cli("cancel", str(borrow_id))
        self.kill9(worker)
        after = {cpu: sorted(h.pid for h in lock_holders([path]).get(path, [])) for cpu, path in REAL_LEASES.items()}
        ours = {proc.pid for proc in self.workers}
        for cpu in (6, 7):
            self.assertFalse(set(after[cpu]) & ours, "the scheduler took a real lease")
        receipt = {"before": before, "after": after, "observations": observations, "both_busy": both_busy,
                   "job_events": self.store.job_events(borrow_id), "worker_log": self.worker_log()}
        out = os.environ.get("SIMQ_TEST_RECEIPTS")
        if out:
            Path(out).mkdir(parents=True, exist_ok=True)
            (Path(out) / "powerflow-lease-observation.json").write_text(json.dumps(receipt, indent=2, default=str))


if __name__ == "__main__":
    unittest.main()
