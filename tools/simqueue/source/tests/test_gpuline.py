"""GPU line (0.6.2): parsing, ordering, and real flock(2) behaviour.

Every test uses private lock files in a temporary directory; no real lease
file is touched.  Raw users are real ``flock(1)`` processes, helpers are real
``gpuline.py lease`` processes and claims real ``gpuline.py claim`` processes
(started directly, not as units; the unit tests are in
test_gpu_line_units.py).
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

from simqueue import gpuline, jobwrap


REPOSITORY = Path(__file__).resolve().parents[1]
GPULINE = REPOSITORY / "simqueue" / "gpuline.py"
PYTHON = "/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen"
FAKE = "GPU-test-line-0000-000000000001"


class LineEnv:
    """Private lease files, a node document and an order log."""

    def __init__(self, base: Path, **line: object) -> None:
        self.base = base
        self.locks = base / "locks"
        self.locks.mkdir()
        self.order = base / "order.txt"
        self.order.write_text("")
        self.doc = {
            "schema": "simqueue-node-config-v1",
            "node": "testnode",
            "default_partition": "pa",
            "gpus": [{"uuid": FAKE, "index": 0, "name": "fake", "memory": "12G", "budget": "10G"}],
            "lock_files": {
                "enabled": True,
                "core_path": str(self.locks / "core{core}.lock"),
                "core_overrides": {"6": str(self.locks / "core.lock")},
                "gpu_path": str(self.locks / "gpu-{uuid}.lock"),
                "gpu_extra": {FAKE: [str(self.locks / "gpu.lock")]},
            },
            "partitions": [{"name": "pa", "cpus": "6-7", "gpus": [FAKE]}],
            "gpu_line": {"poll_s": 0.05, "claim_hold_max_s": 60, **line},
        }
        self.config_path = base / "node.json"
        self.config_path.write_text(json.dumps(self.doc, indent=2))
        self.cfg = gpuline.parse_line_config(self.doc)
        self.gpu = self.cfg.gpus[FAKE]
        self.l1, self.l2 = self.gpu.paths
        for path in (self.l1, self.l2, self.cfg.core_paths[6], self.cfg.core_paths[7]):
            Path(path).touch()
        self.procs: list[subprocess.Popen] = []

    def spawn(self, argv: list[str], **kwargs) -> subprocess.Popen:
        proc = subprocess.Popen(argv, stdin=subprocess.DEVNULL, **kwargs)
        self.procs.append(proc)
        return proc

    def work(self, tag: str, hold: float) -> list[str]:
        return ["sh", "-c", f"echo {tag} >> {self.order}; sleep {hold}"]

    def raw(self, tag: str, hold: float = 0.3) -> subprocess.Popen:
        return self.spawn(["flock", self.l1, "flock", self.l2, *self.work(tag, hold)])

    def raw_loop(self, tag: str, count: int, hold: float = 0.2) -> subprocess.Popen:
        script = (f"for i in $(seq {count}); do flock {self.l1} flock {self.l2} "
                  f"sh -c 'echo {tag} >> {self.order}; sleep {hold}'; done")
        return self.spawn(["bash", "-c", script])

    def helper(self, tag: str, hold: float = 0.3, *options: str) -> subprocess.Popen:
        return self.spawn([PYTHON, str(GPULINE), "lease", "--config", str(self.config_path), "-q", *options, "--",
                           *self.work(tag, hold)], stderr=subprocess.PIPE, text=True)

    def entries(self) -> list[str]:
        return self.order.read_text().split()

    def waiters(self, path: str | None = None) -> list[int]:
        state = gpuline.path_locks([path or self.l1])[path or self.l1]
        return [w.pid for w in state.waiters]

    def holders(self, path: str | None = None) -> list[int]:
        state = gpuline.path_locks([path or self.l1])[path or self.l1]
        return [h.pid for h in state.holders]

    def tickets(self) -> list[dict]:
        return gpuline.load_tickets(self.cfg.ticket_dir, gpu=FAKE)

    def cleanup(self) -> None:
        for proc in self.procs:
            if proc.poll() is None:
                proc.kill()
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                pass


def wait_for(predicate, timeout: float = 20.0, message: str = "condition"):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        value = predicate()
        if value:
            return value
        time.sleep(0.02)
    raise AssertionError(f"timed out waiting for {message}")


class LineTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="simq062-line-")
        self.env = LineEnv(Path(self.temporary.name))

    def tearDown(self) -> None:
        self.env.cleanup()
        self.temporary.cleanup()


# --------------------------------------------------------------------------- pure


class ReadinessTests(unittest.TestCase):
    """Which queued GPU jobs wait only for the GPU lease files (and so get a
    claim): powerflow-like partition h (cores 6-7) sharing GPU-x with g."""

    NOW = 1_000_000.0
    G = 1024**3

    def setUp(self) -> None:
        from simqueue.nodeconfig import parse_node_config

        part = {"default_mem": "1G", "max_mem": "8G", "default_time": "1:00:00", "max_time": "10:00:00"}
        self.config = parse_node_config({
            "schema": "simqueue-node-config-v1", "node": "test", "default_partition": "h",
            "memory_schedulable": "16G",
            "gpus": [{"uuid": "GPU-x", "index": 0, "memory": "12G", "budget": "10G"}],
            "partitions": [{"name": "g", "cpus": "4-5", "gpus": ["GPU-x"], **part},
                           {"name": "h", "cpus": "6-7", "gpus": ["GPU-x"], **part}],
        }, path=Path("/dev/null"), check_topology=False)

    def gpu_job(self, job_id: int, cpu: int, claim_g: float | None):
        from simqueue.scheduler import Allocation, Running

        exclusive = claim_g is None
        allocation = Allocation(cpus=(cpu,), gpus=("GPU-x",), gpu_claim=0 if exclusive else int(claim_g * self.G),
                                gpu_exclusive=exclusive, mem=self.G)
        return Running(job_id, allocation, self.NOW + 600)

    def request(self, job_id: int, gpu_mem_g: float | None, partition: str = "h"):
        from simqueue.scheduler import Request

        return Request(job_id=job_id, partition=partition, ncpus=1, gpus=1,
                       gpu_mem=None if gpu_mem_g is None else int(gpu_mem_g * self.G), mem=self.G, time_limit=120)

    def ready(self, running, owners, *, busy=(), blocked=(), started=()):
        from simqueue.gpuqueue import ready_except_gpu

        return ready_except_gpu(self.config, "GPU-x", running=running, owners=owners, busy_cpus=set(busy),
                                blocked=set(blocked), started=set(started), now=self.NOW)

    def test_shared_job_waiting_for_a_core_of_a_gpu_job_is_not_ready(self) -> None:
        # job 217 on zenobia: core 6 used by a shared GPU job, core 7 leased outside the queue
        self.assertEqual(self.ready([self.gpu_job(215, 6, 3)], [self.request(217, 0.25)], busy={7}), [])

    def test_exclusive_job_behind_queue_gpu_jobs_is_ready(self) -> None:
        # the claim must win the release at the end of the queue's GPU job
        self.assertEqual(self.ready([self.gpu_job(215, 6, 3)], [self.request(301, None)]), [301])
        self.assertEqual(self.ready([self.gpu_job(215, 6, 3)], [self.request(301, None)], busy={7}), [301])

    def test_shared_job_blocked_only_by_an_outside_holder_is_ready(self) -> None:
        self.assertEqual(self.ready([], [self.request(302, 1)], blocked={"GPU-x"}), [302])

    def test_shared_job_waiting_for_gpu_budget_is_not_ready(self) -> None:
        running = [self.gpu_job(211, 4, 6), self.gpu_job(215, 6, 3)]
        self.assertEqual(self.ready(running, [self.request(303, 4)]), [])

    def test_started_jobs_and_dispatch_order(self) -> None:
        owners = [self.request(304, None), self.request(305, None, partition="g")]
        self.assertEqual(self.ready([], owners, blocked={"GPU-x"}), [304])  # one exclusive at a time
        self.assertEqual(self.ready([], owners, blocked={"GPU-x"}, started={304}), [])


class ParsingTests(unittest.TestCase):
    def test_proc_locks_chain_order_and_depth(self) -> None:
        text = (
            "15: FLOCK  ADVISORY  WRITE 100 103:09:77 0 EOF\n"
            "15: -> FLOCK  ADVISORY  WRITE 101 103:09:77 0 EOF\n"
            "15:  -> FLOCK  ADVISORY  WRITE 102 103:09:77 0 EOF\n"
            "15:   -> FLOCK  ADVISORY  WRITE 103 103:09:77 0 EOF\n"
            "16: POSIX  ADVISORY  READ 200 00:3d:12 0 EOF\n"
            "17: OFDLCK ADVISORY  WRITE -1 103:09:5 0 EOF\n"
        )
        records = gpuline.parse_proc_locks(text)
        chain = [r for r in records if r.ino == 77]
        self.assertEqual([r.pid for r in chain], [100, 101, 102, 103])
        self.assertEqual([r.waiter for r in chain], [False, True, True, True])
        self.assertEqual([r.depth for r in chain], [0, 1, 2, 3])
        self.assertIsNone(records[-1].pid)
        self.assertEqual(records[-2].dev, "00:3d")

    def test_line_config_defaults_and_order(self) -> None:
        with tempfile.TemporaryDirectory(prefix="simq062-cfg-") as name:
            env = LineEnv(Path(name))
            self.assertTrue(env.cfg.enabled)
            self.assertEqual(env.gpu.paths, (str(env.locks / "gpu.lock"), str(env.locks / f"gpu-{FAKE}.lock")))
            self.assertEqual(env.cfg.state_path, str(env.locks / "simq-gpu-line.json"))
            self.assertEqual(env.cfg.ticket_dir, str(env.locks / "gpu-line.d"))
            self.assertEqual(env.cfg.core_paths[6], str(env.locks / "core.lock"))
            self.assertEqual(env.cfg.core_paths[7], str(env.locks / "core7.lock"))
            again = gpuline.LineConfig.from_dict(env.cfg.to_dict())
            self.assertEqual(again, env.cfg)
            bad = dict(env.doc, gpu_line={"bogus": 1})
            with self.assertRaises(ValueError):
                gpuline.parse_line_config(bad)
            bad = dict(env.doc, gpu_line={"requeue_commands": "flock"})
            with self.assertRaises(ValueError):
                gpuline.parse_line_config(bad)
            off = dict(env.doc, lock_files={"enabled": False})
            self.assertFalse(gpuline.parse_line_config(off).enabled)

    def test_turn_blockers(self) -> None:
        now = time.time()
        me = {"kind": "helper", "pid": 50, "gpu": FAKE, "since_unix": now - 10}
        tickets = [
            {"kind": "helper", "pid": 40, "gpu": FAKE, "since_unix": now - 20, "state": "queued"},
            {"kind": "helper", "pid": 60, "gpu": FAKE, "since_unix": now - 5, "state": "waiting"},
        ]
        self.assertEqual(gpuline.turn_blockers(me, tickets, None, now), [])
        tickets.append({"kind": "helper", "pid": 41, "gpu": FAKE, "since_unix": now - 15, "state": "waiting"})
        self.assertEqual(len(gpuline.turn_blockers(me, tickets, None, now)), 1)
        tickets.pop()
        state = {"schema": gpuline.STATE_SCHEMA, "heartbeat_unix": now, "stale_after_s": 120,
                 "gpus": {FAKE: {"queue": [{"job_id": 7, "since_unix": now - 30, "state": "waiting"}]}}}
        self.assertIn("queued job 7", gpuline.turn_blockers(me, tickets, state, now)[0])
        # the queue's claim already in the kernel queue covers the queued jobs
        claimed = tickets + [{"kind": "queue", "pid": 70, "gpu": FAKE, "since_unix": now - 30, "state": "queued"}]
        self.assertEqual(gpuline.turn_blockers(me, claimed, state, now), [])
        # a stale or withdrawn state is ignored
        stale = dict(state, heartbeat_unix=now - 1000)
        self.assertEqual(gpuline.turn_blockers(me, tickets, stale, now), [])
        # a later queued job does not block an earlier helper
        late = {"schema": gpuline.STATE_SCHEMA, "heartbeat_unix": now, "stale_after_s": 120,
                "gpus": {FAKE: {"queue": [{"job_id": 8, "since_unix": now - 1}]}}}
        self.assertEqual(gpuline.turn_blockers(me, tickets, late, now), [])
        # claims never look at the published queue; ties go to the queue
        claim = {"kind": "queue", "pid": 80, "gpu": FAKE, "since_unix": now - 10}
        self.assertEqual(gpuline.turn_blockers(claim, [me | {"state": "waiting"}], state, now), [])
        self.assertLess(gpuline.ticket_key(claim), gpuline.ticket_key(me))

    def test_job_control_risk(self) -> None:
        def stat(pid, comm, pgrp, tty=34816, ppid=1):
            return gpuline.ProcStat(pid, comm, "S", ppid, pgrp, pgrp, tty, pgrp, 1)

        table = {
            10: stat(10, "bash", 10),  # interactive shell
            11: stat(11, "flock", 11, ppid=10),  # its job: own process group
            12: stat(12, "flock", 10, ppid=10),  # same group as the shell (no job control)
            20: stat(20, "node", 20, tty=0),
            21: stat(21, "flock", 21, tty=0, ppid=20),  # no terminal
        }
        with mock.patch.object(gpuline, "proc_stat", side_effect=lambda pid: table.get(pid)):
            self.assertTrue(gpuline.job_control_risk(table[11]))
            self.assertFalse(gpuline.job_control_risk(table[12]))
            self.assertFalse(gpuline.job_control_risk(table[21]))


# --------------------------------------------------------------------------- real processes


class RawWaiterTests(LineTestCase):
    def test_kernel_queue_is_fifo_and_requeue_moves_a_waiter_to_the_tail(self) -> None:
        env = self.env
        holder = env.raw("h", 1.0)
        wait_for(lambda: env.holders(), message="holder")
        raws = []
        for index in range(3):
            raws.append(env.raw(f"w{index}", 0.1))
            wait_for(lambda: len(env.waiters()) == index + 1, message=f"waiter {index}")
        member = gpuline.LineMember(env.cfg, env.gpu, kind="helper", since_unix=time.time(), label="test")
        self.assertTrue(member._requeue(raws[0].pid, gpuline.pid_start(raws[0].pid)))
        wait_for(lambda: env.waiters()[-1:] == [raws[0].pid], 5, "requeued waiter at the tail")
        member.close()
        for proc in [holder, *raws]:
            proc.wait(timeout=20)
        self.assertEqual(env.entries(), ["h", "w1", "w2", "w0"])

    def test_helper_wins_the_first_release_over_earlier_raw_waiters(self) -> None:
        env = self.env
        holder = env.raw("r0", 1.5)
        wait_for(lambda: env.holders(), message="raw holder")
        raws = []
        for index in range(3):
            raws.append(env.raw(f"raw{index}", 0.2))
            wait_for(lambda: len(env.waiters()) == index + 1, message=f"raw waiter {index}")
        helper = env.helper("helper", 0.3)
        self.assertEqual(helper.wait(timeout=30), 0, helper.stderr.read() if helper.stderr else "")
        for proc in [holder, *raws]:
            proc.wait(timeout=30)
        entries = env.entries()
        self.assertEqual(entries[:2], ["r0", "helper"], entries)
        self.assertEqual(sorted(entries[2:]), ["raw0", "raw1", "raw2"])
        log = Path(env.cfg.log_path).read_text()
        self.assertIn("moved raw waiter pid", log)

    def test_helpers_are_fifo_by_start_time_against_back_to_back_raw_users(self) -> None:
        env = self.env
        holder = env.raw("r0", 1.5)
        wait_for(lambda: env.holders(), message="raw holder")
        loops = [env.raw_loop("loopA", 4), env.raw_loop("loopB", 4)]
        wait_for(lambda: len(env.waiters()) >= 2, message="raw loop waiters")
        helpers = []
        for index in range(3):
            helpers.append(env.helper(f"h{index}", 0.2))
            wait_for(lambda: len([t for t in env.tickets() if t.get("kind") == "helper"]) == index + 1,
                     message=f"helper {index} ticket")
            time.sleep(0.15)
        for proc in helpers:
            self.assertEqual(proc.wait(timeout=60), 0)
        for proc in [holder, *loops]:
            proc.wait(timeout=60)
        entries = env.entries()
        self.assertEqual(entries[:4], ["r0", "h0", "h1", "h2"], entries)
        self.assertEqual(entries.count("loopA") + entries.count("loopB"), 8)

    def test_helper_yields_to_an_earlier_queued_job_until_it_leaves_the_line(self) -> None:
        env = self.env
        now = time.time()
        state = {"schema": gpuline.STATE_SCHEMA, "heartbeat_unix": now, "stale_after_s": 120, "paused": False,
                 "gpus": {FAKE: {"queue": [{"job_id": 41, "name": "early", "since_unix": now - 60,
                                            "state": "waiting"}]}}}
        gpuline.write_json_atomic(Path(env.cfg.state_path), state)
        helper = env.helper("helper", 0.1)
        wait_for(lambda: any("queued job 41" in str(t.get("note")) for t in env.tickets()), 10, "yield note")
        time.sleep(0.5)
        self.assertEqual(env.entries(), [])
        self.assertIsNone(helper.poll())
        state["gpus"][FAKE]["queue"] = []
        state["heartbeat_unix"] = time.time()
        gpuline.write_json_atomic(Path(env.cfg.state_path), state)
        self.assertEqual(helper.wait(timeout=20), 0)
        self.assertEqual(env.entries(), ["helper"])

    def test_member_in_the_kernel_queue_backs_out_for_an_earlier_member(self) -> None:
        env = self.env
        holder = env.raw("r0", 2.0)
        wait_for(lambda: env.holders(), message="raw holder")
        helper = env.helper("helper", 0.2)
        wait_for(lambda: any(t.get("state") == "queued" for t in env.tickets()), 10, "helper queued in the kernel")
        # an earlier member (a claim on behalf of a job submitted before the helper started)
        script = (
            "import sys, time; sys.path.insert(0, sys.argv[1]); from simqueue import gpuline\n"
            "cfg = gpuline.parse_line_config(__import__('json').load(open(sys.argv[2])))\n"
            "gpu = next(iter(cfg.gpus.values()))\n"
            "m = gpuline.LineMember(cfg, gpu, kind='queue', since_unix=float(sys.argv[3]), label='early claim')\n"
            "m.acquire(); open(sys.argv[4], 'a').write('claim\\n'); time.sleep(0.2); m.close()\n"
        )
        claim = env.spawn([PYTHON, "-c", script, str(REPOSITORY), str(env.config_path), str(time.time() - 30),
                           str(env.order)])
        self.assertEqual(claim.wait(timeout=30), 0)
        self.assertEqual(helper.wait(timeout=30), 0)
        holder.wait(timeout=30)
        self.assertEqual(env.entries(), ["r0", "claim", "helper"])

    def test_reversed_lock_order_user_does_not_deadlock_the_helper(self) -> None:
        env = self.env
        # a holder of the global file only, so that the per-GPU file is free
        blocker = env.spawn(["flock", env.l1, *env.work("r0", 1.0)])
        wait_for(lambda: env.holders(), message="raw holder")
        # pf-launch --wait style: the per-GPU file first, then the global one
        reversed_user = env.spawn(["flock", env.l2, "flock", env.l1, *env.work("reversed", 0.2)])
        wait_for(lambda: len(env.waiters()) == 1, message="reversed user blocked on the global file")
        helper = env.helper("helper", 0.2)
        self.assertEqual(helper.wait(timeout=60), 0)
        self.assertEqual(reversed_user.wait(timeout=60), 0)
        blocker.wait(timeout=30)
        self.assertEqual(sorted(env.entries()), ["helper", "r0", "reversed"])
        self.assertIn("lock-order conflict", Path(env.cfg.log_path).read_text())

    def test_core_holder_waiting_for_the_gpu_does_not_deadlock_the_helper(self) -> None:
        env = self.env
        core = env.cfg.core_paths[6]
        blocker = env.raw("r0", 1.0)
        wait_for(lambda: env.holders(), message="raw holder")
        # core lease first, then the GPU: the order that once deadlocked
        core_first = env.spawn(["flock", core, "flock", env.l1, "flock", env.l2, *env.work("corefirst", 0.2)])
        wait_for(lambda: len(env.waiters()) == 1, message="core-first user blocked on the GPU")
        helper = env.helper("helper", 0.2, "--core", "6")
        self.assertEqual(helper.wait(timeout=60), 0)
        self.assertEqual(core_first.wait(timeout=60), 0)
        blocker.wait(timeout=30)
        self.assertEqual(sorted(env.entries()), ["corefirst", "helper", "r0"])

    def test_helper_core_any_takes_a_free_core_and_sets_the_environment(self) -> None:
        env = self.env
        hold7 = env.spawn(["flock", env.cfg.core_paths[6], "sleep", "3"])
        wait_for(lambda: env.holders(env.cfg.core_paths[6]), message="core 6 holder")
        out = env.base / "env.txt"
        helper = env.spawn([PYTHON, str(GPULINE), "lease", "--config", str(env.config_path), "-q", "--core", "any",
                            "--", "sh", "-c", f"echo $PF_CORE $GPU_LEASE_GPU $OMP_NUM_THREADS > {out}; "
                                              f"grep Cpus_allowed_list /proc/self/status >> {out}"])
        self.assertEqual(helper.wait(timeout=30), 0)
        text = out.read_text()
        self.assertIn(f"7 {FAKE} 1", text)
        self.assertIn("Cpus_allowed_list:\t7", text)
        hold7.kill()

    def test_helper_timeout_and_dead_member_cleanup(self) -> None:
        env = self.env
        holder = env.raw("r0", 3.0)
        wait_for(lambda: env.holders(), message="raw holder")
        helper = env.helper("never", 0.1, "--timeout", "0.5")
        self.assertEqual(helper.wait(timeout=20), gpuline.EXIT_TIMEOUT)
        self.assertEqual([t for t in env.tickets() if t["kind"] == "helper"], [])
        holder.wait(timeout=20)
        self.assertEqual(env.entries(), ["r0"])

    def test_waiter_left_stopped_by_a_dead_member_is_resumed(self) -> None:
        env = self.env
        victim = env.spawn(["sleep", "30"])
        wait_for(lambda: gpuline.proc_stat(victim.pid), message="victim")
        os.kill(victim.pid, signal.SIGSTOP)
        wait_for(lambda: gpuline.proc_stat(victim.pid).state == "T", message="victim stopped")
        dead = env.spawn(["true"])
        dead.wait()
        Path(env.cfg.ticket_dir).mkdir(parents=True, exist_ok=True)
        ticket = {"schema": gpuline.TICKET_SCHEMA, "kind": "helper", "gpu": FAKE, "pid": dead.pid,
                  "pid_start": 1, "since_unix": time.time(), "state": "queued",
                  "stopping": [{"pid": victim.pid, "start": gpuline.pid_start(victim.pid)}]}
        gpuline.write_json_atomic(Path(env.cfg.ticket_dir) / f"helper-{dead.pid}.json", ticket)
        self.assertEqual(env.tickets(), [])
        wait_for(lambda: gpuline.proc_stat(victim.pid).state in ("S", "R"), 5, "victim resumed")
        self.assertFalse((Path(env.cfg.ticket_dir) / f"helper-{dead.pid}.json").exists())

    def test_interactive_job_is_not_requeued(self) -> None:
        stat = gpuline.ProcStat(5, "flock", "S", 4, 5, 4, 34816, 5, 1)
        with mock.patch.object(gpuline, "proc_stat", side_effect=lambda pid: {
                5: stat, 4: gpuline.ProcStat(4, "bash", "S", 1, 4, 4, 34816, 5, 1)}.get(pid)), \
                mock.patch.object(gpuline.os, "stat", return_value=mock.Mock(st_uid=os.getuid())):
            self.assertEqual(gpuline.requeue_refusal(5, ["flock"]), "job of an interactive shell")
        self.assertIn("not in requeue_commands", gpuline.requeue_refusal(os.getpid(), ["flock"]) or "")


class ClaimTests(LineTestCase):
    def start_claim(self, since: float | None = None) -> tuple[subprocess.Popen, dict, Path]:
        env = self.env
        directory = env.base / "claim"
        directory.mkdir(exist_ok=True)
        token = os.urandom(8).hex()
        spec = {
            "schema": gpuline.CLAIM_SPEC_SCHEMA, "token": token, "unit": f"test-claim-{token}.service",
            "gpu": env.gpu.to_dict(), "line": env.cfg.to_dict(), "since_unix": since or time.time(),
            "secret": os.urandom(16).hex(), "socket_name": f"simq-test-claim-{token}",
            "status_path": str(directory / "status.json"), "control_path": str(directory / "control.json"),
            "hold_max_s": 60,
        }
        gpuline.write_json_atomic(directory / "spec.json", spec)
        gpuline.write_json_atomic(directory / "control.json",
                                  {"schema": gpuline.CLAIM_CONTROL_SCHEMA, "since_unix": spec["since_unix"],
                                   "release": False})
        proc = env.spawn([PYTHON, str(GPULINE), "claim", str(directory / "spec.json")],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        return proc, spec, directory

    def status(self, directory: Path) -> dict:
        return gpuline.read_json(directory / "status.json") or {}

    def handoff_spec(self, spec: dict, secret: str | None = None) -> tuple[dict, dict]:
        job = {"job_id": 99, "launch_token": "t" * 32, "unit": "simqt-job-99.service"}
        handoff = {"socket_name": spec["socket_name"], "secret": secret or spec["secret"],
                   "claim_unit": spec["unit"], "timeout_s": 5}
        return job, handoff

    def test_claim_wins_the_release_over_raw_waiters_and_hands_off_without_a_gap(self) -> None:
        env = self.env
        holder = env.raw("r0", 1.5)
        wait_for(lambda: env.holders(), message="raw holder")
        raws = []
        for index in range(3):
            raws.append(env.raw(f"raw{index}", 0.1))
            wait_for(lambda: len(env.waiters()) == index + 1, message=f"raw waiter {index}")
        claim, spec, directory = self.start_claim()
        wait_for(lambda: self.status(directory).get("phase") == "held", 30, "claim held")
        holder.wait(timeout=10)
        self.assertEqual(env.entries(), ["r0"])  # no raw waiter ran
        job, handoff = self.handoff_spec(spec)
        connection, received = jobwrap._receive_gpu_handoff(job, handoff)
        self.assertEqual(list(received), list(env.gpu.paths))
        jobwrap._reply(connection, {"taken": True})
        self.assertEqual(claim.wait(timeout=20), 0)
        self.assertEqual(self.status(directory)["phase"], "handed-off")
        self.assertEqual(self.status(directory)["handed_to"]["job_id"], 99)
        # the claim has exited; the job's descriptors keep both files locked
        time.sleep(0.5)
        self.assertEqual(env.entries(), ["r0"])
        self.assertNotEqual(subprocess.run(["flock", "-n", env.l1, "true"]).returncode, 0)
        for descriptor in received.values():
            os.close(descriptor)
        for proc in raws:
            proc.wait(timeout=20)
        self.assertEqual(sorted(env.entries()[1:]), ["raw0", "raw1", "raw2"])

    def test_wrong_secret_is_refused_and_a_refused_job_leaves_the_claim_holding(self) -> None:
        env = self.env
        claim, spec, directory = self.start_claim()
        wait_for(lambda: self.status(directory).get("phase") == "held", 30, "claim held")
        job, handoff = self.handoff_spec(spec, secret="0" * 32)
        with self.assertRaises(RuntimeError):
            jobwrap._receive_gpu_handoff(job, handoff)
        job, handoff = self.handoff_spec(spec)
        connection, received = jobwrap._receive_gpu_handoff(job, handoff)
        for descriptor in received.values():
            os.close(descriptor)
        jobwrap._reply(connection, {"taken": False, "reason": "core lease busy"})
        time.sleep(0.5)
        self.assertIsNone(claim.poll())
        self.assertEqual(self.status(directory)["phase"], "held")
        self.assertNotEqual(subprocess.run(["flock", "-n", env.l1, "true"]).returncode, 0)
        gpuline.write_json_atomic(directory / "control.json",
                                  {"schema": gpuline.CLAIM_CONTROL_SCHEMA, "release": True})
        self.assertEqual(claim.wait(timeout=20), 0)
        self.assertEqual(self.status(directory)["phase"], "released")
        self.assertEqual(subprocess.run(["flock", "-n", env.l1, "true"]).returncode, 0)

    def test_shared_downgrade_after_handoff_is_atomic(self) -> None:
        env = self.env
        claim, spec, directory = self.start_claim()
        wait_for(lambda: self.status(directory).get("phase") == "held", 30, "claim held")
        job, handoff = self.handoff_spec(spec)
        connection, received = jobwrap._receive_gpu_handoff(job, handoff)
        import fcntl

        for descriptor in received.values():
            fcntl.flock(descriptor, fcntl.LOCK_SH | fcntl.LOCK_NB)
        jobwrap._reply(connection, {"taken": True})
        self.assertEqual(claim.wait(timeout=20), 0)
        self.assertEqual(subprocess.run(["flock", "-n", "-s", env.l1, "true"]).returncode, 0)
        self.assertNotEqual(subprocess.run(["flock", "-n", "-x", env.l1, "true"]).returncode, 0)
        for descriptor in received.values():
            os.close(descriptor)

    def test_stopped_claim_releases(self) -> None:
        env = self.env
        holder = env.raw("r0", 30)
        wait_for(lambda: env.holders(), message="raw holder")
        claim, spec, directory = self.start_claim()
        wait_for(lambda: self.status(directory).get("phase") == "queued", 30, "claim queued in the kernel")
        claim.send_signal(signal.SIGTERM)
        self.assertEqual(claim.wait(timeout=20), 0)
        self.assertEqual(self.status(directory)["phase"], "released")
        self.assertEqual(env.waiters(), [])
        self.assertEqual([t for t in env.tickets() if t["kind"] == "queue"], [])
        holder.kill()


if __name__ == "__main__":
    unittest.main()
