"""upgrade.py install / rollback (0.8.1 to 0.8.2) in a sandbox home.

Nothing real is touched: every path is under a temporary --home, and
--systemctl is a stub that runs the sandbox's simqueue.service ExecStart as a
plain background process (start/stop/show).  The job units and GPU claim
units are real (test prefix, test CPUs, private lease files); GitHub is a
local fake.  A 0.8.1 worker runs a job and queues a second one; the install
keeps the database (schema 10, every row unchanged) and the live node.json
(byte for byte), the 0.8.2 worker adopts the running job and runs the queued
one.  Rollback, while a 0.8.2 --smt job runs on its core and sibling, a
second --smt job waits behind it, a GPU job waits in the GPU line and a
job's GitHub post is due, leaves the schema at 10 and restores every file
byte for byte; 0.8.1 adopts the running --smt job and the GPU claim unit,
makes the due post, and starts the queued --smt job as an exclusive-timing
job on its core only (the safe degradation).
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import signal
import sqlite3
import subprocess
import sys
import time
import unittest
from contextlib import closing
from pathlib import Path

from fake_github import FakeGitHub
from node_fixture import FAKE_GPU, PYTHON_FROZEN, NodeFixture
from simqueue import __version__ as _VERSION

if _VERSION != "0.8.2":
    raise unittest.SkipTest(f"upgrade.py and these tests are for 0.8.1 -> 0.8.2; {_VERSION} needs its own")


RUNTIME = Path(__file__).resolve().parents[2]
UPGRADE = RUNTIME / "upgrade.py"
OLD_RUNTIME = Path("/home/treefor/projects/runtimes/simqueue-0.8.1-20261106")
TOKEN = "ghp_upgradeTestToken_0000"
SHA = "fedcba9876543210fedcba9876543210fedcba98"
STUB = r'''#!%(python)s
import os, shlex, signal, subprocess, sys, time
from pathlib import Path
stub = Path(os.environ["STUB_DIR"]); home = Path(os.environ["STUB_HOME"]); service = "simqueue.service"
args = sys.argv[1:]
with open(stub / "calls.log", "a") as log:
    log.write(" ".join(args) + "\n")
if args and args[0] == "--user":
    args = args[1:]
unit_file = home / ".config" / "systemd" / "user" / service
pidfile = stub / "worker.pid"

def alive():
    try:
        pid = int(pidfile.read_text())
        state = open(f"/proc/{pid}/stat").read().rsplit(")", 1)[1].split()[0]
        return pid if state != "Z" else 0
    except (OSError, ValueError):
        return 0

def directives():
    result = {"Environment": []}
    for line in unit_file.read_text().splitlines():
        key, _, value = line.partition("=")
        if key == "Environment":
            result["Environment"].append(value)
        elif key in ("ExecStart", "WorkingDirectory"):
            result[key] = value
    return result

verb = args[0] if args else ""
if verb == "show":
    unit = next(a for a in args[1:] if not a.startswith("-"))
    prop = next(a.split("=", 1)[1] for a in args if a.startswith("--property="))
    if unit == service:
        if prop == "ActiveState":
            print("active" if alive() else "inactive")
        elif prop == "ExecStart":
            command = directives().get("ExecStart", "")
            print("{ path=" + command.split()[0] + " ; argv[]=" + command + " }")
        else:
            print("")
    else:
        print("idle" if prop == "CPUWeight" else "")
elif verb == "stop":
    pid = alive()
    if pid:
        os.kill(pid, signal.SIGTERM)
        deadline = time.time() + 30
        while alive() and time.time() < deadline:
            time.sleep(0.1)
elif verb == "start":
    if not alive():
        found = directives()
        env = dict(os.environ)
        for item in found["Environment"]:
            key, _, value = item.partition("=")
            env[key] = value
        env.pop("SIMQUEUE_NO_AUTOSTART", None)
        out = open(stub / "worker.out", "ab")
        proc = subprocess.Popen(shlex.split(found["ExecStart"]), env=env, cwd=found.get("WorkingDirectory"),
                                stdin=subprocess.DEVNULL, stdout=out, stderr=subprocess.STDOUT, start_new_session=True)
        pidfile.write_text(str(proc.pid))
sys.exit(0)
'''
UNIT_063 = """[Unit]
Description=SimQueue 0.6.3 node scheduler (jobs run in their own units)
Documentation=file://{runtime}/OPERATIONS.md

[Service]
Type=simple
WorkingDirectory={root}
ExecStart={bin}/simq worker --wait-for-lock
Environment=SIMQUEUE_PROJECT_ROOT={root}
Environment=SIMQUEUE_STATE_DIR={state}
Environment=SIMQUEUE_CONFIG={state}/node.json
# No CPUAffinity: the scheduler places each job unit itself; the worker
# only polls and must not be confined to one partition.
Nice=5
Restart=on-failure
RestartSec=2
# Job units and GPU claim units are separate units: stopping this
# service never stops them.
KillMode=mixed
TimeoutStopSec=30
UMask=0077

[Install]
WantedBy=default.target
"""
SLICE = (
    "[Unit]\n"
    "Description=SimQueue borrow jobs (preemptible; idle CPU priority)\n"
    "\n"
    "[Slice]\n"
    "CPUWeight=idle\n"
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_upgrade_module():
    spec = importlib.util.spec_from_file_location("simq_upgrade", UPGRADE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def config_082(document: dict) -> dict:
    """A node configuration that --config may install: the live one, other notes."""
    new = json.loads(json.dumps(document))
    new["notes"] = ["0.8.2 sandbox"]
    return new


class UpgradeHelperTests(unittest.TestCase):
    def test_config_comparison_ignores_only_notes(self) -> None:
        module = load_upgrade_module()
        old = json.loads((OLD_RUNTIME / "config" / "zenobia-16core-node.json").read_text())
        new = json.loads((RUNTIME / "config" / "zenobia-16core-node.json").read_text())
        # 0.8.2 adds only a note to the template
        self.assertEqual(module.config_differences(module.comparable_config(old), module.comparable_config(new)), [])
        self.assertEqual(new["scheduler"]["reservation_depth"], 8)
        live = Path("/home/treefor/.local/state/simqueue/node.json")
        if live.exists() and module.wrapper_runtime(Path("/home/treefor/.local/bin/simq").read_text()) == OLD_RUNTIME:
            # the live file differs only where it was edited by hand (lend_siblings true, 2026-11-06); the
            # 0.8.2 install keeps it as it is
            self.assertLessEqual(set(module.config_differences(module.comparable_config(json.loads(live.read_text())),
                                                               module.comparable_config(new))),
                                 {"borrow.lend_siblings"})
        for change, key in (
            (lambda d: d["scheduler"].__setitem__("reservation_depth", 4), "scheduler.reservation_depth"),
            (lambda d: d["borrow"].__setitem__("idle_window_s", 5), "borrow.idle_window_s"),
            (lambda d: d.__setitem__("memory_schedulable", "30G"), "memory_schedulable"),
            (lambda d: d["partitions"][0].__setitem__("max_mem", "26G"), "partitions[0].max_mem"),
            (lambda d: d["memory"].__setitem__("auto_default", False), "memory.auto_default"),
            (lambda d: d["borrow"].__setitem__("lend_reserved", False), "borrow.lend_reserved"),
            (lambda d: d["topology"].__setitem__("l3_placement", False), "topology.l3_placement"),
        ):
            changed = json.loads(json.dumps(new))
            change(changed)
            self.assertEqual(module.config_differences(module.comparable_config(new),
                                                       module.comparable_config(changed)), [key])

    def test_wrapper_and_unit_text(self) -> None:
        module = load_upgrade_module()
        live = Path("/home/treefor/.local/bin/simq")
        if live.exists():
            self.assertIn(module.wrapper_runtime(live.read_text()), (*module.OLD_RUNTIMES, RUNTIME))
        self.assertIsNone(module.wrapper_runtime("#!/bin/sh\nexit 0\n"))
        args = type("Args", (), {"execute": False, "home": "/home/x", "runtime": str(RUNTIME),
                                 "old_runtime": str(OLD_RUNTIME), "project_root": "/home/x/projects",
                                 "service": "simqueue.service", "systemctl": "systemctl"})()
        upgrade = module.Upgrade(args)
        unit = upgrade.unit_text(RUNTIME)
        for line in ("Nice=10", "CPUWeight=20", "IOSchedulingClass=idle", "MemoryHigh=96M", "MemoryMax=128M",
                     "Description=SimQueue 0.8.2", "ExecStart=/home/x/.local/bin/simq worker --wait-for-lock"):
            self.assertIn(line, unit)
        self.assertEqual(module.wrapper_runtime(upgrade.wrapper_text(RUNTIME, "0.8.2")), RUNTIME)
        if live.exists() and module.wrapper_runtime(live.read_text()) == OLD_RUNTIME:
            # the live wrapper is what this script writes for 0.8.1
            self.assertEqual(upgrade.wrapper_text(OLD_RUNTIME, "0.8.1").replace("/home/x", "/home/treefor"),
                             live.read_text())

    def test_line_state_path_follows_the_configuration(self) -> None:
        module = load_upgrade_module()
        new = json.loads((RUNTIME / "config" / "zenobia-node.json").read_text())
        self.assertEqual(module.line_state_path(new), new["gpu_line"]["state_path"])


class UpgradeSandbox(NodeFixture):
    """A sandbox home with the 0.8.1 installation as its upgrade.py left it."""

    def setUp(self) -> None:
        super().setUp()
        if not (OLD_RUNTIME / "source" / "simq").exists():
            self.skipTest("0.8.1 runtime not present")
        self.fake = FakeGitHub(TOKEN).start()
        self.home = self.base / "home"
        self.sandbox_state = self.home / ".local" / "state" / "simqueue"
        self.bin = self.home / ".local" / "bin"
        self.units = self.home / ".config" / "systemd" / "user"
        for directory in (self.sandbox_state, self.bin, self.units):
            directory.mkdir(parents=True)
        self.stub_dir = self.base / "stub"
        self.stub_dir.mkdir()
        stub = self.stub_dir / "systemctl"
        stub.write_text(STUB % {"python": PYTHON_FROZEN})
        stub.chmod(0o755)
        token = self.base / "token"
        token.write_text(TOKEN)
        old_doc = self.config_document()
        old_doc["borrow"] = {"slice": self.prefix + "b.slice", "idle_window_s": 3, "busy_window_s": 1}
        old_doc["gpu_line"] = {"poll_s": 0.05, "claim_cpus": list(self.cpus) or [0], "claim_hold_max_s": 60}
        old_doc["launch"]["idle_seconds"] = 5
        old_doc["memory"] = {"sample_seconds": 1, "record_seconds": 2}
        old_doc["github"] = {"api_base": self.fake.api_base, "token_file": str(token), "backoff_seconds": 0.5}
        (self.sandbox_state / "node.json").write_text(json.dumps(old_doc, indent=2))
        new_doc = config_082(old_doc)
        self.new_config = self.base / "node-082.json"
        self.new_config.write_text(json.dumps(new_doc, indent=2))
        self.l1 = self.locks / "gpu.lock"
        self.l2 = self.locks / f"gpu-{FAKE_GPU}.lock"
        for path in (self.l1, self.l2):
            path.touch()
        self.state_file = self.locks / "simq-gpu-line.json"
        module = load_upgrade_module()
        args = type("Args", (), {"execute": False, "home": str(self.home), "runtime": str(OLD_RUNTIME),
                                 "old_runtime": str(OLD_RUNTIME), "project_root": str(self.projects),
                                 "service": "simqueue.service", "systemctl": str(stub)})()
        sandbox = module.Upgrade(args)
        (self.bin / "simq").write_text(sandbox.wrapper_text(OLD_RUNTIME, "0.8.1"))
        (self.bin / "simq").chmod(0o755)
        (self.bin / "simqueue").symlink_to(OLD_RUNTIME / "simq")
        (self.units / "simqueue.service").write_text(sandbox.unit_text(OLD_RUNTIME, "0.8.1"))
        (self.units / "simqborrow.slice").write_text(SLICE)
        self.dropin = self.units / "simq-job-.service.d" / "50-oom-policy-stop.conf"
        self.dropin.parent.mkdir()
        self.dropin.write_text("[Service]\nOOMPolicy=stop\n")
        self.originals = {path: digest(path) for path in (self.bin / "simq", self.units / "simqueue.service",
                                                           self.units / "simqborrow.slice", self.dropin,
                                                           self.sandbox_state / "node.json")}
        self.sandbox_env = {k: v for k, v in os.environ.items() if not k.startswith("SIMQUEUE_")}
        self.sandbox_env.update({"STUB_DIR": str(self.stub_dir), "STUB_HOME": str(self.home)})
        self.pin = ["taskset", "-c", ",".join(str(c) for c in (self.cpus or (0,)))]

    def tearDown(self) -> None:
        self.stub("stop")
        for proc in self.holders:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except OSError:
                pass
        super().tearDown()
        self.fake.stop()
        subprocess.run(["systemctl", "--user", "stop", self.prefix + "b.slice"], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=False)

    def stub(self, *arguments: str) -> None:
        subprocess.run([str(self.stub_dir / "systemctl"), "--user", *arguments], env=self.sandbox_env, check=False)

    def sandbox_simq(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        env = dict(self.sandbox_env, SIMQUEUE_NO_AUTOSTART="1")
        return subprocess.run(["/bin/sh", str(self.bin / "simq"), *arguments], env=env, cwd=self.projects, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120, check=False)

    def upgrade(self, *arguments: str, timeout: str = "60") -> subprocess.CompletedProcess[str]:
        command = [PYTHON_FROZEN, str(UPGRADE), *arguments, "--home", str(self.home), "--project-root",
                   str(self.projects), "--systemctl", str(self.stub_dir / "systemctl"), "--runtime", str(RUNTIME),
                   "--old-runtime", str(OLD_RUNTIME), "--adoption-timeout", timeout]
        if not (RUNTIME / "sealed-files.json").exists() or os.environ.get("SIMQ_UPGRADE_UNSEALED"):
            command.append("--allow-unsealed")
        return subprocess.run(command, env=self.sandbox_env, text=True, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, timeout=600, check=False)

    def db(self, query: str, *params):
        with closing(sqlite3.connect(f"file:{self.sandbox_state / 'queue.sqlite3'}?mode=ro", uri=True)) as c:
            return c.execute(query, params).fetchall()

    def sandbox_state_of(self, job_id: int) -> str:
        return self.db("SELECT state FROM jobs WHERE id = ?", job_id)[0][0]

    def submit(self, *arguments: str) -> int:
        result = self.sandbox_simq("submit", "--json", *arguments)
        self.assertEqual(result.returncode, 0, result.stderr)
        return int(json.loads(result.stdout)["job_id"])

    def worker_log(self) -> str:
        path = self.sandbox_state / "worker.log"
        return path.read_text() if path.exists() else ""

    def claim_statuses(self) -> list[dict]:
        root = self.sandbox_state / "claims"
        found = []
        if root.exists():
            for directory in sorted(root.iterdir()):
                try:
                    found.append(json.loads((directory / "status.json").read_text()))
                except (OSError, ValueError):
                    pass
        return found

    def active_claim_units(self) -> list[str]:
        result = subprocess.run(["systemctl", "--user", "list-units", f"{self.prefix}-gpuclaim-*", "--no-legend",
                                 "--plain", "--state=active,activating,deactivating"], text=True,
                                stdout=subprocess.PIPE, check=False)
        return [line.split()[0] for line in result.stdout.splitlines() if line.strip()]

    def raw_holder(self) -> subprocess.Popen:
        proc = subprocess.Popen([*self.pin, "flock", str(self.l1), "flock", str(self.l2), "sleep", "600"],
                                stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                start_new_session=True)
        self.holders.append(proc)
        return proc


class UpgradePreflightTests(UpgradeSandbox, unittest.TestCase):
    """Dry runs only: nothing is started, no unit and no lease is needed."""

    needs_units = False

    def test_dry_run_changes_nothing(self) -> None:
        self.submit("-p", "pb", "-c", "1", *self.probe_command(1))  # creates the schema-10 database with 0.8.1
        self.submit("-p", "pb", "-c", "1", "--exclusive-timing", *self.probe_command(1))
        before = digest(self.sandbox_state / "queue.sqlite3")
        result = self.upgrade("install")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("[dry-run]", result.stdout)
        self.assertIn(f"node config: the live {self.sandbox_state / 'node.json'} is kept as it is (0.8.2 adds no "
                      "key)", result.stdout)
        self.assertIn("keeping", result.stdout)
        self.assertNotIn("would write " + str(self.sandbox_state / "node.json"), result.stdout)
        self.assertIn("database rehearsal on a copy opened by 0.8.2: schema 10 -> 10 (no schema change), 2 rows",
                      result.stdout)
        self.assertIn("--smt on this node: cores ", result.stdout)
        self.assertIn("queued and active jobs by exclusive_timing: 0: 1, 1: 1 (0 and 1 as under 0.8.1; 0.8.2 reads "
                      "bit 2 as --smt; none has it)", result.stdout)
        self.assertIn("no changes made; re-run with --execute", result.stdout)
        self.assertIn("interim drop-in", result.stdout)
        self.assertIn("would remove the interim drop-in", result.stdout)
        self.assertNotIn("PROBLEM", result.stdout)
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value)
        self.assertEqual(digest(self.sandbox_state / "queue.sqlite3"), before)
        self.assertFalse((self.home / ".local" / "state" / "simqueue-backups").exists())
        # with --config the given file is compared and would be written
        result = self.upgrade("install", "--config", str(self.new_config))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(f"node config {self.new_config}: equal to the live one apart from notes; the install writes "
                      "it", result.stdout)
        self.assertIn("would write " + str(self.sandbox_state / "node.json"), result.stdout)
        self.assertNotIn("PROBLEM", result.stdout)
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value)

    def test_preflight_refuses_other_config_changes(self) -> None:
        self.submit("-p", "pb", "-c", "1", *self.probe_command(1))
        good = json.loads(self.new_config.read_text())
        for change, message in ((lambda d: d["borrow"].__setitem__("idle_window_s", 9), "borrow.idle_window_s"),
                                (lambda d: d["partitions"][0].__setitem__("max_mem", "1G"),
                                 "partitions[0].max_mem")):
            with self.subTest(message=message):
                document = json.loads(json.dumps(good))
                change(document)
                path = self.base / "candidate.json"
                path.write_text(json.dumps(document))
                result = self.upgrade("install", "--config", str(path))
                self.assertEqual(result.returncode, 1, result.stdout)
                self.assertIn(message, result.stdout)
                for file, value in self.originals.items():
                    self.assertEqual(digest(file), value)


class UpgradeTests(UpgradeSandbox, unittest.TestCase):
    def test_install_adopts_and_rollback_restores(self) -> None:
        self.stub("start")
        width_b = len(self.part_b)
        old_job = self.submit("-p", "pb", "-c", str(width_b), *self.probe_command(15))
        self.wait_until(lambda: self.sandbox_state_of(old_job) == "running", message="0.8.1 job running")
        queued_b = self.submit("-p", "pb", "-c", "1", *self.probe_command(1))  # pb is full: it waits
        time.sleep(1.5)
        self.assertEqual(self.sandbox_state_of(queued_b), "queued")
        rows_before = self.db("SELECT id, name, argv_json, env_json, cwd, mem_bytes FROM jobs ORDER BY id")
        result = self.upgrade("install", "--execute")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        backup = Path(next(line.split("--backup ")[1].split()[0] for line in result.stdout.splitlines()
                           if "Rollback:" in line))
        self.assertTrue(backup.name.startswith("0.8.1-to-0.8.2-"))
        receipt = json.loads((backup / "install-receipt.json").read_text())
        verification = receipt["verification"]
        self.assertTrue(verification["ok"])
        self.assertTrue(verification["jobs"][str(old_job)]["adopted"])
        self.assertTrue(verification["wake_fifo"])
        self.assertEqual((receipt["database"]["schema"], receipt["database"]["added_columns"]), (10, []))
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 10)
        self.assertEqual(self.db("SELECT id, name, argv_json, env_json, cwd, mem_bytes FROM jobs WHERE id <= ? "
                                 "ORDER BY id", queued_b), rows_before)
        log = self.worker_log()
        self.assertIn(f"simqueue 0.8.2 from {RUNTIME}", log)
        self.assertIn("SimQueue 0.8.2", (self.units / "simqueue.service").read_text())
        self.assertFalse(self.dropin.exists(), "the interim OOMPolicy drop-in is removed by the install")
        node_json = self.sandbox_state / "node.json"
        self.assertEqual(digest(node_json), self.originals[node_json], "the install keeps node.json byte for byte")
        self.assertTrue(json.loads((backup / "manifest.json").read_text())["node_config_installed"] is False)
        self.wait_until(lambda: self.sandbox_state_of(queued_b) == "succeeded", timeout=60, message="queued job")
        self.wait_until(lambda: self.sandbox_state_of(old_job) == "succeeded", timeout=60, message="old job done")
        peak, source = self.db("SELECT mem_peak_bytes, mem_peak_source FROM jobs WHERE id = ?", old_job)[0]
        self.assertTrue(peak and source == "unit", (peak, source))  # the 0.8.1 wrapper records its unit's peak
        # 0.8.2 work at the rollback: a GitHub post due (GitHub down), a running --smt owner job with a second
        # one queued behind it, a GPU job in the line
        self.fake.down = True
        reported = self.submit("-p", "pa", "-c", "1", "--github-status", f"o/r@{SHA}:simq/upgrade",
                               *self.probe_command(1))
        self.wait_until(lambda: self.sandbox_state_of(reported) == "succeeded", timeout=60, message="reported job")
        self.wait_until(lambda: self.db("SELECT github_attempts FROM jobs WHERE id = ?", reported)[0][0] >= 1,
                        timeout=30, message="GitHub post failed while down")
        from simqueue.nodeconfig import load_node_config

        core = self.part_a[0]
        with_siblings, _ = load_node_config(self.sandbox_state / "node.json").smt_affinity((core,))
        self.assertGreater(len(with_siblings), 1, "the test core has no SMT sibling")
        owner_job = self.submit("-p", "pa", "-c", "1", "--smt", *self.probe_command(60))
        holder = self.raw_holder()
        gpu_job = self.submit("-p", "pb", "--gpus", "1", *self.probe_command(1))
        self.wait_until(lambda: self.sandbox_state_of(owner_job) == "running", message="0.8.2 --smt job running")
        queued_smt = self.submit("-p", "pa", "-c", "1", "--smt", *self.probe_command(1))  # pa is full: it waits
        self.assertEqual(self.db("SELECT exclusive_timing FROM jobs WHERE id = ?", queued_smt)[0][0], 3)
        self.wait_until(lambda: any(s.get("phase") == "queued" for s in self.claim_statuses()), timeout=30,
                        message="GPU claim queued in the kernel")
        claims = set(self.active_claim_units())
        self.assertTrue(claims)
        rolled = self.upgrade("rollback", "--execute", "--backup", str(backup))
        self.assertEqual(rolled.returncode, 0, rolled.stdout + rolled.stderr)
        self.assertIn(f"NOTE: job(s) [{reported}] still have GitHub posts due; 0.8.1 posts them", rolled.stdout)
        self.assertIn("database schema stays 10 (0.8.2 changes no schema)", rolled.stdout)
        self.assertIn("node.json left as it is (the 0.8.2 install did not change it)", rolled.stdout)
        self.assertIn(f"NOTE: job(s) [{owner_job}, {queued_smt}] were submitted with --smt: 0.8.1 reads them as "
                      "--exclusive-timing jobs", rolled.stdout)
        rollback_receipt = json.loads((backup / "rollback-receipt.json").read_text())
        self.assertTrue(rollback_receipt["verification"]["ok"])
        self.assertTrue(rollback_receipt["verification"]["jobs"][str(owner_job)]["adopted"])
        for unit_name in claims:
            self.assertTrue(rollback_receipt["verification"]["claims"][unit_name]["adopted"], unit_name)
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 10)
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value, f"{path} not restored byte for byte")
        self.assertEqual(os.readlink(self.bin / "simqueue"), str(OLD_RUNTIME / "simq"))
        after = self.worker_log().rsplit(f"simqueue 0.8.1 from {OLD_RUNTIME}", 1)[1]
        self.assertIn(f"adopted job {owner_job} ", after)
        self.fake.down = False  # 0.8.1 makes the due post
        self.wait_until(lambda: self.fake.statuses_for("simq/upgrade"), timeout=60, message="due post sent")
        os.killpg(holder.pid, signal.SIGKILL)
        holder.wait(timeout=10)
        self.wait_until(lambda: self.sandbox_state_of(gpu_job) == "succeeded", timeout=90,
                        message="GPU job run by 0.8.1")
        self.wait_until(lambda: self.sandbox_state_of(owner_job) == "succeeded", timeout=120,
                        message="0.8.2 --smt job recorded by 0.8.1")
        self.wait_until(lambda: self.sandbox_state_of(queued_smt) == "succeeded", timeout=120,
                        message="queued --smt job run by 0.8.1")
        records = self.read_intervals()
        self.assertEqual(records[owner_job]["affinity"], list(with_siblings))  # started by 0.8.2: core + sibling
        self.assertEqual(records[queued_smt]["affinity"], [core])  # started by 0.8.1: the core only
        shown = self.sandbox_simq("show", str(queued_smt)).stdout  # 0.8.1's show: an exclusive-timing job
        self.assertIn("exclusive timing", shown)
        self.assertIn("Exclusive timing", shown)

    def test_install_refuses_while_a_job_is_starting(self) -> None:
        job = self.submit("-p", "pa", "-c", "1", *self.probe_command(1))
        with closing(sqlite3.connect(self.sandbox_state / "queue.sqlite3")) as connection:
            connection.execute("UPDATE jobs SET state = 'starting' WHERE id = ?", (job,))
            connection.commit()
        result = subprocess.run(
            [PYTHON_FROZEN, str(UPGRADE), "install", "--execute", "--config", str(self.new_config), "--home",
             str(self.home), "--project-root", str(self.projects), "--systemctl", str(self.stub_dir / "systemctl"),
             "--runtime", str(RUNTIME), "--old-runtime", str(OLD_RUNTIME), "--starting-timeout", "1",
             *(["--allow-unsealed"] if not (RUNTIME / "sealed-files.json").exists()
               or os.environ.get("SIMQ_UPGRADE_UNSEALED") else [])],
            env=self.sandbox_env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=180, check=False)
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn(f"job(s) [{job}] still starting", result.stdout)
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value)
        self.assertFalse((self.home / ".local" / "state" / "simqueue-backups").exists())

    def test_failed_verification_rolls_back_automatically(self) -> None:
        self.stub("start")
        job = self.submit("-p", "pa", "-c", "1", *self.probe_command(20))
        self.wait_until(lambda: self.sandbox_state_of(job) == "running", message="0.8.1 job running")
        result = self.upgrade("install", "--execute", "--config", str(self.new_config), timeout="0")
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn("rolling back automatically", result.stdout)
        self.assertIn("rolled back to SimQueue 0.8.1", result.stdout)
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value)
        self.wait_until(lambda: self.sandbox_state_of(job) == "succeeded", timeout=60, message="job recorded by 0.8.1")
        log = self.worker_log()
        self.assertIn(f"adopted job {job} ", log.rsplit(f"simqueue 0.8.1 from {OLD_RUNTIME}", 1)[1])


if __name__ == "__main__":
    unittest.main()
