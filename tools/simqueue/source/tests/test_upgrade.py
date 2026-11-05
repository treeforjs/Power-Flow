"""upgrade.py install / rollback (0.6.3 to 0.7.0) in a sandbox home.

Nothing real is touched: every path is under a temporary --home, and
--systemctl is a stub that runs the sandbox's simqueue.service ExecStart as a
plain background process (start/stop/show).  The job units and GPU claim
units are real (test prefix, test CPUs, private lease files); GitHub is a
local fake.  A 0.6.3 worker runs a job and queues a second one; the install
migrates the database to schema 9 (every schema-8 column of every row
unchanged), the 0.7.0 worker adopts the running job, samples its memory and
runs the queued one.  Rollback, while a 0.7.0 job runs, a GPU job waits in
the GPU line and a job's GitHub post is due, sets the schema back to 8 in
place, restores every file byte for byte, and 0.6.3 adopts the running job
and the GPU claim unit.
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


RUNTIME = Path(__file__).resolve().parents[2]
UPGRADE = RUNTIME / "upgrade.py"
OLD_RUNTIME = Path("/home/treefor/projects/runtimes/simqueue-0.6.3-20261105")
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


def config_070(document: dict) -> dict:
    """A 0.6.3 node configuration as 0.7.0 ships it: the new keys only."""
    new = json.loads(json.dumps(document))
    new.setdefault("launch", {})["idle_seconds"] = 5
    for part in new["partitions"]:
        part["min_mem"] = "128M"
    new["memory"] = {"sample_seconds": 1, "record_seconds": 2, "notes": "0.7.0"}
    new["notes"] = ["0.7.0 sandbox"]
    return new


class UpgradeHelperTests(unittest.TestCase):
    def test_config_comparison_ignores_only_the_070_keys(self) -> None:
        module = load_upgrade_module()
        old = json.loads((OLD_RUNTIME / "config" / "zenobia-node.json").read_text())
        new = json.loads((RUNTIME / "config" / "zenobia-node.json").read_text())
        self.assertEqual(module.config_differences(module.comparable_config(old), module.comparable_config(new)), [])
        live = Path("/home/treefor/.local/state/simqueue/node.json")
        if live.exists():
            self.assertEqual(module.comparable_config(json.loads(live.read_text())), module.comparable_config(new))
        old16 = json.loads((OLD_RUNTIME / "config" / "zenobia-16core-node.json").read_text())
        new16 = json.loads((RUNTIME / "config" / "zenobia-16core-node.json").read_text())
        self.assertEqual(module.config_differences(module.comparable_config(old16), module.comparable_config(new16)),
                         ["partitions[0].default_mem"])  # 18G in the sealed 0.6.3 file, 4G as on the 8-core node
        for change, key in (
            (lambda d: d["lock_files"].__setitem__("gpu_path", "/elsewhere/gpu-{uuid}.lock"), "lock_files.gpu_path"),
            (lambda d: d["borrow"].__setitem__("idle_window_s", 5), "borrow.idle_window_s"),
            (lambda d: d.__setitem__("memory_schedulable", "30G"), "memory_schedulable"),
            (lambda d: d["partitions"][0].__setitem__("default_mem", "18G"), "partitions[0].default_mem"),
            (lambda d: d["launch"].__setitem__("poll_seconds", 2), "launch.poll_seconds"),
        ):
            changed = json.loads(json.dumps(new))
            change(changed)
            self.assertEqual(module.config_differences(module.comparable_config(old),
                                                       module.comparable_config(changed)), [key])

    def test_wrapper_and_unit_text(self) -> None:
        module = load_upgrade_module()
        live = Path("/home/treefor/.local/bin/simq")
        if live.exists():
            self.assertIn(module.wrapper_runtime(live.read_text()), (OLD_RUNTIME, RUNTIME))
        self.assertIsNone(module.wrapper_runtime("#!/bin/sh\nexit 0\n"))
        args = type("Args", (), {"execute": False, "home": "/home/x", "runtime": str(RUNTIME),
                                 "old_runtime": str(OLD_RUNTIME), "project_root": "/home/x/projects",
                                 "service": "simqueue.service", "systemctl": "systemctl"})()
        upgrade = module.Upgrade(args)
        unit = upgrade.unit_text(RUNTIME)
        for line in ("Nice=10", "CPUWeight=20", "IOSchedulingClass=idle", "MemoryHigh=96M", "MemoryMax=128M",
                     "Description=SimQueue 0.7.0", "ExecStart=/home/x/.local/bin/simq worker --wait-for-lock"):
            self.assertIn(line, unit)
        self.assertEqual(module.wrapper_runtime(upgrade.wrapper_text(RUNTIME, "0.7.0")), RUNTIME)

    def test_line_state_path_follows_the_configuration(self) -> None:
        module = load_upgrade_module()
        new = json.loads((RUNTIME / "config" / "zenobia-node.json").read_text())
        self.assertEqual(module.line_state_path(new), new["gpu_line"]["state_path"])


class UpgradeSandbox(NodeFixture):
    """A sandbox home with the 0.6.3 installation as its upgrade.py left it."""

    def setUp(self) -> None:
        super().setUp()
        if not (OLD_RUNTIME / "source" / "simq").exists():
            self.skipTest("0.6.3 runtime not present")
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
        old_doc = self.config_document()
        old_doc["borrow"] = {"slice": self.prefix + "b.slice", "idle_window_s": 3, "busy_window_s": 1}
        old_doc["gpu_line"] = {"poll_s": 0.05, "claim_cpus": list(self.cpus) or [0], "claim_hold_max_s": 60}
        (self.sandbox_state / "node.json").write_text(json.dumps(old_doc, indent=2))
        token = self.base / "token"
        token.write_text(TOKEN)
        new_doc = config_070(old_doc)
        new_doc["github"] = {"api_base": self.fake.api_base, "token_file": str(token), "backoff_seconds": 0.5}
        self.new_config = self.base / "node-070.json"
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
        (self.bin / "simq").write_text(sandbox.wrapper_text(OLD_RUNTIME, "0.6.3"))
        (self.bin / "simq").chmod(0o755)
        (self.bin / "simqueue").symlink_to(OLD_RUNTIME / "simq")
        (self.units / "simqueue.service").write_text(UNIT_063.format(
            runtime=OLD_RUNTIME, root=self.projects, bin=self.bin, state=self.sandbox_state))
        (self.units / "simqborrow.slice").write_text(SLICE)
        self.originals = {path: digest(path) for path in (self.bin / "simq", self.units / "simqueue.service",
                                                           self.units / "simqborrow.slice",
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

    def test_dry_run_changes_nothing_and_rehearses_the_migration(self) -> None:
        self.submit("-p", "pb", "-c", "1", *self.probe_command(1))  # creates the schema-8 database with 0.6.3
        before = digest(self.sandbox_state / "queue.sqlite3")
        result = self.upgrade("install", "--config", str(self.new_config))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("[dry-run]", result.stdout)
        self.assertIn("node config: equal to the live one apart from the 0.7.0 keys", result.stdout)
        self.assertIn("migration rehearsal on a copy: schema 8 -> 9, 1 rows, every schema-8 column unchanged but "
                      "env_json", result.stdout)
        self.assertIn("environment scrub (counts only): ", result.stdout)
        self.assertIn("partition layout (unchanged by 0.7.0)", result.stdout)
        self.assertIn("migrating the database to schema 9", result.stdout)
        self.assertIn("no changes made; re-run with --execute", result.stdout)
        self.assertNotIn("PROBLEM", result.stdout)
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value)
        self.assertEqual(digest(self.sandbox_state / "queue.sqlite3"), before)
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 8)
        self.assertFalse((self.home / ".local" / "state" / "simqueue-backups").exists())

    def test_preflight_refuses_other_config_changes(self) -> None:
        self.submit("-p", "pb", "-c", "1", *self.probe_command(1))
        good = json.loads(self.new_config.read_text())
        cases = []
        other = json.loads(json.dumps(good))
        other["borrow"]["idle_window_s"] = 9
        cases.append((other, "borrow.idle_window_s"))
        layout = json.loads(json.dumps(good))
        layout["partitions"][0]["default_mem"] = "512M"
        cases.append((layout, "partitions[0].default_mem"))
        for document, message in cases:
            with self.subTest(message=message):
                path = self.base / "candidate.json"
                path.write_text(json.dumps(document))
                result = self.upgrade("install", "--config", str(path))
                self.assertEqual(result.returncode, 1, result.stdout)
                self.assertIn("PROBLEM", result.stdout)
                self.assertIn(message, result.stdout)
                for file, value in self.originals.items():
                    self.assertEqual(digest(file), value)


class UpgradeTests(UpgradeSandbox, unittest.TestCase):
    def test_install_migrates_adopts_and_rollback_downgrades_in_place(self) -> None:
        self.stub("start")
        width_b = len(self.part_b)
        old_job = self.submit("-p", "pb", "-c", str(width_b), *self.probe_command(15))
        self.wait_until(lambda: self.sandbox_state_of(old_job) == "running", message="0.6.3 job running")
        queued_b = self.submit("-p", "pb", "-c", "1", *self.probe_command(1))  # pb is full: it waits
        time.sleep(1.5)
        self.assertEqual(self.sandbox_state_of(queued_b), "queued")
        with closing(sqlite3.connect(f"file:{self.sandbox_state / 'queue.sqlite3'}?mode=ro", uri=True)) as c:
            columns_8 = [row[1] for row in c.execute("PRAGMA table_info(jobs)")]
        rows_before = self.db(f"SELECT {', '.join(columns_8)} FROM jobs ORDER BY id")
        result = self.upgrade("install", "--execute", "--config", str(self.new_config))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        backup = Path(next(line.split("--backup ")[1].split()[0] for line in result.stdout.splitlines()
                           if "Rollback:" in line))
        self.assertTrue(backup.name.startswith("0.6.3-to-0.7.0-"))
        receipt = json.loads((backup / "install-receipt.json").read_text())
        verification = receipt["verification"]
        self.assertTrue(verification["ok"])
        self.assertTrue(verification["jobs"][str(old_job)]["adopted"])
        self.assertTrue(verification["gpu_line_state"]["published"])
        self.assertTrue(verification["wake_fifo"])
        self.assertTrue(receipt["database"]["rows_unchanged"])
        self.assertIn("mem_peak_bytes", receipt["database"]["added_columns"])
        self.assertIn("github_json", receipt["database"]["added_columns"])
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 9)
        rows_after = self.db(f"SELECT {', '.join(columns_8)} FROM jobs WHERE id <= ? ORDER BY id", queued_b)
        for old, new in zip(rows_before, rows_after):
            if old[0] in (old_job, queued_b):
                continue  # their state moves on while the 0.7.0 worker runs them
            self.assertEqual(old, new)
        log = self.worker_log()
        self.assertIn(f"simqueue 0.7.0 from {RUNTIME}", log)
        self.assertIn("event-driven", log)
        self.assertIn(str(RUNTIME), (self.bin / "simq").read_text())
        self.assertEqual(os.readlink(self.bin / "simqueue"), str(RUNTIME / "simq"))
        unit = (self.units / "simqueue.service").read_text()
        self.assertIn("SimQueue 0.7.0", unit)
        self.assertIn("MemoryMax=128M", unit)
        self.assertIn("idle_seconds", (self.sandbox_state / "node.json").read_text())
        self.assertEqual(json.loads(self.state_file.read_text())["version"], "0.7.0")
        manifest = json.loads((backup / "manifest.json").read_text())
        self.assertEqual((manifest["from"], manifest["to"]), ("0.6.3", "0.7.0"))
        # the adopted 0.6.3 job is sampled by 0.7.0; the queued one runs and is measured by its unit
        self.wait_until(lambda: self.sandbox_state_of(queued_b) == "succeeded", timeout=60, message="queued job")
        self.wait_until(lambda: self.sandbox_state_of(old_job) == "succeeded", timeout=60, message="old job done")
        peak, source = self.db("SELECT mem_peak_bytes, mem_peak_source FROM jobs WHERE id = ?", old_job)[0]
        self.assertTrue(peak and source == "sampled", (peak, source))
        peak, source = self.db("SELECT mem_peak_bytes, mem_peak_source FROM jobs WHERE id = ?", queued_b)[0]
        self.assertTrue(peak and source == "unit", (peak, source))
        # 0.7.0 work at the rollback: a GitHub-reported job whose final post is due (GitHub down),
        # a running owner job and a GPU job in the GPU line behind a raw holder
        self.fake.down = True
        reported = self.submit("-p", "pa", "-c", "1", "--github-status", f"o/r@{SHA}:simq/upgrade",
                               *self.probe_command(1))
        self.wait_until(lambda: self.sandbox_state_of(reported) == "succeeded", timeout=60, message="reported job")
        self.wait_until(lambda: self.db("SELECT github_attempts FROM jobs WHERE id = ?", reported)[0][0] >= 1,
                        timeout=30, message="GitHub post failed while down")
        owner_job = self.submit("-p", "pa", "-c", "1", *self.probe_command(60))
        holder = self.raw_holder()
        gpu_job = self.submit("-p", "pb", "--gpus", "1", *self.probe_command(1))
        self.wait_until(lambda: self.sandbox_state_of(owner_job) == "running", message="0.7.0 owner job running")
        self.wait_until(lambda: any(s.get("phase") == "queued" for s in self.claim_statuses()), timeout=30,
                        message="GPU claim queued in the kernel")
        claims = set(self.active_claim_units())
        self.assertTrue(claims)
        self.assertEqual(self.sandbox_state_of(owner_job), "running")
        rolled = self.upgrade("rollback", "--execute", "--backup", str(backup))
        self.assertEqual(rolled.returncode, 0, rolled.stdout + rolled.stderr)
        self.assertIn(f"NOTE: job(s) [{reported}] still have GitHub posts due", rolled.stdout)
        self.assertIn("database schema set back to 8 in place (was 9)", rolled.stdout)
        rollback_receipt = json.loads((backup / "rollback-receipt.json").read_text())
        self.assertTrue(rollback_receipt["verification"]["ok"])
        self.assertTrue(rollback_receipt["verification"]["jobs"][str(owner_job)]["adopted"])
        for unit_name in claims:
            self.assertTrue(rollback_receipt["verification"]["claims"][unit_name]["adopted"], unit_name)
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 8)
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value, f"{path} not restored byte for byte")
        self.assertEqual(os.readlink(self.bin / "simqueue"), str(OLD_RUNTIME / "simq"))
        self.assertTrue((backup / "queue-0.7.0-at-rollback.sqlite3").exists())
        log = self.worker_log()
        after = log.rsplit(f"simqueue 0.6.3 from {OLD_RUNTIME}", 1)[1]
        self.assertIn(f"adopted job {owner_job} ", after)
        os.killpg(holder.pid, signal.SIGKILL)
        holder.wait(timeout=10)
        self.wait_until(lambda: self.sandbox_state_of(gpu_job) == "succeeded", timeout=90,
                        message="GPU job run by 0.6.3")
        self.wait_until(lambda: self.sandbox_state_of(owner_job) == "succeeded", timeout=120,
                        message="0.7.0 job recorded by 0.6.3")
        # installing again picks up where 0.7.0 stopped: the due post goes out
        self.fake.down = False
        again = self.upgrade("install", "--execute", "--config", str(self.new_config))
        self.assertEqual(again.returncode, 0, again.stdout + again.stderr)
        self.wait_until(lambda: self.fake.statuses_for("simq/upgrade"), timeout=30, message="due post sent")
        self.assertEqual(self.fake.statuses_for("simq/upgrade")[-1]["state"], "success")

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
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 8)
        self.assertFalse((self.home / ".local" / "state" / "simqueue-backups").exists())

    def test_failed_verification_rolls_back_automatically(self) -> None:
        self.stub("start")
        job = self.submit("-p", "pa", "-c", "1", *self.probe_command(20))
        self.wait_until(lambda: self.sandbox_state_of(job) == "running", message="0.6.3 job running")
        result = self.upgrade("install", "--execute", "--config", str(self.new_config), timeout="0")
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn("rolling back automatically", result.stdout)
        self.assertIn("rolled back to SimQueue 0.6.3", result.stdout)
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value)
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 8)
        self.wait_until(lambda: self.sandbox_state_of(job) == "succeeded", timeout=60, message="job recorded by 0.6.3")
        log = self.worker_log()
        self.assertIn(f"adopted job {job} ", log.rsplit(f"simqueue 0.6.3 from {OLD_RUNTIME}", 1)[1])


if __name__ == "__main__":
    unittest.main()
