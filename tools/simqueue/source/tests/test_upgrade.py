"""upgrade.py install / rollback in a sandbox home.

Nothing real is touched: every path is under a temporary --home, and
--systemctl is a stub that runs the sandbox's simqueue.service ExecStart as a
plain background process (start/stop/show).  The job units are real (test
prefix, leased test CPUs): a 0.6.0 worker starts a job, the upgrade switches
to 0.6.1, which must adopt it; rollback switches back and 0.6.0 must adopt a
job started by 0.6.1.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sqlite3
import subprocess
import time
import unittest
from contextlib import closing
from pathlib import Path

from node_fixture import PYTHON_FROZEN, NodeFixture


RUNTIME = Path(__file__).resolve().parents[2]
UPGRADE = RUNTIME / "upgrade.py"
OLD_RUNTIME = Path("/home/treefor/projects/runtimes/simqueue-0.6.0-20261030")
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


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_upgrade_module():
    spec = importlib.util.spec_from_file_location("simq_upgrade", UPGRADE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class UpgradeHelperTests(unittest.TestCase):
    def test_config_comparison_ignores_only_the_061_keys(self) -> None:
        module = load_upgrade_module()
        old = json.loads((OLD_RUNTIME / "config" / "zenobia-node.json").read_text())
        new = json.loads((RUNTIME / "config" / "zenobia-node.json").read_text())
        self.assertEqual(module.comparable_config(old), module.comparable_config(new))
        live = Path("/home/treefor/.local/state/simqueue/node.json")
        if live.exists():
            self.assertEqual(module.comparable_config(json.loads(live.read_text())), module.comparable_config(new))
        changed = json.loads(json.dumps(new))
        changed["partitions"][0]["cpus"] = "0-2"
        self.assertNotEqual(module.comparable_config(old), module.comparable_config(changed))

    def test_wrapper_runtime_is_read_from_the_exec_line(self) -> None:
        module = load_upgrade_module()
        live = Path("/home/treefor/.local/bin/simq")
        if live.exists():
            self.assertEqual(module.wrapper_runtime(live.read_text()), OLD_RUNTIME)
        self.assertIsNone(module.wrapper_runtime("#!/bin/sh\nexit 0\n"))


class UpgradeTests(NodeFixture, unittest.TestCase):
    def setUp(self) -> None:
        super().setUp()
        if not (OLD_RUNTIME / "source" / "simq").exists():
            self.skipTest("0.6.0 runtime not present")
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
        (self.sandbox_state / "node.json").write_text(json.dumps(old_doc, indent=2))
        new_doc = json.loads(json.dumps(old_doc))
        new_doc["borrow"] = {"slice": self.prefix + "b.slice", "idle_window_s": 3, "busy_window_s": 1}
        for part in new_doc["partitions"]:
            part["lend"] = True
        self.new_config = self.base / "node-061.json"
        self.new_config.write_text(json.dumps(new_doc, indent=2))
        module = load_upgrade_module()
        args = type("Args", (), {"execute": False, "home": str(self.home), "runtime": str(OLD_RUNTIME),
                                 "old_runtime": str(OLD_RUNTIME), "project_root": str(self.projects),
                                 "service": "simqueue.service", "systemctl": str(stub)})()
        sandbox = module.Upgrade(args)
        # the 0.6.0 installation as switchover.py left it
        old_wrapper = sandbox.wrapper_text(OLD_RUNTIME, "0.6").replace(
            "# borrow jobs run in simqborrow.slice at idle CPU priority.\n", "")
        (self.bin / "simq").write_text(old_wrapper)
        (self.bin / "simq").chmod(0o755)
        (self.bin / "simqueue").symlink_to(OLD_RUNTIME / "simq")
        (self.units / "simqueue.service").write_text(
            sandbox.unit_text(OLD_RUNTIME).replace("SimQueue 0.6.1", "SimQueue 0.6"))
        self.originals = {path: digest(path) for path in (self.bin / "simq", self.units / "simqueue.service",
                                                           self.sandbox_state / "node.json")}
        self.sandbox_env = {k: v for k, v in os.environ.items() if not k.startswith("SIMQUEUE_")}
        self.sandbox_env.update({"STUB_DIR": str(self.stub_dir), "STUB_HOME": str(self.home)})

    def tearDown(self) -> None:
        self.stub("stop")
        super().tearDown()
        subprocess.run(["systemctl", "--user", "revert", self.prefix + "b.slice"], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=False)
        subprocess.run(["systemctl", "--user", "stop", self.prefix + "b.slice"], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=False)

    def stub(self, *arguments: str) -> None:
        subprocess.run([str(self.stub_dir / "systemctl"), "--user", *arguments], env=self.sandbox_env, check=False)

    def sandbox_simq(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        env = dict(self.sandbox_env, SIMQUEUE_NO_AUTOSTART="1")
        return subprocess.run(["/bin/sh", str(self.bin / "simq"), *arguments], env=env, cwd=self.projects, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120, check=False)

    def upgrade(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        command = [PYTHON_FROZEN, str(UPGRADE), *arguments, "--home", str(self.home), "--project-root",
                   str(self.projects), "--systemctl", str(self.stub_dir / "systemctl"), "--runtime", str(RUNTIME),
                   "--old-runtime", str(OLD_RUNTIME), "--adoption-timeout", "60"]
        if not (RUNTIME / "sealed-files.json").exists() or os.environ.get("SIMQ_UPGRADE_UNSEALED"):
            command.append("--allow-unsealed")
        return subprocess.run(command, env=self.sandbox_env, text=True, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, timeout=600, check=False)

    def db(self, query: str, *params):
        with closing(sqlite3.connect(f"file:{self.sandbox_state / 'queue.sqlite3'}?mode=ro", uri=True)) as c:
            return c.execute(query, params).fetchall()

    def sandbox_state_of(self, job_id: int) -> str:
        return self.db("SELECT state FROM jobs WHERE id = ?", job_id)[0][0]

    def submit_old(self, seconds: float) -> int:
        result = self.sandbox_simq("submit", "--json", "-p", "pa", "-c", "1", *self.probe_command(seconds))
        self.assertEqual(result.returncode, 0, result.stderr)
        return int(json.loads(result.stdout)["job_id"])

    def worker_log(self) -> str:
        path = self.sandbox_state / "worker.log"
        return path.read_text() if path.exists() else ""

    def test_dry_run_changes_nothing(self) -> None:
        self.submit_old(1)  # creates the schema-7 database with the 0.6.0 code
        before = digest(self.sandbox_state / "queue.sqlite3")
        self.stub("start")
        result = self.upgrade("install", "--config", str(self.new_config))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("[dry-run]", result.stdout)
        self.stub("stop")
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value)
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 7)
        self.assertNotEqual(before, "")

    def test_install_adopts_running_060_jobs_and_rollback_adopts_061_jobs(self) -> None:
        self.stub("start")
        old_job = self.submit_old(15)
        self.wait_until(lambda: self.sandbox_state_of(old_job) == "running", message="0.6.0 job running")
        result = self.upgrade("install", "--execute", "--config", str(self.new_config))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        backup = Path(next(line.split("--backup ")[1].split()[0] for line in result.stdout.splitlines()
                           if "Rollback:" in line))
        receipt = json.loads((backup / "install-receipt.json").read_text())
        self.assertTrue(receipt["verification"]["ok"])
        self.assertTrue(receipt["verification"]["jobs"][str(old_job)]["adopted"])
        self.assertTrue(receipt["migration"]["rows_unchanged"])
        self.assertIn(f"simqueue 0.6.1 from {RUNTIME}", self.worker_log())
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 8)
        self.assertIn(str(RUNTIME), (self.bin / "simq").read_text())
        self.assertEqual(os.readlink(self.bin / "simqueue"), str(RUNTIME / "simq"))
        self.assertTrue((self.units / "simqborrow.slice").exists())
        self.assertEqual(json.loads((self.sandbox_state / "node.json").read_text())["borrow"]["slice"],
                         self.prefix + "b.slice")
        # the 0.6.0 job keeps running under its 0.6.0 jobwrap and is recorded by 0.6.1
        self.wait_until(lambda: self.sandbox_state_of(old_job) == "succeeded", timeout=60, message="old job done")
        # 0.6.1 work: an owner job and a borrow job
        owner = self.sandbox_simq("submit", "--json", "-p", "pa", "-c", "1", *self.probe_command(20))
        owner_job = int(json.loads(owner.stdout)["job_id"])
        borrow = self.sandbox_simq("submit", "--json", "--borrow", "-p", "pb", "-c", "1", *self.probe_command(60))
        self.assertEqual(borrow.returncode, 0, borrow.stderr)
        borrow_job = int(json.loads(borrow.stdout)["job_id"])
        self.wait_until(lambda: self.sandbox_state_of(owner_job) == "running"
                        and self.sandbox_state_of(borrow_job) == "running", message="0.6.1 jobs running")
        refused = self.upgrade("rollback", "--execute", "--backup", str(backup))
        self.assertEqual(refused.returncode, 1, refused.stdout)
        self.assertIn("borrow job(s) queued or active", refused.stdout)
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 8)
        rolled = self.upgrade("rollback", "--execute", "--backup", str(backup), "--cancel-borrow-jobs")
        self.assertEqual(rolled.returncode, 0, rolled.stdout + rolled.stderr)
        self.assertEqual(self.sandbox_state_of(borrow_job), "canceled")
        rollback_receipt = json.loads((backup / "rollback-receipt.json").read_text())
        self.assertTrue(rollback_receipt["verification"]["jobs"][str(owner_job)]["adopted"])
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 7)
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value, f"{path} not restored byte for byte")
        self.assertEqual(os.readlink(self.bin / "simqueue"), str(OLD_RUNTIME / "simq"))
        self.assertFalse((self.units / "simqborrow.slice").exists())
        log = self.worker_log()
        self.assertIn(f"simqueue 0.6.0 from {OLD_RUNTIME}", log)
        self.assertIn(f"adopted job {owner_job} ", log.rsplit("simqueue 0.6.0 from", 1)[1])
        self.wait_until(lambda: self.sandbox_state_of(owner_job) == "succeeded", timeout=60,
                        message="0.6.1 job recorded by 0.6.0")
        self.assertTrue((backup / "queue-0.6.1-at-rollback.sqlite3").exists())
        # the 0.6.1 columns stay (ignored by 0.6.0)
        self.assertEqual(self.db("SELECT borrow FROM jobs WHERE id = ?", borrow_job)[0][0], 1)

    def test_failed_verification_rolls_back_automatically(self) -> None:
        self.stub("start")
        job = self.submit_old(20)
        self.wait_until(lambda: self.sandbox_state_of(job) == "running", message="0.6.0 job running")
        command_timeout = "0"
        result = subprocess.run(
            [PYTHON_FROZEN, str(UPGRADE), "install", "--execute", "--config", str(self.new_config), "--home",
             str(self.home), "--project-root", str(self.projects), "--systemctl", str(self.stub_dir / "systemctl"),
             "--runtime", str(RUNTIME), "--old-runtime", str(OLD_RUNTIME), "--adoption-timeout", command_timeout,
             *(["--allow-unsealed"] if not (RUNTIME / "sealed-files.json").exists()
               or os.environ.get("SIMQ_UPGRADE_UNSEALED") else [])],
            env=self.sandbox_env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=600, check=False)
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn("rolling back automatically", result.stdout)
        self.assertIn("rolled back to SimQueue 0.6.0", result.stdout)
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value)
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 7)
        self.wait_until(lambda: self.sandbox_state_of(job) == "succeeded", timeout=60, message="job recorded by 0.6.0")
        log = self.worker_log()
        self.assertIn(f"adopted job {job} ", log.rsplit(f"simqueue 0.6.0 from {OLD_RUNTIME}", 1)[1])


if __name__ == "__main__":
    unittest.main()
