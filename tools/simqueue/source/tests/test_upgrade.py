"""upgrade.py install / rollback (0.6.1 to 0.6.2) in a sandbox home.

Nothing real is touched: every path is under a temporary --home, and
--systemctl is a stub that runs the sandbox's simqueue.service ExecStart as a
plain background process (start/stop/show).  The job units and GPU claim
units are real (test prefix, leased test CPUs, private lease files): a 0.6.1
worker starts a job, the upgrade switches to 0.6.2, which must adopt it and
publish the GPU line state; rollback switches back while a 0.6.2 job runs
and a queued GPU job waits in the GPU line: the claim unit is stopped, 0.6.1
adopts the running job and later starts the GPU job itself.
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

from node_fixture import FAKE_GPU, PYTHON_FROZEN, NodeFixture


RUNTIME = Path(__file__).resolve().parents[2]
UPGRADE = RUNTIME / "upgrade.py"
OLD_RUNTIME = Path("/home/treefor/projects/runtimes/simqueue-0.6.1-20261031")
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
SLICE_061 = (
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


class UpgradeHelperTests(unittest.TestCase):
    def test_config_comparison_ignores_only_the_062_keys(self) -> None:
        module = load_upgrade_module()
        old = json.loads((OLD_RUNTIME / "config" / "zenobia-node.json").read_text())
        new = json.loads((RUNTIME / "config" / "zenobia-node.json").read_text())
        self.assertIn("gpu_line", new)
        self.assertEqual(module.comparable_config(old), module.comparable_config(new))
        live = Path("/home/treefor/.local/state/simqueue/node.json")
        if live.exists():
            self.assertEqual(module.comparable_config(json.loads(live.read_text())), module.comparable_config(new))
        changed = json.loads(json.dumps(new))
        changed["partitions"][0]["cpus"] = "0-2"
        self.assertNotEqual(module.comparable_config(old), module.comparable_config(changed))
        changed = json.loads(json.dumps(new))
        changed["lock_files"]["gpu_path"] = "/elsewhere/gpu-{uuid}.lock"
        self.assertNotEqual(module.comparable_config(old), module.comparable_config(changed))

    def test_wrapper_runtime_is_read_from_the_exec_line(self) -> None:
        module = load_upgrade_module()
        live = Path("/home/treefor/.local/bin/simq")
        if live.exists():
            self.assertIn(module.wrapper_runtime(live.read_text()), (OLD_RUNTIME, RUNTIME))
        self.assertIsNone(module.wrapper_runtime("#!/bin/sh\nexit 0\n"))

    def test_line_state_path_follows_the_configuration(self) -> None:
        module = load_upgrade_module()
        new = json.loads((RUNTIME / "config" / "zenobia-node.json").read_text())
        self.assertEqual(module.line_state_path(new), new["gpu_line"]["state_path"])
        old = json.loads((OLD_RUNTIME / "config" / "zenobia-node.json").read_text())
        self.assertEqual(module.line_state_path(old),
                         "/home/treefor/projects/.workspace/outputs/pic-profile/simq-gpu-line.json")


class UpgradeTests(NodeFixture, unittest.TestCase):
    def setUp(self) -> None:
        super().setUp()
        if not (OLD_RUNTIME / "source" / "simq").exists():
            self.skipTest("0.6.1 runtime not present")
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
        for part in old_doc["partitions"]:
            part["lend"] = True
        (self.sandbox_state / "node.json").write_text(json.dumps(old_doc, indent=2))
        new_doc = json.loads(json.dumps(old_doc))
        new_doc["gpu_line"] = {"poll_s": 0.05, "claim_cpus": list(self.cpus), "claim_hold_max_s": 60}
        self.new_config = self.base / "node-062.json"
        self.new_config.write_text(json.dumps(new_doc, indent=2))
        self.l1 = self.locks / "gpu.lock"
        self.l2 = self.locks / f"gpu-{FAKE_GPU}.lock"
        for path in (self.l1, self.l2):
            path.touch()
        self.state_file = self.locks / "simq-gpu-line.json"  # the default next to the private lease files
        module = load_upgrade_module()
        args = type("Args", (), {"execute": False, "home": str(self.home), "runtime": str(OLD_RUNTIME),
                                 "old_runtime": str(OLD_RUNTIME), "project_root": str(self.projects),
                                 "service": "simqueue.service", "systemctl": str(stub)})()
        sandbox = module.Upgrade(args)
        # the 0.6.1 installation as its upgrade.py left it
        old_wrapper = sandbox.wrapper_text(OLD_RUNTIME, "0.6.1").replace(
            "# borrow jobs run in simqborrow.slice at idle CPU priority; queued GPU jobs\n"
            "# win the GPU lease files through the GPU line (gpu-lease for direct use).\n",
            "# borrow jobs run in simqborrow.slice at idle CPU priority.\n")
        (self.bin / "simq").write_text(old_wrapper)
        (self.bin / "simq").chmod(0o755)
        (self.bin / "simqueue").symlink_to(OLD_RUNTIME / "simq")
        (self.units / "simqueue.service").write_text(sandbox.unit_text(OLD_RUNTIME, "0.6.1"))
        (self.units / "simqborrow.slice").write_text(SLICE_061)
        self.originals = {path: digest(path) for path in (self.bin / "simq", self.units / "simqueue.service",
                                                           self.units / "simqborrow.slice",
                                                           self.sandbox_state / "node.json")}
        self.sandbox_env = {k: v for k, v in os.environ.items() if not k.startswith("SIMQUEUE_")}
        self.sandbox_env.update({"STUB_DIR": str(self.stub_dir), "STUB_HOME": str(self.home)})
        self.pin = ["taskset", "-c", ",".join(str(c) for c in self.cpus)]

    def tearDown(self) -> None:
        self.stub("stop")
        for proc in self.holders:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except OSError:
                pass
        super().tearDown()
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

    def test_dry_run_changes_nothing(self) -> None:
        self.submit("-p", "pa", "-c", "1", *self.probe_command(1))  # creates the schema-8 database with 0.6.1
        before = digest(self.sandbox_state / "queue.sqlite3")
        self.stub("start")
        result = self.upgrade("install", "--config", str(self.new_config))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("[dry-run]", result.stdout)
        self.assertIn("gpu_line enabled=True", result.stdout)
        self.stub("stop")
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value)
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 8)
        self.assertNotEqual(before, "")
        self.assertFalse((self.home / ".local" / "state" / "simqueue-backups").exists())

    def test_install_adopts_running_061_jobs_and_rollback_adopts_062_jobs(self) -> None:
        self.stub("start")
        old_job = self.submit("-p", "pa", "-c", "1", *self.probe_command(15))
        self.wait_until(lambda: self.sandbox_state_of(old_job) == "running", message="0.6.1 job running")
        rows_before = self.db("SELECT * FROM jobs ORDER BY id")
        result = self.upgrade("install", "--execute", "--config", str(self.new_config))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        backup = Path(next(line.split("--backup ")[1].split()[0] for line in result.stdout.splitlines()
                           if "Rollback:" in line))
        self.assertTrue(backup.name.startswith("0.6.1-to-0.6.2-"))
        receipt = json.loads((backup / "install-receipt.json").read_text())
        self.assertTrue(receipt["verification"]["ok"])
        self.assertTrue(receipt["verification"]["jobs"][str(old_job)]["adopted"])
        self.assertTrue(receipt["verification"]["gpu_line_state"]["published"])
        self.assertEqual(receipt["verification"]["gpu_line_state"]["path"], str(self.state_file))
        self.assertTrue(receipt["database"]["rows_unchanged"])
        self.assertEqual(receipt["verification"]["partitions_rc"], 0)
        self.assertIn(f"simqueue 0.6.2 from {RUNTIME}", self.worker_log())
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 8)
        self.assertIn(str(RUNTIME), (self.bin / "simq").read_text())
        self.assertEqual(os.readlink(self.bin / "simqueue"), str(RUNTIME / "simq"))
        self.assertIn("SimQueue 0.6.2", (self.units / "simqueue.service").read_text())
        self.assertEqual(digest(self.units / "simqborrow.slice"), self.originals[self.units / "simqborrow.slice"])
        self.assertIn("gpu_line", json.loads((self.sandbox_state / "node.json").read_text()))
        self.assertEqual(json.loads(self.state_file.read_text())["version"], "0.6.2")
        manifest = json.loads((backup / "manifest.json").read_text())
        self.assertEqual(manifest["from"], "0.6.1")
        self.assertEqual(sorted(manifest["files"]), sorted(
            ["simq", "node.json", "worker.log", "units/simqueue.service", "units/simqborrow.slice"]))
        # every row as before (0.6.2 has no schema change)
        rows_after = self.db("SELECT * FROM jobs WHERE id <= ? ORDER BY id", old_job)
        self.assertEqual([r[0] for r in rows_before], [r[0] for r in rows_after])
        # the 0.6.1 job keeps running under its 0.6.1 jobwrap and is recorded by 0.6.2
        self.wait_until(lambda: self.sandbox_state_of(old_job) == "succeeded", timeout=60, message="old job done")
        # 0.6.2 work: an owner job that runs, and a GPU job waiting in the GPU line behind a raw holder
        owner_job = self.submit("-p", "pa", "-c", "1", *self.probe_command(25))
        holder = self.raw_holder()
        gpu_job = self.submit("-p", "pb", "--gpus", "1", *self.probe_command(1))
        self.wait_until(lambda: self.sandbox_state_of(owner_job) == "running", message="0.6.2 owner job running")
        self.wait_until(lambda: any(s.get("phase") == "queued" for s in self.claim_statuses()), timeout=30,
                        message="GPU claim queued in the kernel")
        self.assertTrue(self.active_claim_units())
        rolled = self.upgrade("rollback", "--execute", "--backup", str(backup))
        self.assertEqual(rolled.returncode, 0, rolled.stdout + rolled.stderr)
        rollback_receipt = json.loads((backup / "rollback-receipt.json").read_text())
        self.assertTrue(rollback_receipt["verification"]["jobs"][str(owner_job)]["adopted"])
        self.assertEqual(len(rollback_receipt["stopped_claim_units"]), 1)
        self.wait_until(lambda: not self.active_claim_units(), timeout=30, message="claim unit stopped")
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 8)
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value, f"{path} not restored byte for byte")
        self.assertEqual(os.readlink(self.bin / "simqueue"), str(OLD_RUNTIME / "simq"))
        self.assertTrue((backup / "queue-0.6.2-at-rollback.sqlite3").exists())
        self.assertEqual(json.loads(self.state_file.read_text())["heartbeat_unix"], 0)
        log = self.worker_log()
        self.assertIn(f"simqueue 0.6.1 from {OLD_RUNTIME}", log)
        self.assertIn(f"adopted job {owner_job} ", log.rsplit("simqueue 0.6.1 from", 1)[1])
        # 0.6.1 takes the GPU itself once the raw holder is gone
        self.assertEqual(self.sandbox_state_of(gpu_job), "queued")
        os.killpg(holder.pid, signal.SIGKILL)
        holder.wait(timeout=10)
        self.wait_until(lambda: self.sandbox_state_of(gpu_job) == "succeeded", timeout=60,
                        message="GPU job run by 0.6.1")
        self.wait_until(lambda: self.sandbox_state_of(owner_job) == "succeeded", timeout=60,
                        message="0.6.2 job recorded by 0.6.1")

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
            env=self.sandbox_env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120, check=False)
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn(f"job(s) [{job}] still starting", result.stdout)
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value)
        self.assertFalse((self.home / ".local" / "state" / "simqueue-backups").exists())

    def test_failed_verification_rolls_back_automatically(self) -> None:
        self.stub("start")
        job = self.submit("-p", "pa", "-c", "1", *self.probe_command(20))
        self.wait_until(lambda: self.sandbox_state_of(job) == "running", message="0.6.1 job running")
        result = self.upgrade("install", "--execute", "--config", str(self.new_config), timeout="0")
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn("rolling back automatically", result.stdout)
        self.assertIn("rolled back to SimQueue 0.6.1", result.stdout)
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value)
        self.assertEqual(self.db("PRAGMA user_version")[0][0], 8)
        self.wait_until(lambda: self.sandbox_state_of(job) == "succeeded", timeout=60, message="job recorded by 0.6.1")
        log = self.worker_log()
        self.assertIn(f"adopted job {job} ", log.rsplit(f"simqueue 0.6.1 from {OLD_RUNTIME}", 1)[1])


if __name__ == "__main__":
    unittest.main()
