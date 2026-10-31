"""switchover.py install / rollback in a sandbox home with a stub systemctl.

The real ~/.local/bin/simq, unit files, state and simqueue.service are never
touched: every path is under a temporary --home and systemctl is a stub that
records calls.  Canary jobs run as real units on the leased test CPUs.
"""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import subprocess
import unittest
from contextlib import closing
from pathlib import Path

from node_fixture import PYTHON_FROZEN, NodeFixture

from simqueue.config import format_cpu_list


RUNTIME = Path(__file__).resolve().parents[2]
SWITCHOVER = RUNTIME / "switchover.py"
LEGACY_SOURCE = "/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/source"
STUB = r'''#!/bin/sh
# stub systemctl: records calls, keeps one ActiveState for simqueue.service
log="$STUB_DIR/calls.log"; state="$STUB_DIR/active"
echo "$*" >> "$log"
[ "$1" = "--user" ] && shift
case "$1" in
  show)
    case "$*" in
      *ActiveState*) cat "$state" ;;
      *ExecStart*) echo "{ path=$STUB_HOME/.local/bin/simq ; argv[]=$STUB_HOME/.local/bin/simq worker --wait-for-lock }" ;;
      *) echo "" ;;
    esac ;;
  stop) echo inactive > "$state" ;;
  start) echo active > "$state" ;;
  list-units) : ;;
  *) : ;;
esac
exit 0
'''
LEGACY_IDLE = r'''
import sys
sys.path.insert(0, %(source)r)
from pathlib import Path
import os
from simqueue.store import QueueStore
store = QueueStore()
projects = Path(os.environ["SIMQUEUE_PROJECT_ROOT"])
for name in ("old-a", "old-b"):
    job = store.submit(name=name, cwd=projects, argv=["/bin/true"], environment={"A": "1"},
                       env_override_keys=[], requested_cores=1)
    store.claim_next((0,))
    store.mark_running(job["id"], pid=1, pgid=1, boot_id="b", process_start_ticks=1)
    store.finish(job["id"], state="succeeded", exit_code=0)
%(extra)s
'''


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SwitchoverTests(NodeFixture, unittest.TestCase):
    def setUp(self) -> None:
        super().setUp()
        if not (RUNTIME / "sealed-files.json").exists():
            self.skipTest("runtime not sealed yet")
        self.home = self.base / "home"
        self.sandbox_state = self.home / ".local" / "state" / "simqueue"
        bin_dir = self.home / ".local" / "bin"
        units = self.home / ".config" / "systemd" / "user"
        for directory in (self.sandbox_state, bin_dir, units):
            directory.mkdir(parents=True)
        (bin_dir / "simq").write_text("#!/bin/sh\n# 0.5 wrapper stand-in\nexit 0\n")
        (bin_dir / "simq").chmod(0o755)
        (bin_dir / "simqueue").symlink_to("/nonexistent/simqueue-0.5/simq")
        for name in ("simqueue.service", "simqueue-spheral-dispatch.path", "simqueue-spheral-dispatch.service"):
            (units / name).write_text(f"[Unit]\nDescription=0.5 {name}\n")
        (self.sandbox_state / "worker.log").write_text("0.5 worker log\n")
        self.stub_dir = self.base / "stub"
        self.stub_dir.mkdir()
        (self.stub_dir / "active").write_text("active\n")
        stub = self.stub_dir / "systemctl"
        stub.write_text(STUB)
        stub.chmod(0o755)
        self.originals = {
            path: digest(path)
            for path in [bin_dir / "simq", *(units / n for n in ("simqueue.service",
                         "simqueue-spheral-dispatch.path", "simqueue-spheral-dispatch.service"))]
        }
        self.sandbox_config = self.base / "sandbox-node.json"
        self.sandbox_config.write_text(json.dumps(self.config_document()))

    def legacy_db(self, extra: str = "") -> None:
        env = {"PATH": os.environ["PATH"], "SIMQUEUE_STATE_DIR": str(self.sandbox_state),
               "SIMQUEUE_PROJECT_ROOT": str(self.projects)}
        result = subprocess.run([PYTHON_FROZEN, "-c", LEGACY_IDLE % {"source": LEGACY_SOURCE, "extra": extra}],
                                env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)

    def switchover(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        env = dict(os.environ, STUB_DIR=str(self.stub_dir), STUB_HOME=str(self.home))
        for key in list(env):
            if key.startswith("SIMQUEUE_"):
                del env[key]
        return subprocess.run(
            [PYTHON_FROZEN, str(SWITCHOVER), *arguments, "--home", str(self.home), "--project-root",
             str(self.projects), "--config", str(self.sandbox_config), "--systemctl",
             str(self.stub_dir / "systemctl"), "--no-autostart"],
            env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=600, check=False,
        )

    def user_version(self) -> int:
        with closing(sqlite3.connect(self.sandbox_state / "queue.sqlite3")) as connection:
            return int(connection.execute("PRAGMA user_version").fetchone()[0])

    def test_dry_run_changes_nothing(self) -> None:
        self.legacy_db()
        before = digest(self.sandbox_state / "queue.sqlite3")
        result = self.switchover("install")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("[dry-run]", result.stdout)
        self.assertEqual(digest(self.sandbox_state / "queue.sqlite3"), before)
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value)
        self.assertFalse((self.home / ".local" / "state" / "simqueue-backups").exists())
        calls = (self.stub_dir / "calls.log").read_text()
        self.assertNotIn("stop", calls)

    def test_install_refuses_while_a_05_job_is_active(self) -> None:
        self.legacy_db(
            "j = store.submit(name='busy', cwd=projects, argv=['/bin/true'], environment={}, "
            "env_override_keys=[], requested_cores=1)\nstore.claim_next((0,))\n"
        )
        result = self.switchover("install", "--execute", "--no-canary")
        self.assertEqual(result.returncode, 1)
        self.assertIn("active", result.stdout)
        self.assertEqual(self.user_version(), 6)

    def test_install_canary_and_rollback_restore_everything(self) -> None:
        self.legacy_db()
        result = self.switchover("install", "--execute", "--canary", "pa=run", "--canary", "pb=test-only",
                                 "--foreground-worker", "--canary-timeout", "120")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.user_version(), 7)
        wrapper = (self.home / ".local" / "bin" / "simq").read_text()
        self.assertIn(str(RUNTIME / "source" / "simq"), wrapper)
        self.assertIn("SIMQUEUE_CONFIG=", wrapper)
        unit = (self.home / ".config" / "systemd" / "user" / "simqueue.service").read_text()
        self.assertIn("worker --wait-for-lock", unit)
        self.assertNotIn("CPUAffinity=", unit.replace("# No CPUAffinity", ""))
        self.assertTrue((self.sandbox_state / "node.json").exists())
        self.assertEqual(os.readlink(self.home / ".local" / "bin" / "simqueue"), str(RUNTIME / "simq"))
        backups = list((self.home / ".local" / "state" / "simqueue-backups").iterdir())
        self.assertEqual(len(backups), 1)
        backup = backups[0]
        receipt = json.loads((backup / "install-receipt.json").read_text())
        canaries = {c["partition"]: c for c in receipt["canaries"]}
        self.assertEqual(canaries["pa"]["state"], "succeeded", canaries)
        output = " ".join(canaries["pa"]["output"])
        self.assertIn("partition=pa", output)
        self.assertIn(f"affinity={format_cpu_list(self.part_a)}", output)
        self.assertIn(f"cpu_list={format_cpu_list(self.part_a)}", output)
        self.assertEqual(canaries["pb"]["returncode"], 0)
        self.assertIn('"would_start_now": true', canaries["pb"]["plan"])
        calls = (self.stub_dir / "calls.log").read_text().splitlines()
        self.assertIn("--user stop simqueue.service", calls)
        self.assertIn("--user start simqueue.service", calls)

        result = self.switchover("rollback", "--backup", str(backup), "--execute")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.user_version(), 6)
        for path, value in self.originals.items():
            self.assertEqual(digest(path), value, f"{path} not restored")
        self.assertEqual(os.readlink(self.home / ".local" / "bin" / "simqueue"), "/nonexistent/simqueue-0.5/simq")
        self.assertFalse((self.sandbox_state / "node.json").exists())
        self.assertTrue((backup / "queue-0.6-at-rollback.sqlite3").exists())
        reopened = subprocess.run(
            [PYTHON_FROZEN, "-c",
             f"import sys; sys.path.insert(0, {LEGACY_SOURCE!r}); from simqueue.store import QueueStore;"
             "s = QueueStore(); print(len(s.list_jobs(include_all=True)))"],
            env={"PATH": os.environ["PATH"], "SIMQUEUE_STATE_DIR": str(self.sandbox_state),
                 "SIMQUEUE_PROJECT_ROOT": str(self.projects)},
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        )
        self.assertEqual(reopened.returncode, 0, reopened.stderr)
        self.assertEqual(reopened.stdout.strip(), "2")


if __name__ == "__main__":
    unittest.main()
