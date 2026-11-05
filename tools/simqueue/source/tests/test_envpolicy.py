"""SimQueue 0.7.0: what reaches a job from the submitter's environment.

Tokens in the submitter's environment (GH_TOKEN, RUNNER_TRACKING_ID, ...)
must never reach the database, the job, its log, its spec file or `simq
show`; secrets given explicitly with -e reach the job but are stored
redacted; the schema-9 migration scrubs environments stored by 0.6.3.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import sqlite3
import subprocess
import tempfile
import unittest
from contextlib import closing
from pathlib import Path

from node_fixture import PYTHON_FROZEN, NodeFixture

from simqueue.envpolicy import (
    REDACTED,
    capture,
    is_secret,
    parse_export,
    read_secrets,
    scrub_environment,
    scrub_json,
    scrub_spec,
    secrets_path,
    write_secrets,
)

GH = "ghp_inheritedTokenValue_123456"
RUNNER = "runner-tracking-1f2e3d"
EXPLICIT = "explicit-api-key-987654"
RUNTIME_063 = Path("/home/treefor/projects/runtimes/simqueue-0.6.3-20261105/source")
SUBMITTER = {"HOME": "/home/x", "USER": "x", "PATH": "/usr/bin", "LANG": "C.UTF-8", "TZ": "UTC",
             "GH_TOKEN": GH, "RUNNER_TRACKING_ID": RUNNER, "GITHUB_SHA": "abc", "AWS_SECRET_ACCESS_KEY": "s3",
             "OMP_NUM_THREADS": "4", "SSH_AUTH_SOCK": "/run/agent", "MY_PASSWORD": "pw"}


class PolicyTests(unittest.TestCase):
    def test_patterns(self) -> None:
        for name in ("GH_TOKEN", "github_token", "MY_SECRET", "DB_PASSWORD", "API_KEY", "GITHUB_SHA",
                     "ACTIONS_RUNTIME_URL", "RUNNER_TRACKING_ID", "SSH_KEY_PATH"):
            self.assertTrue(is_secret(name), name)
        for name in ("HOME", "PATH", "OMP_NUM_THREADS", "SSH_AUTH_SOCK", "CUDA_VISIBLE_DEVICES"):
            self.assertFalse(is_secret(name), name)

    def test_allowlist_export_and_explicit_secrets(self) -> None:
        stored, secrets, keys, withheld = capture(SUBMITTER)
        self.assertEqual(stored, {"HOME": "/home/x", "USER": "x", "PATH": "/usr/bin", "LANG": "C.UTF-8", "TZ": "UTC"})
        self.assertEqual((secrets, keys, withheld), ({}, [], []))
        stored, _, _, withheld = capture(SUBMITTER, export=parse_export("OMP_NUM_THREADS,GH_TOKEN"))
        self.assertEqual(stored["OMP_NUM_THREADS"], "4")
        self.assertNotIn("GH_TOKEN", stored)
        self.assertEqual(withheld, ["GH_TOKEN"])
        stored, _, _, withheld = capture(SUBMITTER, export=parse_export("ALL"))
        self.assertEqual(set(stored), {"HOME", "USER", "PATH", "LANG", "TZ", "OMP_NUM_THREADS", "SSH_AUTH_SOCK"})
        self.assertEqual(withheld, sorted(["GH_TOKEN", "RUNNER_TRACKING_ID", "GITHUB_SHA", "AWS_SECRET_ACCESS_KEY",
                                           "MY_PASSWORD"]))
        stored, secrets, keys, _ = capture(SUBMITTER, [("MY_KEY", EXPLICIT), ("MODE", "fast")])
        self.assertEqual((stored["MY_KEY"], stored["MODE"], secrets, keys),
                         (REDACTED, "fast", {"MY_KEY": EXPLICIT}, ["MY_KEY", "MODE"]))
        self.assertEqual(parse_export(None), ())
        with self.assertRaises(ValueError):
            parse_export("A,B C")

    def test_scrub(self) -> None:
        cleaned, names = scrub_environment(SUBMITTER)
        self.assertEqual(sorted(names), sorted(["GH_TOKEN", "RUNNER_TRACKING_ID", "GITHUB_SHA",
                                                "AWS_SECRET_ACCESS_KEY", "MY_PASSWORD"]))
        self.assertNotIn(GH, json.dumps(cleaned))
        text = json.dumps({"PATH": "/bin"})
        self.assertEqual(scrub_json(text), (text, []))  # untouched rows keep their exact text
        base = Path(tempfile.mkdtemp(prefix="simq07-env-"))
        try:
            write_secrets(base, 5, "abcdef0123456789", {"MY_KEY": EXPLICIT})
            path = secrets_path(base, 5, "abcdef0123456789")
            self.assertEqual(path.stat().st_mode & 0o777, 0o600)
            self.assertEqual(read_secrets(base, 5, "abcdef0123456789"), {"MY_KEY": EXPLICIT})
            spec = base / "spec.json"
            spec.write_text(json.dumps({"env": {"GH_TOKEN": GH, "PATH": "/bin"}, "job_id": 5}))
            self.assertEqual(scrub_spec(spec), ["GH_TOKEN"])
            self.assertNotIn(GH, spec.read_text())
            self.assertEqual(json.loads(spec.read_text())["job_id"], 5)
        finally:
            shutil.rmtree(base)


class CliEnvironmentTests(NodeFixture, unittest.TestCase):
    needs_units = False

    def test_submission_stores_no_inherited_token(self) -> None:
        env = dict(self.environment, GH_TOKEN=GH, RUNNER_TRACKING_ID=RUNNER, OMP_NUM_THREADS="4")
        result = self.cli("submit", "--json", "-e", f"MY_KEY={EXPLICIT}", "-e", "MODE=fast", "--export",
                          "OMP_NUM_THREADS", "--", "/bin/true", env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        job = self.store.require(json.loads(result.stdout)["job_id"])
        stored = json.loads(job["env_json"])
        self.assertEqual((stored["MY_KEY"], stored["MODE"], stored["OMP_NUM_THREADS"]), (REDACTED, "fast", "4"))
        self.assertNotIn("GH_TOKEN", stored)
        self.assertNotIn("RUNNER_TRACKING_ID", stored)
        self.assertEqual(read_secrets(self.state, job["id"], job["launch_token"]), {"MY_KEY": EXPLICIT})
        database = (self.state / "queue.sqlite3").read_bytes() + (self.state / "queue.sqlite3-wal").read_bytes() \
            if (self.state / "queue.sqlite3-wal").exists() else (self.state / "queue.sqlite3").read_bytes()
        for value in (GH, RUNNER, EXPLICIT):
            self.assertNotIn(value.encode(), database)
        for args in (("show", str(job["id"])), ("list", "--all")):
            text = self.cli(*args, env=env).stdout
            for value in (GH, RUNNER, EXPLICIT):
                self.assertNotIn(value, text)
        noisy = self.cli("submit", "--export", "ALL", "--", "/bin/true", env=env)
        self.assertIn("not passed to the job", noisy.stderr)
        self.assertIn("GH_TOKEN", noisy.stderr)
        self.assertNotIn(GH, noisy.stderr + noisy.stdout)

    def test_env_i_submission_still_works(self) -> None:
        env = {k: v for k, v in self.environment.items() if k.startswith("SIMQUEUE_")}
        env["PATH"] = "/usr/bin:/bin"
        result = subprocess.run(["env", "-i", *[f"{k}={v}" for k, v in env.items()], PYTHON_FROZEN,
                                 str(self.launcher), "submit", "--json", "--", "/bin/true"],
                                cwd=self.projects, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                check=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        stored = json.loads(self.store.require(json.loads(result.stdout)["job_id"])["env_json"])
        self.assertEqual({k: v for k, v in stored.items() if not k.startswith("SIMQUEUE_")}, {"PATH": "/usr/bin:/bin"})


SETUP_063 = r"""
import sys
sys.path.insert(0, %(source)r)
from pathlib import Path
from simqueue.store import QueueStore
s = QueueStore()
env = {"PATH": "/bin", "GH_TOKEN": %(gh)r, "RUNNER_TRACKING_ID": %(runner)r, "MY_KEY": %(explicit)r}
done = s.submit(name="done", cwd=Path(%(cwd)r), argv=["/bin/true"], environment=env, env_override_keys=[],
                requested_cores=1)
waiting = s.submit(name="waiting", cwd=Path(%(cwd)r), argv=["/bin/true"], environment=env,
                   env_override_keys=["MY_KEY"], requested_cores=1)
plain = s.submit(name="plain", cwd=Path(%(cwd)r), argv=["/bin/true"], environment={"PATH": "/bin"},
                 env_override_keys=[], requested_cores=1)
s.cancel(done["id"])
print(done["id"], waiting["id"], plain["id"])
"""


class ScrubMigrationTests(unittest.TestCase):
    def setUp(self) -> None:
        if not RUNTIME_063.exists():
            self.skipTest("0.6.3 runtime not present")
        self.base = Path(tempfile.mkdtemp(prefix="simq07-scrub-"))
        self.saved = dict(os.environ)
        self.env = {k: v for k, v in os.environ.items() if not k.startswith("SIMQUEUE_")}
        self.env.update({"SIMQUEUE_STATE_DIR": str(self.base / "state"), "SIMQUEUE_PROJECT_ROOT": str(self.base)})
        os.environ.update(self.env)

    def tearDown(self) -> None:
        os.environ.clear()
        os.environ.update(self.saved)
        shutil.rmtree(self.base, ignore_errors=True)

    def test_migration_scrubs_stored_environments_and_keeps_explicit_secrets_of_waiting_jobs(self) -> None:
        output = subprocess.run([PYTHON_FROZEN, "-c", SETUP_063 % {"source": str(RUNTIME_063), "cwd": str(self.base),
                                                                   "gh": GH, "runner": RUNNER, "explicit": EXPLICIT}],
                                env=self.env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        done, waiting, plain = (int(x) for x in output.stdout.split())
        database = self.base / "state" / "queue.sqlite3"
        with closing(sqlite3.connect(database)) as connection:
            plain_before = connection.execute("SELECT env_json FROM jobs WHERE id = ?", (plain,)).fetchone()[0]
        from simqueue.store import QueueStore

        store = QueueStore()
        rows = {job["id"]: json.loads(job["env_json"]) for job in store.list_jobs(include_all=True)}
        for job_id in (done, waiting):
            self.assertEqual((rows[job_id]["GH_TOKEN"], rows[job_id]["RUNNER_TRACKING_ID"], rows[job_id]["MY_KEY"]),
                             (REDACTED, REDACTED, REDACTED))
        self.assertEqual(store.require(plain)["env_json"], plain_before)
        token = store.require(waiting)["launch_token"]
        self.assertEqual(read_secrets(self.base / "state", waiting, token), {"MY_KEY": EXPLICIT})
        self.assertEqual(read_secrets(self.base / "state", done, store.require(done)["launch_token"]), {})
        store.close()
        with closing(sqlite3.connect(database)) as connection:
            connection.execute("PRAGMA wal_checkpoint(TRUNCATE)")
        content = database.read_bytes()
        for value in (GH, RUNNER, EXPLICIT):
            self.assertNotIn(value.encode(), content)
        self.assertEqual(store.get_meta("env_scrubbed_through"), str(plain))


class EnvironmentUnitTests(NodeFixture, unittest.TestCase):
    min_cpus = 1

    def test_job_gets_explicit_secrets_and_never_inherited_tokens(self) -> None:
        self.start_worker()
        script = ("import hashlib, json, os, sys\n"
                  "print(json.dumps({'names': sorted(os.environ), 'key': hashlib.sha256("
                  "os.environ.get('MY_KEY', '').encode()).hexdigest()}))\n")
        (self.projects / "envprobe.py").write_text(script)
        env = dict(self.environment, GH_TOKEN=GH, RUNNER_TRACKING_ID=RUNNER)
        result = self.cli("submit", "--json", "-p", "pa", "-c", "1", "-e", f"MY_KEY={EXPLICIT}", "--",
                          PYTHON_FROZEN, str(self.projects / "envprobe.py"), env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        job_id = int(json.loads(result.stdout)["job_id"])
        self.wait_until(lambda: self.state_of(job_id) == "succeeded", message="env probe")
        job = self.job(job_id)
        log = self.log_text(job_id)
        report = json.loads(next(line for line in log.splitlines() if line.startswith("{")))
        self.assertEqual(report["key"], hashlib.sha256(EXPLICIT.encode()).hexdigest())
        self.assertNotIn("GH_TOKEN", report["names"])
        self.assertNotIn("RUNNER_TRACKING_ID", report["names"])
        self.assertIn("XDG_RUNTIME_DIR", report["names"])  # session plumbing from the user manager
        spec = Path(job["job_dir"]) / "spec.json"
        self.assertEqual(json.loads(spec.read_text())["env"]["MY_KEY"], REDACTED)
        self.assertFalse(secrets_path(self.state, job_id, job["launch_token"]).exists())
        for path in self.state.rglob("*"):
            if path.is_file():
                for value in (GH, RUNNER, EXPLICIT):
                    self.assertNotIn(value.encode(), path.read_bytes(), str(path))
        shown = self.cli("show", str(job_id)).stdout
        for value in (GH, RUNNER, EXPLICIT):
            self.assertNotIn(value, shown)


if __name__ == "__main__":
    unittest.main()
