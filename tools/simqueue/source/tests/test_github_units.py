"""SimQueue 0.7.0 GitHub reporting with real job units and a local fake API.

Every transition (queued, running, success, failure, OOM kill, time limit,
cancellation while running and while queued), the PR comment, retries while
GitHub is down or the token is missing, and that the token never appears in
any output, file, unit or job environment.
"""

from __future__ import annotations

import json
import os
import subprocess
import time
import unittest
from pathlib import Path

from fake_github import FakeGitHub
from node_fixture import PYTHON_FROZEN, NodeFixture

SHA = "89abcdef0123456789abcdef0123456789abcdef"
TOKEN = "ghp_unitTestToken_DO_NOT_LEAK_42"
ALLOC = "import sys\nb = bytearray(int(sys.argv[1]) << 20)\nfor i in range(0, len(b), 4096): b[i] = 1\n"


class GitHubUnitFixture(NodeFixture):
    min_cpus = 1

    def setUp(self) -> None:  # noqa: N802
        self.fake = FakeGitHub(TOKEN).start()
        super().setUp()
        self.token_file = self.base / "github-token"
        self.token_file.write_text(TOKEN + "\n")
        self.token_file.chmod(0o600)
        (self.projects / "alloc.py").write_text(ALLOC)

    def tearDown(self) -> None:  # noqa: N802
        super().tearDown()
        self.fake.stop()

    def config_document(self):
        document = super().config_document()
        cpus = ",".join(str(c) for c in (list(self.cpus) or [0, 1]))
        document["partitions"] = [{"name": "pa", "cpus": cpus, "default_mem": "256M", "max_mem": "2G",
                                   "default_time": "10:00", "max_time": "1:00:00", "default_cpus": 1}]
        document["gpus"] = []
        document["lock_files"].pop("gpu_path", None)
        document["lock_files"].pop("gpu_extra", None)
        document["launch"]["idle_seconds"] = 5
        document["memory"] = {"sample_seconds": 1, "record_seconds": 2}
        document["github"] = {"api_base": self.fake.api_base, "token_file": str(self.base / "github-token"),
                              "backoff_seconds": 0.5, "backoff_max_seconds": 2, "timeout_seconds": 5}
        return document

    def status(self, context: str, *extra: str) -> list[str]:
        return ["--github-status", f"o/r@{SHA}:{context}", *extra]

    def states(self, context: str) -> list[str]:
        return [s["state"] for s in self.fake.statuses_for(context)]

    def wait_status(self, context: str, state: str, timeout: float = 60) -> dict:
        def found():
            items = [s for s in self.fake.statuses_for(context) if s["state"] == state]
            return items[-1] if items else None
        return self.wait_until(found, timeout=timeout, message=f"{context} {state}")

    def wait_done(self, job_id: int, timeout: float = 60) -> dict:
        return self.wait_until(lambda: self.job(job_id) if self.job(job_id).get("github_done") else None,
                               timeout=timeout, message=f"job {job_id} GitHub posts done")


class TransitionTests(GitHubUnitFixture, unittest.TestCase):
    def test_every_transition(self) -> None:
        self.start_worker()
        gate = self.projects / "release"
        hold = self.submit("-p", "pa", "--cpus", ",".join(str(c) for c in self.cpus), "-n", "hold",
                           *self.status("t/hold"), "--", "/bin/sh", "-c", f"while [ ! -e {gate} ]; do sleep 0.1; done")
        self.wait_status("t/hold", "pending")
        queued = self.submit("-p", "pa", "-n", "waits", *self.status("t/waits"), "--", "/bin/true")
        first = self.wait_status("t/waits", "pending")
        self.assertEqual(first["description"], f"SimQueue job {queued} queued (1 cores, 256M)")
        gate.write_text("go")
        success = self.wait_status("t/waits", "success")
        self.assertEqual(self.states("t/waits"), ["pending", "pending", "success"])
        self.assertTrue(self.fake.statuses_for("t/waits")[1]["description"].startswith(
            f"SimQueue job {queued} running on CPUs "))
        self.assertRegex(success["description"], rf"^SimQueue job {queued} succeeded in \d+s, peak \d+M; exit 0$")
        self.assertEqual(self.states("t/hold")[-1], "success")
        failed = self.submit("-p", "pa", "-n", "exit3", *self.status("t/exit3"), "--", "/bin/sh", "-c", "exit 3")
        self.assertTrue(self.wait_status("t/exit3", "failure")["description"].endswith("; exit 3"))
        from node_fixture import real_oom_allowed

        oom = None
        if real_oom_allowed():  # never on zenobia; test_github.py covers the oom-kill status text
            oom = self.submit("-p", "pa", "-n", "oom", "--mem", "64M", *self.status("t/oom"),
                              "--", PYTHON_FROZEN, str(self.projects / "alloc.py"), "300")
            self.assertIn("out of memory at its 64M limit", self.wait_status("t/oom", "error")["description"])
        timeout = self.submit("-p", "pa", "-n", "slow", "-t", "0:02", *self.status("t/time"), "--", "/bin/sleep", "60")
        self.assertIn("time limit exceeded", self.wait_status("t/time", "error")["description"])
        running = self.submit("-p", "pa", "-n", "cancel-me", *self.status("t/cancel"), "--", "/bin/sleep", "60")
        self.wait_until(lambda: self.state_of(running) == "running", message="running")
        self.cli("cancel", str(running))
        self.assertIn("canceled", self.wait_status("t/cancel", "error")["description"])
        self.cli("pause")
        waiting = self.submit("-p", "pa", "-n", "never-runs", *self.status("t/queued-cancel"), "--", "/bin/true")
        self.wait_status("t/queued-cancel", "pending")
        self.cli("cancel", str(waiting))
        self.assertIn("canceled", self.wait_status("t/queued-cancel", "error")["description"])
        self.cli("resume")
        for job_id in (j for j in (hold, queued, failed, oom, timeout, running, waiting) if j is not None):
            done = self.wait_done(job_id)
            self.assertIsNone(done["github_error"], done)
        self.assertTrue(all(r["auth_ok"] for r in self.fake.requests))
        self.assertEqual(len(self.fake.comments), 0, "comments only with --github-comment")
        for status in self.fake.statuses:
            self.assertLessEqual(len(status["description"]), 140)

    def test_pr_comment_with_summary_once_and_sanitized(self) -> None:
        self.start_worker()
        script = ('echo "log line from @trevor about #12"; '
                  'printf "## Pilot result\\n\\nok, cc @team, fixes #5\\n" > "$SIMQ_SUMMARY"')
        job_id = self.submit("-p", "pa", "-n", "pilot", *self.status("t/pilot"), "--github-comment", "o/r#7",
                             "--github-summary", "pilot-summary.md", "--", "/bin/sh", "-c", script)
        self.wait_done(job_id)
        self.assertEqual(self.states("t/pilot")[-1], "success")
        self.assertEqual(len(self.fake.comments), 1)
        comment = next(iter(self.fake.comments.values()))
        body = comment["body"]
        self.assertTrue(body.startswith(f"<!-- simq-job:{job_id} -->"))
        self.assertIn("## Pilot result", body)
        self.assertIn("log line from", body)
        for forbidden in ("@trevor", "@team", "#12", "#5"):
            self.assertNotIn(forbidden, body)
        self.assertEqual(self.job(job_id)["github_comment_id"], comment["id"])
        self.assertTrue((self.projects / "pilot-summary.md").exists())


class RetryTests(GitHubUnitFixture, unittest.TestCase):
    def test_github_down_then_back_retries_with_backoff(self) -> None:
        self.fake.down = True
        self.start_worker()
        job_id = self.submit("-p", "pa", "-n", "down", *self.status("t/down"), "--", "/bin/true")
        self.wait_until(lambda: self.state_of(job_id) == "succeeded", message="the job runs regardless")
        self.wait_until(lambda: int(self.job(job_id).get("github_attempts") or 0) >= 2, timeout=30,
                        message="failed attempts recorded")
        job = self.job(job_id)
        self.assertIn("HTTP 503", job["github_error"])
        self.assertTrue(any(e["kind"] == "github-failed" for e in self.store.job_events(job_id)))
        self.assertIn("retrying with backoff", self.worker_log())
        self.fake.down = False
        done = self.wait_done(job_id, timeout=30)
        self.assertIsNone(done["github_error"])
        self.assertEqual(self.states("t/down")[-1], "success")
        posts = [r for r in self.fake.requests if r["path"].startswith("/repos/o/r/statuses/")]
        self.assertGreaterEqual(len(posts), 3)

    def test_missing_token_is_recorded_then_used_once_present(self) -> None:
        self.token_file.unlink()
        self.start_worker()
        job_id = self.submit("-p", "pa", "-n", "no-token", *self.status("t/token"), "--", "/bin/true")
        self.wait_until(lambda: "cannot read the token file" in str(self.job(job_id).get("github_error") or ""),
                        timeout=30, message="missing token recorded")
        self.token_file.write_text(TOKEN + "\n")
        self.wait_done(job_id, timeout=30)
        self.assertEqual(self.states("t/token")[-1], "success")


class TokenHygieneTests(GitHubUnitFixture, unittest.TestCase):
    def test_token_never_appears_in_any_output(self) -> None:
        self.start_worker()
        job_id = self.submit("-p", "pa", "-n", "env-dump", *self.status("t/env"), "--github-comment", "o/r#9",
                             "--github-summary", "s.md", "--", "/bin/sh", "-c", "env; env > \"$SIMQ_SUMMARY\"; sleep 3")
        self.wait_until(lambda: self.state_of(job_id) == "running", message="running")
        unit = str(self.job(job_id)["unit_name"])
        shown = subprocess.run(["systemctl", "--user", "show", unit], text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.DEVNULL, check=False).stdout
        self.assertIn("SIMQ_SUMMARY=", shown)
        self.assertNotIn(TOKEN, shown)
        self.wait_done(job_id)
        outputs = [self.cli(*args).stdout + self.cli(*args).stderr for args in
                   (("show", str(job_id)), ("list", "--all"), ("partitions", "--json", "--sample", "0"), ("doctor",))]
        for text in outputs:
            self.assertNotIn(TOKEN, text)
        for path in self.state.rglob("*"):
            if path.is_file():
                self.assertNotIn(TOKEN.encode(), path.read_bytes(), str(path))
        for path in self.projects.rglob("*"):
            if path.is_file():
                self.assertNotIn(TOKEN.encode(), path.read_bytes(), str(path))
        for path in self.base.glob("worker-*.out"):
            self.assertNotIn(TOKEN, path.read_text(errors="replace"))
        for request in self.fake.requests:
            self.assertNotIn(TOKEN, json.dumps(request["body"]))
        self.assertTrue(all(r["auth_ok"] for r in self.fake.requests))
        self.assertIn("SIMQ_SUMMARY=", self.log_text(job_id))


if __name__ == "__main__":
    unittest.main()
