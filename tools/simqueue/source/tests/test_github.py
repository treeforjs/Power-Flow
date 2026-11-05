"""SimQueue 0.7.0 GitHub reporting: pure parts and the client against a local
fake API (no unit).  Transitions of real jobs are in test_github_units.py."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import threading
import unittest
from pathlib import Path

from fake_github import FakeGitHub
from node_fixture import NodeFixture

from simqueue.github import (
    FINAL,
    ZWSP,
    Client,
    GitHubConfig,
    PostError,
    Poster,
    TokenSource,
    build_task,
    comment_body,
    comment_marker,
    final_status,
    job_phase,
    parse_comment_target,
    parse_github_config,
    parse_status_target,
    post_task,
    sanitize,
    status_payload,
    targets_from_args,
)

SHA = "0123456789abcdef0123456789abcdef01234567"
TOKEN = "ghp_faketoken_SECRET_0123456789"


def job(**fields):
    base = {"id": 42, "name": "ci-tests", "state": "succeeded", "result": "success", "exit_code": 0,
            "message": None, "started_at": "2026-11-05T01:00:00Z", "finished_at": "2026-11-05T01:12:03Z",
            "cpu_list": "6-7", "requested_cores": 2, "mem_bytes": 4 << 30, "mem_peak_bytes": 1288490188,
            "log_path": None, "preempt_count": 0, "created_at": "2026-11-05T00:59:00Z"}
    base.update(fields)
    return base


class TargetTests(unittest.TestCase):
    def test_parse(self) -> None:
        self.assertEqual(parse_status_target(f"treeforjs/Power-Flow@{SHA.upper()}:simq/tests"),
                         {"repo": "treeforjs/Power-Flow", "sha": SHA, "context": "simq/tests"})
        self.assertEqual(parse_comment_target("o/r#12"), {"repo": "o/r", "pr": 12})
        for bad in ("o/r@abc:ctx", f"o@{SHA}:c", f"o/r@{SHA}", "o/r:ctx"):
            with self.assertRaises(ValueError, msg=bad):
                parse_status_target(bad)
        for bad in ("o/r", "o/r#x", "#3"):
            with self.assertRaises(ValueError, msg=bad):
                parse_comment_target(bad)

    def test_targets_from_args(self) -> None:
        self.assertIsNone(targets_from_args(None, None, None, Path("/w")))
        document = targets_from_args([f"o/r@{SHA}:a", f"o/r@{SHA}:b"], "o/r#3", "out/summary.md", Path("/w"))
        self.assertEqual([s["context"] for s in document["statuses"]], ["a", "b"])
        self.assertEqual(document["summary"], "/w/out/summary.md")
        with self.assertRaises(ValueError):
            targets_from_args([f"o/r@{SHA}:a", f"o/r@{SHA}:a"], None, None, Path("/w"))

    def test_config(self) -> None:
        self.assertEqual(parse_github_config(None), GitHubConfig())
        config = parse_github_config({"api_base": "http://127.0.0.1:9/", "token_file": "/x/t", "max_attempts": 3})
        self.assertEqual((config.api_base, config.token_file, config.max_attempts), ("http://127.0.0.1:9", "/x/t", 3))
        for bad in ({"api_base": "ftp://x"}, {"token_file": "rel"}, {"nope": 1}, {"timeout_seconds": 0}):
            with self.assertRaises(ValueError, msg=str(bad)):
                parse_github_config(bad)


class PayloadTests(unittest.TestCase):
    def test_phases(self) -> None:
        self.assertEqual(job_phase(job(state="queued")), "queued:0")
        self.assertEqual(job_phase(job(state="queued", preempt_count=2)), "queued:2")
        self.assertEqual(job_phase(job(state="running")), "running:2026-11-05T01:00:00Z")
        self.assertEqual(job_phase(job(state="starting")), "running:2026-11-05T01:00:00Z")
        for state in ("succeeded", "failed", "canceled", "interrupted"):
            self.assertEqual(job_phase(job(state=state)), FINAL)

    def test_final_states(self) -> None:
        cases = [
            (job(), ("success", "exit 0")),
            (job(state="failed", exit_code=3, result="exit-code"), ("failure", "exit 3")),
            (job(state="failed", exit_code=None, result="timeout"), ("error", "time limit exceeded")),
            (job(state="failed", exit_code=-9, result="oom-kill"), ("error", "out of memory at its 4.0G limit")),
            (job(state="canceled", exit_code=-15, result="success"), ("error", "canceled")),
            (job(state="canceled", result="preempted"), ("error", "preempted and canceled")),
        ]
        for item, expected in cases:
            self.assertEqual(final_status(item), expected)
        state, reason = final_status(job(state="interrupted", message="machine restarted while the job ran"))
        self.assertEqual(state, "error")
        self.assertTrue(reason.startswith("lost: machine restarted"))

    def test_descriptions_are_short_and_carry_the_facts(self) -> None:
        state, text = status_payload(job(), FINAL)
        self.assertEqual(state, "success")
        self.assertEqual(text, "SimQueue job 42 succeeded in 12m03s, peak 1.2G; exit 0")
        state, text = status_payload(job(state="queued", started_at=None), "queued:0")
        self.assertEqual((state, text), ("pending", "SimQueue job 42 queued (2 cores, 4.0G)"))
        state, text = status_payload(job(state="running"), "running:x")
        self.assertEqual((state, text), ("pending", "SimQueue job 42 running on CPUs 6-7"))
        long = job(state="interrupted", message="x" * 400)
        self.assertLessEqual(len(status_payload(long, FINAL)[1]), 140)

    def test_sanitize_prevents_mentions_and_references(self) -> None:
        text = sanitize("ping @trevor and @org/team about #12 and owner/repo#345, not #x")
        self.assertNotIn("@trevor", text)
        self.assertNotIn("@org", text)
        self.assertNotIn("#12", text)
        self.assertNotIn("#345", text)
        self.assertIn("#x", text)
        self.assertEqual(text.replace(ZWSP, ""), "ping @trevor and @org/team about #12 and owner/repo#345, not #x")
        _, description = status_payload(job(state="failed", exit_code=1, result="exit-code",
                                             message="see #9"), FINAL)
        self.assertNotIn("#9", description)

    def test_comment_body(self) -> None:
        base = Path(tempfile.mkdtemp(prefix="simq07-gh-"))
        try:
            log = base / "job.log"
            log.write_text("".join(f"line {i}\n" for i in range(50)) + "boom @trevor fixes #77\n```\n")
            summary = base / "summary.md"
            summary.write_text("## Pilot\n- result ok, thanks @someone, closes #5\n")
            body = comment_body(job(log_path=str(log)), GitHubConfig(), "zenobia", str(summary))
            self.assertTrue(body.startswith(comment_marker(42)))
            self.assertIn("| State | succeeded (exit 0) |", body)
            self.assertIn("| Peak memory | 1.2G of 4.0G requested |", body)
            self.assertIn("| Cores | 2 (CPUs 6-7) |", body)
            self.assertIn("line 49", body)
            self.assertNotIn("line 29\n", body)  # only the last 20 lines
            self.assertIn("## Pilot", body)
            for forbidden in ("@trevor", "@someone", "#77", "#5"):
                self.assertNotIn(forbidden, body)
            self.assertEqual(body.count("```"), 2)  # the job's own fence cannot break out
        finally:
            shutil.rmtree(base)


class ClientTests(unittest.TestCase):
    def setUp(self) -> None:
        self.fake = FakeGitHub(TOKEN).start()
        self.base = Path(tempfile.mkdtemp(prefix="simq07-gh-"))
        self.token_file = self.base / "token"
        self.token_file.write_text(TOKEN + "\n")
        self.config = GitHubConfig(api_base=self.fake.api_base, token_file=str(self.token_file), timeout_seconds=5)
        self.client = Client(self.config, TokenSource(self.config))

    def tearDown(self) -> None:
        self.fake.stop()
        shutil.rmtree(self.base, ignore_errors=True)

    def task(self, item: dict, phase: str = FINAL, comment: bool = True) -> dict:
        document = {"statuses": [{"repo": "o/r", "sha": SHA, "context": "simq/test"}]}
        if comment:
            document["comment"] = {"repo": "o/r", "pr": 7}
        return build_task(dict(item, github_json=json.dumps(document), github_comment_id=None), phase, "zenobia")

    def test_statuses_and_one_comment_edited_in_place(self) -> None:
        result = post_task(self.client, self.task(job(state="running"), "running:x"))
        self.assertTrue(result["ok"])
        self.assertIsNone(result.get("comment_id"))  # the comment comes only at the end
        result = post_task(self.client, self.task(job()))
        comment_id = result["comment_id"]
        self.assertEqual([s["state"] for s in self.fake.statuses_for("simq/test")], ["pending", "success"])
        self.assertEqual(len(self.fake.comments), 1)
        # posted again (a retry after a lost answer, a restarted worker): found by its marker, edited
        again = post_task(self.client, self.task(job(state="failed", exit_code=2, result="exit-code")))
        self.assertEqual(again["comment_id"], comment_id)
        self.assertEqual(len(self.fake.comments), 1)
        self.assertEqual(self.fake.comments[comment_id]["edits"], 1)
        self.assertIn("failed (exit 2)", self.fake.comments[comment_id]["body"])
        methods = [r["method"] for r in self.fake.requests if "/comments" in r["path"]]
        self.assertEqual(methods, ["GET", "POST", "GET", "PATCH"])
        self.assertTrue(all(r["auth_ok"] for r in self.fake.requests))

    def test_known_comment_id_is_patched_and_recreated_when_deleted(self) -> None:
        task = self.task(job())
        first = post_task(self.client, task)
        task["comment_id"] = first["comment_id"]
        post_task(self.client, task)
        self.assertEqual(self.fake.comments[first["comment_id"]]["edits"], 1)
        del self.fake.comments[first["comment_id"]]
        third = post_task(self.client, task)
        self.assertNotEqual(third["comment_id"], first["comment_id"])
        self.assertEqual(len(self.fake.comments), 1)

    def test_errors_are_classified_and_never_carry_the_token(self) -> None:
        self.fake.down = True
        with self.assertRaises(PostError) as caught:
            post_task(self.client, self.task(job(), comment=False))
        self.assertTrue(caught.exception.retry)
        self.assertIn("HTTP 503", str(caught.exception))
        self.fake.down = False
        self.fake.fail_status = 422
        with self.assertRaises(PostError) as caught:
            post_task(self.client, self.task(job(), comment=False))
        self.assertFalse(caught.exception.retry)
        self.fake.fail_status = None
        self.token_file.write_text("wrong-token\n")
        tokens = TokenSource(GitHubConfig(api_base=self.fake.api_base, token_file=str(self.token_file)))
        with self.assertRaises(PostError) as caught:
            post_task(Client(self.config, tokens), self.task(job(), comment=False))
        self.assertTrue(caught.exception.auth and caught.exception.retry)
        unreachable = GitHubConfig(api_base="http://127.0.0.1:9", token_file=str(self.token_file), timeout_seconds=2)
        with self.assertRaises(PostError) as caught:
            post_task(Client(unreachable, TokenSource(unreachable)), self.task(job(), comment=False))
        self.assertTrue(caught.exception.retry)
        missing = GitHubConfig(api_base=self.fake.api_base, token_file=str(self.base / "absent"))
        with self.assertRaises(PostError) as caught:
            TokenSource(missing).get()
        self.assertTrue(caught.exception.retry)
        for text in (str(caught.exception),):
            self.assertNotIn(TOKEN, text)

    def test_gh_auth_token_is_used_without_a_token_file(self) -> None:
        gh = self.base / "gh"
        gh.write_text(f"#!/bin/sh\n[ \"$1 $2\" = 'auth token' ] && echo {TOKEN}\n")
        gh.chmod(0o700)
        config = GitHubConfig(api_base=self.fake.api_base, gh=str(gh))
        result = post_task(Client(config, TokenSource(config)), self.task(job(), comment=False))
        self.assertTrue(result["ok"])
        failing = self.base / "gh-out"
        failing.write_text("#!/bin/sh\nexit 1\n")
        failing.chmod(0o700)
        with self.assertRaises(PostError):
            TokenSource(GitHubConfig(gh=str(failing))).get()

    def test_poster_thread_hands_back_results(self) -> None:
        woke = threading.Event()
        poster = Poster(self.config, woke.set)
        poster.submit(self.task(job(), comment=False))
        self.assertTrue(woke.wait(10))
        results = poster.results()
        self.assertEqual([(r["job_id"], r["ok"]) for r in results], [(42, True)])
        self.assertFalse(poster.in_flight)
        self.fake.down = True
        woke.clear()
        poster.submit(self.task(job(), comment=False))
        self.assertTrue(woke.wait(10))
        failed = poster.results()[0]
        self.assertEqual((failed["ok"], failed["retry"]), (False, True))
        self.assertNotIn(TOKEN, json.dumps(failed))
        poster.stop()


class CliGitHubTests(NodeFixture, unittest.TestCase):
    needs_units = False

    def test_options_are_validated_recorded_and_fingerprinted(self) -> None:
        result = self.cli("submit", "--json", "-p", "pa", "--github-status", f"o/r@{SHA}:simq/a",
                          "--github-status", f"o/r@{SHA}:simq/b", "--github-comment", "o/r#4",
                          "--github-summary", "out/summary.md", "--idempotency-key", "gh-1", "--", "/bin/true")
        self.assertEqual(result.returncode, 0, result.stderr)
        payload = json.loads(result.stdout)
        stored = json.loads(self.store.require(payload["job_id"])["github_json"])
        self.assertEqual([s["context"] for s in stored["statuses"]], ["simq/a", "simq/b"])
        self.assertEqual(stored["comment"], {"repo": "o/r", "pr": 4})
        self.assertEqual(stored["summary"], str(self.projects / "out" / "summary.md"))
        show = self.cli("show", str(payload["job_id"])).stdout
        self.assertIn("GitHub status", show)
        self.assertIn("GitHub posted          queued due", show)
        different = self.cli("submit", "--json", "-p", "pa", "--github-status", f"o/r@{SHA}:simq/c",
                             "--idempotency-key", "gh-1", "--", "/bin/true")
        self.assertNotEqual(different.returncode, 0)
        self.assertIn("different submission request", different.stderr)
        bad = self.cli("submit", "-p", "pa", "--github-status", "o/r@123:x", "--", "/bin/true")
        self.assertEqual(bad.returncode, 2)
        self.assertIn("full commit SHA", bad.stderr)
        run_help = self.cli("run", "--help").stdout
        for option in ("--github-status", "--github-comment", "--github-summary"):
            self.assertIn(option, run_help)

    def test_disabled_in_the_node_config(self) -> None:
        document = self.config_document()
        document["github"] = {"enabled": False}
        self.write_config(document)
        result = self.cli("submit", "-p", "pa", "--github-status", f"o/r@{SHA}:x", "--", "/bin/true")
        self.assertEqual(result.returncode, 2)
        self.assertIn("disabled", result.stderr)


if __name__ == "__main__":
    unittest.main()
