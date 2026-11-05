"""GitHub reporting of job state (SimQueue 0.7.0).

A job submitted with ``--github-status OWNER/REPO@SHA:CONTEXT`` (repeatable)
gets a commit status per target on every transition; with ``--github-comment
OWNER/REPO#PR`` one pull-request comment at its end, edited in place when
posted again.  The worker posts, not the job: the submission records the
targets in the job row (``github_json``) and the worker compares each job's
*phase* (queued, running, final) with the phase last posted.

Statuses
    ``pending`` while queued and while running, ``success`` / ``failure``
    from the exit code, ``error`` for cancellation, time limit, out-of-memory
    kill, preemption-cancel and lost units.  Descriptions stay within 140
    characters (job id, runtime, peak memory, exit code or reason).
    Commit statuses never send notifications.
Comment
    Only with ``--github-comment``: job id, state, runtime, peak memory,
    cores, the last log lines and the job's summary markdown (``--github-
    summary PATH``, exported to the job as ``SIMQ_SUMMARY``).  One comment per
    job, found again by the hidden marker ``<!-- simq-job:N -->`` and edited
    (PATCH) instead of posted again; edits do not notify.  Text from the job
    (name, log lines, summary) is sanitized: ``@`` and ``#<digits>`` get a
    zero-width space, so nothing in a log can mention a person or reference
    an issue.

Posting never blocks scheduling: one daemon thread of the worker process
makes the HTTP requests (``urllib``, ``timeout_seconds``) and hands results
back through a queue; the worker records them, and retries failures with
exponential backoff from its own event loop.  The token comes from
``github.token_file`` or ``gh auth token`` at post time, is kept only in
that thread's memory (at most ``token_cache_seconds``) and is never written
to the database, logs, job environment or unit files.
"""

from __future__ import annotations

import json
import os
import queue
import re
import subprocess
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from . import __version__

STATUS_PATTERN = re.compile(r"^([A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+)@([0-9a-fA-F]{40}|[0-9a-fA-F]{64}):(.{1,255})$")
COMMENT_PATTERN = re.compile(r"^([A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+)#([0-9]{1,9})$")
DESCRIPTION_MAX = 140
COMMENT_MAX = 60000
ZWSP = "​"
FINAL = "final"


@dataclass(frozen=True)
class GitHubConfig:
    enabled: bool = True
    api_base: str = "https://api.github.com"
    token_file: str | None = None
    gh: str = "gh"
    timeout_seconds: float = 10.0
    token_cache_seconds: float = 600.0
    backoff_seconds: float = 15.0
    backoff_max_seconds: float = 900.0
    max_attempts: int = 12
    log_lines: int = 20
    summary_max_bytes: int = 40000


def parse_github_config(value: Any) -> GitHubConfig:
    if value is None:
        return GitHubConfig()
    if not isinstance(value, dict):
        raise ValueError("github must be a JSON object")
    allowed = {"enabled", "api_base", "token_file", "gh", "timeout_seconds", "token_cache_seconds",
               "backoff_seconds", "backoff_max_seconds", "max_attempts", "log_lines", "summary_max_bytes", "notes"}
    extra = sorted(set(value) - allowed)
    if extra:
        raise ValueError(f"github has unknown key(s): {', '.join(extra)}")
    defaults = GitHubConfig()
    enabled = value.get("enabled", defaults.enabled)
    if not isinstance(enabled, bool):
        raise ValueError("github.enabled must be true or false")
    api_base = str(value.get("api_base", defaults.api_base)).rstrip("/")
    if not api_base.startswith(("https://", "http://")):
        raise ValueError("github.api_base must be an http(s) URL")
    token_file = value.get("token_file")
    if token_file is not None and (not isinstance(token_file, str) or not token_file.startswith("/")):
        raise ValueError("github.token_file must be an absolute path or null")

    def number(key: str, low: float, high: float) -> float:
        item = value.get(key, getattr(defaults, key))
        if isinstance(item, bool) or not isinstance(item, (int, float)) or not low <= float(item) <= high:
            raise ValueError(f"github.{key} must be a number between {low:g} and {high:g}")
        return float(item)

    return GitHubConfig(
        enabled=enabled,
        api_base=api_base,
        token_file=token_file,
        gh=str(value.get("gh", defaults.gh)),
        timeout_seconds=number("timeout_seconds", 1, 120),
        token_cache_seconds=number("token_cache_seconds", 0, 86400),
        backoff_seconds=number("backoff_seconds", 0.1, 3600),
        backoff_max_seconds=number("backoff_max_seconds", 1, 86400),
        max_attempts=int(number("max_attempts", 1, 1000)),
        log_lines=int(number("log_lines", 0, 500)),
        summary_max_bytes=int(number("summary_max_bytes", 0, 60000)),
    )


# ----- submission side -------------------------------------------------------------


def parse_status_target(text: str) -> dict[str, str]:
    match = STATUS_PATTERN.fullmatch(text.strip())
    if match is None:
        raise ValueError(f"--github-status {text!r}: use OWNER/REPO@SHA:CONTEXT with the full commit SHA")
    return {"repo": match.group(1), "sha": match.group(2).lower(), "context": match.group(3)}


def parse_comment_target(text: str) -> dict[str, Any]:
    match = COMMENT_PATTERN.fullmatch(text.strip())
    if match is None:
        raise ValueError(f"--github-comment {text!r}: use OWNER/REPO#PR")
    return {"repo": match.group(1), "pr": int(match.group(2))}


def targets_from_args(statuses: list[str] | None, comment: str | None, summary: str | None,
                      cwd: Path) -> dict[str, Any] | None:
    """The ``github_json`` of a submission, or ``None`` without options."""
    if not statuses and not comment and not summary:
        return None
    document: dict[str, Any] = {"statuses": [parse_status_target(item) for item in statuses or []]}
    contexts = [(s["repo"], s["sha"], s["context"]) for s in document["statuses"]]
    if len(set(contexts)) != len(contexts):
        raise ValueError("--github-status names the same repository, commit and context twice")
    if comment:
        document["comment"] = parse_comment_target(comment)
    if summary:
        path = Path(summary).expanduser()
        document["summary"] = str(path if path.is_absolute() else (cwd / path))
    return document


# ----- what to post ----------------------------------------------------------------


def sanitize(text: str) -> str:
    """No mention and no issue reference can come from job text."""
    text = text.replace("@", "@" + ZWSP)
    return re.sub(r"#(?=\d)", "#" + ZWSP, text)


def job_phase(job: dict[str, Any]) -> str:
    """The transition a job is at, as reported to GitHub."""
    state = str(job.get("state"))
    if state == "queued":
        return f"queued:{int(job.get('preempt_count') or 0)}"
    if state in ("starting", "running", "cancelling"):
        return f"running:{job.get('started_at') or ''}"
    return FINAL


def _seconds(text: str | None) -> float | None:
    if not text:
        return None
    from datetime import datetime

    try:
        return datetime.fromisoformat(str(text).replace("Z", "+00:00")).timestamp()
    except ValueError:
        return None


def runtime_text(job: dict[str, Any]) -> str | None:
    start = _seconds(job.get("started_at"))
    end = _seconds(job.get("finished_at")) or time.time()
    if start is None:
        return None
    seconds = max(0, int(end - start))
    hours, rest = divmod(seconds, 3600)
    minutes, secs = divmod(rest, 60)
    return f"{hours}h{minutes:02d}m" if hours else (f"{minutes}m{secs:02d}s" if minutes else f"{secs}s")


def _size(value: Any) -> str | None:
    if value is None:
        return None
    value = int(value)
    if value >= 1024**3:
        return f"{value / 1024**3:.1f}G"
    return f"{value / 1024**2:.0f}M"


def final_status(job: dict[str, Any]) -> tuple[str, str]:
    """(GitHub state, reason) of a finished job."""
    state = str(job.get("state"))
    result = str(job.get("result") or "")
    code = job.get("exit_code")
    if state == "succeeded":
        return "success", "exit 0"
    if result == "timeout":
        return "error", "time limit exceeded"
    if result == "oom-kill":
        return "error", f"out of memory at its {_size(job.get('mem_bytes')) or '?'} limit"
    if state == "canceled":
        return "error", "preempted and canceled" if result == "preempted" else "canceled"
    if state == "interrupted":
        return "error", "lost: " + str(job.get("message") or "job unit vanished")[:60]
    if state == "failed" and code is not None:
        return "failure", f"exit {code}"
    return "error", str(job.get("message") or state)[:60]


def status_payload(job: dict[str, Any], phase: str) -> tuple[str, str]:
    """(state, description) of a commit status for the job's phase."""
    job_id = int(job["id"])
    if phase.startswith("queued"):
        cores = int(job.get("requested_cores") or 0)
        mem = _size(job.get("mem_bytes"))
        state, text = "pending", f"SimQueue job {job_id} queued ({cores} cores" + (f", {mem}" if mem else "") + ")"
    elif phase.startswith("running"):
        state, text = "pending", f"SimQueue job {job_id} running on CPUs {job.get('cpu_list') or '?'}"
    else:
        state, reason = final_status(job)
        runtime = runtime_text(job)
        peak = _size(job.get("mem_peak_bytes"))
        text = (f"SimQueue job {job_id} {job.get('state')}" + (f" in {runtime}" if runtime else "")
                + (f", peak {peak}" if peak else "") + f"; {reason}")
    text = sanitize(text)
    if len(text) > DESCRIPTION_MAX:
        text = text[: DESCRIPTION_MAX - 1] + "…"
    return state, text


def comment_marker(job_id: int) -> str:
    return f"<!-- simq-job:{int(job_id)} -->"


def _tail(path: str | None, lines: int) -> str:
    if not path or lines <= 0:
        return ""
    try:
        with open(path, "rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(0, size - 65536))
            data = handle.read()
    except OSError:
        return ""
    text = data.decode("utf-8", "replace").splitlines()[-lines:]
    return "\n".join(line[:300] for line in text)


def comment_body(job: dict[str, Any], config: GitHubConfig, node: str, summary_path: str | None) -> str:
    job_id = int(job["id"])
    state, reason = final_status(job) if job_phase(job) == FINAL else ("pending", str(job.get("state")))
    peak = _size(job.get("mem_peak_bytes"))
    request = _size(job.get("mem_bytes"))
    rows = [
        ("State", f"{job.get('state')} ({reason})"),
        ("Runtime", runtime_text(job) or "-"),
        ("Peak memory", (peak or "not measured") + (f" of {request} requested" if request else "")),
        ("Cores", f"{job.get('requested_cores')} (CPUs {job.get('cpu_list') or '-'})"),
        ("Node", node),
    ]
    out = [comment_marker(job_id),
           f"### SimQueue job {job_id}: {sanitize(str(job.get('name') or ''))}",
           "", "| | |", "|---|---|"]
    out += [f"| {key} | {sanitize(str(value)).replace('|', '/')} |" for key, value in rows]
    log = _tail(job.get("log_path"), config.log_lines)
    if log:
        log = sanitize(log).replace("```", "``" + ZWSP + "`")
        out += ["", f"<details><summary>Last {config.log_lines} log lines</summary>", "", "```", log, "```",
                "", "</details>"]
    if summary_path and config.summary_max_bytes > 0:
        try:
            with open(summary_path, "rb") as handle:
                summary = handle.read(config.summary_max_bytes).decode("utf-8", "replace")
        except OSError:
            summary = ""
        if summary.strip():
            out += ["", "---", "", sanitize(summary)]
    out += ["", f"<sub>Posted by SimQueue {__version__} on {node}; edited in place on later posts.</sub>"]
    body = "\n".join(out)
    return body if len(body) <= COMMENT_MAX else body[: COMMENT_MAX - 20] + "\n\n(truncated)"


def build_task(job: dict[str, Any], phase: str, node: str) -> dict[str, Any]:
    """Everything the posting thread needs, copied out of the job row."""
    targets = json.loads(str(job.get("github_json") or "{}"))
    keep = ("id", "name", "state", "result", "exit_code", "message", "started_at", "finished_at", "cpu_list",
            "requested_cores", "mem_bytes", "mem_peak_bytes", "log_path", "preempt_count", "created_at")
    return {
        "job": {key: job.get(key) for key in keep},
        "phase": phase,
        "statuses": list(targets.get("statuses") or []),
        "comment": targets.get("comment") if phase == FINAL else None,
        "summary": targets.get("summary"),
        "comment_id": job.get("github_comment_id"),
        "node": node,
    }


# ----- posting (worker thread) --------------------------------------------------------------


class PostError(Exception):
    def __init__(self, message: str, *, retry: bool, auth: bool = False) -> None:
        super().__init__(message)
        self.retry = retry
        self.auth = auth


class TokenSource:
    """``token_file`` or ``gh auth token``; kept in memory only."""

    def __init__(self, config: GitHubConfig) -> None:
        self.config = config
        self._token: str | None = None
        self._at = 0.0

    def get(self) -> str:
        if self._token and time.monotonic() - self._at < self.config.token_cache_seconds:
            return self._token
        self._token = None
        if self.config.token_file:
            try:
                token = Path(self.config.token_file).read_text(encoding="utf-8").strip()
            except OSError as error:
                raise PostError(f"cannot read the token file ({error.strerror})", retry=True, auth=True) from None
        else:
            try:
                result = subprocess.run([self.config.gh, "auth", "token"], stdin=subprocess.DEVNULL,
                                        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
                                        timeout=30, check=False)
            except (OSError, subprocess.TimeoutExpired) as error:
                raise PostError(f"gh auth token failed ({type(error).__name__})", retry=True, auth=True) from None
            token = result.stdout.strip() if result.returncode == 0 else ""
            if not token:
                raise PostError(f"no GitHub token (gh auth token exited {result.returncode})", retry=True, auth=True)
        if not token or any(c.isspace() for c in token):
            raise PostError("no usable GitHub token", retry=True, auth=True)
        self._token = token
        self._at = time.monotonic()
        return token

    def forget(self) -> None:
        self._token = None


class Client:
    def __init__(self, config: GitHubConfig, tokens: TokenSource) -> None:
        self.config = config
        self.tokens = tokens

    def request(self, method: str, path: str, body: dict[str, Any] | None = None) -> Any:
        # imported on first use: a worker that never posts does not carry them
        import urllib.error
        import urllib.request

        token = self.tokens.get()
        data = None if body is None else json.dumps(body).encode("utf-8")
        request = urllib.request.Request(
            self.config.api_base + path, data=data, method=method,
            headers={"Authorization": f"Bearer {token}", "Accept": "application/vnd.github+json",
                     "X-GitHub-Api-Version": "2022-11-28", "User-Agent": f"simqueue/{__version__}",
                     **({"Content-Type": "application/json"} if data is not None else {})},
        )
        try:
            with urllib.request.urlopen(request, timeout=self.config.timeout_seconds) as response:
                text = response.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as error:
            detail = ""
            try:
                detail = str(json.loads(error.read().decode("utf-8", "replace")).get("message") or "")[:160]
            except (ValueError, OSError, AttributeError):
                pass
            if error.code == 401:
                self.tokens.forget()
            retry = error.code in (401, 408, 429) or error.code >= 500 or (
                error.code == 403 and "rate limit" in detail.lower())
            raise PostError(f"{method} {path}: HTTP {error.code} {detail}".strip(), retry=retry,
                            auth=error.code == 401) from None
        except (urllib.error.URLError, OSError, TimeoutError) as error:
            reason = getattr(error, "reason", error)
            raise PostError(f"{method} {path}: {type(error).__name__}: {str(reason)[:120]}", retry=True) from None
        return json.loads(text) if text.strip() else None


def post_task(client: Client, task: dict[str, Any]) -> dict[str, Any]:
    """Post one job's phase; returns the result for the worker to record."""
    job = task["job"]
    result: dict[str, Any] = {"job_id": int(job["id"]), "phase": task["phase"], "ok": True, "posted": []}
    state, description = status_payload(job, task["phase"])
    for target in task["statuses"]:
        client.request("POST", f"/repos/{target['repo']}/statuses/{target['sha']}",
                       {"state": state, "description": description, "context": target["context"]})
        result["posted"].append(f"status {target['context']}={state}")
    comment = task.get("comment")
    if comment:
        body = comment_body(job, client.config, task["node"], task.get("summary"))
        marker = comment_marker(int(job["id"]))
        comment_id = task.get("comment_id") or find_comment(client, comment["repo"], int(comment["pr"]), marker)
        if comment_id:
            try:
                client.request("PATCH", f"/repos/{comment['repo']}/issues/comments/{int(comment_id)}", {"body": body})
            except PostError as error:
                if "HTTP 404" not in str(error):
                    raise
                comment_id = None  # deleted meanwhile: create it again
        if not comment_id:
            created = client.request("POST", f"/repos/{comment['repo']}/issues/{int(comment['pr'])}/comments",
                                     {"body": body})
            comment_id = int((created or {}).get("id") or 0) or None
        result["comment_id"] = comment_id
        result["posted"].append(f"comment {comment['repo']} PR {comment['pr']}")
    return result


def find_comment(client: Client, repo: str, pr: int, marker: str) -> int | None:
    for page in range(1, 11):
        items = client.request("GET", f"/repos/{repo}/issues/{pr}/comments?per_page=100&page={page}") or []
        for item in items:
            if marker in str(item.get("body") or ""):
                return int(item["id"])
        if len(items) < 100:
            return None
    return None


class Poster:
    """One daemon thread that posts tasks in order; results come back
    through :meth:`results` (the worker records them).  *notify* is called
    from the thread after each result (the worker wakes its loop)."""

    def __init__(self, config: GitHubConfig, notify: Callable[[], None]) -> None:
        self.config = config
        self.notify = notify
        self.client = Client(config, TokenSource(config))
        self._tasks: queue.Queue = queue.Queue()
        self._results: queue.Queue = queue.Queue()
        self.in_flight: set[int] = set()
        self._thread: threading.Thread | None = None

    def submit(self, task: dict[str, Any]) -> None:
        if self._thread is None:
            self._thread = threading.Thread(target=self._run, name="simq-github", daemon=True)
            self._thread.start()
        self.in_flight.add(int(task["job"]["id"]))
        self._tasks.put(task)

    def results(self) -> list[dict[str, Any]]:
        found = []
        while True:
            try:
                item = self._results.get_nowait()
            except queue.Empty:
                return found
            self.in_flight.discard(int(item["job_id"]))
            found.append(item)

    def stop(self) -> None:
        if self._thread is not None:
            self._tasks.put(None)

    def _run(self) -> None:
        while True:
            task = self._tasks.get()
            if task is None:
                return
            try:
                result = post_task(self.client, task)
            except PostError as error:
                result = {"job_id": int(task["job"]["id"]), "phase": task["phase"], "ok": False,
                          "error": str(error), "retry": error.retry}
            except Exception as error:  # noqa: BLE001 - never let the thread die
                result = {"job_id": int(task["job"]["id"]), "phase": task["phase"], "ok": False,
                          "error": f"{type(error).__name__}: {str(error)[:160]}", "retry": True}
            self._results.put(result)
            try:
                self.notify()
            except Exception:  # noqa: BLE001
                pass
