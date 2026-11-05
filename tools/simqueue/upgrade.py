#!/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen
"""Upgrade the node from SimQueue 0.6.3 to 0.7.0, or roll back.  Dry run by default.

    upgrade.py install  [--execute] [--config FILE] [--keep-node-config] [--adoption-timeout S]
                        [--starting-timeout S] [--no-auto-rollback] [--require-doctor]
    upgrade.py rollback --backup DIR [--execute]
    upgrade.py status

0.7.0 makes the worker event-driven and low-priority, tracks every job's
memory and CPU time from its unit's cgroup, adds the memory admission guard,
memory-pressure handling and --mem auto, and reports job state to GitHub
(--github-status, --github-comment).  The database moves from schema 8 to 9:
columns are added (memory, CPU time, GitHub targets and posting state) and an
index; no existing column of any row changes.  The partition layout stays.
Job units are independent of the worker, so running jobs are never touched:
the worker is stopped, its wrapper, unit and node.json are switched, the
database is migrated, and the restarted worker re-adopts every running job
unit by unit name and launch token, and every GPU claim unit.

install (with --execute):
  1. preflight: 0.7.0 runtime sealed; the new node config valid for 0.7.0 on
     this machine and equal to the live node.json apart from the 0.7.0 keys
     (memory, github, launch.idle_seconds, partitions[].min_mem) and the
     notes; every partition name of the live configuration and of every
     queued or active job resolves; wrapper points at 0.6.3; database schema
     8, and a copy of it migrates to 9 with every schema-8 column of every
     row unchanged; records the active jobs, their units and the GPU claim
     units
  2. waits until no job is in state "starting" (--starting-timeout), stops
     simqueue.service (job and claim units keep running) and checks again: a
     job that became "starting" meanwhile makes the install start the 0.6.3
     worker again and stop without changing anything
  3. backs up the database (sqlite3 backup API, source opened read-only),
     the wrapper, the simqueue link, the unit files and node.json and
     worker.log, with sha256 sums in manifest.json
  4. installs the 0.7.0 wrapper, link, simqueue.service (Nice=10,
     CPUWeight=20, IOSchedulingClass=idle, MemoryHigh=96M, MemoryMax=128M)
     and node.json -- from here on every `simq` call runs 0.7.0
  5. migrates the database to schema 9 with the 0.7.0 code and verifies that
     every schema-8 column of every row is unchanged
  6. daemon-reload, enable, start simqueue.service
  7. verifies: the 0.7.0 worker started, every job unit that was active
     before is adopted (worker.log "adopted job N") or has been recorded as
     finished, every GPU claim unit that was active is adopted or has ended,
     the GPU line state file is published by 0.7.0, the wake FIFO exists,
     `simq doctor` and `simq partitions` run
  On a failed verification the install rolls itself back (unless
  --no-auto-rollback).  Receipt: <backup>/install-receipt.json

rollback --backup DIR (with --execute):
  1. checks the backup against its manifest
  2. stops simqueue.service (job units and GPU claim units keep running; the
     0.6.3 worker adopts both)
  3. keeps a copy of the schema-9 database in the backup directory
  4. restores the wrapper, link, unit files and node.json byte for byte
  5. sets the database schema back to 8 in place: rows, including jobs that
     0.7.0 queued or ran, stay; the 0.7.0 columns stay and are ignored by
     0.6.3; checks that the 0.6.3 code opens it
  6. daemon-reload, start simqueue.service (0.6.3) and verifies that it
     adopted every active job unit
  Jobs with GitHub targets whose final post is still due are listed: 0.6.3
  does not post them.  Receipt: <backup>/rollback-receipt.json

Job units started by 0.7.0 run the 0.7.0 jobwrap.py: keep the 0.7.0 runtime
directory after a rollback until they have ended (and the 0.6.3 directory
after the install, for the same reason).

--home, --project-root, --systemctl, --service, --runtime and --old-runtime
exist so that the procedure can be tested in a sandbox; the defaults are the
real paths.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import time
from contextlib import closing
from datetime import UTC, datetime
from pathlib import Path


RUNTIME = Path(__file__).resolve().parent
OLD_RUNTIME = Path("/home/treefor/projects/runtimes/simqueue-0.6.3-20261105")
FROM_VERSION, TO_VERSION = "0.6.3", "0.7.0"
PYTHON = "/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen"
REAL_HOME = "/home/treefor"
UNIT_FILES = ("simqueue.service", "simqueue-spheral-dispatch.path", "simqueue-spheral-dispatch.service",
              "simqborrow.slice")
# What 0.7.0 may change in the node configuration; everything else must match.
NEW_TOP_KEYS = {"memory", "github", "notes"}
NEW_LAUNCH_KEYS = {"idle_seconds"}
NEW_PARTITION_KEYS = {"min_mem"}
NEW_BORROW_KEYS = {"notes"}
NEW_LOCK_KEYS = {"notes"}
NEW_LINE_KEYS = {"notes"}
SCHEMA_OLD, SCHEMA_NEW = 8, 9
SCHEMA = SCHEMA_OLD
ACTIVE_STATES = ("starting", "running", "cancelling")
# Run with the 0.6.3 code (python_in): both layouts as the scheduler reads them.
LAYOUT_PROBE = r"""
import json
from pathlib import Path
from simqueue.config import format_cpu_list
from simqueue.nodeconfig import load_node_config, parse_node_config
live = parse_node_config(json.loads(Path(sys.argv[2]).read_text()), path=Path(sys.argv[2]), check_topology=False)
new = load_node_config(Path(sys.argv[3]))
names = json.loads(sys.argv[4])
def layout(config):
    return [{"name": p.name, "cpus": format_cpu_list(p.cpus), "gpus": len(p.gpus), "aliases": list(p.aliases)}
            for p in config.partitions.values()]
print(json.dumps({
    "live": layout(live), "new": layout(new),
    "live_cpus": list(live.node_cpus()), "new_cpus": list(new.node_cpus()),
    "live_gpus": sorted(live.gpus), "new_gpus": sorted(new.gpus),
    "live_names": sorted(set(live.partitions) | set(live.aliases)),
    "resolves": {name: (new.canonical_partition(name) if new.find_partition(name) is not None else None)
                 for name in sorted(set(names) | set(live.partitions) | set(live.aliases))},
    "default": new.default_partition,
    "backfill_waiting_owners": new.borrow.backfill_waiting_owners,
}))
"""


def now_tag() -> str:
    return datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")


def sha256(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return None


def comparable_config(document: dict) -> dict:
    """The node configuration without what 0.7.0 may change: the memory and
    github sections, launch.idle_seconds, partition min_mem and the notes."""
    result = {k: v for k, v in document.items() if k not in NEW_TOP_KEYS}
    for key, drop in (("lock_files", NEW_LOCK_KEYS), ("borrow", NEW_BORROW_KEYS), ("launch", NEW_LAUNCH_KEYS),
                      ("gpu_line", NEW_LINE_KEYS)):
        if isinstance(result.get(key), dict):
            result[key] = {k: v for k, v in result[key].items() if k not in drop}
    if isinstance(result.get("partitions"), list):
        result["partitions"] = [{k: v for k, v in part.items() if k not in NEW_PARTITION_KEYS}
                                if isinstance(part, dict) else part for part in result["partitions"]]
    return result


def config_differences(live: dict, new: dict) -> list[str]:
    """Keys (dotted) whose values differ between two comparable configs."""
    found: list[str] = []

    def walk(a, b, path: str) -> None:
        if isinstance(a, dict) and isinstance(b, dict):
            for key in sorted(set(a) | set(b)):
                walk(a.get(key), b.get(key), f"{path}.{key}" if path else str(key))
        elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
            for index, (x, y) in enumerate(zip(a, b)):
                walk(x, y, f"{path}[{index}]")
        elif a != b:
            found.append(path or "(document)")

    walk(live, new, "")
    return found


def last_line(text: str) -> str:
    """The last line of a traceback: the error itself."""
    lines = [line for line in text.strip().splitlines() if line.strip()]
    return lines[-1][-400:] if lines else "(no message)"


def layout_text(layout: list[dict]) -> str:
    return ", ".join(
        f"{p['name']}={p['cpus']}" + (f"+{p['gpus']}gpu" if p["gpus"] else "")
        + (f" (aliases {', '.join(p['aliases'])})" if p["aliases"] else "")
        for p in layout
    )


def layout_problems(probe: dict, job_names: dict[str, list[int]]) -> list[str]:
    """Why the new partition layout cannot replace the live one."""
    problems = []
    lost = sorted(set(probe["live_cpus"]) - set(probe["new_cpus"]))
    if lost:
        problems.append(f"the new node config drops CPU(s) {lost} of the live one")
    gone = sorted(set(probe["live_gpus"]) - set(probe["new_gpus"]))
    if gone:
        problems.append(f"the new node config drops GPU(s) {gone}")
    for name in probe["live_names"]:
        if probe["resolves"].get(name) is None:
            problems.append(f"partition name {name!r} of the live node config is neither a partition nor an alias "
                            "in the new one (scripts using -p " + name + " would fail)")
    for name, jobs in sorted(job_names.items()):
        if probe["resolves"].get(name) is None:
            problems.append(f"job(s) {jobs} recorded under partition {name!r} could not be scheduled with the new "
                            "node config")
    return problems


def database_summary(path: Path) -> dict:
    wal = Path(f"{path}-wal")
    uri = f"file:{path}?mode=ro" if wal.exists() and wal.stat().st_size else f"file:{path}?immutable=1"
    with closing(sqlite3.connect(uri, uri=True)) as connection:
        version = connection.execute("PRAGMA user_version").fetchone()[0]
        states = dict(connection.execute("SELECT state, COUNT(*) FROM jobs GROUP BY state").fetchall())
        count = connection.execute("SELECT COUNT(*) FROM jobs").fetchone()[0]
        columns = [row[1] for row in connection.execute("PRAGMA table_info(jobs)").fetchall()]
        active = [
            dict(zip(("id", "name", "state", "unit_name", "launch_token"), row))
            for row in connection.execute(
                "SELECT id, name, state, " + ("unit_name" if "unit_name" in columns else "NULL")
                + ", launch_token FROM jobs WHERE state IN ('starting', 'running', 'cancelling') ORDER BY id"
            )
        ]
        # queued and active jobs by recorded partition name (NULL: a 0.5 job, default partition)
        partition = "partition" if "partition" in columns else "NULL"
        unfinished: dict[str, list[int]] = {}
        for job_id, name in connection.execute(
            f"SELECT id, {partition} FROM jobs WHERE state IN ('queued', 'starting', 'running', 'cancelling') "
            "ORDER BY id"
        ):
            unfinished.setdefault(name or "", []).append(int(job_id))
    return {"user_version": version, "jobs": count, "states": states, "columns": columns, "active": active,
            "starting": [job["id"] for job in active if job["state"] == "starting"],
            "unfinished_by_partition": unfinished}


def row_digests(path: Path, columns: list[str]) -> dict[int, str]:
    """sha256 of the given columns of every job row, by id."""
    with closing(sqlite3.connect(f"file:{path}?mode=ro", uri=True)) as connection:
        rows = connection.execute(f"SELECT id, {', '.join(columns)} FROM jobs ORDER BY id").fetchall()
    return {int(row[0]): hashlib.sha256(json.dumps(row[1:], default=str).encode()).hexdigest() for row in rows}


def _envpolicy(runtime: Path):
    """The 0.7.0 environment policy (secret patterns, scrubbing)."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("simq_envpolicy", runtime / "source" / "simqueue" / "envpolicy.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PATTERN_GROUPS = ("TOKEN", "SECRET", "PASSWORD", "KEY", "GITHUB_", "ACTIONS_", "RUNNER_")


def _group_of(name: str) -> str:
    upper = name.upper()
    for group in PATTERN_GROUPS:
        if (group.endswith("_") and upper.startswith(group)) or (not group.endswith("_") and group in upper):
            return group
    return "other"


def secret_census(database: Path, state: Path, runtime: Path) -> dict:
    """How many stored environments carry secret-pattern values: counts per
    job state and per pattern, and of job spec files.  Never values."""
    policy = _envpolicy(runtime)
    census: dict = {"jobs_scanned": 0, "jobs_with_secrets": 0, "by_state": {}, "values_by_pattern": {},
                    "explicit_secrets_of_unfinished_jobs": 0, "ensemble_members_with_secrets": 0,
                    "spec_files_scanned": 0, "spec_files_with_secrets": 0}
    with closing(sqlite3.connect(f"file:{database}?mode=ro", uri=True)) as connection:
        for job_state, env_json, override_json in connection.execute(
                "SELECT state, env_json, env_override_keys FROM jobs"):
            census["jobs_scanned"] += 1
            try:
                env = json.loads(env_json or "{}")
                explicit = set(json.loads(override_json or "[]"))
            except ValueError:
                continue
            names = [k for k, v in env.items() if policy.is_secret(k) and v != policy.REDACTED]
            if not names:
                continue
            census["jobs_with_secrets"] += 1
            census["by_state"][job_state] = census["by_state"].get(job_state, 0) + 1
            for name in names:
                group = _group_of(name)
                census["values_by_pattern"][group] = census["values_by_pattern"].get(group, 0) + 1
                if name in explicit and job_state not in ("succeeded", "failed", "canceled", "interrupted"):
                    census["explicit_secrets_of_unfinished_jobs"] += 1
        for (env_json,) in connection.execute("SELECT env_json FROM ensemble_members"):
            if policy.scrub_json(env_json)[1]:
                census["ensemble_members_with_secrets"] += 1
    for spec in (state / "jobs").glob("*/spec.json"):
        census["spec_files_scanned"] += 1
        try:
            env = json.loads(spec.read_text()).get("env") or {}
        except (OSError, ValueError, AttributeError):
            continue
        if any(policy.is_secret(k) and v != policy.REDACTED for k, v in env.items()):
            census["spec_files_with_secrets"] += 1
    return census


def census_text(census: dict) -> str:
    return (f"{census['jobs_with_secrets']} of {census['jobs_scanned']} stored job environments carry "
            f"secret-pattern values (by state {census['by_state']}; values by pattern "
            f"{census['values_by_pattern']}; {census['explicit_secrets_of_unfinished_jobs']} given with -e to "
            f"unfinished jobs, kept in their secrets files); ensemble members {census['ensemble_members_with_secrets']}; "
            f"job spec files {census['spec_files_with_secrets']} of {census['spec_files_scanned']}")


def wrapper_runtime(text: str) -> Path | None:
    match = re.search(r"exec \S+ (\S+)/source/simq ", text)
    return Path(match.group(1)) if match else None


def line_state_path(document: dict) -> str | None:
    """The published GPU line state file of a node configuration (0.6.2)."""
    line = document.get("gpu_line") or {}
    if line.get("state_path"):
        return str(line["state_path"])
    locks = document.get("lock_files") or {}
    for paths in (locks.get("gpu_extra") or {}).values():
        for path in paths or []:
            return os.path.join(os.path.dirname(str(path)), "simq-gpu-line.json")
    if locks.get("gpu_path"):
        return os.path.join(os.path.dirname(str(locks["gpu_path"])), "simq-gpu-line.json")
    return None


class Upgrade:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.execute = bool(getattr(args, "execute", False))
        self.home = Path(args.home)
        self.runtime = Path(args.runtime)
        self.old_runtime = Path(args.old_runtime)
        self.project_root = Path(args.project_root)
        self.bin = self.home / ".local" / "bin"
        self.units = self.home / ".config" / "systemd" / "user"
        self.state = self.home / ".local" / "state" / "simqueue"
        self.backups = self.home / ".local" / "state" / "simqueue-backups"
        self.database = self.state / "queue.sqlite3"
        self.service = args.service
        self.problems: list[str] = []
        self.steps: list[dict] = []

    # ----- utilities ------------------------------------------------------

    def say(self, text: str) -> None:
        print(("" if self.execute else "[dry-run] ") + text, flush=True)

    def step(self, name: str, **detail) -> None:
        self.steps.append({"step": name, "at": now_tag(), **detail})

    def problem(self, text: str) -> None:
        self.problems.append(text)
        print("PROBLEM: " + text, flush=True)

    def systemctl(self, *arguments: str, check: bool = True) -> subprocess.CompletedProcess[str]:
        command = [self.args.systemctl, "--user", *arguments]
        read_only = arguments and arguments[0] in ("show", "is-active", "is-enabled", "list-units", "cat")
        if not self.execute and not read_only:
            self.say("would run: " + " ".join(command))
            return subprocess.CompletedProcess(command, 0, "", "")
        result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        if check and result.returncode != 0:
            raise RuntimeError(f"{' '.join(command)} failed: {result.stderr.strip()}")
        return result

    def unit_property(self, unit: str, prop: str) -> str:
        return self.systemctl("show", unit, f"--property={prop}", "--value", check=False).stdout.strip()

    def env_for(self) -> dict[str, str]:
        env = {k: v for k, v in os.environ.items() if not k.startswith("SIMQUEUE_")}
        env.update({
            "SIMQUEUE_PROJECT_ROOT": str(self.project_root),
            "SIMQUEUE_STATE_DIR": str(self.state),
            "SIMQUEUE_CONFIG": str(self.state / "node.json"),
            "SIMQUEUE_NO_AUTOSTART": "1",
        })
        return env

    def probe_env(self) -> dict[str, str]:
        """Environment of read-only probes: a throwaway state directory, so a
        probe can never open (and migrate) the live database."""
        env = self.env_for()
        scratch = Path(os.environ.get("TMPDIR") or "/tmp") / f"simq-upgrade-probe-{os.getpid()}"
        if not scratch.exists():
            scratch.mkdir(mode=0o700)
            import atexit
            atexit.register(shutil.rmtree, scratch, True)
        env["SIMQUEUE_STATE_DIR"] = str(scratch)
        return env

    def python_in(self, runtime: Path, code: str, *arguments: str, probe: bool = False,
                  env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [PYTHON, "-c", "import sys; sys.path.insert(0, sys.argv[1]);" + code, str(runtime / "source"), *arguments],
            env=env or (self.probe_env() if probe else self.env_for()), text=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False, timeout=300,
        )

    def simq(self, runtime: Path, *arguments: str, timeout: float = 120) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [PYTHON, str(runtime / "source" / "simq"), *arguments], env=self.env_for(), cwd=str(self.project_root),
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout, check=False,
        )

    def wrapper_text(self, runtime: Path, version: str) -> str:
        layout = ("# SimQueue {version} node scheduler: one pool of every physical core and the GPU\n"
                  "# (node.json; the old partition names are aliases); each job runs in its own\n"
                  "# systemd user unit; `simq run` waits for a job and returns its exit code;\n")
        if version == TO_VERSION:
            layout += ("# the worker is event-driven and tracks job memory from the units' cgroups;\n"
                       "# --mem auto; --github-status/--github-comment report job state to GitHub;\n")
        return (
            "#!/bin/sh\n"
            + layout.format(version=version)
            + "# borrow jobs run in simqborrow.slice at idle CPU priority; queued GPU jobs\n"
            "# win the GPU lease files through the GPU line (gpu-lease for direct use).\n"
            "set -eu\n"
            "unset SIMQUEUE_PCORES\n"
            f"export SIMQUEUE_PROJECT_ROOT={self.project_root}\n"
            f"export SIMQUEUE_STATE_DIR={self.state}\n"
            f"export SIMQUEUE_CONFIG={self.state / 'node.json'}\n"
            f"exec {PYTHON} {runtime / 'source' / 'simq'} \"$@\"\n"
        )

    def unit_text(self, runtime: Path, version: str = TO_VERSION) -> str:
        return (
            "[Unit]\n"
            f"Description=SimQueue {version} node scheduler (jobs run in their own units)\n"
            f"Documentation=file://{runtime / 'OPERATIONS.md'}\n"
            "\n"
            "[Service]\n"
            "Type=simple\n"
            f"WorkingDirectory={self.project_root}\n"
            f"ExecStart={self.bin / 'simq'} worker --wait-for-lock\n"
            f"Environment=SIMQUEUE_PROJECT_ROOT={self.project_root}\n"
            f"Environment=SIMQUEUE_STATE_DIR={self.state}\n"
            f"Environment=SIMQUEUE_CONFIG={self.state / 'node.json'}\n"
            "# No CPUAffinity: the scheduler places each job unit itself; the worker\n"
            "# must not be confined to one partition.\n"
            "# 0.7.0: the worker is event-driven and runs at low priority; job units\n"
            "# and GPU claim units are started by the user manager and inherit none\n"
            "# of this.  MemoryMax bounds the worker (about 30 MB) only.\n"
            "Nice=10\n"
            "CPUWeight=20\n"
            "IOSchedulingClass=idle\n"
            "MemoryHigh=96M\n"
            "MemoryMax=128M\n"
            "Restart=on-failure\n"
            "RestartSec=2\n"
            "# Job units and GPU claim units are separate units: stopping this\n"
            "# service never stops them.\n"
            "KillMode=mixed\n"
            "TimeoutStopSec=30\n"
            "UMask=0077\n"
            "\n"
            "[Install]\n"
            "WantedBy=default.target\n"
        )

    def write_file(self, path: Path, data: bytes | str, mode: int) -> None:
        if not self.execute:
            self.say(f"would write {path} ({len(data)} bytes, mode {mode:o})")
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".simq-new")
        if isinstance(data, str):
            temporary.write_text(data, encoding="utf-8")
        else:
            temporary.write_bytes(data)
        temporary.chmod(mode)
        os.replace(temporary, path)

    def point_link(self, target: str | None) -> None:
        link = self.bin / "simqueue"
        if not self.execute:
            self.say(f"would point {link} at {target}")
            return
        temporary = self.bin / "simqueue.simq-new"
        if temporary.is_symlink() or temporary.exists():
            temporary.unlink()
        if target is None:
            link.unlink(missing_ok=True)
            return
        temporary.symlink_to(target)
        os.replace(temporary, link)

    def wait_service(self, active: bool, timeout: float = 30) -> None:
        deadline = time.monotonic() + timeout
        wanted = ("active",) if active else ("inactive", "failed")
        while self.unit_property(self.service, "ActiveState") not in wanted:
            if time.monotonic() > deadline:
                raise RuntimeError(f"{self.service} did not become {'active' if active else 'inactive'}")
            time.sleep(0.3)

    def worker_log_since(self, offset: int) -> str:
        path = self.state / "worker.log"
        try:
            with path.open("rb") as handle:
                handle.seek(offset)
                return handle.read().decode("utf-8", "replace")
        except OSError:
            return ""

    def worker_log_size(self) -> int:
        try:
            return (self.state / "worker.log").stat().st_size
        except OSError:
            return 0

    def unit_active(self, unit: str) -> bool:
        """Claim units are real units even in a sandbox (names from its state)."""
        result = subprocess.run(["systemctl", "--user", "is-active", "--", unit], text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False, timeout=60)
        return result.stdout.strip() in ("active", "activating", "deactivating", "reloading")

    def verify_adoption(self, *, active_before: list[dict], log_offset: int, version: str,
                        runtime: Path, timeout: float, state_path: str | None = None,
                        claims_before: list[str] = ()) -> dict:
        """Wait until the new worker has started and adopted (or recorded)
        every job that was active before the switch, and adopted every GPU
        claim unit that was active (or seen it end); with *state_path*, also
        until that worker has published the GPU line state."""
        deadline = time.monotonic() + timeout
        result: dict = {"version": version, "jobs": {}, "claims": {}}
        started_at = time.time()
        while True:
            log = self.worker_log_since(log_offset)
            started = f"node worker started; simqueue {version} from {runtime}" in log
            summary = database_summary(self.database)
            states = {job["id"]: job for job in summary["active"]}
            pending = []
            for job in active_before:
                job_id = job["id"]
                adopted = f"adopted job {job_id} " in log
                finished = job_id not in states
                result["jobs"][job_id] = {"adopted": adopted, "finished_meanwhile": finished and not adopted,
                                          "unit": job.get("unit_name")}
                if not adopted and not finished:
                    pending.append(job_id)
            for unit in claims_before:
                adopted = f"adopted GPU claim {unit} " in log
                ended = not adopted and not self.unit_active(unit)
                result["claims"][unit] = {"adopted": adopted, "ended_meanwhile": ended}
                if not adopted and not ended:
                    pending.append(unit)
            result["worker_started"] = started
            result["service_active"] = self.unit_property(self.service, "ActiveState") == "active"
            published = True
            if state_path is not None:
                document = None
                try:
                    document = json.loads(Path(state_path).read_text())
                except (OSError, ValueError):
                    pass
                published = bool(document and document.get("version") == version
                                 and float(document.get("heartbeat_unix") or 0) >= started_at - 1)
                result["gpu_line_state"] = {"path": state_path, "published": published,
                                            "version": (document or {}).get("version"),
                                            "updated_at": (document or {}).get("updated_at")}
            if started and not pending and result["service_active"] and published:
                result["ok"] = True
                return result
            if time.monotonic() > deadline:
                result["ok"] = False
                result["pending"] = pending
                return result
            time.sleep(0.5)

    def sandbox_guard(self) -> None:
        if str(self.home) != REAL_HOME and self.args.systemctl == "systemctl":
            raise SystemExit("a sandbox --home needs a sandbox --systemctl so that the real service is never touched")

    # ----- install --------------------------------------------------------

    def preflight_install(self) -> dict:
        print(f"SimQueue {TO_VERSION} runtime: {self.runtime}")
        self.layout = None
        sealed = self.runtime / "sealed-files.json"
        if sealed.exists():
            manifest = json.loads(sealed.read_text())
            bad = [name for name, digest in manifest["files"].items() if sha256(self.runtime / name) != digest]
            if bad:
                self.problem(f"sealed runtime files changed: {bad}")
        elif not self.args.allow_unsealed:
            self.problem("sealed-files.json missing from the runtime")
        if not os.access(PYTHON, os.X_OK):
            self.problem(f"pinned interpreter not executable: {PYTHON}")
        wrapper = self.bin / "simq"
        text = wrapper.read_text(errors="replace") if wrapper.exists() else ""
        current = wrapper_runtime(text)
        print(f"wrapper {wrapper} runs {current}")
        if current is None:
            self.problem(f"cannot read the runtime from {wrapper}")
        elif current == self.runtime:
            self.problem(f"{TO_VERSION} already installed")
        elif current != self.old_runtime:
            self.problem(f"wrapper runs {current}, expected {self.old_runtime}")
        if not (self.old_runtime / "source" / "simqueue" / "jobwrap.py").exists():
            self.problem(f"{FROM_VERSION} runtime missing: {self.old_runtime}")
        live_config = self.state / "node.json"
        new_config = Path(self.args.config)
        try:
            live_doc = json.loads(live_config.read_text())
            new_doc = json.loads(new_config.read_text())
        except (OSError, ValueError) as error:
            self.problem(f"cannot read node configs: {error}")
            live_doc = new_doc = {}
        check = self.python_in(
            self.runtime,
            "from pathlib import Path; from simqueue.nodeconfig import load_node_config;"
            "c = load_node_config(Path(sys.argv[2])); l = c.line; m = c.memory; g = c.github;"
            "print(c.node, ' '.join(f'{p.name}={p.cpus} default_mem={p.default_mem >> 20}M' for p in c.partitions.values()),"
            "f'gpu_line enabled={l.enabled}; memory sample={m.sample_seconds}s guard margin={m.admission_margin}'"
            "f' high={m.high_water_available >> 20}M/{m.high_water_percent}% low={m.low_water_available >> 20}M/'"
            "f'{m.low_water_percent}%; github enabled={g.enabled} api={g.api_base}'"
            "f' token={\"file\" if g.token_file else \"gh auth token\"}; launch poll={c.launch.poll_seconds}s'"
            "f' idle={c.launch.idle_seconds}s')",
            str(new_config if not self.args.keep_node_config else live_config), probe=True,
        )
        if check.returncode != 0:
            self.problem(f"node config invalid for {TO_VERSION}: {last_line(check.stderr)}")
        else:
            print(f"{TO_VERSION} node config ok: {check.stdout.strip()}")
        if not self.args.keep_node_config and live_doc and new_doc:
            differences = config_differences(comparable_config(live_doc), comparable_config(new_doc))
            if differences:
                self.problem(f"{new_config} differs from the live {live_config} beyond the {TO_VERSION} keys "
                             f"({', '.join(differences[:8])}); reconcile them or use --keep-node-config")
            else:
                print(f"node config: equal to the live one apart from the {TO_VERSION} keys "
                      "(memory, github, launch.idle_seconds, partitions[].min_mem, notes)")
        if not self.database.exists():
            self.problem(f"database missing: {self.database}")
            return {}
        summary = database_summary(self.database)
        print(f"database: schema {summary['user_version']}, {summary['jobs']} jobs, states {summary['states']}")
        if summary["user_version"] != SCHEMA_OLD:
            self.problem(f"database schema is {summary['user_version']}, expected {SCHEMA_OLD} ({FROM_VERSION})")
        else:
            self.rehearse_migration(summary)
        for job in summary["active"]:
            print(f"active job {job['id']} ({job['name']}): {job['state']}, unit {job['unit_name']}")
        ensembles = self.active_ensembles()
        if ensembles:
            print(f"NOTE: active ensemble job(s) {ensembles}: their {FROM_VERSION} supervisors opened the database "
                  "before the migration and keep working; a new 0.6.3 process could not open schema 9")
        unfinished = summary["unfinished_by_partition"]
        names: dict[str, list[int]] = {}
        for name, jobs in unfinished.items():
            names.setdefault(name or str(live_doc.get("default_partition") or ""), []).extend(jobs)
        print("queued and active jobs by recorded partition: "
              + (", ".join(f"{name}: {jobs}" for name, jobs in sorted(names.items())) or "none"))
        if live_doc:
            probe = self.python_in(self.runtime, LAYOUT_PROBE, str(live_config),
                                   str(new_config if not self.args.keep_node_config else live_config),
                                   json.dumps(sorted(n for n in names if n)), probe=True)
            if probe.returncode != 0:
                self.problem(f"cannot compare the partition layouts: {last_line(probe.stderr)}")
            else:
                self.layout = json.loads(probe.stdout)
                print(f"partition layout (unchanged by {TO_VERSION}): {layout_text(self.layout['new'])}")
                if self.layout["live"] != self.layout["new"]:
                    self.problem(f"the partition layout differs: live {layout_text(self.layout['live'])}")
                for text in layout_problems(self.layout, {n: j for n, j in names.items() if n}):
                    self.problem(text)
        claims = self.claim_units()
        print(f"GPU claim units active: {', '.join(claims) or 'none'} (the {TO_VERSION} worker adopts them)")
        gh = shutil.which("gh")
        print(f"GitHub reporting: gh {'found at ' + gh if gh else 'NOT found (only a token_file would work)'}; "
              "posts happen only for jobs submitted with --github-status/--github-comment")
        exec_start = self.unit_property(self.service, "ExecStart")
        print(f"{self.service} ExecStart: {exec_start[:160]}")
        if str(self.bin / "simq") not in exec_start:
            self.problem(f"{self.service} does not run {self.bin / 'simq'}")
        return summary

    def rehearse_migration(self, summary: dict) -> None:
        """Migrate a copy of the database with the 0.7.0 code and check that
        no schema-8 column of any row changes (the live file is only read)."""
        scratch = Path(self.probe_env()["SIMQUEUE_STATE_DIR"]) / "rehearsal"
        if scratch.exists():
            shutil.rmtree(scratch)
        scratch.mkdir(mode=0o700)
        copy = scratch / "queue.sqlite3"
        with closing(sqlite3.connect(f"file:{self.database}?mode=ro", uri=True)) as source, \
                closing(sqlite3.connect(copy)) as target:
            source.backup(target)
        env = self.env_for()
        env["SIMQUEUE_STATE_DIR"] = str(scratch)
        result = self.python_in(self.runtime, "from simqueue.store import QueueStore; QueueStore(); print('ok')",
                                env=env)
        if result.returncode != 0:
            self.problem(f"the {TO_VERSION} code cannot migrate a copy of the database: {last_line(result.stderr)}")
            return
        after = database_summary(copy)
        columns = [c for c in summary["columns"] if c != "env_json"]
        old = row_digests(self.database, columns)
        new = row_digests(copy, columns)
        changed = [job_id for job_id, digest in old.items() if new.get(job_id) != digest]
        added = [c for c in after["columns"] if c not in summary["columns"]]
        wrong_env = self.env_scrub_problems(self.database, copy)
        census = secret_census(self.database, self.state, self.runtime)
        self.census = census
        if after["user_version"] != SCHEMA_NEW or changed or len(new) != len(old) or wrong_env:
            self.problem(f"migration rehearsal: schema {after['user_version']}, rows changed {changed[:10]}, "
                         f"environments not scrubbed as expected {wrong_env[:10]}")
        else:
            print(f"migration rehearsal on a copy: schema {SCHEMA_OLD} -> {SCHEMA_NEW}, {len(old)} rows, every "
                  f"schema-{SCHEMA_OLD} column unchanged but env_json, whose secret-pattern values are "
                  f"redacted; {len(added)} column(s) added: {', '.join(added)}")
            print("environment scrub (counts only): " + census_text(census))
        shutil.rmtree(scratch, ignore_errors=True)

    def env_scrub_problems(self, before: Path, after: Path) -> list[int]:
        """Jobs whose env_json after the migration is not the scrubbed
        env_json of before (rows added meanwhile are not compared)."""
        policy = _envpolicy(self.runtime)
        with closing(sqlite3.connect(f"file:{before}?mode=ro", uri=True)) as a, \
                closing(sqlite3.connect(f"file:{after}?mode=ro", uri=True)) as b:
            old = dict(a.execute("SELECT id, env_json FROM jobs"))
            new = dict(b.execute("SELECT id, env_json FROM jobs"))
        return [job_id for job_id, text in old.items() if policy.scrub_json(text)[0] != new.get(job_id)]

    def scrub_spec_files(self) -> int:
        """Redact the spec.json of every job that is not active (active
        jobs' specs are redacted by the worker when they end)."""
        policy = _envpolicy(self.runtime)
        with closing(sqlite3.connect(f"file:{self.database}?mode=ro", uri=True)) as connection:
            active = {str(row[0]) for row in connection.execute(
                "SELECT job_dir FROM jobs WHERE state IN ('starting', 'running', 'cancelling') AND job_dir IS NOT NULL")}
        count = 0
        for spec in (self.state / "jobs").glob("*/spec.json"):
            if str(spec.parent) in active:
                continue
            if policy.scrub_spec(spec):
                count += 1
        return count

    def scrub_backup(self, path: Path) -> int:
        """Redact the backup copy's stored environments (its schema stays)."""
        policy = _envpolicy(self.runtime)
        count = 0
        with closing(sqlite3.connect(path)) as connection:
            for job_id, text in connection.execute("SELECT id, env_json FROM jobs").fetchall():
                cleaned, names = policy.scrub_json(text)
                if names:
                    connection.execute("UPDATE jobs SET env_json = ? WHERE id = ?", (cleaned, job_id))
                    count += 1
            for job_id, index, text in connection.execute(
                    "SELECT job_id, member_index, env_json FROM ensemble_members").fetchall():
                cleaned, names = policy.scrub_json(text)
                if names:
                    connection.execute("UPDATE ensemble_members SET env_json = ? WHERE job_id = ? AND member_index = ?",
                                       (cleaned, job_id, index))
            connection.commit()
            connection.execute("VACUUM")
        return count

    def active_ensembles(self) -> list[int]:
        try:
            with closing(sqlite3.connect(f"file:{self.database}?mode=ro", uri=True)) as connection:
                rows = connection.execute(
                    "SELECT j.id FROM jobs j JOIN ensembles e ON e.job_id = j.id "
                    "WHERE j.state IN ('starting', 'running', 'cancelling') ORDER BY j.id").fetchall()
        except sqlite3.Error:
            return []
        return [int(row[0]) for row in rows]

    def github_due(self) -> list[int]:
        """Jobs whose GitHub posts 0.7.0 has not finished (rollback note)."""
        try:
            with closing(sqlite3.connect(f"file:{self.database}?mode=ro", uri=True)) as connection:
                columns = {row[1] for row in connection.execute("PRAGMA table_info(jobs)")}
                if "github_json" not in columns:
                    return []
                rows = connection.execute("SELECT id FROM jobs WHERE github_json IS NOT NULL AND github_done = 0 "
                                          "ORDER BY id").fetchall()
        except sqlite3.Error:
            return []
        return [int(row[0]) for row in rows]

    def wait_no_starting(self, timeout: float) -> list[int]:
        """Jobs still in state "starting" after waiting up to *timeout* s."""
        deadline = time.monotonic() + timeout
        while True:
            starting = database_summary(self.database)["starting"]
            if not starting or time.monotonic() > deadline:
                return starting
            time.sleep(0.5)

    def install(self) -> int:
        self.sandbox_guard()
        summary = self.preflight_install()
        if self.problems:
            print(f"install refused: {len(self.problems)} problem(s)")
            return 1
        tag = now_tag()
        backup = self.backups / f"{FROM_VERSION}-to-{TO_VERSION}-{tag}"
        new_doc = json.loads(Path(self.args.config if not self.args.keep_node_config
                                  else self.state / "node.json").read_text())
        state_path = line_state_path(new_doc)
        self.say(f"backup directory: {backup}")
        # no job may be in "starting" at the switch (a launch in progress)
        starting = self.wait_no_starting(self.args.starting_timeout if self.execute else 0)
        if starting and self.execute:
            print(f"install refused: job(s) {starting} still starting after {self.args.starting_timeout:.0f} s; retry")
            return 1
        if starting:
            self.say(f"job(s) {starting} are starting now; --execute waits up to {self.args.starting_timeout:.0f} s")
        self.say(f"stopping {self.service} (job units keep running)")
        self.systemctl("stop", self.service)
        if self.execute:
            self.wait_service(False)
            starting = database_summary(self.database)["starting"]
            if starting:
                self.systemctl("start", self.service)
                print(f"install refused: job(s) {starting} became starting while the {FROM_VERSION} worker stopped; "
                      f"the {FROM_VERSION} worker was started again and nothing was changed; retry")
                return 1
        self.step(f"stopped-{FROM_VERSION}-worker", starting_at_switch=[])
        # backup
        files: dict[str, dict] = {}
        if self.execute:
            backup.mkdir(parents=True, mode=0o700)
            with closing(sqlite3.connect(f"file:{self.database}?mode=ro", uri=True)) as source, \
                    closing(sqlite3.connect(backup / "queue.sqlite3")) as target:
                source.backup(target)
            for name, path in (
                ("simq", self.bin / "simq"),
                ("node.json", self.state / "node.json"),
                ("worker.log", self.state / "worker.log"),
                *((f"units/{unit}", self.units / unit) for unit in UNIT_FILES),
            ):
                if path.exists():
                    destination = backup / name
                    destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
                    shutil.copy2(path, destination)
                    files[name] = {"source": str(path), "sha256": sha256(destination),
                                   "mode": oct(path.stat().st_mode & 0o777)}
            link = self.bin / "simqueue"
            before = database_summary(backup / "queue.sqlite3")
            manifest = {
                "schema": "simqueue-upgrade-backup-v1",
                "from": FROM_VERSION, "to": TO_VERSION,
                "created": tag,
                "home": str(self.home),
                "runtime": str(self.runtime),
                "old_runtime": str(self.old_runtime),
                "service": self.service,
                "database": {"source": str(self.database), "backup": str(backup / "queue.sqlite3"),
                             **{k: v for k, v in before.items() if k != "columns"}},
                "columns": before["columns"],
                "files": files,
                "simqueue_symlink": os.readlink(link) if link.is_symlink() else None,
                "gpu_line_state_path": state_path,
                "active_before": before["active"],
                "claims_before": self.claim_units(),
                "layout": self.layout,
            }
            (backup / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        else:
            self.say(f"would back up {self.database}, {self.bin / 'simq'}, {self.bin / 'simqueue'}, node.json, "
                     f"worker.log and units {UNIT_FILES}")
            manifest = {"active_before": summary.get("active", []), "columns": summary.get("columns", []),
                        "claims_before": self.claim_units()}
        self.step("backup", directory=str(backup))
        log_offset = self.worker_log_size()
        receipt: dict = {"schema": "simqueue-upgrade-receipt-v1", "action": "install", "from": FROM_VERSION,
                         "to": TO_VERSION, "executed": self.execute, "backup": str(backup),
                         "active_before": manifest["active_before"]}
        try:
            # files first: from here on every simq call runs the new version
            self.write_file(self.bin / "simq", self.wrapper_text(self.runtime, TO_VERSION), 0o755)
            self.point_link(str(self.runtime / "simq"))
            self.write_file(self.units / self.service, self.unit_text(self.runtime), 0o644)
            if not self.args.keep_node_config:
                if not self.execute and self.layout is not None:
                    self.say(f"node.json from {self.args.config}: partitions {layout_text(self.layout['new'])}")
                self.write_file(self.state / "node.json", Path(self.args.config).read_text(), 0o600)
            self.step("installed-wrapper-unit-and-node-config", kept_node_config=bool(self.args.keep_node_config))
            self.say(f"migrating the database to schema {SCHEMA_NEW} with the {TO_VERSION} code")
            if self.execute:
                result = self.python_in(self.runtime, "from simqueue.store import QueueStore; QueueStore(); print('ok')")
                if result.returncode != 0:
                    raise RuntimeError(f"{TO_VERSION} cannot migrate the database: {result.stderr.strip()[-400:]}")
                after = database_summary(self.database)
                if after["user_version"] != SCHEMA_NEW:
                    raise RuntimeError(f"database schema is {after['user_version']} after the migration")
                missing = [c for c in manifest["columns"] if c not in after["columns"]]
                if missing:
                    raise RuntimeError(f"database columns lost: {missing}")
                columns = [c for c in manifest["columns"] if c != "env_json"]
                old = row_digests(backup / "queue.sqlite3", columns)
                new = row_digests(self.database, columns)
                # rows added meanwhile (a submission through the new wrapper) are not compared
                changed = [job_id for job_id, digest in old.items() if new.get(job_id) != digest]
                if changed:
                    raise RuntimeError(f"the migration changed schema-{SCHEMA_OLD} columns of jobs {changed}")
                wrong_env = self.env_scrub_problems(backup / "queue.sqlite3", self.database)
                if wrong_env:
                    raise RuntimeError(f"environments of jobs {wrong_env[:10]} are not the scrubbed ones")
                added = [c for c in after["columns"] if c not in manifest["columns"]]
                specs = self.scrub_spec_files()
                backup_scrubbed = self.scrub_backup(backup / "queue.sqlite3")
                receipt["database"] = {"schema": after["user_version"], "jobs": len(old), "rows_unchanged": True,
                                       "added_columns": added, "env_scrubbed": True,
                                       "spec_files_redacted": specs, "backup_rows_redacted": backup_scrubbed}
                print(f"database migrated: schema {SCHEMA_NEW}, {len(old)} existing rows unchanged but their "
                      f"secret-pattern environment values (redacted); {specs} finished jobs' spec.json redacted; "
                      f"{backup_scrubbed} rows of the backup copy redacted")
            self.step("database-migrated")
            self.systemctl("daemon-reload")
            self.systemctl("enable", self.service)
            self.systemctl("start", self.service)
            self.step(f"started-{TO_VERSION}-worker")
            if self.execute:
                verification = self.verify_adoption(active_before=manifest["active_before"], log_offset=log_offset,
                                                    version=TO_VERSION, runtime=self.runtime,
                                                    timeout=self.args.adoption_timeout, state_path=state_path,
                                                    claims_before=manifest.get("claims_before", []))
                doctor = self.simq(self.runtime, "doctor")
                partitions = self.simq(self.runtime, "partitions", "--sample", "0")
                verification.update({"doctor_rc": doctor.returncode, "doctor": doctor.stdout,
                                     "partitions_rc": partitions.returncode, "partitions": partitions.stdout})
                import stat as _stat
                try:
                    fifo = _stat.S_ISFIFO(os.lstat(self.state / "worker.wake").st_mode)
                except OSError:
                    fifo = False
                verification["wake_fifo"] = fifo
                if not fifo:
                    verification["ok"] = False
                if partitions.returncode != 0:
                    verification["ok"] = False
                    verification["partitions_problem"] = partitions.stderr[-400:]
                if self.args.require_doctor and doctor.returncode != 0:
                    verification["ok"] = False
                    verification["doctor_problem"] = True
                receipt["verification"] = verification
                print(doctor.stdout)
                print(partitions.stdout)
                for job_id, item in verification["jobs"].items():
                    print(f"job {job_id}: " + (f"adopted by {TO_VERSION}" if item["adopted"] else "finished meanwhile"))
                for unit, item in verification["claims"].items():
                    print(f"GPU claim {unit}: " + (f"adopted by {TO_VERSION}" if item["adopted"] else "ended meanwhile"))
                if "gpu_line_state" in verification:
                    print(f"GPU line state: {verification['gpu_line_state']}")
                if not verification["ok"]:
                    raise RuntimeError(f"verification failed: {json.dumps(verification, default=str)[:600]}")
        except Exception as error:
            receipt["error"] = str(error)
            print(f"install failed: {error}", flush=True)
            if self.execute and not self.args.no_auto_rollback:
                print("rolling back automatically", flush=True)
                (backup / "install-receipt.json").write_text(json.dumps(receipt, indent=2, default=str) + "\n")
                self.args.backup = str(backup)
                self.args.adoption_timeout = max(self.args.adoption_timeout, 60.0)
                self.problems = []
                code = self.rollback(automatic=True)
                receipt["auto_rollback_rc"] = code
            if self.execute:
                (backup / "install-receipt.json").write_text(json.dumps(receipt, indent=2, default=str) + "\n")
            return 1
        receipt["steps"] = self.steps
        if self.execute:
            (backup / "install-receipt.json").write_text(json.dumps(receipt, indent=2, default=str) + "\n")
            print(f"installed {TO_VERSION}. Rollback: {self.runtime / 'upgrade.py'} rollback --backup {backup} --execute")
        else:
            self.say("no changes made; re-run with --execute")
        return 0

    # ----- rollback -------------------------------------------------------

    def stranded_jobs(self, backup: Path, summary: dict) -> list[int]:
        """Queued jobs whose recorded partition the backed-up node.json lacks."""
        try:
            document = json.loads((backup / "node.json").read_text())
        except (OSError, ValueError):
            return []
        known = {str(p.get("name")) for p in document.get("partitions") or []}
        for part in document.get("partitions") or []:
            for alias in part.get("aliases") or []:
                known.add(str(alias if isinstance(alias, str) else alias.get("name")))
        default = str(document.get("default_partition") or "")
        stranded = []
        for name, jobs in summary.get("unfinished_by_partition", {}).items():
            if (name or default) not in known:
                stranded.extend(jobs)
        return sorted(stranded)

    def claim_units(self) -> list[str]:
        """GPU claim units recorded in the state directory (0.6.2) that have
        not finished.  The names come from this installation's own state, so
        a sandbox never sees the real ones."""
        root = self.state / "claims"
        names = []
        try:
            directories = sorted(root.iterdir())
        except OSError:
            return names
        for directory in directories:
            if (directory / "finished").exists():
                continue
            try:
                spec = json.loads((directory / "spec.json").read_text())
            except (OSError, ValueError):
                continue
            unit = str(spec.get("unit") or "")
            if re.fullmatch(r"[A-Za-z0-9_.-]+-gpuclaim-[0-9a-f]{12}\.service", unit):
                names.append(unit)
        return names

    def stop_claims(self) -> list[str]:
        names = self.claim_units()
        if not names:
            return []
        if not self.execute:
            self.say(f"would stop GPU claim unit(s) {', '.join(names)}")
            return names
        # Claim units are real transient units even in a sandbox (their names
        # come from the sandbox state); job units are never stopped here.
        stopped = []
        for name in names:
            result = subprocess.run(["systemctl", "--user", "stop", "--", name], text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False, timeout=60)
            if result.returncode == 0:
                stopped.append(name)
            else:
                print(f"could not stop {name}: {result.stderr.strip()}")
        print(f"stopped GPU claim unit(s): {', '.join(stopped) or 'none'}")
        return stopped

    def withdraw_line_state(self, path: str | None) -> None:
        if not path or not self.execute:
            return
        try:
            document = json.loads(Path(path).read_text())
        except (OSError, ValueError):
            return
        if not document.get("heartbeat_unix"):
            return
        document["heartbeat_unix"] = 0
        document["withdrawn_at"] = datetime.now(UTC).isoformat(timespec="seconds")
        temporary = Path(path + ".simq-new")
        temporary.write_text(json.dumps(document, indent=2) + "\n")
        temporary.chmod(0o644)
        os.replace(temporary, path)

    def rollback(self, *, automatic: bool = False) -> int:
        if not automatic:
            self.sandbox_guard()
        backup = Path(self.args.backup)
        manifest_path = backup / "manifest.json"
        if not manifest_path.exists():
            self.problem(f"no manifest in {backup}")
            return 1
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("from") != FROM_VERSION or manifest.get("to") != TO_VERSION:
            self.problem(f"{backup} is a {manifest.get('from')} to {manifest.get('to')} backup, "
                         f"not {FROM_VERSION} to {TO_VERSION}")
        for name, info in manifest["files"].items():
            if sha256(backup / name) != info["sha256"]:
                self.problem(f"backup file {name} does not match its recorded hash")
        old_runtime = Path(manifest.get("old_runtime") or self.old_runtime)
        if not (old_runtime / "source" / "simqueue" / "jobwrap.py").exists():
            self.problem(f"{FROM_VERSION} runtime missing: {old_runtime}")
        summary = database_summary(self.database)
        print(f"database: schema {summary['user_version']}, {summary['jobs']} jobs, states {summary['states']}")
        if summary["user_version"] not in (SCHEMA_OLD, SCHEMA_NEW):
            self.problem(f"database schema is {summary['user_version']}, expected {SCHEMA_NEW} (or {SCHEMA_OLD})")
        if self.problems:
            print(f"rollback refused: {len(self.problems)} problem(s)")
            return 1
        github_due = self.github_due()
        if github_due:
            print(f"NOTE: job(s) {github_due} still have GitHub posts due; {FROM_VERSION} does not post them "
                  f"(their statuses stay as last posted until {TO_VERSION} runs again)")
        receipt: dict = {"schema": "simqueue-upgrade-receipt-v1", "action": "rollback", "from": TO_VERSION,
                         "to": FROM_VERSION, "executed": self.execute, "automatic": automatic,
                         "backup": str(backup), "at": now_tag(), "github_due": github_due}
        state_path = manifest.get("gpu_line_state_path")
        self.say(f"stopping {self.service} (job units keep running)")
        self.systemctl("stop", self.service)
        if self.execute:
            self.wait_service(False)
        active_before = database_summary(self.database)["active"]
        receipt["active_before"] = active_before
        claims_before = self.claim_units()
        receipt["claims_before"] = claims_before
        self.say(f"GPU claim units {', '.join(claims_before) or 'none'} keep running; {FROM_VERSION} adopts them")
        self.withdraw_line_state(state_path)
        log_offset = self.worker_log_size()
        if self.execute:
            with closing(sqlite3.connect(self.database)) as current:
                current.execute("PRAGMA wal_checkpoint(TRUNCATE)")
                with closing(sqlite3.connect(backup / f"queue-{TO_VERSION}-at-rollback.sqlite3")) as keep:
                    current.backup(keep)
        else:
            self.say(f"would keep the database as {backup / f'queue-{TO_VERSION}-at-rollback.sqlite3'}")
        for name, info in manifest["files"].items():
            if name in ("worker.log",):
                continue
            self.write_file(Path(info["source"]), (backup / name).read_bytes(), int(info["mode"], 8))
        self.point_link(manifest.get("simqueue_symlink"))
        self.step("restored-files")
        if self.execute:
            with closing(sqlite3.connect(self.database)) as connection:
                version = connection.execute("PRAGMA user_version").fetchone()[0]
                if version != SCHEMA_OLD:
                    connection.execute(f"PRAGMA user_version = {SCHEMA_OLD}")
                    connection.commit()
            print(f"database schema set back to {SCHEMA_OLD} in place (was {version}); the {TO_VERSION} columns "
                  f"stay and are ignored by {FROM_VERSION}")
            receipt["schema_downgrade"] = {"from": version, "to": SCHEMA_OLD}
            check = self.python_in(
                old_runtime,
                "from simqueue.store import QueueStore, SCHEMA_VERSION; s = QueueStore();"
                "import sqlite3; c = sqlite3.connect(s.path); print(SCHEMA_VERSION, c.execute('PRAGMA user_version')"
                ".fetchone()[0], c.execute('SELECT COUNT(*) FROM jobs').fetchone()[0])",
            )
            if check.returncode != 0:
                print(f"{FROM_VERSION} cannot open the database: {check.stderr.strip()[-400:]}")
                receipt["database_error"] = check.stderr.strip()[-400:]
            else:
                code_version, user_version, jobs = (int(x) for x in check.stdout.split())
                receipt["database"] = {"schema": user_version, "code_schema": code_version, "jobs": jobs,
                                       "opened_by": str(old_runtime)}
                print(f"database schema {user_version}; opened by {old_runtime.name} ({jobs} jobs)")
        else:
            self.say(f"would set the database schema back to {SCHEMA_OLD} in place and check that {old_runtime} "
                     "opens it")
        self.systemctl("daemon-reload")
        self.systemctl("start", self.service)
        if self.execute:
            verification = self.verify_adoption(active_before=active_before, log_offset=log_offset,
                                                version=FROM_VERSION, runtime=old_runtime,
                                                timeout=self.args.adoption_timeout, claims_before=claims_before)
            if "database_error" in receipt:
                verification["ok"] = False
            receipt["verification"] = verification
            for job_id, item in verification["jobs"].items():
                print(f"job {job_id}: " + (f"adopted by {FROM_VERSION}" if item["adopted"] else "finished meanwhile"))
            (backup / ("auto-rollback-receipt.json" if automatic else "rollback-receipt.json")).write_text(
                json.dumps(receipt, indent=2, default=str) + "\n")
            if not verification["ok"]:
                print(f"rollback verification failed: {verification}")
                return 1
            print(f"rolled back to SimQueue {FROM_VERSION}")
        else:
            self.say("no changes made; re-run with --execute")
        return 0

    def status(self) -> int:
        wrapper = self.bin / "simq"
        text = wrapper.read_text() if wrapper.exists() else ""
        print(f"wrapper {wrapper} runs {wrapper_runtime(text)}")
        print(f"{self.service} ExecStart: {self.unit_property(self.service, 'ExecStart')}")
        print(f"{self.service} ActiveState: {self.unit_property(self.service, 'ActiveState')}")
        if self.database.exists():
            summary = database_summary(self.database)
            print(f"database: schema {summary['user_version']}, {summary['jobs']} jobs, states {summary['states']}")
        claims = self.claim_units()
        print(f"GPU claim units recorded as not finished: {', '.join(claims) or 'none'}")
        if self.backups.exists():
            print("backups: " + ", ".join(sorted(p.name for p in self.backups.iterdir())))
        return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=["install", "rollback", "status"])
    parser.add_argument("--execute", action="store_true", help="make the changes (default: dry run)")
    parser.add_argument("--backup", help="rollback: backup directory written by install")
    parser.add_argument("--keep-node-config", action="store_true",
                        help="install: keep the live node.json (the GPU line then uses its defaults)")
    parser.add_argument("--adoption-timeout", type=float, default=60)
    parser.add_argument("--starting-timeout", type=float, default=60,
                        help="install: how long to wait for jobs in state starting before the switch")
    parser.add_argument("--no-auto-rollback", action="store_true")
    parser.add_argument("--require-doctor", action="store_true", help="install: fail when simq doctor reports a problem")
    parser.add_argument("--allow-unsealed", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--home", default=REAL_HOME)
    parser.add_argument("--project-root", default="/home/treefor/projects")
    parser.add_argument("--runtime", default=str(RUNTIME))
    parser.add_argument("--old-runtime", default=str(OLD_RUNTIME))
    parser.add_argument("--config", default=str(RUNTIME / "config" / "zenobia-node.json"))
    parser.add_argument("--systemctl", default="systemctl")
    parser.add_argument("--service", default="simqueue.service")
    args = parser.parse_args()
    upgrade = Upgrade(args)
    if args.action == "install":
        return upgrade.install()
    if args.action == "rollback":
        if not args.backup:
            parser.error("rollback needs --backup DIR")
        return upgrade.rollback()
    return upgrade.status()


if __name__ == "__main__":
    raise SystemExit(main())
