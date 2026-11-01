#!/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen
"""Upgrade the node from SimQueue 0.6.1 to 0.6.2, or roll back.  Dry run by default.

    upgrade.py install  [--execute] [--keep-node-config] [--adoption-timeout S] [--starting-timeout S]
                        [--no-auto-rollback] [--require-doctor]
    upgrade.py rollback --backup DIR [--execute]
    upgrade.py status

Job units are independent of the worker, so running jobs are never touched:
the worker is stopped, its wrapper and unit are switched to the other
runtime, and the restarted worker re-adopts every running job unit.  0.6.2
keeps the 0.6.1 database (schema 8) as it is; nothing is migrated.

install (with --execute):
  1. preflight: 0.6.2 runtime sealed, its node config valid and equal to the
     live node.json apart from the 0.6.2 key (gpu_line) and notes; wrapper
     points at 0.6.1; database schema 8; records the active jobs and units
  2. waits until no job is in state "starting" (--starting-timeout), stops
     simqueue.service (job units keep running) and checks again: a job that
     became "starting" meanwhile makes the install start the 0.6.1 worker
     again and stop without changing anything
  3. backs up the database (sqlite3 backup API, source opened read-only),
     the wrapper, the simqueue link, the unit files (with simqborrow.slice),
     node.json and worker.log, with sha256 sums in manifest.json
  4. installs the 0.6.2 wrapper, link, simqueue.service and node.json --
     from here on every `simq` call runs 0.6.2
  5. opens the database with the 0.6.2 code and checks that its schema and
     every column of every row are unchanged
  6. daemon-reload, enable, start simqueue.service
  7. verifies: the 0.6.2 worker started, every job unit that was active
     before is adopted (worker.log "adopted job N") or has been recorded as
     finished, the GPU line state file is published by 0.6.2, `simq doctor`
     and `simq partitions` run
  On a failed verification the install rolls itself back (unless
  --no-auto-rollback).  Receipt: <backup>/install-receipt.json

rollback --backup DIR (with --execute):
  1. checks the backup against its manifest
  2. stops simqueue.service (job units keep running; the 0.6.2 worker
     withdraws the published GPU line state)
  3. stops the GPU claim units of 0.6.2 (they wait for or hold the GPU lease
     files for queued jobs, and 0.6.1 does not hand them over); job units
     are never stopped
  4. keeps a copy of the database in the backup directory
  5. restores the wrapper, link, unit files and node.json byte for byte
  6. checks that the 0.6.1 code opens the database (schema 8, unchanged)
  7. daemon-reload, start simqueue.service (0.6.1) and verifies that it
     adopted every active job unit
  Receipt: <backup>/rollback-receipt.json

Job units started by 0.6.2 run the 0.6.2 jobwrap.py: keep the 0.6.2 runtime
directory after a rollback until they have ended (and the 0.6.1 directory
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
OLD_RUNTIME = Path("/home/treefor/projects/runtimes/simqueue-0.6.1-20261031")
FROM_VERSION, TO_VERSION = "0.6.1", "0.6.2"
PYTHON = "/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen"
REAL_HOME = "/home/treefor"
UNIT_FILES = ("simqueue.service", "simqueue-spheral-dispatch.path", "simqueue-spheral-dispatch.service",
              "simqborrow.slice")
# Keys added by 0.6.2 to the node configuration; everything else must match.
NEW_CONFIG_KEYS = {"gpu_line"}
NEW_LOCK_KEYS = {"notes"}
SCHEMA = 8  # unchanged by 0.6.2
ACTIVE_STATES = ("starting", "running", "cancelling")


def now_tag() -> str:
    return datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")


def sha256(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return None


def comparable_config(document: dict) -> dict:
    """The node configuration without the key 0.6.2 adds (and notes)."""
    result = {k: v for k, v in document.items() if k not in NEW_CONFIG_KEYS and k != "notes"}
    if "lock_files" in result:
        result["lock_files"] = {k: v for k, v in result["lock_files"].items() if k not in NEW_LOCK_KEYS}
    return result


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
    return {"user_version": version, "jobs": count, "states": states, "columns": columns, "active": active,
            "starting": [job["id"] for job in active if job["state"] == "starting"]}


def row_digests(path: Path, columns: list[str]) -> dict[int, str]:
    """sha256 of the given columns of every job row, by id."""
    with closing(sqlite3.connect(f"file:{path}?mode=ro", uri=True)) as connection:
        rows = connection.execute(f"SELECT id, {', '.join(columns)} FROM jobs ORDER BY id").fetchall()
    return {int(row[0]): hashlib.sha256(json.dumps(row[1:], default=str).encode()).hexdigest() for row in rows}


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

    def python_in(self, runtime: Path, code: str, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [PYTHON, "-c", "import sys; sys.path.insert(0, sys.argv[1]);" + code, str(runtime / "source"), *arguments],
            env=self.env_for(), text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False, timeout=120,
        )

    def simq(self, runtime: Path, *arguments: str, timeout: float = 120) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [PYTHON, str(runtime / "source" / "simq"), *arguments], env=self.env_for(), cwd=str(self.project_root),
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout, check=False,
        )

    def wrapper_text(self, runtime: Path, version: str) -> str:
        return (
            "#!/bin/sh\n"
            f"# SimQueue {version} node scheduler: partitions from node.json; each job runs in its own\n"
            "# systemd user unit (spheral 0-3, raikou-mfem 4-5 + GPU, powerflow 6-7 + GPU);\n"
            "# borrow jobs run in simqborrow.slice at idle CPU priority; queued GPU jobs\n"
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
            "# only polls and must not be confined to one partition.\n"
            "Nice=5\n"
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

    def verify_adoption(self, *, active_before: list[dict], log_offset: int, version: str,
                        runtime: Path, timeout: float, state_path: str | None = None) -> dict:
        """Wait until the new worker has started and adopted (or recorded)
        every job that was active before the switch; with *state_path*, also
        until that worker has published the GPU line state."""
        deadline = time.monotonic() + timeout
        result: dict = {"version": version, "jobs": {}}
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
            "c = load_node_config(Path(sys.argv[2])); l = c.line;"
            "print(c.node, ' '.join(f'{p.name}={p.cpus}' for p in c.partitions.values()),"
            "f'gpu_line enabled={l.enabled} state={l.state_path} tickets={l.ticket_dir}"
            " gpus={list(l.gpus)} requeue={l.requeue_raw_waiters}')",
            str(new_config if not self.args.keep_node_config else live_config),
        )
        if check.returncode != 0:
            self.problem(f"node config invalid for {TO_VERSION}: {check.stderr.strip()[-400:]}")
        else:
            print(f"{TO_VERSION} node config ok: {check.stdout.strip()}")
        if not self.args.keep_node_config and live_doc and new_doc:
            if comparable_config(live_doc) != comparable_config(new_doc):
                self.problem(f"{new_config} differs from the live {live_config} beyond the {TO_VERSION} keys; "
                             "reconcile them or use --keep-node-config")
        if not self.database.exists():
            self.problem(f"database missing: {self.database}")
            return {}
        summary = database_summary(self.database)
        print(f"database: schema {summary['user_version']}, {summary['jobs']} jobs, states {summary['states']}")
        if summary["user_version"] != SCHEMA:
            self.problem(f"database schema is {summary['user_version']}, expected {SCHEMA} ({FROM_VERSION})")
        for job in summary["active"]:
            print(f"active job {job['id']} ({job['name']}): {job['state']}, unit {job['unit_name']}")
        exec_start = self.unit_property(self.service, "ExecStart")
        print(f"{self.service} ExecStart: {exec_start[:160]}")
        if str(self.bin / "simq") not in exec_start:
            self.problem(f"{self.service} does not run {self.bin / 'simq'}")
        return summary

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
            }
            (backup / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        else:
            self.say(f"would back up {self.database}, {self.bin / 'simq'}, {self.bin / 'simqueue'}, node.json, "
                     f"worker.log and units {UNIT_FILES}")
            manifest = {"active_before": summary.get("active", []), "columns": summary.get("columns", [])}
        self.step("backup", directory=str(backup))
        log_offset = self.worker_log_size()
        receipt: dict = {"schema": "simqueue-upgrade-receipt-v1", "action": "install", "from": FROM_VERSION,
                         "to": TO_VERSION, "executed": self.execute, "backup": str(backup),
                         "active_before": manifest["active_before"]}
        try:
            # files first: from here on every simq call runs 0.6.2
            self.write_file(self.bin / "simq", self.wrapper_text(self.runtime, TO_VERSION), 0o755)
            self.point_link(str(self.runtime / "simq"))
            self.write_file(self.units / self.service, self.unit_text(self.runtime), 0o644)
            if not self.args.keep_node_config:
                self.write_file(self.state / "node.json", Path(self.args.config).read_text(), 0o600)
            self.step("installed-wrapper-unit-and-node-config", kept_node_config=bool(self.args.keep_node_config))
            self.say(f"opening the database with the {TO_VERSION} code (no migration: schema {SCHEMA} stays)")
            if self.execute:
                result = self.python_in(self.runtime, "from simqueue.store import QueueStore; QueueStore(); print('ok')")
                if result.returncode != 0:
                    raise RuntimeError(f"{TO_VERSION} cannot open the database: {result.stderr.strip()[-400:]}")
                after = database_summary(self.database)
                if after["user_version"] != SCHEMA:
                    raise RuntimeError(f"database schema is {after['user_version']} after opening it")
                if after["columns"] != manifest["columns"]:
                    raise RuntimeError(f"database columns changed: {after['columns']}")
                old = row_digests(backup / "queue.sqlite3", manifest["columns"])
                new = row_digests(self.database, manifest["columns"])
                changed = [job_id for job_id, digest in old.items() if new.get(job_id) != digest]
                if changed:
                    raise RuntimeError(f"opening the database changed jobs {changed}")
                receipt["database"] = {"schema": after["user_version"], "jobs": len(old), "rows_unchanged": True}
                print(f"database: schema {SCHEMA}, {len(old)} rows unchanged")
            self.step("database-checked")
            self.systemctl("daemon-reload")
            self.systemctl("enable", self.service)
            self.systemctl("start", self.service)
            self.step(f"started-{TO_VERSION}-worker")
            if self.execute:
                verification = self.verify_adoption(active_before=manifest["active_before"], log_offset=log_offset,
                                                    version=TO_VERSION, runtime=self.runtime,
                                                    timeout=self.args.adoption_timeout, state_path=state_path)
                doctor = self.simq(self.runtime, "doctor")
                partitions = self.simq(self.runtime, "partitions", "--sample", "0")
                verification.update({"doctor_rc": doctor.returncode, "doctor": doctor.stdout,
                                     "partitions_rc": partitions.returncode, "partitions": partitions.stdout})
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
        if summary["user_version"] != SCHEMA:
            self.problem(f"database schema is {summary['user_version']}, expected {SCHEMA}")
        if self.problems:
            print(f"rollback refused: {len(self.problems)} problem(s)")
            return 1
        receipt: dict = {"schema": "simqueue-upgrade-receipt-v1", "action": "rollback", "from": TO_VERSION,
                         "to": FROM_VERSION, "executed": self.execute, "automatic": automatic,
                         "backup": str(backup), "at": now_tag()}
        state_path = manifest.get("gpu_line_state_path")
        self.say(f"stopping {self.service} (job units keep running)")
        self.systemctl("stop", self.service)
        if self.execute:
            self.wait_service(False)
        active_before = database_summary(self.database)["active"]
        receipt["active_before"] = active_before
        receipt["stopped_claim_units"] = self.stop_claims()
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
            self.say(f"would check that {old_runtime} opens the database (schema {SCHEMA})")
        self.systemctl("daemon-reload")
        self.systemctl("start", self.service)
        if self.execute:
            verification = self.verify_adoption(active_before=active_before, log_offset=log_offset,
                                                version=FROM_VERSION, runtime=old_runtime,
                                                timeout=self.args.adoption_timeout)
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
