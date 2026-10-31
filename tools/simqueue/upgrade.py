#!/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen
"""Upgrade the node from SimQueue 0.6.0 to 0.6.1, or roll back.  Dry run by default.

    upgrade.py install  [--execute] [--keep-node-config] [--adoption-timeout S] [--no-auto-rollback]
    upgrade.py rollback --backup DIR [--execute] [--cancel-borrow-jobs]
    upgrade.py status

Job units are independent of the worker, so running jobs are never touched:
the worker is stopped, its wrapper and unit are switched to the other
runtime, and the restarted worker re-adopts every running job unit.

install (with --execute):
  1. preflight: 0.6.1 runtime sealed, its node config valid and equal to the
     live node.json apart from the 0.6.1 keys; wrapper points at 0.6.0;
     database schema 7; records the active jobs and their units
  2. stop simqueue.service (job units keep running)
  3. back up the database (sqlite3 backup API, source opened read-only), the
     wrapper, the simqueue link, the unit files, node.json and worker.log
  4. install the 0.6.1 wrapper, link, simqueue.service and simqborrow.slice
     (CPUWeight=idle) -- from here on every `simq` call runs 0.6.1
  5. migrate the database to schema 8 (new columns and the job_events
     table only) and verify that every schema-7 column of every existing row
     is unchanged; install the 0.6.1 node.json
  6. daemon-reload, start simqueue.service
  7. verify: the 0.6.1 worker started, every job unit that was active before
     is adopted (worker.log "adopted job N" and adopted_at) or has been
     recorded as finished, `simq doctor` passes, the borrow slice is idle
  On a failed verification the install rolls itself back (unless
  --no-auto-rollback).  Receipt: <backup>/install-receipt.json

rollback --backup DIR (with --execute):
  1. refuses while borrow jobs are queued or active (0.6.0 would run them as
     owner jobs); --cancel-borrow-jobs cancels them first with the 0.6.1 CLI
  2. stop simqueue.service (job units keep running)
  3. keep a copy of the schema-8 database in the backup directory
  4. restore the wrapper, link, unit file and node.json byte for byte,
     remove simqborrow.slice (if the install created it)
  5. set the database schema back to 7 in place: rows, including jobs that
     0.6.1 started and that still run, are kept; the 0.6.1 columns stay
     and are ignored by 0.6.0; verify that the 0.6.0 code opens it
  6. daemon-reload, start simqueue.service (0.6.0)
  7. verify that the 0.6.0 worker adopted every active job unit
  Receipt: <backup>/rollback-receipt.json

--home, --project-root, --systemctl, --service and --runtime exist so that
the procedure can be tested in a sandbox; the defaults are the real paths.
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
OLD_RUNTIME = Path("/home/treefor/projects/runtimes/simqueue-0.6.0-20261030")
PYTHON = "/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen"
REAL_HOME = "/home/treefor"
UNIT_FILES = ("simqueue.service", "simqueue-spheral-dispatch.path", "simqueue-spheral-dispatch.service")
BORROW_SLICE = "simqborrow.slice"
# Keys added by 0.6.1 to the node configuration; everything else must match.
NEW_CONFIG_KEYS = {"borrow", "smt_siblings"}
NEW_PARTITION_KEYS = {"lend", "lend_leased_idle", "timing_sensitive"}
NEW_LOCK_KEYS = {"borrow_path", "notes"}
SCHEMA_OLD, SCHEMA_NEW = 7, 8


def now_tag() -> str:
    return datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")


def sha256(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return None


def comparable_config(document: dict) -> dict:
    """The node configuration without the keys 0.6.1 adds (and notes)."""
    result = {k: v for k, v in document.items() if k not in NEW_CONFIG_KEYS and k != "notes"}
    if "lock_files" in result:
        result["lock_files"] = {k: v for k, v in result["lock_files"].items() if k not in NEW_LOCK_KEYS}
    result["partitions"] = [
        {k: v for k, v in part.items() if k not in NEW_PARTITION_KEYS} for part in document.get("partitions", [])
    ]
    return result


def database_summary(path: Path) -> dict:
    wal = Path(f"{path}-wal")
    uri = f"file:{path}?mode=ro" if wal.exists() and wal.stat().st_size else f"file:{path}?immutable=1"
    with closing(sqlite3.connect(uri, uri=True)) as connection:
        version = connection.execute("PRAGMA user_version").fetchone()[0]
        states = dict(connection.execute("SELECT state, COUNT(*) FROM jobs GROUP BY state").fetchall())
        count = connection.execute("SELECT COUNT(*) FROM jobs").fetchone()[0]
        columns = [row[1] for row in connection.execute("PRAGMA table_info(jobs)").fetchall()]
        borrow_open = 0
        if "borrow" in columns:
            borrow_open = connection.execute(
                "SELECT COUNT(*) FROM jobs WHERE borrow = 1 AND state IN "
                "('queued', 'starting', 'running', 'cancelling')"
            ).fetchone()[0]
        active = [
            dict(zip(("id", "name", "state", "unit_name", "launch_token"), row))
            for row in connection.execute(
                "SELECT id, name, state, " + ("unit_name" if "unit_name" in columns else "NULL")
                + ", launch_token FROM jobs WHERE state IN ('starting', 'running', 'cancelling') ORDER BY id"
            )
        ]
    return {"user_version": version, "jobs": count, "states": states, "columns": columns,
            "borrow_open": borrow_open, "active": active}


def row_digests(path: Path, columns: list[str]) -> dict[int, str]:
    """sha256 of the given columns of every job row, by id."""
    with closing(sqlite3.connect(f"file:{path}?mode=ro", uri=True)) as connection:
        rows = connection.execute(f"SELECT id, {', '.join(columns)} FROM jobs ORDER BY id").fetchall()
    return {int(row[0]): hashlib.sha256(json.dumps(row[1:], default=str).encode()).hexdigest() for row in rows}


def wrapper_runtime(text: str) -> Path | None:
    match = re.search(r"exec \S+ (\S+)/source/simq ", text)
    return Path(match.group(1)) if match else None


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
            "# borrow jobs run in simqborrow.slice at idle CPU priority.\n"
            "set -eu\n"
            "unset SIMQUEUE_PCORES\n"
            f"export SIMQUEUE_PROJECT_ROOT={self.project_root}\n"
            f"export SIMQUEUE_STATE_DIR={self.state}\n"
            f"export SIMQUEUE_CONFIG={self.state / 'node.json'}\n"
            f"exec {PYTHON} {runtime / 'source' / 'simq'} \"$@\"\n"
        )

    def unit_text(self, runtime: Path) -> str:
        return (
            "[Unit]\n"
            "Description=SimQueue 0.6.1 node scheduler (jobs run in their own units)\n"
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
            "# Job units are separate units: stopping this service never stops a job.\n"
            "KillMode=mixed\n"
            "TimeoutStopSec=30\n"
            "UMask=0077\n"
            "\n"
            "[Install]\n"
            "WantedBy=default.target\n"
        )

    @staticmethod
    def slice_text() -> str:
        return (
            "[Unit]\n"
            "Description=SimQueue borrow jobs (preemptible; idle CPU priority)\n"
            "\n"
            "[Slice]\n"
            "# cpu.idle=1: processes here get CPU time only when no other process of the\n"
            "# user manager wants that CPU (owner job units, terminal sessions, tools).\n"
            "CPUWeight=idle\n"
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
                        runtime: Path, timeout: float) -> dict:
        """Wait until the new worker has started and adopted (or recorded)
        every job that was active before the switch."""
        deadline = time.monotonic() + timeout
        result: dict = {"version": version, "jobs": {}}
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
            if started and not pending and result["service_active"]:
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
        print(f"SimQueue 0.6.1 runtime: {self.runtime}")
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
            self.problem("0.6.1 already installed")
        elif current != self.old_runtime:
            self.problem(f"wrapper runs {current}, expected {self.old_runtime}")
        if not (self.old_runtime / "source" / "simqueue" / "jobwrap.py").exists():
            self.problem(f"0.6.0 runtime missing: {self.old_runtime}")
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
            "c = load_node_config(Path(sys.argv[2]));"
            "print(c.node, ' '.join(f'{p.name}={p.cpus}:lend={p.lend},leased_idle={p.lend_leased_idle},"
            "timing={p.timing_sensitive}' for p in c.partitions.values()), c.borrow.slice)",
            str(new_config if not self.args.keep_node_config else live_config),
        )
        if check.returncode != 0:
            self.problem(f"node config invalid for 0.6.1: {check.stderr.strip()[-400:]}")
        else:
            print(f"0.6.1 node config ok: {check.stdout.strip()}")
        if not self.args.keep_node_config and live_doc and new_doc:
            if comparable_config(live_doc) != comparable_config(new_doc):
                self.problem(f"{new_config} differs from the live {live_config} beyond the 0.6.1 keys; "
                             "reconcile them or use --keep-node-config")
        if not self.database.exists():
            self.problem(f"database missing: {self.database}")
            return {}
        summary = database_summary(self.database)
        print(f"database: schema {summary['user_version']}, {summary['jobs']} jobs, states {summary['states']}")
        if summary["user_version"] != SCHEMA_OLD:
            self.problem(f"database schema is {summary['user_version']}, expected {SCHEMA_OLD} (0.6.0)")
        for job in summary["active"]:
            print(f"active job {job['id']} ({job['name']}): {job['state']}, unit {job['unit_name']}")
        exec_start = self.unit_property(self.service, "ExecStart")
        print(f"{self.service} ExecStart: {exec_start[:160]}")
        if str(self.bin / "simq") not in exec_start:
            self.problem(f"{self.service} does not run {self.bin / 'simq'}")
        if (self.units / BORROW_SLICE).exists():
            self.problem(f"{self.units / BORROW_SLICE} already exists")
        return summary

    def install(self) -> int:
        self.sandbox_guard()
        summary = self.preflight_install()
        if self.problems:
            print(f"install refused: {len(self.problems)} problem(s)")
            return 1
        tag = now_tag()
        backup = self.backups / f"0.6.0-to-0.6.1-{tag}"
        self.say(f"backup directory: {backup}")
        self.say(f"stopping {self.service} (job units keep running)")
        self.systemctl("stop", self.service)
        if self.execute:
            self.wait_service(False)
        self.step("stopped-0.6.0-worker")
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
                "from": "0.6.0", "to": "0.6.1",
                "created": tag,
                "home": str(self.home),
                "runtime": str(self.runtime),
                "old_runtime": str(self.old_runtime),
                "service": self.service,
                "database": {"source": str(self.database), "backup": str(backup / "queue.sqlite3"),
                             **{k: v for k, v in before.items() if k != "columns"}},
                "schema7_columns": before["columns"],
                "files": files,
                "simqueue_symlink": os.readlink(link) if link.is_symlink() else None,
                "borrow_slice_existed": (self.units / BORROW_SLICE).exists(),
                "active_before": before["active"],
            }
            (backup / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        else:
            self.say(f"would back up {self.database}, {self.bin / 'simq'}, {self.bin / 'simqueue'}, node.json, "
                     f"worker.log and units {UNIT_FILES}")
            manifest = {"active_before": summary.get("active", []), "schema7_columns": summary.get("columns", [])}
        self.step("backup", directory=str(backup))
        log_offset = self.worker_log_size()
        receipt: dict = {"schema": "simqueue-upgrade-receipt-v1", "action": "install", "executed": self.execute,
                         "backup": str(backup), "active_before": manifest["active_before"]}
        try:
            # files first: from here on every simq call runs 0.6.1 (which
            # accepts the 0.6.0 node.json and migrates the database itself)
            self.write_file(self.bin / "simq", self.wrapper_text(self.runtime, "0.6.1"), 0o755)
            self.point_link(str(self.runtime / "simq"))
            self.write_file(self.units / self.service, self.unit_text(self.runtime), 0o644)
            self.write_file(self.units / BORROW_SLICE, self.slice_text(), 0o644)
            self.step("installed-wrapper-and-units")
            self.say("migrating the database to schema 8 with the 0.6.1 code")
            if self.execute:
                result = self.python_in(self.runtime, "from simqueue.store import QueueStore; QueueStore(); print('ok')")
                if result.returncode != 0:
                    raise RuntimeError(f"migration failed: {result.stderr.strip()[-400:]}")
                after = database_summary(self.database)
                if after["user_version"] != SCHEMA_NEW:
                    raise RuntimeError(f"database schema is {after['user_version']} after migration")
                old = row_digests(backup / "queue.sqlite3", manifest["schema7_columns"])
                new = row_digests(self.database, manifest["schema7_columns"])
                changed = [job_id for job_id, digest in old.items() if new.get(job_id) != digest]
                if changed:
                    raise RuntimeError(f"migration changed schema-7 columns of jobs {changed}")
                receipt["migration"] = {"schema": after["user_version"], "jobs_before": len(old),
                                        "jobs_after": after["jobs"], "rows_unchanged": True,
                                        "new_columns": [c for c in after["columns"]
                                                        if c not in manifest["schema7_columns"]]}
                print(f"database migrated: schema {SCHEMA_NEW}, {len(old)} existing rows unchanged")
            self.step("migrated")
            if not self.args.keep_node_config:
                self.write_file(self.state / "node.json", Path(self.args.config).read_text(), 0o600)
            self.step("installed-node-config", kept=bool(self.args.keep_node_config))
            self.systemctl("daemon-reload")
            self.systemctl("enable", self.service)
            self.systemctl("start", self.service)
            self.step("started-0.6.1-worker")
            if self.execute:
                verification = self.verify_adoption(active_before=manifest["active_before"], log_offset=log_offset,
                                                    version="0.6.1", runtime=self.runtime,
                                                    timeout=self.args.adoption_timeout)
                doctor = self.simq(self.runtime, "doctor")
                partitions = self.simq(self.runtime, "partitions", "--sample", "0")
                slice_weight = self.unit_property(BORROW_SLICE, "CPUWeight")
                verification.update({"doctor_rc": doctor.returncode, "doctor": doctor.stdout,
                                     "partitions": partitions.stdout, "borrow_slice_cpu_weight": slice_weight})
                if self.args.require_doctor and doctor.returncode != 0:
                    verification["ok"] = False
                    verification["doctor_problem"] = True
                receipt["verification"] = verification
                print(doctor.stdout)
                print(partitions.stdout)
                for job_id, item in verification["jobs"].items():
                    print(f"job {job_id}: " + ("adopted by 0.6.1" if item["adopted"] else "finished meanwhile"))
                if not verification["ok"]:
                    raise RuntimeError(f"verification failed: {json.dumps(verification, default=str)[:600]}")
        except Exception as error:
            receipt["error"] = str(error)
            print(f"install failed: {error}", flush=True)
            if self.execute and not self.args.no_auto_rollback:
                print("rolling back automatically", flush=True)
                (backup / "install-receipt.json").write_text(json.dumps(receipt, indent=2, default=str) + "\n")
                self.args.backup = str(backup)
                self.args.cancel_borrow_jobs = True
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
            print(f"installed 0.6.1. Rollback: {RUNTIME / 'upgrade.py'} rollback --backup {backup} --execute")
        else:
            self.say("no changes made; re-run with --execute")
        return 0

    # ----- rollback -------------------------------------------------------

    def rollback(self, *, automatic: bool = False) -> int:
        if not automatic:
            self.sandbox_guard()
        backup = Path(self.args.backup)
        manifest_path = backup / "manifest.json"
        if not manifest_path.exists():
            self.problem(f"no manifest in {backup}")
            return 1
        manifest = json.loads(manifest_path.read_text())
        for name, info in manifest["files"].items():
            if sha256(backup / name) != info["sha256"]:
                self.problem(f"backup file {name} does not match its recorded hash")
        old_runtime = Path(manifest.get("old_runtime") or self.old_runtime)
        summary = database_summary(self.database)
        print(f"database: schema {summary['user_version']}, {summary['jobs']} jobs, states {summary['states']}, "
              f"{summary['borrow_open']} borrow job(s) queued or active")
        if summary["borrow_open"] and not self.args.cancel_borrow_jobs:
            self.problem(f"{summary['borrow_open']} borrow job(s) queued or active; 0.6.0 would run them as owner "
                         "jobs. Cancel them (simq cancel ID) or re-run with --cancel-borrow-jobs")
        if self.problems:
            print(f"rollback refused: {len(self.problems)} problem(s)")
            return 1
        receipt: dict = {"schema": "simqueue-upgrade-receipt-v1", "action": "rollback", "executed": self.execute,
                         "automatic": automatic, "backup": str(backup), "at": now_tag()}
        if summary["borrow_open"]:
            self.say(f"canceling {summary['borrow_open']} borrow job(s) with the 0.6.1 CLI")
            if self.execute:
                receipt["canceled_borrow_jobs"] = self.cancel_borrow_jobs()
        active_before = database_summary(self.database)["active"]
        self.say(f"stopping {self.service} (job units keep running)")
        self.systemctl("stop", self.service)
        if self.execute:
            self.wait_service(False)
        log_offset = self.worker_log_size()
        if self.execute:
            with closing(sqlite3.connect(self.database)) as current:
                current.execute("PRAGMA wal_checkpoint(TRUNCATE)")
                with closing(sqlite3.connect(backup / "queue-0.6.1-at-rollback.sqlite3")) as keep:
                    current.backup(keep)
        else:
            self.say(f"would keep the schema-8 database as {backup / 'queue-0.6.1-at-rollback.sqlite3'}")
        # restore the 0.6.0 files before the database goes back to schema 7:
        # a 0.6.1 command run meanwhile would migrate it again
        for name, info in manifest["files"].items():
            if name in ("worker.log",):
                continue
            self.write_file(Path(info["source"]), (backup / name).read_bytes(), int(info["mode"], 8))
        self.point_link(manifest.get("simqueue_symlink"))
        if not manifest.get("borrow_slice_existed"):
            slice_path = self.units / BORROW_SLICE
            if self.execute:
                slice_path.unlink(missing_ok=True)
            else:
                self.say(f"would remove {slice_path}")
        self.step("restored-files")
        if self.execute:
            receipt["downgrade"] = self.downgrade_database(old_runtime)
        else:
            self.say(f"would set the database schema back to {SCHEMA_OLD} and open it with {old_runtime}")
        self.systemctl("daemon-reload")
        self.systemctl("start", self.service)
        if self.execute:
            verification = self.verify_adoption(active_before=active_before, log_offset=log_offset,
                                                version="0.6.0", runtime=old_runtime,
                                                timeout=self.args.adoption_timeout)
            receipt["verification"] = verification
            for job_id, item in verification["jobs"].items():
                print(f"job {job_id}: " + ("adopted by 0.6.0" if item["adopted"] else "finished meanwhile"))
            (backup / ("auto-rollback-receipt.json" if automatic else "rollback-receipt.json")).write_text(
                json.dumps(receipt, indent=2, default=str) + "\n")
            if not verification["ok"]:
                print(f"rollback verification failed: {verification}")
                return 1
            print("rolled back to SimQueue 0.6.0")
        else:
            self.say("no changes made; re-run with --execute")
        return 0

    def cancel_borrow_jobs(self) -> list[int]:
        with closing(sqlite3.connect(f"file:{self.database}?mode=ro", uri=True)) as connection:
            ids = [row[0] for row in connection.execute(
                "SELECT id FROM jobs WHERE borrow = 1 AND state IN ('queued', 'starting', 'running', 'cancelling')")]
        for job_id in ids:
            result = self.simq(self.runtime, "cancel", str(job_id))
            print(result.stdout.strip() or result.stderr.strip())
        deadline = time.monotonic() + 120
        while time.monotonic() < deadline:
            if not database_summary(self.database)["borrow_open"]:
                return ids
            if self.unit_property(self.service, "ActiveState") != "active":
                # no worker to record the end: run one pass of the 0.6.1 worker
                self.simq(self.runtime, "worker", "--once", timeout=120)
            time.sleep(1)
        raise RuntimeError("borrow jobs did not end within 120 s; cancel them and retry")

    def downgrade_database(self, old_runtime: Path) -> dict:
        with closing(sqlite3.connect(self.database, timeout=10)) as connection:
            connection.execute("BEGIN IMMEDIATE")
            columns = {row[1] for row in connection.execute("PRAGMA table_info(jobs)")}
            if "borrow" in columns:
                open_borrow = connection.execute(
                    "SELECT COUNT(*) FROM jobs WHERE borrow = 1 AND state IN "
                    "('queued', 'starting', 'running', 'cancelling')").fetchone()[0]
                if open_borrow:
                    connection.rollback()
                    raise RuntimeError(f"{open_borrow} borrow job(s) appeared during the rollback")
            connection.execute(f"PRAGMA user_version = {SCHEMA_OLD}")
            connection.commit()
        check = self.python_in(
            old_runtime,
            "from simqueue.store import QueueStore, SCHEMA_VERSION; s = QueueStore();"
            "import sqlite3; c = sqlite3.connect(s.path); print(SCHEMA_VERSION, c.execute('PRAGMA user_version')"
            ".fetchone()[0], c.execute('SELECT COUNT(*) FROM jobs').fetchone()[0])",
        )
        if check.returncode != 0:
            raise RuntimeError(f"0.6.0 cannot open the downgraded database: {check.stderr.strip()[-400:]}")
        code_version, user_version, jobs = (int(x) for x in check.stdout.split())
        if user_version != SCHEMA_OLD:
            raise RuntimeError(f"database schema is {user_version} after the downgrade")
        print(f"database schema set back to {SCHEMA_OLD}; opened by {old_runtime.name} ({jobs} jobs)")
        return {"schema": user_version, "jobs": jobs, "opened_by": str(old_runtime)}

    def status(self) -> int:
        wrapper = self.bin / "simq"
        text = wrapper.read_text() if wrapper.exists() else ""
        print(f"wrapper {wrapper} runs {wrapper_runtime(text)}")
        print(f"{self.service} ExecStart: {self.unit_property(self.service, 'ExecStart')}")
        print(f"{self.service} ActiveState: {self.unit_property(self.service, 'ActiveState')}")
        print(f"{BORROW_SLICE}: {'installed' if (self.units / BORROW_SLICE).exists() else 'absent'}")
        if self.database.exists():
            summary = database_summary(self.database)
            print(f"database: schema {summary['user_version']}, {summary['jobs']} jobs, states {summary['states']}")
        if self.backups.exists():
            print("backups: " + ", ".join(sorted(p.name for p in self.backups.iterdir())))
        return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=["install", "rollback", "status"])
    parser.add_argument("--execute", action="store_true", help="make the changes (default: dry run)")
    parser.add_argument("--backup", help="rollback: backup directory written by install")
    parser.add_argument("--cancel-borrow-jobs", action="store_true",
                        help="rollback: cancel queued and active borrow jobs first")
    parser.add_argument("--keep-node-config", action="store_true",
                        help="install: keep the live node.json (0.6.1 defaults apply to its partitions)")
    parser.add_argument("--adoption-timeout", type=float, default=60)
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
