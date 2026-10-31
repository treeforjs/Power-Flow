#!/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen
"""Switch the node from SimQueue 0.5 to 0.6, or roll back.  Dry run by default.

    switchover.py install  [--execute] [--allow-queued] [--canary PART=MODE ...] [--no-canary]
    switchover.py rollback --backup DIR [--execute]
    switchover.py status

install (with --execute):
  1. preflight: 0.6 runtime and node config valid; 0.5 queue idle (no active
     jobs, no live 0.5 job processes); database schema 6; units as expected
  2. stop simqueue.service (0.5 worker; idle, so no job is affected)
  3. back up the database (sqlite3 backup API), ~/.local/bin/simq,
     ~/.local/bin/simqueue, the simqueue unit files and worker.log
  4. migrate the database to schema 7 with the 0.6 code; verify every row
  5. install node.json, the 0.6 wrapper and the 0.6 simqueue.service
  6. daemon-reload, start simqueue.service, run `simq doctor`
  7. canaries: powerflow runs a real one-core job that reports its affinity
     and runs `nvidia-smi -L` while holding the GPU leases; spheral runs a
     no-op affinity report on cores 0-3; raikou-mfem is validated with
     `simq submit --test-only` only (no launch on cores 4-5)

rollback (with --execute) restores the backed-up database, wrapper, symlink
and unit files byte for byte and restarts the 0.5 worker.  It refuses while
0.6 job units are active.  The 0.6-era database is kept in the backup
directory as queue-0.6-at-rollback.sqlite3.

--home, --systemctl and --runtime exist so the procedure can be tested in a
sandbox; the defaults are the real paths.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import time
from contextlib import closing
from datetime import UTC, datetime
from pathlib import Path


RUNTIME = Path(__file__).resolve().parent
RUNTIME_05 = Path("/home/treefor/projects/runtimes/simqueue-0.5.0-20261029")
PYTHON = str(RUNTIME_05 / "python-frozen")
DEFAULT_CANARIES = ["powerflow=run-gpu", "spheral=run", "raikou-mfem=test-only"]
UNIT_FILES = ("simqueue.service", "simqueue-spheral-dispatch.path", "simqueue-spheral-dispatch.service")


def now_tag() -> str:
    return datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")


def sha256(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return None


class Plan:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.execute = bool(getattr(args, "execute", False))
        self.home = Path(args.home)
        self.runtime = Path(args.runtime)
        self.bin = self.home / ".local" / "bin"
        self.units = self.home / ".config" / "systemd" / "user"
        self.state = self.home / ".local" / "state" / "simqueue"
        self.backups = self.home / ".local" / "state" / "simqueue-backups"
        self.database = self.state / "queue.sqlite3"
        self.project_root = Path(args.project_root)
        self.steps: list[dict] = []
        self.problems: list[str] = []

    # ----- utilities ------------------------------------------------------

    def say(self, text: str) -> None:
        prefix = "" if self.execute else "[dry-run] "
        print(prefix + text, flush=True)

    def step(self, name: str, **detail) -> None:
        self.steps.append({"step": name, "at": now_tag(), **detail})

    def problem(self, text: str) -> None:
        self.problems.append(text)
        print("PROBLEM: " + text, flush=True)

    def systemctl(self, *arguments: str, check: bool = True) -> subprocess.CompletedProcess[str]:
        command = [self.args.systemctl, "--user", *arguments]
        if not self.execute and arguments and arguments[0] not in ("show", "is-active", "is-enabled", "list-units", "cat"):
            self.say("would run: " + " ".join(command))
            return subprocess.CompletedProcess(command, 0, "", "")
        result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        if check and result.returncode != 0:
            raise RuntimeError(f"{' '.join(command)} failed: {result.stderr.strip()}")
        return result

    def unit_property(self, unit: str, prop: str) -> str:
        result = self.systemctl("show", unit, f"--property={prop}", "--value", check=False)
        return result.stdout.strip()

    def simq_env(self) -> dict[str, str]:
        env = {k: v for k, v in os.environ.items() if not k.startswith("SIMQUEUE_")}
        env.update(
            {
                "SIMQUEUE_PROJECT_ROOT": str(self.project_root),
                "SIMQUEUE_STATE_DIR": str(self.state),
                "SIMQUEUE_CONFIG": str(self.state / "node.json"),
            }
        )
        if self.args.no_autostart:
            env["SIMQUEUE_NO_AUTOSTART"] = "1"
        return env

    def simq(self, *arguments: str, timeout: float = 120) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [PYTHON, str(self.runtime / "source" / "simq"), *arguments],
            env=self.simq_env(), cwd=str(self.project_root), text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout, check=False,
        )

    def wrapper_text(self) -> str:
        return (
            "#!/bin/sh\n"
            "# SimQueue 0.6 node scheduler: partitions from node.json; each job runs in its own\n"
            "# systemd user unit (spheral 0-3, raikou-mfem 4-5 + GPU, powerflow 6-7 + GPU).\n"
            "set -eu\n"
            "unset SIMQUEUE_PCORES\n"
            f"export SIMQUEUE_PROJECT_ROOT={self.project_root}\n"
            f"export SIMQUEUE_STATE_DIR={self.state}\n"
            f"export SIMQUEUE_CONFIG={self.state / 'node.json'}\n"
            f"exec {PYTHON} {self.runtime / 'source' / 'simq'} \"$@\"\n"
        )

    def unit_text(self) -> str:
        return (
            "[Unit]\n"
            "Description=SimQueue 0.6 node scheduler (jobs run in their own units)\n"
            f"Documentation=file://{self.runtime / 'OPERATIONS.md'}\n"
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

    def write_file(self, path: Path, text: str, mode: int) -> None:
        if not self.execute:
            self.say(f"would write {path} ({len(text)} bytes, mode {mode:o})")
            return
        temporary = path.with_name(path.name + ".simq-new")
        temporary.write_text(text, encoding="utf-8")
        temporary.chmod(mode)
        os.replace(temporary, path)

    # ----- checks ---------------------------------------------------------

    def legacy_live_processes(self) -> list[int]:
        found = []
        for entry in Path("/proc").iterdir():
            if not entry.name.isdigit() or int(entry.name) == os.getpid():
                continue
            try:
                environ = (entry / "environ").read_bytes()
            except OSError:
                continue
            if b"SIMQ_LAUNCH_TOKEN=" in environ and b"SIMQ_UNIT=" not in environ:
                found.append(int(entry.name))
        return found

    def database_summary(self, path: Path) -> dict:
        # Without a WAL file the database is fully checkpointed: read it as
        # immutable so that not even SQLite's -shm index file is created next
        # to a live database.  With a WAL file present, a normal read-only
        # connection is needed to see its content.
        wal = Path(f"{path}-wal")
        uri = f"file:{path}?mode=ro" if wal.exists() and wal.stat().st_size else f"file:{path}?immutable=1"
        with closing(sqlite3.connect(uri, uri=True)) as connection:
            version = connection.execute("PRAGMA user_version").fetchone()[0]
            states = dict(connection.execute("SELECT state, COUNT(*) FROM jobs GROUP BY state").fetchall())
            count = connection.execute("SELECT COUNT(*) FROM jobs").fetchone()[0]
            digest = hashlib.sha256()
            columns = [row[1] for row in connection.execute("PRAGMA table_info(jobs)").fetchall()]
            base_columns = [c for c in columns if c in LEGACY_COLUMNS]
            for row in connection.execute(f"SELECT {', '.join(base_columns)} FROM jobs ORDER BY id"):
                digest.update(json.dumps(row, default=str).encode())
        return {"user_version": version, "jobs": count, "states": states, "legacy_columns_sha256": digest.hexdigest()}

    def preflight_install(self) -> dict:
        print(f"SimQueue 0.6 runtime: {self.runtime}")
        sealed = self.runtime / "sealed-files.json"
        if sealed.exists():
            manifest = json.loads(sealed.read_text())
            bad = [name for name, digest in manifest["files"].items() if sha256(self.runtime / name) != digest]
            if bad:
                self.problem(f"sealed runtime files changed: {bad}")
        else:
            self.problem("sealed-files.json missing from the runtime")
        if not os.access(PYTHON, os.X_OK):
            self.problem(f"pinned interpreter not executable: {PYTHON}")
        config_source = Path(self.args.config)
        check = subprocess.run(
            [PYTHON, "-c",
             "import sys; sys.path.insert(0, sys.argv[1]); from pathlib import Path;"
             "from simqueue.nodeconfig import load_node_config; c = load_node_config(Path(sys.argv[2]));"
             "print(c.node, ' '.join(f'{p.name}={p.cpus}' for p in c.partitions.values()))",
             str(self.runtime / "source"), str(config_source)],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        )
        if check.returncode != 0:
            self.problem(f"node config {config_source} invalid: {check.stderr.strip()}")
        else:
            print(f"node config ok: {check.stdout.strip()}")
        if not self.database.exists():
            self.problem(f"database missing: {self.database}")
            return {}
        summary = self.database_summary(self.database)
        print(f"database: schema {summary['user_version']}, {summary['jobs']} jobs, states {summary['states']}")
        if summary["user_version"] != 6:
            self.problem(f"database schema is {summary['user_version']}, expected 6 (0.5)")
        active = sum(summary["states"].get(s, 0) for s in ("starting", "running", "cancelling"))
        queued = summary["states"].get("queued", 0)
        if active:
            self.problem(f"{active} 0.5 job(s) active; wait until the 0.5 queue is idle")
        if queued and not self.args.allow_queued:
            self.problem(f"{queued} job(s) queued in 0.5; they would be dispatched by 0.6 to the default "
                         "partition. Re-run with --allow-queued if that is intended")
        live = self.legacy_live_processes()
        if live:
            self.problem(f"processes from 0.5 jobs still alive: {live[:10]}")
        exec_start = self.unit_property("simqueue.service", "ExecStart")
        print(f"simqueue.service ExecStart: {exec_start[:160]}")
        wrapper = self.bin / "simq"
        wrapper_text = wrapper.read_text(errors="replace") if wrapper.exists() else ""
        if "simqueue-0.6." in exec_start or "simqueue-0.6." in wrapper_text:
            self.problem("0.6 already appears to be installed")
        return summary

    # ----- install --------------------------------------------------------

    def install(self) -> int:
        summary = self.preflight_install()
        if self.problems:
            print(f"install refused: {len(self.problems)} problem(s)")
            return 1
        tag = now_tag()
        backup = self.backups / f"0.5-to-0.6-{tag}"
        self.say(f"backup directory: {backup}")
        # 1. stop the idle 0.5 worker
        self.say("stopping simqueue.service (0.5 worker, idle)")
        self.systemctl("stop", "simqueue.service")
        if self.execute:
            deadline = time.monotonic() + 30
            while self.unit_property("simqueue.service", "ActiveState") not in ("inactive", "failed"):
                if time.monotonic() > deadline:
                    raise RuntimeError("simqueue.service did not stop")
                time.sleep(0.5)
            summary_after_stop = self.database_summary(self.database)
            if summary_after_stop["legacy_columns_sha256"] != summary["legacy_columns_sha256"]:
                raise RuntimeError("the queue changed while stopping the 0.5 worker; aborting before changes")
        self.step("stopped-0.5-worker")
        # 2. backup
        files = {}
        if self.execute:
            backup.mkdir(parents=True, mode=0o700)
            with closing(sqlite3.connect(f"file:{self.database}?mode=ro", uri=True)) as source, \
                    closing(sqlite3.connect(backup / "queue.sqlite3")) as target:
                source.backup(target)
            for name, path in (
                ("simq", self.bin / "simq"),
                ("worker.log", self.state / "worker.log"),
                *((f"units/{unit}", self.units / unit) for unit in UNIT_FILES),
            ):
                if path.exists():
                    destination = backup / name
                    destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
                    shutil.copy2(path, destination)
                    files[name] = {"source": str(path), "sha256": sha256(destination), "mode": oct(path.stat().st_mode & 0o777)}
            link = self.bin / "simqueue"
            manifest = {
                "schema": "simqueue-switchover-backup-v1",
                "created": tag,
                "home": str(self.home),
                "database": {"source": str(self.database), "backup": str(backup / "queue.sqlite3"), **summary},
                "files": files,
                "simqueue_symlink": os.readlink(link) if link.is_symlink() else None,
                "node_json_existed": (self.state / "node.json").exists(),
            }
            (backup / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
            check = self.database_summary(backup / "queue.sqlite3")
            if check["legacy_columns_sha256"] != summary["legacy_columns_sha256"]:
                raise RuntimeError("database backup does not match the live database")
        else:
            self.say(f"would back up {self.database}, {self.bin / 'simq'}, {self.bin / 'simqueue'} and units {UNIT_FILES}")
        self.step("backup", directory=str(backup))
        # 3. migrate
        self.say("migrating the database to schema 7 with the 0.6 code")
        if self.execute:
            result = subprocess.run(
                [PYTHON, "-c",
                 "import sys; sys.path.insert(0, sys.argv[1]);"
                 "from simqueue.store import QueueStore, SCHEMA_VERSION; QueueStore(); print('migrated', SCHEMA_VERSION)",
                 str(self.runtime / "source")],
                env=self.simq_env(), text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
            if result.returncode != 0:
                raise RuntimeError(f"migration failed: {result.stderr.strip()}; run rollback --backup {backup}")
            after = self.database_summary(self.database)
            expected = int(result.stdout.split()[-1])
            if after["user_version"] != expected or after["jobs"] != summary["jobs"] or \
                    after["legacy_columns_sha256"] != summary["legacy_columns_sha256"]:
                raise RuntimeError(f"migration verification failed: {after}; run rollback --backup {backup}")
            print(f"database migrated: schema {expected}, {after['jobs']} jobs, rows unchanged")
        self.step("migrated")
        # 4. install files
        self.write_file(self.state / "node.json", Path(self.args.config).read_text(), 0o600)
        self.write_file(self.bin / "simq", self.wrapper_text(), 0o755)
        if self.execute:
            link = self.bin / "simqueue"
            temporary = self.bin / "simqueue.simq-new"
            if temporary.is_symlink() or temporary.exists():
                temporary.unlink()
            temporary.symlink_to(self.runtime / "simq")
            os.replace(temporary, link)
        else:
            self.say(f"would point {self.bin / 'simqueue'} at {self.runtime / 'simq'}")
        self.write_file(self.units / "simqueue.service", self.unit_text(), 0o644)
        self.step("installed-files")
        # 5. start
        self.systemctl("daemon-reload")
        self.systemctl("enable", "simqueue.service")
        self.systemctl("start", "simqueue.service")
        self.step("started-0.6-worker")
        if self.execute:
            time.sleep(2)
            doctor = self.simq("doctor")
            print(doctor.stdout)
            if doctor.returncode != 0 and not self.args.no_autostart:
                print("simq doctor reported problems (see above)")
            partitions = self.simq("partitions")
            print(partitions.stdout)
        receipt = {"schema": "simqueue-switchover-receipt-v1", "action": "install", "executed": self.execute,
                   "backup": str(backup), "steps": self.steps}
        if self.execute and not self.args.no_canary:
            receipt["canaries"] = self.canaries()
        if self.execute:
            (backup / "install-receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
            print(f"installed. Rollback: {RUNTIME / 'switchover.py'} rollback --backup {backup} --execute")
        else:
            self.say("no changes made; re-run with --execute")
        return 0

    # ----- canaries -------------------------------------------------------

    def canaries(self) -> list[dict]:
        results = []
        specs = [] if self.args.no_canary else (self.args.canary or DEFAULT_CANARIES)
        script = (
            "echo CANARY affinity=$(grep Cpus_allowed_list /proc/self/status | cut -f2) "
            "cpu_list=$SIMQ_CPU_LIST count=$SIMQ_CPU_COUNT partition=$SIMQ_PARTITION "
            "gpus=$SIMQ_GPUS cuda=$CUDA_VISIBLE_DEVICES unit=$SIMQ_UNIT"
        )
        for spec in specs:
            partition, _, mode = spec.partition("=")
            entry: dict = {"partition": partition, "mode": mode}
            if mode == "test-only":
                result = self.simq("submit", "--test-only", "-p", partition, "--", "/bin/true")
                entry.update({"returncode": result.returncode, "plan": result.stdout.strip()})
                results.append(entry)
                print(f"canary {partition}: test-only rc={result.returncode}")
                continue
            command = ["submit", "--json", "-n", f"canary-{partition}", "-p", partition, "--mem", "256M",
                       "-t", "5:00"]
            body = script
            if mode == "run-gpu":
                command += ["-c", "1", "--gpus", "1"]
                body += "; nvidia-smi -L"
            command += ["--", "/bin/sh", "-c", body]
            result = self.simq(*command)
            if result.returncode != 0:
                entry.update({"error": result.stderr.strip()})
                results.append(entry)
                continue
            job_id = json.loads(result.stdout)["job_id"]
            entry["job_id"] = job_id
            if self.args.foreground_worker:
                self.simq("worker", "--drain", timeout=self.args.canary_timeout + 60)
            deadline = time.monotonic() + self.args.canary_timeout
            state = "queued"
            while time.monotonic() < deadline:
                shown = self.simq("show", str(job_id))
                state = next((line.split()[-1] for line in shown.stdout.splitlines() if line.startswith("State ")), "?")
                if state in ("succeeded", "failed", "canceled", "interrupted"):
                    break
                time.sleep(2)
            if state not in ("succeeded", "failed", "canceled", "interrupted"):
                self.simq("cancel", str(job_id))
                entry["note"] = "canary did not finish in time (leases held elsewhere?); canceled"
            log = self.simq("log", "-n", "60", str(job_id)).stdout
            entry.update({"state": state, "output": [l for l in log.splitlines() if l.startswith(("CANARY", "GPU "))]})
            print(f"canary {partition}: {state} {' | '.join(entry['output'])}")
            results.append(entry)
        return results

    # ----- rollback -------------------------------------------------------

    def rollback(self) -> int:
        backup = Path(self.args.backup)
        manifest_path = backup / "manifest.json"
        if not manifest_path.exists():
            self.problem(f"no manifest in {backup}")
            return 1
        manifest = json.loads(manifest_path.read_text())
        listing = self.systemctl("list-units", "simq-job-*", "--no-legend", "--plain", "--state=active", check=False)
        active_units = [line.split()[0] for line in listing.stdout.splitlines() if line.strip()]
        if active_units:
            self.problem(f"0.6 job units are active: {active_units}; wait for them or cancel them first")
        for name, info in manifest["files"].items():
            if sha256(backup / name) != info["sha256"]:
                self.problem(f"backup file {name} does not match its recorded hash")
        if self.problems:
            print(f"rollback refused: {len(self.problems)} problem(s)")
            return 1
        self.say("stopping simqueue.service (0.6 worker)")
        self.systemctl("stop", "simqueue.service")
        if self.execute:
            deadline = time.monotonic() + 30
            while self.unit_property("simqueue.service", "ActiveState") not in ("inactive", "failed"):
                if time.monotonic() > deadline:
                    raise RuntimeError("simqueue.service did not stop")
                time.sleep(0.5)
            # keep the 0.6-era history
            with closing(sqlite3.connect(self.database)) as current:
                current.execute("PRAGMA wal_checkpoint(TRUNCATE)")
                with closing(sqlite3.connect(backup / "queue-0.6-at-rollback.sqlite3")) as keep:
                    current.backup(keep)
            for suffix in ("-wal", "-shm"):
                Path(f"{self.database}{suffix}").unlink(missing_ok=True)
            temporary = self.database.with_name("queue.sqlite3.rollback")
            shutil.copy2(backup / "queue.sqlite3", temporary)
            temporary.chmod(0o600)
            os.replace(temporary, self.database)
            for name, info in manifest["files"].items():
                if name == "worker.log":
                    continue
                target = Path(info["source"])
                temporary = target.with_name(target.name + ".simq-rollback")
                shutil.copy2(backup / name, temporary)
                os.replace(temporary, target)
            link = self.bin / "simqueue"
            if manifest.get("simqueue_symlink"):
                temporary = self.bin / "simqueue.simq-rollback"
                if temporary.is_symlink() or temporary.exists():
                    temporary.unlink()
                temporary.symlink_to(manifest["simqueue_symlink"])
                os.replace(temporary, link)
            node_json = self.state / "node.json"
            if node_json.exists() and not manifest.get("node_json_existed"):
                os.replace(node_json, backup / "node.json.at-rollback")
            summary = self.database_summary(self.database)
            if summary["user_version"] != 6 or summary["legacy_columns_sha256"] != manifest["database"]["legacy_columns_sha256"]:
                raise RuntimeError(f"restored database does not match the backup: {summary}")
            print(f"database restored: schema 6, {summary['jobs']} jobs")
        else:
            self.say(f"would save the 0.6 database to {backup / 'queue-0.6-at-rollback.sqlite3'}")
            self.say(f"would restore {self.database} and {list(manifest['files'])} from {backup}")
        self.systemctl("daemon-reload")
        self.systemctl("start", "simqueue.service")
        if self.execute:
            (backup / "rollback-receipt.json").write_text(json.dumps(
                {"schema": "simqueue-switchover-receipt-v1", "action": "rollback", "at": now_tag()}, indent=2) + "\n")
            print("rolled back to SimQueue 0.5")
        else:
            self.say("no changes made; re-run with --execute")
        return 0

    def status(self) -> int:
        print(f"wrapper {self.bin / 'simq'}:")
        print((self.bin / "simq").read_text() if (self.bin / "simq").exists() else "  missing")
        print(f"simqueue.service ExecStart: {self.unit_property('simqueue.service', 'ExecStart')}")
        print(f"simqueue.service ActiveState: {self.unit_property('simqueue.service', 'ActiveState')}")
        if self.database.exists():
            print(f"database: {self.database_summary(self.database)}")
        if self.backups.exists():
            print("backups: " + ", ".join(sorted(p.name for p in self.backups.iterdir())))
        return 0


LEGACY_COLUMNS = {
    "id", "name", "state", "launch_token", "cwd", "argv_json", "env_json", "env_override_keys",
    "requested_cores", "cpu_list", "log_path", "created_at", "started_at", "finished_at", "exit_code",
    "pid", "pgid", "boot_id", "process_start_ticks", "cancel_requested_at", "message", "retry_of",
    "attempt", "idempotency_key", "request_sha256", "queue_order",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=["install", "rollback", "status"])
    parser.add_argument("--execute", action="store_true", help="make the changes (default: dry run)")
    parser.add_argument("--backup", help="rollback: backup directory written by install")
    parser.add_argument("--allow-queued", action="store_true")
    parser.add_argument("--canary", action="append", metavar="PARTITION=MODE",
                        help="run | run-gpu | test-only (default: " + ", ".join(DEFAULT_CANARIES) + ")")
    parser.add_argument("--no-canary", action="store_true")
    parser.add_argument("--canary-timeout", type=float, default=600)
    parser.add_argument("--foreground-worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--no-autostart", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--home", default="/home/treefor")
    parser.add_argument("--project-root", default="/home/treefor/projects")
    parser.add_argument("--runtime", default=str(RUNTIME))
    parser.add_argument("--config", default=str(RUNTIME / "config" / "zenobia-node.json"))
    parser.add_argument("--systemctl", default="systemctl")
    args = parser.parse_args()
    if args.home != "/home/treefor" and not args.no_autostart:
        parser.error("a sandbox --home needs --no-autostart so that nothing touches the real simqueue.service")
    plan = Plan(args)
    if args.action == "install":
        return plan.install()
    if args.action == "rollback":
        if not args.backup:
            parser.error("rollback needs --backup DIR")
        return plan.rollback()
    return plan.status()


if __name__ == "__main__":
    raise SystemExit(main())
