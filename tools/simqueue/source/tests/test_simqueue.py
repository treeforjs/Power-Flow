from __future__ import annotations

import json
import fcntl
import os
import pwd
import signal
import sqlite3
import subprocess
import sys
import tempfile
import time
import unittest
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
from pathlib import Path
from unittest.mock import patch

from simqueue.config import _capacity_pcores, format_cpu_list, parse_cpu_list
from simqueue.store import SCHEMA_VERSION, QueueStore


REPOSITORY = Path(__file__).resolve().parents[1]
def process_is_live(pid: int) -> bool:
    try:
        stat = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
    except OSError:
        return False
    fields = stat[stat.rfind(")") + 2 :].split()
    return bool(fields) and fields[0] != "Z"


class ConfigTests(unittest.TestCase):
    def test_cpu_list_round_trip(self) -> None:
        self.assertEqual(parse_cpu_list("0-3,7,9-10"), (0, 1, 2, 3, 7, 9, 10))
        self.assertEqual(format_cpu_list((0, 1, 2, 3, 7, 9, 10)), "0-3,7,9-10")

    def test_invalid_cpu_range(self) -> None:
        with self.assertRaises(ValueError):
            parse_cpu_list("7-3")

    def test_capacity_fallback_keeps_the_whole_high_capacity_cluster(self) -> None:
        capacities = {0: 1024, 1: 992, 2: 713, 3: 713}

        def capacity(path: Path) -> int:
            cpu = int(path.parent.name.removeprefix("cpu"))
            return capacities[cpu]

        with patch("simqueue.config._read_int", side_effect=capacity):
            self.assertEqual(_capacity_pcores((0, 1, 2, 3)), (0, 1))


class QueueFixture:
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.projects = self.base / "projects"
        self.projects.mkdir()
        allowed = sorted(os.sched_getaffinity(0))
        self.cpus = tuple(allowed[: min(2, len(allowed))])
        self.environment = {
            "SIMQUEUE_STATE_DIR": str(self.base / "state"),
            "SIMQUEUE_PROJECT_ROOT": str(self.projects),
            "SIMQUEUE_PCORES": format_cpu_list(self.cpus),
            "SIMQUEUE_NO_AUTOSTART": "1",
            "SIMQUEUE_CANCEL_GRACE_SECONDS": "0.5",
        }
        self.patch = patch.dict(os.environ, self.environment, clear=False)
        self.patch.start()
        self.launcher = self.base / "simq-test-launcher"
        self.launcher.write_text(
            "#!/usr/bin/env python3\n"
            "import sys\n"
            "from pathlib import Path\n"
            f"sys.path.insert(0, {str(REPOSITORY)!r})\n"
            "from simqueue import exclusive\n"
            f"exclusive.CONFIG_PATH = Path({str(self.base / 'no-exclusive-config')!r})\n"
            "from simqueue.cli import main\n"
            "raise SystemExit(main())\n",
            encoding="utf-8",
        )
        self.launcher.chmod(0o700)
        self.store = QueueStore()

    def tearDown(self) -> None:
        self.patch.stop()
        self.temporary.cleanup()

    def submit(
        self,
        argv: list[str],
        *,
        name: str = "test",
        cwd: Path | None = None,
    ) -> dict[str, object]:
        return self.store.submit(
            name=name,
            cwd=cwd or self.projects,
            argv=argv,
            environment=dict(os.environ),
            env_override_keys=[],
            requested_cores=len(self.cpus),
        )

    def run_worker(self, *, timeout: float = 15) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(self.launcher), "worker", "--drain"],
            cwd=REPOSITORY,
            env=dict(os.environ),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )

    def run_cli(
        self, *arguments: str, timeout: float = 15
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [str(self.launcher), *arguments],
            cwd=REPOSITORY,
            env=dict(os.environ),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )


class QueueTestCase(QueueFixture, unittest.TestCase):
    def test_waiting_legacy_worker_rechecks_exclusive_marker_after_lock(self) -> None:
        job = self.submit(["/bin/true"], name="must-remain-queued")
        lock_path = self.base / "state" / "worker.lock"
        with lock_path.open("a+", encoding="utf-8") as lock_handle:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
            worker = subprocess.Popen(
                [
                    sys.executable,
                    str(self.launcher),
                    "worker",
                    "--wait-for-lock",
                    "--drain",
                ],
                cwd=REPOSITORY,
                env=dict(os.environ),
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            try:
                time.sleep(0.2)
                self.assertIsNone(worker.poll())
                marker = {
                    "version": 1,
                    "user": pwd.getpwuid(os.getuid()).pw_name,
                    "uid": os.getuid(),
                    "project_root": str(self.projects),
                    "state_dir": str(self.base / "state"),
                    "pcores": format_cpu_list(self.cpus),
                    "ecores": format_cpu_list(self.cpus),
                    "runner_unit": "simqueue-runner.service",
                    "helper": "/usr/local/libexec/simqueue-shield",
                }
                (self.base / "no-exclusive-config").write_text(
                    json.dumps(marker), encoding="utf-8"
                )
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
                self.assertEqual(worker.wait(timeout=5), 2)
            finally:
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
                if worker.poll() is None:
                    worker.kill()
                    worker.wait()
        self.assertEqual(self.store.require(int(job["id"]))["state"], "queued")

    def test_fifo_store_and_pause(self) -> None:
        first = self.submit(["/bin/true"], name="first")
        second = self.submit(["/bin/true"], name="second")
        self.store.set_paused(True)
        self.assertIsNone(self.store.claim_next(self.cpus))
        self.store.set_paused(False)
        claimed = self.store.claim_next(self.cpus)
        self.assertEqual(claimed["id"], first["id"])
        self.store.finish(int(first["id"]), state="succeeded", exit_code=0)
        claimed = self.store.claim_next(self.cpus)
        self.assertEqual(claimed["id"], second["id"])

    def test_promote_and_demote_use_absolute_queued_endpoints_without_touching_active(
        self,
    ) -> None:
        active = self.submit(["/bin/true"], name="active")
        first = self.submit(["/bin/true"], name="first")
        middle = self.submit(["/bin/true"], name="middle")
        last = self.submit(["/bin/true"], name="last")
        claimed = self.store.claim_next(self.cpus)
        self.assertIsNotNone(claimed)
        assert claimed is not None
        self.assertEqual(claimed["id"], active["id"])
        active_snapshot = {
            key: self.store.require(int(active["id"]))[key]
            for key in ("state", "queue_order", "cpu_list", "started_at")
        }

        promoted = self.store.promote(int(last["id"]))
        self.assertEqual(promoted["job_id"], last["id"])
        self.assertTrue(promoted["changed"])
        self.assertEqual(promoted["direction"], "front")
        self.assertEqual(promoted["queue_position"], 1)
        self.assertEqual(promoted["queued_count"], 3)
        self.assertEqual(promoted["job"]["id"], last["id"])
        self.assertEqual(
            [job["id"] for job in self.store.queued_jobs()],
            [last["id"], first["id"], middle["id"]],
        )
        duplicate_promote = self.store.promote(int(last["id"]))
        self.assertFalse(duplicate_promote["changed"])
        self.assertEqual(duplicate_promote["queue_position"], 1)

        demoted = self.store.demote(int(last["id"]))
        self.assertTrue(demoted["changed"])
        self.assertEqual(demoted["direction"], "back")
        self.assertEqual(demoted["queue_position"], 3)
        self.assertEqual(
            [job["id"] for job in self.store.queued_jobs()],
            [first["id"], middle["id"], last["id"]],
        )
        duplicate_demote = self.store.demote(int(last["id"]))
        self.assertFalse(duplicate_demote["changed"])
        self.assertEqual(duplicate_demote["queue_position"], 3)

        self.store.promote(int(middle["id"]))
        self.store.demote(int(first["id"]))
        self.assertEqual(
            [job["id"] for job in self.store.queued_jobs()],
            [middle["id"], last["id"], first["id"]],
            "moves were relative instead of absolute front/back operations",
        )
        self.assertEqual(
            {
                key: self.store.require(int(active["id"]))[key]
                for key in ("state", "queue_order", "cpu_list", "started_at")
            },
            active_snapshot,
            "reordering queued work mutated the active job",
        )

    def test_promote_and_demote_reject_nonqueued_jobs_without_reordering(self) -> None:
        active = self.submit(["/bin/true"], name="active-reject")
        waiting = self.submit(["/bin/true"], name="waiting-preserved")
        claimed = self.store.claim_next(self.cpus)
        self.assertIsNotNone(claimed)
        assert claimed is not None
        self.assertEqual(claimed["id"], active["id"])

        for operation in (self.store.promote, self.store.demote):
            with self.subTest(operation=operation.__name__, state="starting"):
                with self.assertRaisesRegex(ValueError, "only queued jobs"):
                    operation(int(active["id"]))
        self.assertEqual(
            [job["id"] for job in self.store.queued_jobs()], [waiting["id"]]
        )

        self.store.finish(int(active["id"]), state="succeeded", exit_code=0)
        for operation in (self.store.promote, self.store.demote):
            with self.subTest(operation=operation.__name__, state="succeeded"):
                with self.assertRaisesRegex(ValueError, "only queued jobs"):
                    operation(int(active["id"]))
            with self.subTest(operation=operation.__name__, state="missing"):
                with self.assertRaisesRegex(KeyError, "no job with ID"):
                    operation(999_999_999)
        self.assertEqual(
            [job["id"] for job in self.store.queued_jobs()], [waiting["id"]]
        )

    def test_reordered_queue_persists_and_new_submissions_append_after_demoted_tail(
        self,
    ) -> None:
        first = self.submit(["/bin/true"], name="persistent-first")
        middle = self.submit(["/bin/true"], name="persistent-middle")
        last = self.submit(["/bin/true"], name="persistent-last")
        self.store.demote(int(first["id"]))

        reopened = QueueStore()
        self.assertEqual(
            [job["id"] for job in reopened.queued_jobs()],
            [middle["id"], last["id"], first["id"]],
        )

        newcomer = self.submit(["/bin/true"], name="persistent-newcomer")
        self.assertEqual(
            [job["id"] for job in reopened.queued_jobs()],
            [middle["id"], last["id"], first["id"], newcomer["id"]],
            "a normal submission did not append behind the explicitly demoted tail",
        )
        claimed = reopened.claim_next(self.cpus)
        self.assertIsNotNone(claimed)
        assert claimed is not None
        self.assertEqual(claimed["id"], middle["id"])

    def test_concurrent_reorders_leave_one_complete_dispatch_order(self) -> None:
        jobs = [
            self.submit(["/bin/true"], name=f"concurrent-order-{index}")
            for index in range(16)
        ]
        operations = [
            (self.store.promote if index % 2 == 0 else self.store.demote, job)
            for index, job in enumerate(jobs)
        ] * 3
        with ThreadPoolExecutor(max_workers=8) as pool:
            results = list(
                pool.map(
                    lambda operation: operation[0](int(operation[1]["id"])),
                    operations,
                )
            )
        self.assertEqual(len(results), len(operations))
        queued = self.store.queued_jobs()
        self.assertEqual(
            {job["id"] for job in queued},
            {job["id"] for job in jobs},
        )
        self.assertEqual(
            [job["queue_order"] for job in queued],
            list(range(1, len(jobs) + 1)),
            "concurrent transactions left duplicate or missing queue positions",
        )

    def test_reorder_cli_json_persists_while_global_queue_is_paused(self) -> None:
        first = self.submit(["/bin/true"], name="cli-first")
        middle = self.submit(["/bin/true"], name="cli-middle")
        last = self.submit(["/bin/true"], name="cli-last")
        self.store.set_paused(True)

        promote = self.run_cli("promote", str(last["id"]), "--json")
        self.assertEqual(promote.returncode, 0, promote.stderr)
        payload = json.loads(promote.stdout)
        self.assertEqual(payload["schema_version"], "simqueue-reorder-result-v1")
        self.assertEqual(payload["action"], "promote")
        self.assertTrue(payload["changed"])
        self.assertEqual(payload["job_id"], last["id"])
        self.assertEqual(payload["queue_position"], 1)
        self.assertEqual(payload["queued_count"], 3)
        self.assertEqual(payload["worker_mode"], "paused")
        self.assertIsNone(self.store.claim_next(self.cpus))

        demote = self.run_cli("demote", str(last["id"]))
        self.assertEqual(demote.returncode, 0, demote.stderr)
        self.assertIn(
            f"Demoted job {last['id']} (cli-last) to queue position 3 of 3.",
            demote.stdout,
        )
        self.assertIn("global queue remains paused", demote.stdout)
        self.assertEqual(
            [job["id"] for job in self.store.queued_jobs()],
            [first["id"], middle["id"], last["id"]],
        )

        duplicate = self.run_cli("demote", str(last["id"]), "--json")
        self.assertEqual(duplicate.returncode, 0, duplicate.stderr)
        duplicate_payload = json.loads(duplicate.stdout)
        self.assertFalse(duplicate_payload["changed"])
        self.assertEqual(duplicate_payload["queue_position"], 3)

        self.store.promote(int(last["id"]))
        self.store.set_paused(False)
        claimed = self.store.claim_next(self.cpus)
        self.assertIsNotNone(claimed)
        assert claimed is not None
        self.assertEqual(claimed["id"], last["id"])

    def test_migration_gate_blocks_mutations_and_dispatch(self) -> None:
        queued = self.submit(["/bin/true"], name="held")
        sentinel = self.base / "exclusive-installing"
        sentinel.touch()
        with patch("simqueue.store.MIGRATION_SENTINEL", sentinel):
            with self.assertRaisesRegex(RuntimeError, "installation is in progress"):
                self.submit(["/bin/true"], name="too-late")
            with self.assertRaisesRegex(RuntimeError, "installation is in progress"):
                self.store.set_paused(True)
            with self.assertRaisesRegex(RuntimeError, "installation is in progress"):
                self.store.cancel(int(queued["id"]))
            with self.assertRaisesRegex(RuntimeError, "installation is in progress"):
                self.store.promote(int(queued["id"]))
            with self.assertRaisesRegex(RuntimeError, "installation is in progress"):
                self.store.demote(int(queued["id"]))
            self.assertIsNone(self.store.claim_next(self.cpus))
        self.assertEqual(self.store.require(int(queued["id"]))["state"], "queued")

    def test_old_database_is_migrated_with_launch_tokens(self) -> None:
        legacy_path = self.base / "legacy.sqlite3"
        QueueStore(path=legacy_path)
        with closing(sqlite3.connect(legacy_path)) as connection, connection:
            connection.execute("ALTER TABLE jobs DROP COLUMN launch_token")
        migrated = QueueStore(path=legacy_path)
        job = migrated.submit(
            name="migrated",
            cwd=self.projects,
            argv=["/bin/true"],
            environment=dict(os.environ),
            env_override_keys=[],
            requested_cores=len(self.cpus),
        )
        self.assertEqual(len(str(job["launch_token"])), 32)

    def test_old_database_migration_is_concurrency_safe(self) -> None:
        legacy_path = self.base / "concurrent-legacy.sqlite3"
        QueueStore(path=legacy_path)
        with closing(sqlite3.connect(legacy_path)) as connection, connection:
            connection.execute("DROP INDEX jobs_idempotency_key_uidx")
            connection.execute("DROP INDEX jobs_dispatch_order_idx")
            connection.execute("ALTER TABLE jobs DROP COLUMN request_sha256")
            connection.execute("ALTER TABLE jobs DROP COLUMN idempotency_key")
            connection.execute("ALTER TABLE jobs DROP COLUMN launch_token")
            connection.execute("ALTER TABLE jobs DROP COLUMN queue_order")
            for column in (
                "launch_paused",
                "control_generation",
                "control_changed_at",
                "drained_generation",
                "drained_at",
                "supervisor_protocol",
                "supervisor_pid",
            ):
                connection.execute(f"ALTER TABLE ensembles DROP COLUMN {column}")
            connection.execute("PRAGMA user_version = 3")
        with ThreadPoolExecutor(max_workers=8) as pool:
            stores = list(pool.map(lambda _index: QueueStore(path=legacy_path), range(16)))
        self.assertEqual(len(stores), 16)
        with closing(sqlite3.connect(legacy_path)) as connection:
            columns = {
                row[1] for row in connection.execute("PRAGMA table_info(jobs)").fetchall()
            }
            ensemble_columns = {
                row[1]
                for row in connection.execute(
                    "PRAGMA table_info(ensembles)"
                ).fetchall()
            }
            indexes = {
                row[1] for row in connection.execute("PRAGMA index_list(jobs)").fetchall()
            }
            version = connection.execute("PRAGMA user_version").fetchone()[0]
        self.assertTrue(
            {
                "launch_token",
                "idempotency_key",
                "request_sha256",
                "queue_order",
            }.issubset(columns)
        )
        self.assertTrue(
            {
                "launch_paused",
                "control_generation",
                "control_changed_at",
                "drained_generation",
                "drained_at",
                "supervisor_protocol",
                "supervisor_pid",
            }.issubset(ensemble_columns)
        )
        self.assertIn("jobs_idempotency_key_uidx", indexes)
        self.assertIn("jobs_dispatch_order_idx", indexes)
        self.assertEqual(version, SCHEMA_VERSION)

    def test_v5_queue_order_migration_preserves_historical_fifo(self) -> None:
        first = self.submit(["/bin/true"], name="legacy-first")
        second = self.submit(["/bin/true"], name="legacy-second")
        third = self.submit(["/bin/true"], name="legacy-third")
        with closing(sqlite3.connect(self.store.path)) as connection, connection:
            connection.execute("DROP INDEX jobs_dispatch_order_idx")
            connection.execute("ALTER TABLE jobs DROP COLUMN queue_order")
            connection.execute("PRAGMA user_version = 5")

        migrated = QueueStore()
        queued = migrated.queued_jobs()
        self.assertEqual(
            [job["id"] for job in queued],
            [first["id"], second["id"], third["id"]],
        )
        self.assertEqual(
            [job["queue_order"] for job in queued],
            [first["id"], second["id"], third["id"]],
        )
        with closing(sqlite3.connect(self.store.path)) as connection:
            self.assertEqual(
                connection.execute("PRAGMA user_version").fetchone()[0], SCHEMA_VERSION
            )

    def test_idempotent_command_replay_lookup_and_collision(self) -> None:
        key = "msc080-eq:v1:campaign:reconciler:0"
        environment = dict(os.environ)
        first = self.store.submit(
            name="idempotent-command",
            cwd=self.projects,
            argv=["/bin/true"],
            environment=environment,
            env_override_keys=[],
            requested_cores=len(self.cpus),
            idempotency_key=key,
        )
        replay = self.store.submit(
            name="idempotent-command",
            cwd=self.projects,
            argv=["/bin/true"],
            environment=environment,
            env_override_keys=[],
            requested_cores=len(self.cpus),
            idempotency_key=key,
        )
        self.assertTrue(first["submission_created"])
        self.assertFalse(replay["submission_created"])
        self.assertEqual(replay["id"], first["id"])
        self.assertEqual(
            self.store.get_by_idempotency_key(key)["id"], first["id"]
        )
        self.assertRegex(str(first["request_sha256"]), r"^[0-9a-f]{64}$")
        with self.assertRaisesRegex(ValueError, "different submission request"):
            self.store.submit(
                name="idempotent-command",
                cwd=self.projects,
                argv=["/bin/false"],
                environment=environment,
                env_override_keys=[],
                requested_cores=len(self.cpus),
                idempotency_key=key,
            )
        self.assertEqual(len(self.store.list_jobs(include_all=True)), 1)

    def test_idempotent_command_concurrent_replay_creates_one_job(self) -> None:
        key = "concurrent/idempotent-command"
        environment = dict(os.environ)

        def submit(_index: int) -> dict[str, object]:
            return QueueStore().submit(
                name="concurrent-idempotent",
                cwd=self.projects,
                argv=["/bin/true"],
                environment=environment,
                env_override_keys=[],
                requested_cores=len(self.cpus),
                idempotency_key=key,
            )

        with ThreadPoolExecutor(max_workers=16) as pool:
            results = list(pool.map(submit, range(32)))
        self.assertEqual({int(item["id"]) for item in results}, {1})
        self.assertEqual(
            sum(bool(item["submission_created"]) for item in results), 1
        )
        jobs = self.store.list_jobs(include_all=True)
        self.assertEqual(len(jobs), 1)
        self.assertTrue(Path(str(jobs[0]["log_path"])).is_file())

    def test_idempotent_fingerprint_ignores_runtime_queue_environment_only(self) -> None:
        key = "runtime-environment-replay"
        first_environment = dict(os.environ)
        first_environment.update({"PWD": "/old", "SIMQ_JOB_ID": "10"})
        first = self.store.submit(
            name="runtime-normalized",
            cwd=self.projects,
            argv=["/bin/true"],
            environment=first_environment,
            env_override_keys=["PWD", "SIMQ_JOB_ID"],
            requested_cores=len(self.cpus),
            idempotency_key=key,
        )
        replay_environment = dict(first_environment)
        replay_environment.update({"PWD": "/new", "SIMQ_JOB_ID": "99"})
        replay = self.store.submit(
            name="runtime-normalized",
            cwd=self.projects,
            argv=["/bin/true"],
            environment=replay_environment,
            env_override_keys=["SIMQ_JOB_ID", "PWD"],
            requested_cores=len(self.cpus),
            idempotency_key=key,
        )
        self.assertEqual(replay["id"], first["id"])
        changed_environment = dict(replay_environment)
        changed_environment["PATH"] = "/different"
        with self.assertRaisesRegex(ValueError, "different submission request"):
            self.store.submit(
                name="runtime-normalized",
                cwd=self.projects,
                argv=["/bin/true"],
                environment=changed_environment,
                env_override_keys=["SIMQ_JOB_ID", "PWD"],
                requested_cores=len(self.cpus),
                idempotency_key=key,
            )

    def test_terminal_replay_returns_same_job_and_retry_gets_no_key(self) -> None:
        key = "terminal-command-replay"
        original = self.store.submit(
            name="terminal-idempotent",
            cwd=self.projects,
            argv=["/bin/true"],
            environment=dict(os.environ),
            env_override_keys=[],
            requested_cores=len(self.cpus),
            idempotency_key=key,
        )
        self.assertEqual(self.store.cancel(int(original["id"])), "canceled")
        replay = self.store.submit(
            name="terminal-idempotent",
            cwd=self.projects,
            argv=["/bin/true"],
            environment=dict(os.environ),
            env_override_keys=[],
            requested_cores=len(self.cpus),
            idempotency_key=key,
        )
        self.assertEqual(replay["id"], original["id"])
        self.assertEqual(replay["state"], "canceled")
        retry = self.store.retry(int(original["id"]))
        self.assertNotEqual(retry["id"], original["id"])
        self.assertIsNone(retry["idempotency_key"])
        self.assertIsNone(retry["request_sha256"])

    def test_unkeyed_identical_submissions_remain_distinct(self) -> None:
        first = self.submit(["/bin/true"], name="same-unkeyed")
        second = self.submit(["/bin/true"], name="same-unkeyed")
        self.assertNotEqual(first["id"], second["id"])
        self.assertIsNone(first["idempotency_key"])
        self.assertIsNone(second["request_sha256"])

    def test_command_cli_json_replay_returns_stable_job_id(self) -> None:
        arguments = (
            "submit",
            "-n",
            "json-idempotent",
            "-C",
            str(self.projects),
            "-c",
            str(len(self.cpus)),
            "--idempotency-key",
            "json-command-replay",
            "--json",
            "--",
            "/bin/true",
        )
        first = self.run_cli(*arguments)
        replay = self.run_cli(*arguments)
        self.assertEqual(first.returncode, 0, first.stderr)
        self.assertEqual(replay.returncode, 0, replay.stderr)
        first_payload = json.loads(first.stdout)
        replay_payload = json.loads(replay.stdout)
        self.assertTrue(first_payload["created"])
        self.assertFalse(replay_payload["created"])
        self.assertEqual(replay_payload["job_id"], first_payload["job_id"])
        self.assertEqual(first_payload["worker_mode"], "disabled")
        self.assertEqual(len(self.store.list_jobs(include_all=True)), 1)

    def test_queued_cancel_and_retry_preserve_history(self) -> None:
        original = self.submit(["/bin/true"])
        self.assertEqual(self.store.cancel(int(original["id"])), "canceled")
        retried = self.store.retry(int(original["id"]))
        self.assertEqual(retried["retry_of"], original["id"])
        self.assertEqual(retried["attempt"], 2)
        self.assertEqual(self.store.require(int(original["id"]))["state"], "canceled")

    def test_log_creation_failure_does_not_leave_a_queued_job(self) -> None:
        logs = self.base / "state" / "logs"
        logs.chmod(0o500)
        try:
            with self.assertRaises(PermissionError):
                self.submit(["/bin/true"])
        finally:
            logs.chmod(0o700)
        self.assertEqual(self.store.list_jobs(include_all=True), [])

    def test_queued_cancel_never_launches(self) -> None:
        marker = self.projects / "should-not-exist"
        job = self.submit(
            [sys.executable, "-c", "from pathlib import Path; Path('should-not-exist').touch()"]
        )
        self.store.cancel(int(job["id"]))
        result = self.run_worker()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(marker.exists())

    def test_worker_preserves_argv_logs_streams_and_affinity(self) -> None:
        code = (
            "import os,sys; "
            "print(repr(sys.argv[1:])); "
            "print('stderr-line', file=sys.stderr); "
            "print('affinity=' + ','.join(map(str, sorted(os.sched_getaffinity(0)))))"
        )
        job = self.submit(
            [sys.executable, "-c", code, "space value", "$literal", "semi;colon"]
        )
        result = self.run_worker()
        self.assertEqual(result.returncode, 0, result.stderr)
        finished = self.store.require(int(job["id"]))
        self.assertEqual(finished["state"], "succeeded")
        log = Path(str(finished["log_path"])).read_text(encoding="utf-8")
        self.assertIn("['space value', '$literal', 'semi;colon']", log)
        self.assertIn("stderr-line", log)
        self.assertIn(f"affinity={','.join(map(str, self.cpus))}", log)

    def test_failure_does_not_block_next_job(self) -> None:
        failed = self.submit([sys.executable, "-c", "raise SystemExit(7)"], name="bad")
        succeeded = self.submit(["/bin/true"], name="good")
        result = self.run_worker()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.store.require(int(failed["id"]))["state"], "failed")
        self.assertEqual(self.store.require(int(failed["id"]))["exit_code"], 7)
        self.assertEqual(self.store.require(int(succeeded["id"]))["state"], "succeeded")

    def test_reserved_ensemble_drain_exit_does_not_suspend_a_normal_command(self) -> None:
        ordinary = self.submit(
            [sys.executable, "-c", "raise SystemExit(75)"],
            name="ordinary-exit-75",
        )
        following = self.submit(["/bin/true"], name="after-ordinary-exit-75")

        result = self.run_worker()

        self.assertEqual(result.returncode, 0, result.stderr)
        finished = self.store.require(int(ordinary["id"]))
        self.assertEqual(finished["state"], "failed")
        self.assertEqual(finished["exit_code"], 75)
        self.assertIn("status 75", str(finished["message"]))
        self.assertEqual(
            self.store.require(int(following["id"]))["state"], "succeeded"
        )

    def test_cleared_environment_detached_child_cannot_overlap_next_job(self) -> None:
        child_pid_file = self.projects / "cleared-env-child.pid"
        child_done_file = self.projects / "cleared-env-child.done"
        child_code = (
            "import sys,time\n"
            "from pathlib import Path\n"
            "Path(sys.argv[1]).write_text(str(__import__('os').getpid()))\n"
            "time.sleep(0.8)\n"
            "Path(sys.argv[2]).write_text('done')\n"
        )
        parent_code = (
            "import subprocess,sys\n"
            f"subprocess.Popen([sys.executable, '-c', {child_code!r}, "
            "sys.argv[1], sys.argv[2]], env={}, start_new_session=True)\n"
        )
        detached = self.submit(
            [
                sys.executable,
                "-c",
                parent_code,
                str(child_pid_file),
                str(child_done_file),
            ],
            name="cleared-env-detached-child",
        )
        overlap_check = (
            "import sys\n"
            "from pathlib import Path\n"
            "raise SystemExit(0 if Path(sys.argv[1]).is_file() else 23)\n"
        )
        following = self.submit(
            [sys.executable, "-c", overlap_check, str(child_done_file)],
            name="after-cleared-env-child",
        )

        result = self.run_worker()

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            self.store.require(int(detached["id"]))["state"], "succeeded"
        )
        self.assertEqual(
            self.store.require(int(following["id"]))["state"], "succeeded"
        )
        self.assertTrue(child_pid_file.is_file())
        self.assertTrue(child_done_file.is_file())
        child_pid = int(child_pid_file.read_text(encoding="utf-8"))
        self.assertFalse(process_is_live(child_pid))

    def test_deleted_working_directory_fails_cleanly(self) -> None:
        transient = self.projects / "removed"
        transient.mkdir()
        job = self.submit(["/bin/true"], cwd=transient)
        transient.rmdir()
        result = self.run_worker()
        self.assertEqual(result.returncode, 0, result.stderr)
        finished = self.store.require(int(job["id"]))
        self.assertEqual(finished["state"], "failed")
        self.assertIn("working directory", finished["message"])

    def test_replaced_working_directory_cannot_escape_project_root(self) -> None:
        submitted = self.projects / "replaced"
        submitted.mkdir()
        job = self.submit(["/bin/true"], cwd=submitted)
        submitted.rmdir()
        outside = self.base / "outside"
        outside.mkdir()
        submitted.symlink_to(outside, target_is_directory=True)
        result = self.run_worker()
        self.assertEqual(result.returncode, 0, result.stderr)
        finished = self.store.require(int(job["id"]))
        self.assertEqual(finished["state"], "failed")
        self.assertIn("resolves outside", finished["message"])

    def test_running_cancel_then_continue(self) -> None:
        slow = self.submit(
            [sys.executable, "-c", "import time; print('ready', flush=True); time.sleep(30)"],
            name="slow",
        )
        following = self.submit(["/bin/true"], name="following")
        worker = subprocess.Popen(
            [sys.executable, str(self.launcher), "worker", "--drain"],
            cwd=REPOSITORY,
            env=dict(os.environ),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        deadline = time.monotonic() + 8
        while time.monotonic() < deadline:
            if self.store.require(int(slow["id"]))["state"] == "running":
                break
            time.sleep(0.05)
        else:
            worker.kill()
            self.fail("job did not reach running state")
        self.assertEqual(self.store.cancel(int(slow["id"])), "cancelling")
        self.assertEqual(worker.wait(timeout=8), 0)
        self.assertEqual(self.store.require(int(slow["id"]))["state"], "canceled")
        self.assertEqual(self.store.require(int(following["id"]))["state"], "succeeded")

    def test_worker_sigterm_interrupts_current_and_leaves_queue_intact(self) -> None:
        slow = self.submit(
            [sys.executable, "-c", "import time; time.sleep(30)"], name="current"
        )
        following = self.submit(["/bin/true"], name="queued")
        worker = subprocess.Popen(
            [sys.executable, str(self.launcher), "worker"],
            cwd=REPOSITORY,
            env=dict(os.environ),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        deadline = time.monotonic() + 8
        while time.monotonic() < deadline:
            if self.store.require(int(slow["id"]))["state"] == "running":
                break
            time.sleep(0.05)
        else:
            worker.kill()
            self.fail("job did not reach running state")
        worker.terminate()
        self.assertEqual(worker.wait(timeout=8), 0)
        self.assertEqual(
            self.store.require(int(slow["id"]))["state"], "interrupted"
        )
        self.assertEqual(self.store.require(int(following["id"]))["state"], "queued")
        drained = self.run_worker()
        self.assertEqual(drained.returncode, 0, drained.stderr)
        self.assertEqual(
            self.store.require(int(following["id"]))["state"], "succeeded"
        )

    def test_cancel_escalates_and_kills_process_group(self) -> None:
        child_pid_file = self.projects / "child.pid"
        child_code = (
            "import signal,time; "
            "signal.signal(signal.SIGTERM, signal.SIG_IGN); "
            "time.sleep(30)"
        )
        parent_code = (
            "import subprocess,sys,time; "
            f"child=subprocess.Popen([sys.executable, '-c', {child_code!r}]); "
            "open(sys.argv[1], 'w').write(str(child.pid)); "
            "print('ready', flush=True); "
            "time.sleep(30)"
        )
        job = self.submit(
            [sys.executable, "-c", parent_code, str(child_pid_file)], name="stubborn"
        )
        worker = subprocess.Popen(
            [sys.executable, str(self.launcher), "worker", "--drain"],
            cwd=REPOSITORY,
            env=dict(os.environ),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        deadline = time.monotonic() + 8
        while time.monotonic() < deadline and not child_pid_file.exists():
            time.sleep(0.05)
        if not child_pid_file.exists():
            worker.kill()
            self.fail("job did not create its child process")
        child_pid = int(child_pid_file.read_text(encoding="utf-8"))
        self.assertEqual(self.store.cancel(int(job["id"])), "cancelling")
        self.assertEqual(worker.wait(timeout=8), 0)
        finished = self.store.require(int(job["id"]))
        self.assertEqual(finished["state"], "canceled")
        log = Path(str(finished["log_path"])).read_text(encoding="utf-8")
        self.assertIn("sent SIGKILL", log)
        deadline = time.monotonic() + 3
        while time.monotonic() < deadline and Path(f"/proc/{child_pid}").exists():
            time.sleep(0.05)
        self.assertFalse(Path(f"/proc/{child_pid}").exists())

    def test_cancel_kills_nested_sessions_without_touching_unrelated_process(self) -> None:
        child_pid_file = self.projects / "nested-child.pid"
        leaf_pid_file = self.projects / "nested-leaf.pid"
        leaf_code = (
            "import signal,time\n"
            "signal.signal(signal.SIGTERM, signal.SIG_IGN)\n"
            "time.sleep(30)\n"
        )
        child_code = (
            "import os,signal,subprocess,sys,threading,time\n"
            "from pathlib import Path\n"
            "term_seen = threading.Event()\n"
            "signal.signal(signal.SIGTERM, lambda *_: term_seen.set())\n"
            "Path(sys.argv[2]).write_text(str(os.getpid()))\n"
            "term_seen.wait()\n"
            "time.sleep(0.1)\n"
            f"leaf = subprocess.Popen([sys.executable, '-c', {leaf_code!r}], "
            "start_new_session=True)\n"
            "Path(sys.argv[1]).write_text(str(leaf.pid))\n"
            "time.sleep(30)\n"
        )
        parent_code = (
            "import subprocess,sys,time\n"
            f"subprocess.Popen([sys.executable, '-c', {child_code!r}, "
            "sys.argv[1], sys.argv[2]], start_new_session=True)\n"
            "time.sleep(30)\n"
        )
        job = self.submit(
            [
                sys.executable,
                "-c",
                parent_code,
                str(leaf_pid_file),
                str(child_pid_file),
            ],
            name="nested-sessions",
        )
        following_code = (
            "import sys\n"
            "from pathlib import Path\n"
            "def live(path):\n"
            "    if not Path(path).is_file():\n"
            "        return False\n"
            "    pid = int(Path(path).read_text())\n"
            "    try:\n"
            "        stat = Path(f'/proc/{pid}/stat').read_text()\n"
            "    except OSError:\n"
            "        return False\n"
            "    return stat[stat.rfind(')') + 2:].split()[0] != 'Z'\n"
            "raise SystemExit(17 if any(live(path) for path in sys.argv[1:]) else 0)\n"
        )
        following = self.submit(
            [
                sys.executable,
                "-c",
                following_code,
                str(child_pid_file),
                str(leaf_pid_file),
            ],
            name="after-nested-cancel",
        )
        unrelated_environment = dict(os.environ)
        unrelated_environment.pop("SIMQ_JOB_ID", None)
        unrelated_environment.pop("SIMQ_LAUNCH_TOKEN", None)
        unrelated = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(30)"],
            cwd=self.projects,
            env=unrelated_environment,
            start_new_session=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        worker = subprocess.Popen(
            [sys.executable, str(self.launcher), "worker", "--drain"],
            cwd=REPOSITORY,
            env=dict(os.environ),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        try:
            deadline = time.monotonic() + 8
            while time.monotonic() < deadline and not child_pid_file.exists():
                time.sleep(0.05)
            if not child_pid_file.exists():
                self.fail("job did not create its nested session")

            self.assertEqual(self.store.cancel(int(job["id"])), "cancelling")
            worker_stdout, worker_stderr = worker.communicate(timeout=10)
            worker_code = worker.returncode
            diagnostics = worker_stdout + worker_stderr
            if worker_code:
                for path in (child_pid_file, leaf_pid_file):
                    if path.exists():
                        pid = int(path.read_text(encoding="utf-8"))
                        try:
                            diagnostics += f"\n{pid}: " + Path(
                                f"/proc/{pid}/stat"
                            ).read_text(encoding="utf-8")
                        except OSError as error:
                            diagnostics += f"\n{pid}: {error}"
            self.assertEqual(worker_code, 0, diagnostics)
            finished = self.store.require(int(job["id"]))
            self.assertEqual(finished["state"], "canceled")
            self.assertEqual(
                self.store.require(int(following["id"]))["state"], "succeeded"
            )
            self.assertTrue(
                leaf_pid_file.exists(),
                "nested child did not create a new session during TERM grace",
            )
            child_pid = int(child_pid_file.read_text(encoding="utf-8"))
            leaf_pid = int(leaf_pid_file.read_text(encoding="utf-8"))
            self.assertFalse(process_is_live(child_pid))
            self.assertFalse(process_is_live(leaf_pid))
            self.assertIsNone(
                unrelated.poll(), "unrelated process was targeted by cancellation"
            )
            log = Path(str(finished["log_path"])).read_text(encoding="utf-8")
            self.assertIn("sent SIGKILL", log)
        finally:
            if worker.poll() is None:
                worker.kill()
                worker.wait()
            for path in (leaf_pid_file, child_pid_file):
                if not path.exists():
                    continue
                pid = int(path.read_text(encoding="utf-8"))
                if not process_is_live(pid):
                    continue
                try:
                    os.killpg(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            if unrelated.poll() is None:
                unrelated.kill()
            unrelated.wait()

    def test_startup_marks_unlaunched_active_job_interrupted(self) -> None:
        stale = self.submit(["/bin/true"], name="stale")
        claimed = self.store.claim_next(self.cpus)
        self.assertEqual(claimed["id"], stale["id"])
        following = self.submit(["/bin/true"], name="following")
        result = self.run_worker()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.store.require(int(stale["id"]))["state"], "interrupted")
        self.assertEqual(self.store.require(int(following["id"]))["state"], "succeeded")

    def test_restart_recovers_launch_before_pid_was_recorded(self) -> None:
        stale = self.submit(["/bin/true"], name="popen-window")
        claimed = self.store.claim_next(self.cpus)
        self.assertEqual(claimed["id"], stale["id"])
        environment = dict(os.environ)
        environment["SIMQ_LAUNCH_TOKEN"] = str(claimed["launch_token"])
        orphan = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(30)"],
            cwd=self.projects,
            env=environment,
            start_new_session=True,
        )
        following = self.submit(["/bin/true"], name="following")
        try:
            result = self.run_worker()
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                self.store.require(int(stale["id"]))["state"], "interrupted"
            )
            self.assertEqual(
                self.store.require(int(following["id"]))["state"], "succeeded"
            )
            self.assertIsNotNone(orphan.poll())
        finally:
            if orphan.poll() is None:
                orphan.kill()
            orphan.wait()

    def test_restart_recovers_group_after_leader_exits(self) -> None:
        child_pid_file = self.projects / "orphan.pid"
        parent_code = (
            "import subprocess,sys; "
            "child=subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(30)']); "
            "open(sys.argv[1], 'w').write(str(child.pid))"
        )
        interrupted = self.submit(
            [sys.executable, "-c", parent_code, str(child_pid_file)], name="orphan-parent"
        )
        following = self.submit(["/bin/true"], name="following")
        worker = subprocess.Popen(
            [sys.executable, str(self.launcher), "worker", "--drain"],
            cwd=REPOSITORY,
            env=dict(os.environ),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        deadline = time.monotonic() + 8
        while time.monotonic() < deadline:
            current = self.store.require(int(interrupted["id"]))
            leader = current.get("pid")
            if (
                child_pid_file.exists()
                and leader is not None
                and not Path(f"/proc/{leader}").exists()
            ):
                break
            time.sleep(0.05)
        else:
            worker.kill()
            self.fail("job leader did not exit while its child remained")
        child_pid = int(child_pid_file.read_text(encoding="utf-8"))
        self.assertTrue(Path(f"/proc/{child_pid}").exists())
        worker.kill()
        worker.wait(timeout=3)

        restarted = self.run_worker()
        self.assertEqual(restarted.returncode, 0, restarted.stderr)
        self.assertEqual(
            self.store.require(int(interrupted["id"]))["state"], "interrupted"
        )
        self.assertEqual(self.store.require(int(following["id"]))["state"], "succeeded")
        deadline = time.monotonic() + 3
        while time.monotonic() < deadline and Path(f"/proc/{child_pid}").exists():
            time.sleep(0.05)
        self.assertFalse(Path(f"/proc/{child_pid}").exists())

    def test_waiting_worker_takes_over_existing_worker_lock(self) -> None:
        first = subprocess.Popen(
            [sys.executable, str(self.launcher), "worker"],
            cwd=REPOSITORY,
            env=dict(os.environ),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        second: subprocess.Popen[bytes] | None = None
        lock_path = self.base / "state" / "worker.lock"
        try:
            deadline = time.monotonic() + 5
            while time.monotonic() < deadline:
                if lock_path.exists() and lock_path.read_text().strip() == str(first.pid):
                    break
                time.sleep(0.05)
            else:
                self.fail("first worker did not acquire the lock")
            time.sleep(0.1)
            second = subprocess.Popen(
                [sys.executable, str(self.launcher), "worker", "--wait-for-lock"],
                cwd=REPOSITORY,
                env=dict(os.environ),
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            time.sleep(0.2)
            self.assertIsNone(second.poll())
            first.terminate()
            self.assertEqual(first.wait(timeout=5), 0)
            deadline = time.monotonic() + 5
            while time.monotonic() < deadline:
                if lock_path.read_text().strip() == str(second.pid):
                    break
                time.sleep(0.05)
            else:
                self.fail("waiting worker did not take over the lock")
            second.terminate()
            self.assertEqual(second.wait(timeout=5), 0)
        finally:
            if first.poll() is None:
                first.kill()
                first.wait()
            if second is not None and second.poll() is None:
                second.kill()
                second.wait()

class CliTests(QueueFixture, unittest.TestCase):
    def test_active_and_queue_shortcuts_filter_live_state(self) -> None:
        active = self.submit(["/bin/true"], name="active-one")
        waiting = self.submit(["/bin/true"], name="waiting-two")
        waiting_next = self.submit(["/bin/true"], name="waiting-three")
        claimed = self.store.claim_next(self.cpus)
        self.assertIsNotNone(claimed)
        self.assertEqual(claimed["id"], active["id"])

        active_result = subprocess.run(
            [str(self.launcher), "-active"],
            cwd=REPOSITORY,
            env=dict(os.environ),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(active_result.returncode, 0, active_result.stderr)
        self.assertIn("active-one", active_result.stdout)
        self.assertNotIn("waiting-two", active_result.stdout)

        queue_result = subprocess.run(
            [str(self.launcher), "-queue"],
            cwd=REPOSITORY,
            env=dict(os.environ),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(queue_result.returncode, 0, queue_result.stderr)
        self.assertIn("waiting-two", queue_result.stdout)
        self.assertIn("waiting-three", queue_result.stdout)
        self.assertNotIn("active-one", queue_result.stdout)
        self.assertLess(
            queue_result.stdout.index("waiting-two"),
            queue_result.stdout.index("waiting-three"),
        )
        self.assertLess(waiting["id"], waiting_next["id"])

    def test_active_and_queue_shortcuts_report_empty(self) -> None:
        for option, message in (
            ("--active", "No active jobs."),
            ("--queue", "No queued jobs."),
        ):
            with self.subTest(option=option):
                result = subprocess.run(
                    [str(self.launcher), option],
                    cwd=REPOSITORY,
                    env=dict(os.environ),
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=False,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn(message, result.stdout)

    def test_status_shortcuts_require_the_exact_singleton_token(self) -> None:
        for abbreviation in ("-a", "-act"):
            with self.subTest(abbreviation=abbreviation):
                result = subprocess.run(
                    [str(self.launcher), abbreviation],
                    cwd=REPOSITORY,
                    env=dict(os.environ),
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=False,
                )
                self.assertEqual(result.returncode, 2)

        result = subprocess.run(
            [str(self.launcher), "submit", "--", "/bin/echo", "-active"],
            cwd=self.projects,
            env=dict(os.environ),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        submitted = self.store.list_jobs(include_all=True)[0]
        self.assertEqual(json.loads(submitted["argv_json"]), ["/bin/echo", "-active"])

    def test_cli_rejects_work_outside_project_root(self) -> None:
        result = subprocess.run(
            [str(self.launcher), "submit", "-C", str(self.base), "--", "/bin/true"],
            cwd=REPOSITORY,
            env=dict(os.environ),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("must be inside", result.stderr)

    def test_shell_mode_is_explicit(self) -> None:
        result = subprocess.run(
            [
                str(self.launcher),
                "submit",
                "-C",
                str(self.projects),
                "--shell",
                "printf 'shell-ok\\n'",
            ],
            cwd=REPOSITORY,
            env=dict(os.environ),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Queued job", result.stdout)
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stderr)
        jobs = self.store.list_jobs(include_all=True)
        log = Path(str(jobs[0]["log_path"])).read_text(encoding="utf-8")
        self.assertIn("shell-ok", log)


if __name__ == "__main__":
    unittest.main()
