"""Schema 6 (SimQueue 0.5) to schema 7 (0.6) migration.

The 0.5 database is created by the real, read-only 0.5 runtime code.  An
optional check runs against a copy of a production database when
SIMQ_LIVE_DB_COPY names one (never the live file itself).
"""

from __future__ import annotations

import json
import os
import shutil
import sqlite3
import subprocess
import tempfile
import unittest
from contextlib import closing
from pathlib import Path

from simqueue.store import NODE_JOB_COLUMNS, SCHEMA_VERSION, QueueStore


LEGACY_SOURCE = "/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/source"
PYTHON_FROZEN = "/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen"
NEW_COLUMNS = {name for name, _ in NODE_JOB_COLUMNS}

LEGACY_SETUP = r'''
import json, os, sys
sys.path.insert(0, %(source)r)
from pathlib import Path
from simqueue.store import QueueStore
store = QueueStore()
projects = Path(os.environ["SIMQUEUE_PROJECT_ROOT"])
env = {"PATH": "/usr/bin", "A": "1"}
done = store.submit(name="done", cwd=projects, argv=["/bin/true"], environment=env,
                    env_override_keys=[], requested_cores=4, idempotency_key="campaign:0")
store.claim_next((0, 1, 2, 3))
store.mark_running(done["id"], pid=1, pgid=1, boot_id="b", process_start_ticks=1)
store.finish(done["id"], state="succeeded", exit_code=0)
running = store.submit(name="running", cwd=projects, argv=["/bin/sleep", "9"], environment=env,
                       env_override_keys=["A"], requested_cores=4)
store.claim_next((0, 1, 2, 3))
store.mark_running(running["id"], pid=2, pgid=2, boot_id="b", process_start_ticks=2)
queued = store.submit(name="queued", cwd=projects, argv=["/bin/echo", "x"], environment=env,
                      env_override_keys=[], requested_cores=2)
later = store.submit(name="later", cwd=projects, argv=["/bin/echo", "y"], environment=env,
                     env_override_keys=[], requested_cores=1)
store.promote(later["id"])
ensemble = store.submit_ensemble(
    name="ens", cwd=projects, environment=env, env_override_keys=[], requested_cores=2,
    runner_prefix=["/bin/false", "_run-ensemble"], schema_version="simqueue-ensemble-v1",
    canonical_spec={"name": "ens"}, members=[
        {"name": "m0", "cwd": str(projects), "argv": ["/bin/true"], "env": {}, "gpu_admission_bytes": 0},
        {"name": "m1", "cwd": str(projects), "argv": ["/bin/true"], "env": {}, "gpu_admission_bytes": 0},
    ], max_parallel=2, gpu_device="0", gpu_budget_bytes=0, diagnostics_interval=1.0,
    mps_mode="off", mps_pipe_directory=None, fail_fast=False, allow_shared_cwd=True)
print(json.dumps({"ids": [done["id"], running["id"], queued["id"], later["id"], ensemble["id"]]}))
'''

LEGACY_OPEN = r'''
import sys
sys.path.insert(0, %(source)r)
from simqueue.store import QueueStore
try:
    QueueStore()
except RuntimeError as error:
    print("REFUSED", error)
else:
    print("OPENED")
'''


def snapshot(path: Path) -> dict:
    with closing(sqlite3.connect(path)) as connection:
        connection.row_factory = sqlite3.Row
        result = {"user_version": connection.execute("PRAGMA user_version").fetchone()[0]}
        for table in ("jobs", "ensembles", "ensemble_members", "queue_meta"):
            rows = connection.execute(f"SELECT * FROM {table}").fetchall()
            result[table] = [dict(row) for row in rows]
    return result


class MigrationTests(unittest.TestCase):
    def setUp(self) -> None:
        if not Path(LEGACY_SOURCE, "simqueue", "store.py").exists():
            self.skipTest("SimQueue 0.5 runtime not present")
        self.temporary = tempfile.TemporaryDirectory(prefix="simq-migration-")
        self.base = Path(self.temporary.name)
        (self.base / "projects").mkdir()
        self.env = {
            "PATH": os.environ.get("PATH", "/usr/bin"),
            "SIMQUEUE_STATE_DIR": str(self.base / "state"),
            "SIMQUEUE_PROJECT_ROOT": str(self.base / "projects"),
        }
        self.saved = dict(os.environ)
        os.environ.update(self.env)

    def tearDown(self) -> None:
        os.environ.clear()
        os.environ.update(self.saved)
        self.temporary.cleanup()

    def legacy(self, template: str) -> str:
        result = subprocess.run(
            [PYTHON_FROZEN, "-c", template % {"source": LEGACY_SOURCE}],
            env=self.env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout

    def assert_migrated(self, before: dict, after: dict) -> None:
        self.assertEqual(after["user_version"], SCHEMA_VERSION)
        self.assertEqual(len(before["jobs"]), len(after["jobs"]))
        for old, new in zip(before["jobs"], after["jobs"]):
            for column, value in old.items():
                self.assertEqual(new[column], value, f"job {old['id']} column {column} changed")
            self.assertTrue(NEW_COLUMNS <= set(new))
            self.assertIsNone(new["partition"])
            self.assertIsNone(new["launch_backend"])
            self.assertEqual(new["priority"], 0)
            self.assertEqual(new["gpus_requested"], 0)
        self.assertEqual(before["ensembles"], after["ensembles"])
        self.assertEqual(before["ensemble_members"], after["ensemble_members"])
        old_meta = {row["key"]: row["value"] for row in before["queue_meta"]}
        new_meta = {row["key"]: row["value"] for row in after["queue_meta"]}
        for key, value in old_meta.items():
            self.assertEqual(new_meta[key], value)
        self.assertIn("schema7_migrated_at", new_meta)

    def test_05_database_migrates_preserving_rows_and_order(self) -> None:
        ids = json.loads(self.legacy(LEGACY_SETUP))["ids"]
        database = self.base / "state" / "queue.sqlite3"
        before = snapshot(database)
        self.assertEqual(before["user_version"], 6)
        store = QueueStore()
        after = snapshot(database)
        self.assert_migrated(before, after)
        QueueStore()  # idempotent
        self.assertEqual(snapshot(database), after)
        # The explicit 0.5 queue order (promote) survives.
        self.assertEqual([int(j["id"]) for j in store.queued_jobs()], [ids[3], ids[2], ids[4]])
        # The running 0.5 job keeps its CPUs for the node scheduler.
        from simqueue.nodeworker import job_allocation, job_request
        from simqueue.scheduler import Request
        from test_scheduler import zenobia

        config = zenobia()
        running = store.require(ids[1])
        self.assertEqual(job_allocation(running, config).cpus, (0, 1, 2, 3))
        self.assertEqual(job_allocation(running, config).mem, 18 * 1024**3)
        queued = store.require(ids[2])
        self.assertEqual(job_request(queued, config),
                         Request(job_id=ids[2], partition="spheral", ncpus=2, cpu_list=None, gpus=0,
                                 gpu_mem=None, mem=18 * 1024**3, time_limit=172800))
        # An idempotent replay from 0.5 still resolves to the same job.
        replay = store.submit(name="done", cwd=self.base / "projects", argv=["/bin/true"],
                              environment={"PATH": "/usr/bin", "A": "1"}, env_override_keys=[],
                              requested_cores=4, idempotency_key="campaign:0")
        self.assertEqual((replay["id"], replay["submission_created"]), (ids[0], False))

    def test_05_code_refuses_a_migrated_database(self) -> None:
        self.legacy(LEGACY_SETUP)
        QueueStore()
        output = self.legacy(LEGACY_OPEN)
        self.assertIn("REFUSED", output)
        self.assertIn("newer than this SimQueue build", output)

    def test_restored_backup_opens_again_under_05(self) -> None:
        self.legacy(LEGACY_SETUP)
        database = self.base / "state" / "queue.sqlite3"
        backup = self.base / "backup.sqlite3"
        with closing(sqlite3.connect(database)) as source, closing(sqlite3.connect(backup)) as target:
            source.backup(target)
        QueueStore()
        for suffix in ("-wal", "-shm"):
            Path(f"{database}{suffix}").unlink(missing_ok=True)
        shutil.copy2(backup, database)
        self.assertIn("OPENED", self.legacy(LEGACY_OPEN))

    def test_copy_of_live_database(self) -> None:
        copy = os.environ.get("SIMQ_LIVE_DB_COPY")
        if not copy:
            self.skipTest("SIMQ_LIVE_DB_COPY not set")
        live = Path("/home/treefor/.local/state/simqueue/queue.sqlite3").resolve()
        self.assertNotEqual(Path(copy).resolve(), live, "refusing to migrate the live database")
        state = self.base / "state"
        state.mkdir()
        target = state / "queue.sqlite3"
        with closing(sqlite3.connect(f"file:{copy}?mode=ro", uri=True)) as source, \
                closing(sqlite3.connect(target)) as destination:
            source.backup(destination)
        before = snapshot(target)
        QueueStore()
        self.assert_migrated(before, snapshot(target))


if __name__ == "__main__":
    unittest.main()
