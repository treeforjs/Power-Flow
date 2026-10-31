"""Schema 7 (SimQueue 0.6.0) to schema 8 (0.6.1) migration, and back.

The schema-7 database is created by the real, read-only 0.6.0 runtime code.
An optional check runs against a copy of a production database when
SIMQ_LIVE_DB_COPY names one: the file is copied to a temporary directory and
only that copy is migrated.
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

from simqueue.store import BORROW_JOB_COLUMNS, SCHEMA_VERSION, QueueStore


SOURCE_060 = "/home/treefor/projects/runtimes/simqueue-0.6.0-20261030/source"
PYTHON_FROZEN = "/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen"
NEW_COLUMNS = {name for name, _ in BORROW_JOB_COLUMNS}

SETUP_060 = r'''
import json, os, sys
sys.path.insert(0, %(source)r)
from pathlib import Path
from simqueue.store import QueueStore
store = QueueStore()
projects = Path(os.environ["SIMQUEUE_PROJECT_ROOT"])
env = {"PATH": "/usr/bin", "A": "1"}
sched = {"partition": "spheral", "mem_bytes": 1 << 30, "time_limit_s": 600, "priority": 3, "swap_bytes": 0}
done = store.submit(name="done", cwd=projects, argv=["/bin/true"], environment=env, env_override_keys=[],
                    requested_cores=4, idempotency_key="k:0", scheduling=sched, explicit_scheduling={"priority": 3})
store.claim_allocation(done["id"], cpu_list="0-3", gpu_list="", gpu_claim_bytes=0, gpu_exclusive=False,
                       mem_bytes=1 << 30, swap_bytes=0, time_limit_s=600, partition="spheral",
                       unit_prefix="simq-job", job_root=Path(os.environ["SIMQUEUE_STATE_DIR"]) / "jobs",
                       lock_paths=[], runtime_dir="/x", busy_cpus=set())
store.mark_unit_running(done["id"], pid=11, boot_id="b")
store.finish(done["id"], state="succeeded", exit_code=0, result="success")
running = store.submit(name="running", cwd=projects, argv=["/bin/sleep", "9"], environment=env,
                       env_override_keys=["A"], requested_cores=1, scheduling=dict(sched, partition="powerflow"))
store.claim_allocation(running["id"], cpu_list="6", gpu_list="", gpu_claim_bytes=0, gpu_exclusive=False,
                       mem_bytes=1 << 30, swap_bytes=0, time_limit_s=600, partition="powerflow",
                       unit_prefix="simq-job", job_root=Path(os.environ["SIMQUEUE_STATE_DIR"]) / "jobs",
                       lock_paths=[{"path": "/l", "mode": "exclusive", "resource": "cpu:6"}], runtime_dir="/x",
                       busy_cpus=set())
store.mark_unit_running(running["id"], pid=12, boot_id="b")
queued = store.submit(name="queued", cwd=projects, argv=["/bin/echo", "x"], environment=env,
                      env_override_keys=[], requested_cores=2, scheduling=dict(sched, priority=0))
later = store.submit(name="later", cwd=projects, argv=["/bin/echo", "y"], environment=env,
                     env_override_keys=[], requested_cores=1)
store.promote(later["id"])
print(json.dumps({"ids": [done["id"], running["id"], queued["id"], later["id"]]}))
'''

OPEN_060 = r'''
import sys
sys.path.insert(0, %(source)r)
from simqueue.store import QueueStore
try:
    store = QueueStore()
except RuntimeError as error:
    print("REFUSED", error)
else:
    job = store.submit(name="after-rollback", cwd=__import__("pathlib").Path(__import__("os").environ["SIMQUEUE_PROJECT_ROOT"]),
                       argv=["/bin/true"], environment={}, env_override_keys=[], requested_cores=1)
    print("OPENED", job["id"], job["borrow"] if "borrow" in job else "no-borrow-column")
'''


def snapshot(path: Path) -> dict:
    with closing(sqlite3.connect(path)) as connection:
        connection.row_factory = sqlite3.Row
        result = {"user_version": connection.execute("PRAGMA user_version").fetchone()[0]}
        for table in ("jobs", "ensembles", "ensemble_members", "queue_meta"):
            result[table] = [dict(row) for row in connection.execute(f"SELECT * FROM {table}").fetchall()]
    return result


class Migration061Tests(unittest.TestCase):
    def setUp(self) -> None:
        if not Path(SOURCE_060, "simqueue", "store.py").exists():
            self.skipTest("SimQueue 0.6.0 runtime not present")
        self.temporary = tempfile.TemporaryDirectory(prefix="simq-migration-061-")
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

    def run_060(self, template: str) -> str:
        result = subprocess.run([PYTHON_FROZEN, "-c", template % {"source": SOURCE_060}], env=self.env,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout

    @property
    def database(self) -> Path:
        return self.base / "state" / "queue.sqlite3"

    def assert_migrated(self, before: dict, after: dict) -> None:
        self.assertEqual(before["user_version"], 7)
        self.assertEqual(after["user_version"], SCHEMA_VERSION)
        self.assertEqual(SCHEMA_VERSION, 8)
        self.assertEqual(len(before["jobs"]), len(after["jobs"]))
        for old, new in zip(before["jobs"], after["jobs"]):
            for column, value in old.items():
                self.assertEqual(new[column], value, f"job {old['id']} column {column} changed")
            self.assertTrue(NEW_COLUMNS <= set(new))
            self.assertEqual((new["borrow"], new["exclusive_timing"], new["preempt_count"]), (0, 0, 0))
            for column in ("borrow_from", "max_borrow_cores", "on_preempt", "frozen_at", "preempt_pending",
                           "preempt_reason", "placement_json"):
                self.assertIsNone(new[column])
        for table in ("ensembles", "ensemble_members"):
            self.assertEqual(before[table], after[table])
        old_meta = {row["key"]: row["value"] for row in before["queue_meta"]}
        new_meta = {row["key"]: row["value"] for row in after["queue_meta"]}
        for key, value in old_meta.items():
            self.assertEqual(new_meta[key], value)
        self.assertIn("schema8_migrated_at", new_meta)

    def test_060_database_migrates_preserving_rows_order_and_running_jobs(self) -> None:
        ids = json.loads(self.run_060(SETUP_060))["ids"]
        before = snapshot(self.database)
        store = QueueStore()
        after = snapshot(self.database)
        self.assert_migrated(before, after)
        self.assertEqual([j["id"] for j in store.queued_jobs()], [ids[3], ids[2]])
        self.assertEqual(store.job_events(ids[0]), [])
        running = store.get(ids[1])
        self.assertEqual((running["state"], running["cpu_list"], running["borrow"]), ("running", "6", 0))
        # a 0.6.1 borrow job sorts after both 0.6.0 owner jobs
        borrow = store.submit(name="b", cwd=self.base / "projects", argv=["/bin/true"], environment={},
                              env_override_keys=[], requested_cores=1,
                              scheduling={"partition": "spheral", "borrow": 1, "priority": 99})
        self.assertEqual([j["id"] for j in store.queued_jobs()], [ids[3], ids[2], borrow["id"]])

    def test_060_code_refuses_schema_8_and_reopens_after_the_downgrade(self) -> None:
        self.run_060(SETUP_060)
        QueueStore()
        self.assertIn("REFUSED", self.run_060(OPEN_060))
        with closing(sqlite3.connect(self.database)) as connection:
            connection.execute("PRAGMA user_version = 7")
            connection.commit()
        output = self.run_060(OPEN_060)
        self.assertIn("OPENED", output)
        self.assertTrue(output.strip().endswith(" 0"), output)  # new 0.6.0 rows get borrow = 0
        with closing(sqlite3.connect(self.database)) as connection:
            self.assertEqual(connection.execute("PRAGMA user_version").fetchone()[0], 7)
        before = snapshot(self.database)
        QueueStore()  # upgrade again
        after = snapshot(self.database)
        self.assertEqual(after["user_version"], 8)
        self.assertEqual(before["jobs"], after["jobs"])

    def test_copy_of_live_database(self) -> None:
        source = os.environ.get("SIMQ_LIVE_DB_COPY")
        if not source:
            self.skipTest("SIMQ_LIVE_DB_COPY not set")
        live = "/home/treefor/.local/state/simqueue/queue.sqlite3"
        self.assertNotEqual(os.path.realpath(source), os.path.realpath(live), "never migrate the live file")
        self.assertFalse(os.access(source, os.W_OK) and os.stat(source).st_mode & 0o222,
                         "the live copy must be read-only")
        (self.base / "state").mkdir(parents=True, exist_ok=True)
        with closing(sqlite3.connect(f"file:{source}?mode=ro", uri=True)) as origin, \
                closing(sqlite3.connect(self.database)) as target:
            origin.backup(target)
        before = snapshot(self.database)
        QueueStore()
        after = snapshot(self.database)
        self.assert_migrated(before, after)
        with closing(sqlite3.connect(self.database)) as connection:
            connection.execute("PRAGMA user_version = 7")
            connection.commit()
        self.assertIn("OPENED", self.run_060(OPEN_060))
        report = os.environ.get("SIMQ_TEST_RECEIPTS")
        if report:
            Path(report).mkdir(parents=True, exist_ok=True)
            (Path(report) / "live-copy-migration.json").write_text(json.dumps({
                "source": source, "jobs": len(before["jobs"]), "schema_before": before["user_version"],
                "schema_after": after["user_version"], "rows_unchanged": True,
                "downgrade_to_7_opened_by_060": True,
            }, indent=2))


if __name__ == "__main__":
    unittest.main()
