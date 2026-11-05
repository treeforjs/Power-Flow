"""SimQueue 0.7.0 memory tracking: pure parts (no unit), store, schema 9, CLI.

Unit-level tests (real job units allocating memory) are in
test_memory_units.py.
"""

from __future__ import annotations

import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time
import unittest
from contextlib import closing
from pathlib import Path

from node_fixture import PYTHON_FROZEN, REPOSITORY, NodeFixture

from simqueue.memory import (
    GiB,
    MiB,
    Series,
    SpikeState,
    admission,
    admission_margin,
    auto_request,
    describe_basis,
    name_family,
    over_request_warning,
    read_meminfo,
    read_node_memory,
    read_usage,
    round_up,
    slice_path,
    unit_cgroup_path,
    water_marks,
)
from simqueue.nodeconfig import ConfigError, MemoryConfig, parse_node_config
from simqueue.store import SCHEMA_VERSION, QueueStore

RUNTIME_063 = Path("/home/treefor/projects/runtimes/simqueue-0.6.3-20261105/source")


def history_job(job_id: int, name: str, peak: int, *, cores: int = 1, mem: int = 4 * GiB,
                result: str = "success") -> dict:
    return {"id": job_id, "name": name, "requested_cores": cores, "mem_peak_bytes": peak, "mem_bytes": mem,
            "result": result, "state": "succeeded"}


class NameFamilyTests(unittest.TestCase):
    def test_trailing_index_tokens_are_dropped(self) -> None:
        cases = {
            "fe-audit-cold-n4": "fe-audit-cold",
            "lv-interp_300K_xe1e-6": "lv-interp",
            "pore1p5-cpu-a": "pore1p5-cpu",
            "source-attempt-owner-c5": "source-attempt-owner",
            "native-v86-population-cpu": "native-v86-population-cpu",
            "common-v93-checkpoint": "common-v93-checkpoint",
            "rz-1p5mm-pr-tests": "rz-1p5mm-pr-tests",
            "job-1": "job",
            "123": "123",
            "Build-Step-07": "build-step",
        }
        for name, family in cases.items():
            self.assertEqual(name_family(name), family, name)


class SeriesTests(unittest.TestCase):
    def test_minute_buckets_keep_min_and_max(self) -> None:
        series = Series(cap=16)
        base = 1_800_000_000
        self.assertTrue(series.add(base, 100 * MiB))
        self.assertFalse(series.add(base + 10, 300 * MiB))
        self.assertFalse(series.add(base + 20, 50 * MiB))
        self.assertTrue(series.add(base + 60, 200 * MiB))
        self.assertEqual(series.buckets, [[base, 50, 300], [base + 60, 200, 200]])
        self.assertEqual(series.max_since(base), 300 * MiB)
        copy = Series.from_json(series.to_json(), 16)
        self.assertEqual((copy.step, copy.buckets), (series.step, series.buckets))

    def test_bounded_by_coarsening(self) -> None:
        series = Series(cap=16)
        base = 1_800_000_000 // 3600 * 3600
        for minute in range(500):
            series.add(base + 60 * minute, (minute % 7) * MiB)
        self.assertLessEqual(len(series.buckets), 16)
        self.assertGreaterEqual(series.step, 60 * 32)
        self.assertEqual(max(b[2] for b in series.buckets), 6)
        self.assertEqual(min(b[1] for b in series.buckets), 0)
        self.assertLess(len(series.to_json()), 600)

    def test_bad_json_gives_an_empty_series(self) -> None:
        self.assertEqual(Series.from_json("not json").buckets, [])
        self.assertEqual(Series.from_json(None).buckets, [])


class AutoRequestTests(unittest.TestCase):
    config = MemoryConfig()

    def auto(self, history, name="sim-a", cores=1):
        return auto_request(history, name=name, cores=cores, config=self.config, min_mem=256 * MiB,
                            max_mem=20 * GiB, default_mem=4 * GiB)

    def test_no_history_gives_the_default(self) -> None:
        request, basis = self.auto([])
        self.assertEqual(request, 4 * GiB)
        self.assertEqual(basis["source"], "default")
        self.assertIn("partition default", describe_basis(basis))

    def test_largest_recent_peak_times_factor_rounded_up(self) -> None:
        history = [history_job(9, "sim-a", 1000 * MiB), history_job(8, "sim-a", 1500 * MiB)]
        request, basis = self.auto(history)
        self.assertEqual(request, round_up(int(1500 * MiB * 1.25), 256 * MiB))
        self.assertEqual(request, 2048 * MiB)
        self.assertEqual((basis["matched_by"], basis["jobs"], basis["largest_peak"]), ("name", [9, 8], 1500 * MiB))
        self.assertIn("1.25 x largest peak", describe_basis(basis))

    def test_only_the_last_auto_history_jobs_count(self) -> None:
        history = [history_job(i, "sim-a", 100 * MiB) for i in range(10, 5, -1)] + [history_job(1, "sim-a", 9 * GiB)]
        request, basis = self.auto(history)
        self.assertEqual(len(basis["jobs"]), 5)
        self.assertEqual(request, 256 * MiB)  # 125M rounded up, at least min_mem

    def test_family_when_the_exact_name_has_no_history(self) -> None:
        history = [history_job(5, "fe-audit-a0-n2", 2 * GiB, cores=2), history_job(4, "other", 9 * GiB)]
        request, basis = self.auto(history, name="fe-audit-a0-n4", cores=4)
        self.assertEqual(basis["matched_by"], "family")
        self.assertEqual(request, round_up(int(2 * GiB * 1.25), 256 * MiB))

    def test_same_core_count_preferred(self) -> None:
        history = [history_job(3, "mpi", 6 * GiB, cores=8), history_job(2, "mpi", 1 * GiB, cores=2)]
        request, basis = self.auto(history, name="mpi", cores=2)
        self.assertTrue(basis["same_cores"])
        self.assertEqual(request, round_up(int(1 * GiB * 1.25), 256 * MiB))
        request, basis = self.auto(history, name="mpi", cores=4)  # no 4-core run: all of them
        self.assertFalse(basis["same_cores"])
        self.assertEqual(request, round_up(int(6 * GiB * 1.25), 256 * MiB))

    def test_an_oom_kill_counts_with_its_limit(self) -> None:
        history = [history_job(7, "sim-a", 1900 * MiB, mem=2 * GiB, result="oom-kill")]
        request, basis = self.auto(history)
        self.assertTrue(basis["oom_in_history"])
        self.assertEqual(request, round_up(int(2 * GiB * 1.25), 256 * MiB))

    def test_bounded_by_the_partition_limit(self) -> None:
        request, basis = self.auto([history_job(1, "sim-a", 19 * GiB)])
        self.assertEqual(request, 20 * GiB)
        self.assertIn("above the limit", basis["note"])

    def test_over_request_warning(self) -> None:
        history = [history_job(1, "sim-a", 1 * GiB)]
        warning = over_request_warning(18 * GiB, history, name="sim-a", cores=1, config=self.config)
        self.assertIn("18.0x", warning)
        self.assertIsNone(over_request_warning(2 * GiB, history, name="sim-a", cores=1, config=self.config))
        self.assertIsNone(over_request_warning(18 * GiB, [], name="sim-a", cores=1, config=self.config))
        oom = [history_job(1, "sim-a", 1 * GiB, result="oom-kill")]
        self.assertIsNone(over_request_warning(18 * GiB, oom, name="sim-a", cores=1, config=self.config))


class AdmissionTests(unittest.TestCase):
    def test_requests_fit_available_minus_margin_and_ramping(self) -> None:
        result = admission([(1, 4 * GiB), (2, 8 * GiB), (3, 2 * GiB)], available=12 * GiB, margin=2 * GiB,
                           ramping=1 * GiB)
        # 9G to give: job 1 (4G) fits, job 2 (8G) does not, job 3 (2G) still fits
        self.assertEqual(result.admitted, [1, 3])
        self.assertIn(2, result.held)
        self.assertIn("needs 8G", result.held[2])
        self.assertIn("not yet allocated", result.held[2])
        self.assertEqual(result.budget, 3 * GiB)

    def test_water_marks_and_default_margin(self) -> None:
        config = MemoryConfig()
        total = 31 * GiB
        high, low = water_marks(config, total)
        self.assertEqual(high, int(total * 0.07))  # 7% of 31G is above 2G
        self.assertEqual(low, max(4 * GiB, int(total * 0.12)))
        self.assertEqual(admission_margin(config, total), high)
        self.assertEqual(admission_margin(MemoryConfig(admission_margin=1 * GiB), total), 1 * GiB)
        small = MemoryConfig(high_water_available=300 * MiB, high_water_percent=0, low_water_available=500 * MiB,
                             low_water_percent=0)
        self.assertEqual(water_marks(small, 1 * GiB), (300 * MiB, 500 * MiB))

    def test_spike_hysteresis(self) -> None:
        state = SpikeState()
        self.assertIsNone(state.update(900, 300, 500))
        self.assertEqual(state.update(299, 300, 500), "enter")
        self.assertIsNone(state.update(450, 300, 500))  # between the marks: still under pressure
        self.assertTrue(state.active)
        self.assertEqual(state.update(500, 300, 500), "leave")
        self.assertFalse(state.active)


class ReadingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.base = Path(tempfile.mkdtemp(prefix="simq07-mem-"))

    def tearDown(self) -> None:
        shutil.rmtree(self.base, ignore_errors=True)

    def test_meminfo_and_simulated_node(self) -> None:
        path = self.base / "meminfo"
        path.write_text("MemTotal:       32479360 kB\nMemFree:  1 kB\nMemAvailable:   16518920 kB\n")
        self.assertEqual(read_meminfo(path), (32479360 * 1024, 16518920 * 1024))
        node = read_node_memory(1 * GiB, [300 * MiB, 200 * MiB])
        self.assertEqual((node.total, node.available, node.simulated), (1 * GiB, 524 * MiB, True))
        real = read_node_memory()
        self.assertGreater(real.total, real.available)

    def test_cgroup_paths(self) -> None:
        self.assertEqual(slice_path("simq.slice"), Path("simq.slice"))
        self.assertEqual(slice_path("a-b-c.slice"), Path("a.slice/a-b.slice/a-b-c.slice"))
        root = Path("/sys/fs/cgroup/user.slice/user-1000.slice/user@1000.service")
        self.assertEqual(unit_cgroup_path("simq-job-1-ab.service", "simq.slice", root=root),
                         root / "simq.slice" / "simq-job-1-ab.service")

    def test_read_usage_of_a_cgroup_directory(self) -> None:
        cgroup = self.base / "unit.service"
        cgroup.mkdir()
        (cgroup / "memory.current").write_text("123456\n")
        (cgroup / "memory.peak").write_text("654321\n")
        (cgroup / "cpu.stat").write_text("usage_usec 4200\nuser_usec 4000\n")
        (cgroup / "memory.events").write_text("low 0\noom 0\noom_kill 1\n")
        (cgroup / "cgroup.events").write_text("populated 1\nfrozen 0\n")
        usage = read_usage(cgroup)
        self.assertEqual((usage.current, usage.peak, usage.cpu_usec, usage.oom_kill, usage.frozen, usage.populated),
                         (123456, 654321, 4200, 1, False, True))
        self.assertIsNone(read_usage(self.base / "gone.service"))
        own = Path("/proc/self/cgroup").read_text().split("::", 1)[1].strip()
        live = read_usage(Path("/sys/fs/cgroup") / own.lstrip("/"))
        self.assertIsNotNone(live)
        self.assertGreater(live.current, 0)

    def test_sampling_cost_is_tiny(self) -> None:
        """One sample of 8 job cgroups and MemAvailable costs microseconds."""
        own = Path("/sys/fs/cgroup") / Path("/proc/self/cgroup").read_text().split("::", 1)[1].strip().lstrip("/")
        rounds = 300
        started = time.process_time()
        for _ in range(rounds):
            for _job in range(8):
                read_usage(own, full=False)
            read_node_memory()
        per_sample = (time.process_time() - started) / rounds
        record = {"per_sample_s_8_jobs": per_sample, "per_hour_at_3s_s": per_sample * 1200}
        receipts = os.environ.get("SIMQ_TEST_RECEIPTS")
        if receipts:
            Path(receipts).mkdir(parents=True, exist_ok=True)
            (Path(receipts) / "sampling-cost.json").write_text(json.dumps(record, indent=2) + "\n")
        self.assertLess(per_sample, 0.005, record)


class MemoryConfigTests(unittest.TestCase):
    def document(self, memory=None, **partition):
        doc = {"schema": "simqueue-node-config-v1", "node": "t", "default_partition": "p",
               "partitions": [{"name": "p", "cpus": "0", "default_mem": "1G", "max_mem": "4G", **partition}]}
        if memory is not None:
            doc["memory"] = memory
        return doc

    def parse(self, document):
        return parse_node_config(document, path=Path("node.json"), check_topology=False)

    def test_defaults(self) -> None:
        config = self.parse(self.document())
        self.assertEqual(config.memory, MemoryConfig())
        self.assertEqual(config.partitions["p"].min_mem, 256 * MiB)
        self.assertEqual(config.launch.idle_seconds, 30.0)

    def test_values_and_validation(self) -> None:
        config = self.parse(self.document({"sample_seconds": 2, "admission_margin": "1G", "test_node_total": "1G",
                                           "auto_round": "512M"}, min_mem="512M"))
        self.assertEqual((config.memory.sample_seconds, config.memory.admission_margin, config.memory.test_node_total,
                          config.memory.auto_round, config.partitions["p"].min_mem),
                         (2.0, GiB, GiB, 512 * MiB, 512 * MiB))
        for bad in ({"sample_seconds": 0}, {"bogus": 1}, {"low_water_available": "1G"}, {"auto_factor": 9},
                    {"record_seconds": 5, "sample_seconds": 10}):
            with self.assertRaises(ConfigError, msg=str(bad)):
                self.parse(self.document(bad))
        with self.assertRaises(ConfigError):
            self.parse(self.document(min_mem="8G"))


class StoreMemoryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.base = Path(tempfile.mkdtemp(prefix="simq07-store-"))
        self.saved = dict(os.environ)
        os.environ["SIMQUEUE_STATE_DIR"] = str(self.base / "state")
        os.environ["SIMQUEUE_PROJECT_ROOT"] = str(self.base)

    def tearDown(self) -> None:
        os.environ.clear()
        os.environ.update(self.saved)
        shutil.rmtree(self.base, ignore_errors=True)

    def submit(self, store, name="sim", mem=GiB):
        return store.submit(name=name, cwd=self.base, argv=["/bin/true"], environment={}, env_override_keys=[],
                            requested_cores=1, scheduling={"partition": "p", "mem_bytes": mem})

    def test_usage_is_recorded_and_the_peak_only_grows(self) -> None:
        store = QueueStore()
        job = self.submit(store)
        store.record_memory({job["id"]: {"current": 100 * MiB, "peak": 300 * MiB, "source": "sampled",
                                         "sampled_at": "2026-11-05T00:00:00Z", "series": '{"v":1}', "cpu_usec": 5}})
        store.record_memory({job["id"]: {"current": 50 * MiB, "peak": 200 * MiB, "source": "sampled"}})
        row = store.require(job["id"])
        self.assertEqual((row["mem_current_bytes"], row["mem_peak_bytes"], row["mem_peak_source"]),
                         (50 * MiB, 300 * MiB, "sampled"))
        with closing(sqlite3.connect(store.path)) as connection:
            connection.execute("UPDATE jobs SET state = 'running' WHERE id = ?", (job["id"],))
            connection.commit()
        self.assertTrue(store.finish(job["id"], state="succeeded", exit_code=0,
                                     usage={"peak": 310 * MiB, "source": "unit", "cpu_usec": 7_000_000}))
        row = store.require(job["id"])
        self.assertEqual((row["mem_current_bytes"], row["mem_peak_bytes"], row["mem_peak_source"],
                          row["cpu_usage_usec"]), (None, 310 * MiB, "unit", 7_000_000))
        history = store.memory_history()
        self.assertEqual([h["id"] for h in history], [job["id"]])

    def test_data_version_ignores_the_workers_own_writes(self) -> None:
        worker = QueueStore(persistent=True)
        before = worker.data_version()
        job = self.submit(worker)  # the worker's own connection
        worker.set_message(job["id"], "pending: x")
        self.assertEqual(worker.data_version(), before)
        client = QueueStore()
        client.cancel(job["id"])
        self.assertNotEqual(worker.data_version(), before)
        worker.close()

    def test_submit_extra_columns(self) -> None:
        store = QueueStore()
        job = store.submit(name="g", cwd=self.base, argv=["/bin/true"], environment={}, env_override_keys=[],
                           requested_cores=1, extra_columns={"github_json": '{"statuses": []}', "mem_auto_json": None})
        self.assertEqual(store.require(job["id"])["github_json"], '{"statuses": []}')
        with self.assertRaises(ValueError):
            store.submit(name="g", cwd=self.base, argv=["/bin/true"], environment={}, env_override_keys=[],
                         requested_cores=1, extra_columns={"state": "running"})


SETUP_063 = r"""
import sys, json
sys.path.insert(0, %(source)r)
from pathlib import Path
from simqueue.store import QueueStore, SCHEMA_VERSION
s = QueueStore()
ids = []
for i in range(3):
    j = s.submit(name=f"j{i}", cwd=Path(%(cwd)r), argv=["/bin/true"], environment={"A": "1"}, env_override_keys=[],
                 requested_cores=1, scheduling={"partition": "p", "mem_bytes": 1 << 30})
    ids.append(j["id"])
s.cancel(ids[0])
print(json.dumps({"schema": SCHEMA_VERSION, "ids": ids}))
"""
OPEN_063 = r"""
import sys
sys.path.insert(0, %(source)r)
from pathlib import Path
from simqueue.store import QueueStore
try:
    s = QueueStore()
except RuntimeError as error:
    print("REFUSED", error)
else:
    j = s.submit(name="after", cwd=Path(%(cwd)r), argv=["/bin/true"], environment={}, env_override_keys=[],
                 requested_cores=1)
    print("OPENED", len(s.list_jobs(include_all=True)))
"""


class Schema9MigrationTests(unittest.TestCase):
    """0.6.3 (schema 8) -> 0.7.0 (schema 9) and the rollback's in-place downgrade."""

    def setUp(self) -> None:
        if not RUNTIME_063.exists():
            self.skipTest("0.6.3 runtime not present")
        self.base = Path(tempfile.mkdtemp(prefix="simq07-mig-"))
        self.saved = dict(os.environ)
        self.env = {k: v for k, v in os.environ.items() if not k.startswith("SIMQUEUE_")}
        self.env.update({"SIMQUEUE_STATE_DIR": str(self.base / "state"), "SIMQUEUE_PROJECT_ROOT": str(self.base)})
        os.environ.update(self.env)

    def tearDown(self) -> None:
        os.environ.clear()
        os.environ.update(self.saved)
        shutil.rmtree(self.base, ignore_errors=True)

    def run_063(self, template: str) -> str:
        result = subprocess.run([PYTHON_FROZEN, "-c", template % {"source": str(RUNTIME_063), "cwd": str(self.base)}],
                                env=self.env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout

    def rows(self) -> tuple[int, list[dict]]:
        with closing(sqlite3.connect(self.base / "state" / "queue.sqlite3")) as connection:
            connection.row_factory = sqlite3.Row
            version = connection.execute("PRAGMA user_version").fetchone()[0]
            return version, [dict(r) for r in connection.execute("SELECT * FROM jobs ORDER BY id")]

    def test_migration_keeps_rows_and_the_downgrade_reopens_with_063(self) -> None:
        created = json.loads(self.run_063(SETUP_063))
        self.assertEqual(created["schema"], 8)
        version, before = self.rows()
        self.assertEqual(version, 8)
        store = QueueStore()
        version, after = self.rows()
        self.assertEqual((version, SCHEMA_VERSION), (9, 9))
        for old, new in zip(before, after):
            for column, value in old.items():
                self.assertEqual(new[column], value, f"job {old['id']} column {column}")
            for column in ("mem_peak_bytes", "cpu_usage_usec", "github_json", "mem_auto_json"):
                self.assertIsNone(new[column])
            self.assertEqual((new["github_done"], new["github_attempts"]), (0, 0))
        self.assertEqual(store.get_meta("schema9_migrated_at") is not None, True)
        self.assertIn("REFUSED", self.run_063(OPEN_063))
        # rollback: schema back to 8 in place; 0.6.3 opens it and adds a row
        with closing(sqlite3.connect(self.base / "state" / "queue.sqlite3")) as connection:
            connection.execute("PRAGMA user_version = 8")
            connection.commit()
        self.assertIn("OPENED 4", self.run_063(OPEN_063))
        QueueStore()  # install again
        version, again = self.rows()
        self.assertEqual((version, len(again)), (9, 4))

    def test_copy_of_live_database(self) -> None:
        source = os.environ.get("SIMQ_LIVE_DB_COPY")
        if not source:
            self.skipTest("SIMQ_LIVE_DB_COPY not set")
        live = "/home/treefor/.local/state/simqueue/queue.sqlite3"
        self.assertNotEqual(os.path.realpath(source), os.path.realpath(live), "never migrate the live file")
        state = self.base / "state"
        state.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, state / "queue.sqlite3")
        (state / "queue.sqlite3").chmod(0o600)
        version, before = self.rows()
        QueueStore()
        version_after, after = self.rows()
        self.assertEqual((version, version_after, len(before)), (8, 9, len(after)))
        for old, new in zip(before, after):
            self.assertEqual({k: new[k] for k in old}, old)


class CliMemoryTests(NodeFixture, unittest.TestCase):
    """--mem auto, over-request warnings, show and list (no unit runs)."""

    needs_units = False

    def seed(self, name: str, peak: int, *, cores: int = 1, result: str = "success") -> int:
        job = self.store.submit(name=name, cwd=self.projects, argv=["/bin/true"], environment={},
                                env_override_keys=[], requested_cores=cores,
                                scheduling={"partition": "pa", "mem_bytes": 2 * GiB})
        with closing(sqlite3.connect(self.store.path)) as connection:
            connection.execute("UPDATE jobs SET state = 'running', started_at = '2026-11-05T00:00:00Z' WHERE id = ?",
                               (job["id"],))
            connection.commit()
        self.store.finish(job["id"], state="succeeded", exit_code=0, result=result,
                          finished_at="2026-11-05T00:10:00Z",
                          usage={"peak": peak, "source": "unit", "cpu_usec": 90_000_000})
        return int(job["id"])

    def test_mem_auto_uses_the_history_and_a_replay_keeps_its_job(self) -> None:
        seeded = self.seed("sweep-n1", 300 * MiB)
        result = self.cli("submit", "--json", "-p", "pa", "-c", "1", "--mem", "auto", "-n", "sweep-n1",
                          "--idempotency-key", "auto-1", "--", "/bin/true")
        self.assertEqual(result.returncode, 0, result.stderr)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["mem_bytes"], 512 * MiB)  # 375M rounded up to 256M steps
        self.assertEqual(payload["mem_auto"]["jobs"], [seeded])
        job = self.store.require(payload["job_id"])
        self.assertEqual(json.loads(job["mem_auto_json"])["largest_peak"], 300 * MiB)
        self.seed("sweep-n1", 3 * GiB)  # the history changes; the replay still finds its job
        again = json.loads(self.cli("submit", "--json", "-p", "pa", "-c", "1", "--mem", "auto", "-n", "sweep-n1",
                                    "--idempotency-key", "auto-1", "--", "/bin/true").stdout)
        self.assertEqual((again["job_id"], again["created"]), (payload["job_id"], False))
        show = self.cli("show", str(payload["job_id"])).stdout
        self.assertIn("Memory request", show)
        self.assertIn("1.25 x largest peak", show)
        text = self.cli("submit", "-p", "pa", "-c", "1", "--mem", "auto", "-n", "sweep-n9", "--", "/bin/true")
        self.assertIn("(auto)", text.stdout)  # the family sweep-n has history

    def test_no_mem_uses_the_history_else_the_default(self) -> None:
        """0.7.1: memory.auto_default and an alias default_mem "auto"."""
        document = self.config_document()
        document["memory"] = {"auto_default": True}
        document["partitions"][0]["aliases"] = [{"name": "sph", "default_mem": "auto"}]
        self.write_config(document)
        self.seed("known", 300 * MiB)
        known = json.loads(self.cli("submit", "--json", "-p", "pa", "-n", "known", "--", "/bin/true").stdout)
        self.assertEqual(known["mem_bytes"], 512 * MiB)
        self.assertTrue(known["mem_auto"]["implicit"])
        self.assertIn("auto (default without --mem)", self.cli("show", str(known["job_id"])).stdout)
        unknown = json.loads(self.cli("submit", "--json", "-p", "pa", "-n", "fresh", "--", "/bin/true").stdout)
        self.assertEqual(unknown["mem_bytes"], 256 * MiB)  # the partition default, no basis recorded
        self.assertNotIn("mem_auto", unknown)
        explicit = json.loads(self.cli("submit", "--json", "-p", "pa", "-n", "known", "--mem", "1G", "--",
                                       "/bin/true").stdout)
        self.assertEqual(explicit["mem_bytes"], 1 << 30)  # an explicit --mem always wins
        document["memory"] = {"auto_default": False}
        self.write_config(document)
        plain = json.loads(self.cli("submit", "--json", "-p", "pa", "-n", "known", "--", "/bin/true").stdout)
        self.assertEqual(plain["mem_bytes"], 256 * MiB)
        alias = json.loads(self.cli("submit", "--json", "-p", "sph", "-n", "known", "--", "/bin/true").stdout)
        self.assertEqual(alias["mem_bytes"], 512 * MiB)  # the alias asks for auto by itself
        self.assertIn("auto (else 256M)", self.cli("partitions", "--sample", "0").stdout)

    def test_over_request_warning_on_submit(self) -> None:
        self.seed("big-ask", 200 * MiB)
        # pa allows 2G at most in the fixture: write a config with room
        document = self.config_document()
        document["partitions"][0]["max_mem"] = "6G"
        self.write_config(document)
        result = self.cli("submit", "--json", "-p", "pa", "--mem", "5G", "-n", "big-ask", "--", "/bin/true")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("simq: warning: --mem 5G is", result.stderr)
        self.assertIn("warnings", json.loads(result.stdout))
        quiet = self.cli("submit", "--json", "-p", "pa", "--mem", "512M", "-n", "big-ask", "--", "/bin/true")
        self.assertNotIn("warning", quiet.stderr)

    def test_show_and_list_report_measured_use(self) -> None:
        job_id = self.seed("measured", 700 * MiB)
        show = self.cli("show", str(job_id)).stdout
        self.assertIn("Memory peak", show)
        self.assertIn("700M", show)
        self.assertIn("memory.peak recorded by the unit", show)
        self.assertIn("CPU time               90.0 s, 0.15 cores on average", show)
        listing = self.cli("list").stdout
        self.assertIn("MEM", listing.splitlines()[0])
        self.assertIn("700M/2G", listing)
        partitions = json.loads(self.cli("partitions", "--json", "--sample", "0").stdout)
        memory = partitions["node_memory"]
        self.assertGreater(memory["total"], memory["available"])
        self.assertIn("high_water", memory)


if __name__ == "__main__":
    unittest.main()
