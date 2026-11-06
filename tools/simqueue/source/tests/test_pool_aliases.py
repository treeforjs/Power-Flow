"""0.6.3: one pool with partition aliases; ``simq run`` exit statuses.

Pure planning tests use an explicit (fake) clock and launch nothing; the
command-line tests submit into a private state directory and start no unit.
The unit-level tests of aliases and ``simq run`` are in test_run_units.py.
"""

from __future__ import annotations

import json
import tempfile
import os
import unittest
from pathlib import Path

from node_fixture import FAKE_GPU, NodeFixture

from simqueue.borrow import Borrower, CoreState, eligible_partitions, plan_node
from simqueue.cli import job_exit_status
from simqueue.nodeconfig import ConfigError, load_node_config, parse_node_config
from simqueue.nodeworker import borrower_entries, job_request
from simqueue.scheduler import INFINITY, Allocation, Request, Running, build_pool, plan


G = 1024**3
RUNTIME = Path(__file__).resolve().parents[2]
ZENOBIA = RUNTIME / "config" / "zenobia-node.json"
ZENOBIA_16 = RUNTIME / "config" / "zenobia-16core-node.json"
ZENOBIA_062 = Path(__file__).resolve().parent / "fixtures" / "zenobia-node-0.6.2.json"
GPU = "GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9"
OLD_NAMES = ("spheral", "raikou-mfem", "powerflow")
NOW = 5_000_000.0


def pool_config(*, backfill: bool = True, path: Path = ZENOBIA):
    """The zenobia pool for planning (no topology check: tests may run anywhere)."""
    document = json.loads(path.read_text())
    document["borrow"]["backfill_waiting_owners"] = backfill
    return parse_node_config(document, path=path, check_topology=False)


def req(job_id, partition="zenobia", ncpus=1, *, minutes=60, mem=G, borrow=False, gpus=0, gpu_mem=None,
        blocked=False, borrow_from=None, max_cpus=None):
    return Request(job_id=job_id, partition=partition, ncpus=ncpus, gpus=gpus, gpu_mem=gpu_mem, mem=mem,
                   time_limit=None if minutes is None else minutes * 60, blocked=blocked, borrow=borrow,
                   borrow_from=borrow_from, max_cpus=max_cpus)


def owner(job_id, cpus, *, end, mem=G):
    return Running(job_id, Allocation(cpus=tuple(cpus), mem=mem), end)


def borrower(job_id, cpus, *, home="zenobia", frozen=None, deadline=None, mode="requeue"):
    return Borrower(job_id=job_id, home=home, cpus=tuple(cpus), mem=G, on_preempt=mode,
                    placement={cpu: "free" for cpu in cpus}, frozen_since=frozen, deadline=deadline,
                    order=job_id)


def node_plan(cfg, queue, *, running=(), borrowers=(), cores=None):
    running = list(running)
    pool = build_pool(cfg, running)
    owner_of = {cpu: item.job_id for item in running for cpu in item.allocation.cpus}
    return plan_node(cfg, pool, running, list(queue), NOW, borrowers=list(borrowers), cores=cores or {},
                     owner_of=owner_of)


def started(decision):
    return {request.job_id: allocation for request, allocation in decision.starts}


class AliasConfigTests(unittest.TestCase):
    def test_zenobia_is_one_pool_and_the_old_partition_names_are_aliases(self) -> None:
        config = load_node_config(ZENOBIA)  # with the topology check of this machine
        self.assertEqual(list(config.partitions), ["zenobia"])
        self.assertEqual(config.default_partition, "zenobia")
        pool = config.partitions["zenobia"]
        self.assertEqual(pool.cpus, tuple(range(8)))
        self.assertEqual(pool.gpus, (GPU,))
        self.assertEqual(pool.aliases, OLD_NAMES)
        self.assertEqual((pool.lend, pool.lend_leased_idle, pool.timing_sensitive), (True, False, False))
        from simqueue.topology import thread_siblings

        # this machine's sibling of core 6 (14 on the 9850X3D, 22 on the 9950X3D2)
        self.assertEqual(config.siblings(6), thread_siblings([6])[6])
        # the defaults each old name had in 0.6.2, its old cores first; the pool's limits
        old = load_node_config(ZENOBIA_062)
        for name in OLD_NAMES:
            view = config.partition(name)
            before = old.partitions[name]
            self.assertEqual((view.name, view.alias), ("zenobia", name))
            self.assertEqual(view.resolve_default_cpus(), before.resolve_default_cpus(), name)
            if name == "spheral":  # 0.7.1: auto from the history, else the pool's 4G (it was 18G)
                self.assertEqual((view.mem_auto, view.default_mem), (True, 4 * G))
            else:
                self.assertEqual(view.default_mem, before.default_mem, name)
            self.assertEqual(view.default_time, before.default_time, name)
            self.assertEqual(view.prefer_cpus, before.cpus, name)
            self.assertEqual((view.max_mem, view.max_time), (20 * G, 14 * 86400))
            self.assertEqual(view.gpus, (GPU,))
            self.assertEqual(config.canonical_partition(name), "zenobia")
        # without -p: the old default partition's defaults, lowest cores first
        default = config.partition(None)
        spheral = old.partitions["spheral"]
        self.assertEqual((default.name, default.alias, default.prefer_cpus), ("zenobia", None, ()))
        # the pool's own default memory was lowered from spheral's 18G to 4G on 2026-11-05
        self.assertEqual((default.resolve_default_cpus(), default.default_mem, default.default_time),
                         (spheral.resolve_default_cpus(), 4 * G, spheral.default_time))
        self.assertEqual(config.canonical_partition(None), "zenobia")
        self.assertEqual(config.canonical_partition("nope"), "nope")
        self.assertIsNone(config.find_partition("nope"))
        with self.assertRaises(ValueError) as caught:
            config.partition("nope")
        self.assertIn("aliases: powerflow, raikou-mfem, spheral", str(caught.exception))
        self.assertEqual(config.partition_names("powerflow"), ("zenobia", *OLD_NAMES))

    def test_lease_files_are_the_062_ones(self) -> None:
        new, old = load_node_config(ZENOBIA), load_node_config(ZENOBIA_062)
        for cpu in range(8):
            self.assertEqual(new.locks.core_lock_paths(cpu), old.locks.core_lock_paths(cpu))
            self.assertEqual(new.locks.borrow_lock_paths(cpu), old.locks.borrow_lock_paths(cpu))
        self.assertEqual(new.locks.gpu_lock_order(new.gpus[GPU]), old.locks.gpu_lock_order(old.gpus[GPU]))
        self.assertEqual(new.line, old.line)
        self.assertEqual(new.memory_schedulable, old.memory_schedulable)

    def test_16_core_config_for_the_9950x3d2(self) -> None:
        """0.8.0: the 16-core file reads cpus, siblings and L3 domains from sysfs."""
        from test_topology import make_sysfs

        document = json.loads(ZENOBIA_16.read_text())
        self.assertEqual(document["partitions"][0]["cpus"], "auto")
        self.assertNotIn("smt_siblings", document)
        self.assertEqual(document["memory_schedulable"], "28G")
        # off on zenobia: a busy sibling halves an integer-heavy owner (2026-11-06); opt in with true
        self.assertFalse(document["borrow"]["lend_siblings"])
        with tempfile.TemporaryDirectory() as base:
            fake = make_sysfs(Path(base) / "sys", cores=16, domains=2)
            document["topology"]["sysfs_cpu"] = str(fake)
            document["smt_siblings"] = {str(c): [c + 16] for c in range(16)}
            config = parse_node_config(document, path=ZENOBIA_16, check_topology=False)
        self.assertEqual(config.node_cpus(), tuple(range(16)))
        self.assertEqual(config.topology.domains, (tuple(range(8)), tuple(range(8, 16))))
        base = "/home/treefor/projects/.workspace/outputs/pic-profile/"
        self.assertEqual(config.locks.core_lock_paths(6), (base + "testing-core.lock",))
        for cpu in (*range(6), *range(7, 16)):
            self.assertEqual(config.locks.core_lock_paths(cpu), (base + f"testing-core{cpu}.lock",))
        eight = pool_config()
        for name in (None, *OLD_NAMES):
            a, b = config.partition(name), eight.partition(name)
            self.assertEqual((a.resolve_default_cpus(), a.default_mem, a.default_time, a.prefer_cpus,
                              a.max_mem, a.max_time, a.gpus),
                             (b.resolve_default_cpus(), b.default_mem, b.default_time, b.prefer_cpus,
                              b.max_mem, b.max_time, b.gpus))
        # on this machine the same file loads with this machine's cores
        here = load_node_config(ZENOBIA_16)
        from simqueue.topology import physical_cores

        self.assertEqual(here.node_cpus(), physical_cores())

    def test_written_sibling_map_must_match_this_machine(self) -> None:
        from simqueue.topology import thread_siblings

        document = json.loads(ZENOBIA.read_text())
        actual = thread_siblings(range(8))
        document["smt_siblings"] = {str(c): list(actual[c]) for c in range(8)}
        self.assertEqual(parse_node_config(document, path=ZENOBIA).siblings(3), actual[3])
        document["smt_siblings"]["3"] = [actual[3][0] + 1]
        with self.assertRaises(ConfigError):
            parse_node_config(document, path=ZENOBIA)
        parse_node_config(document, path=ZENOBIA, check_topology=False)

    def test_alias_validation(self) -> None:
        def with_aliases(aliases, **partition):
            document = json.loads(ZENOBIA.read_text())
            document["partitions"][0]["aliases"] = aliases
            document["partitions"][0].update(partition)
            return document

        bad = [
            with_aliases(["zenobia"]),
            with_aliases(["a", "a"]),
            with_aliases(["bad name"]),
            with_aliases("spheral"),
            with_aliases([{"name": "a", "default_mem": "21G"}]),
            with_aliases([{"name": "a", "default_cpus": 9}]),
            with_aliases([{"name": "a", "default_time": "15-00:00:00"}]),
            with_aliases([{"name": "a", "default_time": "UNLIMITED"}]),
            with_aliases([{"name": "a", "prefer_cpus": [8]}]),
            with_aliases([{"name": "a", "prefer_cpus": [True]}]),
            with_aliases([{"name": "a", "unknown": 1}]),
            with_aliases([{"default_mem": "1G"}]),
        ]
        two = json.loads(ZENOBIA_062.read_text())
        two["partitions"][0]["aliases"] = ["powerflow"]  # a partition's name
        bad.append(two)
        two = json.loads(ZENOBIA_062.read_text())
        two["partitions"][0]["aliases"] = ["old"]
        two["partitions"][1]["aliases"] = ["old"]
        bad.append(two)
        for document in bad:
            with self.subTest(aliases=[p.get("aliases") for p in document["partitions"]]):
                with self.assertRaises(ConfigError):
                    parse_node_config(document, path=ZENOBIA, check_topology=False)

    def test_plain_alias_and_aliases_of_one_of_several_partitions(self) -> None:
        document = json.loads(ZENOBIA_062.read_text())
        document["partitions"][2]["aliases"] = ["pf", {"name": "pf-short", "default_time": "30"}]
        config = parse_node_config(document, path=ZENOBIA_062, check_topology=False)
        plain = config.partition("pf")
        self.assertEqual(plain.name, "powerflow")
        self.assertEqual((plain.cpus, plain.default_mem, plain.default_time, plain.prefer_cpus),
                         ((6, 7), 4 * G, 43200, ()))
        self.assertEqual(config.partition("pf-short").default_time, 1800)
        self.assertEqual(config.partitions["powerflow"].aliases, ("pf", "pf-short"))
        self.assertEqual(config.partitions["spheral"].aliases, ())
        # an alias spans only its partition
        decision = plan(config, build_pool(config, []), [], [req(1, "pf", 2)], NOW)
        self.assertEqual(started(decision)[1].cpus, (6, 7))
        decision = plan(config, build_pool(config, []), [], [req(1, "pf", 3)], NOW)
        self.assertIn("2 cores", decision.skipped[1])


class PoolPlanningTests(unittest.TestCase):
    def decide(self, queue, running=(), cfg=None, **pool_args):
        cfg = cfg or pool_config()
        running = list(running)
        return plan(cfg, build_pool(cfg, running, **pool_args), running, list(queue), NOW)

    def test_jobs_under_every_name_share_one_pool_in_queue_order(self) -> None:
        queue = [
            req(1, "spheral", 4, mem=18 * G),
            req(2, "powerflow", 2, mem=4 * G),
            req(3, "raikou-mfem", 2, mem=4 * G),
            req(4, "zenobia", 1, mem=G),
        ]
        decision = self.decide(queue)
        allocations = started(decision)
        self.assertEqual(list(allocations), [1, 2, 3])
        self.assertEqual(allocations[1].cpus, (0, 1, 2, 3))  # each old name gets its old cores first
        self.assertEqual(allocations[2].cpus, (6, 7))
        self.assertEqual(allocations[3].cpus, (4, 5))
        self.assertEqual(decision.skipped[4], "resources busy")
        # a GPU request under any name: the pool has the GPU
        decision = self.decide([req(5, "spheral", 1, gpus=1, gpu_mem=G), req(6, "powerflow", 1, gpus=1, gpu_mem=G)])
        self.assertEqual([a.gpus for a in started(decision).values()], [(GPU,), (GPU,)])

    def test_preferred_cores_first_then_the_lowest_free(self) -> None:
        running = [owner(9, [6], end=NOW + 600)]
        decision = self.decide([req(1, "powerflow", 2)], running)
        self.assertEqual(started(decision)[1].cpus, (0, 7))
        decision = self.decide([req(1, "zenobia", 2)], running)
        self.assertEqual(started(decision)[1].cpus, (0, 1))
        decision = self.decide([req(1, "spheral", 6)], running)
        self.assertEqual(started(decision)[1].cpus, (0, 1, 2, 3, 4, 5))
        decision = self.decide([req(1, "powerflow", 1)], running, external_cpus={7})
        self.assertEqual(started(decision)[1].cpus, (0,))

    def test_a_paused_ensemble_holds_the_whole_pool(self) -> None:
        decision = self.decide([req(1, "powerflow", 1, blocked=True), req(2, "spheral", 1), req(3, "zenobia", 1)])
        self.assertEqual(started(decision), {})
        self.assertEqual(decision.skipped[2], "held behind a paused ensemble")
        self.assertEqual(decision.skipped[3], "held behind a paused ensemble")

    def test_eight_core_owner_job_reserves_the_pool_and_short_jobs_backfill(self) -> None:
        running = [owner(9, [0, 1], end=NOW + 3600)]
        queue = [req(1, "zenobia", 8, mem=8 * G), req(2, "powerflow", 1, minutes=30),
                 req(3, "spheral", 1, minutes=120)]
        decision = self.decide(queue, running)
        self.assertEqual(list(started(decision)), [2])
        self.assertEqual(decision.reservation, (1, NOW + 3600))
        self.assertEqual(decision.skipped[3], "backfill would delay job 1")
        decision = self.decide(queue, [])
        self.assertEqual(started(decision)[1].cpus, tuple(range(8)))


class PoolBorrowTests(unittest.TestCase):
    """An 8-core owner job alongside borrow jobs in the single pool."""

    def test_eight_core_owner_starts_at_once_and_preempts_every_borrower(self) -> None:
        cfg = pool_config()
        borrowers = [borrower(20, [6, 7], home="powerflow"), borrower(21, [4, 5], home="spheral")]
        decision = node_plan(cfg, [req(1, "zenobia", 8, mem=8 * G)], borrowers=borrowers)
        self.assertEqual(started(decision)[1].cpus, tuple(range(8)))
        preempted = {p.job_id: p for p in decision.preemptions}
        self.assertEqual(sorted(preempted), [20, 21])
        self.assertTrue(all(p.by_job == 1 and p.mode == "requeue" for p in preempted.values()))

    def test_borrow_jobs_use_the_idle_cores_of_a_running_eight_core_owner(self) -> None:
        cfg = pool_config()
        running = [owner(1, range(8), end=NOW + 7200, mem=8 * G)]
        cores = {6: CoreState(6, owner_idle=True), 7: CoreState(7, owner_idle=True)}
        queue = [req(10, "powerflow", 2, borrow=True), req(11, "zenobia", 1, borrow=True)]
        decision = node_plan(cfg, queue, running=running, cores=cores)
        allocation = started(decision)[10]
        self.assertTrue(allocation.borrow)
        self.assertEqual(allocation.placement, ((6, "zenobia", "owner-idle"), (7, "zenobia", "owner-idle")))
        self.assertNotIn(11, started(decision))
        self.assertIn("lendable", decision.skipped[11])

    def test_free_cores_go_to_borrowers_with_the_highest_cores_first(self) -> None:
        cfg = pool_config()
        running = [owner(1, [0, 1, 2, 3], end=NOW + 7200)]
        decision = node_plan(cfg, [req(10, "spheral", 2, borrow=True, mem=2 * G)], running=running)
        self.assertEqual(started(decision)[10].cpus, (6, 7))

    def test_waiting_eight_core_owner_lends_free_cores_to_borrowers_that_end_first(self) -> None:
        running = [owner(2, [0, 1], end=NOW + 3600)]
        queue = [req(1, "zenobia", 8, mem=8 * G), req(10, "powerflow", 2, borrow=True, minutes=30),
                 req(11, "spheral", 2, borrow=True, minutes=120)]
        decision = node_plan(pool_config(backfill=True), queue, running=running)
        self.assertNotIn(1, started(decision))
        self.assertEqual(decision.reservation, (1, NOW + 3600))
        self.assertEqual(started(decision)[10].cpus, (6, 7))
        self.assertNotIn(11, started(decision))
        self.assertIn("may start before this job's time limit ends", decision.skipped[11])
        # 0.6.1 rule (backfill_waiting_owners false): no free core is lent while the owner waits
        decision = node_plan(pool_config(backfill=False), queue, running=running)
        self.assertEqual([r.job_id for r, _ in decision.starts], [])
        self.assertIn("an owner job of the partition is waiting", decision.skipped[10])
        # an owner waiting behind a job without a time limit has no start time to protect
        running = [owner(2, [0, 1], end=INFINITY)]
        decision = node_plan(pool_config(backfill=True), queue, running=running)
        self.assertEqual(sorted(started(decision)), [10, 11])

    def test_borrow_backfill_does_not_delay_the_owner_when_it_becomes_ready(self) -> None:
        cfg = pool_config(backfill=True)
        # the blocking job ended early: the owner starts and preempts the borrower at once
        queue = [req(1, "zenobia", 8, mem=8 * G)]
        decision = node_plan(cfg, queue, borrowers=[borrower(10, [6, 7], home="powerflow")])
        self.assertEqual(started(decision)[1].cpus, tuple(range(8)))
        self.assertEqual([p.job_id for p in decision.preemptions], [10])

    def test_frozen_borrower_resumes_only_when_it_ends_before_the_waiting_owner(self) -> None:
        cfg = pool_config(backfill=True)
        running = [owner(2, [0, 1], end=NOW + 3600)]
        queue = [req(1, "zenobia", 8, mem=8 * G)]
        early = borrower(20, [6, 7], frozen=NOW - 10, deadline=NOW + 1800, mode="freeze")
        decision = node_plan(cfg, queue, running=running, borrowers=[early])
        self.assertEqual(decision.resumes, [(20, ((6, "zenobia", "free"), (7, "zenobia", "free")))])
        late = borrower(20, [6, 7], frozen=NOW - 10, deadline=NOW + 7200, mode="freeze")
        decision = node_plan(cfg, queue, running=running, borrowers=[late])
        self.assertEqual(decision.resumes, [])

    def test_borrow_from_names_resolve_to_the_pool(self) -> None:
        cfg = pool_config()
        probe = Request(job_id=1, partition="powerflow", ncpus=1, borrow=True, borrow_from=("spheral", "raikou-mfem"))
        self.assertEqual(eligible_partitions(cfg, probe), ["zenobia"])
        decision = node_plan(cfg, [req(10, "powerflow", 3, borrow=True, borrow_from=("spheral",))])
        self.assertEqual(started(decision)[10].cpus, (5, 6, 7))
        old = parse_node_config(json.loads(ZENOBIA_062.read_text()), path=ZENOBIA_062, check_topology=False)
        self.assertEqual(eligible_partitions(old, probe), ["powerflow", "spheral", "raikou-mfem"])

    def test_jobs_recorded_under_the_old_names_are_scheduled_in_the_pool(self) -> None:
        cfg = pool_config()
        queued = {"id": 31, "partition": "powerflow", "requested_cores": 2, "cpu_request": None,
                  "gpus_requested": 0, "gpu_mem_bytes": None, "mem_bytes": None, "time_limit_s": None,
                  "borrow": 0, "state": "queued"}
        request = job_request(queued, cfg)
        self.assertEqual(request.partition, "powerflow")  # the recorded name stays
        self.assertEqual((request.mem, request.time_limit), (4 * G, 43200))  # the alias's defaults
        unrecorded = dict(queued, id=32, partition=None)
        self.assertEqual(job_request(unrecorded, cfg).partition, "zenobia")
        running_borrow = {
            "id": 30, "partition": "spheral", "state": "running", "borrow": 1, "cpu_list": "2", "mem_bytes": G,
            "on_preempt": "requeue", "started_at": None, "time_limit_s": 3600, "queue_order": 30,
            "placement_json": json.dumps({"home": "spheral", "cpus": [{"cpu": 2, "partition": "spheral",
                                                                       "class": "free"}]}),
            "borrow_from": json.dumps(["raikou-mfem"]),
        }
        entries = borrower_entries([running_borrow], cfg)
        self.assertEqual(entries[0].home, "zenobia")
        # nothing conflicts with the old placement; an owner job recorded as spheral avoids the
        # borrower's core while other cores are free (as in 0.6.1), else takes it
        decision = node_plan(cfg, [], borrowers=entries)
        self.assertEqual(decision.preemptions, [])
        old_owner = dict(queued, id=33, partition="spheral", requested_cores=4, mem_bytes=8 * G)
        decision = node_plan(cfg, [job_request(old_owner, cfg)], borrowers=entries)
        self.assertEqual(started(decision)[33].cpus, (0, 1, 3, 4))
        self.assertEqual(decision.preemptions, [])
        old_owner = dict(old_owner, requested_cores=8)
        decision = node_plan(cfg, [job_request(old_owner, cfg)], borrowers=entries)
        self.assertEqual(started(decision)[33].cpus, tuple(range(8)))
        self.assertEqual([(p.job_id, p.by_job) for p in decision.preemptions], [(30, 33)])

    def test_unknown_recorded_name_waits_without_blocking_others(self) -> None:
        cfg = pool_config()
        decision = node_plan(cfg, [req(1, "gone", 1), req(2, "powerflow", 1),
                                   req(3, "gone", 1, borrow=True)])
        self.assertEqual(list(started(decision)), [2])
        self.assertEqual(decision.skipped[1], "unknown partition gone")
        self.assertEqual(decision.skipped[3], "unknown partition gone")


class LiveCopyTests(unittest.TestCase):
    """SIMQ_LIVE_DB_COPY names a read-only copy of the live database: every
    recorded job of it must resolve and plan under the 0.6.3 pool."""

    def test_copy_of_live_database_under_the_pool(self) -> None:
        import shutil
        import sqlite3
        import tempfile
        from contextlib import closing

        from simqueue.nodeworker import job_allocation, running_entries
        from simqueue.store import SCHEMA_VERSION, QueueStore

        source = os.environ.get("SIMQ_LIVE_DB_COPY")
        if not source:
            self.skipTest("SIMQ_LIVE_DB_COPY not set")
        live = "/home/treefor/.local/state/simqueue/queue.sqlite3"
        self.assertNotEqual(os.path.realpath(source), os.path.realpath(live), "never open the live file")
        self.assertFalse(os.stat(source).st_mode & 0o222, "the live copy must be read-only")
        with tempfile.TemporaryDirectory(prefix="simq063-live-") as temporary:
            state = Path(temporary) / "state"
            state.mkdir()
            database = state / "queue.sqlite3"
            with closing(sqlite3.connect(f"file:{source}?mode=ro", uri=True)) as origin, \
                    closing(sqlite3.connect(database)) as target:
                origin.backup(target)

            with closing(sqlite3.connect(database)) as connection:
                columns = [row[1] for row in connection.execute("PRAGMA table_info(jobs)")]

            def rows():
                with closing(sqlite3.connect(database)) as connection:
                    return connection.execute(f"SELECT {', '.join(columns)} FROM jobs ORDER BY id").fetchall()

            before = rows()
            saved = dict(os.environ)
            os.environ.update({"SIMQUEUE_STATE_DIR": str(state), "SIMQUEUE_CONFIG": str(ZENOBIA),
                               "SIMQUEUE_NO_AUTOSTART": "1"})
            try:
                store = QueueStore()
                with closing(sqlite3.connect(database)) as connection:
                    self.assertEqual(connection.execute("PRAGMA user_version").fetchone()[0], SCHEMA_VERSION)
                self.assertEqual(rows(), before)  # 0.7.0 and 0.8.0 add columns and change no existing one
                config = pool_config()
                jobs = store.node_jobs() + store.list_jobs(include_all=True, limit=100000)
                names = {job.get("partition") for job in jobs}
                # jobs recorded before 0.6.3 under the old names, since then under the pool's own
                self.assertTrue(names <= {None, "spheral", "raikou-mfem", "powerflow", "zenobia"}, names)
                for job in jobs:
                    request = job_request(job, config)
                    self.assertEqual(config.canonical_partition(request.partition), "zenobia")
                    self.assertIsNotNone(config.find_partition(request.partition))
                    job_allocation(job, config)
                active = [job for job in jobs if job["state"] in ("starting", "running", "cancelling")]
                running_entries(active, config, NOW)
                borrower_entries(active, config)
                queued = [job_request(job, config) for job in jobs if job["state"] == "queued"]
                decision = plan(config, build_pool(config, []), [], [r for r in queued if not r.borrow], NOW)
                self.assertFalse([reason for reason in decision.skipped.values() if "unknown partition" in reason])
            finally:
                os.environ.clear()
                os.environ.update(saved)


class RunExitStatusTests(unittest.TestCase):
    def test_exit_status_of_a_finished_job(self) -> None:
        cases = [
            ({"state": "succeeded", "exit_code": 0}, 0),
            ({"state": "succeeded", "exit_code": None}, 0),
            ({"state": "failed", "exit_code": 3}, 3),
            ({"state": "failed", "exit_code": None}, 1),
            ({"state": "failed", "exit_code": 0}, 1),  # never success for a job that did not succeed
            ({"state": "canceled", "exit_code": -15}, 143),
            ({"state": "failed", "exit_code": -9}, 137),
            ({"state": "interrupted", "exit_code": None}, 1),
            ({"state": "canceled", "exit_code": None}, 1),
        ]
        for job, status in cases:
            with self.subTest(job=job):
                self.assertEqual(job_exit_status(job), status)


class PoolCliTests(NodeFixture, unittest.TestCase):
    """Submissions under aliases (private state; no unit is started)."""

    needs_units = False

    def config_document(self):
        document = super().config_document()
        cpus = ",".join(str(c) for c in sorted(set(self.part_a) | set(self.part_b)))
        document["default_partition"] = "pool"
        document["partitions"] = [{
            "name": "pool", "cpus": cpus, "gpus": [FAKE_GPU], "default_cpus": 1, "default_mem": "256M",
            "max_mem": "2G", "default_time": "10:00", "max_time": "1:00:00",
            "aliases": [
                {"name": "pa", "default_mem": "300M", "default_time": "5:00", "prefer_cpus": list(self.part_a)},
                {"name": "pb", "default_mem": "200M", "default_time": "7:00", "prefer_cpus": list(self.part_b)},
            ],
        }]
        return document

    def submit_json(self, *arguments):
        result = self.cli("submit", "--json", *arguments, "--", "/bin/true")
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads(result.stdout)

    def test_submit_under_an_alias_records_it_with_the_alias_defaults(self) -> None:
        payload = self.submit_json("-p", "pa", "-c", "1")
        self.assertEqual(payload["partition"], "pa")
        self.assertEqual((payload["mem_bytes"], payload["time_limit_s"]), (300 * 1024**2, 300))
        first = payload["job_id"]
        self.assertEqual(self.job(first)["partition"], "pa")
        payload = self.submit_json("-p", "pb", "-c", "1", "--gpus", "1", "--gpu-mem", "1G")
        self.assertEqual((payload["partition"], payload["mem_bytes"], payload["gpus_requested"]),
                         ("pb", 200 * 1024**2, 1))
        payload = self.submit_json("-c", "1")
        self.assertEqual((payload["partition"], payload["mem_bytes"], payload["time_limit_s"]),
                         ("pool", 256 * 1024**2, 600))
        result = self.cli("submit", "-p", "pa", "-c", "1", "--", "/bin/true")
        self.assertIn("in partition pa (pool pool)", result.stdout)
        # limits are the pool's, whatever the name
        result = self.cli("submit", "-p", "pa", "--mem", "3G", "--", "/bin/true")
        self.assertEqual(result.returncode, 2)
        self.assertIn("--mem must be between 1M and 2G in partition pool (as pa)", result.stderr)
        result = self.cli("submit", "-p", "nope", "--", "/bin/true")
        self.assertEqual(result.returncode, 2)
        self.assertIn("aliases: pa, pb", result.stderr)
        show = self.cli("show", str(first))
        self.assertIn("pa (pool pool)", show.stdout)

    def test_test_only_partitions_and_borrow_from_with_aliases(self) -> None:
        result = self.cli("submit", "--test-only", "-p", "pb", "-c", "1", "--", "/bin/true")
        self.assertEqual(result.returncode, 0, result.stderr)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["partition"], "pb")
        self.assertTrue(payload["would_start_now"])
        self.assertEqual(payload["cpu_list"], str(self.part_b[0]))
        self.submit_json("-p", "pa", "-c", "1")
        self.submit_json("-p", "pb", "-c", "1")
        result = self.cli("partitions", "--json", "--sample", "0")
        self.assertEqual(result.returncode, 0, result.stderr)
        status = json.loads(result.stdout)
        self.assertEqual([p["name"] for p in status["partitions"]], ["pool"])
        self.assertEqual(status["partitions"][0]["queued"], 2)
        self.assertEqual([a["name"] for a in status["partitions"][0]["aliases"]], ["pa", "pb"])
        self.assertEqual(status["partitions"][0]["aliases"][0]["default_mem"], 300 * 1024**2)
        text = self.cli("partitions", "--sample", "0").stdout
        self.assertIn("alias pa: defaults 1 core(s), 300M, 00:05:00", text)
        payload = self.submit_json("-p", "pa", "--borrow-from", "pb", "--mem", "256M")
        self.assertEqual(payload["borrow_from"], ["pb"])
        self.assertEqual(self.job(payload["job_id"])["partition"], "pa")

    def test_alias_replay_keeps_its_fingerprint(self) -> None:
        first = self.submit_json("-p", "pa", "-c", "1", "--idempotency-key", "ci:1")
        again = self.submit_json("-p", "pa", "-c", "1", "--idempotency-key", "ci:1")
        self.assertEqual((again["job_id"], again["created"]), (first["job_id"], False))
        other = self.cli("submit", "--json", "-p", "pb", "-c", "1", "--idempotency-key", "ci:1", "--", "/bin/true")
        self.assertEqual(other.returncode, 2)
        self.assertIn("different submission request", other.stderr)


if __name__ == "__main__":
    unittest.main()
