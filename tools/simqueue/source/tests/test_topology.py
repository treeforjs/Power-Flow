"""SimQueue 0.8.0 L3-domain placement on fake sysfs topologies (no unit).

Topologies: one domain (Ryzen 7 9850X3D: 8 cores, 16 CPUs), two (Ryzen 9
9950X3D2: two 8-core chiplets, CPUs 0-15 and their siblings 16-31) and four
(a station-like 32-core part: four 8-core domains).
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

from node_fixture import PYTHON_FROZEN, NodeFixture

from simqueue.config import format_cpu_list
from simqueue.nodeconfig import ConfigError, parse_node_config
from simqueue.nodeworker import l3_groups, placement_environment
from simqueue.scheduler import Allocation, Request, Running, build_pool, plan, request_is_feasible
from simqueue.store import SCHEMA_VERSION, QueueStore
from simqueue.topology import core_domains, discover_l3, thread_siblings

RUNTIME_070 = Path("/home/treefor/projects/runtimes/simqueue-0.7.0-20261105/source")


def make_sysfs(root: Path, *, cores: int, domains: int, smt: bool = True, l3: bool = True) -> Path:
    """A fake /sys/devices/system/cpu: *cores* physical cores (CPU c and its
    sibling c + cores), split evenly into *domains* L3 domains."""
    total = cores * (2 if smt else 1)
    (root).mkdir(parents=True, exist_ok=True)
    (root / "online").write_text(f"0-{total - 1}\n")
    per = cores // domains
    for cpu in range(total):
        core = cpu % cores
        domain = core // per
        base = root / f"cpu{cpu}"
        (base / "topology").mkdir(parents=True)
        siblings = [core, core + cores] if smt else [core]
        (base / "topology" / "thread_siblings_list").write_text(",".join(map(str, siblings)) + "\n")
        for index, (level, kind, size) in enumerate(((1, "Data", "48K"), (1, "Instruction", "32K"),
                                                    (2, "Unified", "1024K"), (3, "Unified", "98304K"))):
            if level == 3 and not l3:
                continue
            entry = base / "cache" / f"index{index}"
            entry.mkdir(parents=True)
            (entry / "level").write_text(f"{level}\n")
            (entry / "type").write_text(kind + "\n")
            (entry / "size").write_text(size + "\n")
            if level == 3:
                shared = list(range(domain * per, (domain + 1) * per))
                if smt:
                    shared += [c + cores for c in shared]
                (entry / "shared_cpu_list").write_text(format_cpu_list(shared) + "\n")
                (entry / "id").write_text(f"{domain}\n")
            else:
                (entry / "shared_cpu_list").write_text(",".join(map(str, siblings)) + "\n")
    return root


def config_for(sysfs: Path, cores: int, **topology) -> object:
    document = {
        "schema": "simqueue-node-config-v1", "node": "t", "default_partition": "pool", "memory_schedulable": "64G",
        "allow_smt_siblings": True,
        "topology": {"sysfs_cpu": str(sysfs), **topology},
        "partitions": [{"name": "pool", "cpus": f"0-{cores - 1}", "default_mem": "1G", "max_mem": "8G",
                        "default_time": "1:00:00", "max_time": "2:00:00",
                        "aliases": [{"name": "fe", "prefer_cpus": [8, 9]}] if cores > 9 else []}],
    }
    return parse_node_config(document, path=Path("node.json"), check_topology=False)


def probe_cpus(payload: dict) -> list[int]:
    from simqueue.config import parse_cpu_list

    return list(parse_cpu_list(payload["cpu_list"]))


def req(job_id: int, ncpus: int, **fields) -> Request:
    return Request(job_id=job_id, partition=fields.pop("partition", "pool"), ncpus=ncpus, mem=1 << 30,
                   time_limit=3600, **fields)


class TopologyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.base = Path(tempfile.mkdtemp(prefix="simq08-topo-"))

    def tearDown(self) -> None:
        shutil.rmtree(self.base, ignore_errors=True)

    def sysfs(self, cores: int, domains: int, **kw) -> Path:
        root = self.base / f"sys-{cores}-{domains}-{len(kw)}"
        return root if root.exists() else make_sysfs(root, cores=cores, domains=domains, **kw)

    def start(self, config, requests, running=()):
        decision = plan(config, build_pool(config, list(running)), list(running), requests, 0.0)
        return {r.job_id: a.cpus for r, a in decision.starts}, decision


class DiscoveryTests(TopologyTestCase):
    def test_one_two_and_four_domains(self) -> None:
        one = discover_l3(self.sysfs(8, 1))
        self.assertEqual([d.cpus for d in one], [tuple(range(16))])
        self.assertEqual(one[0].size_bytes, 96 * 1024**2)
        two = discover_l3(self.sysfs(16, 2))
        self.assertEqual([format_cpu_list(d.cpus) for d in two], ["0-7,16-23", "8-15,24-31"])
        self.assertEqual([d.cache_id for d in two], [0, 1])
        four = discover_l3(self.sysfs(32, 4))
        self.assertEqual(len(four), 4)
        self.assertEqual(format_cpu_list(four[3].cpus), "24-31,56-63")
        self.assertEqual(core_domains(range(32), four),
                         tuple(tuple(range(8 * k, 8 * k + 8)) for k in range(4)))
        self.assertEqual(core_domains(range(16), two), (tuple(range(8)), tuple(range(8, 16))))

    def test_siblings_and_missing_l3(self) -> None:
        root = self.sysfs(16, 2)
        self.assertEqual(thread_siblings([0, 15, 31], root), {0: (16,), 15: (31,), 31: (15,)})
        bare = discover_l3(self.sysfs(4, 1, l3=False))
        self.assertEqual([d.cpus for d in bare], [tuple(range(8))])  # no L3 information: one domain
        self.assertIsNone(bare[0].size_bytes)

    def test_this_machine(self) -> None:
        domains = discover_l3()
        self.assertGreaterEqual(len(domains), 1)
        cpus = [c for d in domains for c in d.cpus]
        self.assertEqual(len(cpus), len(set(cpus)))

    def test_config(self) -> None:
        config = config_for(self.sysfs(16, 2), 16)
        self.assertEqual(config.topology.domains, (tuple(range(8)), tuple(range(8, 16))))
        self.assertEqual((config.topology.heavy(), config.topology.domain_of(9)), (4, 1))
        explicit = parse_node_config({
            "schema": "simqueue-node-config-v1", "node": "t", "default_partition": "p",
            "topology": {"l3_domains": ["0-3", "4-7"], "heavy_cores": 3},
            "partitions": [{"name": "p", "cpus": "0-7", "default_mem": "1G", "max_mem": "2G"}]},
            path=Path("n"), check_topology=False)
        self.assertEqual((explicit.topology.domains, explicit.topology.heavy()), (((0, 1, 2, 3), (4, 5, 6, 7)), 3))
        for bad in ({"l3_domains": ["0-3", "3-7"]}, {"heavy_cores": 0}, {"sysfs_cpu": "rel"}, {"x": 1}):
            with self.assertRaises(ConfigError, msg=str(bad)):
                parse_node_config({"schema": "simqueue-node-config-v1", "node": "t", "default_partition": "p",
                                   "topology": bad, "partitions": [{"name": "p", "cpus": "0-7"}]},
                                  path=Path("n"), check_topology=False)


class PlacementTests(TopologyTestCase):
    def test_one_domain_keeps_the_old_placement(self) -> None:
        config = config_for(self.sysfs(8, 1), 8)
        starts, _ = self.start(config, [req(1, 2), req(2, 4), req(3, 1)])
        self.assertEqual(starts, {1: (0, 1), 2: (2, 3, 4, 5), 3: (6,)})
        self.assertIn("1 L3 domain", request_is_feasible(req(4, 4, l3_domains=2), config.partition("pool"), config))

    def test_two_domains_keep_jobs_inside_one_and_spread_heavy_jobs(self) -> None:
        config = config_for(self.sysfs(16, 2), 16)
        starts, _ = self.start(config, [req(1, 4), req(2, 4), req(3, 6)])
        self.assertEqual(starts[1], (0, 1, 2, 3))
        self.assertEqual(starts[2], (8, 9, 10, 11))  # heavy: the domain with the most free cores
        self.assertEqual(len({c // 8 for c in starts[3]}), 2)  # 6 cores fit no single domain now
        # light jobs pack into the fuller domain and leave a whole domain for a heavy job
        busy = [Running(10, Allocation(cpus=(0, 1, 2, 3, 4), mem=1 << 30), 7200.0)]
        starts, _ = self.start(config, [req(4, 1), req(5, 2), req(6, 8)], running=busy)
        self.assertEqual(starts[4], (5,))
        self.assertEqual(starts[5], (6, 7))
        self.assertEqual(starts[6], tuple(range(8, 16)))

    def test_an_8_core_job_takes_a_whole_chiplet_when_one_is_free(self) -> None:
        config = config_for(self.sysfs(16, 2), 16)
        busy = [Running(10, Allocation(cpus=(8,), mem=1 << 30), 7200.0)]
        starts, _ = self.start(config, [req(1, 8)], running=busy)
        self.assertEqual(starts[1], tuple(range(8)))

    def test_l3_domains_and_split_domains(self) -> None:
        config = config_for(self.sysfs(16, 2), 16)
        starts, _ = self.start(config, [req(1, 8, l3_domains=2)])
        self.assertEqual(starts[1], (0, 1, 2, 3, 8, 9, 10, 11))
        starts, _ = self.start(config, [req(2, 8, split_domains=(6, 2))])
        cpus = starts[2]
        self.assertEqual(sorted(len([c for c in cpus if c // 8 == d]) for d in (0, 1)), [2, 6])
        groups = l3_groups(config, cpus, (6, 2))
        self.assertEqual([len(g) for g in groups], [6, 2])
        groups = l3_groups(config, cpus, (2, 6))
        self.assertEqual([len(g) for g in groups], [2, 6])
        self.assertIsNotNone(request_is_feasible(req(3, 10, split_domains=(9, 1)), config.partition("pool"), config))
        self.assertIsNotNone(request_is_feasible(req(3, 9, split_domains=(3, 3, 3)), config.partition("pool"), config))
        # a split waits until both domains have room
        busy = [Running(10, Allocation(cpus=tuple(range(8)), mem=1 << 30), 7200.0)]
        starts, decision = self.start(config, [req(4, 4, split_domains=(2, 2))], running=busy)
        self.assertEqual(starts, {})
        self.assertEqual(decision.reservation[0], 4)

    def test_four_domains(self) -> None:
        config = config_for(self.sysfs(32, 4), 32)
        starts, _ = self.start(config, [req(i, 6) for i in range(1, 5)])
        self.assertEqual(sorted({c // 8 for c in starts[i]}.pop() for i in range(1, 5)), [0, 1, 2, 3])
        starts, _ = self.start(config, [req(5, 4, l3_domains=4)])
        self.assertEqual(sorted(c // 8 for c in starts[5]), [0, 1, 2, 3])
        starts, _ = self.start(config, [req(6, 24, split_domains=(8, 8, 8))])
        self.assertEqual(len(starts[6]), 24)
        self.assertEqual(len({c // 8 for c in starts[6]}), 3)

    def test_placement_can_be_turned_off(self) -> None:
        config = config_for(self.sysfs(16, 2), 16, l3_placement=False)
        starts, _ = self.start(config, [req(1, 4), req(2, 4)])
        self.assertEqual(starts, {1: (0, 1, 2, 3), 2: (4, 5, 6, 7)})

    def test_backfill_with_domains(self) -> None:
        config = config_for(self.sysfs(16, 2), 16)
        busy = [Running(10, Allocation(cpus=tuple(range(12)), mem=1 << 30), 600.0)]
        # job 1 needs every core (waits for job 10); job 2 (1 h) would delay it; job 3 (5 min) backfills
        starts, decision = self.start(config, [req(1, 16), req(2, 2), Request(job_id=3, partition="pool", ncpus=2,
                                                                              mem=1 << 30, time_limit=300)],
                                      running=busy)
        self.assertEqual(set(starts), {3})
        self.assertEqual(decision.reservation[0], 1)

    def test_environment_lists_the_jobs_cores_per_domain(self) -> None:
        config = config_for(self.sysfs(16, 2), 16)
        job = {"id": 7, "launch_token": "abc", "log_path": "/x", "split_domains": "[2, 6]"}
        env = placement_environment(job=job, config=config, partition="pool",
                                    allocation=Allocation(cpus=(0, 1, 2, 3, 4, 5, 8, 9), mem=1 << 30),
                                    time_limit=60, unit="u")
        self.assertEqual(env["SIMQ_L3_DOMAINS"], "8-9;0-5")
        self.assertEqual(env["SIMQ_L3_DOMAIN_IDS"], "1,0")
        plain = placement_environment(job={"id": 8, "launch_token": "abc", "log_path": "/x"}, config=config,
                                      partition="pool", allocation=Allocation(cpus=(2, 3), mem=1), time_limit=None,
                                      unit="u")
        self.assertEqual((plain["SIMQ_L3_DOMAINS"], plain["SIMQ_L3_DOMAIN_IDS"]), ("2-3", "0"))


class SiblingLendingTests(TopologyTestCase):
    """0.8.0: owner jobs hold whole physical cores; their SMT siblings go to
    borrow jobs unless the owner runs with --exclusive-timing."""

    def config(self, lend: bool = True):
        from simqueue.borrow import plan_node  # noqa: F401 - imported for the tests below

        document = {
            "schema": "simqueue-node-config-v1", "node": "t", "default_partition": "pool",
            "memory_schedulable": "64G", "smt_siblings": {str(c): [c + 16] for c in range(16)},
            "borrow": {"lend_siblings": lend},
            "topology": {"sysfs_cpu": str(self.sysfs(16, 2))},
            "partitions": [{"name": "pool", "cpus": "auto", "default_mem": "1G", "max_mem": "8G",
                            "default_time": "1:00:00", "max_time": "2:00:00"}],
        }
        return parse_node_config(document, path=Path("node.json"), check_topology=False)

    def plan(self, config, running, queue, **kw):
        from simqueue.borrow import plan_node

        pool = build_pool(config, running)
        owner_of = {cpu: r.job_id for r in running for cpu in r.allocation.cpus}
        return plan_node(config, pool, running, queue, 0.0, owner_of=owner_of, **kw)

    def test_siblings_of_a_running_sim_are_lent(self) -> None:
        config = self.config()
        self.assertEqual(config.partitions["pool"].cpus, tuple(range(16)))
        sims = [Running(1, Allocation(cpus=tuple(range(8)), mem=1 << 30), 7200.0),
                Running(2, Allocation(cpus=tuple(range(8, 16)), mem=1 << 30), 7200.0)]
        build = req(3, 4, borrow=True)
        decision = self.plan(config, sims, [build])
        (request, allocation), = decision.starts
        self.assertEqual(allocation.cpus, (28, 29, 30, 31))  # highest first, as for every borrow job
        self.assertEqual({cls for _, _, cls in allocation.placement}, {"sibling"})
        # without lend_siblings the borrow job waits
        self.assertEqual(self.plan(self.config(lend=False), sims, [build]).starts, [])

    def test_exclusive_timing_keeps_its_siblings(self) -> None:
        config = self.config()
        sims = [Running(1, Allocation(cpus=tuple(range(8)), mem=1 << 30), 7200.0),
                Running(2, Allocation(cpus=tuple(range(8, 16)), mem=1 << 30), 7200.0)]
        decision = self.plan(config, sims, [req(3, 4, borrow=True)], exclusive_owner_cpus=range(8))
        (_, allocation), = decision.starts
        self.assertTrue(set(allocation.cpus) <= set(range(24, 32)), allocation.cpus)  # never 16-23

    def test_free_cores_go_first_and_an_exclusive_start_preempts_a_sibling_borrower(self) -> None:
        from simqueue.borrow import Borrower

        config = self.config()
        sims = [Running(1, Allocation(cpus=tuple(range(8)), mem=1 << 30), 7200.0)]
        decision = self.plan(config, sims, [req(3, 2, borrow=True)])
        (_, allocation), = decision.starts
        self.assertEqual({cls for _, _, cls in allocation.placement}, {"free"})  # a free core beats a sibling
        borrower = Borrower(job_id=5, home="pool", cpus=(16,), mem=1 << 30, placement={16: "sibling"})
        timing = req(6, 8, exclusive_timing=True)
        decision = self.plan(config, [], [timing], borrowers=[borrower])
        self.assertEqual([p.job_id for p in decision.preemptions], [5])


class UtilizationAndLimitTests(TopologyTestCase):
    def test_proc_stat(self) -> None:
        from simqueue.nodeworker import read_proc_stat

        path = self.base / "stat"
        path.write_text("cpu  10 0 10 80 0 0 0 0 0 0\ncpu0 5 0 5 40 0 0 0 0\ncpu1 5 0 5 40 0 0 0 0\nintr 1\n")
        stat = read_proc_stat(str(path))
        self.assertEqual(stat["cpus"], {0: (10, 50), 1: (10, 50)})
        self.assertEqual(len(read_proc_stat()["cpus"]), os.cpu_count())

    def test_priority_max_mem(self) -> None:
        document = {"schema": "simqueue-node-config-v1", "node": "t", "default_partition": "p",
                    "memory_schedulable": "28G",
                    "partitions": [{"name": "p", "cpus": "0-1", "default_mem": "4G", "max_mem": "20G",
                                    "priority_max_mem": {"min_priority": 10, "max_mem": "26G"}}]}
        part = parse_node_config(document, path=Path("n"), check_topology=False).partitions["p"]
        self.assertEqual((part.mem_limit(0), part.mem_limit(10), part.mem_limit(10, borrow=True)),
                         (20 << 30, 26 << 30, 20 << 30))
        document["partitions"][0]["priority_max_mem"]["max_mem"] = "30G"  # above memory_schedulable
        with self.assertRaises(ConfigError):
            parse_node_config(document, path=Path("n"), check_topology=False)


class CliTopologyTests(NodeFixture, unittest.TestCase):
    needs_units = False

    def setUp(self) -> None:
        super().setUp()
        self.fake = make_sysfs(self.base / "sysfs", cores=16, domains=2)
        document = self.config_document()
        document["allow_smt_siblings"] = True
        document["partitions"] = [{"name": "pool", "cpus": "0-15", "default_mem": "256M", "max_mem": "2G",
                                   "default_time": "10:00", "max_time": "1:00:00", "default_cpus": 1}]
        document["default_partition"] = "pool"
        document["gpus"] = []
        document["lock_files"].pop("gpu_path", None)
        document["lock_files"].pop("gpu_extra", None)
        document["topology"] = {"sysfs_cpu": str(self.fake)}
        self.write_config(document)

    def test_options_validation_storage_and_display(self) -> None:
        probe = json.loads(self.cli("submit", "--test-only", "--split-domains", "6,2", "--", "/bin/true").stdout)
        self.assertTrue(probe["would_start_now"])
        self.assertEqual(sorted(len([c for c in probe_cpus(probe) if c // 8 == d]) for d in (0, 1)), [2, 6])
        result = self.cli("submit", "--json", "--split-domains", "4,4", "-n", "pilot", "--", "/bin/true")
        self.assertEqual(result.returncode, 0, result.stderr)
        job = self.store.require(json.loads(result.stdout)["job_id"])
        self.assertEqual((job["requested_cores"], json.loads(job["split_domains"])), (8, [4, 4]))
        self.assertIn("split over L3 domains 4+4", self.cli("show", str(job["id"])).stdout)
        result = self.cli("submit", "--json", "-c", "6", "--l3-domains", "2", "--", "/bin/true")
        self.assertEqual(self.store.require(json.loads(result.stdout)["job_id"])["l3_domains"], 2)
        for args, message in ((["--split-domains", "9,1"], "L3 domain"), (["--split-domains", "4"], "two or more"),
                              (["--l3-domains", "2", "--cpus", "0-1"], "cannot be combined"),
                              (["--l3-domains", "3", "-c", "6"], "L3 domain"),
                              (["-c", "3", "--split-domains", "2,2"], "does not match")):
            bad = self.cli("submit", *args, "--", "/bin/true")
            self.assertEqual(bad.returncode, 2, args)
            self.assertIn(message, bad.stderr, args)
        # 0.8.0 utilization: a finished job's average against its cores, live use per core
        seeded = self.store.submit(name="lazy", cwd=self.projects, argv=["/bin/true"], environment={},
                                   env_override_keys=[], requested_cores=4, scheduling={"partition": "pool"})
        with closing(sqlite3.connect(self.store.path)) as connection:
            connection.execute("UPDATE jobs SET state = 'running', cpu_list = '0-3', "
                               "started_at = '2026-11-05T00:00:00Z' WHERE id = ?", (seeded["id"],))
            connection.commit()
        self.store.finish(seeded["id"], state="succeeded", exit_code=0, finished_at="2026-11-05T00:10:00Z",
                          usage={"peak": 1 << 28, "source": "unit", "cpu_usec": 600_000_000, "utilization": 0.25})
        shown = self.cli("show", str(seeded["id"])).stdout
        self.assertIn("1.00 cores on average = 25% of its 4 allocated core(s) (UNDER-USED", shown)
        listing = self.cli("list").stdout
        self.assertIn("USE", listing.splitlines()[0])
        self.assertIn(" 25% ", listing)
        text = self.cli("partitions", "--sample", "0").stdout
        self.assertIn("CPU use over 0.5 s: cores", text)
        self.assertIn("L3 domain 0: cores 0-7 (CPUs 0-7,16-23, 96 MiB): 0/8 busy", text)
        self.assertIn("L3 domain 1: cores 8-15 (CPUs 8-15,24-31, 96 MiB): 0/8 busy", text)
        status = json.loads(self.cli("partitions", "--json", "--sample", "0").stdout)
        self.assertEqual(len(status["partitions"][0]["l3_domains"]), 2)


class Schema10Tests(unittest.TestCase):
    def setUp(self) -> None:
        if not RUNTIME_070.exists():
            self.skipTest("0.7.0 runtime not present")
        self.base = Path(tempfile.mkdtemp(prefix="simq08-mig-"))
        self.saved = dict(os.environ)
        self.env = {k: v for k, v in os.environ.items() if not k.startswith("SIMQUEUE_")}
        self.env.update({"SIMQUEUE_STATE_DIR": str(self.base / "state"), "SIMQUEUE_PROJECT_ROOT": str(self.base)})
        os.environ.update(self.env)

    def tearDown(self) -> None:
        os.environ.clear()
        os.environ.update(self.saved)
        shutil.rmtree(self.base, ignore_errors=True)

    def test_schema_9_database_migrates_to_10_with_rows_unchanged(self) -> None:
        code = ("import sys; sys.path.insert(0, sys.argv[1]); from pathlib import Path;"
                "from simqueue.store import QueueStore; s = QueueStore();"
                "[s.submit(name=f'j{i}', cwd=Path(sys.argv[2]), argv=['/bin/true'], environment={}, "
                "env_override_keys=[], requested_cores=1) for i in range(3)]")
        subprocess.run([PYTHON_FROZEN, "-c", code, str(RUNTIME_070), str(self.base)], env=self.env, check=True)
        database = self.base / "state" / "queue.sqlite3"
        with closing(sqlite3.connect(database)) as connection:
            columns = [row[1] for row in connection.execute("PRAGMA table_info(jobs)")]
            before = connection.execute(f"SELECT {', '.join(columns)} FROM jobs ORDER BY id").fetchall()
            self.assertEqual(connection.execute("PRAGMA user_version").fetchone()[0], 9)
        QueueStore()
        with closing(sqlite3.connect(database)) as connection:
            self.assertEqual(connection.execute("PRAGMA user_version").fetchone()[0], SCHEMA_VERSION)
            after = connection.execute(f"SELECT {', '.join(columns)} FROM jobs ORDER BY id").fetchall()
            new = [row[1] for row in connection.execute("PRAGMA table_info(jobs)")]
        self.assertEqual((SCHEMA_VERSION, before), (10, after))
        self.assertIn("split_domains", new)


if __name__ == "__main__":
    unittest.main()
