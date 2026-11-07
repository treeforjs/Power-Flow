"""SimQueue 0.8.2: ``submit --smt`` (a job's cores and their SMT siblings).

Pure tests: the sibling map, the stored request (bit 2 of exclusive_timing,
no schema change), the placement environment, the borrow planner's
protection, submit/show/doctor, and the wrapper's affinity check (run
directly, no unit).  StartPathTests run ``simq worker --once`` with stub
systemd-run/systemctl (nothing is started) under 0.8.2 and under the 0.8.1
runtime: a job without --smt gets the same unit command and spec from both.
SmtUnitTests (SIMQ_TEST_CPUS) start real --smt units on the leased cores;
one runs Open MPI's --report-bindings inside one (skipped without the MOOSE
Open MPI launcher).
"""

from __future__ import annotations

import json
import os
import re
import sqlite3
import subprocess
import unittest
from contextlib import closing
from pathlib import Path

from node_fixture import FAKE_GPU, PYTHON_FROZEN, NodeFixture

from simqueue.config import format_cpu_list, parse_cpu_list
from simqueue.nodeconfig import SMT_FLAG, parse_node_config, smt_job


OLD_RUNTIME = Path("/home/treefor/projects/runtimes/simqueue-0.8.1-20261106/source")
MOOSE_MPIEXEC = Path("/home/treefor/projects/runtimes/rz-spheral-20261029/bin/moose-mpiexec")
ENV_PROBE = (
    "import json,os\n"
    "keys=sorted(k for k in os.environ if k.startswith('SIMQ_'))\n"
    "print('PROBE '+json.dumps({'env':{k:os.environ[k] for k in keys},'affinity':sorted(os.sched_getaffinity(0))}))\n"
)
SYSTEMD_RUN_STUB = r'''#!%(python)s
import json, os, sys
with open(os.environ.get("STUB_RECORD") or %(record)r, "a") as handle:
    handle.write(json.dumps(sys.argv[1:]) + "\n")
'''
SYSTEMCTL_STUB = r'''#!%(python)s
import sys
args = sys.argv[1:]
if "show" in args:
    if "--value" in args:
        print("idle")
    else:
        names = args[args.index("--") + 1:] if "--" in args else []
        for name in names:
            print(f"Id={name}\nLoadState=not-found\nActiveState=inactive\nSubState=dead\n")
'''


def probe_payload(text: str) -> dict:
    for line in text.splitlines():
        if line.startswith("PROBE "):
            return json.loads(line[len("PROBE "):])
    raise AssertionError("probe output missing:\n" + text)


def sixteen_core_config(*, lend_siblings: bool = True, extra_partition: str | None = None):
    """The 9950X3D2 layout: cores 0-15, siblings c + 16 (no sysfs read)."""
    partitions = [{"name": "pool", "cpus": "0-15", "default_mem": "1G", "max_mem": "8G",
                   "default_time": "1:00:00", "max_time": "2:00:00"}]
    if extra_partition:
        partitions.append({"name": "sib", "cpus": extra_partition, "default_mem": "1G", "max_mem": "8G",
                           "default_time": "1:00:00", "max_time": "2:00:00"})
    document = {
        "schema": "simqueue-node-config-v1", "node": "t", "default_partition": "pool",
        "memory_schedulable": "64G", "allow_smt_siblings": bool(extra_partition),
        "smt_siblings": {str(c): [c + 16] for c in range(16)},
        "borrow": {"lend_siblings": lend_siblings},
        "partitions": partitions,
    }
    return parse_node_config(document, path=Path("node.json"), check_topology=False)


class FlagAndAffinityTests(unittest.TestCase):
    def test_bit_two_of_exclusive_timing_marks_smt(self) -> None:
        self.assertEqual(SMT_FLAG, 2)
        for value, expected in ((None, False), (0, False), (1, False), (2, True), (3, True), ("3", True), ("x", False)):
            self.assertEqual(smt_job({"exclusive_timing": value}), expected, value)
        self.assertFalse(smt_job({}))

    def test_affinity_is_the_cores_and_their_siblings(self) -> None:
        config = sixteen_core_config()
        self.assertEqual(config.smt_affinity(tuple(range(16))), (tuple(range(32)), 2))
        self.assertEqual(config.smt_affinity((2, 3)), ((2, 3, 18, 19), 2))
        self.assertEqual(config.smt_affinity((8, 9, 10, 11)), ((8, 9, 10, 11, 24, 25, 26, 27), 2))

    def test_siblings_that_are_partition_cpus_or_missing_are_not_added(self) -> None:
        config = sixteen_core_config(extra_partition="16")  # allow_smt_siblings: CPU 16 is allocatable
        self.assertEqual(config.smt_affinity((0, 1)), ((0, 1, 17), 1))
        self.assertEqual(config.smt_siblings_of((0, 1)), {0: (), 1: (17,)})
        plain = parse_node_config({
            "schema": "simqueue-node-config-v1", "node": "t", "default_partition": "p", "smt_siblings": {},
            "partitions": [{"name": "p", "cpus": "0-3", "default_mem": "1G", "max_mem": "2G"}]},
            path=Path("node.json"), check_topology=False)
        self.assertEqual(plain.smt_affinity((0, 1)), ((0, 1), 1))

    def test_placement_environment(self) -> None:
        from simqueue.nodeworker import placement_environment
        from simqueue.scheduler import Allocation

        config = sixteen_core_config()

        def env(flag, borrow=0):
            return placement_environment(job={"id": 7, "launch_token": "ab" * 16, "log_path": "/x",
                                              "exclusive_timing": flag, "borrow": borrow},
                                         config=config, partition="pool",
                                         allocation=Allocation(cpus=(2, 3), mem=1 << 30), time_limit=60, unit="u")

        base = env(0)
        self.assertEqual(env(1), base)  # --exclusive-timing: no new variable
        self.assertFalse({k for k in base if "SMT" in k or "THREADS" in k})
        self.assertEqual(env(3, borrow=1), base)  # a borrow row never gets the SMT variables
        smt = env(3)
        self.assertEqual({k: v for k, v in smt.items() if k not in base},
                         {"SIMQ_SMT": "1", "SIMQ_THREADS_PER_CORE": "2", "SIMQ_CPU_LIST_SMT": "2-3,18-19"})
        self.assertEqual({k: v for k, v in smt.items() if k in base}, base)  # SIMQ_CPU_LIST etc. unchanged
        self.assertEqual((smt["SIMQ_CPU_LIST"], smt["SIMQ_NCPUS"], smt["SIMQ_CPU_COUNT"]), ("2-3", "2", "2"))

    def test_an_smt_job_is_protected_like_an_exclusive_timing_job(self) -> None:
        from simqueue.borrow import Borrower, plan_node
        from simqueue.nodeworker import job_request
        from simqueue.scheduler import Allocation, Running, build_pool

        config = sixteen_core_config()
        row = {"id": 6, "requested_cores": 8, "partition": "pool", "mem_bytes": 1 << 30, "time_limit_s": 600,
               "exclusive_timing": 3}
        request = job_request(row, config)
        self.assertTrue(request.exclusive_timing)
        # a starting SMT owner preempts a borrower on one of its siblings
        borrower = Borrower(job_id=5, home="pool", cpus=(16,), mem=1 << 30, placement={16: "sibling"})
        decision = plan_node(config, build_pool(config, []), [], [request], 0.0, borrowers=[borrower])
        self.assertEqual([p.job_id for p in decision.preemptions], [5])
        # a running SMT owner on cores 0-7 (the worker passes its cores as exclusive_owner_cpus): borrow jobs
        # never get 16-23, only the siblings of the other job's cores
        running = [Running(1, Allocation(cpus=tuple(range(8)), mem=1 << 30), 7200.0),
                   Running(2, Allocation(cpus=tuple(range(8, 16)), mem=1 << 30), 7200.0)]
        owners = [{"id": 1, "cpu_list": "0-7", "exclusive_timing": 3}, {"id": 2, "cpu_list": "8-15",
                                                                        "exclusive_timing": 0}]
        exclusive = {cpu for job in owners if job.get("exclusive_timing")
                     for cpu in parse_cpu_list(job["cpu_list"])}  # as NodeWorker._full_pass computes it
        build = job_request({"id": 3, "requested_cores": 4, "partition": "pool", "borrow": 1, "mem_bytes": 1 << 30,
                             "time_limit_s": 600}, config)
        pool = build_pool(config, running)
        owner_of = {cpu: r.job_id for r in running for cpu in r.allocation.cpus}
        decision = plan_node(config, pool, running, [build], 0.0, owner_of=owner_of, exclusive_owner_cpus=exclusive)
        (_, allocation), = decision.starts
        self.assertTrue(set(allocation.cpus) <= set(range(24, 32)), allocation.cpus)


class WrapperCheckTests(NodeFixture, unittest.TestCase):
    """jobwrap.py run directly (no unit, require_unit false) on this process's CPUs."""

    needs_units = False

    def run_wrapper(self, spec_extra: dict) -> tuple[int, dict, str]:
        from simqueue.units import JOBWRAP

        marker = self.base / "ran"
        marker.unlink(missing_ok=True)
        own = sorted(os.sched_getaffinity(0))
        spec = {
            "schema": "simqueue-job-spec-v1", "job_id": 1, "name": "x", "launch_token": "ab" * 16,
            "unit": "none.service", "argv": ["/bin/touch", str(marker)], "env": {"PATH": "/usr/bin"},
            "cwd": str(self.projects), "project_root": str(self.projects),
            "log_path": str(self.base / "x.log"), "status_path": str(self.base / "status.json"),
            "stop_post_path": str(self.base / "stop-post.json"),
            "cpu_list": format_cpu_list(own), "locks": [], "require_unit": False, **spec_extra,
        }
        path = self.base / "spec.json"
        path.write_text(json.dumps(spec))
        result = subprocess.run([PYTHON_FROZEN, str(JOBWRAP), "run", str(path)], stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, check=False)
        return result.returncode, json.loads((self.base / "status.json").read_text()), result.stdout

    def test_smt_spec_is_checked_against_the_whole_affinity(self) -> None:
        own = sorted(os.sched_getaffinity(0))
        if len(own) < 2:
            self.skipTest("needs a process affinity of two or more CPUs")
        # cpu_list (the cores) is a strict subset of the affinity: an --smt spec passes...
        code, status, out = self.run_wrapper({"cpu_list": str(own[0]), "affinity_cpu_list": format_cpu_list(own),
                                              "smt_text": "cores with siblings"})
        self.assertEqual(code, 0, out)
        self.assertTrue((self.base / "ran").exists())
        self.assertIn("SMT:       cores with siblings", out)
        # ...a spec without affinity_cpu_list (every non-SMT job) is checked against cpu_list as in 0.8.1
        code, status, out = self.run_wrapper({"cpu_list": str(own[0])})
        self.assertEqual((code, status["refusal"]), (254, "placement"), out)
        self.assertFalse((self.base / "ran").exists())
        self.assertNotIn("SMT:", out)
        # and an --smt spec whose siblings the unit lacks is refused
        missing = max(own) + 1
        code, status, out = self.run_wrapper({"affinity_cpu_list": format_cpu_list([*own, missing])})
        self.assertEqual((code, status["refusal"]), (254, "placement"), out)
        self.assertIn("does not match the allocation", status["reason"])


class SubmitShowDoctorTests(NodeFixture, unittest.TestCase):
    needs_units = False

    def setUp(self) -> None:
        super().setUp()
        from simqueue.nodeconfig import load_node_config

        config = load_node_config(self.config_path)
        if not all(config.smt_siblings_of(config.partitions["pa"].cpus).values()):
            self.skipTest("the test CPUs have no SMT siblings on this machine")

    def stored(self, job_id: int) -> int:
        with closing(sqlite3.connect(f"file:{self.state / 'queue.sqlite3'}?mode=ro", uri=True)) as connection:
            return connection.execute("SELECT exclusive_timing FROM jobs WHERE id = ?", (job_id,)).fetchone()[0]

    def test_submit_stores_bits_and_reports_smt(self) -> None:
        result = self.cli("submit", "--json", "-p", "pa", "-c", "1", "--smt", "--", "/bin/true")
        self.assertEqual(result.returncode, 0, result.stderr)
        payload = json.loads(result.stdout)
        self.assertEqual((payload["smt"], payload["exclusive_timing"]), (True, True))
        smt = int(payload["job_id"])
        both = self.submit("-p", "pa", "-c", "1", "--smt", "--exclusive-timing", "--", "/bin/true")
        timing = self.submit("-p", "pa", "-c", "1", "--exclusive-timing", "--", "/bin/true")
        result = self.cli("submit", "--json", "-p", "pa", "-c", "1", "--", "/bin/true")
        plain = int(json.loads(result.stdout)["job_id"])
        self.assertNotIn("smt", json.loads(result.stdout))
        self.assertEqual([self.stored(j) for j in (smt, both, timing, plain)], [3, 3, 1, 0])
        shown = self.cli("show", str(smt)).stdout
        self.assertIn("smt: cores with their SMT siblings, exclusive timing", shown)
        self.assertRegex(shown, r"SMT +yes: once started, its CPU affinity is its cores and their SMT siblings, 2 "
                                r"threads per core")
        self.assertIn("Exclusive timing       yes (--smt)", shown)
        shown = self.cli("show", str(timing)).stdout
        self.assertNotRegex(shown, r"(?m)^SMT ")
        self.assertIn("exclusive timing", shown)
        self.assertIn("Exclusive timing       yes: its cores", shown)
        self.assertNotIn("Exclusive timing", self.cli("show", str(plain)).stdout)

    def test_refusals(self) -> None:
        result = self.cli("submit", "-p", "pa", "-c", "1", "--smt", "--borrow", "--", "/bin/true")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--smt is for owner jobs", result.stderr)
        # every sibling of the test cores made a partition CPU (allow_smt_siblings): none is free
        from simqueue.nodeconfig import load_node_config

        config = load_node_config(self.config_path)
        siblings = sorted({s for c in self.part_a + self.part_b for s in config.siblings(c)})
        document = self.config_document()
        document["partitions"].append({"name": "ps", "cpus": format_cpu_list(siblings), "default_mem": "256M",
                                       "max_mem": "2G"})
        self.write_config(document)
        result = self.cli("submit", "-p", "pa", "-c", "1", "--smt", "--", "/bin/true")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("have no free SMT sibling", result.stderr)
        self.assertIn("--smt unavailable", self.cli("doctor").stdout)

    def test_idempotent_replay_distinguishes_smt(self) -> None:
        first = self.submit("--idempotency-key", "k1", "-p", "pa", "-c", "1", "--smt", "--", "/bin/true")
        self.assertEqual(self.submit("--idempotency-key", "k1", "-p", "pa", "-c", "1", "--smt", "--", "/bin/true"),
                         first)
        result = self.cli("submit", "--idempotency-key", "k1", "-p", "pa", "-c", "1", "--", "/bin/true")
        self.assertNotEqual(result.returncode, 0, "a replay without --smt is a different request")

    def test_test_only_and_doctor(self) -> None:
        result = self.cli("submit", "--test-only", "-p", "pa", "-c", "1", "--smt", "--", "/bin/true")
        self.assertEqual(result.returncode, 0, result.stderr)
        payload = json.loads(result.stdout)
        core = self.part_a[0]
        from simqueue.nodeconfig import load_node_config

        config = load_node_config(self.config_path)
        affinity, threads = config.smt_affinity((core,))
        self.assertEqual((payload["smt"], payload["cpu_list"], payload["cpu_list_smt"], payload["threads_per_core"]),
                         (True, str(core), format_cpu_list(affinity), threads))
        plain = json.loads(self.cli("submit", "--test-only", "-p", "pa", "-c", "1", "--", "/bin/true").stdout)
        self.assertNotIn("smt", plain)
        self.assertNotIn("cpu_list_smt", plain)
        job = self.submit("-p", "pa", "-c", "1", "--smt", "--", "/bin/true")
        doctor = self.cli("doctor").stdout
        line = next(item for item in doctor.splitlines() if item.startswith("SMT jobs"))
        self.assertIn("submit --smt: CPU affinity = the job's cores and their SMT siblings", line)
        self.assertIn(f"{threads} threads per core", line)
        self.assertIn(f"queued: {job}", line)
        self.assertIn("none running", line)


class StartPathTests(NodeFixture, unittest.TestCase):
    """``simq worker --once`` with stub systemd-run/systemctl: the unit
    command and spec of each start, compared between 0.8.1 and 0.8.2."""

    needs_units = False
    maxDiff = None
    SUBMISSIONS = (
        ("plain", ("-p", "pa", "-c", "1")),
        ("timing", ("-p", "pa", "-c", "1", "--exclusive-timing")),
        ("gpu", ("-p", "pb", "-c", "1", "--gpus", "1", "--gpu-mem", "1G", "--mem", "300M", "-t", "7:00")),
    )

    def setUp(self) -> None:
        super().setUp()
        if not (OLD_RUNTIME / "simq").exists():
            self.skipTest("0.8.1 runtime not present")

    def stub_config(self, record: Path) -> dict:
        stubs = self.base / f"stubs-{record.stem}"
        stubs.mkdir()
        run, ctl = stubs / "systemd-run", stubs / "systemctl"
        run.write_text(SYSTEMD_RUN_STUB % {"python": PYTHON_FROZEN, "record": str(record)})
        ctl.write_text(SYSTEMCTL_STUB % {"python": PYTHON_FROZEN})
        for path in (run, ctl):
            path.chmod(0o755)
        document = self.config_document()
        # nothing runs: two cores in each partition, whatever SIMQ_TEST_CPUS holds
        document["partitions"][0]["cpus"] = "0-1"
        document["partitions"][1]["cpus"] = "2-3"
        document["launch"].update({"systemd_run": str(run), "systemctl": str(ctl)})
        document["memory"] = {"enabled": False}
        document["borrow"] = {"enabled": False}
        document["gpu_line"] = {"enabled": False}
        return document

    def start_all(self, repository: Path, tag: str, submissions) -> dict[str, dict]:
        """Submit and start every job in a fresh state with *repository*'s code."""
        state = self.base / f"state-{tag}"
        record = self.base / f"systemd-run-{tag}.jsonl"
        (self.base / f"cfg-{tag}").mkdir()
        config = self.base / f"cfg-{tag}" / "node.json"
        config.write_text(json.dumps(self.stub_config(record)))
        env = dict(self.environment, SIMQUEUE_STATE_DIR=str(state), SIMQUEUE_CONFIG=str(config))
        launcher = self.make_launcher(repository)
        ids = {}
        for name, arguments in submissions:
            result = self.cli("submit", "--json", "-n", name, *arguments, "--", "/bin/true", launcher=launcher,
                              env=env)
            self.assertEqual(result.returncode, 0, result.stderr)
            ids[name] = int(json.loads(result.stdout)["job_id"])
        result = self.cli("worker", "--once", launcher=launcher, env=env)
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        commands = [json.loads(line) for line in record.read_text().splitlines()]
        self.assertEqual(len(commands), len(submissions), (state / "worker.log").read_text())
        with closing(sqlite3.connect(f"file:{state / 'queue.sqlite3'}?mode=ro", uri=True)) as connection:
            connection.row_factory = sqlite3.Row
            rows = {int(r["id"]): dict(r) for r in connection.execute("SELECT * FROM jobs")}
        found = {}
        for name, job_id in ids.items():
            row = rows[job_id]
            unit = str(row["unit_name"])
            command = next(c for c in commands if f"--unit={unit}" in c)
            spec = json.loads((Path(str(row["job_dir"])) / "spec.json").read_text())

            def normal(value, token=str(row["launch_token"])):
                text = json.dumps(value, sort_keys=True)
                text = text.replace(token, "<token>").replace(token[:12], "<token12>")
                text = text.replace(str(state), "<state>").replace(str(config.parent), "<cfg>")
                text = text.replace(str(self.base), "<base>")
                text = text.replace(str(repository), "<runtime>")
                return json.loads(re.sub(r"\d{4}-\d\d-\d\dT[0-9:.]+(Z|\+00:00)", "<time>", text))

            spec = normal(spec)
            for key in ("scheduler_version", "scheduler_runtime"):
                spec.pop(key, None)
            found[name] = {"command": normal(command), "spec": spec, "row": row}
        return found

    def test_jobs_without_smt_start_exactly_as_under_081(self) -> None:
        from simqueue.units import JOBWRAP

        new = self.start_all(JOBWRAP.parents[1], "082", self.SUBMISSIONS)
        old = self.start_all(OLD_RUNTIME, "081", self.SUBMISSIONS)
        for name, _ in self.SUBMISSIONS:
            with self.subTest(job=name):
                self.assertEqual(new[name]["command"], old[name]["command"])
                self.assertEqual(new[name]["spec"], old[name]["spec"])
                self.assertEqual(new[name]["row"]["exclusive_timing"], old[name]["row"]["exclusive_timing"])
                self.assertNotIn("affinity_cpu_list", new[name]["spec"])
                self.assertFalse([k for k in new[name]["spec"]["env"] if "SMT" in k or "THREADS_PER" in k])
        self.assertIn("--property=CPUAffinity=" + new["plain"]["row"]["cpu_list"], new["plain"]["command"])

    def test_smt_job_gets_its_siblings(self) -> None:
        from simqueue.nodeconfig import load_node_config
        from simqueue.units import JOBWRAP

        (self.base / "cfg-probe").mkdir()
        probe = self.base / "cfg-probe" / "node.json"
        probe.write_text(json.dumps(self.stub_config(self.base / "probe.jsonl")))
        if load_node_config(probe).smt_affinity((0, 1, 2, 3))[1] < 2:
            self.skipTest("CPUs 0-3 have no SMT siblings on this machine")
        found = self.start_all(JOBWRAP.parents[1], "smt", (("smt", ("-p", "pa", "-c", "1", "--smt")),
                                                           ("plain", ("-p", "pb", "-c", "1"))))
        smt, plain = found["smt"], found["plain"]
        core = int(smt["row"]["cpu_list"])
        affinity, threads = load_node_config(probe).smt_affinity((core,))
        self.assertIn("--property=CPUAffinity=" + " ".join(map(str, affinity)), smt["command"])
        self.assertIn("--property=CPUAffinity=" + plain["row"]["cpu_list"], plain["command"])
        text = format_cpu_list(affinity)
        for key, value in (("SIMQ_SMT", "1"), ("SIMQ_THREADS_PER_CORE", str(threads)), ("SIMQ_CPU_LIST_SMT", text),
                           ("SIMQ_CPU_LIST", str(core)), ("SIMQ_NCPUS", "1")):
            self.assertIn(f"--setenv={key}={value}", smt["command"])
            self.assertEqual(smt["spec"]["env"][key], value)
        self.assertEqual((smt["spec"]["cpu_list"], smt["spec"]["affinity_cpu_list"]), (str(core), text))
        self.assertIn(f"CPU affinity {text}", smt["spec"]["smt_text"])
        self.assertFalse([a for a in plain["command"] if "SMT" in a])
        log = (self.base / "state-smt" / "worker.log").read_text()
        self.assertIn(f"CPUs {core} (SMT: affinity {text})", log)
        self.assertEqual(smt["row"]["exclusive_timing"], 3)
        self.assertEqual(smt["row"]["cpu_list"], str(core))  # the allocation stays the physical core


class SmtUnitTests(NodeFixture, unittest.TestCase):
    """Real --smt units on the leased test cores (their siblings are idle:
    run-tests.py --in-job runs inside a queue job that holds those cores)."""

    min_cpus = 2

    def config_document(self):
        document = super().config_document()
        cpus = ",".join(str(c) for c in self.cpus) or "0,1"
        document["partitions"] = [{"name": "pa", "cpus": cpus, "default_mem": "256M", "max_mem": "2G",
                                   "default_time": "10:00", "max_time": "1:00:00"}]
        return document

    def setUp(self) -> None:
        super().setUp()
        from simqueue.nodeconfig import load_node_config

        config = load_node_config(self.config_path)
        self.affinity, self.threads = config.smt_affinity(self.cpus)
        if self.threads < 2:
            self.skipTest("the test cores have no SMT siblings")

    def test_smt_unit_runs_on_its_cores_and_siblings(self) -> None:
        cores = format_cpu_list(self.cpus)
        smt = self.submit("-c", str(len(self.cpus)), "--smt", "-n", "smt", "--", PYTHON_FROZEN, "-c", ENV_PROBE)
        plain = self.submit("-c", str(len(self.cpus)), "-n", "plain", "--", PYTHON_FROZEN, "-c", ENV_PROBE)
        self.run_drain()
        for job_id in (smt, plain):
            self.assertEqual(self.state_of(job_id), "succeeded", self.log_text(job_id))
        payload = probe_payload(self.log_text(smt))
        self.assertEqual(payload["affinity"], list(self.affinity))
        env = payload["env"]
        self.assertEqual((env["SIMQ_SMT"], env["SIMQ_THREADS_PER_CORE"], env["SIMQ_CPU_LIST_SMT"]),
                         ("1", str(self.threads), format_cpu_list(self.affinity)))
        self.assertEqual((env["SIMQ_CPU_LIST"], env["SIMQ_NCPUS"]), (cores, str(len(self.cpus))))
        self.assertIn(f"SMT:       cores {cores} with their SMT siblings: CPU affinity "
                      f"{format_cpu_list(self.affinity)}", self.log_text(smt))
        shown = self.cli("show", str(smt)).stdout
        self.assertIn(f"CPU affinity {format_cpu_list(self.affinity)}", shown)
        self.assertEqual(self.job(smt)["cpu_list"], cores)
        # the same request without --smt: the cores only, as 0.8.1
        payload = probe_payload(self.log_text(plain))
        self.assertEqual(payload["affinity"], list(self.cpus))
        self.assertFalse([k for k in payload["env"] if "SMT" in k or "THREADS" in k])
        self.assertNotIn("SMT:", self.log_text(plain))

    def test_open_mpi_binds_ranks_to_whole_cores(self) -> None:
        if not MOOSE_MPIEXEC.exists():
            self.skipTest(f"{MOOSE_MPIEXEC} not present")
        ranks = len(self.cpus)
        show = "grep Cpus_allowed_list /proc/self/status"
        script = (f"{MOOSE_MPIEXEC} -n {ranks} --map-by core --bind-to core --report-bindings "
                  f"/bin/sh -c 'echo CORE $({show} | cut -f2)'; "
                  f"{MOOSE_MPIEXEC} -n {ranks} --bind-to none /bin/sh -c 'echo NONE $({show} | cut -f2)'")
        job = self.submit("-c", str(ranks), "--smt", "-n", "smt-mpi", "-t", "5:00", "--shell", script)
        self.run_drain(timeout=300)
        log = self.log_text(job)
        self.assertEqual(self.state_of(job), "succeeded", log)
        bound = sorted(tuple(parse_cpu_list(line.split()[1])) for line in log.splitlines() if line.startswith("CORE "))
        free = [tuple(parse_cpu_list(line.split()[1])) for line in log.splitlines() if line.startswith("NONE ")]
        siblings = self.sibling_map()
        expected = sorted(tuple(sorted({cpu, *siblings[cpu]})) for cpu in self.cpus)
        self.assertEqual(bound, expected, log)  # each rank: one core and its sibling
        self.assertEqual(free, [tuple(self.affinity)] * ranks, log)  # unbound ranks float over all of them
        out = os.environ.get("SIMQ_TEST_RECEIPTS")
        if out:
            Path(out).mkdir(parents=True, exist_ok=True)
            (Path(out) / "smt-mpi-bindings.json").write_text(json.dumps(
                {"cores": list(self.cpus), "affinity": list(self.affinity), "bound": bound, "unbound": free,
                 "log": log.splitlines()}, indent=2))

    def sibling_map(self) -> dict[int, tuple[int, ...]]:
        from simqueue.nodeconfig import load_node_config

        return load_node_config(self.config_path).smt_siblings_of(self.cpus)


if __name__ == "__main__":
    unittest.main()
