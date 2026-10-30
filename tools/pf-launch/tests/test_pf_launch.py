"""Unit tests for pf_launch with fake topologies, fake GPU inventories and fake probes.

No program is launched.  Lease tests use real flock(2) on files in a temporary directory.
Run:  python3 -m unittest discover -s tests -v     (from the pf-launch directory)
"""

import io
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import pf_launch as pf  # noqa: E402

UUIDS = [f"GPU-0000000{i}-aaaa-bbbb-cccc-00000000000{i}" for i in range(4)]
SIMQ = ("/opt/simq", "active")
IDLE = (0, "No active jobs.\n", "")
ACTIVE = (0, "   ID  STATE       CPU       TIME  NAME                     DIRECTORY\n"
             "   14  running     0-5      3m23s  spheral-run              runs/x\n", "")


def zen_topology():
    """8 cores / 16 threads, siblings (k, k+8), one NUMA node (like the workstation)."""
    return pf.Topology.from_siblings([(k, k + 8) for k in range(8)], cpu_model="Fake 8-core")


def station_topology():
    """64 cores / 128 threads, siblings (k, k+64), two NUMA nodes."""
    return pf.Topology.from_siblings([(k, k + 64) for k in range(64)], numa_of=lambda c: 0 if c < 32 else 1,
                                     cpu_model="Fake 64-core")


def intel_topology():
    """Adjacent siblings (0,1), (2,3), ..."""
    return pf.Topology.from_siblings([(2 * k, 2 * k + 1) for k in range(4)])


def fake_gpus(n):
    return [{"index": i, "uuid": UUIDS[i], "name": "Fake GPU", "compute_cap": "12.0", "memory_total": "24463 MiB",
             "driver_version": "580.00", "pci_bus_id": f"00000000:0{i}:00.0", "numa": 0} for i in range(n)]


class FakeSystem(pf.System):
    def __init__(self, topology=None, gpus=None, env=None, probes=None, busy=(), nice=0, files=None,
                 existing=(), host="testhost", cgroup=None, inventory_note=None):
        self.sysfs = "/fake/sys"
        self.procfs = "/fake/proc"
        self._topology = topology or zen_topology()
        self._gpus = fake_gpus(1) if gpus is None else gpus
        self.environ = dict(env or {"PATH": "/usr/bin"})
        self.probes = dict(probes or {})
        self.busy = set(busy)
        self._nice = nice
        self.files = {"/fake/sys/module/nvidia/version": "580.00", "/fake/sys/module/nvidia/srcversion": "ABC123"}
        self.files.update(files or {})
        self.existing = set(existing) | {"/frozen/run-frozen", "/var/simq/queue.sqlite3"}
        self.host = host
        self.cgroup = cgroup
        self.inventory_note = inventory_note
        self.calls = []

    def hostname(self):
        return self.host

    def topology(self):
        return self._topology

    def cpu_model(self):
        return self._topology.cpu_model

    def kernel_release(self):
        return self.files.get("kernel", "6.1.0-test")

    def read_text(self, path):
        return self.files.get(path)

    def exists(self, path):
        return path in self.existing

    def sha256_file(self, path):
        return self.files.get(("sha", path))

    def run(self, argv, timeout):
        self.calls.append(tuple(argv))
        result = self.probes.get(tuple(argv))
        if isinstance(result, BaseException):
            raise result
        if result is None:
            raise FileNotFoundError(argv[0])
        return result

    def gpu_inventory(self, nvidia_smi="nvidia-smi"):
        return [dict(g) for g in self._gpus], self.inventory_note

    def lease_busy(self, path):
        return path in self.busy

    def current_nice(self):
        return self._nice

    def cgroup_cpus(self):
        return self.cgroup

    def same_dir(self, a, b):
        return False

    def now(self):
        return "2026-10-30T00:00:00.000000Z"


def zen_site(**overrides):
    site = {
        "schema": pf.SITE_SCHEMA, "site": "testhost", "hostname": "testhost",
        "nice": 15, "default_gpus": "auto", "core_preference": [6, 7],
        "reserved_cores": {"cores": [0, 1, 2, 3, 4, 5], "reason": "allocated elsewhere"},
        "reserved_gpus": {"gpus": [], "reason": ""},
        "reservations": [{"name": "queue-active", "cores": [0, 1, 2, 3],
                          "probe": {"argv": list(SIMQ), "timeout_s": 5, "require_paths": ["/var/simq/queue.sqlite3"],
                                    "busy_regex": r"^\s*\d+\s+\S+", "free_regex": r"No active jobs\."}}],
        "leases": {"core_path": "/leases/testing-core{core}.lock",
                   "core_overrides": {"6": "/leases/testing-core.lock"},
                   "gpu_path": "/leases/testing-gpu-{uuid}.lock",
                   "gpu_extra": {UUIDS[0]: ["/leases/testing-gpu.lock"]}},
        "launcher": {"mpirun": "mpirun", "mpirun_args": [], "srun": "srun",
                     "srun_args": ["--cpu-bind=cores", "--ntasks-per-core=1", "--kill-on-bad-exit=1"],
                     "srun_in_allocation": True, "nvidia_smi": "nvidia-smi", "wrapper_python": None},
        "execution": {"kind": "frozen-root", "prefix": ["/frozen/run-frozen"], "mpi_inside_prefix": True,
                      "wrapper_python": "/usr/bin/python3", "root": "/frozen/root", "manifest": "/frozen/root-files.json"},
        "environment": {"set": {"OMP_NUM_THREADS": "4", "SITE_VAR": "x"}, "unset": ["UNWANTED"]},
        "identity": {"policy": "strict", "record": {"kernel_release": "6.1.0-test", "nvidia_driver_version": "580.00",
                                                    "nvidia_srcversion": "ABC123", "gpu_uuids": [UUIDS[0]],
                                                    "frozen_root_manifest_sha256": "f" * 64}},
        "receipt_dir": None,
    }
    site.update(overrides)
    return site


def zen_system(**kwargs):
    files = {("sha", "/frozen/root-files.json"): "f" * 64}
    files.update(kwargs.pop("files", {}))
    kwargs.setdefault("probes", {SIMQ: IDLE})
    return FakeSystem(files=files, **kwargs)


def options(*args):
    opts = pf.make_parser().parse_args(list(args))
    opts.argv = ["pf-launch"] + list(args)
    return opts


SITE_INFO = {"path": "/sites/testhost.json", "sha256": "0" * 64, "schema": pf.SITE_SCHEMA, "site": "testhost"}


def plan_for(*args, site=None, system=None):
    return pf.build_plan(options(*args), site or zen_site(), SITE_INFO, system or zen_system())


def codes(result):
    return [r["code"] for r in result["refusals"]]


# --------------------------------------------------------------------------- CPU lists and topology


class CpuListTests(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(pf.parse_cpu_list("6-7"), [6, 7])
        self.assertEqual(pf.parse_cpu_list("0,2,5"), [0, 2, 5])
        self.assertEqual(pf.parse_cpu_list("4-7"), [4, 5, 6, 7])
        self.assertEqual(pf.parse_cpu_list(" 0-2, 8 "), [0, 1, 2, 8])
        self.assertEqual(pf.parse_cpu_list("7,6"), [7, 6])  # order kept: rank 0 on CPU 7

    def test_parse_errors(self):
        for bad in ("", "7-6", "a", "1,,2", "3,3", "0-3,2", "-1", "1-", "1.5"):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                pf.parse_cpu_list(bad)

    def test_format(self):
        self.assertEqual(pf.format_cpu_list([0, 1, 2, 3, 5, 7, 8]), "0-3,5,7-8")
        self.assertEqual(pf.format_cpu_list([]), "")

    def test_sysfs_topology(self):
        with tempfile.TemporaryDirectory() as root:
            sysfs, procfs = os.path.join(root, "sys"), os.path.join(root, "proc")
            base = os.path.join(sysfs, "devices/system/cpu")
            os.makedirs(base)
            with open(os.path.join(base, "online"), "w") as f:
                f.write("0-6\n")  # CPU 7 offline
            for cpu in range(8):
                topo = os.path.join(base, f"cpu{cpu}", "topology")
                os.makedirs(topo)
                for name, value in (("thread_siblings_list", f"{cpu % 4},{cpu % 4 + 4}"), ("core_id", cpu % 4),
                                    ("physical_package_id", 0)):
                    with open(os.path.join(topo, name), "w") as f:
                        f.write(f"{value}\n")
            for node, cpus in ((0, "0-1,4-5"), (1, "2-3,6-7")):
                os.makedirs(os.path.join(sysfs, f"devices/system/node/node{node}"))
                with open(os.path.join(sysfs, f"devices/system/node/node{node}/cpulist"), "w") as f:
                    f.write(cpus + "\n")
            os.makedirs(procfs)
            with open(os.path.join(procfs, "cpuinfo"), "w") as f:
                f.write("processor\t: 0\nmodel name\t: Fake CPU 9000\n")
            topology = pf.read_topology(sysfs, procfs)
        self.assertEqual(topology.primaries(), [0, 1, 2, 3])
        self.assertEqual(topology.cores[0]["cpus"], [0, 4])
        self.assertEqual(topology.cores[3]["cpus"], [3])  # sibling 7 offline
        self.assertEqual(topology.cores[2]["numa"], 1)
        self.assertEqual(topology.cpu_model, "Fake CPU 9000")
        self.assertEqual(topology.cpu_core[6], 2)
        self.assertEqual(topology.to_json()["numa_nodes"], {"0": "0-1,4-5", "1": "2-3,6"})


# --------------------------------------------------------------------------- core selection


class SmtSiblingTests(unittest.TestCase):
    def test_siblings_refused(self):
        result = pf.select_cores(zen_topology(), ("list", [6, 14]))
        self.assertEqual(codes(result), [pf.EXIT_USAGE])
        self.assertIn("SMT siblings of physical core 6", result["refusals"][0]["message"])
        self.assertIn("thread_siblings_list 6,14", result["refusals"][0]["message"])
        self.assertEqual(result["ranks"], [])

    def test_adjacent_siblings_refused(self):
        result = pf.select_cores(intel_topology(), ("list", [0, 1]))
        self.assertEqual(codes(result), [pf.EXIT_USAGE])
        ok = pf.select_cores(intel_topology(), ("list", [0, 2]))
        self.assertEqual([r["cpu"] for r in ok["ranks"]], [0, 2])

    def test_secondary_thread_alone_is_one_core(self):
        result = pf.select_cores(zen_topology(), ("list", [14]))
        self.assertEqual(result["refusals"], [])
        self.assertEqual(result["ranks"], [{"rank": 0, "cpu": 14, "core": 6, "siblings": [6, 14]}])

    def test_count_never_uses_siblings(self):
        result = pf.select_cores(zen_topology(), ("count", 8))
        cpus = [r["cpu"] for r in result["ranks"]]
        self.assertEqual(sorted(cpus), list(range(8)))
        self.assertEqual(len({r["core"] for r in result["ranks"]}), 8)


class CoreCountTests(unittest.TestCase):
    def test_counts_on_station(self):
        topology = station_topology()
        for n in (1, 2, 3, 5, 8):
            with self.subTest(n=n):
                result = pf.select_cores(topology, ("count", n))
                self.assertEqual(result["refusals"], [])
                cores = [r["core"] for r in result["ranks"]]
                self.assertEqual(len(cores), n)
                self.assertEqual(len(set(cores)), n)
                self.assertTrue(all(r["cpu"] == r["core"] < 64 for r in result["ranks"]))
                self.assertEqual(cores, list(range(n)))  # NUMA node 0 first, ascending
                self.assertEqual([r["rank"] for r in result["ranks"]], list(range(n)))

    def test_counts_with_reservation(self):
        reserved = {c: ["static: elsewhere"] for c in range(6)}
        topology = zen_topology()
        one = pf.select_cores(topology, ("count", 1), reserved, preference=[6, 7])
        self.assertEqual([r["cpu"] for r in one["ranks"]], [6])
        two = pf.select_cores(topology, ("count", 2), reserved, preference=[6, 7])
        self.assertEqual([r["cpu"] for r in two["ranks"]], [6, 7])
        for n in (3, 5, 8):
            with self.subTest(n=n):
                refused = pf.select_cores(topology, ("count", n), reserved, preference=[6, 7])
                self.assertEqual(codes(refused), [pf.EXIT_BUSY])
                message = refused["refusals"][0]["message"]
                self.assertIn(f"--cores {n}: only 2 of 8 physical cores available (6-7)", message)
                self.assertIn("cores 0-5 (static: elsewhere)", message)
                self.assertIn("at most 2 physical core(s) are unreserved on this host (6-7)", message)
                self.assertIn("--cores-available", message)
                self.assertFalse(refused["refusals"][0]["retryable"])
        fewer = pf.select_cores(topology, ("count", 8), reserved, preference=[6, 7], allow_fewer=True)
        self.assertEqual([r["cpu"] for r in fewer["ranks"]], [6, 7])
        self.assertIn("--cores-available", fewer["notes"][0])

    def test_more_than_exists(self):
        refused = pf.select_cores(zen_topology(), ("count", 9))
        self.assertEqual(codes(refused), [pf.EXIT_ABSENT])
        shrunk = pf.select_cores(zen_topology(), ("count", 9), allow_fewer=True)
        self.assertEqual(len(shrunk["ranks"]), 8)

    def test_auto(self):
        reserved = {c: ["r"] for c in range(6)}
        result = pf.select_cores(zen_topology(), ("auto", None), reserved, preference=[6, 7])
        self.assertEqual([r["cpu"] for r in result["ranks"]], [6, 7])
        none_free = pf.select_cores(zen_topology(), ("auto", None), {c: ["r"] for c in range(8)})
        self.assertEqual(codes(none_free), [pf.EXIT_BUSY])

    def test_count_one_picks_any_free_core(self):
        reserved = {c: ["r"] for c in range(6)}
        result = pf.select_cores(zen_topology(), ("count", 1), reserved, busy={6: "/l/core6"}, preference=[6, 7])
        self.assertEqual([r["cpu"] for r in result["ranks"]], [7])
        both = pf.select_cores(zen_topology(), ("count", 1), reserved, busy={6: "/a", 7: "/b"}, preference=[6, 7])
        self.assertEqual(codes(both), [pf.EXIT_BUSY])
        self.assertTrue(both["refusals"][0]["retryable"])  # --wait polls

    def test_leased_and_reserved_lists(self):
        reserved = {c: ["static: elsewhere"] for c in range(6)}
        busy = {7: "/leases/testing-core7.lock"}
        topology = zen_topology()
        result = pf.select_cores(topology, ("count", 2), reserved, busy, preference=[6, 7])
        self.assertEqual(codes(result), [pf.EXIT_BUSY])
        self.assertIn("leased: core 7 (/leases/testing-core7.lock)", result["refusals"][0]["message"])
        self.assertTrue(result["refusals"][0]["retryable"])
        auto = pf.select_cores(topology, ("auto", None), reserved, busy, preference=[6, 7])
        self.assertEqual([r["cpu"] for r in auto["ranks"]], [6])
        listed = pf.select_cores(topology, ("list", [7]), reserved, busy)
        self.assertEqual(codes(listed), [pf.EXIT_BUSY])
        self.assertIn("is leased: lock held on /leases/testing-core7.lock", listed["refusals"][0]["message"])
        waited = pf.select_cores(topology, ("list", [6, 7]), reserved, busy, wait=True)
        self.assertEqual(waited["refusals"], [])
        self.assertIn("--wait queues for it", waited["notes"][0])
        reserved_list = pf.select_cores(topology, ("list", [3, 6]), reserved, busy, wait=True)
        self.assertEqual(codes(reserved_list), [pf.EXIT_BUSY])
        self.assertIn("CPU 3 (physical core 3) is reserved: static: elsewhere", reserved_list["refusals"][0]["message"])
        self.assertFalse(reserved_list["refusals"][0]["retryable"])

    def test_absent_cpu(self):
        result = pf.select_cores(zen_topology(), ("list", [6, 99]))
        self.assertEqual(codes(result), [pf.EXIT_ABSENT])
        self.assertIn("CPU(s) 99 not present", result["refusals"][0]["message"])

    def test_preference_order(self):
        order = pf.core_order(zen_topology(), [14, 7])
        self.assertEqual(order, [6, 7, 0, 1, 2, 3, 4, 5])


# --------------------------------------------------------------------------- GPUs


class GpuTests(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(pf.parse_gpu_request("none"), ("none", None))
        self.assertEqual(pf.parse_gpu_request("auto"), ("auto", None))
        self.assertEqual(pf.parse_gpu_request("0"), ("list", [0]))
        self.assertEqual(pf.parse_gpu_request("0,2"), ("list", [0, 2]))
        self.assertEqual(pf.parse_gpu_request(UUIDS[1]), ("list", [UUIDS[1]]))
        self.assertEqual(pf.parse_gpu_request(f"1,{UUIDS[3]}"), ("list", [1, UUIDS[3]]))
        self.assertEqual(pf.parse_gpu_request("count:2"), ("count", 2))
        for bad in ("", "x", "0,,1", "0,0", "-1", "gpu0"):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                pf.parse_gpu_request(bad)

    def test_auto(self):
        inventory = fake_gpus(4)
        two = pf.select_gpus(inventory, ("auto", None), nranks=2)
        self.assertEqual([g["index"] for g in two["gpus"]], [0, 1])
        self.assertIn("4 free GPU(s); 2 used by 2 rank(s)", two["notes"][0])
        skip_busy = pf.select_gpus(inventory, ("auto", None), nranks=2, busy={UUIDS[0]: "/l/gpu0"})
        self.assertEqual([g["index"] for g in skip_busy["gpus"]], [1, 2])
        skip_reserved = pf.select_gpus(inventory, ("auto", None), nranks=8, reserved={UUIDS[2]: ["r"]})
        self.assertEqual([g["index"] for g in skip_reserved["gpus"]], [0, 1, 3])
        none_free = pf.select_gpus(fake_gpus(1), ("auto", None), nranks=2, busy={UUIDS[0]: "/l/gpu"})
        self.assertEqual(none_free["gpus"], [])
        self.assertIn("ranks run CPU-only", none_free["notes"][0])
        no_gpu = pf.select_gpus([], ("auto", None), nranks=2, inventory_note="nvidia-smi not found")
        self.assertEqual(no_gpu["gpus"], [])
        self.assertIn("nvidia-smi not found", no_gpu["notes"][0])

    def test_none(self):
        result = pf.select_gpus(fake_gpus(4), ("none", None), nranks=4)
        self.assertEqual(result["gpus"], [])
        self.assertEqual(result["refusals"], [])

    def test_index_and_uuid(self):
        inventory = fake_gpus(4)
        by_index = pf.select_gpus(inventory, ("list", [1, 3]), nranks=2)
        self.assertEqual([g["uuid"] for g in by_index["gpus"]], [UUIDS[1], UUIDS[3]])
        by_uuid = pf.select_gpus(inventory, ("list", [UUIDS[2].lower()]), nranks=1)
        self.assertEqual([g["index"] for g in by_uuid["gpus"]], [2])
        mixed = pf.select_gpus(inventory, ("list", [UUIDS[3], 0]), nranks=2)
        self.assertEqual([g["index"] for g in mixed["gpus"]], [3, 0])

    def test_list_refusals(self):
        inventory = fake_gpus(2)
        absent = pf.select_gpus(inventory, ("list", [7]), nranks=1)
        self.assertEqual(codes(absent), [pf.EXIT_ABSENT])
        self.assertIn("GPU 7 not present", absent["refusals"][0]["message"])
        absent_uuid = pf.select_gpus(inventory, ("list", ["GPU-ffff"]), nranks=1)
        self.assertEqual(codes(absent_uuid), [pf.EXIT_ABSENT])
        too_many = pf.select_gpus(inventory, ("list", [0, 1]), nranks=1)
        self.assertEqual(codes(too_many), [pf.EXIT_USAGE])
        self.assertIn("use at most 1", too_many["refusals"][0]["message"])
        leased = pf.select_gpus(inventory, ("list", [0]), nranks=1, busy={UUIDS[0]: "/l/gpu.lock"})
        self.assertEqual(codes(leased), [pf.EXIT_BUSY])
        self.assertTrue(leased["refusals"][0]["retryable"])
        reserved = pf.select_gpus(inventory, ("list", [1]), nranks=1, reserved={UUIDS[1]: ["static: display"]})
        self.assertEqual(codes(reserved), [pf.EXIT_BUSY])
        waited = pf.select_gpus(inventory, ("list", [0]), nranks=1, busy={UUIDS[0]: "/l/gpu.lock"}, wait=True)
        self.assertEqual(waited["refusals"], [])

    def test_count(self):
        inventory = fake_gpus(4)
        self.assertEqual([g["index"] for g in pf.select_gpus(inventory, ("count", 2), nranks=2)["gpus"]], [0, 1])
        self.assertEqual(codes(pf.select_gpus(inventory, ("count", 5), nranks=8)), [pf.EXIT_ABSENT])
        self.assertEqual(codes(pf.select_gpus(inventory, ("count", 3), nranks=2)), [pf.EXIT_USAGE])
        busy = {UUIDS[0]: "a", UUIDS[1]: "b", UUIDS[2]: "c"}
        self.assertEqual(codes(pf.select_gpus(inventory, ("count", 2), nranks=2, busy=busy)), [pf.EXIT_BUSY])

    def test_rank_gpu_map(self):
        gpus = fake_gpus(2)
        self.assertEqual(pf.rank_gpu_map(5, gpus, 2), [UUIDS[0], UUIDS[0], UUIDS[1], UUIDS[1], None])
        self.assertEqual(pf.rank_gpu_map(3, gpus, 1), [UUIDS[0], UUIDS[1], None])
        self.assertEqual(pf.rank_gpu_map(2, [], 1), [None, None])
        self.assertEqual(pf.rank_gpu_map(1, gpus[:1], 1), [UUIDS[0]])

    def test_ranks_per_gpu_usable(self):
        result = pf.select_gpus(fake_gpus(4), ("auto", None), nranks=5, ranks_per_gpu=2)
        self.assertEqual(len(result["gpus"]), 3)  # ceil(5 / 2)

    def test_parse_nvidia_smi(self):
        text = ("0, GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9, NVIDIA GeForce RTX 3080 Ti, 8.6, 12288 MiB, 615.78.08\n"
                "1, GPU-11111111-2222-3333-4444-555555555555, NVIDIA RTX PRO 4000 Blackwell, 12.0, 24467 MiB, 615.78.08\n")
        gpus = pf.parse_gpu_csv(text)
        self.assertEqual([g["index"] for g in gpus], [0, 1])
        self.assertEqual(gpus[1]["compute_cap"], "12.0")
        self.assertEqual(gpus[0]["name"], "NVIDIA GeForce RTX 3080 Ti")
        with self.assertRaises(ValueError):
            pf.parse_gpu_csv("0, GPU-x, name\n")


# --------------------------------------------------------------------------- reservations, identity, leases


class ReservationProbeTests(unittest.TestCase):
    def probe(self, result, existing=True):
        system = zen_system(probes={SIMQ: result})
        if not existing:
            system.existing.discard("/var/simq/queue.sqlite3")
        return pf.run_reservation_probe(zen_site()["reservations"][0], system)

    def test_idle_and_active(self):
        self.assertEqual(self.probe(IDLE), (False, "idle"))
        paused = (0, "Queue is paused (a running job is allowed to finish).\nNo active jobs.\n", "")
        self.assertEqual(self.probe(paused)[0], False)
        reserved, detail = self.probe(ACTIVE)
        self.assertTrue(reserved)
        self.assertIn("14  running", detail)

    def test_errors_count_as_reserved(self):
        for result in ((1, "", "boom"), (0, "something else\n", ""), (0, "", ""),
                       subprocess.TimeoutExpired("simq", 5), OSError("no exec")):
            with self.subTest(result=result):
                reserved, detail = self.probe(result)
                self.assertTrue(reserved)
                self.assertIn("reserved", detail)
        reserved, detail = self.probe(IDLE, existing=False)
        self.assertTrue(reserved)
        self.assertIn("required path missing", detail)

    def test_evaluate(self):
        system = zen_system(probes={SIMQ: ACTIVE})
        result = pf.evaluate_reservations(zen_site(), zen_topology(), system, fake_gpus(1))
        self.assertEqual(sorted(result["cores"]), [0, 1, 2, 3, 4, 5])
        self.assertEqual(len(result["cores"][0]), 2)
        self.assertEqual(len(result["cores"][4]), 1)
        self.assertTrue(result["probes"][0]["reserved"])
        idle = pf.evaluate_reservations(zen_site(), zen_topology(), zen_system(), fake_gpus(1))
        self.assertEqual(len(idle["cores"][0]), 1)


class IdentityTests(unittest.TestCase):
    def test_strict(self):
        site = zen_site()
        observed = pf.observe_identity(zen_system(), site, fake_gpus(1))
        self.assertTrue(pf.admit_identity(site["identity"], observed)["admitted"])
        changed = dict(observed, kernel_release="6.2.0-other")
        result = pf.admit_identity(site["identity"], changed)
        self.assertFalse(result["admitted"])
        self.assertEqual([c["key"] for c in result["checks"] if not c["match"]], ["kernel_release"])
        for key, value in (("nvidia_srcversion", "OTHER"), ("gpu_uuids", [UUIDS[1]]),
                           ("frozen_root_manifest_sha256", "0" * 64)):
            with self.subTest(key=key):
                self.assertFalse(pf.admit_identity(site["identity"], dict(observed, **{key: value}))["admitted"])

    def test_record_only_and_order(self):
        observed = {"kernel_release": "x", "gpu_uuids": [UUIDS[1], UUIDS[0]]}
        record = {"policy": "record-only", "record": {"kernel_release": "y", "gpu_uuids": [UUIDS[0], UUIDS[1]]}}
        result = pf.admit_identity(record, observed)
        self.assertTrue(result["admitted"])
        self.assertEqual([c["match"] for c in result["checks"]], [True, False])  # gpu_uuids order-free
        with self.assertRaises(pf.PfError):
            pf.admit_identity({"policy": "lenient"}, observed)


class LeaseTests(unittest.TestCase):
    def test_names(self):
        site = zen_site()
        self.assertEqual(pf.core_lease_path(site, 6), "/leases/testing-core.lock")
        self.assertEqual(pf.core_lease_path(site, 7), "/leases/testing-core7.lock")
        self.assertEqual(pf.core_lease_path(site, 4), "/leases/testing-core4.lock")
        self.assertEqual(pf.gpu_lease_paths(site, fake_gpus(1)[0]),
                         [f"/leases/testing-gpu-{UUIDS[0]}.lock", "/leases/testing-gpu.lock"])
        self.assertIsNone(pf.core_lease_path({"leases": {}}, 3))

    def test_order_ascending(self):
        site = zen_site()
        order = pf.lease_order(site, [7, 6], list(reversed(fake_gpus(2))))
        self.assertEqual(order, ["/leases/testing-core.lock", "/leases/testing-core7.lock",
                                 f"/leases/testing-gpu-{UUIDS[0]}.lock", f"/leases/testing-gpu-{UUIDS[1]}.lock",
                                 "/leases/testing-gpu.lock"])

    def test_flock(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "sub", "core6.lock")
            self.assertFalse(pf.lease_probe_busy(path))  # missing file is free and not created
            self.assertFalse(os.path.exists(path))
            first = pf.LeaseSet()
            self.assertTrue(first.acquire(path, block=False))
            self.assertTrue(pf.lease_probe_busy(path))
            second = pf.LeaseSet()
            self.assertFalse(second.acquire(path, block=False))
            first.release()
            self.assertFalse(pf.lease_probe_busy(path))
            self.assertTrue(second.acquire(path, block=False))
            second.release()

    def test_acquire_plan_leases_releases_on_busy(self):
        with tempfile.TemporaryDirectory() as directory:
            a, b = os.path.join(directory, "a.lock"), os.path.join(directory, "b.lock")
            holder = pf.LeaseSet()
            holder.acquire(b, block=False)
            leases, busy = pf.acquire_plan_leases({"leases": [a, b]}, block=False, err=io.StringIO())
            self.assertIsNone(leases)
            self.assertEqual(busy, b)
            self.assertFalse(pf.lease_probe_busy(a))  # partial acquisition released
            holder.release()
            leases, busy = pf.acquire_plan_leases({"leases": [a, b]}, block=False, err=io.StringIO())
            self.assertEqual(leases.paths(), [a, b])
            leases.release()


# --------------------------------------------------------------------------- plans, launch argv, Slurm


class PlanTests(unittest.TestCase):
    def test_two_ranks_gpu0(self):
        plan = plan_for("--cores", "2", "--gpus", "0")
        self.assertEqual(plan["refusals"], [])
        self.assertEqual(plan["launcher"], "mpirun")
        self.assertEqual([p["cpus"] for p in plan["processes"]], [[6], [7]])
        self.assertEqual([p["cuda_visible_devices"] for p in plan["processes"]], [UUIDS[0], ""])
        self.assertEqual(plan["leases"], ["/leases/testing-core.lock", "/leases/testing-core7.lock",
                                          f"/leases/testing-gpu-{UUIDS[0]}.lock", "/leases/testing-gpu.lock"])
        self.assertEqual(plan["nice"]["effective"], 15)
        argv = pf.build_child_argv(plan, ["./prog", "-x"])
        self.assertEqual(argv[:7], ["/frozen/run-frozen", "mpirun", "--bind-to", "none", "-np", "2", "/usr/bin/python3"])
        self.assertEqual(argv[7:9], [pf.SELF_PATH, "--rank-exec"])
        rank_plan = json.loads(argv[9])
        self.assertEqual(rank_plan["cpus"], [6, 7])
        self.assertEqual(rank_plan["gpus"], [UUIDS[0], ""])
        self.assertEqual(argv[10:], ["--", "./prog", "-x"])

    def test_single_rank_exec(self):
        plan = plan_for("--cpus", "6", "--gpus", "none")
        self.assertEqual(plan["launcher"], "exec")
        self.assertEqual(pf.build_child_argv(plan, ["./prog"]), ["/frozen/run-frozen", "./prog"])
        env = pf.process_env(plan, plan["processes"][0])
        self.assertEqual(env["CUDA_VISIBLE_DEVICES"], "")
        self.assertEqual(env["PF_CORE"], "6")
        for key in pf.THREAD_VARS:
            self.assertEqual(env[key], "1")  # site OMP_NUM_THREADS=4 is overridden
        self.assertEqual(env["SITE_VAR"], "x")

    def test_environment(self):
        env = pf.build_base_env({"UNWANTED": "1", "KEEP": "2", "OMP_NUM_THREADS": "16"}, zen_site())
        self.assertNotIn("UNWANTED", env)
        self.assertEqual(env["KEEP"], "2")
        self.assertEqual(env["OMP_NUM_THREADS"], "1")
        self.assertEqual(env["CUDA_DEVICE_ORDER"], "PCI_BUS_ID")
        self.assertEqual(pf.env_digest(env), pf.env_digest(dict(reversed(list(env.items())))))

    def test_host_and_container_execution(self):
        host = plan_for("--cores", "2", "--gpus", "none", "--execution", "host")
        argv = pf.build_child_argv(host, ["./prog"])
        self.assertEqual(argv[:5], ["mpirun", "--bind-to", "none", "-np", "2"])
        self.assertEqual(argv[5], sys.executable)
        container = plan_for("--cores", "2", "--gpus", "none", "--prefix", "apptainer exec --nv img.sif")
        argv = pf.build_child_argv(container, ["./prog"])
        self.assertEqual(argv[0], "mpirun")
        self.assertEqual(argv[argv.index("--"):], ["--", "apptainer", "exec", "--nv", "img.sif", "./prog"])
        self.assertEqual(container["execution"]["kind"], "container")

    def test_group_mode(self):
        plan = plan_for("--cpus", "6-7", "--gpus", "none", "--launcher", "exec")
        self.assertTrue(plan["group"])
        self.assertEqual(len(plan["processes"]), 1)
        self.assertEqual(plan["processes"][0]["cpus"], [6, 7])
        env = pf.process_env(plan, plan["processes"][0])
        self.assertEqual(env["PF_CORE"], "6,7")
        self.assertEqual(env["MAKEFLAGS"], "-j2")
        self.assertEqual(env["CMAKE_BUILD_PARALLEL_LEVEL"], "2")

    def test_wait_count_is_retryable(self):
        system = zen_system(busy={"/leases/testing-core.lock", "/leases/testing-core7.lock"})
        plan = plan_for("--cores", "1", "--gpus", "none", "--wait", system=system)
        self.assertEqual(pf.plan_exit_code(plan), pf.EXIT_BUSY)
        self.assertTrue(all(r["retryable"] for r in plan["refusals"]))
        one_busy = zen_system(busy={"/leases/testing-core.lock"})
        plan = plan_for("--cores", "1", "--gpus", "none", system=one_busy)
        self.assertEqual([p["cpus"] for p in plan["processes"]], [[7]])
        self.assertEqual(plan["leases"], ["/leases/testing-core7.lock"])

    def test_nice_not_lowered(self):
        plan = plan_for("--cpus", "6", "--gpus", "none", system=zen_system(nice=19))
        self.assertEqual(plan["nice"]["effective"], 19)
        self.assertTrue(any("cannot be lowered" in n for n in plan["notes"]))

    def test_cgroup_restriction(self):
        plan = plan_for("--cpus", "6", "--gpus", "none", system=zen_system(cgroup={0, 1, 2, 3, 4, 5, 7, 8}))
        self.assertEqual(pf.plan_exit_code(plan), pf.EXIT_BUSY)
        self.assertIn("cgroup cpuset", plan["refusals"][0]["message"])

    def test_identity_refusal(self):
        system = zen_system(files={"/fake/sys/module/nvidia/srcversion": "CHANGED"})
        plan = plan_for("--cpus", "6", "--gpus", "none", system=system)
        self.assertEqual(pf.plan_exit_code(plan), pf.EXIT_CONFIG)
        self.assertIn("nvidia_srcversion: recorded 'ABC123', observed 'CHANGED'", plan["refusals"][0]["message"])

    def test_station_auto(self):
        site = zen_site(reserved_cores={"cores": [], "reason": ""}, reservations=[], core_preference=[],
                        identity={"policy": "record-only", "record": {}})
        system = zen_system(topology=station_topology(), gpus=fake_gpus(1))
        plan = plan_for("--gpus", "auto", site=site, system=system)
        self.assertEqual(len(plan["processes"]), 64)
        self.assertEqual(plan["processes"][0]["cuda_visible_devices"], UUIDS[0])
        self.assertTrue(all(p["cuda_visible_devices"] == "" for p in plan["processes"][1:]))
        plan5 = plan_for("--cores", "5", "--gpus", "0", "--ranks-per-gpu", "2", site=site, system=system)
        self.assertEqual([p["cuda_visible_devices"] for p in plan5["processes"]], [UUIDS[0], UUIDS[0], "", "", ""])


class SlurmTests(unittest.TestCase):
    ENV = {"PATH": "/usr/bin", "SLURM_JOB_ID": "4242", "SLURM_NTASKS": "4", "SLURM_JOB_NODELIST": "n[01-02]"}

    def test_srun_plan(self):
        system = zen_system(env=self.ENV)
        plan = plan_for("--gpus", "auto", system=system)
        self.assertEqual(plan["refusals"], [])
        self.assertEqual(plan["launcher"], "srun")
        self.assertEqual(len(plan["processes"]), 4)
        self.assertEqual(plan["leases"], [])
        self.assertEqual(system.calls, [])  # reservation probe not run inside an allocation
        argv = pf.build_child_argv(plan, ["./prog"])
        self.assertEqual(argv[:6], ["srun", "--ntasks=4", "--cpus-per-task=1", "--cpu-bind=cores",
                                    "--ntasks-per-core=1", "--kill-on-bad-exit=1"])
        rank_plan = json.loads(argv[argv.index("--rank-exec") + 1])
        self.assertEqual(rank_plan["mode"], "srun")
        self.assertEqual(rank_plan["gpu_request"], {"kind": "auto", "value": None})
        self.assertEqual(argv[-3:], ["--", "/frozen/run-frozen", "./prog"])
        self.assertEqual(plan["slurm"]["SLURM_JOB_ID"], "4242")

    def test_srun_counts(self):
        system = zen_system(env=self.ENV)
        self.assertEqual(len(plan_for("--cores", "3", system=system)["processes"]), 3)
        refused = plan_for("--cores", "8", system=zen_system(env=self.ENV))
        self.assertEqual(pf.plan_exit_code(refused), pf.EXIT_BUSY)
        self.assertIn("SLURM_NTASKS", refused["refusals"][0]["message"])
        shrunk = plan_for("--cores", "8", "--cores-available", system=zen_system(env=self.ENV))
        self.assertEqual(len(shrunk["processes"]), 4)
        listed = plan_for("--cpus", "6-7", system=zen_system(env=self.ENV))
        self.assertEqual(pf.plan_exit_code(listed), pf.EXIT_USAGE)
        self.assertIn("Slurm allocation", listed["refusals"][0]["message"])

    def test_srun_rank_assignment(self):
        base = {"schema": pf.RANK_PLAN_SCHEMA, "mode": "srun", "nranks": 4, "nice": 10, "ranks_per_gpu": 1,
                "gpu_request": {"kind": "auto", "value": None}}
        env = {"SLURM_PROCID": "3", "SLURM_LOCALID": "1", "CUDA_VISIBLE_DEVICES": "0,1"}
        a = pf.rank_assignment(base, env)
        self.assertEqual((a["rank"], a["local_rank"], a["cpu"], a["cuda_visible_devices"]), (3, 1, None, "1"))
        self.assertEqual(pf.rank_assignment(dict(base, ranks_per_gpu=2), env)["cuda_visible_devices"], "0")
        self.assertEqual(pf.rank_assignment(dict(base, gpu_request={"kind": "none", "value": None}), env)
                         ["cuda_visible_devices"], "")
        third = dict(env, SLURM_LOCALID="2")
        self.assertEqual(pf.rank_assignment(base, third)["cuda_visible_devices"], "")
        listed = dict(base, gpu_request={"kind": "list", "value": [1]})
        self.assertEqual(pf.rank_assignment(listed, dict(env, SLURM_LOCALID="0"))["cuda_visible_devices"], "1")
        self.assertEqual(pf.slurm_visible_gpus({"SLURM_STEP_GPUS": "2,3"}), ["2", "3"])

    def test_mpirun_rank_assignment(self):
        plan = {"schema": pf.RANK_PLAN_SCHEMA, "mode": "mpirun", "nranks": 2, "nice": 15, "cpus": [6, 7],
                "gpus": [UUIDS[0], ""], "ranks_per_gpu": 1}
        a = pf.rank_assignment(plan, {"OMPI_COMM_WORLD_RANK": "1", "OMPI_COMM_WORLD_LOCAL_RANK": "1"})
        self.assertEqual((a["cpu"], a["cuda_visible_devices"], a["env"]["PF_CORE"]), (7, "", "7"))
        self.assertEqual(pf.rank_assignment(plan, {"PMIX_RANK": "0"})["cuda_visible_devices"], UUIDS[0])
        with self.assertRaises(pf.PfError):
            pf.rank_assignment(plan, {})
        with self.assertRaises(pf.PfError):
            pf.rank_assignment(plan, {"OMPI_COMM_WORLD_RANK": "2"})
        single = dict(plan, nranks=1, cpus=[6], gpus=[""])
        self.assertEqual(pf.rank_assignment(single, {})["rank"], 0)


# --------------------------------------------------------------------------- CLI dry-run and receipts


class CliTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.site_path = os.path.join(self.tmp.name, "testhost.json")
        self.write_site(zen_site())

    def tearDown(self):
        self.tmp.cleanup()

    def write_site(self, site):
        with open(self.site_path, "w") as f:
            json.dump(site, f)

    def run_cli(self, *args, system=None):
        out, err = io.StringIO(), io.StringIO()
        code = pf.main(["--site", self.site_path] + list(args), system=system or zen_system(), out=out, err=err)
        return code, out.getvalue(), err.getvalue()

    def test_ok(self):
        code, out, _ = self.run_cli("--dry-run", "--cores", "2", "--gpus", "0", "--", "./prog")
        self.assertEqual(code, 0)
        self.assertIn("result     OK", out)
        self.assertIn("rank 0     cpu 6 (physical core 6; siblings 6,14)  CUDA_VISIBLE_DEVICES=" + UUIDS[0], out)
        self.assertIn("rank 1     cpu 7 (physical core 7; siblings 7,15)  CUDA_VISIBLE_DEVICES=\n", out)

    def test_refusal_messages(self):
        cases = [
            (("--cores", "8"), pf.EXIT_BUSY, "--cores 8: only 2 of 8 physical cores available (6-7)"),
            (("--cores", "3"), pf.EXIT_BUSY, "use --cores 2 or --cores-available to accept fewer"),
            (("--cores", "9"), pf.EXIT_ABSENT, "this host has 8 physical cores"),
            (("--cpus", "6,14"), pf.EXIT_USAGE, "SMT siblings of physical core 6"),
            (("--cpus", "2"), pf.EXIT_BUSY, "CPU 2 (physical core 2) is reserved: static: allocated elsewhere"),
            (("--cpus", "6", "--gpus", "3"), pf.EXIT_ABSENT, "GPU 3 not present"),
            (("--cpus", "6", "--gpus", "x"), pf.EXIT_USAGE, "--gpus: invalid element 'x'"),
            (("--cores", "zero"), pf.EXIT_USAGE, "--cores expects a positive integer or 'auto'"),
            (("--cpus", "6", "--ranks-per-gpu", "0"), pf.EXIT_USAGE, "--ranks-per-gpu must be at least 1"),
        ]
        for args, expected, text in cases:
            with self.subTest(args=args):
                code, out, _ = self.run_cli("--dry-run", *args, "--", "./prog")
                self.assertEqual(code, expected)
                self.assertIn("REFUSED", out)
                self.assertIn(text, out)

    def test_leased_refusals(self):
        system = zen_system(busy={"/leases/testing-gpu.lock"})
        code, out, _ = self.run_cli("--dry-run", "--cpus", "6", "--gpus", "0", "--", "./prog", system=system)
        self.assertEqual(code, pf.EXIT_BUSY)
        self.assertIn("is leased: lock held on /leases/testing-gpu.lock", out)
        code, out, _ = self.run_cli("--dry-run", "--cpus", "6", "--", "./prog", system=system)
        self.assertEqual(code, 0)  # auto skips the leased GPU
        self.assertIn("--gpus auto: no free GPU", out)
        system = zen_system(busy={"/leases/testing-core7.lock"})
        code, out, _ = self.run_cli("--dry-run", "--cpus", "6-7", "--gpus", "none", "--", "./prog", system=system)
        self.assertEqual(code, pf.EXIT_BUSY)
        self.assertIn("CPU 7 (physical core 7) is leased", out)

    def test_dry_run_wait(self):
        system = zen_system(busy={"/leases/testing-core7.lock"})
        receipt_path = os.path.join(self.tmp.name, "wait.json")
        code, out, _ = self.run_cli("--dry-run", "--cpus", "6-7", "--gpus", "none", "--wait", "--receipt",
                                    receipt_path, "--", "x", system=system)
        self.assertEqual(code, pf.EXIT_BUSY)
        self.assertIn("result     WAIT (exit 75): with --wait a real run blocks on /leases/testing-core7.lock", out)
        with open(receipt_path) as f:
            self.assertEqual(json.load(f)["exit_codes"]["status"], "dry-run-wait")
        system = zen_system(busy={"/leases/testing-core.lock", "/leases/testing-core7.lock"})
        code, out, _ = self.run_cli("--dry-run", "--cores", "1", "--gpus", "none", "--wait", "--", "x", system=system)
        self.assertEqual(code, pf.EXIT_BUSY)
        self.assertIn("--wait queues until the leased cores are released", out)
        self.assertIn("with --wait a real run re-plans every 3 s", out)
        code, out, _ = self.run_cli("--dry-run", "--cores", "2", "--gpus", "none", "--wait", "--", "x")
        self.assertEqual(code, 0)
        self.assertIn("result     OK", out)

    def test_cores_available(self):
        code, out, _ = self.run_cli("--dry-run", "--cores", "8", "--cores-available", "--gpus", "none", "--", "x")
        self.assertEqual(code, 0)
        self.assertIn("--cores 8: 2 physical core(s) available (6-7); running 2 rank(s) (--cores-available)", out)

    def test_identity_refusal(self):
        system = zen_system(files={"kernel": "6.9.9-new"})
        code, out, _ = self.run_cli("--dry-run", "--cpus", "6", "--gpus", "none", "--", "x", system=system)
        self.assertEqual(code, pf.EXIT_CONFIG)
        self.assertIn("identity admission (strict) failed: kernel_release", out)

    def test_missing_and_bad_site(self):
        out, err = io.StringIO(), io.StringIO()
        code = pf.main(["--site", os.path.join(self.tmp.name, "absent.json"), "--dry-run"], system=zen_system(),
                       out=out, err=err)
        self.assertEqual(code, pf.EXIT_CONFIG)
        self.assertIn("pf-launch --detect", err.getvalue())
        self.write_site({"schema": "other"})
        code, _, err_text = self.run_cli("--dry-run")
        self.assertEqual(code, pf.EXIT_CONFIG)
        self.assertIn("schema", err_text)

    def test_json_receipt_schema(self):
        receipt_path = os.path.join(self.tmp.name, "r.json")
        code, out, _ = self.run_cli("--dry-run", "--json", "--cores", "2", "--gpus", "0", "--receipt", receipt_path,
                                    "--", "./prog")
        self.assertEqual(code, 0)
        receipt = json.loads(out)
        with open(receipt_path) as f:
            self.assertEqual(json.load(f), receipt)
        self.assertEqual(tuple(receipt), pf.RECEIPT_FIELDS)
        self.assertEqual(receipt["rank_count"], 2)
        self.assertEqual([r["cpus"] for r in receipt["rank_to_core"]], [[6], [7]])
        self.assertEqual([r["gpus"] for r in receipt["rank_to_gpu"]], [[UUIDS[0]], []])
        self.assertEqual(receipt["gpus"][0]["uuid"], UUIDS[0])
        self.assertEqual(set(receipt["gpus"][0]), {"index", "uuid", "name", "compute_cap", "memory_total",
                                                   "driver_version", "pci_bus_id", "numa"})
        self.assertEqual(receipt["site_file_digest"]["sha256"], pf.sha256_file(self.site_path))
        self.assertEqual(receipt["exit_codes"]["status"], "dry-run")
        self.assertTrue(receipt["launcher_argv"]["dry_run"])
        self.assertEqual(receipt["root_or_image_digest"]["kind"], "frozen-root")
        self.assertEqual(receipt["root_or_image_digest"]["manifest_sha256"], "f" * 64)
        self.assertEqual(len(receipt["environment_digest"]["sha256"]), 64)
        self.assertEqual(receipt["selected_cores"]["leases"],
                         ["/leases/testing-core.lock", "/leases/testing-core7.lock"])
        self.assertEqual(receipt["selected_gpus"]["leases"],
                         [f"/leases/testing-gpu-{UUIDS[0]}.lock", "/leases/testing-gpu.lock"])
        self.assertTrue(receipt["host"]["identity"]["admitted"])

    def test_refused_receipt_schema(self):
        receipt_path = os.path.join(self.tmp.name, "refused.json")
        code, _, _ = self.run_cli("--dry-run", "--cores", "8", "--receipt", receipt_path, "--", "x")
        self.assertEqual(code, pf.EXIT_BUSY)
        with open(receipt_path) as f:
            receipt = json.load(f)
        self.assertEqual(tuple(receipt), pf.RECEIPT_FIELDS)
        self.assertEqual(receipt["exit_codes"]["status"], "dry-run-refused")
        self.assertEqual(receipt["exit_codes"]["pf_launch"], pf.EXIT_BUSY)
        self.assertIn("--cores 8", receipt["exit_codes"]["messages"][0])
        self.assertEqual(receipt["rank_count"], 0)

    def test_receipt_with_rank_reports(self):
        plan = plan_for("--cores", "2", "--gpus", "0")
        reports = {0: {"host": "testhost", "affinity": [6], "nice": 15, "cuda_visible_devices": UUIDS[0], "exit_code": 0},
                   1: {"host": "testhost", "affinity": [7], "nice": 15, "cuda_visible_devices": "", "exit_code": 3}}
        receipt = pf.build_receipt(plan, status="completed", pf_code=3, child_argv=["mpirun"], start_time="s",
                                   end_time="e", launcher_rc=3, rank_reports=reports, leases_held=plan["leases"],
                                   system=zen_system())
        self.assertEqual(tuple(receipt), pf.RECEIPT_FIELDS)
        self.assertEqual(receipt["exit_codes"]["ranks"], [{"rank": 0, "host": "testhost", "exit_code": 0},
                                                          {"rank": 1, "host": "testhost", "exit_code": 3}])
        self.assertEqual(receipt["rank_to_core"][1]["observed"]["affinity"], [7])
        self.assertEqual(receipt["rank_to_gpu"][0]["observed_cuda_visible_devices"], UUIDS[0])
        self.assertEqual(receipt["selected_cores"]["leases_held"], plan["core_leases"])
        json.dumps(receipt)
        exec_plan = plan_for("--cpus", "6", "--gpus", "none")
        exec_receipt = pf.build_receipt(exec_plan, status="completed", pf_code=0, child_argv=["x"], start_time="s",
                                        end_time="e", launcher_rc=0, system=zen_system())
        self.assertEqual(exec_receipt["exit_codes"]["ranks"], [{"rank": 0, "host": "testhost", "exit_code": 0}])

    def test_slurm_receipt_schema(self):
        system = zen_system(env=SlurmTests.ENV)
        plan = plan_for("--cores", "2", system=system)
        receipt = pf.build_receipt(plan, status="dry-run", pf_code=0, child_argv=None, start_time="s",
                                   end_time=None, system=system)
        self.assertEqual(tuple(receipt), pf.RECEIPT_FIELDS)
        self.assertEqual(receipt["rank_to_core"][0]["binding"], "srun --cpu-bind (allocation)")
        self.assertEqual(receipt["selected_cores"]["slurm"]["SLURM_JOB_ID"], "4242")

    def test_overlap_check(self):
        reports = {0: {"host": "testhost", "affinity": [6, 14]}, 1: {"host": "testhost", "affinity": [14]},
                   2: {"host": "other", "affinity": [6]}}
        self.assertEqual(pf.check_rank_overlap(reports, zen_topology(), "testhost"),
                         ["ranks 0 and 1 share physical core 6"])

    def test_detect(self):
        path = os.path.join(self.tmp.name, "sites", "detected.json")
        system = zen_system(topology=station_topology(), gpus=fake_gpus(2))
        out, err = io.StringIO(), io.StringIO()
        self.assertEqual(pf.main(["--detect", "--site", path], system=system, out=out, err=err), 0)
        with open(path) as f:
            site = json.load(f)
        self.assertEqual(site["schema"], pf.SITE_SCHEMA)
        self.assertEqual(site["topology"]["physical_core_count"], 64)
        self.assertEqual(site["topology"]["numa_nodes"], {"0": "0-31,64-95", "1": "32-63,96-127"})
        self.assertEqual([g["uuid"] for g in site["gpus"]], UUIDS[:2])
        self.assertEqual(site["identity"]["policy"], "record-only")
        self.assertEqual(site["identity"]["record"]["gpu_uuids"], UUIDS[:2])
        self.assertEqual(pf.main(["--detect", "--site", path], system=system, out=out, err=err), pf.EXIT_CONFIG)
        self.assertIn("--force", err.getvalue())
        self.assertEqual(pf.main(["--detect", "--site", path, "--force"], system=system, out=out, err=err), 0)
        loaded, info = pf.load_site(path)  # a detected site file loads and plans without edits
        plan = pf.build_plan(options("--cores", "3"), loaded, info, system)
        self.assertEqual(plan["refusals"], [])
        self.assertEqual([p["cpus"] for p in plan["processes"]], [[0], [1], [2]])


if __name__ == "__main__":
    unittest.main()
