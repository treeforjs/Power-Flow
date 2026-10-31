"""pf-launch 1.1 inside a SimQueue allocation (SIMQ_JOB_ID + SIMQ_CPU_LIST / SIMQ_GPUS).

Like a Slurm allocation, the queue owns the placement: pf-launch picks no
cores of its own, evaluates no site reservations and takes no leases (the
queue's job unit already holds them).  Fake topology, GPUs and leases; no
program is launched.
"""

import io
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import pf_launch as pf  # noqa: E402
from test_pf_launch import UUIDS, codes, plan_for, zen_site, zen_system  # noqa: E402

POWERFLOW = {"PATH": "/usr/bin", "SIMQ_JOB_ID": "12", "SIMQ_PARTITION": "powerflow", "SIMQ_CPU_LIST": "6-7",
             "SIMQ_CPU_COUNT": "2", "SIMQ_NCPUS": "2", "SIMQ_GPUS": UUIDS[0], "SIMQ_GPU_MODE": "exclusive",
             "SIMQ_UNIT": "simq-job-12-abc.service", "CUDA_VISIBLE_DEVICES": UUIDS[0]}
ALL_LEASES = {"/leases/testing-core.lock", "/leases/testing-core7.lock", "/leases/testing-gpu.lock",
              f"/leases/testing-gpu-{UUIDS[0]}.lock"}


def env(**changes):
    value = dict(POWERFLOW)
    for key, item in changes.items():
        if item is None:
            value.pop(key, None)
        else:
            value[key] = item
    return value


class DetectionTests(unittest.TestCase):
    def test_allocation_detection(self):
        self.assertIsNone(pf.simqueue_allocation({"PATH": "/usr/bin"}))
        self.assertIsNone(pf.simqueue_allocation({"SIMQ_JOB_ID": "3"}))
        alloc = pf.simqueue_allocation(POWERFLOW)
        self.assertEqual((alloc["cpus"], alloc["gpus"], alloc["gpus_managed"]), ([6, 7], [UUIDS[0]], True))
        legacy = pf.simqueue_allocation({"SIMQ_JOB_ID": "3", "SIMQ_CPU_LIST": "0-3", "SIMQ_CPU_COUNT": "4"})
        self.assertEqual((legacy["cpus"], legacy["gpus_managed"]), ([0, 1, 2, 3], False))
        empty = pf.simqueue_allocation(env(SIMQ_GPUS=""))
        self.assertEqual((empty["gpus"], empty["gpus_managed"]), ([], True))


class AllocationPlanTests(unittest.TestCase):
    def test_uses_allocation_takes_no_leases_runs_no_probe(self):
        system = zen_system(env=POWERFLOW, busy=ALL_LEASES)  # the job unit holds every lease
        plan = plan_for("--gpus", "auto", system=system)
        self.assertEqual(plan["refusals"], [])
        self.assertEqual(plan["launcher"], "mpirun")
        self.assertEqual([p["cpus"] for p in plan["processes"]], [[6], [7]])
        self.assertEqual([p["cuda_visible_devices"] for p in plan["processes"]], [UUIDS[0], ""])
        self.assertEqual(plan["leases"], [])
        self.assertEqual(plan["waits"], [])
        self.assertEqual(system.calls, [])  # the queue-active reservation probe is not run
        self.assertEqual(plan["simqueue"]["SIMQ_JOB_ID"], "12")
        self.assertIsNone(plan["slurm"])

    def test_allocation_outside_site_core_preference_is_used_as_given(self):
        # The site reserves cores 0-5 for other work; a queue allocation there is authoritative.
        system = zen_system(env=env(SIMQ_PARTITION="spheral", SIMQ_CPU_LIST="0-3", SIMQ_GPUS=""))
        plan = plan_for("--gpus", "none", system=system)
        self.assertEqual(plan["refusals"], [])
        self.assertEqual([r["cpu"] for r in plan["ranks"]], [0, 1, 2, 3])
        self.assertEqual(plan["leases"], [])

    def test_core_requests_within_the_allocation(self):
        self.assertEqual([r["cpu"] for r in plan_for("--cores", "1", system=zen_system(env=POWERFLOW))["ranks"]], [6])
        too_many = plan_for("--cores", "3", system=zen_system(env=POWERFLOW))
        self.assertEqual(pf.plan_exit_code(too_many), pf.EXIT_BUSY)
        self.assertIn("exceeds SimQueue job 12's allocation of 2 core(s)", too_many["refusals"][0]["message"])
        shrunk = plan_for("--cores", "3", "--cores-available", system=zen_system(env=POWERFLOW))
        self.assertEqual(len(shrunk["ranks"]), 2)
        listed = plan_for("--cpus", "7", system=zen_system(env=POWERFLOW))
        self.assertEqual([r["cpu"] for r in listed["ranks"]], [7])
        outside = plan_for("--cpus", "5", system=zen_system(env=POWERFLOW))
        self.assertEqual(codes(outside), [pf.EXIT_BUSY])
        self.assertIn("not in SimQueue job 12's allocation", outside["refusals"][0]["message"])
        siblings = plan_for("--cpus", "6,14", system=zen_system(env=env(SIMQ_CPU_LIST="6,14")))
        self.assertIn(pf.EXIT_USAGE, codes(siblings))

    def test_gpu_requests_within_the_allocation(self):
        cpu_only = plan_for("--gpus", "auto", system=zen_system(env=env(SIMQ_GPUS="")))
        self.assertEqual(cpu_only["refusals"], [])
        self.assertEqual([p["cuda_visible_devices"] for p in cpu_only["processes"]], ["", ""])
        self.assertTrue(any("no GPU allocated" in n for n in cpu_only["notes"]))
        not_ours = plan_for("--gpus", "0", system=zen_system(env=env(SIMQ_GPUS="")))
        self.assertEqual(codes(not_ours), [pf.EXIT_BUSY])
        self.assertIn("is not in SimQueue job 12's allocation", not_ours["refusals"][0]["message"])
        by_index = plan_for("--cpus", "6", "--gpus", "0", system=zen_system(env=POWERFLOW, busy=ALL_LEASES))
        self.assertEqual(by_index["refusals"], [])
        self.assertEqual(by_index["processes"][0]["cuda_visible_devices"], UUIDS[0])
        counted = plan_for("--gpus", "count:2", system=zen_system(env=POWERFLOW, gpus=None))
        self.assertEqual(codes(counted), [pf.EXIT_BUSY])
        none = plan_for("--gpus", "none", system=zen_system(env=POWERFLOW))
        self.assertEqual([p["cuda_visible_devices"] for p in none["processes"]], ["", ""])

    def test_slurm_takes_precedence(self):
        system = zen_system(env=env(SLURM_JOB_ID="77", SLURM_NTASKS="2"))
        plan = plan_for("--gpus", "auto", system=system)
        self.assertEqual(plan["launcher"], "srun")
        self.assertIsNone(plan["simqueue"])

    def test_05_queue_job_cores_from_allocation_gpus_leased_locally(self):
        legacy = {"PATH": "/usr/bin", "SIMQ_JOB_ID": "9", "SIMQ_CPU_LIST": "6-7", "SIMQ_CPU_COUNT": "2"}
        plan = plan_for("--gpus", "auto", system=zen_system(env=legacy))
        self.assertEqual(plan["refusals"], [])
        self.assertEqual([r["cpu"] for r in plan["ranks"]], [6, 7])
        self.assertEqual(plan["core_leases"], [])
        self.assertEqual(plan["gpu_leases"], [f"/leases/testing-gpu-{UUIDS[0]}.lock", "/leases/testing-gpu.lock"])

    def test_srun_launcher_is_refused_inside_a_simqueue_allocation(self):
        plan = plan_for("--launcher", "srun", system=zen_system(env=POWERFLOW))
        self.assertIsNone(plan["simqueue"])  # explicit srun asks Slurm for an allocation (unchanged behaviour)


class AllocationCliTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.site_path = os.path.join(self.tmp.name, "testhost.json")
        with open(self.site_path, "w") as handle:
            json.dump(zen_site(), handle)

    def tearDown(self):
        self.tmp.cleanup()

    def test_dry_run_text_and_receipt(self):
        out, err = io.StringIO(), io.StringIO()
        receipt = os.path.join(self.tmp.name, "receipt.json")
        code = pf.main(["--site", self.site_path, "--dry-run", "--receipt", receipt, "--gpus", "auto", "--", "./prog"],
                       system=zen_system(env=POWERFLOW, busy=ALL_LEASES), out=out, err=err)
        self.assertEqual(code, 0, out.getvalue())
        self.assertIn("SimQueue job 12 partition powerflow: CPUs 6-7", out.getvalue())
        self.assertIn("result     OK", out.getvalue())
        with open(receipt) as handle:
            data = json.load(handle)
        self.assertEqual(data["selected_cores"]["simqueue"]["SIMQ_JOB_ID"], "12")
        self.assertEqual(data["selected_cores"]["leases"], [])
        self.assertEqual(data["selected_gpus"]["leases"], [])
        self.assertEqual(data["selected_gpus"]["uuids"], [UUIDS[0]])
        self.assertEqual(sorted(data), sorted(pf.RECEIPT_FIELDS))


if __name__ == "__main__":
    unittest.main()
