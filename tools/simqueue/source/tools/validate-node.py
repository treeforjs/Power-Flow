#!/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen
"""Real-configuration validation of the 0.6 node scheduler on zenobia.

Uses the production node config (real partitions and the real lease files)
with a scratch state directory and a distinct unit prefix, so it never
touches the live queue.  Jobs: powerflow probes on cores 6-7 (whichever are
not leased elsewhere), one powerflow GPU canary that only runs
`nvidia-smi -L` while holding the GPU leases, and a spheral affinity no-op on
cores 0-3.  The worker is killed with SIGKILL while jobs run and restarted
from a relocated copy of the code; the receipt records adoption, results,
lease holders and per-job CPU intervals.

Usage: validate-node.py OUTPUT_DIR
"""

from __future__ import annotations

import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[1]
RUNTIME = SOURCE.parent
PYTHON = "/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen"
CONTROLLER_CPUS = {8, 9, 10, 11}  # SMT siblings of the spheral cores 0-3
PROBE = r'''
import json, os, sys, time
path, seconds, code = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
r = {"job": os.environ["SIMQ_JOB_ID"], "partition": os.environ["SIMQ_PARTITION"],
     "cpu_list": os.environ["SIMQ_CPU_LIST"], "cpu_count": os.environ["SIMQ_CPU_COUNT"],
     "gpus": os.environ["SIMQ_GPUS"], "cuda": os.environ.get("CUDA_VISIBLE_DEVICES"),
     "affinity": sorted(os.sched_getaffinity(0)), "start": time.time(),
     "cgroup": open("/proc/self/cgroup").read().strip()}
time.sleep(seconds)
r["end"] = time.time()
open(path, "a").write(json.dumps(r) + "\n")
print("PROBE " + json.dumps(r), flush=True)
sys.exit(code)
'''


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    os.sched_setaffinity(0, CONTROLLER_CPUS)
    os.nice(15)
    projects = out / "projects"
    state = out / "state"
    projects.mkdir()
    config = json.loads((RUNTIME / "config" / "zenobia-node.json").read_text())
    config["launch"]["unit_prefix"] = "simqv"
    config["launch"]["slice"] = "simqv.slice"
    config["launch"]["poll_seconds"] = 0.5
    config_path = out / "node.json"
    config_path.write_text(json.dumps(config, indent=2))
    probe = projects / "probe.py"
    probe.write_text(PROBE)
    intervals = projects / "intervals.jsonl"
    env = {k: v for k, v in os.environ.items() if not k.startswith(("SIMQUEUE_", "SIMQ_", "SLURM_"))}
    env.update({"SIMQUEUE_STATE_DIR": str(state), "SIMQUEUE_PROJECT_ROOT": str(projects),
                "SIMQUEUE_CONFIG": str(config_path), "SIMQUEUE_NO_AUTOSTART": "1"})
    receipt: dict = {"schema": "simqueue-node-validation-v1", "runtime": str(RUNTIME), "events": []}

    def event(text: str, **data) -> None:
        receipt["events"].append({"t": round(time.time(), 3), "event": text, **data})
        print(f"[{time.strftime('%H:%M:%S')}] {text} {json.dumps(data) if data else ''}", flush=True)

    def simq(*args: str, source: Path = SOURCE) -> subprocess.CompletedProcess[str]:
        return subprocess.run([PYTHON, str(source / "simq"), *args], env=env, cwd=projects, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)

    def submit(*args: str) -> int:
        result = simq("submit", "--json", *args)
        if result.returncode != 0:
            raise RuntimeError(result.stderr)
        return int(json.loads(result.stdout)["job_id"])

    def jobs() -> dict[int, dict]:
        import sqlite3
        connection = sqlite3.connect(f"file:{state / 'queue.sqlite3'}?mode=ro", uri=True)
        connection.row_factory = sqlite3.Row
        rows = {int(r["id"]): dict(r) for r in connection.execute("SELECT * FROM jobs")}
        connection.close()
        return rows

    def lease_holders() -> dict[str, list]:
        sys.path.insert(0, str(SOURCE))
        from simqueue.locks import lock_holders
        paths = [config["lock_files"]["core_overrides"]["6"],
                 config["lock_files"]["core_path"].replace("{core}", "7"),
                 *(config["lock_files"]["core_path"].replace("{core}", str(c)) for c in range(4))]
        return {p: [(h.mode, h.pid) for h in hs] for p, hs in lock_holders(paths).items()}

    def start_worker(source: Path, *extra: str) -> subprocess.Popen:
        log = open(out / f"worker-{source.name}.out", "ab")
        return subprocess.Popen([PYTHON, str(source / "simq"), "worker", *extra], env=env, cwd=projects,
                                stdout=log, stderr=subprocess.STDOUT, start_new_session=True)

    event("leases before start", holders=lease_holders())
    pf_long = submit("-n", "pf-long", "-p", "powerflow", "-c", "1", "--mem", "256M", "-t", "10:00",
                     "--", PYTHON, str(probe), str(intervals), "15", "0")
    pf_fail = submit("-n", "pf-exit3", "-p", "powerflow", "-c", "1", "--mem", "256M", "-t", "10:00",
                     "--", PYTHON, str(probe), str(intervals), "3", "3")
    gpu = submit("-n", "pf-gpu-canary", "-p", "powerflow", "-c", "1", "--gpus", "1", "--mem", "256M",
                 "-t", "5:00", "--", "/bin/sh", "-c",
                 "echo CANARY cpu=$SIMQ_CPU_LIST gpus=$SIMQ_GPUS cuda=$CUDA_VISIBLE_DEVICES; nvidia-smi -L")
    spheral = submit("-n", "spheral-noop", "-p", "spheral", "-c", "4", "--mem", "256M", "-t", "5:00",
                     "--", PYTHON, str(probe), str(intervals), "3", "0")
    event("submitted", pf_long=pf_long, pf_fail=pf_fail, gpu=gpu, spheral=spheral)
    receipt["test_only_raikou_mfem"] = json.loads(simq("submit", "--test-only", "-p", "raikou-mfem", "--gpus", "1",
                                                       "--", "/bin/true").stdout)
    worker = start_worker(SOURCE)
    deadline = time.time() + 3 * 3600
    while time.time() < deadline:
        rows = jobs()
        if rows[pf_long]["state"] == "running" and rows[spheral]["state"] == "running":
            break
        time.sleep(0.5)
    else:
        raise RuntimeError("jobs did not start")
    event("running", pf_long_cpus=jobs()[pf_long]["cpu_list"], spheral_cpus=jobs()[spheral]["cpu_list"],
          holders=lease_holders(), partitions=simq("partitions").stdout)
    os.kill(worker.pid, signal.SIGKILL)
    worker.wait()
    event("worker killed with SIGKILL", pid=worker.pid)
    units = {j: jobs()[j]["unit_name"] for j in (pf_long, spheral)}
    active = {j: subprocess.run(["systemctl", "--user", "is-active", u], stdout=subprocess.PIPE, text=True).stdout.strip()
              for j, u in units.items()}
    event("units after worker death", active=active)
    time.sleep(6)
    status = json.loads((Path(jobs()[spheral]["job_dir"]) / "status.json").read_text())
    event("spheral job finished while worker down", phase=status["phase"], db_state=jobs()[spheral]["state"],
          holders=lease_holders())
    relocated = out / "relocated-runtime"
    shutil.copytree(SOURCE / "simqueue", relocated / "simqueue", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copy2(SOURCE / "simq", relocated / "simq")
    init = relocated / "simqueue" / "__init__.py"
    init.write_text(init.read_text().replace('"0.6.0"', '"0.6.0+relocated"'))
    restarted = start_worker(relocated, "--drain")
    restarted.wait(timeout=3 * 3600)
    event("relocated worker drained", returncode=restarted.returncode)
    rows = jobs()
    receipt["jobs"] = {
        str(j): {k: rows[j][k] for k in ("name", "state", "exit_code", "partition", "cpu_list", "gpu_list",
                                        "unit_name", "result", "adopted_at", "started_at", "finished_at", "message")}
        for j in (pf_long, pf_fail, gpu, spheral)
    }
    receipt["intervals"] = [json.loads(l) for l in intervals.read_text().splitlines()]
    receipt["gpu_canary_log"] = simq("log", "-n", "30", str(gpu)).stdout
    receipt["worker_log"] = (state / "worker.log").read_text()
    receipt["leases_after"] = lease_holders()
    checks = {
        "pf_long_succeeded_and_adopted": rows[pf_long]["state"] == "succeeded" and bool(rows[pf_long]["adopted_at"]),
        "spheral_finished_while_down_reconciled": rows[spheral]["state"] == "succeeded",
        "spheral_affinity_0_3": any(i["job"] == str(spheral) and i["affinity"] == [0, 1, 2, 3]
                                    for i in receipt["intervals"]),
        "pf_exit3_failed_with_3": rows[pf_fail]["state"] == "failed" and rows[pf_fail]["exit_code"] == 3,
        "gpu_canary_succeeded": rows[gpu]["state"] == "succeeded" and "RTX 3080 Ti" in receipt["gpu_canary_log"],
        "relocated_code_adopted": "0.6.0+relocated" in receipt["worker_log"] and f"adopted job {pf_long}" in receipt["worker_log"],
    }
    by_job = {i["job"]: i for i in receipt["intervals"]}
    overlap = []
    items = list(by_job.values())
    for a in range(len(items)):
        for b in range(a + 1, len(items)):
            shared = set(items[a]["affinity"]) & set(items[b]["affinity"])
            if shared and min(items[a]["end"], items[b]["end"]) > max(items[a]["start"], items[b]["start"]):
                overlap.append((items[a]["job"], items[b]["job"], sorted(shared)))
    checks["no_cpu_double_allocation"] = not overlap
    receipt["checks"] = checks
    receipt["passed"] = all(checks.values())
    (out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    event("done", checks=checks)
    return 0 if receipt["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
