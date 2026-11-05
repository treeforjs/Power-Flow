#!/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen
"""Measure the worker's own CPU time and memory (SimQueue 0.6.3 against 0.7.0).

Usage: bench-overhead.py --out DIR [--minutes M] [--cpus LIST] [--runtimes R063,R070]

For each runtime and each scenario a sandbox worker runs as a transient
systemd user unit (simqbench-<version>-<scenario>.service, with the unit
properties that version is installed with), so its cgroup's cpu.stat
(usage_usec, children included) and memory.current/memory.peak are exact:

idle  empty queue
busy  a pool of the given CPUs (default 2), 12 one-core ``sleep`` jobs: two
      running, ten queued behind them (the jobs sleep: they use no CPU)

The sandbox node config mirrors zenobia's (lease files, a GPU with the GPU
line, borrowing on, poll_seconds 1) with private paths.  All four workers
run at once for the same window after a warm-up.  Receipt:
DIR/bench-overhead.json.
"""

from __future__ import annotations

import argparse
import json
import os
import secrets
import shutil
import subprocess
import sys
import time
from pathlib import Path

PYTHON = "/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen"
RUNTIMES = {
    "0.6.3": Path("/home/treefor/projects/runtimes/simqueue-0.6.3-20261105"),
    "0.7.0": Path(__file__).resolve().parents[2],
}
UNIT_PROPERTIES = {
    "0.6.3": ["Nice=5"],
    "0.7.0": ["Nice=10", "CPUWeight=20", "IOSchedulingClass=idle", "MemoryHigh=96M", "MemoryMax=128M"],
}
FAKE_GPU = "GPU-bench-0000-0000-0000-000000000001"
USER_CGROUP = Path(f"/sys/fs/cgroup/user.slice/user-{os.getuid()}.slice/user@{os.getuid()}.service")


def node_config(base: Path, prefix: str, cpus: str) -> dict:
    locks = base / "locks"
    return {
        "schema": "simqueue-node-config-v1",
        "node": "bench",
        "default_partition": "pool",
        "memory_schedulable": "28G",
        "allow_smt_siblings": True,
        "gpus": [{"uuid": FAKE_GPU, "index": 0, "name": "fake", "memory": "12G", "budget": "11G"}],
        "lock_files": {"enabled": True, "core_path": str(locks / "core{core}.lock"),
                       "gpu_path": str(locks / "gpu-{uuid}.lock"), "gpu_extra": {FAKE_GPU: [str(locks / "gpu.lock")]},
                       "borrow_path": str(locks / "borrow-core{core}.lock")},
        "launch": {"unit_prefix": prefix, "slice": "simqbench.slice", "kill_grace_seconds": 2, "poll_seconds": 1},
        "borrow": {"enabled": True, "slice": "simqbenchborrow.slice"},
        "gpu_line": {"enabled": True, "state_path": str(locks / "simq-gpu-line.json"),
                     "ticket_dir": str(locks / "gpu-line.d"), "log_path": str(locks / "simq-gpu-line.log")},
        "partitions": [{"name": "pool", "cpus": cpus, "gpus": [FAKE_GPU], "default_cpus": 1, "default_mem": "256M",
                        "max_mem": "2G", "default_time": "2-00:00:00", "max_time": "14-00:00:00", "nice": 15}],
    }


class Sandbox:
    def __init__(self, root: Path, version: str, scenario: str, cpus: str) -> None:
        self.version = version
        self.scenario = scenario
        self.runtime = RUNTIMES[version]
        self.base = root / f"{version}-{scenario}"
        if self.base.exists():
            shutil.rmtree(self.base)
        (self.base / "projects").mkdir(parents=True)
        (self.base / "locks" / "gpu-line.d").mkdir(parents=True)
        self.prefix = "simqb" + secrets.token_hex(3)
        self.unit = f"simqbench-{version.replace('.', '')}-{scenario}-{secrets.token_hex(2)}.service"
        self.config = self.base / "node.json"
        self.config.write_text(json.dumps(node_config(self.base, self.prefix, cpus), indent=2))
        self.env = {
            "SIMQUEUE_STATE_DIR": str(self.base / "state"),
            "SIMQUEUE_PROJECT_ROOT": str(self.base / "projects"),
            "SIMQUEUE_CONFIG": str(self.config),
            "SIMQUEUE_NO_AUTOSTART": "1",
        }

    def cli(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        env = dict(os.environ, **self.env)
        return subprocess.run([PYTHON, str(self.runtime / "source" / "simq"), *arguments], env=env,
                              cwd=self.base / "projects", text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              check=False, timeout=120)

    def start(self) -> None:
        command = ["systemd-run", "--user", f"--unit={self.unit}", "--collect", "--quiet",
                   f"--description=SimQueue overhead bench {self.version} {self.scenario}",
                   f"--working-directory={self.base / 'projects'}"]
        command += [f"--property={item}" for item in UNIT_PROPERTIES[self.version]]
        command += [f"--setenv={key}={value}" for key, value in self.env.items()]
        command += ["--", PYTHON, str(self.runtime / "source" / "simq"), "worker"]
        subprocess.run(command, check=True)
        if self.scenario == "busy":
            for index in range(12):
                result = self.cli("submit", "-p", "pool", "-c", "1", "-n", f"sleeper-{index}", "--", "/bin/sleep",
                                  "7200")
                if result.returncode != 0:
                    raise RuntimeError(result.stderr)

    @property
    def cgroup(self) -> Path:
        return USER_CGROUP / "app.slice" / self.unit

    def read(self) -> dict:
        def number(name: str) -> int | None:
            try:
                return int((self.cgroup / name).read_text().split()[0])
            except (OSError, ValueError):
                return None

        usage = None
        try:
            for line in (self.cgroup / "cpu.stat").read_text().splitlines():
                if line.startswith("usage_usec "):
                    usage = int(line.split()[1])
        except OSError:
            pass
        pids = []
        try:
            pids = [int(p) for p in (self.cgroup / "cgroup.procs").read_text().split()]
        except OSError:
            pass
        rss = 0
        for pid in pids:
            try:
                for line in Path(f"/proc/{pid}/status").read_text().splitlines():
                    if line.startswith("VmRSS:"):
                        rss += int(line.split()[1]) * 1024
            except OSError:
                pass
        return {"cpu_usec": usage, "memory_current": number("memory.current"), "memory_peak": number("memory.peak"),
                "rss": rss, "unix": time.time()}

    def stop(self) -> dict:
        state = self.cli("list", "--all").stdout
        subprocess.run(["systemctl", "--user", "stop", self.unit], check=False)
        subprocess.run(["systemctl", "--user", "stop", f"{self.prefix}-*"], check=False)
        subprocess.run(["systemctl", "--user", "reset-failed", f"{self.prefix}-*"], check=False,
                       stderr=subprocess.DEVNULL)
        log = (self.base / "state" / "worker.log")
        text = log.read_text(errors="replace") if log.exists() else ""
        return {"listing": state.splitlines()[:6], "worker_log_tail": text.splitlines()[-3:]}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--minutes", type=float, default=10.0)
    parser.add_argument("--warmup", type=float, default=30.0)
    parser.add_argument("--cpus", default="14,15")
    parser.add_argument("--versions", default="0.6.3,0.7.0")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    root = out / "sandboxes"
    boxes = [Sandbox(root, version, scenario, args.cpus)
             for version in args.versions.split(",") for scenario in ("idle", "busy")]
    started = time.time()
    try:
        for box in boxes:
            box.start()
        time.sleep(args.warmup)
        before = {box.unit: box.read() for box in boxes}
        time.sleep(args.minutes * 60)
        after = {box.unit: box.read() for box in boxes}
    finally:
        details = {box.unit: box.stop() for box in boxes}
    results = []
    for box in boxes:
        a, b = before[box.unit], after[box.unit]
        span = b["unix"] - a["unix"]
        cpu = (b["cpu_usec"] - a["cpu_usec"]) / 1e6 if a["cpu_usec"] is not None and b["cpu_usec"] is not None else None
        results.append({
            "version": box.version, "scenario": box.scenario, "unit": box.unit, "window_s": round(span, 1),
            "cpu_s": None if cpu is None else round(cpu, 3),
            "cpu_s_per_hour": None if cpu is None else round(cpu / span * 3600, 3),
            "percent_of_one_core": None if cpu is None else round(100 * cpu / span, 4),
            "rss_bytes": b["rss"], "memory_current": b["memory_current"], "memory_peak": b["memory_peak"],
            **details[box.unit],
        })
    receipt = {"schema": "simqueue-bench-overhead-v1", "started_unix": started, "minutes": args.minutes,
               "warmup_s": args.warmup, "cpus": args.cpus, "results": results}
    (out / "bench-overhead.json").write_text(json.dumps(receipt, indent=2) + "\n")
    for item in results:
        print(f"{item['version']:<6} {item['scenario']:<5} {item['cpu_s_per_hour']} CPU-s/h "
              f"({item['percent_of_one_core']}% of a core), RSS {item['rss_bytes'] / 2**20:.1f} MiB, cgroup memory "
              f"{(item['memory_current'] or 0) / 2**20:.1f} MiB (peak {(item['memory_peak'] or 0) / 2**20:.1f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
