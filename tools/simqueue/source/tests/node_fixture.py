"""Shared fixture for node-scheduler tests that start real systemd user units.

Real units run only on the CPUs in ``SIMQ_TEST_CPUS``, which ``run-tests.py``
sets after taking those cores' lease files.  Without it these tests skip.
Each test uses its own unit-name prefix, state directory and temporary lease
files, and stops its units in tearDown.
"""

from __future__ import annotations

import fcntl
import json
import os
import secrets
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Callable

from simqueue.config import parse_cpu_list


REPOSITORY = Path(__file__).resolve().parents[1]
PYTHON_FROZEN = "/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen"
FAKE_GPU = "GPU-test-0000-0000-0000-000000000001"
PROBE = r'''
import json, os, sys, time
path, seconds, code = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
record = {"job": os.environ.get("SIMQ_JOB_ID"), "cpus": os.environ.get("SIMQ_CPU_LIST"),
          "affinity": sorted(os.sched_getaffinity(0)), "start": time.time()}
time.sleep(seconds)
record["end"] = time.time()
with open(path, "a") as handle:
    handle.write(json.dumps(record) + "\n")
sys.exit(code)
'''


def systemd_user_available() -> bool:
    try:
        result = subprocess.run(
            ["systemctl", "--user", "is-system-running"],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, timeout=10, check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.stdout.strip() in ("running", "degraded")


def test_cpus() -> tuple[int, ...]:
    text = os.environ.get("SIMQ_TEST_CPUS", "")
    if not text:
        return ()
    return tuple(parse_cpu_list(text))


class NodeFixture:
    min_cpus = 2
    needs_units = True

    def setUp(self) -> None:  # noqa: N802
        self.cpus = test_cpus()
        if self.needs_units:
            if len(self.cpus) < self.min_cpus:
                self.skipTest(f"needs SIMQ_TEST_CPUS with {self.min_cpus}+ leased CPUs (use run-tests.py)")
            if not systemd_user_available():
                self.skipTest("systemd user manager unavailable")
        self.temporary = tempfile.TemporaryDirectory(prefix="simq06-")
        self.base = Path(self.temporary.name)
        self.projects = self.base / "projects"
        self.projects.mkdir()
        self.state = self.base / "state"
        self.locks = self.base / "locks"
        self.locks.mkdir()
        self.prefix = "simqt" + secrets.token_hex(4)
        self.probe = self.projects / "probe.py"
        self.probe.write_text(PROBE, encoding="utf-8")
        self.intervals = self.projects / "intervals.jsonl"
        self.config_path = self.base / "node.json"
        cpus = list(self.cpus) or [0, 1]
        half = max(1, len(cpus) // 2)
        self.part_a = cpus[:half] if len(cpus) > 1 else cpus
        self.part_b = cpus[half:] if len(cpus) > 1 else cpus
        self.write_config()
        self.environment = dict(os.environ)
        self.environment.update(
            {
                "SIMQUEUE_STATE_DIR": str(self.state),
                "SIMQUEUE_PROJECT_ROOT": str(self.projects),
                "SIMQUEUE_CONFIG": str(self.config_path),
                "SIMQUEUE_NO_AUTOSTART": "1",
                "SIMQUEUE_CANCEL_GRACE_SECONDS": "1",
            }
        )
        self._saved_environ = dict(os.environ)
        os.environ.update(self.environment)
        self.launcher = self.make_launcher(REPOSITORY)
        self.workers: list[subprocess.Popen] = []
        self.holders: list[subprocess.Popen] = []
        from simqueue.store import QueueStore

        self.store = QueueStore()

    def tearDown(self) -> None:  # noqa: N802
        for proc in self.workers + self.holders:
            if proc.poll() is None:
                proc.kill()
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                pass
        if self.needs_units:
            subprocess.run(["systemctl", "--user", "stop", f"{self.prefix}-*"],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
            subprocess.run(["systemctl", "--user", "reset-failed", f"{self.prefix}-*"],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        os.environ.clear()
        os.environ.update(self._saved_environ)
        self.temporary.cleanup()

    # ----- configuration ---------------------------------------------------

    def config_document(self) -> dict[str, Any]:
        def fmt(values):
            return ",".join(str(v) for v in values)

        return {
            "schema": "simqueue-node-config-v1",
            "node": "testnode",
            "default_partition": "pa",
            "memory_schedulable": "8G",
            "allow_smt_siblings": True,
            "gpus": [{"uuid": FAKE_GPU, "index": 0, "name": "fake", "memory": "12G", "budget": "10G"}],
            "lock_files": {
                "enabled": True,
                "core_path": str(self.locks / "core{core}.lock"),
                "gpu_path": str(self.locks / "gpu-{uuid}.lock"),
                "gpu_extra": {FAKE_GPU: [str(self.locks / "gpu.lock")]},
            },
            "launch": {
                "unit_prefix": self.prefix,
                "slice": "simqtest.slice",
                "wrapper_python": PYTHON_FROZEN,
                "kill_grace_seconds": 2,
                "poll_seconds": 0.2,
            },
            "partitions": [
                {"name": "pa", "cpus": fmt(self.part_a), "default_mem": "256M", "max_mem": "2G",
                 "default_time": "10:00", "max_time": "1:00:00"},
                {"name": "pb", "cpus": fmt(self.part_b), "gpus": [FAKE_GPU], "default_mem": "256M",
                 "max_mem": "2G", "default_time": "10:00", "max_time": "1:00:00", "default_cpus": 1},
            ],
        }

    def write_config(self, document: dict[str, Any] | None = None) -> None:
        self.config_path.write_text(json.dumps(document or self.config_document(), indent=2), encoding="utf-8")

    def make_launcher(self, repository: Path) -> Path:
        launcher = self.base / f"simq-launcher-{secrets.token_hex(3)}"
        launcher.write_text(
            "import sys\n"
            "from pathlib import Path\n"
            f"sys.path.insert(0, {str(repository)!r})\n"
            f"sys.executable = {PYTHON_FROZEN!r}\n"
            "from simqueue import exclusive\n"
            f"exclusive.CONFIG_PATH = Path({str(self.base / 'no-exclusive-config')!r})\n"
            "from simqueue.cli import main\n"
            "raise SystemExit(main())\n",
            encoding="utf-8",
        )
        return launcher

    # ----- helpers -----------------------------------------------------------

    def cli(self, *arguments: str, launcher: Path | None = None, env: dict[str, str] | None = None,
            timeout: float = 60) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [PYTHON_FROZEN, str(launcher or self.launcher), *arguments],
            cwd=self.projects,
            env=env or self.environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )

    def submit(self, *arguments: str, env: dict[str, str] | None = None) -> int:
        result = self.cli("submit", "--json", *arguments, env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        return int(json.loads(result.stdout)["job_id"])

    def probe_command(self, seconds: float, code: int = 0) -> list[str]:
        return ["--", PYTHON_FROZEN, str(self.probe), str(self.intervals), str(seconds), str(code)]

    def start_worker(self, *arguments: str, launcher: Path | None = None) -> subprocess.Popen:
        log = open(self.base / f"worker-{len(self.workers)}.out", "ab")
        proc = subprocess.Popen(
            [PYTHON_FROZEN, str(launcher or self.launcher), "worker", *arguments],
            cwd=self.projects,
            env=self.environment,
            stdin=subprocess.DEVNULL,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        log.close()
        self.workers.append(proc)
        return proc

    def run_drain(self, launcher: Path | None = None, timeout: float = 120) -> subprocess.CompletedProcess[str]:
        result = self.cli("worker", "--drain", launcher=launcher, timeout=timeout)
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        return result

    def wait_until(self, predicate: Callable[[], Any], timeout: float = 60, message: str = "condition") -> Any:
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            value = predicate()
            if value:
                return value
            time.sleep(0.1)
        self.fail(f"timed out waiting for {message}")

    def job(self, job_id: int) -> dict[str, Any]:
        return self.store.require(job_id)

    def state_of(self, job_id: int) -> str:
        return str(self.job(job_id)["state"])

    def unit_active(self, unit: str) -> bool:
        result = subprocess.run(["systemctl", "--user", "is-active", unit], stdout=subprocess.PIPE,
                                stderr=subprocess.DEVNULL, text=True, check=False)
        return result.stdout.strip() in ("active", "activating", "deactivating")

    def log_text(self, job_id: int) -> str:
        return Path(str(self.job(job_id)["log_path"])).read_text(encoding="utf-8", errors="replace")

    def worker_log(self) -> str:
        path = self.state / "worker.log"
        return path.read_text(encoding="utf-8") if path.exists() else ""

    def read_intervals(self) -> dict[int, dict[str, Any]]:
        records: dict[int, dict[str, Any]] = {}
        if self.intervals.exists():
            for line in self.intervals.read_text(encoding="utf-8").splitlines():
                record = json.loads(line)
                records[int(record["job"])] = record
        return records

    def assert_no_cpu_overlap(self) -> None:
        records = list(self.read_intervals().values())
        for i, first in enumerate(records):
            for second in records[i + 1 :]:
                shared = set(first["affinity"]) & set(second["affinity"])
                if not shared:
                    continue
                overlap = min(first["end"], second["end"]) - max(first["start"], second["start"])
                self.assertLessEqual(
                    overlap, 0.0,
                    f"jobs {first['job']} and {second['job']} overlapped on CPUs {sorted(shared)}",
                )

    def hold_lease(self, path: Path, *, shared: bool = False) -> subprocess.Popen:
        """Hold a lease file from a separate (non-queue) process."""
        mode = "LOCK_SH" if shared else "LOCK_EX"
        proc = subprocess.Popen(
            [sys.executable, "-c",
             "import fcntl,sys,time\n"
             f"h=open({str(path)!r},'a'); fcntl.flock(h, fcntl.{mode}); print('held', flush=True)\n"
             "sys.stdin.read()\n"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True,
        )
        self.assertEqual(proc.stdout.readline().strip(), "held")
        self.holders.append(proc)
        return proc

    @staticmethod
    def release_lease(proc: subprocess.Popen) -> None:
        proc.stdin.close()
        proc.wait(timeout=10)

    def lease_holder_pids(self, path: Path) -> list[int]:
        from simqueue.locks import lock_holders

        return [h.pid for h in lock_holders([str(path)]).get(str(path), []) if h.pid]

    def kill9(self, proc: subprocess.Popen) -> None:
        os.kill(proc.pid, signal.SIGKILL)
        proc.wait(timeout=10)

    def copy_runtime(self, version_suffix: str) -> Path:
        target = self.base / f"runtime-{version_suffix}"
        shutil.copytree(REPOSITORY / "simqueue", target / "simqueue",
                        ignore=shutil.ignore_patterns("__pycache__"))
        shutil.copy2(REPOSITORY / "simq", target / "simq")
        init = target / "simqueue" / "__init__.py"
        text = init.read_text(encoding="utf-8").replace('__version__ = "0.6.0"', f'__version__ = "0.6.0+{version_suffix}"')
        init.write_text(text, encoding="utf-8")
        return target


def flock_free(path: Path) -> bool:
    handle = open(path, "a")
    try:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        handle.close()
        return False
    fcntl.flock(handle, fcntl.LOCK_UN)
    handle.close()
    return True
