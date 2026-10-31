"""SimQueue 0.6 job-unit wrapper.

This file is executed by systemd as the main process of each job's transient
user service (``simq-job-<id>-<token>.service``).  It deliberately imports
nothing from the ``simqueue`` package so that a running job never depends on
the scheduler's code, process or version.

``jobwrap.py run SPEC``
    Verify the unit's CPU affinity and working directory, take the job's
    lease files (inherited by the job's processes, so they stay held while
    any of them lives), write the status file, run the command with the
    recorded environment, wait until every other process in the unit's
    cgroup has exited, record the exit status and exit with it.

A borrow job (spec ``"borrow": true``, SimQueue 0.6.1) holds the shared
borrow leases of its cores instead of the exclusive core leases, and runs
the command under ``SCHED_IDLE`` (set by the unit; checked and, if needed,
set here before the command starts, so every process of the job inherits it).

``jobwrap.py stop-post SPEC``
    ExecStopPost hook.  Records systemd's $SERVICE_RESULT, $EXIT_CODE and
    $EXIT_STATUS so that the final result is recoverable even if the wrapper
    itself was killed (out-of-memory kill, time limit, SIGKILL).
"""

from __future__ import annotations

import ctypes
import fcntl
import json
import os
import shlex
import signal
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path


SPEC_SCHEMA = "simqueue-job-spec-v1"
STATUS_SCHEMA = "simqueue-job-status-v1"
STOP_POST_SCHEMA = "simqueue-job-stop-post-v1"
EXIT_REFUSED = 254
PR_SET_CHILD_SUBREAPER = 36
POLL_SECONDS = 0.25


def _now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z")


def _boot_id() -> str | None:
    try:
        return Path("/proc/sys/kernel/random/boot_id").read_text(encoding="utf-8").strip()
    except OSError:
        return None


def _atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(path.name + ".tmp")
    data = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")
    descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    try:
        os.write(descriptor, data)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.replace(temporary, path)


def _log(text: str) -> None:
    # stdout/stderr of the unit are appended to the job log by systemd.
    sys.stdout.write(text)
    sys.stdout.flush()


def _own_cgroup() -> str | None:
    try:
        for line in Path("/proc/self/cgroup").read_text(encoding="utf-8").splitlines():
            if line.startswith("0::"):
                return line[3:]
    except OSError:
        return None
    return None


def _live_cgroup_pids(cgroup: str) -> set[int]:
    pids: set[int] = set()
    try:
        text = (Path("/sys/fs/cgroup") / cgroup.lstrip("/") / "cgroup.procs").read_text(encoding="ascii")
    except OSError:
        return pids
    for value in text.split():
        try:
            pid = int(value)
        except ValueError:
            continue
        try:
            stat = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
            state = stat[stat.rfind(")") + 2 :].split()[0]
        except (OSError, IndexError):
            continue
        if state != "Z":
            pids.add(pid)
    return pids


def _reap() -> None:
    while True:
        try:
            pid, _status = os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            return
        except InterruptedError:
            continue
        if pid == 0:
            return


def _cpu_set(text: str) -> set[int]:
    cpus: set[int] = set()
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        low, _, high = part.partition("-")
        cpus.update(range(int(low), int(high or low) + 1))
    return cpus


def _ensure_sched_idle() -> str:
    """Make this process SCHED_IDLE (inherited by the command); return the
    resulting policy name."""
    names = {getattr(os, "SCHED_OTHER", 0): "other", getattr(os, "SCHED_BATCH", 3): "batch",
             getattr(os, "SCHED_IDLE", 5): "idle", getattr(os, "SCHED_FIFO", 1): "fifo",
             getattr(os, "SCHED_RR", 2): "rr"}
    try:
        if os.sched_getscheduler(0) != os.SCHED_IDLE:
            os.sched_setscheduler(0, os.SCHED_IDLE, os.sched_param(0))
        return names.get(os.sched_getscheduler(0), "unknown")
    except (AttributeError, OSError) as error:
        return f"error: {error}"


class Wrapper:
    def __init__(self, spec_path: Path) -> None:
        self.spec_path = spec_path
        self.spec = json.loads(spec_path.read_text(encoding="utf-8"))
        if self.spec.get("schema") != SPEC_SCHEMA:
            raise RuntimeError(f"unsupported job spec schema {self.spec.get('schema')!r}")
        self.status_path = Path(self.spec["status_path"])
        self.status: dict = {
            "schema": STATUS_SCHEMA,
            "job_id": self.spec["job_id"],
            "launch_token": self.spec["launch_token"],
            "unit": self.spec["unit"],
            "phase": "starting",
            "wrapper_pid": os.getpid(),
            "boot_id": _boot_id(),
            "invocation_id": os.environ.get("INVOCATION_ID"),
            "wrapper_started_at": _now(),
            "wrapper_started_unix": time.time(),
        }
        self.stop_signal: int | None = None

    def save(self) -> None:
        _atomic_write(self.status_path, self.status)

    def refuse(self, reason: str, *, kind: str) -> int:
        self.status.update(
            {
                "phase": "refused",
                "refusal": kind,
                "reason": reason,
                "command_started": False,
                "finished_at": _now(),
                "finished_unix": time.time(),
            }
        )
        self.save()
        _log(f"simq: launch refused ({kind}): {reason}\n")
        return EXIT_REFUSED

    def _handle_signal(self, signum: int, _frame: object) -> None:
        # systemd signals every process of the unit's cgroup itself; the
        # wrapper keeps waiting so that it can record the command's status.
        if self.stop_signal is None:
            self.stop_signal = signum

    def run(self) -> int:
        spec = self.spec
        signal.signal(signal.SIGTERM, self._handle_signal)
        signal.signal(signal.SIGINT, self._handle_signal)
        signal.signal(signal.SIGHUP, self._handle_signal)
        self.save()
        argv = spec["argv"]
        _log(f"===== simq job {spec['job_id']}: {spec['name']} =====\n")
        _log(f"submitted: {spec.get('submitted_at', '-')}\n")
        _log(f"started:   {_now()}\n")
        _log(f"cwd:       {spec['cwd']}\n")
        _log(f"partition: {spec.get('partition', '-')}\n")
        _log(f"P-cores:   {spec['cpu_list']}\n")
        _log(f"GPUs:      {spec.get('gpu_text', 'none')}\n")
        _log(f"memory:    {spec.get('mem_text', '-')}\n")
        _log(f"time:      {spec.get('time_text', '-')}\n")
        if spec.get("borrow"):
            _log(f"borrowed:  {spec.get('borrow_text', '-')}\n")
        _log(f"unit:      {spec['unit']}\n")
        _log(f"command:   {shlex.join(argv)}\n")
        _log("========================================\n")

        cgroup = _own_cgroup()
        in_unit = cgroup is not None and cgroup.rstrip("/").endswith("/" + spec["unit"])
        if spec.get("require_unit", True) and not in_unit:
            return self.refuse(f"not running inside unit {spec['unit']} (cgroup {cgroup})", kind="placement")
        expected = _cpu_set(spec["cpu_list"])
        observed = set(os.sched_getaffinity(0))
        if observed != expected:
            return self.refuse(
                f"CPU affinity {sorted(observed)} does not match the allocation {sorted(expected)}",
                kind="placement",
            )
        cwd = Path(spec["cwd"])
        try:
            resolved = cwd.resolve(strict=True)
        except OSError:
            return self.refuse(f"working directory is unavailable: {cwd}", kind="cwd")
        root = Path(spec["project_root"])
        if not resolved.is_dir():
            return self.refuse(f"working directory is unavailable: {cwd}", kind="cwd")
        if not resolved.is_relative_to(root):
            return self.refuse(f"working directory now resolves outside {root}: {resolved}", kind="cwd")

        if spec.get("borrow"):
            policy = _ensure_sched_idle()
            self.status["sched_policy"] = policy
            if policy != "idle":
                return self.refuse(f"cannot run a borrow job under SCHED_IDLE (policy {policy})", kind="placement")

        lock_fds: list[int] = []
        held: list[dict] = []
        for lease in spec.get("locks", []):
            path = lease["path"]
            try:
                descriptor = os.open(path, os.O_RDONLY | os.O_CREAT | os.O_APPEND, 0o644)
            except OSError as error:
                return self.refuse(f"cannot open lease file {path}: {error}", kind="lock-error")
            operation = fcntl.LOCK_SH if lease.get("mode") == "shared" else fcntl.LOCK_EX
            try:
                fcntl.flock(descriptor, operation | fcntl.LOCK_NB)
            except BlockingIOError:
                os.close(descriptor)
                return self.refuse(f"lease file {path} is held by another process", kind="lock-busy")
            os.set_inheritable(descriptor, True)
            lock_fds.append(descriptor)
            held.append({"path": path, "mode": lease.get("mode", "exclusive"), "fd": descriptor})
        if held:
            _log("leases:    " + ", ".join(f"{item['path']} ({item['mode']})" for item in held) + "\n")
        self.status["leases"] = held

        try:
            ctypes.CDLL(None, use_errno=True).prctl(PR_SET_CHILD_SUBREAPER, 1, 0, 0, 0)
        except (AttributeError, OSError):
            pass
        environment = {str(k): str(v) for k, v in spec["env"].items()}
        environment["PWD"] = str(resolved)
        try:
            child = subprocess.Popen(
                argv,
                cwd=resolved,
                env=environment,
                stdin=subprocess.DEVNULL,
                stdout=sys.stdout.fileno(),
                stderr=subprocess.STDOUT,
                pass_fds=tuple(lock_fds),
                start_new_session=True,
            )
        except (OSError, ValueError) as error:
            self.status.update(
                {
                    "phase": "finished",
                    "command_started": False,
                    "launch_error": str(error),
                    "returncode": None,
                    "exit_status": EXIT_REFUSED,
                    "finished_at": _now(),
                    "finished_unix": time.time(),
                }
            )
            self.save()
            _log(f"simq: launch failed: {error}\n")
            return EXIT_REFUSED
        self.status.update(
            {
                "phase": "running",
                "command_started": True,
                "child_pid": child.pid,
                "started_at": _now(),
                "started_unix": time.time(),
            }
        )
        self.save()

        returncode: int | None = None
        while returncode is None:
            try:
                pid, wait_status = os.waitpid(child.pid, 0)
            except InterruptedError:
                continue
            except ChildProcessError:
                returncode = child.poll()
                break
            if pid == child.pid:
                returncode = os.waitstatus_to_exitcode(wait_status)
        child.returncode = returncode
        self.status.update({"phase": "draining", "returncode": returncode, "leader_exit_at": _now()})
        self.save()

        # A job is finished when every process it started has exited, as in
        # SimQueue 0.5.  Two consecutive empty scans close fork races.
        quiet = 0
        while quiet < 2:
            _reap()
            others = (_live_cgroup_pids(cgroup) - {os.getpid()}) if in_unit and cgroup else set()
            quiet = 0 if others else quiet + 1
            time.sleep(POLL_SECONDS if quiet < 2 else 0)
        if returncode is None:
            exit_status = EXIT_REFUSED
        elif returncode < 0:
            exit_status = 128 - returncode
        else:
            exit_status = returncode & 0xFF
        self.status.update(
            {
                "phase": "finished",
                "returncode": returncode,
                "exit_status": exit_status,
                "stop_signal": self.stop_signal,
                "finished_at": _now(),
                "finished_unix": time.time(),
            }
        )
        self.save()
        _log("========================================\n")
        _log(f"finished:  {_now()}\n")
        _log(f"exit code: {returncode}\n")
        return exit_status


def stop_post(spec_path: Path) -> int:
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    record = {
        "schema": STOP_POST_SCHEMA,
        "job_id": spec.get("job_id"),
        "launch_token": spec.get("launch_token"),
        "unit": spec.get("unit"),
        "service_result": os.environ.get("SERVICE_RESULT"),
        "exit_code": os.environ.get("EXIT_CODE"),
        "exit_status": os.environ.get("EXIT_STATUS"),
        "invocation_id": os.environ.get("INVOCATION_ID"),
        "boot_id": _boot_id(),
        "recorded_at": _now(),
        "recorded_unix": time.time(),
    }
    _atomic_write(Path(spec["stop_post_path"]), record)
    status_path = Path(spec["status_path"])
    try:
        phase = json.loads(status_path.read_text(encoding="utf-8")).get("phase")
    except (OSError, ValueError):
        phase = None
    if phase not in ("finished", "refused"):
        _log(
            f"simq: unit ended before the wrapper recorded a result: "
            f"result={record['service_result']} code={record['exit_code']} "
            f"status={record['exit_status']}\n"
        )
    return 0


def main(argv: list[str]) -> int:
    if len(argv) != 3 or argv[1] not in ("run", "stop-post"):
        print("usage: jobwrap.py run|stop-post SPEC.json", file=sys.stderr)
        return 2
    os.umask(0o077)
    spec_path = Path(argv[2])
    if argv[1] == "stop-post":
        return stop_post(spec_path)
    try:
        wrapper = Wrapper(spec_path)
    except Exception as error:  # the spec is unreadable: nothing was started
        print(f"simq: cannot read job spec {spec_path}: {error}", file=sys.stderr)
        return EXIT_REFUSED
    return wrapper.run()


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
