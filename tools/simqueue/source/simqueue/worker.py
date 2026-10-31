from __future__ import annotations

import fcntl
import ctypes
import json
import os
import shlex
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path
from types import FrameType
from typing import Any

from .config import detect_pcores, format_cpu_list, project_root, state_dir
from .exclusive import current_cgroup, load_config, verify_worker_context
from .store import (
    ENSEMBLE_CONTROL_PROTOCOL,
    ENSEMBLE_DRAIN_EXIT_CODE,
    ENSEMBLE_DRAIN_PROTOCOL_ENV,
    QueueStore,
    utc_now,
)


POLL_SECONDS = 0.25
KILL_SETTLE_SECONDS = 2.0
PR_SET_CHILD_SUBREAPER = 36
EXCLUSIVE_RUNNER_CGROUP = "/simqueue.slice/simqueue-runner.service"
try:
    CANCEL_GRACE_SECONDS = max(
        0.1, float(os.environ.get("SIMQUEUE_CANCEL_GRACE_SECONDS", "10"))
    )
except ValueError:
    CANCEL_GRACE_SECONDS = 10.0


def _boot_id() -> str | None:
    try:
        return Path("/proc/sys/kernel/random/boot_id").read_text(
            encoding="utf-8"
        ).strip()
    except OSError:
        return None


def _process_identity(pid: int) -> tuple[str | None, int, int] | None:
    """Return boot ID, process start ticks, and process group for *pid*."""
    record = _process_record(pid)
    if record is None:
        return None
    _state, process_group, start_ticks = record
    return _boot_id(), start_ticks, process_group


def _process_record(pid: int) -> tuple[str, int, int] | None:
    """Return process state, process group, and start ticks from procfs."""
    try:
        stat = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
        remainder = stat[stat.rfind(")") + 2 :].split()
        state = remainder[0]  # proc(5) field 3
        process_group = int(remainder[2])  # proc(5) field 5
        start_ticks = int(remainder[19])  # proc(5) field 22
        return state, process_group, start_ticks
    except (OSError, ValueError, IndexError):
        return None


def _process_parent_and_group(pid: int) -> tuple[str, int, int] | None:
    """Return process state, parent PID, and process group from procfs."""
    try:
        stat_line = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
        remainder = stat_line[stat_line.rfind(")") + 2 :].split()
        return remainder[0], int(remainder[1]), int(remainder[2])
    except (OSError, ValueError, IndexError):
        return None


def _descendant_process_groups(root_pid: int) -> set[int]:
    """Find live process groups below *root_pid*, independent of environment."""
    children: dict[int, list[int]] = {}
    records: dict[int, tuple[str, int]] = {}
    try:
        entries = Path("/proc").iterdir()
    except OSError:
        return set()
    for entry in entries:
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        record = _process_parent_and_group(pid)
        if record is None:
            continue
        state, parent, process_group = record
        children.setdefault(parent, []).append(pid)
        records[pid] = (state, process_group)

    descendants: set[int] = set()
    pending = list(children.get(root_pid, ()))
    while pending:
        pid = pending.pop()
        if pid in descendants:
            continue
        descendants.add(pid)
        pending.extend(children.get(pid, ()))
    return {
        records[pid][1]
        for pid in descendants
        if pid in records and records[pid][0] != "Z"
    }


def _enable_child_subreaper() -> bool:
    """Make orphaned daemon descendants reparent to this worker on Linux."""
    try:
        libc = ctypes.CDLL(None, use_errno=True)
        return libc.prctl(PR_SET_CHILD_SUBREAPER, 1, 0, 0, 0) == 0
    except (AttributeError, OSError):
        return False


def _reap_exited_children() -> None:
    """Reap exited orphan descendants adopted by the worker."""
    while True:
        try:
            pid, _status = os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            return
        except InterruptedError:
            continue
        if pid == 0:
            return


def _live_process_group_members(
    process_group: int, *, minimum_start_ticks: int | None = None
) -> set[int]:
    """Find live (non-zombie) processes still belonging to a job's group."""
    members: set[int] = set()
    try:
        entries = Path("/proc").iterdir()
    except OSError:
        return members
    for entry in entries:
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        record = _process_record(pid)
        if record is None:
            continue
        state, group, start_ticks = record
        if (
            group == process_group
            and state != "Z"
            and (minimum_start_ticks is None or start_ticks >= minimum_start_ticks)
        ):
            members.add(pid)
    return members


def _process_has_launch_token(pid: int, token: str) -> bool:
    expected = f"SIMQ_LAUNCH_TOKEN={token}".encode("utf-8")
    try:
        environment = Path(f"/proc/{pid}/environ").read_bytes().split(b"\0")
    except OSError:
        return False
    return expected in environment


def _verified_job_process_groups(job: dict[str, Any]) -> set[int]:
    """Resolve a job's live groups without risking a recycled PID/PGID.

    The random launch token covers the small Popen-to-database window in
    which a hard worker death can leave no recorded PID yet. Exact procfs
    identity remains a compatibility fallback for older queue records.
    """
    token = job.get("launch_token")
    groups: set[int] = set()
    if token:
        try:
            entries = Path("/proc").iterdir()
        except OSError:
            entries = iter(())
        for entry in entries:
            if not entry.name.isdigit():
                continue
            pid = int(entry.name)
            record = _process_record(pid)
            if (
                record is not None
                and record[0] != "Z"
                and _process_has_launch_token(pid, str(token))
            ):
                groups.add(record[1])
    if groups:
        groups.discard(os.getpgrp())
        return groups
    if _matches_recorded_process(job):
        return {int(job["pgid"])}
    return set()


def _matches_recorded_process(job: dict[str, Any]) -> bool:
    pid = job.get("pid")
    recorded_ticks = job.get("process_start_ticks")
    recorded_group = job.get("pgid")
    recorded_boot = job.get("boot_id")
    if not all(value is not None for value in (pid, recorded_ticks, recorded_group, recorded_boot)):
        return False
    record = _process_record(int(pid))
    if record is None or record[0] == "Z":
        return False
    _state, process_group, ticks = record
    return (
        _boot_id() == recorded_boot
        and ticks == int(recorded_ticks)
        and process_group == int(recorded_group)
        and process_group == int(pid)
    )


class Worker:
    def __init__(self, *, drain: bool = False, wait_for_lock: bool = False) -> None:
        self.store = QueueStore()
        self.drain = drain
        self.wait_for_lock = wait_for_lock
        self.stop_requested = False
        self.stop_signal: int | None = None
        self.current_process: subprocess.Popen[bytes] | None = None
        self.current_job: dict[str, Any] | None = None
        self.exclusive_mode = False
        self.subreaper_enabled = False
        self._lock_file: Any = None

    def _worker_log(self, message: str) -> None:
        path = state_dir() / "worker.log"
        line = f"[{utc_now()}] {message}\n"
        try:
            with path.open("a", encoding="utf-8") as handle:
                handle.write(line)
            path.chmod(0o600)
        except OSError:
            print(line, file=sys.stderr, end="")

    def _acquire_lock(self) -> bool:
        path = state_dir() / "worker.lock"
        handle = path.open("a+", encoding="utf-8")
        path.chmod(0o600)
        try:
            operation = fcntl.LOCK_EX
            if not self.wait_for_lock:
                operation |= fcntl.LOCK_NB
            fcntl.flock(handle.fileno(), operation)
        except BlockingIOError:
            handle.close()
            return False
        handle.seek(0)
        handle.truncate()
        handle.write(f"{os.getpid()}\n")
        handle.flush()
        self._lock_file = handle
        return True

    def _release_lock(self) -> None:
        if self._lock_file is not None:
            self._lock_file.close()
            self._lock_file = None

    def _handle_stop(self, signum: int, _frame: FrameType | None) -> None:
        self.stop_requested = True
        self.stop_signal = signum

    def _install_signal_handlers(self) -> None:
        signal.signal(signal.SIGTERM, self._handle_stop)
        signal.signal(signal.SIGINT, self._handle_stop)

    @staticmethod
    def _signal_group(process: subprocess.Popen[bytes], sig: signal.Signals) -> None:
        try:
            os.killpg(process.pid, sig)
        except ProcessLookupError:
            pass

    @staticmethod
    def _signal_groups(process_groups: set[int], sig: signal.Signals) -> None:
        for process_group in sorted(process_groups, reverse=True):
            try:
                os.killpg(process_group, sig)
            except ProcessLookupError:
                pass

    def _exclusive_process_groups(self) -> set[int]:
        if current_cgroup() != EXCLUSIVE_RUNNER_CGROUP:
            raise RuntimeError("exclusive worker moved outside its protected cgroup")
        root = Path("/sys/fs/cgroup") / EXCLUSIVE_RUNNER_CGROUP.lstrip("/")
        try:
            directories = [root]
            directories.extend(
                path for path in root.rglob("*") if path.is_dir()
            )
        except OSError as error:
            raise RuntimeError(f"cannot enumerate exclusive runner cgroup: {error}") from error
        groups: set[int] = set()
        for directory in directories:
            try:
                values = (directory / "cgroup.procs").read_text(
                    encoding="ascii"
                ).split()
            except FileNotFoundError:
                continue
            except OSError as error:
                raise RuntimeError(
                    f"cannot inspect exclusive runner processes: {error}"
                ) from error
            for value in values:
                try:
                    pid = int(value)
                except ValueError as error:
                    raise RuntimeError(
                        f"invalid PID in exclusive runner cgroup: {value!r}"
                    ) from error
                if pid == os.getpid():
                    continue
                record = _process_record(pid)
                if record is not None and record[0] != "Z":
                    groups.add(record[1])
        groups.discard(os.getpgrp())
        return groups

    def _job_process_groups(self, job: dict[str, Any]) -> set[int]:
        groups = _verified_job_process_groups(job)
        if (
            self.subreaper_enabled
            and self.current_job is not None
            and self.current_job.get("launch_token") == job.get("launch_token")
        ):
            groups.update(_descendant_process_groups(os.getpid()))
        if self.exclusive_mode:
            groups.update(self._exclusive_process_groups())
        groups.discard(os.getpgrp())
        return groups

    def _terminate_verified_job_processes(
        self,
        job: dict[str, Any],
        *,
        log_handle: Any | None = None,
        term_message: str | None = None,
    ) -> bool:
        """Terminate every process group carrying *job*'s launch token.

        Nested launchers are allowed to create new sessions, so a job cannot
        be represented by only its original process group.  Rescanning during
        both the TERM grace period and KILL settle period also catches a child
        which creates another session while cancellation is in progress.
        """
        saw_process = False
        term_signaled_groups: set[int] = set()
        quiet_scans = 0
        term_deadline = time.monotonic() + CANCEL_GRACE_SECONDS
        while time.monotonic() < term_deadline:
            process_groups = self._job_process_groups(job)
            if process_groups:
                saw_process = True
                quiet_scans = 0
                new_groups = process_groups - term_signaled_groups
                if new_groups:
                    self._signal_groups(new_groups, signal.SIGTERM)
                    term_signaled_groups.update(new_groups)
                    if log_handle is not None and term_message is not None:
                        self._write_log(log_handle, term_message)
                        term_message = None
            else:
                quiet_scans += 1
                # A second scan closes the fork/setsid race at the edge of a
                # procfs walk without imposing the full grace period on a
                # cooperative cancellation.
                if quiet_scans >= 2:
                    return saw_process
            time.sleep(POLL_SECONDS)

        quiet_scans = 0
        kill_logged = False
        kill_deadline = time.monotonic() + KILL_SETTLE_SECONDS
        while time.monotonic() < kill_deadline:
            process_groups = self._job_process_groups(job)
            if process_groups:
                saw_process = True
                quiet_scans = 0
                self._signal_groups(process_groups, signal.SIGKILL)
                if log_handle is not None and not kill_logged:
                    self._write_log(
                        log_handle, "simq: grace period expired; sent SIGKILL\n"
                    )
                    kill_logged = True
            else:
                quiet_scans += 1
                if quiet_scans >= 2:
                    return saw_process
            time.sleep(POLL_SECONDS)

        survivors = self._job_process_groups(job)
        if survivors:
            rendered = ",".join(str(group) for group in sorted(survivors))
            raise RuntimeError(
                f"job process group(s) {rendered} survived SIGKILL; "
                "queue cannot continue safely"
            )
        return saw_process

    def _clean_up_current_process_after_error(self) -> None:
        process = self.current_process
        if process is None:
            return
        if self.current_job is not None:
            self._terminate_verified_job_processes(self.current_job)
        elif process.poll() is None or _live_process_group_members(process.pid):
            # The launch failed before an identity-bearing current_job could
            # be recorded.  At that point only the original session exists.
            self._signal_group(process, signal.SIGKILL)
        try:
            process.wait(timeout=KILL_SETTLE_SECONDS)
        except subprocess.TimeoutExpired as error:
            raise RuntimeError(
                f"job leader {process.pid} did not exit; queue cannot continue safely"
            ) from error
        self.current_process = None
        self.current_job = None

    def _terminate_verified_stale_process(self, job: dict[str, Any]) -> bool:
        return self._terminate_verified_job_processes(job)

    def _recover_stale_jobs(self) -> None:
        for job in self.store.active_jobs():
            terminated = self._terminate_verified_stale_process(job)
            detail = "verified stale process group terminated" if terminated else "no matching live process"
            if self.store.requeue_drained_ensemble(int(job["id"])):
                self._worker_log(
                    f"recovered drained ensemble {job['id']} as the same queued job "
                    f"({detail})"
                )
                continue
            self.store.mark_interrupted(
                int(job["id"]),
                f"worker restarted; {detail}; retry explicitly if appropriate",
            )
            self._worker_log(f"marked stale job {job['id']} interrupted ({detail})")

    @staticmethod
    def _write_log(handle: Any, text: str) -> None:
        handle.write(text.encode("utf-8", errors="replace"))
        handle.flush()

    def _finish_launch_error(
        self, job: dict[str, Any], log_handle: Any, message: str
    ) -> None:
        self._write_log(log_handle, f"simq: launch failed: {message}\n")
        self._write_log(log_handle, f"===== simq job finished {utc_now()} =====\n")
        self.store.finish(
            int(job["id"]), state="failed", exit_code=None, message=message
        )

    def _run_job(self, job: dict[str, Any], pcores: tuple[int, ...]) -> None:
        job_id = int(job["id"])
        requested = int(job["requested_cores"])
        cpus = pcores[:requested]
        cpu_list = format_cpu_list(cpus)
        log_path = Path(str(job["log_path"]))
        log_path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
        with log_path.open("ab", buffering=0) as log_handle:
            log_path.chmod(0o600)
            try:
                argv = json.loads(str(job["argv_json"]))
            except json.JSONDecodeError as error:
                self._finish_launch_error(job, log_handle, f"invalid stored argv: {error}")
                return
            if not isinstance(argv, list) or not argv or not all(
                isinstance(argument, str) for argument in argv
            ):
                self._finish_launch_error(job, log_handle, "stored argv is not a nonempty string list")
                return
            if requested > len(pcores):
                self._finish_launch_error(
                    job,
                    log_handle,
                    f"job requests {requested} cores but only {len(pcores)} P-cores are available",
                )
                return
            submitted_cwd = Path(str(job["cwd"]))
            try:
                cwd = submitted_cwd.resolve(strict=True)
            except OSError:
                cwd = submitted_cwd
            if not cwd.is_dir():
                self._finish_launch_error(job, log_handle, f"working directory is unavailable: {submitted_cwd}")
                return
            if not cwd.is_relative_to(project_root()):
                self._finish_launch_error(
                    job,
                    log_handle,
                    f"working directory now resolves outside {project_root()}: {cwd}",
                )
                return
            if self.store.cancellation_requested(job_id):
                self._write_log(log_handle, "simq: canceled before process launch\n")
                self.store.finish(
                    job_id,
                    state="canceled",
                    exit_code=None,
                    message="canceled before process launch",
                )
                return
            taskset = shutil.which("taskset")
            if taskset is None:
                self._finish_launch_error(job, log_handle, "taskset is not installed")
                return
            command = [taskset, "--cpu-list", cpu_list, *argv]
            try:
                environment = dict(json.loads(str(job["env_json"])))
            except (json.JSONDecodeError, TypeError, ValueError) as error:
                self._finish_launch_error(job, log_handle, f"invalid stored environment: {error}")
                return
            environment.update(
                {
                    "PWD": str(cwd),
                    "SIMQ_JOB_ID": str(job_id),
                    "SIMQ_CPU_LIST": cpu_list,
                    "SIMQ_CPU_COUNT": str(len(cpus)),
                    "SIMQ_LOG": str(log_path),
                    "SIMQ_LAUNCH_TOKEN": str(job["launch_token"]),
                    # A new ensemble supervisor may return the reserved drain
                    # exit only when its outer worker advertises that it knows
                    # how to suspend/requeue the parent instead of failing it.
                    ENSEMBLE_DRAIN_PROTOCOL_ENV: (
                        f"{ENSEMBLE_CONTROL_PROTOCOL}:{job['launch_token']}"
                    ),
                }
            )
            self._write_log(log_handle, f"===== simq job {job_id}: {job['name']} =====\n")
            self._write_log(log_handle, f"submitted: {job['created_at']}\n")
            self._write_log(log_handle, f"started:   {utc_now()}\n")
            self._write_log(log_handle, f"cwd:       {cwd}\n")
            self._write_log(log_handle, f"P-cores:   {cpu_list}\n")
            self._write_log(log_handle, f"command:   {shlex.join(argv)}\n")
            self._write_log(log_handle, "========================================\n")
            try:
                process = subprocess.Popen(
                    command,
                    cwd=cwd,
                    env=environment,
                    stdin=subprocess.DEVNULL,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )
            except (OSError, ValueError) as error:
                self._finish_launch_error(job, log_handle, str(error))
                return

            self.current_process = process
            identity = _process_identity(process.pid)
            boot = _boot_id()
            ticks = None
            process_group = process.pid
            if identity is not None:
                boot, ticks, process_group = identity
            running_job = dict(job)
            running_job.update(
                {
                    "pid": process.pid,
                    "pgid": process_group,
                    "boot_id": boot,
                    "process_start_ticks": ticks,
                }
            )
            self.current_job = running_job
            self.store.mark_running(
                job_id,
                pid=process.pid,
                pgid=process_group,
                boot_id=boot,
                process_start_ticks=ticks,
            )
            reason: str | None = None
            completion_quiet_scans = 0
            while True:
                if reason is None and self.stop_requested:
                    reason = "worker stopped"
                elif reason is None and self.store.cancellation_requested(job_id):
                    reason = "canceled by user"

                if reason is not None:
                    action = (
                        "simq: cancellation requested; sent SIGTERM\n"
                        if reason == "canceled by user"
                        else "simq: worker stopping; sent SIGTERM\n"
                    )
                    self._terminate_verified_job_processes(
                        running_job,
                        log_handle=log_handle,
                        term_message=action,
                    )
                    break

                leader_exit_code = process.poll()
                if leader_exit_code is not None:
                    if self._job_process_groups(running_job):
                        completion_quiet_scans = 0
                    else:
                        completion_quiet_scans += 1
                        # Require two empty procfs/cgroup scans so a final
                        # fork+setsid at leader exit cannot slip into the next
                        # serial job between enumeration and reparenting.
                        if completion_quiet_scans >= 2:
                            break
                time.sleep(POLL_SECONDS)

            exit_code = process.wait()
            _reap_exited_children()
            self.current_process = None
            self.current_job = None
            if (
                reason is None
                and exit_code == ENSEMBLE_DRAIN_EXIT_CODE
                and self.store.requeue_drained_ensemble(job_id)
            ):
                control = self.store.ensemble_control_snapshot(job_id)
                suspended_state = (
                    "paused" if control["launch_paused"] else "resuming"
                )
                self._write_log(log_handle, "========================================\n")
                self._write_log(log_handle, f"suspended: {utc_now()}\n")
                self._write_log(log_handle, f"state:     {suspended_state}\n")
                self._write_log(
                    log_handle,
                    f"remaining: {control['queued']} queued member(s)\n",
                )
                self._write_log(
                    log_handle,
                    f"===== simq ensemble {job_id} allocation released =====\n",
                )
                self._worker_log(
                    f"suspended ensemble {job_id}; same job ID retained with "
                    f"{control['queued']} queued member(s)"
                )
                return
            if reason == "canceled by user":
                final_state = "canceled"
            elif reason == "worker stopped":
                final_state = "interrupted"
            elif exit_code == 0:
                final_state = "succeeded"
            else:
                final_state = "failed"
            message = reason
            if final_state == "failed" and message is None:
                message = f"command exited with status {exit_code}"
            self._write_log(log_handle, "========================================\n")
            self._write_log(log_handle, f"finished:  {utc_now()}\n")
            self._write_log(log_handle, f"state:     {final_state}\n")
            self._write_log(log_handle, f"exit code: {exit_code}\n")
            self._write_log(log_handle, f"===== simq job {job_id} finished =====\n")
            self.store.finish(
                job_id,
                state=final_state,
                exit_code=exit_code,
                message=message,
            )

    def run(self) -> int:
        if not self._acquire_lock():
            self._worker_log("another worker already owns the queue lock")
            return 0
        # This check deliberately happens after the lock is acquired.  A
        # legacy --wait-for-lock worker can begin before exclusive migration,
        # then wake after the root-owned marker is published.  Revalidating
        # here prevents that process from ever dispatching outside the shield.
        try:
            exclusive = load_config()
            if exclusive is not None:
                verify_worker_context(exclusive)
        except RuntimeError as error:
            self._worker_log(f"cannot start protected worker: {error}")
            self._release_lock()
            return 2
        self.exclusive_mode = exclusive is not None
        self.subreaper_enabled = _enable_child_subreaper()
        if not self.exclusive_mode and not self.subreaper_enabled:
            self._worker_log(
                "cannot start: child-subreaper support is required to preserve "
                "strict job serialization"
            )
            self._release_lock()
            return 2
        self._install_signal_handlers()
        try:
            pcores = detect_pcores()
        except (RuntimeError, ValueError) as error:
            self._worker_log(f"cannot start: {error}")
            self._release_lock()
            return 2
        self._worker_log(
            f"worker started; P-cores={format_cpu_list(pcores)}; drain={self.drain}"
        )
        self._recover_stale_jobs()
        worker_exit_code = 0
        while not self.stop_requested:
            if exclusive is not None:
                try:
                    # Recheck between jobs so a privileged cgroup change or
                    # CPU-topology event cannot silently expose a later job.
                    verify_worker_context(exclusive)
                except RuntimeError as error:
                    self._worker_log(
                        f"protected worker lost its P-core shield: {error}"
                    )
                    worker_exit_code = 2
                    break
            job = self.store.claim_next(pcores)
            if job is None:
                if self.drain:
                    break
                time.sleep(1.0)
                continue
            self._worker_log(f"starting job {job['id']}: {job['name']}")
            try:
                self._run_job(job, pcores)
            except Exception as error:  # preserve the queue after an unexpected job error
                self._worker_log(f"unexpected error in job {job['id']}: {error!r}")
                self._clean_up_current_process_after_error()
                self.store.finish(
                    int(job["id"]),
                    state="failed",
                    exit_code=None,
                    message=f"worker error: {error!r}",
                )
        suffix = f" after signal {self.stop_signal}" if self.stop_signal else ""
        self._worker_log(f"worker stopped{suffix}")
        self._release_lock()
        return worker_exit_code


def run_worker(*, drain: bool = False, wait_for_lock: bool = False, once: bool = False) -> int:
    """Run the node scheduler when a node config exists, else the 0.5 worker."""
    from .nodeconfig import ConfigError, load_node_config

    try:
        node_config = load_node_config()
    except ConfigError as error:
        print(f"simq: {error}", file=sys.stderr)
        return 2
    if node_config is not None:
        from .nodeworker import NodeWorker

        return NodeWorker(
            node_config, drain=drain, wait_for_lock=wait_for_lock, once=once
        ).run()
    return Worker(drain=drain, wait_for_lock=wait_for_lock).run()
