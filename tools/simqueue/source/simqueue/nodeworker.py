"""SimQueue 0.6 node scheduler.

One worker per node.  It decides which queued jobs start (priority/FIFO with
EASY backfill over partitions, GPUs and memory), starts each job as its own
transient systemd user service, and reconciles the database with those units.
It never owns job processes: stopping, crashing or upgrading the worker does
not affect running jobs, and a restarted worker (of any 0.6 build) re-adopts
them by unit name and launch token.

0.6.1 adds borrow jobs (``simq submit --borrow``): they start on lendable
cores of any partition after every owner job has been considered, run at idle
CPU priority, and are frozen, requeued or canceled as soon as an owner needs
their cores (see ``borrow.py``).  With no borrow job in the queue the worker
behaves exactly as 0.6.0.
"""

from __future__ import annotations

import fcntl
import json
import os
import signal
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from types import FrameType
from typing import Any, Callable

from . import __version__
from .config import format_cpu_list, parse_cpu_list, project_root, state_dir
from .borrow import Borrower, Preemption, plan_node
from .locks import borrow_lock_plan, busy_resources, lock_plan_for, process_description
from .nodeconfig import NodeConfig, format_size, format_time_limit
from .scheduler import INFINITY, Allocation, Request, Running, build_pool, plan
from .telemetry import NodeTelemetry
from .store import (
    ACTIVE_STATES,
    ENSEMBLE_CONTROL_PROTOCOL,
    ENSEMBLE_DRAIN_EXIT_CODE,
    ENSEMBLE_DRAIN_PROTOCOL_ENV,
    QueueStore,
    utc_now,
)
from .units import (
    JOBWRAP,
    UnitError,
    UnitState,
    build_run_command,
    ensure_borrow_slice,
    freeze_unit,
    read_json,
    show_units,
    start_unit,
    stop_unit,
    thaw_unit,
)


SPEC_SCHEMA = "simqueue-job-spec-v1"
LOCK_BUSY_REFUSAL = "lock-busy"
STOP_RESEND_SECONDS = 30.0


def _boot_id() -> str | None:
    try:
        return Path("/proc/sys/kernel/random/boot_id").read_text(encoding="utf-8").strip()
    except OSError:
        return None


def _parse_time(value: str | None) -> float | None:
    if not value:
        return None
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00")).timestamp()
    except ValueError:
        return None


def _iso(unix: float | None) -> str | None:
    if unix is None:
        return None
    return datetime.fromtimestamp(unix, UTC).isoformat(timespec="seconds").replace("+00:00", "Z")


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    try:
        os.write(descriptor, (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8"))
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.replace(temporary, path)


def job_allocation(job: dict[str, Any], config: NodeConfig) -> Allocation:
    try:
        cpus = parse_cpu_list(str(job.get("cpu_list") or ""))
    except ValueError:
        cpus = ()
    gpus = tuple(g for g in str(job.get("gpu_list") or "").split(",") if g)
    mem = job.get("mem_bytes")
    if mem is None:
        # A 0.5 job has no recorded memory request; charge the default of
        # the partition that owns its first core.
        owner = next(
            (part for part in config.partitions.values() if cpus and cpus[0] in part.cpus),
            None,
        )
        mem = owner.default_mem if owner is not None else 0
    return Allocation(
        cpus=tuple(cpus),
        gpus=gpus,
        gpu_claim=int(job.get("gpu_claim_bytes") or 0),
        gpu_exclusive=bool(job.get("gpu_exclusive")),
        mem=int(mem),
    )


def job_request(job: dict[str, Any], config: NodeConfig) -> Request:
    partition_name = job.get("partition") or config.default_partition
    partition = config.partitions.get(partition_name)
    cpu_request = None
    if job.get("cpu_request"):
        try:
            cpu_request = parse_cpu_list(str(job["cpu_request"]))
        except ValueError:
            cpu_request = ()
    mem = job.get("mem_bytes")
    time_limit = job.get("time_limit_s")
    if partition is not None:
        if mem is None:
            mem = partition.default_mem
        if time_limit is None:
            time_limit = partition.default_time
    max_cpus = job.get("max_borrow_cores")
    return Request(
        job_id=int(job["id"]),
        partition=str(partition_name),
        ncpus=int(job["requested_cores"]),
        cpu_list=cpu_request,
        gpus=int(job.get("gpus_requested") or 0),
        gpu_mem=job.get("gpu_mem_bytes"),
        mem=int(mem or 0),
        time_limit=None if time_limit is None else int(time_limit),
        blocked=bool(job.get("ensemble_launch_paused")),
        borrow=bool(job.get("borrow")),
        max_cpus=None if max_cpus is None else int(max_cpus),
        borrow_from=_borrow_from(job),
        exclusive_timing=bool(job.get("exclusive_timing")),
        on_preempt=str(job.get("on_preempt") or config.borrow.default_on_preempt),
    )


def _borrow_from(job: dict[str, Any]) -> tuple[str, ...] | None:
    raw = job.get("borrow_from")
    if not raw:
        return None
    try:
        value = json.loads(str(raw))
    except ValueError:
        return None
    if not isinstance(value, list):
        return None
    return tuple(str(item) for item in value)


def job_placement(job: dict[str, Any]) -> dict[str, Any]:
    try:
        value = json.loads(str(job.get("placement_json") or "{}"))
    except ValueError:
        return {}
    return value if isinstance(value, dict) else {}


def _job_cpus(job: dict[str, Any]) -> tuple[int, ...]:
    try:
        return parse_cpu_list(str(job.get("cpu_list") or ""))
    except ValueError:
        return ()


def borrower_entries(jobs: list[dict[str, Any]], config: NodeConfig) -> list[Borrower]:
    """Active borrow jobs (running, frozen or being stopped)."""
    entries: list[Borrower] = []
    for job in jobs:
        if job["state"] not in ACTIVE_STATES or not job.get("borrow"):
            continue
        placement = job_placement(job)
        items = [item for item in placement.get("cpus", []) if isinstance(item, dict) and "cpu" in item]
        started = _parse_time(job.get("started_at"))
        limit = job.get("time_limit_s")
        entries.append(
            Borrower(
                job_id=int(job["id"]),
                home=str(job.get("partition") or config.default_partition),
                cpus=_job_cpus(job),
                mem=int(job.get("mem_bytes") or 0),
                on_preempt=str(job.get("on_preempt") or config.borrow.default_on_preempt),
                placement={int(item["cpu"]): str(item.get("class") or "free") for item in items},
                owners={int(item["cpu"]): int(item["owner_job"]) for item in items
                        if item.get("owner_job") is not None},
                holder_keys={int(item["cpu"]): str(item["holder_key"]) for item in items if item.get("holder_key")},
                frozen_since=_parse_time(job.get("frozen_at")),
                pending=job.get("preempt_pending") or None,
                borrow_from=_borrow_from(job),
                deadline=None if started is None or limit is None else started + int(limit),
                order=int(job.get("queue_order") or job["id"]),
            )
        )
    return entries


def running_entries(jobs: list[dict[str, Any]], config: NodeConfig, now: float) -> list[Running]:
    """Owner jobs only: borrow jobs hold no reservation against owners."""
    entries: list[Running] = []
    for job in jobs:
        if job["state"] not in ACTIVE_STATES or job.get("borrow"):
            continue
        allocation = job_allocation(job, config)
        started = _parse_time(job.get("started_at"))
        limit = job.get("time_limit_s")
        if started is None or limit is None:
            end = INFINITY
        else:
            end = max(now, started + int(limit))
        entries.append(Running(int(job["id"]), allocation, end))
    return entries


def placement_environment(
    *,
    job: dict[str, Any],
    config: NodeConfig,
    partition: str,
    allocation: Allocation,
    time_limit: int | None,
    unit: str,
) -> dict[str, str]:
    cpu_list = format_cpu_list(allocation.cpus)
    gpus = ",".join(allocation.gpus)
    if not allocation.gpus:
        gpu_mode = "none"
    elif allocation.gpu_exclusive:
        gpu_mode = "exclusive"
    else:
        gpu_mode = "shared"
    return {
        "SIMQ_JOB_ID": str(job["id"]),
        "SIMQ_LAUNCH_TOKEN": str(job["launch_token"]),
        "SIMQ_PARTITION": partition,
        "SIMQ_NODE": config.node,
        "SIMQ_CPU_LIST": cpu_list,
        "SIMQ_CPU_COUNT": str(len(allocation.cpus)),
        "SIMQ_NCPUS": str(len(allocation.cpus)),
        "SIMQ_GPUS": gpus,
        "SIMQ_GPU_MODE": gpu_mode,
        "SIMQ_GPU_MEM": str(allocation.gpu_claim if allocation.gpus else 0),
        "SIMQ_MEM": str(allocation.mem),
        "SIMQ_TIME_LIMIT": "" if time_limit is None else str(time_limit),
        "SIMQ_LOG": str(job["log_path"]),
        "SIMQ_UNIT": unit,
        "CUDA_VISIBLE_DEVICES": gpus,
        "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
        ENSEMBLE_DRAIN_PROTOCOL_ENV: f"{ENSEMBLE_CONTROL_PROTOCOL}:{job['launch_token']}",
    }


def job_environment(job: dict[str, Any], placement: dict[str, str]) -> dict[str, str]:
    """The captured submission environment plus the placement variables.

    Inherited ``SLURM_*`` variables are removed (unless passed with ``-e``)
    because a stale ``SLURM_JOB_ID`` makes launchers such as pf-launch,
    Open MPI and MPICH take their Slurm code paths.
    """
    environment = {str(k): str(v) for k, v in dict(json.loads(str(job["env_json"]))).items()}
    try:
        explicit = set(json.loads(str(job.get("env_override_keys") or "[]")))
    except ValueError:
        explicit = set()
    for key in list(environment):
        if key.startswith("SLURM_") and key not in explicit:
            del environment[key]
        elif key.startswith("SIMQ_"):
            del environment[key]
    environment.update(placement)
    return environment


class NodeWorker:
    def __init__(
        self,
        config: NodeConfig,
        *,
        drain: bool = False,
        wait_for_lock: bool = False,
        once: bool = False,
        clock: Callable[[], float] = time.time,
    ) -> None:
        self.config = config
        self.store = QueueStore()
        self.drain = drain
        self.wait_for_lock = wait_for_lock
        self.once = once
        self.clock = clock
        self.stop_requested = False
        self.stop_signal: int | None = None
        self._lock_file: Any = None
        self._stop_sent: dict[int, float] = {}
        self._reasons: dict[int, str] = {}
        self._noted: set[tuple[int, str]] = set()
        self.job_root = state_dir() / "jobs"
        self.runtime_dir = str(Path(__file__).resolve().parents[2])
        self.last_started = 0
        self.telemetry = NodeTelemetry(config)
        self._telemetry_used = False
        self._slice_checked_at: float | None = None
        self._slice_error: str | None = None

    # ----- infrastructure ------------------------------------------------

    def log(self, message: str) -> None:
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

    def _note_once(self, job_id: int, key: str, message: str) -> None:
        if (job_id, key) not in self._noted:
            self._noted.add((job_id, key))
            self.log(message)

    # ----- reconciliation ------------------------------------------------

    def reconcile(self, *, startup: bool = False) -> None:
        jobs = self.store.active_jobs()
        unit_jobs = [j for j in jobs if j.get("launch_backend") == "systemd-user" and j.get("unit_name")]
        other_jobs = [j for j in jobs if j not in unit_jobs]
        try:
            states = show_units(self.config.launch, [str(j["unit_name"]) for j in unit_jobs])
        except UnitError as error:
            # Never finalize anything while systemd cannot be queried.
            self.log(f"cannot query job units: {error}")
            states = None
        if states is not None:
            for job in unit_jobs:
                try:
                    self._reconcile_unit_job(job, states[str(job["unit_name"])], startup=startup)
                except Exception as error:  # keep reconciling the others
                    self.log(f"reconcile error for job {job['id']}: {error!r}")
        for job in other_jobs:
            self._reconcile_legacy_job(job)

    def _job_files(self, job: dict[str, Any]) -> tuple[dict | None, dict | None]:
        job_dir = Path(str(job.get("job_dir") or ""))
        token = str(job["launch_token"])
        status = read_json(job_dir / "status.json") if job.get("job_dir") else None
        stop_post = read_json(job_dir / "stop-post.json") if job.get("job_dir") else None
        if status is not None and status.get("launch_token") != token:
            status = None
        if stop_post is not None and stop_post.get("launch_token") != token:
            stop_post = None
        return status, stop_post

    def _reconcile_unit_job(self, job: dict[str, Any], unit: UnitState, *, startup: bool) -> None:
        job_id = int(job["id"])
        status, stop_post = self._job_files(job)
        token = str(job["launch_token"])
        if unit.found and unit.active:
            if not unit.has_token(token):
                self._note_once(
                    job_id,
                    "foreign-unit",
                    f"job {job_id}: unit {unit.name} is active without this job's launch token; "
                    "not adopting it and keeping the allocation reserved",
                )
                return
            if job["state"] == "starting" or startup:
                self.store.mark_unit_running(job_id, pid=unit.main_pid or None, boot_id=_boot_id())
            if startup:
                self.store.mark_adopted(job_id)
                borrow_note = ""
                if job.get("borrow"):
                    borrow_note = ", borrow job" + (", frozen" if unit.frozen else "")
                self.log(
                    f"adopted job {job_id} (unit {unit.name}, main pid {unit.main_pid}, "
                    f"CPUs {job.get('cpu_list')}, GPUs {job.get('gpu_list') or 'none'}{borrow_note})"
                )
            if job.get("cancel_requested_at"):
                last = self._stop_sent.get(job_id)
                now = time.monotonic()
                if last is None or now - last > STOP_RESEND_SECONDS:
                    stop_unit(self.config.launch, unit.name)
                    self._stop_sent[job_id] = now
                    if last is None:
                        self.log(f"stopping unit {unit.name} for canceled job {job_id}")
                return
            if job.get("borrow"):
                self._reconcile_borrow_unit(job, unit)
            return
        self._finalize_unit_job(job, unit, status, stop_post)

    def _reconcile_borrow_unit(self, job: dict[str, Any], unit: UnitState) -> None:
        job_id = int(job["id"])
        if job.get("preempt_pending"):
            last = self._stop_sent.get(job_id)
            now = time.monotonic()
            if last is None or now - last > STOP_RESEND_SECONDS:
                stop_unit(self.config.launch, unit.name)
                self._stop_sent[job_id] = now
            return
        if not unit.freezer_state:
            return
        if job.get("frozen_at") and unit.freezer_state in ("running", "thawing"):
            # The recorded intent is frozen (its cores were given to an owner).
            if freeze_unit(self.config.launch, unit.name):
                self.log(f"job {job_id}: unit {unit.name} was running although recorded frozen; froze it again")
        elif not job.get("frozen_at") and unit.frozen:
            if self.store.mark_frozen_outside(job_id):
                self.log(f"job {job_id}: unit {unit.name} was frozen outside the scheduler; recorded as frozen")

    def _append_job_log(self, job: dict[str, Any], text: str) -> None:
        try:
            with open(str(job["log_path"]), "ab") as handle:
                handle.write(text.encode("utf-8", errors="replace"))
        except OSError:
            pass

    def _finish(
        self,
        job: dict[str, Any],
        *,
        state: str,
        exit_code: int | None,
        message: str | None,
        result: str | None,
        finished_at: str | None = None,
    ) -> None:
        job_id = int(job["id"])
        if self.store.finish(
            job_id,
            state=state,
            exit_code=exit_code,
            message=message,
            result=result,
            finished_at=finished_at,
        ):
            self._append_job_log(job, f"state:     {state}\n===== simq job {job_id} finished =====\n")
            self.log(
                f"job {job_id} {state}"
                + (f" (exit {exit_code})" if exit_code is not None else "")
                + (f": {message}" if message else "")
            )
        self._stop_sent.pop(job_id, None)

    def _finalize_unit_job(
        self,
        job: dict[str, Any],
        unit: UnitState,
        status: dict | None,
        stop_post: dict | None,
    ) -> None:
        job_id = int(job["id"])
        cancel = bool(job.get("cancel_requested_at"))
        service_result = (stop_post or {}).get("service_result") or (unit.result if unit.found else None)
        pending = job.get("preempt_pending")
        if pending and not cancel:
            done = (
                status is not None
                and status.get("phase") == "finished"
                and status.get("returncode") == 0
                and status.get("command_started", True)
            )
            if not done:
                self._finish_preempted(job, str(pending), status, stop_post, service_result)
                return
        if status is not None and status.get("phase") == "refused":
            reason = str(status.get("reason") or "launch refused")
            if status.get("refusal") == LOCK_BUSY_REFUSAL and not cancel:
                if self.store.requeue_unstarted(job_id, f"lease busy at launch; requeued ({reason})"):
                    self.log(f"job {job_id}: {reason}; command not started, requeued")
                return
            if cancel:
                self._finish(job, state="canceled", exit_code=None, message="canceled before the command started",
                             result=service_result)
            else:
                self._finish(job, state="failed", exit_code=None, message=f"launch refused: {reason}",
                             result=service_result)
            return
        if status is None and stop_post is None and not unit.found:
            if job["state"] == "starting" or (job["state"] == "cancelling" and not job.get("pid")):
                # The wrapper never wrote anything: the command never ran.
                if self.store.requeue_unstarted(job_id, "unit did not start; requeued"):
                    self.log(f"job {job_id}: unit {unit.name} never started; requeued or canceled")
                return
            self._finish(job, state="interrupted", exit_code=None,
                         message="job unit vanished without a status record; check outputs and retry explicitly",
                         result=None)
            return
        finished_at = None
        if status is not None and status.get("finished_unix"):
            finished_at = _iso(float(status["finished_unix"]))
        elif stop_post is not None and stop_post.get("recorded_unix"):
            finished_at = _iso(float(stop_post["recorded_unix"]))
        phase = (status or {}).get("phase")
        returncode = (status or {}).get("returncode")
        if phase == "finished" and not status.get("command_started", True):
            self._finish(job, state="canceled" if cancel else "failed", exit_code=None,
                         message=f"launch failed: {status.get('launch_error')}", result=service_result,
                         finished_at=finished_at)
            return
        if service_result == "timeout":
            limit = job.get("time_limit_s")
            self._finish(job, state="failed", exit_code=returncode,
                         message=f"time limit exceeded ({format_time_limit(limit)})", result="timeout",
                         finished_at=finished_at)
            return
        if service_result == "oom-kill":
            self._finish(job, state="failed", exit_code=returncode,
                         message="out of memory: the job unit's memory limit was reached", result="oom-kill",
                         finished_at=finished_at)
            return
        if phase == "finished":
            if cancel:
                self._finish(job, state="canceled", exit_code=returncode, message="canceled by user",
                             result=service_result, finished_at=finished_at)
            elif returncode == 0:
                self._finish(job, state="succeeded", exit_code=0, message=None, result=service_result,
                             finished_at=finished_at)
            elif returncode == ENSEMBLE_DRAIN_EXIT_CODE and self.store.requeue_drained_ensemble(job_id):
                control = self.store.ensemble_control_snapshot(job_id)
                self._append_job_log(
                    job,
                    f"suspended: {utc_now()}\nremaining: {control['queued']} queued member(s)\n"
                    f"===== simq ensemble {job_id} allocation released =====\n",
                )
                self.log(f"suspended ensemble {job_id}; {control['queued']} queued member(s) remain")
            else:
                self._finish(job, state="failed", exit_code=returncode,
                             message=f"command exited with status {returncode}", result=service_result,
                             finished_at=finished_at)
            return
        # The wrapper did not record a result: it was killed with its unit.
        if status is not None and status.get("boot_id") not in (None, _boot_id()) and stop_post is None:
            self._finish(job, state="interrupted", exit_code=None,
                         message="machine restarted while the job ran; retry explicitly if appropriate",
                         result=None)
            return
        if cancel:
            self._finish(job, state="canceled", exit_code=returncode, message="canceled by user",
                         result=service_result, finished_at=finished_at)
            return
        detail = "unit ended before the wrapper recorded a result"
        if stop_post is not None:
            detail += (
                f" (result={stop_post.get('service_result')}, code={stop_post.get('exit_code')}, "
                f"status={stop_post.get('exit_status')})"
            )
            self._finish(job, state="failed", exit_code=returncode, message=detail, result=service_result,
                         finished_at=finished_at)
        else:
            self._finish(job, state="interrupted", exit_code=returncode,
                         message=detail + "; retry explicitly if appropriate", result=service_result)

    def _finish_preempted(
        self,
        job: dict[str, Any],
        mode: str,
        status: dict | None,
        stop_post: dict | None,
        service_result: str | None,
    ) -> None:
        job_id = int(job["id"])
        reason = str(job.get("preempt_reason") or "preempted")
        returncode = (status or {}).get("returncode")
        self._stop_sent.pop(job_id, None)
        if mode == "cancel":
            finished_at = None
            if status is not None and status.get("finished_unix"):
                finished_at = _iso(float(status["finished_unix"]))
            elif stop_post is not None and stop_post.get("recorded_unix"):
                finished_at = _iso(float(stop_post["recorded_unix"]))
            self._finish(job, state="canceled", exit_code=returncode,
                         message=f"preempted: {reason}; canceled (--on-preempt cancel)",
                         result="preempted", finished_at=finished_at)
            return
        detail = {"returncode": returncode, "service_result": service_result}
        if self.store.requeue_preempted(job_id, f"pending: preempted ({reason}); requeued", detail=detail):
            self._append_job_log(
                job,
                f"simq: preempted ({reason}); requeued at {utc_now()}\n"
                f"===== simq job {job_id} preempted; requeued under the same ID =====\n",
            )
            self.log(f"job {job_id}: preempted ({reason}); requeued under the same ID")

    def _reconcile_legacy_job(self, job: dict[str, Any]) -> None:
        """A job started by a 0.5 worker as its child process.

        It cannot be adopted (its exit status belongs to the dead 0.5
        worker), but it is never killed implicitly: its cores stay reserved
        until its processes are gone, then it is marked interrupted.
        """
        from .worker import _verified_job_process_groups

        job_id = int(job["id"])
        groups = _verified_job_process_groups(job)
        groups.discard(os.getpgrp())
        if groups:
            self._note_once(
                job_id,
                "legacy",
                f"job {job_id} was started by a 0.5 worker and is still running "
                f"(process groups {sorted(groups)}); observing it and keeping CPUs "
                f"{job.get('cpu_list')} reserved",
            )
            if job.get("cancel_requested_at"):
                now = time.monotonic()
                first = self._stop_sent.setdefault(job_id, now)
                sig = signal.SIGKILL if now - first > self.config.launch.kill_grace_seconds else signal.SIGTERM
                for group in groups:
                    try:
                        os.killpg(group, sig)
                    except ProcessLookupError:
                        pass
            return
        if job.get("cancel_requested_at"):
            self._finish(job, state="canceled", exit_code=None, message="canceled by user", result=None)
            return
        self.store.mark_interrupted(
            job_id,
            "started by a 0.5 worker that is no longer supervising it; exit status unknown; "
            "check outputs and retry explicitly if appropriate",
        )
        self.log(f"marked 0.5 job {job_id} interrupted (no live processes)")

    # ----- dispatch --------------------------------------------------------

    def dispatch(self) -> int:
        jobs = self.store.node_jobs()
        if any(job.get("borrow") for job in jobs):
            return self._dispatch_with_borrowers(jobs)
        if self._telemetry_used:
            self.telemetry.reset()
            self._telemetry_used = False
        # Without borrow jobs: exactly the 0.6.0 dispatch.
        if self.store.is_paused():
            return 0
        queued = [j for j in jobs if j["state"] == "queued"]
        if not queued:
            return 0
        now = self.clock()
        active = [j for j in jobs if j["state"] in ACTIVE_STATES]
        running = running_entries(active, self.config, now)
        allocated_cpus = {cpu for item in running for cpu in item.allocation.cpus}
        allocated_gpus = {gpu for item in running for gpu in item.allocation.gpus}
        busy_cpus, busy_gpus = busy_resources(
            self.config, allocated_cpus=allocated_cpus, allocated_gpus=allocated_gpus
        )
        pool = build_pool(self.config, running, external_cpus=busy_cpus, external_gpus=busy_gpus)
        requests = [job_request(j, self.config) for j in queued]
        decision = plan(self.config, pool, running, requests, now)
        by_id = {int(j["id"]): j for j in queued}
        started = 0
        for request, allocation in decision.starts:
            if self.stop_requested:
                break
            if self._launch(by_id[request.job_id], request, allocation, busy_cpus=set(busy_cpus)):
                started += 1
        self._record_reasons(decision.skipped, decision.reservation, busy_cpus)
        self.last_started = started
        return started

    def _record_reasons(
        self,
        skipped: dict[int, str],
        reservation: tuple[int, float] | None,
        busy_cpus: dict[int, list],
    ) -> None:
        for job_id, reason in skipped.items():
            if reason == "resources busy":
                holders = [
                    f"CPU {cpu}: {process_description(h.pid)}"
                    for cpu, found in busy_cpus.items()
                    for h in found[:1]
                ]
                if holders:
                    reason = "resources busy; lease held outside the queue (" + "; ".join(holders) + ")"
            if reservation and reservation[0] == job_id:
                shadow = reservation[1]
                if shadow != INFINITY:
                    reason += f"; reserved to start by {_iso(shadow)}"
            if self._reasons.get(job_id) != reason:
                self._reasons[job_id] = reason
                self.store.set_message(job_id, f"pending: {reason}")

    # ----- borrowing (0.6.1) ---------------------------------------------

    def _dispatch_with_borrowers(self, jobs: list[dict[str, Any]]) -> int:
        config = self.config
        paused = self.store.is_paused()
        now = self.clock()
        active = [j for j in jobs if j["state"] in ACTIVE_STATES]
        queued = [j for j in jobs if j["state"] == "queued"]
        owners = [j for j in active if not j.get("borrow")]
        borrow_active = [j for j in active if j.get("borrow")]
        running = running_entries(owners, config, now)
        owner_of = {cpu: int(item.job_id) for item in running if item.job_id is not None
                    for cpu in item.allocation.cpus}
        allocated_gpus = {gpu for item in running for gpu in item.allocation.gpus}
        busy_cpus, busy_gpus = busy_resources(config, allocated_cpus=set(owner_of), allocated_gpus=allocated_gpus)
        pool = build_pool(config, running, external_cpus=busy_cpus, external_gpus=busy_gpus)
        borrowers = borrower_entries(borrow_active, config)
        self._telemetry_used = True
        cores = self.telemetry.sample(owners=owners, external=busy_cpus, borrowers=borrowers)
        exclusive = {cpu for job in owners if job.get("exclusive_timing") for cpu in _job_cpus(job)}
        requests = [job_request(j, config) for j in queued]
        decision = plan_node(
            config, pool, running, requests, now,
            borrowers=borrowers, cores=cores, owner_of=owner_of,
            exclusive_owner_cpus=exclusive, starts_allowed=not paused,
        )
        active_by_id = {int(j["id"]): j for j in borrow_active}
        for item in decision.preemptions:
            if item.job_id in active_by_id:
                self._preempt(active_by_id[item.job_id], item)
        for job_id, placement in decision.resumes:
            if job_id in active_by_id:
                self._resume(active_by_id[job_id], placement, owner_of, cores)
        by_id = {int(j["id"]): j for j in queued}
        started = 0
        skipped = dict(decision.skipped)
        for request, allocation in decision.starts:
            if self.stop_requested:
                break
            busy = set(busy_cpus)
            if allocation.borrow:
                error = self._borrow_slice_ready()
                if error is not None:
                    skipped[request.job_id] = f"borrow slice unavailable: {error}"
                    continue
                busy -= {cpu for cpu, _, cls in allocation.placement if cls == "leased-idle"}
            try:
                launched = self._launch(by_id[request.job_id], request, allocation, busy_cpus=busy,
                                        owner_of=owner_of, cores=cores)
            except RuntimeError as error:
                # The claim re-checks every active job: a borrow job whose
                # preemption did not take effect keeps its cores until the
                # next pass.
                self.log(f"job {request.job_id}: not started this pass: {error}")
                skipped[request.job_id] = "waiting for a preempted borrow job"
                continue
            if launched:
                started += 1
        self._record_reasons(skipped, decision.reservation, busy_cpus)
        self.last_started = started
        return started

    def _borrow_slice_ready(self) -> str | None:
        now = time.monotonic()
        if self._slice_checked_at is not None and (self._slice_error is None or now - self._slice_checked_at < 60):
            return self._slice_error
        self._slice_checked_at = now
        self._slice_error = ensure_borrow_slice(self.config.launch, self.config.borrow)
        if self._slice_error is None:
            self.log(f"borrow slice {self.config.borrow.slice}: CPUWeight={self.config.borrow.slice_cpu_weight}")
        else:
            self.log(f"borrow slice not ready: {self._slice_error}")
        return self._slice_error

    def _placement_document(
        self,
        home: str,
        on_preempt: str,
        placement: tuple[tuple[int, str, str], ...],
        owner_of: dict[int, int],
        cores: dict[int, Any],
    ) -> dict[str, Any]:
        items = []
        for cpu, partition, cls in placement:
            entry: dict[str, Any] = {"cpu": cpu, "partition": partition, "class": cls}
            if cls == "owner-idle":
                entry["owner_job"] = owner_of.get(cpu)
            if cls == "leased-idle" and cpu in cores:
                entry["holder_key"] = cores[cpu].holder_key
                entry["holder"] = cores[cpu].holder
            items.append(entry)
        return {"home": home, "on_preempt": on_preempt, "cpus": items}

    def _preempt(self, job: dict[str, Any], item: Preemption) -> None:
        job_id = int(job["id"])
        unit = str(job.get("unit_name") or "")
        reason = item.reason
        mode = item.mode
        escalated = item.escalated
        if mode == "freeze":
            if not self.store.mark_preempted(job_id, mode="freeze", reason=reason, by_job=item.by_job):
                return
            if unit and freeze_unit(self.config.launch, unit):
                self.log(f"froze borrow job {job_id} ({unit}): {reason}")
                self._append_job_log(job, f"simq: frozen at {utc_now()}: {reason}\n")
                return
            reason += "; freezing failed, stopping instead"
            mode, escalated = "requeue", True
        if not self.store.mark_preempted(job_id, mode=mode, reason=reason, by_job=item.by_job, escalated=escalated):
            return
        if unit:
            stop_unit(self.config.launch, unit)
            self._stop_sent[job_id] = time.monotonic()
        self.log(f"preempting borrow job {job_id} ({mode}{', escalated from freeze' if escalated else ''}): {reason}")
        self._append_job_log(job, f"simq: preempted at {utc_now()} ({mode}): {reason}\n")

    def _resume(
        self,
        job: dict[str, Any],
        placement: tuple[tuple[int, str, str], ...],
        owner_of: dict[int, int],
        cores: dict[int, Any],
    ) -> None:
        job_id = int(job["id"])
        unit = str(job.get("unit_name") or "")
        if not unit or not thaw_unit(self.config.launch, unit):
            self._note_once(job_id, f"thaw-failed-{job.get('frozen_at')}", f"job {job_id}: cannot thaw {unit}")
            return
        home = str(job.get("partition") or self.config.default_partition)
        on_preempt = str(job.get("on_preempt") or self.config.borrow.default_on_preempt)
        document = self._placement_document(home, on_preempt, placement, owner_of, cores)
        self.store.mark_resumed(job_id, document)
        self.log(f"resumed borrow job {job_id} ({unit}) on CPUs {job.get('cpu_list')}")
        self._append_job_log(job, f"simq: resumed at {utc_now()}\n")

    def _launch(
        self,
        job: dict[str, Any],
        request: Request,
        allocation: Allocation,
        *,
        busy_cpus: set[int],
        owner_of: dict[int, int] | None = None,
        cores: dict[int, Any] | None = None,
    ) -> bool:
        config = self.config
        partition = config.partitions[request.partition]
        borrow = allocation.borrow
        placement_doc: dict[str, Any] | None = None
        owner_idle: set[int] = set()
        if borrow:
            locks = borrow_lock_plan(config, allocation.cpus)
            placement_doc = self._placement_document(
                partition.name, request.on_preempt, allocation.placement, owner_of or {}, cores or {}
            )
            owner_idle = {cpu for cpu, _, cls in allocation.placement if cls == "owner-idle"}
        else:
            locks = lock_plan_for(config, allocation.cpus, allocation.gpus, gpu_shared=not allocation.gpu_exclusive)
        mem = allocation.mem or partition.default_mem
        swap = partition.swap_max
        claimed = self.store.claim_allocation(
            int(job["id"]),
            cpu_list=format_cpu_list(allocation.cpus),
            gpu_list=",".join(allocation.gpus),
            gpu_claim_bytes=allocation.gpu_claim,
            gpu_exclusive=allocation.gpu_exclusive,
            mem_bytes=mem,
            swap_bytes=swap,
            time_limit_s=request.time_limit,
            partition=partition.name,
            unit_prefix=config.launch.unit_prefix,
            job_root=self.job_root,
            lock_paths=[{"path": lock.path, "mode": lock.mode, "resource": lock.resource} for lock in locks],
            runtime_dir=self.runtime_dir,
            busy_cpus=busy_cpus,
            **({"borrow": True, "placement": placement_doc, "owner_idle_cpus": owner_idle} if borrow else {}),
        )
        if claimed is None:
            return False
        job_id = int(claimed["id"])
        self._reasons.pop(job_id, None)
        unit = str(claimed["unit_name"])
        job_dir = Path(str(claimed["job_dir"]))
        try:
            job_dir.mkdir(mode=0o700, parents=True, exist_ok=False)
            self.job_root.chmod(0o700)
            placement = placement_environment(
                job=claimed,
                config=config,
                partition=partition.name,
                allocation=Allocation(allocation.cpus, allocation.gpus, allocation.gpu_claim,
                                      allocation.gpu_exclusive, mem),
                time_limit=request.time_limit,
                unit=unit,
            )
            if borrow:
                placement["SIMQ_BORROW"] = "1"
                placement["SIMQ_ON_PREEMPT"] = request.on_preempt
            environment = job_environment(claimed, placement)
            argv = json.loads(str(claimed["argv_json"]))
            if not isinstance(argv, list) or not argv or not all(isinstance(a, str) for a in argv):
                raise ValueError("stored argv is not a nonempty string list")
            gpu_text = "none"
            if allocation.gpus:
                gpu_text = (
                    f"{','.join(allocation.gpus)} ("
                    + ("exclusive" if allocation.gpu_exclusive else f"shared, claim {format_size(allocation.gpu_claim)}")
                    + ")"
                )
            spec = {
                "schema": SPEC_SCHEMA,
                "job_id": job_id,
                "name": str(claimed["name"]),
                "launch_token": str(claimed["launch_token"]),
                "unit": unit,
                "argv": argv,
                "env": environment,
                "cwd": str(claimed["cwd"]),
                "project_root": str(project_root()),
                "log_path": str(claimed["log_path"]),
                "status_path": str(job_dir / "status.json"),
                "stop_post_path": str(job_dir / "stop-post.json"),
                "cpu_list": placement["SIMQ_CPU_LIST"],
                "partition": partition.name,
                "gpu_text": gpu_text,
                "mem_text": f"{format_size(mem)} (swap {format_size(swap)})",
                "time_text": format_time_limit(request.time_limit),
                "locks": [{"path": lock.path, "mode": lock.mode, "resource": lock.resource} for lock in locks],
                "submitted_at": str(claimed["created_at"]),
                "scheduler_version": __version__,
                "scheduler_runtime": self.runtime_dir,
                "require_unit": True,
            }
            if borrow:
                spec["borrow"] = True
                spec["borrow_text"] = (
                    "; ".join(f"CPU {cpu} from {part} ({cls})" for cpu, part, cls in allocation.placement)
                    + f"; on preempt: {request.on_preempt}"
                )
            spec_path = job_dir / "spec.json"
            _atomic_json(spec_path, spec)
            unit_env = {key: value for key, value in placement.items() if key.startswith("SIMQ_")}
            unit_env["CUDA_VISIBLE_DEVICES"] = placement["CUDA_VISIBLE_DEVICES"]
            runtime_max = request.time_limit if partition.enforce_time_limit else None
            command = build_run_command(
                config.launch,
                unit=unit,
                description=f"simq job {job_id}: {claimed['name']}",
                cwd=str(claimed["cwd"]),
                cpus=allocation.cpus,
                mem_bytes=mem,
                swap_bytes=swap,
                runtime_max_s=runtime_max,
                nice=partition.nice,
                log_path=str(claimed["log_path"]),
                spec_path=str(spec_path),
                setenv=unit_env,
                wrapper=JOBWRAP,
                borrow=config.borrow if borrow else None,
            )
            start_unit(command)
        except (UnitError, OSError, ValueError) as error:
            # The unit may exist if systemd-run failed after starting it.
            try:
                state = show_units(config.launch, [unit])[unit]
            except UnitError:
                state = None
            if state is not None and state.active:
                self.log(f"job {job_id}: systemd-run reported '{error}' but {unit} is active; adopting")
                return True
            self._finish(claimed, state="failed", exit_code=None, message=f"launch failed: {error}",
                         result=None)
            return False
        self.log(
            f"started {'borrow ' if borrow else ''}job {job_id} ({claimed['name']}) on partition {partition.name}: "
            f"CPUs {placement['SIMQ_CPU_LIST']}, GPUs {placement['SIMQ_GPUS'] or 'none'}, "
            f"mem {format_size(mem)}, time {format_time_limit(request.time_limit)}, unit {unit}"
            + (f"; {spec.get('borrow_text')}" if borrow else "")
        )
        return True

    # ----- main loop -------------------------------------------------------

    def run(self) -> int:
        if not self._acquire_lock():
            self.log("another worker already owns the queue lock")
            return 0
        signal.signal(signal.SIGTERM, self._handle_stop)
        signal.signal(signal.SIGINT, self._handle_stop)
        self.job_root.mkdir(mode=0o700, parents=True, exist_ok=True)
        self.log(
            f"node worker started; simqueue {__version__} from {self.runtime_dir}; "
            f"config {self.config.path}; partitions "
            + ", ".join(
                f"{p.name}={format_cpu_list(p.cpus)}" + (f"+{len(p.gpus)}gpu" if p.gpus else "")
                for p in self.config.partitions.values()
            )
            + f"; drain={self.drain}"
        )
        exit_code = 0
        failures = 0
        startup = True
        while not self.stop_requested:
            try:
                self.reconcile(startup=startup)
                startup = False
                started = self.dispatch()
                failures = 0
            except Exception as error:  # one bad pass must not end the scheduler
                failures += 1
                self.log(f"node worker pass failed ({failures} in a row): {error!r}")
                if failures >= 30:
                    exit_code = 2
                    break
                started = 0
            if self.once:
                break
            if self.drain and failures == 0 and started == 0 and not self.store.active_jobs():
                break
            deadline = time.monotonic() + self.config.launch.poll_seconds
            while not self.stop_requested and time.monotonic() < deadline:
                time.sleep(min(0.1, self.config.launch.poll_seconds))
        remaining = len(self.store.active_jobs())
        suffix = f" after signal {self.stop_signal}" if self.stop_signal else ""
        self.log(f"node worker stopped{suffix}; {remaining} job unit(s) continue independently")
        self._release_lock()
        return exit_code
