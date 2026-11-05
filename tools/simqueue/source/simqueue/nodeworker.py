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

0.6.2 adds the GPU line (``gpuqueue.py``, ``gpuline.py``): a queued GPU job
that is ready except for the GPU joins the line through a claim unit that
wins the GPU lease files in the kernel queue and hands them to the job unit,
so raw ``flock`` users can no longer take the GPU from it at every release.

0.6.3 adds partition aliases (``nodeconfig.py``): a job recorded under an
alias (for example a partition name from before a node became one pool) is
scheduled, adopted and lent under the partition the alias resolves to, and
keeps its recorded name in the database and in ``SIMQ_PARTITION``.

0.7.0 makes the worker event-driven and memory-aware (``memory.py``):

* It sleeps until something happens: a client writes a byte to the wake
  FIFO ``<state>/worker.wake`` after a submission or cancellation, a job
  unit's stop-post hook does the same when the job ends, and the queue
  database's ``data_version`` changes (any other writer).  Between passes it
  only reads files: each running job's cgroup (``memory.current``,
  ``memory.peak``) every ``memory.sample_seconds`` and its existence, which
  is how the end of a job started by an older wrapper is noticed.  A pass
  queries systemd only for jobs that need it (starting, being stopped,
  unit gone, freezer state differs); every ``SAFETY_RECONCILE_IDLE``
  idle intervals it queries every unit.  Passes run every
  ``launch.poll_seconds`` only while something must be watched (a job
  starting or being stopped, borrow jobs, the GPU line, queued jobs held by
  leases outside the queue or by the memory guard), else at most every
  ``launch.idle_seconds``.
* Memory: samples are kept per job as min/max buckets and written to the
  database every ``memory.record_seconds``; at a job's end the unit's own
  record (``memory.peak``, CPU time) is stored.  A start needs its request
  below MemAvailable minus the margin (admission guard).  When MemAvailable
  falls below the high-water mark, dispatch stops, borrow jobs are
  preempted by their ``--on-preempt`` policy and the jobs that grew are
  recorded; owner jobs are never stopped for it (their own MemoryMax bounds
  them).  Dispatch resumes at the low-water mark.
"""

from __future__ import annotations

import fcntl
import json
import os
import re
import select
import signal
import stat
import sys
import time
from collections import deque
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from types import FrameType
from typing import Any, Callable

from . import __version__
from .config import format_cpu_list, parse_cpu_list, project_root, state_dir
from .borrow import Borrower, Preemption, plan_node
from .gpuqueue import GpuQueue, GpuView
from .locks import borrow_lock_plan, busy_resources, lock_plan_for, process_description
from .nodeconfig import NodeConfig, format_size, format_time_limit
from .scheduler import INFINITY, Allocation, Request, Running, build_pool, plan
from .telemetry import NodeTelemetry
from .github import FINAL, Poster, build_task, job_phase
from .memory import (
    NodeMemory,
    Series,
    SpikeState,
    admission,
    admission_margin,
    cgroup_from_control_group,
    format_gib,
    read_node_memory,
    read_usage,
    unit_cgroup_path,
    user_manager_cgroup,
    water_marks,
)
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
    reset_failed,
    show_units,
    start_unit,
    stop_unit,
    thaw_unit,
)

SHM_DIR = Path("/dev/shm")
# Shared-memory segments that MPI libraries leave behind when a job is killed
# (Open MPI sm/vader BTL, PSM2, MPICH/Hydra, UCX); charged to the dead job's
# cgroup as shmem until deleted (job units have no swap).
SHM_PATTERNS = re.compile(r"^(sm_segment\.|vader_segment\.|psm2_shm\.|psm_shm\.|ompi\.|open_mpi|openmpi|"
                          r"mpich_shm|mpich2_shm|ucx_shm|ucx-|shm_mpi)")


SPEC_SCHEMA = "simqueue-job-spec-v1"
LOCK_BUSY_REFUSAL = "lock-busy"
STOP_RESEND_SECONDS = 30.0
WAKE_FIFO = "worker.wake"
EXIT_WATCH_SECONDS = 5.0  # after a wake-up, look for ended job units this long ...
EXIT_WATCH_STEP = 0.2  # ... this often (a stop-post hook wakes the worker just before its unit ends)
IDLE_CHECK_SECONDS = 5.0  # with no job running: how often data_version is compared
SAFETY_RECONCILE_IDLE = 10  # every N idle intervals every job unit is queried from systemd
GROWTH_WINDOW_SECONDS = 60.0  # memory growth is measured over this window
PRESSURE_META = "memory_pressure"


@dataclass
class JobTrack:
    """Memory of one running job unit, sampled from its cgroup."""

    job_id: int
    name: str
    unit: str
    cgroup: Path
    request: int
    borrow: bool
    series: Series
    peak: int = 0
    current: int = 0
    cpu_usec: int | None = None
    recent: deque = field(default_factory=lambda: deque(maxlen=64))  # (unix, current)
    near_limit_noted: bool = False
    dirty: bool = False
    gone: bool = False

    def observe(self, unix: float, current: int, peak: int | None) -> None:
        self.current = current
        self.peak = max(self.peak, current, peak or 0)
        self.series.add(unix, current)
        self.recent.append((unix, current))
        self.dirty = True

    def growth(self, unix: float, window: float = GROWTH_WINDOW_SECONDS) -> int:
        """How much the job grew over the last *window* seconds."""
        values = [value for at, value in self.recent if at >= unix - window]
        return self.current - min(values) if values else 0


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


def _shm_in_use(directory: str) -> set[str]:
    """Paths under *directory* that some process maps or holds open."""
    found: set[str] = set()
    prefix = directory.rstrip("/") + "/"
    for pid in os.listdir("/proc"):
        if not pid.isdigit():
            continue
        try:
            with open(f"/proc/{pid}/maps", "rb") as handle:
                for line in handle.read().decode("utf-8", "replace").splitlines():
                    at = line.find(prefix)
                    if at >= 0:
                        found.add(line[at:].split(" (deleted)")[0].strip())
        except OSError:
            pass
        try:
            for fd in os.listdir(f"/proc/{pid}/fd"):
                try:
                    target = os.readlink(f"/proc/{pid}/fd/{fd}")
                except OSError:
                    continue
                if target.startswith(prefix):
                    found.add(target)
        except OSError:
            pass
    return found


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
    # The recorded name may be an alias (0.6.3); its defaults apply to a job
    # that has no recorded memory or time limit.
    partition_name = job.get("partition") or config.default_partition
    partition = config.find_partition(partition_name)
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
                home=config.canonical_partition(job.get("partition")),
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


def running_entries(jobs: list[dict[str, Any]], config: NodeConfig, now: float,
                    currents: dict[int, int] | None = None) -> list[Running]:
    """Owner jobs only: borrow jobs hold no reservation against owners.
    *currents* (0.7.1): what each job uses now (its end frees that much)."""
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
        entries.append(Running(int(job["id"]), allocation, end,
                               None if currents is None else currents.get(int(job["id"]))))
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
    from .envpolicy import REDACTED, SESSION, read_secrets

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
        elif environment[key] == REDACTED:
            del environment[key]  # 0.7.0: a redacted value never reaches the job
    # 0.7.0: secrets given with -e, from their private file (not the database)
    environment.update(read_secrets(state_dir(), int(job["id"]), str(job["launch_token"])))
    # session plumbing of the user manager, so that `systemctl --user` works in a job
    for key in SESSION:
        if key not in environment and os.environ.get(key):
            environment[key] = os.environ[key]
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
        # one persistent connection: its own writes leave data_version alone
        self.store = QueueStore(persistent=True)
        self.drain = drain
        self.wait_for_lock = wait_for_lock
        self.once = once
        self.clock = clock
        self.stop_requested = False
        self.stop_signal: int | None = None
        self._lock_file: Any = None
        self._stop_sent: dict[int, float] = {}
        self._reasons: dict[int, str] = {}
        self._reason_written: dict[int, float] = {}
        self._noted: set[tuple[int, str]] = set()
        self.job_root = state_dir() / "jobs"
        self.runtime_dir = str(Path(__file__).resolve().parents[2])
        self.last_started = 0
        self.telemetry = NodeTelemetry(config)
        self._telemetry_used = False
        self._slice_checked_at: float | None = None
        self._slice_error: str | None = None
        self.gpu = GpuQueue(config, self.store, self.log, runtime_dir=self.runtime_dir)
        # 0.7.0: event loop and memory
        self.memory = config.memory
        self.wake_path = state_dir() / WAKE_FIFO
        self._wake_fd: int | None = None
        self._wake_full = True
        self._next_full = 0.0
        self._next_sample = 0.0
        self._next_record = 0.0
        self._next_check = 0.0
        self._watch_until = 0.0
        self._next_safety = 0.0
        self._data_version: int | None = None
        self._fast_reasons: list[str] = []
        self._external_busy = False
        self._memory_held: dict[int, str] = {}
        self._verified: dict[int, Path] = {}  # job id -> cgroup of its unit (verified by systemd)
        self._tracks: dict[int, JobTrack] = {}
        self._pressure = SpikeState()
        self._pressure_preempted: dict[int, float] = {}  # borrow job -> when it was frozen for pressure
        self._node: NodeMemory | None = None
        self._manager_cgroup = user_manager_cgroup()
        self._idle_seen: deque = deque()  # (unix, MemAvailable + the queue jobs' use), the last hour
        self._last_unit: dict[int, UnitState] = {}  # job id -> its unit's state when it ended
        self.shm_dir = SHM_DIR
        self._admissible: int | None = None
        self.passes = 0
        self.unit_queries = 0
        # 0.7.0: GitHub reporting (posting thread started on first use)
        self.github = Poster(config.github, self._github_notify) if config.github.enabled else None
        self._github_next: float | None = None  # wall clock of the earliest retry
        self.github_posts = 0

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
        if self._wake_fd is not None:
            try:
                os.write(self._wake_fd, b"s")
            except OSError:
                pass

    def _note_once(self, job_id: int, key: str, message: str) -> None:
        if (job_id, key) not in self._noted:
            self._noted.add((job_id, key))
            self.log(message)

    # ----- reconciliation ------------------------------------------------

    def reconcile(self, *, startup: bool = False, full: bool = True) -> None:
        """Bring the database in line with the job units.

        With *full* false (0.7.0 event loop) only the units that need it are
        queried from systemd: jobs starting or being stopped, jobs whose
        cgroup is gone or was never verified, borrow jobs whose freezer
        state differs from the record.  A running job whose cgroup exists
        needs nothing."""
        jobs = self.store.active_jobs()
        unit_jobs = [j for j in jobs if j.get("launch_backend") == "systemd-user" and j.get("unit_name")]
        other_jobs = [j for j in jobs if j not in unit_jobs]
        active_ids = {int(j["id"]) for j in jobs}
        for job_id in list(self._verified):
            if job_id not in active_ids:
                self._verified.pop(job_id, None)
        query = unit_jobs if (startup or full) else [j for j in unit_jobs if self._needs_query(j)]
        states = None
        if query:
            try:
                self.unit_queries += 1
                states = show_units(self.config.launch, [str(j["unit_name"]) for j in query])
            except UnitError as error:
                # Never finalize anything while systemd cannot be queried.
                self.log(f"cannot query job units: {error}")
                states = None
        if states is not None:
            for job in query:
                unit = states[str(job["unit_name"])]
                try:
                    if unit.found and unit.active and unit.has_token(str(job["launch_token"])):
                        self._verified[int(job["id"])] = (
                            cgroup_from_control_group(unit.control_group) or self._computed_cgroup(job)
                        )
                    self._reconcile_unit_job(job, unit, startup=startup)
                except Exception as error:  # keep reconciling the others
                    self.log(f"reconcile error for job {job['id']}: {error!r}")
        for job in other_jobs:
            self._reconcile_legacy_job(job)

    def _computed_cgroup(self, job: dict[str, Any]) -> Path:
        slice_name = self.config.borrow.slice if job.get("borrow") else self.config.launch.slice
        return unit_cgroup_path(str(job["unit_name"]), slice_name, root=self._manager_cgroup)

    def _needs_query(self, job: dict[str, Any]) -> bool:
        if job["state"] != "running" or job.get("cancel_requested_at") or job.get("preempt_pending"):
            return True
        cgroup = self._verified.get(int(job["id"]))
        if cgroup is None or not os.path.isdir(cgroup):
            return True
        if job.get("borrow"):
            usage = read_usage(cgroup)
            if usage is None or usage.frozen is None or usage.frozen != bool(job.get("frozen_at")):
                return True
        return False

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
        self._last_unit[job_id] = unit
        self._finalize_unit_job(job, unit, status, stop_post)
        self._last_unit.pop(job_id, None)
        if unit.found and not unit.active and self.store.get(job_id) and \
                self.store.require(job_id)["state"] not in ACTIVE_STATES:
            # 0.7.1: failed job units are no longer collected by systemd; the
            # end is recorded, unload the unit
            reset_failed(self.config.launch, unit.name)

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
        usage = self._final_usage(job)
        if self.store.finish(
            job_id,
            state=state,
            exit_code=exit_code,
            message=message,
            result=result,
            finished_at=finished_at,
            usage=usage,
        ):
            self._append_job_log(job, f"state:     {state}\n===== simq job {job_id} finished =====\n")
            self.log(
                f"job {job_id} {state}"
                + (f" (exit {exit_code})" if exit_code is not None else "")
                + (f": {message}" if message else "")
            )
        self._stop_sent.pop(job_id, None)
        self._tracks.pop(job_id, None)
        self._verified.pop(job_id, None)
        self._forget_secrets(job)
        try:
            self._clean_shm(job)
        except Exception as error:  # noqa: BLE001 - never block a job's end
            self.log(f"job {job_id}: shared-memory cleanup failed: {error!r}")

    def _clean_shm(self, job: dict[str, Any]) -> None:
        """0.7.1: remove the shared-memory segments an MPI job left in /dev/shm
        (killed ranks never unlink them; their shmem stays charged and cannot
        be reclaimed).  Only files of this user, matching MPI's names, made
        after the job started and mapped or held open by no process."""
        started = _parse_time(job.get("started_at"))
        if started is None:
            return
        try:
            entries = [entry for entry in os.scandir(self.shm_dir)
                       if SHM_PATTERNS.match(entry.name) and entry.is_file(follow_symlinks=False)]
        except OSError:
            return
        candidates = []
        uid = os.getuid()
        for entry in entries:
            try:
                info = entry.stat(follow_symlinks=False)
            except OSError:
                continue
            if info.st_uid == uid and info.st_mtime >= started - 5:
                candidates.append((entry.path, info.st_size))
        if not candidates:
            return
        in_use = _shm_in_use(str(self.shm_dir))
        removed, freed = 0, 0
        for path, size in candidates:
            if path in in_use:
                continue
            try:
                os.unlink(path)
            except OSError:
                continue
            removed += 1
            freed += size
        if removed:
            self.log(f"job {job['id']}: removed {removed} leftover shared-memory segment(s) from {self.shm_dir} "
                     f"({format_size(freed) if freed else '0B'})")

    def _forget_secrets(self, job: dict[str, Any]) -> None:
        """0.7.0: a finished job's secrets file goes, its spec.json is redacted."""
        from .envpolicy import delete_secrets, scrub_spec

        try:
            delete_secrets(state_dir(), int(job["id"]), str(job["launch_token"]))
            if job.get("job_dir"):
                scrub_spec(Path(str(job["job_dir"])) / "spec.json")
        except Exception as error:  # noqa: BLE001 - never block a job's end
            self.log(f"job {job['id']}: could not remove its secrets: {error!r}")

    def _final_usage(self, job: dict[str, Any]) -> dict[str, Any] | None:
        """What the job used: the unit's own record at its end (0.7.0 jobwrap:
        stop-post, else the finished status), else the worker's samples."""
        job_id = int(job["id"])
        recorded: dict[str, Any] = {}
        if job.get("job_dir"):
            status, stop_post = self._job_files(job)
            for record in (status, stop_post):  # stop-post is later: it wins
                if record and isinstance(record.get("usage"), dict):
                    recorded = record["usage"]
        track = self._tracks.get(job_id)
        unit = self._last_unit.get(job_id)
        if unit is not None and not recorded.get("memory_peak") and unit.memory_peak:
            # 0.7.1: a killed job (no wrapper record): systemd's own accounting
            recorded = dict(recorded, memory_peak=unit.memory_peak)
        if unit is not None and recorded.get("cpu_usage_usec") is None and unit.cpu_usage_nsec is not None:
            recorded = dict(recorded, cpu_usage_usec=unit.cpu_usage_nsec // 1000)
        peak = recorded.get("memory_peak")
        source = "unit" if peak is not None else None
        if track is not None and track.peak and (peak is None or track.peak > peak):
            peak, source = track.peak, source or "sampled"
        if peak is None and job.get("mem_peak_bytes") is not None:
            peak, source = int(job["mem_peak_bytes"]), job.get("mem_peak_source") or "sampled"
        cpu = recorded.get("cpu_usage_usec")
        if cpu is None and track is not None:
            cpu = track.cpu_usec
        if peak is None and cpu is None:
            return None
        usage: dict[str, Any] = {"peak": peak, "source": source, "cpu_usec": cpu, "sampled_at": utc_now()}
        if track is not None and track.series.buckets:
            usage["series"] = track.series.to_json()
        return usage

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
        if service_result is None and status is not None and status.get("phase") == "finished":
            # 0.7.1: no stop-post hook; a unit that ended normally is unloaded at
            # once, so systemd's result follows from the wrapper's exit status
            code = status.get("returncode")
            service_result = "success" if code == 0 else ("signal" if isinstance(code, int) and code < 0
                                                          else "exit-code")
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
        self.gpu.refresh(jobs)
        if any(job.get("borrow") for job in jobs):
            return self._dispatch_with_borrowers(jobs)
        if self._telemetry_used:
            self.telemetry.reset()
            self._telemetry_used = False
        # Without borrow jobs: the 0.6.0 dispatch plus the GPU line (0.6.2).
        paused = self.store.is_paused()
        queued = [j for j in jobs if j["state"] == "queued"]
        now = self.clock()
        active = [j for j in jobs if j["state"] in ACTIVE_STATES]
        running = running_entries(active, self.config, now, self._currents())
        self._external_busy = False
        self._memory_held = {}
        if paused or not queued:
            self.gpu.idle_pass(jobs=jobs, running=running, now=now, paused=paused)
            return 0
        allocated_cpus = {cpu for item in running for cpu in item.allocation.cpus}
        allocated_gpus = {gpu for item in running for gpu in item.allocation.gpus}
        busy_cpus, busy_gpus = busy_resources(
            self.config, allocated_cpus=allocated_cpus, allocated_gpus=allocated_gpus
        )
        self._external_busy = bool(busy_cpus or busy_gpus)
        view = self.gpu.prepare(jobs=jobs, running=running, busy_gpus=busy_gpus, now=now)
        actual_free, parked = self._memory_budget(active, queued)
        pool = build_pool(self.config, running, external_cpus=busy_cpus, external_gpus=set(view.blocked),
                          actual_free=actual_free)
        requests = [job_request(j, self.config) for j in queued if int(j["id"]) not in parked]
        decision = plan(self.config, pool, running, requests, now)
        starts, deferred = self.gpu.filter_starts(decision.starts, view)
        starts, held = self._admit(starts, active)
        deferred.update(held)
        deferred.update(parked)
        deferred.update(self._memory_reasons(decision.skipped, requests, pool, actual_free))
        self.gpu.after_plan(view, jobs=jobs, running=running, busy_cpus=busy_cpus, requests=requests,
                            starts=starts, now=now, paused=False)
        by_id = {int(j["id"]): j for j in queued}
        started = 0
        for request, allocation in starts:
            if self.stop_requested:
                break
            if self._launch(by_id[request.job_id], request, allocation, busy_cpus=set(busy_cpus), view=view):
                started += 1
        skipped = dict(decision.skipped)
        skipped.update(deferred)
        self._record_reasons(skipped, decision.reservation, busy_cpus)
        self.last_started = started
        return started

    def _record_reasons(
        self,
        skipped: dict[int, str],
        reservation: tuple[int, float] | None,
        busy_cpus: dict[int, list],
    ) -> None:
        for job_id, reason in skipped.items():
            if job_id in self.gpu.reasons:
                reason = self.gpu.reasons[job_id]
            elif reason == "resources busy":
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
            previous = self._reasons.get(job_id)
            if previous != reason:
                # a memory reason carries live numbers: rewrite it at most every 30 s
                if (previous is not None and previous.split(":", 1)[0] == reason.split(":", 1)[0]
                        and reason.startswith(("memory:", "memory pressure", "memory busy", "parked:"))
                        and time.monotonic() - self._reason_written.get(job_id, 0.0) < 30.0):
                    continue
                self._reasons[job_id] = reason
                self._reason_written[job_id] = time.monotonic()
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
        running = running_entries(owners, config, now, self._currents())
        owner_of = {cpu: int(item.job_id) for item in running if item.job_id is not None
                    for cpu in item.allocation.cpus}
        allocated_gpus = {gpu for item in running for gpu in item.allocation.gpus}
        busy_cpus, busy_gpus = busy_resources(config, allocated_cpus=set(owner_of), allocated_gpus=allocated_gpus)
        self._external_busy = bool(busy_cpus or busy_gpus)
        self._memory_held = {}
        pressure = self._refresh_pressure()
        view = self.gpu.prepare(jobs=jobs, running=running, busy_gpus=busy_gpus, now=now)
        actual_free, parked = self._memory_budget(active, queued)
        pool = build_pool(config, running, external_cpus=busy_cpus, external_gpus=set(view.blocked),
                          actual_free=actual_free)
        borrowers = borrower_entries(borrow_active, config)
        self._telemetry_used = True
        cores = self.telemetry.sample(owners=owners, external=busy_cpus, borrowers=borrowers)
        exclusive = {cpu for job in owners if job.get("exclusive_timing") for cpu in _job_cpus(job)}
        requests = [job_request(j, config) for j in queued if int(j["id"]) not in parked]
        decision = plan_node(
            config, pool, running, requests, now,
            borrowers=borrowers, cores=cores, owner_of=owner_of,
            exclusive_owner_cpus=exclusive, starts_allowed=not paused and not pressure,
        )
        starts, deferred = self.gpu.filter_starts(decision.starts, view)
        starts, held = self._admit(starts, active)
        deferred.update(held)
        deferred.update(parked)
        deferred.update(self._memory_reasons(decision.skipped, requests, pool, actual_free))
        if pressure:
            for job_id in [r.job_id for r in requests if r.job_id not in deferred]:
                deferred[job_id] = self._pressure_reason()
        if paused or not queued:
            self.gpu.idle_pass(jobs=jobs, running=running, now=now, paused=paused)
        else:
            self.gpu.after_plan(view, jobs=jobs, running=running, busy_cpus=busy_cpus, requests=requests,
                                starts=starts, now=now, paused=False)
        active_by_id = {int(j["id"]): j for j in borrow_active}
        preemptions = list(decision.preemptions)
        if pressure:
            preemptions = self._pressure_preemptions(borrowers, preemptions)
        for item in preemptions:
            if item.job_id in active_by_id:
                self._preempt(active_by_id[item.job_id], item)
        for job_id, placement in decision.resumes:
            if job_id in active_by_id:
                self._resume(active_by_id[job_id], placement, owner_of, cores)
        by_id = {int(j["id"]): j for j in queued}
        started = 0
        skipped = dict(decision.skipped)
        skipped.update(deferred)
        for request, allocation in starts:
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
                                        owner_of=owner_of, cores=cores, view=view)
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
        home = self.config.canonical_partition(job.get("partition"))
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
        view: GpuView | None = None,
    ) -> bool:
        config = self.config
        handoff = None if view is None else view.handoffs.get(request.job_id)
        # request.partition is the job's recorded name (a partition or an
        # alias, 0.6.3); the job keeps it, the view resolves the partition.
        partition = config.partition(request.partition)
        recorded = request.partition
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
            partition=recorded,
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
                partition=recorded,
                allocation=Allocation(allocation.cpus, allocation.gpus, allocation.gpu_claim,
                                      allocation.gpu_exclusive, mem),
                time_limit=request.time_limit,
                unit=unit,
            )
            if borrow:
                placement["SIMQ_BORROW"] = "1"
                placement["SIMQ_ON_PREEMPT"] = request.on_preempt
            if claimed.get("github_json"):
                # 0.7.0: the summary markdown the job may write for its PR comment
                summary = (json.loads(str(claimed["github_json"])) or {}).get("summary")
                if summary:
                    placement["SIMQ_SUMMARY"] = str(summary)
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
                "partition": recorded if recorded == partition.name else f"{recorded} (pool {partition.name})",
                "gpu_text": gpu_text,
                "mem_text": f"{format_size(mem)} (swap {format_size(swap)})",
                "time_text": format_time_limit(request.time_limit),
                "locks": [{"path": lock.path, "mode": lock.mode, "resource": lock.resource} for lock in locks],
                "submitted_at": str(claimed["created_at"]),
                "scheduler_version": __version__,
                "scheduler_runtime": self.runtime_dir,
                "require_unit": True,
                "wake_path": str(self.wake_path),
            }
            if handoff is not None:
                handoff_paths = set(config.locks.gpu_lock_order(config.gpus[handoff["gpu"]]))
                spec["gpu_handoff"] = {k: v for k, v in handoff.items()}
                for lease in spec["locks"]:
                    if lease["path"] in handoff_paths:
                        lease["handoff"] = True
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
            f"started {'borrow ' if borrow else ''}job {job_id} ({claimed['name']}) on partition "
            f"{partition.display_name()}: "
            f"CPUs {placement['SIMQ_CPU_LIST']}, GPUs {placement['SIMQ_GPUS'] or 'none'}, "
            f"mem {format_size(mem)}, time {format_time_limit(request.time_limit)}, unit {unit}"
            + (f"; {spec.get('borrow_text')}" if borrow else "")
            + (f"; GPU handed over by {handoff['claim_unit']}" if handoff is not None else "")
        )
        if view is not None:
            self.gpu.launched(claimed, view, self.clock())
        return True

    # ----- memory (0.7.0) ------------------------------------------------------

    def _job_currents(self) -> list[int]:
        return [track.current for track in self._tracks.values() if not track.gone]

    def _read_node(self) -> NodeMemory | None:
        node = read_node_memory(self.memory.test_node_total, self._job_currents())
        if node is not None:
            self._node = node
        return node

    def _refresh_pressure(self) -> bool:
        """Read MemAvailable now and update the memory-pressure state."""
        if not self.memory.enabled:
            return False
        node = self._read_node()
        if node is None:
            return self._pressure.active
        high, low = water_marks(self.memory, node.total)
        transition = self._pressure.update(node.available, high, low)
        if transition == "enter":
            self._enter_pressure(node, high, low)
        elif transition == "leave":
            self._leave_pressure(node, high, low)
        return self._pressure.active

    def _pressure_reason(self) -> str:
        record = self._pressure.record
        return (f"memory pressure since {record.get('since', '?')}: MemAvailable "
                f"{format_gib(self._node.available if self._node else None)} below the high-water mark "
                f"{format_gib(record.get('high_water'))}; dispatch resumes at {format_gib(record.get('low_water'))}")

    def _job_rows(self, unix: float) -> list[dict[str, Any]]:
        rows = []
        for track in self._tracks.values():
            if track.gone:
                continue
            rows.append({
                "job_id": track.job_id, "name": track.name, "borrow": track.borrow,
                "current": track.current, "request": track.request, "peak": track.peak,
                "percent_of_request": round(100.0 * track.current / track.request, 1) if track.request else None,
                "grew": track.growth(unix),
            })
        return rows

    def _enter_pressure(self, node: NodeMemory, high: int, low: int) -> None:
        unix = time.time()
        self._pressure.since = unix
        rows = self._job_rows(unix)
        jobs_total = sum(row["current"] for row in rows)
        grew = sorted((row for row in rows if row["grew"] > 0), key=lambda row: -row["grew"])
        culprit = grew[0] if grew else (max(rows, key=lambda row: row["current"]) if rows else None)
        record = {
            "since": _iso(unix),
            "since_unix": unix,
            "available": node.available,
            "total": node.total,
            "simulated": node.simulated,
            "high_water": high,
            "low_water": low,
            "queue_jobs_current": jobs_total,
            "outside_queue": max(0, node.used - jobs_total),
            "culprit": culprit,
            "jobs": rows,
        }
        self._pressure.record = record
        self.store.set_meta(PRESSURE_META, json.dumps(record, sort_keys=True, default=str))
        culprit_text = "no queue job grew"
        if culprit is not None:
            culprit_text = (f"job {culprit['job_id']} ({culprit['name']}) uses {format_gib(culprit['current'])} of its "
                            f"{format_gib(culprit['request'])} request (grew {format_gib(culprit['grew'])} in "
                            f"{GROWTH_WINDOW_SECONDS:.0f} s)")
        self.log(f"memory pressure: MemAvailable {format_gib(node.available)} below {format_gib(high)}; "
                 f"stopping dispatch, preempting borrow jobs; {culprit_text}; queue jobs use {format_gib(jobs_total)}, "
                 f"outside the queue {format_gib(record['outside_queue'])}")
        for row in rows:
            if row["grew"] > 0 or row is culprit:
                self.store.add_event(row["job_id"], "memory-pressure", {
                    "available": format_gib(node.available), "high_water": format_gib(high),
                    "current": format_gib(row["current"]), "request": format_gib(row["request"]),
                    "grew": format_gib(row["grew"]), "window_s": GROWTH_WINDOW_SECONDS,
                    "culprit": row is culprit, "action": "preempted per --on-preempt" if row["borrow"]
                    else "recorded (owner jobs are not stopped)",
                })
                if not row["borrow"]:
                    self.log(f"job {row['job_id']}: {format_gib(row['current'])} of {format_gib(row['request'])} "
                             f"requested (grew {format_gib(row['grew'])}); owner job left running")
        self._record_samples(force=True)
        self._wake_full = True

    def _leave_pressure(self, node: NodeMemory, high: int, low: int) -> None:
        since = self._pressure.record.get("since_unix")
        lasted = "" if since is None else f" after {time.time() - float(since):.0f} s"
        self.log(f"memory pressure over{lasted}: MemAvailable {format_gib(node.available)} at or above "
                 f"{format_gib(low)}; dispatch resumes")
        for job_id in list(self._pressure_preempted):
            self.store.add_event(job_id, "memory-recovered", {"available": format_gib(node.available)})
        self._pressure_preempted.clear()
        self._pressure.record = {}
        self.store.set_meta(PRESSURE_META, None)
        self._wake_full = True

    def _pressure_preemptions(self, borrowers: list[Borrower], planned: list[Preemption]) -> list[Preemption]:
        """Memory pressure: preempt every running borrow job by its policy; a
        frozen one is stopped when the pressure lasts longer than
        ``pressure_escalate_seconds`` (a frozen job keeps its memory)."""
        result = {item.job_id: item for item in planned}
        reason = (f"memory pressure (MemAvailable {format_gib(self._node.available if self._node else None)} "
                  f"below {format_gib(self._pressure.record.get('high_water'))})")
        now = time.time()
        for borrower in borrowers:
            if borrower.pending:
                continue
            mode = borrower.on_preempt
            if borrower.frozen_since is not None:
                since = self._pressure_preempted.get(borrower.job_id, borrower.frozen_since)
                if now - since < self.memory.pressure_escalate_seconds:
                    continue
                result[borrower.job_id] = Preemption(borrower.job_id, "requeue",
                                                     reason + f"; frozen {now - since:.0f} s keeps its memory",
                                                     None, True, borrower.cpus)
                continue
            if borrower.job_id in result and result[borrower.job_id].mode != "freeze":
                continue
            result[borrower.job_id] = Preemption(borrower.job_id, mode, reason, None, False, borrower.cpus)
            self._pressure_preempted.setdefault(borrower.job_id, now)
        return [result[key] for key in sorted(result)]

    def _memory_reasons(self, skipped: dict[int, str], requests: list[Request], pool: Any,
                        actual_free: int | None) -> dict[int, str]:
        """Jobs the planner held because actual memory is short (0.7.1)."""
        if actual_free is None:
            return {}
        found: dict[int, str] = {}
        for request in requests:
            reason = skipped.get(request.job_id, "")
            if reason.startswith("resources busy") and request.mem > actual_free:
                found[request.job_id] = (f"memory busy: needs {format_size(request.mem)}, "
                                         f"{format_gib(max(0, actual_free))} actually free now")
        if found:
            held = dict(self._memory_held)
            held.update(found)
            self._memory_held = held
        return found

    def _currents(self) -> dict[int, int]:
        return {job_id: track.current for job_id, track in self._tracks.items() if not track.gone}

    def _memory_budget(self, active: list[dict[str, Any]],
                       queued: list[dict[str, Any]]) -> tuple[int | None, dict[int, str]]:
        """0.7.1: the actual memory starts may take now (MemAvailable - margin -
        what young jobs have not allocated), for the planner, and the queued
        jobs to park: a request above what this node can ever admit -- the
        best MemAvailable seen with the queue's jobs taken away, minus the
        margin, or memory_schedulable minus the margin -- would hold the head
        reservation forever; parked, it waits without blocking the queue and
        is looked at again on every pass."""
        if not self.memory.enabled:
            return None, {}
        node = self._read_node()
        if node is None:
            return None, {}
        margin = admission_margin(self.memory, node.total)
        idle = node.available + sum(self._currents().values())
        now = time.time()
        self._idle_seen.append((now, idle))
        while self._idle_seen and self._idle_seen[0][0] < now - 3600:
            self._idle_seen.popleft()
        best = max(value for _, value in self._idle_seen)
        ceiling = best - margin
        if self.config.memory_schedulable is not None:
            ceiling = min(ceiling, self.config.memory_schedulable - margin)
        self._admissible = ceiling
        parked: dict[int, str] = {}
        for job in queued:
            request = int(job.get("mem_bytes") or job_request(job, self.config).mem or 0)
            if request > ceiling:
                parked[int(job["id"])] = (
                    f"parked: mem {format_size(request)} exceeds admissible ~{format_gib(ceiling)} on "
                    f"{self.config.node} (MemAvailable with the queue idle ~{format_gib(best)} - margin "
                    f"{format_gib(margin)}); it starts when that much is free, other jobs go ahead")
        return node.available - margin - self._ramping(active), parked

    def _ramping(self, active: list[dict[str, Any]]) -> int:
        """Requests of jobs started less than ramp_seconds ago that they have
        not allocated yet (counted as used by the admission guard)."""
        now = time.time()
        total = 0
        for job in active:
            started = _parse_time(job.get("started_at"))
            if job["state"] != "starting" and (started is None or now - started > self.memory.ramp_seconds):
                continue
            request = int(job.get("mem_bytes") or 0)
            track = self._tracks.get(int(job["id"]))
            current = track.current if track is not None else 0
            if track is None and job["state"] != "starting":
                cgroup = self._verified.get(int(job["id"])) or (
                    self._computed_cgroup(job) if job.get("unit_name") else None)
                usage = read_usage(cgroup, full=False) if cgroup is not None else None
                current = usage.current if usage is not None else 0
            total += max(0, request - current)
        return total

    def _admit(self, starts: list[tuple[Request, Allocation]],
               active: list[dict[str, Any]]) -> tuple[list[tuple[Request, Allocation]], dict[int, str]]:
        """The actual-memory admission guard (0.7.0).  Reservations stay the
        primary constraint (the planner); this only holds starts back."""
        if not starts or not self.memory.enabled:
            return starts, {}
        if self._refresh_pressure():
            held = {request.job_id: self._pressure_reason() for request, _ in starts}
            self._memory_held = held
            return [], held
        node = self._node
        if node is None:
            return starts, {}
        needs = []
        for request, allocation in starts:
            mem = allocation.mem or self.config.partition(request.partition).default_mem
            needs.append((request.job_id, int(mem)))
        result = admission(needs, available=node.available, margin=admission_margin(self.memory, node.total),
                           ramping=self._ramping(active))
        if result.held:
            for job_id, reason in result.held.items():
                self._note_once(job_id, "memory-held", f"job {job_id}: held by the memory guard: {reason}")
        self._memory_held = dict(result.held)
        admitted = set(result.admitted)
        return [(r, a) for r, a in starts if r.job_id in admitted], dict(result.held)

    def _sync_tracks(self) -> None:
        """Follow every running job unit whose cgroup is known."""
        if not self.memory.enabled:
            return
        jobs = {int(j["id"]): j for j in self.store.active_jobs()
                if j.get("unit_name") and j["state"] in ("running", "cancelling")}
        for job_id in list(self._tracks):
            if job_id not in jobs:
                track = self._tracks.pop(job_id)
                # requeued or recorded elsewhere: keep what it used
                if track.peak:
                    self.store.record_final_usage(job_id, {
                        "peak": track.peak, "source": "sampled", "cpu_usec": track.cpu_usec,
                        "sampled_at": utc_now(), "series": track.series.to_json()})
        for job_id, job in jobs.items():
            cgroup = self._verified.get(job_id)
            if cgroup is None:
                continue
            track = self._tracks.get(job_id)
            if track is not None and track.cgroup == cgroup:
                continue
            track = JobTrack(
                job_id=job_id,
                name=str(job.get("name") or ""),
                unit=str(job["unit_name"]),
                cgroup=cgroup,
                request=int(job.get("mem_bytes") or 0),
                borrow=bool(job.get("borrow")),
                series=Series.from_json(job.get("mem_series_json"), self.memory.series_buckets),
                peak=int(job.get("mem_peak_bytes") or 0),
                cpu_usec=job.get("cpu_usage_usec"),
            )
            usage = read_usage(cgroup, full=False)
            if usage is not None:
                track.observe(time.time(), usage.current, usage.peak)
            self._tracks[job_id] = track

    def _sample(self) -> None:
        """Read every running job's cgroup and MemAvailable (no process is
        spawned); note ended units, near-limit jobs and memory pressure."""
        self._next_sample = time.monotonic() + self.memory.sample_seconds
        unix = time.time()
        for track in self._tracks.values():
            usage = read_usage(track.cgroup, full=False)
            if usage is None:
                if not track.gone:
                    track.gone = True
                    self._wake_full = True  # the unit ended: reconcile it
                continue
            track.observe(unix, usage.current, usage.peak)
            if (not track.near_limit_noted and track.request
                    and usage.current >= self.memory.near_limit_fraction * track.request):
                track.near_limit_noted = True
                self.log(f"job {track.job_id} ({track.name}) uses {format_gib(usage.current)} of its "
                         f"{format_gib(track.request)} memory limit ({100.0 * usage.current / track.request:.0f}%)")
                self.store.add_event(track.job_id, "memory-near-limit", {
                    "current": format_gib(usage.current), "limit": format_gib(track.request),
                    "percent": round(100.0 * usage.current / track.request, 1)})
        self._refresh_pressure()
        if time.monotonic() >= self._next_record:
            self._record_samples()

    def _record_samples(self, *, force: bool = False) -> None:
        self._next_record = time.monotonic() + self.memory.record_seconds
        updates: dict[int, dict[str, Any]] = {}
        now = utc_now()
        for track in self._tracks.values():
            if track.gone or not (track.dirty or force):
                continue
            usage = read_usage(track.cgroup)
            if usage is not None and usage.cpu_usec is not None:
                track.cpu_usec = usage.cpu_usec
            updates[track.job_id] = {
                "current": track.current, "peak": track.peak, "source": "sampled", "sampled_at": now,
                "series": track.series.to_json(), "cpu_usec": track.cpu_usec,
            }
            track.dirty = False
        if updates:
            self.store.record_memory(updates)

    def _check_exits(self) -> None:
        """A verified job unit whose cgroup is gone has ended: reconcile it."""
        for job_id, cgroup in self._verified.items():
            if not os.path.isdir(cgroup):
                track = self._tracks.get(job_id)
                if track is not None:
                    track.gone = True
                self._wake_full = True
                return

    def _tick(self) -> None:
        """Between passes, every sample_seconds while job units run: unit ends
        (cgroup existence) and, with memory tracking, the samples."""
        self._next_sample = time.monotonic() + self.memory.sample_seconds
        self._check_exits()
        if self.memory.enabled and self._tracks:
            self._sample()

    def _check_clients(self) -> None:
        """A client wrote to the queue database (data_version changed)."""
        self._next_check = time.monotonic() + IDLE_CHECK_SECONDS
        try:
            version = self.store.data_version()
        except Exception:  # noqa: BLE001 - the next pass reports database errors
            self._wake_full = True
            return
        if self._data_version is not None and version != self._data_version:
            self._wake_full = True

    # ----- GitHub reporting (0.7.0) ---------------------------------------------

    def _github_notify(self) -> None:
        """Called by the posting thread: wake the loop to record results."""
        if self._wake_fd is not None:
            try:
                os.write(self._wake_fd, b"g")
            except OSError:
                pass

    def _github_check(self) -> None:
        """Hand every job whose phase changed since its last post to the
        posting thread; never waits for GitHub."""
        if self.github is None:
            return
        self._github_results()
        now = time.time()
        self._github_next = None
        for job in self.store.github_pending():
            job_id = int(job["id"])
            phase = job_phase(job)
            if phase == job.get("github_phase") or job_id in self.github.in_flight:
                continue
            retry_at = job.get("github_next_unix")
            if retry_at is not None and float(retry_at) > now and int(job.get("github_attempts") or 0) > 0:
                self._github_next = min(self._github_next or float(retry_at), float(retry_at))
                continue
            try:
                task = build_task(job, phase, self.config.node)
            except ValueError as error:
                self.store.record_github(job_id, phase=phase, ok=False, error=f"bad GitHub targets: {error}",
                                         next_unix=None, give_up=True, final=phase == FINAL, comment_id=None)
                continue
            self.github.submit(task)

    def _github_results(self) -> None:
        if self.github is None:
            return
        for result in self.github.results():
            job_id = int(result["job_id"])
            phase = str(result["phase"])
            final = phase == FINAL
            if result.get("ok"):
                self.github_posts += 1
                self.store.record_github(job_id, phase=phase, ok=True, error=None, next_unix=None, give_up=False,
                                         final=final, comment_id=result.get("comment_id"))
                continue
            job = self.store.get(job_id) or {}
            attempts = int(job.get("github_attempts") or 0) + 1
            error = str(result.get("error") or "unknown error")
            give_up = not result.get("retry") or attempts >= self.config.github.max_attempts
            delay = min(self.config.github.backoff_max_seconds,
                        self.config.github.backoff_seconds * (2 ** (attempts - 1)))
            self.store.record_github(job_id, phase=phase, ok=False, error=error,
                                     next_unix=None if give_up else time.time() + delay, give_up=give_up,
                                     final=final, comment_id=None)
            self.store.add_event(job_id, "github-failed", {
                "phase": phase.split(":", 1)[0], "attempt": attempts, "error": error,
                "next": "given up" if give_up else f"retry in {delay:.0f} s"})
            self._note_once(job_id, f"github-{phase}-{give_up}",
                            f"job {job_id}: GitHub post ({phase.split(':', 1)[0]}) failed: {error}; "
                            + ("given up" if give_up else f"retrying with backoff (attempt {attempts})"))
            if not give_up:
                retry_at = time.time() + delay
                self._github_next = min(self._github_next or retry_at, retry_at)

    # ----- event loop (0.7.0) -----------------------------------------------------

    def _open_wake(self) -> None:
        path = self.wake_path
        try:
            info = os.lstat(path)
            if not stat.S_ISFIFO(info.st_mode):
                path.unlink()
                raise FileNotFoundError
        except FileNotFoundError:
            try:
                os.mkfifo(path, 0o600)
            except FileExistsError:
                pass
        except OSError as error:
            self.log(f"wake FIFO {path} unusable ({error}); polling instead")
            return
        try:
            # read and write: the FIFO never reports end-of-file, writers never block
            self._wake_fd = os.open(path, os.O_RDWR | os.O_NONBLOCK)
        except OSError as error:
            self.log(f"cannot open the wake FIFO {path}: {error}; polling instead")
            self._wake_fd = None

    def _close_wake(self) -> None:
        if self._wake_fd is not None:
            os.close(self._wake_fd)
            self._wake_fd = None

    def _drain_wake(self) -> bytes:
        received = b""
        while self._wake_fd is not None:
            try:
                data = os.read(self._wake_fd, 4096)
            except (BlockingIOError, InterruptedError):
                break
            except OSError:
                break
            if not data:
                break
            received += data
        return received

    def _wait(self, timeout: float) -> None:
        if timeout <= 0:
            return
        if self._wake_fd is None:
            time.sleep(min(timeout, 0.5))
            return
        try:
            ready, _, _ = select.select([self._wake_fd], [], [], timeout)
        except InterruptedError:
            ready = []
        if not ready:
            return
        received = self._drain_wake()
        if received.strip(b"g"):  # anything but the posting thread's results
            self._wake_full = True
            self._watch_until = time.monotonic() + EXIT_WATCH_SECONDS

    def _full_pass(self, startup: bool) -> int:
        self._wake_full = False
        self.passes += 1
        try:
            self._data_version = self.store.data_version()
        except Exception:  # noqa: BLE001
            self._data_version = None
        now = time.monotonic()
        safety = startup or now >= self._next_safety
        if safety:
            self._next_safety = now + self.config.launch.idle_seconds * SAFETY_RECONCILE_IDLE
        self.reconcile(startup=startup, full=safety)
        if self._pressure.active:
            self._refresh_pressure()  # also when no job runs any more (no samples then)
        started = self.dispatch()
        self._sync_tracks()
        try:
            self._github_check()
        except Exception as error:  # noqa: BLE001 - reporting never stops scheduling
            self.log(f"GitHub reporting pass failed: {error!r}")
        self._fast_reasons = self._watch_reasons()
        interval = self.config.launch.poll_seconds if self._fast_reasons else self.config.launch.idle_seconds
        self._next_full = time.monotonic() + interval
        return started

    def _watch_reasons(self) -> list[str]:
        """Why the next pass must come after poll_seconds rather than idle_seconds."""
        reasons = []
        jobs = self.store.node_jobs()
        queued = [j for j in jobs if j["state"] == "queued"]
        if self.drain:
            reasons.append("drain")
        if any(j["state"] in ("starting", "cancelling") or j.get("cancel_requested_at") or j.get("preempt_pending")
               for j in jobs):
            reasons.append("jobs starting or stopping")
        if any(j.get("borrow") for j in jobs):
            reasons.append("borrow jobs")
        if self.gpu.claims or self.gpu.ready_since:
            reasons.append("GPU line")
        if queued and self._external_busy:
            reasons.append("leases held outside the queue")
        if queued and self._memory_held:
            reasons.append("memory guard")
        if self._pressure.active:
            reasons.append("memory pressure")
        if any(job["state"] in ACTIVE_STATES and job.get("unit_name") and int(job["id"]) not in self._verified
               for job in jobs):
            reasons.append("job units not yet verified")
        return reasons

    def _timeout(self) -> float:
        now = time.monotonic()
        deadlines = [self._next_full]
        if self._verified:
            deadlines.append(self._next_sample)
        deadlines.append(self._next_check)
        if self._watch_until > now:
            deadlines.append(now + EXIT_WATCH_STEP)
        if self._github_next is not None:
            deadlines.append(now + max(0.0, self._github_next - time.time()))
        return max(0.0, min(deadlines) - now)

    # ----- main loop -------------------------------------------------------

    def run(self) -> int:
        if not self._acquire_lock():
            self.log("another worker already owns the queue lock")
            return 0
        signal.signal(signal.SIGTERM, self._handle_stop)
        signal.signal(signal.SIGINT, self._handle_stop)
        self.job_root.mkdir(mode=0o700, parents=True, exist_ok=True)
        if not self.once:
            self._open_wake()
        mem = self.memory
        self.log(
            f"node worker started; simqueue {__version__} from {self.runtime_dir}; "
            f"config {self.config.path}; partitions "
            + ", ".join(
                f"{p.name}={format_cpu_list(p.cpus)}" + (f"+{len(p.gpus)}gpu" if p.gpus else "")
                + (f" (aliases {','.join(p.aliases)})" if p.aliases else "")
                for p in self.config.partitions.values()
            )
            + f"; drain={self.drain}"
            + f"; event-driven (wake FIFO {'on' if self._wake_fd is not None else 'off'}, poll "
              f"{self.config.launch.poll_seconds:g} s while watching, else {self.config.launch.idle_seconds:g} s)"
            + (f"; memory sampling every {mem.sample_seconds:g} s"
               + (f" (simulated node of {format_gib(mem.test_node_total)})" if mem.test_node_total else "")
               if mem.enabled else "; memory tracking off")
        )
        try:
            self.gpu.adopt()
        except Exception as error:  # noqa: BLE001 - the line is re-established on later passes
            self.log(f"GPU claim adoption failed: {error!r}")
        if mem.enabled:
            # a previous worker's pressure record ends with it; it is re-detected if it persists
            try:
                if self.store.get_meta(PRESSURE_META) is not None:
                    self.store.set_meta(PRESSURE_META, None)
            except Exception:  # noqa: BLE001
                pass
        exit_code = 0
        failures = 0
        startup = True
        started = 0
        while not self.stop_requested:
            passed = False
            try:
                if self._wake_full or time.monotonic() >= self._next_full:
                    started = self._full_pass(startup)
                    startup = False
                    passed = True
                if self._verified and time.monotonic() >= self._next_sample:
                    self._tick()
                elif self._watch_until > time.monotonic():
                    self._check_exits()
                if time.monotonic() >= self._next_check or self._verified:
                    self._check_clients()
                if self.github is not None and (self.github.in_flight or (
                        self._github_next is not None and time.time() >= self._github_next)):
                    self._github_check()
                failures = 0
            except Exception as error:  # one bad pass must not end the scheduler
                failures += 1
                self.log(f"node worker pass failed ({failures} in a row): {error!r}")
                if failures >= 30:
                    exit_code = 2
                    break
                started = 0
                passed = True
                self._wake_full = False
                self._next_full = time.monotonic() + self.config.launch.poll_seconds
            if self.once:
                break
            if self.drain and passed and failures == 0 and started == 0 and not self.store.active_jobs():
                break
            if self._wake_full and failures == 0:
                continue
            self._wait(self._timeout())
        remaining = len(self.store.active_jobs())
        if not self.once:
            self.gpu.withdraw()
            if self._pressure.active:
                self.store.set_meta(PRESSURE_META, None)
        if self._tracks:
            try:
                self._record_samples(force=True)
            except Exception as error:  # noqa: BLE001
                self.log(f"could not record the last memory samples: {error!r}")
        suffix = f" after signal {self.stop_signal}" if self.stop_signal else ""
        self.log(f"node worker stopped{suffix}; {remaining} job unit(s) continue independently; "
                 f"{self.passes} pass(es), {self.unit_queries} systemd unit quer(ies)")
        if self.github is not None:
            self._github_flush(timeout=min(5.0, self.config.github.timeout_seconds))
            self.github.stop()
        self._close_wake()
        self._release_lock()
        self.store.close()
        return exit_code

    def _github_flush(self, timeout: float) -> None:
        """At exit (and for --once/--drain): post what is due, waiting at most
        *timeout* seconds; the rest is posted by the next worker."""
        if self.github is None:
            return
        try:
            self._github_check()
            deadline = time.monotonic() + timeout
            while self.github.in_flight and time.monotonic() < deadline:
                time.sleep(0.05)
                self._github_results()
            self._github_results()
        except Exception as error:  # noqa: BLE001
            self.log(f"GitHub reporting at exit failed: {error!r}")
