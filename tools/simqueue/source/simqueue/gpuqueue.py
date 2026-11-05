"""Worker side of the GPU line (SimQueue 0.6.2).

Each pass the node worker

1. decides, per GPU, whether queued jobs may take it directly now (it is
   free, nobody earlier is in the GPU line) or not (held outside the queue,
   the queue's claim is still waiting, or an earlier ``gpu-lease`` user);
2. plans as before (EASY backfill), the GPU counting as busy when queued
   jobs may not take it;
3. finds the jobs that are *ready except for the GPU*: an exclusive job
   that the same planner would start if the GPU were released (by an
   outside holder or by the queue jobs on it), a shared (``--gpu-mem``) job
   that it would start if outside holders released it while the queue's
   own jobs on it keep running (see ``ready_except_gpu``).  When there are
   any, the queue joins the GPU line
   with one *claim* per GPU: a small systemd unit (``gpuline.py claim``)
   that waits for the GPU lease files in the kernel queue (raw ``flock``
   waiters are moved behind it), holds them, and hands them to the job unit
   the worker starts next;
4. publishes the line (``<pic-profile>/simq-gpu-line.json``): the queued GPU
   jobs in dispatch order with their since-times, the claim, the holders and
   the raw waiters, so that ``gpu-lease`` and other tools can yield.

A claim's line position (its ``since``) is the submit time of the oldest
ready job; ``gpu-lease`` users that started earlier go first, later ones
after the queue.  Claims are units of their own: stopping, restarting or
upgrading the worker leaves them waiting or holding, and the next worker
adopts them by unit name from ``<state>/claims/<token>/``.
"""

from __future__ import annotations

import ctypes
import json
import os
import secrets
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable

from . import __version__, gpuline
from .config import state_dir
from .nodeconfig import NodeConfig, format_size
from .scheduler import Allocation, Request, Running, build_pool, plan
from .units import (
    UnitError,
    UnitState,
    build_claim_command,
    claim_unit_name,
    read_json,
    show_units,
    start_unit,
    stop_unit,
)


CLAIM_IDLE_GRACE_S = 5.0  # a held claim no queued job can use is released after this
HEARTBEAT_S = 10.0  # the published state is rewritten at least this often
CLAIM_DIRS_KEPT = 40  # finished claim directories kept for diagnosis
HANDOFF_STUCK_S = 60.0
STOP_RESEND_S = 30.0


def ready_except_gpu(
    config: NodeConfig,
    uuid: str,
    *,
    running: list[Running],
    owners: list[Request],
    busy_cpus: Iterable[int],
    blocked: Iterable[str],
    started: set[int],
    now: float,
) -> list[int]:
    """Queued jobs that wait only for the lease files of GPU *uuid*, in
    dispatch order (*owners* is the planner's queue).

    An exclusive request is ready when the planner would start it once
    everybody holding the GPU released it: outside holders and the queue's
    own jobs on it (whose cores free at the same moment).  A shared
    (``--gpu-mem``) request is ready only when the planner would start it
    once outside holders released the GPU while the queue's jobs on it keep
    running: a shared job that waits for cores or for GPU budget held by
    queue jobs is not waiting for the lease files (it joins their shared
    lock when it starts), and a claim for it would hold back every queued
    GPU job until the GPU drained.
    """
    external = set(blocked) - {uuid}
    cpus = set(busy_cpus)
    without = [r for r in running if uuid not in r.allocation.gpus]
    whole = plan(config, build_pool(config, without, external_cpus=cpus, external_gpus=external), without,
                 owners, now)
    kept = plan(config, build_pool(config, running, external_cpus=cpus, external_gpus=external), running,
                owners, now)
    ready = {request.job_id for request, allocation in whole.starts
             if allocation.gpus == (uuid,) and allocation.gpu_exclusive}
    ready |= {request.job_id for request, allocation in kept.starts
              if allocation.gpus == (uuid,) and not allocation.gpu_exclusive}
    return [request.job_id for request in owners if request.job_id in ready and request.job_id not in started]


def _iso(unix: float | None) -> str | None:
    return None if unix is None else gpuline.now_iso(unix)


class _StatxTimestamp(ctypes.Structure):
    _fields_ = [("tv_sec", ctypes.c_int64), ("tv_nsec", ctypes.c_uint32), ("reserved", ctypes.c_int32)]


class _Statx(ctypes.Structure):
    # struct statx (linux/stat.h); padded to its full 256 bytes
    _fields_ = [
        ("stx_mask", ctypes.c_uint32), ("stx_blksize", ctypes.c_uint32), ("stx_attributes", ctypes.c_uint64),
        ("stx_nlink", ctypes.c_uint32), ("stx_uid", ctypes.c_uint32), ("stx_gid", ctypes.c_uint32),
        ("stx_mode", ctypes.c_uint16), ("spare0", ctypes.c_uint16), ("stx_ino", ctypes.c_uint64),
        ("stx_size", ctypes.c_uint64), ("stx_blocks", ctypes.c_uint64), ("stx_attributes_mask", ctypes.c_uint64),
        ("stx_atime", _StatxTimestamp), ("stx_btime", _StatxTimestamp), ("stx_ctime", _StatxTimestamp),
        ("stx_mtime", _StatxTimestamp), ("rest", ctypes.c_uint8 * 128),
    ]


_STATX_BTIME = 0x800
_AT_FDCWD = -100
_libc: Any = None
_SUBMIT_CACHE: dict[tuple[int, str], float | None] = {}


def birth_time(path: str) -> float | None:
    """A file's creation time (statx STATX_BTIME), or ``None``."""
    global _libc
    try:
        if _libc is None:
            _libc = ctypes.CDLL(None, use_errno=True)
        buffer = _Statx()
        if _libc.statx(_AT_FDCWD, os.fsencode(path), 0, _STATX_BTIME, ctypes.byref(buffer)) != 0:
            return None
        if not buffer.stx_mask & _STATX_BTIME or buffer.stx_btime.tv_sec <= 0:
            return None
        return buffer.stx_btime.tv_sec + buffer.stx_btime.tv_nsec / 1e9
    except (AttributeError, OSError, ValueError):
        return None


def created_unix(job: dict[str, Any]) -> float:
    """A job's place in the GPU line: its submit time.

    The database records ``created_at`` in whole seconds, so the exact time
    comes from the job's log file, which ``simq submit`` creates (empty) in
    the same transaction: its birth time, or, without one, its modification
    time while that still lies within the ``created_at`` second (nothing has
    been written yet).  Without either, the end of the ``created_at`` second
    is used (a gpu-lease user that started within that second goes first)."""
    value = gpuline.parse_iso(job.get("created_at"))
    if value is None:
        return float(job.get("id") or 0)
    base = float(value)
    log_path = str(job.get("log_path") or "")
    key = (int(job.get("id") or 0), log_path)
    if key not in _SUBMIT_CACHE:
        exact: float | None = None
        if log_path:
            born = birth_time(log_path)
            if born is not None and base - 1.0 <= born <= base + 60.0:
                exact = born
            else:
                try:
                    modified = os.stat(log_path).st_mtime
                except OSError:
                    modified = None
                if modified is not None and base <= modified < base + 1.0:
                    exact = modified
        if exact is not None or len(_SUBMIT_CACHE) < 100000:
            _SUBMIT_CACHE[key] = exact
    exact = _SUBMIT_CACHE.get(key)
    return exact if exact is not None else base + 0.999


@dataclass
class ClaimRecord:
    token: str
    uuid: str
    unit: str
    directory: Path
    spec: dict[str, Any]
    created: float
    status: dict[str, Any] | None = None
    unit_state: UnitState | None = None
    handoff_job: int | None = None
    handoff_at: float | None = None
    idle_since: float | None = None
    control_since: float | None = None
    stop_sent: float | None = None

    @property
    def phase(self) -> str:
        return str((self.status or {}).get("phase") or "starting")

    @property
    def active(self) -> bool:
        return self.unit_state is None or self.unit_state.active

    @property
    def held(self) -> bool:
        return self.active and self.phase in ("held", "handing-off") and self.stop_sent is None

    def brief(self) -> dict[str, Any]:
        status = self.status or {}
        return {
            "unit": self.unit,
            "phase": self.phase,
            "pid": status.get("pid"),
            "waiter_pid": status.get("waiter_pid"),
            "since_unix": self.control_since or self.spec.get("since_unix"),
            "held_since_unix": status.get("held_since_unix"),
            "handoff_job": self.handoff_job,
            "note": status.get("note"),
            "requeued": status.get("requeued"),
        }


@dataclass
class GpuView:
    """Per-pass GPU line view."""

    blocked: dict[str, str] = field(default_factory=dict)  # uuid -> why queue jobs cannot take it directly
    reserved: set[str] = field(default_factory=set)  # uuids held by the queue's claim, free for a hand-off
    tickets: dict[str, list[dict[str, Any]]] = field(default_factory=dict)
    snapshots: dict[str, dict[str, Any]] = field(default_factory=dict)
    ready: dict[str, list[int]] = field(default_factory=dict)
    handoffs: dict[int, dict[str, Any]] = field(default_factory=dict)  # job id -> hand-off spec


class GpuQueue:
    def __init__(self, config: NodeConfig, store: Any, log: Callable[[str], None], *, runtime_dir: str) -> None:
        self.config = config
        self.cfg = config.line
        self.store = store
        self.log = log
        self.runtime_dir = runtime_dir
        self.enabled = bool(config.line.enabled and config.locks.enabled and config.gpus)
        self.root = state_dir() / "claims"
        self.claims: dict[str, ClaimRecord] = {}
        self.ready_since: dict[int, float] = {}
        self.reasons: dict[int, str] = {}
        self._noted: set[tuple[int, str]] = set()
        self._published_body: str | None = None
        self._published_at = 0.0
        self._publish_error: str | None = None

    # ----- helpers -------------------------------------------------------------

    def gpu_locks(self, uuid: str) -> gpuline.GpuLocks:
        gpu = self.config.gpus[uuid]
        return gpuline.GpuLocks(uuid, gpu.index, gpu.name, self.config.locks.gpu_lock_order(gpu))

    def _note_once(self, job_id: int, key: str, message: str) -> None:
        if (job_id, key) not in self._noted:
            self._noted.add((job_id, key))
            self.log(message)

    # ----- claims: adoption and refresh ---------------------------------------

    def adopt(self) -> None:
        """Re-adopt the claims of a previous worker (they run in their own units)."""
        if not self.enabled:
            return
        self.root.mkdir(mode=0o700, parents=True, exist_ok=True)
        state = gpuline.load_queue_state(self.cfg.state_path)
        for item in ((state or {}).get("gpus") or {}).values():
            for entry in item.get("queue") or []:
                if entry.get("job_id") is not None and entry.get("ready_since_unix"):
                    self.ready_since[int(entry["job_id"])] = float(entry["ready_since_unix"])
        records: list[ClaimRecord] = []
        for directory in sorted(self.root.iterdir()):
            spec = read_json(directory / "spec.json")
            if spec is None or spec.get("schema") != gpuline.CLAIM_SPEC_SCHEMA:
                continue
            records.append(ClaimRecord(token=str(spec["token"]), uuid=str(spec["gpu"]["uuid"]), unit=str(spec["unit"]),
                                       directory=directory, spec=spec, created=float(spec.get("created_unix") or 0),
                                       control_since=spec.get("since_unix")))
        if not records:
            return
        try:
            states = show_units(self.config.launch, [r.unit for r in records])
        except UnitError as error:
            self.log(f"cannot query GPU claim units: {error}; adopting them on a later pass")
            for record in records:
                self.claims.setdefault(record.uuid, record)
            return
        for record in sorted(records, key=lambda r: r.created):
            record.unit_state = states[record.unit]
            record.status = read_json(record.directory / "status.json")
            control = read_json(record.directory / "control.json") or {}
            if control.get("since_unix") is not None:
                record.control_since = float(control["since_unix"])
            if record.unit_state.active:
                older = self.claims.get(record.uuid)
                if older is not None:
                    self.stop_claim(older, "a newer claim exists for the same GPU")
                self.claims[record.uuid] = record
                self.log(f"adopted GPU claim {record.unit} for {record.uuid} ({record.phase})")
            else:
                self._retire(record, log=False)
        self._prune()

    def refresh(self, jobs: list[dict[str, Any]]) -> None:
        if not self.enabled or not self.claims:
            return
        try:
            states = show_units(self.config.launch, [c.unit for c in self.claims.values()])
        except UnitError as error:
            self.log(f"cannot query GPU claim units: {error}")
            return
        by_id = {int(j["id"]): j for j in jobs}
        for uuid, record in list(self.claims.items()):
            record.unit_state = states[record.unit]
            record.status = read_json(record.directory / "status.json") or record.status
            if record.handoff_job is not None:
                job = by_id.get(record.handoff_job)
                if job is None or job["state"] not in ("starting", "running", "cancelling"):
                    record.handoff_job = None  # refused before the hand-off: the claim keeps the GPU
                elif record.held and record.handoff_at and time.time() - record.handoff_at > HANDOFF_STUCK_S:
                    self.log(f"GPU claim {record.unit}: job {record.handoff_job} did not take the GPU; "
                             "the claim is available again")
                    record.handoff_job = None
            if not record.active:
                del self.claims[uuid]
                self._retire(record)
            elif record.stop_sent is not None and time.monotonic() - record.stop_sent > STOP_RESEND_S:
                stop_unit(self.config.launch, record.unit)
                record.stop_sent = time.monotonic()

    def _retire(self, record: ClaimRecord, *, log: bool = True) -> None:
        status = read_json(record.directory / "status.json") or record.status or {}
        if log:
            detail = status.get("phase") or "ended"
            if status.get("handed_to"):
                handed = status["handed_to"]
                detail += f" to job {handed.get('job_id')}"
            if status.get("reason"):
                detail += f": {status['reason']}"
            if status.get("requeued"):
                detail += f"; requeued {status['requeued']} raw waiter(s)"
            self.log(f"GPU claim {record.unit} ended ({detail})")
        try:
            (record.directory / "finished").write_text(json.dumps(status, default=str) + "\n", encoding="utf-8")
        except OSError:
            pass

    def _prune(self) -> None:
        try:
            finished = sorted((d for d in self.root.iterdir() if (d / "finished").exists()),
                              key=lambda d: d.stat().st_mtime)
        except OSError:
            return
        for directory in finished[:-CLAIM_DIRS_KEPT]:
            shutil.rmtree(directory, ignore_errors=True)

    def stop_claim(self, record: ClaimRecord, reason: str) -> None:
        if record.stop_sent is not None:
            return
        try:
            gpuline.write_json_atomic(record.directory / "control.json",
                                      {"schema": gpuline.CLAIM_CONTROL_SCHEMA, "since_unix": record.control_since,
                                       "release": True, "reason": reason, "updated_unix": time.time()}, mode=0o600)
        except OSError:
            pass
        stop_unit(self.config.launch, record.unit)
        record.stop_sent = time.monotonic()
        self.log(f"stopping GPU claim {record.unit} ({record.phase}): {reason}")

    def stop_all(self, reason: str) -> None:
        for record in list(self.claims.values()):
            self.stop_claim(record, reason)

    def start_claim(self, uuid: str, since: float, job_ids: list[int]) -> ClaimRecord | None:
        token = secrets.token_hex(16)
        unit = claim_unit_name(self.config.launch.unit_prefix, token)
        directory = self.root / token
        gpu = self.gpu_locks(uuid)
        try:
            self.root.mkdir(mode=0o700, parents=True, exist_ok=True)
            directory.mkdir(mode=0o700)
            Path(self.cfg.ticket_dir).mkdir(parents=True, exist_ok=True)
            spec = {
                "schema": gpuline.CLAIM_SPEC_SCHEMA,
                "token": token,
                "unit": unit,
                "gpu": gpu.to_dict(),
                "line": self.cfg.to_dict(),
                "since_unix": since,
                "secret": secrets.token_hex(32),
                "socket_name": f"simq-gpuclaim-{token}",
                "status_path": str(directory / "status.json"),
                "control_path": str(directory / "control.json"),
                "hold_max_s": self.cfg.claim_hold_max_s,
                "created_unix": time.time(),
                "for_jobs": job_ids,
                "worker_version": __version__,
                "worker_runtime": self.runtime_dir,
            }
            gpuline.write_json_atomic(directory / "control.json",
                                      {"schema": gpuline.CLAIM_CONTROL_SCHEMA, "since_unix": since, "release": False,
                                       "updated_unix": time.time()}, mode=0o600)
            gpuline.write_json_atomic(directory / "spec.json", spec, mode=0o600)
            command = build_claim_command(self.config.launch, unit=unit,
                                          description=f"simq GPU claim for {gpu.name or uuid}",
                                          spec_path=str(directory / "spec.json"),
                                          log_path=str(directory / "claim.log"),
                                          cpus=self.cfg.claim_cpus)
            start_unit(command)
        except (OSError, UnitError) as error:
            self.log(f"cannot start a GPU claim for {uuid}: {error}")
            shutil.rmtree(directory, ignore_errors=True)
            return None
        record = ClaimRecord(token=token, uuid=uuid, unit=unit, directory=directory, spec=spec, created=time.time(),
                             control_since=since)
        self.claims[uuid] = record
        self.log(f"started GPU claim {unit} for {gpu.name or uuid}: line since {_iso(since)}, "
                 f"queued GPU job(s) {', '.join(map(str, job_ids))}")
        return record

    def _set_since(self, record: ClaimRecord, since: float) -> None:
        if record.control_since is not None and abs(record.control_since - since) < 1e-6:
            return
        try:
            gpuline.write_json_atomic(record.directory / "control.json",
                                      {"schema": gpuline.CLAIM_CONTROL_SCHEMA, "since_unix": since, "release": False,
                                       "updated_unix": time.time()}, mode=0o600)
            record.control_since = since
        except OSError as error:
            self.log(f"cannot update GPU claim {record.unit}: {error}")

    # ----- planning --------------------------------------------------------------

    def prepare(
        self,
        *,
        jobs: list[dict[str, Any]],
        running: list[Running],
        busy_gpus: dict[str, list],
        now: float,
    ) -> GpuView:
        """Which GPUs queued jobs may not take directly this pass, and why."""
        view = GpuView()
        if not self.enabled:
            for uuid, found in busy_gpus.items():
                view.blocked[uuid] = "held outside the queue"
            return view
        queued_gpu = [j for j in jobs if j["state"] == "queued" and int(j.get("gpus_requested") or 0) > 0]
        for uuid in self.config.gpus:
            tickets = gpuline.load_tickets(self.cfg.ticket_dir, gpu=uuid, log_path=self.cfg.log_path)
            view.tickets[uuid] = tickets
            claim = self.claims.get(uuid)
            allocated = any(uuid in r.allocation.gpus for r in running)
            if claim is not None and claim.handoff_job is not None:
                continue  # the GPU is allocated to the job being started
            if claim is not None and claim.held:
                view.reserved.add(uuid)
                continue
            if claim is not None and claim.active:
                view.blocked[uuid] = f"the queue's GPU claim is {claim.phase}"
                continue
            if uuid in busy_gpus and not allocated:
                view.blocked[uuid] = "held outside the queue"
                continue
            parts = {p.name for p in self.config.partitions.values() if uuid in p.gpus}
            candidates = [created_unix(j) for j in queued_gpu
                          if self.config.canonical_partition(j.get("partition")) in parts]
            if not candidates:
                continue
            oldest = min(candidates)
            early = [t for t in tickets if t.get("kind") == "helper"
                     and t.get("state") not in ("held", "running", "handing-off")
                     and float(t.get("since_unix") or 0) < oldest]
            if early:
                first = early[0]
                view.blocked[uuid] = (f"yielding to gpu-lease pid {first.get('pid')} "
                                      f"(waiting since {_iso(float(first.get('since_unix') or 0))})")
        return view

    def filter_starts(
        self, starts: list[tuple[Request, Allocation]], view: GpuView
    ) -> tuple[list[tuple[Request, Allocation]], dict[int, str]]:
        """At most one start per claimed GPU per pass (a claim hands off once)."""
        kept: list[tuple[Request, Allocation]] = []
        deferred: dict[int, str] = {}
        used: set[str] = set()
        for request, allocation in starts:
            claimed = [g for g in allocation.gpus if g in view.reserved]
            if claimed:
                if len(allocation.gpus) != 1 or claimed[0] in used:
                    deferred[request.job_id] = "waiting for the GPU hand-off of the queue's claim"
                    continue
                used.add(claimed[0])
                record = self.claims[claimed[0]]
                view.handoffs[request.job_id] = {
                    "claim_unit": record.unit,
                    "claim_token": record.token,
                    "socket_name": record.spec["socket_name"],
                    "secret": record.spec["secret"],
                    "gpu": claimed[0],
                    "timeout_s": 20,
                }
            kept.append((request, allocation))
        return kept, deferred

    def launched(self, job: dict[str, Any], view: GpuView, now: float) -> None:
        job_id = int(job["id"])
        handoff = view.handoffs.get(job_id)
        waited = self.ready_since.pop(job_id, None)
        self.reasons.pop(job_id, None)
        if handoff is None:
            return
        record = self.claims.get(str(handoff["gpu"]))
        if record is not None:
            record.handoff_job = job_id
            record.handoff_at = time.time()
        detail = {"claim": handoff["claim_unit"]}
        if waited is not None:
            detail["waited_s"] = round(now - waited, 1)
        try:
            self.store.add_event(job_id, "gpu-handoff", detail)
        except Exception:  # noqa: BLE001 - history only
            pass

    def after_plan(
        self,
        view: GpuView,
        *,
        jobs: list[dict[str, Any]],
        running: list[Running],
        busy_cpus: dict[int, list],
        requests: list[Request],
        starts: list[tuple[Request, Allocation]],
        now: float,
        paused: bool,
    ) -> None:
        """Readiness, claims, published state (before the starts are launched)."""
        if not self.enabled:
            return
        by_id = {int(j["id"]): j for j in jobs}
        started = {request.job_id for request, _ in starts}
        owners = [r for r in requests if not r.borrow]
        for uuid in self.config.gpus:
            if paused or not owners:
                view.ready[uuid] = []
                continue
            view.ready[uuid] = ready_except_gpu(self.config, uuid, running=running, owners=owners,
                                                busy_cpus=busy_cpus, blocked=view.blocked, started=started, now=now)
        ready_all = {job_id for ids in view.ready.values() for job_id in ids}
        queued_ids = {int(j["id"]) for j in jobs if j["state"] == "queued"}
        for job_id in list(self.ready_since):
            if job_id not in queued_ids:
                self.ready_since.pop(job_id, None)
                self.reasons.pop(job_id, None)
        for job_id in ready_all:
            if job_id not in self.ready_since:
                self.ready_since[job_id] = now
        for job_id in list(self.reasons):
            if job_id not in ready_all:
                self.reasons.pop(job_id, None)
        self._manage(view, by_id, now, paused)
        self._snapshots(view, now)
        for uuid, ids in view.ready.items():
            snapshot = view.snapshots.get(uuid) or {}
            holder = self._holder_summary(uuid, snapshot, running, by_id)
            claim = self.claims.get(uuid)
            for position, job_id in enumerate(ids, start=1):
                waited = now - self.ready_since.get(job_id, now)
                text = (f"waiting for the GPU (line position {position} of the queue, "
                        f"{gpuline.duration_text(waited)}"
                        + (f"; claim {claim.phase}" if claim is not None else "")
                        + (f"; GPU held by {holder}" if holder else "") + ")")
                self.reasons[job_id] = text
                if (job_id, "gpu-wait") not in self._noted:
                    self._note_once(job_id, "gpu-wait", f"job {job_id}: {text}")
                    if not any(e.get("kind") == "gpu-wait" for e in self._events(job_id)):
                        self._event(job_id, "gpu-wait", {"holder": holder or "-",
                                                         "claim": claim.unit if claim else None})
                if waited > self.cfg.starvation_warn_s:
                    self._note_once(job_id, "gpu-starving",
                                    f"WARNING: job {job_id} has waited {gpuline.duration_text(waited)} for the GPU "
                                    f"(threshold {gpuline.duration_text(self.cfg.starvation_warn_s)}); "
                                    f"GPU held by {holder or 'nobody'}")
        self.publish(view, jobs=jobs, running=running, starts=starts, now=now, paused=paused)

    def _events(self, job_id: int) -> list[dict[str, Any]]:
        try:
            return self.store.job_events(job_id)
        except Exception:  # noqa: BLE001
            return []

    def _event(self, job_id: int, kind: str, detail: dict[str, Any]) -> None:
        try:
            self.store.add_event(job_id, kind, detail)
        except Exception:  # noqa: BLE001 - history only
            pass

    def _manage(self, view: GpuView, by_id: dict[int, dict[str, Any]], now: float, paused: bool) -> None:
        for uuid in self.config.gpus:
            ids = [] if paused else view.ready.get(uuid, [])
            record = self.claims.get(uuid)
            if record is not None and (record.handoff_job is not None or record.stop_sent is not None):
                continue
            if not ids:
                if record is None:
                    continue
                if paused:
                    self.stop_claim(record, "the queue is paused")
                elif record.held:
                    record.idle_since = record.idle_since or now
                    if now - record.idle_since >= CLAIM_IDLE_GRACE_S:
                        self.stop_claim(record, "no queued job can use the GPU now")
                else:
                    self.stop_claim(record, "no queued GPU job is waiting for the GPU")
                continue
            since = min(created_unix(by_id[job_id]) for job_id in ids if job_id in by_id)
            if record is None:
                record = self.start_claim(uuid, since, ids)
                if record is not None:
                    for job_id in ids:
                        self._event(job_id, "gpu-claim", {"claim": record.unit, "line_since": _iso(since)})
                continue
            record.idle_since = None
            self._set_since(record, since)

    def idle_pass(self, *, jobs: list[dict[str, Any]], running: list[Running], now: float, paused: bool) -> None:
        """A pass without planning (nothing queued, or paused)."""
        if not self.enabled:
            return
        view = GpuView()
        for uuid in self.config.gpus:
            view.tickets[uuid] = gpuline.load_tickets(self.cfg.ticket_dir, gpu=uuid, log_path=self.cfg.log_path)
            view.ready[uuid] = []
        by_id = {int(j["id"]): j for j in jobs}
        queued = {int(j["id"]) for j in jobs if j["state"] == "queued"}
        for job_id in list(self.ready_since):
            if job_id not in queued:
                self.ready_since.pop(job_id, None)
        if paused:
            # paused: keep the since-times of waiting jobs, publish nothing to yield to
            pass
        self.reasons.clear()
        self._manage(view, by_id, now, paused)
        self._snapshots(view, now)
        self.publish(view, jobs=jobs, running=running, starts=[], now=now, paused=paused)

    # ----- published state -----------------------------------------------------------

    def _snapshots(self, view: GpuView, now: float) -> None:
        for uuid in self.config.gpus:
            try:
                view.snapshots[uuid] = gpuline.gpu_snapshot(
                    self.cfg, self.gpu_locks(uuid), tickets=view.tickets.get(uuid), queue_state=None, now=now,
                    scan_fds=False)
            except Exception as error:  # noqa: BLE001 - diagnostics only
                view.snapshots[uuid] = {"error": str(error), "holders": [], "line": [], "raw_waiters": []}

    def _holder_summary(self, uuid: str, snapshot: dict[str, Any], running: list[Running],
                        by_id: dict[int, dict[str, Any]]) -> str:
        jobs = [r.job_id for r in running if uuid in r.allocation.gpus]
        if jobs:
            return "queue job " + ", ".join(str(j) for j in jobs)
        texts = [gpuline.holder_text(h) for h in snapshot.get("holders") or []]
        return "; ".join(texts[:2])

    def publish(
        self,
        view: GpuView,
        *,
        jobs: list[dict[str, Any]],
        running: list[Running],
        starts: list[tuple[Request, Allocation]],
        now: float,
        paused: bool,
    ) -> None:
        if not self.enabled or not self.cfg.state_path:
            return
        by_id = {int(j["id"]): j for j in jobs}
        gpus: dict[str, Any] = {}
        for uuid, gpu in self.config.gpus.items():
            entries: list[dict[str, Any]] = []
            starting_ids = [request.job_id for request, allocation in starts
                            if uuid in allocation.gpus and request.job_id not in view.handoffs]
            starting_ids += [int(j["id"]) for j in jobs if j["state"] == "starting"
                             and uuid in str(j.get("gpu_list") or "").split(",") and int(j["id"]) not in starting_ids]
            ordered = starting_ids + [i for i in view.ready.get(uuid, []) if i not in starting_ids]
            claim = self.claims.get(uuid)
            for job_id in ordered:
                job = by_id.get(job_id)
                if job is None:
                    continue
                gpu_mem = job.get("gpu_mem_bytes")
                entries.append({
                    "job_id": job_id,
                    "name": job.get("name"),
                    "partition": job.get("partition") or self.config.default_partition,
                    "mode": "exclusive" if gpu_mem is None else "shared",
                    "gpu_mem": None if gpu_mem is None else format_size(int(gpu_mem)),
                    "since_unix": created_unix(job),
                    "since": job.get("created_at"),
                    "ready_since_unix": self.ready_since.get(job_id),
                    "waited_s": None if job_id not in self.ready_since else round(now - self.ready_since[job_id], 1),
                    "state": "starting" if job_id in starting_ids else "waiting",
                })
            snapshot = view.snapshots.get(uuid) or {}
            on_gpu = []
            for item in running:
                if uuid in item.allocation.gpus:
                    job = by_id.get(int(item.job_id or 0)) or {}
                    on_gpu.append({"job_id": item.job_id, "partition": job.get("partition"),
                                   "mode": "exclusive" if item.allocation.gpu_exclusive else "shared",
                                   "claim": format_size(item.allocation.gpu_claim), "started_at": job.get("started_at")})
            holders = snapshot.get("holders") or []
            if on_gpu:
                holders = [{"kind": "job", "job_id": j["job_id"], "mode": j["mode"]} for j in on_gpu]
            warnings = []
            for entry in entries:
                waited = entry.get("waited_s")
                if waited is not None and waited > self.cfg.starvation_warn_s:
                    warnings.append(f"queue job {entry['job_id']} has waited {gpuline.duration_text(waited)} for the "
                                    f"GPU (threshold {gpuline.duration_text(self.cfg.starvation_warn_s)})")
            helpers = [{k: t.get(k) for k in ("pid", "state", "since_unix", "label", "command", "note",
                                              "held_since_unix", "cores")}
                       for t in view.tickets.get(uuid, []) if t.get("kind") == "helper"]
            gpus[uuid] = {
                "index": gpu.index,
                "name": gpu.name,
                "lock_order": list(self.config.locks.gpu_lock_order(gpu)),
                "jobs": on_gpu,
                "claim": None if claim is None else claim.brief(),
                "queue": entries,
                "line_since_unix": min((e["since_unix"] for e in entries), default=None),
                "holders": holders,
                "helpers": helpers,
                "raw_waiters": snapshot.get("raw_waiters") or [],
                "blocked": view.blocked.get(uuid),
                "warnings": warnings,
            }
        body = {
            "schema": gpuline.STATE_SCHEMA,
            "version": __version__,
            "node": self.config.node,
            "writer": {"pid": os.getpid(), "runtime": self.runtime_dir},
            "stale_after_s": self.cfg.stale_after_s,
            "starvation_warn_s": self.cfg.starvation_warn_s,
            "paused": paused,
            "ticket_dir": self.cfg.ticket_dir,
            "log_path": self.cfg.log_path,
            "gpus": gpus,
        }
        text = json.dumps(body, sort_keys=True, default=str)
        if text == self._published_body and now - self._published_at < HEARTBEAT_S:
            return
        document = dict(body, heartbeat_unix=now, updated_at=_iso(now))
        try:
            gpuline.write_json_atomic(Path(self.cfg.state_path), document, mode=0o644)
            self._published_body = text
            self._published_at = now
            self._publish_error = None
        except OSError as error:
            if str(error) != self._publish_error:
                self._publish_error = str(error)
                self.log(f"cannot publish the GPU line state {self.cfg.state_path}: {error}")

    def withdraw(self) -> None:
        """Mark the published state stale (worker stopping): heartbeat 0 tells
        helpers not to wait for queue jobs while no worker runs."""
        if not self.enabled or not self.cfg.state_path:
            return
        document = gpuline.load_queue_state(self.cfg.state_path)
        if document is None:
            return
        document["heartbeat_unix"] = 0
        document["withdrawn_at"] = _iso(time.time())
        try:
            gpuline.write_json_atomic(Path(self.cfg.state_path), document, mode=0o644)
        except OSError:
            pass


def describe_job(state: dict[str, Any] | None, job_id: int, now: float | None = None) -> str | None:
    """``simq show``: the job's place in the published GPU line."""
    if state is None:
        return None
    now = time.time() if now is None else now
    for uuid, item in (state.get("gpus") or {}).items():
        entries = item.get("queue") or []
        for position, entry in enumerate(entries, start=1):
            if int(entry.get("job_id") or -1) != job_id:
                continue
            waited = entry.get("ready_since_unix")
            claim = item.get("claim") or {}
            text = (f"{entry.get('state')}; position {position} of {len(entries)} queued GPU job(s) on {uuid}; "
                    f"waiting for the GPU {gpuline.duration_text(None if waited is None else now - float(waited))}")
            if claim:
                text += f"; claim {claim.get('unit')} {claim.get('phase')}"
            holders = item.get("holders") or []
            if holders:
                text += "; GPU held by " + "; ".join(gpuline.holder_text(h, now) for h in holders[:2])
            if not gpuline.queue_state_fresh(state, now):
                text += " (published state is stale)"
            return text
    return None


def iter_waiting(state: dict[str, Any] | None) -> Iterable[tuple[str, dict[str, Any]]]:
    for uuid, item in ((state or {}).get("gpus") or {}).items():
        for entry in item.get("queue") or []:
            yield uuid, entry
