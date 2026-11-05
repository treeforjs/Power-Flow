"""SimQueue 0.6.2 GPU line: fair hand-off of the GPU lease files.

The GPU is guarded by flock(2) lease files (on zenobia ``testing-gpu.lock``
and ``testing-gpu-<UUID>.lock``, taken in that order).  Before 0.6.2 the
scheduler only *observed* them, so a queued GPU job lost every release to
processes blocked in ``flock``: the kernel hands a released flock to the
first blocked waiter at once, and the scheduler looked again a second later.

The GPU line fixes this with three parts that share this module:

* **Members.**  Everyone who wants the GPU in order is a *member* of the
  line and publishes a ticket (``<ticket_dir>/<kind>-<pid>.json``): the
  queue's *claim* (a small systemd unit started by the worker on behalf of
  the queued GPU jobs) and ``gpu-lease`` helper users.  Members are ordered
  by their ``since`` time (the claim's is the submit time of the oldest
  queued job that is ready except for the GPU; a helper's is its start
  time; ties go to the queue).  A member enters the kernel's wait queue on
  the first GPU lock only when every earlier member is already in it (or
  holds the GPU), so the kernel's FIFO order equals the line order.  It
  leaves the kernel queue again when an earlier member appears.
* **Raw waiters go last.**  Processes blocked in plain ``flock`` on the GPU
  lock (not members) that sit ahead of a member in the kernel queue are
  moved to its tail: SIGSTOP interrupts their blocked flock(2), SIGCONT
  restarts it, and the restarted request queues behind every member.  Only
  processes whose command is listed in ``requeue_commands`` (default
  ``flock``) and that are not job-control jobs of an interactive shell are
  touched; holders are never touched.  Every requeue is logged.
* **Hand-off.**  A claim that holds the GPU locks passes the open file
  descriptions to the job unit over an abstract unix socket (SCM_RIGHTS),
  so the lock is never released between the claim and the job; the job's
  time limit starts when its unit starts, after the GPU was won.

Lock order: the GPU lock files in configured order (global file first, then
the per-GPU file), then core leases.  Only the first GPU lock is waited for
in the kernel; later locks are polled without blocking, and a member that
finds the holder of a lock it needs blocked on a lock the member holds
(reversed order elsewhere) backs off: it releases everything and re-enters
the line with its original ``since``.

This module imports nothing from the simqueue package: the claim units and
the helper keep working whatever happens to the scheduler.

    python gpuline.py claim SPEC.json        (main process of a claim unit)
    python gpuline.py lease [options] -- CMD (the gpu-lease helper)
    python gpuline.py status [--json]        (holders and line of each GPU)
"""

from __future__ import annotations

import argparse
import ctypes
import fcntl
import json
import os
import random
import select
import signal
import socket
import struct
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable, Iterable


VERSION = "0.6.3"
TICKET_SCHEMA = "simqueue-gpu-ticket-v1"
STATE_SCHEMA = "simqueue-gpu-line-v1"
CLAIM_SPEC_SCHEMA = "simqueue-gpu-claim-spec-v1"
CLAIM_STATUS_SCHEMA = "simqueue-gpu-claim-status-v1"
CLAIM_CONTROL_SCHEMA = "simqueue-gpu-claim-control-v1"
HANDOFF_PROTOCOL = "simqueue-gpu-handoff-v1"
EXIT_TIMEOUT = 75  # as any-core.sh when it gives up
EXIT_USAGE = 2
PR_SET_PDEATHSIG = 1
KIND_RANK = {"queue": 0, "helper": 1}
# Ticket states.  A member in one of KERNEL_STATES is in the kernel's wait
# queue of the first GPU lock or holds the GPU; later members may then queue
# behind it.
WAITING_STATES = ("waiting", "entering")
KERNEL_STATES = ("queued", "acquiring", "held", "running", "handing-off")
HOLDING_STATES = ("acquiring", "held", "running", "handing-off")
SHELLS = {"bash", "sh", "zsh", "fish", "dash", "ksh", "mksh", "tcsh", "csh", "yash", "xonsh", "nu", "elvish"}
DEFAULT_CONFIG_CANDIDATES = (
    "/home/treefor/.local/state/simqueue/node.json",
)
LOG_LIMIT_BYTES = 2 * 1024 * 1024


# --------------------------------------------------------------------------- small utilities


def now_iso(unix: float | None = None) -> str:
    moment = datetime.fromtimestamp(time.time() if unix is None else unix, UTC)
    return moment.isoformat(timespec="seconds").replace("+00:00", "Z")


def parse_iso(value: Any) -> float | None:
    if not value:
        return None
    try:
        return datetime.fromisoformat(str(value).replace("Z", "+00:00")).timestamp()
    except ValueError:
        return None


def duration_text(seconds: float | None) -> str:
    if seconds is None:
        return "-"
    seconds = max(0, int(seconds))
    if seconds < 60:
        return f"{seconds}s"
    minutes, secs = divmod(seconds, 60)
    if minutes < 60:
        return f"{minutes}m{secs:02d}s"
    hours, minutes = divmod(minutes, 60)
    if hours < 48:
        return f"{hours}h{minutes:02d}m"
    days, hours = divmod(hours, 24)
    return f"{days}d{hours:02d}h"


def write_json_atomic(path: Path, payload: dict[str, Any], mode: int = 0o644) -> None:
    path = Path(path)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    data = (json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n").encode("utf-8")
    descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_TRUNC | os.O_CLOEXEC, mode)
    try:
        os.write(descriptor, data)
        os.fchmod(descriptor, mode)
    finally:
        os.close(descriptor)
    os.replace(temporary, path)


def read_json(path: Path | str) -> dict[str, Any] | None:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return value if isinstance(value, dict) else None


def append_log(path: str | None, text: str) -> None:
    """Append one line to a shared log (O_APPEND; rotated at 2 MiB)."""
    if not path:
        return
    try:
        target = Path(path)
        try:
            if target.stat().st_size > LOG_LIMIT_BYTES:
                os.replace(target, target.with_name(target.name + ".1"))
        except OSError:
            pass
        descriptor = os.open(target, os.O_WRONLY | os.O_CREAT | os.O_APPEND | os.O_CLOEXEC, 0o644)
        try:
            os.write(descriptor, f"{now_iso()} [{os.getpid()}] {text}\n".encode("utf-8", "replace"))
        finally:
            os.close(descriptor)
    except OSError:
        pass


# --------------------------------------------------------------------------- configuration


@dataclass(frozen=True)
class GpuLocks:
    uuid: str
    index: int | None
    name: str
    paths: tuple[str, ...]  # acquisition order: global file(s) first, then the per-GPU file

    def to_dict(self) -> dict[str, Any]:
        return {"uuid": self.uuid, "index": self.index, "name": self.name, "paths": list(self.paths)}

    @staticmethod
    def from_dict(value: dict[str, Any]) -> "GpuLocks":
        return GpuLocks(str(value["uuid"]), value.get("index"), str(value.get("name", "")),
                        tuple(str(p) for p in value["paths"]))


@dataclass(frozen=True)
class LineConfig:
    enabled: bool = True
    state_path: str = ""
    ticket_dir: str = ""
    log_path: str = ""
    starvation_warn_s: float = 600.0
    claim_hold_max_s: float = 300.0
    stale_after_s: float = 120.0
    requeue_raw_waiters: bool = True
    requeue_commands: tuple[str, ...] = ("flock",)
    poll_s: float = 0.1
    helper_cores: tuple[int, ...] = (6, 7)  # gpu-lease --core any|both
    claim_cpus: tuple[int, ...] = ()  # CPUAffinity of claim units (empty: unpinned, like the worker)
    core_paths: dict[int, str] = field(default_factory=dict)
    gpus: dict[str, GpuLocks] = field(default_factory=dict)

    def gpu(self, selector: str | None = None) -> GpuLocks:
        if not self.gpus:
            raise ValueError("no GPU is configured")
        if selector in (None, ""):
            return next(iter(self.gpus.values()))
        if selector in self.gpus:
            return self.gpus[selector]
        for gpu in self.gpus.values():
            if gpu.index is not None and str(gpu.index) == str(selector):
                return gpu
        raise ValueError(f"unknown GPU {selector!r} (configured: {', '.join(self.gpus)})")

    def to_dict(self) -> dict[str, Any]:
        return {
            "enabled": self.enabled,
            "state_path": self.state_path,
            "ticket_dir": self.ticket_dir,
            "log_path": self.log_path,
            "starvation_warn_s": self.starvation_warn_s,
            "claim_hold_max_s": self.claim_hold_max_s,
            "stale_after_s": self.stale_after_s,
            "requeue_raw_waiters": self.requeue_raw_waiters,
            "requeue_commands": list(self.requeue_commands),
            "poll_s": self.poll_s,
            "helper_cores": list(self.helper_cores),
            "claim_cpus": list(self.claim_cpus),
            "core_paths": {str(k): v for k, v in self.core_paths.items()},
            "gpus": {k: v.to_dict() for k, v in self.gpus.items()},
        }

    @staticmethod
    def from_dict(value: dict[str, Any]) -> "LineConfig":
        return LineConfig(
            enabled=bool(value.get("enabled", True)),
            state_path=str(value.get("state_path", "")),
            ticket_dir=str(value.get("ticket_dir", "")),
            log_path=str(value.get("log_path", "")),
            starvation_warn_s=float(value.get("starvation_warn_s", 600.0)),
            claim_hold_max_s=float(value.get("claim_hold_max_s", 300.0)),
            stale_after_s=float(value.get("stale_after_s", 120.0)),
            requeue_raw_waiters=bool(value.get("requeue_raw_waiters", True)),
            requeue_commands=tuple(str(c) for c in value.get("requeue_commands", ["flock"])),
            poll_s=float(value.get("poll_s", 0.1)),
            helper_cores=tuple(int(c) for c in value.get("helper_cores", [6, 7])),
            claim_cpus=tuple(int(c) for c in value.get("claim_cpus", [])),
            core_paths={int(k): str(v) for k, v in dict(value.get("core_paths", {})).items()},
            gpus={k: GpuLocks.from_dict(v) for k, v in dict(value.get("gpus", {})).items()},
        )


LINE_KEYS = {
    "enabled", "state_path", "ticket_dir", "log_path", "starvation_warn_s", "claim_hold_max_s",
    "stale_after_s", "requeue_raw_waiters", "requeue_commands", "poll_s", "helper_cores", "claim_cpus", "notes",
}


def parse_line_config(document: dict[str, Any]) -> LineConfig:
    """The ``gpu_line`` section of a node configuration (raises ValueError)."""
    raw = document.get("gpu_line", {})
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise ValueError("gpu_line must be an object")
    unknown = set(raw) - LINE_KEYS
    if unknown:
        raise ValueError(f"gpu_line: unknown key(s) {', '.join(sorted(unknown))}")
    locks = document.get("lock_files") or {}
    if not isinstance(locks, dict):
        raise ValueError("lock_files must be an object")
    gpu_path = locks.get("gpu_path")
    extras = locks.get("gpu_extra") or {}
    gpus: dict[str, GpuLocks] = {}
    for item in document.get("gpus") or []:
        if not isinstance(item, dict) or "uuid" not in item:
            continue
        uuid = str(item["uuid"])
        index = item.get("index")
        paths: list[str] = []
        for extra in extras.get(uuid, []) or []:
            if extra not in paths:
                paths.append(str(extra))
        if gpu_path:
            own = str(gpu_path).replace("{uuid}", uuid).replace("{index}", "" if index is None else str(index))
            if own not in paths:
                paths.append(own)
        gpus[uuid] = GpuLocks(uuid, index if isinstance(index, int) else None, str(item.get("name", "")),
                              tuple(paths))
    core_paths: dict[int, str] = {}
    overrides = locks.get("core_overrides") or {}
    core_template = locks.get("core_path")
    cpus: set[int] = set()
    for part in document.get("partitions") or []:
        value = part.get("cpus", part.get("cores")) if isinstance(part, dict) else None
        if isinstance(value, str):
            for chunk in value.split(","):
                chunk = chunk.strip()
                if not chunk:
                    continue
                low, _, high = chunk.partition("-")
                if low.isdigit() and (not high or high.isdigit()):
                    cpus.update(range(int(low), int(high or low) + 1))
        elif isinstance(value, list):
            cpus.update(int(v) for v in value if isinstance(v, int))
    for cpu in sorted(cpus):
        if str(cpu) in overrides:
            core_paths[cpu] = str(overrides[str(cpu)])
        elif core_template:
            core_paths[cpu] = str(core_template).replace("{core}", str(cpu))
    base: str | None = None
    for candidate in [p for g in gpus.values() for p in g.paths] + list(core_paths.values()):
        base = os.path.dirname(candidate)
        break
    enabled_default = bool(locks.get("enabled")) and bool(gpus) and base is not None
    enabled = raw.get("enabled", enabled_default)
    if not isinstance(enabled, bool):
        raise ValueError("gpu_line.enabled must be true or false")
    if enabled and not enabled_default:
        raise ValueError("gpu_line.enabled needs lock_files.enabled and at least one GPU with lease files")

    def path_value(key: str, default: str | None) -> str:
        value = raw.get(key, default)
        if value is None:
            return ""
        if not isinstance(value, str) or not value.startswith("/"):
            raise ValueError(f"gpu_line.{key} must be an absolute path")
        return value

    def number(key: str, default: float, low: float, high: float) -> float:
        value = raw.get(key, default)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not low <= float(value) <= high:
            raise ValueError(f"gpu_line.{key} must be a number in [{low}, {high}]")
        return float(value)

    requeue = raw.get("requeue_raw_waiters", True)
    if not isinstance(requeue, bool):
        raise ValueError("gpu_line.requeue_raw_waiters must be true or false")
    commands = raw.get("requeue_commands", ["flock"])
    if not isinstance(commands, list) or not all(isinstance(c, str) and c and len(c) <= 15 for c in commands):
        raise ValueError("gpu_line.requeue_commands must be an array of command names (as in /proc/PID/comm)")
    helper_cores = raw.get("helper_cores", [6, 7])
    if not isinstance(helper_cores, list) or not all(isinstance(c, int) and not isinstance(c, bool)
                                                     for c in helper_cores):
        raise ValueError("gpu_line.helper_cores must be an array of CPU numbers")
    claim_cpus = raw.get("claim_cpus", [])
    if not isinstance(claim_cpus, list) or not all(isinstance(c, int) and not isinstance(c, bool) and c >= 0
                                                   for c in claim_cpus):
        raise ValueError("gpu_line.claim_cpus must be an array of CPU numbers")
    return LineConfig(
        enabled=enabled,
        state_path=path_value("state_path", None if base is None else os.path.join(base, "simq-gpu-line.json")),
        ticket_dir=path_value("ticket_dir", None if base is None else os.path.join(base, "gpu-line.d")),
        log_path=path_value("log_path", None if base is None else os.path.join(base, "simq-gpu-line.log")),
        starvation_warn_s=number("starvation_warn_s", 600.0, 1.0, 7 * 86400.0),
        claim_hold_max_s=number("claim_hold_max_s", 300.0, 5.0, 86400.0),
        stale_after_s=number("stale_after_s", 120.0, 5.0, 86400.0),
        requeue_raw_waiters=requeue,
        requeue_commands=tuple(commands),
        poll_s=number("poll_s", 0.1, 0.02, 2.0),
        helper_cores=tuple(helper_cores),
        claim_cpus=tuple(claim_cpus),
        core_paths=core_paths,
        gpus=gpus,
    )


def load_line_config(path: str | None = None) -> tuple[LineConfig, str]:
    """The line configuration from a node.json (helper and status use)."""
    candidates: list[str] = []
    if path:
        candidates.append(path)
    else:
        for variable in ("GPU_LEASE_CONFIG", "SIMQUEUE_CONFIG"):
            if os.environ.get(variable):
                candidates.append(os.environ[variable])
        candidates.extend(DEFAULT_CONFIG_CANDIDATES)
        candidates.append(str(Path(__file__).resolve().parents[2] / "config" / "zenobia-node.json"))
    errors = []
    for candidate in candidates:
        document = read_json(candidate)
        if document is None:
            errors.append(f"{candidate}: unreadable")
            continue
        try:
            return parse_line_config(document), candidate
        except ValueError as error:
            errors.append(f"{candidate}: {error}")
    raise ValueError("no usable node configuration: " + "; ".join(errors))


# --------------------------------------------------------------------------- /proc/locks


@dataclass(frozen=True)
class LockRecord:
    seq: int
    lock_id: int
    kind: str  # FLOCK / POSIX / OFDLCK / LEASE
    mode: str  # WRITE / READ
    pid: int | None
    dev: str  # as printed by /proc/locks, e.g. "103:09"
    ino: int
    waiter: bool
    depth: int


def parse_proc_locks(text: str) -> list[LockRecord]:
    records: list[LockRecord] = []
    for seq, line in enumerate(text.splitlines()):
        head, sep, rest = line.partition(":")
        if not sep or not head.strip().isdigit():
            continue
        stripped = rest.lstrip(" ")
        waiter = stripped.startswith("->")
        depth = 0
        if waiter:
            depth = max(1, len(rest) - len(stripped))
            stripped = stripped[2:]
        fields = stripped.split()
        if len(fields) < 5:
            continue
        kind, mode = fields[0], fields[2]
        try:
            pid: int | None = int(fields[3])
        except ValueError:
            pid = None
        if pid is not None and pid <= 0:
            pid = None
        parts = fields[4].split(":")
        if len(parts) != 3:
            continue
        try:
            ino = int(parts[2])
        except ValueError:
            continue
        records.append(LockRecord(seq, int(head), kind, mode, pid, f"{parts[0]}:{parts[1]}".lower(), ino,
                                  waiter, depth))
    return records


def read_proc_locks() -> list[LockRecord]:
    try:
        return parse_proc_locks(Path("/proc/locks").read_text(encoding="ascii", errors="replace"))
    except OSError:
        return []


def file_identity(path: str) -> tuple[str, int] | None:
    try:
        stat = os.stat(path)
    except OSError:
        return None
    return f"{os.major(stat.st_dev):02x}:{os.minor(stat.st_dev):02x}", stat.st_ino


@dataclass
class PathLocks:
    path: str
    identity: tuple[str, int] | None
    holders: list[LockRecord] = field(default_factory=list)
    waiters: list[LockRecord] = field(default_factory=list)  # kernel queue order (chain order)


def path_locks(paths: Iterable[str], records: list[LockRecord] | None = None) -> dict[str, PathLocks]:
    paths = list(paths)
    result = {path: PathLocks(path, file_identity(path)) for path in paths}
    by_identity: dict[tuple[str, int], list[str]] = {}
    for path, item in result.items():
        if item.identity is not None:
            by_identity.setdefault(item.identity, []).append(path)
    if not by_identity:
        return result
    for record in read_proc_locks() if records is None else records:
        for path in by_identity.get((record.dev, record.ino), ()):
            (result[path].waiters if record.waiter else result[path].holders).append(record)
    return result


# --------------------------------------------------------------------------- processes


@dataclass(frozen=True)
class ProcStat:
    pid: int
    comm: str
    state: str
    ppid: int
    pgrp: int
    session: int
    tty_nr: int
    tpgid: int
    start: int


def proc_stat(pid: int) -> ProcStat | None:
    try:
        text = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8", errors="replace")
    except (OSError, ValueError):
        return None
    try:
        comm = text[text.index("(") + 1 : text.rindex(")")]
        rest = text[text.rindex(")") + 2 :].split()
        return ProcStat(pid, comm, rest[0], int(rest[1]), int(rest[2]), int(rest[3]), int(rest[4]),
                        int(rest[5]), int(rest[19]))
    except (IndexError, ValueError):
        return None


def pid_start(pid: int | None) -> int | None:
    if not pid:
        return None
    stat = proc_stat(pid)
    return None if stat is None else stat.start


def same_process(pid: Any, start: Any) -> bool:
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return False
    stat = proc_stat(pid)
    if stat is None or stat.state in ("Z", "X"):
        return False
    return start is None or stat.start == int(start)


def proc_cmdline(pid: int | None, width: int = 160) -> str:
    if not pid:
        return "?"
    try:
        raw = Path(f"/proc/{pid}/cmdline").read_bytes()
    except OSError:
        return "(exited)"
    command = " ".join(part.decode("utf-8", "replace") for part in raw.split(b"\0") if part)
    if not command:
        stat = proc_stat(pid)
        command = f"[{stat.comm}]" if stat else "?"
    return command if len(command) <= width else command[: width - 1] + "…"


def process_table() -> dict[int, ProcStat]:
    table: dict[int, ProcStat] = {}
    try:
        names = os.listdir("/proc")
    except OSError:
        return table
    for name in names:
        if name.isdigit():
            stat = proc_stat(int(name))
            if stat is not None:
                table[stat.pid] = stat
    return table


def descendants(table: dict[int, ProcStat], roots: Iterable[int]) -> set[int]:
    children: dict[int, list[int]] = {}
    for stat in table.values():
        children.setdefault(stat.ppid, []).append(stat.pid)
    found: set[int] = set()
    stack = [pid for pid in roots if pid in table]
    while stack:
        pid = stack.pop()
        if pid in found:
            continue
        found.add(pid)
        stack.extend(children.get(pid, ()))
    return found


def fd_holders(identities: Iterable[tuple[str, int]]) -> dict[tuple[str, int], set[int]]:
    """Processes of this account that have one of the files open."""
    wanted = set(identities)
    found: dict[tuple[str, int], set[int]] = {key: set() for key in wanted}
    if not wanted:
        return found
    uid = os.getuid()
    try:
        names = os.listdir("/proc")
    except OSError:
        return found
    for name in names:
        if not name.isdigit():
            continue
        base = f"/proc/{name}"
        try:
            if os.stat(base).st_uid != uid:
                continue
            fds = os.listdir(base + "/fd")
        except OSError:
            continue
        for fd in fds:
            try:
                stat = os.stat(f"{base}/fd/{fd}")
            except OSError:
                continue
            key = (f"{os.major(stat.st_dev):02x}:{os.minor(stat.st_dev):02x}", stat.st_ino)
            if key in wanted:
                found[key].add(int(name))
    return found


def send_signal(pid: int, start: int | None, sig: int) -> bool:
    """Signal *pid* only if it is still the process that started at *start*."""
    try:
        handle = os.pidfd_open(pid)
    except OSError:
        return False
    try:
        if start is not None and pid_start(pid) != start:
            return False
        signal.pidfd_send_signal(handle, sig)
        return True
    except OSError:
        return False
    finally:
        os.close(handle)


def job_control_risk(stat: ProcStat) -> bool:
    """True for a job of an interactive job-control shell: stopping it would
    make the shell report it stopped and take back the terminal."""
    if stat.tty_nr == 0:
        return False
    parent = proc_stat(stat.ppid)
    if parent is None:
        return False
    return stat.pgrp != parent.pgrp and parent.comm in SHELLS


def requeue_refusal(pid: int, commands: Iterable[str]) -> str | None:
    """Why a blocked raw waiter must not be requeued, or ``None``."""
    stat = proc_stat(pid)
    if stat is None:
        return "exited"
    try:
        if os.stat(f"/proc/{pid}").st_uid != os.getuid():
            return "another user's process"
    except OSError:
        return "exited"
    if stat.comm not in set(commands):
        return f"command {stat.comm!r} is not in requeue_commands"
    if stat.state not in ("S", "D"):
        return f"not blocked (state {stat.state})"
    if job_control_risk(stat):
        return "job of an interactive shell"
    return None


# --------------------------------------------------------------------------- tickets


def ticket_key(ticket: dict[str, Any]) -> tuple[float, int, int]:
    return (float(ticket.get("since_unix") or 0.0), KIND_RANK.get(str(ticket.get("kind")), 9),
            int(ticket.get("pid") or 0))


def member_pids(tickets: Iterable[dict[str, Any]]) -> set[int]:
    pids: set[int] = set()
    for ticket in tickets:
        for key in ("pid", "waiter_pid", "child_pid"):
            if ticket.get(key):
                pids.add(int(ticket[key]))
    return pids


def recover_stopped(entries: Iterable[dict[str, Any]]) -> list[int]:
    """SIGCONT raw waiters a crashed member left stopped while requeueing them."""
    resumed = []
    for entry in entries or ():
        try:
            pid, start = int(entry["pid"]), entry.get("start")
        except (KeyError, TypeError, ValueError):
            continue
        stat = proc_stat(pid)
        if stat is not None and stat.state in ("T", "t") and (start is None or stat.start == int(start)):
            if send_signal(pid, stat.start, signal.SIGCONT):
                resumed.append(pid)
    return resumed


def load_tickets(ticket_dir: str, *, gpu: str | None = None, cleanup: bool = True,
                 log_path: str | None = None) -> list[dict[str, Any]]:
    """Live tickets in line order.  Tickets of exited members are removed
    (after resuming any raw waiter they left stopped)."""
    tickets: list[dict[str, Any]] = []
    try:
        names = os.listdir(ticket_dir)
    except OSError:
        return tickets
    for name in names:
        if not name.endswith(".json") or name.startswith("."):
            continue
        path = os.path.join(ticket_dir, name)
        data = read_json(path)
        if data is None or data.get("schema") != TICKET_SCHEMA:
            continue
        if not same_process(data.get("pid"), data.get("pid_start")):
            if cleanup:
                resumed = recover_stopped(data.get("stopping") or [])
                if resumed:
                    append_log(log_path, f"resumed raw waiter(s) {resumed} left stopped by exited member "
                                         f"pid {data.get('pid')}")
                try:
                    os.unlink(path)
                except OSError:
                    pass
            continue
        if gpu is not None and data.get("gpu") != gpu:
            continue
        data["_path"] = path
        tickets.append(data)
    tickets.sort(key=ticket_key)
    return tickets


def load_queue_state(path: str) -> dict[str, Any] | None:
    data = read_json(path) if path else None
    if data is None or data.get("schema") != STATE_SCHEMA:
        return None
    return data


def queue_state_fresh(state: dict[str, Any] | None, now: float | None = None) -> bool:
    if state is None:
        return False
    now = time.time() if now is None else now
    try:
        heartbeat = float(state.get("heartbeat_unix") or 0.0)
        stale = float(state.get("stale_after_s") or 120.0)
    except (TypeError, ValueError):
        return False
    return now - heartbeat <= stale


def queue_entries(state: dict[str, Any] | None, gpu: str, now: float | None = None) -> list[dict[str, Any]]:
    """Queue jobs waiting for *gpu* according to a fresh published state."""
    if not queue_state_fresh(state, now) or state.get("paused"):
        return []
    entries = ((state.get("gpus") or {}).get(gpu) or {}).get("queue") or []
    return [e for e in entries if isinstance(e, dict) and e.get("since_unix") is not None]


def turn_blockers(me: dict[str, Any], tickets: Iterable[dict[str, Any]], queue_state: dict[str, Any] | None,
                  now: float | None = None, *, positions: dict[int, int] | None = None) -> list[str]:
    """Earlier line members that are not yet in the kernel queue (empty: it
    is *me*'s turn to enter it).

    With *positions* (waiter pid -> place in the kernel queue of the first
    GPU file; *me* is queued there), an earlier member queued *behind* me
    is a blocker too: I must leave the kernel queue and re-enter behind it.
    """
    key = ticket_key(me)
    gpu = me.get("gpu")
    reasons: list[str] = []
    claim_in_kernel = False
    mine = None if positions is None else positions.get(int(me.get("waiter_pid") or 0))
    for ticket in tickets:
        if ticket.get("gpu") != gpu or int(ticket.get("pid") or 0) == int(me.get("pid") or 0):
            continue
        state = ticket.get("state")
        if ticket.get("kind") == "queue" and state in KERNEL_STATES:
            claim_in_kernel = True
        if ticket_key(ticket) >= key:
            continue
        behind = False
        if mine is not None and state == "queued":
            theirs = positions.get(int(ticket.get("waiter_pid") or 0))
            behind = theirs is not None and theirs > mine
        if state not in KERNEL_STATES or behind:
            reasons.append(f"{ticket.get('kind')} pid {ticket.get('pid')} ({ticket.get('label') or '-'}, "
                           f"since {now_iso(float(ticket.get('since_unix') or 0))}"
                           + (", queued behind me)" if behind else ")"))
    if me.get("kind") == "helper" and not claim_in_kernel:
        for entry in queue_entries(queue_state, str(gpu), now):
            if float(entry["since_unix"]) < float(me.get("since_unix") or 0.0):
                reasons.append(f"queued job {entry.get('job_id')} ({entry.get('name') or '-'}, "
                               f"waiting since {now_iso(float(entry['since_unix']))})")
    return reasons


# --------------------------------------------------------------------------- the member engine


class Stopped(Exception):
    pass


class TimedOut(Exception):
    pass


def _set_pdeathsig() -> None:
    try:
        ctypes.CDLL(None, use_errno=True).prctl(PR_SET_PDEATHSIG, signal.SIGKILL, 0, 0, 0)
    except (AttributeError, OSError):
        pass


class LineMember:
    """One member of the GPU line (a claim or a helper)."""

    def __init__(
        self,
        cfg: LineConfig,
        gpu: GpuLocks,
        *,
        kind: str,
        since_unix: float,
        label: str,
        command: str = "",
        extra: dict[str, Any] | None = None,
        should_stop: Callable[[], bool] = lambda: False,
        since_getter: Callable[[], float | None] | None = None,
        deadline: float | None = None,
        on_state: Callable[[str, dict[str, Any]], None] | None = None,
        log: Callable[[str], None] | None = None,
        requeue: bool | None = None,
    ) -> None:
        if kind not in KIND_RANK:
            raise ValueError(f"unknown member kind {kind!r}")
        self.cfg = cfg
        self.gpu = gpu
        self.kind = kind
        self.should_stop = should_stop
        self.since_getter = since_getter
        self.deadline = deadline
        self.on_state = on_state
        self._log = log
        self.requeue_enabled = cfg.requeue_raw_waiters if requeue is None else requeue
        self.requeued: dict[tuple[int, int | None], float] = {}
        self.refused_noted: set[tuple[int, str]] = set()
        self.requeue_count = 0
        self.fds: list[int] = []
        Path(cfg.ticket_dir).mkdir(parents=True, exist_ok=True)
        self.ticket_path = Path(cfg.ticket_dir) / f"{kind}-{os.getpid()}.json"
        self.ticket: dict[str, Any] = {
            "schema": TICKET_SCHEMA,
            "version": VERSION,
            "kind": kind,
            "gpu": gpu.uuid,
            "lock_paths": list(gpu.paths),
            "pid": os.getpid(),
            "pid_start": pid_start(os.getpid()),
            "since_unix": float(since_unix),
            "since": now_iso(since_unix),
            "label": label,
            "command": command,
            "state": "new",
            "stopping": [],
            **(extra or {}),
        }

    # ----- bookkeeping

    def log(self, text: str) -> None:
        append_log(self.cfg.log_path, f"{self.kind} pid {os.getpid()} ({self.ticket.get('label')}): {text}")
        if self._log is not None:
            self._log(text)

    def write(self) -> None:
        self.ticket["updated_unix"] = time.time()
        write_json_atomic(self.ticket_path, {k: v for k, v in self.ticket.items() if not k.startswith("_")})

    def set_state(self, state: str, **fields: Any) -> None:
        changed = self.ticket.get("state") != state or any(self.ticket.get(k) != v for k, v in fields.items())
        self.ticket["state"] = state
        self.ticket.update(fields)
        if changed:
            self.write()
            if self.on_state is not None:
                self.on_state(state, dict(fields))

    def remove_ticket(self) -> None:
        try:
            self.ticket_path.unlink()
        except OSError:
            pass

    def release(self) -> None:
        for descriptor in self.fds:
            try:
                os.close(descriptor)
            except OSError:
                pass
        self.fds = []

    def close(self) -> None:
        self.release()
        self.remove_ticket()

    def _refresh_since(self) -> None:
        if self.since_getter is None:
            return
        value = self.since_getter()
        if value is not None and abs(float(value) - float(self.ticket["since_unix"])) > 1e-6:
            self.ticket["since_unix"] = float(value)
            self.ticket["since"] = now_iso(float(value))
            self.write()

    def _check(self) -> None:
        if self.should_stop():
            raise Stopped()
        if self.deadline is not None and time.time() > self.deadline:
            raise TimedOut()

    def blockers(self, positions: dict[int, int] | None = None) -> list[str]:
        tickets = load_tickets(self.cfg.ticket_dir, gpu=self.gpu.uuid, log_path=self.cfg.log_path)
        state = load_queue_state(self.cfg.state_path) if self.kind == "helper" else None
        return turn_blockers(self.ticket, tickets, state, positions=positions)

    # ----- acquisition

    def acquire(self) -> list[int]:
        """Wait for the turn, win the first GPU lock in the kernel queue, take
        the others; returns the descriptors (in order).  Raises Stopped or
        TimedOut."""
        self.set_state("waiting")
        while True:
            self._wait_turn()
            first = self._kernel_wait(self.gpu.paths[0])
            if first is None:
                continue
            fds = [first]
            self.set_state("acquiring", waiter_pid=None)
            ok = True
            try:
                for index, path in enumerate(self.gpu.paths[1:], start=1):
                    descriptor = self.acquire_nonblocking(path, self.gpu.paths[:index])
                    if descriptor is None:
                        ok = False
                        break
                    fds.append(descriptor)
            except BaseException:
                for descriptor in fds:
                    os.close(descriptor)
                raise
            if ok:
                self.fds = fds
                self.set_state("held", held_since_unix=time.time(), note=None)
                return list(fds)
            for descriptor in fds:
                os.close(descriptor)
            self.set_state("waiting", note="backed off: lock-order conflict")
            self._sleep(1.0 + random.random())

    def _sleep(self, seconds: float) -> None:
        end = time.monotonic() + seconds
        while time.monotonic() < end:
            self._check()
            time.sleep(min(self.cfg.poll_s, max(0.0, end - time.monotonic())))

    def _wait_turn(self) -> None:
        last_note = None
        while True:
            self._check()
            self._refresh_since()
            reasons = self.blockers()
            if not reasons:
                return
            note = "behind " + "; ".join(reasons[:3])
            if note != last_note:
                self.set_state("waiting", note=note)
                last_note = note
            time.sleep(self.cfg.poll_s)

    def _kernel_wait(self, path: str) -> int | None:
        descriptor = os.open(path, os.O_RDONLY | os.O_CREAT | os.O_CLOEXEC, 0o644)
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            return descriptor
        except BlockingIOError:
            pass
        except BaseException:
            os.close(descriptor)
            raise
        parent = os.getpid()
        child = os.fork()
        if child == 0:  # the waiter: blocks in the kernel queue on the shared open file description
            try:
                _set_pdeathsig()
                if os.getppid() != parent:
                    os._exit(3)
                for number in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP, signal.SIGUSR1):
                    signal.signal(number, signal.SIG_DFL)
                fcntl.flock(descriptor, fcntl.LOCK_EX)
                os._exit(0)
            except BaseException:
                os._exit(2)
        confirmed = False
        waiter: int | None = child
        try:
            self.set_state("entering", waiter_pid=child, waiter_start=pid_start(child))
            while True:
                done, status = os.waitpid(child, os.WNOHANG)
                if done == child:
                    waiter = None
                    if os.waitstatus_to_exitcode(status) == 0:
                        return descriptor
                    os.close(descriptor)
                    time.sleep(0.2)
                    return None
                state = path_locks([path])[path]
                if not confirmed and any(w.pid == child for w in state.waiters):
                    confirmed = True
                    self.set_state("queued", note=None)
                if confirmed and self.requeue_enabled:
                    self._requeue_ahead(path, state, child)
                reason: str | None = None
                stop: type[Exception] | None = None
                if self.should_stop():
                    reason, stop = "stop requested", Stopped
                elif self.deadline is not None and time.time() > self.deadline:
                    reason, stop = "timed out", TimedOut
                elif confirmed:
                    self._refresh_since()
                    positions = {w.pid: index for index, w in enumerate(state.waiters) if w.pid}
                    blockers = self.blockers(positions)
                    if blockers:
                        reason = "an earlier member arrived: " + blockers[0]
                if reason is not None:
                    os.kill(child, signal.SIGKILL)
                    os.waitpid(child, 0)
                    waiter = None
                    try:
                        # the waiter may have won just before it died: the lock is then ours
                        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        if stop is None:
                            return descriptor
                    except BlockingIOError:
                        pass
                    os.close(descriptor)
                    if stop is not None:
                        raise stop()
                    self.set_state("waiting", waiter_pid=None, note=reason)
                    return None
                time.sleep(self.cfg.poll_s)
        except BaseException:
            if waiter is not None:
                try:
                    os.kill(waiter, signal.SIGKILL)
                    os.waitpid(waiter, 0)
                except OSError:
                    pass
            raise

    def _requeue_ahead(self, path: str, state: PathLocks, my_waiter: int) -> None:
        tickets = load_tickets(self.cfg.ticket_dir, log_path=self.cfg.log_path)
        members = member_pids(tickets) | {os.getpid(), my_waiter}
        now = time.monotonic()
        for record in state.waiters:
            if record.pid == my_waiter:
                break
            if record.pid is None or record.pid in members:
                continue
            start = pid_start(record.pid)
            key = (record.pid, start)
            if now - self.requeued.get(key, -1e9) < 0.5:
                continue
            refusal = requeue_refusal(record.pid, self.cfg.requeue_commands)
            if refusal is not None:
                if (record.pid, refusal) not in self.refused_noted:
                    self.refused_noted.add((record.pid, refusal))
                    self.log(f"not requeueing waiter pid {record.pid} ({proc_cmdline(record.pid, 100)}) "
                             f"on {os.path.basename(path)}: {refusal}")
                continue
            if self._requeue(record.pid, start):
                self.requeued[key] = now
                self.requeue_count += 1
                self.log(f"moved raw waiter pid {record.pid} ({proc_cmdline(record.pid, 120)}) "
                         f"behind this member on {os.path.basename(path)}")

    def _requeue(self, pid: int, start: int | None) -> bool:
        """SIGSTOP + SIGCONT: the blocked flock(2) is abandoned and restarted,
        so the request joins the tail of the kernel queue."""
        entry = {"pid": pid, "start": start}
        self.ticket.setdefault("stopping", []).append(entry)
        self.write()
        try:
            if not send_signal(pid, start, signal.SIGSTOP):
                return False
            end = time.monotonic() + 0.5
            while time.monotonic() < end:
                stat = proc_stat(pid)
                if stat is None or (start is not None and stat.start != start) or stat.state in ("T", "t"):
                    break
                time.sleep(0.002)
            # Let the requests chained behind the stopped one (this member's
            # among them) re-queue before it restarts its flock(2).
            time.sleep(0.03)
            return True
        finally:
            send_signal(pid, start, signal.SIGCONT)
            try:
                self.ticket["stopping"].remove(entry)
            except ValueError:
                pass
            self.write()

    def acquire_nonblocking(self, path: str, holding: Iterable[str], *, shared: bool = False,
                            note: str | None = None) -> int | None:
        """Poll *path* while holding *holding*; ``None`` when the holder of
        *path* waits for one of *holding* (a lock-order conflict)."""
        holding = list(holding)
        descriptor = os.open(path, os.O_RDONLY | os.O_CREAT | os.O_CLOEXEC, 0o644)
        operation = (fcntl.LOCK_SH if shared else fcntl.LOCK_EX) | fcntl.LOCK_NB
        noted = False
        last_check = 0.0
        try:
            while True:
                try:
                    fcntl.flock(descriptor, operation)
                    return descriptor
                except BlockingIOError:
                    pass
                self._check()
                if not noted:
                    self.set_state(str(self.ticket.get("state")), note=note or f"waiting for {path}")
                    noted = True
                if time.monotonic() - last_check > 0.5:
                    last_check = time.monotonic()
                    conflict = lock_order_conflict(path, holding)
                    if conflict:
                        self.log(f"lock-order conflict on {path}: {conflict}; releasing and re-entering the line")
                        os.close(descriptor)
                        return None
                time.sleep(max(self.cfg.poll_s, 0.1))
        except BaseException:
            os.close(descriptor)
            raise


def lock_order_conflict(path: str, holding: Iterable[str]) -> str | None:
    """Describe a holder of *path* that is blocked on a lock in *holding*."""
    holding = list(holding)
    if not holding:
        return None
    records = read_proc_locks()
    locks = path_locks([path, *holding], records)
    waiting = {w.pid: p for p in holding for w in locks[p].waiters if w.pid}
    if not waiting:
        return None
    table = process_table()
    roots = {r.pid for r in locks[path].holders if r.pid and r.pid in table}
    if any(r.pid and r.pid not in table for r in locks[path].holders) and locks[path].identity is not None:
        roots |= fd_holders([locks[path].identity]).get(locks[path].identity, set())
    tree = descendants(table, roots)
    for pid in sorted(tree & set(waiting)):
        return f"pid {pid} ({proc_cmdline(pid, 80)}) holds it and waits for {waiting[pid]}"
    return None


# --------------------------------------------------------------------------- describing a GPU


def gpu_snapshot(
    cfg: LineConfig,
    gpu: GpuLocks,
    *,
    tickets: list[dict[str, Any]] | None = None,
    queue_state: dict[str, Any] | None = None,
    records: list[LockRecord] | None = None,
    now: float | None = None,
    cleanup: bool = False,
    scan_fds: bool = True,
) -> dict[str, Any]:
    """Who holds *gpu*'s lease files and who waits for them, from
    /proc/locks, the tickets and (for queue jobs) the published state.
    With *scan_fds* a lock whose locker exited is attributed to the
    processes that have the file open (a scan of /proc/*/fd)."""
    now = time.time() if now is None else now
    if tickets is None:
        tickets = load_tickets(cfg.ticket_dir, gpu=gpu.uuid, cleanup=cleanup, log_path=cfg.log_path)
    tickets = [t for t in tickets if t.get("gpu") == gpu.uuid]
    locks = path_locks(gpu.paths, read_proc_locks() if records is None else records)
    by_pid: dict[int, dict[str, Any]] = {}
    for ticket in tickets:
        for key in ("pid", "waiter_pid", "child_pid"):
            if ticket.get(key):
                by_pid[int(ticket[key])] = ticket
    holders: list[dict[str, Any]] = []
    seen: set[int] = set()
    dead_paths: list[str] = []
    for path in gpu.paths:
        for record in locks[path].holders:
            if record.pid is None or record.pid in seen:
                continue
            seen.add(record.pid)
            alive = proc_stat(record.pid) is not None
            if not alive:
                dead_paths.append(path)
                continue
            ticket = by_pid.get(record.pid)
            holders.append({
                "pid": record.pid,
                "kind": ("claim" if ticket and ticket.get("kind") == "queue" else
                         "helper" if ticket else "raw"),
                "mode": "exclusive" if record.mode == "WRITE" else "shared",
                "command": proc_cmdline(record.pid, 160),
                "ticket": _ticket_brief(ticket) if ticket else None,
                "path": path,
            })
    if dead_paths and not scan_fds:
        if not any(t.get("state") in HOLDING_STATES for t in tickets):
            holders.append({"pid": None, "kind": "unknown", "mode": "exclusive",
                            "command": "lock held; the process that took it has exited", "ticket": None,
                            "path": dead_paths[0]})
    elif dead_paths:
        # the process that took the lock exited; whoever has the file open keeps it
        identities = [locks[p].identity for p in dead_paths if locks[p].identity is not None]
        waiting = {w.pid for p in gpu.paths for w in locks[p].waiters if w.pid}
        # members that merely have the file open while waiting do not hold it
        not_holding = {int(t[k]) for t in tickets if t.get("state") not in HOLDING_STATES
                       for k in ("pid", "waiter_pid") if t.get(k)}
        open_by = fd_holders(identities)
        keepers = sorted({pid for pids in open_by.values() for pid in pids} - waiting - seen - not_holding)
        table = process_table() if keepers else {}
        holding_members = {int(t["pid"]) for t in tickets if t.get("state") in HOLDING_STATES}
        covered = descendants(table, holding_members) - holding_members
        kept = 0
        for pid in keepers:
            ticket = by_pid.get(pid)
            if ticket is None and pid in covered:
                continue  # a child of a member that holds the GPU
            if ticket is None:
                kept += 1
                if kept > 3:
                    continue
            holders.append({
                "pid": pid,
                "kind": ("claim" if ticket and ticket.get("kind") == "queue" else
                         "helper" if ticket else "kept"),
                "mode": "exclusive",
                "command": proc_cmdline(pid, 160),
                "ticket": _ticket_brief(ticket) if ticket else None,
                "path": dead_paths[0],
            })
        if kept > 3:
            holders.append({"pid": None, "kind": "unknown", "mode": "exclusive",
                            "command": f"and {kept - 3} more process(es) with the file open", "ticket": None,
                            "path": dead_paths[0]})
        if not keepers:
            holders.append({"pid": None, "kind": "unknown", "mode": "exclusive",
                            "command": "lock held; the process that took it has exited", "ticket": None,
                            "path": dead_paths[0]})
    for ticket in tickets:  # members that hold the GPU but took the lock in a child that exited
        if ticket.get("state") in ("held", "running", "handing-off") and int(ticket["pid"]) not in seen:
            if not any(h.get("ticket") and h["ticket"]["pid"] == ticket["pid"] for h in holders):
                holders = [h for h in holders if h["kind"] not in ("unknown", "kept") or h["pid"] in by_pid]
                holders.append({"pid": int(ticket["pid"]),
                                "kind": "claim" if ticket.get("kind") == "queue" else "helper",
                                "mode": "exclusive", "command": ticket.get("command") or ticket.get("label") or "",
                                "ticket": _ticket_brief(ticket), "path": gpu.paths[0]})
    # one command usually holds both files (``flock L1 flock L2 cmd``): show the outermost process only
    holder_pids = {h["pid"] for h in holders if h.get("pid")}
    if len(holder_pids) > 1:
        def ancestors(pid: int) -> set[int]:
            found: set[int] = set()
            for _ in range(32):
                stat = proc_stat(pid)
                if stat is None or stat.ppid <= 1 or stat.ppid in found:
                    break
                found.add(stat.ppid)
                pid = stat.ppid
            return found

        inner = {pid for pid in holder_pids if ancestors(pid) & (holder_pids - {pid})}
        holders = [h for h in holders if h.get("pid") not in inner]
    members = member_pids(tickets)
    raw_waiters = []
    for path in gpu.paths:
        for record in locks[path].waiters:
            if record.pid is None or record.pid in members:
                continue
            raw_waiters.append({"pid": record.pid, "path": path, "command": proc_cmdline(record.pid, 160)})
    line: list[dict[str, Any]] = []
    queue = queue_entries(queue_state, gpu.uuid, now)
    claim = next((t for t in tickets if t.get("kind") == "queue"), None)
    for entry in queue:
        line.append({"kind": "job", "job_id": entry.get("job_id"), "name": entry.get("name"),
                     "partition": entry.get("partition"), "since_unix": float(entry["since_unix"]),
                     "ready_since_unix": entry.get("ready_since_unix"), "state": entry.get("state"),
                     "claim_state": claim.get("state") if claim else None,
                     "key": (float(entry["since_unix"]), 0, int(entry.get("job_id") or 0))})
    for ticket in tickets:
        if ticket.get("kind") == "queue":
            if not queue and ticket.get("state") not in HOLDING_STATES:
                line.append({"kind": "claim", "pid": ticket.get("pid"), "since_unix": ticket.get("since_unix"),
                             "state": ticket.get("state"), "label": ticket.get("label"), "key": ticket_key(ticket)})
            continue
        if ticket.get("state") in ("held", "running", "handing-off"):
            continue
        line.append({"kind": "helper", "pid": ticket.get("pid"), "since_unix": ticket.get("since_unix"),
                     "state": ticket.get("state"), "label": ticket.get("label"),
                     "command": ticket.get("command"), "note": ticket.get("note"), "key": ticket_key(ticket)})
    line.sort(key=lambda item: item["key"])
    for item in line:
        item.pop("key", None)
    return {
        "gpu": gpu.uuid,
        "index": gpu.index,
        "name": gpu.name,
        "lock_order": list(gpu.paths),
        "holders": holders,
        "line": line,
        "raw_waiters": raw_waiters,
        "queue_state_fresh": queue_state_fresh(queue_state, now),
        "queue_state_age_s": (None if queue_state is None
                              else round(now - float(queue_state.get("heartbeat_unix") or 0.0), 1)),
    }


def _ticket_brief(ticket: dict[str, Any] | None) -> dict[str, Any] | None:
    if ticket is None:
        return None
    return {k: ticket.get(k) for k in ("kind", "pid", "state", "since_unix", "label", "command", "held_since_unix",
                                       "claim_unit", "note", "cores")}


def holder_text(holder: dict[str, Any], now: float | None = None) -> str:
    now = time.time() if now is None else now
    kind = holder.get("kind")
    ticket = holder.get("ticket") or {}
    if kind == "job":
        return f"queue job {holder.get('job_id')} ({holder.get('mode')})"
    if kind == "claim":
        since = ticket.get("held_since_unix")
        extra = f" since {now_iso(since)} ({duration_text(now - since)})" if since else ""
        return f"SimQueue GPU claim pid {holder.get('pid')} (held for the queue{extra})"
    if kind == "helper":
        since = ticket.get("held_since_unix")
        extra = f", held {duration_text(now - since)}" if since else ""
        return f"gpu-lease helper pid {holder.get('pid')}{extra}: {ticket.get('command') or holder.get('command')}"
    if kind == "kept":
        return f"raw lease kept by pid {holder.get('pid')} (the locker exited): {holder.get('command')}"
    if kind == "unknown":
        return f"raw lease: {holder.get('command')}"
    return f"raw flock holder pid {holder.get('pid')}: {holder.get('command')}"


def render_snapshot(snapshot: dict[str, Any], *, indent: str = "    ", now: float | None = None,
                    warn_s: float | None = None, queue_jobs: list[dict[str, Any]] | None = None) -> list[str]:
    """Text lines for ``simq partitions`` and ``gpu-lease --status``."""
    now = time.time() if now is None else now
    lines: list[str] = []
    holders = snapshot.get("holders") or []
    if not holders:
        lines.append(f"{indent}holder: none (free)")
    for holder in holders:
        lines.append(f"{indent}holder: {holder_text(holder, now)}")
    line = snapshot.get("line") or []
    if line:
        lines.append(f"{indent}waiting (in order):")
    for position, item in enumerate(line, start=1):
        if item["kind"] == "job":
            ready = item.get("ready_since_unix") or item.get("since_unix")
            state = item.get("state") or "waiting"
            claim = f", claim {item['claim_state']}" if item.get("claim_state") else ""
            lines.append(f"{indent}  {position}. queue job {item.get('job_id')} ({item.get('partition') or '-'}, "
                         f"{item.get('name') or '-'}): {state}{claim}; submitted {now_iso(item['since_unix'])}, "
                         f"waiting for the GPU {duration_text(now - float(ready))}")
        elif item["kind"] == "claim":
            lines.append(f"{indent}  {position}. SimQueue GPU claim pid {item.get('pid')} ({item.get('state')}), "
                         f"since {now_iso(float(item.get('since_unix') or 0))}")
        else:
            note = f" [{item['note']}]" if item.get("note") else ""
            lines.append(f"{indent}  {position}. gpu-lease pid {item.get('pid')} ({item.get('state')}) since "
                         f"{now_iso(float(item.get('since_unix') or 0))} "
                         f"({duration_text(now - float(item.get('since_unix') or now))}): "
                         f"{item.get('command') or item.get('label')}{note}")
    raw = snapshot.get("raw_waiters") or []
    if raw:
        lines.append(f"{indent}raw flock waiters (served after the line):")
        for item in raw:
            lines.append(f"{indent}  pid {item['pid']} on {os.path.basename(item['path'])}: {item['command']}")
    if not snapshot.get("queue_state_fresh"):
        age = snapshot.get("queue_state_age_s")
        lines.append(f"{indent}(queue state " + ("missing" if age is None else f"{age:.0f} s old: worker down?")
                     + "; queue jobs not shown)")
    if warn_s is not None:
        for item in line:
            if item["kind"] != "job":
                continue
            ready = float(item.get("ready_since_unix") or item.get("since_unix") or now)
            if now - ready > warn_s:
                lines.append(f"{indent}WARNING: queue job {item.get('job_id')} has waited "
                             f"{duration_text(now - ready)} for the GPU (threshold {duration_text(warn_s)})")
    return lines


# --------------------------------------------------------------------------- claim units (the queue's member)


def _peer_uid(connection: socket.socket) -> int | None:
    try:
        data = connection.getsockopt(socket.SOL_SOCKET, socket.SO_PEERCRED, struct.calcsize("3i"))
    except OSError:
        return None
    _pid, uid, _gid = struct.unpack("3i", data)
    return uid


class Claim:
    """Main process of ``<prefix>-gpuclaim-<token>.service``."""

    def __init__(self, spec_path: Path) -> None:
        self.spec_path = spec_path
        self.spec = json.loads(spec_path.read_text(encoding="utf-8"))
        if self.spec.get("schema") != CLAIM_SPEC_SCHEMA:
            raise RuntimeError(f"unsupported claim spec schema {self.spec.get('schema')!r}")
        self.cfg = LineConfig.from_dict(self.spec["line"])
        self.gpu = GpuLocks.from_dict(self.spec["gpu"])
        self.status_path = Path(self.spec["status_path"])
        self.control_path = Path(self.spec["control_path"])
        self.stop_signal: int | None = None
        self._control_mtime: float | None = None
        self.control: dict[str, Any] = {}
        self.status: dict[str, Any] = {
            "schema": CLAIM_STATUS_SCHEMA,
            "token": self.spec["token"],
            "unit": self.spec["unit"],
            "gpu": self.gpu.uuid,
            "pid": os.getpid(),
            "pid_start": pid_start(os.getpid()),
            "phase": "starting",
            "started_unix": time.time(),
            "boot_id": _boot_id(),
            "invocation_id": os.environ.get("INVOCATION_ID"),
        }

    def save(self, **fields: Any) -> None:
        self.status.update(fields)
        self.status["updated_unix"] = time.time()
        write_json_atomic(self.status_path, self.status, mode=0o600)

    def print(self, text: str) -> None:
        sys.stdout.write(f"{now_iso()} {text}\n")
        sys.stdout.flush()

    def _handle(self, signum: int, _frame: Any) -> None:
        self.stop_signal = signum

    def read_control(self) -> dict[str, Any]:
        try:
            mtime = self.control_path.stat().st_mtime
        except OSError:
            return self.control
        if mtime != self._control_mtime:
            data = read_json(self.control_path)
            if data is not None and data.get("schema") == CLAIM_CONTROL_SCHEMA:
                self.control = data
                self._control_mtime = mtime
        return self.control

    def since(self) -> float | None:
        value = self.read_control().get("since_unix")
        return None if value is None else float(value)

    def should_stop(self) -> bool:
        return self.stop_signal is not None or bool(self.read_control().get("release"))

    def run(self) -> int:
        for number in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
            signal.signal(number, self._handle)
        self.read_control()
        listener = socket.socket(socket.AF_UNIX, socket.SOCK_SEQPACKET | socket.SOCK_CLOEXEC)
        try:
            listener.bind("\0" + self.spec["socket_name"])
        except OSError as error:
            self.save(phase="failed", error=f"cannot bind the hand-off socket: {error}", finished_unix=time.time())
            return 1
        listener.listen(8)
        since = self.since() or float(self.spec["since_unix"])
        member = LineMember(
            self.cfg, self.gpu, kind="queue", since_unix=since,
            label=f"SimQueue claim {self.spec['unit']}",
            command=f"SimQueue GPU claim for the queue ({self.spec['unit']})",
            extra={"claim_unit": self.spec["unit"], "claim_token": self.spec["token"]},
            should_stop=self.should_stop, since_getter=self.since,
            on_state=lambda state, fields: self.save(phase=state, **{k: v for k, v in fields.items()
                                                                     if k in ("waiter_pid", "note",
                                                                              "held_since_unix")}),
            log=self.print,
        )
        self.save(phase="waiting", since_unix=since)
        self.print(f"claim for {self.gpu.uuid} ({', '.join(self.gpu.paths)}), line since {now_iso(since)}")
        try:
            try:
                fds = member.acquire()
            except Stopped:
                self.print("stopped before the GPU was won")
                self.save(phase="released", reason="stopped before the GPU was won", finished_unix=time.time(),
                          requeued=member.requeue_count)
                return 0
            held_at = time.time()
            self.print(f"holding the GPU locks (requeued {member.requeue_count} raw waiter(s))")
            self.save(phase="held", held_since_unix=held_at, requeued=member.requeue_count)
            return self.serve(listener, member, fds, held_at)
        finally:
            member.close()
            listener.close()

    def serve(self, listener: socket.socket, member: LineMember, fds: list[int], held_at: float) -> int:
        hold_max = float(self.spec.get("hold_max_s") or self.cfg.claim_hold_max_s)
        while True:
            if self.should_stop():
                reason = "stopped by the worker" if self.stop_signal else "released on request"
                self.print(f"releasing the GPU locks: {reason}")
                self.save(phase="released", reason=reason, finished_unix=time.time())
                return 0
            if time.time() - held_at > hold_max:
                reason = f"no hand-off within {int(hold_max)} s"
                self.print(f"releasing the GPU locks: {reason}")
                self.save(phase="released", reason=reason, finished_unix=time.time())
                return 0
            ready, _, _ = select.select([listener], [], [], 0.2)
            if not ready:
                continue
            try:
                connection, _ = listener.accept()
            except OSError:
                continue
            with connection:
                result = self.handle(connection, member, fds)
            if result is not None:
                self.print(f"handed the GPU locks to job {result.get('job_id')} ({result.get('unit')})")
                member.fds = []  # the descriptors are closed below; the job keeps the locks
                for descriptor in fds:
                    os.close(descriptor)
                self.save(phase="handed-off", handed_to=result, finished_unix=time.time())
                return 0

    def handle(self, connection: socket.socket, member: LineMember, fds: list[int]) -> dict[str, Any] | None:
        connection.settimeout(30)
        if _peer_uid(connection) != os.getuid():
            return None
        try:
            request = json.loads(connection.recv(65536).decode("utf-8"))
        except (OSError, ValueError):
            return None
        if request.get("protocol") != HANDOFF_PROTOCOL or request.get("secret") != self.spec["secret"]:
            try:
                connection.send(json.dumps({"ok": False, "error": "refused"}).encode())
            except OSError:
                pass
            return None
        member.set_state("handing-off", note=f"handing off to job {request.get('job_id')}")
        header = {"ok": True, "protocol": HANDOFF_PROTOCOL, "paths": list(self.gpu.paths),
                  "identities": [list(file_identity(p) or ("", 0)) for p in self.gpu.paths]}
        try:
            socket.send_fds(connection, [json.dumps(header).encode()], fds)
            reply = json.loads(connection.recv(65536).decode("utf-8") or "{}")
        except (OSError, ValueError) as error:
            reply = {"taken": False, "reason": f"no reply: {error}"}
        if reply.get("taken"):
            return {"job_id": request.get("job_id"), "unit": request.get("unit"),
                    "launch_token": request.get("launch_token"), "at": now_iso()}
        self.print(f"hand-off to job {request.get('job_id')} not completed: {reply.get('reason')}; still holding")
        member.set_state("held", note=f"hand-off to job {request.get('job_id')} not completed")
        self.save(last_failed_handoff={"job_id": request.get("job_id"), "reason": reply.get("reason"),
                                       "at": now_iso()})
        return None


def _boot_id() -> str | None:
    try:
        return Path("/proc/sys/kernel/random/boot_id").read_text(encoding="utf-8").strip()
    except OSError:
        return None


def claim_main(spec_path: str) -> int:
    os.umask(0o077)
    try:
        claim = Claim(Path(spec_path))
    except Exception as error:  # noqa: BLE001 - nothing was taken
        print(f"gpu claim: cannot read spec {spec_path}: {error}", file=sys.stderr)
        return 1
    return claim.run()


# --------------------------------------------------------------------------- gpu-lease (direct GPU use)


LEASE_USAGE = """\
gpu-lease [--core any|both|N|none] [--gpu UUID|INDEX] [--timeout S] [--label TEXT] [-q] [--] COMMAND [ARGS...]
gpu-lease --status [--json]

Run COMMAND holding the GPU lease files (and optionally a core lease), in a
fair line: queued SimQueue GPU jobs that were waiting before this helper
started go first, other gpu-lease users are served in start order, and
processes blocked in a raw `flock testing-gpu.lock` are served after the line.
Locks are taken GPU first (testing-gpu.lock, then testing-gpu-<UUID>.lock),
then the core lease; everything is released when COMMAND and its children
have closed them.
"""


def _lease_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="gpu-lease", usage=LEASE_USAGE, add_help=True)
    parser.add_argument("--core", default="none",
                        help="any: whichever helper core (6 or 7) is free first, as any-core.sh --wait; "
                             "both: cores 6 and 7, as both-cores.sh; N: that core; none (default)")
    parser.add_argument("--gpu", default=None, help="GPU UUID or index (default: the node's GPU)")
    parser.add_argument("--timeout", type=float, default=None, help="give up after S seconds of waiting (exit 75)")
    parser.add_argument("--label", default=None, help="text shown in simq partitions (default: the command)")
    parser.add_argument("--nice", type=int, default=None, help="nice level for COMMAND (default 15 with --core)")
    parser.add_argument("--config", default=None, help="node configuration (default: the live node.json)")
    parser.add_argument("--no-requeue-raw", action="store_true",
                        help="do not move raw flock waiters behind this request")
    parser.add_argument("--status", action="store_true", help="show the holders and the line of each GPU")
    parser.add_argument("--json", action="store_true", help="with --status: JSON")
    parser.add_argument("-q", "--quiet", action="store_true", help="no progress messages on stderr")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    return parser


def _core_choice(text: str, cfg: LineConfig) -> tuple[str, list[int]]:
    value = (text or "none").strip().lower()
    if value in ("none", "no", ""):
        return "none", []
    if value == "any":
        return "any", list(cfg.helper_cores)
    if value == "both":
        return "both", sorted(cfg.helper_cores)
    if value.isdigit():
        return "one", [int(value)]
    raise ValueError(f"--core must be any, both, none or a CPU number (got {text!r})")


def status_main(cfg: LineConfig, *, as_json: bool) -> int:
    now = time.time()
    state = load_queue_state(cfg.state_path)
    snapshots = [gpu_snapshot(cfg, gpu, queue_state=state, now=now) for gpu in cfg.gpus.values()]
    if as_json:
        print(json.dumps({"gpus": snapshots, "state_path": cfg.state_path, "ticket_dir": cfg.ticket_dir},
                         indent=2, sort_keys=True, default=str))
        return 0
    for snapshot in snapshots:
        print(f"GPU {snapshot['index'] if snapshot['index'] is not None else '-'} {snapshot['gpu']} "
              f"{snapshot['name']}  (locks: {', '.join(os.path.basename(p) for p in snapshot['lock_order'])})")
        jobs = ((state or {}).get("gpus") or {}).get(snapshot["gpu"], {}).get("jobs") or []
        for job in jobs:
            print(f"    queue job {job.get('job_id')} runs on it ({job.get('mode')})")
        for line in render_snapshot(snapshot, now=now, warn_s=cfg.starvation_warn_s):
            print(line)
    return 0


def lease_main(argv: list[str]) -> int:
    parser = _lease_parser()
    args = parser.parse_args(argv)
    try:
        cfg, source = load_line_config(args.config)
    except ValueError as error:
        print(f"gpu-lease: {error}", file=sys.stderr)
        return EXIT_USAGE
    if args.status:
        return status_main(cfg, as_json=args.json)
    command = list(args.command)
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        parser.print_usage(sys.stderr)
        print("gpu-lease: no command given", file=sys.stderr)
        return EXIT_USAGE
    if not cfg.enabled:
        print(f"gpu-lease: the GPU line is disabled in {source}", file=sys.stderr)
        return EXIT_USAGE
    try:
        gpu = cfg.gpu(args.gpu)
        mode, cores = _core_choice(args.core, cfg)
    except ValueError as error:
        print(f"gpu-lease: {error}", file=sys.stderr)
        return EXIT_USAGE
    for core in cores:
        if core not in cfg.core_paths:
            print(f"gpu-lease: no lease file is configured for core {core}", file=sys.stderr)
            return EXIT_USAGE

    def say(text: str) -> None:
        if not args.quiet:
            print(f"gpu-lease: {text}", file=sys.stderr, flush=True)

    stop: dict[str, int | None] = {"signal": None}

    def on_signal(signum: int, _frame: Any) -> None:
        stop["signal"] = signum

    previous = {number: signal.signal(number, on_signal) for number in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)}
    since = time.time()
    label = args.label or " ".join(command)
    member = LineMember(
        cfg, gpu, kind="helper", since_unix=since, label=label[:160], command=" ".join(command)[:300],
        extra={"cores_requested": args.core, "cwd": os.getcwd()},
        should_stop=lambda: stop["signal"] is not None,
        deadline=None if args.timeout is None else since + max(0.0, args.timeout),
        requeue=False if args.no_requeue_raw else None,
    )
    core_fds: list[int] = []
    held_cores: list[int] = []
    try:
        while True:
            reasons = member.blockers()
            if reasons:
                say("waiting behind " + "; ".join(reasons[:3]))
            fds = member.acquire()
            say(f"holding {', '.join(os.path.basename(p) for p in gpu.paths)}")
            if mode == "none":
                break
            got = _acquire_cores(member, gpu, mode, cores, cfg, say)
            if got is None:  # lock-order conflict: release the GPU and line up again
                member.release()
                member.set_state("waiting", note="backed off: a core holder waits for the GPU")
                time.sleep(1.0 + random.random())
                continue
            core_fds, held_cores = got
            break
    except TimedOut:
        say(f"gave up after {args.timeout:g} s")
        member.close()
        return EXIT_TIMEOUT
    except Stopped:
        member.close()
        return 128 + int(stop["signal"] or signal.SIGTERM)
    env = dict(os.environ)
    env.update({"GPU_LEASE_GPU": gpu.uuid, "GPU_LEASE_LOCKS": ":".join(gpu.paths)})
    nice = args.nice
    if held_cores:
        env.update({"PF_CORE": ",".join(str(c) for c in held_cores), "GPU_LEASE_CORES": ",".join(map(str, held_cores)),
                    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
        if len(held_cores) > 1:
            env.update({"MAKEFLAGS": f"-j{len(held_cores)}", "CMAKE_BUILD_PARALLEL_LEVEL": str(len(held_cores))})
        if nice is None:
            nice = 15

    def prepare() -> None:
        if held_cores:
            os.sched_setaffinity(0, set(held_cores))
        if nice is not None:
            current = os.nice(0)
            if nice > current:
                os.nice(nice - current)

    inherited = list(member.fds) + list(core_fds)
    try:
        child = subprocess.Popen(command, env=env, pass_fds=inherited, preexec_fn=prepare)
    except OSError as error:
        print(f"gpu-lease: cannot run {command[0]}: {error}", file=sys.stderr)
        member.close()
        for descriptor in core_fds:
            os.close(descriptor)
        return 127
    member.set_state("running", child_pid=child.pid, cores=held_cores, started_unix=time.time(), note=None)
    forwarded = {signal.SIGTERM, signal.SIGHUP}

    def forward(signum: int, _frame: Any) -> None:
        if signum in forwarded:
            try:
                child.send_signal(signum)
            except OSError:
                pass

    for number in (signal.SIGTERM, signal.SIGHUP, signal.SIGINT):
        signal.signal(number, forward)
    while True:
        try:
            returncode = child.wait()
            break
        except InterruptedError:
            continue
    for number, handler in previous.items():
        signal.signal(number, handler)
    member.close()
    for descriptor in core_fds:
        os.close(descriptor)
    if returncode < 0:
        return 128 - returncode
    return returncode


def _acquire_cores(member: LineMember, gpu: GpuLocks, mode: str, cores: list[int], cfg: LineConfig,
                   say: Callable[[str], None]) -> tuple[list[int], list[int]] | None:
    """Core leases after the GPU (the documented order).  ``None`` when a
    core's holder waits for the GPU (lock-order conflict)."""
    holding = list(gpu.paths)
    if mode in ("both", "one"):
        fds: list[int] = []
        for core in cores:
            path = cfg.core_paths[core]
            descriptor = member.acquire_nonblocking(path, holding, note=f"holding the GPU; waiting for core {core}")
            if descriptor is None:
                for item in fds:
                    os.close(item)
                return None
            fds.append(descriptor)
            holding.append(path)
        say(f"holding core lease(s) {', '.join(map(str, cores))}")
        return fds, list(cores)
    # any: the first of the helper cores that is free
    handles = {core: os.open(cfg.core_paths[core], os.O_RDONLY | os.O_CREAT | os.O_CLOEXEC, 0o644) for core in cores}
    noted = False
    last_check = 0.0
    try:
        while True:
            for core, descriptor in handles.items():
                try:
                    fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    continue
                for other, handle in handles.items():
                    if other != core:
                        os.close(handle)
                say(f"holding core lease {core}")
                handles = {}
                return [descriptor], [core]
            member._check()
            if not noted:
                member.set_state(str(member.ticket.get("state")),
                                 note=f"holding the GPU; waiting for one of cores {', '.join(map(str, cores))}")
                say(f"waiting for one of cores {', '.join(map(str, cores))}")
                noted = True
            if time.monotonic() - last_check > 0.5:
                last_check = time.monotonic()
                conflicts = [lock_order_conflict(cfg.core_paths[core], holding) for core in cores]
                if all(conflicts):
                    member.log("every candidate core is held by a process waiting for the GPU: " + conflicts[0])
                    return None
            time.sleep(0.2)
    finally:
        for handle in handles.values():
            os.close(handle)


# --------------------------------------------------------------------------- entry point


def main(argv: list[str]) -> int:
    if len(argv) >= 3 and argv[1] == "claim":
        return claim_main(argv[2])
    if len(argv) >= 2 and argv[1] == "lease":
        return lease_main(argv[2:])
    if len(argv) >= 2 and argv[1] == "status":
        try:
            cfg, _source = load_line_config(None)
        except ValueError as error:
            print(error, file=sys.stderr)
            return EXIT_USAGE
        return status_main(cfg, as_json="--json" in argv[2:])
    print("usage: gpuline.py claim SPEC.json | lease [options] -- COMMAND | status [--json]", file=sys.stderr)
    return EXIT_USAGE


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
