"""Interoperation with the ad-hoc flock(2) lease files used outside the queue.

The scheduler never takes these locks itself.  Each job unit's wrapper takes
the lease files of its allocated cores and GPUs for the job's lifetime, and
the scheduler only *observes* them through ``/proc/locks`` so that a probe can
never make a manual ``flock -n`` fail.
"""

from __future__ import annotations

import fcntl
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .nodeconfig import NodeConfig


PROC_LOCKS = Path("/proc/locks")


@dataclass(frozen=True)
class LockHolder:
    path: str
    kind: str  # FLOCK / POSIX / OFDLCK
    mode: str  # WRITE / READ
    pid: int | None


@dataclass(frozen=True)
class LockPlan:
    path: str
    mode: str  # "exclusive" | "shared"
    resource: str  # "cpu:6" / "gpu:GPU-..."


def lock_plan_for(
    config: NodeConfig,
    cpus: Iterable[int],
    gpus: Iterable[str],
    *,
    gpu_shared: bool,
) -> list[LockPlan]:
    """Lease files a job unit must hold, in a global ascending order."""
    plans: dict[str, LockPlan] = {}
    for cpu in cpus:
        for path in config.locks.core_lock_paths(cpu):
            plans[path] = LockPlan(path, "exclusive", f"cpu:{cpu}")
    for uuid in gpus:
        gpu = config.gpus.get(uuid)
        if gpu is None:
            continue
        for path in config.locks.gpu_lock_paths(gpu):
            mode = "shared" if gpu_shared else "exclusive"
            previous = plans.get(path)
            if previous is not None and previous.mode == "exclusive":
                continue
            plans[path] = LockPlan(path, mode, f"gpu:{uuid}")
    return [plans[key] for key in sorted(plans)]


def _device_key(path: str) -> tuple[str, int] | None:
    try:
        stat = os.stat(path)
    except OSError:
        return None
    return f"{os.major(stat.st_dev):02x}:{os.minor(stat.st_dev):02x}", stat.st_ino


def read_proc_locks(text: str | None = None) -> list[tuple[str, str, int | None, str, int]]:
    """Parse ``/proc/locks`` into (kind, mode, pid, device, inode)."""
    if text is None:
        try:
            text = PROC_LOCKS.read_text(encoding="ascii", errors="replace")
        except OSError:
            return []
    records: list[tuple[str, str, int | None, str, int]] = []
    for line in text.splitlines():
        fields = line.split()
        if "->" in fields:
            continue  # a blocked waiter, not a holder
        if len(fields) < 6:
            continue
        kind, mode = fields[1], fields[3]
        try:
            pid: int | None = int(fields[4])
        except ValueError:
            pid = None
        if pid is not None and pid <= 0:
            pid = None
        parts = fields[5].split(":")
        if len(parts) != 3:
            continue
        try:
            inode = int(parts[2])
        except ValueError:
            continue
        records.append((kind, mode, pid, f"{parts[0]}:{parts[1]}".lower(), inode))
    return records


def lock_holders(paths: Iterable[str], *, proc_locks_text: str | None = None) -> dict[str, list[LockHolder]]:
    """Map each lease path to the locks currently held on it."""
    keys: dict[tuple[str, int], list[str]] = {}
    for path in paths:
        key = _device_key(path)
        if key is not None:
            keys.setdefault(key, []).append(path)
    result: dict[str, list[LockHolder]] = {}
    if not keys:
        return result
    for kind, mode, pid, device, inode in read_proc_locks(proc_locks_text):
        for path in keys.get((device, inode), ()):
            result.setdefault(path, []).append(LockHolder(path, kind, mode, pid))
    return result


def probe_free(path: str) -> bool:
    """Fallback probe used only when /proc/locks is unreadable."""
    try:
        handle = open(path, "a")
    except OSError:
        return False
    try:
        fcntl.flock(handle.fileno(), fcntl.LOCK_SH | fcntl.LOCK_NB)
    except BlockingIOError:
        handle.close()
        return False
    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    handle.close()
    return True


def busy_resources(
    config: NodeConfig,
    *,
    allocated_cpus: set[int],
    allocated_gpus: set[str],
    proc_locks_text: str | None = None,
) -> tuple[dict[int, list[LockHolder]], dict[str, list[LockHolder]]]:
    """Resources whose lease files are held although no queue job owns them.

    A lease held on a resource that an active queue job owns is that job's own
    lease (or a lingering descendant of it) and is already accounted for.
    """
    if not config.locks.enabled:
        return {}, {}
    cpu_paths: dict[str, list[int]] = {}
    for cpu in config.node_cpus():
        for path in config.locks.core_lock_paths(cpu):
            cpu_paths.setdefault(path, []).append(cpu)
    gpu_paths: dict[str, list[str]] = {}
    for uuid, gpu in config.gpus.items():
        for path in config.locks.gpu_lock_paths(gpu):
            gpu_paths.setdefault(path, []).append(uuid)
    all_paths = list(cpu_paths) + [p for p in gpu_paths if p not in cpu_paths]
    use_proc = proc_locks_text is not None or PROC_LOCKS.exists()
    holders = lock_holders(all_paths, proc_locks_text=proc_locks_text) if use_proc else {}
    if not use_proc:
        for path in all_paths:
            if Path(path).exists() and not probe_free(path):
                holders[path] = [LockHolder(path, "FLOCK", "WRITE", None)]
    busy_cpus: dict[int, list[LockHolder]] = {}
    busy_gpus: dict[str, list[LockHolder]] = {}
    for path, found in holders.items():
        for cpu in cpu_paths.get(path, ()):
            if cpu not in allocated_cpus:
                busy_cpus.setdefault(cpu, []).extend(found)
        for uuid in gpu_paths.get(path, ()):
            if uuid not in allocated_gpus:
                busy_gpus.setdefault(uuid, []).extend(found)
    return busy_cpus, busy_gpus


def process_description(pid: int | None) -> str:
    if pid is None:
        return "unknown process"
    try:
        raw = Path(f"/proc/{pid}/cmdline").read_bytes()
    except OSError:
        return f"pid {pid} (exited; lease kept by a descendant)"
    command = " ".join(part.decode("utf-8", "replace") for part in raw.split(b"\0") if part)
    return f"pid {pid} {command[:120]}"
