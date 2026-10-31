"""CPU usage of lease holders, owner jobs and logical CPUs (SimQueue 0.6.1).

The node worker samples ``/proc`` once per pass while borrow jobs exist (a
full pass reads every ``/proc/<pid>/stat`` and, when a lease is held outside
the queue, the ``fd`` links of this account's processes: a few milliseconds).

* A *holder tree* is every process that has a core's lease file open (found
  through ``/proc/<pid>/fd``, because the pid in ``/proc/locks`` is the
  process that took the lock and may have exited, leaving the lock to its
  descendants), plus all their descendants.  For an owner queue job it is
  the job unit's main process (the wrapper) and its descendants.
* Usage of a tree is accumulated from per-process deltas of
  ``utime + stime + cutime + cstime`` (so children reaped inside the tree are
  counted), in cores averaged over a window.
* Per-CPU busy fractions come from ``/proc/stat``.

:func:`classify_external` and :func:`classify_owner` turn usage figures into
:class:`borrow.CoreState` flags; they are pure and unit tested with fake
numbers.
"""

from __future__ import annotations

import os
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping

from .borrow import in_use_cores
from .nodeconfig import BorrowConfig


try:
    CLK_TCK = os.sysconf("SC_CLK_TCK")
except (AttributeError, ValueError, OSError):  # pragma: no cover
    CLK_TCK = 100
PROC = Path("/proc")


@dataclass(frozen=True)
class ProcInfo:
    pid: int
    ppid: int
    start: int  # start time in clock ticks since boot (identity across samples)
    ticks: int  # utime + stime + cutime + cstime
    state: str


def parse_stat(pid: int, text: str) -> ProcInfo | None:
    try:
        rest = text[text.rfind(")") + 2 :].split()
        return ProcInfo(
            pid=pid,
            ppid=int(rest[1]),
            start=int(rest[19]),
            ticks=int(rest[11]) + int(rest[12]) + int(rest[13]) + int(rest[14]),
            state=rest[0],
        )
    except (IndexError, ValueError):
        return None


def read_process_table(proc: Path = PROC) -> dict[int, ProcInfo]:
    table: dict[int, ProcInfo] = {}
    try:
        names = os.listdir(proc)
    except OSError:
        return table
    for name in names:
        if not name.isdigit():
            continue
        try:
            with open(proc / name / "stat", "r", encoding="utf-8", errors="replace") as handle:
                text = handle.read()
        except OSError:
            continue
        info = parse_stat(int(name), text)
        if info is not None:
            table[info.pid] = info
    return table


def descendants(table: Mapping[int, ProcInfo], roots: Iterable[int]) -> set[int]:
    children: dict[int, list[int]] = {}
    for info in table.values():
        children.setdefault(info.ppid, []).append(info.pid)
    found: set[int] = set()
    stack = [pid for pid in roots if pid in table]
    while stack:
        pid = stack.pop()
        if pid in found:
            continue
        found.add(pid)
        stack.extend(children.get(pid, ()))
    return found


def file_key(path: str) -> tuple[int, int] | None:
    try:
        stat = os.stat(path)
    except OSError:
        return None
    return stat.st_dev, stat.st_ino


def find_fd_holders(keys: Iterable[tuple[int, int]], proc: Path = PROC) -> dict[tuple[int, int], set[int]]:
    """Processes (of this user) that have any of the files *keys* open."""
    wanted = set(keys)
    found: dict[tuple[int, int], set[int]] = {key: set() for key in wanted}
    if not wanted:
        return found
    uid = os.getuid()
    try:
        names = os.listdir(proc)
    except OSError:
        return found
    for name in names:
        if not name.isdigit():
            continue
        base = proc / name
        try:
            if os.stat(base).st_uid != uid:
                continue
            fds = os.listdir(base / "fd")
        except OSError:
            continue
        for fd in fds:
            try:
                stat = os.stat(base / "fd" / fd)
            except OSError:
                continue
            key = (stat.st_dev, stat.st_ino)
            if key in wanted:
                found[key].add(int(name))
    return found


def command_summary(pid: int, width: int = 100, proc: Path = PROC) -> str:
    try:
        raw = (proc / str(pid) / "cmdline").read_bytes()
    except OSError:
        return "(exited)"
    command = " ".join(part.decode("utf-8", "replace") for part in raw.split(b"\0") if part)
    if not command:
        try:
            command = "[" + (proc / str(pid) / "comm").read_text(encoding="utf-8").strip() + "]"
        except OSError:
            command = "?"
    return command if len(command) <= width else command[: width - 1] + "…"


def read_cpu_times(text: str | None = None) -> dict[int, tuple[int, int]]:
    """Per logical CPU: (busy ticks, total ticks) from /proc/stat."""
    if text is None:
        try:
            text = (PROC / "stat").read_text(encoding="ascii", errors="replace")
        except OSError:
            return {}
    result: dict[int, tuple[int, int]] = {}
    for line in text.splitlines():
        if not line.startswith("cpu") or line.startswith("cpu "):
            continue
        fields = line.split()
        try:
            cpu = int(fields[0][3:])
            values = [int(v) for v in fields[1:9]]
        except (ValueError, IndexError):
            continue
        total = sum(values)
        idle = values[3] + values[4]  # idle + iowait
        result[cpu] = (total - idle, total)
    return result


@dataclass
class _TreeSeries:
    identity: str
    prev: dict[int, tuple[int, int]] = field(default_factory=dict)
    cumulative: float = 0.0
    history: deque = field(default_factory=deque)  # (time, cumulative cpu seconds)
    initialized: bool = False
    last_seen: float = 0.0


class UsageTracker:
    """Rolling CPU usage of process trees and logical CPUs."""

    def __init__(self, *, hz: int = CLK_TCK, keep_s: float = 900.0) -> None:
        self.hz = hz
        self.keep_s = keep_s
        self.trees: dict[str, _TreeSeries] = {}
        self.cpu_history: dict[int, deque] = {}

    def reset(self) -> None:
        self.trees.clear()
        self.cpu_history.clear()

    def observe_tree(
        self,
        key: str,
        identity: str,
        pids: Iterable[int],
        table: Mapping[int, ProcInfo],
        now: float,
    ) -> None:
        series = self.trees.get(key)
        if series is None or series.identity != identity:
            series = _TreeSeries(identity=identity)
            self.trees[key] = series
        delta = 0
        current: dict[int, tuple[int, int]] = {}
        for pid in pids:
            info = table.get(pid)
            if info is None:
                continue
            previous = series.prev.get(pid)
            if previous is not None and previous[0] == info.start:
                delta += max(0, info.ticks - previous[1])
            elif series.initialized:
                delta += info.ticks  # joined the tree since the last sample
            current[pid] = (info.start, info.ticks)
        series.prev = current
        if series.initialized:
            series.cumulative += delta / self.hz
        series.initialized = True
        series.last_seen = now
        series.history.append((now, series.cumulative))
        while series.history and now - series.history[0][0] > self.keep_s:
            series.history.popleft()

    def observe_cpus(self, times: Mapping[int, tuple[int, int]], now: float) -> None:
        for cpu, (busy, total) in times.items():
            history = self.cpu_history.setdefault(cpu, deque())
            history.append((now, busy, total))
            while history and now - history[0][0] > self.keep_s:
                history.popleft()

    def forget_trees(self, now: float, max_age: float = 60.0) -> None:
        for key in [k for k, s in self.trees.items() if now - s.last_seen > max_age]:
            del self.trees[key]

    @staticmethod
    def _window_start(history: deque, now: float, window: float, full: bool) -> int | None:
        """Index of the newest sample at least *window* before the last one."""
        if len(history) < 2:
            return None
        last_time = history[-1][0]
        chosen = None
        for index in range(len(history) - 2, -1, -1):
            if last_time - history[index][0] >= window - 1e-6:
                chosen = index
                break
        if chosen is None:
            if full:
                return None
            chosen = 0
            if last_time - history[0][0] < min(0.2, window):
                return None
        return chosen

    def tree_usage(self, key: str, window: float, now: float, *, full: bool = True) -> float | None:
        """Average cores used by tree *key* over the last *window* seconds;
        ``None`` while fewer than *window* seconds have been observed
        (unless *full* is false)."""
        series = self.trees.get(key)
        if series is None:
            return None
        start = self._window_start(series.history, now, window, full)
        if start is None:
            return None
        then, before = series.history[start]
        last_time, last = series.history[-1]
        span = last_time - then
        if span <= 0:
            return None
        return max(0.0, (last - before) / span)

    def cpu_busy(self, cpu: int, window: float, now: float, *, full: bool = True) -> float | None:
        """Busy fraction of logical CPU *cpu* over the last *window* seconds."""
        history = self.cpu_history.get(cpu)
        if not history:
            return None
        start = self._window_start(history, now, window, full)
        if start is None:
            return None
        _, busy0, total0 = history[start]
        _, busy1, total1 = history[-1]
        if total1 <= total0:
            return None
        return min(1.0, max(0.0, (busy1 - busy0) / (total1 - total0)))


@dataclass(frozen=True)
class HolderGroup:
    """Cores whose leases are held by one process tree outside the queue."""

    key: str  # identity: the root pids
    cpus: tuple[int, ...]
    description: str
    usage_idle: float | None  # cores, averaged over the idle window
    usage_busy: float | None  # cores, averaged over the busy window
    stale: bool = False


@dataclass(frozen=True)
class CoreFlags:
    external: bool = False
    external_idle: bool = False
    external_busy: bool = False
    stale: bool = False
    holder: str = ""
    holder_key: str = ""
    owner_idle: bool = False
    owner_busy: bool = False


def _pick_idle_cpus(
    cpus: Iterable[int],
    count: int,
    cpu_busy: Mapping[int, float | None],
    occupied: set[int],
    threshold: float,
) -> list[int]:
    if count <= 0:
        return []
    candidates = []
    for cpu in cpus:
        if cpu in occupied:
            continue
        busy = cpu_busy.get(cpu)
        if busy is None or busy > threshold:
            continue
        candidates.append((busy, -cpu, cpu))
    candidates.sort()
    return [cpu for _, _, cpu in candidates[:count]]


def classify_external(
    group: HolderGroup,
    *,
    cpu_busy: Mapping[int, float | None],
    running_borrowed: set[int],
    occupied: set[int],
    borrow: BorrowConfig,
) -> dict[int, CoreFlags]:
    """Lease-held cores of one holder tree.

    *running_borrowed*: cores of the group that running (not frozen)
    borrowers use; *occupied*: cores with any borrower.  The tree is busy
    when it uses more than the cores left to it plus ``busy_threshold``;
    idle cores are those beyond its usage (``in_use_cores``) whose CPU is
    also measured idle.
    """
    flags: dict[int, CoreFlags] = {}
    if group.stale:
        for cpu in group.cpus:
            flags[cpu] = CoreFlags(external=True, stale=True, holder=group.description, holder_key=group.key)
        return flags
    busy, idle_cpus = _idle_and_busy(group.cpus, group.usage_idle, group.usage_busy, cpu_busy,
                                     running_borrowed, occupied, borrow)
    for cpu in group.cpus:
        flags[cpu] = CoreFlags(
            external=True,
            external_idle=cpu in idle_cpus,
            external_busy=busy,
            holder=group.description,
            holder_key=group.key,
        )
    return flags


def _idle_and_busy(
    cpus: tuple[int, ...],
    usage_idle: float | None,
    usage_busy: float | None,
    cpu_busy: Mapping[int, float | None],
    running_borrowed: set[int],
    occupied: set[int],
    borrow: BorrowConfig,
) -> tuple[bool, set[int]]:
    """(busy, idle cores) of a tree holding *cpus*.

    Busy: over the busy window it uses more than the cores not lent to
    running borrowers, plus ``busy_threshold``.  Idle cores: the tree's
    usage over the idle window (rounded with ``idle_threshold``) leaves
    room; cores already carrying a borrower stay idle while that room
    covers them (this is what lets a frozen borrower resume), and further
    free cores are added only if their CPU is measured idle too.
    """
    lent = len(set(cpus) & running_borrowed)
    busy = usage_busy is not None and usage_busy > (len(cpus) - lent) + borrow.busy_threshold
    idle: set[int] = set()
    if usage_idle is None:
        return busy, idle
    room = len(cpus) - in_use_cores(usage_idle, borrow.idle_threshold)
    taken = [cpu for cpu in cpus if cpu in occupied]
    if room >= len(taken):
        idle.update(taken)
        idle.update(_pick_idle_cpus(cpus, room - len(taken), cpu_busy, occupied, borrow.cpu_idle_threshold))
    return busy, idle


def classify_owner(
    cpus: tuple[int, ...],
    *,
    usage_idle: float | None,
    usage_busy: float | None,
    cpu_busy: Mapping[int, float | None],
    running_borrowed: set[int],
    occupied: set[int],
    borrow: BorrowConfig,
) -> dict[int, CoreFlags]:
    """Cores of one running owner job: which are idle (lendable) and whether
    the job needs its lent cores back."""
    busy, idle_cpus = _idle_and_busy(cpus, usage_idle, usage_busy, cpu_busy, running_borrowed, occupied, borrow)
    return {cpu: CoreFlags(owner_idle=cpu in idle_cpus, owner_busy=busy) for cpu in cpus}


class NodeTelemetry:
    """Per-pass sampling for the node worker and ``simq partitions``.

    ``sample()`` returns a :class:`borrow.CoreState` for every core that is
    held outside the queue, allocated to an owner job, or marked
    not-lendable through its borrow lease.  Readers are injectable for tests.
    """

    def __init__(
        self,
        config,
        *,
        clock=None,
        read_table=read_process_table,
        read_cpus=read_cpu_times,
        fd_holders=find_fd_holders,
        lock_holders=None,
    ) -> None:
        import time as _time

        from .locks import lock_holders as _lock_holders

        self.config = config
        self.clock = clock or _time.monotonic
        self.read_table = read_table
        self.read_cpus = read_cpus
        self.fd_holders = fd_holders
        self.lock_holders = lock_holders or _lock_holders
        self.tracker = UsageTracker()
        self.groups: list[HolderGroup] = []
        self.group_trees: dict[str, set[int]] = {}
        self.owner_usage: dict[int, tuple[float | None, float | None]] = {}

    def reset(self) -> None:
        self.tracker.reset()
        self.groups = []
        self.group_trees = {}
        self.owner_usage = {}

    def _external_groups(self, external: Mapping[int, list], table: Mapping[int, ProcInfo], now: float,
                         idle_window: float, busy_window: float) -> list[HolderGroup]:
        config = self.config
        paths_by_cpu = {cpu: config.locks.core_lock_paths(cpu) for cpu in external}
        keys_by_cpu = {cpu: [k for k in (file_key(p) for p in paths) if k] for cpu, paths in paths_by_cpu.items()}
        holders = self.fd_holders({key for keys in keys_by_cpu.values() for key in keys})
        roots_by_cpu: dict[int, frozenset[int]] = {}
        lock_pids: dict[int, set[int]] = {}
        for cpu, found in external.items():
            roots = set()
            for key in keys_by_cpu.get(cpu, ()):
                roots |= holders.get(key, set())
            pids = {h.pid for h in found if h.pid}
            lock_pids[cpu] = pids
            roots |= {pid for pid in pids if pid in table}
            roots_by_cpu[cpu] = frozenset(pid for pid in roots if pid in table and table[pid].state != "Z")
        by_roots: dict[frozenset[int], list[int]] = {}
        for cpu, roots in roots_by_cpu.items():
            # a lease whose holders are not visible: one group per core
            by_roots.setdefault(roots if roots else frozenset({-cpu - 1}), []).append(cpu)
        groups = []
        for roots, cpus in by_roots.items():
            visible = sorted(pid for pid in roots if pid >= 0)
            if not visible:
                cpu = cpus[0]
                pids = sorted(lock_pids.get(cpu, ()))
                description = "holder gone (" + (", ".join(f"locked by pid {p}" for p in pids) or "pid unknown") + ")"
                groups.append(HolderGroup(key=f"stale:{cpu}:{pids}", cpus=tuple(sorted(cpus)),
                                          description=description, usage_idle=None, usage_busy=None, stale=True))
                continue
            identity = ",".join(str(pid) + ":" + str(table[pid].start) for pid in visible)
            tree = descendants(table, visible)
            key = "lease:" + identity
            self.tracker.observe_tree(key, identity, tree, table, now)
            self.group_trees[identity] = tree
            description = "; ".join(f"pid {pid} {command_summary(pid, 70)}" for pid in visible[:4])
            if len(visible) > 4:
                description += f"; +{len(visible) - 4} more"
            groups.append(HolderGroup(
                key=identity,
                cpus=tuple(sorted(cpus)),
                description=description,
                usage_idle=self.tracker.tree_usage(key, idle_window, now),
                usage_busy=self.tracker.tree_usage(key, busy_window, now, full=False),
            ))
        return groups

    def sample(self, *, owners: list[dict], external: Mapping[int, list], borrowers: list) -> dict[int, "CoreState"]:
        from .borrow import CoreState
        from .config import parse_cpu_list

        config = self.config
        cfg = config.borrow
        now = self.clock()
        table = self.read_table()
        self.tracker.observe_cpus(self.read_cpus(), now)
        idle_window, busy_window = cfg.idle_window_s, cfg.busy_window_s
        cpu_busy = {cpu: self.tracker.cpu_busy(cpu, idle_window, now) for cpu in config.node_cpus()}
        occupied = {cpu for b in borrowers for cpu in b.cpus}
        running_borrowed = {cpu for b in borrowers if b.frozen_since is None and not b.pending for cpu in b.cpus}
        flags: dict[int, CoreFlags] = {}
        self.groups = self._external_groups(external, table, now, idle_window, busy_window) if external else []
        for group in self.groups:
            flags.update(classify_external(group, cpu_busy=cpu_busy, running_borrowed=running_borrowed,
                                           occupied=occupied, borrow=cfg))
        self.owner_usage = {}
        for job in owners:
            try:
                cpus = parse_cpu_list(str(job.get("cpu_list") or ""))
            except ValueError:
                continue
            pid = job.get("pid")
            if not cpus or not pid or int(pid) not in table:
                continue
            key = f"job:{job['id']}:{job.get('launch_token')}"
            tree = descendants(table, [int(pid)])
            self.tracker.observe_tree(key, key, tree, table, now)
            usage_idle = self.tracker.tree_usage(key, idle_window, now)
            usage_busy = self.tracker.tree_usage(key, busy_window, now, full=False)
            self.owner_usage[int(job["id"])] = (usage_idle, usage_busy)
            if job.get("exclusive_timing"):
                usage_idle = None  # never lent
            for cpu, item in classify_owner(cpus, usage_idle=usage_idle, usage_busy=usage_busy, cpu_busy=cpu_busy,
                                            running_borrowed=running_borrowed, occupied=occupied,
                                            borrow=cfg).items():
                flags[cpu] = item
        no_borrow: set[int] = set()
        borrow_paths = {path: cpu for cpu in config.node_cpus() for path in config.locks.borrow_lock_paths(cpu)}
        if borrow_paths:
            for path, found in self.lock_holders(list(borrow_paths)).items():
                if any(h.mode == "WRITE" for h in found):
                    no_borrow.add(borrow_paths[path])
        self.tracker.forget_trees(now)
        # A free core is lent only while its CPU is idle: processes pinned
        # there without a lease file are invisible to the queue otherwise.
        # Measured over the idle window, or the history so far (at least two
        # samples); unknown counts as busy.
        cpu_recent = {cpu: self.tracker.cpu_busy(cpu, idle_window, now, full=False) for cpu in config.node_cpus()}
        cpu_busy_flags = {
            cpu: value is None or value > cfg.cpu_idle_threshold for cpu, value in cpu_recent.items()
        }
        states: dict[int, CoreState] = {}
        for cpu in set(flags) | no_borrow | set(config.node_cpus()):
            item = flags.get(cpu, CoreFlags())
            states[cpu] = CoreState(
                cpu=cpu,
                external=item.external,
                external_idle=item.external_idle,
                external_busy=item.external_busy,
                stale=item.stale,
                holder=item.holder,
                holder_key=item.holder_key,
                owner_idle=item.owner_idle,
                owner_busy=item.owner_busy,
                no_borrow=cpu in no_borrow,
                cpu_busy=cpu_busy_flags.get(cpu, False),
            )
        return states
