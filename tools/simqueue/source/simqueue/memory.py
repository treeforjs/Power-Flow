"""Memory tracking and memory-aware admission (SimQueue 0.7.0).

Everything here reads files only: cgroup v2 files of job units under the
user manager (``memory.current``, ``memory.peak``, ``cpu.stat``,
``memory.events``, ``cgroup.events``) and ``/proc/meminfo``.  No process is
spawned, so the worker can sample every few seconds at negligible cost.

Pieces
------
``UnitCgroup`` / ``read_usage``
    Where a job unit's cgroup lives and what it uses now (current memory,
    high-water mark, CPU time, OOM kills, frozen).
``read_node_memory``
    MemTotal and MemAvailable of the node (or of a simulated node, for
    tests: ``memory.test_node_total``).
``Series``
    A job's memory over time as compact min/max buckets (one per minute,
    coarsened when a job runs long), so the database never holds more than
    ``series_buckets`` entries per job.
``SpikeState``
    High/low-water hysteresis on MemAvailable.
``admission``
    The actual-memory admission guard: a start needs its request below
    MemAvailable minus the margin minus what recently started jobs have
    not allocated yet.
``auto_request`` / ``name_family``
    ``--mem auto``: the request derived from the recorded peaks of earlier
    jobs with the same name (else the same name family).
"""

from __future__ import annotations

import json
import math
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from .nodeconfig import format_size

MiB = 1024**2
GiB = 1024**3
CGROUP_ROOT = Path("/sys/fs/cgroup")
MEMINFO = Path("/proc/meminfo")
SERIES_SCHEMA = 1


# ----- cgroup paths -------------------------------------------------------------


def user_manager_cgroup() -> Path:
    """The cgroup of this user's systemd manager (``user@UID.service``).

    Taken from this process's own cgroup when it runs below the manager (the
    worker unit, a terminal scope, a job unit), else the standard path."""
    uid = os.getuid()
    marker = f"user@{uid}.service"
    try:
        for line in Path("/proc/self/cgroup").read_text(encoding="ascii").splitlines():
            if line.startswith("0::"):
                parts = line[3:].strip("/").split("/")
                if marker in parts:
                    return CGROUP_ROOT.joinpath(*parts[: parts.index(marker) + 1])
    except OSError:
        pass
    return CGROUP_ROOT / "user.slice" / f"user-{uid}.slice" / marker


def slice_path(name: str) -> Path:
    """Relative cgroup path of a slice: ``a-b.slice`` is ``a.slice/a-b.slice``."""
    stem = name[: -len(".slice")] if name.endswith(".slice") else name
    if not stem or stem == "-":
        return Path()
    parts = stem.split("-")
    return Path(*(("-".join(parts[: index + 1]) + ".slice") for index in range(len(parts))))


def unit_cgroup_path(unit: str, slice_name: str, *, root: Path | None = None) -> Path:
    return (root or user_manager_cgroup()) / slice_path(slice_name) / unit


def cgroup_from_control_group(value: str) -> Path | None:
    """systemd's ``ControlGroup`` property as a sysfs path."""
    value = (value or "").strip()
    if not value.startswith("/"):
        return None
    return CGROUP_ROOT / value.lstrip("/")


def _read_int(path: Path) -> int | None:
    try:
        with open(path, "rb") as handle:
            return int(handle.read().split()[0])
    except (OSError, ValueError, IndexError):
        return None


def _read_keyed(path: Path) -> dict[str, int] | None:
    try:
        with open(path, "rb") as handle:
            text = handle.read().decode("ascii", "replace")
    except OSError:
        return None
    values: dict[str, int] = {}
    for line in text.splitlines():
        key, _, number = line.partition(" ")
        try:
            values[key] = int(number)
        except ValueError:
            continue
    return values


@dataclass(frozen=True)
class Usage:
    """One reading of a job unit's cgroup."""

    current: int
    peak: int | None  # memory.peak (kernel 5.19+): high-water mark since the cgroup was created
    cpu_usec: int | None  # cpu.stat usage_usec
    oom_kill: int | None  # memory.events oom_kill
    frozen: bool | None  # cgroup.events frozen
    populated: bool | None  # cgroup.events populated


def read_usage(cgroup: Path, *, full: bool = True) -> Usage | None:
    """Read a job cgroup; ``None`` when it does not exist (the unit ended).

    With *full* false only ``memory.current`` and ``memory.peak`` are read
    (the sampling tick); CPU time, OOM count and freezer state are read at
    a job's end and when the worker needs them."""
    current = _read_int(cgroup / "memory.current")
    if current is None:
        return None
    peak = _read_int(cgroup / "memory.peak")
    if not full:
        return Usage(current, peak, None, None, None, None)
    cpu = _read_keyed(cgroup / "cpu.stat") or {}
    events = _read_keyed(cgroup / "memory.events") or {}
    state = _read_keyed(cgroup / "cgroup.events") or {}
    return Usage(
        current=current,
        peak=peak,
        cpu_usec=cpu.get("usage_usec"),
        oom_kill=events.get("oom_kill"),
        frozen=None if "frozen" not in state else bool(state["frozen"]),
        populated=None if "populated" not in state else bool(state["populated"]),
    )


def own_cgroup_usage() -> dict[str, Any] | None:
    """Usage of the cgroup this process runs in, as recorded by jobwrap."""
    try:
        for line in Path("/proc/self/cgroup").read_text(encoding="ascii").splitlines():
            if line.startswith("0::"):
                usage = read_usage(CGROUP_ROOT / line[3:].lstrip("/"))
                return None if usage is None else usage_record(usage)
    except OSError:
        return None
    return None


def usage_record(usage: Usage) -> dict[str, Any]:
    return {"memory_current": usage.current, "memory_peak": usage.peak, "cpu_usage_usec": usage.cpu_usec,
            "oom_kill": usage.oom_kill}


# ----- node memory ----------------------------------------------------------------


@dataclass(frozen=True)
class NodeMemory:
    total: int
    available: int
    simulated: bool = False

    @property
    def used(self) -> int:
        return max(0, self.total - self.available)


def read_meminfo(path: Path = MEMINFO) -> tuple[int, int] | None:
    """(MemTotal, MemAvailable) in bytes."""
    try:
        with open(path, "rb") as handle:
            text = handle.read(4096).decode("ascii", "replace")
    except OSError:
        return None
    total = available = None
    for line in text.splitlines():
        if line.startswith("MemTotal:"):
            total = int(line.split()[1]) * 1024
        elif line.startswith("MemAvailable:"):
            available = int(line.split()[1]) * 1024
        if total is not None and available is not None:
            return total, available
    return None


def read_node_memory(test_total: int | None = None, job_currents: Iterable[int] = ()) -> NodeMemory | None:
    """The node's memory.  With *test_total* (``memory.test_node_total``) the
    node is simulated: its total is that size and only the queue's own job
    units use it, so tests can drive the guard and spike handling with real
    allocations of a few hundred MiB."""
    if test_total is not None:
        used = sum(job_currents)
        return NodeMemory(test_total, max(0, test_total - used), simulated=True)
    found = read_meminfo()
    if found is None:
        return None
    return NodeMemory(found[0], found[1])


def water_marks(config: Any, total: int) -> tuple[int, int]:
    """(high-water, low-water) as MemAvailable thresholds: memory pressure
    starts when MemAvailable falls below the first and ends when it is back
    at or above the second."""
    high = max(config.high_water_available, int(total * config.high_water_percent / 100.0))
    low = max(config.low_water_available, int(total * config.low_water_percent / 100.0))
    return high, max(low, high)


def admission_margin(config: Any, total: int) -> int:
    """MemAvailable that must remain after a start: ``admission_margin``, or
    by default the high-water mark (a job that grows to its request then
    does not by itself cause memory pressure)."""
    if config.admission_margin is not None:
        return config.admission_margin
    return water_marks(config, total)[0]


# ----- series -------------------------------------------------------------------------


class Series:
    """Min/max memory per bucket (``step`` seconds, 60 at first).

    When more than *cap* buckets exist, neighbouring buckets are merged
    and the step doubles, so a job's whole history stays within *cap*
    entries.  Values are stored in MiB."""

    def __init__(self, cap: int = 240, step: int = 60) -> None:
        self.cap = max(8, int(cap))
        self.step = int(step)
        self.buckets: list[list[int]] = []  # [start_unix, min_mib, max_mib]

    def add(self, unix: float, value: int) -> bool:
        """Add a sample; returns True when a new bucket was started."""
        mib = int(math.ceil(value / MiB))
        start = int(unix // self.step * self.step)
        if self.buckets and self.buckets[-1][0] == start:
            bucket = self.buckets[-1]
            bucket[1] = min(bucket[1], mib)
            bucket[2] = max(bucket[2], mib)
            return False
        if self.buckets and start < self.buckets[-1][0]:
            bucket = self.buckets[-1]  # clock went back: fold into the last bucket
            bucket[1] = min(bucket[1], mib)
            bucket[2] = max(bucket[2], mib)
            return False
        self.buckets.append([start, mib, mib])
        while len(self.buckets) > self.cap:
            self._coarsen()
        return True

    def _coarsen(self) -> None:
        self.step *= 2
        merged: list[list[int]] = []
        for bucket in self.buckets:
            start = bucket[0] // self.step * self.step
            if merged and merged[-1][0] == start:
                merged[-1][1] = min(merged[-1][1], bucket[1])
                merged[-1][2] = max(merged[-1][2], bucket[2])
            else:
                merged.append([start, bucket[1], bucket[2]])
        self.buckets = merged

    def to_json(self) -> str:
        return json.dumps({"v": SERIES_SCHEMA, "step": self.step, "unit": "MiB", "b": self.buckets},
                          separators=(",", ":"))

    @classmethod
    def from_json(cls, text: str | None, cap: int = 240) -> "Series":
        series = cls(cap)
        if not text:
            return series
        try:
            value = json.loads(text)
            series.step = int(value.get("step") or 60)
            series.buckets = [[int(a), int(b), int(c)] for a, b, c in value.get("b") or []]
        except (ValueError, TypeError, AttributeError):
            return cls(cap)
        return series

    def max_since(self, unix: float) -> int | None:
        values = [b[2] for b in self.buckets if b[0] + self.step > unix]
        return max(values) * MiB if values else None

    def min_since(self, unix: float) -> int | None:
        values = [b[1] for b in self.buckets if b[0] + self.step > unix]
        return min(values) * MiB if values else None


def sparkline(series: Series, width: int = 40) -> str:
    """The maxima of the last *width* buckets as a text sparkline."""
    blocks = " ▁▂▃▄▅▆▇█"
    values = [b[2] for b in series.buckets[-width:]]
    if not values:
        return ""
    top = max(values) or 1
    return "".join(blocks[min(8, int(round(v / top * 8)))] for v in values)


# ----- spike state ----------------------------------------------------------------


@dataclass
class SpikeState:
    """Memory pressure with hysteresis."""

    active: bool = False
    since: float | None = None
    record: dict[str, Any] = field(default_factory=dict)

    def update(self, available: int, high: int, low: int) -> str | None:
        """Returns "enter", "leave" or None."""
        if not self.active and available < high:
            self.active = True
            return "enter"
        if self.active and available >= low:
            self.active = False
            self.since = None
            return "leave"
        return None


# ----- admission guard --------------------------------------------------------------


@dataclass(frozen=True)
class Admission:
    admitted: list[int]
    held: dict[int, str]
    budget: int  # MemAvailable left for further starts after the admitted ones


def admission(
    requests: list[tuple[int, int]],
    *,
    available: int,
    margin: int,
    ramping: int,
) -> Admission:
    """Admit starts (job id, request bytes) in order while their requests fit
    MemAvailable minus *margin* minus *ramping* (the requests of recently
    started jobs that they have not allocated yet).  A job that does not fit
    is held; later (smaller) ones may still start: the reservations already
    keep its place."""
    budget = available - margin - ramping
    admitted: list[int] = []
    held: dict[int, str] = {}
    for job_id, need in requests:
        if need <= budget:
            admitted.append(job_id)
            budget -= need
            continue
        left = max(0, budget)
        held[job_id] = (
            f"memory: needs {format_size(need)} but MemAvailable {format_gib(available)} "
            f"- margin {format_gib(margin)}"
            + (f" - {format_gib(ramping)} not yet allocated by starting jobs" if ramping else "")
            + f" leaves {format_gib(left)}"
        )
    return Admission(admitted, held, budget)


def format_gib(value: int | None) -> str:
    if value is None:
        return "-"
    if abs(value) >= GiB:
        return f"{value / GiB:.1f}G"
    return f"{value / MiB:.0f}M"


# ----- --mem auto -------------------------------------------------------------------

_INDEX_TOKEN = re.compile(r"^(?:[a-z]{0,6}\d+(?:\.\d+)?[a-z]{0,2}|[a-z])$")


def name_family(name: str) -> str:
    """A job name without its trailing index tokens: ``fe-audit-cold-n4`` ->
    ``fe-audit-cold``, ``lv-interp_300K_xe1e-6`` -> ``lv-interp``,
    ``pore1p5-cpu-a`` -> ``pore1p5-cpu``.  At least one token stays."""
    text = name.strip().lower()
    pieces = re.split(r"([-_.:/ ]+)", text)
    tokens = pieces[0::2]
    separators = pieces[1::2]
    while len(tokens) > 1 and (not tokens[-1] or _INDEX_TOKEN.fullmatch(tokens[-1])):
        tokens.pop()
        separators = separators[: len(tokens) - 1]
    result = tokens[0]
    for separator, token in zip(separators, tokens[1:]):
        result += separator + token
    return result or text


def round_up(value: int, step: int) -> int:
    step = max(1, int(step))
    return int(math.ceil(value / step) * step)


def auto_request(
    history: list[dict[str, Any]],
    *,
    name: str,
    cores: int,
    config: Any,
    min_mem: int,
    max_mem: int,
    default_mem: int,
) -> tuple[int, dict[str, Any]]:
    """The ``--mem auto`` request and its basis.

    *history* holds finished jobs (newest first) with ``name``,
    ``requested_cores``, ``mem_peak_bytes``, ``mem_bytes``, ``result`` and
    ``id``.  Jobs with the same name are used, else jobs of the same name
    family; among those, jobs with the same core count when there are any.
    The request is ``auto_factor`` times the largest of the last
    ``auto_history`` peaks (a job killed at its memory limit counts with its
    limit: it needed more), rounded up to ``auto_round`` and bounded by the
    partition's ``min_mem`` and ``max_mem``.  Without history: the default."""
    exact = [job for job in history if job.get("name") == name]
    family = name_family(name)
    matched_by = "name"
    pool = exact
    if not pool:
        pool = [job for job in history if name_family(str(job.get("name") or "")) == family]
        matched_by = "family"
    same_cores = [job for job in pool if int(job.get("requested_cores") or 0) == cores]
    if same_cores:
        pool = same_cores
    pool = pool[: max(1, int(config.auto_history))]
    if not pool:
        return default_mem, {"source": "default", "family": family, "request": default_mem,
                             "note": "no recorded peak for this name or family; partition default"}
    peaks = []
    for job in pool:
        peak = int(job.get("mem_peak_bytes") or 0)
        if job.get("result") == "oom-kill" and job.get("mem_bytes"):
            peak = max(peak, int(job["mem_bytes"]))
        peaks.append(peak)
    largest = max(peaks)
    wanted = round_up(int(largest * float(config.auto_factor)), int(config.auto_round))
    request = min(max(wanted, min_mem), max_mem)
    basis = {
        "source": "history",
        "matched_by": matched_by,
        "family": family,
        "jobs": [int(job["id"]) for job in pool],
        "peaks": peaks,
        "largest_peak": largest,
        "factor": float(config.auto_factor),
        "wanted": wanted,
        "request": request,
        "same_cores": bool(same_cores),
        "oom_in_history": any(job.get("result") == "oom-kill" for job in pool),
    }
    if wanted > max_mem:
        basis["note"] = f"history suggests {format_size(wanted)}, above the limit {format_size(max_mem)}"
    return request, basis


def auto_time_limit(history: list[dict[str, Any]], *, name: str, default: int | None,
                    maximum: int | None, factor: float = 1.5, count: int = 20, floor: int = 600,
                    step: int = 300) -> tuple[int | None, dict[str, Any]]:
    """0.8.0 ``-t auto``: the 95th percentile of the last *count* runtimes of
    succeeded jobs with the same name (else name family) times *factor*,
    rounded up to *step*, at least *floor*, at most the partition's limit;
    without history the default.  Accurate limits let backfill fill gaps."""
    pool = [job for job in history if job.get("name") == name]
    matched = "name"
    if not pool:
        family = name_family(name)
        pool = [job for job in history if name_family(str(job.get("name") or "")) == family]
        matched = "family"
    runtimes = sorted(float(job["runtime_s"]) for job in pool[:count] if job.get("runtime_s"))
    if not runtimes:
        return default, {"source": "default"}
    p95 = runtimes[min(len(runtimes) - 1, int(math.ceil(0.95 * len(runtimes))) - 1)]
    limit = max(floor, round_up(int(p95 * factor), step))
    if maximum is not None:
        limit = min(limit, maximum)
    return limit, {"source": "history", "matched_by": matched, "runs": len(runtimes), "p95_s": round(p95),
                   "factor": factor, "limit_s": limit}


def describe_basis(basis: dict[str, Any] | None) -> str:
    if not basis:
        return ""
    if basis.get("source") != "history":
        return "auto: " + str(basis.get("note") or "partition default")
    prefix = "auto (default without --mem)" if basis.get("implicit") else "auto"
    jobs = basis.get("jobs") or []
    text = (f"{prefix}: {basis['factor']:g} x largest peak {format_gib(basis['largest_peak'])} of "
            f"{len(jobs)} earlier job(s) with the same {basis.get('matched_by')}"
            + (" and core count" if basis.get("same_cores") else "")
            + f" ({', '.join(str(j) for j in jobs[:6])}{', ...' if len(jobs) > 6 else ''})")
    if basis.get("oom_in_history"):
        text += "; an out-of-memory kill counts with its limit"
    if basis.get("note"):
        text += "; " + str(basis["note"])
    return text


def over_request_warning(
    request: int, history: list[dict[str, Any]], *, name: str, cores: int, config: Any,
) -> str | None:
    """A warning when *request* is far above what this job's name (or family)
    used before, else ``None``."""
    exact = [job for job in history if job.get("name") == name]
    pool = exact or [job for job in history if name_family(str(job.get("name") or "")) == name_family(name)]
    same = [job for job in pool if int(job.get("requested_cores") or 0) == cores]
    pool = (same or pool)[: max(1, int(config.auto_history))]
    if not pool or any(job.get("result") == "oom-kill" for job in pool):
        return None
    largest = max(int(job.get("mem_peak_bytes") or 0) for job in pool)
    if largest <= 0:
        return None
    if request >= largest * float(config.warn_factor) and request - largest >= int(config.warn_min_excess):
        return (f"--mem {format_size(request)} is {request / largest:.1f}x the largest recent peak "
                f"({format_gib(largest)}) of {len(pool)} earlier job(s) named like this one; the unused "
                "reservation keeps other jobs waiting (--mem auto would derive the request from that history)")
    return None
