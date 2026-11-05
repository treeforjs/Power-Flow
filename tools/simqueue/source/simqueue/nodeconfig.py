"""Node and partition configuration for the SimQueue 0.6 node scheduler.

A node configuration file switches SimQueue from the 0.5 serial queue
(``legacy`` mode, driven by ``SIMQUEUE_PCORES``) to the concurrent node
scheduler.  The file is JSON with schema ``simqueue-node-config-v1``.  Its
path is ``$SIMQUEUE_CONFIG`` when set, otherwise ``<state dir>/node.json``.

0.6.3: a partition may carry ``aliases``, further names under which jobs
request it.  An alias is a plain name or an object that also sets the
defaults a job gets when it names the alias (``default_cpus``,
``default_mem``, ``default_time``) and the cores tried first
(``prefer_cpus``); limits, GPUs, lending and the cores themselves stay the
partition's.  Every partition-by-name lookup resolves aliases, so a single
pool can keep the names of the partitions it replaced.  A job keeps the name
it was submitted with (``SIMQ_PARTITION``, ``simq`` listings) and is
scheduled under the partition the name resolves to.

0.7.0: a ``memory`` section (``MemoryConfig``) sets how job memory is
sampled and used: the sampling interval, the actual-memory admission guard,
the high/low-water marks of memory pressure and ``--mem auto``.  A partition
may set ``min_mem`` (the smallest request ``--mem auto`` gives).  ``launch``
gains ``idle_seconds``: the longest the event-driven worker sleeps when
nothing needs attention.
"""

from __future__ import annotations

import dataclasses
import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .config import format_cpu_list, parse_cpu_list, state_dir
from .gpuline import LineConfig, parse_line_config
from .github import GitHubConfig, parse_github_config


CONFIG_SCHEMA = "simqueue-node-config-v1"
PARTITION_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
GPU_UUID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$")
DEFAULT_WRAPPER_PYTHON = "/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen"
UNLIMITED_WORDS = {"unlimited", "infinite", "none"}
ON_PREEMPT_MODES = ("freeze", "requeue", "cancel")
DEFAULT_BORROW_SLICE = "simqborrow.slice"
SIZE_PATTERN = re.compile(r"^([0-9]+)\s*([A-Za-z]*)$")
SIZE_UNITS = {
    # Slurm convention: a bare number is MiB and single-letter suffixes are binary.
    "": 1024**2,
    "k": 1024,
    "m": 1024**2,
    "g": 1024**3,
    "t": 1024**4,
    "b": 1,
    "kib": 1024,
    "mib": 1024**2,
    "gib": 1024**3,
    "tib": 1024**4,
    "kb": 1000,
    "mb": 1000**2,
    "gb": 1000**3,
    "tb": 1000**4,
}


class ConfigError(ValueError):
    """The node configuration file is missing required data or invalid."""


def config_path() -> Path:
    configured = os.environ.get("SIMQUEUE_CONFIG")
    if configured:
        return Path(configured).expanduser().resolve()
    return state_dir() / "node.json"


def parse_size(value: Any, *, field_name: str = "size", bare_unit: int | None = None) -> int:
    """Parse a memory size.

    ``18G``/``18GiB`` are binary (Slurm style), ``18GB`` is decimal, ``123B``
    is bytes and a bare number is MiB (as for ``sbatch --mem``) unless
    *bare_unit* overrides it.  Integers given in JSON are bytes.
    """
    if isinstance(value, bool):
        raise ValueError(f"{field_name} must be a size")
    if isinstance(value, int):
        if value < 0:
            raise ValueError(f"{field_name} must be nonnegative")
        return value
    if not isinstance(value, str):
        raise ValueError(f"{field_name} must be a size string such as 4G")
    match = SIZE_PATTERN.fullmatch(value.strip())
    if match is None:
        raise ValueError(f"invalid {field_name} {value!r}; use e.g. 512M, 4G, 18GiB")
    unit = match.group(2).lower()
    if unit not in SIZE_UNITS:
        raise ValueError(f"invalid {field_name} unit {match.group(2)!r}")
    multiplier = SIZE_UNITS[unit]
    if unit == "" and bare_unit is not None:
        multiplier = bare_unit
    return int(match.group(1)) * multiplier


def format_size(value: int | None) -> str:
    if value is None:
        return "-"
    for suffix, multiplier in (("T", 1024**4), ("G", 1024**3), ("M", 1024**2), ("K", 1024)):
        if value >= multiplier and value % multiplier == 0:
            return f"{value // multiplier}{suffix}"
    if value >= 1024**3:
        return f"{value / 1024**3:.1f}G"
    if value >= 1024**2:
        return f"{value / 1024**2:.1f}M"
    return f"{value}B"


def parse_time_limit(value: Any, *, field_name: str = "time limit") -> int | None:
    """Parse a Slurm time string into seconds; ``None`` means unlimited.

    Accepted: ``MM``, ``MM:SS``, ``HH:MM:SS``, ``D-HH``, ``D-HH:MM``,
    ``D-HH:MM:SS``, plus ``UNLIMITED``.  JSON integers are seconds.
    """
    if value is None:
        return None
    if isinstance(value, bool):
        raise ValueError(f"invalid {field_name}")
    if isinstance(value, int):
        if value <= 0:
            raise ValueError(f"{field_name} must be positive")
        return value
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"invalid {field_name} {value!r}")
    text = value.strip()
    if text.lower() in UNLIMITED_WORDS:
        return None
    days = 0
    rest = text
    if "-" in text:
        day_text, rest = text.split("-", 1)
        if not day_text.isdigit():
            raise ValueError(f"invalid {field_name} {value!r}")
        days = int(day_text)
        parts = rest.split(":")
        if not 1 <= len(parts) <= 3 or not all(part.isdigit() for part in parts):
            raise ValueError(f"invalid {field_name} {value!r}")
        numbers = [int(part) for part in parts] + [0] * (3 - len(parts))
        hours, minutes, seconds = numbers
    else:
        parts = rest.split(":")
        if not 1 <= len(parts) <= 3 or not all(part.isdigit() for part in parts):
            raise ValueError(f"invalid {field_name} {value!r}")
        numbers = [int(part) for part in parts]
        if len(numbers) == 1:
            hours, minutes, seconds = 0, numbers[0], 0
        elif len(numbers) == 2:
            hours, minutes, seconds = 0, numbers[0], numbers[1]
        else:
            hours, minutes, seconds = numbers
    total = ((days * 24 + hours) * 60 + minutes) * 60 + seconds
    if total <= 0:
        raise ValueError(f"{field_name} must be positive")
    return total


def format_time_limit(seconds: int | None) -> str:
    if seconds is None:
        return "UNLIMITED"
    days, remainder = divmod(int(seconds), 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, secs = divmod(remainder, 60)
    if days:
        return f"{days}-{hours:02d}:{minutes:02d}:{secs:02d}"
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


@dataclass(frozen=True)
class GpuSpec:
    uuid: str
    index: int | None
    name: str
    memory_bytes: int
    budget_bytes: int


@dataclass(frozen=True)
class LockConfig:
    enabled: bool = False
    core_path: str | None = None
    core_overrides: dict[int, str] = field(default_factory=dict)
    gpu_path: str | None = None
    gpu_extra: dict[str, tuple[str, ...]] = field(default_factory=dict)
    borrow_path: str | None = None

    def borrow_lock_paths(self, cpu: int) -> tuple[str, ...]:
        """Shared "borrow" lease of a core (0.6.1).

        A borrower holds this file shared instead of the core's exclusive
        lease, so the core's own lease stays free for its owner.  Holding it
        exclusively from outside the queue marks the core as not lendable.
        Default: ``borrow-core<K>.lock`` next to the core's lease file.
        """
        if not self.enabled:
            return ()
        if self.borrow_path:
            return (self.borrow_path.replace("{core}", str(cpu)),)
        core = self.core_lock_paths(cpu)
        if not core:
            return ()
        return (str(Path(core[0]).with_name(f"borrow-core{cpu}.lock")),)

    def core_lock_paths(self, cpu: int) -> tuple[str, ...]:
        if not self.enabled:
            return ()
        if cpu in self.core_overrides:
            return (self.core_overrides[cpu],)
        if self.core_path:
            return (self.core_path.replace("{core}", str(cpu)),)
        return ()

    def gpu_lock_paths(self, gpu: GpuSpec) -> tuple[str, ...]:
        if not self.enabled:
            return ()
        paths: list[str] = []
        if self.gpu_path:
            paths.append(
                self.gpu_path.replace("{uuid}", gpu.uuid).replace(
                    "{index}", "" if gpu.index is None else str(gpu.index)
                )
            )
        paths.extend(self.gpu_extra.get(gpu.uuid, ()))
        return tuple(paths)

    def gpu_lock_order(self, gpu: GpuSpec) -> tuple[str, ...]:
        """The GPU's lease files in acquisition order (0.6.2): the global
        files (``gpu_extra``, e.g. testing-gpu.lock) first, then the per-GPU
        file (testing-gpu-<UUID>.lock), as raw flock users take them."""
        if not self.enabled:
            return ()
        paths: list[str] = []
        for extra in self.gpu_extra.get(gpu.uuid, ()):
            if extra not in paths:
                paths.append(extra)
        if self.gpu_path:
            own = self.gpu_path.replace("{uuid}", gpu.uuid).replace(
                "{index}", "" if gpu.index is None else str(gpu.index)
            )
            if own not in paths:
                paths.append(own)
        return tuple(paths)


@dataclass(frozen=True)
class Partition:
    name: str
    cpus: tuple[int, ...]
    gpus: tuple[str, ...]
    default_mem: int
    max_mem: int
    default_time: int | None
    max_time: int | None
    default_cpus: int | None  # None means every core of the partition
    swap_max: int
    enforce_time_limit: bool
    nice: int
    description: str = ""
    # 0.6.1 lending policy (see OPERATIONS.md, "Borrowing idle cores").
    lend: bool = True
    lend_leased_idle: bool = False
    timing_sensitive: bool = False
    # 0.6.3: further names of this partition, and (in the view of an alias,
    # see NodeConfig.find_partition) the alias used and its preferred cores.
    aliases: tuple[str, ...] = ()
    alias: str | None = None
    prefer_cpus: tuple[int, ...] = ()
    # 0.7.0: the smallest request --mem auto gives
    min_mem: int = 256 * 1024**2
    # 0.7.1: a job without --mem gets --mem auto (an alias with default_mem "auto")
    mem_auto: bool = False

    def resolve_default_cpus(self) -> int:
        return len(self.cpus) if self.default_cpus is None else self.default_cpus

    def lends(self) -> bool:
        """Whether jobs of other partitions may borrow this partition's cores."""
        return self.lend and not self.timing_sensitive

    def display_name(self) -> str:
        """``pool`` or ``pool (as alias)`` for messages."""
        return self.name if self.alias is None else f"{self.name} (as {self.alias})"


@dataclass(frozen=True)
class PartitionAlias:
    """Another name of a partition (0.6.3), with optional own defaults."""

    name: str
    partition: str
    description: str = ""
    # Only the keys given in the configuration: default_cpus (int, or None
    # for every core), default_mem (bytes), default_time (seconds or None).
    defaults: tuple[tuple[str, Any], ...] = ()
    prefer_cpus: tuple[int, ...] = ()


@dataclass(frozen=True)
class BorrowConfig:
    """Node-wide settings for opportunistic (``--borrow``) jobs."""

    enabled: bool = True
    slice: str = DEFAULT_BORROW_SLICE
    slice_cpu_weight: str = "idle"
    unit_cpu_weight: int = 1
    idle_window_s: float = 30.0
    busy_window_s: float = 2.0
    idle_threshold: float = 0.10
    cpu_idle_threshold: float = 0.15
    busy_threshold: float = 0.25
    freeze_max_s: float | None = 1800.0
    default_on_preempt: str = "requeue"
    # 0.6.3: a free core that a waiting owner job needs may still be lent to
    # a borrow job whose time limit ends no later than that owner job's
    # earliest possible start (EASY backfill for borrow jobs).  False keeps
    # the 0.6.1 rule: such cores are never lent while the owner job waits.
    backfill_waiting_owners: bool = False


@dataclass(frozen=True)
class LaunchConfig:
    backend: str = "systemd-user"
    unit_prefix: str = "simq-job"
    slice: str = "simq.slice"
    systemd_run: str = "systemd-run"
    systemctl: str = "systemctl"
    wrapper_python: str = DEFAULT_WRAPPER_PYTHON
    kill_grace_seconds: float = 10.0
    # 0.7.0: the worker is event-driven.  poll_seconds is its pass interval
    # while something needs watching (a job starting or being stopped,
    # borrow jobs, the GPU line, queued jobs held by leases outside the
    # queue or by the memory guard); otherwise it sleeps up to idle_seconds
    # between passes, woken at once by submissions, cancellations and job
    # ends.
    poll_seconds: float = 1.0
    idle_seconds: float = 30.0


@dataclass(frozen=True)
class MemoryConfig:
    """Node-wide memory tracking and memory-aware dispatch (0.7.0)."""

    enabled: bool = True
    # running jobs' cgroups and MemAvailable are read this often (no process is spawned)
    sample_seconds: float = 3.0
    # samples are written to the database at most this often (min/max buckets)
    record_seconds: float = 60.0
    series_buckets: int = 240
    # admission guard: a start needs MemAvailable - margin - (unallocated
    # requests of jobs started less than ramp_seconds ago) >= its request;
    # margin None = the high-water mark
    admission_margin: int | None = None
    ramp_seconds: float = 60.0
    # memory pressure: MemAvailable below max(high_water_available,
    # high_water_percent of MemTotal) stops dispatch and preempts borrow
    # jobs; at or above max(low_water_*) it ends
    high_water_available: int = 2 * 1024**3
    high_water_percent: float = 7.0
    low_water_available: int = 4 * 1024**3
    low_water_percent: float = 12.0
    # a freeze-policy borrow job frozen for memory pressure is stopped
    # (requeued) when the pressure lasts this long: a frozen job keeps its memory
    pressure_escalate_seconds: float = 120.0
    # a job using at least this fraction of its request is warned about once
    near_limit_fraction: float = 0.9
    # 0.7.1: a job submitted without --mem gets --mem auto when its name (or
    # family) has a recorded peak, else the partition's or alias's default
    auto_default: bool = False
    # --mem auto: auto_factor x the largest of the last auto_history peaks,
    # rounded up to auto_round
    auto_factor: float = 1.25
    auto_history: int = 5
    auto_round: int = 256 * 1024**2
    # submit warns when --mem is at least warn_factor x the recent peak and
    # at least warn_min_excess above it
    warn_factor: float = 2.0
    warn_min_excess: int = 2 * 1024**3
    # tests only: simulate a node of this size used only by queue jobs
    test_node_total: int | None = None


@dataclass(frozen=True)
class NodeConfig:
    path: Path
    node: str
    default_partition: str
    partitions: dict[str, Partition]
    gpus: dict[str, GpuSpec]
    locks: LockConfig
    launch: LaunchConfig
    memory_schedulable: int | None
    allow_smt_siblings: bool
    raw: dict[str, Any]
    borrow: BorrowConfig = field(default_factory=BorrowConfig)
    smt_siblings: dict[int, tuple[int, ...]] = field(default_factory=dict)
    # 0.6.2: the GPU line (fair hand-off of the GPU lease files; gpuline.py)
    line: LineConfig = field(default_factory=lambda: LineConfig(enabled=False))
    # 0.6.3: alias name -> alias (see PartitionAlias)
    aliases: dict[str, PartitionAlias] = field(default_factory=dict)
    # 0.7.0: memory tracking, admission guard, pressure, --mem auto
    memory: MemoryConfig = field(default_factory=MemoryConfig)
    # 0.7.0: GitHub reporting of job state (github.py)
    github: GitHubConfig = field(default_factory=GitHubConfig)

    def canonical_partition(self, name: str | None) -> str:
        """The partition a requested or recorded name stands for: the name
        itself, the partition of an alias, or the default partition for no
        name.  An unknown name is returned unchanged."""
        key = name or self.default_partition
        alias = self.aliases.get(key)
        return alias.partition if alias is not None else key

    def find_partition(self, name: str | None) -> Partition | None:
        """The partition for a name, as seen through an alias: the alias's
        defaults and preferred cores replace the partition's own.  ``None``
        when the name is unknown."""
        key = name or self.default_partition
        alias = self.aliases.get(key)
        if alias is None:
            return self.partitions.get(key)
        base = self.partitions.get(alias.partition)
        if base is None:
            return None
        return dataclasses.replace(base, alias=alias.name, prefer_cpus=alias.prefer_cpus, **dict(alias.defaults))

    def partition(self, name: str | None) -> Partition:
        found = self.find_partition(name)
        if found is None:
            key = name or self.default_partition
            known = ", ".join(sorted(self.partitions))
            if self.aliases:
                known += "; aliases: " + ", ".join(sorted(self.aliases))
            raise ValueError(f"unknown partition {key!r} (configured: {known})")
        return found

    def partition_names(self, name: str) -> tuple[str, ...]:
        """Every name that resolves to the same partition as *name*."""
        canonical = self.canonical_partition(name)
        return (canonical, *(a for a, alias in self.aliases.items() if alias.partition == canonical))

    def node_cpus(self) -> tuple[int, ...]:
        cpus: set[int] = set()
        for partition in self.partitions.values():
            cpus.update(partition.cpus)
        return tuple(sorted(cpus))

    def partitions_for_cpu(self, cpu: int) -> list[str]:
        return [name for name, part in self.partitions.items() if cpu in part.cpus]

    def partition_of(self, cpu: int) -> Partition | None:
        for part in self.partitions.values():
            if cpu in part.cpus:
                return part
        return None

    def siblings(self, cpu: int) -> tuple[int, ...]:
        """Other logical CPUs on the same physical core as *cpu*."""
        return tuple(c for c in self.smt_siblings.get(cpu, ()) if c != cpu)

    def with_siblings(self, cpus: Any) -> set[int]:
        result: set[int] = set()
        for cpu in cpus:
            result.add(cpu)
            result.update(self.siblings(cpu))
        return result


def _read_text(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8").strip()
    except OSError:
        return None


def physical_core_of(cpu: int) -> int | None:
    """Return the lowest logical CPU sharing *cpu*'s physical core."""
    siblings = _read_text(
        Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")
    )
    if siblings is None:
        return None
    try:
        values = parse_cpu_list(siblings)
    except ValueError:
        return None
    return min(values) if values else None


def sibling_map(cpus: Any) -> dict[int, tuple[int, ...]]:
    """SMT siblings of each CPU in *cpus*, read from sysfs."""
    result: dict[int, tuple[int, ...]] = {}
    for cpu in cpus:
        text = _read_text(Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list"))
        if text is None:
            continue
        try:
            values = parse_cpu_list(text)
        except ValueError:
            continue
        result[cpu] = tuple(v for v in values if v != cpu)
    return result


def online_cpus() -> tuple[int, ...] | None:
    text = _read_text(Path("/sys/devices/system/cpu/online"))
    if text is None:
        return None
    try:
        return parse_cpu_list(text)
    except ValueError:
        return None


def _object(value: Any, name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ConfigError(f"{name} must be a JSON object")
    return value


def _unknown(value: dict[str, Any], allowed: set[str], name: str) -> None:
    extra = sorted(set(value) - allowed)
    if extra:
        raise ConfigError(f"{name} has unknown key(s): {', '.join(extra)}")


def _cpu_list(value: Any, name: str) -> tuple[int, ...]:
    if isinstance(value, str):
        try:
            return parse_cpu_list(value)
        except ValueError as error:
            raise ConfigError(f"{name}: {error}") from error
    if isinstance(value, list) and all(isinstance(v, int) and not isinstance(v, bool) and v >= 0 for v in value):
        return tuple(sorted(set(value)))
    raise ConfigError(f"{name} must be a CPU list string or an array of CPU numbers")


def parse_node_config(document: dict[str, Any], *, path: Path, check_topology: bool = True) -> NodeConfig:
    _object(document, "node config")
    _unknown(
        document,
        {
            "schema",
            "node",
            "default_partition",
            "memory_schedulable",
            "allow_smt_siblings",
            "gpus",
            "lock_files",
            "launch",
            "partitions",
            "notes",
            "borrow",
            "smt_siblings",
            "gpu_line",
            "memory",
            "github",
        },
        "node config",
    )
    if document.get("schema") != CONFIG_SCHEMA:
        raise ConfigError(f"node config schema must be {CONFIG_SCHEMA!r}")
    node = document.get("node", "localhost")
    if not isinstance(node, str) or not node:
        raise ConfigError("node must be a nonempty string")
    allow_smt = document.get("allow_smt_siblings", False)
    if not isinstance(allow_smt, bool):
        raise ConfigError("allow_smt_siblings must be true or false")

    gpus: dict[str, GpuSpec] = {}
    raw_gpus = document.get("gpus", [])
    if not isinstance(raw_gpus, list):
        raise ConfigError("gpus must be an array")
    for index, raw in enumerate(raw_gpus):
        item = _object(raw, f"gpus[{index}]")
        _unknown(item, {"uuid", "index", "name", "memory", "budget"}, f"gpus[{index}]")
        uuid = item.get("uuid")
        if not isinstance(uuid, str) or GPU_UUID_PATTERN.fullmatch(uuid) is None:
            raise ConfigError(f"gpus[{index}].uuid must be a GPU UUID")
        if uuid in gpus:
            raise ConfigError(f"duplicate GPU {uuid}")
        gpu_index = item.get("index")
        if gpu_index is not None and (isinstance(gpu_index, bool) or not isinstance(gpu_index, int)):
            raise ConfigError(f"gpus[{index}].index must be an integer")
        try:
            memory = parse_size(item.get("memory", 0), field_name=f"gpus[{index}].memory")
            budget = parse_size(item.get("budget", item.get("memory", 0)), field_name=f"gpus[{index}].budget")
        except ValueError as error:
            raise ConfigError(str(error)) from error
        if budget <= 0:
            raise ConfigError(f"gpus[{index}].budget must be positive")
        if memory and budget > memory:
            raise ConfigError(f"gpus[{index}].budget exceeds its memory")
        gpus[uuid] = GpuSpec(uuid, gpu_index, str(item.get("name", "")), memory, budget)

    lock_raw = _object(document.get("lock_files", {"enabled": False}), "lock_files")
    _unknown(lock_raw, {"enabled", "core_path", "core_overrides", "gpu_path", "gpu_extra", "borrow_path", "notes"},
             "lock_files")
    enabled = lock_raw.get("enabled", False)
    if not isinstance(enabled, bool):
        raise ConfigError("lock_files.enabled must be true or false")
    overrides_raw = _object(lock_raw.get("core_overrides", {}), "lock_files.core_overrides")
    overrides: dict[int, str] = {}
    for key, value in overrides_raw.items():
        if not str(key).isdigit() or not isinstance(value, str) or not value.startswith("/"):
            raise ConfigError("lock_files.core_overrides maps a CPU number to an absolute path")
        overrides[int(key)] = value
    extra_raw = _object(lock_raw.get("gpu_extra", {}), "lock_files.gpu_extra")
    extras: dict[str, tuple[str, ...]] = {}
    for key, value in extra_raw.items():
        if not isinstance(value, list) or not all(isinstance(v, str) and v.startswith("/") for v in value):
            raise ConfigError("lock_files.gpu_extra maps a GPU UUID to an array of absolute paths")
        extras[str(key)] = tuple(value)
    for key in ("core_path", "gpu_path", "borrow_path"):
        value = lock_raw.get(key)
        if value is not None and (not isinstance(value, str) or not value.startswith("/")):
            raise ConfigError(f"lock_files.{key} must be an absolute path template")
    if enabled and not lock_raw.get("core_path") and not overrides:
        raise ConfigError("lock_files.enabled requires core_path or core_overrides")
    locks = LockConfig(
        enabled=enabled,
        core_path=lock_raw.get("core_path"),
        core_overrides=overrides,
        gpu_path=lock_raw.get("gpu_path"),
        gpu_extra=extras,
        borrow_path=lock_raw.get("borrow_path"),
    )

    launch_raw = _object(document.get("launch", {}), "launch")
    _unknown(
        launch_raw,
        {"backend", "unit_prefix", "slice", "systemd_run", "systemctl", "wrapper_python",
         "kill_grace_seconds", "poll_seconds", "idle_seconds"},
        "launch",
    )
    backend = launch_raw.get("backend", "systemd-user")
    if backend != "systemd-user":
        raise ConfigError("launch.backend must be 'systemd-user'")
    unit_prefix = launch_raw.get("unit_prefix", "simq-job")
    if not isinstance(unit_prefix, str) or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,40}", unit_prefix) is None:
        raise ConfigError("launch.unit_prefix must be a short unit-name prefix")
    slice_name = launch_raw.get("slice", "simq.slice")
    if not isinstance(slice_name, str) or not slice_name.endswith(".slice"):
        raise ConfigError("launch.slice must name a .slice unit")
    try:
        grace = float(launch_raw.get("kill_grace_seconds", 10.0))
        poll = float(launch_raw.get("poll_seconds", 1.0))
        idle = float(launch_raw.get("idle_seconds", max(30.0, poll)))
    except (TypeError, ValueError) as error:
        raise ConfigError("launch timing values must be numbers") from error
    if not 0.1 <= grace <= 3600 or not 0.05 <= poll <= 60 or not poll <= idle <= 300:
        raise ConfigError("launch timing values are out of range")
    launch = LaunchConfig(
        backend=backend,
        unit_prefix=unit_prefix,
        slice=slice_name,
        systemd_run=str(launch_raw.get("systemd_run", "systemd-run")),
        systemctl=str(launch_raw.get("systemctl", "systemctl")),
        wrapper_python=str(launch_raw.get("wrapper_python", DEFAULT_WRAPPER_PYTHON)),
        kill_grace_seconds=grace,
        poll_seconds=poll,
        idle_seconds=idle,
    )

    try:
        memory_schedulable = (
            None
            if document.get("memory_schedulable") is None
            else parse_size(document["memory_schedulable"], field_name="memory_schedulable")
        )
    except ValueError as error:
        raise ConfigError(str(error)) from error

    raw_partitions = document.get("partitions")
    if not isinstance(raw_partitions, list) or not raw_partitions:
        raise ConfigError("partitions must be a nonempty array")
    partitions: dict[str, Partition] = {}
    aliases: dict[str, PartitionAlias] = {}
    online = online_cpus() if check_topology else None
    for index, raw in enumerate(raw_partitions):
        item = _object(raw, f"partitions[{index}]")
        _unknown(
            item,
            {"name", "cpus", "cores", "gpus", "default_mem", "max_mem", "default_time", "max_time",
             "default_cpus", "swap_max", "enforce_time_limit", "nice", "description",
             "lend", "lend_leased_idle", "timing_sensitive", "aliases", "min_mem"},
            f"partitions[{index}]",
        )
        name = item.get("name")
        if not isinstance(name, str) or PARTITION_NAME_PATTERN.fullmatch(name) is None:
            raise ConfigError(f"partitions[{index}].name is invalid")
        if name in partitions:
            raise ConfigError(f"duplicate partition {name}")
        cpus = _cpu_list(item.get("cpus", item.get("cores")), f"partition {name} cpus")
        if not cpus:
            raise ConfigError(f"partition {name} has no CPUs")
        if check_topology:
            if online is not None and not set(cpus) <= set(online):
                raise ConfigError(f"partition {name} names offline CPUs")
            if not allow_smt:
                cores_seen: dict[int, int] = {}
                for cpu in cpus:
                    core = physical_core_of(cpu)
                    if core is None:
                        continue
                    if core != cpu:
                        raise ConfigError(
                            f"partition {name}: CPU {cpu} is an SMT sibling of CPU {core}; "
                            "name each physical core by its lowest logical CPU"
                        )
                    if core in cores_seen:
                        raise ConfigError(f"partition {name}: CPUs share physical core {core}")
                    cores_seen[core] = cpu
        part_gpus = item.get("gpus", [])
        if not isinstance(part_gpus, list) or not all(isinstance(g, str) for g in part_gpus):
            raise ConfigError(f"partition {name} gpus must be an array of GPU UUIDs")
        for gpu in part_gpus:
            if gpu not in gpus:
                raise ConfigError(f"partition {name} names unknown GPU {gpu}")
        try:
            default_mem = parse_size(item.get("default_mem", "1G"), field_name=f"partition {name} default_mem")
            max_mem = parse_size(item.get("max_mem", item.get("default_mem", "1G")), field_name=f"partition {name} max_mem")
            swap_max = parse_size(item.get("swap_max", "0"), field_name=f"partition {name} swap_max", bare_unit=1024**2)
            min_mem = parse_size(item.get("min_mem", "256M"), field_name=f"partition {name} min_mem")
            default_time = parse_time_limit(item.get("default_time"), field_name=f"partition {name} default_time")
            max_time = parse_time_limit(item.get("max_time"), field_name=f"partition {name} max_time")
        except ValueError as error:
            raise ConfigError(str(error)) from error
        if default_mem <= 0 or max_mem < default_mem:
            raise ConfigError(f"partition {name}: need 0 < default_mem <= max_mem")
        if min_mem <= 0 or min_mem > max_mem:
            raise ConfigError(f"partition {name}: need 0 < min_mem <= max_mem")
        if memory_schedulable is not None and max_mem > memory_schedulable:
            raise ConfigError(f"partition {name}: max_mem exceeds memory_schedulable")
        if max_time is not None and (default_time is None or default_time > max_time):
            raise ConfigError(f"partition {name}: default_time must be set and <= max_time")
        default_cpus = item.get("default_cpus", "all")
        if default_cpus == "all":
            resolved_default: int | None = None
        elif isinstance(default_cpus, int) and not isinstance(default_cpus, bool) and 1 <= default_cpus <= len(cpus):
            resolved_default = default_cpus
        else:
            raise ConfigError(f"partition {name}: default_cpus must be 'all' or 1..{len(cpus)}")
        enforce = item.get("enforce_time_limit", True)
        nice = item.get("nice", 0)
        if not isinstance(enforce, bool):
            raise ConfigError(f"partition {name}: enforce_time_limit must be true or false")
        if isinstance(nice, bool) or not isinstance(nice, int) or not 0 <= nice <= 19:
            raise ConfigError(f"partition {name}: nice must be 0..19")
        policy: dict[str, bool] = {}
        for key, default in (("lend", True), ("lend_leased_idle", False), ("timing_sensitive", False)):
            value = item.get(key, default)
            if not isinstance(value, bool):
                raise ConfigError(f"partition {name}: {key} must be true or false")
            policy[key] = value
        partitions[name] = Partition(
            name=name,
            cpus=cpus,
            gpus=tuple(part_gpus),
            default_mem=default_mem,
            max_mem=max_mem,
            default_time=default_time,
            max_time=max_time,
            default_cpus=resolved_default,
            swap_max=swap_max,
            enforce_time_limit=enforce,
            nice=nice,
            description=str(item.get("description", "")),
            lend=policy["lend"],
            lend_leased_idle=policy["lend_leased_idle"],
            timing_sensitive=policy["timing_sensitive"],
            min_mem=min_mem,
        )
        for alias in _parse_aliases(item.get("aliases", []), partitions[name]):
            if alias.name in aliases:
                raise ConfigError(f"duplicate partition alias {alias.name}")
            aliases[alias.name] = alias
        partitions[name] = dataclasses.replace(
            partitions[name], aliases=tuple(a.name for a in aliases.values() if a.partition == name)
        )
    clash = sorted(set(aliases) & set(partitions))
    if clash:
        raise ConfigError(f"partition alias(es) {', '.join(clash)} also name a partition")
    default_partition = document.get("default_partition")
    if default_partition not in partitions:
        raise ConfigError("default_partition must name a configured partition")
    borrow = _parse_borrow(document.get("borrow", {}))
    memory = _parse_memory(document.get("memory", {}))
    try:
        github = parse_github_config(document.get("github"))
    except ValueError as error:
        raise ConfigError(str(error)) from error
    node_cpus = sorted({cpu for part in partitions.values() for cpu in part.cpus})
    raw_siblings = document.get("smt_siblings")
    if raw_siblings is not None:
        sibling_raw = _object(raw_siblings, "smt_siblings")
        siblings: dict[int, tuple[int, ...]] = {}
        for key, value in sibling_raw.items():
            if not str(key).isdigit():
                raise ConfigError("smt_siblings maps a CPU number to a CPU list")
            siblings[int(key)] = tuple(c for c in _cpu_list(value, f"smt_siblings[{key}]") if c != int(key))
        if check_topology:
            # 0.6.3: a written sibling map must describe this machine's CPU
            measured = sibling_map(siblings)
            wrong = sorted(cpu for cpu, found in measured.items() if set(found) != set(siblings[cpu]))
            if wrong:
                raise ConfigError(
                    f"smt_siblings disagrees with this machine's topology for CPU(s) {format_cpu_list(wrong)}"
                )
    elif check_topology:
        siblings = sibling_map(node_cpus)
    else:
        siblings = {}
    try:
        line = parse_line_config(document)
    except ValueError as error:
        raise ConfigError(str(error)) from error
    return NodeConfig(
        path=path,
        node=node,
        default_partition=str(default_partition),
        partitions=partitions,
        gpus=gpus,
        locks=locks,
        launch=launch,
        memory_schedulable=memory_schedulable,
        allow_smt_siblings=allow_smt,
        raw=document,
        borrow=borrow,
        smt_siblings=siblings,
        line=line,
        aliases=aliases,
        memory=memory,
        github=github,
    )


def _parse_aliases(value: Any, partition: Partition) -> list[PartitionAlias]:
    """The ``aliases`` of one partition: names, or objects with defaults."""
    where = f"partition {partition.name} aliases"
    if not isinstance(value, list):
        raise ConfigError(f"{where} must be an array of names or objects")
    result: list[PartitionAlias] = []
    for index, raw in enumerate(value):
        item = {"name": raw} if isinstance(raw, str) else _object(raw, f"{where}[{index}]")
        _unknown(item, {"name", "description", "default_cpus", "default_mem", "default_time", "prefer_cpus"},
                 f"{where}[{index}]")
        name = item.get("name")
        if not isinstance(name, str) or PARTITION_NAME_PATTERN.fullmatch(name) is None:
            raise ConfigError(f"{where}[{index}].name is invalid")
        label = f"partition {partition.name} alias {name}"
        defaults: list[tuple[str, Any]] = []
        if "default_cpus" in item:
            cpus_value = item["default_cpus"]
            if cpus_value == "all":
                defaults.append(("default_cpus", None))
            elif (isinstance(cpus_value, int) and not isinstance(cpus_value, bool)
                  and 1 <= cpus_value <= len(partition.cpus)):
                defaults.append(("default_cpus", cpus_value))
            else:
                raise ConfigError(f"{label}: default_cpus must be 'all' or 1..{len(partition.cpus)}")
        try:
            if "default_mem" in item and str(item["default_mem"]).strip().lower() == "auto":
                # 0.7.1: --mem auto when -p names this alias and --mem is not given
                # (the partition's default_mem without history)
                defaults.append(("mem_auto", True))
            elif "default_mem" in item:
                mem = parse_size(item["default_mem"], field_name=f"{label} default_mem")
                if mem <= 0 or mem > partition.max_mem:
                    raise ConfigError(f"{label}: need 0 < default_mem <= the partition's max_mem")
                defaults.append(("default_mem", mem))
            if "default_time" in item:
                limit = parse_time_limit(item["default_time"], field_name=f"{label} default_time")
                if partition.max_time is not None and (limit is None or limit > partition.max_time):
                    raise ConfigError(f"{label}: default_time must be set and <= the partition's max_time")
                defaults.append(("default_time", limit))
        except ValueError as error:
            if isinstance(error, ConfigError):
                raise
            raise ConfigError(str(error)) from error
        prefer: tuple[int, ...] = ()
        if "prefer_cpus" in item:
            raw_prefer = item["prefer_cpus"]
            if isinstance(raw_prefer, list):
                # an array keeps its order: the first listed core is tried first
                if not all(isinstance(c, int) and not isinstance(c, bool) and c >= 0 for c in raw_prefer):
                    raise ConfigError(f"{label}: prefer_cpus must be a CPU list or an array of CPU numbers")
                prefer = tuple(dict.fromkeys(raw_prefer))
            else:
                prefer = _cpu_list(raw_prefer, f"{label} prefer_cpus")
            outside = sorted(set(prefer) - set(partition.cpus))
            if outside:
                raise ConfigError(f"{label}: prefer_cpus {format_cpu_list(outside)} are not in the partition")
        result.append(PartitionAlias(
            name=name,
            partition=partition.name,
            description=str(item.get("description", "")),
            defaults=tuple(defaults),
            prefer_cpus=prefer,
        ))
    return result


def _parse_memory(value: Any) -> MemoryConfig:
    """The ``memory`` section (0.7.0); every key is optional."""
    raw = _object(value, "memory")
    defaults = MemoryConfig()
    sizes = ("admission_margin", "high_water_available", "low_water_available", "auto_round", "warn_min_excess",
             "test_node_total")
    numbers = {
        "sample_seconds": (0.5, 60.0),
        "record_seconds": (1.0, 3600.0),
        "ramp_seconds": (0.0, 3600.0),
        "high_water_percent": (0.0, 50.0),
        "low_water_percent": (0.0, 60.0),
        "pressure_escalate_seconds": (0.0, 86400.0),
        "near_limit_fraction": (0.1, 1.0),
        "auto_factor": (1.0, 4.0),
        "warn_factor": (1.0, 100.0),
    }
    integers = {"series_buckets": (16, 10000), "auto_history": (1, 100)}
    _unknown(raw, {"enabled", "auto_default", "notes", *sizes, *numbers, *integers}, "memory")
    values: dict[str, Any] = {}
    enabled = raw.get("enabled", True)
    if not isinstance(enabled, bool):
        raise ConfigError("memory.enabled must be true or false")
    values["enabled"] = enabled
    auto_default = raw.get("auto_default", False)
    if not isinstance(auto_default, bool):
        raise ConfigError("memory.auto_default must be true or false")
    values["auto_default"] = auto_default
    for key in sizes:
        if key not in raw:
            continue
        if raw[key] is None and key in ("admission_margin", "test_node_total"):
            values[key] = None
            continue
        try:
            values[key] = parse_size(raw[key], field_name=f"memory.{key}")
        except ValueError as error:
            raise ConfigError(str(error)) from error
    for key, (low, high) in numbers.items():
        if key not in raw:
            continue
        item = raw[key]
        if isinstance(item, bool) or not isinstance(item, (int, float)) or not low <= float(item) <= high:
            raise ConfigError(f"memory.{key} must be a number between {low:g} and {high:g}")
        values[key] = float(item)
    for key, (low, high) in integers.items():
        if key not in raw:
            continue
        item = raw[key]
        if isinstance(item, bool) or not isinstance(item, int) or not low <= item <= high:
            raise ConfigError(f"memory.{key} must be an integer between {low} and {high}")
        values[key] = item
    config = dataclasses.replace(defaults, **values)
    if config.auto_round <= 0:
        raise ConfigError("memory.auto_round must be positive")
    if config.record_seconds < config.sample_seconds:
        raise ConfigError("memory.record_seconds must be at least memory.sample_seconds")
    if (config.low_water_available < config.high_water_available
            or config.low_water_percent < config.high_water_percent):
        raise ConfigError("memory: the low-water mark must not be below the high-water mark")
    return config


def describe_memory(config: MemoryConfig) -> dict[str, Any]:
    return {
        "enabled": config.enabled,
        "sample_seconds": config.sample_seconds,
        "record_seconds": config.record_seconds,
        "admission_margin": config.admission_margin,
        "ramp_seconds": config.ramp_seconds,
        "high_water_available": config.high_water_available,
        "high_water_percent": config.high_water_percent,
        "low_water_available": config.low_water_available,
        "low_water_percent": config.low_water_percent,
        "pressure_escalate_seconds": config.pressure_escalate_seconds,
        "near_limit_fraction": config.near_limit_fraction,
        "auto_default": config.auto_default,
        "auto_factor": config.auto_factor,
        "auto_history": config.auto_history,
        "auto_round": config.auto_round,
        "warn_factor": config.warn_factor,
        "warn_min_excess": config.warn_min_excess,
        "test_node_total": config.test_node_total,
    }


def _parse_borrow(value: Any) -> BorrowConfig:
    raw = _object(value, "borrow")
    _unknown(
        raw,
        {"enabled", "slice", "slice_cpu_weight", "unit_cpu_weight", "idle_window_s", "busy_window_s",
         "idle_threshold", "cpu_idle_threshold", "busy_threshold", "freeze_max_s", "default_on_preempt",
         "backfill_waiting_owners", "notes"},
        "borrow",
    )
    defaults = BorrowConfig()
    enabled = raw.get("enabled", defaults.enabled)
    if not isinstance(enabled, bool):
        raise ConfigError("borrow.enabled must be true or false")
    slice_name = raw.get("slice", defaults.slice)
    # A dash nests a slice under its prefix ("simq-borrow.slice" is a child of
    # simq.slice).  The borrow slice must sit directly under the user manager
    # so that its idle weight applies against every other user slice.
    if (not isinstance(slice_name, str) or not slice_name.endswith(".slice") or "-" in slice_name
            or re.fullmatch(r"[A-Za-z0-9_.]+\.slice", slice_name) is None):
        raise ConfigError("borrow.slice must name a top-level slice without '-', e.g. simqborrow.slice")
    slice_weight = raw.get("slice_cpu_weight", defaults.slice_cpu_weight)
    if slice_weight != "idle" and not (isinstance(slice_weight, int) and not isinstance(slice_weight, bool)
                                       and 1 <= slice_weight <= 10000):
        raise ConfigError("borrow.slice_cpu_weight must be 'idle' or 1..10000")
    unit_weight = raw.get("unit_cpu_weight", defaults.unit_cpu_weight)
    if isinstance(unit_weight, bool) or not isinstance(unit_weight, int) or not 1 <= unit_weight <= 10000:
        raise ConfigError("borrow.unit_cpu_weight must be 1..10000")

    def number(key: str, low: float, high: float) -> float:
        item = raw.get(key, getattr(defaults, key))
        if isinstance(item, bool) or not isinstance(item, (int, float)) or not low <= float(item) <= high:
            raise ConfigError(f"borrow.{key} must be a number in [{low}, {high}]")
        return float(item)

    freeze_max = raw.get("freeze_max_s", defaults.freeze_max_s)
    if freeze_max is not None and (isinstance(freeze_max, bool) or not isinstance(freeze_max, (int, float))
                                   or freeze_max <= 0):
        raise ConfigError("borrow.freeze_max_s must be a positive number or null")
    mode = raw.get("default_on_preempt", defaults.default_on_preempt)
    if mode not in ON_PREEMPT_MODES:
        raise ConfigError(f"borrow.default_on_preempt must be one of {', '.join(ON_PREEMPT_MODES)}")
    backfill = raw.get("backfill_waiting_owners", defaults.backfill_waiting_owners)
    if not isinstance(backfill, bool):
        raise ConfigError("borrow.backfill_waiting_owners must be true or false")
    return BorrowConfig(
        enabled=enabled,
        slice=slice_name,
        slice_cpu_weight=str(slice_weight),
        unit_cpu_weight=unit_weight,
        idle_window_s=number("idle_window_s", 0.5, 3600),
        busy_window_s=number("busy_window_s", 0.2, 600),
        idle_threshold=number("idle_threshold", 0.0, 1.0),
        cpu_idle_threshold=number("cpu_idle_threshold", 0.0, 1.0),
        busy_threshold=number("busy_threshold", 0.0, 1.0),
        freeze_max_s=None if freeze_max is None else float(freeze_max),
        default_on_preempt=mode,
        backfill_waiting_owners=backfill,
    )


def load_node_config(path: Path | None = None) -> NodeConfig | None:
    """Return the node configuration, or ``None`` for 0.5 legacy mode."""
    target = path or config_path()
    try:
        raw = target.read_text(encoding="utf-8")
    except FileNotFoundError:
        if os.environ.get("SIMQUEUE_CONFIG"):
            raise ConfigError(f"SIMQUEUE_CONFIG names a missing file: {target}") from None
        return None
    except OSError as error:
        raise ConfigError(f"cannot read node config {target}: {error}") from error
    try:
        document = json.loads(raw)
    except json.JSONDecodeError as error:
        raise ConfigError(f"invalid JSON in node config {target}: {error}") from error
    return parse_node_config(document, path=target)


def describe_partition(partition: Partition) -> dict[str, Any]:
    return {
        "name": partition.name,
        "cpus": format_cpu_list(partition.cpus),
        "gpus": list(partition.gpus),
        "default_cpus": partition.resolve_default_cpus(),
        "default_mem": partition.default_mem,
        "max_mem": partition.max_mem,
        "min_mem": partition.min_mem,
        "default_time_s": partition.default_time,
        "max_time_s": partition.max_time,
        "swap_max": partition.swap_max,
        "enforce_time_limit": partition.enforce_time_limit,
        "nice": partition.nice,
        "description": partition.description,
        "lend": partition.lend,
        "lend_leased_idle": partition.lend_leased_idle,
        "timing_sensitive": partition.timing_sensitive,
        "aliases": list(partition.aliases),
    }


def describe_alias(config: NodeConfig, name: str) -> dict[str, Any]:
    """An alias as ``simq partitions --json`` shows it: its effective defaults."""
    alias = config.aliases[name]
    view = config.partition(name)
    return {
        "name": name,
        "partition": alias.partition,
        "description": alias.description,
        "default_cpus": view.resolve_default_cpus(),
        "default_mem": view.default_mem,
        "default_mem_auto": view.mem_auto,
        "default_time_s": view.default_time,
        "prefer_cpus": list(alias.prefer_cpus),
    }
