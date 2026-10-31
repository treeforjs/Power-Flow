"""Node and partition configuration for the SimQueue 0.6 node scheduler.

A node configuration file switches SimQueue from the 0.5 serial queue
(``legacy`` mode, driven by ``SIMQUEUE_PCORES``) to the concurrent node
scheduler.  The file is JSON with schema ``simqueue-node-config-v1``.  Its
path is ``$SIMQUEUE_CONFIG`` when set, otherwise ``<state dir>/node.json``.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .config import format_cpu_list, parse_cpu_list, state_dir


CONFIG_SCHEMA = "simqueue-node-config-v1"
PARTITION_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
GPU_UUID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$")
DEFAULT_WRAPPER_PYTHON = "/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen"
UNLIMITED_WORDS = {"unlimited", "infinite", "none"}
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

    def resolve_default_cpus(self) -> int:
        return len(self.cpus) if self.default_cpus is None else self.default_cpus


@dataclass(frozen=True)
class LaunchConfig:
    backend: str = "systemd-user"
    unit_prefix: str = "simq-job"
    slice: str = "simq.slice"
    systemd_run: str = "systemd-run"
    systemctl: str = "systemctl"
    wrapper_python: str = DEFAULT_WRAPPER_PYTHON
    kill_grace_seconds: float = 10.0
    poll_seconds: float = 1.0


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

    def partition(self, name: str | None) -> Partition:
        key = name or self.default_partition
        try:
            return self.partitions[key]
        except KeyError:
            known = ", ".join(sorted(self.partitions))
            raise ValueError(f"unknown partition {key!r} (configured: {known})") from None

    def node_cpus(self) -> tuple[int, ...]:
        cpus: set[int] = set()
        for partition in self.partitions.values():
            cpus.update(partition.cpus)
        return tuple(sorted(cpus))

    def partitions_for_cpu(self, cpu: int) -> list[str]:
        return [name for name, part in self.partitions.items() if cpu in part.cpus]


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
    _unknown(lock_raw, {"enabled", "core_path", "core_overrides", "gpu_path", "gpu_extra", "notes"}, "lock_files")
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
    for key in ("core_path", "gpu_path"):
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
    )

    launch_raw = _object(document.get("launch", {}), "launch")
    _unknown(
        launch_raw,
        {"backend", "unit_prefix", "slice", "systemd_run", "systemctl", "wrapper_python",
         "kill_grace_seconds", "poll_seconds"},
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
    except (TypeError, ValueError) as error:
        raise ConfigError("launch timing values must be numbers") from error
    if not 0.1 <= grace <= 3600 or not 0.05 <= poll <= 60:
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
    online = online_cpus() if check_topology else None
    for index, raw in enumerate(raw_partitions):
        item = _object(raw, f"partitions[{index}]")
        _unknown(
            item,
            {"name", "cpus", "cores", "gpus", "default_mem", "max_mem", "default_time", "max_time",
             "default_cpus", "swap_max", "enforce_time_limit", "nice", "description"},
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
            default_time = parse_time_limit(item.get("default_time"), field_name=f"partition {name} default_time")
            max_time = parse_time_limit(item.get("max_time"), field_name=f"partition {name} max_time")
        except ValueError as error:
            raise ConfigError(str(error)) from error
        if default_mem <= 0 or max_mem < default_mem:
            raise ConfigError(f"partition {name}: need 0 < default_mem <= max_mem")
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
        )
    default_partition = document.get("default_partition")
    if default_partition not in partitions:
        raise ConfigError("default_partition must name a configured partition")
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
        "default_time_s": partition.default_time,
        "max_time_s": partition.max_time,
        "swap_max": partition.swap_max,
        "enforce_time_limit": partition.enforce_time_limit,
        "nice": partition.nice,
        "description": partition.description,
    }
