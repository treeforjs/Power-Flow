"""CPU topology: L3 cache domains and SMT siblings, read from sysfs (0.8.0).

Nothing is hard-coded.  At start the worker (and every CLI that loads the
node configuration) reads, for each online CPU,
``cpu*/cache/index*/{level,type,shared_cpu_list,size,id}`` (the level-3
cache) and ``cpu*/topology/thread_siblings_list``.  CPUs that share one L3
form a *domain*: a Ryzen 7 9850X3D has one (16 logical CPUs, 96 MiB), a
Ryzen 9 9950X3D2 two (one per chiplet, each with stacked cache), a 64-core
station as many as it has CCDs.  A machine without L3 information is one
domain.

The scheduler works on physical cores named by their lowest logical CPU (a
partition's ``cpus``); :func:`core_domains` maps those onto the domains.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .config import format_cpu_list, parse_cpu_list

SYSFS_CPU = Path("/sys/devices/system/cpu")
SYSFS_NODE = Path("/sys/devices/system/node")


@dataclass(frozen=True)
class NumaNode:
    node: int  # the kernel's node number
    cpus: tuple[int, ...]


def discover_numa(root: Path = SYSFS_NODE, online: tuple[int, ...] = ()) -> tuple[NumaNode, ...]:
    """NUMA nodes with CPUs (``node*/cpulist``; memory-only nodes skipped).
    NPS1/2/4 on EPYC gives 1/2/4 nodes per socket; without the directory the
    machine is one node."""
    nodes = []
    try:
        entries = sorted(root.glob("node[0-9]*"), key=lambda p: int(p.name[4:]))
    except (OSError, ValueError):
        entries = []
    for entry in entries:
        text = _read(entry / "cpulist")
        if not text:
            continue
        try:
            cpus = parse_cpu_list(text)
        except ValueError:
            continue
        if cpus:
            nodes.append(NumaNode(int(entry.name[4:]), cpus))
    return tuple(nodes) or (NumaNode(0, tuple(online)),)


def node_meminfo(root: Path = SYSFS_NODE) -> dict[int, dict[str, int]]:
    """Per NUMA node: MemTotal, MemFree, FilePages (bytes)."""
    result: dict[int, dict[str, int]] = {}
    for entry in root.glob("node[0-9]*"):
        text = _read(entry / "meminfo")
        if not text:
            continue
        values: dict[str, int] = {}
        for line in text.splitlines():
            parts = line.split()
            if len(parts) >= 4 and parts[2].rstrip(":") in ("MemTotal", "MemFree", "FilePages"):
                values[parts[2].rstrip(":")] = int(parts[3]) * 1024
        result[int(entry.name[4:])] = values
    return result


@dataclass(frozen=True)
class L3Domain:
    index: int
    cpus: tuple[int, ...]  # every logical CPU sharing this L3 (SMT siblings included)
    size_bytes: int | None = None
    cache_id: int | None = None

    def describe(self) -> str:
        size = f", {self.size_bytes // 1024**2} MiB L3" if self.size_bytes else ""
        return f"L3 domain {self.index}: CPUs {format_cpu_list(self.cpus)}{size}"


def _read(path: Path) -> str | None:
    try:
        return path.read_text(encoding="ascii").strip()
    except OSError:
        return None


def _size(text: str | None) -> int | None:
    if not text:
        return None
    units = {"K": 1024, "M": 1024**2, "G": 1024**3}
    try:
        if text[-1].upper() in units:
            return int(text[:-1]) * units[text[-1].upper()]
        return int(text)
    except ValueError:
        return None


def online_cpus(root: Path = SYSFS_CPU) -> tuple[int, ...]:
    text = _read(root / "online")
    if text:
        try:
            return parse_cpu_list(text)
        except ValueError:
            pass
    found = []
    for entry in root.glob("cpu[0-9]*"):
        try:
            found.append(int(entry.name[3:]))
        except ValueError:
            continue
    return tuple(sorted(found))


def l3_of(cpu: int, root: Path = SYSFS_CPU) -> tuple[tuple[int, ...], int | None, int | None] | None:
    """(shared CPUs, size, id) of the level-3 cache of *cpu*, if sysfs has it."""
    cache = root / f"cpu{cpu}" / "cache"
    try:
        entries = sorted(cache.glob("index*"))
    except OSError:
        return None
    for entry in entries:
        if _read(entry / "level") != "3" or _read(entry / "type") not in ("Unified", "Data"):
            continue
        shared = _read(entry / "shared_cpu_list")
        if not shared:
            continue
        try:
            cpus = parse_cpu_list(shared)
        except ValueError:
            continue
        cache_id = _read(entry / "id")
        return cpus, _size(_read(entry / "size")), int(cache_id) if cache_id and cache_id.isdigit() else None
    return None


def discover_l3(root: Path = SYSFS_CPU) -> tuple[L3Domain, ...]:
    """Every L3 domain of the machine, ordered by its lowest CPU."""
    online = online_cpus(root)
    groups: dict[tuple[int, ...], tuple[int | None, int | None]] = {}
    unknown: list[int] = []
    for cpu in online:
        found = l3_of(cpu, root)
        if found is None:
            unknown.append(cpu)
            continue
        shared, size, cache_id = found
        groups.setdefault(tuple(c for c in shared if c in online) or shared, (size, cache_id))
    covered = {c for group in groups for c in group}
    unknown = [c for c in unknown if c not in covered]
    if unknown:
        groups[tuple(unknown)] = (None, None)
    if not groups:
        return (L3Domain(0, online),)
    ordered = sorted(groups.items(), key=lambda item: min(item[0]))
    return tuple(L3Domain(index, cpus, size, cache_id) for index, (cpus, (size, cache_id)) in enumerate(ordered))


def thread_siblings(cpus: Iterable[int], root: Path = SYSFS_CPU) -> dict[int, tuple[int, ...]]:
    """SMT siblings of each CPU (``thread_siblings_list``), the CPU itself excluded."""
    result: dict[int, tuple[int, ...]] = {}
    for cpu in cpus:
        text = _read(root / f"cpu{cpu}" / "topology" / "thread_siblings_list")
        if text is None:
            continue
        try:
            result[cpu] = tuple(c for c in parse_cpu_list(text) if c != cpu)
        except ValueError:
            continue
    return result


def physical_cores(root: Path = SYSFS_CPU) -> tuple[int, ...]:
    """Every physical core of the machine, named by its lowest logical CPU."""
    online = online_cpus(root)
    siblings = thread_siblings(online, root)
    return tuple(sorted({min((cpu, *siblings.get(cpu, ()))) for cpu in online}))


def core_domains(cores: Iterable[int], domains: tuple[L3Domain, ...]) -> tuple[tuple[int, ...], ...]:
    """The scheduler's cores (lowest logical CPU of each physical core)
    grouped by L3 domain, in domain order; cores outside every domain form a
    last group.  Empty groups are dropped."""
    cores = sorted(set(cores))
    groups = [tuple(c for c in cores if c in domain.cpus) for domain in domains]
    seen = {c for group in groups for c in group}
    rest = tuple(c for c in cores if c not in seen)
    return tuple(group for group in [*groups, rest] if group)


def sysfs_root(configured: str | None = None) -> Path:
    """The sysfs CPU directory: a configured one (tests use fake trees) or
    ``$SIMQUEUE_SYSFS_CPU`` or the real one."""
    return Path(configured or os.environ.get("SIMQUEUE_SYSFS_CPU") or SYSFS_CPU)
