from __future__ import annotations

import os
from pathlib import Path


def project_root() -> Path:
    configured = os.environ.get("SIMQUEUE_PROJECT_ROOT", "~/projects")
    return Path(configured).expanduser().resolve()


def state_dir() -> Path:
    configured = os.environ.get("SIMQUEUE_STATE_DIR")
    if configured:
        return Path(configured).expanduser().resolve()
    base = Path(os.environ.get("XDG_STATE_HOME", "~/.local/state")).expanduser()
    return (base / "simqueue").resolve()


def database_path() -> Path:
    return state_dir() / "queue.sqlite3"


def log_dir() -> Path:
    return state_dir() / "logs"


def parse_cpu_list(value: str) -> tuple[int, ...]:
    """Parse Linux cpulist syntax such as ``0-3,8,10-11``."""
    cpus: set[int] = set()
    if not value.strip():
        return ()
    for raw_part in value.split(","):
        part = raw_part.strip()
        if not part:
            raise ValueError(f"invalid empty CPU-list component in {value!r}")
        if "-" in part:
            bounds = part.split("-", 1)
            if len(bounds) != 2:
                raise ValueError(f"invalid CPU range {part!r}")
            start, end = (int(item) for item in bounds)
            if start < 0 or end < start:
                raise ValueError(f"invalid CPU range {part!r}")
            cpus.update(range(start, end + 1))
        else:
            cpu = int(part)
            if cpu < 0:
                raise ValueError(f"invalid CPU {part!r}")
            cpus.add(cpu)
    return tuple(sorted(cpus))


def format_cpu_list(cpus: tuple[int, ...] | list[int] | set[int]) -> str:
    ordered = sorted(set(cpus))
    if not ordered:
        return ""
    ranges: list[str] = []
    start = previous = ordered[0]
    for cpu in ordered[1:]:
        if cpu == previous + 1:
            previous = cpu
            continue
        ranges.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = cpu
    ranges.append(str(start) if start == previous else f"{start}-{previous}")
    return ",".join(ranges)


def _read_int(path: Path) -> int | None:
    try:
        return int(path.read_text(encoding="utf-8").strip())
    except (OSError, ValueError):
        return None


def _allowed_cpus() -> tuple[int, ...]:
    try:
        return tuple(sorted(os.sched_getaffinity(0)))
    except AttributeError:
        return tuple(range(os.cpu_count() or 1))


def _physical_representatives(cpus: tuple[int, ...]) -> tuple[int, ...]:
    representatives: dict[tuple[int, int], int] = {}
    for cpu in cpus:
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        package = _read_int(topology / "physical_package_id")
        core = _read_int(topology / "core_id")
        if package is None or core is None:
            return cpus
        key = (package, core)
        representatives[key] = min(cpu, representatives.get(key, cpu))
    return tuple(sorted(representatives.values()))


def _capacity_pcores(cpus: tuple[int, ...]) -> tuple[int, ...]:
    capacities = {
        cpu: _read_int(Path(f"/sys/devices/system/cpu/cpu{cpu}/cpu_capacity"))
        for cpu in cpus
    }
    values = {value for value in capacities.values() if value is not None}
    if len(values) < 2:
        return ()
    ordered = sorted(values)
    gap, split_index = max(
        (upper - lower, index)
        for index, (lower, upper) in enumerate(zip(ordered, ordered[1:]))
    )
    if gap < ordered[-1] * 0.05:
        return ()
    threshold = ordered[split_index + 1]
    return tuple(
        cpu
        for cpu in cpus
        if capacities[cpu] is not None and int(capacities[cpu]) >= threshold
    )


def detect_pcores() -> tuple[int, ...]:
    """Return the P-core CPU IDs available to this process.

    An explicit ``SIMQUEUE_PCORES`` setting wins. Linux's hybrid-core sysfs
    interface is preferred, followed by scheduler capacity data. The final
    fallback is one logical CPU per physical core.
    """
    allowed = set(_allowed_cpus())
    configured = os.environ.get("SIMQUEUE_PCORES")
    if configured is not None:
        detected = parse_cpu_list(configured)
    else:
        hybrid_path = Path("/sys/devices/cpu_core/cpus")
        try:
            detected = parse_cpu_list(hybrid_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            detected = _capacity_pcores(tuple(sorted(allowed)))
            if not detected:
                detected = _physical_representatives(tuple(sorted(allowed)))

    available = tuple(cpu for cpu in detected if cpu in allowed)
    if not available:
        raise RuntimeError(
            "no configured P-cores are available in this process's CPU affinity"
        )
    return available


def ensure_state_directories() -> None:
    for directory in (state_dir(), log_dir()):
        directory.mkdir(mode=0o700, parents=True, exist_ok=True)
        try:
            directory.chmod(0o700)
        except OSError:
            pass
