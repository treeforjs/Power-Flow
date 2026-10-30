from __future__ import annotations

import csv
import ctypes
import json
import math
import os
import re
import shlex
import shutil
import signal
import stat
import subprocess
import time
from collections import Counter, deque
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from types import FrameType
from typing import Any, BinaryIO

from .config import parse_cpu_list, project_root
from .store import (
    ENSEMBLE_CONTROL_PROTOCOL,
    ENSEMBLE_DRAIN_EXIT_CODE,
    ENSEMBLE_DRAIN_PROTOCOL_ENV,
    MEMBER_TERMINAL_STATES,
    QueueStore,
    utc_now,
)


SCHEMA_VERSION = "simqueue-ensemble-v1"
POLL_SECONDS = 0.1
FAIL_FAST_GRACE_SECONDS = 5.0
CLEANUP_KILL_SETTLE_SECONDS = 5.0
PR_SET_CHILD_SUBREAPER = 36
SQLITE_INTEGER_MAX = 2**63 - 1
NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
ENV_NAME_PATTERN = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
SIZE_PATTERN = re.compile(r"^([0-9]+)\s*([A-Za-z]*)$")
BYTE_UNITS = {
    "": 1,
    "b": 1,
    "kb": 1000,
    "mb": 1000**2,
    "gb": 1000**3,
    "tb": 1000**4,
    "kib": 1024,
    "mib": 1024**2,
    "gib": 1024**3,
    "tib": 1024**4,
}
THREAD_LIMIT_ENVIRONMENT = {
    "OMP_NUM_THREADS": "1",
    "OMP_THREAD_LIMIT": "1",
    "OMP_DYNAMIC": "FALSE",
    "OMP_PROC_BIND": "FALSE",
    "KOKKOS_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "BLIS_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}
# MPI launchers commonly discover all host CPUs and replace inherited
# sched_setaffinity masks. Disable their own binding so the enclosing
# taskset lane remains authoritative for rank-1 members. PRTE covers Open MPI
# 5, while OMPI covers older Open MPI releases.
MPI_AFFINITY_ENVIRONMENT = {
    "OMPI_MCA_hwloc_base_binding_policy": "none",
    "PRTE_MCA_bindto": "none",
    "I_MPI_PIN": "0",
    "MV2_ENABLE_AFFINITY": "0",
}
AFFINITY_ENVIRONMENT_TO_CLEAR = {
    "GOMP_CPU_AFFINITY",
    "KMP_AFFINITY",
    "KMP_HW_SUBSET",
    "KMP_PLACE_THREADS",
    "LIKWID_PIN",
    "MV2_CPU_BINDING_POLICY",
    "MV2_ENABLE_AFFINITY",
    "OMP_DISPLAY_AFFINITY",
    "OMP_PLACES",
    "OMPI_MCA_hwloc_base_binding_policy",
    "OMPI_MCA_rmaps_base_mapping_policy",
    "PRTE_MCA_bind_to",
    "PRTE_MCA_bindto",
    "PRTE_MCA_hwloc_default_binding_policy",
    "SLURM_CPU_BIND",
    "SLURM_CPU_BIND_LIST",
    "SUNW_MP_PROCBIND",
    "I_MPI_PIN",
    "I_MPI_PIN_PROCESSOR_LIST",
}
DIAGNOSTIC_ENVIRONMENT_NAMES = {
    "CUDA_VISIBLE_DEVICES",
    "CUDA_LAUNCH_BLOCKING",
    "CUDA_MANAGED_FORCE_DEVICE_ALLOC",
    "CUDA_MPS_PIPE_DIRECTORY",
    "CUDA_MPS_ACTIVE_THREAD_PERCENTAGE",
    "CUDA_MPS_PINNED_DEVICE_MEM_LIMIT",
    "OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT",
    "OMP_DYNAMIC",
    "OMP_PLACES",
    "OMP_PROC_BIND",
    "OMPI_MCA_hwloc_base_binding_policy",
    "PRTE_MCA_bindto",
    "I_MPI_PIN",
    "MV2_ENABLE_AFFINITY",
    "KOKKOS_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "TPETRA_ASSUME_CUDA_AWARE_MPI",
}
GPU_QUERY_FIELDS = (
    "timestamp",
    "index",
    "uuid",
    "memory.total",
    "memory.used",
    "memory.free",
    "utilization.gpu",
    "utilization.memory",
    "pstate",
    "temperature.gpu",
    "power.draw",
)
GPU_INTEGER_FIELDS = {
    "index",
    "memory.total",
    "memory.used",
    "memory.free",
    "utilization.gpu",
    "utilization.memory",
    "temperature.gpu",
}
UVM_INTEGER_FIELDS = {
    "limit_bytes",
    "live_allocations",
    "live_bytes",
    "peak_live_bytes",
    "live_reserved_bytes",
    "peak_reserved_bytes",
    "allocation_requests",
    "raw_allocations",
    "pool_reclaims",
    "teardown_ok",
}


def parse_byte_size(value: Any, *, field: str = "size") -> int:
    if isinstance(value, bool):
        raise ValueError(f"{field} must be a nonnegative byte size")
    if isinstance(value, int):
        if value < 0:
            raise ValueError(f"{field} must be a nonnegative byte size")
        if value > SQLITE_INTEGER_MAX:
            raise ValueError(
                f"{field} exceeds the maximum supported size "
                f"({SQLITE_INTEGER_MAX} bytes)"
            )
        return value
    if not isinstance(value, str):
        raise ValueError(f"{field} must be an integer or a byte-size string")
    match = SIZE_PATTERN.fullmatch(value.strip())
    if match is None:
        raise ValueError(
            f"invalid {field} size {value!r}; use bytes or B/KiB/MiB/GiB"
        )
    unit = match.group(2).lower()
    multiplier = BYTE_UNITS.get(unit)
    if multiplier is None:
        raise ValueError(
            f"invalid {field} size unit {match.group(2)!r}; use B/KiB/MiB/GiB"
        )
    size = int(match.group(1)) * multiplier
    if size > SQLITE_INTEGER_MAX:
        raise ValueError(
            f"{field} exceeds the maximum supported size "
            f"({SQLITE_INTEGER_MAX} bytes)"
        )
    return size


def format_byte_size(value: int) -> str:
    for suffix, multiplier in (
        ("TiB", 1024**4),
        ("GiB", 1024**3),
        ("MiB", 1024**2),
        ("KiB", 1024),
    ):
        if value >= multiplier and value % multiplier == 0:
            return f"{value // multiplier}{suffix}"
    return f"{value}B"


def _elapsed_seconds(started_at: str | None, finished_at: str | None) -> float | None:
    if not started_at:
        return None
    try:
        started = datetime.fromisoformat(started_at.replace("Z", "+00:00"))
        finished = (
            datetime.fromisoformat(finished_at.replace("Z", "+00:00"))
            if finished_at
            else datetime.now(UTC)
        )
    except ValueError:
        return None
    return round(max(0.0, (finished - started).total_seconds()), 3)


def _object(value: Any, *, field: str) -> dict[str, Any]:
    if not isinstance(value, dict) or not all(isinstance(key, str) for key in value):
        raise ValueError(f"{field} must be a JSON object")
    return dict(value)


def _reject_unknown(
    document: dict[str, Any], allowed: set[str], *, field: str
) -> None:
    unknown = sorted(set(document) - allowed)
    if unknown:
        raise ValueError(f"unknown {field} field(s): {', '.join(unknown)}")


def _environment(value: Any, *, field: str) -> dict[str, str]:
    document = _object(value, field=field)
    environment: dict[str, str] = {}
    for key, raw_value in document.items():
        if ENV_NAME_PATTERN.fullmatch(key) is None or key.startswith("SIMQ_"):
            raise ValueError(f"invalid or reserved environment name {key!r} in {field}")
        if not isinstance(raw_value, str) or "\x00" in raw_value:
            raise ValueError(f"{field}.{key} must be a string without NUL bytes")
        environment[key] = raw_value
    return environment


def _canonical_launch_gate(
    value: Any, *, field: str = "launch_gate"
) -> dict[str, Any] | None:
    if value is None:
        return None
    gate = _object(value, field=field)
    _reject_unknown(
        gate,
        {"kind", "contains", "timeout_seconds"},
        field=field,
    )
    kind = gate.get("kind")
    if kind != "member-log-substring":
        raise ValueError(f"{field}.kind must be 'member-log-substring'")
    contains = gate.get("contains")
    if (
        not isinstance(contains, str)
        or not contains
        or any(character in contains for character in ("\x00", "\n", "\r"))
        or len(contains.encode("utf-8")) > 1024
    ):
        raise ValueError(
            f"{field}.contains must be a nonempty, single-line UTF-8 "
            "substring of at most 1024 bytes"
        )
    timeout_seconds = gate.get("timeout_seconds", 300.0)
    if (
        isinstance(timeout_seconds, bool)
        or not isinstance(timeout_seconds, (int, float))
        or not math.isfinite(float(timeout_seconds))
        or float(timeout_seconds) <= 0
    ):
        raise ValueError(f"{field}.timeout_seconds must be a positive number")
    return {
        "kind": kind,
        "contains": contains,
        "timeout_seconds": float(timeout_seconds),
    }


def _member_cwd(value: Any, *, manifest_directory: Path, field: str) -> Path:
    if not isinstance(value, str) or not value.strip() or "\x00" in value:
        raise ValueError(f"{field} must be a nonempty path string")
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = manifest_directory / path
    try:
        path = path.resolve(strict=True)
    except OSError as error:
        raise ValueError(f"{field} is unavailable: {path}: {error}") from error
    if not path.is_dir():
        raise ValueError(f"{field} is not a directory: {path}")
    root = project_root()
    if not path.is_relative_to(root):
        raise ValueError(f"{field} must resolve inside {root}: {path}")
    return path


def _member_argv(value: Any, *, field: str) -> list[str]:
    if (
        not isinstance(value, list)
        or not value
        or not all(isinstance(argument, str) and "\x00" not in argument for argument in value)
        or not value[0]
    ):
        raise ValueError(
            f"{field} must be a nonempty string array with a nonempty executable "
            "and no NUL bytes"
        )
    return list(value)


def _uvm_arena(environment: dict[str, str], *, member_name: str) -> int:
    backing = environment.get("MSIC_KOKKOS_CUDAUVM_BACKING")
    raw_limit = environment.get("MSIC_KOKKOS_CUDAUVM_LIMIT_BYTES")
    if backing is None:
        if raw_limit is not None:
            raise ValueError(
                f"member {member_name!r} sets a UVM limit without "
                "MSIC_KOKKOS_CUDAUVM_BACKING"
            )
        return 0
    if backing not in ("legacy", "managed-pool", "mapped-host"):
        raise ValueError(
            f"member {member_name!r} has unsupported UVM backing {backing!r}"
        )
    if backing == "legacy":
        if raw_limit is not None:
            raise ValueError(
                f"member {member_name!r} uses legacy UVM backing with a forbidden limit"
            )
        return 0
    if raw_limit is None or re.fullmatch(r"[1-9][0-9]*", raw_limit) is None:
        raise ValueError(
            f"member {member_name!r} UVM limit must be a canonical positive decimal integer"
        )
    limit = int(raw_limit)
    if limit > SQLITE_INTEGER_MAX:
        raise ValueError(
            f"member {member_name!r} UVM limit exceeds the maximum supported "
            f"size ({SQLITE_INTEGER_MAX} bytes)"
        )
    if backing == "managed-pool" and limit % 256:
        raise ValueError(
            f"member {member_name!r} managed-pool UVM limit must be a multiple of 256"
        )
    return limit


def load_manifest(
    path: Path,
    *,
    available_cores: int,
    name_override: str | None = None,
    max_parallel_override: int | None = None,
    gpu_budget_override: Any | None = None,
    gpu_device_override: str | None = None,
    diagnostics_interval_override: float | None = None,
    fail_fast_override: bool | None = None,
) -> dict[str, Any]:
    """Load, validate, and canonicalize an immutable ensemble manifest."""
    manifest_path = path.expanduser().resolve()
    try:
        document = _object(
            json.loads(manifest_path.read_text(encoding="utf-8")),
            field="manifest",
        )
    except OSError as error:
        raise ValueError(f"cannot read ensemble manifest {manifest_path}: {error}") from error
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid JSON in ensemble manifest {manifest_path}: {error}") from error
    _reject_unknown(
        document,
        {
            "schema_version",
            "name",
            "max_parallel",
            "gpu",
            "defaults",
            "members",
            "diagnostics_interval",
            "fail_fast",
            "allow_shared_cwd",
            "launch_gate",
        },
        field="manifest",
    )
    if document.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            f"unsupported ensemble schema {document.get('schema_version')!r}; "
            f"expected {SCHEMA_VERSION!r}"
        )

    raw_name = name_override if name_override is not None else document.get("name")
    if not isinstance(raw_name, str) or NAME_PATTERN.fullmatch(raw_name) is None:
        raise ValueError(
            "ensemble name must start with an alphanumeric and contain only "
            "letters, digits, dot, underscore, or dash"
        )
    defaults = _object(document.get("defaults", {}), field="defaults")
    _reject_unknown(defaults, {"cwd", "env", "gpu_memory"}, field="defaults")
    default_environment = _environment(defaults.get("env", {}), field="defaults.env")

    raw_members = document.get("members")
    if not isinstance(raw_members, list) or not raw_members:
        raise ValueError("manifest.members must be a nonempty array")
    manifest_directory = manifest_path.parent
    members: list[dict[str, Any]] = []
    names: set[str] = set()
    working_directories: dict[Path, str] = {}
    allow_shared_cwd = document.get("allow_shared_cwd", False)
    if not isinstance(allow_shared_cwd, bool):
        raise ValueError("allow_shared_cwd must be true or false")
    for index, raw_member in enumerate(raw_members):
        member = _object(raw_member, field=f"members[{index}]")
        _reject_unknown(
            member,
            {"name", "cwd", "argv", "env", "gpu_memory"},
            field=f"members[{index}]",
        )
        member_name = member.get("name")
        if (
            not isinstance(member_name, str)
            or NAME_PATTERN.fullmatch(member_name) is None
            or member_name.isdecimal()
        ):
            raise ValueError(
                f"members[{index}].name must be a nonnumeric safe identifier"
            )
        if member_name in names:
            raise ValueError(f"duplicate ensemble member name {member_name!r}")
        names.add(member_name)
        cwd_value = member.get("cwd", defaults.get("cwd", "."))
        cwd = _member_cwd(
            cwd_value,
            manifest_directory=manifest_directory,
            field=f"members[{index}].cwd",
        )
        if not allow_shared_cwd and cwd in working_directories:
            raise ValueError(
                f"members {working_directories[cwd]!r} and {member_name!r} share "
                f"working directory {cwd}; set allow_shared_cwd only when outputs are isolated"
            )
        working_directories[cwd] = member_name
        environment = dict(default_environment)
        environment.update(
            _environment(member.get("env", {}), field=f"members[{index}].env")
        )
        claim = parse_byte_size(
            member.get("gpu_memory", defaults.get("gpu_memory", 0)),
            field=f"members[{index}].gpu_memory",
        )
        uvm_arena = _uvm_arena(environment, member_name=member_name)
        if environment.get("MSIC_KOKKOS_CUDAUVM_BACKING") == "managed-pool" and uvm_arena > claim:
            raise ValueError(
                f"member {member_name!r} UVM arena ({uvm_arena} bytes) exceeds "
                f"its gpu_memory admission claim ({claim} bytes)"
            )
        members.append(
            {
                "name": member_name,
                "cwd": str(cwd),
                "argv": _member_argv(
                    member.get("argv"), field=f"members[{index}].argv"
                ),
                "env": environment,
                "gpu_admission_bytes": claim,
                "uvm_arena_bytes": uvm_arena,
            }
        )

    if available_cores < 1:
        raise ValueError("no P-cores are available for ensemble lanes")
    configured_parallel = (
        max_parallel_override
        if max_parallel_override is not None
        else document.get("max_parallel", min(len(members), available_cores))
    )
    if isinstance(configured_parallel, bool) or not isinstance(configured_parallel, int):
        raise ValueError("max_parallel must be a positive integer")
    if configured_parallel < 1 or configured_parallel > available_cores:
        raise ValueError(
            f"max_parallel must be between 1 and {available_cores} on this host"
        )
    max_parallel = min(configured_parallel, len(members))

    gpu = _object(document.get("gpu", {}), field="gpu")
    _reject_unknown(gpu, {"device", "budget", "mps"}, field="gpu")
    budget_value = gpu_budget_override if gpu_budget_override is not None else gpu.get("budget")
    any_gpu_claim = any(member["gpu_admission_bytes"] for member in members)
    if budget_value is None:
        if any_gpu_claim:
            raise ValueError(
                "gpu.budget (or --gpu-budget) is required when members claim GPU memory"
            )
        gpu_budget_bytes = 0
    else:
        gpu_budget_bytes = parse_byte_size(budget_value, field="gpu budget")
    for member in members:
        if member["gpu_admission_bytes"] > gpu_budget_bytes:
            raise ValueError(
                f"member {member['name']!r} gpu_memory claim exceeds the GPU budget"
            )
    gpu_device = gpu_device_override if gpu_device_override is not None else gpu.get("device", "0")
    if not isinstance(gpu_device, str) or not gpu_device or "\x00" in gpu_device:
        raise ValueError("gpu.device must be a nonempty string")
    mps_mode = gpu.get("mps", "off")
    if mps_mode != "off":
        raise ValueError(
            "gpu.mps currently supports only 'off'; external MPS is disabled "
            "until SimQueue can guarantee controller-mediated client teardown"
        )
    mps_pipe = None

    diagnostics_interval = (
        diagnostics_interval_override
        if diagnostics_interval_override is not None
        else document.get("diagnostics_interval", 1.0)
    )
    if (
        isinstance(diagnostics_interval, bool)
        or not isinstance(diagnostics_interval, (int, float))
        or not math.isfinite(float(diagnostics_interval))
        or float(diagnostics_interval) <= 0
    ):
        raise ValueError("diagnostics_interval must be a positive number of seconds")
    fail_fast = (
        fail_fast_override
        if fail_fast_override is not None
        else document.get("fail_fast", False)
    )
    if not isinstance(fail_fast, bool):
        raise ValueError("fail_fast must be true or false")

    launch_gate = _canonical_launch_gate(document.get("launch_gate"))

    canonical = {
        "schema_version": SCHEMA_VERSION,
        "name": raw_name,
        "max_parallel": max_parallel,
        "gpu": {
            "device": gpu_device,
            "budget_bytes": gpu_budget_bytes,
            "mps": mps_mode,
            "mps_pipe_directory": mps_pipe,
        },
        "diagnostics_interval": float(diagnostics_interval),
        "fail_fast": fail_fast,
        "allow_shared_cwd": allow_shared_cwd,
        "launch_gate": launch_gate,
        "members": members,
    }
    return canonical


def parse_uvm_final_marker(path: Path) -> dict[str, Any] | None:
    marker: dict[str, Any] | None = None
    try:
        handle = path.open("r", encoding="utf-8", errors="replace")
    except OSError:
        return None
    with handle:
        for line in handle:
            if (
                "MSIC_KOKKOS_CUDAUVM_BACKING_V1" not in line
                or "event=finalize" not in line
            ):
                continue
            parsed: dict[str, Any] = {}
            for token in line.split():
                if "=" not in token:
                    continue
                key, value = token.split("=", 1)
                if key in UVM_INTEGER_FIELDS:
                    try:
                        parsed[key] = int(value)
                    except ValueError:
                        parsed[key] = value
                else:
                    parsed[key] = value
            parsed["raw"] = line.strip()
            marker = parsed
    return marker


def _json_append(path: Path, payload: dict[str, Any]) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, ensure_ascii=False, sort_keys=True) + "\n")
    path.chmod(0o600)


def _number(value: str) -> int | float | str | None:
    cleaned = value.strip()
    if not cleaned or cleaned in ("N/A", "[N/A]", "Not Supported"):
        return None
    try:
        return int(cleaned)
    except ValueError:
        try:
            return float(cleaned)
        except ValueError:
            return cleaned


def _gpu_sample(device: str) -> dict[str, Any]:
    sample: dict[str, Any] = {
        "timestamp": utc_now(),
        "requested_device": device,
        "units": {
            "memory.total": "MiB",
            "memory.used": "MiB",
            "memory.free": "MiB",
            "utilization.gpu": "percent",
            "utilization.memory": "percent",
            "temperature.gpu": "C",
            "power.draw": "W",
            "processes.used_memory_mib": "MiB",
        },
    }
    executable = shutil.which("nvidia-smi")
    if executable is None:
        sample["error"] = "nvidia-smi is not installed"
        return sample
    command = [
        executable,
        f"--query-gpu={','.join(GPU_QUERY_FIELDS)}",
        "--format=csv,noheader,nounits",
        "--id",
        device,
    ]
    try:
        result = subprocess.run(
            command,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=3,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        sample["error"] = str(error)
        return sample
    if result.returncode != 0:
        sample["error"] = result.stderr.strip() or f"nvidia-smi exited {result.returncode}"
        return sample
    rows = list(csv.reader(line for line in result.stdout.splitlines() if line.strip()))
    devices: list[dict[str, Any]] = []
    for row in rows:
        if len(row) != len(GPU_QUERY_FIELDS):
            continue
        record: dict[str, Any] = {}
        for field, raw_value in zip(GPU_QUERY_FIELDS, row):
            value = raw_value.strip()
            record[field] = (
                _number(value)
                if field in GPU_INTEGER_FIELDS or field == "power.draw"
                else value
            )
        devices.append(record)
    sample["devices"] = devices
    process_command = [
        executable,
        "--query-compute-apps=pid,gpu_uuid,used_memory,process_name",
        "--format=csv,noheader,nounits",
        "--id",
        device,
    ]
    try:
        process_result = subprocess.run(
            process_command,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=3,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        sample["process_error"] = str(error)
        return sample
    processes: list[dict[str, Any]] = []
    if process_result.returncode == 0:
        for row in csv.reader(
            line for line in process_result.stdout.splitlines() if line.strip()
        ):
            if len(row) != 4:
                continue
            processes.append(
                {
                    "pid": _number(row[0]),
                    "gpu_uuid": row[1].strip(),
                    "used_memory_mib": _number(row[2]),
                    "process_name": row[3].strip(),
                }
            )
    else:
        sample["process_error"] = (
            process_result.stderr.strip()
            or f"nvidia-smi compute query exited {process_result.returncode}"
        )
    sample["processes"] = processes
    return sample


def _process_environment_value(pid: int, key: str) -> str | None:
    prefix = f"{key}=".encode("utf-8")
    try:
        entries = Path(f"/proc/{pid}/environ").read_bytes().split(b"\0")
    except OSError:
        return None
    for entry in entries:
        if entry.startswith(prefix):
            return entry[len(prefix) :].decode("utf-8", errors="replace")
    return None


def _enable_child_subreaper() -> bool:
    try:
        libc = ctypes.CDLL(None, use_errno=True)
        return libc.prctl(PR_SET_CHILD_SUBREAPER, 1, 0, 0, 0) == 0
    except (AttributeError, OSError):
        return False


def _process_records() -> dict[int, tuple[str, int, int]]:
    """Return PID -> (state, parent PID, process group) from procfs."""
    records: dict[int, tuple[str, int, int]] = {}
    try:
        entries = list(Path("/proc").iterdir())
    except OSError:
        return records
    for entry in entries:
        if not entry.name.isdecimal():
            continue
        try:
            stat = (entry / "stat").read_text(encoding="utf-8")
            fields = stat[stat.rfind(")") + 2 :].split()
            records[int(entry.name)] = (fields[0], int(fields[1]), int(fields[2]))
        except (OSError, ValueError, IndexError):
            continue
    return records


def _descendants_from_records(
    root_pid: int, records: dict[int, tuple[str, int, int]]
) -> set[int]:
    children: dict[int, list[int]] = {}
    for pid, (_state, parent, _process_group) in records.items():
        children.setdefault(parent, []).append(pid)
    descendants: set[int] = set()
    pending = list(children.get(root_pid, ()))
    while pending:
        pid = pending.pop()
        if pid in descendants:
            continue
        descendants.add(pid)
        pending.extend(children.get(pid, ()))
    return {
        pid
        for pid in descendants
        if pid in records and records[pid][0] != "Z"
    }


def _member_process_metrics_by_token(
    tokens: set[str],
) -> dict[str, dict[str, Any]]:
    metrics: dict[str, dict[str, Any]] = {
        token: {
            "pids": [],
            "process_groups": set(),
            "cpu_affinity": {},
            "rss_bytes": 0,
            "threads": 0,
            "cpu_ticks": 0,
        }
        for token in tokens
    }
    try:
        entries = list(Path("/proc").iterdir())
    except OSError:
        entries = []
    for entry in entries:
        if not entry.name.isdecimal():
            continue
        pid = int(entry.name)
        token = _process_environment_value(pid, "SIMQ_MEMBER_TOKEN")
        if token not in metrics:
            continue
        record = metrics[token]
        try:
            status = (entry / "status").read_text(encoding="utf-8")
            stat = (entry / "stat").read_text(encoding="utf-8")
        except OSError:
            continue
        record["pids"].append(pid)
        for line in status.splitlines():
            if line.startswith("VmRSS:"):
                try:
                    record["rss_bytes"] += int(line.split()[1]) * 1024
                except (ValueError, IndexError):
                    pass
            elif line.startswith("Threads:"):
                try:
                    record["threads"] += int(line.split()[1])
                except (ValueError, IndexError):
                    pass
            elif line.startswith("Cpus_allowed_list:"):
                record["cpu_affinity"][str(pid)] = line.split(":", 1)[1].strip()
        fields = stat[stat.rfind(")") + 2 :].split()
        try:
            if fields[0] != "Z":
                record["process_groups"].add(int(fields[2]))
            record["cpu_ticks"] += int(fields[11]) + int(fields[12])
        except (ValueError, IndexError):
            pass
    for record in metrics.values():
        record["pids"] = sorted(record["pids"])
        record["process_groups"] = sorted(record["process_groups"])
    return metrics


def _process_snapshot(pid: int) -> dict[str, Any]:
    snapshot: dict[str, Any] = {"pid": pid}
    try:
        stat = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
        fields = stat[stat.rfind(")") + 2 :].split()
        snapshot.update(
            {
                "state": fields[0],
                "parent_pid": int(fields[1]),
                "process_group": int(fields[2]),
            }
        )
    except (OSError, ValueError, IndexError):
        pass
    try:
        status = Path(f"/proc/{pid}/status").read_text(encoding="utf-8")
        for line in status.splitlines():
            if line.startswith("VmRSS:"):
                snapshot["rss_bytes"] = int(line.split()[1]) * 1024
            elif line.startswith("Threads:"):
                snapshot["threads"] = int(line.split()[1])
            elif line.startswith("Cpus_allowed_list:"):
                snapshot["cpu_affinity"] = line.split(":", 1)[1].strip()
    except (OSError, ValueError, IndexError):
        pass
    try:
        arguments = Path(f"/proc/{pid}/cmdline").read_bytes().split(b"\0")
        snapshot["command"] = " ".join(
            argument.decode("utf-8", errors="replace") for argument in arguments if argument
        )[:1000]
    except OSError:
        pass
    return snapshot


@dataclass
class ActiveMember:
    row: dict[str, Any]
    process: subprocess.Popen[bytes]
    log_handle: BinaryIO
    lane: int
    cpu: int
    started_monotonic: float
    known_process_groups: set[int] = field(default_factory=set)
    descendant_wait_logged: bool = False
    completion_quiet_scans: int = 0
    affinity_warning: tuple[str, ...] | None = None
    launch_gate_ready: bool = False
    launch_gate_log_handle: BinaryIO | None = None
    launch_gate_scan_offset: int = 0
    launch_gate_scan_tail: bytes = b""
    launch_gate_timeout_logged: bool = False
    launch_gate_scan_failed: bool = False


class EnsembleSupervisor:
    def __init__(self, job_id: int, *, store: QueueStore | None = None) -> None:
        self.store = store or QueueStore()
        self.job_id = job_id
        self.job = self.store.require(job_id)
        self.ensemble = self.store.require_ensemble(job_id)
        canonical_spec = dict(json.loads(str(self.ensemble["spec_json"])))
        try:
            self.launch_gate = _canonical_launch_gate(
                canonical_spec.get("launch_gate"), field="stored launch_gate"
            )
        except ValueError as error:
            raise RuntimeError(
                f"invalid stored ensemble launch gate: {error}"
            ) from error
        self.members = self.store.list_ensemble_members(job_id)
        self.artifact_dir = Path(str(self.ensemble["artifact_dir"]))
        self.events_path = self.artifact_dir / "events.jsonl"
        self.gpu_samples_path = self.artifact_dir / "gpu_samples.jsonl"
        self.process_samples_path = self.artifact_dir / "process_samples.jsonl"
        self.stop_signal: int | None = None
        self.outer_worker_supports_drain = (
            os.environ.get(ENSEMBLE_DRAIN_PROTOCOL_ENV)
            == f"{ENSEMBLE_CONTROL_PROTOCOL}:{self.job['launch_token']}"
        )
        self.active: dict[int, ActiveMember] = {}
        self.parent_environment = dict(json.loads(str(self.job["env_json"])))
        # These protected values are assigned by the outer worker at launch,
        # after the submission environment was snapshotted in the database.
        # Every member must retain them so outer-job cancellation and stale
        # process recovery continue to cover the entire cohort.
        for key in (
            "SIMQ_JOB_ID",
            "SIMQ_LAUNCH_TOKEN",
            "SIMQ_CPU_LIST",
            "SIMQ_CPU_COUNT",
            "SIMQ_LOG",
        ):
            if key in os.environ:
                self.parent_environment[key] = os.environ[key]
        cpu_list = os.environ.get("SIMQ_CPU_LIST", "")
        try:
            self.cpus = parse_cpu_list(cpu_list)
        except ValueError as error:
            raise RuntimeError(f"invalid outer SIMQ_CPU_LIST {cpu_list!r}: {error}") from error

    def _handle_stop(self, signum: int, _frame: FrameType | None) -> None:
        self.stop_signal = signum

    def _install_signal_handlers(self) -> None:
        signal.signal(signal.SIGTERM, self._handle_stop)
        signal.signal(signal.SIGINT, self._handle_stop)

    def _event(self, event: str, **fields: Any) -> None:
        payload = {"timestamp": utc_now(), "event": event, "job_id": self.job_id}
        payload.update(fields)
        _json_append(self.events_path, payload)
        detail = " ".join(f"{key}={value}" for key, value in fields.items())
        print(
            f"[{payload['timestamp']}] ensemble {event}"
            f"{' ' if detail else ''}{detail}",
            flush=True,
        )

    def _verify_outer_contract(self) -> None:
        if os.environ.get("SIMQ_JOB_ID") != str(self.job_id):
            raise RuntimeError("ensemble supervisor must run as its queued outer job")
        if os.environ.get("SIMQ_LAUNCH_TOKEN") != str(self.job["launch_token"]):
            raise RuntimeError("ensemble supervisor launch token does not match its job")
        expected = int(self.ensemble["max_parallel"])
        if len(self.cpus) < expected:
            raise RuntimeError(
                f"ensemble needs {expected} lane CPUs but outer allocation exposes {len(self.cpus)}"
            )
        if shutil.which("taskset") is None:
            raise RuntimeError("taskset is not installed")
        if not _enable_child_subreaper():
            raise RuntimeError(
                "Linux child-subreaper support is required for safe ensemble lanes"
            )

    def _verify_mps(self) -> None:
        mode = str(self.ensemble["mps_mode"])
        if mode != "off":
            raise RuntimeError(
                "stored ensemble requests unsupported MPS mode; direct signal "
                "cleanup is unsafe for active MPS clients"
            )
        disabled_path = self.artifact_dir / "mps-disabled"
        if disabled_path.exists():
            raise RuntimeError(
                f"reserved MPS bypass path unexpectedly exists: {disabled_path}"
            )
        self._event("mps-disabled", pipe_directory=str(disabled_path))

    @staticmethod
    def _resolve_runtime_cwd(member: dict[str, Any]) -> Path:
        submitted = Path(str(member["cwd"]))
        try:
            cwd = submitted.resolve(strict=True)
        except OSError as error:
            raise RuntimeError(f"working directory is unavailable: {submitted}: {error}") from error
        if not cwd.is_dir() or not cwd.is_relative_to(project_root()):
            raise RuntimeError(
                f"working directory resolves outside {project_root()}: {cwd}"
            )
        return cwd

    def _member_environment(
        self, member: dict[str, Any], *, lane: int, cpu: int
    ) -> dict[str, str]:
        environment = dict(self.parent_environment)
        environment.update(dict(json.loads(str(member["env_json"]))))
        for key in AFFINITY_ENVIRONMENT_TO_CLEAR:
            environment.pop(key, None)
        environment.update(THREAD_LIMIT_ENVIRONMENT)
        environment.update(MPI_AFFINITY_ENVIRONMENT)
        environment.update(
            {
                "CUDA_VISIBLE_DEVICES": str(self.ensemble["gpu_device"]),
                "SIMQ_ENSEMBLE_JOB_ID": str(self.job_id),
                "SIMQ_ENSEMBLE_MEMBER_ID": str(member["member_index"]),
                "SIMQ_ENSEMBLE_MEMBER_NAME": str(member["name"]),
                "SIMQ_MEMBER_TOKEN": str(member["member_token"]),
                "SIMQ_LANE": str(lane),
                "SIMQ_LANE_CPU": str(cpu),
                "SIMQ_CPU_LIST": str(cpu),
                "SIMQ_CPU_COUNT": "1",
                "SIMQ_LOG": str(member["log_path"]),
                "PWD": str(member["cwd"]),
            }
        )
        environment["CUDA_MPS_PIPE_DIRECTORY"] = str(
            self.artifact_dir / "mps-disabled"
        )
        return environment

    @staticmethod
    def _diagnostic_environment(environment: dict[str, str]) -> dict[str, str]:
        return {
            key: value
            for key, value in sorted(environment.items())
            if key in DIAGNOSTIC_ENVIRONMENT_NAMES
            or key.startswith("MSIC_KOKKOS_CUDAUVM")
            or key.startswith("MSIC_CUDA")
            or key.startswith("MSIC_PANZER")
        }

    def _launch(self, member: dict[str, Any], *, lane: int, cpu: int) -> bool:
        index = int(member["member_index"])
        log_path = Path(str(member["log_path"]))
        log_handle: BinaryIO | None = None
        launch_gate_log_handle: BinaryIO | None = None
        process: subprocess.Popen[bytes] | None = None
        try:
            cwd = self._resolve_runtime_cwd(member)
            argv = list(json.loads(str(member["argv_json"])))
            if not argv or not all(isinstance(argument, str) for argument in argv):
                raise RuntimeError("stored argv is not a nonempty string array")
            environment = self._member_environment(member, lane=lane, cpu=cpu)
            taskset = shutil.which("taskset")
            if taskset is None:
                raise RuntimeError("taskset is not installed")
            command = [taskset, "--cpu-list", str(cpu), *argv]
            log_handle = log_path.open("ab", buffering=0)
            log_path.chmod(0o600)
            header = [
                f"===== simq ensemble {self.job_id} member {index}: {member['name']} =====\n",
                f"started:   {utc_now()}\n",
                f"cwd:       {cwd}\n",
                f"lane:      {lane}\n",
                f"P-core:    {cpu}\n",
                f"GPU claim: {member['gpu_admission_bytes']} bytes\n",
                f"UVM arena: {member['uvm_arena_bytes']} bytes\n",
                f"command:   {shlex.join(argv)}\n",
            ]
            for key, value in self._diagnostic_environment(environment).items():
                header.append(f"environment: {key}={value}\n")
            header.append("========================================\n")
            log_handle.write("".join(header).encode("utf-8", errors="replace"))
            launch_gate_scan_offset = log_handle.tell()
            if self.launch_gate is not None:
                launch_gate_log_handle = log_path.open("rb", buffering=0)
                writer_stat = os.fstat(log_handle.fileno())
                reader_stat = os.fstat(launch_gate_log_handle.fileno())
                if (
                    not stat.S_ISREG(writer_stat.st_mode)
                    or not stat.S_ISREG(reader_stat.st_mode)
                    or (writer_stat.st_dev, writer_stat.st_ino)
                    != (reader_stat.st_dev, reader_stat.st_ino)
                ):
                    raise RuntimeError(
                        "launch gate writer and pinned reader do not reference "
                        "the same regular member log"
                    )
            process = subprocess.Popen(
                command,
                cwd=cwd,
                env=environment,
                stdin=subprocess.DEVNULL,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
            # The database launch reservation was committed before Popen.
            # Track the process locally before filling its PID so run()'s
            # finally block still owns the just-launched tree if bookkeeping
            # itself fails.
            self.active[lane] = ActiveMember(
                row=member,
                process=process,
                log_handle=log_handle,
                lane=lane,
                cpu=cpu,
                started_monotonic=time.monotonic(),
                known_process_groups={process.pid},
                launch_gate_ready=self.launch_gate is None,
                launch_gate_log_handle=launch_gate_log_handle,
                launch_gate_scan_offset=launch_gate_scan_offset,
            )
            self.store.record_member_process(
                self.job_id,
                index,
                pid=process.pid,
                pgid=process.pid,
            )
            self._event(
                "member-launch",
                member_index=index,
                member_name=member["name"],
                lane=lane,
                cpu=cpu,
                pid=process.pid,
                gpu_admission_bytes=member["gpu_admission_bytes"],
                uvm_arena_bytes=member["uvm_arena_bytes"],
            )
            if self.launch_gate is not None:
                self._event(
                    "launch-gate-armed",
                    member_index=index,
                    member_name=member["name"],
                    lane=lane,
                    timeout_seconds=self.launch_gate["timeout_seconds"],
                )
            return True
        except (OSError, RuntimeError, ValueError, json.JSONDecodeError) as error:
            if process is not None:
                # A post-Popen failure cannot be reduced to an ordinary
                # member launch failure: descendants may already exist. Keep
                # the process registered and unwind through the supervisor's
                # cohort-wide cleanup.
                raise RuntimeError(
                    f"post-launch tracking failed for member {member['name']!r}: {error}"
                ) from error
            if log_handle is not None:
                log_handle.close()
            if launch_gate_log_handle is not None:
                launch_gate_log_handle.close()
            try:
                log_path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
                with log_path.open("ab") as handle:
                    handle.write(f"simq: member launch failed: {error}\n".encode())
                log_path.chmod(0o600)
            except OSError:
                pass
            self.store.finish_member(
                self.job_id,
                index,
                state="failed",
                exit_code=None,
                message=f"launch failed: {error}",
            )
            self._event(
                "member-launch-failed",
                member_index=index,
                member_name=member["name"],
                error=str(error),
            )
            return False

    def _refresh_launch_gate(self) -> None:
        """Release a follower only after its predecessor emits the gate marker.

        The scan starts after SimQueue's own member header so a marker embedded
        in argv or environment diagnostics cannot release the gate. Any read
        failure or truncation is fail-closed and aborts the ensemble before a
        follower can launch.
        """
        if self.launch_gate is None:
            return
        marker = str(self.launch_gate["contains"]).encode("utf-8")
        timeout_seconds = float(self.launch_gate["timeout_seconds"])
        for active in self.active.values():
            if active.launch_gate_ready:
                continue
            handle = active.launch_gate_log_handle
            if handle is None:
                raise RuntimeError(
                    f"launch gate has no pinned log reader for member "
                    f"{active.row['name']!r}"
                )
            try:
                snapshot_size = os.fstat(handle.fileno()).st_size
                if snapshot_size < active.launch_gate_scan_offset:
                    raise RuntimeError("member log was truncated after launch")
                while active.launch_gate_scan_offset < snapshot_size:
                    chunk = os.pread(
                        handle.fileno(),
                        min(65536, snapshot_size - active.launch_gate_scan_offset),
                        active.launch_gate_scan_offset,
                    )
                    if not chunk:
                        break
                    combined = active.launch_gate_scan_tail + chunk
                    active.launch_gate_scan_offset += len(chunk)
                    overlap = len(marker) - 1
                    active.launch_gate_scan_tail = (
                        combined[-overlap:] if overlap > 0 else b""
                    )
                    if marker in combined:
                        active.launch_gate_ready = True
                        self._event(
                            "launch-gate-ready",
                            member_index=active.row["member_index"],
                            member_name=active.row["name"],
                            lane=active.lane,
                            elapsed_seconds=round(
                                time.monotonic() - active.started_monotonic, 3
                            ),
                            after_timeout=int(active.launch_gate_timeout_logged),
                        )
                        break
            except (OSError, RuntimeError) as error:
                active.launch_gate_scan_failed = True
                self._event(
                    "launch-gate-read-failed",
                    member_index=active.row["member_index"],
                    member_name=active.row["name"],
                    lane=active.lane,
                    error=str(error),
                    policy="fail-ensemble-no-follower-launch",
                )
                raise RuntimeError(
                    f"launch gate log scan failed for member "
                    f"{active.row['name']!r}: {error}"
                ) from error

            if active.launch_gate_ready:
                continue

            elapsed = time.monotonic() - active.started_monotonic
            if elapsed >= timeout_seconds and not active.launch_gate_timeout_logged:
                active.launch_gate_timeout_logged = True
                self._event(
                    "launch-gate-timeout",
                    member_index=active.row["member_index"],
                    member_name=active.row["name"],
                    lane=active.lane,
                    elapsed_seconds=round(elapsed, 3),
                    timeout_seconds=timeout_seconds,
                    policy="fail-ensemble-no-follower-launch",
                )
                raise RuntimeError(
                    f"launch gate timed out after {timeout_seconds:g}s for "
                    f"member {active.row['name']!r}"
                )

    def _cohort_processes(self) -> tuple[dict[int, set[int]], set[int]]:
        """Classify live supervisor descendants by active lane.

        Because this supervisor is a subreaper, a member which daemonizes,
        starts a new session, and clears its environment is reparented here
        instead of escaping to the outer worker. Such an orphan cannot be
        assigned to a lane reliably, so it conservatively keeps every exited
        member's lane and GPU claim occupied until it is gone.
        """
        records = _process_records()
        leader_pids = {active.process.pid for active in self.active.values()}
        reaped: list[int] = []
        for pid, (state, parent, _process_group) in tuple(records.items()):
            if state != "Z" or parent != os.getpid() or pid in leader_pids:
                continue
            try:
                waited, _status = os.waitpid(pid, os.WNOHANG)
            except (ChildProcessError, OSError):
                continue
            if waited == pid:
                reaped.append(pid)
                records.pop(pid, None)
        if reaped:
            self._event("adopted-descendant-reap", pids=sorted(reaped))
        descendants = _descendants_from_records(os.getpid(), records)
        owned: dict[int, set[int]] = {lane: set() for lane in self.active}
        recognized: set[int] = set()
        token_lanes = {
            str(active.row["member_token"]): lane
            for lane, active in self.active.items()
        }
        for lane, active in self.active.items():
            leader = active.process.pid
            leader_record = records.get(leader)
            if leader_record is not None and leader_record[0] != "Z":
                tree = {leader, *_descendants_from_records(leader, records)}
                owned[lane].update(tree)
                recognized.update(tree)
                active.known_process_groups.update(
                    records[pid][2] for pid in tree if pid in records
                )

        # A token attributes new sessions even after their launcher exits.
        # Only cohort descendants are inspected, avoiding a full procfs scan
        # once per member on every scheduler tick.
        for pid in descendants:
            lane = token_lanes.get(
                _process_environment_value(pid, "SIMQ_MEMBER_TOKEN") or ""
            )
            if lane is None:
                continue
            owned[lane].add(pid)
            recognized.add(pid)
            active = self.active[lane]
            if pid in records:
                active.known_process_groups.add(records[pid][2])

        for lane, active in self.active.items():
            live_groups = {
                records[pid][2]
                for pid in descendants
                if pid in records
                and records[pid][0] != "Z"
                and records[pid][2] in active.known_process_groups
            }
            active.known_process_groups = live_groups
            for pid, (state, _parent, process_group) in records.items():
                if (
                    pid in descendants
                    and state != "Z"
                    and process_group in active.known_process_groups
                ):
                    owned[lane].add(pid)
                    recognized.add(pid)
        return owned, descendants - recognized

    def _unattributed_descendants(self) -> list[int]:
        _owned, unattributed = self._cohort_processes()
        return sorted(unattributed)

    def _finish_active(
        self,
        lane: int,
        *,
        forced_state: str | None = None,
        forced_message: str | None = None,
    ) -> str:
        active = self.active.pop(lane)
        return_code = active.process.wait()
        active.log_handle.write("========================================\n".encode())
        active.log_handle.write(f"finished:  {utc_now()}\n".encode())
        active.log_handle.write(f"exit code: {return_code}\n".encode())
        active.log_handle.close()
        if active.launch_gate_log_handle is not None:
            active.launch_gate_log_handle.close()
        marker = parse_uvm_final_marker(Path(str(active.row["log_path"])))
        if forced_state is not None:
            state = forced_state
        else:
            state = "succeeded" if return_code == 0 else "failed"
        message: str | None = None
        if state == "failed":
            message = f"command exited with status {return_code}"
        elif state == "canceled":
            message = "canceled with ensemble"
        elif state == "interrupted":
            message = "ensemble supervisor interrupted"
        if forced_message is not None:
            message = forced_message
        if (
            marker is None
            and active.row.get("uvm_arena_bytes")
            and state == "succeeded"
        ):
            message = "configured UVM arena emitted no final allocator marker"
        self.store.finish_member(
            self.job_id,
            int(active.row["member_index"]),
            state=state,
            exit_code=return_code,
            message=message,
            uvm_marker=marker,
        )
        self._event(
            "member-finish",
            member_index=active.row["member_index"],
            member_name=active.row["name"],
            lane=lane,
            cpu=active.cpu,
            state=state,
            exit_code=return_code,
            elapsed_seconds=round(time.monotonic() - active.started_monotonic, 3),
            uvm_peak_reserved_bytes=(
                marker.get("peak_reserved_bytes") if marker else None
            ),
            uvm_teardown_ok=(marker.get("teardown_ok") if marker else None),
        )
        return state

    def _sample(self) -> None:
        gpu = _gpu_sample(str(self.ensemble["gpu_device"]))
        gpu.update(
            {
                "job_id": self.job_id,
                "gpu_budget_bytes": int(self.ensemble["gpu_budget_bytes"]),
                "active_gpu_admission_bytes": sum(
                    int(active.row["gpu_admission_bytes"])
                    for active in self.active.values()
                ),
            }
        )
        gpu["mps"] = {"mode": "off"}
        _json_append(self.gpu_samples_path, gpu)
        process_sample: dict[str, Any] = {
            "timestamp": utc_now(),
            "job_id": self.job_id,
            "members": [],
            "unattributed_descendants": [
                _process_snapshot(pid) for pid in self._unattributed_descendants()
            ],
        }
        metrics_by_token = _member_process_metrics_by_token(
            {str(active.row["member_token"]) for active in self.active.values()}
        )
        for lane, active in sorted(self.active.items()):
            metrics = metrics_by_token[str(active.row["member_token"])]
            record = {
                "member_index": active.row["member_index"],
                "member_name": active.row["name"],
                "lane": lane,
                "cpu": active.cpu,
                "leader_pid": active.process.pid,
                "elapsed_seconds": round(
                    time.monotonic() - active.started_monotonic, 3
                ),
                **metrics,
            }
            violations = tuple(
                f"{pid}:{affinity}"
                for pid, affinity in sorted(metrics["cpu_affinity"].items())
                if affinity != str(active.cpu)
            )
            if violations and violations != active.affinity_warning:
                self._event(
                    "member-affinity-warning",
                    member_index=active.row["member_index"],
                    member_name=active.row["name"],
                    expected_cpu=active.cpu,
                    observed=violations,
                )
            active.affinity_warning = violations or None
            process_sample["members"].append(record)
            self.store.update_member_sample(
                self.job_id,
                int(active.row["member_index"]),
                rss_bytes=int(metrics["rss_bytes"]),
                threads=int(metrics["threads"]),
            )
        _json_append(self.process_samples_path, process_sample)

    def _signal_active(self, sig: signal.Signals) -> None:
        records = _process_records()
        process_groups = {
            records[pid][2]
            for pid in _descendants_from_records(os.getpid(), records)
            if pid in records
        }
        for active in list(self.active.values()):
            process_groups.update(active.known_process_groups)
        process_groups.discard(os.getpgrp())
        for process_group in sorted(process_groups, reverse=True):
            try:
                os.killpg(process_group, sig)
            except ProcessLookupError:
                pass

    def _kill_active_until_quiet(self) -> None:
        """Kill the cohort and close the final fork/setsid race.

        This is used only while unwinding an exceptional supervisor exit.
        Two empty scans are required after all leaders have exited because a
        just-forked child may not yet have been reparented to this subreaper.
        """
        quiet_scans = 0
        deadline = time.monotonic() + CLEANUP_KILL_SETTLE_SECONDS
        survivors: set[int] = set()
        while quiet_scans < 2:
            self._signal_active(signal.SIGKILL)
            owned, unattributed = self._cohort_processes()
            survivors = set(unattributed)
            for pids in owned.values():
                survivors.update(pids)
            survivors.update(
                active.process.pid
                for active in self.active.values()
                if active.process.poll() is None
            )
            if survivors:
                quiet_scans = 0
            else:
                quiet_scans += 1
            if quiet_scans >= 2:
                return
            if time.monotonic() >= deadline:
                rendered = ",".join(str(pid) for pid in sorted(survivors))
                raise RuntimeError(
                    "ensemble process(es) survived repeated SIGKILL scans"
                    f"{': ' + rendered if rendered else ''}"
                )
            time.sleep(POLL_SECONDS)

    def run(self) -> int:
        self._verify_outer_contract()
        self.store.register_ensemble_supervisor(
            self.job_id,
            pid=os.getpid(),
            launch_token=str(os.environ["SIMQ_LAUNCH_TOKEN"]),
            protocol=ENSEMBLE_CONTROL_PROTOCOL,
        )
        self._install_signal_handlers()
        self._verify_mps()
        pending = deque(member for member in self.members if member["state"] == "queued")
        invalid_members = [
            member
            for member in self.members
            if member["state"] not in ("queued", *MEMBER_TERMINAL_STATES)
        ]
        if invalid_members:
            rendered = ", ".join(
                f"{member['member_index']}:{member['state']}"
                for member in invalid_members
            )
            raise RuntimeError(
                "ensemble resume found nonterminal members already admitted: "
                f"{rendered}"
            )
        max_parallel = int(self.ensemble["max_parallel"])
        gpu_budget = int(self.ensemble["gpu_budget_bytes"])
        fail_fast = bool(self.ensemble["fail_fast"])
        prior_counts = Counter(str(member["state"]) for member in self.members)
        self._event(
            "start",
            members=len(self.members),
            pending_members=len(pending),
            previously_succeeded=prior_counts.get("succeeded", 0),
            max_parallel=max_parallel,
            cpus=",".join(str(cpu) for cpu in self.cpus[:max_parallel]),
            gpu_device=self.ensemble["gpu_device"],
            gpu_budget_bytes=gpu_budget,
            mps_mode=self.ensemble["mps_mode"],
            fail_fast=fail_fast,
            launch_gate=(
                self.launch_gate["kind"] if self.launch_gate is not None else "off"
            ),
        )
        next_sample = 0.0
        stop_kind: str | None = None
        termination_started: float | None = None
        kill_escalated = False
        saw_failure = bool(prior_counts.get("failed", 0))
        admission_blocked: tuple[int, int] | None = None
        launch_gate_blocked: tuple[int, int] | None = None
        observed_pause_generation: int | None = None
        held_pause_generation: int | None = None
        if fail_fast and saw_failure and pending:
            stop_kind = "fail-fast"
            self._event(
                "fail-fast-triggered",
                reason="failure preserved from an earlier drained execution",
            )
        try:
            while pending or self.active:
                control = self.store.ensemble_control_snapshot(self.job_id)
                launch_paused = bool(control["launch_paused"])
                control_generation = int(control["control_generation"])
                if launch_paused and observed_pause_generation != control_generation:
                    observed_pause_generation = control_generation
                    self._event(
                        "pause-observed",
                        control_generation=control_generation,
                        active_members=len(self.active),
                        queued_members=len(pending),
                    )
                elif not launch_paused and observed_pause_generation is not None:
                    self._event(
                        "resume-observed",
                        control_generation=control_generation,
                        queued_members=len(pending),
                    )
                    observed_pause_generation = None
                    held_pause_generation = None

                if stop_kind is None and control["cancel_requested_at"] is not None:
                    stop_kind = "canceled"
                    self._event("cancellation-observed")
                elif stop_kind is None and self.stop_signal is not None:
                    stop_kind = "interrupted"
                    self._event("stop-signal", signal=self.stop_signal)

                # Cancellation and supervisor signals win over a gate timeout
                # observed in the same polling cycle. Once the cohort is
                # stopping there is no follower left to admit.
                if stop_kind is None:
                    self._refresh_launch_gate()

                owned_processes, unattributed_descendants = self._cohort_processes()
                completion_stop_kind = stop_kind
                trigger_fail_fast = False
                for lane in list(self.active):
                    active = self.active[lane]
                    if active.process.poll() is None:
                        active.completion_quiet_scans = 0
                        continue
                    descendants = sorted(
                        (owned_processes.get(lane, set()) - {active.process.pid})
                        | set(unattributed_descendants)
                    )
                    if descendants:
                        active.completion_quiet_scans = 0
                        if not active.descendant_wait_logged:
                            self._event(
                                "member-leader-exited-waiting-descendants",
                                member_index=active.row["member_index"],
                                member_name=active.row["name"],
                                lane=lane,
                                leader_exit_code=active.process.returncode,
                                descendant_pids=descendants,
                            )
                            active.descendant_wait_logged = True
                        continue
                    active.completion_quiet_scans += 1
                    # The leader can fork and exit between procfs enumeration
                    # and poll(). A second empty scan after its observed exit
                    # lets subreaper reparenting settle before this lane and
                    # GPU claim are released.
                    if active.completion_quiet_scans < 2:
                        continue
                    gate_marker_missing = (
                        self.launch_gate is not None
                        and not active.launch_gate_ready
                        and completion_stop_kind is None
                    )
                    if gate_marker_missing:
                        self._event(
                            "launch-gate-marker-missing",
                            member_index=active.row["member_index"],
                            member_name=active.row["name"],
                            lane=lane,
                            exit_code=active.process.returncode,
                            policy="fail-ensemble-no-follower-launch",
                        )
                    forced = (
                        "canceled"
                        if completion_stop_kind in ("canceled", "fail-fast")
                        else "interrupted"
                        if completion_stop_kind == "interrupted"
                        else "failed"
                        if gate_marker_missing
                        else None
                    )
                    state = self._finish_active(
                        lane,
                        forced_state=forced,
                        forced_message=(
                            "member exited before its required launch-gate marker"
                            if gate_marker_missing
                            else None
                        ),
                    )
                    if state == "failed":
                        saw_failure = True
                        if (
                            (fail_fast or gate_marker_missing)
                            and completion_stop_kind is None
                        ):
                            trigger_fail_fast = True

                if trigger_fail_fast and stop_kind is None:
                    stop_kind = "fail-fast"
                    self._event("fail-fast-triggered")

                # A pause closes only future launch admission. Already admitted
                # members keep running and are reaped normally. Once the cohort
                # is quiet, persist a generation-specific drain marker before
                # asking the outer worker to suspend this same parent job.
                if (
                    stop_kind is None
                    and launch_paused
                    and pending
                    and not self.active
                ):
                    if not self.outer_worker_supports_drain:
                        # A supervisor launched by a pre-upgrade worker must
                        # stay alive while paused. Returning the reserved code
                        # would make that old worker fail the parent and
                        # terminalize every pending member. This compatibility
                        # mode still stops backfill; a one-time worker restart
                        # enables allocation-releasing suspension.
                        if held_pause_generation != control_generation:
                            held_pause_generation = control_generation
                            self._event(
                                "pause-held-for-legacy-worker",
                                control_generation=control_generation,
                                queued_members=len(pending),
                            )
                    else:
                        if self.store.mark_ensemble_drained(
                            self.job_id, control_generation
                        ):
                            # Nothing after the durable marker may turn this
                            # into an ordinary supervisor failure: the outer
                            # worker uses the reserved exit code to preserve
                            # the pending rows. Diagnostics are best-effort once
                            # suspension is committed.
                            try:
                                self._sample()
                                self._event(
                                    "pause-drained",
                                    control_generation=control_generation,
                                    queued_members=len(pending),
                                    exit_code=ENSEMBLE_DRAIN_EXIT_CODE,
                                )
                            except Exception as error:
                                try:
                                    print(
                                        "simq: ensemble drain committed but final "
                                        f"diagnostics failed: {error}",
                                        flush=True,
                                    )
                                except Exception:
                                    pass
                            return ENSEMBLE_DRAIN_EXIT_CODE
                        # Resume or cancellation won the transaction. Re-read
                        # control before either launching or stopping.
                        continue

                if stop_kind is not None:
                    pending_state = "interrupted" if stop_kind == "interrupted" else "canceled"
                    if pending:
                        self.store.finish_unstarted_members(
                            self.job_id,
                            state=pending_state,
                            message=(
                                "not launched because ensemble was interrupted"
                                if stop_kind == "interrupted"
                                else "not launched after ensemble stop"
                            ),
                        )
                        pending.clear()
                    if self.active and termination_started is None:
                        termination_started = time.monotonic()
                        self._signal_active(signal.SIGTERM)
                        self._event("member-termination", signal="SIGTERM")
                    elif (
                        self.active
                        and termination_started is not None
                        and time.monotonic() - termination_started >= FAIL_FAST_GRACE_SECONDS
                    ):
                        self._signal_active(signal.SIGKILL)
                        if not kill_escalated:
                            self._event("member-termination", signal="SIGKILL")
                            kill_escalated = True
                elif not launch_paused:
                    free_lanes = [
                        lane for lane in range(max_parallel) if lane not in self.active
                    ]
                    while free_lanes and pending:
                        gate_holders = [
                            active
                            for active in self.active.values()
                            if not active.launch_gate_ready
                        ]
                        if gate_holders:
                            holder = min(
                                gate_holders,
                                key=lambda active: int(active.row["member_index"]),
                            )
                            key = (int(holder.row["member_index"]), len(pending))
                            if launch_gate_blocked != key:
                                self._event(
                                    "launch-gate-wait",
                                    holder_member_index=holder.row["member_index"],
                                    holder_member_name=holder.row["name"],
                                    holder_lane=holder.lane,
                                    queued_members=len(pending),
                                )
                                launch_gate_blocked = key
                            break
                        launch_gate_blocked = None
                        active_gpu = sum(
                            int(active.row["gpu_admission_bytes"])
                            for active in self.active.values()
                        )
                        eligible_offset: int | None = None
                        for offset, member in enumerate(pending):
                            if active_gpu + int(member["gpu_admission_bytes"]) <= gpu_budget:
                                eligible_offset = offset
                                break
                        if eligible_offset is None:
                            key = (active_gpu, len(pending))
                            if admission_blocked != key:
                                self._event(
                                    "admission-wait",
                                    active_gpu_admission_bytes=active_gpu,
                                    gpu_budget_bytes=gpu_budget,
                                    queued_members=len(pending),
                                )
                                admission_blocked = key
                            break
                        member = pending[eligible_offset]
                        # Narrow the cancellation/admission race: a gate scan
                        # or admission search may take long enough for the
                        # parent to be canceled after the top-of-loop check.
                        # Recheck immediately before dequeue and Popen.
                        if self.store.cancellation_requested(self.job_id):
                            stop_kind = "canceled"
                            self._event("cancellation-observed-before-launch")
                            break
                        if self.stop_signal is not None:
                            stop_kind = "interrupted"
                            self._event(
                                "stop-signal-before-launch", signal=self.stop_signal
                            )
                            break
                        lane = free_lanes[0]
                        if not self.store.reserve_member_launch(
                            self.job_id,
                            int(member["member_index"]),
                            lane=lane,
                            cpu=self.cpus[lane],
                        ):
                            # Pause and cancellation are linearized by the
                            # reservation transaction. Leave the member in the
                            # local pending deque and observe the winner on the
                            # next scheduler cycle.
                            break
                        del pending[eligible_offset]
                        free_lanes.pop(0)
                        launched = self._launch(
                            member, lane=lane, cpu=self.cpus[lane]
                        )
                        admission_blocked = None
                        if not launched:
                            saw_failure = True
                            if fail_fast:
                                stop_kind = "fail-fast"
                                self._event("fail-fast-triggered")
                                break

                now = time.monotonic()
                if now >= next_sample:
                    self._sample()
                    next_sample = time.monotonic() + float(
                        self.ensemble["diagnostics_interval"]
                    )
                if pending or self.active:
                    time.sleep(POLL_SECONDS)

            self._sample()
            members = self.store.list_ensemble_members(self.job_id)
            counts = Counter(str(member["state"]) for member in members)
            self._event("finish", **dict(sorted(counts.items())))
            if stop_kind in ("canceled", "interrupted"):
                return 130
            return 1 if saw_failure or counts.get("failed", 0) else 0
        finally:
            if self.active:
                self._kill_active_until_quiet()
                for lane in list(self.active):
                    try:
                        self._finish_active(lane, forced_state="interrupted")
                    except Exception:
                        pass


def run_ensemble(job_id: int) -> int:
    try:
        return EnsembleSupervisor(job_id).run()
    except Exception as error:
        print(f"simq: ensemble supervisor failed: {error}", flush=True)
        try:
            store = QueueStore()
            store.finish_unstarted_members(
                job_id,
                state="interrupted",
                message=f"ensemble supervisor failed: {error}",
            )
        except Exception:
            pass
        return 2


def _launch_gate_runtime_status(
    *, artifact_dir: Path, job_state: str
) -> dict[str, Any]:
    holder_index: int | None = None
    holder_name: str | None = None
    armed_at: str | None = None
    passed: set[int] = set()
    failure_event: str | None = None
    events_path = artifact_dir / "events.jsonl"
    try:
        lines = events_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        lines = []
    for line in lines:
        try:
            event = json.loads(line)
        except (json.JSONDecodeError, TypeError):
            # A status read can race the single append which is currently
            # extending the final JSONL line. The next status call will see it.
            continue
        if not isinstance(event, dict):
            continue
        kind = event.get("event")
        if kind in ("member-launch", "launch-gate-armed"):
            try:
                holder_index = int(event["member_index"])
            except (KeyError, TypeError, ValueError):
                continue
            holder_name = str(event.get("member_name", holder_index))
            armed_at = str(event.get("timestamp")) if event.get("timestamp") else None
        elif kind == "launch-gate-ready":
            try:
                ready_index = int(event["member_index"])
            except (KeyError, TypeError, ValueError):
                continue
            passed.add(ready_index)
            if holder_index == ready_index:
                holder_index = None
                holder_name = None
                armed_at = None
        elif kind in (
            "launch-gate-timeout",
            "launch-gate-read-failed",
            "launch-gate-marker-missing",
        ):
            failure_event = str(kind)

    if failure_event is not None:
        state = "failed"
    elif job_state == "succeeded":
        state = "complete"
    elif job_state in ("failed", "canceled", "interrupted"):
        state = job_state
    elif holder_index is not None:
        state = "waiting-marker"
    elif passed:
        state = "released"
    elif job_state in ("starting", "running", "cancelling"):
        state = "awaiting-first-launch"
    else:
        state = "queued"
    return {
        "state": state,
        "holder_member_index": holder_index,
        "holder_member_name": holder_name,
        "holder_armed_at": armed_at,
        "holder_elapsed_seconds": (
            _elapsed_seconds(armed_at, None) if armed_at is not None else None
        ),
        "passed_member_indices": sorted(passed),
        "passed_count": len(passed),
        "failure_event": failure_event,
    }


def ensemble_status_payload(job_id: int, store: QueueStore) -> dict[str, Any]:
    job = store.require(job_id)
    ensemble = store.require_ensemble(job_id)
    members = store.list_ensemble_members(job_id)
    counts = Counter(str(member["state"]) for member in members)
    launch_paused = bool(ensemble["launch_paused"])
    job_state = str(job["state"])
    if job_state in ("succeeded", "failed", "canceled", "interrupted"):
        control_state = "terminal"
    elif job_state == "cancelling":
        control_state = "cancelling"
    elif launch_paused and job_state in ("starting", "running"):
        control_state = "draining"
    elif launch_paused and job_state == "queued":
        control_state = "paused"
    elif (
        not launch_paused
        and job_state == "queued"
        and ensemble["drained_generation"] is not None
    ):
        control_state = "resuming"
    elif job_state in ("starting", "running"):
        control_state = "running"
    else:
        control_state = "queued"
    public_job = {
        key: job.get(key)
        for key in (
            "id",
            "name",
            "state",
            "requested_cores",
            "cpu_list",
            "created_at",
            "started_at",
            "finished_at",
            "exit_code",
            "message",
            "retry_of",
            "attempt",
            "log_path",
        )
    }
    canonical_spec = dict(json.loads(str(ensemble["spec_json"])))
    public_ensemble = {
        "schema_version": ensemble["schema_version"],
        "max_parallel": int(ensemble["max_parallel"]),
        "gpu_device": ensemble["gpu_device"],
        "gpu_budget_bytes": int(ensemble["gpu_budget_bytes"]),
        "diagnostics_interval": float(ensemble["diagnostics_interval"]),
        "mps_mode": ensemble["mps_mode"],
        "mps_pipe_directory": ensemble["mps_pipe_directory"],
        "fail_fast": bool(ensemble["fail_fast"]),
        "allow_shared_cwd": bool(ensemble["allow_shared_cwd"]),
        "launch_gate": canonical_spec.get("launch_gate"),
        "artifact_dir": ensemble["artifact_dir"],
        "launch_paused": launch_paused,
        "control_state": control_state,
        "control_generation": int(ensemble["control_generation"]),
        "control_changed_at": ensemble["control_changed_at"],
        "drained_generation": (
            int(ensemble["drained_generation"])
            if ensemble["drained_generation"] is not None
            else None
        ),
        "drained_at": ensemble["drained_at"],
        "allocation_released": job_state
        not in ("starting", "running", "cancelling"),
    }
    if public_ensemble["launch_gate"] is not None:
        public_ensemble["launch_gate_runtime"] = _launch_gate_runtime_status(
            artifact_dir=Path(str(ensemble["artifact_dir"])),
            job_state=str(job["state"]),
        )
    else:
        public_ensemble["launch_gate_runtime"] = None
    integer_fields = {
        "member_index",
        "gpu_admission_bytes",
        "uvm_arena_bytes",
        "lane",
        "cpu",
        "exit_code",
        "pid",
        "pgid",
        "last_rss_bytes",
        "peak_rss_bytes",
        "last_threads",
        "uvm_limit_bytes",
        "uvm_peak_live_bytes",
        "uvm_peak_reserved_bytes",
        "uvm_pool_reclaims",
        "uvm_teardown_ok",
    }
    public_members: list[dict[str, Any]] = []
    for member in members:
        public = {
            key: value
            for key, value in member.items()
            if key not in ("argv_json", "env_json", "member_token", "uvm_marker_json")
        }
        public["argv"] = list(json.loads(str(member["argv_json"])))
        public["elapsed_seconds"] = _elapsed_seconds(
            member.get("started_at"), member.get("finished_at")
        )
        for field in integer_fields:
            if public.get(field) is not None:
                public[field] = int(public[field])
        public_members.append(public)
    active_gpu = sum(
        int(member["gpu_admission_bytes"])
        for member in members
        if member["state"] == "running"
    )
    summary: dict[str, Any] = {
        "total": len(members),
        "active_gpu_admission_bytes": active_gpu,
    }
    for state in ("queued", "running", *MEMBER_TERMINAL_STATES):
        summary[state] = counts.get(state, 0)
    return {
        "job": public_job,
        "ensemble": public_ensemble,
        "summary": summary,
        "members": public_members,
    }


def print_ensemble_status(payload: dict[str, Any]) -> None:
    job = payload["job"]
    ensemble = payload["ensemble"]
    summary = payload["summary"]
    print(
        f"Ensemble {job['id']}: {job['name']} [{job['state']}]  "
        f"members={summary['succeeded']}/{summary['total']} succeeded, "
        f"active={summary['running']}, queued={summary['queued']}, "
        f"failed={summary['failed']}, canceled={summary['canceled']}, "
        f"interrupted={summary['interrupted']}"
    )
    print(
        f"Lanes: {ensemble['max_parallel']}  CPUs: {job['cpu_list'] or 'not assigned'}  "
        f"GPU: {ensemble['gpu_device']}  admission: "
        f"{format_byte_size(summary['active_gpu_admission_bytes'])}/"
        f"{format_byte_size(ensemble['gpu_budget_bytes'])}  MPS: {ensemble['mps_mode']}"
    )
    print(
        f"Launch control: {ensemble['control_state']}; "
        f"generation={ensemble['control_generation']}; "
        f"allocation={'released' if ensemble['allocation_released'] else 'held'}"
    )
    if ensemble.get("launch_gate") is not None:
        gate = ensemble["launch_gate"]
        runtime = ensemble["launch_gate_runtime"]
        print(
            "Launch gate: one pre-marker member per GPU; "
            f"kind={gate['kind']} timeout={gate['timeout_seconds']:g}s "
            f"contains={gate['contains']!r}"
        )
        holder = (
            f"{runtime['holder_member_index']}:{runtime['holder_member_name']}"
            if runtime["holder_member_index"] is not None
            else "none"
        )
        elapsed = (
            f"{runtime['holder_elapsed_seconds']:.1f}s"
            if runtime["holder_elapsed_seconds"] is not None
            else "-"
        )
        print(
            f"Gate runtime: {runtime['state']}; holder={holder}; "
            f"wait={elapsed}; passed={runtime['passed_count']}"
            + (
                f"; failure={runtime['failure_event']}"
                if runtime["failure_event"] is not None
                else ""
            )
        )
    print(f"Artifacts: {ensemble['artifact_dir']}")
    print(
        f"{'IDX':>4}  {'STATE':<11} {'LANE':>4} {'CPU':>4} {'PID':>7} "
        f"{'ELAPSED':>8} {'GPU CLAIM':>10} {'UVM PEAK':>10}  NAME"
    )
    for member in payload["members"]:
        elapsed = (
            f"{member['elapsed_seconds']:.1f}s"
            if member["elapsed_seconds"] is not None
            else "-"
        )
        uvm_peak = (
            format_byte_size(member["uvm_peak_reserved_bytes"])
            if member["uvm_peak_reserved_bytes"] is not None
            else "-"
        )
        print(
            f"{member['member_index']:>4}  {member['state']:<11} "
            f"{str(member['lane'] if member['lane'] is not None else '-'):>4} "
            f"{str(member['cpu'] if member['cpu'] is not None else '-'):>4} "
            f"{str(member['pid'] if member['pid'] is not None else '-'):>7} "
            f"{elapsed:>8} {format_byte_size(member['gpu_admission_bytes']):>10} "
            f"{uvm_peak:>10}  {member['name']}"
        )
    detailed = [member for member in payload["members"] if member.get("message")]
    if detailed:
        print("Details:")
        for member in detailed:
            print(
                f"  {member['member_index']} {member['name']}: "
                f"{member['message']}  log={member['log_path']}"
            )


def tail_member_log(
    job_id: int,
    selector: str,
    store: QueueStore,
    *,
    lines: int,
    follow: bool,
) -> int:
    member = store.require_ensemble_member(job_id, selector)
    path = Path(str(member["log_path"]))
    command = ["tail", "-n", str(lines)]
    if follow:
        command.append("-f")
    command.append(str(path))
    try:
        return subprocess.call(command)
    except KeyboardInterrupt:
        return 130
