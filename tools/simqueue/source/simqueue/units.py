"""systemd user-unit backend: each job is its own transient service.

Jobs never run as children of the scheduler.  ``systemd-run --user`` starts
``jobwrap.py`` as the main process of ``<prefix>-<id>-<token12>.service``
with the CPU affinity, memory limits, time limit, working directory and
placement environment set on the unit.  The scheduler can therefore stop,
crash or be replaced without affecting running jobs; on restart it re-adopts
them by unit name and launch token.

Borrow jobs (0.6.1) run in a separate top-level slice whose CPU weight is
``idle`` (cgroup ``cpu.idle``), with ``CPUWeight=1``, ``SCHED_IDLE`` and idle
I/O class on the unit itself, so every other process of the user manager
(owner job units, terminal sessions, tools holding leases) wins the CPU.
They are preempted with the cgroup freezer (``systemctl --user freeze``) or
by stopping the unit; a frozen unit must be thawed before systemd accepts a
stop request.
"""

from __future__ import annotations

import json
import re
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from .nodeconfig import BorrowConfig, LaunchConfig


JOBWRAP = Path(__file__).resolve().with_name("jobwrap.py")
SHOW_PROPERTIES = (
    "Id",
    "LoadState",
    "ActiveState",
    "SubState",
    "Result",
    "MainPID",
    "ExecMainStatus",
    "ExecMainCode",
    "InvocationID",
    "Environment",
    "FreezerState",
)
ACTIVE_STATES = {"active", "activating", "deactivating", "reloading", "refreshing"}


class UnitError(RuntimeError):
    pass


def unit_name(prefix: str, job_id: int, launch_token: str) -> str:
    token = re.sub(r"[^0-9a-fA-F]", "", launch_token)[:12]
    return f"{prefix}-{job_id}-{token}.service"


@dataclass(frozen=True)
class UnitState:
    name: str
    load_state: str
    active_state: str
    sub_state: str
    result: str
    main_pid: int
    exec_main_status: int | None
    exec_main_code: int | None
    invocation_id: str
    environment: str
    freezer_state: str = ""

    @property
    def active(self) -> bool:
        return self.active_state in ACTIVE_STATES

    @property
    def found(self) -> bool:
        return self.load_state not in ("not-found", "")

    @property
    def frozen(self) -> bool:
        return self.freezer_state in ("frozen", "freezing")

    def has_token(self, token: str) -> bool:
        return f"SIMQ_LAUNCH_TOKEN={token}" in self.environment


def _systemd_value(value: str) -> str:
    return value.replace("%", "%%")


def build_run_command(
    launch: LaunchConfig,
    *,
    unit: str,
    description: str,
    cwd: str,
    cpus: Iterable[int],
    mem_bytes: int,
    swap_bytes: int,
    runtime_max_s: int | None,
    nice: int,
    log_path: str,
    spec_path: str,
    setenv: dict[str, str],
    wrapper: Path = JOBWRAP,
    borrow: BorrowConfig | None = None,
) -> list[str]:
    python = launch.wrapper_python
    stop_post = _systemd_value(
        " ".join(shlex.quote(part) for part in (python, str(wrapper), "stop-post", spec_path))
    )
    command = [
        launch.systemd_run,
        "--user",
        f"--unit={unit}",
        f"--slice={launch.slice if borrow is None else borrow.slice}",
        "--service-type=exec",
        "--collect",
        "--quiet",
        f"--description={_systemd_value(description)}",
        f"--working-directory={cwd}",
        "--property=CPUAffinity=" + " ".join(str(cpu) for cpu in sorted(cpus)),
        f"--property=MemoryMax={mem_bytes}",
        f"--property=MemorySwapMax={swap_bytes}",
        "--property=OOMPolicy=kill",
        "--property=KillMode=control-group",
        f"--property=TimeoutStopSec={max(1, int(round(launch.kill_grace_seconds)))}s",
        "--property=SendSIGKILL=yes",
        "--property=UMask=0077",
        f"--property=StandardOutput=append:{log_path}",
        f"--property=StandardError=append:{log_path}",
        f"--property=ExecStopPost={stop_post}",
    ]
    if runtime_max_s is not None:
        command.append(f"--property=RuntimeMaxSec={int(runtime_max_s)}s")
    if borrow is not None:
        command.extend(
            [
                f"--property=CPUWeight={int(borrow.unit_cpu_weight)}",
                "--property=CPUSchedulingPolicy=idle",
                "--property=IOSchedulingClass=idle",
            ]
        )
    elif nice:
        command.append(f"--property=Nice={int(nice)}")
    for key, value in sorted(setenv.items()):
        command.append(f"--setenv={key}={value}")
    command.extend(["--", python, str(wrapper), "run", spec_path])
    return command


def start_unit(command: list[str], *, timeout: float = 60.0) -> None:
    try:
        result = subprocess.run(
            command,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise UnitError(f"systemd-run failed: {error}") from error
    if result.returncode != 0:
        detail = (result.stderr or result.stdout).strip().splitlines()
        raise UnitError("systemd-run failed: " + (detail[-1] if detail else f"exit {result.returncode}"))


def _parse_show(text: str, names: list[str]) -> dict[str, UnitState]:
    blocks: list[dict[str, str]] = []
    current: dict[str, str] = {}
    for line in text.splitlines():
        if not line.strip():
            if current:
                blocks.append(current)
                current = {}
            continue
        key, _, value = line.partition("=")
        current[key] = value
    if current:
        blocks.append(current)
    states: dict[str, UnitState] = {}
    for name, block in zip(names, blocks):
        def number(key: str) -> int | None:
            try:
                return int(block.get(key, ""))
            except ValueError:
                return None

        states[name] = UnitState(
            name=name,
            load_state=block.get("LoadState", ""),
            active_state=block.get("ActiveState", ""),
            sub_state=block.get("SubState", ""),
            result=block.get("Result", ""),
            main_pid=number("MainPID") or 0,
            exec_main_status=number("ExecMainStatus"),
            exec_main_code=number("ExecMainCode"),
            invocation_id=block.get("InvocationID", ""),
            environment=block.get("Environment", ""),
            freezer_state=block.get("FreezerState", ""),
        )
    return states


def show_units(launch: LaunchConfig, names: list[str]) -> dict[str, UnitState]:
    if not names:
        return {}
    command = [launch.systemctl, "--user", "show"]
    for prop in SHOW_PROPERTIES:
        command.append(f"--property={prop}")
    command.append("--")
    command.extend(names)
    try:
        result = subprocess.run(
            command,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise UnitError(f"systemctl show failed: {error}") from error
    if result.returncode != 0:
        raise UnitError("systemctl show failed: " + result.stderr.strip())
    states = _parse_show(result.stdout, names)
    if len(states) != len(names):
        raise UnitError("systemctl show returned an unexpected number of units")
    return states


def _systemctl(launch: LaunchConfig, *arguments: str, timeout: float = 30) -> subprocess.CompletedProcess[str] | None:
    try:
        return subprocess.run(
            [launch.systemctl, "--user", *arguments],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None


def stop_unit(launch: LaunchConfig, name: str) -> bool:
    # systemd refuses to stop a frozen unit ("Cannot perform operation on
    # frozen unit"); thawing a unit that is not frozen is a no-op.
    _systemctl(launch, "thaw", "--", name)
    result = _systemctl(launch, "stop", "--no-block", "--", name)
    return result is not None and result.returncode == 0


def freeze_unit(launch: LaunchConfig, name: str) -> bool:
    """Freeze every process of the unit (cgroup v2 freezer); returns once frozen."""
    result = _systemctl(launch, "freeze", "--", name)
    return result is not None and result.returncode == 0


def thaw_unit(launch: LaunchConfig, name: str) -> bool:
    result = _systemctl(launch, "thaw", "--", name)
    return result is not None and result.returncode == 0


def ensure_borrow_slice(launch: LaunchConfig, borrow: BorrowConfig) -> str | None:
    """Give the borrow slice its CPU weight (runtime drop-in); returns an
    error text, or ``None`` when the slice is configured."""
    current = _systemctl(launch, "show", "--property=CPUWeight", "--value", "--", borrow.slice)
    wanted = str(borrow.slice_cpu_weight)
    if current is not None and current.returncode == 0 and current.stdout.strip() == wanted:
        return None
    result = _systemctl(launch, "set-property", "--runtime", "--", borrow.slice, f"CPUWeight={wanted}")
    if result is None or result.returncode != 0:
        detail = "" if result is None else (result.stderr or result.stdout).strip()
        return f"cannot set CPUWeight={wanted} on {borrow.slice}: {detail or 'systemctl failed'}"
    check = _systemctl(launch, "show", "--property=CPUWeight", "--value", "--", borrow.slice)
    if check is None or check.stdout.strip() != wanted:
        shown = "?" if check is None else check.stdout.strip()
        return f"{borrow.slice} reports CPUWeight={shown} after setting {wanted}"
    return None


def read_json(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return value if isinstance(value, dict) else None
