from __future__ import annotations

import argparse
import json
import os
import pwd
import shlex
import signal
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Sequence

from . import __version__
from .config import (
    database_path,
    detect_pcores,
    format_cpu_list,
    parse_cpu_list,
    project_root,
    state_dir,
)
from .exclusive import helper_status, load_config
from .ensemble import (
    ensemble_status_payload,
    load_manifest,
    print_ensemble_status,
    run_ensemble,
    tail_member_log,
)
from .nodeconfig import (
    EXCLUSIVE_TIMING_FLAG,
    SMT_FLAG,
    ConfigError,
    NodeConfig,
    config_path,
    describe_alias,
    format_size,
    format_time_limit,
    load_node_config,
    parse_size,
    parse_time_limit,
    smt_job,
)
from .memory import (
    Series,
    auto_request,
    describe_basis,
    format_gib,
    over_request_warning,
    read_node_memory,
    read_usage,
    sparkline,
    unit_cgroup_path,
    water_marks,
)
from .github import job_phase, targets_from_args
from .store import ACTIVE_STATES, TERMINAL_STATES, QueueStore
from .units import stop_unit
from .worker import run_worker


EXCLUSIVE_START_TIMEOUT_SECONDS = 10.0


class WorkerStartError(RuntimeError):
    """The queue entry is durable, but its configured worker could not start."""


def _cli_pcores(exclusive: dict[str, Any] | None = None) -> tuple[int, ...]:
    """Return queue P-cores even when this outside CLI is shielded onto E-cores."""
    config = load_config() if exclusive is None else exclusive
    if config is not None:
        return parse_cpu_list(str(config["pcores"]))
    return detect_pcores()


def _parse_time(value: str | None) -> datetime | None:
    if not value:
        return None
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None


def _compact_duration(seconds: float) -> str:
    seconds = max(0, int(seconds))
    if seconds < 60:
        return f"{seconds}s"
    minutes, seconds = divmod(seconds, 60)
    if minutes < 60:
        return f"{minutes}m{seconds:02d}s"
    hours, minutes = divmod(minutes, 60)
    if hours < 24:
        return f"{hours}h{minutes:02d}m"
    days, hours = divmod(hours, 24)
    return f"{days}d{hours:02d}h"


def _job_age(job: dict[str, Any]) -> str:
    now = datetime.now(UTC)
    start = _parse_time(job.get("started_at"))
    end = _parse_time(job.get("finished_at"))
    created = _parse_time(job.get("created_at"))
    if start is not None:
        return _compact_duration(((end or now) - start).total_seconds())
    if created is not None:
        return _compact_duration((now - created).total_seconds())
    return "-"


def _abbreviate(value: str, width: int) -> str:
    if len(value) <= width:
        return value
    if width <= 1:
        return value[:width]
    return value[: width - 1] + "…"


def _job_cgroup(job: dict[str, Any], config: NodeConfig | None) -> Path | None:
    if config is None or not job.get("unit_name") or job.get("launch_backend") != "systemd-user":
        return None
    slice_name = config.borrow.slice if job.get("borrow") else config.launch.slice
    return unit_cgroup_path(str(job["unit_name"]), slice_name)


def _live_usage(job: dict[str, Any], config: NodeConfig | None) -> Any:
    """The running job unit's cgroup reading now (0.7.0), else ``None``."""
    if job.get("state") not in ACTIVE_STATES:
        return None
    cgroup = _job_cgroup(job, config)
    return None if cgroup is None else read_usage(cgroup)


def _mem_column(job: dict[str, Any], config: NodeConfig | None) -> str:
    """MEM of ``simq list``: now/request while running, peak/request after."""
    request = job.get("mem_bytes")
    request_text = format_size(int(request)) if request is not None else "-"
    if job["state"] in ACTIVE_STATES:
        usage = _live_usage(job, config)
        current = usage.current if usage is not None else job.get("mem_current_bytes")
        return f"{format_gib(current)}/{request_text}" if current is not None else request_text
    if job["state"] in TERMINAL_STATES and job.get("mem_peak_bytes") is not None:
        return f"{format_gib(int(job['mem_peak_bytes']))}/{request_text}"
    return request_text


def _use_column(job: dict[str, Any], config: NodeConfig | None) -> str:
    """USE of ``simq list`` (0.8.0): CPU time / (runtime x allocated cores)."""
    if job["state"] in TERMINAL_STATES:
        value = job.get("cpu_utilization")
        return "-" if value is None else f"{100 * float(value):.0f}%"
    if job["state"] not in ACTIVE_STATES or not job.get("cpu_list"):
        return "-"
    usage = _live_usage(job, config)
    start = _parse_time(job.get("started_at"))
    if usage is None or usage.cpu_usec is None or start is None:
        return "-"
    seconds = (datetime.now(UTC) - start).total_seconds()
    cores = len(parse_cpu_list(str(job["cpu_list"])))
    return "-" if seconds <= 0 or not cores else f"{100 * usage.cpu_usec / 1e6 / seconds / cores:.0f}%"


def _print_node_job_rows(jobs: Sequence[dict[str, Any]]) -> None:
    root = project_root()
    borrow_shown = False
    try:
        config = load_node_config()
    except ConfigError:
        config = None
    print(
        f"{'ID':>5}  {'STATE':<11} {'PARTITION':<11} {'CPU':<7} {'GPU':<4} {'MEM':<11} {'USE':>4} "
        f"{'TIME':>8}  {'NAME':<24} DIRECTORY / REASON"
    )
    for job in jobs:
        cwd = Path(str(job["cwd"]))
        try:
            directory = str(cwd.relative_to(root)) or "."
        except ValueError:
            directory = str(cwd)
        cpu = str(job.get("cpu_list") or f"{job['requested_cores']}c")
        gpus = job.get("gpu_list")
        gpu = (str(len(gpus.split(","))) if gpus else "-") if job["state"] != "queued" else (
            str(job.get("gpus_requested") or "-")
        )
        detail = directory
        message = str(job.get("message") or "")
        if job["state"] == "queued" and message.startswith("pending: "):
            detail = message[len("pending: "):]
        partition = str(job.get("partition") or "-")
        state = str(job["state"])
        if job.get("borrow"):
            borrow_shown = True
            partition = "~" + partition
            if job.get("frozen_at") and state in ACTIVE_STATES:
                state = "frozen"
        print(
            f"{job['id']:>5}  {state:<11} {_abbreviate(partition, 11):<11} "
            f"{_abbreviate(cpu, 7):<7} {gpu:<4} {_abbreviate(_mem_column(job, config), 11):<11} "
            f"{_use_column(job, config):>4} {_job_age(job):>8}  "
            f"{_abbreviate(str(job['name']), 24):<24} {_abbreviate(detail, 60)}"
        )
    if borrow_shown:
        print("~ borrow job: preemptible, idle CPU priority (simq show ID for its placement and history)")


def _print_job_rows(
    jobs: Sequence[dict[str, Any]], *, empty_message: str
) -> None:
    if not jobs:
        print(empty_message)
        return
    try:
        node_mode = load_node_config() is not None
    except ConfigError:
        node_mode = False
    if node_mode:
        _print_node_job_rows(jobs)
        return
    root = project_root()
    print(f"{'ID':>5}  {'STATE':<11} {'CPU':<5} {'TIME':>8}  {'NAME':<24} DIRECTORY")
    for job in jobs:
        cwd = Path(str(job["cwd"]))
        try:
            directory = str(cwd.relative_to(root)) or "."
        except ValueError:
            directory = str(cwd)
        cpu = str(job.get("cpu_list") or f"{job['requested_cores']}c")
        print(
            f"{job['id']:>5}  {job['state']:<11} {_abbreviate(cpu, 5):<5} "
            f"{_job_age(job):>8}  {_abbreviate(str(job['name']), 24):<24} "
            f"{_abbreviate(directory, 45)}"
        )


def _print_jobs(store: QueueStore, *, include_all: bool, limit: int) -> None:
    if store.is_paused():
        print("Queue is paused (a running job is allowed to finish).")
    _print_job_rows(
        store.list_jobs(include_all=include_all, limit=limit),
        empty_message="No jobs have been submitted.",
    )


def _print_active_jobs(store: QueueStore) -> None:
    if store.is_paused():
        print("Queue is paused (a running job is allowed to finish).")
    _print_job_rows(store.active_jobs(), empty_message="No active jobs.")


def _print_queued_jobs(store: QueueStore) -> None:
    if store.is_paused():
        print("Queue is paused.")
    _print_job_rows(store.queued_jobs(), empty_message="No queued jobs.")


def _parse_environment(items: Sequence[str], export: str | None = None) -> tuple[dict[str, str], list[str]]:
    """The job's environment (0.7.0: an allowlist, see envpolicy.py); the
    secrets given with -e are kept aside (``_parse_environment.secrets``)."""
    environment, _secrets, override_keys, _withheld = _capture_environment(items, export)
    return environment, override_keys


def _capture_environment(items: Sequence[str], export: str | None = None):
    from .envpolicy import capture, parse_export

    overrides: list[tuple[str, str]] = []
    for item in items:
        if "=" not in item:
            raise ValueError(f"environment override must be KEY=VALUE, got {item!r}")
        key, value = item.split("=", 1)
        if not key or "\x00" in key or "=" in key:
            raise ValueError(f"invalid environment variable name {key!r}")
        overrides.append((key, value))
    return capture(os.environ, overrides, parse_export(export))


def _resolve_cwd(value: str | None) -> Path:
    cwd = Path(value).expanduser() if value else Path.cwd()
    cwd = cwd.resolve()
    root = project_root()
    if not cwd.is_dir():
        raise ValueError(f"working directory does not exist: {cwd}")
    if not cwd.is_relative_to(root):
        raise ValueError(f"working directory must be inside {root}: {cwd}")
    return cwd


def _start_system_unit(unit: str) -> None:
    try:
        result = subprocess.run(
            ["systemctl", "start", unit],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=15,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise WorkerStartError(f"cannot start {unit}: {error}") from error
    if result.returncode != 0:
        detail = result.stderr.strip().splitlines()
        message = detail[-1] if detail else f"systemctl exited {result.returncode}"
        raise WorkerStartError(f"cannot start {unit}: {message}")


def _start_exclusive_runner(store: QueueStore, config: dict[str, Any]) -> str:
    blocker = store.dispatch_blocker()
    if blocker is not None and not store.active_jobs():
        if blocker["kind"] == "queue-paused":
            return "paused"
        if blocker["kind"] == "ensemble-paused":
            return "ensemble-paused"
    unit = str(config["runner_unit"])
    deadline = time.monotonic() + EXCLUSIVE_START_TIMEOUT_SECONDS
    starts = 0
    while time.monotonic() < deadline:
        blocker = store.dispatch_blocker()
        if blocker is not None and not store.active_jobs():
            if blocker["kind"] == "queue-paused":
                return "paused"
            if blocker["kind"] == "ensemble-paused":
                return "ensemble-paused"
        active_state = _systemctl_property("ActiveState", unit=unit, user=False)
        if active_state not in ("active", "activating", "reloading"):
            if starts >= 3:
                break
            _start_system_unit(unit)
            starts += 1
            active_state = _systemctl_property("ActiveState", unit=unit, user=False)

        # A running/starting database row proves the drain worker has crossed
        # the submit-vs-exit race.  If a command was very short, an empty
        # unfinished set proves it already completed.
        if active_state in ("active", "activating", "reloading"):
            if store.active_jobs() or not store.unfinished_ids():
                return "exclusive-systemd"
        if store.is_paused():
            return "paused"
        if not store.unfinished_ids():
            return "exclusive-systemd"
        time.sleep(0.1)
    state = _systemctl_property("ActiveState", unit=unit, user=False)
    raise WorkerStartError(
        f"{unit} did not take ownership of the queued work (state: {state})"
    )


WAKE_FIFO = "worker.wake"
MUTATING_ACTIONS = {"cancel", "retry", "promote", "demote", "pause", "resume", "priority", "ensemble"}


def _wake_worker() -> None:
    """Wake the 0.7.0 event-driven worker at once (a byte into its FIFO).
    Never blocks; without a listening worker nothing happens."""
    try:
        descriptor = os.open(state_dir() / WAKE_FIFO, os.O_WRONLY | os.O_NONBLOCK)
    except OSError:
        return
    try:
        os.write(descriptor, b"c")
    except OSError:
        pass
    finally:
        os.close(descriptor)


def _nudge_worker(store: QueueStore) -> str:
    _wake_worker()
    config = load_config()
    blocker = store.dispatch_blocker()
    if blocker is not None and not store.active_jobs():
        if blocker["kind"] == "queue-paused":
            return "paused"
        if blocker["kind"] == "ensemble-paused":
            return "ensemble-paused"
    if os.environ.get("SIMQUEUE_NO_AUTOSTART") == "1":
        return "exclusive-disabled" if config is not None else "disabled"
    if config is not None:
        return _start_exclusive_runner(store, config)
    try:
        result = subprocess.run(
            ["systemctl", "--user", "start", "simqueue.service"],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if result.returncode == 0:
            return "systemd"
    except OSError:
        pass

    launcher = Path(__file__).resolve().parents[1] / "simq"
    worker_log = state_dir() / "worker.log"
    with worker_log.open("ab") as handle:
        subprocess.Popen(
            [sys.executable, str(launcher), "worker"],
            cwd=project_root(),
            stdin=subprocess.DEVNULL,
            stdout=handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            close_fds=True,
        )
    worker_log.chmod(0o600)
    return "detached"


NODE_OPTIONS = (
    "partition", "cpus", "gpus", "gpu_mem", "mem", "time", "priority",
    "borrow_from", "max_borrow_cores", "on_preempt", "l3_domains", "split_domains",
)
NODE_FLAGS = ("borrow", "exclusive_timing", "smt")


def _submission_name(args: argparse.Namespace) -> str | None:
    """The job name a command submission will get (as ``_enqueue`` sets it)."""
    if getattr(args, "name", None):
        return str(args.name)
    if getattr(args, "shell", None):
        try:
            words = shlex.split(args.shell)
        except ValueError:
            words = []
        return Path(words[0]).name if words else "shell-job"
    command = [part for part in (getattr(args, "command", None) or [])]
    if command and command[0] == "--":
        command = command[1:]
    if command:
        return Path(command[0]).name or command[0]
    return None


def _node_request(
    args: argparse.Namespace, config: NodeConfig, *, cores: int | None = None,
    store: QueueStore | None = None, name: str | None = None,
) -> tuple[int, dict[str, Any], dict[str, Any]]:
    """Resolve a submission's resource request against the node config.

    Returns (requested cores, stored scheduling fields, explicitly given
    fields).  Only the explicit fields enter the idempotency fingerprint, so
    a 0.5-style replay keeps its 0.5 fingerprint.

    0.6.3: ``-p`` may name an alias of a partition.  The job is recorded
    under the name given (so ``SIMQ_PARTITION`` and fingerprints keep it)
    and gets the alias's defaults; limits are the partition's.
    """
    partition = config.partition(args.partition)
    recorded = args.partition or config.default_partition
    where = partition.display_name()
    explicit: dict[str, Any] = {}
    cpu_request = None
    cpu_list_arg = getattr(args, "cpus", None)
    borrow = _borrow_options(args, config, recorded)
    if cpu_list_arg is not None:
        cpus = parse_cpu_list(cpu_list_arg)
        if not cpus:
            raise ValueError("--cpus needs at least one CPU")
        outside = sorted(set(cpus) - set(partition.cpus))
        if outside:
            raise ValueError(
                f"--cpus {cpu_list_arg}: CPU(s) {format_cpu_list(outside)} are not in partition "
                f"{where} ({format_cpu_list(partition.cpus)}); name each physical core by "
                "its partition CPU"
            )
        cpu_request = format_cpu_list(cpus)
        explicit["cpu_request"] = cpu_request
    # 0.8.0: L3 domains
    l3_count = getattr(args, "l3_domains", None)
    split_text = getattr(args, "split_domains", None)
    split: list[int] | None = None
    if split_text is not None:
        try:
            split = [int(part) for part in split_text.split(",") if part.strip()]
        except ValueError:
            raise ValueError("--split-domains takes core counts, e.g. 4,4") from None
        if len(split) < 2 or min(split) < 1:
            raise ValueError("--split-domains needs two or more positive core counts, e.g. 4,4")
    if (l3_count is not None or split is not None) and (cpu_list_arg is not None or borrow is not None):
        raise ValueError("--l3-domains and --split-domains cannot be combined with --cpus or --borrow")
    if l3_count is not None and split is not None:
        raise ValueError("use either --l3-domains or --split-domains")
    if l3_count is not None and l3_count < 1:
        raise ValueError("--l3-domains must be positive")
    requested = cores
    if split is not None:
        if requested is not None and requested != sum(split):
            raise ValueError(f"--cores {requested} does not match --split-domains (total {sum(split)})")
        requested = sum(split)
    if requested is None:
        requested = len(parse_cpu_list(cpu_request)) if cpu_request else partition.resolve_default_cpus()
        if borrow is not None and borrow["max_borrow_cores"] is not None:
            requested = 1  # elastic: 1 .. --max-borrow-cores
    if l3_count is not None and requested < l3_count:
        raise ValueError(f"--l3-domains {l3_count} needs at least {l3_count} cores")
    if cpu_request is not None and requested != len(parse_cpu_list(cpu_request)):
        raise ValueError("--cores and --cpus disagree on the number of cores")
    if borrow is not None:
        reachable = borrow["reachable"]
        if requested < 1 or requested > reachable:
            raise ValueError(
                f"--cores must be between 1 and {reachable} for a borrow job "
                f"(cores of {', '.join(borrow['eligible'])})"
            )
        maximum = borrow["max_borrow_cores"]
        if maximum is not None and not requested <= maximum <= reachable:
            raise ValueError(f"--max-borrow-cores must be between {requested} (--cores) and {reachable}")
    elif requested < 1 or requested > len(partition.cpus):
        raise ValueError(
            f"--cores must be between 1 and {len(partition.cpus)} in partition {where}"
        )
    gpus = getattr(args, "gpus", None)
    gpu_mem = None
    if getattr(args, "gpu_mem", None) is not None:
        gpu_mem = parse_size(args.gpu_mem, field_name="--gpu-mem")
        if gpu_mem <= 0:
            raise ValueError("--gpu-mem must be positive")
        explicit["gpu_mem_bytes"] = gpu_mem
        if gpus is None:
            gpus = 1
    if gpus is not None:
        if gpus < 0:
            raise ValueError("--gpus cannot be negative")
        explicit["gpus_requested"] = gpus
    gpus = gpus or 0
    if gpus > len(partition.gpus):
        raise ValueError(f"partition {where} has {len(partition.gpus)} GPU(s)")
    if gpus == 0 and gpu_mem is not None:
        raise ValueError("--gpu-mem needs --gpus 1 or more")
    if gpu_mem is not None:
        budgets = [config.gpus[uuid].budget_bytes for uuid in partition.gpus]
        if not any(gpu_mem <= budget for budget in budgets):
            raise ValueError(
                f"--gpu-mem {format_size(gpu_mem)} exceeds the GPU budget "
                f"({', '.join(format_size(b) for b in budgets)})"
            )
    mem = partition.default_mem
    args.mem_basis = None
    args.mem_warning = None
    mem_arg = getattr(args, "mem", None)
    history: list[dict[str, Any]] | None = None
    if mem_arg is not None and str(mem_arg).strip().lower() == "auto":
        # 0.7.0: from the recorded peaks of earlier jobs of this name (family)
        history = store.memory_history() if store is not None else []
        mem, args.mem_basis = auto_request(
            history, name=name or "", cores=requested, config=config.memory,
            min_mem=partition.min_mem,
            max_mem=partition.mem_limit(int(getattr(args, "priority", None) or 0), borrow is not None),
            default_mem=partition.default_mem,
        )
        explicit["mem_bytes"] = "auto"  # a replay keeps its fingerprint while the history grows
    elif mem_arg is None and (config.memory.auto_default or partition.mem_auto):
        # 0.7.1: no --mem: the history when there is one, else the default
        history = store.memory_history() if store is not None else []
        derived, basis = auto_request(
            history, name=name or "", cores=requested, config=config.memory,
            min_mem=partition.min_mem, max_mem=partition.max_mem, default_mem=partition.default_mem,
        )
        if basis.get("source") == "history":
            mem = derived
            args.mem_basis = dict(basis, implicit=True)
    elif mem_arg is not None:
        mem = parse_size(args.mem, field_name="--mem")
        limit = partition.mem_limit(int(getattr(args, "priority", None) or 0), borrow is not None)
        if mem <= 0 or mem > limit:
            raise ValueError(
                f"--mem must be between 1M and {format_size(limit)} in partition {where}"
                + ("" if limit != partition.max_mem or partition.priority_max_mem is None else
                   f" ({format_size(partition.priority_max_mem[1])} for owner jobs of priority "
                   f"{partition.priority_max_mem[0]} or more)")
            )
        explicit["mem_bytes"] = mem
        if store is not None and name:
            history = store.memory_history()
            args.mem_warning = over_request_warning(mem, history, name=name, cores=requested,
                                                    config=config.memory)
    time_limit = partition.default_time
    args.time_basis = None
    if getattr(args, "time", None) is not None and str(args.time).strip().lower() == "auto":
        # 0.8.0: from the runtimes of earlier jobs of this name (family)
        from .memory import auto_time_limit

        history = store.runtime_history() if store is not None else []
        time_limit, args.time_basis = auto_time_limit(history, name=name or "", default=partition.default_time,
                                                      maximum=partition.max_time)
        explicit["time_limit_s"] = "auto"
    elif getattr(args, "time", None) is not None:
        time_limit = parse_time_limit(args.time, field_name="--time")
        explicit["time_limit_s"] = time_limit
        if partition.max_time is not None and (time_limit is None or time_limit > partition.max_time):
            raise ValueError(
                f"--time exceeds partition {where}'s limit {format_time_limit(partition.max_time)}"
            )
    priority = getattr(args, "priority", None)
    if priority is not None:
        explicit["priority"] = priority
    if (borrow is None and requested >= 0.75 * len(partition.cpus) and len(partition.cpus) > 1
            and (priority or 0) < 10):
        # 0.8.0: a near-whole-pool job at default priority holds the queue head
        print(f"simq: warning: {requested} of {len(partition.cpus)} cores at priority {priority or 0}: while it "
              "waits for them, the reservation keeps every longer job behind it waiting; ask for fewer cores, give "
              "an accurate -t (or -t auto), or use --priority 10 for a pilot", file=sys.stderr)
    if args.partition is not None:
        explicit["partition"] = recorded
    scheduling = {
        "partition": recorded,
        "cpu_request": cpu_request,
        "gpus_requested": gpus,
        "gpu_mem_bytes": gpu_mem,
        "mem_bytes": mem,
        "swap_bytes": partition.swap_max,
        "time_limit_s": time_limit,
        "priority": priority or 0,
    }
    # 0.6.1 fields enter the idempotency fingerprint only when given.
    if borrow is not None:
        scheduling.update(
            {
                "borrow": 1,
                "borrow_from": borrow["borrow_from"],
                "max_borrow_cores": borrow["max_borrow_cores"],
                "on_preempt": borrow["on_preempt"],
            }
        )
        explicit["borrow"] = 1
        for key in ("borrow_from", "max_borrow_cores"):
            if borrow[key] is not None:
                explicit[key] = borrow[key]
        if borrow["on_preempt_given"]:
            explicit["on_preempt"] = borrow["on_preempt"]
    if getattr(args, "exclusive_timing", False):
        scheduling["exclusive_timing"] = 1
        explicit["exclusive_timing"] = 1
    if getattr(args, "smt", False):
        # 0.8.2: the job's CPU affinity adds its cores' SMT siblings.  Stored
        # as bit 2 of exclusive_timing together with bit 1: an SMT job's cores
        # and siblings are never lent (no schema change; 0.8.1 reads it as an
        # exclusive-timing job).
        candidates = parse_cpu_list(cpu_request) if cpu_request else partition.cpus
        lacking = [cpu for cpu, found in config.smt_siblings_of(candidates).items() if not found]
        if lacking:
            raise ValueError(f"--smt: CPU(s) {format_cpu_list(lacking)} of partition {where} have no free SMT "
                             "sibling on this node (SMT off, or the siblings are partition CPUs themselves)")
        scheduling["exclusive_timing"] = explicit["exclusive_timing"] = EXCLUSIVE_TIMING_FLAG | SMT_FLAG
    if l3_count is not None or split is not None:
        from .scheduler import Request, request_is_feasible

        probe = Request(job_id=0, partition=recorded, ncpus=requested, l3_domains=l3_count,
                        split_domains=None if split is None else tuple(split))
        reason = request_is_feasible(probe, partition, config)
        if reason is not None:
            raise ValueError(reason)
        if l3_count is not None:
            scheduling["l3_domains"] = explicit["l3_domains"] = l3_count
        if split is not None:
            scheduling["split_domains"] = explicit["split_domains"] = json.dumps(split)
    return requested, scheduling, explicit


def _borrow_options(args: argparse.Namespace, config: NodeConfig, home: str) -> dict[str, Any] | None:
    """Validate the 0.6.1 borrow options; ``None`` for an owner job.

    *home* and the ``--borrow-from`` names may be aliases (0.6.3); they are
    recorded as given and resolve to their partitions."""
    from .borrow import eligible_partitions, protected_cpus
    from .scheduler import Request

    borrow_from_text = getattr(args, "borrow_from", None)
    maximum = getattr(args, "max_borrow_cores", None)
    on_preempt = getattr(args, "on_preempt", None)
    borrow = bool(getattr(args, "borrow", False)) or borrow_from_text is not None or maximum is not None
    if getattr(args, "exclusive_timing", False) and borrow:
        raise ValueError("--exclusive-timing is for owner jobs; a borrow job cannot request it")
    if getattr(args, "smt", False) and borrow:
        raise ValueError("--smt is for owner jobs; a borrow job cannot request it")
    if not borrow:
        if on_preempt is not None:
            raise ValueError("--on-preempt needs --borrow")
        return None
    if not config.borrow.enabled:
        raise ValueError(f"borrowing is disabled in the node configuration ({config.path})")
    if getattr(args, "cpus", None) is not None:
        raise ValueError("--cpus cannot be combined with --borrow; the scheduler places borrow jobs")
    if getattr(args, "gpus", None) or getattr(args, "gpu_mem", None) is not None:
        raise ValueError("borrow jobs cannot use GPUs")
    names: list[str] | None = None
    if borrow_from_text is not None:
        names = [name.strip() for name in borrow_from_text.split(",") if name.strip()]
        if not names:
            raise ValueError("--borrow-from needs at least one partition")
        for name in names:
            part = config.partition(name)
            if part.name != config.canonical_partition(home) and not part.lends():
                why = "is timing-sensitive" if part.timing_sensitive else "does not lend (lend: false)"
                raise ValueError(f"--borrow-from: partition {name} {why}")
    if maximum is not None and maximum < 1:
        raise ValueError("--max-borrow-cores must be positive")
    probe = Request(job_id=0, partition=home, ncpus=1, borrow=True,
                    borrow_from=None if names is None else tuple(names))
    eligible = eligible_partitions(config, probe)
    protected = protected_cpus(config)
    reachable = sum(len([c for c in config.partitions[name].cpus if c not in protected]) for name in eligible)
    if reachable == 0:
        raise ValueError(f"no partition lends cores to a borrow job of partition {home}")
    return {
        "borrow_from": None if names is None else json.dumps(names),
        "max_borrow_cores": maximum,
        "on_preempt": on_preempt or config.borrow.default_on_preempt,
        "on_preempt_given": on_preempt is not None,
        "eligible": eligible,
        "reachable": reachable,
    }


def _reject_node_options_without_config(args: argparse.Namespace) -> None:
    given = [name for name in NODE_OPTIONS if getattr(args, name, None) is not None]
    given += [name for name in NODE_FLAGS if getattr(args, name, False)]
    if given:
        rendered = ", ".join("--" + name.replace("_", "-") for name in given)
        raise ValueError(
            f"{rendered} need a node configuration ({config_path()}); this queue runs in 0.5 mode"
        )


def _describe_request(job: dict[str, Any]) -> str:
    parts = [f"{job['requested_cores']} core(s)"]
    if job.get("cpu_request"):
        parts[0] = f"CPUs {job['cpu_request']}"
    gpus = int(job.get("gpus_requested") or 0)
    if gpus:
        claim = job.get("gpu_mem_bytes")
        parts.append(f"{gpus} GPU" + (f" (shared, {format_size(int(claim))})" if claim else " (exclusive)"))
    if job.get("mem_bytes") is not None:
        parts.append(f"mem {format_size(int(job['mem_bytes']))}" + (" (auto)" if job.get("mem_auto_json") else ""))
    parts.append(f"time {format_time_limit(job.get('time_limit_s'))}")
    if int(job.get("priority") or 0):
        parts.append(f"priority {job['priority']}")
    if job.get("borrow"):
        text = "borrow"
        if job.get("max_borrow_cores") is not None:
            text += f" {job['requested_cores']}-{job['max_borrow_cores']} cores"
        if job.get("borrow_from"):
            try:
                text += " from " + ",".join(json.loads(str(job["borrow_from"])))
            except ValueError:
                pass
        text += f", on preempt {job.get('on_preempt') or 'requeue'}"
        parts.append(text)
    if smt_job(job):
        parts.append("smt: cores with their SMT siblings, exclusive timing")
    elif job.get("exclusive_timing"):
        parts.append("exclusive timing")
    if job.get("l3_domains"):
        parts.append(f"over {job['l3_domains']} L3 domains")
    if job.get("split_domains"):
        try:
            parts.append("split over L3 domains " + "+".join(str(p) for p in json.loads(str(job["split_domains"]))))
        except ValueError:
            pass
    return ", ".join(parts)


def _test_only(config: NodeConfig, store: QueueStore, requested: int, scheduling: dict[str, Any]) -> int:
    """Validate a request and report where it would run now, like
    ``sbatch --test-only``.  Nothing is queued or launched."""
    from .borrow import plan_node
    from .locks import borrow_lock_plan, busy_resources, lock_plan_for
    from .nodeworker import borrower_entries, job_request, running_entries
    from .scheduler import INFINITY, build_pool

    now = time.time()
    jobs = store.node_jobs()
    running = running_entries(jobs, config, now)
    owner_of = {cpu: int(item.job_id) for item in running if item.job_id is not None for cpu in item.allocation.cpus}
    allocated_gpus = {gpu for item in running for gpu in item.allocation.gpus}
    busy_cpus, busy_gpus = busy_resources(config, allocated_cpus=set(owner_of), allocated_gpus=allocated_gpus)
    pool = build_pool(config, running, external_cpus=busy_cpus, external_gpus=busy_gpus)
    waiting = [j for j in jobs if j["state"] == "queued"]
    probe = dict(scheduling, id=0, requested_cores=requested)
    # A new submission goes behind every queued job of its class (owner jobs
    # before borrow jobs) with equal or higher priority.
    probe_class = int(scheduling.get("borrow") or 0)
    position = sum(
        1 for j in waiting
        if int(j.get("borrow") or 0) < probe_class
        or (int(j.get("borrow") or 0) == probe_class
            and int(j.get("priority") or 0) >= int(scheduling["priority"] or 0))
    )
    waiting.insert(position, probe)
    queue = [job_request(j, config) for j in waiting]
    borrowers = borrower_entries([j for j in jobs if j["state"] in ACTIVE_STATES], config)
    # Without the worker's usage history, leased and allocated cores count
    # as busy here; only free cores are shown as lendable.
    decision = plan_node(config, pool, running, queue, now, borrowers=borrowers, owner_of=owner_of)
    started = {request.job_id: allocation for request, allocation in decision.starts}
    payload: dict[str, Any] = {
        "schema_version": "simqueue-test-only-v1",
        "partition": scheduling["partition"],
        "requested_cores": requested,
        "request": {k: scheduling[k] for k in scheduling},
        "queued_ahead": position,
    }
    if smt_job(probe):
        payload["smt"] = True
    if 0 in started:
        allocation = started[0]
        locks = (
            borrow_lock_plan(config, allocation.cpus)
            if allocation.borrow
            else lock_plan_for(config, allocation.cpus, allocation.gpus, gpu_shared=not allocation.gpu_exclusive)
        )
        payload.update(
            {
                "would_start_now": True,
                "cpu_list": format_cpu_list(allocation.cpus),
                "gpus": list(allocation.gpus),
                "gpu_mode": "none" if not allocation.gpus else ("exclusive" if allocation.gpu_exclusive else "shared"),
                "lease_files": [{"path": lock.path, "mode": lock.mode} for lock in locks],
            }
        )
        if smt_job(probe):
            affinity, threads = config.smt_affinity(allocation.cpus)
            payload.update({"cpu_list_smt": format_cpu_list(affinity), "threads_per_core": threads})
        if allocation.borrow:
            payload["borrow_placement"] = [
                {"cpu": cpu, "partition": part, "class": cls} for cpu, part, cls in allocation.placement
            ]
        preempted = [item.job_id for item in decision.preemptions]
        if preempted:
            payload["would_preempt_borrow_jobs"] = preempted
    else:
        payload["would_start_now"] = False
        payload["reason"] = decision.skipped.get(0, "waiting behind earlier jobs")
        if busy_cpus:
            payload["leases_held_outside_queue"] = sorted(busy_cpus)
    print(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


def _enqueue(args: argparse.Namespace, store: QueueStore) -> int | tuple[dict[str, Any], dict[str, Any], str]:
    """Validate and queue a command submission (``submit`` and ``run``).

    Returns the job, the submission payload and the placement text, or the
    exit status of ``--test-only`` (which queues nothing)."""
    node_config = load_node_config()
    scheduling: dict[str, Any] | None = None
    explicit: dict[str, Any] | None = None
    if node_config is not None:
        requested, scheduling, explicit = _node_request(args, node_config, cores=args.cores, store=store,
                                                        name=_submission_name(args))
        pcores: tuple[int, ...] = ()
        if getattr(args, "test_only", False):
            _resolve_cwd(args.cwd)
            return _test_only(node_config, store, requested, scheduling)
    else:
        _reject_node_options_without_config(args)
        if getattr(args, "test_only", False):
            raise ValueError("--test-only needs a node configuration")
        pcores = _cli_pcores()
        requested = args.cores if args.cores is not None else len(pcores)
        if requested < 1 or requested > len(pcores):
            raise ValueError(
                f"--cores must be between 1 and {len(pcores)} on this host"
            )
    command = list(args.command)
    if command and command[0] == "--":
        command.pop(0)
    if args.shell is not None:
        if command:
            raise ValueError("use either --shell COMMAND or -- COMMAND, not both")
        argv = ["/bin/bash", "-lc", args.shell]
        try:
            shell_words = shlex.split(args.shell)
        except ValueError:
            shell_words = []
        default_name = Path(shell_words[0]).name if shell_words else "shell-job"
    else:
        if not command:
            raise ValueError("missing command; use -- COMMAND [ARGS...] or --shell COMMAND")
        argv = command
        default_name = Path(argv[0]).name or argv[0]
    name = args.name or default_name
    if not name.strip():
        raise ValueError("job name cannot be empty")
    cwd = _resolve_cwd(args.cwd)
    environment, secret_env, override_keys, withheld = _capture_environment(args.env, getattr(args, "export", None))
    if withheld:
        print(f"simq: not passed to the job (secret-looking names; use -e NAME=VALUE): {', '.join(withheld)}",
              file=sys.stderr)
    # 0.7.0: GitHub reporting targets (posted by the worker on transitions)
    github = targets_from_args(getattr(args, "github_status", None), getattr(args, "github_comment", None),
                               getattr(args, "github_summary", None), cwd)
    if github is not None:
        if node_config is None:
            raise ValueError("--github-status, --github-comment and --github-summary need a node configuration")
        if not node_config.github.enabled:
            raise ValueError(f"GitHub reporting is disabled in the node configuration ({node_config.path})")
        explicit = dict(explicit or {}, github=github)
    basis = getattr(args, "mem_basis", None)
    job = store.submit(
        name=name,
        cwd=cwd,
        argv=argv,
        environment=environment,
        env_override_keys=override_keys,
        requested_cores=requested,
        idempotency_key=args.idempotency_key,
        scheduling=scheduling,
        explicit_scheduling=explicit,
        secret_env=secret_env,
        extra_columns={
            "github_json": None if github is None else json.dumps(github, sort_keys=True),
            "mem_auto_json": None if basis is None else json.dumps(basis, sort_keys=True, default=str),
        },
    )
    created = bool(job.get("submission_created", True))
    cpu_list = format_cpu_list(pcores[:requested])
    placement_text = f"on P-cores {cpu_list}"
    if node_config is not None:
        placement_text = f"in partition {_partition_text(node_config, job.get('partition'))}: {_describe_request(job)}"
    payload = {
        "schema_version": "simqueue-submit-result-v1",
        "kind": "command",
        "job_id": int(job["id"]),
        "created": created,
        "name": str(job["name"]),
        "state": str(job["state"]),
        "requested_cores": int(job["requested_cores"]),
        "cpu_list": job.get("cpu_list"),
        "log_path": str(job["log_path"]),
        "idempotency_key": job.get("idempotency_key"),
        "request_sha256": job.get("request_sha256"),
    }
    if node_config is not None:
        payload.update(
            {
                "partition": job.get("partition"),
                "cpu_request": job.get("cpu_request"),
                "gpus_requested": int(job.get("gpus_requested") or 0),
                "gpu_mem_bytes": job.get("gpu_mem_bytes"),
                "mem_bytes": job.get("mem_bytes"),
                "time_limit_s": job.get("time_limit_s"),
                "priority": int(job.get("priority") or 0),
            }
        )
        if job.get("borrow"):
            payload.update(
                {
                    "borrow": True,
                    "borrow_from": None if not job.get("borrow_from") else json.loads(str(job["borrow_from"])),
                    "max_borrow_cores": job.get("max_borrow_cores"),
                    "on_preempt": job.get("on_preempt"),
                }
            )
        if job.get("exclusive_timing"):
            payload["exclusive_timing"] = True
        if smt_job(job):
            payload["smt"] = True
        if basis is not None:
            payload["mem_auto"] = basis
        ceiling = _admissible_memory(node_config, store)
        mem_now = job.get("mem_bytes")
        if ceiling is not None and mem_now is not None and int(mem_now) > ceiling[0]:
            note = (f"mem {format_size(int(mem_now))} exceeds what {node_config.node} can admit now "
                    f"(~{format_gib(ceiling[0])}: {ceiling[1]}); the job is parked until that much is free and "
                    "does not hold up other jobs")
            payload.setdefault("warnings", []).append(note)
            print(f"simq: warning: {note}", file=sys.stderr)
        if github is not None:
            payload["github"] = github
        warning = getattr(args, "mem_warning", None)
        if warning:
            payload["warnings"] = [warning]
            print(f"simq: warning: {warning}", file=sys.stderr)
    return job, payload, placement_text


def _admissible_memory(config: NodeConfig, store: QueueStore) -> tuple[int, str] | None:
    """0.7.1: the largest request the memory guard can admit on this node
    now: MemAvailable with the queue's running jobs taken away, minus the
    margin, and memory_schedulable minus the margin.  (The worker also
    remembers the best value over the last hour.)"""
    from .memory import admission_margin

    if not config.memory.enabled:
        return None
    currents = []
    for job in store.active_jobs():
        usage = _live_usage(job, config)
        if usage is not None:
            currents.append(usage.current)
    node = read_node_memory(config.memory.test_node_total, currents)
    if node is None:
        return None
    margin = admission_margin(config.memory, node.total)
    idle = node.available + sum(currents)
    ceiling = idle - margin
    text = f"MemAvailable with the queue idle ~{format_gib(idle)} - margin {format_gib(margin)}"
    if config.memory_schedulable is not None and config.memory_schedulable - margin < ceiling:
        ceiling = config.memory_schedulable - margin
        text = f"memory_schedulable {format_gib(config.memory_schedulable)} - margin {format_gib(margin)}"
    return ceiling, text


def _partition_text(config: NodeConfig, name: str | None) -> str:
    """``name``, or ``alias (pool partition)`` for an alias (0.6.3)."""
    recorded = name or config.default_partition
    canonical = config.canonical_partition(recorded)
    return recorded if canonical == recorded else f"{recorded} (pool {canonical})"


def _submit(args: argparse.Namespace, store: QueueStore) -> int:
    submitted = _enqueue(args, store)
    if isinstance(submitted, int):
        return submitted
    job, payload, placement_text = submitted
    created = bool(payload["created"])
    try:
        mode = _nudge_worker(store) if job["state"] not in TERMINAL_STATES else "terminal"
    except WorkerStartError as error:
        payload["worker_mode"] = "error"
        payload["worker_error"] = str(error)
        if args.json:
            print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
        else:
            action = "Queued" if created else "Existing"
            print(
                f"{action} job {job['id']} ({job['name']}) {placement_text}.\n"
                f"Log: {job['log_path']}"
            )
        print(
            f"simq: job {job['id']} remains queued, but the protected runner "
            f"did not start: {error}",
            file=sys.stderr,
        )
        return 1
    payload["worker_mode"] = mode
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
    else:
        action = "Queued" if created else "Existing"
        print(
            f"{action} job {job['id']} ({job['name']}) {placement_text}.\n"
            f"Log: {job['log_path']}"
        )
        if mode == "disabled":
            print("Worker autostart is disabled; run `simq worker --drain` to execute it.")
        elif mode == "exclusive-disabled":
            print(
                "Worker autostart is disabled; start the protected runner with "
                "`systemctl start simqueue-runner.service`."
            )
        elif mode == "paused":
            print("Queue is paused; the protected runner will start on `simq resume`.")
        elif mode == "ensemble-paused":
            print(
                "An ensemble is paused at the dispatch head; this job remains "
                "queued until that ensemble is resumed or reordered."
            )
    return 0


# ----- simq run (0.6.3) ---------------------------------------------------

RUN_POLL_SECONDS = 0.2
RUN_FINISH_GRACE_SECONDS = 2.0  # wait this long for the worker's last log lines
RUN_CANCEL_WAIT_SECONDS = 20.0  # after a signal: wait this long for the job to end
RUN_SIGNALS = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)


def _say(text: str) -> None:
    print(f"simq: {text}", file=sys.stderr, flush=True)


def _request_cancel(store: QueueStore, job_id: int) -> str:
    """Cancel a job as ``simq cancel`` does; returns the state afterwards."""
    state = store.cancel(job_id)
    if state == "cancelling":
        job = store.require(job_id)
        node_config = load_node_config()
        if node_config is not None and job.get("launch_backend") == "systemd-user" and job.get("unit_name"):
            # the job unit is independent of the worker: stop it directly
            stop_unit(node_config.launch, str(job["unit_name"]))
        try:
            _nudge_worker(store)
        except WorkerStartError as error:
            _say(f"cancellation is recorded, but the protected runner did not start: {error}")
    return state


def job_exit_status(job: dict[str, Any]) -> int:
    """The exit status ``simq run`` returns for a finished job: the command's
    own exit code when recorded (a signal N as 128+N, as a shell reports
    it), else 0 for ``succeeded`` and 1 otherwise.  A job that did not
    succeed never returns 0."""
    code = job.get("exit_code")
    succeeded = job.get("state") == "succeeded"
    if code is None:
        return 0 if succeeded else 1
    code = int(code)
    status = min(128 - code if code < 0 else code, 255)
    if status == 0 and not succeeded:
        return 1
    return status


class _LogFollower:
    """Copies a job log to stdout as it grows (bytes, unchanged)."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.offset = 0
        self.tail = b""
        self.stdout_open = True

    def pump(self) -> None:
        try:
            with self.path.open("rb") as handle:
                handle.seek(self.offset)
                data = handle.read()
        except OSError:
            return
        if not data:
            return
        self.offset += len(data)
        self.tail = (self.tail + data)[-512:]
        if not self.stdout_open:
            return
        try:
            sys.stdout.buffer.write(data)
            sys.stdout.buffer.flush()
        except (BrokenPipeError, ValueError, OSError):
            # nobody reads the output any more; keep waiting for the job
            self.stdout_open = False


def _follow_job(store: QueueStore, job_id: int, received: list[int]) -> int:
    """Stream a job's log until it ends; cancel it once a signal is in
    *received* (filled by the SIGINT/SIGTERM/SIGHUP handlers of ``_run``).

    Returns the job's exit status (``job_exit_status``), or 128+signal when
    a signal canceled it."""
    follower = _LogFollower(Path(str(store.require(job_id)["log_path"])))
    marker = f"===== simq job {job_id} finished =====".encode()
    canceled_by: int | None = None
    cancel_deadline = 0.0
    terminal_since: float | None = None
    last_state = None
    while True:
        follower.pump()
        if received and canceled_by is None:
            canceled_by = received[0]
            name = signal.Signals(canceled_by).name
            state = _request_cancel(store, job_id)
            _say(f"{name}: canceling job {job_id} ({state})")
            cancel_deadline = time.monotonic() + RUN_CANCEL_WAIT_SECONDS
        elif canceled_by is not None and len(received) > 1:
            _say(f"second signal: not waiting for job {job_id} to end (its cancellation is recorded)")
            break
        job = store.require(job_id)
        state = str(job["state"])
        if state != last_state and last_state is not None and state == "queued":
            _say(f"job {job_id} was requeued (see `simq show {job_id}`)")
        last_state = state
        if state in TERMINAL_STATES:
            if terminal_since is None:
                terminal_since = time.monotonic()
            if marker in follower.tail or time.monotonic() - terminal_since >= RUN_FINISH_GRACE_SECONDS:
                follower.pump()
                break
        elif canceled_by is not None and time.monotonic() > cancel_deadline:
            _say(f"job {job_id} is still {state}; its cancellation is recorded")
            break
        time.sleep(RUN_POLL_SECONDS)
    job = store.require(job_id)
    if canceled_by is not None:
        return 128 + int(canceled_by)
    status = job_exit_status(job)
    detail = f"exit {job['exit_code']}" if job.get("exit_code") is not None else "no exit code"
    _say(f"job {job_id} {job['state']} ({detail}); exit status {status}")
    return status


def _run(args: argparse.Namespace, store: QueueStore) -> int:
    """``simq run``: submit, stream the job's output, wait, return its status
    (Slurm ``srun``-shaped; the command runs in its own job unit as for
    ``submit``, stdin is not forwarded).  The signal handlers are in place
    before the job is queued, so a signal during the submission cancels the
    job right after it."""
    received: list[int] = []

    def on_signal(signum: int, _frame: Any) -> None:
        received.append(signum)

    previous = {sig: signal.signal(sig, on_signal) for sig in RUN_SIGNALS}
    try:
        return _run_with_handlers(args, store, received)
    finally:
        for sig, handler in previous.items():
            signal.signal(sig, handler)


def _run_with_handlers(args: argparse.Namespace, store: QueueStore, received: list[int]) -> int:
    submitted = _enqueue(args, store)
    if isinstance(submitted, int):
        return submitted  # --test-only
    job, payload, placement_text = submitted
    job_id = int(job["id"])
    created = bool(payload["created"])
    try:
        mode = _nudge_worker(store) if job["state"] not in TERMINAL_STATES else "terminal"
    except WorkerStartError as error:
        payload["worker_mode"] = "error"
        payload["worker_error"] = str(error)
        _say(json.dumps(payload, ensure_ascii=False, sort_keys=True) if args.json
             else f"job {job_id} ({job['name']}) remains queued, but the protected runner did not start: {error}")
        return 1
    payload["worker_mode"] = mode
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, sort_keys=True), file=sys.stderr, flush=True)
    elif created:
        _say(f"job {job_id} ({job['name']}) queued {placement_text}")
    else:
        _say(f"attached to existing job {job_id} ({job['name']}, {job['state']}) {placement_text}")
    if mode == "disabled":
        _say("worker autostart is disabled; waiting for a worker to run the job")
    elif mode == "paused":
        _say("the queue is paused; waiting until it is resumed")
    elif mode == "ensemble-paused":
        _say("an ensemble is paused at the dispatch head; waiting until it is resumed or reordered")
    return _follow_job(store, job_id, received)


def _ensemble_manifest(args: argparse.Namespace) -> dict[str, Any]:
    node_config = load_node_config()
    if node_config is not None:
        available = len(node_config.partition(getattr(args, "partition", None)).cpus)
    else:
        _reject_node_options_without_config(args)
        available = len(_cli_pcores())
    return load_manifest(
        Path(args.manifest),
        available_cores=available,
        name_override=args.name,
        max_parallel_override=args.max_parallel,
        gpu_budget_override=args.gpu_budget,
        gpu_device_override=args.gpu_device,
        diagnostics_interval_override=args.diagnostic_interval,
        fail_fast_override=args.fail_fast,
    )


def _ensemble_submit(args: argparse.Namespace, store: QueueStore) -> int:
    canonical = _ensemble_manifest(args)
    launcher = Path(sys.argv[0]).resolve()
    if not launcher.is_file():
        launcher = Path(__file__).resolve().parents[1] / "simq"
    members = list(canonical["members"])
    ensemble = canonical["gpu"]
    node_config = load_node_config()
    scheduling: dict[str, Any] | None = None
    explicit: dict[str, Any] | None = None
    if node_config is not None:
        # The ensemble runs inside its own allocation: max_parallel cores of
        # the partition and, when members claim GPU memory, one GPU shared
        # up to the manifest's budget.
        budget = int(ensemble["budget_bytes"])
        args.gpu_mem = f"{budget}B" if budget > 0 else None
        args.gpus = None
        args.cpus = None
        _requested, scheduling, explicit = _node_request(
            args, node_config, cores=int(canonical["max_parallel"]), store=store, name=str(canonical["name"])
        )
        # The GPU request follows from the manifest, which is already part
        # of the ensemble fingerprint.
        explicit.pop("gpu_mem_bytes", None)
        explicit.pop("gpus_requested", None)
    job = store.submit_ensemble(
        name=str(canonical["name"]),
        cwd=project_root(),
        environment=_capture_environment([])[0],
        env_override_keys=[],
        requested_cores=int(canonical["max_parallel"]),
        runner_prefix=[sys.executable, str(launcher), "_run-ensemble"],
        schema_version=str(canonical["schema_version"]),
        canonical_spec=canonical,
        members=members,
        max_parallel=int(canonical["max_parallel"]),
        gpu_device=str(ensemble["device"]),
        gpu_budget_bytes=int(ensemble["budget_bytes"]),
        diagnostics_interval=float(canonical["diagnostics_interval"]),
        mps_mode=str(ensemble["mps"]),
        mps_pipe_directory=ensemble["mps_pipe_directory"],
        fail_fast=bool(canonical["fail_fast"]),
        allow_shared_cwd=bool(canonical["allow_shared_cwd"]),
        idempotency_key=args.idempotency_key,
        scheduling=scheduling,
        explicit_scheduling=explicit,
    )
    created = bool(job.get("submission_created", True))
    if getattr(args, "mem_basis", None) is not None and created:
        store.set_memory_basis(int(job["id"]), args.mem_basis)
    stored = store.require_ensemble(int(job["id"]))
    payload = {
        "schema_version": "simqueue-submit-result-v1",
        "kind": "ensemble",
        "job_id": int(job["id"]),
        "created": created,
        "name": str(job["name"]),
        "state": str(job["state"]),
        "requested_cores": int(job["requested_cores"]),
        "cpu_list": job.get("cpu_list"),
        "log_path": str(job["log_path"]),
        "artifact_dir": str(stored["artifact_dir"]),
        "member_count": len(members),
        "max_parallel": int(canonical["max_parallel"]),
        "idempotency_key": job.get("idempotency_key"),
        "request_sha256": job.get("request_sha256"),
    }
    try:
        mode = _nudge_worker(store) if job["state"] not in TERMINAL_STATES else "terminal"
    except WorkerStartError as error:
        payload["worker_mode"] = "error"
        payload["worker_error"] = str(error)
        if args.json:
            print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
        else:
            action = "Queued" if created else "Existing"
            print(
                f"{action} ensemble {job['id']} ({job['name']}) with {len(members)} "
                f"members, up to {canonical['max_parallel']} active.\n"
                f"Log: {job['log_path']}\nArtifacts: {stored['artifact_dir']}"
            )
        print(
            f"simq: ensemble {job['id']} remains queued, but the protected runner "
            f"did not start: {error}",
            file=sys.stderr,
        )
        return 1
    payload["worker_mode"] = mode
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
    else:
        action = "Queued" if created else "Existing"
        print(
            f"{action} ensemble {job['id']} ({job['name']}) with {len(members)} members, "
            f"up to {canonical['max_parallel']} active.\n"
            f"Log: {job['log_path']}\n"
            f"Artifacts: {stored['artifact_dir']}"
        )
        if mode == "disabled":
            print("Worker autostart is disabled; run `simq worker --drain` to execute it.")
        elif mode == "exclusive-disabled":
            print(
                "Worker autostart is disabled; start the protected runner with "
                "`systemctl start simqueue-runner.service`."
            )
        elif mode == "paused":
            print("Queue is paused; the protected runner will start on `simq resume`.")
        elif mode == "ensemble-paused":
            print(
                "An ensemble is paused at the dispatch head; this ensemble remains "
                "queued behind it."
            )
    return 0


def _ensemble_validate(args: argparse.Namespace) -> int:
    canonical = _ensemble_manifest(args)
    if args.json:
        print(json.dumps(canonical, ensure_ascii=False, indent=2, sort_keys=True))
    else:
        print(
            f"Valid {canonical['schema_version']} manifest: {canonical['name']}  "
            f"members={len(canonical['members'])}  "
            f"max_parallel={canonical['max_parallel']}  "
            f"gpu_budget_bytes={canonical['gpu']['budget_bytes']}  "
            f"mps={canonical['gpu']['mps']}"
        )
    return 0


def _ensemble_status(args: argparse.Namespace, store: QueueStore) -> int:
    payload = ensemble_status_payload(args.job_id, store)
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True))
    else:
        print_ensemble_status(payload)
    return 0


def _ensemble_pause(args: argparse.Namespace, store: QueueStore) -> int:
    before = store.ensemble_control_snapshot(args.job_id)
    control = store.set_ensemble_paused(args.job_id, True)
    changed = int(control["control_generation"]) != int(
        before["control_generation"]
    )
    payload = {
        "schema_version": "simqueue-ensemble-control-v1",
        "action": "pause",
        "changed": changed,
        "job_id": args.job_id,
        **control,
    }
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
    elif control["job_state"] in ("starting", "running"):
        prefix = "Pause requested for" if changed else "Pause already active for"
        print(
            f"{prefix} ensemble {args.job_id}; {control['running']} admitted "
            f"member(s) may finish and {control['queued']} queued member(s) "
            "will not start."
        )
        print(
            "The allocation will be released automatically after the admitted "
            "members drain."
        )
    else:
        prefix = "Paused" if changed else "Already paused"
        print(
            f"{prefix} ensemble {args.job_id}; {control['queued']} member(s) "
            "remain queued under the same job ID."
        )
    return 0


def _ensemble_resume(args: argparse.Namespace, store: QueueStore) -> int:
    before = store.ensemble_control_snapshot(args.job_id)
    control = store.set_ensemble_paused(args.job_id, False)
    changed = int(control["control_generation"]) != int(
        before["control_generation"]
    )
    # An active supervisor polls the durable control state itself, while a
    # process-tree-quiet drained parent needs a worker to reclaim it. Nudging is
    # safe in both cases: the worker lock and protected systemd unit serialize
    # ownership, and it also closes the supervisor-exit/worker-crash window.
    mode = _nudge_worker(store)
    payload = {
        "schema_version": "simqueue-ensemble-control-v1",
        "action": "resume",
        "changed": changed,
        "job_id": args.job_id,
        "worker_mode": mode,
        **control,
    }
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
    else:
        prefix = "Resumed" if changed else "Already enabled"
        print(
            f"{prefix} ensemble {args.job_id}; {control['queued']} queued "
            "member(s) will continue under the same job ID."
        )
        if mode == "paused":
            print(
                "The global queue is still paused; run `simq resume` when you "
                "want outer-job dispatch to continue."
            )
        elif mode == "ensemble-paused":
            print(
                "An ensemble is paused at the dispatch head; resume or reorder "
                "it before this one can dispatch."
            )
        elif mode in ("disabled", "exclusive-disabled"):
            print("Worker autostart is disabled; start the queue worker to continue.")
    return 0


def _reorder_job(
    args: argparse.Namespace, store: QueueStore, *, direction: str
) -> int:
    result = (
        store.promote(args.job_id)
        if direction == "front"
        else store.demote(args.job_id)
    )
    try:
        mode = _nudge_worker(store)
        worker_error: str | None = None
    except WorkerStartError as error:
        mode = "error"
        worker_error = str(error)
    action = "promote" if direction == "front" else "demote"
    payload = {
        "schema_version": "simqueue-reorder-result-v1",
        "action": action,
        "changed": bool(result["changed"]),
        "job_id": int(result["job_id"]),
        "name": str(result["job"]["name"]),
        "state": str(result["job"]["state"]),
        "queue_position": int(result["queue_position"]),
        "queue_order": int(result["job"]["queue_order"]),
        "queued_count": int(result["queued_count"]),
        "worker_mode": mode,
    }
    if worker_error is not None:
        payload["worker_error"] = worker_error
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
    else:
        if not result["changed"]:
            endpoint = "front" if direction == "front" else "back"
            print(
                f"Job {args.job_id} ({result['job']['name']}) is already at the "
                f"{endpoint}; "
                f"queue position {result['queue_position']} of {result['queued_count']}."
            )
        else:
            verb = "Promoted" if direction == "front" else "Demoted"
            print(
                f"{verb} job {args.job_id} ({result['job']['name']}) to queue "
                f"position {result['queue_position']} of {result['queued_count']}."
            )
        if mode == "paused":
            print("The global queue remains paused; ordering is saved for later dispatch.")
        elif mode == "ensemble-paused":
            print("A paused ensemble remains at the dispatch front.")
    if worker_error is not None:
        print(
            f"simq: queue order was saved, but the protected runner did not start: "
            f"{worker_error}",
            file=sys.stderr,
        )
        return 1
    return 0


def _show(job_id: int, store: QueueStore) -> int:
    job = store.require(job_id)
    ensemble = store.get_ensemble(job_id)
    argv = list(json.loads(job["argv_json"]))
    keys = list(json.loads(job["env_override_keys"]))
    fields = [
        ("ID", job["id"]),
        ("Type", "ensemble" if ensemble is not None else "command"),
        ("Name", job["name"]),
        ("State", job["state"]),
        ("Attempt", job["attempt"]),
        ("Retry of", job["retry_of"] or "-"),
        ("Directory", job["cwd"]),
        ("Command", shlex.join(argv)),
        ("Requested cores", job["requested_cores"]),
        ("CPU list", job["cpu_list"] or "not assigned"),
        ("Created", job["created_at"]),
        ("Started", job["started_at"] or "-"),
        ("Finished", job["finished_at"] or "-"),
        ("Exit code", job["exit_code"] if job["exit_code"] is not None else "-"),
        ("PID", job["pid"] or "-"),
        ("Environment overrides", ", ".join(keys) if keys else "none"),
        ("Log", job["log_path"]),
        ("Message", job["message"] or "-"),
    ]
    if job.get("partition") or job.get("launch_backend"):
        partition_text = str(job.get("partition") or "-")
        try:
            show_config = load_node_config()
        except ConfigError:
            show_config = None
        if show_config is not None and job.get("partition"):
            partition_text = _partition_text(show_config, str(job["partition"]))
        fields.extend(
            [
                ("Partition", partition_text),
                ("Priority", job.get("priority") or 0),
                ("Request", _describe_request(job)),
                ("GPUs", job.get("gpu_list") or "none"),
                ("GPU claim", format_size(job["gpu_claim_bytes"]) if job.get("gpu_list") else "-"),
                ("Memory limit", format_size(job.get("mem_bytes"))),
                *_memory_fields(job, show_config),
                ("Time limit", format_time_limit(job.get("time_limit_s"))),
                ("Unit", job.get("unit_name") or "-"),
                ("Result", job.get("result") or "-"),
                ("Adopted", job.get("adopted_at") or "-"),
                ("Job files", job.get("job_dir") or "-"),
                *_containment_fields(job),
            ]
        )
    if job.get("borrow") or job.get("exclusive_timing"):
        placement_text = "-"
        if job.get("placement_json"):
            try:
                placement = json.loads(str(job["placement_json"]))
                placement_text = "; ".join(
                    f"CPU {item['cpu']} {item.get('partition')} ({item.get('class')}"
                    + (f", job {item['owner_job']}" if item.get("owner_job") is not None else "")
                    + ")"
                    for item in placement.get("cpus", [])
                ) or "-"
            except (ValueError, KeyError, TypeError):
                placement_text = str(job["placement_json"])
        if job.get("borrow"):
            borrow_from = "every lending partition"
            if job.get("borrow_from"):
                try:
                    borrow_from = ", ".join(json.loads(str(job["borrow_from"])))
                except ValueError:
                    borrow_from = str(job["borrow_from"])
            fields.extend(
                [
                    ("Class", "borrow job (preemptible; idle CPU priority)"),
                    ("Borrow from", f"{job.get('partition')} and {borrow_from}"),
                    ("Cores", f"{job['requested_cores']}" + (
                        f" to {job['max_borrow_cores']}" if job.get("max_borrow_cores") is not None else "")),
                    ("On preempt", job.get("on_preempt") or "requeue"),
                    ("Placement", placement_text),
                    ("Frozen since", job.get("frozen_at") or "-"),
                    ("Preemption pending", job.get("preempt_pending") or "-"),
                    ("Preemptions", job.get("preempt_count") or 0),
                ]
            )
        elif smt_job(job):
            fields.append(("Exclusive timing", "yes (--smt): its cores and their SMT siblings are never lent"))
            fields.append(("SMT", _smt_text(job)))
        else:
            fields.append(("Exclusive timing", "yes: its cores and their SMT siblings are never lent"))
    if job["state"] == "queued":
        queued = store.queued_jobs()
        position = next(
            index
            for index, queued_job in enumerate(queued, start=1)
            if int(queued_job["id"]) == job_id
        )
        fields.append(("Queue position", f"{position} of {len(queued)}"))
        if int(job.get("gpus_requested") or 0) > 0:
            fields.append(("GPU line", _gpu_line_text(job_id) or "not waiting for the GPU alone "
                                                                  "(other resources busy, or it starts directly)"))
    if job.get("github_json"):
        fields.extend(_github_fields(job))
    if ensemble is not None:
        members = store.list_ensemble_members(job_id)
        fields.extend(
            [
                ("Ensemble members", len(members)),
                ("Max parallel", ensemble["max_parallel"]),
                ("GPU budget bytes", ensemble["gpu_budget_bytes"]),
                ("MPS mode", ensemble["mps_mode"]),
                ("Member launches paused", bool(ensemble["launch_paused"])),
                ("Control generation", ensemble["control_generation"]),
                ("Drained generation", ensemble["drained_generation"] or "-"),
                ("Artifacts", ensemble["artifact_dir"]),
            ]
        )
    width = max(len(label) for label, _ in fields)
    for label, value in fields:
        print(f"{label:<{width}}  {value}")
    events = store.job_events(job_id)
    if events:
        print("History:")
        for event in events:
            detail = event.get("detail") or {}
            text = ", ".join(f"{key}={value}" for key, value in sorted(detail.items()) if value not in (None, "", []))
            print(f"  {event['at']}  {event['kind']:<14} {text}")
    return 0


def _memory_fields(job: dict[str, Any], config: NodeConfig | None) -> list[tuple[str, Any]]:
    """Memory and CPU use of a job for ``simq show`` (0.7.0)."""
    fields: list[tuple[str, Any]] = []
    request = job.get("mem_bytes")
    basis = None
    if job.get("mem_auto_json"):
        try:
            basis = json.loads(str(job["mem_auto_json"]))
        except ValueError:
            basis = None
    if basis is not None:
        fields.append(("Memory request", describe_basis(basis)))
    usage = _live_usage(job, config)
    peak = job.get("mem_peak_bytes")
    if usage is not None:
        percent = f" ({100.0 * usage.current / int(request):.0f}% of the limit)" if request else ""
        fields.append(("Memory now", f"{format_gib(usage.current)}{percent}, read from the unit's cgroup"))
        if usage.peak is not None:
            peak = max(int(peak or 0), int(usage.peak))
    elif job["state"] in ACTIVE_STATES and job.get("mem_current_bytes") is not None:
        fields.append(("Memory now", f"{format_gib(int(job['mem_current_bytes']))} "
                                     f"(sampled {job.get('mem_sampled_at') or '?'})"))
    if peak is not None:
        source = {"unit": "memory.peak recorded by the unit at its end",
                  "sampled": "highest memory.peak the worker sampled"}.get(str(job.get("mem_peak_source") or ""), "")
        if usage is not None:
            source = "memory.peak of the running unit"
        percent = f", {100.0 * int(peak) / int(request):.0f}% of the limit" if request else ""
        fields.append(("Memory peak", f"{format_gib(int(peak))}{percent}" + (f" ({source})" if source else "")))
    elif job["state"] in TERMINAL_STATES and job.get("launch_backend") == "systemd-user":
        fields.append(("Memory peak", "not measured (job ended before SimQueue 0.7.0)"))
    cpu = usage.cpu_usec if usage is not None and usage.cpu_usec is not None else job.get("cpu_usage_usec")
    if cpu is not None:
        start = _parse_time(job.get("started_at"))
        end = _parse_time(job.get("finished_at")) or datetime.now(UTC)
        average = ""
        if start is not None and (end - start).total_seconds() > 0:
            cores_used = int(cpu) / 1e6 / (end - start).total_seconds()
            average = f", {cores_used:.2f} cores on average"
            allocated = len(parse_cpu_list(str(job.get("cpu_list") or ""))) if job.get("cpu_list") else 0
            if allocated:
                share = cores_used / allocated
                average += f" = {100 * share:.0f}% of its {allocated} allocated core(s)"
                if smt_job(job):
                    average += " (with their SMT siblings: can exceed 100%)"
                if share < 0.5 and (end - start).total_seconds() >= 60:
                    average += " (UNDER-USED: fewer cores would do)"
        fields.append(("CPU time", f"{int(cpu) / 1e6:.1f} s{average}"))
    series = Series.from_json(job.get("mem_series_json"))
    if series.buckets:
        low = min(b[1] for b in series.buckets)
        high = max(b[2] for b in series.buckets)
        fields.append(("Memory history",
                       f"{sparkline(series)}  {low}-{high} MiB over {len(series.buckets)} x {series.step} s"))
    return fields


def _smt_text(job: dict[str, Any]) -> str:
    """0.8.2: an --smt job's CPU affinity for ``simq show``: the unit's, from
    its spec, once it has started."""
    spec = None
    if job["state"] != "queued" and job.get("job_dir"):
        spec = _read_json_file(Path(str(job["job_dir"])) / "spec.json")
    variables = "SIMQ_SMT=1, SIMQ_THREADS_PER_CORE, SIMQ_CPU_LIST_SMT"
    if isinstance(spec, dict) and spec.get("affinity_cpu_list"):
        return f"yes: {spec.get('smt_text') or 'CPU affinity ' + str(spec['affinity_cpu_list'])} ({variables})"
    if isinstance(spec, dict):
        # started by a SimQueue without --smt (0.8.1 after a rollback)
        return (f"requested, but started by SimQueue {spec.get('scheduler_version') or '?'}, which has no --smt: "
                f"CPU affinity {job.get('cpu_list') or '?'} (its cores only); its siblings stay unlent")
    try:
        config = load_node_config()
    except ConfigError:
        config = None
    threads = ""
    if config is not None:
        partition = config.find_partition(job.get("partition"))
        if partition is not None:
            threads = f", {config.smt_affinity(partition.cpus)[1]} threads per core on {config.node}"
    return f"yes: once started, its CPU affinity is its cores and their SMT siblings{threads} ({variables})"


def _containment_fields(job: dict[str, Any]) -> list[tuple[str, Any]]:
    """0.7.2: the job's rendezvous directory and processes seen outside its cgroup."""
    from .containment import describe_escaped, read_escaped

    fields: list[tuple[str, Any]] = []
    job_dir = job.get("job_dir")
    if not job_dir:
        return fields
    spec = _read_json_file(Path(str(job_dir)) / "spec.json")
    if isinstance(spec, dict) and spec.get("tmpdir"):
        live = job["state"] in ("starting", "running", "cancelling")
        fields.append(("Temporary dir", f"{spec['tmpdir']} (TMPDIR and the MPI rendezvous; "
                       + ("removed when the job ends)" if live else "removed at the job's end)")))
    escaped = read_escaped(str(job_dir))
    if escaped:
        fields.append(("ESCAPED", f"{len(escaped)} process(es) carrying this job's launch token ran outside its "
                       f"cgroup: {describe_escaped(escaped)} (escaped.json)"))
    return fields


def _read_json_file(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _github_fields(job: dict[str, Any]) -> list[tuple[str, Any]]:
    try:
        targets = json.loads(str(job["github_json"]))
    except ValueError:
        return [("GitHub", "unreadable targets")]
    fields: list[tuple[str, Any]] = []
    for item in targets.get("statuses") or []:
        fields.append(("GitHub status", f"{item['repo']}@{item['sha'][:12]} {item['context']}"))
    if targets.get("comment"):
        comment = targets["comment"]
        fields.append(("GitHub comment", f"{comment['repo']} PR {comment['pr']}"
                       + (f" (comment {job['github_comment_id']})" if job.get("github_comment_id") else "")))
    if targets.get("summary"):
        fields.append(("Summary file", targets["summary"]))
    phase = str(job.get("github_phase") or "")
    current = job_phase(job)
    if phase == current:
        posted = f"{current.split(':', 1)[0]} posted at {job.get('github_posted_at') or '-'}"
    elif job.get("github_attempts"):
        posted = f"{current.split(':', 1)[0]} not posted yet ({job['github_attempts']} failed attempt(s))"
    else:
        posted = f"{current.split(':', 1)[0]} due"
    fields.append(("GitHub posted", posted))
    if job.get("github_error"):
        fields.append(("GitHub error", job["github_error"]))
    return fields


def _gpu_line_text(job_id: int) -> str | None:
    try:
        config = load_node_config()
    except ConfigError:
        return None
    if config is None or not config.line.enabled:
        return None
    from . import gpuline
    from .gpuqueue import describe_job

    state = gpuline.load_queue_state(config.line.state_path)
    text = describe_job(state, job_id)
    if text is None:
        return None
    for uuid, gpu in config.gpus.items():
        for item in (((state or {}).get("gpus") or {}).get(uuid) or {}).get("queue") or []:
            if int(item.get("job_id") or -1) == job_id:
                waited = item.get("ready_since_unix")
                if waited is not None and time.time() - float(waited) > config.line.starvation_warn_s:
                    text += (f"; WARNING: waiting longer than {gpuline.duration_text(config.line.starvation_warn_s)}")
    return text


def _tail(job_id: int, store: QueueStore, *, lines: int, follow: bool) -> int:
    job = store.require(job_id)
    path = Path(str(job["log_path"]))
    if follow:
        try:
            return subprocess.call(["tail", "-n", str(lines), "-f", str(path)])
        except KeyboardInterrupt:
            return 130
    result = subprocess.run(["tail", "-n", str(lines), str(path)], check=False)
    return result.returncode


def _wait(job_ids: list[int], store: QueueStore) -> int:
    ids = job_ids or store.unfinished_ids()
    if not ids:
        print("No unfinished jobs.")
        return 0
    for job_id in ids:
        store.require(job_id)
    remaining = set(ids)
    final: dict[int, str] = {}
    try:
        while remaining:
            for job_id in list(remaining):
                state = str(store.require(job_id)["state"])
                if state in TERMINAL_STATES:
                    final[job_id] = state
                    remaining.remove(job_id)
                    print(f"Job {job_id}: {state}")
            if remaining:
                time.sleep(1.0)
    except KeyboardInterrupt:
        print("Wait interrupted; jobs were not canceled.", file=sys.stderr)
        return 130
    return 0 if all(state == "succeeded" for state in final.values()) else 1


def _systemctl_property(
    property_name: str,
    *,
    unit: str = "simqueue.service",
    user: bool = True,
) -> str:
    command = ["systemctl"]
    if user:
        command.append("--user")
    command.extend(
        ["show", unit, f"--property={property_name}", "--value"]
    )
    try:
        result = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    except OSError:
        return "unavailable"
    if result.returncode != 0:
        return "not installed"
    return result.stdout.strip() or "unknown"


def _doctor(store: QueueStore) -> int:
    node_config = load_node_config()
    if node_config is not None:
        return _node_doctor(node_config, store)
    try:
        exclusive = load_config()
        config_error: str | None = None
    except RuntimeError as error:
        exclusive = None
        config_error = str(error)
    try:
        detected = (
            parse_cpu_list(str(exclusive["pcores"]))
            if exclusive is not None
            else detect_pcores()
        )
        pcores = format_cpu_list(detected)
    except (RuntimeError, ValueError) as error:
        pcores = f"ERROR: {error}"

    if exclusive is not None:
        unit = str(exclusive["runner_unit"])
        load_state = _systemctl_property("LoadState", unit=unit, user=False)
        active_state = _systemctl_property("ActiveState", unit=unit, user=False)
        unit_file_state = _systemctl_property("UnitFileState", unit=unit, user=False)
        active_jobs = store.active_jobs()
        queued_jobs = store.queued_jobs()
        blocker = store.dispatch_blocker()
        runner_required = bool(active_jobs) or (bool(queued_jobs) and blocker is None)
        active_ok = active_state in ("active", "activating", "reloading")
        idle_ok = active_state == "inactive" and not runner_required
        try:
            shield = helper_status(exclusive)
            shield_state = str(shield.get("state", "unknown"))
            shield_valid = shield.get("valid") is True
            shield_detail = shield_state if shield_valid else f"ERROR: {shield_state}"
        except RuntimeError as error:
            shield_state = "error"
            shield_valid = False
            shield_detail = f"ERROR: {error}"
        if active_state in ("active", "activating", "reloading"):
            shield_matches_runner = shield_state == "root"
        else:
            shield_matches_runner = shield_state in ("member", "absent")
        service_ok = (
            load_state == "loaded"
            and unit_file_state in ("enabled", "static")
            and (active_ok or idle_ok)
            and shield_valid
            and shield_matches_runner
        )
        mode = "exclusive cgroup-v2 P-core partition"
        linger = "not required (system service)"
    else:
        unit = "simqueue.service"
        user = pwd.getpwuid(os.getuid()).pw_name
        linger = "unknown"
        try:
            result = subprocess.run(
                ["loginctl", "show-user", user, "--property=Linger", "--value"],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                check=False,
            )
            if result.returncode == 0:
                linger = result.stdout.strip()
        except OSError:
            pass
        load_state = _systemctl_property("LoadState")
        active_state = _systemctl_property("ActiveState")
        unit_file_state = _systemctl_property("UnitFileState")
        service_ok = (
            config_error is None
            and load_state == "loaded"
            and active_state == "active"
            and unit_file_state in ("enabled", "static")
        )
        mode = (
            f"ERROR: {config_error}"
            if config_error is not None
            else "cooperative affinity (not exclusive)"
        )
        shield_detail = "not installed"
    fields = [
        ("Project root", project_root()),
        ("State directory", state_dir()),
        ("Database", database_path()),
        ("Detected P-cores", pcores),
        ("Queue", "paused" if store.is_paused() else "accepting jobs"),
        ("CPU protection", mode),
        ("P-core shield", shield_detail),
        ("Runner unit", unit),
        ("Service installed", load_state),
        ("Service active", active_state),
        ("Service enabled", unit_file_state),
        ("Survives logout (linger)", linger),
    ]
    width = max(len(label) for label, _ in fields)
    for label, value in fields:
        print(f"{label:<{width}}  {value}")
    return 0 if not pcores.startswith("ERROR") and service_ok else 1


def _node_status(config: NodeConfig, store: QueueStore, *, sample_s: float = 0.0) -> dict[str, Any]:
    """Node state for ``simq partitions``.

    Every core is classified as ``free``, ``job`` (an owner queue job),
    ``borrowed`` (only borrow jobs on it), ``external-busy`` /
    ``external-idle`` (core lease held outside the queue; CPU use of the
    holder's process tree over *sample_s* seconds), ``external`` (not
    sampled) or ``stale`` (lease held but no holder process visible).
    """
    import dataclasses

    from .borrow import in_use_cores
    from .locks import busy_resources, process_description
    from .nodeworker import borrower_entries, running_entries
    from .telemetry import NodeTelemetry, command_summary

    jobs = store.node_jobs()
    now = time.time()
    running = running_entries(jobs, config, now)
    borrowers = borrower_entries([j for j in jobs if j["state"] in ACTIVE_STATES], config)
    owner: dict[int, int] = {}
    gpu_users: dict[str, list[dict[str, Any]]] = {}
    mem_used = 0
    for item in running:
        for cpu in item.allocation.cpus:
            owner[cpu] = int(item.job_id or 0)
        for uuid in item.allocation.gpus:
            gpu_users.setdefault(uuid, []).append(
                {
                    "job_id": item.job_id,
                    "claim": item.allocation.gpu_claim,
                    "exclusive": item.allocation.gpu_exclusive,
                }
            )
        mem_used += item.allocation.mem
    borrowed_mem = sum(b.mem for b in borrowers if not b.pending)
    on_cpu: dict[int, list[dict[str, Any]]] = {}
    for borrower in borrowers:
        for cpu in borrower.cpus:
            on_cpu.setdefault(cpu, []).append(
                {
                    "job_id": borrower.job_id,
                    "class": borrower.placement.get(cpu, "free"),
                    "frozen": borrower.frozen_since is not None,
                    "stopping": borrower.pending,
                }
            )
    busy_cpus, busy_gpus = busy_resources(
        config, allocated_cpus=set(owner), allocated_gpus=set(gpu_users)
    )
    groups_by_cpu: dict[int, Any] = {}
    if busy_cpus:
        window = max(0.2, sample_s) if sample_s > 0 else config.borrow.idle_window_s
        probe_config = dataclasses.replace(
            config, borrow=dataclasses.replace(config.borrow, idle_window_s=window, busy_window_s=window)
        )
        telemetry = NodeTelemetry(probe_config)
        telemetry.sample(owners=[], external=busy_cpus, borrowers=borrowers)
        if sample_s > 0:
            time.sleep(sample_s)
            telemetry.sample(owners=[], external=busy_cpus, borrowers=borrowers)
        for group in telemetry.groups:
            for cpu in group.cpus:
                groups_by_cpu[cpu] = (group, telemetry.group_trees.get(group.key, set()))
    # 0.8.0: per-CPU utilization over a short window (/proc/stat)
    from .nodeworker import read_proc_stat

    first = read_proc_stat()
    time.sleep(max(0.5, min(sample_s, 2.0)) if sample_s >= 0 else 0.5)
    second = read_proc_stat()
    utilization: dict[int, float] = {}
    if first and second:
        for cpu, (busy, total) in second["cpus"].items():
            old = first["cpus"].get(cpu)
            if old and total > old[1]:
                utilization[cpu] = max(0.0, min(1.0, (busy - old[0]) / (total - old[1])))
    no_borrow: set[int] = set()
    from .locks import lock_holders

    borrow_paths = {path: cpu for cpu in config.node_cpus() for path in config.locks.borrow_lock_paths(cpu)}
    for path, found in lock_holders(list(borrow_paths)).items():
        if any(h.mode == "WRITE" for h in found):
            no_borrow.add(borrow_paths[path])
    partitions = []
    for part in config.partitions.values():
        cpus = []
        for cpu in part.cpus:
            entry: dict[str, Any] = {"cpu": cpu, "borrowers": on_cpu.get(cpu, []), "no_borrow": cpu in no_borrow,
                                     "utilization": utilization.get(cpu),
                                     "siblings": [{"cpu": s, "utilization": utilization.get(s),
                                                   "borrowers": on_cpu.get(s, [])} for s in config.siblings(cpu)]}
            if cpu in owner:
                entry.update({"state": f"job {owner[cpu]}", "class": "job", "job_id": owner[cpu]})
            elif cpu in busy_cpus:
                entry["state"] = "lease held: " + process_description(busy_cpus[cpu][0].pid)
                entry["lock_pids"] = sorted({h.pid for h in busy_cpus[cpu] if h.pid})
                group_tree = groups_by_cpu.get(cpu)
                if group_tree is None:
                    entry["class"] = "external"
                else:
                    group, tree = group_tree
                    roots = [int(x.split(":")[0]) for x in group.key.split(",") if x and not group.stale]
                    entry["holders"] = [{"pid": pid, "command": command_summary(pid)} for pid in roots]
                    entry["tree_pids"] = len(tree)
                    entry["holder_cpus"] = list(group.cpus)
                    if group.stale:
                        entry["class"] = "stale"
                        entry["holder"] = group.description
                    elif group.usage_idle is None:
                        entry["class"] = "external"
                    else:
                        entry["usage_cores"] = round(group.usage_idle, 3)
                        entry["usage_window_s"] = sample_s
                        in_use = in_use_cores(group.usage_idle, config.borrow.idle_threshold)
                        entry["class"] = "external-idle" if in_use < len(group.cpus) else "external-busy"
            elif on_cpu.get(cpu):
                entry.update({"state": "free", "class": "borrowed"})
            else:
                entry.update({"state": "free", "class": "free"})
            cpus.append(entry)
        queued = [j for j in jobs if j["state"] == "queued" and config.canonical_partition(j.get("partition")) == part.name]
        domains = []
        for index, group in enumerate(config.topology.domains):
            cores = [cpu for cpu in group if cpu in part.cpus]
            if not cores:
                continue
            info = config.topology.info[index] if index < len(config.topology.info) else None
            users: dict[int, list[int]] = {}
            for cpu in cores:
                if cpu in owner:
                    users.setdefault(owner[cpu], []).append(cpu)
            domains.append({
                "index": index,
                "cores": cores,
                "logical_cpus": list(info.cpus) if info is not None else cores,
                "l3_bytes": None if info is None else info.size_bytes,
                "busy": sum(len(v) for v in users.values()),
                "borrowed": sum(1 for cpu in cores if cpu not in owner and on_cpu.get(cpu)),
                "jobs": {str(job_id): format_cpu_list(cpus) for job_id, cpus in sorted(users.items())},
            })
        partitions.append(
            {
                "name": part.name,
                "aliases": [describe_alias(config, name) for name in part.aliases],
                "default": part.name == config.default_partition,
                "cpus": cpus,
                "gpus": list(part.gpus),
                "default_mem": part.default_mem,
                "max_mem": part.max_mem,
                "default_time_s": part.default_time,
                "max_time_s": part.max_time,
                "queued": len(queued),
                "queued_borrow": sum(1 for j in queued if j.get("borrow")),
                "lend": part.lend,
                "lend_leased_idle": part.lend_leased_idle,
                "timing_sensitive": part.timing_sensitive,
                "l3_domains": domains,
                "l3_placement": config.topology.l3_placement,
            }
        )
    gpus = []
    from . import gpuline

    line_state = gpuline.load_queue_state(config.line.state_path) if config.line.enabled else None
    for uuid, gpu in config.gpus.items():
        users = gpu_users.get(uuid, [])
        entry = {
            "uuid": uuid,
            "index": gpu.index,
            "name": gpu.name,
            "budget": gpu.budget_bytes,
            "claimed": sum(u["claim"] for u in users),
            "jobs": users,
            "lease_held_outside_queue": [
                process_description(h.pid) for h in busy_gpus.get(uuid, [])
            ],
        }
        if config.line.enabled:
            locks = gpuline.GpuLocks(uuid, gpu.index, gpu.name, config.locks.gpu_lock_order(gpu))
            try:
                snapshot = gpuline.gpu_snapshot(config.line, locks, queue_state=line_state, now=now,
                                                scan_fds=not users)
            except Exception as error:  # noqa: BLE001 - diagnostics only
                snapshot = {"error": str(error), "holders": [], "line": [], "raw_waiters": []}
            if users:
                snapshot["holders"] = [
                    {"kind": "job", "job_id": u["job_id"], "mode": "exclusive" if u["exclusive"] else "shared"}
                    for u in users
                ]
            published = ((line_state or {}).get("gpus") or {}).get(uuid) or {}
            snapshot["claim"] = published.get("claim") if gpuline.queue_state_fresh(line_state, now) else None
            snapshot["warnings"] = []
            for item in snapshot.get("line") or []:
                if item.get("kind") != "job":
                    continue
                ready = item.get("ready_since_unix") or item.get("since_unix")
                if ready is not None and now - float(ready) > config.line.starvation_warn_s:
                    snapshot["warnings"].append(
                        f"queue job {item.get('job_id')} has waited {gpuline.duration_text(now - float(ready))} "
                        f"for the GPU (threshold {gpuline.duration_text(config.line.starvation_warn_s)})")
            entry["line"] = snapshot
        gpus.append(entry)
    return {
        "schema_version": "simqueue-node-status-v1",
        "node": config.node,
        "config": str(config.path),
        "paused": store.is_paused(),
        "memory_schedulable": config.memory_schedulable,
        "memory_allocated": mem_used,
        "memory_borrowed": borrowed_mem,
        "node_memory": _node_memory_status(config, store, jobs),
        "cpu_utilization": {"window_s": 0.5 if sample_s <= 0 else min(sample_s, 2.0),
                            "cpus": {str(k): round(v, 3) for k, v in sorted(utilization.items())},
                            "recorded_by_worker": _meta_json(store, "cpu_utilization")},
        "lock_files": config.locks.enabled,
        "reservation_depth": config.scheduler.reservation_depth,
        "sample_s": sample_s,
        "borrow": {
            "enabled": config.borrow.enabled,
            "slice": config.borrow.slice,
            "jobs": len(borrowers),
            "frozen": sum(1 for b in borrowers if b.frozen_since is not None),
        },
        "partitions": partitions,
        "gpus": gpus,
        "gpu_line": {
            "enabled": config.line.enabled,
            "state_path": config.line.state_path,
            "ticket_dir": config.line.ticket_dir,
            "starvation_warn_s": config.line.starvation_warn_s,
        },
    }


def _node_memory_status(config: NodeConfig, store: QueueStore, jobs: list[dict[str, Any]]) -> dict[str, Any]:
    """Requested against actual memory per running job, and the node's free
    memory (0.7.0), read from the job units' cgroups and /proc/meminfo."""
    rows = []
    for job in jobs:
        if job["state"] not in ACTIVE_STATES:
            continue
        usage = _live_usage(job, config)
        rows.append({
            "job_id": int(job["id"]),
            "name": job.get("name"),
            "borrow": bool(job.get("borrow")),
            "frozen": bool(job.get("frozen_at")),
            "request": job.get("mem_bytes"),
            "current": None if usage is None else usage.current,
            "peak": None if usage is None else usage.peak,
        })
    currents = [row["current"] or 0 for row in rows]
    node = read_node_memory(config.memory.test_node_total, currents)
    status: dict[str, Any] = {"jobs": rows, "queue_jobs_current": sum(currents),
                              "queue_jobs_requested": sum(int(row["request"] or 0) for row in rows)}
    if node is not None:
        high, low = water_marks(config.memory, node.total)
        status.update({"total": node.total, "available": node.available, "simulated": node.simulated,
                       "high_water": high, "low_water": low,
                       "outside_queue": max(0, node.used - sum(currents))})
    pressure = store.get_meta("memory_pressure")
    if pressure:
        try:
            status["pressure"] = json.loads(pressure)
        except ValueError:
            status["pressure"] = {"raw": pressure}
    # 0.8.0: the admission ceiling -- the worker's (best over the last hour)
    # and the one this moment's memory gives
    admission: dict[str, Any] = {}
    worker = _read_json_file(state_dir() / "admission.json")
    if isinstance(worker, dict) and worker.get("ceiling") is not None:
        admission["worker"] = worker
    if node is not None and config.memory.enabled:
        from .memory import admission_margin

        margin = admission_margin(config.memory, node.total)
        now = node.available + sum(currents) - margin
        if config.memory_schedulable is not None:
            now = min(now, config.memory_schedulable - margin)
        admission["now"] = now
    if admission:
        status["admission"] = admission
    return status


def _print_node_utilization(status: dict[str, Any]) -> None:
    """Node-level used against idle cores (0.8.0)."""
    cores = [cpu for part in status["partitions"] for cpu in part["cpus"]]
    measured = [cpu for cpu in cores if cpu.get("utilization") is not None]
    if not measured:
        return
    used = sum(cpu["utilization"] for cpu in measured)
    idle = [cpu["cpu"] for cpu in measured if cpu["utilization"] < 0.1]
    siblings = [s for cpu in cores for s in cpu.get("siblings") or [] if s.get("utilization") is not None]
    sibling_text = ""
    if siblings:
        sibling_text = (f"; SMT siblings {sum(s['utilization'] for s in siblings):.1f} of {len(siblings)} busy")
    print(f"CPU use over {status['cpu_utilization']['window_s']:g} s: cores {used:.1f} of {len(measured)} busy "
          f"({100 * used / len(measured):.0f}%), idle cores: {format_cpu_list(idle) or 'none'}{sibling_text}")


def _print_node_memory(status: dict[str, Any]) -> None:
    memory = status.get("node_memory") or {}
    if "available" in memory:
        print(f"Node memory: {format_gib(memory['available'])} available of {format_gib(memory['total'])}"
              + (" (simulated node)" if memory.get("simulated") else "")
              + f"; high-water {format_gib(memory['high_water'])}, low-water {format_gib(memory['low_water'])}; "
              f"queue jobs use {format_gib(memory['queue_jobs_current'])} of "
              f"{format_gib(memory['queue_jobs_requested'])} requested, outside the queue "
              f"{format_gib(memory['outside_queue'])}")
    admission = memory.get("admission") or {}
    worker = admission.get("worker")
    if worker:
        cap = worker.get("cap")
        print(f"Admission ceiling: ~{format_gib(worker['ceiling'])} (a request above it is parked, not reserved): "
              f"best MemAvailable with the queue idle over the last hour ~{format_gib(worker['best_idle_last_hour'])}"
              f" - margin {format_gib(worker['margin'])}"
              + (f", at most memory_schedulable - margin = {format_gib(cap)}" if cap is not None else "")
              + (f"; this moment's memory alone gives ~{format_gib(admission['now'])}" if "now" in admission else "")
              + f" [worker, {worker.get('at')}]")
    elif "now" in admission:
        print(f"Admission ceiling: ~{format_gib(admission['now'])} from this moment's memory "
              "(the worker has not published its value; it uses the best over the last hour)")
    pressure = memory.get("pressure")
    if pressure:
        culprit = pressure.get("culprit") or {}
        print(f"  MEMORY PRESSURE since {pressure.get('since')}: dispatch stopped, borrow jobs preempted"
              + (f"; job {culprit.get('job_id')} grew {format_gib(culprit.get('grew'))}" if culprit else ""))
    for row in memory.get("jobs") or []:
        request = row.get("request")
        current = row.get("current")
        percent = f" ({100.0 * current / request:.0f}%)" if current is not None and request else ""
        flags = (" borrow" if row.get("borrow") else "") + (" frozen" if row.get("frozen") else "")
        print(f"  job {row['job_id']:>5} {_abbreviate(str(row.get('name') or ''), 24):<24} "
              f"now {format_gib(current):>6}{percent:<7} peak {format_gib(row.get('peak')):>6} "
              f"request {format_size(request) if request is not None else '-':>5}{flags}")


def _meta_json(store: QueueStore, key: str) -> Any:
    raw = store.get_meta(key)
    try:
        return None if raw is None else json.loads(raw)
    except ValueError:
        return None


def _policy_text(part: dict[str, Any]) -> str:
    if part["timing_sensitive"]:
        return "timing-sensitive"
    if not part["lend"]:
        return "no-lend"
    return "lend" + (",leased-idle" if part["lend_leased_idle"] else "")


def _core_text(cpu: dict[str, Any]) -> str:
    cls = cpu.get("class")
    if cls == "job":
        text = f"job {cpu['job_id']}"
    elif cls == "borrowed":
        text = "free"
    elif cls in ("external-busy", "external-idle"):
        word = "busy" if cls == "external-busy" else "idle"
        text = (f"external lease {word} ({cpu['usage_cores']:.2f} cores over {cpu['usage_window_s']:g} s, "
                f"{cpu.get('tree_pids', 0)} process(es))")
    elif cls == "stale":
        text = f"stale lease: {cpu.get('holder')}"
    elif cls == "external":
        text = "external lease (not sampled)"
    else:
        text = "free"
    for holder in cpu.get("holders", [])[:4]:
        text += f"\n        holder pid {holder['pid']}: {holder['command']}"
    if cpu.get("utilization") is not None:
        text += f" [busy {100 * cpu['utilization']:.0f}%"
        siblings = [s for s in cpu.get("siblings") or [] if s.get("utilization") is not None]
        if siblings:
            text += "; sibling " + ", ".join(f"{s['cpu']} {100 * s['utilization']:.0f}%" for s in siblings)
        text += "]"
    for sibling in cpu.get("siblings") or []:
        for borrower in sibling.get("borrowers") or []:
            text += f"\n        sibling CPU {sibling['cpu']}: borrow job {borrower['job_id']}"
    for borrower in cpu.get("borrowers", []):
        if borrower["frozen"]:
            state = "frozen"
        elif borrower["stopping"]:
            state = f"being stopped ({borrower['stopping']})"
        else:
            state = "running"
        text += f"\n        borrow job {borrower['job_id']} ({borrower['class']}, {state})"
    if cpu.get("no_borrow"):
        text += "\n        borrow lease held exclusively: not lendable"
    return text


def _partitions(args: argparse.Namespace, store: QueueStore) -> int:
    config = load_node_config()
    if config is None:
        raise ValueError(f"no node configuration ({config_path()}); this queue runs in 0.5 mode")
    if args.sample < 0 or args.sample > 60:
        raise ValueError("--sample must be between 0 and 60 seconds")
    status = _node_status(config, store, sample_s=args.sample)
    if args.json:
        print(json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True))
        return 0
    print(f"Node {status['node']}  config {status['config']}" + ("  [queue paused]" if status["paused"] else ""))
    if status["memory_schedulable"] is not None:
        borrowed = f" (+{format_size(status['memory_borrowed'])} borrow jobs)" if status["memory_borrowed"] else ""
        print(
            f"Memory: {format_size(status['memory_allocated'])} allocated{borrowed} of "
            f"{format_size(status['memory_schedulable'])} schedulable"
        )
    _print_node_memory(status)
    for line in _outside_daemon_lines():
        print(line)
    _print_node_utilization(status)
    print(f"{'PARTITION':<13} {'CPUS':<8} {'FREE':<5} {'GPUS':<5} {'MEM(def/max)':<13} {'TIME(def/max)':<27} "
          f"{'QUEUED':<7} POLICY")
    for part in status["partitions"]:
        cpus = [c["cpu"] for c in part["cpus"]]
        free = sum(1 for c in part["cpus"] if c["class"] in ("free", "borrowed"))
        name = part["name"] + ("*" if part["default"] else "")
        times = f"{format_time_limit(part['default_time_s'])}/{format_time_limit(part['max_time_s'])}"
        mems = f"{format_size(part['default_mem'])}/{format_size(part['max_mem'])}"
        queued = str(part["queued"]) + (f" ({part['queued_borrow']}~)" if part["queued_borrow"] else "")
        print(
            f"{name:<13} {format_cpu_list(cpus):<8} {free:<5} {len(part['gpus']):<5} {mems:<13} "
            f"{times:<27} {queued:<7} {_policy_text(part)}"
        )
        for alias in part.get("aliases") or []:
            prefer = f", prefers CPUs {format_cpu_list(alias['prefer_cpus'])}" if alias["prefer_cpus"] else ""
            mem_text = (f"auto (else {format_size(alias['default_mem'])})" if alias.get("default_mem_auto")
                        else format_size(alias["default_mem"]))
            print(
                f"  alias {alias['name']}: defaults {alias['default_cpus']} core(s), "
                f"{mem_text}, {format_time_limit(alias['default_time_s'])}{prefer}"
            )
        for domain in part.get("l3_domains") or []:
            size = f", {domain['l3_bytes'] // 1024**2} MiB" if domain.get("l3_bytes") else ""
            jobs = "; ".join(f"job {job_id}: {cpus}" for job_id, cpus in domain["jobs"].items())
            print(f"  L3 domain {domain['index']}: cores {format_cpu_list(domain['cores'])} (CPUs "
                  f"{format_cpu_list(domain['logical_cpus'])}{size}): {domain['busy']}/{len(domain['cores'])} busy"
                  + (f", {domain['borrowed']} borrowed" if domain["borrowed"] else "")
                  + (f" ({jobs})" if jobs else "")
                  + ("" if part.get("l3_placement", True) else " [L3 placement off]"))
        for cpu in part["cpus"]:
            print(f"    CPU {cpu['cpu']}: {_core_text(cpu)}")
    for gpu in status["gpus"]:
        print(
            f"GPU {gpu['index'] if gpu['index'] is not None else '-'} {gpu['uuid']} {gpu['name']}: "
            f"claimed {format_size(gpu['claimed'])} of {format_size(gpu['budget'])} budget"
        )
        for user in gpu["jobs"]:
            mode = "exclusive" if user["exclusive"] else f"shared {format_size(user['claim'])}"
            print(f"    job {user['job_id']}: {mode}")
        line = gpu.get("line")
        if line is None:
            for holder in gpu["lease_held_outside_queue"]:
                print(f"    lease held outside the queue: {holder}")
            continue
        from . import gpuline

        claim = line.get("claim")
        for text in gpuline.render_snapshot(line, now=time.time()):
            print(text)
        if claim:
            held = claim.get("held_since_unix")
            print(f"    queue claim: {claim.get('unit')} {claim.get('phase')}"
                  + (f" (holding since {gpuline.now_iso(float(held))})" if held else "")
                  + (f", handing off to job {claim['handoff_job']}" if claim.get("handoff_job") else ""))
        for warning in line.get("warnings") or []:
            print(f"    WARNING: {warning}")
    return 0


def _node_doctor(config: NodeConfig, store: QueueStore) -> int:
    import shutil as _shutil

    from .units import JOBWRAP

    problems: list[str] = []
    load_state = _systemctl_property("LoadState")
    active_state = _systemctl_property("ActiveState")
    unit_file_state = _systemctl_property("UnitFileState")
    exec_start = _systemctl_property("ExecStart")
    runtime = str(Path(__file__).resolve().parents[2])
    if active_state != "active":
        problems.append("simqueue.service is not active")
    service_build = "unknown"
    if "path=" in exec_start:
        launcher = exec_start.split("path=", 1)[1].split(";", 1)[0].strip()
        try:
            text = Path(launcher).read_text(encoding="utf-8", errors="replace")
        except OSError:
            text = ""
        service_build = "this runtime" if runtime in text else f"another build ({launcher})"
        if runtime not in text:
            problems.append(f"simqueue.service does not run this runtime ({runtime})")
    wrapper_python = config.launch.wrapper_python
    if not os.access(wrapper_python, os.X_OK):
        problems.append(f"wrapper python is not executable: {wrapper_python}")
    if _shutil.which(config.launch.systemd_run) is None:
        problems.append(f"{config.launch.systemd_run} not found")
    manager = "unknown"
    try:
        result = subprocess.run(
            [config.launch.systemctl, "--user", "is-system-running"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False, timeout=10,
        )
        manager = result.stdout.strip() or "unknown"
    except (OSError, subprocess.TimeoutExpired):
        manager = "unreachable"
    if manager not in ("running", "degraded"):
        problems.append(f"user service manager state: {manager}")
    units = "unknown"
    try:
        result = subprocess.run(
            [config.launch.systemctl, "--user", "list-units", f"{config.launch.unit_prefix}-*",
             "--all", "--no-legend", "--plain"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False, timeout=10,
        )
        names = [line.split()[0] for line in result.stdout.splitlines() if line.strip()]
        units = ", ".join(names) if names else "none"
    except (OSError, subprocess.TimeoutExpired):
        pass
    lock_detail = "disabled"
    if config.locks.enabled:
        missing = [
            path
            for cpu in config.node_cpus()
            for path in config.locks.core_lock_paths(cpu)
            if not Path(path).parent.is_dir()
        ]
        lock_detail = "enabled" + (f"; missing directories for {missing}" if missing else "")
        if missing:
            problems.append("lease directory missing")
    user = pwd.getpwuid(os.getuid()).pw_name
    linger = "unknown"
    try:
        result = subprocess.run(
            ["loginctl", "show-user", user, "--property=Linger", "--value"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False,
        )
        if result.returncode == 0:
            linger = result.stdout.strip()
    except OSError:
        pass
    controllers_path = Path(
        f"/sys/fs/cgroup/user.slice/user-{os.getuid()}.slice/user@{os.getuid()}.service/cgroup.controllers"
    )
    try:
        controllers = controllers_path.read_text(encoding="ascii").split()
    except OSError:
        controllers = []
    freezer = controllers_path.with_name("cgroup.freeze").exists()
    borrow_detail = "disabled"
    if config.borrow.enabled:
        weight = _systemctl_property("CPUWeight", unit=config.borrow.slice)
        borrow_detail = (
            f"slice {config.borrow.slice} CPUWeight={weight or '?'} (wanted {config.borrow.slice_cpu_weight}); "
            f"units CPUWeight={config.borrow.unit_cpu_weight}, SCHED_IDLE, idle I/O; "
            f"cpu controller {'delegated' if 'cpu' in controllers else 'NOT delegated'}; "
            f"freezer {'available' if freezer else 'unavailable'}"
        )
        if "cpu" not in controllers:
            problems.append("the cpu controller is not delegated to the user manager (borrow jobs need cpu.weight)")
        if not freezer:
            problems.append("cgroup freezer unavailable (--on-preempt freeze cannot work)")
    line_detail = "disabled"
    if config.line.enabled:
        from . import gpuline

        state = gpuline.load_queue_state(config.line.state_path)
        if state is None:
            age_text = "not published yet"
        else:
            age = time.time() - float(state.get("heartbeat_unix") or 0.0)
            age_text = f"published {age:.0f} s ago" if state.get("heartbeat_unix") else "withdrawn (worker stopped)"
        ticket_dir = Path(config.line.ticket_dir)
        writable = os.access(ticket_dir if ticket_dir.exists() else ticket_dir.parent, os.W_OK)
        line_detail = (
            f"state {config.line.state_path} ({age_text}); tickets {config.line.ticket_dir}"
            f"{'' if writable else ' (NOT writable)'}; claim hold limit {config.line.claim_hold_max_s:g} s; "
            f"starvation warning after {gpuline.duration_text(config.line.starvation_warn_s)}; raw waiters "
            + ("requeued (" + ", ".join(config.line.requeue_commands) + ")" if config.line.requeue_raw_waiters
               else "not requeued")
        )
        if not writable:
            problems.append(f"GPU line ticket directory not writable: {config.line.ticket_dir}")
    lending = "; ".join(
        f"{p.name}=" + ("timing-sensitive" if p.timing_sensitive else
                        ("lend" + (",leased-idle" if p.lend_leased_idle else "") if p.lend else "no-lend"))
        for p in config.partitions.values()
    )
    fields = [
        ("SimQueue", f"{__version__} ({runtime})"),
        ("Mode", "node scheduler (partitions, concurrent jobs, one systemd unit per job)"),
        ("Node config", config.path),
        ("Project root", project_root()),
        ("State directory", state_dir()),
        ("Database", database_path()),
        ("Default partition", config.default_partition),
        ("Partitions", "; ".join(
            f"{p.name}={format_cpu_list(p.cpus)}" + (f"+{len(p.gpus)}GPU" if p.gpus else "")
            + (f" (aliases {', '.join(p.aliases)})" if p.aliases else "")
            for p in config.partitions.values()
        )),
        ("Lease files", lock_detail),
        ("Lending", lending),
        ("Borrow jobs", borrow_detail),
        ("GPU line", line_detail),
        ("Memory", _memory_doctor_text(config, store)),
        ("GitHub reporting", _github_doctor_text(config, store)),
        ("Containment", _containment_doctor_text(config)),
        ("Backfill", f"reservations for the first {config.scheduler.reservation_depth} blocked owner job(s) in "
                     "dispatch order (scheduler.reservation_depth; 1 = EASY, the head only)"),
        ("SMT jobs", _smt_doctor_text(config, store)),
        ("Worker wake FIFO", _wake_fifo_text()),
        ("Queue", "paused" if store.is_paused() else "accepting jobs"),
        ("Worker unit", "simqueue.service"),
        ("Service installed", load_state),
        ("Service active", active_state),
        ("Service enabled", unit_file_state),
        ("Service ExecStart", exec_start.split(" ; ignore_errors", 1)[0]),
        ("Service runs", service_build),
        ("User manager", manager),
        ("Job wrapper", f"{wrapper_python} {JOBWRAP}"),
        ("Job units", units),
        ("Survives logout (linger)", linger),
    ]
    width = max(len(label) for label, _ in fields)
    for label, value in fields:
        print(f"{label:<{width}}  {value}")
    for problem in problems:
        print(f"PROBLEM: {problem}")
    return 0 if not problems else 1


def _smt_doctor_text(config: NodeConfig, store: QueueStore) -> str:
    """0.8.2: what ``submit --smt`` gives on this node, and the SMT jobs now."""
    cores = config.node_cpus()
    siblings = config.smt_siblings_of(cores)
    added = sorted({cpu for found in siblings.values() for cpu in found})
    if not added:
        return "--smt unavailable: no core has a free SMT sibling on this node (SMT off?)"
    _, threads = config.smt_affinity(cores)
    lacking = [cpu for cpu, found in siblings.items() if not found]
    jobs = [job for job in store.node_jobs() if smt_job(job)]
    running = [f"{job['id']} on {job.get('cpu_list')}" for job in jobs if job["state"] in ACTIVE_STATES]
    queued = [str(job["id"]) for job in jobs if job["state"] == "queued"]
    now = ("running: " + ", ".join(running) if running else "none running") + (
        f"; queued: {', '.join(queued)}" if queued else "")
    return (f"submit --smt: CPU affinity = the job's cores and their SMT siblings (siblings {format_cpu_list(added)} "
            f"of cores {format_cpu_list(sorted(set(cores) - set(lacking)))}, {threads} threads per core"
            + (f"; cores {format_cpu_list(lacking)} have none" if lacking else "")
            + "); SIMQ_SMT=1, SIMQ_THREADS_PER_CORE, SIMQ_CPU_LIST_SMT; never lent (as --exclusive-timing); "
            + now)


def _outside_daemon_lines() -> list[str]:
    """0.7.2: MPI launcher daemons the worker saw outside every job's cgroup."""
    from .containment import read_daemons

    document = read_daemons(state_dir() / "containment.json")
    if not document:
        return []
    lines = []
    for daemon in document.get("outside_mpi_daemons") or []:
        children = daemon.get("children") or []
        lines.append(f"OUTSIDE MPI DAEMON: {daemon.get('comm')} pid {daemon.get('pid')} in "
                     f"{daemon.get('cgroup')} ({len(children)} child process(es)); jobs that attach to it "
                     f"(prun, --dvm) run their ranks outside the queue [scan {document.get('scanned_at')}]")
    return lines


def _containment_doctor_text(config: NodeConfig) -> str:
    containment = config.containment
    if not containment.enabled:
        return "off (jobs share /tmp; no escape scan)"
    scan = (f"escape scan every {containment.scan_seconds:g} s while jobs run and at each job's end"
            if containment.scan_seconds else "escape scan at each job's end")
    outside = _outside_daemon_lines()
    return (f"per-job TMPDIR/PMIx/PRRTE rendezvous under {containment.tmp_base}/simq.<job>.<token>; {scan}; "
            + (f"{len(outside)} outside MPI daemon(s) seen" if outside else "no outside MPI daemon seen"))


def _memory_doctor_text(config: NodeConfig, store: QueueStore) -> str:
    memory = config.memory
    if not memory.enabled:
        return "tracking off (memory.enabled false)"
    node = read_node_memory(memory.test_node_total)
    marks = ""
    if node is not None:
        high, low = water_marks(memory, node.total)
        marks = (f"; MemAvailable {format_gib(node.available)} of {format_gib(node.total)}; pressure below "
                 f"{format_gib(high)}, over at {format_gib(low)}")
    margin = "the high-water mark" if memory.admission_margin is None else format_gib(memory.admission_margin)
    pressure = store.get_meta("memory_pressure")
    return (f"cgroups sampled every {memory.sample_seconds:g} s, recorded every {memory.record_seconds:g} s; "
            f"admission guard margin {margin}{marks}; --mem auto {memory.auto_factor:g} x the largest of the "
            f"last {memory.auto_history} peaks" + ("; MEMORY PRESSURE now" if pressure else ""))


def _github_doctor_text(config: NodeConfig, store: QueueStore) -> str:
    github = config.github
    if not github.enabled:
        return "disabled (github.enabled false)"
    import shutil as _shutil

    if github.token_file:
        source = f"token file {github.token_file}" + ("" if os.path.exists(github.token_file) else " (MISSING)")
    else:
        found = _shutil.which(github.gh)
        source = f"`{github.gh} auth token`" + (f" ({found})" if found else " (gh NOT found)")
    due = store.github_pending()
    failing = [int(job["id"]) for job in due if job.get("github_error")]
    return (f"{github.api_base}; token from {source}; {len(due)} job(s) with posts due"
            + (f", failing: {failing[:10]}" if failing else ""))


def _wake_fifo_text() -> str:
    import stat as _stat

    path = state_dir() / WAKE_FIFO
    try:
        info = os.lstat(path)
    except OSError:
        return f"{path} missing (worker not started by 0.7.0)"
    return f"{path}" + ("" if _stat.S_ISFIFO(info.st_mode) else " (NOT a FIFO)")


def _add_ensemble_manifest_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("manifest", metavar="MANIFEST.json")
    parser.add_argument("-n", "--name", help="override the manifest ensemble name")
    parser.add_argument(
        "--max-parallel", type=int, help="maximum simultaneously active one-core members"
    )
    parser.add_argument(
        "--gpu-budget", metavar="SIZE", help="override gpu.budget (for example 8GiB)"
    )
    parser.add_argument(
        "--gpu-device", metavar="DEVICE", help="override gpu.device / CUDA_VISIBLE_DEVICES"
    )
    parser.add_argument(
        "--diagnostic-interval",
        type=float,
        metavar="SECONDS",
        help="GPU/process sample interval (default: 1.0)",
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        default=None,
        help="stop the cohort after the first member failure",
    )


def _add_node_request_options(parser: argparse.ArgumentParser, *, cpus: bool, gpus: bool) -> None:
    parser.add_argument(
        "-p", "--partition", metavar="PARTITION",
        help="node mode: partition or partition alias (default: the configured default partition)",
    )
    if cpus:
        parser.add_argument(
            "--cpus", metavar="LIST",
            help="node mode: explicit cores by partition CPU number, e.g. 6 or 6-7",
        )
    if gpus:
        parser.add_argument("--gpus", type=int, metavar="N", help="node mode: GPUs (default 0)")
        parser.add_argument(
            "--gpu-mem", metavar="SIZE",
            help="node mode: per-GPU memory claim; shares the GPU within its budget "
            "(default with --gpus: exclusive GPU)",
        )
    else:
        parser.set_defaults(gpus=None, gpu_mem=None, cpus=None)
    parser.add_argument(
        "--mem", metavar="SIZE",
        help="node mode: memory limit, e.g. 4G (bare number = MiB); 'auto': 1.25 x the largest recent peak "
        "of earlier jobs with this name (or name family), bounded by the partition's min_mem and max_mem",
    )
    parser.add_argument(
        "-t", "--time", metavar="TIME",
        help="node mode: time limit, e.g. 30, 2:00:00, 1-12:00:00; 'auto' (0.8): 1.5 x the 95th percentile runtime "
        "of earlier jobs with this name (or name family)",
    )
    parser.add_argument(
        "--l3-domains", type=int, metavar="N",
        help="node mode (0.8): spread the cores evenly over N L3 cache domains (default: one domain when it fits)",
    )
    parser.add_argument(
        "--split-domains", metavar="A,B",
        help="node mode (0.8): A cores in one L3 domain, B in another (sets the core count); the job finds them "
        "in SIMQ_L3_DOMAINS in that order",
    )
    parser.add_argument("--priority", type=int, metavar="N", help="node mode: priority (higher first; default 0)")


def _add_borrow_options(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("borrowing (0.6.1)")
    group.add_argument(
        "--borrow", "--opportunistic", dest="borrow", action="store_true",
        help="restartable work (builds, tests): may start on idle cores of any lending partition; "
        "runs at idle CPU priority and is preempted whenever an owner job needs the cores",
    )
    group.add_argument(
        "--borrow-from", metavar="P1,P2",
        help="borrow only from these partitions (besides the job's own); implies --borrow",
    )
    group.add_argument(
        "--max-borrow-cores", type=int, metavar="N",
        help="take up to N cores when that many are idle (at least -c, default 1); implies --borrow",
    )
    group.add_argument(
        "--on-preempt", choices=("freeze", "requeue", "cancel"),
        help="what preemption does to a borrow job (default requeue; freeze resumes it later)",
    )
    group.add_argument(
        "--exclusive-timing", action="store_true",
        help="owner job: its cores and their SMT siblings are never lent while it runs",
    )
    group.add_argument(
        "--smt", action="store_true",
        help="owner job (0.8.2): CPU affinity = its cores and their SMT siblings (e.g. -c 16: CPUs 0-15 and "
        "16-31), for 2 threads per core; exports SIMQ_SMT=1, SIMQ_THREADS_PER_CORE and SIMQ_CPU_LIST_SMT "
        "(SIMQ_CPU_LIST stays the cores); implies --exclusive-timing",
    )


def _add_submit_arguments(parser: argparse.ArgumentParser) -> None:
    """Options of ``simq submit`` (and ``simq run``, which takes the same)."""
    parser.add_argument("-n", "-J", "--name", "--job-name", dest="name", help="short job name")
    parser.add_argument("-C", "--cwd", help="working directory (default: current directory)")
    parser.add_argument(
        "-c",
        "--cores",
        "--cpus-per-task",
        dest="cores",
        type=int,
        help="physical cores (one logical CPU each; default: the default of the partition "
        "or alias named by -p)",
    )
    _add_node_request_options(parser, cpus=True, gpus=True)
    parser.add_argument(
        "-e",
        "--env",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="set a job environment variable (repeatable); a secret-looking name (*TOKEN*, *SECRET*, *PASSWORD*, "
        "*KEY*, GITHUB_*, ACTIONS_*, RUNNER_*) reaches the job but is stored redacted",
    )
    parser.add_argument(
        "--export",
        metavar="NAMES|ALL|NONE",
        help="pass these variables of your environment to the job besides HOME, USER, LOGNAME, LANG, PATH and TZ "
        "(ALL: every one); secret-looking names are never taken from the environment (0.7.0)",
    )
    parser.add_argument(
        "--shell",
        metavar="COMMAND",
        help="run COMMAND with /bin/bash -lc (explicit shell mode)",
    )
    parser.add_argument(
        "--idempotency-key",
        metavar="KEY",
        help="return the existing job when KEY is replayed with the same request",
    )
    parser.add_argument(
        "--json", action="store_true", help="print a structured submission result"
    )
    parser.add_argument(
        "--test-only",
        action="store_true",
        help="node mode: validate the request and show where it would run now; queue nothing",
    )
    _add_borrow_options(parser)
    group = parser.add_argument_group("GitHub reporting (0.7.0; posted by the worker, never by the job)")
    group.add_argument(
        "--github-status", action="append", metavar="OWNER/REPO@SHA:CONTEXT",
        help="commit status CONTEXT on SHA: pending when queued and running, success/failure from the exit "
        "code, error on cancel, time limit, out-of-memory kill or a lost unit (repeatable)",
    )
    group.add_argument(
        "--github-comment", metavar="OWNER/REPO#PR",
        help="one comment on the pull request at the job's end (state, runtime, peak memory, cores, last log "
        "lines, summary), edited in place on later posts; opt-in",
    )
    group.add_argument(
        "--github-summary", metavar="PATH",
        help="markdown file the job may write (exported as SIMQ_SUMMARY); appended to the comment",
    )
    parser.add_argument("command", nargs=argparse.REMAINDER, help="command after --")


def build_parser() -> argparse.ArgumentParser:
    invoked_as = Path(sys.argv[0]).name
    program = invoked_as if invoked_as in {"simq", "simqueue"} else "simq"
    parser = argparse.ArgumentParser(
        prog=program,
        description=(
            "Persistent simulation queue with serial outer jobs and optional "
            "multi-simulation ensembles."
        ),
        epilog=(
            f"shortcuts: {program} -active\n"
            f"           {program} -queue\n\n"
            "example: simq submit -n aging -C ~/projects/my-sim -- ./run.sh\n"
            "         simq submit -n setup-and-run --shell 'source env.sh && exec ./run.sh'"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"simq {__version__}")
    subparsers = parser.add_subparsers(dest="action")

    submit = subparsers.add_parser("submit", aliases=["add"], help="enqueue a command")
    _add_submit_arguments(submit)
    run = subparsers.add_parser(
        "run", aliases=["srun"],
        help="submit a command, stream its output and wait for it (srun-shaped)",
        description=(
            "Submit COMMAND with the options of `simq submit`, print the job ID, stream the job's "
            "combined log while it grows, wait until the job ends and exit with its exit code "
            "(0 for succeeded, 1 for other ends without one). SIGINT, SIGTERM or SIGHUP cancel the "
            "job and exit with 128+signal. With --idempotency-key, a replay of the same request "
            "reattaches to the existing job instead of submitting a second one. Messages go to "
            "stderr; stdout carries only the job log."
        ),
    )
    _add_submit_arguments(run)

    listing = subparsers.add_parser("list", aliases=["ls", "status"], help="list jobs")
    listing.add_argument(
        "--all", action="store_true", help="show the entire history (ignores --limit)"
    )
    listing.add_argument(
        "--limit", type=int, default=20, help="maximum rows in the default view (default: 20)"
    )

    subparsers.add_parser("active", help="show only active jobs")
    subparsers.add_parser("queue", help="show only jobs waiting in the queue")

    show = subparsers.add_parser("show", help="show one job")
    show.add_argument("job_id", type=int, metavar="ID")

    log = subparsers.add_parser("log", help="show a job's combined output")
    log.add_argument("job_id", type=int, metavar="ID")
    log.add_argument("-f", "--follow", action="store_true")
    log.add_argument("-n", "--lines", type=int, default=50)

    cancel = subparsers.add_parser("cancel", help="cancel a queued or running job")
    cancel.add_argument("job_id", type=int, metavar="ID")

    retry = subparsers.add_parser("retry", help="enqueue a fresh copy of a finished job")
    retry.add_argument("job_id", type=int, metavar="ID")

    promote = subparsers.add_parser(
        "promote", help="move a queued outer job to the dispatch front"
    )
    promote.add_argument("job_id", type=int, metavar="ID")
    promote.add_argument("--json", action="store_true")

    demote = subparsers.add_parser(
        "demote", help="move a queued outer job to the dispatch back"
    )
    demote.add_argument("job_id", type=int, metavar="ID")
    demote.add_argument("--json", action="store_true")

    wait = subparsers.add_parser("wait", help="wait without canceling on Ctrl-C")
    wait.add_argument("job_ids", type=int, nargs="*", metavar="ID")

    subparsers.add_parser("pause", help="stop launching new jobs")
    subparsers.add_parser("resume", help="resume launching jobs")
    subparsers.add_parser("doctor", help="show queue, CPU, and service configuration")

    ensemble = subparsers.add_parser(
        "ensemble", help="validate, submit, and inspect multi-simulation cohorts"
    )
    ensemble_commands = ensemble.add_subparsers(
        dest="ensemble_action", required=True
    )
    ensemble_submit = ensemble_commands.add_parser(
        "submit", help="enqueue a manifest as one outer job"
    )
    _add_ensemble_manifest_options(ensemble_submit)
    _add_node_request_options(ensemble_submit, cpus=False, gpus=False)
    ensemble_submit.add_argument(
        "--idempotency-key",
        metavar="KEY",
        help="return the existing ensemble when KEY is replayed with the same request",
    )
    ensemble_submit.add_argument(
        "--json", action="store_true", help="print a structured submission result"
    )
    ensemble_validate = ensemble_commands.add_parser(
        "validate", help="validate and canonicalize a manifest without submitting"
    )
    _add_ensemble_manifest_options(ensemble_validate)
    ensemble_validate.add_argument("-p", "--partition", metavar="PARTITION")
    ensemble_validate.add_argument(
        "--json", action="store_true", help="print the canonical manifest"
    )
    ensemble_status = ensemble_commands.add_parser(
        "status", help="show parent, lane, member, GPU, and UVM state"
    )
    ensemble_status.add_argument("job_id", type=int, metavar="ID")
    ensemble_status.add_argument("--json", action="store_true")
    ensemble_pause = ensemble_commands.add_parser(
        "pause", help="stop member backfill after admitted members finish"
    )
    ensemble_pause.add_argument("job_id", type=int, metavar="ID")
    ensemble_pause.add_argument("--json", action="store_true")
    ensemble_resume = ensemble_commands.add_parser(
        "resume", help="continue queued members under the same ensemble ID"
    )
    ensemble_resume.add_argument("job_id", type=int, metavar="ID")
    ensemble_resume.add_argument("--json", action="store_true")
    ensemble_log = ensemble_commands.add_parser(
        "log", help="show one ensemble member's output"
    )
    ensemble_log.add_argument("job_id", type=int, metavar="ID")
    ensemble_log.add_argument("member", metavar="MEMBER")
    ensemble_log.add_argument("-f", "--follow", action="store_true")
    ensemble_log.add_argument("-n", "--lines", type=int, default=50)

    worker = subparsers.add_parser(
        "worker", help="run the worker (normally managed by systemd)"
    )
    worker.add_argument(
        "--drain", action="store_true", help="process queued jobs and exit when empty"
    )
    worker.add_argument(
        "--wait-for-lock",
        action="store_true",
        help="wait behind an existing worker (used by the systemd service)",
    )
    worker.add_argument(
        "--once", action="store_true", help="node mode: one reconcile and dispatch pass, then exit"
    )
    partitions = subparsers.add_parser(
        "partitions", aliases=["sinfo"], help="show partitions, free cores, GPUs and lease holders"
    )
    partitions.add_argument("--json", action="store_true")
    partitions.add_argument(
        "--sample", type=float, default=1.0, metavar="SECONDS",
        help="CPU sampling window for lease holders (default 1.0; 0 skips sampling)",
    )
    priority = subparsers.add_parser("priority", help="set a queued job's priority (higher runs first)")
    priority.add_argument("job_id", type=int, metavar="ID")
    priority.add_argument("value", type=int, metavar="PRIORITY")
    internal_ensemble = subparsers.add_parser("_run-ensemble", help=argparse.SUPPRESS)
    internal_ensemble.add_argument("job_id", type=int, metavar="ID")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    raw_arguments = list(sys.argv[1:] if argv is None else argv)
    shortcut = {
        ("-active",): ["active"],
        ("--active",): ["active"],
        ("-queue",): ["queue"],
        ("--queue",): ["queue"],
    }
    raw_arguments = shortcut.get(tuple(raw_arguments), raw_arguments)
    parser = build_parser()
    args = parser.parse_args(raw_arguments)
    if args.action == "worker":
        return run_worker(drain=args.drain, wait_for_lock=args.wait_for_lock, once=args.once)
    try:
        # Validate the root-owned policy before QueueStore can create or mutate
        # state. This keeps XDG_STATE_HOME/SIMQUEUE_STATE_DIR overrides from
        # splitting retry, cancel, or resume operations into an orphaned DB.
        load_config()
        store = QueueStore()
        if args.action in MUTATING_ACTIONS:
            try:
                return _dispatch_action(args, store, parser)
            finally:
                _wake_worker()
        return _dispatch_action(args, store, parser)
    except (KeyError, RuntimeError, ValueError, ConfigError) as error:
        print(f"simq: {error}", file=sys.stderr)
        return 2


def _dispatch_action(args: argparse.Namespace, store: QueueStore, parser: argparse.ArgumentParser) -> int:
    if True:
        if args.action == "_run-ensemble":
            return run_ensemble(args.job_id)
        if args.action is None:
            _print_jobs(store, include_all=False, limit=20)
            return 0
        if args.action in ("submit", "add"):
            return _submit(args, store)
        if args.action in ("run", "srun"):
            return _run(args, store)
        if args.action in ("list", "ls", "status"):
            if args.limit < 1:
                raise ValueError("--limit must be positive")
            _print_jobs(store, include_all=args.all, limit=args.limit)
            return 0
        if args.action == "active":
            _print_active_jobs(store)
            return 0
        if args.action == "queue":
            _print_queued_jobs(store)
            return 0
        if args.action == "show":
            return _show(args.job_id, store)
        if args.action == "log":
            if args.lines < 0:
                raise ValueError("--lines cannot be negative")
            return _tail(args.job_id, store, lines=args.lines, follow=args.follow)
        if args.action == "cancel":
            state = store.cancel(args.job_id)
            if state == "cancelling":
                print(f"Cancellation requested for running job {args.job_id}.")
                job = store.require(args.job_id)
                node_config = load_node_config()
                if (
                    node_config is not None
                    and job.get("launch_backend") == "systemd-user"
                    and job.get("unit_name")
                ):
                    # The job unit is independent of the worker, so stop it
                    # directly; the worker records the final state.
                    stop_unit(node_config.launch, str(job["unit_name"]))
                try:
                    _nudge_worker(store)
                except WorkerStartError as error:
                    print(
                        f"simq: cancellation is recorded, but the protected runner "
                        f"did not start: {error}",
                        file=sys.stderr,
                    )
                    return 1
            elif state == "canceled":
                print(f"Canceled queued job {args.job_id}.")
            else:
                print(f"Job {args.job_id} is already {state}.")
            return 0
        if args.action == "retry":
            job = store.retry(args.job_id)
            print(f"Queued job {job['id']} as retry of job {args.job_id}.")
            try:
                _nudge_worker(store)
            except WorkerStartError as error:
                print(
                    f"simq: job {job['id']} remains queued, but the protected runner "
                    f"did not start: {error}",
                    file=sys.stderr,
                )
                return 1
            return 0
        if args.action == "promote":
            return _reorder_job(args, store, direction="front")
        if args.action == "demote":
            return _reorder_job(args, store, direction="back")
        if args.action == "wait":
            return _wait(args.job_ids, store)
        if args.action == "pause":
            store.set_paused(True)
            print("Queue paused; a running job will be allowed to finish.")
            return 0
        if args.action == "resume":
            store.set_paused(False)
            _nudge_worker(store)
            print("Queue resumed.")
            return 0
        if args.action == "doctor":
            return _doctor(store)
        if args.action in ("partitions", "sinfo"):
            return _partitions(args, store)
        if args.action == "priority":
            job = store.set_priority(args.job_id, args.value)
            print(f"Job {args.job_id} ({job['name']}) priority set to {job['priority']}.")
            _nudge_worker(store)
            return 0
        if args.action == "ensemble":
            if args.ensemble_action == "submit":
                return _ensemble_submit(args, store)
            if args.ensemble_action == "validate":
                return _ensemble_validate(args)
            if args.ensemble_action == "status":
                return _ensemble_status(args, store)
            if args.ensemble_action == "pause":
                return _ensemble_pause(args, store)
            if args.ensemble_action == "resume":
                return _ensemble_resume(args, store)
            if args.ensemble_action == "log":
                if args.lines < 0:
                    raise ValueError("--lines cannot be negative")
                return tail_member_log(
                    args.job_id,
                    args.member,
                    store,
                    lines=args.lines,
                    follow=args.follow,
                )
    parser.error(f"unhandled action: {args.action}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
