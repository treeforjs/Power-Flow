from __future__ import annotations

import argparse
import json
import os
import pwd
import shlex
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
    ConfigError,
    NodeConfig,
    config_path,
    format_size,
    format_time_limit,
    load_node_config,
    parse_size,
    parse_time_limit,
)
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


def _print_node_job_rows(jobs: Sequence[dict[str, Any]]) -> None:
    root = project_root()
    print(
        f"{'ID':>5}  {'STATE':<11} {'PARTITION':<11} {'CPU':<7} {'GPU':<4} "
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
        print(
            f"{job['id']:>5}  {job['state']:<11} {_abbreviate(str(job.get('partition') or '-'), 11):<11} "
            f"{_abbreviate(cpu, 7):<7} {gpu:<4} {_job_age(job):>8}  "
            f"{_abbreviate(str(job['name']), 24):<24} {_abbreviate(detail, 60)}"
        )


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


def _parse_environment(items: Sequence[str]) -> tuple[dict[str, str], list[str]]:
    environment = dict(os.environ)
    override_keys: list[str] = []
    for item in items:
        if "=" not in item:
            raise ValueError(f"environment override must be KEY=VALUE, got {item!r}")
        key, value = item.split("=", 1)
        if not key or "\x00" in key or "=" in key:
            raise ValueError(f"invalid environment variable name {key!r}")
        environment[key] = value
        override_keys.append(key)
    return environment, override_keys


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


def _nudge_worker(store: QueueStore) -> str:
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


NODE_OPTIONS = ("partition", "cpus", "gpus", "gpu_mem", "mem", "time", "priority")


def _node_request(
    args: argparse.Namespace, config: NodeConfig, *, cores: int | None = None
) -> tuple[int, dict[str, Any], dict[str, Any]]:
    """Resolve a submission's resource request against the node config.

    Returns (requested cores, stored scheduling fields, explicitly given
    fields).  Only the explicit fields enter the idempotency fingerprint, so
    a 0.5-style replay keeps its 0.5 fingerprint.
    """
    partition = config.partition(args.partition)
    explicit: dict[str, Any] = {}
    cpu_request = None
    cpu_list_arg = getattr(args, "cpus", None)
    if cpu_list_arg is not None:
        cpus = parse_cpu_list(cpu_list_arg)
        if not cpus:
            raise ValueError("--cpus needs at least one CPU")
        outside = sorted(set(cpus) - set(partition.cpus))
        if outside:
            raise ValueError(
                f"--cpus {cpu_list_arg}: CPU(s) {format_cpu_list(outside)} are not in partition "
                f"{partition.name} ({format_cpu_list(partition.cpus)}); name each physical core by "
                "its partition CPU"
            )
        cpu_request = format_cpu_list(cpus)
        explicit["cpu_request"] = cpu_request
    requested = cores
    if requested is None:
        requested = len(parse_cpu_list(cpu_request)) if cpu_request else partition.resolve_default_cpus()
    if cpu_request is not None and requested != len(parse_cpu_list(cpu_request)):
        raise ValueError("--cores and --cpus disagree on the number of cores")
    if requested < 1 or requested > len(partition.cpus):
        raise ValueError(
            f"--cores must be between 1 and {len(partition.cpus)} in partition {partition.name}"
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
        raise ValueError(f"partition {partition.name} has {len(partition.gpus)} GPU(s)")
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
    if getattr(args, "mem", None) is not None:
        mem = parse_size(args.mem, field_name="--mem")
        if mem <= 0 or mem > partition.max_mem:
            raise ValueError(
                f"--mem must be between 1M and {format_size(partition.max_mem)} in partition {partition.name}"
            )
        explicit["mem_bytes"] = mem
    time_limit = partition.default_time
    if getattr(args, "time", None) is not None:
        time_limit = parse_time_limit(args.time, field_name="--time")
        explicit["time_limit_s"] = time_limit
        if partition.max_time is not None and (time_limit is None or time_limit > partition.max_time):
            raise ValueError(
                f"--time exceeds partition {partition.name}'s limit {format_time_limit(partition.max_time)}"
            )
    priority = getattr(args, "priority", None)
    if priority is not None:
        explicit["priority"] = priority
    if args.partition is not None:
        explicit["partition"] = partition.name
    scheduling = {
        "partition": partition.name,
        "cpu_request": cpu_request,
        "gpus_requested": gpus,
        "gpu_mem_bytes": gpu_mem,
        "mem_bytes": mem,
        "swap_bytes": partition.swap_max,
        "time_limit_s": time_limit,
        "priority": priority or 0,
    }
    return requested, scheduling, explicit


def _reject_node_options_without_config(args: argparse.Namespace) -> None:
    given = [name for name in NODE_OPTIONS if getattr(args, name, None) is not None]
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
        parts.append(f"mem {format_size(int(job['mem_bytes']))}")
    parts.append(f"time {format_time_limit(job.get('time_limit_s'))}")
    if int(job.get("priority") or 0):
        parts.append(f"priority {job['priority']}")
    return ", ".join(parts)


def _test_only(config: NodeConfig, store: QueueStore, requested: int, scheduling: dict[str, Any]) -> int:
    """Validate a request and report where it would run now, like
    ``sbatch --test-only``.  Nothing is queued or launched."""
    from .locks import busy_resources, lock_plan_for
    from .nodeworker import job_request, running_entries
    from .scheduler import INFINITY, build_pool, plan

    now = time.time()
    jobs = store.node_jobs()
    running = running_entries(jobs, config, now)
    allocated_cpus = {cpu for item in running for cpu in item.allocation.cpus}
    allocated_gpus = {gpu for item in running for gpu in item.allocation.gpus}
    busy_cpus, busy_gpus = busy_resources(config, allocated_cpus=allocated_cpus, allocated_gpus=allocated_gpus)
    pool = build_pool(config, running, external_cpus=busy_cpus, external_gpus=busy_gpus)
    waiting = [j for j in jobs if j["state"] == "queued"]
    probe = dict(scheduling, id=0, requested_cores=requested)
    # A new submission goes behind every queued job of equal or higher priority.
    position = sum(1 for j in waiting if int(j.get("priority") or 0) >= int(scheduling["priority"] or 0))
    waiting.insert(position, probe)
    queue = [job_request(j, config) for j in waiting]
    decision = plan(config, pool, running, queue, now)
    started = {request.job_id: allocation for request, allocation in decision.starts}
    payload: dict[str, Any] = {
        "schema_version": "simqueue-test-only-v1",
        "partition": scheduling["partition"],
        "requested_cores": requested,
        "request": {k: scheduling[k] for k in scheduling},
        "queued_ahead": position,
    }
    if 0 in started:
        allocation = started[0]
        payload.update(
            {
                "would_start_now": True,
                "cpu_list": format_cpu_list(allocation.cpus),
                "gpus": list(allocation.gpus),
                "gpu_mode": "none" if not allocation.gpus else ("exclusive" if allocation.gpu_exclusive else "shared"),
                "lease_files": [
                    {"path": lock.path, "mode": lock.mode}
                    for lock in lock_plan_for(config, allocation.cpus, allocation.gpus,
                                              gpu_shared=not allocation.gpu_exclusive)
                ],
            }
        )
    else:
        payload["would_start_now"] = False
        payload["reason"] = decision.skipped.get(0, "waiting behind earlier jobs")
        if busy_cpus:
            payload["leases_held_outside_queue"] = sorted(busy_cpus)
    print(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


def _submit(args: argparse.Namespace, store: QueueStore) -> int:
    node_config = load_node_config()
    scheduling: dict[str, Any] | None = None
    explicit: dict[str, Any] | None = None
    if node_config is not None:
        requested, scheduling, explicit = _node_request(args, node_config, cores=args.cores)
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
    environment, override_keys = _parse_environment(args.env)
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
    )
    created = bool(job.get("submission_created", True))
    cpu_list = format_cpu_list(pcores[:requested])
    placement_text = f"on P-cores {cpu_list}"
    if node_config is not None:
        placement_text = f"in partition {job.get('partition')}: {_describe_request(job)}"
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
            args, node_config, cores=int(canonical["max_parallel"])
        )
        # The GPU request follows from the manifest, which is already part
        # of the ensemble fingerprint.
        explicit.pop("gpu_mem_bytes", None)
        explicit.pop("gpus_requested", None)
    job = store.submit_ensemble(
        name=str(canonical["name"]),
        cwd=project_root(),
        environment=dict(os.environ),
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
        fields.extend(
            [
                ("Partition", job.get("partition") or "-"),
                ("Priority", job.get("priority") or 0),
                ("Request", _describe_request(job)),
                ("GPUs", job.get("gpu_list") or "none"),
                ("GPU claim", format_size(job["gpu_claim_bytes"]) if job.get("gpu_list") else "-"),
                ("Memory limit", format_size(job.get("mem_bytes"))),
                ("Time limit", format_time_limit(job.get("time_limit_s"))),
                ("Unit", job.get("unit_name") or "-"),
                ("Result", job.get("result") or "-"),
                ("Adopted", job.get("adopted_at") or "-"),
                ("Job files", job.get("job_dir") or "-"),
            ]
        )
    if job["state"] == "queued":
        queued = store.queued_jobs()
        position = next(
            index
            for index, queued_job in enumerate(queued, start=1)
            if int(queued_job["id"]) == job_id
        )
        fields.append(("Queue position", f"{position} of {len(queued)}"))
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
    return 0


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


def _node_status(config: NodeConfig, store: QueueStore) -> dict[str, Any]:
    from .locks import busy_resources, process_description
    from .nodeworker import running_entries

    jobs = store.node_jobs()
    now = time.time()
    running = running_entries(jobs, config, now)
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
    busy_cpus, busy_gpus = busy_resources(
        config, allocated_cpus=set(owner), allocated_gpus=set(gpu_users)
    )
    partitions = []
    for part in config.partitions.values():
        cpus = []
        for cpu in part.cpus:
            if cpu in owner:
                state = f"job {owner[cpu]}"
            elif cpu in busy_cpus:
                state = "lease held: " + process_description(busy_cpus[cpu][0].pid)
            else:
                state = "free"
            cpus.append({"cpu": cpu, "state": state})
        queued = [j for j in jobs if j["state"] == "queued" and (j.get("partition") or config.default_partition) == part.name]
        partitions.append(
            {
                "name": part.name,
                "default": part.name == config.default_partition,
                "cpus": cpus,
                "gpus": list(part.gpus),
                "default_mem": part.default_mem,
                "max_mem": part.max_mem,
                "default_time_s": part.default_time,
                "max_time_s": part.max_time,
                "queued": len(queued),
            }
        )
    gpus = []
    for uuid, gpu in config.gpus.items():
        users = gpu_users.get(uuid, [])
        gpus.append(
            {
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
        )
    return {
        "schema_version": "simqueue-node-status-v1",
        "node": config.node,
        "config": str(config.path),
        "paused": store.is_paused(),
        "memory_schedulable": config.memory_schedulable,
        "memory_allocated": mem_used,
        "lock_files": config.locks.enabled,
        "partitions": partitions,
        "gpus": gpus,
    }


def _partitions(args: argparse.Namespace, store: QueueStore) -> int:
    config = load_node_config()
    if config is None:
        raise ValueError(f"no node configuration ({config_path()}); this queue runs in 0.5 mode")
    status = _node_status(config, store)
    if args.json:
        print(json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True))
        return 0
    print(f"Node {status['node']}  config {status['config']}" + ("  [queue paused]" if status["paused"] else ""))
    if status["memory_schedulable"] is not None:
        print(
            f"Memory: {format_size(status['memory_allocated'])} allocated of "
            f"{format_size(status['memory_schedulable'])} schedulable"
        )
    print(f"{'PARTITION':<13} {'CPUS':<8} {'FREE':<5} {'GPUS':<5} {'MEM(def/max)':<13} {'TIME(def/max)':<27} QUEUED")
    for part in status["partitions"]:
        cpus = [c["cpu"] for c in part["cpus"]]
        free = sum(1 for c in part["cpus"] if c["state"] == "free")
        name = part["name"] + ("*" if part["default"] else "")
        times = f"{format_time_limit(part['default_time_s'])}/{format_time_limit(part['max_time_s'])}"
        mems = f"{format_size(part['default_mem'])}/{format_size(part['max_mem'])}"
        print(
            f"{name:<13} {format_cpu_list(cpus):<8} {free:<5} {len(part['gpus']):<5} {mems:<13} "
            f"{times:<27} {part['queued']}"
        )
        for cpu in part["cpus"]:
            if cpu["state"] != "free":
                print(f"    CPU {cpu['cpu']}: {cpu['state']}")
    for gpu in status["gpus"]:
        print(
            f"GPU {gpu['index'] if gpu['index'] is not None else '-'} {gpu['uuid']} {gpu['name']}: "
            f"claimed {format_size(gpu['claimed'])} of {format_size(gpu['budget'])} budget"
        )
        for user in gpu["jobs"]:
            mode = "exclusive" if user["exclusive"] else f"shared {format_size(user['claim'])}"
            print(f"    job {user['job_id']}: {mode}")
        for holder in gpu["lease_held_outside_queue"]:
            print(f"    lease held outside the queue: {holder}")
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
            for p in config.partitions.values()
        )),
        ("Lease files", lock_detail),
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
        help="node mode: partition (default: the configured default partition)",
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
    parser.add_argument("--mem", metavar="SIZE", help="node mode: memory limit, e.g. 4G (bare number = MiB)")
    parser.add_argument("-t", "--time", metavar="TIME", help="node mode: time limit, e.g. 30, 2:00:00, 1-12:00:00")
    parser.add_argument("--priority", type=int, metavar="N", help="node mode: priority (higher first; default 0)")


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
    submit.add_argument("-n", "-J", "--name", "--job-name", dest="name", help="short job name")
    submit.add_argument("-C", "--cwd", help="working directory (default: current directory)")
    submit.add_argument(
        "-c",
        "--cores",
        "--cpus-per-task",
        dest="cores",
        type=int,
        help="physical cores (one logical CPU each; default: the partition default, "
        "normally all of its cores, as in 0.5)",
    )
    _add_node_request_options(submit, cpus=True, gpus=True)
    submit.add_argument(
        "-e",
        "--env",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="set a job environment variable (repeatable)",
    )
    submit.add_argument(
        "--shell",
        metavar="COMMAND",
        help="run COMMAND with /bin/bash -lc (explicit shell mode)",
    )
    submit.add_argument(
        "--idempotency-key",
        metavar="KEY",
        help="return the existing job when KEY is replayed with the same request",
    )
    submit.add_argument(
        "--json", action="store_true", help="print a structured submission result"
    )
    submit.add_argument(
        "--test-only",
        action="store_true",
        help="node mode: validate the request and show where it would run now; queue nothing",
    )
    submit.add_argument("command", nargs=argparse.REMAINDER, help="command after --")

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
        if args.action == "_run-ensemble":
            return run_ensemble(args.job_id)
        if args.action is None:
            _print_jobs(store, include_all=False, limit=20)
            return 0
        if args.action in ("submit", "add"):
            return _submit(args, store)
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
    except (KeyError, RuntimeError, ValueError, ConfigError) as error:
        print(f"simq: {error}", file=sys.stderr)
        return 2
    parser.error(f"unhandled action: {args.action}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
