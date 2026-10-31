"""Resource accounting and EASY-backfill planning for the node scheduler.

Everything here is pure: callers pass the current time, the free resources
and the queue, and receive the jobs that may start now.  This keeps the
dispatch decisions testable with a fake clock.

Policy
------
Waiting jobs are considered in priority order (higher ``priority`` first,
then queue order).  A job that fits the free resources starts, unless an
earlier job is already blocked: the first blocked job receives a reservation
at its *shadow time* (the earliest time it can start, assuming every running
job ends at its time limit).  A later job may then start only if starting it
does not move that shadow time later (EASY backfill).  Jobs without a time
limit are treated as never ending, so they can backfill only into resources
the reserved job will not need.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable

from .nodeconfig import NodeConfig, Partition


INFINITY = math.inf


@dataclass(frozen=True)
class Request:
    job_id: int
    partition: str
    ncpus: int
    cpu_list: tuple[int, ...] | None = None  # explicit cores
    gpus: int = 0
    gpu_mem: int | None = None  # per-GPU claim; None = exclusive GPU
    mem: int = 0
    time_limit: int | None = None
    blocked: bool = False  # a paused ensemble holds its queue position


@dataclass(frozen=True)
class Allocation:
    cpus: tuple[int, ...]
    gpus: tuple[str, ...] = ()
    gpu_claim: int = 0  # bytes per GPU
    gpu_exclusive: bool = False
    mem: int = 0


@dataclass
class Running:
    job_id: int | None
    allocation: Allocation
    end_time: float  # INFINITY when unknown or unlimited


@dataclass
class Pool:
    """Free node resources."""

    free_cpus: set[int]
    gpu_budget: dict[str, int]
    gpu_used: dict[str, int] = field(default_factory=dict)
    gpu_exclusive: set[str] = field(default_factory=set)
    gpu_blocked: set[str] = field(default_factory=set)  # held outside the queue
    mem_free: int | None = None

    def copy(self) -> "Pool":
        return Pool(
            free_cpus=set(self.free_cpus),
            gpu_budget=dict(self.gpu_budget),
            gpu_used=dict(self.gpu_used),
            gpu_exclusive=set(self.gpu_exclusive),
            gpu_blocked=set(self.gpu_blocked),
            mem_free=self.mem_free,
        )

    def _gpu_available(self, uuid: str, claim: int | None) -> bool:
        if uuid in self.gpu_blocked or uuid in self.gpu_exclusive or uuid not in self.gpu_budget:
            return False
        used = self.gpu_used.get(uuid, 0)
        if claim is None:
            return used == 0
        return used + claim <= self.gpu_budget[uuid]

    def try_allocate(self, request: Request, partition: Partition) -> Allocation | None:
        if request.cpu_list is not None:
            if not set(request.cpu_list) <= self.free_cpus:
                return None
            cpus = tuple(sorted(request.cpu_list))
        else:
            candidates = sorted(cpu for cpu in partition.cpus if cpu in self.free_cpus)
            if len(candidates) < request.ncpus:
                return None
            cpus = tuple(candidates[: request.ncpus])
        if self.mem_free is not None and request.mem > self.mem_free:
            return None
        chosen: list[str] = []
        if request.gpus:
            for uuid in partition.gpus:
                if self._gpu_available(uuid, request.gpu_mem):
                    chosen.append(uuid)
                    if len(chosen) == request.gpus:
                        break
            if len(chosen) < request.gpus:
                return None
        exclusive = bool(request.gpus) and request.gpu_mem is None
        claim = 0
        if request.gpus:
            claim = (
                max(self.gpu_budget[uuid] for uuid in chosen) if exclusive else int(request.gpu_mem or 0)
            )
        return Allocation(cpus=cpus, gpus=tuple(chosen), gpu_claim=claim, gpu_exclusive=exclusive, mem=request.mem)

    def take(self, allocation: Allocation) -> None:
        self.free_cpus -= set(allocation.cpus)
        for uuid in allocation.gpus:
            if allocation.gpu_exclusive:
                self.gpu_exclusive.add(uuid)
            self.gpu_used[uuid] = self.gpu_used.get(uuid, 0) + allocation.gpu_claim
        if self.mem_free is not None:
            self.mem_free -= allocation.mem

    def give(self, allocation: Allocation) -> None:
        self.free_cpus |= set(allocation.cpus)
        for uuid in allocation.gpus:
            if allocation.gpu_exclusive:
                self.gpu_exclusive.discard(uuid)
            self.gpu_used[uuid] = max(0, self.gpu_used.get(uuid, 0) - allocation.gpu_claim)
        if self.mem_free is not None:
            self.mem_free += allocation.mem


def build_pool(
    config: NodeConfig,
    running: Iterable[Running],
    *,
    external_cpus: Iterable[int] = (),
    external_gpus: Iterable[str] = (),
) -> Pool:
    """Free resources after subtracting running allocations and outside holders."""
    pool = Pool(
        free_cpus=set(config.node_cpus()),
        gpu_budget={uuid: gpu.budget_bytes for uuid, gpu in config.gpus.items()},
        mem_free=config.memory_schedulable,
    )
    pool.free_cpus -= set(external_cpus)
    pool.gpu_blocked |= set(external_gpus)
    for item in running:
        pool.take(item.allocation)
    return pool


def request_is_feasible(request: Request, partition: Partition, config: NodeConfig) -> str | None:
    """Return why *request* can never run on an empty node, else ``None``."""
    if request.cpu_list is not None:
        if not set(request.cpu_list) <= set(partition.cpus):
            return "requested CPUs are outside the partition"
    elif request.ncpus < 1 or request.ncpus > len(partition.cpus):
        return f"partition {partition.name} has {len(partition.cpus)} cores"
    if request.gpus > len(partition.gpus):
        return f"partition {partition.name} has {len(partition.gpus)} GPU(s)"
    if request.gpus and request.gpu_mem is not None:
        if any(request.gpu_mem > config.gpus[uuid].budget_bytes for uuid in partition.gpus[: request.gpus]):
            return "GPU memory claim exceeds the GPU budget"
    if config.memory_schedulable is not None and request.mem > config.memory_schedulable:
        return "memory request exceeds the node's schedulable memory"
    return None


def earliest_start(
    request: Request,
    partition: Partition,
    pool: Pool,
    running: list[Running],
    now: float,
) -> float:
    """Shadow time: earliest start assuming running jobs end at their limits."""
    trial = pool.copy()
    if trial.try_allocate(request, partition) is not None:
        return now
    for item in sorted(running, key=lambda r: r.end_time):
        if item.end_time == INFINITY:
            break
        trial.give(item.allocation)
        if trial.try_allocate(request, partition) is not None:
            return max(now, item.end_time)
    return INFINITY


@dataclass
class Decision:
    starts: list[tuple[Request, Allocation]]
    reservation: tuple[int, float] | None  # (job id, shadow time)
    skipped: dict[int, str]  # job id -> reason (diagnostic)


def plan(
    config: NodeConfig,
    pool: Pool,
    running: list[Running],
    queue: list[Request],
    now: float,
    *,
    backfill: bool = True,
) -> Decision:
    """Choose the queued jobs that start now.

    *queue* must already be in dispatch order.  *pool* and *running* are not
    modified.
    """
    pool = pool.copy()
    running = list(running)
    starts: list[tuple[Request, Allocation]] = []
    skipped: dict[int, str] = {}
    reservation: tuple[Request, Partition, float] | None = None
    held_partitions: set[str] = set()
    for request in queue:
        partition = config.partitions.get(request.partition)
        if partition is None:
            skipped[request.job_id] = f"unknown partition {request.partition}"
            continue
        if request.partition in held_partitions:
            skipped[request.job_id] = "held behind a paused ensemble"
            continue
        if request.blocked:
            held_partitions.add(request.partition)
            skipped[request.job_id] = "paused ensemble"
            continue
        reason = request_is_feasible(request, partition, config)
        if reason is not None:
            skipped[request.job_id] = reason
            continue
        allocation = pool.try_allocate(request, partition)
        if allocation is not None:
            end = INFINITY if request.time_limit is None else now + request.time_limit
            if reservation is None:
                pool.take(allocation)
                running.append(Running(request.job_id, allocation, end))
                starts.append((request, allocation))
                continue
            if not backfill:
                skipped[request.job_id] = "waiting behind a higher-priority job"
                continue
            reserved_request, reserved_partition, shadow = reservation
            trial = pool.copy()
            trial.take(allocation)
            trial_running = [*running, Running(request.job_id, allocation, end)]
            new_shadow = earliest_start(reserved_request, reserved_partition, trial, trial_running, now)
            if new_shadow <= shadow:
                pool.take(allocation)
                running.append(Running(request.job_id, allocation, end))
                starts.append((request, allocation))
            else:
                skipped[request.job_id] = (
                    f"backfill would delay job {reserved_request.job_id}"
                )
            continue
        skipped[request.job_id] = "resources busy"
        if reservation is None:
            shadow = earliest_start(request, partition, pool, running, now)
            reservation = (request, partition, shadow)
    return Decision(
        starts=starts,
        reservation=None if reservation is None else (reservation[0].job_id, reservation[2]),
        skipped=skipped,
    )
