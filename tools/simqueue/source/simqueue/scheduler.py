"""Resource accounting and backfill planning for the node scheduler.

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

0.8.1: the first ``scheduler.reservation_depth`` blocked jobs (default 8)
hold reservations, not only the first.  Reservations are placed one after
another in a resource-time profile: each later reservation gets the
earliest time at which it fits beside the running jobs (until their limits)
and beside the reservations ahead of it (each holding, from its reserved
time until its time limit, the cores, memory and GPU it would take then).
A job that fits now starts only if every reserved job ahead of it can still
start at its reserved time once the job runs (its own and the earlier
reservations' resources re-placed around it).  Depth 1 is exactly the EASY
backfill of 0.8.0.  Without the deeper reservations, narrow low-priority
jobs could take every core that frees up one by one while a wide job behind
a far-off head reservation never saw enough free cores at once (zenobia,
2026-11-06).  Jobs refused because they would delay a reservation, parked
jobs (memory) and held or infeasible jobs hold no reservation.

Borrowing (0.6.1) is layered on top in ``borrow.py``: jobs submitted with
``--borrow`` are planned after every owner job, and cores they occupy are
"soft": free for owners (who preempt the borrower) but chosen only after
truly free cores.  With no borrow jobs the functions here behave exactly as
in 0.6.0.
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
    partition: str  # as submitted: a partition or an alias of one (0.6.3)
    ncpus: int
    cpu_list: tuple[int, ...] | None = None  # explicit cores
    gpus: int = 0
    gpu_mem: int | None = None  # per-GPU claim; None = exclusive GPU
    mem: int = 0
    time_limit: int | None = None
    blocked: bool = False  # a paused ensemble holds its queue position
    # 0.6.1 borrowing (ignored by plan(); see borrow.py)
    borrow: bool = False
    max_cpus: int | None = None  # borrow: take up to this many cores (>= ncpus)
    borrow_from: tuple[str, ...] | None = None  # borrow: lenders; None = every lending partition
    exclusive_timing: bool = False  # owner: its cores and SMT siblings are never lent
    on_preempt: str = "requeue"
    # 0.8.0 L3 domains: spread the cores evenly over N domains, or give each
    # listed part (core counts, in order) a domain of its own
    l3_domains: int | None = None
    split_domains: tuple[int, ...] | None = None


@dataclass(frozen=True)
class Allocation:
    cpus: tuple[int, ...]
    gpus: tuple[str, ...] = ()
    gpu_claim: int = 0  # bytes per GPU
    gpu_exclusive: bool = False
    mem: int = 0
    borrow: bool = False
    # borrow: (cpu, partition, class) per core; class is free | leased-idle | owner-idle
    placement: tuple[tuple[int, str, str], ...] = ()


@dataclass
class Running:
    job_id: int | None
    allocation: Allocation
    end_time: float  # INFINITY when unknown or unlimited
    # 0.7.1: what the job uses now; its end gives back this much actual memory
    # (None: its whole request, as for a job the plan itself starts)
    mem_actual: int | None = None


@dataclass
class Pool:
    """Free node resources."""

    free_cpus: set[int]
    gpu_budget: dict[str, int]
    gpu_used: dict[str, int] = field(default_factory=dict)
    gpu_exclusive: set[str] = field(default_factory=set)
    gpu_blocked: set[str] = field(default_factory=set)  # held outside the queue
    mem_free: int | None = None
    # 0.6.1: free for owners but occupied by a borrower; chosen last
    soft_cpus: set[int] = field(default_factory=set)
    # 0.7.1: actual memory a start may take (MemAvailable - margin - what
    # young jobs have not allocated); None: not accounted
    actual_free: int | None = None
    # 0.8.0: cores grouped by L3 domain (one group: no domain placement)
    domains: tuple[tuple[int, ...], ...] = ()
    heavy: int = 2
    domain_nodes: tuple[int, ...] = ()  # NUMA node of each domain group

    def copy(self) -> "Pool":
        return Pool(
            free_cpus=set(self.free_cpus),
            gpu_budget=dict(self.gpu_budget),
            gpu_used=dict(self.gpu_used),
            gpu_exclusive=set(self.gpu_exclusive),
            gpu_blocked=set(self.gpu_blocked),
            mem_free=self.mem_free,
            soft_cpus=set(self.soft_cpus),
            actual_free=self.actual_free,
            domains=self.domains,
            heavy=self.heavy,
            domain_nodes=self.domain_nodes,
        )

    def _choose_cpus(self, request: Request, partition: Partition) -> tuple[int, ...] | None:
        """Free cores of the partition for *request* (0.8.0 placement by L3
        domain and NUMA node).

        Within a domain: truly free cores before cores a borrower occupies,
        then the alias's preferred cores, then the lowest.  A job goes into
        as few NUMA nodes and, within them, as few L3 domains as possible: one
        domain when one has room -- a heavy job (at least ``heavy`` cores) the
        domain with the most free cores (ties: the node with the most free
        cores), so concurrent heavy jobs spread over caches and nodes; a light
        one the fullest domain that still fits, so whole domains stay free for
        heavy jobs.  Otherwise one NUMA node filled domain by domain (most
        free first), else nodes in turn.  ``l3_domains`` / ``split_domains``
        ask for several domains, kept in as few nodes as possible."""
        rank = {cpu: index for index, cpu in enumerate(partition.prefer_cpus)}

        def order(cpus) -> list[int]:
            return sorted(cpus, key=lambda cpu: (cpu in self.soft_cpus, rank.get(cpu, len(rank)), cpu))

        free = [cpu for cpu in partition.cpus if cpu in self.free_cpus]
        indexed = [(i, group) for i, group in enumerate(self.domains or (tuple(partition.cpus),))
                   if set(group) & set(partition.cpus)]
        groups = [group for _, group in indexed]
        nodes_of = [self.domain_nodes[i] if i < len(self.domain_nodes) else 0 for i, _ in indexed]
        if len(groups) <= 1 and not request.l3_domains and not request.split_domains:
            chosen = order(free)
            return tuple(sorted(chosen[: request.ncpus])) if len(chosen) >= request.ncpus else None
        by_domain = [order(cpu for cpu in free if cpu in group) for group in groups]
        node_free: dict[int, int] = {}
        for index, node in enumerate(nodes_of):
            node_free[node] = node_free.get(node, 0) + len(by_domain[index])

        def preference(index: int) -> tuple:
            cpus = by_domain[index]
            hard = sum(1 for cpu in cpus if cpu not in self.soft_cpus)
            preferred = sum(1 for cpu in cpus if cpu in rank)
            return hard, preferred, node_free[nodes_of[index]], len(cpus)

        parts: list[int] | None = None
        if request.split_domains:
            parts = list(request.split_domains)
        elif request.l3_domains:
            count = request.l3_domains
            base, extra = divmod(request.ncpus, count)
            parts = [base + (1 if i < extra else 0) for i in range(count)]
        if parts is not None:
            if len(parts) > len(groups) or any(part < 1 for part in parts):
                return None
            # the node (or fewest nodes) that can hold every part, most free first
            best_plan: list[int] | None = None
            for home in sorted(node_free, key=lambda n: (-node_free[n], n)):
                used: list[int] = []
                for part in sorted(parts, reverse=True):
                    fitting = [i for i in range(len(groups)) if i not in used and len(by_domain[i]) >= part]
                    if not fitting:
                        used = []
                        break
                    used.append(max(fitting, key=lambda i: (nodes_of[i] == home,
                                                            any(nodes_of[i] == nodes_of[u] for u in used),
                                                            preference(i), -i)))
                if used and (best_plan is None or len({nodes_of[i] for i in used})
                             < len({nodes_of[i] for i in best_plan})):
                    best_plan = used
            if best_plan is None:
                return None
            chosen = []
            for part, index in zip(sorted(parts, reverse=True), best_plan):
                chosen.extend(by_domain[index][:part])
            return tuple(sorted(chosen))
        fitting = [i for i in range(len(groups)) if len(by_domain[i]) >= request.ncpus]
        if fitting:
            if request.ncpus >= self.heavy:
                best = max(fitting, key=lambda i: (preference(i), -i))
            else:
                best = min(fitting, key=lambda i: (len(by_domain[i]) - request.ncpus,
                                                  -preference(i)[1], i))
            return tuple(sorted(by_domain[best][: request.ncpus]))
        if len(free) < request.ncpus:
            return None

        def fill(indices: list[int]) -> list[int]:
            chosen: list[int] = []
            for index in sorted(indices, key=lambda i: (-len(by_domain[i]), i)):
                chosen.extend(by_domain[index][: request.ncpus - len(chosen)])
                if len(chosen) == request.ncpus:
                    break
            return chosen

        homes = [n for n in node_free if node_free[n] >= request.ncpus]
        if homes:
            # one NUMA node: the one needing the fewest domains, then the most free
            def domains_needed(node: int) -> int:
                sizes = sorted((len(by_domain[i]) for i in range(len(groups)) if nodes_of[i] == node), reverse=True)
                total = 0
                for count, size in enumerate(sizes, 1):
                    total += size
                    if total >= request.ncpus:
                        return count
                return len(sizes) + 1

            home = min(homes, key=lambda n: (domains_needed(n), -node_free[n], n))
            return tuple(sorted(fill([i for i in range(len(groups)) if nodes_of[i] == home])))
        chosen = []
        for node in sorted(node_free, key=lambda n: (-node_free[n], n)):
            remaining = request.ncpus - len(chosen)
            part = fill([i for i in range(len(groups)) if nodes_of[i] == node])[:remaining]
            chosen.extend(part)
            if len(chosen) == request.ncpus:
                break
        return tuple(sorted(chosen))

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
            # Truly free cores before cores a borrower occupies; then the
            # cores the requested alias prefers (0.6.3); then the lowest --
            # inside one L3 domain where possible (0.8.0, _choose_cpus).
            found = self._choose_cpus(request, partition)
            if found is None:
                return None
            cpus = found
        if self.mem_free is not None and request.mem > self.mem_free:
            return None
        if self.actual_free is not None and request.mem > self.actual_free:
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

    def take(self, allocation: Allocation, *, actual: bool = True) -> None:
        self.free_cpus -= set(allocation.cpus)
        for uuid in allocation.gpus:
            if allocation.gpu_exclusive:
                self.gpu_exclusive.add(uuid)
            self.gpu_used[uuid] = self.gpu_used.get(uuid, 0) + allocation.gpu_claim
        if self.mem_free is not None:
            self.mem_free -= allocation.mem
        if self.actual_free is not None and actual:
            self.actual_free -= allocation.mem

    def give(self, allocation: Allocation, mem_actual: int | None = None) -> None:
        self.free_cpus |= set(allocation.cpus)
        for uuid in allocation.gpus:
            if allocation.gpu_exclusive:
                self.gpu_exclusive.discard(uuid)
            self.gpu_used[uuid] = max(0, self.gpu_used.get(uuid, 0) - allocation.gpu_claim)
        if self.mem_free is not None:
            self.mem_free += allocation.mem
        if self.actual_free is not None:
            self.actual_free += allocation.mem if mem_actual is None else mem_actual


def build_pool(
    config: NodeConfig,
    running: Iterable[Running],
    *,
    external_cpus: Iterable[int] = (),
    external_gpus: Iterable[str] = (),
    actual_free: int | None = None,
) -> Pool:
    """Free resources after subtracting running allocations and outside holders.

    *actual_free* (0.7.1): the actual memory starts may take now; running
    jobs already use theirs, so it is not reduced by their requests."""
    topology = config.topology
    pool = Pool(
        free_cpus=set(config.node_cpus()),
        gpu_budget={uuid: gpu.budget_bytes for uuid, gpu in config.gpus.items()},
        mem_free=config.memory_schedulable,
        domains=topology.domains if topology.l3_placement else (),
        heavy=topology.heavy(),
        domain_nodes=topology.domain_nodes() if topology.l3_placement else (),
    )
    pool.free_cpus -= set(external_cpus)
    pool.gpu_blocked |= set(external_gpus)
    for item in running:
        pool.take(item.allocation, actual=False)
    pool.actual_free = actual_free
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
    if request.l3_domains or request.split_domains:
        groups = [g for g in config.topology.domains if set(g) & set(partition.cpus)] or [partition.cpus]
        sizes = sorted((len(set(g) & set(partition.cpus)) for g in groups), reverse=True)
        parts = list(request.split_domains or ())
        if request.l3_domains:
            base, extra = divmod(request.ncpus, request.l3_domains)
            parts = [base + (1 if i < extra else 0) for i in range(request.l3_domains)]
        parts.sort(reverse=True)
        if len(parts) > len(sizes) or any(part > size for part, size in zip(parts, sizes)) or min(parts) < 1:
            return (f"partition {partition.name} has {len(groups)} L3 domain(s) of "
                    f"{', '.join(str(s) for s in sizes)} core(s): cannot place {'+'.join(map(str, parts))}")
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
        trial.give(item.allocation, item.mem_actual)
        if trial.try_allocate(request, partition) is not None:
            return max(now, item.end_time)
    return INFINITY


@dataclass
class Decision:
    starts: list[tuple[Request, Allocation]]
    reservation: tuple[int, float] | None  # (job id, shadow time) of the first reservation
    skipped: dict[int, str]  # job id -> reason (diagnostic)
    # 0.8.1: every reservation in dispatch order, the first included
    reservations: tuple[tuple[int, float], ...] = ()
    # 0.8.1: job refused as backfill -> the reserved job it would delay
    delays: dict[int, int] = field(default_factory=dict)


@dataclass(frozen=True)
class Reservation:
    """A blocked job's reservation (0.8.1): the time it is promised to start
    by (*shadow*) and what it would take then (*allocation*; ``None`` while
    the shadow time is unknown), which later reservations leave alone."""

    request: Request
    partition: Partition
    shadow: float
    allocation: Allocation | None

    @property
    def end(self) -> float:
        if self.shadow == INFINITY or self.request.time_limit is None:
            return INFINITY
        return self.shadow + self.request.time_limit


def _shape(request: Request) -> tuple:
    """What decides whether a request fits and what it would delay: two
    queued jobs of the same shape get the same answer from the same pool."""
    return (request.partition, request.ncpus, request.cpu_list, request.gpus, request.gpu_mem, request.mem,
            request.time_limit, request.l3_domains, request.split_domains)


def _window(free: Pool, later: list[Running], start: float, end: float, ahead: list[Reservation]) -> Pool:
    """What a job may hold from *start* until *end*: *free* (the free
    resources at *start*, every running job that ended by then given back)
    less what the reservations *ahead* hold at any moment of that window.
    *later* lists the running jobs that end after *start*; the memory and GPU
    budget they give back during the window count from their end on."""
    overlapping = [r for r in ahead if r.allocation is not None and r.shadow < end and r.end > start]
    if not overlapping:
        return free
    window = free.copy()
    for reservation in overlapping:
        allocation = reservation.allocation
        assert allocation is not None
        window.free_cpus -= set(allocation.cpus)
        if allocation.gpu_exclusive:
            window.gpu_exclusive |= set(allocation.gpus)
    # Counted resources: the least that is left at any moment of the window;
    # it changes only where a reservation starts.
    points = [start, *sorted({r.shadow for r in overlapping if r.shadow > start})]
    mem_low: int | None = None
    actual_low: int | None = None
    gpu_high: dict[str, int] = {}
    for point in points:
        ended = [item for item in later if item.end_time <= point]
        active = [r.allocation for r in overlapping if r.shadow <= point < r.end]
        held = sum(a.mem for a in active if a is not None)
        if free.mem_free is not None:
            value = free.mem_free + sum(item.allocation.mem for item in ended) - held
            mem_low = value if mem_low is None else min(mem_low, value)
        if free.actual_free is not None:
            value = (free.actual_free - held + sum(item.allocation.mem if item.mem_actual is None
                                                   else item.mem_actual for item in ended))
            actual_low = value if actual_low is None else min(actual_low, value)
        for uuid in {g for a in active if a is not None and not a.gpu_exclusive for g in a.gpus}:
            used = (free.gpu_used.get(uuid, 0)
                    - sum(item.allocation.gpu_claim for item in ended if uuid in item.allocation.gpus)
                    + sum(a.gpu_claim for a in active if a is not None and not a.gpu_exclusive and uuid in a.gpus))
            gpu_high[uuid] = max(gpu_high.get(uuid, 0), used)
    if mem_low is not None:
        window.mem_free = mem_low
    if actual_low is not None:
        window.actual_free = actual_low
    for uuid, used in gpu_high.items():
        window.gpu_used[uuid] = max(window.gpu_used.get(uuid, 0), used)
    return window


def _fits_at(request: Request, partition: Partition, pool: Pool, running: list[Running], start: float,
             ahead: list[Reservation]) -> Allocation | None:
    """What *request* would take when it starts at *start* (every running job
    that ends by then has given its resources back; the reservations *ahead*
    keep theirs), or ``None`` when it does not fit then."""
    free = pool.copy()
    later: list[Running] = []
    for item in running:
        if item.end_time <= start:
            free.give(item.allocation, item.mem_actual)
        else:
            later.append(item)
    end = INFINITY if request.time_limit is None else start + request.time_limit
    return _window(free, later, start, end, ahead).try_allocate(request, partition)


def _reserve(request: Request, partition: Partition, pool: Pool, running: list[Running], now: float,
             ahead: list[Reservation]) -> Reservation:
    """The reservation of a blocked request behind the reservations *ahead*.

    With none ahead this is the 0.8.0 shadow time (:func:`earliest_start`).
    Otherwise the request starts at the earliest moment (now, the end of a
    running job, or the end of a reservation ahead) at which it fits for its
    whole time limit beside the running jobs and the reservations ahead."""
    placed = [r for r in ahead if r.allocation is not None]
    if not placed:
        # the shadow time exactly as 0.8.0 computes it; what the job takes
        # then is what it gets once every job ending by then has ended
        shadow = earliest_start(request, partition, pool, running, now)
        allocation = None if shadow == INFINITY else _fits_at(request, partition, pool, running, shadow, [])
        return Reservation(request, partition, shadow, allocation)
    order = sorted(running, key=lambda r: r.end_time)
    moments = sorted({now}
                     | {max(now, r.end_time) for r in order if r.end_time != INFINITY}
                     | {max(now, r.end) for r in placed if r.end != INFINITY})
    free = pool.copy()
    index = 0
    limit = request.time_limit
    for moment in moments:
        while index < len(order) and order[index].end_time <= moment:
            free.give(order[index].allocation, order[index].mem_actual)
            index += 1
        end = INFINITY if limit is None else moment + limit
        allocation = _window(free, order[index:], moment, end, placed).try_allocate(request, partition)
        if allocation is not None:
            return Reservation(request, partition, moment, allocation)
    return Reservation(request, partition, INFINITY, None)


def _first_delayed(request: Request, allocation: Allocation, end: float, pool: Pool, running: list[Running],
                   reservations: list[Reservation]) -> tuple[int | None, list[Reservation]]:
    """Whether starting *request* now (until *end*) would delay a reservation.

    Returns the first reserved job that could no longer start at its reserved
    time, or ``None`` and the reservations re-placed around the started job.
    A reservation that starts after the job ends -- and after it every
    earlier reservation starts as well -- sees the same resources as before
    and keeps its place."""
    trial = pool.copy()
    trial.take(allocation)
    trial_running = [*running, Running(request.job_id, allocation, end)]
    placed: list[Reservation] = []
    horizon = INFINITY
    for reservation in reservations:
        if reservation.shadow == INFINITY:
            placed.append(reservation)  # no start time to keep (0.8.0: anything <= infinity)
            continue
        horizon = min(horizon, reservation.shadow)
        if end < horizon:
            placed.append(reservation)
            continue
        moved = _fits_at(reservation.request, reservation.partition, trial, trial_running, reservation.shadow,
                         placed)
        if moved is None:
            return reservation.request.job_id, reservations
        placed.append(Reservation(reservation.request, reservation.partition, reservation.shadow, moved))
    return None, placed


def plan(
    config: NodeConfig,
    pool: Pool,
    running: list[Running],
    queue: list[Request],
    now: float,
    *,
    backfill: bool = True,
    depth: int | None = None,
) -> Decision:
    """Choose the queued jobs that start now.

    *queue* must already be in dispatch order.  *pool* and *running* are not
    modified.  *depth* (0.8.1): how many blocked jobs hold a reservation;
    ``None`` takes ``scheduler.reservation_depth`` of the node configuration.
    """
    depth = config.scheduler.reservation_depth if depth is None else max(1, int(depth))
    pool = pool.copy()
    running = list(running)
    starts: list[tuple[Request, Allocation]] = []
    skipped: dict[int, str] = {}
    delays: dict[int, int] = {}
    reservations: list[Reservation] = []
    # shape -> the reserved job a request of that shape would delay; valid
    # until the next start (more reservations do not change the answer)
    refused: dict[tuple, int] = {}
    held_partitions: set[str] = set()
    for request in queue:
        # The request names a partition or an alias (0.6.3); the alias view
        # carries the alias's preferred cores, the hold is per partition.
        partition = config.find_partition(request.partition)
        if partition is None:
            skipped[request.job_id] = f"unknown partition {request.partition}"
            continue
        if partition.name in held_partitions:
            skipped[request.job_id] = "held behind a paused ensemble"
            continue
        if request.blocked:
            held_partitions.add(partition.name)
            skipped[request.job_id] = "paused ensemble"
            continue
        reason = request_is_feasible(request, partition, config)
        if reason is not None:
            skipped[request.job_id] = reason
            continue
        allocation = pool.try_allocate(request, partition)
        if allocation is not None:
            end = INFINITY if request.time_limit is None else now + request.time_limit
            if not reservations:
                pool.take(allocation)
                running.append(Running(request.job_id, allocation, end))
                starts.append((request, allocation))
                continue
            if not backfill:
                skipped[request.job_id] = "waiting behind a higher-priority job"
                continue
            key = _shape(request)
            delayed = refused.get(key)
            if delayed is None:
                delayed, placed = _first_delayed(request, allocation, end, pool, running, reservations)
                if delayed is None:
                    pool.take(allocation)
                    running.append(Running(request.job_id, allocation, end))
                    starts.append((request, allocation))
                    reservations = placed
                    refused.clear()
                    continue
                refused[key] = delayed
            skipped[request.job_id] = f"backfill would delay job {delayed}"
            delays[request.job_id] = delayed
            continue
        skipped[request.job_id] = "resources busy"
        if len(reservations) < depth:
            reservations.append(_reserve(request, partition, pool, running, now, reservations))
    return Decision(
        starts=starts,
        reservation=None if not reservations else (reservations[0].request.job_id, reservations[0].shadow),
        skipped=skipped,
        reservations=tuple((r.request.job_id, r.shadow) for r in reservations),
        delays=delays,
    )
