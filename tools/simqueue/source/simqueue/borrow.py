"""Opportunistic borrowing of idle cores (SimQueue 0.6.1): pure planning.

Terms
-----
owner job
    A job submitted without ``--borrow``.  It owns the cores the scheduler
    gives it (its unit holds their exclusive lease files) and is never
    preempted.
borrow job (borrower)
    A job submitted with ``--borrow`` (restartable work such as builds and
    tests).  It runs in the idle-weight borrow slice with ``CPUWeight=1`` and
    ``SCHED_IDLE``, holds only the shared borrow lease of its cores, and is
    preempted (freeze, requeue or cancel) as soon as an owner needs a core.
lendable core
    A physical core a borrow job may use now, by class:

    ``free``         no queue job and no lease holder, and the CPU measured
                     idle (work pinned there without a lease stays alone);
                     never lent while an owner job of its partition is waiting
    ``leased-idle``  core lease held outside the queue and the holder's
                     process tree measured idle (partition ``lend_leased_idle``)
    ``owner-idle``   allocated to a running owner job that leaves it unused
                     (partition ``lend``; never for ``--exclusive-timing`` jobs)

Everything here is pure: the worker passes the measured core states and the
current jobs, and receives the starts, preemptions and resumes to carry out.
With no borrow jobs, :func:`plan_node` returns exactly what
:func:`scheduler.plan` returns.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, Mapping

from .config import format_cpu_list
from .nodeconfig import NodeConfig
from .scheduler import Allocation, Decision, Pool, Request, Running, plan


CLASS_RANK = {"free": 0, "leased-idle": 1, "owner-idle": 2}
# Skip reasons of owner jobs that are waiting for cores (as opposed to being
# held by a paused ensemble or infeasible).
WAITING_REASONS = ("resources busy", "waiting behind a higher-priority job", "backfill would delay")
# A frozen borrower is requeued instead when its time limit is this close,
# so that systemd never has to end a frozen unit at its time limit.
FROZEN_DEADLINE_MARGIN_S = 60.0


@dataclass(frozen=True)
class CoreState:
    """Measured state of one physical core (by its partition CPU number)."""

    cpu: int
    external: bool = False  # core lease held outside the queue
    external_idle: bool = False  # ... and the holder tree is idle over the idle window
    external_busy: bool = False  # ... and the holder tree uses more than it leaves to borrowers
    stale: bool = False  # lease held but no holder process is visible
    holder: str = ""  # description of the holder, for messages
    holder_key: str = ""  # identity of the holder (changes when another process takes the lease)
    owner_idle: bool = False  # the owner job allocated this core leaves it unused
    owner_busy: bool = False  # the owner job needs its lent cores back
    no_borrow: bool = False  # borrow lease held exclusively outside the queue
    cpu_busy: bool = False  # the CPU itself is busy (or not measured yet), e.g. unleased work pinned there


@dataclass
class Borrower:
    """A running, frozen or stopping borrow job."""

    job_id: int
    home: str
    cpus: tuple[int, ...]
    mem: int
    on_preempt: str = "requeue"
    placement: dict[int, str] = field(default_factory=dict)  # cpu -> class at placement
    owners: dict[int, int] = field(default_factory=dict)  # cpu -> owner job (owner-idle)
    holder_keys: dict[int, str] = field(default_factory=dict)  # cpu -> holder (leased-idle)
    frozen_since: float | None = None
    pending: str | None = None  # "requeue" / "cancel" while its unit is being stopped
    borrow_from: tuple[str, ...] | None = None
    deadline: float | None = None  # time-limit end (unix seconds)
    order: int = 0  # dispatch order; lower resumes first and is stopped last


@dataclass(frozen=True)
class Preemption:
    job_id: int
    mode: str  # freeze | requeue | cancel
    reason: str
    by_job: int | None = None
    escalated: bool = False  # freeze turned into requeue
    cpus: tuple[int, ...] = ()


@dataclass
class NodeDecision:
    starts: list[tuple[Request, Allocation]]
    reservation: tuple[int, float] | None
    skipped: dict[int, str]
    preemptions: list[Preemption] = field(default_factory=list)
    resumes: list[tuple[int, tuple[tuple[int, str, str], ...]]] = field(default_factory=list)
    lendable: dict[int, str] = field(default_factory=dict)  # diagnostic: cpu -> class, after this pass


def protected_cpus(config: NodeConfig, exclusive_owner_cpus: Iterable[int] = ()) -> set[int]:
    """Logical CPUs no borrower may use: timing-sensitive partitions and the
    cores of running ``--exclusive-timing`` jobs, each with its SMT siblings."""
    protected: set[int] = set()
    for part in config.partitions.values():
        if part.timing_sensitive:
            protected |= config.with_siblings(part.cpus)
    protected |= config.with_siblings(exclusive_owner_cpus)
    return protected


def eligible_partitions(config: NodeConfig, request: Request) -> list[str]:
    """Partitions whose cores a borrow job may use, in preference order."""
    home = request.partition
    names = [home]
    others = request.borrow_from if request.borrow_from is not None else tuple(config.partitions)
    for name in others:
        part = config.partitions.get(name)
        if part is None or name in names or not part.lends():
            continue
        names.append(name)
    home_part = config.partitions.get(home)
    if home_part is not None and home_part.timing_sensitive:
        names.remove(home)
    return names


def borrow_feasibility(config: NodeConfig, request: Request) -> str | None:
    """Why a borrow request can never run on this node, else ``None``."""
    if request.partition not in config.partitions:
        return f"unknown partition {request.partition}"
    if request.gpus:
        return "borrow jobs cannot use GPUs"
    if request.cpu_list is not None:
        return "borrow jobs cannot name explicit CPUs"
    reachable = sum(
        len([c for c in config.partitions[name].cpus if c not in protected_cpus(config)])
        for name in eligible_partitions(config, request)
    )
    if request.ncpus < 1 or request.ncpus > reachable:
        return f"needs {request.ncpus} core(s); the partitions it may use have {reachable}"
    if config.memory_schedulable is not None and request.mem > config.memory_schedulable:
        return "memory request exceeds the node's schedulable memory"
    return None


def plan_node(
    config: NodeConfig,
    pool: Pool,
    running: list[Running],
    queue: list[Request],
    now: float,
    *,
    borrowers: Iterable[Borrower] = (),
    cores: Mapping[int, CoreState] | None = None,
    owner_of: Mapping[int, int] | None = None,
    exclusive_owner_cpus: Iterable[int] = (),
    backfill: bool = True,
    starts_allowed: bool = True,
) -> NodeDecision:
    """Owner starts (EASY backfill, as in 0.6.0), then borrower preemptions,
    resumes and starts.

    *pool* holds the free resources after owner allocations and outside
    lease holders; cores and memory of borrowers are *not* subtracted (they
    are reclaimable for owners).  *running* lists owner jobs only.
    *owner_of* maps each core of a running owner job to that job.
    """
    borrowers = list(borrowers)
    cores = dict(cores or {})
    owner_of = dict(owner_of or {})
    owner_queue = [r for r in queue if not r.borrow]
    borrow_queue = [r for r in queue if r.borrow]
    by_id = {b.job_id: b for b in borrowers}
    occupied: dict[int, int] = {}
    for borrower in borrowers:
        for cpu in borrower.cpus:
            occupied[cpu] = borrower.job_id

    owner_pool = pool.copy()
    owner_pool.soft_cpus = set(occupied) & owner_pool.free_cpus
    if starts_allowed and owner_queue:
        base = plan(config, owner_pool, running, owner_queue, now, backfill=backfill)
    else:
        base = Decision(starts=[], reservation=None, skipped={})
    if not borrowers and not borrow_queue:
        return NodeDecision(starts=base.starts, reservation=base.reservation, skipped=base.skipped)

    preemptions: dict[int, Preemption] = {}

    def preempt(borrower: Borrower, reason: str, *, by_job: int | None = None, stop: bool = False) -> None:
        """Record a preemption; *stop* turns a freeze into a requeue."""
        if borrower.pending:
            return
        mode = borrower.on_preempt
        escalated = False
        if stop and mode == "freeze":
            mode, escalated = "requeue", True
        previous = preemptions.get(borrower.job_id)
        if previous is not None:
            if previous.mode == "freeze" and mode != "freeze":
                preemptions[borrower.job_id] = Preemption(
                    borrower.job_id, mode, reason, by_job, escalated, borrower.cpus
                )
            return
        if mode == "freeze" and borrower.frozen_since is not None:
            return  # already frozen
        preemptions[borrower.job_id] = Preemption(borrower.job_id, mode, reason, by_job, escalated, borrower.cpus)

    exclusive = set(exclusive_owner_cpus)
    started_owner: dict[int, int] = {}
    for request, allocation in base.starts:
        for cpu in allocation.cpus:
            started_owner[cpu] = request.job_id
        if request.exclusive_timing:
            exclusive |= set(allocation.cpus)
    protected = protected_cpus(config, exclusive)

    # 1. Owners never wait: borrowers on the cores (or, for --exclusive-timing,
    #    the cores and SMT siblings) of a starting owner are preempted now.
    for request, allocation in base.starts:
        guard = config.with_siblings(allocation.cpus) if request.exclusive_timing else set(allocation.cpus)
        victims = sorted({occupied[cpu] for cpu in guard if cpu in occupied})
        for job_id in victims:
            borrower = by_id[job_id]
            shared = sorted(set(borrower.cpus) & guard)
            preempt(borrower, f"CPU {format_cpu_list(shared)} needed by job {request.job_id}", by_job=request.job_id)

    # 2. Memory: owners count borrower memory as free.  When starting owners
    #    needs it, stop borrowers (newest first); freezing keeps memory, so a
    #    freeze becomes a requeue here.
    if pool.mem_free is not None and base.starts:
        owner_new = sum(allocation.mem for _, allocation in base.starts)

        def staying() -> list[Borrower]:
            return [
                b for b in borrowers
                if not b.pending and not (b.job_id in preemptions and preemptions[b.job_id].mode != "freeze")
            ]

        free_after = pool.mem_free - owner_new - sum(b.mem for b in staying())
        if free_after < 0:
            first_owner = base.starts[0][0].job_id
            for borrower in sorted(staying(), key=lambda b: -b.order):
                if free_after >= 0:
                    break
                preempt(borrower, f"memory needed by job {first_owner} (a frozen job keeps its memory)",
                        by_job=first_owner, stop=True)
                free_after += borrower.mem

    # 3. Measured state: a lease holder that becomes busy, a lease taken on a
    #    lent core, an owner that needs its idle cores back, a policy change.
    for borrower in borrowers:
        if borrower.pending:
            continue
        if borrower.frozen_since is not None:
            if config.borrow.freeze_max_s is not None and now - borrower.frozen_since >= config.borrow.freeze_max_s:
                preempt(borrower, f"frozen for more than {int(config.borrow.freeze_max_s)} s", stop=True)
                continue
            if borrower.deadline is not None and now >= borrower.deadline - FROZEN_DEADLINE_MARGIN_S:
                preempt(borrower, "time limit reached while frozen", stop=True)
                continue
        reason = _state_conflict(config, borrower, cores, owner_of, protected)
        if reason is not None:
            preempt(borrower, reason)

    # Cores of borrowers that stay (running or frozen) or are being stopped
    # remain occupied for other borrowers; owners may take them (above).
    owner_cpus = set(owner_of) | set(started_owner)
    waiting: set[int] = set()
    started_ids = {request.job_id for request, _ in base.starts}
    for request in owner_queue:
        if request.job_id in started_ids:
            continue
        reason = base.skipped.get(request.job_id, "")
        if not reason.startswith(WAITING_REASONS):
            continue
        part = config.partitions.get(request.partition)
        if part is not None:
            waiting |= set(request.cpu_list) if request.cpu_list is not None else set(part.cpus)

    def lendable_class(cpu: int, home: str, *, occupant: int | None = None) -> str | None:
        part = config.partition_of(cpu)
        if part is None or cpu in protected or part.timing_sensitive:
            return None
        if part.name != home and not part.lends():
            return None
        if occupied.get(cpu, occupant) != occupant:
            return None
        state = cores.get(cpu, CoreState(cpu))
        if state.no_borrow:
            return None
        if cpu in owner_cpus:
            if cpu in started_owner or cpu in exclusive:
                return None
            if state.owner_idle and part.lend:
                return "owner-idle"
            return None
        if state.external:
            if state.stale or not state.external_idle:
                return None
            if part.lend and part.lend_leased_idle:
                return "leased-idle"
            return None
        if cpu in waiting or cpu not in pool.free_cpus or state.cpu_busy:
            return None
        return "free"

    def why_not(cpu: int, home: str) -> str:
        """Short reason a core is not lendable to a borrow job of *home*."""
        part = config.partition_of(cpu)
        if part is None:
            return "not in a partition"
        if cpu in owner_cpus and cpu in exclusive:
            return f"job {owner_of.get(cpu, started_owner.get(cpu))} (exclusive timing)"
        if part.timing_sensitive or cpu in protected_cpus(config):
            return "timing-sensitive"
        if cpu in protected:
            return "SMT sibling of an exclusive-timing job"
        if part.name != home and not part.lends():
            return "not lent"
        if cpu in occupied:
            return f"used by borrow job {occupied[cpu]}"
        state = cores.get(cpu, CoreState(cpu))
        if state.no_borrow:
            return "borrow lease held"
        if cpu in owner_cpus:
            owner = owner_of.get(cpu, started_owner.get(cpu))
            return f"job {owner}" + (" (exclusive timing)" if cpu in exclusive else "")
        if state.external:
            if state.stale:
                return "lease held, holder not visible"
            if not (part.lend and part.lend_leased_idle):
                return "lease held"
            return "lease held, " + ("idle" if state.external_idle else "busy")
        if cpu in waiting:
            return "an owner job of the partition is waiting"
        if state.cpu_busy:
            return "CPU busy outside the queue"
        return "unavailable"

    # 4. Resume frozen borrowers whose cores are lendable again (oldest first).
    resumes: list[tuple[int, tuple[tuple[int, str, str], ...]]] = []
    if starts_allowed:
        for borrower in sorted(borrowers, key=lambda b: b.order):
            if borrower.frozen_since is None or borrower.pending or borrower.job_id in preemptions:
                continue
            placement = []
            for cpu in borrower.cpus:
                cls = lendable_class(cpu, borrower.home, occupant=borrower.job_id)
                part = config.partition_of(cpu)
                allowed = part is not None and (
                    part.name == borrower.home or borrower.borrow_from is None or part.name in borrower.borrow_from
                )
                if cls is None or not allowed:
                    placement = []
                    break
                placement.append((cpu, part.name, cls))
            if placement:
                resumes.append((borrower.job_id, tuple(placement)))

    # 5. Start borrow jobs on lendable cores, in dispatch order, greedily.
    starts = list(base.starts)
    skipped = dict(base.skipped)
    mem_free: int | None = None
    if pool.mem_free is not None:
        stopping = {job_id for job_id, p in preemptions.items() if p.mode != "freeze"}
        mem_free = (
            pool.mem_free
            - sum(allocation.mem for _, allocation in base.starts)
            - sum(b.mem for b in borrowers if not b.pending and b.job_id not in stopping)
        )
    for request in borrow_queue:
        if not starts_allowed:
            skipped[request.job_id] = "queue paused"
            continue
        if not config.borrow.enabled:
            skipped[request.job_id] = "borrowing is disabled in the node configuration"
            continue
        reason = borrow_feasibility(config, request)
        if reason is not None:
            skipped[request.job_id] = reason
            continue
        names = eligible_partitions(config, request)
        candidates: list[tuple[int, int, int, str, str]] = []  # (rank, partition index, -cpu, ...)
        for index, name in enumerate(names):
            for cpu in config.partitions[name].cpus:
                if config.partition_of(cpu) is not config.partitions[name]:
                    continue  # a CPU listed in two partitions is lent by its first
                cls = lendable_class(cpu, request.partition)
                if cls is not None:
                    candidates.append((CLASS_RANK[cls], index, -cpu, name, cls))
        candidates.sort()
        want = max(request.ncpus, request.max_cpus or request.ncpus)
        chosen = candidates[:want]
        if len(chosen) < request.ncpus:
            details = []
            for name in names:
                reasons: dict[str, list[int]] = {}
                for cpu in config.partitions[name].cpus:
                    if lendable_class(cpu, request.partition) is None:
                        reasons.setdefault(why_not(cpu, request.partition), []).append(cpu)
                if reasons:
                    details.append(name + ": " + "; ".join(
                        f"CPU {format_cpu_list(cpus)} {why}" for why, cpus in reasons.items()))
            skipped[request.job_id] = (
                f"borrow: {len(candidates)} of {request.ncpus} core(s) lendable now"
                + (" (" + " | ".join(details) + ")" if details else "")
            )
            continue
        if mem_free is not None and request.mem > mem_free:
            skipped[request.job_id] = "borrow: memory busy"
            continue
        placement = tuple(sorted((-negative, name, cls) for _, _, negative, name, cls in chosen))
        cpus = tuple(cpu for cpu, _, _ in placement)
        assert not set(cpus) & protected, "borrower placed on protected CPUs"
        allocation = Allocation(cpus=cpus, mem=request.mem, borrow=True, placement=placement)
        for cpu in cpus:
            occupied[cpu] = request.job_id
        if mem_free is not None:
            mem_free -= request.mem
        starts.append((request, allocation))

    lendable = {}
    for cpu in config.node_cpus():
        part = config.partition_of(cpu)
        if part is None:
            continue
        cls = lendable_class(cpu, "")
        if cls is not None:
            lendable[cpu] = cls
    return NodeDecision(
        starts=starts,
        reservation=base.reservation,
        skipped=skipped,
        preemptions=[preemptions[key] for key in sorted(preemptions)],
        resumes=resumes,
        lendable=lendable,
    )


def _state_conflict(
    config: NodeConfig,
    borrower: Borrower,
    cores: Mapping[int, CoreState],
    owner_of: Mapping[int, int],
    protected: set[int],
) -> str | None:
    for cpu in borrower.cpus:
        state = cores.get(cpu, CoreState(cpu))
        cls = borrower.placement.get(cpu, "free")
        part = config.partition_of(cpu)
        if cpu in protected or (part is not None and part.timing_sensitive):
            return f"CPU {cpu} is reserved for timing-sensitive work"
        if part is not None and part.name != borrower.home and not part.lends():
            return f"partition {part.name} no longer lends"
        if state.no_borrow:
            return f"borrow lease of CPU {cpu} is held exclusively outside the queue"
        if state.external:
            if cls != "leased-idle":
                return f"lease of CPU {cpu} taken outside the queue ({state.holder or 'unknown holder'})"
            if borrower.holder_keys.get(cpu) not in (None, "", state.holder_key):
                return f"lease of CPU {cpu} changed hands ({state.holder or 'unknown holder'})"
            if state.external_busy:
                return f"lease holder of CPU {cpu} is busy ({state.holder or 'unknown holder'})"
            if part is not None and not (part.lend and part.lend_leased_idle):
                return f"partition {part.name} no longer lends leased idle cores"
        if cls == "owner-idle":
            owner = borrower.owners.get(cpu)
            if owner is not None and owner_of.get(cpu) == owner and state.owner_busy:
                return f"job {owner} needs CPU {cpu} again"
        current_owner = owner_of.get(cpu)
        if current_owner is not None and borrower.owners.get(cpu) != current_owner:
            return f"CPU {cpu} is allocated to job {current_owner}"
    return None


def in_use_cores(usage: float, threshold: float) -> int:
    """Cores a process tree occupies, from its average usage in cores.

    Usage up to *threshold* above a whole number still counts as that number
    (1.08 cores with threshold 0.10 is one core in use)."""
    return max(0, math.ceil(usage - threshold - 1e-9))
