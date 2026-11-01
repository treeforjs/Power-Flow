Transient SimQueue timing shield — source and simulated tests, 2026-11-01

Live canary: UNEXECUTED. No privileged execution, cgroup/system configuration
changes, queue submissions or numerical benchmarks were performed.

The helper accepts only a fixed 75-second waiting canary in a normally allocated
CPU6 SimQueue owner job (4 GiB, no GPU, exclusive-timing). It preserves the job
unit, lease, CPU affinity and memory cap while a separate watchdog would hold
a transient remote isolated cpuset for CPU6/14, then compare-and-restore only
its own matching writes. CANARY.txt describes staging, prerequisites and
acceptance; this artifact does not establish kernel behavior or a speedup.

The frozen-interpreter metadata suite passed 32 tests. Tests simulate all cgroup
writes in memory, including every setup-write failure before/after its effect,
audit storage failure, allocation/conflict/inode changes and cleanup. One test
sends TERM solely to its newly forked fake-watchdog child. No other process
is signaled. SIGKILL of the watchdog is not independently recoverable.

Persistent audit failure does not stop safely owned restoration. A changed
value/inode or failed undo stops ancestor cleanup and requires inspection.
The existing timing load/memory thresholds and IRQ affinities are unchanged.
Public output redacts full launch tokens; raw root journals must not be
published. The administrative credential is currently unavailable, so live
qualification remains pending.

Runtime source identity: 3597743d14030e4b512ce035740bb9c3af8aeea4fd68daa104a7308b6d07d95a
Helper identity: b06fc68a0c513e1fdfb27e80a164c742d6cae547bb8db4b0ef8642377d8ffd9c
