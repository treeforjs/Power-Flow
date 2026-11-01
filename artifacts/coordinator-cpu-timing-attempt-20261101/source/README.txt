Coordinator CPU performance candidate: bounded timing attempt, 2026-11-01

Functional correctness passed; timing is NOT qualified.
SimQueue job 351 executed six three-step manufactured million-carrier runs on owner
core 6, 4 GiB cap, no GPU. All three alternating baseline/candidate comparisons were
bit-identical; all recorded source/dependency/binary integrity checks passed.
All six runs exceeded the unchanged 1% SMT sibling busy-time limit: CPU 14 activity
was 1.21%-1.83% over the complete run windows. Zero pairs qualify. Raw timing
values are retained for audit only; no speedup or production runtime claim is
made. No automatic repeat was submitted.

Numerical candidate c2 is unchanged from the separately sealed functional artifact
coordinator-v8.0.2-cpu-20261031. Baseline is qualified coordinator-v8.0.1/c2.
The intended timing scope was the manufactured 1e6-carrier, 200 um load fixture,
steps 2 and 3, all-CPU backends; it would not establish full physical-run speed.
Physical-memory admission passed before every run: minimum MemAvailable was
24055660544 bytes, above the 8 GiB required reserve+job ceiling. Study cooperative
sibling lock was held; actual /proc/stat activity still exceeded the limit.
The frozen PID namespace task-list diagnostic cannot identify external sibling
workers. Host CPU counters do include that load and conservatively reject it.

Launcher corrections and evidence
  Job 338 was canceled while queued because SQLite scheduler metadata is hidden
  by the frozen runtime's temporary ~/.local mount. Host-side read-only scheduler
  verification now precedes the unchanged frozen numerical launch.
  Job 350 failed before numerical work: the command can start while SimQueue's row
  is still starting. Store allocation occurs before later running reconciliation.
  The corrected launcher waits at most 10 s for starting->running while rechecking
  stable request, unit/cgroup, owner core/lease/affinity and actual memory.max on
  every poll. Running is still required before frozen entry. Job 351 recorded this
  transition at 1.008 s. Cancellation/identity changes and timeout refuse.
  Frozen phase rechecks actual unit/cgroup/affinity/memory and exact launcher and
  comparison pins. Sibling-lock acquisition remains bounded 30 s; actual quiet-SMT
  and physical-memory gates are unchanged. Fifteen metadata regressions passed.
  Failed 350 source/log and sanitized queued-job cancellation metadata are retained.

No precision, physics, schedule, particles, guards, source model or numerical
binary was changed for the launcher corrections. Source/dependency pins of c2
remain those of its original functional seal. This artifact contains only timing
harness source and selected evidence, not a promotion to physical plasma results.
