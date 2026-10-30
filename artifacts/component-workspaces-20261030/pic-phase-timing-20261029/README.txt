Additive per-step phase timing, 2026-10-29

Ready for isolated checkpoint replay. No GPU workload or production signal
was executed while preparing these builds.

Executables:
  baseline/build/heating_history_pic_gpu
  cpu-v2/build/heating_history_pic_gpu
Both paths are relative to:
  /home/treefor/projects/runtimes/pic-phase-timing-20261029

Use the identical five application arguments, saved checkpoint copy, cudss
backend, sealed CUDA LD_LIBRARY_PATH, and one-thread OpenMP/BLAS environment
used for each uninstrumented variant. Run through the existing frozen
launcher. GPU correctness/performance replay remains unrun.

Each successful run writes OUTPUT/phase_profile.json alongside the original
benchmark_profile.json. Existing benchmark timers remain intact.
The new file contains definitions, aggregate counts/times, every recorded
timestep, actual step indices, physical start times, and accepted dt values.

Only seconds_exclusive form an additive partition. Their sum equals each
step's wall_seconds; additive_error_seconds exposes floating-point summation
roundoff. seconds_inclusive describes nested elapsed phases and overlaps
with child phases. Do not add inclusive values or subtract overlapping
benchmark_profile totals to infer CPU-only time.

The timing window begins before each successful driver while-condition and
ends after that iteration's periodic output checks. The loop condition calls
diagnostics(), which scans particles; its cost belongs to the step.
The final false loop-condition is reported separately under
terminal_loop_condition, with no phantom timestep. Startup, checkpoint load,
unit-bias setup, initial progress, final snapshot/save/progress, solver
release, and JSON serialization are outside step totals.

Phases:
  driver_selection: field/boundary updates and initial timestep constraints.
  thermal_source_reservoirs: thermal interpolation, source timestep limits,
    reservoir accumulation, and injection.
  diagnostics: coarse diagnostics() calls, including kinetic-energy scans,
    even when nested under selection or periodic output.
  drift_first / drift_second: the full half-drift calls.
  solve_fields_first / solve_fields_second: complete field solves, with
    nested CUDA, linear CPU residual guard, and any factor setup separated.
    Their exclusive values retain deposition, field assembly, recovery,
    reaction/energy work, and field-level diagnostic residual work.
  force_gather_boris: one scope around the entire particle force loop,
    containing self/external gathers, imposed fields, and Boris updates.
  cuda_solve_first / cuda_solve_second: RHS transfer, unchanged full cuDSS
    SOLVE (including both refinement iterations), and answer transfer.
  cpu_residual_first / cpu_residual_second: LinearSystem residual vector and
    norms, solver bookkeeping, and unchanged numerical guards.
  advance: remaining advance validation, transfer reset, and clock update.
  advance_empty: whole empty-vacuum advance.
  periodic_output: conditional frame/checkpoint/progress work inside a step.
  driver_step: inclusive whole step; exclusive remaining driver glue.
  cuda_solve_other_step / cpu_residual_other_step: unexpected step-time
    linear solves outside the two advance field scopes, tracked explicitly.
  cuda_factor_setup_in_step: factorization/setup only if it actually occurs
    in a driver iteration. Startup/setup calls cannot leak into these totals.

No timer was inserted inside a per-particle loop, field callback, or Boris
function. Timing state is separate from physical/checkpoint data. Per-step
records remain in memory and are serialized after completion; use bounded
replays. Clock calls and scope accounting have small instrumentation cost.
The row-copy bookkeeping between iterations is outside physical step totals.

Validation:
  Both Release builds passed with existing frozen compiler/dependencies.
  Builds used CPU6, nice15, one compile job, and exclusive testing-core.lock.
  The CPU-only accounting test passed nested exclusive sums, first/second
    solve routing, idempotent stop, phase counts, startup exclusion, and the
    separate terminal condition. Its JSON parsed and matched those checks.
  Removing marked timing scaffolding and the while-condition comma hook
    restores each original source apart from whitespace. Physics/guards,
    callback bodies, invocation order, checkpoint layout, tolerances, and
    refinement configuration were not changed.

Review evidence:
  evidence/baseline.diff
  evidence/cpu-v2.diff
  evidence/instrumentation-validation.json
  evidence/original-sources.json
  evidence/accounting-test.json
  validate_phase_json.py
  manifest.json

Validate any completed replay phase file without rerunning the simulation:
  taskset -c 6 nice -n 15 python \
    /home/treefor/projects/runtimes/pic-phase-timing-20261029/validate_phase_json.py \
    /absolute/replay-output/phase_profile.json
The validator checks per-step and aggregate exclusive sums, finite/nonnegative
times, phase counts, step ordering, accepted dt, and the separate terminal check.

Build invocation:
  taskset -c 6 nice -n 15 \
    flock <workspace>/outputs/pic-profile/testing-core.lock \
    /bin/bash /home/treefor/projects/runtimes/pic-phase-timing-20261029/build.sh
