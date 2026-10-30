Bounded observational RHS diagnostic, 2026-10-29

Ready for a separately controlled copied-state qualification replay.
No PIC or GPU execution occurred during construction. Production remains
untouched. This is an observer on the stable CPU-v2 source; no solver,
precision, refinement, tolerance, field/particle callback, timestep rule,
guard, or checkpoint-layout change was implemented.

Diagnostic executable:
  /home/treefor/projects/runtimes/pic-rhs-observer-20261029/build/heating_history_pic_gpu
Original source:
  /home/treefor/projects/runtimes/pic-cpu-candidate-v2-20261029/source

The environment variable RAIKOU_PIC_RHS_OBSERVE_DIR enables observation.
If unset, the observer produces no files. Its directory must be absolute and
new or empty; existing capture data is never overwritten. The helper checks
this before starting the application. Normal PIC output remains separate.

Parent-controlled launch with the same five PIC arguments and copied state:
  /home/treefor/projects/runtimes/pic-rhs-observer-20261029/bin/run-observer \
    /absolute/new-capture-directory \
    FIELD_DIR THERMAL_BINARY COPIED_PIC_OUTPUT END_NS VOLTAGE_CSV

The wrapper supplies the sealed CUDA libraries, cudss backend, and one-thread
OpenMP/BLAS environment. It opens the frozen runtime itself. Use the ordinary
direct binary plus RAIKOU_PIC_RHS_OBSERVE_DIR when the surrounding execution
harness already supplies that runtime/environment. No process/cgroup control
is built into the observer or wrapper.

Selection is the first eight valid populated GapPic::advance entries, based
on the actual particle count at advance entry. Both original half drifts and
both original solve_fields calls execute in their original order. Startup,
checkpoint-load reconstruction, unit-bias solves, and empty-vacuum advances
are outside observation. Later advances execute normally without more copies.

Each selected call captures exactly six actual FP64 vectors:
  projected_previous_phi: actual trial immediately BEFORE FormLinearSystem.
  full_charge_rhs: existing full_rhs assembled-source copy before elimination.
  constrained_rhs: actual b after FormLinearSystem, before the backend.
  incoming_zeroed_X: actual X after the default copy_interior=0 operation.
  accepted_solution: actual answer AFTER the existing CPU guards passed.
  guard_residual: existing CPU error vector from that guard.

The projected prior is the actual full trial, not a reconstructed vector and
not the raw unprojected potential_ field. Original ProjectBdrCoefficient runs
before the finite-foil code replaces the final essential list. It can modify
DOFs outside that final mask. Consequently second-half projected prior values
need not equal the previous answer on every nonessential DOF. The files retain
those differences and the validator reports their counts without inventing
or silently correcting them. No seventh vector is exported.

One exact native constrained MFEM CSR is captured in its original row/column
order, including stored zeros, together with the exact int32 essential list
and a uint8 essential mask. The matrix/list are checked bit for bit on later
selected calls. Metadata records actual step index, populated-step ordinal,
half, driver absolute start time and dt, advance ledger time, current lower/
upper potentials, particle count at entry, real linear-system ID and addresses,
field-state/matrix identities, and the existing residual/RHS norms.
No extra solver getter, solve, residual norm, or matrix-vector evaluation runs.

The capture buffers are bounded. Both calls' six vectors are reserved before
a populated step is recorded; matrix/list/mask storage is reserved separately.
Payload must fit within 126 MiB and the JSON manifest within a reserved 2 MiB,
so total capture output cannot exceed 128 MiB (134217728 bytes).
Capacity refusal disables further observation, preserves scientific execution,
and produces status=capacity_limited rather than claiming a complete fixture.
A complete fixture requires eight populated steps and sixteen checked calls.
Incomplete steps are removed from the capture.

For n=130321, six vectors per call across sixteen solves require 100086528
bytes. Native CSR/list/mask add about 14-15 MB. The runtime checks actual
storage sizes and the final exact byte count; this estimate is not the cap.
Memory capture occurs at the hooks; files are written only after the original
final snapshots/save/progress and original benchmark/resource reports. This
is an observational diagnostic and must not be used as a timing benchmark.

Output:
  capture.json
  matrix.rows.i32, matrix.cols.i32, matrix.values.f64
  essential.i32, essential_mask.u8
  call_0001_*.f64 through call_0016_*.f64

Post-run CPU-only validation and SHA256 receipt:
  taskset -c 6 nice -n 15 \
    flock -s -w 1 /home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock \
    python /home/treefor/projects/runtimes/pic-rhs-observer-20261029/validate_capture.py \
    /absolute/new-capture-directory \
    --audit-prefix /home/treefor/sim-workspace/powerflow/outputs/pic-profile/cudss-tuning-validation/v2-default-audit/poisson.0 \
    --write-receipt

The validator checks complete selection, finite vectors, exact essential
mask/current essential seed, zero incoming interiors, step/half/time metadata,
original guard values, payload sizes, and the directory's total byte cap.
The optional audit comparison filters/sorts only an upper-storage view of
the captured native coefficients and compares it bit for bit with the existing
native cuDSS audit. It never assembles a matrix, invents an RHS/prior, solves,
or invokes a GPU. Its bounded validation.json includes SHA256 for every raw
file and capture.json; it is also included in the total byte-cap check.

Verification performed here:
  Frozen Release diagnostic build with original compiler/dependencies.
  One compile job, CPU6/nice15, exclusive shared testing-core.lock lease.
  Standard-library CPU I/O contracts under CPU/address-space limits:
    first-eight selection, startup/empty exclusion, deep pre-form copies,
    exact CSR zeros/mask, cap refusal, oversized-dimension refusal, no overwrite,
    and native matrix mutation detection.
  CPU-only native MFEM FormLinearSystem contract under resource limits:
    real zero-interior overwrite and projection outside the final essential
    mask, with a known zero-RHS answer; no solver/PIC/GPU was called.
  Captures in evidence/contracts-* are explicitly marked synthetic contract
    fixtures and rejected by the validator unless that test mode is requested.
  Removing marked observation scaffolding restores all original source tokens
    and operation ordering. CudssPoisson source/configuration is byte unchanged.

Evidence:
  evidence/source.diff
  evidence/source-removal-proof.json
  evidence/original-sources.json
  evidence/build-manifest.json
  evidence/latest-contract-root.txt

The real 1 ps copied-state GPU replay and unchanged full saved-state comparison
remain pending. No solver skipping, deployment, production resume/unfreeze,
new solver package, or precision/refinement change is authorized by this build.
