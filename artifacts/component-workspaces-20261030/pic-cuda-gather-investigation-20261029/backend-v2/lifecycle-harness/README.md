Preparation is complete only when the fresh harness hashes and CPU build are
published. GPU dispatch remains on HOLD until independent review accepts this
harness and the parent dispatches it. The accepted v2 core is unchanged.

The shared probe performs exactly:

1. Original cuDSS solve with the original finite/residual guard.
2. Force on all 46,726 real captured midpoint particles, fresh original CPU
   E/B/u reference, original componentwise/weighted gates, transactional commit.
3. The same cuDSS system solves again while the backend remains alive.
4. A second fully checked force on the immutable midpoint fixture.
5. Explicit backend close, empty cleanup issue list, own primary reference
   release, unchanged current shared-primary identity and active state.
6. The cuDSS system solves again after close.
7. Separate wrong-current and lost-current probes reject before asynchronous
   transfers or kernel dispatch. Closing each backend establishes the required
   shared primary temporarily, drains/releases its own resources/reference and
   restores the caller's wrong or null context. Only the deliberately created
   private context is destroyed. A cuDSS solve follows each probe.

The Poisson probe uses a manufactured Q1 RHS/operator at the actual 360x360
mesh resolution, original epsilon0, homogeneous horizontal Dirichlet boundaries
and original frozen CudssPoisson implementation. Its RHS is A times the captured
phi with those boundary nodes zeroed. This tests real cuDSS/context lifetime;
it does not claim a real charge-deposition solve or an evolving PIC timestep.
Every solve enforces residual <= max(1e-30, 2e-13*||rhs||2) and finite outputs.

A separate process tests explicit standalone-primary selection, two fully
checked actual-midpoint forces, drain/release, restoration of null current
context and the original primary activity state. Link-time wrappers observe
successful backend retain/release and asynchronous dispatch counts. They do
not interpose cuDSS/runtime shared-library internals. A shared reset/destroy
attempt is blocked and disqualifies the run; no such attempt is intentional.
There is no real device fault/illegal-memory injection.

The third, separately bounded process uses the existing compiled
reference_inclusive_gpu_timing adapter. Its 3 warmups plus 20 measured samples
include fresh original CPU reference generation, GPU input checks/packing,
H2D/kernel/D2H/synchronization, every existing reference gate and transactional
velocity commit in one measured force-pipeline total. Setup/refresh, resetting
the immutable fixture and extra post-timer golden checks are excluded and
labeled. All original API gates are still executed inside the measured total.
This is a fixed-state force microbenchmark, not an integrated application
speedup. The separate native 210-step CPU7 application mean (23.421 ms force)
and the matched fixed-fixture CPU6 timings are not interchangeable baselines.

All cleanup sinks live outside the outer try/catch and are persisted on success
or error. Any cleanup issue, retained resource, failed reference/residual gate,
missing receipt, context identity mismatch or abnormal child exit disqualifies
the run. Terminate the isolated candidate process; never reset/destroy the
borrowed primary or signal/resume production.

The parent alone may run parent-gpu-request.json after fresh harness review.
It must acquire testing-gpu.lock THEN testing-core7.lock, retain both across
the serial command list, and stop on the first failure. Each command has an
external 60-second timeout and runs CPU7/nice15/one thread through the existing
host pidfd guard, owned frozen launcher and read-only production-hold check.
The host guard verifies production before, throughout and after each process;
only its owned child group may be terminated. No command acquires a lock in
reverse order or modifies production/reservations. Existing host guard and
launcher files are hash-pinned, unchanged dependencies.

Before dispatch, run verify_bundle.py under the CPU6 testing lease and compare
the accepted independent review with artifact-sha256.json. The verifier is
read-only and grants no execution authorization. Request phase output paths
must be fresh; retain all stdout/stderr, host-guard, preflight, cleanup,
comparison, result and timing receipts. GPU qualification remains distinct
from later evolving full-step integration review.

FUTURE-RUNTIME-PROPOSAL.md documents a possible explicit shadow/runtime API
split. It changes no current reference policy or pre-existing scientific guard.
