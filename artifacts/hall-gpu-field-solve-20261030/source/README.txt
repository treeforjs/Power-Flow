Hall field GPU backend (gpu1), v1 - 2026-10-30

Scope: the corrected Hall implicit edge solve
  A = M + dt/2 (S+Y) + dt^2/4 C^T H C        (field_model.hpp:79-88, commit 8dddec0)
in the scaled variable v = sqrt(M) dE, solved on the GPU with deterministic
FP64 Jacobi CG, with every acceptance check evaluated by the unchanged CPU model.
Sealed sources are only read (snapshotted per run from git and verified blob by
blob). Results and receipts:
  /home/treefor/sim-workspace/powerflow/outputs/hall-gpu-resident-v1-20261030/field/<label>/

Layout
  src/csr_operator.hpp  operator container, host products, fixed-tree reduction
                        emulation, host emulation of the device CG, plain CPU
                        PCG timing baseline, Lanczos estimate, fixture file I/O.
                        No Hall/MFEM headers.
  src/assemble.hpp      assembly from axial::Mesh, dt, cell sigma.
  src/model_seam.hpp    solver seam into Model::prepare (see below).
  src/cpu_main.cpp      hall_field_cpu prep|guard.
  src/cu_driver.hpp     minimal CUDA driver API + NVRTC wrapper.
  src/kernels.hpp       device code (NVRTC at run time).
  src/gpu_cg.hpp        device-resident CG.
  src/gpu_main.cpp      hall_field_gpu.
  tests/run_field.py    receipted runner (build, prep, gpu, guard, report).
  tests/run_all.sh      runs all phases of one label under the leases.

Operator representations (both on the device, both FP64)
  assembled  one scaled CSR matrix, rows I + D^-1/2 (dt/2 (S+Y) + dt^2/4 C^T H C) D^-1/2
             on free edges, identity rows on PEC edges, PEC columns dropped.
             Exactly symmetric (entry formula symmetric in e,f). ~12 entries/row.
             Not bitwise equal to the matrix-free apply: the sealed apply never forms
             A_ef; it rounds C E, then H (C E), then C^T(.), then the 1/sqrt(M)
             scaling, and FP multiplication does not distribute. Differences are at
             rounding level (normwise ~2e-16, componentwise <= ~1e-15 of |A||v|).
  factored   curl rows (4 per face) and transposed curl rows (faces ascending per
             edge), sqrt(M), a = (.5 dt)(S+Y), q = (.25 dt) dt, PEC mask: two kernels
             (face, edge) that reproduce the operation order of the sealed apply
             (field_model.hpp:84-85, axial_mesh.hpp:84-87). Bitwise equal to the
             sealed matrix-free apply on every probe; also less memory traffic, so
             it is the faster device form on refined meshes. Suggested default.
  The preconditioner diagonal, a, q and sqrt(M) use the same expression text as
  Model::prepare; the diagonal is checked bitwise against the one prepare hands
  to the solver.

Solver seam (no source change)
  Model::prepare calls foil::solve(lambda, rhs, diagonal, 200); the sealed solver
  takes std::function, so the lambda needs a user-defined conversion. model_seam.hpp
  declares a constrained function template foil::solve after the headers whose
  solves must stay sealed (DC seed in axial_fields.hpp) and before field_model.hpp;
  for that single call the template is an exact match and wins overload resolution.
  It forwards (sealed apply closure, rhs, diagonal, maxit) to an installed backend.
  Everything before and after the solve inside Model::prepare - RHS in compensated
  words, update, midpoint, Gauss increment, energy bound, admit (absolute Gauss and
  divB) - is the sealed code. Model construction runs with no backend installed,
  which shows that the DC seed solves did not go through the seam. model_seam.hpp
  must be the only way field_model.hpp enters a program (one definition of prepare).

Device CG (gpu_cg.hpp, kernels.hpp)
  Stopping rule of foil::solve: inner loop while it < maxit and ||r|| > 0.2 target,
  target = max(1e-30, 2e-13 ||b||), then true residual ||b - A x|| on the device,
  up to three correction cycles. Determinism: row-per-thread products, every
  reduction over a fixed grid (512 blocks x 256 threads, grid-stride rows), shared
  memory halving trees, one finishing block of 512 threads, no atomics; NVRTC with
  --fmad=false --prec-div=true --prec-sqrt=true --ftz=false; architecture chosen at
  run time (SASS when NVRTC supports the device, else newest lower compute_XX PTX).
  The launch geometry does not depend on the device, so the summation tree is the
  same on every architecture; a host emulation of the same algorithm
  (emulate_device_cg) reproduces the device solution bit for bit.
  Launch structure: one head (target, init, BATCH gated iterations, true residual
  into a separate buffer) captured as a CUDA graph, one 112-byte control readback;
  further gated batches only if needed. Gated kernels are no-ops after the stop,
  so BATCH and graph use do not change the result (checked by the batch sweep).
  sigma(T) refresh per step: only conduction rows change (CSR diagonal, Jacobi
  diagonal, a); refresh_sigma_rows() uploads those values and scatters them.

Acceptance of a GPU increment (hall_field_cpu guard)
  1. RHS and diagonal handed over by Model::prepare equal the stored fixture bitwise.
  2. Independent FP64 residual with the sealed apply closure and foil::norm:
     ||b - A x||_2 <= max(1e-30, 2e-13 ||b||_2), else refused.
  3. Model::prepare continues on the GPU increment: its Gauss-increment and energy
     guards and admit() (absolute Gauss, divB) run; failures throw.
  4. The v2 thermal owner prepares heat from the candidate receipt (Joule operand
     identity check inside Thermal::prepare); the FieldHeatProposal energy bound of
     owner.hpp:20-24 is evaluated on the result.
  Negative controls: a 1e-6 relative change of one component must be refused by
  step 2; whole-increment scalings with step 2 bypassed are recorded to show what
  the model guards alone detect.

Lifecycle rules of the driver wrapper (closing the earlier review findings)
  Ordering: every copy, memset and launch on one created non-blocking stream; no
  legacy default stream; uploads from session-owned pinned buffers; static uploads
  staged and the stream synchronized before the operator is used.
  Teardown: one Session owns all memory, module, graphs, events and the stream;
  close() makes its context current, drains the stream (context synchronize as
  fallback) and only then frees and releases; if completion cannot be shown the
  session is quarantined (nothing freed, process must exit). Cleanup never throws
  and never replaces the original error; issues go to a caller-owned report that
  is written into every gpu.json.

Candidates for common/ (owned elsewhere; not edited here)
  cu_driver.hpp: Session (primary-context retain/push/pop/release, stream, pinned
  and device allocation tracking, ordered quarantining teardown, graph capture),
  compile() (NVRTC runtime architecture choice, fixed FP options, image hashing),
  DeviceInfo. The fixed-tree reduction constants and block_tree/final_tree device
  helpers, and the ArrayFile fixture format, are also shareable.

Launch (outputs, leases and core policy as of 2026-10-30)
  build:  core-lease-20261030/both-cores.sh /usr/bin/python -I tests/run_field.py build LABEL
  phases: tests/run_all.sh LABEL MESH...   (prep, gpu, guard per mesh, then report)
  GPU phases hold testing-gpu.lock, then one core lease (any-core.sh); CPU phases
  one core lease. Every compile and run step runs inside the pinned frozen launcher
  arch-base-20261028/run-frozen, which binds /dev; its root carries the same
  user-space driver as the host module (615.78.08), so no step ran outside it.

Results of record: label qual-02 (SUMMARY.txt there). tests/verify_label.py LABEL
recomputes every sha256 recorded in a label's receipts (read-only).
common/ now has its own device primitives (owned elsewhere); reconciling
cu_driver.hpp and the reduction helpers with it is left to that owner's interface.

Open items
  - Cross-architecture bit identity (sm_89, sm_120) not yet shown; the launch
    geometry and summation tree are device-independent by construction.
  - Model::prepare's own guards do not detect a 0.1 % solver error on
    near-equilibrium steps; every backend must keep the independent residual gate.
  - The RHS is still built on the CPU in compensated words; a fully resident step
    needs it on the device (b/x transfer costs are reported separately).
  - Small-mesh latency floor ~0.13-0.2 ms (about 25 launches in one graph plus one
    control readback); kernel fusion (redundant per-block final reductions) not done.
  - Optional FP32 cuDSS preconditioner for Poisson-type solves not attempted: the
    Hall DC seed solve runs once per run inside Model construction.
