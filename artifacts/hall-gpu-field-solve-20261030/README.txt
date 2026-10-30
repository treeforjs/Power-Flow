GPU backend for the corrected Hall field solve (gpu1), qualification qual-02 — 2026-10-30

source/ is the backend (runtimes/hall-gpu-resident-v1-20261030/field); qual-02/ holds the summary, report, build and per-mesh
phase receipts (workspace paths written as /home/treefor/sim-workspace/powerflow). The corrected Hall model (commit 8dddec0) is
unchanged: the GPU result enters Model::prepare through a foil::solve overload that replaces only that call, so all of the model's
post-solve computation and guards run unchanged on the GPU increment.

Results (RTX 3080 Ti, sm_86, NVRTC 12.8 with --fmad=false --prec-div=true --prec-sqrt=true --ftz=false, all inside run-frozen):
- Operator: the factored device form is bitwise equal to the sealed matrix-free apply on every probe of all 24 fixtures (dc and
  heated states, dt 10 fs and 1 fs, meshes from 12,300 to 1,904,588 edges); the assembled CSR form is exactly symmetric with
  normwise relative difference <= 5.7e-16.
- Guards on every fixture: independent CPU FP64 residual at most 0.114 of the 2e-13 relative bound; absolute Gauss <= 3.3e-25 C;
  divB <= 7.9e-12 T/m; Gauss increment <= 1.2e-26 C; field energy balance <= 1.2e-21 J; Joule transfer accepted by the thermal
  owner. A 1e-6 change of one component is refused by the residual check.
- Determinism: bit-identical solutions across repeat solves and fresh processes, and equal to a host emulation of the same
  fixed-tree algorithm. Iteration counts equal the CPU solver's (2-5); Lanczos kappa 1.000001-1.013.
- Timing, median ms per solve, one core: 41,676 edges (default mesh) CPU 1.2-2.0 vs GPU 0.18-0.28 resident (0.26-0.36 with b/x
  transfers and sigma refresh); 478,448 edges CPU 21-29 vs GPU 1.0-1.7; 1,904,588 edges CPU 112-175 vs GPU 4.2-7.2. The GPU is
  faster at every accepted mesh size (break-even ~5e3 edges vs the current CPU solver, ~1.2e4 vs a single-thread CSR solver).
- Both review findings of the earlier GPU backend (stream ordering, async error teardown) are closed.

Open: the model's own guards miss a 0.1 % scaling of the increment, so every backend must keep the independent residual check;
bit identity on sm_89 and sm_120 not yet shown; the right-hand side is still assembled on the CPU (needed on the device for a fully
resident step); device context and reduction helpers to be merged into the shared GPU common layer; FP32 cuDSS for the one-time
DC seed solve not pursued.
