cuDSS CUDA Graph candidate v2 — 2026-10-29

Prefix: /home/treefor/projects/runtimes/pic-gpu-candidate-v2-20261029
State: built; seven CPU record tests passed; no GPU execution.
v1 at pic-gpu-candidate-20261029 and the sealed production runtime remain unchanged.

Changes from v1
1. Graph-enabled construction explicitly completes default-stream initialization
   with checked cudaDeviceSynchronize BEFORE cuDSS analysis/factorization on its
   nonblocking stream. This adds one setup-only synchronization per solver.
2. A process-local run record retains every solver, off/fallback/qualification
   outcome, late solve/construction/cleanup failure and normal closure. Atomic
   replacement keeps each JSON snapshot complete. There is no per-solve file I/O.
   The main program begins a fresh run record before validating arguments and
   completes it after output, solver release and parallel execution bookkeeping.
   Exception paths mark the run failed before preserving the existing error exit.

Source diffs
  source-from-v1.patch       Directly reviewable incremental changes
  source-from-sealed.patch   Complete changes relative to the sealed source
  candidate-manifest.json    Observed source/binary hashes and verification
  READY.json                Handoff for isolated GPU qualification

Frozen model/runtime
FP64 SPD matrix, cached factorization, two iterative-refinement steps, original
bitwise first-RHS gate and the outer residual/conservation/timestep guards are
unchanged. Particle/options/checkpoint structs and ABI headers are byte-identical.
Dependencies remain GNU16.2.1, MFEM4.10, cuDSS0.6.0.5, cuBLAS12.8.4.1 and CUDA
runtime12.8.90 from the same existing prefixes. No package was installed.

Graph environment
  RAIKOU_PIC_CUDSS_GRAPH=off|required|auto
    off is the default; required is mandatory for performance acceptance.
    auto still has the v1 experimental unsupported-capture recovery contract;
    this follow-up does not establish opaque cuDSS state recovery after failure.
  RAIKOU_PIC_CUDSS_GRAPH_RECORD=/absolute/output/run-record.json
    Optional diagnostic JSON. The parent directory must already exist.
    The process replaces this file at run start, including off mode.
    Use a unique output path for each process; concurrent writers are unsupported.
  RAIKOU_PIC_FIELD_BACKEND=cudss
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1

bin/run-candidate MODE ABSOLUTE_RECORD_PATH PIC_ARGUMENT...
sets these variables and the original sealed CUDA library paths. It targets v2.
The integration chat owns isolated GPU execution and a successful outer process
receipt. No GPU command was run during this implementation.

Record schema version2
  run_id = process PID plus start timestamp in nanoseconds; solver_id = sequential
  process-local integer. Each solver retains matrix_rows, requested_mode, status,
  reason, node counts, first_rhs_bitwise_qualified, direct_solve_calls,
  graph_replays and closed. Direct counts include the warmup reference solve.
  graph_replays counts completed launches after successful stream synchronization.
  Counts are flushed on state changes/closure, not on every steady solve.

  run_status: running | completed | failed
    Qualified solvers do not imply run completion. Exceptions and cleanup errors
    make failure sticky. A successful complete() requires every solver closed.
  graph_status: pending | off | mixed | fallback | qualified | failed
    Aggregate priority: failure, fallback, pending, then mixed/off/qualified.
    One fallback cannot be hidden by another solver qualifying later.
  all_solvers_graph_qualified: false for empty, failed, off, pending or fallback
    runs; true only if every registered solver has qualified and no run failure.
  all_solvers_requested_required: true only for nonempty all-required solver sets.
  outer_process_receipt_required: always true.

Performance acceptance requires ALL of:
  - independent successful outer process receipt for this exact run/binary;
  - run_status=completed, graph_status=qualified;
  - all_solvers_graph_qualified=true and all_solvers_requested_required=true;
  - expected solver count (two for this checkpoint replay), every solver closed,
    and graph_replays>0 for each;
  - unchanged numerical guards and full field/history/ledger/particle comparisons;
  - repeated unprofiled end-to-end timing against identical baseline/off replays.
A SIGKILL, abort, pre-main failure or failure to write the record cannot be inferred
as success. The outer receipt is mandatory even if the last snapshot says completed.
The small qualification executable deliberately includes an off reference solver,
so its aggregate graph status is mixed; its successful exit is only a small test.

Build and CPU-only verification
Commands used CPU6, nice15, one job/thread, holding exclusive flock:
  /home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock
Frozen launcher:
  /home/treefor/projects/runtimes/arch-base-20261028/run-frozen
Configuration:
  /usr/bin/cmake -S PREFIX/source -B PREFIX/build -DCMAKE_BUILD_TYPE=Release
  -DMFEM_DIR=/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/cmake/mfem
Build:
  /usr/bin/cmake --build PREFIX/build --parallel 1
CPU-only executable:
  PREFIX/build/cudss_graph_record_tests SCENARIO ABSOLUTE_OUTPUT_JSON
Scenarios: off, fallback, qualified, late_failure, outer_failure, unclosed,
invalid_path. It includes only the record header and does not link CUDA/MFEM.
JSON verification: evidence/verify-record-tests.py.
Evidence: logs/configure.log, logs/build.log, evidence/record-test-results.json
and evidence/record-tests/*.json. Native/GPU foundation and graph qualification
executables were compiled but not executed here.

Complexity and benefit
See complexity-audit.md. The baseline offers limited submission-only headroom;
keep graphs optional/default-off until an application gain exceeds timing noise
and justifies the capture, lifetime and error-recovery code. No speedup is claimed.
