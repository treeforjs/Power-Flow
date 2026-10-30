cuDSS0.6 direct-solve tuning candidate — 2026-10-29

This isolated source starts from the sealed GPU source. It does not include the
experimental CUDA Graph implementation or change particle/gather code. V1, v2,
shared repositories and production sources remain unchanged. GPU execution is
parent-controlled. The full run is under a manual hold; do not resume/restart it.

Knobs, all set before analysis; defaults preserve the sealed configuration:
  RAIKOU_PIC_CUDSS_REORDER=default|amd             default: default
    default: cuDSS customized METIS nested dissection; amd: CUDSS_ALG_3.
  RAIKOU_PIC_CUDSS_ND_NLEVELS=default|128|256|512|1024
    default: leave library default (documented512). Explicit ND with AMD rejects.
  RAIKOU_PIC_CUDSS_HYBRID_EXECUTE=0|1             default:0
    Optional CPU/GPU placement heuristic, one host thread; hybrid memory stays0.
There is no solve-algorithm sweep:0.6 supports only CUDSS_ALG_DEFAULT for solve.
FP64 SPD, matrix values, two refinement steps, factorization reuse, all outer
residual/timestep/conservation guards, checkpoint structures and ABI are unchanged.
No factorization/pivot/matching/precision or superpanel settings are exposed.

Diagnostics (use only for an initial inspection run, then disable for timings):
  RAIKOU_PIC_CUDSS_AUDIT_PREFIX=/absolute/new-output/poisson
    Writes poisson.0.json and poisson.1.json for the two cached systems. Files
    contain effective queried knobs, matrix structure, raw LU_NNZ and memory
    estimates with query return status/byte counts. These are diagnostics only,
    never correctness or successful-run receipts. Use a fresh prefix each run.
  RAIKOU_PIC_CUDSS_DUMP_MATRIX=1
    With audit enabled, also writes native int32 CSR upper offsets/columns,
    FP64 values and the row permutation if retrieval/validation succeeds.
    Native little-endian raw arrays; n and nnz are in the JSON. Compare their
    SHA256 with evidence/cpu-matrix.* to establish exact CPU reconstruction parity.
    No matrix or permutation is changed by dumping it. No files are overwritten.
  cuda_analysis and cuda_factorization are separate nested setup timers.
    cuda_factor_setup includes both and optional diagnostics. Existing solve and
    outer timers remain unchanged; do not add nested timers.

bin/run-candidate PIC_ARGUMENT... sets frozen libraries and one thread. The parent
must acquire testing-core.lock, use CPU6/nice15 and the copied-state replay harness
that preserves a pre-existing stop/manual hold. No GPU command was run here.

Build: frozen CMake, Release -O3 -DNDEBUG, MFEM4.10, --parallel1, CPU6/nice15,
one job/thread under exclusive outputs/pic-profile/testing-core.lock. No dependencies
were installed. Logs, binary/source hashes and ready marker are local to this prefix.

CPU-only utilities:
  build/cudss_tuning_options_tests validates supported settings and rejects bad
    combinations; it does not link or call CUDA/MFEM.
  build/poisson_matrix_audit ABSOLUTE_FRESH_PREFIX reconstructs the sealed serial
    Q1 mesh/foil-eliminated matrix without CUDA, and exports its structure/CSR.
    It does not factor or solve the matrix. Zero boundary values do not change
    the matrix; a native audit array comparison remains the direct parity witness.
  Only these CPU utilities were executed here. Foundation/GPU executables were
    compiled, not executed. See benchmark-plan.md for the bounded parent sweep.

Official0.6 docs already cached under pic-gpu-candidate-20261029/docs; primary links:
  https://docs.nvidia.com/cuda/archive/12.9.1/cudss/types.html
  https://docs.nvidia.com/cuda/archive/12.9.1/cudss/advanced_features.html
  https://docs.nvidia.com/cuda/archive/12.9.1/cudss/functions.html
