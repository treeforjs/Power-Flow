cuDSS tuning v2: setup diagnostics only — 2026-10-29

Base: /home/treefor/projects/runtimes/pic-cudss-tuning-20261029
V1 source, binaries and saved screening evidence remain unchanged.
No GPU execution or production signal occurred here. Full production stays frozen.

The tuning variables and permitted values are unchanged from v1. This revision
fixes only setup metadata handling: every ConfigSet/Get is labelled with parameter,
native scalar type, return code, supplied buffer bytes, raw sizeWritten and raw
bytes. ConfigSet has no sizeWritten API argument, so that field is null for setters.
DataGet diagnostics are also logged, with at most32 raw bytes for larger arrays.

RAIKOU_PIC_CUDSS_AUDIT_PREFIX=/absolute/fresh/output/poisson
creates poisson.0.config.jsonl and poisson.0.json, then .1 files for the next solver.
The JSONL flushes each call so a later setter/runtime failure retains its evidence.
Get buffers start with0xa5. The setup JSON reports a getter value only if return
status is SUCCESS and sizeWritten equals the exact requested scalar/array size.
Otherwise it reports null. A successful query with unexpected byte count, or an
unsupported/failed diagnostic query, is recorded without aborting the solve.
Raw bytes are evidence, not a substitute for the API's size/return contract.
Diagnostic I/O errors still fail explicitly; use a fresh writable output prefix.
RAIKOU_PIC_CUDSS_DUMP_MATRIX=1 retains optional CSR/permutation dumps from v1.
Unset both audit/dump variables for timing runs; no logging occurs per solve.

Setter failures remain mandatory errors. No fallback, precision change, refinement
reduction or scientific guard change was introduced. The v1 metadata readback
assertion was diagnostic-only and has been removed; the actual Set calls remain
checked, with the same two-refinement/one-thread/hybrid-memory-off values.

Exact pinned0.6 scalar contracts, verified against installed cudss.h and cached
official0.6 types/functions docs:
  HOST_NTHREADS, HYBRID_EXECUTE_MODE, HYBRID_MODE, IR_N_STEPS, ND_NLEVELS: int.
  REORDERING_ALG, SOLVE_ALG, FACTORIZATION_ALG: cudssAlgType_t.
  LU_NNZ: int64_t; MEMORY_ESTIMATES: int64_t[16]; row permutation: int32 indices.
The native int and enum are both4 bytes in this frozen build, enforced at compile
time. No int64 scalar workaround or guessed reinterpretation is attempted.
Get sizeWritten is valid only on SUCCESS. The nonfatal getter path never changes
config, input matrix, RHS, solution or numerical acceptance bounds.

Build: frozen CMake Release -O3 -DNDEBUG, CPU6/nice15, --parallel1 and one thread
under exclusive outputs/pic-profile/testing-core.lock; dependencies unchanged.
CPU tests passed for parser cases and mocked cuDSS diagnostic returns, including
SUCCESS with sizeWritten0, NOT_SUPPORTED, known enum values and failed DataGet.
The mock executable links no cuDSS/CUDA runtime and does not exercise a GPU.
Saved logs/mock JSON prove reporting semantics, not real cuDSS query behavior.

bin/run-candidate targets this prefix. READY.json and candidate-manifest.json
identify the build and protected-source hashes. source-from-v1.patch is the delta.
See followup-plan.md for the bounded hybrid trial and ND1024 decision.
