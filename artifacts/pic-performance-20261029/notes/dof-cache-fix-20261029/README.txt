Isolated DOF cache invalidation fix

Status: prepared and CPU-verified. Shared candidate and active GPU run untouched.

Source copy:
/home/treefor/sim-workspace/powerflow/outputs/pic-profile/cpu-independent-review/dof-cache-fix-20261029/source
Base source:
/home/treefor/projects/runtimes/pic-cpu-candidate-20261029/source
Patch:
/home/treefor/sim-workspace/powerflow/outputs/pic-profile/cpu-independent-review/dof-cache-fix-20261029/dof-cache-invalidation.patch
Patch SHA-256: aa423ea8f3f221dcc1060d2cdb1bdcfa2498577d365461aa93509804c125f387

Minimal production change
Q1ScalarGather now exposes invalidate(), which clears the cached connectivity and permanently selects the existing MFEM fallback. The header requires invalidating before operations that change numbering without updating a sequence counter, such as FiniteElementSpace::ReorderElementToDofTable(). Construct a fresh Q1ScalarGather to enable caching again. Ordinary versioned mesh/space updates retain their existing automatic fallback. The value/gradient fast paths and all physics code remain unchanged.

Usage for an unversioned connectivity operation:
    gather.invalidate();
    space.ReorderElementToDofTable();
    // gather.value()/gradient() now use MFEM directly.
    Q1ScalarGather rebuilt(field); // optional: cache the current numbering.

This is an explicit invalidation contract, not automatic detection of unversioned MFEM mutations. No existing PIC call site performs such a mutation, so no call-site change is required. All gathers referencing a space must be invalidated before an unversioned numbering change. No shared candidate change or deployment was performed.

Only two source files differ from the shared candidate:
- include/raikou/q1_scalar_gather.h: documentation and invalidate().
- tests/test_q1_scalar_gather.cpp: regression for unchanged mesh/space sequence counters, repeated invalidation, fallback values/gradients and reconstructed cache values/gradients on all six elements at corners/interiors.

CPU validation
The original reproducer was rerun on CPU 6 at nice 15 and confirmed the original wrong values on 6/6 elements. The isolated fix then passed:
- Existing foundation tests with serial CG.
- Gather tests, including the new renumbering regression, Q1 geometry cases, Q2 and triangle fallback, and refinement invalidation.
- Independent edge cases: changing coefficients/geometry, signed DOFs, nonidentity DOF transforms, Q2 rebinding, invalidation and cache reconstruction. The renumbered values/gradients match MFEM bit for bit; zero mismatches.
- Supplied serial trace: 455537 bytes, byte-identical to sealed-baseline evidence; SHA-256 79947a84f5e45e8011b046d584fde84c0c731320df4fc37cb067df17d4bbc226.
- Complete serial state trace, including every diagnostic double and Poisson residual: 328651 bytes, byte-identical; SHA-256 1624f743214d33a3986e8ca5b7d6271a6a012bc105bcff7db9fe9f6a69955cd3.
- Native field/time-cache trace: 31386 bytes, byte-identical; SHA-256 64fe3cfda218d0b1fea64b578187145544532abdba6db160ae2a646f21a794c9. The two-million-sample checksum also matches.

All builds, utilities and tests used CPU 6, nice 15 and one local job; OMP/BLAS thread counts were 1 and CUDA_VISIBLE_DEVICES was empty for tests. The initial build/full checks completed before the lease update was observed. A final six-check verification was then rerun beneath exclusive flock on:
/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock
The final verification confirmed lock exclusion, CPU affinity [6] and nice 15. See leased-verification.json and leased-logs/. The rerun entry point now acquires the same lease before any build or test.

Reproduce all CPU checks (the script takes the lease internally):
taskset -c 6 nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= python3 -u /home/treefor/sim-workspace/powerflow/outputs/pic-profile/cpu-independent-review/dof-cache-fix-20261029/run_cpu_checks.py

Evidence
results.json: complete first-run result and trace comparisons.
leased-verification.json: final checks under the testing-core lease.
commands.json: exact frozen-runtime configure/build/compiler/test commands.
shared-source-before.json and source-changes.json: complete input hashes and the two changed source hashes.
logs/ and leased-logs/: raw compiler/test outputs.

Limits
Unversioned numbering mutations require the explicit invalidate() call. No automatic mutation tracking, dependency changes, physics/guard changes, GPU execution, combined checkpoint replay, end-to-end performance claim, reservation changes or deployment are included. Dependencies/compiler came from the frozen runtime. The shared source was hash-checked after the checks and remains unchanged.
