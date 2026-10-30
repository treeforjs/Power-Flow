CPU hot-path candidate v2

V2 copies the completed v1 source and applies the independently reviewed minimal DOF-cache invalidation patch, SHA-256 aa423ea8f3f221dcc1060d2cdb1bdcfa2498577d365461aa93509804c125f387. Exactly two files differ from v1: include/raikou/q1_scalar_gather.h and tests/test_q1_scalar_gather.cpp. Both patched files match the independently verified fix source. V1 is unchanged.

Contract: invalidate every Q1ScalarGather that references a space before an unversioned numbering mutation such as FiniteElementSpace::ReorderElementToDofTable(). invalidate() permanently selects the existing MFEM fallback; construct a new gather to cache again. Ordinary versioned mesh/space updates keep their existing automatic fallback. This is explicit invalidation, not automatic detection of unversioned changes. Current PIC runtime sources contain no such numbering operations, so there is no call-site change.

The unchanged Release CMake project builds the same heating_history_pic_gpu target and all existing utilities/tests against frozen MFEM 4.10 and frozen CUDA SDK/runtime libraries. readelf confirms cuDSS and CUDA runtime dynamic dependencies. The GPU-capable executable was built but never run.

The original PIC foundation tests and the new DOF invalidation regression pass. Serial PIC output is byte-identical to the sealed-source baseline: 455,537 bytes, SHA-256 79947a84f5e45e8011b046d584fde84c0c731320df4fc37cb067df17d4bbc226.

Compiled ABI probes compare sizeof, alignof, every member offset and member size against active, v1 and v2 headers. All three match, and the checkpoint structures remain standard-layout and trivially copyable: GapPicOptions 168 bytes, GapPicDiagnostics 200 bytes, GapParticle 80 bytes and GapFoil 32 bytes; alignment is 8 for each. ABI probe outputs and their source are included.

Every build and test command runs with logical CPU 6, nice 15, one compiler job, OMP/BLAS thread counts 1, and an exclusive testing-core.lock lease; lock exclusion is verified. CUDA_VISIBLE_DEVICES is empty for this CPU validation. commands.json records the exact frozen-runtime commands and their exit statuses. No dependency, physics, tolerance, source particle count, checkpoint structure, GPU implementation or production process was changed. No GPU execution, timing claim or deployment is included. Complete checkpoint validation and timing remain separate steps.

Reproduce the CPU build and checks:
flock /home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock taskset -c 6 nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= python3 -u /home/treefor/projects/runtimes/pic-cpu-candidate-v2-20261029/evidence/build_and_check.py

Use handoff.json and READY as the stable review entry point. diff-from-v1.patch contains the two-file change; combined-diff-from-sealed-source.patch contains the complete candidate changes. source-manifest.json and build-manifest.json bind sources and build outputs.
