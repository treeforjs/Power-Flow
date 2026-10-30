Comparator v2: explicit replay CPU policy, unchanged scientific checks

This is a new isolated copy. The original checkpoint-validation tool, original evidence, CPU candidate and reviewed DOF patch remain untouched. No production operations or GPU execution are performed. Comparator execution is still required on CPU6 at nice15, one thread/job, beneath testing-core.lock.

New CLI:
    --expected-replay-cpu 6|7       default 6
    --allocation-receipt PATH      default outputs/pic-testing-allocation.json

Both baseline and candidate run receipts must match the selected expected replay CPU. This option describes saved replay execution; it never changes comparator CPU affinity. CPU8 is rejected by the API and CLI regardless of receipt contents.

Default CPU6 behavior retains the v1 receipt rule and needs no new allocation proof. Explicit CPU7 requires:
- Recorded status paused_job_frozen_core7_released_for_testing.
- native_stopped=true, production_cgroup_frozen=true and gpu_benchmarks_serial=true.
- Authorized testing CPUs containing6 and7, with CPU8 excluded.
- A recorded production identity and finite allocation timestamp.
- Each replay receipt recording profile_cpu=7, production_left_stopped=true, a stopped native state T/t, and production identity keys matching the allocation (pid, executable, start_ticks and arguments where recorded).
- Exactly one taskset -c7 or taskset --cpu-list7 command record, represented as the ordinary separate tokens taskset,-c,7 or taskset,--cpu-list,7. The expected CPU argument, receipt label and recorded command must agree.

The allocation SHA-256 and proof fields are included in the comparison. Allocation replacement during the comparison is rejected. This is recorded allocation/receipt validation, not a live-process/cgroup inspection or a new authorization to run a GPU workload. No metadata or numerical input is rewritten.

Example for the completed CPU7 phase pair:
taskset -c 6 nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= python3 /home/treefor/sim-workspace/powerflow/outputs/pic-profile/cpu-independent-review/checkpoint-validation-v2/compare_saved_pic.py /home/treefor/sim-workspace/powerflow/outputs/pic-profile/phase-baseline-late-cpu7 /home/treefor/sim-workspace/powerflow/outputs/pic-profile/phase-cpu-v2-late-cpu7 /home/treefor/sim-workspace/powerflow/outputs/pic-profile/seed-late-phase --expected-replay-cpu 7 --output /home/treefor/sim-workspace/powerflow/outputs/pic-profile/cpu-independent-review/checkpoint-validation-v2/another-cpu7-comparison.json

Do not add an outer flock: the CLI acquires the CPU6 lease itself. Always check the process exit status and report policy fields; an invalid CLI invocation does not make an older report a new successful validation.

Scientific behavior preserved
POLICY and all finite/ABI/state/kinetic/norm/tolerance/RHS/density reconstruction function ASTs are unchanged. The ABI schema, CPU-only RHS probe source/binary, and all numerical constants/production residual guard are unchanged. Only replay receipt-policy plumbing and proof metadata are added. V1-README.txt describes the full scientific coverage and tolerances. The existing corruption harness is included unchanged; original v1 146-fixture evidence stays at its original location.

Evidence in this directory
- science-preservation.json: unchanged scientific function ASTs, numerical policy, probe and ABI hashes.
- v1-preservation-before.json: original tool/evidence hashes, rechecked after v2 tests.
- policy-tests.json: short CPU6/7/8 and allocation/identity/affinity proof fixtures, real default6 comparison and actual CPU7 phase-pair CLI test.
- comparison-default6.json: original CPU6 pair, with scientific result sections exactly matching v1.
- comparison-phase-cpu7.json: completed CPU7 phase baseline/candidate compared with explicit selection and recorded allocation proof, while this comparator runs onCPU6.
- v2-policy.patch: complete comparator-source difference from v1.
- cpu7-cli.log and cpu8-cli.log: CLI pass and CPU8 rejection.

Exit codes are unchanged:0 pass,1 finite numerical comparison failure,2 invalid/unsupported/proof/guard failure. No acceptance tolerances are widened. The original limitations remain: no nodal potential/error vector/crossing logs in saved files, unmatched intermediate history rows are not inferred, and short replay agreement is not physical/grid/time convergence.
