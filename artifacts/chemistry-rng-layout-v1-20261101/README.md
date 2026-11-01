# Disjoint chemistry random-word layout

This isolated successor fixes an inherited overlap: the four target-Maxwellian rejection attempts used words 4–19, while ordinary collision kinematics started at word 8. Rejected target samples could therefore condition the same random words later used for scattering angles. The sealed v3.1.1 runtime is unchanged.

Layout version 2 preserves candidate time, target-species and channel draws. It reserves half-open ranges [4,20) for the four target attempts, [32,96) for ordinary kinematics, and [96,352) for a future bounded H2 partition sampler, inside the existing 512-word event domain. Compile-time interval checks and pre-draw ordinary-kinematics bounds prevent overlap and overflow. The H2 sampler is not implemented by this artifact. Blob and species/process ABI are unchanged; numerical kernel identity changes. Old fixed-seed trajectories are intentionally not claimed identical.

CPU SimQueue job 503 used one allocated core, frozen userspace and 2 GiB. An actual manufactured MCC fixture with a one-sigma target truncation reproduced the old defect: among 1,678 events accepted on the second target attempt, the mean reconstructed CM scattering cosine was 0.67517 and correlation with the reused uniform draw was effectively one. With the corrected independent range the corresponding mean was 0.000744 and correlation -0.01873. The first-event target identity was identical before and after the correction. This strongly exercised fixture does not quantify the bias frequency under the production seven-sigma truncation.

The CPU checks cover 10,297 populated first events, 277 direct interval and boundary controls, 32 actual late-rejection deterministic replays, energy/momentum closure, and all 11 inherited chemistry audit regressions. Corrected normal and UBSan reports were identical. Maximum relative event energy defect was 4.41e-16.

GPU SimQueue job 506 used one allocated core and a shared 3 GiB GPU reservation. Seven cases across three launch layouts had zero CPU/GPU carrier, event-record or ledger mismatches. The added rejection fixture contained 1,258 events, including 500 actual late target rejections. The math fixture also matched 4,194,304 values. Timings in raw receipts are diagnostic execution records, not performance benchmarks. CPU and GPU evidence have separate input manifests; the later GPU addendum did not change the tested CPU headers.

A read-only peer review found no concrete defect in the two changed headers. This qualifies the scoped implementation checks, not the physical chemistry network. Native adoption requires a new explicit kernel/configuration/checkpoint identity and separately qualified actual native fixtures; no existing native binding or production run was changed.

Reproduce only through SimQueue using the exact runtime and dependencies in INPUTS.json and GPU-INPUTS.json:

```sh
simq submit -p powerflow -c 1 --mem 2G -t 00:05:00 -n chemistry-rng-layout-cpu -e CUDA_VISIBLE_DEVICES= -- /home/treefor/projects/runtimes/arch-base-20261028/run-frozen /usr/bin/python3 -I /home/treefor/projects/runtimes/chemistry-rng-layout-v1-20261101/tools/check.py
simq submit -p powerflow -c 1 --mem 3G --gpus 1 --gpu-mem 3G -t 00:05:00 -n chemistry-rng-layout-gpu -- /home/treefor/projects/runtimes/arch-base-20261028/run-frozen /usr/bin/python3 -I /home/treefor/projects/runtimes/chemistry-rng-layout-v1-20261101/tools/check_gpu.py
```

The shared source is derived from the project's sealed chemistry implementation; no primary raw cross-section arrays, table blob or binary is included here. Exact dependencies remain local and pinned. `production_qualified=false`, `physical_model_validated=false`, `network_complete=false`, `native_adopted=false`.
