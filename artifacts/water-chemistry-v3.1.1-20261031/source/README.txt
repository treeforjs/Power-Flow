Water chemistry v3.1.1 audit correction
======================================
Isolated successor to sealed v3.1. Active v3.2 is unchanged.
Read INTEGRATION.txt before linking. Tests and evidence are in tests/ and evidence/.
No data tables derived from LXCat may be placed in the repository artifact.
The runtime retains locally frozen source data and fixed runtime/compiler paths.

Submit tools/audit_cpu.sh, tools/audit_gpu_build.sh, tools/audit_gpu.sh through
SimQueue, always with explicit --mem. Use one allocated core. GPU verification
is short and shared (--gpus 1 --gpu-mem 3G); timing output is not a benchmark.
No raw GPU/core locks or affinity overrides are used by these audit scripts.
Standalone checks passed; see SEAL.json and RESULTS.txt. Full coupled physics remains unqualified.
Historical SWARM-VALIDATION.txt describes v3.1 only and is not a new validation.
