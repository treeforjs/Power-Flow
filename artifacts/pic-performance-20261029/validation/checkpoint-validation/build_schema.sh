#!/bin/bash
set -eu
out=/home/treefor/sim-workspace/powerflow/outputs/pic-profile/cpu-independent-review/checkpoint-validation
runtime=/home/treefor/projects/runtimes
exec 9> /home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock
flock -x 9
"$runtime/arch-base-20261028/run-frozen" /usr/sbin/c++ -std=c++20 -O2 -ffunction-sections -fdata-sections -Wl,--gc-sections -I"$runtime/mfem-gpu-pic-20261029/source/include" -isystem "$runtime/mfem-4.10-pic-20261029/install-mpi/include" -isystem "$runtime/mfem-4.10-pic-20261029/tpl/hypre/include" -isystem "$runtime/mfem-4.10-pic-20261029/tpl/metis/include" "$out/schema_probe.cpp" -o "$out/schema_probe" > "$out/schema-build.log" 2>&1
"$runtime/arch-base-20261028/run-frozen" "$out/schema_probe" > "$out/abi-schema.json"
