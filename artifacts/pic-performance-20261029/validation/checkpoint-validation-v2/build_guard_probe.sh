#!/bin/bash
set -eu
out=/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/cpu-independent-review/checkpoint-validation-v2
runtime=/home/treefor/projects/runtimes
mf="$runtime/mfem-4.10-pic-20261029"
exec 9> /home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core.lock
flock -x 9
"$runtime/arch-base-20261028/run-frozen" /usr/sbin/c++ -std=c++20 -O2 -I"$runtime/mfem-gpu-pic-20261029/source/include" -isystem "$mf/install-mpi/include" -isystem "$mf/tpl/hypre/include" -isystem "$mf/tpl/metis/include" "$out/residual_guard_probe.cpp" -o "$out/residual_guard_probe" "$mf/install-mpi/lib/libmfem.a" "$mf/tpl/hypre/lib/libHYPRE.a" "$mf/tpl/metis/lib/libmetis.a" "$mf/tpl/GKlib/lib/libGKlib.a" /usr/lib/libmpi.so -Wl,-rpath,/usr/lib > "$out/guard-build.log" 2>&1
