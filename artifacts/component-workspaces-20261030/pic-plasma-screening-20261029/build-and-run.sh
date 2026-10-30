#!/bin/bash
set -euo pipefail
out=/home/treefor/sim-workspace/powerflow/outputs/pic-plasma-screening-20261029
root=/home/treefor/projects/runtimes/pic-plasma-screening-20261029
evidence=${1:-"$out/evidence-review"}
if [ -e "$evidence" ]; then printf '%s\n' 'Use a fresh evidence directory.' >&2; exit 2; fi
exec 9>/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core7.lock
flock -x 9
cd "$root"
taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= cmake -S "$root" -B "$root/build" -DMFEM_DIR=/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/cmake/mfem -DCMAKE_BUILD_TYPE=Release >"$out/configure.log" 2>&1
taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= cmake --build "$root/build" -j1 >"$out/build.log" 2>&1
taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= "$root/build/plasma_screening" "$evidence" >"$out/run.log" 2>&1
