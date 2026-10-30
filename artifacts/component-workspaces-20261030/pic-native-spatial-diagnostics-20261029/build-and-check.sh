#!/bin/sh
set -eu
prefix=/home/treefor/projects/runtimes/pic-native-spatial-diagnostics-20261029
lease=<workspace>/outputs/pic-profile/testing-core7.lock
frozen=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
mfem=/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/cmake/mfem
# Fresh output directory is required. No simulation or GPU operation is present.
flock "$lease" taskset -c 7 nice -n 15 "$frozen" env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 cmake -S "$prefix" -B "$prefix/build" -G Ninja -DMFEM_DIR="$mfem"
flock "$lease" taskset -c 7 nice -n 15 "$frozen" env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 cmake --build "$prefix/build" --parallel 1
flock "$lease" taskset -c 7 nice -n 15 "$frozen" env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python "$prefix/tools/run_contracts.py" --out-dir "$1"
