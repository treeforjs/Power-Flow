#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
build_dir="${1:?explicit isolated build directory required}"
exec flock -w 30 /home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core7.lock \
 taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen \
 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= \
 bash -c 'set -euo pipefail
 cmake -S native/mfem_thermoelastic_v3 -B "$1" -DCMAKE_BUILD_TYPE=Release -DMFEM_DIR=/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/cmake/mfem
 cmake --build "$1" -j1' bash "$build_dir"
