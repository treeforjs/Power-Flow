#!/bin/sh
set -eu
prefix=/home/treefor/projects/runtimes/pic-collision-grid-design-20261029
lease=/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core7.lock
frozen=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
mkdir "$1"
flock "$lease" taskset -c 7 nice -n 15 "$frozen" env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python "$prefix/microcase.py" --out "$1/geometry"
flock "$lease" taskset -c 7 nice -n 15 "$frozen" env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python "$prefix/prepare_sparse_case.py" --out "$1/inputs"
flock "$lease" taskset -c 7 nice -n 15 "$frozen" env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python "$prefix/run_sparse_test.py" --input "$1/inputs" --out "$1/parity"
