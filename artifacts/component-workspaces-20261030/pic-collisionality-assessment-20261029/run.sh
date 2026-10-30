#!/bin/sh
set -eu
prefix=/home/treefor/projects/runtimes/pic-collisionality-assessment-20261029
lease=/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core7.lock
frozen=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
python=/home/treefor/projects/Power-Flow/.venv/bin/python
flock "$lease" taskset -c 7 nice -n 15 "$frozen" env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 "$python" "$prefix/assess.py" --out "$1"
flock "$lease" taskset -c 7 nice -n 15 "$frozen" env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 "$python" "$prefix/plot_assessment.py" --assessment-dir "$1"
