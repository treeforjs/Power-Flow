#!/usr/bin/env bash
set -eu
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
runtime=/home/treefor/projects/runtimes/pic-source-sampling-validation-20261029/coupled-validation
output=/home/treefor/sim-workspace/powerflow/outputs/pic-profile/source-sampling-validation/coupled-preflight
exec 9>>/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock
if ! taskset -c 6 nice -n 15 flock -n 9; then
  printf '%s\n' 'Testing CPU6 lease busy; no build or preflight started.' >&2
  exit 75
fi
taskset -c 6 nice -n 15 timeout 180s python3 "$runtime/preflight.py" --build --output "$output"
