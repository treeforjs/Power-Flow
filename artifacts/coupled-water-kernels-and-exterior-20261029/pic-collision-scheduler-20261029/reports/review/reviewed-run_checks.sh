#!/bin/bash
set -euo pipefail
base=/home/treefor/projects/runtimes/pic-collision-scheduler-20261029
report=/home/treefor/sim-workspace/powerflow/outputs/pic-collision-scheduler-20261029
mkdir -p "$base/build" "$report"
exec flock /home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core7.lock \
  taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen \
  /bin/bash -c '
    set -euo pipefail
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    ulimit -v 1048576
    ulimit -t 60
    /usr/bin/g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror "$1/source/test_collision_scheduler.cpp" -o "$1/build/test_collision_scheduler"
    "$1/build/test_collision_scheduler" "$2/checks.json"
  ' -- "$base" "$report"
