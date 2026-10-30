#!/bin/bash
set -euo pipefail
exec flock /home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core7.lock taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= bash -c '
set -e; ulimit -v 1048576; ulimit -t 60
p=/home/treefor/projects/runtimes/pic-collision-exposure-clock-review-20261029
o=/home/treefor/sim-workspace/powerflow/outputs/pic-collision-exposure-clock-review-20261029
/usr/bin/time -v g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror "$p/source/probe.cpp" -o "$p/build/probe" 2> "$o/build.log"
/usr/bin/time -v "$p/build/probe" "$o/probes.json" > "$o/probes.log" 2> "$o/run-resource.log"
'
