#!/bin/bash
set -euo pipefail
p=/home/treefor/sim-workspace/powerflow/outputs/pic-water-event-bridge-20261029/review
exec flock /home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core7.lock \
 taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen /bin/bash -c '
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=
ulimit -v 1048576
ulimit -t 60
p=$1
mkdir -p "$p/build"
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -I"$p" -I"$p/inputs" -I"$p/inputs/include" "$p/probe.cpp" -o "$p/build/probe"
read -r kernel_sha _ < <(sha256sum "$p/inputs/include/raikou/h2o_conditional_parent.h")
read -r data_sha _ < <(sha256sum "$p/inputs/h2o-gas.json")
"$p/build/probe" "$p/inputs/h2o-gas.json" "$data_sha" "$kernel_sha" "$p/probes.jsonl"
python3 "$p/oracle.py"
' -- "$p"
