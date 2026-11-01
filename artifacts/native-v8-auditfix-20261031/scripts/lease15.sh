#!/bin/bash
# Core loan mode (Power-Flow cores 6-7 lent to Spheral): one job on logical CPU 15 only, nice 19,
# holding the shared testing-cpu15.lock (all Power-Flow agents share it). Waits for the lock.
# usage: scripts/lease15.sh COMMAND [ARGS...]
set -u
P=/home/treefor/projects/.workspace/outputs/pic-profile
exec {fd}>>"$P/testing-cpu15.lock"
flock "$fd"
PF_CORE=15 taskset -c 15 nice -n 19 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "$@"
rc=$?; exec {fd}>&-; exit "$rc"
