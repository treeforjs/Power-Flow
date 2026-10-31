#!/bin/sh
# Run one CPU-only command inside the frozen root on a lent core K (0, 1 or 3): flock the core lease, taskset, nice 15.
# Refuses when RETURN-CORES-0-3 exists. Usage: lent-core.sh K COMMAND [ARGS...]
set -eu
K=$1; shift
P=/home/treefor/projects/.workspace/outputs/pic-profile
case "$K" in 0|1|3) ;; *) echo "only cores 0, 1, 3 are lent" >&2; exit 64;; esac
[ -e "$P/RETURN-CORES-0-3" ] && { echo "cores 0-3 have been recalled" >&2; exit 75; }
exec flock "$P/testing-core$K.lock" taskset -c "$K" nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    /home/treefor/projects/runtimes/arch-base-20261028/run-frozen "$@"
