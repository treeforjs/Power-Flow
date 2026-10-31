#!/bin/bash
# Run one single-threaded CPU-only job on a lent core (0, 1 or 3) holding its lease: lent.sh K COMMAND [ARGS...]
# Refuses when pic-profile/RETURN-CORES-0-3 exists. Never used for GPU work.
set -u
P=/home/treefor/projects/.workspace/outputs/pic-profile
K=$1; shift
case "$K" in 0|1|3) ;; *) echo "core $K is not lent" >&2; exit 64;; esac
[ -e "$P/RETURN-CORES-0-3" ] && { echo "lent cores recalled" >&2; exit 75; }
exec flock "$P/testing-core$K.lock" taskset -c "$K" nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PF_CORE=$K "$@"
