#!/bin/bash
# One single-threaded job on a free core lease: Power-Flow cores 6 and 7 (testing-core.lock,
# testing-core7.lock), then the lent cores 0 and 1 (testing-core<K>.lock) while RETURN-CORES-0-3 does
# not exist (never cores 2-5). Waits for a free core. nice 15, one thread.
# usage: scripts/lease.sh COMMAND [ARGS...]
set -u
P=/home/treefor/projects/.workspace/outputs/pic-profile
declare -A LOCK=([6]=$P/testing-core.lock [7]=$P/testing-core7.lock [0]=$P/testing-core0.lock [1]=$P/testing-core1.lock)
while :; do
  for c in 6 7 0 1; do
    if [ "$c" -le 1 ]; then [ -e "$P/RETURN-CORES-0-3" ] && continue; fi
    exec {fd}>>"${LOCK[$c]}"
    if flock -n "$fd"; then
      PF_CORE=$c taskset -c "$c" nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "$@"
      rc=$?; exec {fd}>&-; exit "$rc"
    fi
    exec {fd}>&-
  done
  sleep 2
done
