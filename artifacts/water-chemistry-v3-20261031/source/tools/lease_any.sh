#!/bin/bash
# Run one single-threaded CPU job on the first free core among 6, 7 (Power-Flow) and 0, 1, 3 (lent), holding that
# core's lease file for the whole job (taskset, nice 15, one thread). Cores 0-3 are skipped once
# RETURN-CORES-0-3 exists. Cores 2, 4, 5 are never used.  usage: lease_any.sh COMMAND [ARGS...]
set -u
P=/home/treefor/projects/.workspace/outputs/pic-profile
declare -A LOCK=([6]=$P/testing-core.lock [7]=$P/testing-core7.lock [0]=$P/testing-core0.lock [1]=$P/testing-core1.lock [3]=$P/testing-core3.lock)
while :; do
  for c in 6 7 0 1 3; do
    if [ "$c" -lt 4 ] && [ -e "$P/RETURN-CORES-0-3" ]; then continue; fi
    exec {fd}>>"${LOCK[$c]}"
    if flock -n "$fd"; then
      PF_CORE=$c taskset -c "$c" nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "$@"
      rc=$?; exec {fd}>&-; exit "$rc"
    fi
    exec {fd}>&-
  done
  sleep 2
done
