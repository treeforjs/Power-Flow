#!/bin/bash
# One single-threaded job on a free core lease: Power-Flow cores 6 and 7 (testing-core.lock,
# testing-core7.lock), then the lent cores 0 and 1 (testing-core<K>.lock) while RETURN-CORES-0-3 does
# not exist (never cores 2-5). Waits for a free core. nice 15, one thread.
# PF_LEASE_SMT=1 (non-timing work only: comparisons, test compiles) also tries the SMT siblings 15 and 14
# of cores 7 and 6 (testing-core15.lock, testing-core14.lock) after 6 and 7.
# usage: scripts/lease.sh COMMAND [ARGS...]
set -u
# PF_CPU15_ONLY=1 (core loan): route to scripts/lease15.sh (logical CPU 15, nice 19, testing-cpu15.lock).
[ "${PF_CPU15_ONLY:-0}" = 1 ] && exec "$(dirname "$0")/lease15.sh" "$@"
P=/home/treefor/projects/.workspace/outputs/pic-profile
declare -A LOCK=([6]=$P/testing-core.lock [7]=$P/testing-core7.lock [15]=$P/testing-core15.lock [14]=$P/testing-core14.lock [0]=$P/testing-core0.lock [1]=$P/testing-core1.lock)
CPUS="6 7 0 1"; [ "${PF_LEASE_SMT:-0}" = 1 ] && CPUS="6 7 15 14 0 1"
while :; do
  for c in $CPUS; do
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
