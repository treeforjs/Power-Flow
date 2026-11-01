#!/bin/bash
# Run one single-threaded CPU job on the first free leased core, holding that core's lease file for the whole job
# (taskset, nice, one thread). Cores and lease files: PFC_LEASE (default "6 7", the Power-Flow cores); during a core
# loan PFC_LEASE=15 selects logical CPU 15 with testing-cpu15.lock at nice 19 (shared by the Power-Flow agents);
# PFC_LEASE=none inside a SimQueue job (the queue holds the cores).
# Cores 0-3 are skipped once RETURN-CORES-0-3 exists; cores 2-5 are never used.  usage: lease_any.sh COMMAND [ARGS...]
set -u
P=/home/treefor/projects/.workspace/outputs/pic-profile
declare -A LOCK=([6]=$P/testing-core.lock [7]=$P/testing-core7.lock [0]=$P/testing-core0.lock [1]=$P/testing-core1.lock [15]=$P/testing-cpu15.lock)
declare -A NICE=([6]=15 [7]=15 [0]=15 [1]=15 [15]=19)
CORES=${PFC_LEASE:-6 7}
if [ "$CORES" = none ]; then   # inside a SimQueue job: the queue holds the cores; run on the job's CPUs as given
  exec nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "$@"
fi
while :; do
  for c in $CORES; do
    case $c in 2|3|4|5) continue ;; esac
    if [ "$c" -lt 4 ] && [ -e "$P/RETURN-CORES-0-3" ]; then continue; fi
    exec {fd}>>"${LOCK[$c]}"
    if flock -n "$fd"; then
      PF_CORE=$c taskset -c "$c" nice -n "${NICE[$c]}" env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "$@"
      rc=$?; exec {fd}>&-; exit "$rc"
    fi
    exec {fd}>&-
  done
  sleep 2
done
