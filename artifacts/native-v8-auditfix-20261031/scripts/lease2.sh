#!/bin/bash
# Two single-threaded CPUs for a build: any two free leases among Power-Flow cores 6 and 7
# (testing-core.lock, testing-core7.lock), their SMT siblings 15 and 14 (testing-core15.lock,
# testing-core14.lock; builds only) and the lent cores 0 and 1 (testing-core<K>.lock, only while
# RETURN-CORES-0-3 does not exist; never cores 2-5). Non-blocking attempts, released and retried when
# fewer than two are free, so waiting jobs cannot deadlock. nice 15, one thread per process.
# usage: scripts/lease2.sh COMMAND [ARGS...]
set -u
# PF_CPU15_ONLY=1 (core loan): route to scripts/lease15.sh (logical CPU 15, nice 19, testing-cpu15.lock).
[ "${PF_CPU15_ONLY:-0}" = 1 ] && exec "$(dirname "$0")/lease15.sh" "$@"
P=/home/treefor/projects/.workspace/outputs/pic-profile
declare -A LOCK=([6]=$P/testing-core.lock [7]=$P/testing-core7.lock [15]=$P/testing-core15.lock [14]=$P/testing-core14.lock [0]=$P/testing-core0.lock [1]=$P/testing-core1.lock)
while :; do
  got=(); fds=()
  for c in 7 6 15 14 0 1; do
    if [ "$c" -le 1 ]; then [ -e "$P/RETURN-CORES-0-3" ] && continue; fi
    exec {fd}>>"${LOCK[$c]}"
    if flock -n "$fd"; then got+=("$c"); fds+=("$fd"); else exec {fd}>&-; fi
    [ "${#got[@]}" -eq 2 ] && break
  done
  if [ "${#got[@]}" -eq 2 ]; then
    cpus="${got[0]},${got[1]}"
    PF_CORE=$cpus taskset -c "$cpus" nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "$@"
    rc=$?; for f in "${fds[@]}"; do exec {f}>&-; done; exit "$rc"
  fi
  for f in "${fds[@]}"; do exec {f}>&-; done
  sleep 3
done
