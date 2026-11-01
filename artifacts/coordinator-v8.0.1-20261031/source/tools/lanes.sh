#!/bin/bash
# Run a list of run_v8.py commands on CPU lanes: lent core 1 (flock + taskset + nice 15; skipped once
# RETURN-CORES-0-3 exists) and the Power-Flow cores 6/7 (core-lease any-core.sh --wait). One job per lane at a time.
# Usage: lanes.sh LABEL JOBFILE   (JOBFILE: one line per job: NAME MODE ARGS...)
set -u
LABEL=$1; JOBS=$2
D=/home/treefor/projects/runtimes/pic-water-common-commit-v8.0.1-20261031
L=/home/treefor/projects/.workspace/outputs/pic-profile
RF=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
O=/home/treefor/projects/.workspace/outputs/pic-water-common-commit-v8.0.1-20261031/dev
ENVS="env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1"
mkdir -p $O
cd $D
mapfile -t LINES < <(grep -v '^\s*#' "$JOBS" | grep -v '^\s*$')
next=0
lock=$O/lanes-$$.lock; : > $lock
take() { flock $lock bash -c "n=\$(cat $O/lanes-$$.next 2>/dev/null || echo 0); echo \$((n+1)) > $O/lanes-$$.next; echo \$n"; }
echo 0 > $O/lanes-$$.next
lane() {
  local kind=$1
  while :; do
    i=$(take)
    [ "$i" -ge "${#LINES[@]}" ] && return
    set -- ${LINES[$i]}
    name=$1; mode=$2; shift 2
    if [ "$kind" = lent ]; then
      if [ -e $L/RETURN-CORES-0-3 ]; then
        /home/treefor/projects/runtimes/core-lease-20261030/any-core.sh --wait $RF /usr/bin/python -I run_v8.py run $LABEL $mode $name "$@" > $O/$name.out 2>&1
      else
        flock $L/testing-core1.lock taskset -c 1 nice -n 15 $ENVS $RF /usr/bin/python -I run_v8.py run $LABEL $mode $name "$@" > $O/$name.out 2>&1
      fi
    else
      /home/treefor/projects/runtimes/core-lease-20261030/any-core.sh --wait $RF /usr/bin/python -I run_v8.py run $LABEL $mode $name "$@" > $O/$name.out 2>&1
    fi
    echo "$(date +%T) $name rc=$? $(tail -n 1 $O/$name.out)"
  done
}
lane lent & lane pf & lane pf &
wait
rm -f $O/lanes-$$.next $lock
