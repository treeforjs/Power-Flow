#!/bin/sh
# Native v8 CPU controls (host mode), each on one leased core (scripts/lease.sh), two lanes in parallel.
# usage: scripts/cpu_suite.sh OUTDIR
O=$1; mkdir -p $O
D=/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-full-time-candidate008-auditfix
F=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
B=$D/build
run() { name=$1; shift; mkdir -p $O/$name; $D/scripts/lease.sh env CUDA_VISIBLE_DEVICES= $F prlimit --as=4294967296 -- "$@" > $O/$name.log 2>&1; echo "$name exit $? $(tail -n 1 $O/$name.log | cut -c1-160)" >> $O/progress.txt; }
lane1() {
  run containers $B/test_v8_containers
  run chem-v3 $B/test_v8_chemistry $O/chem-v3 full
  PF_V8_CHEM=v3.1 run chem-v3.1 $B/test_v8_chemistry $O/chem-v3.1 full
  run sinks $B/test_v8_sinks $O/sinks
}
lane2() {
  for K in 2 3 5; do for m in unbroken save restore tamper; do
    [ $m = tamper ] && [ $K != 3 ] && continue
    run restore-k$K-$m $B/test_v8_restore $O/restore-k$K $m $K
  done; done
}
lane1 & lane2 & wait
echo "done" >> $O/progress.txt
