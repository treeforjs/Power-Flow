#!/bin/sh
# CPU sampling cost per collision, v3.1 against v3.2, on identical inputs and one core (run under tools/lease_any.sh).
# usage: tools/bench_cost.sh OUTDIR   (OUTDIR holds pop-1000.bin, pop-4711.bin, pop-99.bin from bench_gen)
# Cases as tests/gpu_parity.cpp: (pop seed, cells, seed, step, dt). Each line: configuration, case, JSON of bench_step.
# The node's cores are shared (SMT siblings busy): ROUNDS rounds (default 4) interleave the configurations, each run
# taking the minimum of 5 repetitions; the summary keeps the minimum per configuration and case over the rounds.
set -eu
R=/home/treefor/projects/runtimes/pic-water-chemistry-v3.2.2-auditprep-20261101
V31=/home/treefor/projects/runtimes/pic-water-chemistry-v3.1-20261031/data/tables/chemistry-v3.1.blob
V32=$R/data/tables/chemistry-v3.2.blob
FR=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
O=$1
ROUNDS=${2:-4}
r=0
while [ "$r" -lt "$ROUNDS" ]; do
r=$((r + 1))
for cfg in "v31code-v31tables $R/build/bin/bench_v31code $V31" "v32code-v31tables $R/build/bin/bench_v32 $V31" "v32code-v32tables $R/build/bin/bench_v32 $V32"; do
  set -- $cfg
  name=$1; bin=$2; blob=$3
  for c in "1000 6840227781351331448 4242 4e-11" "4711 12379813738877118345 77 5e-12" "99 3 549755813893 2e-10"; do
    set -- $c
    printf '%s %s ' "$name" "$1"
    $FR "$bin" run "$blob" "$O/pop-$1.bin" "$2" "$3" "$4" 5
  done
done
done
