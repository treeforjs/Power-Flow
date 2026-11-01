#!/bin/sh
set -eu
base=/home/treefor/projects/runtimes/pic-magnetized-emission-20261029
out=/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-magnetized-emission-20261029/tangential-onepass
lock=/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core7.lock
frozen=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
flock "$lock" taskset -c 7 nice -n 15 "$frozen" /bin/sh -c 'cmake -S "$1" -B "$1/build" -DCMAKE_BUILD_TYPE=Release && cmake --build "$1/build" --parallel 1' shell "$base"
for T in 298.2 300 773.15 1490; do
 for B in 0 30; do
  for refinement in 1 2; do
   flock "$lock" taskset -c 7 nice -n 15 "$frozen" "$base/build/joint_qualification" "$out" "$T" "$B" "$refinement"
  done
 done
done
flock "$lock" taskset -c 7 nice -n 15 "$frozen" "$base/build/joint_admission" "$out"
