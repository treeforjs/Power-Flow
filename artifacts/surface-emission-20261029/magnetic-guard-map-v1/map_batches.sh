#!/bin/sh
set -eu
base=/home/treefor/projects/runtimes/pic-magnetized-emission-20261029
out=/home/treefor/sim-workspace/powerflow/outputs/pic-magnetized-emission-20261029
probes=/home/treefor/sim-workspace/powerflow/outputs/cathode-emission-audit-20261029/fields/probes.csv
lock=/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core7.lock
frozen=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
flock "$lock" taskset -c 7 nice -n 15 "$frozen" cmake --build "$base/build" --parallel 1
first=0
while [ "$first" -lt 1056 ]; do
  count=64
  if [ $((first+count)) -gt 1056 ]; then count=$((1056-first)); fi
  flock "$lock" taskset -c 7 nice -n 15 "$frozen" "$base/build/guard_assessment" map "$probes" "$out/rate-map" "$first" "$count"
  first=$((first+count))
done
