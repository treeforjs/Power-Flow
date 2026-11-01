#!/bin/bash
# Load-fixture host profiles of the coordinator v7.1 copy with the Hall v3.4 resident pair (the v7.1 prof6 settings:
# 8 / 4 / 2 steps at 1e4 / 1e5 / 1e6 carriers on 200 and 50 um, GPU push, host-work v7.1). One run per GPU lock (both
# lock files, then a core lease). Usage: run_profile.sh BUILD LOG
set -u
B=${1:?build}; LOG=${2:?log}
cd /home/treefor/projects/runtimes/hall-corrective-v3.4-20261031/coordinator-v7.1-copy
L=/home/treefor/projects/.workspace/outputs/pic-profile
for N in 10000 100000 1000000; do
  case $N in 10000) S=8;; 100000) S=4;; *) S=2;; esac
  for M in 200 50; do
    MESH=(); [ $M = 50 ] && MESH=(--mesh /home/treefor/projects/runtimes/resolution-series-prep-20261031/meshes/mesh-50um.json)
    echo "$(date +%H:%M:%S) START $N $M" >> $LOG
    flock $L/testing-gpu.lock flock $L/testing-gpu-GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9.lock /home/treefor/projects/runtimes/core-lease-20261030/any-core.sh --wait \
      /home/treefor/projects/runtimes/arch-base-20261028/run-frozen /usr/bin/python -I run_v7_1.py run $B run $B-prof-n$N-${M}um \
      --fixture load --load-carriers $N --steps $S --progress-every 1 --host-work v7.1 "${MESH[@]}" --field gpu --push gpu >> $LOG 2>&1
    echo "$(date +%H:%M:%S) END rc=$? $N $M" >> $LOG
  done
done
echo "$(date +%H:%M:%S) DONE" >> $LOG
