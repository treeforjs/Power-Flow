#!/bin/bash
# Hall v3.4 host/device timing matrix (test_v34 bench): modes owner (v3.3 Owner, host arithmetic, GPU push and gpu1 solve),
# v33 (resident owner, v3.3 host path), device (v3.4 device path with host record arrays for native v7) and lane (v3.4
# device path without host record arrays: native v7.1 device lane); 1e4/1e5/1e6 carriers on 200 um (45) and 50 um (180).
# One run per GPU lock (both lock files, then a core lease); 1e6 runs one at a time. Usage: run_bench.sh LABEL LOG
set -u
LABEL=${1:?label}; LOG=${2:?log}
D=/home/treefor/projects/runtimes/hall-corrective-v3.4-20261031
L=/home/treefor/projects/.workspace/outputs/pic-profile
cd $D
for N in 10000 100000 1000000; do
  case $N in 10000) S=6;; 100000) S=4;; *) S=3;; esac
  for M in 45 180; do
    for MODE in ${MODES:-owner v33 device lane}; do
      echo "$(date +%H:%M:%S) START $MODE $M $N" >> $LOG
      flock $L/testing-gpu.lock flock $L/testing-gpu-GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9.lock /home/treefor/projects/runtimes/core-lease-20261030/any-core.sh --wait \
        /usr/bin/python -I run_hall_v34.py gpu $LABEL bench-$MODE-$M-$N test_v34 bench {TRACE} {OUT} {NVRTC} {UUID} $M $N $S $MODE >> $LOG 2>&1
      echo "$(date +%H:%M:%S) END rc=$? $MODE $M $N" >> $LOG
    done
  done
done
echo "$(date +%H:%M:%S) DONE" >> $LOG
