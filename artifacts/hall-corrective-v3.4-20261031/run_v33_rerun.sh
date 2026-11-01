#!/bin/bash
# Reruns, on the sealed Hall v3.3 build s1 (runtimes/hall-corrective-v3.3-20261031, read in place), the acceptance stages
# that were stopped to free the GPU: fresh-process restore-continue, the dense blobs and the v7 sink runs. Each stage
# takes both GPU lock files, then a core lease (any-core.sh), and runs through run_hall_v33.py (frozen launcher, receipts
# in the v3.3 outputs s1/). Usage: run_v33_rerun.sh LOGFILE
set -u
LOG=${1:?log}
D=/home/treefor/projects/runtimes/hall-corrective-v3.3-20261031
L=/home/treefor/projects/.workspace/outputs/pic-profile
CL=/home/treefor/projects/runtimes/core-lease-20261030
O=/home/treefor/projects/.workspace/outputs/hall-corrective-v3.3-20261031
cd $D
gpu() { echo "$(date +%H:%M:%S) START $*" >> $LOG
  flock $L/testing-gpu.lock flock $L/testing-gpu-GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9.lock $CL/any-core.sh --wait \
    /usr/bin/python -I run_hall_v33.py gpu "$@" >> $LOG 2>&1
  echo "$(date +%H:%M:%S) END rc=$? $*" >> $LOG; }
LR=$O/s1/gpu-ledger-run-122842
gpu s1 restore-continue test_v32 restore-continue {TRACE} $LR/checkpoint-end $LR/continuation-reference.json {OUT} {NVRTC} {UUID}
gpu s1 dense-light test_v32 dense {TRACE} {OUT} 3.45e10 3 1e-14 400 gpu {NVRTC} {UUID}
gpu s1 dense-medium test_v32 dense {TRACE} {OUT} 8e10 3 1e-14 400 gpu {NVRTC} {UUID}
gpu s1 blob10 test_v32 dense {TRACE} {OUT} 1.6e11 10 5.4e-15 150 gpu {NVRTC} {UUID}
gpu s1 v7-sinks-gpu test_v32 dense {TRACE} {OUT} 3.45e10 3 1e-14 150 gpu {NVRTC} {UUID} v7 sinks
gpu s1 v7-sinks-cpu test_v32 dense {TRACE} {OUT} 3.45e10 3 1e-14 150 cpu {NVRTC} {UUID} v7 sinks
echo "$(date +%H:%M:%S) DONE" >> $LOG
