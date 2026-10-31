#!/bin/bash
# Field-only scaling runs of the sealed Hall v3.1 resident path on the T-series meshes (dt = 10 fs).
cd /home/treefor/projects/runtimes/resolution-series-prep-20261031
O=/home/treefor/projects/.workspace/outputs/resolution-series-prep-20261031
for spec in "45 4 1 1000" "90 4 1 600" "180 4 1 400" "360 4 1 300"; do
  set -- $spec
  /home/treefor/projects/runtimes/core-lease-20261030/any-core.sh --wait /usr/bin/python -I build_run.py run b1 gpu-vacuum-$1x$2 field_scaling \
    gpu-vacuum {OUT} {TRACE} {NVRTC} {UUID} $1 $2 $3 $4 1e-14 1e-14 100 > $O/gpu-vacuum-$1x$2.out 2>&1
done
