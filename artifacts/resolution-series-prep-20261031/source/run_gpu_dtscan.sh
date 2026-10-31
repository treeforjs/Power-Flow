#!/bin/bash
# CG iterations and time per step versus dt (field-only, sealed Hall v3.1 resident path, dt cap raised for the scan only).
cd /home/treefor/projects/runtimes/resolution-series-prep-20261031
O=/home/treefor/projects/.workspace/outputs/resolution-series-prep-20261031
while pgrep -f run_gpu_series.sh > /dev/null; do sleep 10; done
for n in 45 90 180 360; do
  /home/treefor/projects/runtimes/core-lease-20261030/any-core.sh --wait /usr/bin/python -I build_run.py run b1 gpu-dtscan-${n}x4 field_scaling \
    gpu-dtscan {OUT} {TRACE} {NVRTC} {UUID} $n 4 1 20 1e-14 2e-14 5e-14 1e-13 2e-13 5e-13 1e-12 > $O/gpu-dtscan-${n}x4.out 2>&1
done
