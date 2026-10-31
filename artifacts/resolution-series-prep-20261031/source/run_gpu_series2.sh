#!/bin/bash
# Second GPU series: axial-refined proxies of the 25 um size (sealed Hall v3.1 admits 180x16+4 and 90x8+2) and the closed
# coupled Hall step (GPU push) with 2, 1e4 and 1e5 carriers per mesh.
cd /home/treefor/projects/runtimes/resolution-series-prep-20261031
O=/home/treefor/projects/.workspace/outputs/resolution-series-prep-20261031
W=/home/treefor/projects/runtimes/core-lease-20261030/any-core.sh
while pgrep -f gmg_solve > /dev/null; do sleep 10; done
$W --wait /usr/bin/python -I build_run.py run b3 gpu-vacuum-180x16 field_scaling gpu-vacuum {OUT} {TRACE} {NVRTC} {UUID} 180 16 4 300 1e-14 1e-14 100 > $O/gpu-vacuum-180x16.out 2>&1
$W --wait /usr/bin/python -I build_run.py run b3 gpu-vacuum-90x8 field_scaling gpu-vacuum {OUT} {TRACE} {NVRTC} {UUID} 90 8 2 600 1e-14 1e-14 100 > $O/gpu-vacuum-90x8.out 2>&1
$W --wait /usr/bin/python -I build_run.py run b3 gpu-dtscan-180x16 field_scaling gpu-dtscan {OUT} {TRACE} {NVRTC} {UUID} 180 16 4 10 1e-14 1e-13 5e-13 1e-12 > $O/gpu-dtscan-180x16.out 2>&1
for n in 45 90 180; do
  for spec in "2 200" "10000 40" "100000 8"; do set -- $spec
    $W --wait /usr/bin/python -I build_run.py run b3 gpu-coupled-${n}x4-c$1 field_scaling gpu-coupled {OUT} {TRACE} {NVRTC} {UUID} $n 4 1 $2 $1 1e-14 > $O/gpu-coupled-${n}x4-c$1.out 2>&1
  done
done
for spec in "2 100" "10000 30"; do set -- $spec
  $W --wait /usr/bin/python -I build_run.py run b3 gpu-coupled-180x16-c$1 field_scaling gpu-coupled {OUT} {TRACE} {NVRTC} {UUID} 180 16 4 $2 $1 1e-14 > $O/gpu-coupled-180x16-c$1.out 2>&1
done
