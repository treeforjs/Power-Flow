#!/bin/bash
# GPU series via SimQueue (SimQueue orders the cores and the GPU itself; one job per run, each well under 10 minutes).
# Usage: submit_gpu_series.sh LABEL
set -u
L=${1:-b9}
D=/home/treefor/projects/runtimes/resolution-series-prep-20261031
sub() {   # name gpu-mem time-min args...
  local name=$1 mem=$2 t=$3; shift 3
  simq submit -n rs-$name -p powerflow -c 1 --mem 8G --gpus 1 --gpu-mem $mem -t $t -C $D -- /usr/bin/python -I build_run.py run $L $name "$@" | tail -1
}
sub gmg-360x4-smooth 3G 15 gmg_solve {OUT} {TRACE} {NVRTC} {UUID} 360 4 1 smooth 0,1,4 1e-14 1e-13 5e-13 1e-12 2e-12
sub gpu-vacuum-180x16 4G 12 field_scaling gpu-vacuum {OUT} {TRACE} {NVRTC} {UUID} 180 16 4 300 1e-14 1e-14 100
sub gpu-vacuum-90x8 2G 8 field_scaling gpu-vacuum {OUT} {TRACE} {NVRTC} {UUID} 90 8 2 600 1e-14 1e-14 100
sub gpu-dtscan-180x16 4G 15 field_scaling gpu-dtscan {OUT} {TRACE} {NVRTC} {UUID} 180 16 4 10 1e-14 1e-13 5e-13 1e-12
for n in 45 90 180; do
  sub gpu-coupled-${n}x4-c2 2G 10 field_scaling gpu-coupled {OUT} {TRACE} {NVRTC} {UUID} $n 4 1 200 2 1e-14
  sub gpu-coupled-${n}x4-c10000 2G 10 field_scaling gpu-coupled {OUT} {TRACE} {NVRTC} {UUID} $n 4 1 40 10000 1e-14
  sub gpu-coupled-${n}x4-c100000 3G 12 field_scaling gpu-coupled {OUT} {TRACE} {NVRTC} {UUID} $n 4 1 8 100000 1e-14
done
sub gpu-coupled-180x16-c2 4G 12 field_scaling gpu-coupled {OUT} {TRACE} {NVRTC} {UUID} 180 16 4 100 2 1e-14
sub gpu-coupled-180x16-c10000 4G 12 field_scaling gpu-coupled {OUT} {TRACE} {NVRTC} {UUID} 180 16 4 30 10000 1e-14
