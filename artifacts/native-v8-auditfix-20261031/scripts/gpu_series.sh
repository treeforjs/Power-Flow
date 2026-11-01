#!/bin/sh
# Native v8 GPU series (run through SimQueue: one core, one GPU share): device chemistry
# controls, the device_verified regression of the v7.1 dangling host-record defect, and the
# device-mode step timings. usage: scripts/gpu_series.sh OUTDIR "N..." "MESH..." "WEIGHT..."
set -u
O=$1; mkdir -p $O
NV=/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/lib/libnvrtc.so.12
U=GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9
D=/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-full-time-candidate008-auditfix
F=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
export CUDA_VISIBLE_DEVICES=$U OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
nice -n 15 $F $D/build/test_v8_chemistry_gpu $O/gpu-controls $NV $U > $O/gpu-controls.log 2>&1; echo "controls exit $?" >> $O/progress.txt
PF_V8_CHEM=v3.1 nice -n 15 $F $D/build/test_v8_chemistry_gpu $O/gpu-controls-v31 $NV $U > $O/gpu-controls-v31.log 2>&1; echo "controls v3.1 exit $?" >> $O/progress.txt
export PF_V8_CHEM=v3.1
nice -n 15 $F $D/build/test_step_timing $O/regression-device_verified-45-10000.json device_verified 10000 4 45 $NV $U > $O/regression.log 2>&1; echo "regression exit $?" >> $O/progress.txt
for w in $4; do for mesh in $3; do for n in $2; do
  steps=5; [ $n -ge 1000000 ] && steps=4
  nice -n 15 $F $D/build/test_step_timing_v8 $O/v8-device-$mesh-$n-w$w.json device $n $steps $mesh $w $NV $U > $O/v8-device-$mesh-$n-w$w.log 2>&1
  echo "$mesh $n $w exit $?" >> $O/progress.txt
done; done; done
