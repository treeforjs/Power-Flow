#!/bin/sh
# Host-mode native v8 step timings (CPU only). usage: scripts/timing_host.sh OUTDIR "N..." "MESH..." "WEIGHT..."
O=$1; mkdir -p $O
export PF_V8_CHEM=${PF_V8_CHEM:-v3.1}
NV=/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/lib/libnvrtc.so.12
U=GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9
D=/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-full-time-candidate008
for w in $4; do for mesh in $3; do for n in $2; do
  steps=5; [ $n -ge 1000000 ] && steps=3
  $D/scripts/lease.sh env CUDA_VISIBLE_DEVICES= /home/treefor/projects/runtimes/arch-base-20261028/run-frozen prlimit --as=8589934592 -- \
    $D/build/test_step_timing_v8 $O/v8-host-$mesh-$n-w$w.json host $n $steps $mesh $w $NV $U > $O/v8-host-$mesh-$n-w$w.log 2>&1
  echo "$mesh $n $w exit $?" >> $O/progress.txt
done; done; done
