#!/bin/bash
# measure_modes.sh BUILD: per-kernel device time per coupled step (tools/kernel_timer stand-in libcuda.so.1, built in ktimer/)
# for the timing table against v3.3: test_v34 bench in modes v33 (v3.3 host path of the resident owner), lane (v3.4,
# native v7.1 configuration) and owner (v3.3 host Owner), 200 um (45) and 50 um (180), 1e5 and 1e6 carriers, 2 and 5
# steps each (steady step = difference / 3). Lane at 200 um is job 189 (run-N-S). Priority order; one SimQueue powerflow
# job (1 core, shared GPU). Measurement only.
set -u
B=${1:?build}
O=/home/treefor/projects/.workspace/outputs/hall-corrective-v3.4-20261031/gpu-kernel-measure
R=/home/treefor/projects/runtimes
T=$R/hall-corrective-v3.4-20261031/build/$B/test_v34
TRACE=$R/hall-coupled-transport-20261029/source/development/hall-finite-foil-em-20261029/data/shot_17133_current.csv
NVRTC=$R/pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/lib/libnvrtc.so.12
UUID=GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9
cd /home/treefor/projects
echo "binary $(sha256sum $T)"
one() {  # MODE TRANSVERSE CARRIERS
  for S in 2 5; do
    W=$O/m2/run-$1-$2-$3-$S; mkdir -p $W
    echo "$(date +%H:%M:%S) START $1 $2 $3 $S"
    nice -n 15 $R/arch-base-20261028/run-frozen env LD_LIBRARY_PATH=$O/ktimer HALL_KTIMER_OUT=$W/kernels.json \
      OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 $T bench $TRACE $W $NVRTC $UUID $2 $3 $S $1 > $W/bench.log 2>&1
    echo "$(date +%H:%M:%S) END rc=$? $1 $2 $3 $S"
  done
}
one v33 45 100000;   one v33 45 1000000
one v33 180 100000;  one v33 180 1000000
one lane 180 100000; one lane 180 1000000
one owner 45 100000; one owner 45 1000000
one owner 180 100000; one owner 180 1000000
echo DONE
