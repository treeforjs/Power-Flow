#!/bin/bash
# measure.sh BUILD OUTDIR: kernel breakdown of the Hall v3.4 coupled step (test_v34 bench, native v7.1 lane configuration,
# 200 um) at 1e5 and 1e6 carriers, each with 2 and 5 steps (the difference gives three steady steps without setup and
# first-step work). Runs inside a SimQueue powerflow job (1 core, shared GPU). Measurement only (README.txt).
set -eu
B=${1:?build}; O=${2:?outdir}
R=/home/treefor/projects/runtimes
T=$R/hall-corrective-v3.4-20261031/build/$B/test_v34
TRACE=$R/hall-coupled-transport-20261029/source/development/hall-finite-foil-em-20261029/data/shot_17133_current.csv
NVRTC=$R/pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/lib/libnvrtc.so.12
UUID=GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9
cd /home/treefor/projects
for N in 100000 1000000; do
  for S in 2 5; do
    W=$O/run-$N-$S; mkdir -p $W
    echo "$(date +%H:%M:%S) START $N $S"
    nice -n 15 $R/arch-base-20261028/run-frozen env LD_LIBRARY_PATH=$O/ktimer HALL_KTIMER_OUT=$W/kernels.json \
      OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 $T bench $TRACE $W $NVRTC $UUID 45 $N $S lane > $W/bench.log 2>&1
    echo "$(date +%H:%M:%S) END rc=$? $N $S"
  done
done
echo DONE
