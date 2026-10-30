#!/bin/bash
set -euo pipefail
exec flock /home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core7.lock taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen env PIC_GENERIC_REBUILD="${PIC_GENERIC_REBUILD:-0}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= bash -c '
set -euo pipefail
ulimit -v 1572864
ulimit -t 180
p=/home/treefor/projects/runtimes/pic-generic-water-bridge-20261029
if [ "$PIC_GENERIC_REBUILD" = 1 ]; then python "$p/rebuild_foundation.py"; fi
python "$p/build_and_check.py"
python "$p/test_native.py"
python "$p/independent_oracle.py"
'
