#!/bin/bash
set -euo pipefail
: "${SIMQ_JOB_ID:?submit this check through SimQueue}"
test -n "${SIMQ_GPUS:-}" || { echo "GPU allocation required" >&2; exit 2; }
cd /home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031
for process in 1 2; do
 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen build/bin/gpu_parity \
 data/tables/chemistry-v3.1.blob "evidence/gpu-process${process}" \
 /home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/lib/libnvrtc.so.12 \
 GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9 2000 1
done
