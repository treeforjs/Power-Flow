#!/bin/sh
# Development build of the GPU parity driver (qualification builds: run_checks.py).
set -eu
R=/home/treefor/projects/runtimes/pic-water-chemistry-v3-20261031
ARR=/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include/pf_arrival
GPU=/home/treefor/projects/runtimes/hall-gpu-resident-v1-20261030
NVRTC_DIR=/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc
CUDA_INC=/home/treefor/projects/runtimes/mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include
exec /usr/bin/c++ -std=c++20 -O2 -Wall -Wextra -Wno-c++26-extensions -ffp-contract=off -fno-fast-math \
  -I "$R/include/pf_chem_v3" -I "$ARR" -I "$GPU" -isystem "$CUDA_INC" -isystem "$NVRTC_DIR/include" \
  "$R/tests/gpu_parity.cpp" "$GPU/common/hall_gpu_common.cpp" -ldl -lcrypto -o "$1"
