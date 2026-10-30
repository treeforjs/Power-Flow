#!/bin/bash
# Serial bounded execution; failure to acquire CPU7 is a queue event, not permission to use another CPU.
set -euo pipefail
cd "$(dirname "$0")/.."
exec flock -n /home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core7.lock \
  taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen \
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  CUDA_VISIBLE_DEVICES= PYTHONDONTWRITEBYTECODE=1 \
  /home/treefor/projects/Power-Flow/.venv/bin/python "$@"
