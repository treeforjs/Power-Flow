#!/bin/bash
set -euo pipefail
exec 9>/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core7.lock
flock -x 9
cd /home/treefor/projects/Power-Flow
taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= /home/treefor/projects/Power-Flow/.venv/bin/python /home/treefor/Documents/Codex/2026-10-28/we/outputs/cathode-emission-audit-20261029/preview/render_source_cloud.py
