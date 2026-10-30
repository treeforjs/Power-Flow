#!/bin/bash
set -euo pipefail
p=/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-collision-scheduler-20261029/review-v3
exec 9>/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core7.lock
flock -x 9
cd "$p"
taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= python3 generate_weights.py
taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= g++ -std=c++20 -O2 -ffp-contract=off -Wall -Wextra adversarial.cpp -o adversarial >compile-adversarial.log 2>&1
taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= ./adversarial weights.csv >adversarial.txt
