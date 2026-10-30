#!/bin/bash
set -euo pipefail
p=/home/treefor/sim-workspace/powerflow/outputs/pic-collision-scheduler-20261029/review-v2
exec 9>/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core7.lock
flock -x 9
cd "$p"
taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= python3 generate_cases.py
taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= g++ -std=c++20 -O2 -ffp-contract=off kinematics.cpp -o kinematics >compile-kinematics.log 2>&1
taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= ./kinematics cases.csv >native-kinematics.csv
