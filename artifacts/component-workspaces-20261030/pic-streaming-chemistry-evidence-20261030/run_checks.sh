#!/bin/bash
set -euo pipefail
exec flock <workspace>/outputs/pic-profile/testing-core7.lock taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen /usr/bin/env -u PYTHONOPTIMIZE -u PYTHONSTARTUP OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= /usr/bin/python /home/treefor/projects/runtimes/pic-streaming-chemistry-evidence-20261030/build_and_check.py
