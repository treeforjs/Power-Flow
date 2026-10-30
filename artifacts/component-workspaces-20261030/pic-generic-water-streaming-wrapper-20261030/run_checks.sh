#!/bin/bash
set -euo pipefail
runtime_prefix="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec flock /home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core7.lock taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen /usr/bin/env -u PYTHONOPTIMIZE -u PYTHONSTARTUP OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= /usr/bin/python "$runtime_prefix/build_and_check.py"
