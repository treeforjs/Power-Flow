#!/bin/bash
set -euo pipefail
export CUDA_VISIBLE_DEVICES= OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/usr/bin/python build-components.py normal v8-chemistry
/usr/bin/python audit-tests.py cpu final-cpu
