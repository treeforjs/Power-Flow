#!/bin/bash
set -euo pipefail
export CUDA_VISIBLE_DEVICES= OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/usr/bin/python build-components.py normal core hall v8-sinks step-timing audit-target-guard audit-parcels-codec v8-chemistry v8-restore flows work-gross v8-stage-b v8-chemistry-gpu v8-flows-gpu
/usr/bin/python audit-tests.py cpu final-cpu-v4
