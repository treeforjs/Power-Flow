#!/bin/bash
set -euo pipefail
: "${SIMQ_JOB_ID:?submit this check through SimQueue}"
cd /home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031
/home/treefor/projects/runtimes/arch-base-20261028/run-frozen tools/dev_build_gpu.sh build/bin/gpu_parity
