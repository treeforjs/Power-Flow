#!/bin/bash
set -euo pipefail
exec flock /home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core7.lock taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen /bin/bash -c '
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=
ulimit -v 1048576
ulimit -t 60
p=/home/treefor/projects/runtimes/pic-collision-exposure-clock-20261029
o=/home/treefor/sim-workspace/powerflow/outputs/pic-collision-exposure-clock-20261029
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -I"$p/source" "$p/source/test_existing_contract.cpp" -o "$p/build/test_existing_contract"
"$p/build/test_existing_contract" "$o/existing-contract.json" > "$o/existing-contract.log"
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -I"$p/source" "$p/source/test_exposure_clock.cpp" -o "$p/build/test_exposure_clock"
"$p/build/test_exposure_clock" "$o/exposure-results.json" > "$o/exposure-results.log"
cat "$o/exposure-results.log"
'
