#!/usr/bin/env bash
set -eu
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
runtime=/home/treefor/projects/runtimes/pic-source-sampling-validation-20261029
output=<workspace>/outputs/pic-profile/source-sampling-validation
exec 9>><workspace>/outputs/pic-profile/testing-core.lock
if ! taskset -c 6 nice -n 15 flock -n 9; then
  printf '%s\n' 'Testing CPU6 lease is busy; no build or replay started.' >&2
  exit 75
fi
taskset -c 6 nice -n 15 python3 "$runtime/extract_native.py"
taskset -c 6 nice -n 15 /usr/sbin/c++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Werror=return-type "$runtime/validate.cpp" -o "$runtime/validate" > "$output/build.log" 2>&1
taskset -c 6 nice -n 15 "$runtime/validate" \
  <workspace>/outputs/heating-history-full.bin \
  <workspace>/outputs/pic-profile/particle-audit/copied-checkpoint \
  "$output" > "$output/run.log" 2>&1
taskset -c 6 nice -n 15 python3 "$runtime/analyze.py" > "$output/analysis.log" 2>&1
