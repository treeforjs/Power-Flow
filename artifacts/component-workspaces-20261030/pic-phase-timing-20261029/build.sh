#!/bin/bash
set -euo pipefail
phase_prefix=/home/treefor/projects/runtimes/pic-phase-timing-20261029
phase_frozen=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
for phase_variant in baseline cpu-v2; do
    "$phase_frozen" /usr/bin/cmake \
        -S "$phase_prefix/$phase_variant/source" \
        -B "$phase_prefix/$phase_variant/build" \
        -DCMAKE_BUILD_TYPE=Release \
        -DMFEM_DIR=/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/cmake/mfem \
        > "$phase_prefix/logs/$phase_variant-configure.log" 2>&1
    "$phase_frozen" /usr/bin/cmake --build "$phase_prefix/$phase_variant/build" \
        --target heating_history_pic_gpu --parallel 1 \
        > "$phase_prefix/logs/$phase_variant-build.log" 2>&1
    printf '%s\n' "$phase_variant build passed"
done
"$phase_frozen" /usr/bin/c++ -std=c++20 -O2 -Wall -Wextra -pedantic \
    "$phase_prefix/phase_accounting_test.cpp" -o "$phase_prefix/phase_accounting_test" \
    > "$phase_prefix/logs/accounting-build.log" 2>&1
"$phase_frozen" "$phase_prefix/phase_accounting_test" "$phase_prefix/evidence/accounting-test.json" \
    > "$phase_prefix/logs/accounting-test.log" 2>&1
cat "$phase_prefix/logs/accounting-test.log"
