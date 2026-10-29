#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
b_prefix=/home/treefor/projects/runtimes/pic-native-b-export-20261029
b_frozen=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
"$b_frozen" /usr/bin/cmake -S "$b_prefix/source" -B "$b_prefix/build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DMFEM_DIR=/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/cmake/mfem \
    > "$b_prefix/logs/configure.log" 2>&1
"$b_frozen" /usr/bin/cmake --build "$b_prefix/build" --parallel 1 \
    > "$b_prefix/logs/build.log" 2>&1
"$b_frozen" /usr/bin/ldd "$b_prefix/build/native_b_export" > "$b_prefix/logs/ldd.txt"
if /usr/bin/rg 'libcuda|libcudss|not found' "$b_prefix/logs/ldd.txt"; then
    exit 1
fi
"$b_frozen" /usr/bin/nm -C "$b_prefix/build/native_b_export" > "$b_prefix/logs/symbols.txt"
if /usr/bin/rg ' [Tt] sealed_pic_driver_main_unused\(|GapPic::advance|GapPic::solve_fields|CudssPoisson' "$b_prefix/logs/symbols.txt"; then
    exit 1
fi
printf '%s\n' 'Build passed; unused driver/PIC solve symbols removed; no CUDA/cuDSS library linked.'
