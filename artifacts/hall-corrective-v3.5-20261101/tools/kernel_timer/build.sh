#!/bin/sh
# build.sh OUTDIR: the stand-in libcuda.so.1 (frozen toolchain; run inside run-frozen).
set -eu
OUT=${1:?outdir}; D=$(dirname "$0")
mkdir -p "$OUT"
nm -D --defined-only /usr/lib/libcuda.so.615.78.08 | awk '$2=="T"{print $3}' > "$OUT/symbols.txt"
/usr/bin/python -I "$D/gen_stubs.py" "$OUT/symbols.txt" "$OUT/stubs.S" "$OUT/stubs.h"
gcc -O2 -fPIC -shared -Wl,-soname,libhallktimer.so -I"$OUT" -o "$OUT/libcuda.so.1" "$D/shim.c" "$OUT/stubs.S" -ldl
echo built "$OUT/libcuda.so.1"
