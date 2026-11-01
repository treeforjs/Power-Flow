#!/bin/bash
set -euo pipefail
: "${SIMQ_JOB_ID:?submit through SimQueue}"
[[ -z "${PF_BUILD_ISA_FLAGS:-}" ]] || { echo 'unqualified build flags refused' >&2; exit 2; }
cd /home/treefor/projects/runtimes/water-network-scope-guard-v1-20261101
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=
mode=$1
case "$mode" in
  cpu) flags=(-O2) ;;
  sanitized) flags=(-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined -fno-sanitize-recover=all) ;;
  *) exit 2 ;;
esac
python3 -I tools/verify_pins.py
/usr/bin/c++ -std=c++20 "${flags[@]}" -Wall -Wextra -Wpedantic -Werror -ffp-contract=off -fno-fast-math \
  -I include src/water_network_scope.cpp tests/test_scope.cpp -lcrypto -o "build/scope-$mode"
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
 "build/scope-$mode" \
 /home/treefor/projects/runtimes/pic-water-chemistry-v3.2.2-auditprep-20261101/data/tables/chemistry-v3.2.blob \
 4ebf093eaf5a7fed50504cf4410b04b6183463faa59d67736da16777f348cafc \
 "evidence/synthetic-contract-$mode.txt" | tee "evidence/check-$mode.log"
python3 -I tools/verify_pins.py
sha256sum "build/scope-$mode" "evidence/check-$mode.log" "evidence/synthetic-contract-$mode.txt" > "evidence/check-$mode.sha256"
