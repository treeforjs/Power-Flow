#!/bin/sh
# Development build of the CPU tests (qualification builds: run_checks.py).
# usage: tools/dev_build.sh [normal|sanitized] TARGET.cpp OUTPUT
# PF_BUILD_ISA_FLAGS (whitespace-split, default empty) is appended to the C++ compile flags: hook for a per-site ISA flavour
# (e.g. -march=x86-64-v4); the qualification receipt records its value.
set -eu
R=/home/treefor/projects/runtimes/pic-water-chemistry-v3.2.2-auditprep-20261101
ARR=/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include/pf_arrival
mode=$1; src=$2; out=$3
if [ "$mode" = sanitized ]; then
  OPT="-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined -fno-sanitize-recover=all"
else
  OPT="-O2"
fi
exec /usr/bin/c++ -std=c++20 $OPT -Wall -Wextra -Wpedantic -Werror -ffp-contract=off -fno-fast-math \
  ${PF_BUILD_ISA_FLAGS:-} -I "$R/include/pf_chem_v3" -I "$ARR" "$src" -lcrypto -o "$out"
