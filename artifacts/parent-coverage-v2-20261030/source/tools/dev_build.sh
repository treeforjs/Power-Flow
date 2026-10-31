#!/bin/bash
# Development build helper (the qualified build is run_coverage_v2_checks.py).
# Usage: dev_build.sh SOURCE OUTPUT [normal|ubsan]
set -eu
R=/home/treefor/projects/runtimes
P=$R/pic-water-parent-coverage-v2-20261030
N=$R/pic-reaction-transactions-20261029/native-full-time-candidate003
I=$R/pic-water-full-time-subflow-20261029/inputs
M=$R/mfem-4.10-pic-20261029
ID=c045c42d2f57082e3fb5e12172c440d2470e681a38df9ac406bd7e063609a08d
opt="-O2"; [ "${3:-normal}" = ubsan ] && opt="-O1 -g -fsanitize=undefined -fno-sanitize-recover=all"
cd "$P"
exec $R/arch-base-20261028/run-frozen /usr/bin/mpic++ -std=c++20 -Wall -Wextra -Wpedantic -Werror -ffp-contract=off -fno-fast-math \
  -ffunction-sections -fdata-sections -I$P/source -I$R/pic-water-parent-coverage-v1-20261030/source -I$R/pic-water-common-commit-v3-20261030/source \
  -DPIC_REACTION_NATIVE_SOURCE_SHA256=\"$ID\" -isystem $N/source/include -isystem $I -isystem $I/include -isystem $N/inputs \
  -isystem $M/install-mpi/include -isystem $M/tpl/hypre/include -isystem $M/tpl/metis/include $opt "$1" \
  $N/build/libpic_native_candidate003.a $M/install-mpi/lib/libmfem.a $M/tpl/hypre/lib/libHYPRE.a $M/tpl/metis/lib/libmetis.a \
  $M/tpl/GKlib/lib/libGKlib.a -Wl,--gc-sections -lcrypto -o "$2"
