#!/bin/bash
# Syntax-only probe of the v4 translation units (no artifacts produced).
R=/home/treefor/projects/runtimes; P=$R/pic-water-common-commit-v4-20261030
N=$R/pic-reaction-transactions-20261029/native-full-time-candidate003; I=$R/pic-water-full-time-subflow-20261029/inputs
H=$R/hall-coupled-transport-20261029/source; M=$R/mfem-4.10-pic-20261029; V3=$R/pic-water-common-commit-v3-20261030
FL="-std=c++20 -O2 -ffp-contract=off -fno-fast-math -Wall -Wextra -Wpedantic -Werror -I$P/source -I$V3/source -DPIC_REACTION_NATIVE_SOURCE_SHA256=\"c045c42d2f57082e3fb5e12172c440d2470e681a38df9ac406bd7e063609a08d\""
for d in $N/source/include $I $I/include $N/inputs $H/development/hall-field-correction-20261030 $H/development $M/install-mpi/include $M/tpl/hypre/include $M/tpl/metis/include; do FL="$FL -isystem $d"; done
ulimit -v 3145728
$R/arch-base-20261028/run-frozen /usr/bin/mpic++ $FL -fsyntax-only "$@"
