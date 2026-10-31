#!/bin/sh
P=/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-full-time-candidate006
MFEM=/home/treefor/projects/runtimes/mfem-4.10-pic-20261029
VENDOR=/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/source/vendor
mkdir -p $P/build
/usr/bin/mpic++ -std=c++20 -O0 -g0 -ffp-contract=off -fno-fast-math -I$P/source/include -isystem $P/inputs -I$VENDOR -isystem $MFEM/install-mpi/include -isystem $MFEM/tpl/hypre/include -isystem $MFEM/tpl/metis/include -DPIC_REACTION_NATIVE_SOURCE_SHA256=\"syntaxcheck\" -Wall -Wextra -Wpedantic -Werror -fsyntax-only "$@"
