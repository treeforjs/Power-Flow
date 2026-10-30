#!/bin/bash
# usage: syntax.sh <source> [extra flags]: -fsyntax-only compile with the component flags
. /home/treefor/projects/runtimes/pic-electrode-heat-source-v5-20261030/scripts/flags.sh
src=$1; shift
eval /usr/sbin/c++ $FLAGS "$@" -fsyntax-only "$src"
