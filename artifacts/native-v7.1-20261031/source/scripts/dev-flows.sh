#!/bin/sh
# Development cycle: rebuild selected normal targets, then run one test.
# usage: scripts/dev-flows.sh TARGET BINARY OUTDIR [ARGS...]
set -e
cd /home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-full-time-candidate007.1
T=$1; BIN=$2; D=$3; shift 3
/home/treefor/projects/runtimes/core-lease-20261030/both-cores.sh env CUDA_VISIBLE_DEVICES= python3 build-components.py normal $T > $D.build.log 2>&1 || { grep -E "error" -A3 $D.build.log $(ls -t /home/treefor/projects/.workspace/outputs/native-v7.1-20261031/build-normal/*.log | head -3) | head -40; exit 3; }
mkdir -p $D
/home/treefor/projects/runtimes/core-lease-20261030/any-core.sh --wait env CUDA_VISIBLE_DEVICES= /home/treefor/projects/runtimes/arch-base-20261028/run-frozen prlimit --as=3221225472 -- $BIN $D "$@" > $D.log 2>&1
