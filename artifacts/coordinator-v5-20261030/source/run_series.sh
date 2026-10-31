#!/bin/bash
# Acceptance (b) GPU step and (c) series for one build label, sequentially, each under its leases.
# Usage: run_series.sh LABEL
set -u
L=/home/treefor/projects/.workspace/outputs/pic-profile
R=/home/treefor/projects/runtimes
P=$R/pic-water-common-commit-v5-20261030
LABEL=$1
cd $P
gpu() { flock $L/testing-gpu.lock flock $L/testing-gpu-GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9.lock $R/core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_v5.py "$@"; }
cpu() { $R/core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_v5.py "$@"; }
BIG="--carrier-cap 1000000 --hall-budget 8e9 --native-budget 4e9"
echo "== gpu step"; gpu gpu $LABEL
echo "== multistep 2"; gpu multistep $LABEL gpu 2 2000 1; cpu multistep $LABEL cpu 2 2000 0
echo "== multistep 1e4"; gpu multistep $LABEL gpu 10000 200 20; cpu multistep $LABEL cpu 10000 200 0
echo "== hall-loop 1e4"; gpu hall-loop $LABEL gpu 10000 50 10; cpu hall-loop $LABEL cpu 10000 50 0
echo "== multistep 1e5"; gpu multistep $LABEL gpu 100000 30 10; cpu multistep $LABEL cpu 100000 30 0
echo "== hall-loop 1e5"; gpu hall-loop $LABEL gpu 100000 20 10; cpu hall-loop $LABEL cpu 100000 20 0
echo "== 1e6 (GPU, configured caps)"; gpu hall-loop $LABEL gpu 1000000 10 0 $BIG; gpu multistep $LABEL gpu 1000000 5 0 $BIG
echo "== summary"; /usr/bin/python -I series_summary.py /home/treefor/projects/.workspace/outputs/pic-water-common-commit-v5-20261030/series-$LABEL
echo "== series done"
