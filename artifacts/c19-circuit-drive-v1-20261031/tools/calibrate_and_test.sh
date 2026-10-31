#!/bin/sh
# Final calibration then the test batch (2026-10-31; lent cores 1 and 3, see tools/lent-core.sh).
P=/home/treefor/projects/runtimes/c19-circuit-drive-v1-20261031
OUT=/home/treefor/projects/.workspace/outputs/c19-circuit-drive-v1-20261031
R=/home/treefor/projects/runtimes
cd $P
tools/lent-core.sh 1 /home/treefor/projects/Power-Flow/.venv/bin/python -I calibration/calibrate.py $OUT/build/c19_sim calibration/hall_ss304_law_table.csv calibration $OUT/calib-scratch > $OUT/calibrate.log 2>&1
tools/lent-core.sh 3 sh tools/run_tests.sh > $OUT/tests/run.log 2>&1
echo done >> $OUT/calibrate.log
