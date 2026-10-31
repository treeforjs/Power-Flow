#!/bin/sh
# Run the c19-circuit-drive-v1 test batch inside the frozen root on one leased core (any-core.sh wraps this script).
# Logs: /home/treefor/projects/.workspace/outputs/c19-circuit-drive-v1-20261031/tests/
set -u
P=/home/treefor/projects/runtimes/c19-circuit-drive-v1-20261031
OUT=/home/treefor/projects/.workspace/outputs/c19-circuit-drive-v1-20261031
B=$OUT/build; T=$OUT/tests; mkdir -p $T
TRACE=/home/treefor/projects/Power-Flow/data/dissertation/shot_17133_current.csv
LAW=$P/calibration/hall_ss304_law_table.csv
rc=0
for v in A_supplied C_traces C3_traces_45mm; do
  PAR=$P/calibration/c19_circuit_v1_parameters_$v.txt
  $B/test_circuit $PAR $LAW > $T/test_circuit_$v.log 2>&1 || rc=1
  $B/test_hall_port $TRACE $PAR $LAW 2000 2000 1000 $T/hall_${v}_ > $T/test_hall_port_$v.log 2>&1 || rc=1
done
grep -h -E "^(PASS|FAIL|ALL|FAILED|INFO)" $T/test_circuit_*.log $T/test_hall_port_*.log > $T/summary.txt
exit $rc
