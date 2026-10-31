#!/bin/sh
# Build the c19-circuit-drive-v1 binaries inside the frozen root (run through any-core.sh or simq).
set -eu
R=/home/treefor/projects/runtimes; P=$R/c19-circuit-drive-v1-20261031
OUT=/home/treefor/projects/.workspace/outputs/c19-circuit-drive-v1-20261031/build; mkdir -p $OUT
N=$R/pic-reaction-transactions-20261029/native-full-time-candidate003; I=$R/pic-water-full-time-subflow-20261029/inputs
H=$R/hall-coupled-transport-20261029/source; M=$R/mfem-4.10-pic-20261029; V31=$R/hall-corrective-v3.1-20261031/source
MV=$H/development/hall-canonical-driver-20261029/bindings/heat-v2/material-v2
COMMON="-std=c++20 -O2 -DNDEBUG -ffp-contract=off -fno-fast-math"
ISYS="-isystem $N/source/include -isystem $I -isystem $I/include -isystem $N/inputs -isystem $H/development -isystem $M/install-mpi/include -isystem $M/tpl/hypre/include -isystem $M/tpl/metis/include"
LIBS="$M/install-mpi/lib/libmfem.a $M/tpl/hypre/lib/libHYPRE.a $M/tpl/metis/lib/libmetis.a $M/tpl/GKlib/lib/libGKlib.a -lcrypto"
cd $P
/usr/bin/c++ $COMMON -Wall -Wextra tools/c19_sim.cpp -o $OUT/c19_sim
/usr/bin/c++ $COMMON -Wall -Wextra tests/test_circuit.cpp -o $OUT/test_circuit
/usr/bin/c++ $COMMON -I$MV tools/ss304_law_table.cpp -o $OUT/ss304_law_table
/usr/bin/mpic++ $COMMON $ISYS tools/hall_port_load.cpp -Wl,--gc-sections $LIBS -o $OUT/hall_port_load
/usr/bin/mpic++ $COMMON -I$V31 $ISYS tests/test_hall_port.cpp -Wl,--gc-sections $LIBS -o $OUT/test_hall_port
sha256sum $OUT/c19_sim $OUT/test_circuit $OUT/ss304_law_table $OUT/hall_port_load $OUT/test_hall_port > $OUT/binaries.sha256
cat $OUT/binaries.sha256
