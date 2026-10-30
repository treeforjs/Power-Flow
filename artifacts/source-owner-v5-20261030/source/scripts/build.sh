#!/bin/bash
# Configure and build the owner v5 targets, then run the compile-reject cases.
# Runs inside run-frozen under both-cores.sh (cores 6,7; -j2; nice 15).
# usage: build.sh <attempt-directory>
set -u
V5=/home/treefor/projects/runtimes/pic-electrode-heat-source-v5-20261030
A=$1; BUILD=$V5/build
COUNT=/home/treefor/projects/runtimes/pic-water-exact-count-subflow-20261029/native/source/include/raikou/exact_particle_count.h
MFEM_DIR=/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/cmake/mfem
TARGETS="electrode-heat-source-v5-tests electrode-heat-source-v5-sanitized source_v5_port_normal source_v5_port_sanitized min_repro_a08_gaps_v5_normal min_repro_a08_gaps_v5_sanitized"
t0=$(date +%s.%N)
/usr/bin/cmake -S $V5 -B $BUILD -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/sbin/c++ \
  -DEXACT_PARTICLE_COUNT_HEADER=$COUNT -DMFEM_DIR=$MFEM_DIR -DCMAKE_EXPORT_COMPILE_COMMANDS=ON > $A/logs/configure.log 2>&1 || { echo "configure failed"; exit 10; }
t1=$(date +%s.%N)
/usr/bin/cmake --build $BUILD -j 2 --target $TARGETS > $A/logs/build.log 2>&1 || { echo "build failed"; tail -50 $A/logs/build.log; exit 11; }
t2=$(date +%s.%N)
/usr/bin/python3 -I $V5/scripts/reject_cases.py $BUILD $V5/audit-port/compile_reject.cpp $A/reject || { echo "compile-reject driver failed"; exit 12; }
t3=$(date +%s.%N)
/usr/bin/python3 -I -c "import sys;a=[float(x) for x in sys.argv[1:]];print('configure_s=%.2f build_s=%.2f reject_s=%.2f'%(a[1]-a[0],a[2]-a[1],a[3]-a[2]))" $t0 $t1 $t2 $t3
