#!/bin/bash
# Hall corrective v3.1 acceptance batch. Every stage takes its leases through the sealed wrappers
# (both-cores.sh for builds, any-core.sh --wait for CPU runs, both GPU lock files + any-core.sh --wait
# for GPU runs); every compile and run executes inside arch-base-20261028/run-frozen (run_hall_v31.py,
# run_coordinator_parity.py). Usage: run_acceptance.sh LABEL
set -u
LABEL=${1:?label}
D=/home/treefor/projects/runtimes/hall-corrective-v3.1-20261031
L=/home/treefor/projects/.workspace/outputs/pic-profile
CL=/home/treefor/projects/runtimes/core-lease-20261030
LOG=/home/treefor/projects/.workspace/outputs/hall-corrective-v3.1-20261031/acceptance-$LABEL.log
cd $D
step() { echo "$(date +%H:%M:%S) START $*" >> $LOG; "$@" >> $LOG 2>&1; echo "$(date +%H:%M:%S) END rc=$? $*" >> $LOG; }
build() { step $CL/both-cores.sh /usr/bin/python -I run_hall_v31.py build "$@"; }
cpu() { step $CL/any-core.sh --wait /usr/bin/python -I run_hall_v31.py run "$@"; }
gpu() { step flock $L/testing-gpu.lock flock $L/testing-gpu-GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9.lock $CL/any-core.sh --wait /usr/bin/python -I run_hall_v31.py gpu "$@"; }

# SKIP_BUILD=1: LABEL was built beforehand (run_hall_v31.py build LABEL ... on another admitted core)
[ "${SKIP_BUILD:-0}" = 1 ] || build $LABEL test_hall_v3 test_resident profile_v3_baseline profile_v3_baseline_gpu
[ "${WITH_CPU_ACCEPTANCE:-1}" = 1 ] && step $CL/both-cores.sh /usr/bin/python -I run_coordinator_parity.py build $LABEL
# device-resident qualification
gpu $LABEL units test_resident units {TRACE} {OUT} {NVRTC} {UUID}
gpu $LABEL negative test_resident negative {TRACE} {OUT} {NVRTC} {UUID}
gpu $LABEL parity-vacuum test_resident parity-vacuum {TRACE} {OUT} 200 {NVRTC} {UUID}
gpu $LABEL parity-coupled-2 test_resident parity-coupled {TRACE} {OUT} 50 2 {NVRTC} {UUID}
gpu $LABEL parity-coupled-1e4 test_resident parity-coupled {TRACE} {OUT} 3 10000 {NVRTC} {UUID}
# timings and run-to-run identity
gpu $LABEL vacuum10k-a test_resident vacuum-run {TRACE} {OUT} 10000 {NVRTC} {UUID} sample 500
gpu $LABEL vacuum10k-b test_resident vacuum-run {TRACE} {OUT} 10000 {NVRTC} {UUID} sample 500
gpu $LABEL vacuum1k-verify-all test_resident vacuum-run {TRACE} {OUT} 1000 {NVRTC} {UUID} verify-all 100
gpu $LABEL coupled1k-gpu-a test_resident coupled-run {TRACE} {OUT} 1000 2 {NVRTC} {UUID} gpu
gpu $LABEL coupled1k-gpu-b test_resident coupled-run {TRACE} {OUT} 1000 2 {NVRTC} {UUID} gpu
gpu $LABEL coupled1k-cpu test_resident coupled-run {TRACE} {OUT} 1000 2 {NVRTC} {UUID} cpu
gpu $LABEL coupled-1e4 test_resident coupled-run {TRACE} {OUT} 20 10000 {NVRTC} {UUID} gpu
gpu $LABEL coupled-1e5 test_resident coupled-run {TRACE} {OUT} 5 100000 {NVRTC} {UUID} gpu
# before: Hall v3 (sealed sources) with the gpu1 backend, and with cpu-pcg
gpu $LABEL base-vacuum-gpu1 profile_v3_baseline_gpu vacuum {TRACE} {OUT} 200 gpu {NVRTC} {UUID}
gpu $LABEL base-coupled2-gpu1 profile_v3_baseline_gpu coupled {TRACE} {OUT} 100 2 gpu {NVRTC} {UUID}
gpu $LABEL base-coupled1e4-gpu1 profile_v3_baseline_gpu coupled {TRACE} {OUT} 20 10000 gpu {NVRTC} {UUID}
gpu $LABEL base-coupled1e5-gpu1 profile_v3_baseline_gpu coupled {TRACE} {OUT} 6 100000 gpu {NVRTC} {UUID}
cpu $LABEL base-vacuum-cpu profile_v3_baseline vacuum {TRACE} {OUT} 200
cpu $LABEL base-coupled2-cpu profile_v3_baseline coupled {TRACE} {OUT} 100 2
# Hall v3 component tests on the v3.1 sources (CPU backends) and coordinator v5 T1-T7 (skipped with
# WITH_CPU_ACCEPTANCE=0 when they ran under another label)
if [ "${WITH_CPU_ACCEPTANCE:-1}" = 1 ]; then
cpu $LABEL components test_hall_v3 components {TRACE} {OUT}/components
cpu $LABEL parity-dump test_hall_v3 parity-dump {TRACE} {OUT}/parity-v3.json
cpu $LABEL motion-dump test_hall_v3 motion-dump {TRACE} {OUT}/motion-dump-v3.json
cpu $LABEL sink-dump test_hall_v3 sink-dump {TRACE} {OUT}/sink-dump-cpu.json
step $CL/any-core.sh --wait /usr/bin/python -I run_coordinator_parity.py run $LABEL
fi
# diagnostic step-size study (non-default dt upper bound in this study only)
gpu $LABEL dt-study-100ps test_resident dt-study {TRACE} {OUT} {NVRTC} {UUID} 1e-10 2.5e-15 5e-15 1e-14 2e-14 5e-14 1e-13 2e-13 5e-13 1e-12
gpu $LABEL dt-study-1ns test_resident dt-study {TRACE} {OUT} {NVRTC} {UUID} 1e-9 1e-14 1e-13 1e-12
echo "$(date +%H:%M:%S) DONE" >> $LOG
