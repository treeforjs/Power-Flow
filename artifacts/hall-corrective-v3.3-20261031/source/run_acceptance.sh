#!/bin/bash
# Hall corrective v3.2 acceptance batch. Builds and CPU-only stages take loaned cores 0, 1, 3 (flock
# testing-core<K>.lock, taskset, nice 15; not started once pic-profile/RETURN-CORES-0-3 exists) or cores 6/7
# (core-lease any-core.sh / both-cores.sh); GPU stages hold both GPU lock files and one of cores 6/7
# (any-core.sh --wait). Every compile and run executes inside arch-base-20261028/run-frozen (run_hall_v33.py,
# run_coordinator_parity.py). Usage: run_acceptance.sh LABEL
set -u
LABEL=${1:?label}
D=/home/treefor/projects/runtimes/hall-corrective-v3.3-20261031
L=/home/treefor/projects/.workspace/outputs/pic-profile
CL=/home/treefor/projects/runtimes/core-lease-20261030
O=/home/treefor/projects/.workspace/outputs/hall-corrective-v3.3-20261031
LOG=$O/acceptance-$LABEL.log
cd $D
step() { echo "$(date +%H:%M:%S) START $*" >> $LOG; "$@" >> $LOG 2>&1; echo "$(date +%H:%M:%S) END rc=$? $*" >> $LOG; }
# one loaned core (0, 1, 3) when free, else cores 6/7 through any-core.sh
loan1() {
  while :; do
    if [ ! -e $L/RETURN-CORES-0-3 ]; then
      for k in 0 1 3; do
        exec {fd}>>"$L/testing-core$k.lock"
        if flock -n "$fd"; then taskset -c $k nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "$@"; rc=$?; exec {fd}>&-; return $rc; fi
        exec {fd}>&-
      done
    fi
    exec {fd}>>"$L/testing-core.lock"; if flock -n "$fd"; then exec {fd}>&-; $CL/any-core.sh --wait "$@"; return $?; fi; exec {fd}>&-
    sleep 5
  done
}
build() { step flock $L/testing-core0.lock flock $L/testing-core1.lock taskset -c 0,1 nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
            /usr/bin/python -I run_hall_v33.py build $LABEL test_hall_v3 test_resident test_v32; }
cpu() { step loan1 /usr/bin/python -I run_hall_v33.py run "$@"; }
gpu() { step flock $L/testing-gpu.lock flock $L/testing-gpu-GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9.lock $CL/any-core.sh --wait /usr/bin/python -I run_hall_v33.py gpu "$@"; }

[ "${SKIP_BUILD:-0}" = 1 ] || build
[ "${SKIP_BUILD:-0}" = 1 ] || step flock $L/testing-core0.lock flock $L/testing-core1.lock taskset -c 0,1 nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 /usr/bin/python -I run_coordinator_parity.py build $LABEL
# device-resident qualification (v3.1 tests with the ledger in the state comparison)
gpu $LABEL units test_resident units {TRACE} {OUT} {NVRTC} {UUID}
gpu $LABEL negative test_resident negative {TRACE} {OUT} {NVRTC} {UUID}
gpu $LABEL parity-vacuum test_resident parity-vacuum {TRACE} {OUT} 200 {NVRTC} {UUID}
gpu $LABEL parity-coupled-2 test_resident parity-coupled {TRACE} {OUT} 50 2 {NVRTC} {UUID}
gpu $LABEL parity-coupled-1e4 test_resident parity-coupled {TRACE} {OUT} 3 10000 {NVRTC} {UUID}
gpu $LABEL vacuum10k test_resident vacuum-run {TRACE} {OUT} 10000 {NVRTC} {UUID} sample 500
# v3.2: admitted Gauss ledger, A09, injection, tamper, checkpoint and continuation
gpu $LABEL a09 test_v32 a09 {TRACE} {OUT} {NVRTC} {UUID}
gpu $LABEL ledger-run test_v32 ledger-run {TRACE} {OUT} 55000 250 500 {NVRTC} {UUID}
LR=$(ls -d $O/$LABEL/gpu-ledger-run-* | tail -1)
gpu $LABEL restore-continue test_v32 restore-continue {TRACE} $LR/checkpoint-end $LR/continuation-reference.json {OUT} {NVRTC} {UUID}
# v3.2: motion-work-gross-v1 on the coordinator v7 dense fixtures (v3 forms and native v7 forms; CPU and GPU push)
gpu $LABEL dense-light test_v32 dense {TRACE} {OUT} 3.45e10 3 1e-14 400 gpu {NVRTC} {UUID}
gpu $LABEL dense-medium test_v32 dense {TRACE} {OUT} 8e10 3 1e-14 400 gpu {NVRTC} {UUID}
gpu $LABEL blob10 test_v32 dense {TRACE} {OUT} 1.6e11 10 5.4e-15 150 gpu {NVRTC} {UUID}
gpu $LABEL v7-sinks-gpu test_v32 dense {TRACE} {OUT} 3.45e10 3 1e-14 150 gpu {NVRTC} {UUID} v7 sinks
gpu $LABEL v7-sinks-cpu test_v32 dense {TRACE} {OUT} 3.45e10 3 1e-14 150 cpu {NVRTC} {UUID} v7 sinks
gpu $LABEL v6-sinks-gpu test_v32 dense {TRACE} {OUT} 3.45e10 3 1e-14 150 gpu {NVRTC} {UUID} v6 sinks
gpu $LABEL v6-sinks-cpu test_v32 dense {TRACE} {OUT} 3.45e10 3 1e-14 150 cpu {NVRTC} {UUID} v6 sinks
# Hall v3 component tests on the v3.2 sources (CPU backends) and coordinator v5 T1-T7
cpu $LABEL components test_hall_v3 components {TRACE} {OUT}/components
cpu $LABEL parity-dump test_hall_v3 parity-dump {TRACE} {OUT}/parity-v3.json
cpu $LABEL motion-dump test_hall_v3 motion-dump {TRACE} {OUT}/motion-dump-v3.json
cpu $LABEL sink-dump test_hall_v3 sink-dump {TRACE} {OUT}/sink-dump-cpu.json
step /usr/bin/python -I tests/compare_v3_evidence.py $O/$LABEL/v3-evidence-comparison.json $(ls -d $O/$LABEL/cpu-components-* | tail -1)/components \
    $(ls -d $O/$LABEL/cpu-parity-dump-* | tail -1)/parity-v3.json $(ls -d $O/$LABEL/cpu-motion-dump-* | tail -1)/motion-dump-v3.json \
    $(ls -d $O/$LABEL/cpu-sink-dump-* | tail -1)/sink-dump-cpu.json
step loan1 /usr/bin/python -I run_coordinator_parity.py run $LABEL
echo "$(date +%H:%M:%S) DONE" >> $LOG
