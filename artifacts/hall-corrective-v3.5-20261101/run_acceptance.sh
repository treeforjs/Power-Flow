#!/bin/bash
# Hall corrective v3.5 acceptance batch. CPU-only stages take one of cores 6/7 (core-lease any-core.sh); GPU stages run
# through SimQueue's gpu-lease helper (GPU locks, then one of cores 6/7). Every compile and run executes inside
# arch-base-20261028/run-frozen (run_hall_v35.py).
#   v3.5 qualification: restore_as_genesis (test_v35 genesis: the v3.3 and v3.4 20 ns checkpoints, transfers, lockstep
#   against the plain restore and a fresh genesis, refusals, host Owner); regression of the v3.4 suites that the v3.5
#   changes touch (identities and arithmetic of bootstraps, steps, checkpoints and restores unchanged): test_v34 parity
#   200 um, tamper, divB, block digests; the v3.3 suites in device mode; Hall v3 component tests and dumps (CPU).
# Usage: run_acceptance.sh LABEL        (the build LABEL must exist; SKIP_V35=1 / SKIP_REG=1 skip the groups)
set -u
LABEL=${1:?label}
D=/home/treefor/projects/runtimes/hall-corrective-v3.5-20261101
CL=/home/treefor/projects/runtimes/core-lease-20261030
G=/home/treefor/projects/runtimes/simqueue-0.6.2-20261031/bin/gpu-lease
O=/home/treefor/projects/.workspace/outputs/hall-corrective-v3.5-20261101
LOG=$O/acceptance-$LABEL.log
V33=/home/treefor/projects/.workspace/outputs/hall-corrective-v3.3-20261031/s1/gpu-prerun-20ns-130711/checkpoint-20.000000ns
V34=/home/treefor/projects/.workspace/outputs/hall-corrective-v3.4-20261031/b13/gpu-prerun-20ns-200um-212515/checkpoint-20.000000ns
cd $D
step() { echo "$(date +%H:%M:%S) START $*" >> $LOG; "$@" >> $LOG 2>&1; echo "$(date +%H:%M:%S) END rc=$? $*" >> $LOG; }
skip() { if [ -n "${START_AT:-}" ] && [ "$2" != "$START_AT" ]; then echo "$(date +%H:%M:%S) SKIP $2" >> $LOG; return 0; fi; START_AT=; return 1; }
cpu() { skip "$@" && return; step $CL/any-core.sh --wait /usr/bin/python -I run_hall_v35.py run "$@"; }
gpu() { skip "$@" && return; sleep ${GPU_GAP:-20}; step $G --core any -q --label "hall-v35-$2" -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /usr/bin/python -I run_hall_v35.py gpu "$@"; }

if [ "${SKIP_V35:-0}" != 1 ]; then
  gpu $LABEL genesis test_v35 genesis {TRACE} {OUT} {NVRTC} {UUID} $V33 $V34
  gpu $LABEL genesis-run test_v35 genesis-run {TRACE} {OUT} {NVRTC} {UUID} $V34
  GR=$(ls -d $O/$LABEL/gpu-genesis-run-* | tail -1)
  gpu $LABEL genesis-continue test_v35 genesis-continue {TRACE} {OUT} {NVRTC} {UUID} $GR/rebased-50 $GR/continuation-reference.json
fi
if [ "${SKIP_REG:-0}" != 1 ]; then
  cpu $LABEL x87-support test_x87 support {OUT} 60
  gpu $LABEL v34-parity-200um test_v34 parity {TRACE} {OUT} {NVRTC} {UUID} 45
  gpu $LABEL v34-tamper test_v34 tamper {TRACE} {OUT} {NVRTC} {UUID}
  gpu $LABEL v34-divb test_v34 divb {TRACE} {OUT} {NVRTC} {UUID} 2000
  gpu $LABEL v34-blocks test_v34 blocks {TRACE} {OUT} {NVRTC} {UUID}
  for MODE in device; do
    E="--env HALL_HOST_WORK=$MODE"
    gpu $LABEL $MODE-units $E test_resident units {TRACE} {OUT} {NVRTC} {UUID}
    gpu $LABEL $MODE-negative $E test_resident negative {TRACE} {OUT} {NVRTC} {UUID}
    gpu $LABEL $MODE-parity-vacuum $E test_resident parity-vacuum {TRACE} {OUT} 200 {NVRTC} {UUID}
    gpu $LABEL $MODE-parity-coupled-2 $E test_resident parity-coupled {TRACE} {OUT} 50 2 {NVRTC} {UUID}
    gpu $LABEL $MODE-vacuum10k $E test_resident vacuum-run {TRACE} {OUT} 10000 {NVRTC} {UUID} sample 500
    gpu $LABEL $MODE-a09 $E test_v32 a09 {TRACE} {OUT} {NVRTC} {UUID}
    gpu $LABEL $MODE-ledger-run $E test_v32 ledger-run {TRACE} {OUT} 55000 250 500 {NVRTC} {UUID}
    LR=$(ls -d $O/$LABEL/gpu-$MODE-ledger-run-* | tail -1)
    gpu $LABEL $MODE-restore-continue $E test_v32 restore-continue {TRACE} $LR/checkpoint-end $LR/continuation-reference.json {OUT} {NVRTC} {UUID}
    gpu $LABEL $MODE-dense-light $E test_v32 dense {TRACE} {OUT} 3.45e10 3 1e-14 400 gpu {NVRTC} {UUID}
    gpu $LABEL $MODE-v7-sinks-gpu $E test_v32 dense {TRACE} {OUT} 3.45e10 3 1e-14 150 gpu {NVRTC} {UUID} v7 sinks
  done
  cpu $LABEL components test_hall_v3 components {TRACE} {OUT}/components
  cpu $LABEL parity-dump test_hall_v3 parity-dump {TRACE} {OUT}/parity-v3.json
  cpu $LABEL motion-dump test_hall_v3 motion-dump {TRACE} {OUT}/motion-dump-v3.json
  cpu $LABEL sink-dump test_hall_v3 sink-dump {TRACE} {OUT}/sink-dump-cpu.json
fi
echo "$(date +%H:%M:%S) DONE" >> $LOG
