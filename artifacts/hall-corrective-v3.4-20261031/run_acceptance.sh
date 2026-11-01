#!/bin/bash
# Hall corrective v3.4 acceptance batch. CPU-only stages take a loaned core 0 or 1 (flock testing-core<K>.lock, taskset,
# nice 15; not started once pic-profile/RETURN-CORES-0-3 exists) or cores 6/7 (core-lease any-core.sh); GPU stages hold
# both GPU lock files, then one of cores 6/7 (any-core.sh --wait). Every compile and run executes inside
# arch-base-20261028/run-frozen (run_hall_v34.py, run_coordinator_parity.py).
#   v3.4 qualification: x87 emulation (CPU), device/host-path parity (200 and 50 um), tamper, divB ledger, block digests,
#   host-loop timing; the v3.3 suites in both host-work modes (HALL_HOST_WORK=v33 | device); the native v7 forms dense
#   and sink runs (device path in device mode); Hall v3 component tests (CPU).
# Usage: run_acceptance.sh LABEL        (the build LABEL must exist; SKIP_V34=1 / SKIP_V33=1 skip the groups)
set -u
LABEL=${1:?label}
D=/home/treefor/projects/runtimes/hall-corrective-v3.4-20261031
L=/home/treefor/projects/.workspace/outputs/pic-profile
CL=/home/treefor/projects/runtimes/core-lease-20261030
O=/home/treefor/projects/.workspace/outputs/hall-corrective-v3.4-20261031
LOG=$O/acceptance-$LABEL.log
cd $D
step() { echo "$(date +%H:%M:%S) START $*" >> $LOG; "$@" >> $LOG 2>&1; echo "$(date +%H:%M:%S) END rc=$? $*" >> $LOG; }
loan1() {
  while :; do
    if [ ! -e $L/RETURN-CORES-0-3 ]; then
      for k in 0 1; do
        exec {fd}>>"$L/testing-core$k.lock"
        if flock -n "$fd"; then taskset -c $k nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "$@"; rc=$?; exec {fd}>&-; return $rc; fi
        exec {fd}>&-
      done
    fi
    exec {fd}>>"$L/testing-core.lock"; if flock -n "$fd"; then exec {fd}>&-; $CL/any-core.sh --wait "$@"; return $?; fi; exec {fd}>&-
    sleep 5
  done
}
# START_AT=NAME: skip the stages before NAME (resume a batch). Between GPU stages: a gap of GPU_GAP seconds (default
# 45) and, while SimQueue holds queued GPU jobs, up to 30 min of waiting so that queued jobs get the GPU.
skip() { if [ -n "${START_AT:-}" ] && [ "$2" != "$START_AT" ]; then echo "$(date +%H:%M:%S) SKIP $2" >> $LOG; return 0; fi; START_AT=; return 1; }
# The GPU locks are taken only when a Power-Flow core (6 or 7) is free, so that a stage does not hold the GPU while it
# waits for a core.
polite() {
  sleep ${GPU_GAP:-45}
  for i in $(seq 1 60); do simq list 2>/dev/null | awk '$2=="queued" && $5=="1"' | grep -q . || break; sleep 30; done
  for i in $(seq 1 240); do
    for c in testing-core.lock testing-core7.lock; do flock -n "$L/$c" true && return 0; done
    sleep 15
  done
}
cpu() { skip "$@" && return; step loan1 /usr/bin/python -I run_hall_v34.py run "$@"; }
gpu() { skip "$@" && return; polite; step flock $L/testing-gpu.lock flock $L/testing-gpu-GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9.lock $CL/any-core.sh --wait /usr/bin/python -I run_hall_v34.py gpu "$@"; }

if [ "${SKIP_V34:-0}" != 1 ]; then
  cpu $LABEL x87-ops test_x87 ops {OUT} 4000000
  cpu $LABEL x87-compare test_x87 compare {OUT} 4000000
  cpu $LABEL x87-support test_x87 support {OUT} 60
  gpu $LABEL v34-parity-200um test_v34 parity {TRACE} {OUT} {NVRTC} {UUID} 45
  gpu $LABEL v34-parity-50um test_v34 parity {TRACE} {OUT} {NVRTC} {UUID} 180
  gpu $LABEL v34-tamper test_v34 tamper {TRACE} {OUT} {NVRTC} {UUID}
  gpu $LABEL v34-divb test_v34 divb {TRACE} {OUT} {NVRTC} {UUID} 2000
  gpu $LABEL v34-blocks test_v34 blocks {TRACE} {OUT} {NVRTC} {UUID}
  gpu $LABEL v34-hostloop test_v34 hostloop {TRACE} {OUT} {NVRTC} {UUID} 10000 30
fi
if [ "${SKIP_V33:-0}" != 1 ]; then
  for MODE in v33 device; do
    E="--env HALL_HOST_WORK=$MODE"
    gpu $LABEL $MODE-units $E test_resident units {TRACE} {OUT} {NVRTC} {UUID}
    gpu $LABEL $MODE-negative $E test_resident negative {TRACE} {OUT} {NVRTC} {UUID}
    gpu $LABEL $MODE-parity-vacuum $E test_resident parity-vacuum {TRACE} {OUT} 200 {NVRTC} {UUID}
    gpu $LABEL $MODE-parity-coupled-2 $E test_resident parity-coupled {TRACE} {OUT} 50 2 {NVRTC} {UUID}
    gpu $LABEL $MODE-parity-coupled-1e4 $E test_resident parity-coupled {TRACE} {OUT} 3 10000 {NVRTC} {UUID}
    gpu $LABEL $MODE-vacuum10k $E test_resident vacuum-run {TRACE} {OUT} 10000 {NVRTC} {UUID} sample 500
    gpu $LABEL $MODE-a09 $E test_v32 a09 {TRACE} {OUT} {NVRTC} {UUID}
    gpu $LABEL $MODE-ledger-run $E test_v32 ledger-run {TRACE} {OUT} 55000 250 500 {NVRTC} {UUID}
    LR=$(ls -d $O/$LABEL/gpu-$MODE-ledger-run-* | tail -1)
    gpu $LABEL $MODE-restore-continue $E test_v32 restore-continue {TRACE} $LR/checkpoint-end $LR/continuation-reference.json {OUT} {NVRTC} {UUID}
    gpu $LABEL $MODE-dense-light $E test_v32 dense {TRACE} {OUT} 3.45e10 3 1e-14 400 gpu {NVRTC} {UUID}
    gpu $LABEL $MODE-dense-medium $E test_v32 dense {TRACE} {OUT} 8e10 3 1e-14 400 gpu {NVRTC} {UUID}
    gpu $LABEL $MODE-blob10 $E test_v32 dense {TRACE} {OUT} 1.6e11 10 5.4e-15 150 gpu {NVRTC} {UUID}
    gpu $LABEL $MODE-v7-dense-light $E test_v32 dense {TRACE} {OUT} 3.45e10 3 1e-14 400 gpu {NVRTC} {UUID} v7
    gpu $LABEL $MODE-v7-dense-medium $E test_v32 dense {TRACE} {OUT} 8e10 3 1e-14 400 gpu {NVRTC} {UUID} v7
    gpu $LABEL $MODE-v7-blob10 $E test_v32 dense {TRACE} {OUT} 1.6e11 10 5.4e-15 150 gpu {NVRTC} {UUID} v7
    gpu $LABEL $MODE-v7-sinks-gpu $E test_v32 dense {TRACE} {OUT} 3.45e10 3 1e-14 150 gpu {NVRTC} {UUID} v7 sinks
    gpu $LABEL $MODE-v7-sinks-cpu $E test_v32 dense {TRACE} {OUT} 3.45e10 3 1e-14 150 cpu {NVRTC} {UUID} v7 sinks
    gpu $LABEL $MODE-v6-sinks-gpu $E test_v32 dense {TRACE} {OUT} 3.45e10 3 1e-14 150 gpu {NVRTC} {UUID} v6 sinks
  done
  cpu $LABEL components test_hall_v3 components {TRACE} {OUT}/components
  cpu $LABEL parity-dump test_hall_v3 parity-dump {TRACE} {OUT}/parity-v3.json
  cpu $LABEL motion-dump test_hall_v3 motion-dump {TRACE} {OUT}/motion-dump-v3.json
  cpu $LABEL sink-dump test_hall_v3 sink-dump {TRACE} {OUT}/sink-dump-cpu.json
  step /usr/bin/python -I tests/compare_v3_evidence.py $O/$LABEL/v3-evidence-comparison.json $(ls -d $O/$LABEL/cpu-components-* | tail -1)/components \
      $(ls -d $O/$LABEL/cpu-parity-dump-* | tail -1)/parity-v3.json $(ls -d $O/$LABEL/cpu-motion-dump-* | tail -1)/motion-dump-v3.json \
      $(ls -d $O/$LABEL/cpu-sink-dump-* | tail -1)/sink-dump-cpu.json
fi
echo "$(date +%H:%M:%S) DONE" >> $LOG
