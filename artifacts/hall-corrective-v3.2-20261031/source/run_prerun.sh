#!/bin/bash
# Hall v3.2 field-only (vacuum) runs inside a SimQueue powerflow job (1 core, 1 GPU; the job unit holds the core and GPU
# lease files): the 10 fs / 100 fs cross-check from the supported start to 2.75 ns, or the pre-run to 20 ns at 100 fs
# with checkpoints every 1 ns and at 20 ns. Model: resident backend, dt upper bound 1e-13, native v7 motion forms.
# Every run goes through run_hall_v32.py (frozen launcher, receipts). Usage: run_prerun.sh LABEL crosscheck|prerun
set -eu
LABEL=${1:?label}; MODE=${2:?mode}
D=/home/treefor/projects/runtimes/hall-corrective-v3.2-20261031
O=/home/treefor/projects/.workspace/outputs/hall-corrective-v3.2-20261031/$LABEL
cd $D
G="nice -n 15 /usr/bin/python -I run_hall_v32.py gpu $LABEL"
last() { ls -d $O/gpu-$1-* | tail -1; }
case $MODE in
crosscheck)
  $G xc-100fs test_v32 prerun {TRACE} {OUT} 1e-13 2.75e-9 1e-11 0 1e-13 {NVRTC} {UUID}
  $G xc-10fs test_v32 prerun {TRACE} {OUT} 1e-14 2.75e-9 1e-11 0 1e-13 {NVRTC} {UUID}
  A=$(last xc-100fs); B=$(last xc-10fs)
  $G xc-compare test_v32 compare {TRACE} $A/checkpoint-2.750000ns $B/checkpoint-2.750000ns {OUT}
  C=$(last xc-compare)
  nice -n 15 /usr/bin/python -I tests/compare_series.py $A/series.jsonl $B/series.jsonl $C/compare.json $O/crosscheck.json
  nice -n 15 /usr/bin/python -I tests/plot_series.py $A/series.jsonl $A/plots "100 fs 0.75-2.75 ns"
  nice -n 15 /usr/bin/python -I tests/plot_series.py $B/series.jsonl $B/plots "10 fs 0.75-2.75 ns" ;;
prerun)
  $G prerun-20ns test_v32 prerun {TRACE} {OUT} 1e-13 2e-8 1e-11 1e-9 1e-13 {NVRTC} {UUID}
  P=$(last prerun-20ns)
  nice -n 15 /usr/bin/python -I tests/plot_series.py $P/series.jsonl $P/plots "Hall v3.2 pre-run 100 fs" ;;
*) echo "mode crosscheck|prerun"; exit 2 ;;
esac
