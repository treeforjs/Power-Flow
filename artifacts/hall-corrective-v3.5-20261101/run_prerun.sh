#!/bin/bash
# Hall v3.4 field-only (vacuum) runs inside a SimQueue powerflow job (1 core, 1 GPU; the job unit holds the core and GPU
# lease files): the divB rule divb-coefficient-gross-v1 over the 20 ns pre-run on 200 um (100 fs, checkpoints every 1 ns),
# and the longest affordable windows on 100 um and 50 um. Model: resident backend, dt upper bound 1e-13, native v7 motion
# forms. Every run goes through run_hall_v34.py (frozen launcher, receipts). Usage: run_prerun.sh LABEL prerun200|window100|window50
set -eu
LABEL=${1:?label}; MODE=${2:?mode}
D=/home/treefor/projects/runtimes/hall-corrective-v3.4-20261031
O=/home/treefor/projects/.workspace/outputs/hall-corrective-v3.4-20261031/$LABEL
cd $D
G="nice -n 15 /usr/bin/python -I run_hall_v34.py gpu $LABEL"
case $MODE in
prerun200) $G prerun-20ns-200um test_v32 prerun {TRACE} {OUT} 1e-13 2e-8 1e-11 1e-9 1e-13 {NVRTC} {UUID} 45 ;;
window100) $G window-100um test_v32 prerun {TRACE} {OUT} 1e-13 5e-9 1e-11 0 1e-13 {NVRTC} {UUID} 90 ;;
window50)  $G window-50um test_v32 prerun {TRACE} {OUT} 1e-13 3e-9 1e-11 0 1e-13 {NVRTC} {UUID} 180 ;;
*) echo "mode prerun200|window100|window50"; exit 2 ;;
esac
