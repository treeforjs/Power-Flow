#!/bin/bash
# Coordinator v8 stage-1 GPU work, run inside one SimQueue job (1 core, shared GPU): the job unit holds the core and GPU
# lease files; pf-launch (through pf-coupled-v8) uses the job's allocation.
#   1. pre-run checkpoint (coordinator Model) vs the Hall v3.3 20 ns checkpoint, word for word
#   2. owner-type parity (host Owner vs ResidentOwner)
#   3. coupled runs started from the 20 ns pre-run checkpoint (pre-run then genesis) through the launch interface
#   4. restart tests on the resident Hall owner
set -u
LABEL=${1:-s3}
D=/home/treefor/projects/runtimes/pic-water-common-commit-v8-20261031
O=/home/treefor/projects/.workspace/outputs/pic-water-common-commit-v8-20261031
RF=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
REF=/home/treefor/projects/.workspace/outputs/hall-corrective-v3.3-20261031/s1/gpu-prerun-20ns-130711/checkpoint-20.000000ns
OURS=$O/prerun-20ns/checkpoint-20.000000ns
cd $D
echo "== 1 hallcompare $(date +%T)"
nice -n 15 $RF /usr/bin/python -I run_v8.py run $LABEL hallcompare hallcompare-20ns --hall-owner resident --field gpu --hall-dt-max 1e-12 \
  --source-dt-max 1e-12 --dt-cap 1e-12 --hall-checkpoint $OURS --compare-hall-checkpoint $REF | tail -n 1
echo "== 2 owner parity $(date +%T)"
tools/gpu_owner_parity.sh $LABEL op
echo "== 3 start from the 20 ns checkpoint $(date +%T)"
tools/pf-coupled-v8 --cores 1 --gpus 0 --label $LABEL --out start20-pilot --hall-checkpoint $OURS --dt-cap 1e-12 --fixture pilot --window 5e-12 \
  --checkpoint-every 100 -- --progress-every 50 --frame-every 100 | tail -n 3
tools/pf-coupled-v8 --cores 1 --gpus 0 --label $LABEL --out start20-qual --hall-checkpoint $OURS --dt-cap 1e-12 --fixture qual --steps 50 \
  --checkpoint-every 25 -- --progress-every 10 --frame-every 25 | tail -n 3
echo "== 4 resident restart tests $(date +%T)"
python3 -I tools/restart_tests.py $LABEL rs3r --resident --lane nice -n 15 $RF /usr/bin/python -I | tail -n 30
echo "== done $(date +%T)"
