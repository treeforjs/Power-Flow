#!/bin/bash
# Owner-type parity (host Owner vs ResidentOwner) on the same binary: CPU motion, sealed CPU PCG on both sides
# (the resident owner's host_cpu_pcg parity solve), fixtures qual (20 steps), dense adaptive (30 steps) and pilot
# (200 steps). The host-Owner runs need no GPU; the resident runs need the GPU locks (held by the SimQueue job).
# Usage (inside a SimQueue job or with the GPU locks and a core lease held): gpu_owner_parity.sh LABEL TAG
set -u
LABEL=$1; TAG=$2
D=/home/treefor/projects/runtimes/pic-water-common-commit-v8-20261031
RF=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
cd $D
run() { nice -n 15 $RF /usr/bin/python -I run_v8.py run $LABEL run "$@" | tail -n 1; }
for owner in owner resident; do
  run $TAG-$owner-qual --hall-owner $owner --fixture qual --steps 20 --frame-every 10
  run $TAG-$owner-dense --hall-owner $owner --fixture dense --blob-weight 8e10 --dt-mode adaptive --steps 30 --frame-every 10
  run $TAG-$owner-pilot --hall-owner $owner --fixture pilot --steps 200 --frame-every 100
done
# resident owner with the device CG and the GPU push (the production configuration), for the record
run $TAG-resident-gpu-qual --hall-owner resident --field gpu --push gpu --parity-every 1 --fixture qual --steps 20 --frame-every 10
