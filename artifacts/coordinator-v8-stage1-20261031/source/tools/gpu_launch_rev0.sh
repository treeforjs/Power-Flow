#!/bin/bash
# Start-from-checkpoint path through the launch interface (pf-coupled-v8 -> pf-launch -> run_v8.py) at revision-0 Hall
# checkpoints (the hall/ part of a run's genesis checkpoint, resident owner): the started run must reproduce the
# fresh-genesis run's physics. Inside a SimQueue job (1 core, shared GPU).
set -u
LABEL=${1:-s3}
D=/home/treefor/projects/runtimes/pic-water-common-commit-v8-20261031
O=/home/treefor/projects/.workspace/outputs/pic-water-common-commit-v8-20261031
cd $D
tools/pf-coupled-v8 --cores 1 --gpus 0 --label $LABEL --out start0-pilot --hall-checkpoint $O/rs3r-pilot-ref/checkpoints/ckpt-000000000000/hall \
  --dt-cap 1e-12 --field cpu --fixture pilot --steps 60 --checkpoint-every 20 | tail -n 2
tools/pf-coupled-v8 --cores 1 --gpus 0 --label $LABEL --out start0-qual --hall-checkpoint $O/rs3r-qual-ref/checkpoints/ckpt-000000000000/hall \
  --dt-mode fixed --preset benchmark --field cpu --fixture qual --steps 30 --checkpoint-every 10 | tail -n 2
# production configuration (device CG) from the same revision-0 checkpoint; restart from its checkpoint 20 (replay)
tools/pf-coupled-v8 --cores 1 --gpus 0 --label $LABEL --out start0-pilot-gpu --hall-checkpoint $O/rs3r-pilot-ref/checkpoints/ckpt-000000000000/hall \
  --dt-cap 1e-12 --fixture pilot --steps 60 --checkpoint-every 20 | tail -n 2
