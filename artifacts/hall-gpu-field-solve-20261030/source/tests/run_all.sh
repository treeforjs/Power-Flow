#!/bin/sh
# Runs every phase of one qualification label in order, each phase under its
# own leases (core lease for CPU phases; GPU lease, then a core lease, for GPU).
set -eu
LABEL=$1; shift
T=/home/treefor/projects/runtimes/hall-gpu-resident-v1-20261030/field/tests
L=/home/treefor/projects/runtimes/core-lease-20261030
G=/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-gpu.lock
for M in "$@"; do
  $L/any-core.sh --wait /usr/bin/python -I $T/run_field.py prep "$LABEL" "$M"
  flock "$G" $L/any-core.sh --wait /usr/bin/python -I $T/run_field.py gpu "$LABEL" "$M"
  $L/any-core.sh --wait /usr/bin/python -I $T/run_field.py guard "$LABEL" "$M"
done
$L/any-core.sh --wait /usr/bin/python -I $T/run_field.py report "$LABEL"
