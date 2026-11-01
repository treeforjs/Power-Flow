#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
reports=${1:?Pass a report directory}
probes=${2:-/home/treefor/sim-workspace/powerflow/outputs/cathode-emission-audit-20261029/fields/probes.csv}
# Caller must hold testing-core7.lock, use CPU7/nice15/one job and frozen userspace.
cmake -S "$root" -B "$root/build" -DCMAKE_BUILD_TYPE=Release
cmake --build "$root/build" --parallel 1
"$root/build/emission_qualification" "$reports" "$probes"
"$root/build/emission_micro_qualification" "$reports/micro"
