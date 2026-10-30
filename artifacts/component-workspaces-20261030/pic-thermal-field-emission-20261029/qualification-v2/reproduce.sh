#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
reports=${1:?Pass report directory}
probes=${2:?Pass original probes.csv}
snapshot=${3:-"$root/../snapshots/usable-api-v1"}
if [ ! -f "$snapshot/include/surface_emission.hpp" ] && [ -f "$root/../include/surface_emission.hpp" ]; then
  snapshot="$root/.."
fi
if [ ! -f "$snapshot/lib/libsurface_emission.a" ]; then
  cmake -S "$snapshot" -B "$snapshot/build" -DCMAKE_BUILD_TYPE=Release
  cmake --build "$snapshot/build" --target surface_emission --parallel 1
fi
cmake -S "$root" -B "$root/build" -DCMAKE_BUILD_TYPE=Release -Dsnapshot="$snapshot"
cmake --build "$root/build" --parallel 1
"$root/build/integration_qualification" "$reports" "$probes"
"$root/build/canonical_barrier_diagnostic" "$reports"
