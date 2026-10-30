#!/bin/bash
# Serial post-reboot re-runs of the corrective component suites on CPU7.
# Every step writes to a new directory or backs up the fixed evidence directory first.
set -u
OUT=/home/treefor/sim-workspace/powerflow/outputs
RT=/home/treefor/projects/runtimes
LOCK=$OUT/pic-profile/testing-core7.lock
FROZEN=$RT/arch-base-20261028/run-frozen
TAG=post-reboot-rerun-001
SUMMARY=$1
ENV1="env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES="
: > "$SUMMARY"
cd /home/treefor/projects

step() {  # name, log, command...
    local name=$1 log=$2; shift 2
    local t0=$(date +%s)
    "$@" > "$log" 2>&1
    local rc=$?
    printf '%s\trc=%s\t%ss\t%s\n' "$name" "$rc" "$(( $(date +%s) - t0 ))" "$(tail -n 2 "$log" | tr '\n' ' ' | cut -c1-240)" >> "$SUMMARY"
}

L=$OUT/post-reboot-rerun-20261030
mkdir -p "$L"

# 1-2. Common coordinator v3 (runners take the CPU7 lease themselves; each makes a new numbered directory)
step v3-declaration "$L/v3-declaration.log" taskset -c 7 nice -n 15 python3 $RT/pic-water-common-commit-v3-20261030/run_declaration_checks.py
step v3-parent-regression "$L/v3-parent-regression.log" taskset -c 7 nice -n 15 python3 $RT/pic-water-common-commit-v3-20261030/run_parent_regression.py

# 3. Source v4: confirm the build is current, then run normal and sanitized binaries into new directories
B=$RT/pic-electrode-heat-feedback-20261029/build-heat-source-v4
V4=$OUT/pic-electrode-heat-source-v4-20261030
step v4-build "$L/v4-build.log" flock "$LOCK" taskset -c 7 nice -n 15 $ENV1 "$FROZEN" cmake --build "$B" --target electrode-heat-source-v4-tests electrode-heat-source-v4-sanitized -j1
step v4-checks "$L/v4-checks.log" flock "$LOCK" taskset -c 7 nice -n 15 $ENV1 "$FROZEN" "$B/electrode-heat-source-v4-tests" "$V4/checks-$TAG"
step v4-sanitized "$L/v4-sanitized.log" flock "$LOCK" taskset -c 7 nice -n 15 $ENV1 ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 "$FROZEN" "$B/electrode-heat-source-v4-sanitized" "$V4/checks-sanitized-$TAG"

# 4. Native candidate003: back up the fixed evidence directory, then run verify-final under an external lease
E=$OUT/pic-profile/reaction-transactions-20261029/native-full-time-candidate003
C=$RT/pic-reaction-transactions-20261029/native-full-time-candidate003
cp -a "$E" "$L/native-candidate003-evidence-before-rerun"
cp -a "$C/build" "$L/native-candidate003-build-before-rerun"
step native003-verify-final "$L/native003-verify-final.log" flock "$LOCK" taskset -c 7 nice -n 15 $ENV1 python3 $C/verify-final.py
cp -a "$E" "$L/native-candidate003-evidence-after-rerun"
sha256sum "$C"/build/*.a > "$L/native-candidate003-rebuilt-libraries.sha256" 2>/dev/null
rm -rf "$E" && cp -a "$L/native-candidate003-evidence-before-rerun" "$E"
rm -rf "$C/build" && cp -a "$L/native-candidate003-build-before-rerun" "$C/build"

# 5. Streaming primitive (fixed output: back up and restore) and wrapper (new output directory)
SP=$OUT/pic-streaming-chemistry-evidence-20261030
cp -a "$SP" "$L/streaming-primitive-evidence-before-rerun"
step streaming-primitive "$L/streaming-primitive.log" $RT/pic-streaming-chemistry-evidence-20261030/run_checks.sh
cp -a "$SP" "$L/streaming-primitive-evidence-after-rerun"
rm -rf "$SP" && cp -a "$L/streaming-primitive-evidence-before-rerun" "$SP"
step streaming-wrapper "$L/streaming-wrapper.log" env PIC_STREAMING_WRAPPER_OUTPUT=$OUT/pic-generic-water-streaming-wrapper-20261030/$TAG $RT/pic-generic-water-streaming-wrapper-20261030/run_checks.sh

# 6. Dispatcher injected primitives
step dispatcher "$L/dispatcher.log" flock "$LOCK" taskset -c 7 nice -n 15 $ENV1 "$FROZEN" /usr/bin/python3 $RT/pic-test-dispatch-v1-20261030/tests/test_dispatch.py "$L/dispatcher-primitive-results.json"

echo done >> "$SUMMARY"
