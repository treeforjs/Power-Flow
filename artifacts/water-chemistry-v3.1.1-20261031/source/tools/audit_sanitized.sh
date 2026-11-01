#!/bin/bash
set -euo pipefail
: "${SIMQ_JOB_ID:?submit this check through SimQueue}"
cd /home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031
frozen=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
"$frozen" tools/dev_build.sh sanitized tests/test_audit_regressions.cpp build/bin/audit_regressions_sanitized
"$frozen" env ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
 build/bin/audit_regressions_sanitized data/tables/chemistry-v3.1.blob | tee evidence/audit-sanitized.log
"$frozen" tools/dev_build.sh sanitized tests/test_core.cpp build/bin/test_core_sanitized
"$frozen" env ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
 build/bin/test_core_sanitized data/tables/chemistry-v3.1.blob data/tables/chemistry-v3.1-checks.json \
 evidence/core-sanitized 10000 | tee evidence/core-sanitized.log
