#!/bin/bash
set -euo pipefail
: "${SIMQ_JOB_ID:?submit this check through SimQueue}"
cd /home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031
mkdir -p evidence build/bin
frozen=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
bash tools/audit_baseline.sh
"$frozen" python3 -I tests/test_source_bindings.py | tee evidence/input-bindings.log
"$frozen" python3 -I tools/build_tables.py data/tables/chemistry-v3.1-spec.json data/tables/chemistry-v3.1 > evidence/tables-build.log
"$frozen" python3 -I tools/make_coverage_policy.py data/tables/chemistry-v3.1-manifest.json data/tables/chemistry-v3.1-coverage-policy.json > evidence/policy-build.log
"$frozen" python3 -I tests/test_audit_policy.py | tee evidence/audit-policy.log
"$frozen" tools/dev_build.sh normal tests/test_audit_regressions.cpp build/bin/audit_regressions
"$frozen" build/bin/audit_regressions data/tables/chemistry-v3.1.blob | tee evidence/audit-corrected.log
"$frozen" tools/dev_build.sh normal tests/test_core.cpp build/bin/test_core
"$frozen" build/bin/test_core data/tables/chemistry-v3.1.blob data/tables/chemistry-v3.1-checks.json evidence/core-normal 100000 | tee evidence/core-normal.log
