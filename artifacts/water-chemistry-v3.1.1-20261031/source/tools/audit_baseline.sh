#!/bin/bash
set -euo pipefail
: "${SIMQ_JOB_ID:?submit this check through SimQueue}"
cd /home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031
mkdir -p evidence build/bin
/home/treefor/projects/runtimes/arch-base-20261028/run-frozen \
 /usr/bin/c++ -std=c++20 -O2 -Wall -Wextra -Werror -ffp-contract=off -fno-fast-math \
 -I /home/treefor/projects/runtimes/pic-water-chemistry-v3.1-20261031/include/pf_chem_v3 \
 -I /home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include/pf_arrival \
 tests/test_audit_regressions.cpp -lcrypto -o build/bin/audit_baseline
set +e
/home/treefor/projects/runtimes/arch-base-20261028/run-frozen build/bin/audit_baseline data/tables/chemistry-v3.1.blob > evidence/audit-baseline.log 2>&1
rc=$?
set -e
cat evidence/audit-baseline.log
printf '%s\n' "$rc" > evidence/audit-baseline.exit
test "$rc" = 1
