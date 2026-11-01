#!/bin/bash
set -euo pipefail
: "${SIMQ_JOB_ID:?submit this check through SimQueue}"
cd /home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031
frozen=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
"$frozen" python3 -I adoption/v3.2/verify_base.py adoption/v3.2/base | tee evidence/adoption-base-check.log
"$frozen" /usr/bin/c++ -std=c++20 -Wall -Wextra -Werror -ffp-contract=off -fno-fast-math -fsyntax-only \
 -I adoption/v3.2/patched/include/pf_chem_v3 \
 -I /home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include/pf_arrival \
 tests/test_audit_regressions.cpp
echo 'PASS patched v3.2 headers compile with the audit regression interface (syntax only)' | tee evidence/adoption-syntax.log
