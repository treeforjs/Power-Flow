#!/bin/bash
set -euo pipefail
label=${1:?build label}
for name in baseline_defects transactions law owner v11 v12 regression; do
  /usr/bin/python -I run_owner_v1_3.py run "$label" "$name" "test_$name" '{OUT}'
done
