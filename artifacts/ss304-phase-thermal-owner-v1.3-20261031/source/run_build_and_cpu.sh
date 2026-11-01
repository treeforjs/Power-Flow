#!/bin/bash
set -euo pipefail
label=${1:?build label}
/usr/bin/python -I run_owner_v1_3.py build "$label"
/bin/bash run_cpu_regressions.sh "$label"
