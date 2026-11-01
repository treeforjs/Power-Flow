#!/bin/bash
# Run inside a one-core, frozen-runtime SimQueue allocation.
set -eu
label=${1:?build label required}
binding=${2:?corrected or parity required}
case "$binding" in
  corrected) hall=v3.3; native=v8-auditfix ;;
  parity) hall=v3.2; native=v7 ;;
  *) echo 'binding must be corrected or parity' >&2; exit 2 ;;
esac
/usr/bin/python -I run_v8.py build "$label" --hall "$hall" --native "$native"
/usr/bin/python -I run_v8.py run "$label" checks "$label-checks"
# Both bindings deliberately retain the v5 native configuration: genesis
# restores directly, while nonempty coordinator checkpoints require replay.
/usr/bin/python -I tools/audit_regressions.py "$label" "$label-continuation" replay
/usr/bin/python -I tools/audit_run_integrity.py "/home/treefor/projects/.workspace/outputs/pic-water-common-commit-v8.0.1-20261031/$label-run-integrity-results.json"
