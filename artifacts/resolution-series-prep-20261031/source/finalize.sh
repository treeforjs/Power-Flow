#!/bin/bash
# Regenerate the derived results from the recorded runs: dt limits, run-matrix inputs and table, acceptance checks.
# CPU only: run on a lent core, e.g. ./lent.sh 3 ./finalize.sh
set -eu
cd /home/treefor/projects/runtimes/resolution-series-prep-20261031
O=/home/treefor/projects/.workspace/outputs/resolution-series-prep-20261031
REP=$(ls -d $O/b1/mesh-report-*/mesh-report.json | tail -1)
/usr/bin/python -I runmatrix/collect_inputs.py $O runmatrix/inputs.json > /dev/null
VAC=$(/usr/bin/python -I -c "import json;print(' '.join(json.load(open('runmatrix/inputs.json'))['gpu_vacuum']))")
SCN=$(/usr/bin/python -I -c "import json;print(' '.join(json.load(open('runmatrix/inputs.json'))['gpu_dtscan']))")
/usr/bin/python -I timestep/dt_limits.py $O/dt-limits $REP $VAC --dtscan $SCN > /dev/null
cp $O/dt-limits/dt-limits.txt timestep/DT-LIMITS.txt
/usr/bin/python -I runmatrix/run_matrix.py $O/runmatrix runmatrix/inputs.json > /dev/null
cp $O/runmatrix/run-matrix.txt runmatrix/RUN-MATRIX.txt
/usr/bin/python -I tests/check_results.py $O $O/tests-final.json || true
cp $O/tests-final.json tests/TESTS-RESULT.json
