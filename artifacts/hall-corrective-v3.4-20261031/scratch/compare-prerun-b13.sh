#!/bin/bash
# b13 20 ns pre-run against the sealed v3.3 s1 pre-run and against the b4 pre-run (tools/compare_prerun_v33.py).
set -u
R=/home/treefor/projects/runtimes/hall-corrective-v3.4-20261031
O=/home/treefor/projects/.workspace/outputs/hall-corrective-v3.4-20261031
N=$O/b13/gpu-prerun-20ns-200um-212515
cd /home/treefor/projects
/usr/bin/python -I $R/tools/compare_prerun_v33.py /home/treefor/projects/.workspace/outputs/hall-corrective-v3.3-20261031/s1/gpu-prerun-20ns-130711 $N $N/v33-comparison.json; echo "rc v33 $?"
/usr/bin/python -I $R/tools/compare_prerun_v33.py $O/b4/gpu-prerun-20ns-200um-150433 $N $N/b4-comparison.json; echo "rc b4 $?"
