#!/bin/bash
# Graded-box study (CPU only): sizes/memory of graded meshes and the DC port inductance versus box size and grading.
# Run as one SimQueue CPU job, e.g.
#   simq submit -n rs-graded -p powerflow -c 1 --mem 6G -t 90 --borrow -C <this dir> -- graded/run_graded_study.sh b11
set -u
L=${1:-b11}
cd /home/treefor/projects/runtimes/resolution-series-prep-20261031
R="/usr/bin/python -I build_run.py run $L"
$R graded-sizes field_scaling graded-sizes {OUT} 2e-4,1e-4,5e-5,2.5e-5 .0045,.0175,.0225 1.1,1.2 2e-3
for b in .0063 .009 .0135 .0173 .0175 .0225 .0227 .030 .045; do
  $R graded-dc-200-q1.1-b$b field_scaling graded-dc {OUT} {TRACE} 2e-4 .0045 $b 1.1 2e-3
done
for b in .0175 .0225; do
  $R graded-dc-200-q1.2-b$b field_scaling graded-dc {OUT} {TRACE} 2e-4 .0045 $b 1.2 2e-3
  $R graded-dc-200-q1.0-b$b field_scaling graded-dc {OUT} {TRACE} 2e-4 .0045 $b 1.0 2e-4
done
$R graded-dc-200-q1.05-b.0225 field_scaling graded-dc {OUT} {TRACE} 2e-4 .0045 .0225 1.05 2e-3
$R graded-dc-200-q1.1-h1mm-b.0225 field_scaling graded-dc {OUT} {TRACE} 2e-4 .0045 .0225 1.1 1e-3
$R graded-dc-100-uniform-b.0045 field_scaling graded-dc {OUT} {TRACE} 1e-4 .0045 .0045 1.1 2e-3
for b in .0175 .0225; do
  $R graded-dc-100-q1.1-b$b field_scaling graded-dc {OUT} {TRACE} 1e-4 .0045 $b 1.1 2e-3
done
$R graded-dc-50-q1.1-b.0225 field_scaling graded-dc {OUT} {TRACE} 5e-5 .0045 .0225 1.1 2e-3
$R admit-225x4-w45 field_scaling admit {OUT} {TRACE} 225 4 1 .045
$R admit-175x4-w35 field_scaling admit {OUT} {TRACE} 175 4 1 .035
