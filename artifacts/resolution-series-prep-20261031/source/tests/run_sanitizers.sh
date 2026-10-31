#!/bin/bash
# AddressSanitizer + UndefinedBehaviorSanitizer runs of the CPU paths (build s1 = build_run.py build s1 --sanitize):
# mesh construction/report, Hall admission with the instrumented DC seed, the DC accuracy diagnostic and the GMG CPU reference
# with the Hall independent residual check. CPU only (lent core), no GPU. Sanitizer defaults (halt on error, leak check at exit).
set -u
cd /home/treefor/projects/runtimes/resolution-series-prep-20261031
K=${1:-1}
R="./lent.sh $K /usr/bin/python -I build_run.py run s1"
$R san-mesh-report field_scaling mesh-report {OUT} 45 4 1 90 4 1
$R san-admit-45 field_scaling admit {OUT} {TRACE} 45 4 1
$R san-dc-floor-45 field_scaling dc-floor {OUT} {TRACE} 45 4 1 1
$R san-gmg-45-random gmg_solve {OUT} {TRACE} none none 45 4 1 random 0,1,4 1e-14 1e-12
$R san-gmg-90-physical gmg_solve {OUT} {TRACE} none none 90 4 1 physical 0,1,4 1e-13 1e-12
