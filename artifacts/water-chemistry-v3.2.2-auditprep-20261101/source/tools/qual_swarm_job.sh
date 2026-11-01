#!/bin/sh
# SimQueue job: the swarm stage of a qualification run (the other stages run under the core leases first).
# usage (submitted once): simq submit -p powerflow -c 1 --borrow --max-borrow-cores 4 --mem 4G -t 3:00:00 \
#            -n pfchem-v32-swarm -C <component> -- sh tools/qual_swarm_job.sh LABEL
# The queue holds the job's CPUs, so the swarm job runner uses them directly (PFC_LEASE=none, no lease files).
set -eu
cd /home/treefor/projects/runtimes/pic-water-chemistry-v3.2.2-auditprep-20261101
export PFC_LEASE=none
exec python3 -I run_checks.py "$1" swarm inputs   # inputs again: the final component file hashes
