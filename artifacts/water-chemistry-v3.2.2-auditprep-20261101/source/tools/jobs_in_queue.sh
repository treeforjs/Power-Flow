#!/bin/sh
# SimQueue job: run a tools/run_jobs.py job list on the job's CPUs (the queue holds them; no lease files).
# usage (submitted once): simq submit -p powerflow -c 1 --borrow --max-borrow-cores N --mem 4G -t 2:00:00 \
#            -n NAME -C <component> -- sh tools/jobs_in_queue.sh JOBS.json RECEIPT.json
set -eu
cd /home/treefor/projects/runtimes/pic-water-chemistry-v3.2.2-auditprep-20261101
export PFC_LEASE=none
exec python3 -I tools/run_jobs.py "$1" "$2"
