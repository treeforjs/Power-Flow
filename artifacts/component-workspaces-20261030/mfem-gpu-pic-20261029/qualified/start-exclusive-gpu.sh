#!/bin/sh
set -eu
exec /usr/bin/pkexec /usr/bin/systemd-run --unit=powerflow-pic-gpu-reserved '--description=Powerflow PIC GPU with reserved physical core' --property=AllowedCPUs=0-6,8-14 --property=KillMode=mixed --property=TimeoutStopSec=90 /usr/bin/python /home/treefor/projects/runtimes/mfem-gpu-pic-20261029/qualified/reserve-core-and-run.py
