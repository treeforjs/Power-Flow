"""Queue bounded preparation checks; no swarm matrix or production promotion."""
import argparse
from pathlib import Path
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('stage', choices=('regenerate', 'cpu', 'sanitized', 'gpu-build', 'gpu'))
p.add_argument('--borrow', action='store_true')
a = p.parse_args()
r = Path(__file__).resolve().parent
cmd = ['simq', 'submit', '-p', 'powerflow', '-c', '1', '--mem', '3G', '-t',
       '00:05:00' if a.stage in ('regenerate', 'gpu') else '00:20:00',
       '-n', 'chem-auditprep-' + a.stage, '-C', str(r)]
if a.borrow:
    if a.stage == 'gpu': p.error('GPU parity uses an explicit GPU allocation, not borrowed CPU timing')
    cmd += ['--borrow', '--max-borrow-cores', '1', '--on-preempt', 'freeze']
if a.stage == 'gpu': cmd += ['--gpus', '1', '--gpu-mem', '3G']
cmd += ['--', '/home/treefor/projects/runtimes/arch-base-20261028/run-frozen',
        '/usr/bin/python3', '-I', 'tools/prep_stage.py', a.stage]
raise SystemExit(subprocess.run(cmd).returncode)
