#!/usr/bin/env python3
"""Submit a short audit check through SimQueue; no local resource bypass.

Usage: run_checks.py cpu|gpu-build|gpu|sanitized
Run cpu and gpu-build successfully before gpu. Inspect simq show/log receipts.
"""
import argparse
import subprocess
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('stage', choices=['cpu', 'gpu-build', 'gpu', 'sanitized'])
p.add_argument('--borrow', action='store_true', help='optional preemptible use of idle lending cores')
a = p.parse_args()
r = Path(__file__).resolve().parent
scripts = {'cpu': 'audit_cpu.sh', 'gpu-build': 'audit_gpu_build.sh',
           'gpu': 'audit_gpu.sh', 'sanitized': 'audit_sanitized.sh'}
cmd = ['simq', 'submit', '-p', 'powerflow', '-c', '1', '--mem', '3G', '-t', '00:05:00' if a.stage == 'gpu' else '00:20:00',
       '-n', 'chemistry-v311-' + a.stage, '-C', str(r)]
if a.stage == 'gpu':
    cmd += ['--gpus', '1', '--gpu-mem', '3G']
if a.borrow:
    cmd += ['--borrow-from', 'raikou-mfem', '--on-preempt', 'freeze']
cmd += ['--', 'bash', 'tools/' + scripts[a.stage]]
raise SystemExit(subprocess.run(cmd).returncode)
