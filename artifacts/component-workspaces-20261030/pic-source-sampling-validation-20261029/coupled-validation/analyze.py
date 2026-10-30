#!/usr/bin/env python3
"""Analyze already-exported observable bundles; no simulation or activation path."""
import argparse
import json
import fcntl
import math
import os
from pathlib import Path
import metrics

ROOT=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/source-sampling-validation/coupled-preflight')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('bundle',type=Path);p.add_argument('output',type=Path);args=p.parse_args()
    if os.sched_getaffinity(0)!={6} or os.getpriority(os.PRIO_PROCESS,0)!=15:raise RuntimeError('CPU6/nice15 required')
    out=args.output.resolve();out.relative_to(ROOT.resolve())
    if out.exists():raise RuntimeError('Analysis output must be new')
    lock=(ROOT.parent.parent/'testing-core.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    def number(s):
        value=float(s)
        if not math.isfinite(value):raise ValueError('Nonfinite observable')
        return value
    def constant(s):raise ValueError('Nonfinite JSON constant')
    data=json.loads(args.bundle.read_text(),parse_float=number,parse_constant=constant);result=metrics.evaluate_bundle(data)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
