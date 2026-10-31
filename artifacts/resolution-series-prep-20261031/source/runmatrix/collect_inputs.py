#!/usr/bin/python -I
"""Collect the measured files used by run_matrix.py (latest run per mesh and kind) into an inputs JSON.

  collect_inputs.py OUTPUT_ROOT INPUTS.json
"""
import glob, json, sys
from pathlib import Path


def latest(paths, key):
    best = {}
    for p in sorted(paths):
        try:
            d = json.load(open(p))
        except Exception:
            continue
        if d.get('steps_done', 1) == 0 or ('rows' not in d and 'schema' in d and 'dtscan' in d['schema']):
            continue
        k = key(d)
        best[k] = p   # sorted by label/time: later wins
    return sorted(best.values())


def main(argv):
    root = Path(argv[0])
    inp = {
        'mesh_report': sorted(glob.glob(str(root / 'b1/mesh-report-*/mesh-report.json')))[-1],
        'gpu_vacuum': latest(glob.glob(str(root / '*/gpu-vacuum-*/gpu-vacuum.json')), lambda d: tuple(d['spec'][:3])),
        'gpu_dtscan': latest(glob.glob(str(root / '*/gpu-dtscan-*/gpu-dtscan.json')), lambda d: tuple(d['spec'][:3])),
        'gpu_coupled': latest(glob.glob(str(root / '*/gpu-coupled-*/gpu-coupled.json')), lambda d: (tuple(d['spec'][:3]), d['carriers'])),
        'gmg': latest(glob.glob(str(root / '*/gmg-*/gmg.json')), lambda d: (tuple(d['spec'][:3]), d['rhs'], d.get('options', {}).get('nu', 1))),
        'dt_limits': str(root / 'dt-limits/dt-limits.json'),
        'graded_sizes': (sorted(glob.glob(str(root / '*/graded-sizes-*/graded-sizes.json'))) or [None])[-1],
    }
    Path(argv[1]).write_text(json.dumps(inp, indent=1) + '\n')
    print(json.dumps(inp, indent=1))


if __name__ == '__main__':
    main(sys.argv[1:])
