#!/usr/bin/env python3
"""Print per-step ledger columns (median and last value) of a coordinator v6 run directory."""
import json, sys, struct
from pathlib import Path

def load(d):
    d = Path(d)
    cols = json.loads((d / 'ledger-schema.json').read_text())['columns']
    raw = (d / 'ledger.bin').read_bytes()
    n = len(cols)
    rows = [struct.unpack_from('<%dd' % n, raw, k * 8 * n) for k in range(len(raw) // (8 * n))]
    return cols, rows

if __name__ == '__main__':
    cols, rows = load(sys.argv[1])
    pat = sys.argv[2] if len(sys.argv) > 2 else ''
    for i, c in enumerate(cols):
        if pat and pat not in c:
            continue
        v = sorted(r[i] for r in rows)
        print('%-48s median %-14.6g last %-14.6g max %-14.6g' % (c, v[len(v) // 2], rows[-1][i], v[-1]))
