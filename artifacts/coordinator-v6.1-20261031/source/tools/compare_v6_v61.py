#!/usr/bin/python -I
"""Coordinator v6 -> v6.1 comparisons (post-processing only; not part of the binary).

  identity A B [--steps N] [--out JSON]
      A: a v6 run directory (reference), B: a v6.1 run directory. Byte identity of genesis.json,
      ledger-schema.json, the first N rows of deterministic.bin and the first N lines of
      step-digests.txt (N = the smaller step count, or --steps); every deterministic ledger
      column of ledger.bin compared bit for bit by name; frames.bin frames present in both.
  timings DIR [DIR...] [--out JSON]
      Median and mean of every t_* ledger column (ms) and of the v6.1 sub-phase columns.
  checks A B [--out JSON]
      Case-by-case pass/fail of two checks-result.json files.
"""
import array, hashlib, json, statistics, sys
from pathlib import Path


def sha(b):
    return hashlib.sha256(b).hexdigest()


def ledger(d):
    cols = json.loads((d / 'ledger-schema.json').read_text())['columns']
    raw = (d / 'ledger.bin').read_bytes()
    n = len(cols)
    rows = len(raw) // (8 * n)
    a = array.array('d')
    a.frombytes(raw[:rows * n * 8])
    return cols, rows, a


def deterministic(name):
    return not name.startswith('t_') and name not in ('alloc_swap_window', 'rss_kib')


def identity(A, B, steps=None):
    A, B = Path(A), Path(B)
    ca, ra, la = ledger(A)
    cb, rb, lb = ledger(B)
    n = min(ra, rb) if steps is None else steps
    if ra < n or rb < n:
        raise SystemExit('fewer rows than requested: %d %d < %d' % (ra, rb, n))
    out = {'mode': 'identity', 'reference_v6': str(A), 'candidate_v61': str(B), 'steps_compared': n, 'rows': [ra, rb], 'files': {}}
    for f in ('genesis.json', 'ledger-schema.json'):
        x, y = (A / f).read_bytes(), (B / f).read_bytes()
        out['files'][f] = {'v6_sha256': sha(x), 'v61_sha256': sha(y), 'identical': x == y}
    ndet = sum(1 for c in ca if deterministic(c))
    x = (A / 'deterministic.bin').read_bytes()[:n * ndet * 8]
    y = (B / 'deterministic.bin').read_bytes()[:n * ndet * 8]
    out['files']['deterministic.bin'] = {'rows': n, 'bytes': len(x), 'v6_sha256': sha(x), 'v61_sha256': sha(y), 'identical': x == y and len(x) == n * ndet * 8}
    da = (A / 'step-digests.txt').read_text().splitlines()[:n]
    db = (B / 'step-digests.txt').read_text().splitlines()[:n]
    out['files']['step-digests.txt'] = {'lines': n, 'v6_sha256': sha('\n'.join(da).encode()), 'v61_sha256': sha('\n'.join(db).encode()),
                                        'identical': da == db and len(da) == n,
                                        'first_difference': next((i + 1 for i, (p, q) in enumerate(zip(da, db)) if p != q), None),
                                        'last_line': db[-1] if db else ''}
    differing = []
    for c in ca:
        if not deterministic(c):
            continue
        if c not in cb:
            differing.append(c + ' (missing)')
            continue
        i, j = ca.index(c), cb.index(c)
        for r in range(n):
            u, v = la[r * len(ca) + i], lb[r * len(cb) + j]
            if array.array('d', [u]).tobytes() != array.array('d', [v]).tobytes():
                differing.append(c)
                break
    out['ledger_deterministic_columns'] = sum(1 for c in ca if deterministic(c))
    out['ledger_columns_differing'] = differing
    fa, fb = A / 'frames.bin', B / 'frames.bin'
    if fa.exists() and fb.exists():
        lay = json.loads((A / 'frames-layout.json').read_text())
        nx, ny = lay['nx'], lay['ny']
        size = 8 * (2 + 6 * nx * ny + (nx + 1) * (ny + 1))
        ba, bb = fa.read_bytes(), fb.read_bytes()
        k = min(len(ba), len(bb)) // size
        same = [ba[i * size:(i + 1) * size] == bb[i * size:(i + 1) * size] for i in range(k)]
        out['files']['frames.bin'] = {'frames_compared': k, 'identical_frames': sum(same),
                                      'layout_identical': (A / 'frames-layout.json').read_bytes() == (B / 'frames-layout.json').read_bytes()}
    out['bit_identical'] = all(v['identical'] for k, v in out['files'].items() if 'identical' in v) and not differing and \
        out['files'].get('frames.bin', {}).get('identical_frames', 0) == out['files'].get('frames.bin', {}).get('frames_compared', 0)
    return out


def timings(dirs):
    out = {'mode': 'timings', 'runs': {}}
    for d in map(Path, dirs):
        cols, rows, a = ledger(d)
        rec = {'steps': rows}
        for c in cols:
            if c.startswith('t_'):
                v = [a[r * len(cols) + cols.index(c)] for r in range(rows)]
                rec[c] = {'median': statistics.median(v), 'mean': statistics.fmean(v)}
        sp = d / 'subphase-schema.json'
        if sp.exists():
            sc = json.loads(sp.read_text())['columns']
            raw = (d / 'subphase-timings.bin').read_bytes()
            b = array.array('d')
            b.frombytes(raw[:len(raw) // (8 * len(sc)) * 8 * len(sc)])
            m = len(b) // len(sc)
            for i, c in enumerate(sc):
                if i == 0:
                    continue
                v = [b[r * len(sc) + i] for r in range(m)]
                rec['sub:' + c] = {'median': statistics.median(v), 'mean': statistics.fmean(v)}
        out['runs'][str(d)] = rec
    return out


def checks(A, B):
    a = json.loads((Path(A) / 'checks-result.json').read_text())['cases']
    b = json.loads((Path(B) / 'checks-result.json').read_text())['cases']
    pa = [(x['test'], x['name'], x['passed']) for x in a]
    pb = [(x['test'], x['name'], x['passed']) for x in b]
    return {'mode': 'checks', 'v6': str(A), 'v61': str(B), 'cases_v6': len(pa), 'cases_v61': len(pb), 'passed_v6': sum(x[2] for x in pa),
            'passed_v61': sum(x[2] for x in pb), 'same_cases_and_outcomes': pa == pb,
            'refusal_points': [[x['test'], x['name'], x['evidence'].get('refusal_point'), y['evidence'].get('refusal_point')]
                               for x, y in zip(a, b) if 'refusal_point' in x['evidence']]}


if __name__ == '__main__':
    args = sys.argv[1:]
    out_path = None
    if '--out' in args:
        i = args.index('--out'); out_path = args[i + 1]; del args[i:i + 2]
    steps = None
    if '--steps' in args:
        i = args.index('--steps'); steps = int(args[i + 1]); del args[i:i + 2]
    if args and args[0] == 'identity':
        r = identity(args[1], args[2], steps)
    elif args and args[0] == 'timings':
        r = timings(args[1:])
    elif args and args[0] == 'checks':
        r = checks(args[1], args[2])
    else:
        raise SystemExit(__doc__)
    text = json.dumps(r, indent=1)
    if out_path:
        Path(out_path).write_text(text + '\n')
    print(text if len(text) < 6000 else text[:6000] + '\n...')
