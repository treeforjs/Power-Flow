#!/usr/bin/python -I
"""Coordinator v8 comparisons (post-processing only; not part of the binary).

  identity A B [--steps N] [--first-step K] [--out JSON]
      A: reference run directory (coordinator v7.1 or v8), B: candidate. Byte identity of genesis.json (unless
      --first-step > 1), ledger-schema.json, deterministic.bin rows, step-digests.txt lines and step-receipts.jsonl lines of the
      compared steps; every deterministic ledger column by name; the control ledger without its timing columns; probes.bin
      rows (step >= 1); frames present in both. --first-step K compares the reference's steps K.. with the candidate's rows
      from its first row on (a restored run that starts at step K).
  physics A B [--steps N] [--out JSON]
      Owner-type parity (host Owner vs ResidentOwner): every deterministic ledger column, frames, probes and the control
      ledger compared bit for bit; step digests compared column by column (time, native prefix, field, heat, kinetic
      ledgers, source seal, RNG, arrival), the field-state digests are expected to differ by design (different
      field-state digest schemes), so differences are reported per column.
  checks A B [--out JSON]
      Case-by-case pass/fail of two checks-result.json files (the reference cases first, then the new ones).
"""
import array, hashlib, json, sys
from pathlib import Path


def sha(b):
    return hashlib.sha256(b).hexdigest()


def table(d, data, schema):
    cols = json.loads((d / schema).read_text())['columns']
    raw = (d / data).read_bytes()
    n = len(cols)
    rows = len(raw) // (8 * n)
    a = array.array('d')
    a.frombytes(raw[:rows * n * 8])
    return cols, rows, a


def deterministic(name):
    return not name.startswith('t_') and name not in ('alloc_swap_window', 'rss_kib')


def rows_of(cols, a, rows):
    n = len(cols)
    return [a[r * n:(r + 1) * n] for r in range(rows)]


def bits(x):
    return array.array('d', [x]).tobytes()


def compare_tables(A, B, data, schema, first_step, steps, keep):
    ca, ra, la = table(A, data, schema)
    cb, rb, lb = table(B, data, schema)
    RA, RB = rows_of(ca, la, ra), rows_of(cb, lb, rb)
    if 'step' in ca:
        si = ca.index('step')
        RA = [r for r in RA if r[si] >= first_step]
    n = min(len(RA), len(RB)) if steps is None else steps
    differing = []
    for c in ca:
        if not keep(c):
            continue
        if c not in cb:
            differing.append(c + ' (missing)')
            continue
        i, j = ca.index(c), cb.index(c)
        if any(bits(RA[r][i]) != bits(RB[r][j]) for r in range(n)):
            differing.append(c)
    return {'rows_compared': n, 'columns_compared': sum(1 for c in ca if keep(c)), 'columns_differing': differing, 'identical': not differing and n > 0}


def lines(p, first_step, key):
    out = []
    if not p.exists():
        return out
    for line in p.read_text().splitlines():
        s = key(line)
        if s is not None and s >= first_step:
            out.append(line)
    return out


def identity(A, B, steps=None, first_step=1):
    A, B = Path(A), Path(B)
    out = {'mode': 'identity', 'reference': str(A), 'candidate': str(B), 'first_step': first_step, 'files': {}}
    files = ['ledger-schema.json'] + (['genesis.json'] if first_step <= 1 and (B / 'genesis.json').exists() else [])
    for f in files:
        x, y = (A / f).read_bytes(), (B / f).read_bytes()
        out['files'][f] = {'reference_sha256': sha(x), 'candidate_sha256': sha(y), 'identical': x == y}
    out['ledger'] = compare_tables(A, B, 'ledger.bin', 'ledger-schema.json', first_step, steps, deterministic)
    n = out['ledger']['rows_compared']
    ca = json.loads((A / 'ledger-schema.json').read_text())['columns']
    ndet = sum(1 for c in ca if deterministic(c))
    x = (A / 'deterministic.bin').read_bytes()[(first_step - 1) * ndet * 8:(first_step - 1 + n) * ndet * 8]
    y = (B / 'deterministic.bin').read_bytes()[:n * ndet * 8]
    out['files']['deterministic.bin'] = {'rows': n, 'identical': x == y and len(x) == n * ndet * 8}
    da = lines(A / 'step-digests.txt', first_step, lambda l: int(l.split()[0]))[:n]
    db = lines(B / 'step-digests.txt', first_step, lambda l: int(l.split()[0]))[:n]
    out['files']['step-digests.txt'] = {'lines': len(da), 'identical': da == db and len(da) == n,
                                        'first_difference': next((i + first_step for i, (p, q) in enumerate(zip(da, db)) if p != q), None)}
    ra = lines(A / 'step-receipts.jsonl', first_step, lambda l: json.loads(l).get('step'))[:n]
    rb = lines(B / 'step-receipts.jsonl', first_step, lambda l: json.loads(l).get('step'))[:n]
    out['files']['step-receipts.jsonl'] = {'lines': len(ra), 'identical': ra == rb and len(ra) == n,
                                           'first_difference': next((i + first_step for i, (p, q) in enumerate(zip(ra, rb)) if p != q), None)}
    out['control'] = compare_tables(A, B, 'control-ledger.bin', 'control-schema.json', first_step, n, lambda c: not c.startswith('t_'))
    pa, pb = A / 'probes.bin', B / 'probes.bin'
    if pa.exists() and pb.exists():
        cols = json.loads((A / 'probes-schema.json').read_text())['columns']
        w = len(cols) * 8
        xa, xb = pa.read_bytes(), pb.read_bytes()
        RA = [xa[i:i + w] for i in range(0, len(xa) - w + 1, w)]
        RB = [xb[i:i + w] for i in range(0, len(xb) - w + 1, w)]
        RA = [r for r in RA if array.array('d', r[:8])[0] >= first_step]
        RB = [r for r in RB if array.array('d', r[:8])[0] >= first_step]
        k = min(len(RA), len(RB), n)
        out['files']['probes.bin'] = {'rows': k, 'identical': RA[:k] == RB[:k] and k > 0}
    fa, fb = A / 'frames.bin', B / 'frames.bin'
    if fa.exists() and fb.exists() and first_step <= 1:
        lay = json.loads((A / 'frames-layout.json').read_text())
        nx, ny = lay['nx'], lay['ny']
        size = 8 * (2 + 6 * nx * ny + (nx + 1) * (ny + 1))
        ba, bb = fa.read_bytes(), fb.read_bytes()
        k = min(len(ba), len(bb)) // size
        same = [ba[i * size:(i + 1) * size] == bb[i * size:(i + 1) * size] for i in range(k)]
        out['files']['frames.bin'] = {'frames_compared': k, 'identical_frames': sum(same), 'identical': k > 0 and all(same)}
    out['bit_identical'] = all(v.get('identical', True) for v in out['files'].values()) and out['ledger']['identical'] and out['control']['identical']
    return out


def physics(A, B, steps=None):
    A, B = Path(A), Path(B)
    out = {'mode': 'physics', 'reference': str(A), 'candidate': str(B)}
    out['ledger'] = compare_tables(A, B, 'ledger.bin', 'ledger-schema.json', 1, steps, deterministic)
    out['control'] = compare_tables(A, B, 'control-ledger.bin', 'control-schema.json', 1, steps, lambda c: not c.startswith('t_'))
    n = out['ledger']['rows_compared']
    names = ['step', 'time', 'native_prefix', 'field', 'heat', 'kinetic_ledgers', 'source_seal', 'rng', 'arrival']
    da = [l.split() for l in (A / 'step-digests.txt').read_text().splitlines()[:n]]
    db = [l.split() for l in (B / 'step-digests.txt').read_text().splitlines()[:n]]
    out['step_digest_columns_identical'] = {nm: all(x[i] == y[i] for x, y in zip(da, db)) for i, nm in enumerate(names)}
    for f in ('probes.bin', 'frames.bin'):
        if (A / f).exists() and (B / f).exists():
            out[f] = (A / f).read_bytes() == (B / f).read_bytes()
    return out


def checks(A, B):
    a = json.loads((Path(A) / 'checks-result.json').read_text())['cases']
    b = json.loads((Path(B) / 'checks-result.json').read_text())['cases']
    head = b[:len(a)]
    same = [(x['test'], x['name'], x['passed']) == (y['test'], y['name'], y['passed']) for x, y in zip(a, head)]
    evidence = []
    for x, y in zip(a, head):
        ex, ey = x['evidence'], y['evidence']
        if ex != ey:
            evidence.append({'test': x['test'], 'name': x['name'], 'differences': {k: [ex.get(k), ey.get(k)] for k in sorted(set(ex) | set(ey)) if ex.get(k) != ey.get(k)}})
    return {'mode': 'checks', 'reference': str(A), 'candidate': str(B), 'cases_reference': len(a), 'cases_candidate': len(b),
            'passed_reference': sum(x['passed'] for x in a), 'passed_candidate': sum(x['passed'] for x in b),
            'reference_cases_first_same_names_and_outcomes': len(head) == len(a) and all(same), 'evidence_differences': evidence,
            'new_cases': [(x['test'], x['name'], x['passed']) for x in b[len(a):]]}


if __name__ == '__main__':
    args = sys.argv[1:]
    opts = {}
    for k in ('--out', '--steps', '--first-step'):
        if k in args:
            i = args.index(k)
            opts[k] = args[i + 1]
            del args[i:i + 2]
    steps = int(opts['--steps']) if '--steps' in opts else None
    if args and args[0] == 'identity':
        r = identity(args[1], args[2], steps, int(opts.get('--first-step', 1)))
    elif args and args[0] == 'physics':
        r = physics(args[1], args[2], steps)
    elif args and args[0] == 'checks':
        r = checks(args[1], args[2])
    else:
        raise SystemExit(__doc__)
    text = json.dumps(r, indent=1)
    if '--out' in opts:
        Path(opts['--out']).parent.mkdir(parents=True, exist_ok=True)
        Path(opts['--out']).write_text(text + '\n')
    print(text if len(text) < 5000 else text[:5000] + '\n...')
