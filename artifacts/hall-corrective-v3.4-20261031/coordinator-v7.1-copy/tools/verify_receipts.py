#!/usr/bin/env python3
"""Verify the sealed step receipts of a coordinator v7 run directory (post-processing only).

  verify_receipts.py RUN_DIR [--out JSON]

Checks: the chain (receipt k: sha256(previous digest + canonical receipt without its own digest)) from the recorded
genesis link of run-summary.json to the recorded head; each decision digest = sha256 of the canonical decision
record; the dt bits of every receipt equal the dt_s column of the step ledger and of the control ledger; the
receipt's after-digests and end time equal the step's line of step-digests.txt; every attempt but the last is a
refusal and the last attempt's decision digest is the receipt's.
"""
import hashlib, json, struct, sys
from pathlib import Path


def sha(b):
    return hashlib.sha256(b).hexdigest()


def raw_object(text, key):
    """The exact recorded text of the JSON object value of a top-level key (the digests are over these bytes)."""
    i = text.index('"%s":{' % key) + len(key) + 3
    depth, k, string = 0, i, False
    while True:
        c = text[k]
        if string:
            if c == '\\':
                k += 1
            elif c == '"':
                string = False
        elif c == '"':
            string = True
        elif c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                return text[i:k + 1]
        k += 1


def column(d, binname, schema, name):
    cols = json.loads((d / schema).read_text())['columns']
    raw = (d / binname).read_bytes()
    n = len(cols)
    i = cols.index(name)
    return [struct.unpack_from('<d', raw, 8 * (r * n + i))[0] for r in range(len(raw) // (8 * n))]


def main():
    d = Path(sys.argv[1])
    summary = json.loads((d / 'run-summary.json').read_text())
    lines = (d / 'step-receipts.jsonl').read_text().splitlines()
    digests = (d / 'step-digests.txt').read_text().splitlines()
    dt_ledger = column(d, 'ledger.bin', 'ledger-schema.json', 'dt_s')
    dt_control = column(d, 'control-ledger.bin', 'control-schema.json', 'dt_s')
    head = summary['receipt_chain_genesis']
    bad = []
    steps = 0
    for k, line in enumerate(lines):
        rec = json.loads(line)
        if rec['schema'] != 'coordinator-v7-step-receipt-v1':
            if rec['previous_sha256'] != head:
                bad.append((k + 1, 'refusal record not linked to the head'))
            continue
        steps += 1
        own = rec.pop('receipt_sha256')
        if rec['previous_sha256'] != head:
            bad.append((k + 1, 'previous digest'))
        text = line.replace(',"receipt_sha256":"%s"' % own, '', 1)
        if json.loads(text) != rec:
            bad.append((k + 1, 'receipt text differs from its record'))
        head = sha((head + text).encode())
        if head != own:
            bad.append((k + 1, 'receipt digest'))
        if sha(raw_object(text, 'decision').encode()) != rec['decision_sha256']:
            bad.append((k + 1, 'decision digest'))
        at = rec['attempts']
        if at[-1]['outcome'] != 'prepared' or any(a['outcome'] != 'refused' for a in at[:-1]) or at[-1]['decision_sha256'] != rec['decision_sha256']:
            bad.append((k + 1, 'attempt record'))
        bits = int(rec['dt_bits'], 16)
        s = rec['step'] - 1
        for name, col in (('ledger', dt_ledger), ('control', dt_control)):
            if struct.unpack('<Q', struct.pack('<d', col[s]))[0] != bits:
                bad.append((k + 1, 'dt bits differ from the %s ledger' % name))
        if float.fromhex(rec['dt_s']) != col[s]:
            bad.append((k + 1, 'dt text'))
        f = digests[s].split()
        a = rec['after']
        if int(f[0]) != rec['step'] or float.fromhex(f[1]) != float.fromhex(a['time_s']) or \
                f[2:] != [a['native_prefix'], a['field'], a['heat'], a['kinetic_ledgers'], a['source_seal'], a['rng'], a['arrival']]:
            bad.append((k + 1, 'after-digests differ from step-digests.txt'))
    if head != summary['receipt_chain_head']:
        bad.append((0, 'final head differs from run-summary.json'))
    out = {'run': str(d), 'receipts': steps, 'records': len(lines), 'head': head, 'head_matches_summary': head == summary['receipt_chain_head'],
           'failures': bad[:50], 'verified': not bad and steps == summary['steps_completed']}
    text = json.dumps(out, indent=1)
    if '--out' in sys.argv:
        Path(sys.argv[sys.argv.index('--out') + 1]).write_text(text + '\n')
    print(text)
    return 0 if out['verified'] else 1


if __name__ == '__main__':
    sys.exit(main())
