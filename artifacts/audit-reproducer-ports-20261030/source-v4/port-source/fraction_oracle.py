#!/usr/bin/env python3
"""Exact rational re-check of the source-v4 port's count operands.

Independent of ExactParticleCount: every binary64 operand is read from its
hexadecimal form with float.fromhex and summed with fractions.Fraction.
The audit's own recorded values (reproducer.log, exact-allocation-evidence.json)
are re-read so the constants embedded in the port are checked against them.

usage: fraction_oracle.py <exact-counts.json> <port-results.json> <audit-source-review-dir> <out.json>
"""
import json
import re
import sys
from fractions import Fraction
from pathlib import Path


def F(h):
    return Fraction(float.fromhex(h))


def total(xs):
    return sum((F(x) for x in xs), Fraction(0))


def main():
    counts = json.loads(Path(sys.argv[1]).read_text())
    results = json.loads(Path(sys.argv[2]).read_text())
    audit = Path(sys.argv[3])
    probes = {p['id']: p for p in results['probes']}
    checks = []

    def check(name, ok, detail=None):
        checks.append({'name': name, 'passed': bool(ok), 'detail': detail})

    log = (audit / 'reproducer.log').read_text()
    evidence = json.loads((audit / 'exact-allocation-evidence.json').read_text())

    # Embedded audit constants match the audit's own records.
    v3_weight = float(re.search(r'actual_positive_packet_weight=(\S+)', log).group(1))
    v3_log = float(re.search(r'log_count=(\S+)', log).group(1))
    check('c1 port weight equals reproducer.log packet weight',
          F(probes['A07-C1-actual-release']['v3_packet_weight_hex']) == Fraction(v3_weight))
    check('l1 port log constant equals reproducer.log log_count',
          F(probes['A08-L1-exposure-underflow']['v3_log_count_hex']) == Fraction(v3_log))
    check('c2 port weights equal exact-allocation-evidence.json',
          [F(h) for h in counts['c2_audit_weights']] == [F(h) for h in evidence['packet_weights_hex']])

    # A07/C1: exact source + returned leaf conservation, independent arithmetic.
    initial = total(counts['c1_initial'])
    count2 = total(counts['c1_count2_input'])
    check('c1 initial exact stock equals Count2 input hi+lo exactly', initial == count2, str(initial))
    if 'c1_accepted' in counts:
        accepted = total(counts['c1_accepted'])
        returned = total(counts['c1_returned'])
        check('c1 returned leaves equal the caller leaves', returned == total(counts['c1_leaves']))
        check('c1 accepted stock + returned leaves == initial stock (exact)', accepted + returned == initial,
              str(initial - accepted - returned))
        check('c1 accepted stock differs from initial (debit not lost)', accepted != initial)
    else:
        check('c1 accepted stock present', False)
    arith_after = total(counts['c1_arith_after'])
    check('c1 arithmetic: debit result + 1e-30 == {1e15,-0.01} (exact)',
          arith_after + total(counts['c1_arith_amount']) == count2)

    # A07/C2: audit excess reproduced; exact partition admitted exactly.
    audit_sum = total(counts['c2_audit_weights'])
    check('c2 audit ten weights exceed the unit debit by exactly 5/2^55',
          audit_sum - 1 == Fraction(5, 2**55), str(audit_sum - 1))
    check('c2 recorded rational excess matches evidence',
          audit_sum - 1 == Fraction(evidence['packet_sum_minus_debit_exact']))
    check('c2 exact partition sums to exactly 1', total(counts['c2_exact_partition']) == 1)
    if 'c2_exact_returned' in counts:
        check('c2 returned exact-partition leaves sum to exactly 1', total(counts['c2_exact_returned']) == 1)
        check('c2 returned leaves are the caller leaves verbatim',
              [F(h) for h in counts['c2_exact_returned']] == [F(h) for h in counts['c2_exact_partition']])
    else:
        check('c2 exact partition returned', False)

    out = {'schema': 'source-v4-port-fraction-oracle-v1', 'checks': checks,
           'passed': sum(c['passed'] for c in checks), 'failed': sum(not c['passed'] for c in checks)}
    Path(sys.argv[4]).write_text(json.dumps(out, indent=1) + '\n')
    print(f"fraction oracle: {out['passed']} passed, {out['failed']} failed")
    return 0 if out['failed'] == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
