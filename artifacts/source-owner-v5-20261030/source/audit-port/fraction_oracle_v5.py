#!/usr/bin/env python3
"""Exact rational re-check of the source-v5 port's count and A08 precision operands.

Derived from the source-v4 port oracle (same 13 count checks), extended with the
A07 stock-drain refusal and the A08 G1/G2 rows and sweeps of the v5 port.

Independent of ExactParticleCount: every binary64 operand is read from its
hexadecimal form with float.fromhex and summed with fractions.Fraction.
The audit's own recorded values (reproducer.log, exact-allocation-evidence.json)
are re-read so the constants embedded in the port are checked against them.

usage: fraction_oracle_v5.py <exact-counts.json> <port-results.json> <audit-source-review-dir> <out.json>
"""
import json
import re
import sys
from decimal import Decimal, getcontext
from fractions import Fraction

getcontext().prec = 80
T = Fraction(1, 2**40)
EB = Fraction(1, 2**58)
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

    # ---- v5 additions
    def finite(h):
        x = float.fromhex(h)
        return x == x and abs(x) != float('inf')

    def representable(x):
        try:
            y = float(x)
        except OverflowError:
            return False
        return abs(y) != float('inf') and Fraction(y) == x

    def justified(r):
        """observed is not the exact finite binary64 value reference*2^e."""
        if not finite(r['reference']):
            return False
        target = F(r['reference']) * Fraction(2) ** r['binary_exponent']
        if not representable(target) or not finite(r['observed']):
            return True
        return F(r['observed']) != target

    # A07 extension: the rounded full release exceeds the exact stock.
    stock = total(counts['drain_stock'])
    check('drain: exact stock equals Count2 input {1e15,-0.01}', stock == count2)
    check('drain: refused rounded release exceeds the exact stock (Fraction)', F(counts['drain_released']) > stock,
          str(F(counts['drain_released']) - stock))
    check('drain: port row is a typed stock-drain refusal',
          probes['EXT-A07-full-release-of-multiscale-stock']['v5_verdict'].startswith('refused (typed stock-drain'))

    # A08 G1: the w=1e-301 birth refusal is exactly justified.
    g1 = counts['g1_refusal']
    check('G1: 1e-301 birth refused (typed weighted-representation) and justified exactly',
          g1['kind'] == 'weighted-representation' and justified(g1), g1)
    # A08 G2 and the exposure sweep against the exact finite-log count.
    D = lambda h: Decimal(float.fromhex(h))
    current = D(counts['exposure_log_current']).exp()
    area, dt, qe, mean = D(counts['exposure_area']), D(counts['exposure_dt']), D(counts['qe']), D(counts['exposure_mean_eV'])
    Td = Decimal(T.numerator) / Decimal(T.denominator)
    floor = Td - 2 * Decimal(EB.numerator) / Decimal(EB.denominator)

    def exact(frac_hex):
        a = D(frac_hex)
        return current * area * a * dt / qe, -current * mean * area * a * dt
    g2 = counts['g2_refusal']
    cx, _ = exact(g2['area_fraction'])
    check('G2: area 3.16e-245 refused (typed exposure-scalar); refused scalar outside T exactly',
          g2['kind'] == 'exposure-scalar' and abs(D(g2['observed']) - cx) / cx > floor, float(abs(D(g2['observed']) - cx) / cx))
    adm = ref = 0
    bad = []
    for row in counts['exposure_sweep']:
        cx, hx = exact(row['area_fraction'])
        if row['admitted']:
            adm += 1
            rc = abs(D(row['count']) - cx) / cx
            rh = abs(D(row['diagnostic']) - hx) / abs(hx)
            if not (rc <= Td and rh <= Td):
                bad.append((row['area_fraction'], float(rc), float(rh)))
        else:
            ref += 1
            obs = D(row['observed']) if finite(row['observed']) else Decimal(0)
            rel = (abs(obs - cx) / cx) if row['quantity'] == 'electron expected exposure count' else (abs(obs - hx) / abs(hx))
            if not (obs == 0 or rel > floor):
                bad.append((row['area_fraction'], row['quantity'], float(rel)))
    check(f'exposure sweep: {adm} admitted within T, {ref} refused outside T - 2*evaluation_bound (exact)', not bad and adm > 0 and ref > 0, bad[:5])
    # A08 weight sweep: admitted values within T exactly (gross for arrays); refusals justified exactly.
    adm = ref = 0
    bad = []
    for row in counts['weight_sweep']:
        w = F(row['weight'])
        if not row['admitted']:
            ref += 1
            if row['kind'] == 'weighted-representation' and not justified(row):
                bad.append((row['weight'], row['quantity']))
            continue
        adm += 1
        items = [(f'exchange[{i}]', F(u), F(a), abs(F(u))) for i, (u, a) in enumerate(zip(row['unit_exchange'], row['exchange']))]
        for key in ('handoff_current_load', 'handoff_endpoint_charge_impulse'):
            for i, u, a, g in row.get(key, []):
                items.append((f'{key}[{i}]', F(u), F(a), F(g)))
        for name, u, a, g in items:
            ok = (a == 0) if g == 0 else abs(a - w * u) <= T * w * g
            if not ok:
                bad.append((row['weight'], name))
    check(f'weight sweep: {adm} admitted weights within T exactly; {ref} refusals justified', not bad and adm > 0 and ref > 0, bad[:5])

    out = {'schema': 'source-v5-port-fraction-oracle-v1', 'checks': checks,
           'passed': sum(c['passed'] for c in checks), 'failed': sum(not c['passed'] for c in checks)}
    Path(sys.argv[4]).write_text(json.dumps(out, indent=1) + '\n')
    print(f"fraction oracle: {out['passed']} passed, {out['failed']} failed")
    return 0 if out['failed'] == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
