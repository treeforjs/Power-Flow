#!/usr/bin/env python3
"""Independent exact-arithmetic oracle for the source owner v5 precision checks.

Test oracle only; never used by the owner. Every binary64 operand is read from its
hexadecimal form with float.fromhex and compared with fractions.Fraction. The single
non-rational reference, exp(log current), is evaluated with decimal at 80 digits.

Verified for precision-dump.json written by the v5 focused tests:
  * every ADMITTED weighted value a at weight w satisfies |a - w*u| <= T*w*g exactly,
    where u is the unit-weight value and g = |u| (scalars) or the unit gross
    magnitude (deposited arrays); g == 0 requires a == 0;
  * every weighted_representation REFUSAL is justified exactly: its observed value is
    not the exact finite binary64 rescaling reference*2^e;
  * every ADMITTED exposure count/diagnostic is within T of exp(log current)*area*
    fraction*dt/qe (resp. -exp(log current)*mean*area*fraction*dt); every exposure
    REFUSAL has a relative error above T - 2*evaluation_bound (or a zero scalar);
  * the stock-drain refusal is justified: released > exact stock (Fraction);
  * sealed-guard conversions are typed; all normal-range weights were admitted.

usage: precision_oracle.py <precision-dump.json> <out.json>
"""
import json
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

getcontext().prec = 80


def F(h):
    return Fraction(float.fromhex(h))


def finite_hex(h):
    x = float.fromhex(h)
    return x == x and abs(x) != float('inf')


def representable(x):
    """True iff the rational x is exactly a finite binary64 number."""
    try:
        y = float(x)
    except OverflowError:
        return False
    return abs(y) != float('inf') and Fraction(y) == x


def D(h):
    return Decimal(float.fromhex(h))


def main():
    dump = json.loads(Path(sys.argv[1]).read_text())
    checks = []
    stats = {}

    def check(name, ok, detail=None):
        checks.append({'name': name, 'passed': bool(ok), 'detail': detail})
        return bool(ok)

    T = F(dump['relative_tolerance'])
    eb = F(dump['evaluation_bound'])
    check('tolerance is 2^-40', T == Fraction(1, 2**40), str(T))
    check('evaluation bound is 2^-58', eb == Fraction(1, 2**58), str(eb))

    def within(u, a, w, g):
        """Exact |a - w*u| <= T*w*g; returns (ok, exact relative error vs w*g)."""
        if g == 0:
            return a == 0, (Fraction(0) if a == 0 else None)
        err = abs(a - w * u) / (w * g)
        return err <= T, err

    def justified_weighted(r):
        """A weighted_representation refusal: observed != reference*2^e exactly, or not representable."""
        if r['kind'] != 'weighted-representation':
            return None
        if not (finite_hex(r['reference'])):
            return None
        target = F(r['reference']) * Fraction(2) ** r['binary_exponent']
        if not representable(target):
            return True
        if not finite_hex(r['observed']):
            return True
        return F(r['observed']) != target

    def typed(r):
        return bool(r['quantity']) and bool(r['identity']) and finite_hex(r['log_abs_scope'])

    # ---- G1 and G2 single cases
    g1 = dump['g1_refusal']
    check('G1 w=1e-301 birth refused with typed scope', typed(g1), g1['quantity'])
    check('G1 refusal justified exactly', justified_weighted(g1) is True,
          {'quantity': g1['quantity'], 'observed': g1['observed'], 'reference': g1['reference'], 'e': g1['binary_exponent']})

    ex = dump['exposure_sweep']
    area, dt, qe = D(ex['area_m2']), D(ex['duration_s']), D(ex['qe'])
    current = D(ex['log_current_A_m2']).exp()
    mean = D(ex['mean_metal_energy_eV'])
    Td = Decimal(T.numerator) / Decimal(T.denominator)
    refuse_floor = Td - 2 * Decimal(eb.numerator) / Decimal(eb.denominator)

    def exposure_exact(fraction_hex):
        a = D(fraction_hex)
        return current * area * a * dt / qe, -current * mean * area * a * dt

    g2 = dump['g2_refusal']
    count_exact, _ = exposure_exact(float.hex(3.1622776601683792e-245))
    rel_g2 = abs(D(g2['observed']) - count_exact) / count_exact
    check('G2 area 3.1622776601683792e-245 refused as exposure-scalar with typed scope',
          typed(g2) and g2['kind'] == 'exposure-scalar' and g2['quantity'] == 'electron expected exposure count')
    check('G2 refused scalar has exact relative error above tolerance', rel_g2 > refuse_floor, str(rel_g2))

    # ---- birth sweeps
    for sweep in dump['birth_sweeps']:
        name = sweep['fixture']
        n_adm = n_ref = n_val = 0
        worst = Fraction(0)
        bad = []
        unjust = []
        untyped = []
        for row in sweep['rows']:
            w = F(row['weight'])
            if row['admitted']:
                n_adm += 1
                items = [(v['name'], F(v['unit']), F(v['actual']), F(v['gross'])) for v in row['values']]
                for arr in row['arrays']:
                    for i, u, a, g in arr['entries']:
                        items.append((arr['name'] + f'[{i}]', F(u), F(a), F(g)))
                for qn, u, a, g in items:
                    n_val += 1
                    ok, err = within(u, a, w, g)
                    if not ok:
                        bad.append((row['weight'], qn))
                    elif err is not None and err > worst:
                        worst = err
            else:
                n_ref += 1
                if 'refusal' not in row:
                    untyped.append(row['weight'])
                    continue
                r = row['refusal']
                if not typed(r):
                    untyped.append(row['weight'])
                j = justified_weighted(r)
                if j is False:
                    unjust.append((row['weight'], r['quantity']))
        stats[name] = {'admitted_weights': n_adm, 'refused_weights': n_ref, 'admitted_values_checked': n_val,
                       'worst_exact_relative_error': float(worst)}
        check(f'{name}: every admitted weighted value within T exactly ({n_val} values, {n_adm} weights)', not bad, bad[:5])
        check(f'{name}: every refusal typed with quantity, finite log scope and identity ({n_ref})', not untyped, untyped[:5])
        check(f'{name}: every weighted-representation refusal justified exactly', not unjust, unjust[:5])
        check(f'{name}: sweep covers admitted and refused weights', n_adm > 0 and n_ref > 0)

    # ---- exposure sweep
    n_adm = n_ref = 0
    bad = []
    unjust = []
    worst = Decimal(0)
    for row in ex['rows']:
        cx, hx_ = exposure_exact(row['area_fraction'])
        if row['admitted']:
            n_adm += 1
            rc = abs(D(row['count']) - cx) / cx
            rh = abs(D(row['diagnostic_J']) - hx_) / abs(hx_) if hx_ != 0 else (Decimal(0) if D(row['diagnostic_J']) == 0 else Decimal(1))
            worst = max(worst, rc, rh)
            if not (rc <= Td and rh <= Td):
                bad.append((row['area_fraction'], float(rc), float(rh)))
        else:
            n_ref += 1
            r = row['refusal']
            obs = D(r['observed']) if finite_hex(r['observed']) else None
            if r['quantity'] == 'electron expected exposure count':
                rel = Decimal(1) if obs in (None, 0) else abs(obs - cx) / cx
            else:
                rel = Decimal(1) if obs is None else (abs(obs - hx_) / abs(hx_) if hx_ != 0 else Decimal(1))
            if not (typed(r) and r['kind'] == 'exposure-scalar' and rel > refuse_floor):
                unjust.append((row['area_fraction'], r['quantity'], float(rel)))
    stats['exposure'] = {'admitted': n_adm, 'refused': n_ref, 'worst_exact_relative_error': float(worst)}
    check(f'exposure: every admitted count and diagnostic within T of the exact finite-log value ({n_adm})', not bad, bad[:5])
    check(f'exposure: every refusal typed and its scalar outside T - 2*evaluation_bound ({n_ref})', not unjust, unjust[:5])
    check('exposure: sweep covers admitted and refused areas', n_adm > 0 and n_ref > 0)

    # ---- neutral leaves
    nl = dump['neutral_leaves']
    outcomes = []
    for case in nl['cases']:
        total = sum((F(x) for x in case['leaves']), Fraction(0))
        check(f"neutral '{case['partition']}': leaves are an exact partition of 1", total == 1 and case['exact_partition_of_one'])
        if case['admitted']:
            bad = []
            for v in case['values']:
                ok, _ = within(F(v['unit']), F(v['actual']), F(v['weight']), F(v['gross']))
                if not ok:
                    bad.append((v['weight'], v['name']))
            check(f"neutral '{case['partition']}': admitted values within T exactly ({len(case['values'])})", not bad, bad[:5])
        else:
            r = case['refusal']
            check(f"neutral '{case['partition']}': refusal typed and justified exactly", typed(r) and justified_weighted(r) is True,
                  {'quantity': r['quantity'], 'identity': r['identity']})
        outcomes.append((case['partition'], case['admitted']))
    stats['neutral'] = outcomes
    check('neutral: the subnormal-reaching partition is refused and the normal-range one admitted',
          dict(outcomes).get('chain to 2^-1060') is False and dict(outcomes).get('normal-range leaves') is True)

    # ---- stock drain
    sd = dump['stock_drain']
    stock = sum((F(x) for x in sd['stock_carriers']), Fraction(0))
    qr = sd['quote_refusal']
    released = F(qr['observed'])
    check('stock drain: quote refusal typed as stock-drain with stock identity',
          typed(qr) and qr['kind'] == 'stock-drain' and 'exact-stock-sha256=' in qr['identity'])
    check('stock drain: rounded release exceeds the exact stock (Fraction)', released > stock, str(released - stock))
    check('stock drain: exact stock equals Count2 input {1e15,-0.01}', stock == Fraction(1e15) + Fraction(-0.01))
    check('stock drain: prepare refusal typed as stock-drain, owner unchanged',
          sd['prepare_refusal']['kind'] == 'stock-drain' and sd['owner_unchanged'])

    # ---- sealed guard conversion and normal range
    rows = dump['sealed_guard_conversion']
    check('sealed energy-closure guards convert to typed sealed-guard-at-weight refusals',
          len(rows) == 3 and all(r['kind'] == 'sealed-guard-at-weight' and typed(r) for r in rows),
          [r['quantity'] for r in rows])
    nr = dump['normal_range_weights']
    check(f'normal-range weights all admitted ({len(nr)})', all(r['admitted'] for r in nr))

    out = {'schema': 'source-v5-precision-oracle-v1', 'oracle': 'Python fractions.Fraction; decimal (80 digits) for exp only; owner not imported',
           'tolerance': str(T), 'stats': stats, 'checks': checks,
           'passed': sum(c['passed'] for c in checks), 'failed': sum(not c['passed'] for c in checks)}
    Path(sys.argv[2]).write_text(json.dumps(out, indent=1) + '\n')
    print(f"precision oracle: {out['passed']} passed, {out['failed']} failed")
    return 0 if out['failed'] == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
