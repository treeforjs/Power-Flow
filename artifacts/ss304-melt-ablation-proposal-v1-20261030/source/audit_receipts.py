#!/usr/bin/env python3
"""Independent Fraction/Decimal reconstruction of this manufactured fixture.

No production state or library is loaded. The sealed caloric equation is
re-expressed here only to check threshold/projection wiring, not re-qualified.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys

checks = 0


def check(ok, message):
    global checks
    if not ok:
        raise AssertionError(message)
    checks += 1


def dy(j):
    return F(int(j['hex'], 16)) * F(2) ** j['exponent']


def number(x):
    return F(float(x))


def count(stock):
    return sum(map(number, stock['count_positive_carriers']), F(0))


def metal_u(t):
    def cp(x):
        return 472.0 + 0.136*x - 2820000.0/(x*x)
    def hs(x):
        return (x-300.0)*(472.0+0.068*(x+300.0)-2820000.0/(x*300.0))
    def mush(dt):
        cs = cp(1672.0)
        return hs(1672.0)+cs*dt+0.5*(800.0-cs)*dt*dt/55.0
    if t < 1672.0:
        return hs(t)
    if t < 1727.0:
        return mush(t-1672.0)+276320.0*((t-1672.0)/55.0)
    return mush(55.0)+276320.0+800.0*(t-1727.0)


def u(stock_id, t):
    return number(metal_u(t) if stock_id == 'metal' else 2.0*(t-300.0))


def thermal(stocks, t):
    return sum((dy(s['mass_kg'])*u(s['id'], t) for s in stocks.values()
                if s['id'] != 'water'), F(0))


def reference(stock):
    if stock['id'] == 'metal':
        return F(0)
    if stock['id'] == 'oxide':
        return -20*dy(stock['mass_kg'])
    return -F(1, 4)*count(stock)


def mass(stocks, models):
    return sum((dy(s['mass_kg']) if s['id'] != 'water' else
                count(s)*number(models[s['id']]['unit_mass_kg'])
                for s in stocks.values()), F(0))


def projected(stocks, e, t):
    check(298.2 <= t <= 1800., 'projection domain')
    check(thermal(stocks, 298.2) <= e <= thermal(stocks, 1800.), 'exact local range')
    err = abs(thermal(stocks, t)-e)
    for adjacent in (math.nextafter(t, -math.inf), math.nextafter(t, math.inf)):
        if 298.2 <= adjacent <= 1800.:
            check(err <= abs(thermal(stocks, adjacent)-e), 'closest binary64 projection')


def audit(path):
    data = json.loads(path.read_text())
    check(data['status'] == 'passed_stateless_synthetic_component_only', 'test status')
    check(all(c['passed'] for c in data['checks']), 'all focused checks')
    models = {m['id']: m for m in data['fixture_models']}
    for m in models.values():
        check(m['parent_elements'] == m['product_elements'], 'element identity')
        check(m['product'][-1] == 0, 'neutral product')
        check(number(m['product'][2]) == number(m['unit_mass_kg']), 'product mass')
        if m['kind'] == 0:
            check(number(m['u_at_trigger_J_kg']) == u(m['id'], m['trigger_K']), 'sealed/fixture threshold law')
            actual = number(m['product'][3])+number(m['product'][4])-number(m['unit_mass_kg'])*(number(m['condensed_reference_J_kg'])+u(m['id'], m['trigger_K']))
            declared = number(m['unit_mass_kg'])*number(m['removal_J_kg_excluding_actual_K'])
            check(actual > 0 and abs(actual-declared) <= F(8)*number(sys.float_info.epsilon)*(abs(number(m['product'][3]))+abs(number(m['product'][4]))+abs(number(m['unit_mass_kg'])*number(m['condensed_reference_J_kg']))+abs(number(m['unit_mass_kg'])*u(m['id'], m['trigger_K']))+abs(declared)), 'incremental cost excludes already paid threshold')
    for p in data['proofs']:
        before = {s['id']: s for s in p['stocks_before']}
        stocks = dict(before)
        e0 = dy(p['thermal_before_J'])
        e = e0
        impulse = [F(0), F(0), F(0)]
        gas = extra = removed = F(0)
        projected(stocks, e, p['initial_projected_T_K'])
        ref0 = sum(map(reference, stocks.values()), F(0))
        check(ref0 == dy(p['reference_before_J']), 'initial derived reference')
        for row in p['actions']:
            m = models[row['stock']]
            old = stocks[row['stock']]
            n = sum((number(c['weight']) for c in row['carriers']), F(0))
            k = sum((number(c['weight'])*number(c['kinetic_per_unit_J']) for c in row['carriers']), F(0))
            dm = n*number(m['unit_mass_kg'])
            packet = n*(number(row['product_R_J'])+number(row['product_I_J']))+k
            check(n == dy(row['count']) and dm == dy(row['mass_removed_kg']), 'actual leaf sum and removed mass')
            check(packet == dy(row['gas_energy_J']), 'actual R/I/K energy')
            check(number(row['product_R_J']) == number(m['product'][3]) and number(row['product_I_J']) == number(m['product'][4]), 'fixed product reference no event retuning')
            new = row['stock_after']
            if m['kind'] == 0:
                check(dy(old['mass_kg'])-dm == dy(new['mass_kg']), 'condensed exact debit')
                check(count(new) == 0 and dy(new['mass_kg']) >= 0, 'one condensed encoding')
            else:
                check(count(old)-n == count(new) and count(new) >= 0, 'source count exact debit')
                check(dy(new['mass_kg']) == 0, 'one adsorbate encoding')
            dr = reference(new)-reference(old)
            check(dr == dy(row['reference_change_J']), 'derived bound-reference delta')
            x = dy(row['other_energy_out_J'])
            check(x >= 0 and (x == 0 or bool(row['other_energy_owner'])), 'named other-energy receiver')
            if m['kind'] == 0:
                threshold = thermal(stocks, m['trigger_K'])
                allowance = e-threshold
                cost = packet+dr+x-dm*u(row['stock'], m['trigger_K'])
                check(threshold == dy(row['trigger_thermal_energy_J']), 'same-local-inventory threshold')
                check(allowance == dy(row['available_above_trigger_J']), 'present net local allowance')
                check(cost == dy(row['incremental_cost_including_actual_K_J']), 'actual incremental debit')
                check(0 < cost <= allowance, 'no release beyond allowance')
                ratio = row['plan_scale_limit_ratio']
                check(dy(ratio['numerator']) == allowance and dy(ratio['denominator']) == cost, 'exact allowance ratio')
                volume = row['removed_reference_volume_ratio']
                check(dy(volume['mass_numerator_kg']) == dm and number(volume['density_denominator_kg_m3']) == number(m['rho']), 'volume quote from explicit rho only')
            else:
                check(row['trigger'] == 'independent_desorption_no_melt_gate' and row['plan_scale_limit_ratio'] is None, 'water has no melt budget or gate')
            for c in row['carriers']:
                for axis in range(3):
                    impulse[axis] -= number(c['weight'])*number(m['unit_mass_kg'])*number(c['proper_velocity_m_s'][axis])
                with localcontext() as ctx:
                    ctx.prec = 90
                    u2 = sum((Decimal.from_float(float(q))**2 for q in c['proper_velocity_m_s']), Decimal(0))
                    kp = Decimal.from_float(float(m['unit_mass_kg']))*u2/((1+u2/Decimal(299792458)**2).sqrt()+1)
                    sampled = Decimal.from_float(float(c['kinetic_per_unit_J']))
                    check(abs(kp-sampled) <= Decimal(8)*Decimal.from_float(sys.float_info.epsilon)*(abs(kp)+abs(sampled)), '90-digit proper-velocity kinetic energy')
            stocks[row['stock']] = new
            e -= packet+dr+x
            gas += packet
            extra += x
            removed += dm
        after = {s['id']: s for s in p['stocks_after']}
        check(stocks == after, 'ordered actions reproduce full candidate inventory')
        ref1 = sum(map(reference, after.values()), F(0))
        check(ref1 == dy(p['reference_after_J']), 'final derived reference')
        check(e == dy(p['thermal_after_J']), 'single external thermal debit')
        check(gas == dy(p['outgoing_energy_J']) and extra == dy(p['other_energy_out_J']), 'outgoing aggregate')
        check(removed == dy(p['outgoing_mass_kg']), 'removed aggregate')
        check(mass(after, models)-mass(before, models)+removed == 0 == dy(p['exact_mass_residual']), 'independent exact mass closure')
        check(e-e0+ref1-ref0+gas+extra == 0 == dy(p['exact_energy_residual']), 'independent exact total energy closure')
        check(impulse == list(map(dy, p['support_impulse_N_s'])), 'exact proper-momentum support quote')
        check(p['fixed_support_work_J'] == p['charge_change_C'] == 0, 'stationary work and neutral charge')
        outer = 'oxide' if dy(after['oxide']['mass_kg']) > 0 else 'metal'
        check(p['exposure_after']['outer_structural_material'] == outer, 'retained structural exposure')
        check(bool(p['exposure_after']['bound_adsorbates']) == (outer == 'oxide' and count(after['water']) > 0), 'adsorbate support exposure')
        projected(after, e, p['final_projected_T_K'])
    for r in data['budget_refusals']:
        available, required = dy(r['available_J']), dy(r['requested_incremental_J'])
        check(0 <= available < required, 'refused insufficient local allowance')
        check(dy(r['plan_scale_limit_ratio']['numerator']) == available and dy(r['plan_scale_limit_ratio']['denominator']) == required, 'refusal ratio is exact')
        check(r['proposed_weights_unchanged'], 'refusal does not resize source plan')
    return data


if __name__ == '__main__':
    normal, sanitized, target = map(Path, sys.argv[1:])
    a = audit(normal)
    first = checks
    b = audit(sanitized)
    check(a == b, 'normal and UBSan deterministic identical receipts')
    result = dict(status='passed_independent_fraction_decimal_fixture_audit',
                  checks=checks, normal_checks=first, ubsan_checks=checks-first-1,
                  deterministic_comparison_checks=1,
                  focused_checks=len(a['checks']), compile_time_contract_checks=3,
                  proofs=len(a['proofs']), budget_refusals=len(a['budget_refusals']),
                  source_identity=a['source'],
                  inputs={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (normal, sanitized)},
                  scope='synthetic proposal arithmetic and refusal evidence; no integrated physics, source rates, native birth, full kinetics or foil calibration')
    target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
