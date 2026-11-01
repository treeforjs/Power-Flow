#!/usr/bin/env python3
"""Checks of the v3.2 model estimates against measured analogues (documentation; RESULTS.txt and SCOPE.md section 12).

usage: validation_checks.py OUT.json
  1. free-collision breakup estimate (source type free_collision, as sampled for O+, OH+, H2O+, H3O+, H2+, O-, OH-):
     applied to H- detachment (I = EA(H)) against the Simon Wedlund 2019 sigma-10 table, and to He+ electron loss
     (I = 54.4 eV) against the Simon Wedlund 2019 He+ sigma12 fit;
  2. target ionization at 6.6 keV O+: proton sigma11 and He+ sigma11 at the same speed against the measured total
     positive fragment-ion production of Juhasz et al. 2019.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_tables as bt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RK = '/home/treefor/projects/runtimes/h2o-collision-data-20261029/data/collisions/h2o-gas-20261029.json'
U_EV = 931494103.72
ME = 510998.95069


def ext(t):
    return json.load(open(os.path.join(ROOT, 'data/extracted', t + '.json')))


def wfit(t, v_m_s):
    p = [float(x) for x in ext(t)['y']]
    lv = math.log10(v_m_s)
    return 10 ** sum(c * lv ** j for j, c in enumerate(p))


def main(out):
    common = json.load(open(os.path.join(ROOT, 'tests/specs/common-species.json')))
    sp = bt.species_table(common, {})
    mass = {x['name']: x['mass_eV'] for x in sp}
    mass['He+'] = 4.002602 * U_EV - ME   # He+ (rest energy; binding energies negligible here)
    ref = {x['name']: x['ref_eV'] for x in sp}
    ctx = {'files': {}, 'mass_of': lambda n: mass[n]}
    efree = {"type": "sum_of", "hold": True, "parts": [{"type": "raikou", "file": RK, "channel": "H2O:munoz2007:elastic_integral"},
                                                       {"type": "raikou", "file": RK, "channel": "H2O:munoz2007:inelastic_integral"}]}
    res = {'free_collision_vs_H_minus_sigma_m10': [], 'free_collision_vs_He_plus_sigma12': [], 'ionization_6p6keV_O_plus': {}}
    I_H = ref['H'] - ref['H-']
    fc_h = bt.make_source({"type": "free_collision", "projectile": "H-", "ref": efree, "I_eV": I_H, "lo": 1e3, "hi": 1e7}, ctx)
    t = ext('ion-wedlund2019-hminus-h2o-one-electron-loss-sigmam10')
    for x, y in zip(t['x'], t['y']):
        E = x * mass['H-'] / U_EV
        if E < 1e3:
            continue
        est = fc_h(E)
        res['free_collision_vs_H_minus_sigma_m10'].append({'E_eV': E, 'measured_m2': y, 'estimate_m2': est, 'ratio': est / y})
    fc_he = bt.make_source({"type": "free_collision", "projectile": "He+", "ref": efree, "I_eV": 54.4178, "lo": 1e3, "hi": 1e8}, ctx)
    for v_km in (3000, 4000, 5000, 6000, 8000, 10000):
        v = v_km * 1e3
        b = v / bt.C_M_S
        E = mass['He+'] * (1 / math.sqrt(1 - b * b) - 1)
        meas = wfit('ion-wedlund2019-fitcoeffs-heplus-h2o-sigma12', v)
        est = fc_he(E)
        res['free_collision_vs_He_plus_sigma12'].append({'v_km_s': v_km, 'E_eV': E, 'measured_fit_m2': meas, 'estimate_m2': est,
                                                         'ratio': est / meas if meas > 0 else None})
    E = 6600.0
    v = bt.speed_tf(E, mass['O+'])
    j = ext('ion-juhasz2019-oplus-h2o-cation-fragment-total-6p6kev')
    res['ionization_6p6keV_O_plus'] = {'v_m_s': v, 'measured_total_positive_fragment_ion_production_m2': j['y'][0], 'measured_uncertainty_rel': 0.5,
                                       'proton_sigma11_equal_speed_m2': wfit('ion-wedlund2019-fitcoeffs-hplus-h2o-sigma11', v),
                                       'He_plus_sigma11_equal_speed_m2': wfit('ion-wedlund2019-fitcoeffs-heplus-h2o-sigma11', v)}
    json.dump(res, open(out, 'w'), indent=1)
    r1 = [x['ratio'] for x in res['free_collision_vs_H_minus_sigma_m10']]
    r2 = [x['ratio'] for x in res['free_collision_vs_He_plus_sigma12'] if x['ratio']]
    print('free collision / H- sigma-10 (1-130 keV): %.2f - %.2f' % (min(r1), max(r1)))
    print('free collision / He+ sigma12 (3000-10000 km/s): %.2f - %.2f' % (min(r2), max(r2)))
    o = res['ionization_6p6keV_O_plus']
    print('6.6 keV O+: measured %.2g, proton %.2g, He+ %.2g m2' % (o['measured_total_positive_fragment_ion_production_m2'],
                                                                  o['proton_sigma11_equal_speed_m2'], o['He_plus_sigma11_equal_speed_m2']))


if __name__ == '__main__':
    main(sys.argv[1])
