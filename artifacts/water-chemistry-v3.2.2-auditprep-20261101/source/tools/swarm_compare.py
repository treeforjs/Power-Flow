#!/usr/bin/env python3
"""Swarm re-validation of the v3.2 electron set against the measured data, at the agreement bands v3 stated.

usage: swarm_compare.py NEW_SWARM-VALIDATION.json V3_SWARM-VALIDATION.json OUT_PREFIX
writes OUT_PREFIX.json and OUT_PREFIX.txt; exit 0 iff every banded value lies inside its band.

Bands (v3 SCOPE.md section 8, receipt qual-002; deviations MC/measured - 1 of the production set 'v3'):
  (alpha - eta)/N              198-3000 Td   within -11 % / +3 %
  flux drift                    60-800 Td    within -10 % / +9 %
  bulk drift                       <= 350 Td within  -2 % / +15 %
  bulk drift                   400-800 Td    within +19 % / +38 % (the range v3 reported there)
  D_T/mu 64-188 Td and (alpha - eta)/N at 152.5 Td (next to the zero crossing) are reported, not banded, as in v3.
v3's text rounds its extremes (e.g. flux drift -10.4 % at 301.7 Td is summarized as "-10 %"), so each band is widened to
v3's own extreme deviation in its range where needed and rounded outward to a whole percent: the band v3 itself met.
"""
import math
import json
import sys

BANDS = [
    ('alpha_minus_eta_over_N', 'v3', 198.0, 3000.0, -0.11, 0.03, '(alpha - eta)/N'),
    ('drift_velocity', 'v3_flux', 60.0, 800.0, -0.10, 0.09, 'flux drift'),
    ('drift_velocity', 'v3_bulk', 0.0, 350.0, -0.02, 0.15, 'bulk drift'),
    ('drift_velocity', 'v3_bulk', 390.0, 800.0, 0.19, 0.38, 'bulk drift'),
]


def rel_of(row, key):
    k = key + '_rel'
    return row.get(k)


def main(new_p, old_p, out):
    new = json.load(open(new_p))
    old = json.load(open(old_p))
    rows, bad = [], 0
    for table, key, lo, hi, blo, bhi, label in BANDS:
        nt = {round(r['E_over_N_Td'], 3): r for r in new['comparison'][table]}
        ot = {round(r['E_over_N_Td'], 3): r for r in old['comparison'][table]}
        v3vals = [rel_of(r, key) for en, r in ot.items() if lo <= en <= hi and rel_of(r, key) is not None]
        if v3vals and min(v3vals) < blo:
            blo = math.floor(round(100 * min(v3vals), 6)) / 100.0
        if v3vals and max(v3vals) > bhi:
            bhi = math.ceil(round(100 * max(v3vals), 6)) / 100.0
        for en in sorted(nt):
            if not (lo <= en <= hi):
                continue
            a, b = rel_of(nt[en], key), rel_of(ot.get(en, {}), key)
            ok = a is not None and blo <= a <= bhi
            bad += 0 if ok else 1
            rows.append({'quantity': label, 'E_over_N_Td': en, 'band': [blo, bhi], 'v3_rel': b, 'v32_rel': a,
                         'change_points': (a - b) if (a is not None and b is not None) else None, 'inside_band': ok})
    info = []
    for table, key, label in (('DT_over_mu', 'v3', 'D_T/mu (not banded)'), ('alpha_minus_eta_over_N', 'v3', '(alpha - eta)/N (not banded below 198 Td)')):
        nt = {round(r['E_over_N_Td'], 3): r for r in new['comparison'].get(table, [])}
        ot = {round(r['E_over_N_Td'], 3): r for r in old['comparison'].get(table, [])}
        for en in sorted(nt):
            if table == 'alpha_minus_eta_over_N' and en >= 198.0:
                continue
            a, b = rel_of(nt[en], key), rel_of(ot.get(en, {}), key)
            info.append({'quantity': label, 'E_over_N_Td': en, 'v3_rel': b, 'v32_rel': a})
    mc_rows = []
    for en, r in sorted(new['monte_carlo']['v3'].items(), key=lambda x: float(x[0])):
        o = old['monte_carlo']['v3'].get(en, {})
        mc_rows.append({'E_over_N_Td': float(en), 'W_bulk_v3': o.get('W_bulk_m_s'), 'W_bulk_v32': r.get('W_bulk_m_s'),
                        'W_flux_v3': o.get('W_flux_m_s'), 'W_flux_v32': r.get('W_flux_m_s'),
                        'mean_energy_v3': o.get('mean_energy_eV'), 'mean_energy_v32': r.get('mean_energy_eV'),
                        'alpha_minus_eta_v3': o.get('alpha_minus_eta_over_N_m2'), 'alpha_minus_eta_v32': r.get('alpha_minus_eta_over_N_m2')})
    res = {'rule': 'every banded value inside the band v3 stated (SCOPE.md section 8)', 'passed': bad == 0, 'outside_band': bad,
           'banded': rows, 'reported': info, 'monte_carlo_old_new': mc_rows}
    json.dump(res, open(out + '.json', 'w'), indent=1)
    L = ['Swarm re-validation of the v3.2 electron set (production tables) against the measured data, at the bands v3 stated.',
         'Deviation = MC/measured - 1; change = v3.2 - v3 in percentage points. Result: %s (%d values outside their band).' % (
             'PASS' if bad == 0 else 'FAIL', bad), '']
    L.append('  %-18s %8s %16s %9s %9s %9s  %s' % ('quantity', 'E/N(Td)', 'band', 'v3', 'v3.2', 'change', 'inside'))
    for r in rows:
        L.append('  %-18s %8.1f %7.0f%%/%+4.0f%% %+8.1f%% %+8.1f%% %+8.1f  %s' % (r['quantity'], r['E_over_N_Td'], 100 * r['band'][0], 100 * r['band'][1],
                                                                          100 * (r['v3_rel'] or 0), 100 * (r['v32_rel'] or 0), 100 * (r['change_points'] or 0),
                                                                          'yes' if r['inside_band'] else 'NO'))
    L.append('')
    L.append('Reported, not banded (as in v3):')
    for r in info:
        L.append('  %-42s %8.1f %+8.1f%% %+8.1f%%' % (r['quantity'], r['E_over_N_Td'], 100 * (r['v3_rel'] or 0), 100 * (r['v32_rel'] or 0)))
    L.append('')
    L.append('Monte Carlo values, v3 against v3.2 (bulk W, flux W in m/s, mean energy in eV, (alpha - eta)/N in m^2):')
    for r in mc_rows:
        L.append('  %6.0f Td  Wb %.4g / %.4g  Wf %.4g / %.4g  <e> %.4g / %.4g  (a-e)/N %.4g / %.4g' % (
            r['E_over_N_Td'], r['W_bulk_v3'] or 0, r['W_bulk_v32'] or 0, r['W_flux_v3'] or 0, r['W_flux_v32'] or 0,
            r['mean_energy_v3'] or 0, r['mean_energy_v32'] or 0, r['alpha_minus_eta_v3'] or 0, r['alpha_minus_eta_v32'] or 0))
    open(out + '.txt', 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L[:3]))
    return 0 if bad == 0 else 1


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:4]))
