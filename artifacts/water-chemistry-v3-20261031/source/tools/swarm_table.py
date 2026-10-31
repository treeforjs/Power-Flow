#!/usr/bin/env python3
"""Swarm validation table: Monte Carlo swarm parameters of each cross-section set vs measured water-vapour data.

usage: swarm_table.py SWARM_DIR OUT_PREFIX     (writes OUT_PREFIX.json and OUT_PREFIX.txt)
Measured data (data/extracted, LXCat LAPLACE): Hasegawa et al. mobility x N (-> W = muN E/N) and
(alpha - eta)/N; Crompton et al. 1965 D_T/mu. Monte Carlo values are interpolated to the measured
E/N log-linearly in E/N between computed points (no extrapolation).
"""
import glob
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_sets(d):
    res = {}
    for f in sorted(glob.glob(os.path.join(d, '*-[0-9].json'))):
        name = os.path.basename(f).rsplit('-', 1)[0]
        try:
            j = json.load(open(f))
        except Exception:
            continue
        for p in j['points']:
            res.setdefault(name, {})[p['E_over_N_Td']] = p
    return res


def interp(points, en, key):
    xs = sorted(points)
    for a, b in zip(xs, xs[1:]):
        if a <= en <= b:
            ya, yb = points[a][key], points[b][key]
            t = (math.log(en) - math.log(a)) / (math.log(b) - math.log(a))
            return ya + (yb - ya) * t
    return None


def main(d, out):
    res = load_sets(d)
    ext = lambda t: json.load(open(os.path.join(ROOT, 'data/extracted', t + '.json')))
    mu = ext('swarm-laplace-mobility-N-hasegawa')
    ae = ext('swarm-laplace-alpha-minus-eta-over-N-hasegawa')
    ce = ext('swarm-laplace-characteristic-energy-crompton1965')
    order = [n for n in ('v3', 'v3-iso', 'lxcat-muroranit', 'lxcat-triniti', 'lxcat-hayashi', 'lxcat-phelps', 'lxcat-morgan') if n in res]
    rows = {'drift_velocity': [], 'alpha_minus_eta_over_N': [], 'DT_over_mu': []}
    for en, m in zip(mu['x'], mu['y']):
        W = m * en * 1e-21
        r = {'E_over_N_Td': en, 'measured_W_m_s': W}
        for n in order:
            for key, tag in (('W_bulk_m_s', 'bulk'), ('W_flux_m_s', 'flux')):
                v = interp(res[n], en, key)
                r['%s_%s' % (n, tag)] = v
                r['%s_%s_rel' % (n, tag)] = (v / W - 1) if v else None
        rows['drift_velocity'].append(r)
    for en, m in zip(ae['x'], ae['y']):
        r = {'E_over_N_Td': en, 'measured_m2': m}
        for n in order:
            v = interp(res[n], en, 'alpha_minus_eta_over_N_m2')
            r[n] = v
            r[n + '_rel'] = (v / m - 1) if v is not None else None
            r[n + '_alpha'] = interp(res[n], en, 'alpha_over_N_m2')
            r[n + '_eta'] = interp(res[n], en, 'eta_over_N_m2')
        rows['alpha_minus_eta_over_N'].append(r)
    for en, m in zip(ce['x'], ce['y']):
        r = {'E_over_N_Td': en, 'measured_eV': m}
        for n in order:
            v = interp(res[n], en, 'DT_over_mu_eV')
            r[n] = v
            r[n + '_rel'] = (v / m - 1) if v is not None else None
        rows['DT_over_mu'].append(r)
    mc = {n: {str(k): v for k, v in sorted(res[n].items())} for n in order}
    json.dump({'sets': order, 'measured_sources': {'W': mu['locator'], 'alpha_minus_eta': ae['locator'], 'DT_over_mu': ce['locator']},
               'comparison': rows, 'monte_carlo': mc}, open(out + '.json', 'w'), indent=1)
    L = []
    pct = lambda x: '   -   ' if x is None else '%+6.1f%%' % (100 * x)
    L.append('Swarm validation of the chemistry v3 electron + H2O set (0D Monte Carlo, tests/swarm_mc.cpp; N = 3.22e22 m^-3, T_gas = 300 K)')
    L.append('Sets: v3 = production tables; v3-iso = same with isotropic elastic (sigma_m); lxcat-* = complete LXCat sets run through the same')
    L.append('Monte Carlo (energy loss only; isotropic; no data above their last tabulated energy). Deviations are MC/measured - 1.')
    L.append('')
    L.append('Drift velocity, measured W = (mu N)(E/N), Hasegawa et al. (LXCat LAPLACE). bulk = d<z>/dt, flux = <v_z>.')
    L.append('  E/N(Td)  W_meas(m/s)  ' + ''.join('%-20s' % (n.replace('lxcat-', '') + ' bulk/flux') for n in order))
    for r in rows['drift_velocity']:
        L.append('  %7.1f  %10.4g  ' % (r['E_over_N_Td'], r['measured_W_m_s']) + ''.join(
            '%-20s' % ('%s/%s' % (pct(r['%s_bulk_rel' % n]).strip(), pct(r['%s_flux_rel' % n]).strip())) for n in order))
    L.append('')
    L.append('Effective ionization (alpha - eta)/N (m^2), Hasegawa et al. (LXCat LAPLACE).')
    L.append('  E/N(Td)  measured    ' + ''.join('%14s' % n.replace('lxcat-', '') for n in order) + '   v3 alpha/N   v3 eta/N')
    for r in rows['alpha_minus_eta_over_N']:
        L.append('  %7.1f  %10.4g  ' % (r['E_over_N_Td'], r['measured_m2']) + ''.join('%14s' % pct(r[n + '_rel']) for n in order) +
                 '   %10.3g %10.3g' % (r.get('v3_alpha') or 0, r.get('v3_eta') or 0))
    L.append('')
    L.append('Transverse characteristic energy D_T/mu (eV), Crompton et al. 1965 (LXCat LAPLACE).')
    L.append('  E/N(Td)  measured    ' + ''.join('%14s' % n.replace('lxcat-', '') for n in order))
    for r in rows['DT_over_mu']:
        L.append('  %7.1f  %10.4g  ' % (r['E_over_N_Td'], r['measured_eV']) + ''.join('%14s' % pct(r[n + '_rel']) for n in order))
    L.append('')
    L.append('Monte Carlo values of v3 (bulk W, flux W, mean energy, alpha/N, eta/N):')
    for k, p in sorted(res.get('v3', {}).items()):
        L.append('  %6g Td  Wb %.4g  Wf %.4g m/s  <e> %.4g eV  alpha/N %.4g  eta/N %.4g m^2  real events %d' % (
            k, p['W_bulk_m_s'], p['W_flux_m_s'], p['mean_energy_eV'], p['alpha_over_N_m2'], p['eta_over_N_m2'], p['real_events']))
    open(out + '.txt', 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
