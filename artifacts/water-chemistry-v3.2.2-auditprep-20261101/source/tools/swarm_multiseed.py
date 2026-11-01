#!/usr/bin/env python3
"""Multi-seed swarm records at the measured E/N (no interpolation) and the v3.2 re-validation criterion.

usage: swarm_multiseed.py jobs SET BLOB OUTDIR SEED[,SEED...] JOBS.json   one swarm_mc job per seed over every measured E/N
       swarm_multiseed.py summary DIR SET OUT.json                       mean, SD, n of the deviations MC/measured - 1
       swarm_multiseed.py compare NEW.json REF.json OLD_SINGLE_SEED.json OUT_PREFIX
Measured data (data/extracted, LXCat LAPLACE): Hasegawa et al. mobility x N (W = muN E/N) and (alpha - eta)/N,
Crompton et al. 1965 D_T/mu. Monte Carlo: tests/swarm_mc.cpp, N = 3.22e22 m^-3, 300 K, 2000 electrons per seed.
Criterion (observed worker criterion of 2026-10-31; authorization not established here): each v3.2 mean within the v3 multi-seed mean +- max(2 SD_v3, 1 percentage point),
for the drift velocities (bulk and flux) at every measured E/N and (alpha - eta)/N from 198 Td up (as v3 banded them);
D_T/mu and (alpha - eta)/N at 152.5 Td (next to the zero crossing) are reported, not judged, as in v3.
"""
import glob
import json
import math
import os
import sys
import importlib.util
from pathlib import Path

R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location('swarm_validation', Path(__file__).with_name('swarm_validation.py'))
_validation = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_validation)


def contract():
    return _validation.Contract(Path(__file__).with_name('swarm-contract.json'), R)


def measured():
    ext = lambda t: json.load(open(os.path.join(R, 'data/extracted', t + '.json')))
    mu = ext('swarm-laplace-mobility-N-hasegawa')
    ae = ext('swarm-laplace-alpha-minus-eta-over-N-hasegawa')
    ce = ext('swarm-laplace-characteristic-energy-crompton1965')
    m = {'W_bulk_m_s': {en: v * en * 1e-21 for en, v in zip(mu['x'], mu['y'])}}
    m['W_flux_m_s'] = dict(m['W_bulk_m_s'])
    m['alpha_minus_eta_over_N_m2'] = dict(zip(ae['x'], ae['y']))
    m['DT_over_mu_eV'] = dict(zip(ce['x'], ce['y']))
    return m


def judged(q, en):
    if q in ('W_bulk_m_s', 'W_flux_m_s'):
        return True
    if q == 'alpha_minus_eta_over_N_m2':
        return en >= 198.0
    return False


def jobs(st, blob, outdir, seeds, jp):
    c = contract()
    _validation.need(st in c.data['sets'], 'unknown set')
    _validation.need([int(x) for x in seeds.split(',')] == c.data['seeds'], 'unexpected seed plan')
    _validation.need(__import__('hashlib').sha256(Path(blob).read_bytes()).hexdigest() == c.data['sets'][st]['blob_sha256'], 'job blob mismatch')
    m = c.measured
    ens = c.energies
    js = []
    for seed in seeds.split(','):
        out = os.path.join(outdir, '%s-%s.json' % (st, seed))
        js.append({'name': '%s-%s' % (st, seed), 'cwd': R, 'log': out[:-5] + '.log', 'output': out,
                   'argv': [os.path.join(R, 'build/bin/swarm_mc'), blob, out, ','.join(repr(e) for e in ens), '3.22e22', '300', '2000', seed]})
    json.dump(js, open(jp, 'w'), indent=1)
    print(len(js), 'jobs,', len(ens), 'E/N points each')


def summary(d, st, out):
    result = _validation.summarize(d, st, contract())
    _validation.atomic_json(out, result)
    print(st, len(result['seeds']), 'complete runs,', len(result['rows']), 'rows')


LABEL = {'W_bulk_m_s': 'bulk drift', 'W_flux_m_s': 'flux drift', 'alpha_minus_eta_over_N_m2': '(alpha - eta)/N', 'DT_over_mu_eV': 'D_T/mu'}
OLDKEY = {'W_bulk_m_s': ('drift_velocity', 'v3_bulk_rel'), 'W_flux_m_s': ('drift_velocity', 'v3_flux_rel'),
          'alpha_minus_eta_over_N_m2': ('alpha_minus_eta_over_N', 'v3_rel'), 'DT_over_mu_eV': ('DT_over_mu', 'v3_rel')}


def compare(new_p, ref_p, old_p, outp):
    c = contract()
    new, _ = _validation.stable_json(new_p)
    ref, _ = _validation.stable_json(ref_p)
    old, old_hash = _validation.stable_json(old_p)
    _validation.need(old_hash == c.data['old_reference_sha256'], 'historical reference source mismatch')
    _validation.validate_summary(new, 'v32', c)
    _validation.validate_summary(ref, 'v3', c)
    rk = {(r['quantity'], r['E_over_N_Td']): r for r in ref['rows']}
    rows, bad = [], 0
    for r in new['rows']:
        k = (r['quantity'], r['E_over_N_Td'])
        f = rk.get(k)
        _validation.need(f is not None, 'missing reference row')
        tol = max(2 * f['sd_rel'], 0.01)
        diff = r['mean_rel'] - f['mean_rel']
        j = judged(*k)
        ok = abs(diff) <= tol
        if j and not ok:
            bad += 1
        tbl, key = OLDKEY[r['quantity']]
        o = next((x.get(key) for x in old['comparison'][tbl] if abs(x['E_over_N_Td'] - r['E_over_N_Td']) < 1e-6), None)
        rows.append({'quantity': LABEL[r['quantity']], 'E_over_N_Td': r['E_over_N_Td'], 'measured': r['measured'], 'judged': j,
                     'v3_single_seed_interpolated_rel': o, 'v3_mean_rel': f['mean_rel'], 'v3_sd_rel': f['sd_rel'], 'v3_n': f['n'],
                     'v32_mean_rel': r['mean_rel'], 'v32_sd_rel': r['sd_rel'], 'v32_n': r['n'], 'tolerance': tol, 'difference': diff,
                     'pass': ok if j else None})
    _validation.need(len(rows) == len(c.rows) and sum(1 for r in rows if r['judged']) == c.judged and c.judged > 0, 'incomplete or zero judged rows')
    res = {'raw_input_files': {'v3': ref['input_files'], 'v32': new['input_files']}, 'contract_sha256': c.sha256, 'complete': True, 'scientific_criterion_status': c.data['scientific_criterion_status'], 'criterion': 'each v3.2 mean within the v3 multi-seed mean +- max(2 SD_v3, 1 percentage point) (observed worker criterion of 2026-10-31; authorization not established here)',
           'v3_reference': {'blobs': ref['blobs'], 'seeds': ref['seeds']}, 'v32': {'blobs': new['blobs'], 'seeds': new['seeds']},
           'judged': sum(1 for r in rows if r['judged']), 'failing': bad, 'passed': bad == 0, 'rows': rows}
    _validation.validate_comparison(res, new, ref, c)
    _validation.atomic_json(outp + '.json', res)
    L = ['Swarm re-validation of the v3.2 electron set: multi-seed Monte Carlo directly at the measured E/N.',
         'Criterion (observed worker criterion of 2026-10-31; authorization not established here): each v3.2 mean within the v3 multi-seed mean +- max(2 SD_v3, 1 percentage point).',
         'Result: %s (%d of %d judged values outside the tolerance). v3 = v3.1 tables %s; v3.2 tables %s; seeds %s.' % (
             'PASS' if bad == 0 else 'FAIL', bad, res['judged'], ','.join(b[:12] for b in ref['blobs']), ','.join(b[:12] for b in new['blobs']),
             ','.join(str(s) for s in new['seeds'])),
         'Deviation = MC/measured - 1. "old v3 record" = the single-seed value v3 reported (interpolated between computed E/N points).', '']
    L.append('  %-16s %8s %11s %13s %18s %18s %8s %7s' % ('quantity', 'E/N(Td)', 'measured', 'old v3 record', 'v3 mean +- SD', 'v3.2 mean +- SD',
                                                        'diff', 'result'))
    for r in rows:
        L.append('  %-16s %8.1f %11.4g %12s %+9.1f%% +-%4.1f%% %+9.1f%% +-%4.1f%% %+7.2f  %s' % (
            r['quantity'], r['E_over_N_Td'], r['measured'],
            ('%+.1f%%' % (100 * r['v3_single_seed_interpolated_rel'])) if r['v3_single_seed_interpolated_rel'] is not None else '-',
            100 * r['v3_mean_rel'], 100 * r['v3_sd_rel'], 100 * r['v32_mean_rel'], 100 * r['v32_sd_rel'], 100 * r['difference'],
            ('pass' if r['pass'] else 'FAIL') if r['judged'] else 'report'))
    Path(outp + '.txt').write_text('\n'.join(L) + '\n')
    print('\n'.join(L[:3]))
    return 0 if bad == 0 else 1


if __name__ == '__main__':
    a = sys.argv[1:]
    if a[0] == 'jobs':
        jobs(*a[1:6])
    elif a[0] == 'summary':
        summary(*a[1:4])
    elif a[0] == 'compare':
        sys.exit(compare(*a[1:5]))
