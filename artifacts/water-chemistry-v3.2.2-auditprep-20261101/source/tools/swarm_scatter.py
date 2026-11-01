#!/usr/bin/env python3
"""Monte Carlo scatter of the swarm values at chosen E/N (diagnostic for the re-validation; no data or band change).

usage: swarm_scatter.py DIR OUT.json [REFERENCE_SWARM-VALIDATION.json]
DIR (comma-separated directories) holds swarm_mc outputs named <set>-<seed>.json (sets v31, v32). For every E/N point: mean and standard deviation
over the seeds of bulk drift, flux drift and (alpha - eta)/N per set, and the v3.2 - v3.1 shift of the means in units of
the combined standard error. With the reference file: per banded point the measurement, the band v3 met (as in
tools/swarm_compare.py) and both sets' deviations mean +- spread (sample standard deviation over the seeds).
"""
import glob
import json
import math
import os
import sys


def main(d, out, ref=None):
    data = {}
    seen = set()
    for f in sorted(f_ for dd in d.split(',') for f_ in glob.glob(os.path.join(dd, 'v3[12]-*.json'))):
        name = os.path.basename(f)[:-5]
        st, seed = name.split('-', 1)
        if (st, seed) in seen:   # the same seed run in two places is one sample
            continue
        try:
            j = json.load(open(f))
        except ValueError:   # a run still in progress
            continue
        seen.add((st, seed))
        for p in j['points']:
            key = round(p['E_over_N_Td'], 3)
            data.setdefault(key, {}).setdefault(st, []).append(p)
    res = []
    for en in sorted(data):
        row = {'E_over_N_Td': en}
        for q in ('W_bulk_m_s', 'W_flux_m_s', 'alpha_minus_eta_over_N_m2'):
            for st in ('v31', 'v32'):
                v = [p[q] for p in data[en].get(st, [])]
                if not v:
                    continue
                m = sum(v) / len(v)
                sd = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else None
                row['%s_%s' % (st, q)] = {'n': len(v), 'mean': m, 'sd': sd, 'sd_rel': (sd / abs(m)) if (sd is not None and m) else None}
            a, b = row.get('v31_%s' % q), row.get('v32_%s' % q)
            if a and b and a['sd'] is not None and b['sd'] is not None:
                se = math.sqrt(a['sd'] ** 2 / a['n'] + b['sd'] ** 2 / b['n'])
                row['shift_%s' % q] = {'rel': b['mean'] / a['mean'] - 1, 'in_standard_errors': (b['mean'] - a['mean']) / se if se > 0 else None}
        res.append(row)
    report = []
    if ref:
        import swarm_compare as sc
        rj = json.load(open(ref))
        qmap = {'drift_velocity': {'v3_bulk': 'W_bulk_m_s', 'v3_flux': 'W_flux_m_s'}, 'alpha_minus_eta_over_N': {'v3': 'alpha_minus_eta_over_N_m2'}}
        meas_key = {'drift_velocity': 'measured_W_m_s', 'alpha_minus_eta_over_N': 'measured_m2'}
        for table, key, lo, hi, blo, bhi, label in sc.BANDS:
            ot = {round(r['E_over_N_Td'], 3): r for r in rj['comparison'][table]}
            v3vals = [r[key + '_rel'] for en, r in ot.items() if lo <= en <= hi and r.get(key + '_rel') is not None]
            if v3vals and min(v3vals) < blo:
                blo = math.floor(round(100 * min(v3vals), 6)) / 100.0
            if v3vals and max(v3vals) > bhi:
                bhi = math.ceil(round(100 * max(v3vals), 6)) / 100.0
            q = qmap[table][key]
            for en in sorted(data):
                near = [e_ for e_ in ot if abs(e_ - en) < 0.06]   # the MC runs at the printed (rounded) measured E/N
                if not near or not (lo <= en <= hi):
                    continue
                meas = ot[near[0]][meas_key[table]]
                row = {'quantity': label, 'E_over_N_Td': en, 'measured': meas, 'band': [blo, bhi]}
                for st in ('v31', 'v32'):
                    v = [p[q] / meas - 1 for p in data[en].get(st, [])]
                    if v:
                        m = sum(v) / len(v)
                        sd = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else 0.0
                        row[st] = {'n': len(v), 'mean_rel': m, 'sd_rel': sd, 'inside_band_mean': blo <= m <= bhi}
                report.append(row)
    json.dump({'points': res, 'banded_report': report}, open(out, 'w'), indent=1)
    if report:
        print('%-16s %8s %11s %13s %22s %22s' % ('quantity', 'E/N(Td)', 'measured', 'band', 'v3 (v3.1 tables)', 'v3.2'))
        for r in report:
            a, b = r.get('v31', {}), r.get('v32', {})
            print('%-16s %8.1f %11.4g %6.0f%%/%+3.0f%% %+8.1f%% +- %4.1f%% (n %d) %+8.1f%% +- %4.1f%% (n %d)' % (
                r['quantity'], r['E_over_N_Td'], r['measured'], 100 * r['band'][0], 100 * r['band'][1], 100 * a.get('mean_rel', 0), 100 * a.get('sd_rel', 0),
                a.get('n', 0), 100 * b.get('mean_rel', 0), 100 * b.get('sd_rel', 0), b.get('n', 0)))
    for r in res:
        parts = []
        for q, lab in (('W_bulk_m_s', 'bulk'), ('W_flux_m_s', 'flux'), ('alpha_minus_eta_over_N_m2', '(a-e)/N')):
            a, b, s = r.get('v31_' + q), r.get('v32_' + q), r.get('shift_' + q)
            if a and b and s:
                parts.append('%s: sd %.2f%%/%.2f%% shift %+.2f%% (%.1f se)' % (lab, 100 * (a['sd_rel'] or 0), 100 * (b['sd_rel'] or 0), 100 * s['rel'],
                                                                         s['in_standard_errors'] or 0))
        print('%6.0f Td  %s' % (r['E_over_N_Td'], ';  '.join(parts)))


if __name__ == '__main__':
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    main(*sys.argv[1:4])
