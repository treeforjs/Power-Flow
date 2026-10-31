#!/usr/bin/env python3
"""Independent 50-digit Decimal evaluation of BEB and RBEB for the H2O orbitals.

Reads the orbital constants from the frozen outer3 dataset (bytes pinned by
SHA-256), evaluates
  BEB  (Kim & Rudd 1994, Q = 1) and
  RBEB (Kim, Santos & Parente 2000, alpha^4 = (2R/mc^2)^2)
per orbital and the parent-weighted sums, and writes them as JSON. The C++
test values are compared against this file by run_coverage_v2_checks.py.
Inputs are the exact binary64 values the C++ code receives (constants parsed
from the dataset JSON, mc^2 = 510998.95069 as binary64), converted exactly to
Decimal; near threshold the cross sections are proportional to T - B, so the
decimal spelling of B would otherwise differ by the binary64 rounding of B
amplified by T/(T-B).
Usage: oracle_beb_rbeb.py OUTER3_JSON OUTPUT_JSON
"""
import hashlib
import json
import sys
from decimal import Decimal as D, getcontext

getcontext().prec = 50
OUTER3_SHA = 'c436020d3f2d07902fb059cfc0a269b2e914a3cbb2789300093d3631df896551'
MC2 = D(510998.95069)  # exact binary64 value
PI = D('3.14159265358979323846264338327950288419716939937510582097494')
ENERGIES = ['13', '20', '50', '100', '1000', '1500', '2000', '3000', '5000', '8000', '10000', '20000', '50000',
            '100000', '200000', '500000', '1000000', '2000000', '5000000', '10000000', '100000000', '1000000000']


def beta2(x):
    return x * (x + 2) / ((1 + x) * (1 + x))


def beb(o, a0, R, T):
    B, U, N = D(float(o['B_eV'])), D(float(o['U_eV'])), D(float(o['N']))
    if T <= B:
        return D(0)
    t, u = T / B, U / B
    S = 4 * PI * a0 * a0 * N * (R / B) ** 2
    lt = t.ln()
    return S / (t + u + 1) * (lt / 2 * (1 - 1 / (t * t)) + (1 - 1 / t - lt / (t + 1)))


def rbeb(o, a0, R, T):
    B, U, N = D(float(o['B_eV'])), D(float(o['U_eV'])), D(float(o['N']))
    if T <= B:
        return D(0)
    t = T / B
    tp, bp, up = T / MC2, B / MC2, U / MC2
    bt2, bb2, bu2 = beta2(tp), beta2(bp), beta2(up)
    alpha4 = (2 * R / MC2) ** 2
    P = 4 * PI * a0 * a0 * alpha4 * N / ((bt2 + bu2 + bb2) * 2 * bp)
    L = (bt2 / (1 - bt2)).ln() - bt2 - (2 * bp).ln()
    h = 1 + tp / 2
    return P * (L / 2 * (1 - 1 / (t * t)) + 1 - 1 / t - t.ln() / (t + 1) * (1 + 2 * tp) / (h * h)
                + bp * bp / (h * h) * (t - 1) / 2)


def main():
    raw = open(sys.argv[1], 'rb').read()
    if hashlib.sha256(raw).hexdigest() != OUTER3_SHA:
        raise SystemExit('outer3 dataset differs')
    data = json.loads(raw)
    ch = next(c for c in data['channels'] if c['id'] == 'H2O:conditional_outer3:parent')
    a0, R = D(float(ch['bohr_radius_m'])), D(float(ch['rydberg_eV']))
    rows = []
    for e in ENERGIES:
        T = D(e)
        per = []
        for o, wgt in zip(ch['orbitals'], ch['orbital_weights']):
            per.append({'label': o['label'], 'weight': wgt, 'beb_m2': str(beb(o, a0, R, T)), 'rbeb_m2': str(rbeb(o, a0, R, T))})
        wsum = lambda key: sum(D(float(p['weight'])) * D(p[key]) for p in per)
        tsum = lambda key: sum(D(p[key]) for p in per)
        rows.append({'energy_eV': float(e), 'orbitals': per,
                     'beb_parent_m2': str(wsum('beb_m2')), 'rbeb_parent_m2': str(wsum('rbeb_m2')),
                     'beb_valence_m2': str(tsum('beb_m2')), 'rbeb_valence_m2': str(tsum('rbeb_m2'))})
    json.dump({'schema': 'beb-rbeb-decimal-oracle-v1', 'precision_digits': 50, 'electron_rest_energy_eV_binary64_exact': str(MC2),
               'outer3_sha256': OUTER3_SHA, 'rows': rows}, open(sys.argv[2], 'w'), indent=1)


if __name__ == '__main__':
    main()
