#!/usr/bin/env python3
"""Independent oracle for the arrival owner v1 (exact integer/rational arithmetic).

Re-derives from API-v1.txt, without the C++ code: Philox4x32-10 (Random123
known-answer vectors), the stream layout v1, the exact u01_open conversion, the
exact Bernoulli as the definition U < p with U the 1088-bit prefix of the word
sequence, stochastic rounding, the neutral exact partition (Fraction sums) and the
electron count, then compares with the values dumped by test_core.
usage: oracle.py ORACLE_INPUT_JSON RESULT_JSON
"""
import json
import math
import sys
from fractions import Fraction

M32 = 0xFFFFFFFF
STEP_LIMIT = 1 << 40


def philox4x32_10(ctr, key):
    c0, c1, c2, c3 = ctr
    k0, k1 = key
    for r in range(10):
        if r:
            k0 = (k0 + 0x9E3779B9) & M32
            k1 = (k1 + 0xBB67AE85) & M32
        p0 = 0xD2511F53 * c0
        p1 = 0xCD9E8D57 * c2
        hi0, lo0 = p0 >> 32, p0 & M32
        hi1, lo1 = p1 >> 32, p1 & M32
        c0, c1, c2, c3 = hi1 ^ c1 ^ k0, lo1, hi0 ^ c3 ^ k1, lo0
    return [c0, c1, c2, c3]


def block(seed, step, stream, candidate, purpose, b):
    assert 0 <= step < STEP_LIMIT and 0 <= stream < 1 << 16 and 1 <= purpose <= 255 and 0 <= candidate < 1 << 32 and 0 <= b < 1 << 32
    ctr = [b, candidate, stream | (purpose << 16) | ((step >> 32) << 24), step & M32]
    return philox4x32_10(ctr, [seed & M32, seed >> 32])


def word(key, j):
    w = block(*key, j // 2)
    return (w[1] << 32 | w[0]) if j % 2 == 0 else (w[3] << 32 | w[2])


def u01_open(r):
    return float(Fraction(2 * (r >> 12) + 1, 1 << 53))


def bernoulli(key, p):
    """1 iff U < p, U = sum_j R(j) 2^(-64(j+1)), truncated to 17 words (p has no bits beyond 2^-1088)."""
    P = Fraction(p) * (1 << 1088)
    assert P.denominator == 1
    U = 0
    for j in range(17):
        U = (U << 64) | word(key, j)
    return 1 if U < P.numerator else 0


def stochastic_round(key, x):
    n = math.floor(x)
    return n + bernoulli(key, float(Fraction(x) - n))


def main():
    data = json.load(open(sys.argv[1]))
    checks, failures = 0, []

    def check(ok, what):
        nonlocal checks
        checks += 1
        if not ok and len(failures) < 20:
            failures.append(what)

    kat = [([0, 0, 0, 0], [0, 0], [0x6627e8d5, 0xe169c58d, 0xbc57ac4c, 0x9b00dbd8]),
           ([M32] * 4, [M32] * 2, [0x408f276d, 0x41c83b0e, 0xa20bc7c6, 0x6d5451fd]),
           ([0x243f6a88, 0x85a308d3, 0x13198a2e, 0x03707344], [0xa4093822, 0x299f31d0], [0xd16cfe09, 0x94fdcceb, 0x5001e420, 0x24126ea1])]
    for c, k, o in kat:
        check(philox4x32_10(c, k) == o, 'known-answer vector')
    for rec in data['blocks']:
        key = (int(rec['seed'], 16), rec['step'], rec['stream'], rec['candidate'], rec['purpose'])
        check(block(*key, rec['block']) == rec['out'], 'block words %s' % (key,))
        em = [u01_open(word(key, j)) for j in range(4)]
        nm = [u01_open(word(key, j)) for j in range(6)]
        check(em == [float.fromhex(h) for h in rec['electron_marks']], 'electron marks')
        check(nm == [float.fromhex(h) for h in rec['neutral_marks']], 'neutral marks')
        check(all(0 < u < 1 for u in em + nm), 'open interval')
        check(stochastic_round(key, float.fromhex(rec['x'])) == rec['sr'], 'stochastic rounding')
    seeds = data['seeds']
    for c in data['partitions']:
        D, w0 = float.fromhex(c['D']), float.fromhex(c['w0'])
        x = D / w0
        check(x == float.fromhex(c['x']), 'ratio D/w0')
        rounded = stochastic_round((seeds['neutral'], c['step'], seeds['neutral_stream'], 0, 1), x)
        check(rounded == c['rounded'], 'neutral SR count')
        m, e = math.frexp(D)   # D = m 2^e, m in [0.5, 1)
        if D >= 2.0 ** -1022:
            M, E = int(m * 2 ** 53), e - 53
        else:
            M, E = int(Fraction(D) * (1 << 1074)), -1074
        check(Fraction(M) * Fraction(2) ** E == Fraction(D), 'quantum decomposition')
        births = min(max(1, rounded), M)
        check(births == c['births'], 'neutral births')
        q, r = divmod(M, births)
        weights = [float(Fraction(q + (1 if k < r else 0)) * Fraction(2) ** E) for k in range(births)]
        check([float.fromhex(h) for h in c['weights']] == weights[:c['weights_listed']], 'leaf weights')
        check(sum(Fraction(w) for w in weights) == Fraction(D), 'exact leaf sum equals the debit')
    for c in data['electrons']:
        N, w0 = float.fromhex(c['N']), float.fromhex(c['w0'])
        n = stochastic_round((seeds['electron'], c['step'], seeds['electron_stream'], 0, 1), N / w0)
        check(n == c['births'], 'electron births')
    result = {'oracle': 'arrival-owner-v1 independent exact oracle', 'checks': checks, 'failures': failures, 'passed': not failures,
              'blocks': len(data['blocks']), 'partitions': len(data['partitions']), 'electrons': len(data['electrons'])}
    json.dump(result, open(sys.argv[2], 'w'), indent=1)
    print('oracle checks=%d %s' % (checks, 'PASS' if not failures else 'FAIL %s' % failures[:3]))
    return 0 if not failures else 1


if __name__ == '__main__':
    sys.exit(main())
