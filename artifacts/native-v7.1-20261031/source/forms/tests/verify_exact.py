"""Independent check of the v7 exact sums: every set in adversarial.txt is
summed with exact rational arithmetic (fractions.Fraction of the binary64
factors) and rounded by CPython's correctly rounded integer division.
Usage: python -I verify_exact.py OUT_DIR/adversarial.txt"""
import struct
import sys
from fractions import Fraction


def f64(h):
    return struct.unpack('<d', struct.pack('<Q', int(h, 16)))[0]


def main(path):
    sets, bad = 0, 0
    with open(path) as f:
        lines = f.read().split('\n')
    i = 0
    while i < len(lines):
        line = lines[i]
        i += 1
        if not line.startswith('set '):
            continue
        _, name, n, value, status = line.split()
        n, status = int(n), int(status)
        total = Fraction(0)
        for k in range(n):
            a, b = lines[i + k].split()
            total += Fraction(f64(a)) * Fraction(f64(b))
        i += n
        try:
            expect = float(total)
            expect_overflow = False
        except OverflowError:
            expect_overflow = True
        got = f64(value)
        if expect_overflow:
            ok = (status & 2) != 0 and got in (float('inf'), float('-inf'))
        else:
            ok = struct.pack('<d', got) == struct.pack('<d', expect) and (status & 2) == 0
        sets += 1
        bad += not ok
        print('%s %s terms=%d value=%r expected=%s' % ('OK ' if ok else 'BAD', name, n, got,
                                                       'overflow' if expect_overflow else repr(expect)))
    print('RESULT exact_rational_verification sets=%d mismatches=%d' % (sets, bad))
    return 1 if bad or not sets else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
