#!/usr/bin/env python3
"""Build a chemistry v3 table blob ("PFCHEMV3") and its manifest from a spec.

usage: build_tables.py SPEC.json OUT_PREFIX
writes OUT_PREFIX.blob, OUT_PREFIX-manifest.json, OUT_PREFIX-checks.json

The spec names every process, its kind, products, provenance tag and the
cross-section source (tabulated sources are read verbatim from the retained
files; analytic models are sampled and their interpolation error is recorded).
Species masses and reference energies come from ATcT 0 K formation enthalpies
plus CODATA/NIST atomic masses (see species_table). Every number written to the
blob is a binary64 produced here by documented arithmetic; the blob SHA-256 is
the table identity.
"""
import hashlib
import json
import math
import os
import struct
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EV_PER_KJ_MOL = 1.0 / 96.48533212331001   # e * N_A / 1000 exactly (SI 2019)
C_M_S = 299792458.0

KIND = {'ELASTIC': 1, 'EXCITE': 2, 'DISSOC': 3, 'IONIZE': 4, 'CAPTURE': 5, 'CHARGE_X': 6, 'ION_IONIZE': 7, 'BREAKUP': 8}
ANG = {'ISO': 0, 'SCREENED': 1, 'FORWARD': 2}
DEST = {'SURVIVOR': 1, 'BIRTH': 2, 'RESERVOIR': 3, 'GAS': 4}
SINK = {'NONE': 0, 'ROT': 1, 'VIB': 2, 'ELEC': 3, 'ION_INTERNAL': 4, 'FRAGMENT': 5}
TAG = {'measured': 1, 'recommended': 2, 'calculated': 3, 'model_derived': 4}
COVER = {'ZERO': 0, 'DATA': 1, 'MODEL': 2, 'UNRESOLVED': 3}
SPKIND = {'electron': 0, 'cation': 1, 'anion': 2, 'neutral': 3}
NUCLEI = ['H', 'O', 'Fe', 'Cr', 'Ni']
NONE32 = 0xFFFFFFFF
MAX_OUT, MAX_ORB, MAX_COVER = 6, 4, 6


def sha_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def rel(p):
    return p if os.path.isabs(p) else os.path.join(ROOT, p)


# ------------------------------------------------------------ species
def species_table(spec, files):
    th = spec['thermo']
    atct = json.load(open(rel(th['atct_file'])))
    files[th['atct_file']] = sha_file(rel(th['atct_file']))
    u_eV = th['atomic_mass_constant_eV']
    me = th['electron_rest_energy_eV']
    atom_R = {}
    atom_m = {}
    for a, info in th['atoms'].items():
        atom_m[a] = info['mass_u'] * u_eV
        if 'atct' in info:
            atom_R[a] = float(atct[info['atct']]['DHf0_kJ_mol']) * EV_PER_KJ_MOL
        else:
            atom_R[a] = info['ref_eV']
    out = []
    for s in spec['species']:
        rec = dict(s)
        if s['kind'] == 'electron':
            rec['mass_eV'] = me
            rec['ref_eV'] = 0.0
        else:
            if 'atct' in s:
                R = float(atct[s['atct']]['DHf0_kJ_mol']) * EV_PER_KJ_MOL
                rec['ref_source'] = 'ATcT v1.202 %s DHf(0 K) = %s kJ/mol' % (s['atct'], atct[s['atct']]['DHf0_kJ_mol'])
            else:
                R = s['ref_eV']
                rec['ref_source'] = s.get('ref_source', 'spec')
            base = 0.0
            for a, n in s['atoms'].items():
                base += n * (atom_m[a] - atom_R[a])
            rec['mass_eV'] = base - s['charge'] * me + R
            rec['ref_eV'] = R
        out.append(rec)
    return out


# ------------------------------------------------------------ sources
def parse_lxcat(path):
    """Blocks of an LXCat text download: list of dict(keyword, lines, process, comments, table)."""
    lines = open(path, encoding='utf-8', errors='replace').read().splitlines()
    kws = ('ELASTIC', 'EFFECTIVE', 'EXCITATION', 'IONIZATION', 'ATTACHMENT')
    blocks = []
    i = 0
    while i < len(lines):
        if lines[i].strip() in kws:
            b = {'keyword': lines[i].strip(), 'header': [], 'table': [], 'line': i + 1}
            j = i + 1
            while j < len(lines) and not lines[j].startswith('-----'):
                b['header'].append(lines[j])
                j += 1
            j += 1
            while j < len(lines) and not lines[j].startswith('-----'):
                parts = lines[j].split()
                if len(parts) >= 2:
                    b['table'].append((float(parts[0]), float(parts[1])))
                j += 1
            b['process'] = next((h[len('PROCESS:'):].strip() for h in b['header'] if h.startswith('PROCESS:')), '')
            b['comment'] = ' '.join(h[len('COMMENT:'):].strip() for h in b['header'] if h.startswith('COMMENT:'))
            b['param'] = next((h[len('PARAM.:'):].strip() for h in b['header'] if h.startswith('PARAM.:')), '')
            blocks.append(b)
            i = j + 1
        else:
            i += 1
    return blocks


def nudge_duplicates(pts):
    """Strictly ascending energies: a repeated energy (a step) moves to the next binary64 up."""
    out = []
    for E, s in pts:
        if out and E <= out[-1][0]:
            E = math.nextafter(out[-1][0], math.inf)
        out.append((E, s))
    return out


class Tab:
    """Piecewise-linear tabulated function (linear in E and sigma)."""

    def __init__(self, pts, desc):
        self.pts = nudge_duplicates(pts)
        self.desc = desc

    def lo(self):
        return self.pts[0][0]

    def hi(self):
        return self.pts[-1][0]

    def __call__(self, E):
        p = self.pts
        if E < p[0][0] or E > p[-1][0]:
            raise ValueError('outside table')
        lo, hi = 0, len(p) - 1
        if E == p[-1][0]:
            return p[-1][1]
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if E < p[mid][0]:
                hi = mid
            else:
                lo = mid
        x0, y0 = p[lo]
        x1, y1 = p[lo + 1]
        return y0 + (y1 - y0) * ((E - x0) / (x1 - x0))

    def breakpoints(self):
        return [E for E, _ in self.pts]


class Model:
    def __init__(self, fn, lo, hi, per_decade, desc):
        self.fn, self._lo, self._hi, self.per_decade, self.desc = fn, lo, hi, per_decade, desc

    def lo(self):
        return self._lo

    def hi(self):
        return self._hi

    def __call__(self, E):
        return self.fn(E)

    def breakpoints(self):
        n = max(2, int(math.ceil(math.log10(self._hi / self._lo) * self.per_decade)) + 1)
        r = [self._lo * (self._hi / self._lo) ** (k / (n - 1)) for k in range(n)]
        r[0], r[-1] = self._lo, self._hi
        return r


def beb_orbital(E, B, U, N, a0, R, mc2, rel_, Q=1.0):
    if E <= B:
        return 0.0
    t = E / B
    if not rel_:
        u = U / B
        S = 4 * math.pi * a0 * a0 * N * (R / B) ** 2
        lt = math.log(t)
        return S / (t + u + 1) * (Q * lt / 2 * (1 - 1 / t ** 2) + (2 - Q) * (1 - 1 / t - lt / (t + 1)))
    if Q != 1.0:
        raise ValueError('RBEB is defined with Q = 1')
    tp, bp, up = E / mc2, B / mc2, U / mc2
    b2 = lambda x: x * (x + 2) / (1 + x) ** 2
    bt2, bb2, bu2 = b2(tp), b2(bp), b2(up)
    alpha2 = 2 * R / mc2
    P = 4 * math.pi * a0 * a0 * alpha2 ** 2 * N / ((bt2 + bu2 + bb2) * 2 * bp)
    L = math.log(tp * (tp + 2)) - bt2 - math.log(2 * bp)
    h = 1 + tp / 2
    phi = (1 + 2 * tp) / h ** 2
    cb = bp ** 2 / h ** 2
    return P * (L / 2 * (1 - 1 / t ** 2) + 1 - 1 / t - phi * math.log(t) / (t + 1) + cb * (t - 1) / 2)


def rudd_shape_py(P, B, N, Te, a0, R):
    v2 = Te / B
    v = math.sqrt(v2)
    L1 = P[2] * v ** P[3] / (1 + P[4] * v ** (P[3] + 4))
    H1 = P[0] * math.log1p(v2) / (v2 + P[1] / v2)
    L2 = P[7] * v ** P[8]
    H2 = P[5] / v2 + P[6] / v2 ** 2
    F1, F2 = L1 + H1, L2 * H2 / (L2 + H2)
    wc = 4 * v2 - 2 * v - R / (4 * B)
    a = P[9] / v
    S = 4 * math.pi * a0 * a0 * N * (R / B) ** 2
    return F1, F2, wc, a, S


def rudd_sdcs(W, B, N, P, Te, a0, R):
    """dsigma/dW (m^2/eV) of one orbital (Rudd et al. RMP 64, 441 (1992))."""
    F1, F2, wc, a, S = rudd_shape_py(P, B, N, Te, a0, R)
    w = W / B
    z = a * (w - wc)
    c = 0.0 if z > 700 else 1.0 / (1.0 + math.exp(z))
    return S / B * (F1 + F2 * w) / (1 + w) ** 3 * c


def rudd_inner_params():
    """Universal inner-shell Rudd parameters (Rudd et al. 1992 Table V), read from the extracted table and checked against
    the kernel constant PFC_RUDD_INNER_SHELL_PARAMS in chem_event.h (one set of numbers for the builder and the kernel)."""
    d = json.load(open(os.path.join(ROOT, 'data/extracted/ion-rudd1992-table5-sdcs-parameters-inner-shells.json')))
    vals = [float(x) for x in d['y']]
    hdr = open(os.path.join(ROOT, 'include/pf_chem_v3/chem_event.h')).read()
    line = next(l for l in hdr.splitlines() if l.startswith('#define PFC_RUDD_INNER_SHELL_PARAMS'))
    kv = [float(x) for x in line.split(None, 2)[2].split(',')]
    if kv != vals:
        raise ValueError('inner-shell Rudd parameters differ between the extracted table and chem_event.h')
    return vals


def rudd_params_for(o, orbs, P, inner):
    """Target set P, or the inner-shell set for an orbital with B > 2 B_min (Rudd et al. 1992, Sec. VI.D)."""
    bmin = min(x['B'] for x in orbs)
    return inner if o['B'] > 2.0 * bmin else P


def rudd_sigma(E, mP, mT, me, orbs, P, a0, R, q=0.0, nquad=400):
    """Total of the Rudd SDCS over the orbitals (weights w), W in [0, K_cm - B_eff]; Gauss-Legendre in
    t = ln(1 + W/B) on panels (exact enough for the builder: relative error < 1e-8)."""
    b2 = E * (E + 2 * mP) / (E + mP) ** 2
    Te = 0.5 * me * b2
    Kcm = (math.sqrt(mP * mP + mT * mT + 2 * mT * (E + mP)) - mP - mT)
    tot = 0.0
    xs, ws = _gl(16)
    inner = rudd_inner_params()
    P0 = P
    for o in orbs:
        P = rudd_params_for(o, orbs, P0, inner)
        B, N, wt = o['B'], o['N'], o.get('w', 1.0)
        Be = max(B, q)
        Wm = Kcm - Be
        if wt <= 0 or Wm <= 0:
            continue
        tmax = math.log1p(Wm / B)
        npan = 64
        s_ = 0.0
        for k_ in range(npan):
            t0, t1 = tmax * k_ / npan, tmax * (k_ + 1) / npan
            for x, w_ in zip(xs, ws):
                t = 0.5 * (t0 + t1) + 0.5 * (t1 - t0) * x
                W = B * math.expm1(t)
                s_ += 0.5 * (t1 - t0) * w_ * rudd_sdcs(W, B, N, P, Te, a0, R) * B * math.exp(t)
        tot += wt * s_
    return tot


def rudd_bin_fracs(E, mP, mT, me, orbs, P, a0, R, q, edges):
    """Fractions of the summed Rudd SDCS (orbital weights w) in W bins [edges[k], edges[k+1]) (W in eV),
    each orbital truncated at K_cm - B_eff; independent of the kernel (Python quadrature)."""
    b2 = E * (E + 2 * mP) / (E + mP) ** 2
    Te = 0.5 * me * b2
    Kcm = (math.sqrt(mP * mP + mT * mT + 2 * mT * (E + mP)) - mP - mT)
    xs, ws = _gl(24)
    acc = [0.0] * (len(edges) - 1)
    inner = rudd_inner_params()
    P0 = P
    for o in orbs:
        P = rudd_params_for(o, orbs, P0, inner)
        B, N, wt = o['B'], o['N'], o.get('w', 1.0)
        Wm = Kcm - max(B, q)
        if wt <= 0 or Wm <= 0:
            continue
        for k in range(len(edges) - 1):
            a, b = edges[k], min(edges[k + 1], Wm)
            if b <= a:
                continue
            t0, t1 = math.log1p(a / B), math.log1p(b / B)
            npan = 16
            s_ = 0.0
            for j in range(npan):
                u0, u1 = t0 + (t1 - t0) * j / npan, t0 + (t1 - t0) * (j + 1) / npan
                for x, w_ in zip(xs, ws):
                    t = 0.5 * (u0 + u1) + 0.5 * (u1 - u0) * x
                    W = B * math.expm1(t)
                    s_ += 0.5 * (u1 - u0) * w_ * rudd_sdcs(W, B, N, P, Te, a0, R) * B * math.exp(t)
            acc[k] += wt * s_
    tot = sum(acc)
    return [x / tot for x in acc], tot


_GL = {}


def _gl(n):
    if n not in _GL:
        xs, ws = [], []
        for i in range(1, n + 1):
            x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
            for _ in range(100):
                p0, p1 = 1.0, x
                for k in range(2, n + 1):
                    p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
                dp = n * (x * p1 - p0) / (x * x - 1)
                dx = p1 / dp
                x -= dx
                if abs(dx) < 1e-16:
                    break
            xs.append(x)
            ws.append(2 / ((1 - x * x) * dp * dp))
        _GL[n] = (xs, ws)
    return _GL[n]


def make_source(src, ctx):
    """Return a Tab or Model for a source spec; records provenance in ctx['files']."""
    t = src['type']
    if t == 'lxcat':
        path = rel(src['file'])
        ctx['files'][src['file']] = sha_file(path)
        blocks = parse_lxcat(path)
        cand = [b for b in blocks if b['keyword'] == src.get('keyword', b['keyword']) and src['process'] == b['process']
                and src.get('comment', '') in b['comment']]
        k = src.get('occurrence', 0)
        if len(cand) <= k:
            raise ValueError('LXCat block not found: %s' % src)
        b = cand[k]
        scale = src.get('scale', 1.0)
        return Tab([(E, s * scale) for E, s in b['table']], 'LXCat %s line %d: %s | %s' % (src['file'], b['line'], b['process'], b['comment'][:120]))
    if t == 'raikou':
        path = rel(src['file'])
        ctx['files'][src['file']] = sha_file(path)
        d = json.load(open(path))
        ch = next(c for c in d['channels'] if c['id'] == src['channel'])
        sc = ch['ordinate_scale_to_SI'] * src.get('scale', 1.0)
        return Tab([(E, v * sc) for E, v in zip(ch['energy_eV'], ch['values'])], 'dataset %s channel %s (%s %s)' % (d['id'], ch['id'], ch['source'], ch['source_locator']))
    if t == 'extracted':
        path = rel(os.path.join('data/extracted', src['table_id'] + '.json'))
        ctx['files'][os.path.relpath(path, ROOT)] = sha_file(path)
        d = json.load(open(path))
        sc = src.get('scale', 1.0)
        xs = d['x']
        if src.get('x_scale'):
            xs = [x * src['x_scale'] for x in xs]
        return Tab([(x, y * sc) for x, y in zip(xs, d['y'])], 'extracted %s (%s)' % (d['table_id'], d.get('locator', '')[:160]))
    if t == 'inline':
        return Tab(list(zip(src['E'], src['sigma'])), 'inline: ' + src['locator'])
    if t == 'zero':
        return Model(lambda E: 0.0, src['lo'], src['hi'], 1, 'zero')
    if t == 'beb':
        k = src
        orbs = k['orbitals']
        rel_from = k.get('relativistic_from_eV', 1e300)

        def f(E, orbs=orbs):
            s = 0.0
            for o in orbs:
                s += o.get('w', 1.0) * beb_orbital(E, o['B'], o['U'], o['N'], k['bohr_m'], k['rydberg_eV'], k['mc2_eV'], E > rel_from, o.get('Q', 1.0))
            return s * k.get('scale', 1.0)
        return Model(f, k['lo'], k['hi'], k.get('per_decade', 64), 'BEB/RBEB orbital model: ' + k.get('locator', ''))
    if t == 'lennon1988':
        # Lennon et al. JPCRD 17, 1285 (1988) Eq. (1) with exponent i:
        # sigma[cm^2] = 1e-13/(I E) [A ln(E/I) + sum_i B_i (1 - I/E)^i]
        k = src
        I, A, Bs = k['I_eV'], k['A'], k['B']

        def f(E):
            if E <= I:
                return 0.0
            x = 1 - I / E
            return 1e-13 / (I * E) * (A * math.log(E / I) + sum(b * x ** (i + 1) for i, b in enumerate(Bs))) * 1e-4
        return Model(f, k['lo'], k['hi'], k.get('per_decade', 64), 'Lennon 1988 Eq. (1) fit: ' + k.get('locator', ''))
    if t == 'bethe':
        # sigma = (4 pi a0^2 alpha^2 / beta^2) [M2 (ln(beta^2/(1-beta^2)) - beta^2) + C]
        k = src
        a0, alpha, mc2 = k['bohr_m'], k['alpha'], k['mc2_eV']
        pref = 4 * math.pi * a0 * a0 * alpha * alpha

        def f(E):
            g = 1 + E / mc2
            b2 = 1 - 1 / (g * g)
            return pref / b2 * (k['M2'] * (math.log(b2 / (1 - b2)) - b2) + k['C']) * k.get('scale', 1.0)
        return Model(f, k['lo'], k['hi'], k.get('per_decade', 64), 'Bethe fit: ' + k.get('locator', ''))
    if t == 'product':
        a = make_source(src['a'], ctx)
        f = src['factor']
        if isinstance(a, Tab):
            return Tab([(E, s * f) for E, s in a.pts], '%r x (%s)' % (f, a.desc))
        return Model(lambda E: a(E) * f, a.lo(), a.hi(), a.per_decade, '%r x (%s)' % (f, a.desc))
    if t == 'ratio_tab':
        # numerator(E)/denominator(E) * base(E) on the base breakpoints (model join)
        raise ValueError('ratio_tab not supported')
    if t == 'rate_k':
        k = src

        def f(E):
            m = ctx['mass_of'](k['projectile'])
            M = ctx['mass_of'](k['target'])
            v = speed_tf(E, m) if not k.get('cm_energy') else None
            return k['k_m3_s'] / v
        return Model(f, k['lo'], k['hi'], k.get('per_decade', 64), 'constant rate coefficient k/v: ' + k.get('locator', ''))
    if t == 'scaled_shape':
        # anchor(E0) * shape(E) / shape(E0): the energy dependence of `shape`, normalized to `anchor` at E0
        sh = make_source(src['shape'], ctx)
        an = make_source(src['anchor'], ctx)
        E0 = src['E0']
        fac = an(E0) / sh(E0)
        lo, hi = src.get('lo', sh.lo()), src.get('hi', sh.hi())
        desc = 'shape of [%s] scaled to [%s] at %g eV (factor %.10g)' % (sh.desc, an.desc, E0, fac)
        if isinstance(sh, Tab):
            pts = [(E, v * fac) for E, v in sh.pts if lo <= E <= hi]
            if pts[0][0] > lo:
                pts.insert(0, (lo, sh(lo) * fac))
            if pts[-1][0] < hi:
                pts.append((hi, sh(hi) * fac))
            t_ = Tab(pts, desc)
            t_.factor = fac
            return t_
        m_ = Model(lambda E: sh(E) * fac, lo, hi, sh.per_decade, desc)
        m_.factor = fac
        return m_
    if t == 'extend_constant':
        base = make_source(src['base'], ctx)
        pts = [(E, base(E)) for E in base.breakpoints()]
        pts.append((src['hi'], pts[-1][1]))
        return Tab(pts, '[%s] held constant from %g to %g eV (upper bound of a decreasing function)' % (base.desc, pts[-2][0], src['hi']))
    if t == 'restrict':
        base = make_source(src['base'], ctx)
        lo, hi = src['lo'], src['hi']
        pts = [(E, base(E)) for E in base.breakpoints() if lo <= E <= hi]
        if not pts or pts[0][0] > lo:
            pts.insert(0, (lo, base(lo)))
        if pts[-1][0] < hi:
            pts.append((hi, base(hi)))
        return Tab(pts, '[%s] restricted to [%g, %g] eV' % (base.desc, lo, hi))
    if t == 'joined_table':
        # power law (or any source) up to E0, then a measured table: the table is preceded by the anchor value
        # at E0 so the junction is continuous; table points at or below E0 are not used
        an = make_source(src['anchor'], ctx)
        tb = make_source(src['table'], ctx)
        E0 = src['E0']
        pts = [(E0, an(E0))] + [(E, v) for E, v in tb.pts if E > E0]
        t_ = Tab(pts, 'join at %g eV of [%s] (value %.6g m2) and [%s]' % (E0, an.desc, an(E0), tb.desc))
        t_.join_value = an(E0)
        return t_
    if t == 'rate_anchor':
        # constant rate below E1: sigma(E) = sigma_ref(E1) v(E1) / v(E), v the target-frame projectile speed
        ref = make_source(src['ref'], ctx)
        E1 = src['E1']
        m = ctx['mass_of'](src['projectile'])
        k = ref(E1) * speed_tf(E1, m)
        mod = Model(lambda E: k / speed_tf(E, m), src['lo'], E1, src.get('per_decade', 64),
                    'constant rate %.6g m3/s anchored at %g eV of [%s]' % (k, E1, ref.desc))
        mod.k = k
        return mod
    if t == 'loglog_bridge':
        a_ = make_source(src['a'], ctx)
        b_ = make_source(src['b'], ctx)
        Ea, Eb = src['Ea'], src['Eb']
        sa, sb = a_(Ea), b_(Eb)
        pw = math.log(sb / sa) / math.log(Eb / Ea)
        return Model(lambda E: sa * (E / Ea) ** pw, Ea, Eb, src.get('per_decade', 64),
                     'log-log interpolation from %.6g m2 at %g eV [%s] to %.6g m2 at %g eV [%s] (exponent %.4f)' % (sa, Ea, a_.desc, sb, Eb, b_.desc, pw))
    if t == 'none':
        return Model(lambda E: 0.0, src['lo'], src['hi'], 1, 'no data: ' + src.get('locator', ''))
    if t == 'langevin':
        # sigma_m = f * k_L / v, k_L = 2 pi |q| e sqrt(alpha_vol / (4 pi eps0 mu)) (classical r^-4 polarization capture);
        # f = 1.105 is the momentum-transfer / capture ratio of the polarization potential
        k = src
        e, eps0, u = 1.602176634e-19, 8.8541878188e-12, 1.66053906892e-27
        mp = ctx['mass_of'](k['projectile']) / 931494103.72 * u
        mt = ctx['mass_of'](k['target']) / 931494103.72 * u
        mu = mp * mt / (mp + mt)
        kL = 2 * math.pi * abs(k['charge']) * e * math.sqrt(k['alpha_A3'] * 1e-30 / (4 * math.pi * eps0 * mu))
        m_eV = ctx['mass_of'](k['projectile'])
        f = k.get('factor', 1.105)
        mod = Model(lambda E: f * kL / speed_tf(E, m_eV), k['lo'], k['hi'], k.get('per_decade', 64),
                    'Langevin polarization model, alpha = %g A^3 (%s), k_L = %.6g m3/s, sigma_m = %g k_L / v' % (k['alpha_A3'], k.get('locator', ''), kL, f))
        mod.kL = kL
        return mod
    if t == 'dr_powerlaw':
        # sigma(E) = A E^-p_i on consecutive segments (continuous), A fixed so that the Maxwellian rate
        # coefficient at T equals the measured alpha (electron on an ion at rest)
        k = src
        segs = k['segments']   # [[E_lo, E_hi, p], ...] contiguous

        cs = [1.0]
        for i_ in range(1, len(segs)):
            Eb = segs[i_ - 1][1]
            cs.append(cs[-1] * Eb ** -segs[i_ - 1][2] / Eb ** -segs[i_][2])

        def shape(E):
            for (lo_, hi_, p_), c_ in zip(segs, cs):
                if E <= hi_:
                    return c_ * E ** -p_
            return cs[-1] * E ** -segs[-1][2]
        kT = 8.617333262e-5 * k['T_K']
        me = 510998.95069
        # Maxwellian <sigma v> = int shape(E) v(E) f(E) dE, f = 2 sqrt(E/pi) kT^-1.5 exp(-E/kT); log-spaced quadrature
        n = 20000
        lo_, hi_ = 1e-9, 60 * kT
        tot = 0.0
        prev = None
        for i in range(n + 1):
            E = lo_ * (hi_ / lo_) ** (i / n)
            g = shape(E) * speed_tf(E, me) * 2 * math.sqrt(E / math.pi) * kT ** -1.5 * math.exp(-E / kT) * E   # dE = E dlnE
            if prev is not None:
                tot += 0.5 * (g + prev) * math.log(hi_ / lo_) / n
            prev = g
        A = k['alpha_m3_s'] / tot
        mod = Model(lambda E: A * shape(E) * k.get('fraction', 1.0), k['lo'], k['hi'], k.get('per_decade', 64),
                    'power law %s normalized to alpha(%g K) = %g m3/s (%s), branch fraction %g' % (segs, k['T_K'], k['alpha_m3_s'], k.get('locator', ''), k.get('fraction', 1.0)))
        return mod
    if t == 'piecewise_bound':
        parts = [make_source(x, ctx) for x in src['parts']]
        bps = sorted(set(E for p_ in parts for E in p_.breakpoints()))
        pts = []
        for E in bps:
            vals = [p_(E) for p_ in parts if p_.lo() <= E <= p_.hi()]
            pts.append((E, max(vals) if vals else 0.0))
        return Tab(pts, 'max of [' + '; '.join(p_.desc for p_ in parts) + ']')
    if t == 'iam_elastic':
        path = rel(src['file'])
        ctx['files'][src['file']] = sha_file(path)
        a = json.load(open(path))['conditional_IAM_elastic_angular_data']
        mu = a['mu']
        a0sq = a['a0_squared_m2']
        pts = []
        for E, row in zip(a['energy_eV'], a['dcs_a0_squared_per_sr']):
            tot = 0.0
            for i in range(len(mu) - 1):
                m0, m1, d0, d1 = mu[i], mu[i + 1], row[i], row[i + 1]
                h = m1 - m0
                if src['quantity'] == 'integral':
                    tot += 0.5 * h * (d0 + d1)
                else:   # (1 - mu) * dcs, dcs linear in mu: Simpson is exact for the quadratic integrand
                    mm = 0.5 * (m0 + m1)
                    dm = 0.5 * (d0 + d1)
                    tot += h / 6.0 * ((1 - m0) * d0 + 4 * (1 - mm) * dm + (1 - m1) * d1)
            pts.append((float(E), 2 * math.pi * tot * a0sq))
        return Tab(pts, 'IAM (NIST SRD 64 O + 2H, raikou %s) %s cross section integrated from the piecewise-linear DCS' % (a['model'], src['quantity']))
    if t == 'zbl_momentum_transfer':
        # v3.2: equivalent momentum-transfer cross section of projectile P on target T from the ZBL universal
        # nuclear stopping (Ziegler, Biersack & Littmark 1985), Bragg-additive over the atoms of P (each at the
        # common speed, energy E m_a / m_P) and of T: S_n,tot(E) = sum_ab S_n(a -> b, E m_a / m_P) and
        # sigma_m = 2 S_n,tot / (gamma E), gamma = 4 m_P m_T / (m_P + m_T)^2, so that isotropic CM scattering with
        # sigma_m reproduces the summed nuclear stopping (mean energy transfer gamma E / 2 per event)
        k = src
        U = 931494103.72
        Zof = {'H': 1, 'O': 8}
        mP, mT = ctx['mass_of'](k['projectile']), ctx['mass_of'](k['target'])
        ma = {a: ctx['mass_of'](a) / U for a in set(k['projectile_atoms'] + k['target_atoms'])}
        gam = 4 * mP * mT / (mP + mT) ** 2

        def sn_red(eps):
            if eps <= 30.0:
                return math.log(1 + 1.1383 * eps) / (2 * (eps + 0.01321 * eps ** 0.21226 + 0.19593 * eps ** 0.5))
            return math.log(eps) / (2 * eps)

        def f(E):
            tot = 0.0
            for a in k['projectile_atoms']:
                Ea_keV = E * (ma[a] * U / mP) / 1e3
                for b in k['target_atoms']:
                    z1, z2, m1, m2 = Zof[a], Zof[b], ma[a], ma[b]
                    au = z1 ** 0.23 + z2 ** 0.23
                    eps = 32.53 * m2 * Ea_keV / (z1 * z2 * (m1 + m2) * au)
                    tot += 8.462 * z1 * z2 * m1 * sn_red(eps) / ((m1 + m2) * au) * 1e-19   # eV m^2 (1 eV cm^2/1e15 = 1e-19 eV m^2)
            return 2 * tot / (gam * E) * k.get('factor', 1.0)
        return Model(f, k['lo'], k['hi'], k.get('per_decade', 32),
                     'ZBL universal nuclear stopping, Bragg-additive over the atoms of %s %s on %s %s, sigma_m = 2 S_n / (gamma E) (%s)'
                     % (k['projectile'], k['projectile_atoms'], k['target'], k['target_atoms'], k.get('locator', '')))
    if t == 'wedlund_fit':
        # v3.2: Simon Wedlund et al. 2019 recommended polynomial fit log10(sigma[m2]) = sum_j p_j (log10 v[m/s])^j of the
        # extracted coefficient table, v the projectile speed in the target frame
        path = rel(os.path.join('data/extracted', src['table_id'] + '.json'))
        ctx['files'][os.path.relpath(path, ROOT)] = sha_file(path)
        d = json.load(open(path))
        pj = [float(x) for x in d['y']]
        m = ctx['mass_of'](src['projectile'])

        def f(E):
            lv = math.log10(speed_tf(E, m))
            return 10.0 ** sum(c * lv ** j for j, c in enumerate(pj)) * src.get('factor', 1.0)
        return Model(f, src['lo'], src['hi'], src.get('per_decade', 64),
                     'Simon Wedlund 2019 polynomial fit %s (%s; validity %s km/s, confidence %s)' % (d['table_id'], d['locator'][:60], d['validity_speed_km_s'], d['confidence_level']))
    if t == 'nuclear_stopping_mt':
        # v3.2: tabulated nuclear stopping S_n (eV cm2 per 1e15 molecules, energies in MeV) -> equivalent
        # momentum-transfer cross section sigma_m = 2 S_n / (gamma E) (isotropic CM scattering reproduces S_n)
        path = rel(os.path.join('data/extracted', src['table_id'] + '.json'))
        ctx['files'][os.path.relpath(path, ROOT)] = sha_file(path)
        d = json.load(open(path))
        mP, mT = ctx['mass_of'](src['projectile']), ctx['mass_of'](src['target'])
        gam = 4 * mP * mT / (mP + mT) ** 2
        pts = []
        for x, y in zip(d['x'], d['y_eV_cm2_per_1e15_molecules']):
            E = x * 1e6
            if src.get('lo', 0.0) <= E <= src.get('hi', 1e300):
                pts.append((E, 2 * y * 1e-19 / (gam * E)))
        return Tab(pts, 'nuclear stopping %s (%s) as sigma_m = 2 S_n / (gamma E)' % (d['table_id'], d['locator'][:80]))
    if t == 'close_distant_bound':
        # v3.2 upper bound for an inelastic process of projectile P (net charge q) at speed v:
        #   sigma <= pi R^2 + q^2 sigma_ref(v),
        # close collisions (impact parameters below R = R_P + R_T, probability <= 1) plus distant collisions, in
        # which P acts as its net charge and the first-order (Bethe) cross section scales as q^2 relative to a
        # reference projectile of unit charge at the same speed (sigma_ref: measured proton data, or the free-electron
        # total cross section at the electron energy of that speed for detachment, see the spec).
        k = src
        mP = ctx['mass_of'](k['projectile'])
        parts = [make_source(x, ctx) for x in k['ref_parts']]
        geo = math.pi * (k['R_sum_A'] * 1e-10) ** 2
        mode = k.get('ref_mode', 'equal_velocity')
        if mode == 'equal_velocity':
            r = ctx['mass_of'](k['ref_projectile']) / mP
            conv = lambda E: E * r
        else:   # electron with the projectile speed: E_e = m_e (gamma - 1) for the projectile's gamma
            me = ctx['mass_of']('e')
            conv = lambda E: me * (E / mP)

        def ref_at(Eref):
            for p_ in parts:
                if p_.lo() <= Eref <= p_.hi():
                    return p_(Eref)
            if Eref < parts[0].lo():
                return parts[0](parts[0].lo()) if k.get('hold_below') else 0.0
            return parts[-1](parts[-1].hi())   # held constant above (decreasing reference)
        q2 = float(k['q']) ** 2 * k.get('ref_factor', 1.0)
        emin = k.get('ref_min_eV', 0.0)   # reference term only where the reference energy reaches emin

        def f(E):
            Er = conv(E)
            return geo + (q2 * ref_at(Er) if Er >= emin else 0.0)
        lo, hi = k['lo'], k['hi']
        nd = k.get('per_decade', 32)
        nn = max(2, int(math.ceil(math.log10(hi / lo) * nd)) + 1)
        pts = set(lo * (hi / lo) ** (i_ / (nn - 1)) for i_ in range(nn)) | {lo, hi}
        for p_ in parts:   # reference breakpoints mapped to projectile energy
            for Er in p_.breakpoints():
                E = Er / conv(1.0)
                if lo <= E <= hi:
                    pts.add(E)
        if emin > 0:
            Eon = emin / conv(1.0)   # conv is linear in E
            if lo < Eon < hi:
                pts.add(Eon)
                pts.add(math.nextafter(Eon, -math.inf))
        pts = sorted(E for E in pts if lo <= E <= hi)
        return Tab([(E, f(E)) for E in pts],
                   'bound pi (%.3g A)^2 + %g x [%s] (%s%s)' % (k['R_sum_A'], q2, '; '.join(p_.desc for p_ in parts), mode,
                                                               ', reference term from %g eV reference energy' % emin if emin else ''))
    if t == 'equal_velocity':
        # v3.2: sigma(E) = factor * ref(E m_ref / m_P): a reference projectile's cross section at the same speed
        k = src
        ref_ = make_source(k['ref'], ctx)
        r = ctx['mass_of'](k['ref_projectile']) / ctx['mass_of'](k['projectile'])
        lo, hi = max(k.get('lo', 0.0), ref_.lo() / r), min(k.get('hi', 1e300), ref_.hi() / r)
        fac = k.get('factor', 1.0)
        if isinstance(ref_, Tab):
            pts = [(E / r, v * fac) for E, v in ref_.pts if lo <= E / r <= hi]
            return Tab(pts, '%g x [%s] at equal speed (E x %.6g)' % (fac, ref_.desc, r))
        return Model(lambda E: fac * ref_(E * r), lo, hi, ref_.per_decade, '%g x [%s] at equal speed (E x %.6g)' % (fac, ref_.desc, r))
    if t == 'constant':
        return Model(lambda E, v=src['value']: v, src['lo'], src['hi'], 1, 'constant %g m2: %s' % (src['value'], src.get('locator', '')))
    if t == 'powerlaw_from':
        # sigma(E0) of a reference source continued as (E/E0)^p
        k = src
        ref_ = make_source(k['ref'], ctx)
        E0 = k['E0']
        s0 = ref_(E0)
        return Model(lambda E: s0 * (E / E0) ** k['p'], k.get('lo', E0), k['hi'], k.get('per_decade', 16),
                     '%.6g m2 at %g eV of [%s] continued as (E/E0)^%g' % (s0, E0, ref_.desc, k['p']))
    if t == 'product_of':
        # pointwise product a(E) b(E) sampled on both breakpoint sets plus `per_decade` points (linear between)
        a_ = make_source(src['a'], ctx)
        b_ = make_source(src['b'], ctx)
        lo, hi = max(a_.lo(), b_.lo(), src.get('lo', 0.0)), min(a_.hi(), b_.hi(), src.get('hi', 1e300))
        pts = set([lo, hi])
        for E in a_.breakpoints() + b_.breakpoints():
            if lo <= E <= hi:
                pts.add(E)
        nd = src.get('per_decade', 32)
        nn = max(2, int(math.ceil(math.log10(hi / lo) * nd)) + 1)
        for i_ in range(nn):
            pts.add(lo * (hi / lo) ** (i_ / (nn - 1)))
        pts = sorted(E for E in pts if lo <= E <= hi)
        return Tab([(E, a_(E) * b_(E)) for E in pts], '[%s] x [%s]' % (a_.desc, b_.desc))
    if t == 'min_of':
        parts = [make_source(x, ctx) for x in src['parts']]
        lo, hi = max(p_.lo() for p_ in parts), min(p_.hi() for p_ in parts)
        bps = sorted(set(E for p_ in parts for E in p_.breakpoints() if lo <= E <= hi) | {lo, hi})
        return Tab([(E, min(p_(E) for p_ in parts)) for E in bps], 'min of [' + '; '.join(p_.desc for p_ in parts) + ']')
    if t == 'rudd_total':
        # total ionization cross section of the Rudd SDCS (sum over orbitals, integrated over W to K_cm - B_eff)
        k = src
        mP, mT = ctx['mass_of'](k['projectile']), ctx['mass_of'](k['target'])
        me = ctx['mass_of']('e')

        def f(E):
            return rudd_sigma(E, mP, mT, me, k['orbitals'], k['rudd'], k['bohr_m'], k['rydberg_eV'], k.get('q_eV', 0.0)) * k.get('factor', 1.0)
        return Model(f, k['lo'], k['hi'], k.get('per_decade', 32), 'Rudd SDCS integrated (%s)' % k.get('locator', ''))
    if t == 'sum_of':
        # v3.2: pointwise sum of sources on the union of their breakpoints, each held at its end value outside its range
        # when 'hold' is set (otherwise zero outside)
        parts = [make_source(x, ctx) for x in src['parts']]
        lo, hi = src.get('lo', min(p_.lo() for p_ in parts)), src.get('hi', max(p_.hi() for p_ in parts))

        def val(p_, E):
            if p_.lo() <= E <= p_.hi():
                return p_(E)
            if not src.get('hold'):
                return 0.0
            return p_(p_.lo()) if E < p_.lo() else p_(p_.hi())
        bps = sorted(set(E for p_ in parts for E in p_.breakpoints() if lo <= E <= hi) | {lo, hi})
        return Tab([(E, sum(val(p_, E) for p_ in parts)) for E in bps], 'sum of [' + '; '.join(p_.desc for p_ in parts) + ']')
    if t == 'threshold_powerlaw':
        # v3.2 rise from zero at the threshold E_th to the first point (E0, s0) of a reference table, as a power law in
        # (E - E_th) whose log-log slope equals that of the table's first interval (log-log in E - E_th, continuous value
        # and slope at E0): sigma = s0 ((E - E_th)/(E0 - E_th))^p
        k = src
        ref_ = make_source(k['ref'], ctx)
        Eth, E0, E1 = k['E_th'], k['E0'], k['E1']
        s0, s1 = ref_(E0), ref_(E1)
        p_ = math.log(s1 / s0) / math.log((E1 - Eth) / (E0 - Eth))
        mod = Model(lambda E: s0 * ((E - Eth) / (E0 - Eth)) ** p_ if E > Eth else 0.0, k['lo'], E0, k.get('per_decade', 64),
                    'rise from 0 at %.6g eV to %.6g m2 at %.6g eV of [%s] as ((E - E_th)/(E0 - E_th))^%.4f' % (Eth, s0, E0, ref_.desc, p_))
        bps = mod.breakpoints()
        # refine just above the threshold, where the power law is steep
        mod.breakpoints = lambda: sorted(set(bps) | {Eth + (E0 - Eth) * 10 ** (-k_ / 4.0) for k_ in range(1, 17) if k['lo'] <= Eth + (E0 - Eth) * 10 ** (-k_ / 4.0) <= E0})
        return mod
    if t == 'adiabatic_scaled':
        # v3.2 energy-defect scaling of a near-resonant reference capture cross section at the same speed (Massey adiabatic
        # criterion; H. S. W. Massey, Rep. Prog. Phys. 12, 248 (1949); J. B. Hasted, Proc. R. Soc. A 205, 421 (1951)):
        # sigma(E) = ref(E) exp(-a (|dE| - |dE_ref|) / (h v)), a the adiabatic length, v the projectile speed
        k = src
        ref_ = make_source(k['ref'], ctx)
        m = ctx['mass_of'](k['projectile'])
        a = k['a_A'] * 1e-10
        dd = (abs(k['dE_eV']) - abs(k['dE_ref_eV'])) * 1.602176634e-19
        h = 6.62607015e-34
        lo, hi = max(k['lo'], ref_.lo()), min(k['hi'], ref_.hi())
        return Model(lambda E: ref_(E) * math.exp(-a * dd / (h * speed_tf(E, m))), lo, hi, k.get('per_decade', 64),
                     '[%s] x exp(-a (|dE| - |dE_ref|) / (h v)), a = %g A, dE = %.4g eV, dE_ref = %.4g eV' % (ref_.desc, k['a_A'], k['dE_eV'], k['dE_ref_eV']))
    if t == 'free_collision':
        # v3.2 free-collision estimate of projectile electron loss / breakup: the most weakly bound projectile electron
        # (binding or breakup threshold I) scattered as a free electron of the projectile speed by the target; for isotropic
        # scattering the share of collisions transferring W = 2 E_e (1 - cos theta) >= I in the projectile frame is
        # 1 - I / (4 E_e):  sigma = n x ref(E_e) x max(0, 1 - I / (4 E_e)),  E_e = (m_e / m_P) E
        k = src
        mP, me = ctx['mass_of'](k['projectile']), ctx['mass_of']('e')
        ref_ = make_source(k['ref'], ctx)
        I, n_ = k['I_eV'], k.get('n', 1.0)

        def refv(Ee):
            if Ee < ref_.lo():
                return ref_(ref_.lo())
            if Ee > ref_.hi():
                return ref_(ref_.hi())
            return ref_(Ee)

        def f(E):
            Ee = me * E / mP
            fac = 1.0 - I / (4.0 * Ee)
            return n_ * refv(Ee) * fac if fac > 0 else 0.0
        lo, hi = k['lo'], k['hi']
        nd = k.get('per_decade', 32)
        nn = max(2, int(math.ceil(math.log10(hi / lo) * nd)) + 1)
        pts = set(lo * (hi / lo) ** (i_ / (nn - 1)) for i_ in range(nn)) | {lo, hi}
        for Er in ref_.breakpoints():
            E = Er * mP / me
            if lo <= E <= hi:
                pts.add(E)
        Eon = I / 4.0 * mP / me
        if lo < Eon < hi:
            pts.add(Eon)
            for k_ in range(1, 9):   # the onset factor rises steeply: refine the first decade above the onset
                pts.add(Eon * 10 ** (k_ / 8.0)) if Eon * 10 ** (k_ / 8.0) < hi else None
        pts = sorted(E for E in pts if lo <= E <= hi)
        t_ = Tab([(E, f(E)) for E in pts], 'free-collision estimate %g x [%s] x max(0, 1 - %g eV / (4 E_e)), E_e = (m_e/m_P) E' % (n_, ref_.desc, I))
        return t_
    if t == 'rudd_analytic':
        # v3.2: Rudd-type analytic total cross section of the extracted fit table (Gobet et al. 2004 Eq. (1)):
        # sigma = 4 pi a0^2 / (x / (A ln(1 + x) + B) + 1 / (C x^D)), x = (E_p / 1836) / 13.6, for the equal-speed proton
        # energy E_p = E m_p / m_P of projectile P
        path = rel(os.path.join('data/extracted', src['table_id'] + '.json'))
        ctx['files'][os.path.relpath(path, ROOT)] = sha_file(path)
        d = json.load(open(path))
        A, B, C, D = d['y'][d['x'].index(src['row'])]
        a0 = src['bohr_m']
        r = ctx['mass_of']('H+') / ctx['mass_of'](src['projectile'])

        def f(E):
            x = (E * r / 1836.0) / 13.6
            den = x / (A * math.log1p(x) + B) + 1.0 / (C * x ** D)
            return 4 * math.pi * a0 * a0 / den
        return Model(f, src['lo'], src['hi'], src.get('per_decade', 32),
                     '%s row %s (A %g, B %g, C %g, D %g), Rudd-type form' % (d['table_id'], src['row'], A, B, C, D))
    if t == 'branch_fraction':
        # v3.2: fraction parts[index] / sum(parts) evaluated on [lo, hi] (the validity range of the parts) and held at its
        # end values outside, on [range_lo, range_hi]
        parts = [make_source(x, ctx) for x in src['parts']]
        lo, hi = src['lo'], src['hi']
        k_ = src['index']

        def frac(E):
            Ec = min(max(E, lo), hi)
            vals = [p_(Ec) for p_ in parts]
            return vals[k_] / sum(vals)
        mod = Model(frac, src['range_lo'], src['range_hi'], src.get('per_decade', 32),
                    'fraction %d of [%s] on %g-%g eV, held outside' % (k_, '; '.join(p_.desc for p_ in parts), lo, hi))
        bps = mod.breakpoints()
        mod.breakpoints = lambda: sorted(set(bps) | {E for E in (lo, hi) if src['range_lo'] <= E <= src['range_hi']})
        return mod
    if t == 'piecewise':
        raise ValueError('piecewise handled by caller')
    raise ValueError('unknown source type %s' % t)


def speed_tf(E, m):
    """Speed (m/s) of kinetic energy E, rest energy m (eV)."""
    return C_M_S * math.sqrt(E * (E + 2 * m)) / (E + m)


# ------------------------------------------------------------ assembly
def kin_threshold(m, M, Q):
    return Q * (2 * (m + M) + Q) / (2 * M) if Q > 0 else 0.0


def build(spec_path, out_prefix):
    spec = json.load(open(spec_path))
    files = {os.path.relpath(spec_path, ROOT): sha_file(spec_path)}
    species = species_table(spec, files)
    sidx = {s['name']: i for i, s in enumerate(species)}
    mass = {s['name']: s['mass_eV'] for s in species}
    ref = {s['name']: s['ref_eV'] for s in species}
    ctx = {'files': files, 'mass_of': lambda n: mass[n]}
    bottom, top = spec['energy_bottom_eV'], spec['energy_top_eV']
    k_const = spec['constants']

    groups_out, procs_out, grid_all, sigma_all, bound_all, rmq_all, xi_all = [], [], [], [], [], [], []
    bin_region_all, bin_row_all, bin_off = [], [], []
    regions = []
    region_info = []
    manifest_procs, checks = [], []
    rudd_specs = {}
    seams = []
    top_global, bottom_global = top, bottom
    for gspec in spec['groups']:
        P, T = gspec['projectile'], gspec['target']
        bottom = gspec.get('energy_bottom_eV', bottom_global)   # per-group grid bottom (v3.1: 1e-9 eV for electron-ion DR)
        mp, mt = mass[P], mass[T]
        gname = '%s+%s' % (P, T)
        proc_items = []   # (pspec, segments[(lo, hi, src_obj, class, tag)], cover list, bound src)
        pts = set([bottom, top])
        for ps in gspec['processes']:
            if '|' in ps['id'] or (',' in ps['id'] and ps.get('placeholder_below_eV')):
                # region names list process ids separated by ',' between '|' fields; placeholder ids are parsed by the tests
                raise ValueError('process id contains a region-name separator: %s' % ps['id'])
            segs = []
            sig = ps['sigma']
            seglist = sig['segments'] if sig['type'] == 'piecewise' else [{'source': sig, 'class': ps.get('class', 'DATA')}]
            for sg in seglist:
                so = make_source(sg['source'], ctx)
                lo = sg.get('lo', so.lo())
                hi = sg.get('hi', so.hi())
                segs.append({'lo': lo, 'hi': hi, 'src': so, 'class': sg.get('class', 'DATA'), 'tag': sg.get('tag', ps['tag']),
                             'desc': so.desc})
                for E in so.breakpoints():
                    if lo <= E <= hi:
                        pts.add(E)
                pts.add(lo)
                pts.add(hi)
            # junction points: the right segment starts at the next binary64 above a shared edge
            for a, b in zip(segs, segs[1:]):
                if b['lo'] <= a['hi']:
                    b['lo_eff'] = math.nextafter(a['hi'], math.inf)
                    pts.add(b['lo_eff'])
            # support edges: exact zero one binary64 outside every free segment edge, so
            # no linear ramp ever crosses into an unsupported interval (no gap filling)
            for k, s_ in enumerate(segs):
                lo_e = s_.get('lo_eff', s_['lo'])
                if k == 0 or 'lo_eff' not in s_:
                    pts.add(math.nextafter(lo_e, -math.inf))
                joined_right = k + 1 < len(segs) and 'lo_eff' in segs[k + 1]
                if not joined_right:
                    pts.add(math.nextafter(s_['hi'], math.inf))
            # reference energy change and kinematic threshold
            prods = ps['products']
            q = sum(ref[n] for n, _ in prods) - ref[P] - ref[T]
            if ps.get('loss_eV') == 'threshold':
                ps['loss_eV'] = sum(ref[n] for n, _ in prods[1:]) - ref[T]
            if ps['kind'] in ('EXCITE', 'DISSOC'):
                Qk = ps['loss_eV']
            elif ps['kind'] == 'CHARGE_X' and len(prods) > 2:
                Qk = q + ps['loss_eV']   # v3.2: capture with fragmentation of the target-side ion
            elif ps['kind'] == 'BREAKUP':
                Qk = q + ps['loss_eV']   # v3.2: projectile breakup (detachment, stripping, dissociation)
            elif ps['kind'] in ('CAPTURE', 'CHARGE_X', 'IONIZE', 'ION_IONIZE'):
                Qk = q
            else:
                Qk = 0.0
            thr_k = kin_threshold(mp, mt, Qk)
            ps['_q'] = q
            ps['_thr_k'] = thr_k
            if thr_k > 0:
                pts.add(thr_k)
                pts.add(math.nextafter(thr_k, -math.inf))
            for extra in ps.get('extra_grid_eV', []):
                pts.add(extra)
            if ps.get('loss_eV') == 'threshold':
                ps['loss_eV'] = sum(ref[n] for n, _ in prods[1:]) - ref[T]
            proc_items.append((ps, segs))
        xsegs = gspec.get('xi', {}).get('segments', [])
        for si, seg in enumerate(xsegs):
            bps = set([seg['lo'], seg['hi']])
            for key in ('sigma_m', 'sigma_el'):
                so_ = make_source(seg[key], ctx)
                for E in so_.breakpoints():
                    if seg['lo'] <= E <= seg['hi']:
                        bps.add(E)
            bps = sorted(bps)
            # refinement: the ratio of two piecewise-linear functions is not linear; 16 interior points
            # per interval keep the momentum-transfer reproduction error small (recorded per process)
            for a_, b_ in zip(bps, bps[1:]):
                for k_ in range(1, 16):
                    pts.add(a_ * (b_ / a_) ** (k_ / 16.0))
            pts.update(bps)
            # one-ulp steps at the segment edges (xi is discontinuous there)
            pts.add(math.nextafter(seg['lo'], -math.inf))
            if si + 1 == len(xsegs):
                pts.add(math.nextafter(seg['hi'], math.inf))
        grid = sorted(E for E in pts if bottom <= E <= top)
        # strictly ascending, no duplicates
        g2 = []
        for E in grid:
            if not g2 or E > g2[-1]:
                g2.append(E)
        grid = g2
        n = len(grid)
        # xi row
        xi_row = None
        if 'xi' in gspec:
            xi_row = [0.0] * n
            xsegs = gspec['xi']['segments']
            for si, seg in enumerate(xsegs):
                sm = make_source(seg['sigma_m'], ctx)
                se = make_source(seg['sigma_el'], ctx)
                last = si + 1 == len(xsegs)
                for i, E in enumerate(grid):
                    inside = seg['lo'] <= E < seg['hi'] or (last and E == seg['hi'])
                    if inside and sm.lo() <= E <= sm.hi() and se.lo() <= E <= se.hi():
                        r = sm(E) / se(E)
                        xi_row[i] = solve_xi(r)
        sig_rows, cover_rows, bound_rows, bound_names = [], [], [], []
        for ps, segs in proc_items:
            row = []
            cls = []
            thr_k = ps['_thr_k']
            for E in grid:
                val, c = 0.0, None
                for s in segs:
                    lo = s.get('lo_eff', s['lo'])
                    if lo <= E <= s['hi']:
                        val = s['src'](E)
                        c = s['class']
                        break
                if c is None:
                    first_lo = segs[0].get('lo_eff', segs[0]['lo'])
                    if E < first_lo:
                        c = 'ZERO' if E < ps.get('known_zero_below_eV', thr_k) else 'UNRESOLVED'
                    else:
                        c = 'UNRESOLVED'
                if E < thr_k:
                    if val > 0:
                        ps.setdefault('_cut', []).append([E, val])
                    val = 0.0
                if val < 0:
                    raise ValueError('negative cross section %s at %g' % (ps['id'], E))
                row.append(val)
                cls.append(c)
            if ps.get('sigma_is_momentum_transfer'):
                # sampled elastic cross section = sigma_m / (1 - <mu>(xi)), so that the screened
                # angular law reproduces the momentum-transfer cross section exactly at grid points
                if xi_row is None:
                    raise ValueError('momentum-transfer elastic needs a xi row')
                ps['_mt_row'] = list(row)
                row = [r / (1.0 - mean_mu(x)) for r, x in zip(row, xi_row)]
            sig_rows.append(row)
            cover_rows.append(cls)
        # coverage segments per process (energy intervals of constant class)
        proc_cover = []
        for (ps, segs), cls in zip(proc_items, cover_rows):
            # class of interval [grid[i], grid[i+1]): the segment class if both ends lie in
            # one segment, otherwise the class outside the segments
            ivs = []
            # known zero only below the physical threshold (kinematic threshold of the
            # reference energy, or an explicit source statement); between that threshold
            # and the first data point the process is unresolved (no gap filling)
            zb = ps.get('known_zero_below_eV', ps['_thr_k'])
            ps['_zb'] = zb
            if ps.get('constant_rate_below'):
                # sigma v constant below the grid (kernel continuation): same class as the first segment
                ivs.append([0.0, grid[0], segs[0]['class']])
            else:
                ivs.append([0.0, grid[0], 'ZERO' if grid[0] <= zb else 'UNRESOLVED'])
            first_lo = segs[0].get('lo_eff', segs[0]['lo'])

            def seg_of(E):
                for si, s_ in enumerate(segs):
                    if s_.get('lo_eff', s_['lo']) <= E <= s_['hi']:
                        return si
                return None
            for i in range(n - 1):
                a_, b_ = seg_of(grid[i]), seg_of(grid[i + 1])
                if a_ is not None and a_ == b_:
                    c = segs[a_]['class']
                elif a_ is not None and b_ is not None and b_ == a_ + 1 and 'lo_eff' in segs[b_]:
                    c = segs[b_]['class']   # one-ulp junction interval of joined segments
                elif grid[i + 1] <= first_lo:
                    c = 'ZERO' if grid[i + 1] <= zb else 'UNRESOLVED'
                else:
                    c = 'UNRESOLVED'
                if grid[i + 1] <= ps['_thr_k']:
                    c = 'ZERO'
                if ivs[-1][2] == c and ivs[-1][1] == grid[i]:
                    ivs[-1][1] = grid[i + 1]
                else:
                    ivs.append([grid[i], grid[i + 1], c])
            above = ps.get('above_top_class', 'UNRESOLVED')
            if ivs[-1][2] == above:
                ivs[-1][1] = 1e300
            else:
                ivs.append([grid[-1], 1e300, above])
            # merge to at most MAX_COVER
            merged = []
            for iv in ivs:
                if merged and merged[-1][2] == iv[2]:
                    merged[-1][1] = iv[1]
                else:
                    merged.append(iv)
            if len(merged) > MAX_COVER:
                raise ValueError('too many coverage segments for %s: %s' % (ps['id'], merged))
            proc_cover.append(merged)
        # bound rows per process (for unresolved intervals)
        proc_bound_row = []
        for (ps, segs), cls in zip(proc_items, cover_rows):
            b = ps.get('unresolved_bound')
            if b is None:
                proc_bound_row.append(None)
                continue
            so = make_source(b['source'], ctx)
            row = []
            phb = ps.get('placeholder_below_eV', 0.0)
            for E, c in zip(grid, cls):
                if c == 'UNRESOLVED' and so.lo() <= E <= so.hi() and E >= phb:
                    row.append(so(E) * b.get('factor', 1.0))
                else:
                    row.append(0.0)
            proc_bound_row.append(row)
        # coverage bins: intervals (below, 0..n-2, above) -> set of unresolved processes
        bins = []
        def unresolved_set(where):
            s = []
            for j, cov in enumerate(proc_cover):
                E = where
                for lo, hi, c in cov:
                    if lo <= E < hi:
                        if c == 'UNRESOLVED':
                            s.append(j)
                        break
            return tuple(s)
        probe = [grid[0] * 0.5] + [grid[i] for i in range(n - 1)] + [grid[-1] * 2]
        g_region_of = {}
        brow_list = []

        def is_ph(j, E):
            # placeholder (never sampled, exposure recorded, no bound) at energy E: a no-data process, or (v3.2) a process
            # below its placeholder_below_eV edge
            p_ = proc_items[j][0]
            return bool(p_.get('placeholder')) or E < p_.get('placeholder_below_eV', 0.0)
        for bi, E in enumerate(probe):
            us = unresolved_set(E)
            if not us:
                bins.append((NONE32, NONE32))
                continue
            dec = int(math.floor(math.log10(E))) if E > 0 else -99
            ph = tuple(is_ph(j, E) for j in us)
            key = (us, dec, ph)
            if key not in g_region_of:
                g_region_of[key] = len(regions)
                # a process that is a placeholder only in this energy range (below its placeholder_below_eV edge) is marked
                # ' [placeholder here]' in the region name; no-data placeholders carry '(no data)' in their id
                nm = '%s|%s|decade 1e%d eV' % (gname, ','.join(proc_items[j][0]['id'] + (' [placeholder here]' if x and not proc_items[j][0].get('placeholder') else '')
                                                               for j, x in zip(us, ph)), dec)
                if any(nm == r_ for r_ in regions):
                    nm += ' (placeholders: %s)' % ','.join(str(int(x)) for x in ph)
                regions.append(nm)
                region_info.append({'index': len(regions) - 1, 'name': regions[-1], 'group': gname,
                                    'processes': [proc_items[j][0]['id'] for j in us], 'decade_eV': dec,
                                    'placeholder_processes': [proc_items[j][0]['id'] for j, x in zip(us, ph) if x],
                                    'below_grid': bi == 0, 'above_grid': bi == len(probe) - 1,
                                    'placeholders_only': all(ph),
                                    'bounded': (not all(ph)) and all(proc_bound_row[j] is not None for j, x in zip(us, ph) if not x)})
            reg = g_region_of[key]
            # bound row = sum of the bound rows of the non-placeholder unresolved processes, if all of them have one
            # (placeholders are never sampled; their exposure is recorded in the same region)
            real = tuple(j for j, x in zip(us, ph) if not x)
            if real and all(proc_bound_row[j] is not None for j in real):
                bkey = real
                if bkey not in [x[0] for x in brow_list]:
                    brow_list.append((bkey, [sum(proc_bound_row[j][i] for j in real) for i in range(n)]))
                row_idx = [x[0] for x in brow_list].index(bkey)
            else:
                row_idx = NONE32
            bins.append((reg, row_idx))
        # rmq
        S = [sum(sig_rows[j][i] for j in range(len(sig_rows))) for i in range(n)]
        iv_max = []
        for i in range(n - 1):
            iv_max.append(interval_max(grid[i], grid[i + 1], S[i], S[i + 1], mp) * (1 + 2.0 ** -20))
        levels = [iv_max]
        L = 1
        while (1 << L) <= n - 1:
            prev = levels[-1]
            half = 1 << (L - 1)
            cur = [max(prev[i], prev[i + half]) if i + half < len(prev) else prev[i] for i in range(n - 1)]
            levels.append(cur)
            L += 1
        for (ps, segs), row in zip(proc_items, sig_rows):
            if ps.get('constant_rate_below') and not row[0] > 0:
                raise ValueError('constant-rate continuation needs a positive cross section at the grid bottom: %s' % ps['id'])
        low_rate = 1 if any(ps.get('constant_rate_below') for ps, _ in proc_items) else 0
        g = {'low_rate': low_rate, 'projectile': sidx[P], 'target': sidx[T], 'grid_off': len(grid_all), 'grid_n': n, 'proc_off': len(procs_out),
             'proc_n': len(proc_items), 'sigma_off': len(sigma_all), 'bound_off': len(bound_all), 'bound_n': len(brow_list),
             'rmq_off': len(rmq_all), 'rmq_levels': len(levels), 'xi_off': NONE32 if xi_row is None else len(xi_all), 'name': gname}
        grid_all.extend(grid)
        for row in sig_rows:
            sigma_all.extend(row)
        for _, row in brow_list:
            bound_all.extend(row)
        for lv in levels:
            rmq_all.extend(lv)
        if xi_row is not None:
            xi_all.extend(xi_row)
        bin_off.append(len(bin_region_all))
        for reg, row_idx in bins:
            bin_region_all.append(reg)
            bin_row_all.append(row_idx)
        groups_out.append(g)
        # process records
        for j, ((ps, segs), cov) in enumerate(zip(proc_items, proc_cover)):
            prods = ps['products']
            if len(prods) > MAX_OUT:
                raise ValueError('too many products')
            # charge and nuclei balance
            qin = species[sidx[P]]['charge'] + species[sidx[T]]['charge']
            qout = sum(species[sidx[nm]]['charge'] for nm, _ in prods)
            if qin != qout:
                raise ValueError('charge balance %s' % ps['id'])
            for a in NUCLEI:
                nin = species[sidx[P]].get('atoms', {}).get(a, 0) + species[sidx[T]].get('atoms', {}).get(a, 0)
                nout = sum(species[sidx[nm]].get('atoms', {}).get(a, 0) for nm, _ in prods)
                if nin != nout:
                    raise ValueError('nuclei balance %s %s' % (ps['id'], a))
            chain_q = 0.0
            if ps['kind'] == 'DISSOC':
                chain_q = sum(ref[nm] for nm, _ in prods[1:]) - ref[T]
                if ps['loss_eV'] < chain_q:
                    raise ValueError('dissociation loss below the fragment reference energy %s' % ps['id'])
            orbs = ps.get('orbitals', [])
            if ps['kind'] in ('BREAKUP',) or (ps['kind'] == 'CHARGE_X' and len(prods) > 2):
                if not (ps.get('loss_eV', -1.0) >= 0.0):
                    raise ValueError('%s needs loss_eV >= 0 (energy released inside the fragmenting system)' % ps['id'])
            if ps['kind'] == 'BREAKUP' and (len(prods) < 3 or prods[-1][1] not in ('GAS', 'BIRTH')):
                raise ValueError('BREAKUP needs >= 2 projectile fragments and the target last: %s' % ps['id'])
            if ps['kind'] == 'ION_IONIZE':
                if prods[0][1] != 'SURVIVOR' or prods[1][0] != 'e' or not orbs or len(ps.get('rudd', [])) != 10:
                    raise ValueError('ION_IONIZE needs (projectile SURVIVOR, e, ion system...), orbitals and 10 Rudd parameters: %s' % ps['id'])
            rec = {
                'kind': KIND[ps['kind']], 'group': len(groups_out) - 1, 'angular': ANG[ps.get('angular', 'ISO')], 'tag': TAG[ps['tag']],
                'n_out': len(prods), 'internal_sink': SINK[ps.get('sink', 'NONE')], 'n_chain': 0, 'n_orb': len(orbs),
                'out_species': [sidx[nm] for nm, _ in prods], 'out_dest': [DEST[d] for _, d in prods],
                'n_cover': len(cov), 'flags': 1 if ps.get('constant_rate_below') else 0, 'loss_eV': ps.get('loss_eV', 0.0), 'q_eV': ps['_q'],
                'kinetic_fraction': ps.get('kinetic_fraction', 1.0),
                'orb_B': [o['B'] for o in orbs], 'orb_U': [o.get('U', 0.0) for o in orbs], 'orb_N': [o['N'] for o in orbs], 'orb_w': [o['w'] for o in orbs],
                'rudd': list(ps.get('rudd', [0.0] * 10)),
                'cover_lo': [c[0] for c in cov], 'cover_hi': [c[1] for c in cov], 'cover_class': [COVER[c[2]] for c in cov],
                'cover_region': [NONE32] * len(cov), 'chain_q_eV': chain_q,
                'bohr_m': k_const['bohr_m'], 'rydberg_eV': k_const['rydberg_eV'], 'relativistic_from_eV': k_const['relativistic_from_eV'],
            }
            procs_out.append((ps['id'], ps['tag'], rec))
            if ps['kind'] == 'ION_IONIZE':
                rudd_specs[ps['id']] = {'orbitals': orbs, 'rudd': ps['rudd']}
            # manifest + checks
            mp_entry = {'id': ps['id'], 'group': gname, 'kind': ps['kind'], 'tag': ps['tag'], 'products': prods,
                        'q_eV': ps['_q'], 'kinematic_threshold_eV': ps['_thr_k'], 'loss_eV': ps.get('loss_eV'),
                        'angular': ps.get('angular', 'ISO'), 'sink': ps.get('sink', 'NONE'),
                        'kinetic_fraction': ps.get('kinetic_fraction', 1.0), 'orbitals': orbs, 'rudd': ps.get('rudd'),
                        'segments': [{'lo': s.get('lo_eff', s['lo']), 'hi': s['hi'], 'class': s['class'], 'tag': s['tag'], 'source': s['desc']} for s in segs],
                        'coverage': [{'lo': c[0], 'hi': c[1], 'class': c[2]} for c in cov],
                        'known_zero_below_eV': ps['_zb'], 'constant_rate_below_grid': bool(ps.get('constant_rate_below')), 'source_values_below_threshold_set_to_zero': ps.get('_cut', []),
                        'notes': ps.get('notes', ''), 'unresolved_bound': ps.get('unresolved_bound', {}).get('desc') if ps.get('unresolved_bound') else None,
                        'placeholder': bool(ps.get('placeholder')), 'placeholder_below_eV': ps.get('placeholder_below_eV', 0.0)}
            if ps.get('documented_bound'):
                # v3.2: the physics bound of a model estimate (documentation only; not in the blob): bound / sampled sigma
                db = ps['documented_bound']
                so_b = make_source(db['source'], ctx)
                rows_ = []
                for Ed in db.get('energies_eV', [1e3, 1e4, 5e4, 1e5, 1e6, 1e7]):
                    i_ = max(i for i in range(n) if grid[i] <= Ed) if Ed >= grid[0] else None
                    if i_ is None or i_ >= n - 1:
                        continue
                    sv = sig_rows[j][i_] + (sig_rows[j][i_ + 1] - sig_rows[j][i_]) * ((Ed - grid[i_]) / (grid[i_ + 1] - grid[i_]))
                    bv = so_b(Ed) if so_b.lo() <= Ed <= so_b.hi() else None
                    rows_.append({'E_eV': Ed, 'sigma_m2': sv, 'bound_m2': bv, 'bound_over_sigma': (bv / sv if (bv is not None and sv > 0) else None)})
                mp_entry['documented_bound'] = {'desc': db['desc'], 'source': so_b.desc, 'rule': db.get('rule', ''), 'uncertainty': db.get('uncertainty', ''),
                                                'points': rows_}
            # interpolation error of model segments; exactness of tabulated ones
            errs = []
            for s in segs:
                so = s['src']
                lo, hi = s.get('lo_eff', s['lo']), s['hi']
                gi = [i for i, E in enumerate(grid) if lo <= E <= hi]
                worst = 0.0
                for a, b in zip(gi, gi[1:]):
                    Em = 0.5 * (grid[a] + grid[b])
                    if Em < ps['_thr_k']:
                        continue
                    exact = so(Em)
                    lin = sig_rows[j][a] + (sig_rows[j][b] - sig_rows[j][a]) * ((Em - grid[a]) / (grid[b] - grid[a]))
                    if '_mt_row' in ps:
                        # momentum transfer reproduced by (sampled sigma) x (1 - <mu>(xi)) with xi as the kernel interpolates it
                        xm = xi_row[a] + (xi_row[b] - xi_row[a]) * ((Em - grid[a]) / (grid[b] - grid[a]))
                        lin = lin * (1.0 - mean_mu(xm))
                    if exact > 0:
                        worst = max(worst, abs(lin - exact) / exact)
                errs.append({'segment_lo': lo, 'segment_hi': hi, 'type': 'model' if isinstance(so, Model) else 'table',
                             'max_rel_midpoint_error': worst})
                if isinstance(so, Tab):
                    for E, v in so.pts:
                        if lo <= E <= hi and E >= ps['_thr_k'] and grid[0] <= E <= grid[-1]:
                            if '_mt_row' in ps:
                                v = v / (1.0 - mean_mu(xi_row[grid.index(E)]))
                            checks.append({'process': len(procs_out) - 1, 'E': E, 'sigma': v, 'kind': 'breakpoint'})
            mp_entry['interpolation'] = errs
            # seams at junctions
            for a, b in zip(segs, segs[1:]):
                Ej = a['hi']
                left = a['src'](Ej) if a['src'].lo() <= Ej <= a['src'].hi() else None
                bl = b.get('lo_eff', b['lo'])
                right = b['src'](bl) if b['src'].lo() <= bl <= b['src'].hi() else None
                seam = {'process': ps['id'], 'energy_eV': Ej, 'left': left, 'right': right,
                        'relative_jump': (right / left - 1) if left and right else None, 'left_source': a['desc'], 'right_source': b['desc']}
                seams.append(seam)
                mp_entry.setdefault('seams', []).append(seam)
            manifest_procs.append(mp_entry)
        # region ids into process cover records (first matching region per segment for reporting)
        for j in range(len(proc_items)):
            pid = g['proc_off'] + j
            rec = procs_out[pid][2]
            for c in range(rec['n_cover']):
                if rec['cover_class'][c] == COVER['UNRESOLVED']:
                    # region of the bin at the segment start
                    E = max(rec['cover_lo'][c], grid[0] * 0.5)
                    bi = 0 if E < grid[0] else (n if E >= grid[-1] else max(i for i in range(n - 1) if grid[i] <= E) + 1)
                    rec['cover_region'][c] = bins[bi][0]
    # ------------------------------------------------------ write blob
    ns = len(species)
    group_of = [NONE32] * (ns * ns)
    for gi, g in enumerate(groups_out):
        group_of[g['projectile'] * ns + g['target']] = gi
    sec = {}
    sec['meta'] = struct.pack('<4I', ns, len(groups_out), len(procs_out), len(regions))
    sp = b''
    for s in species:
        nuc = [s.get('atoms', {}).get(a, 0) for a in NUCLEI] + [0]
        sp += struct.pack('<ddiI6I4I', s['mass_eV'], s['ref_eV'], s['charge'], SPKIND[s['kind']], *nuc, 0, 0, 0, 0)
    sec['species'] = sp
    gb = b''
    for g in groups_out:
        gb += struct.pack('<16I', g['projectile'], g['target'], g['grid_off'], g['grid_n'], g['proc_off'], g['proc_n'], g['sigma_off'],
                          g['bound_off'], g['bound_n'], g['rmq_off'], g['rmq_levels'], g['xi_off'], g['low_rate'], 0, 0, 0)
    sec['groups'] = gb
    pb = b''
    for _, _, r in procs_out:
        pad = lambda v, k, f: list(v) + [f] * (k - len(v))
        b = struct.pack('<4I', r['kind'], r['group'], r['angular'], r['tag'])
        b += struct.pack('<4I', r['n_out'], r['internal_sink'], r['n_chain'], r['n_orb'])
        b += struct.pack('<6I', *pad(r['out_species'], MAX_OUT, NONE32))
        b += struct.pack('<6I', *pad(r['out_dest'], MAX_OUT, 0))
        b += struct.pack('<4I', r['n_cover'], r['flags'], 0, 0)
        b += struct.pack('<3d', r['loss_eV'], r['q_eV'], r['kinetic_fraction'])
        for key in ('orb_B', 'orb_U', 'orb_N', 'orb_w'):
            b += struct.pack('<4d', *pad(r[key], MAX_ORB, 0.0))
        b += struct.pack('<6d', *pad(r['cover_lo'], MAX_COVER, 0.0))
        b += struct.pack('<6d', *pad(r['cover_hi'], MAX_COVER, 0.0))
        b += struct.pack('<6I', *pad(r['cover_class'], MAX_COVER, 0))
        b += struct.pack('<6I', *pad(r['cover_region'], MAX_COVER, NONE32))
        b += struct.pack('<4d', r['chain_q_eV'], r['bohr_m'], r['rydberg_eV'], r['relativistic_from_eV'])
        b += struct.pack('<10d', *r['rudd'])
        b += struct.pack('<d', 0.0)
        assert len(b) == 512, len(b)
        pb += b
    sec['procs'] = pb
    sec['grid'] = struct.pack('<%dd' % len(grid_all), *grid_all)
    sec['sigma'] = struct.pack('<%dd' % len(sigma_all), *sigma_all)
    sec['bound'] = struct.pack('<%dd' % len(bound_all), *bound_all) if bound_all else struct.pack('<d', 0.0)
    sec['rmq'] = struct.pack('<%dd' % len(rmq_all), *rmq_all)
    sec['xi'] = struct.pack('<%dd' % len(xi_all), *xi_all) if xi_all else struct.pack('<d', 0.0)
    sec['group_of'] = struct.pack('<%dI' % len(group_of), *group_of)
    sec['bin_region'] = struct.pack('<%dI' % len(bin_region_all), *bin_region_all)
    sec['bin_row'] = struct.pack('<%dI' % len(bin_row_all), *bin_row_all)
    sec['bin_off'] = struct.pack('<%dI' % len(bin_off), *bin_off)
    names = [s['name'] for s in species] + [p[0] for p in procs_out] + [p[1] for p in procs_out] + regions
    for nm in names:
        assert '\n' not in nm
    sec['strings'] = ''.join(nm + '\n' for nm in names).encode()
    order = ['meta', 'species', 'groups', 'procs', 'grid', 'sigma', 'bound', 'rmq', 'xi', 'group_of', 'bin_region', 'bin_row', 'bin_off', 'strings']
    header = b'PFCHEMV3' + struct.pack('<II', 1, len(order))
    off = 16 + 40 * len(order)
    off = (off + 63) // 64 * 64
    toc, body = b'', b''
    pos = off
    for nm in order:
        data = sec[nm]
        toc += nm.encode().ljust(24, b'\0') + struct.pack('<QQ', pos, len(data))
        body += data
        padn = (-len(data)) % 64
        body += b'\0' * padn
        pos += len(data) + padn
    blob = header + toc
    blob += b'\0' * (off - len(blob))
    blob += body
    with open(out_prefix + '.blob', 'wb') as f:
        f.write(blob)
    blob_sha = hashlib.sha256(blob).hexdigest()
    manifest = {
        'schema': 'pfc-chem-v3-manifest-1', 'spec_id': spec['id'], 'blob_sha256': blob_sha, 'blob_bytes': len(blob),
        'input_files_sha256': files, 'constants': k_const, 'thermo': spec['thermo'],
        'species': [{k: v for k, v in s.items()} for s in species],
        'groups': [{k: v for k, v in g.items()} for g in groups_out],
        'processes': manifest_procs, 'regions': regions, 'region_info': region_info, 'seams': seams,
        'majorant_safety_factor': 1 + 2.0 ** -20,
    }
    json.dump(manifest, open(out_prefix + '-manifest.json', 'w'), indent=1)
    json.dump({'blob_sha256': blob_sha, 'checks': checks}, open(out_prefix + '-checks.json', 'w'))
    # v3.2: reference ejected-energy distributions of the ION_IONIZE processes (independent quadrature of the
    # Rudd SDCS) at the decade energies 1e3..1e7 eV where the process has a positive cross section
    rref = []
    edges = [0.0, 0.5, 1.0, 2.0, 4.0, 7.0, 11.0, 16.0, 25.0, 40.0, 70.0, 120.0, 250.0, 600.0, 2000.0, 1e7]
    for pi_, (pid, _, r) in enumerate(procs_out):
        if r['kind'] != KIND['ION_IONIZE']:
            continue
        ps = rudd_specs[pid]
        g = groups_out[r['group']]
        gg = grid_all[g['grid_off']:g['grid_off'] + g['grid_n']]
        row = sigma_all[g['sigma_off'] + (pi_ - g['proc_off']) * g['grid_n']:g['sigma_off'] + (pi_ - g['proc_off'] + 1) * g['grid_n']]
        mP, mT = species[g['projectile']]['mass_eV'], species[g['target']]['mass_eV']
        for E in (1e3, 1e4, 1e5, 1e6, 1e7):
            i = min(range(len(gg)), key=lambda k_: abs(math.log(gg[k_] / E)))
            if not row[i] > 0:
                continue
            fr, tot = rudd_bin_fracs(gg[i], mP, mT, mass['e'], ps['orbitals'], ps['rudd'], k_const['bohr_m'], k_const['rydberg_eV'], r['q_eV'], edges)
            rref.append({'proc': pi_, 'E': gg[i], 'edges': edges, 'frac': fr})
    json.dump({'blob_sha256': blob_sha, 'cases': rref}, open(out_prefix + '-rudd-reference.json', 'w'))
    print(blob_sha, len(blob), 'groups', len(groups_out), 'processes', len(procs_out), 'regions', len(regions))


def mean_mu(xi):
    if abs(xi) < 1e-6:
        return xi / 3.0
    return 1.0 / xi - (1 - xi * xi) / (2 * xi * xi) * math.log((1 + xi) / (1 - xi))


def solve_xi(ratio):
    """xi with 1 - <mu>(xi) = ratio (sigma_m / sigma_el), ratio in (0, 2)."""
    target = 1.0 - ratio
    target = max(-1 + 1e-12, min(1 - 1e-12, target))
    lo, hi = -1 + 1e-15, 1 - 1e-15
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mean_mu(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def interval_max(E0, E1, S0, S1, m):
    """Upper bound of (S0 + (S1-S0)(E-E0)/(E1-E0)) v(E) on [E0, E1] (log-concave: ternary search)."""
    f = lambda E: (S0 + (S1 - S0) * ((E - E0) / (E1 - E0))) * speed_tf(E, m)
    best = max(f(E0), f(E1))
    a, b = E0, E1
    for _ in range(200):
        x1 = a + (b - a) / 3
        x2 = b - (b - a) / 3
        if f(x1) < f(x2):
            a = x1
        else:
            b = x2
    best = max(best, f(0.5 * (a + b)))
    return best


if __name__ == '__main__':
    build(sys.argv[1], sys.argv[2])
