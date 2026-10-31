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

KIND = {'ELASTIC': 1, 'EXCITE': 2, 'DISSOC': 3, 'IONIZE': 4, 'CAPTURE': 5, 'CHARGE_X': 6}
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
    manifest_procs, checks = [], []
    seams = []
    for gspec in spec['groups']:
        P, T = gspec['projectile'], gspec['target']
        mp, mt = mass[P], mass[T]
        gname = '%s+%s' % (P, T)
        proc_items = []   # (pspec, segments[(lo, hi, src_obj, class, tag)], cover list, bound src)
        pts = set([bottom, top])
        for ps in gspec['processes']:
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
            elif ps['kind'] in ('CAPTURE', 'CHARGE_X', 'IONIZE'):
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
            for E, c in zip(grid, cls):
                if c == 'UNRESOLVED' and so.lo() <= E <= so.hi():
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
        for bi, E in enumerate(probe):
            us = unresolved_set(E)
            if not us:
                bins.append((NONE32, NONE32))
                continue
            dec = int(math.floor(math.log10(E))) if E > 0 else -99
            key = (us, dec)
            if key not in g_region_of:
                g_region_of[key] = len(regions)
                regions.append('%s|%s|decade 1e%d eV' % (gname, ','.join(proc_items[j][0]['id'] for j in us), dec))
            reg = g_region_of[key]
            # bound row = sum of the processes' bound rows if all exist
            if all(proc_bound_row[j] is not None for j in us):
                bkey = us
                if bkey not in [x[0] for x in brow_list]:
                    brow_list.append((bkey, [sum(proc_bound_row[j][i] for j in us) for i in range(n)]))
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
        g = {'projectile': sidx[P], 'target': sidx[T], 'grid_off': len(grid_all), 'grid_n': n, 'proc_off': len(procs_out),
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
            rec = {
                'kind': KIND[ps['kind']], 'group': len(groups_out) - 1, 'angular': ANG[ps.get('angular', 'ISO')], 'tag': TAG[ps['tag']],
                'n_out': len(prods), 'internal_sink': SINK[ps.get('sink', 'NONE')], 'n_chain': 0, 'n_orb': len(orbs),
                'out_species': [sidx[nm] for nm, _ in prods], 'out_dest': [DEST[d] for _, d in prods],
                'n_cover': len(cov), 'flags': 0, 'loss_eV': ps.get('loss_eV', 0.0), 'q_eV': ps['_q'],
                'kinetic_fraction': ps.get('kinetic_fraction', 1.0),
                'orb_B': [o['B'] for o in orbs], 'orb_U': [o['U'] for o in orbs], 'orb_N': [o['N'] for o in orbs], 'orb_w': [o['w'] for o in orbs],
                'cover_lo': [c[0] for c in cov], 'cover_hi': [c[1] for c in cov], 'cover_class': [COVER[c[2]] for c in cov],
                'cover_region': [NONE32] * len(cov), 'chain_q_eV': chain_q,
                'bohr_m': k_const['bohr_m'], 'rydberg_eV': k_const['rydberg_eV'], 'relativistic_from_eV': k_const['relativistic_from_eV'],
            }
            procs_out.append((ps['id'], ps['tag'], rec))
            # manifest + checks
            mp_entry = {'id': ps['id'], 'group': gname, 'kind': ps['kind'], 'tag': ps['tag'], 'products': prods,
                        'q_eV': ps['_q'], 'kinematic_threshold_eV': ps['_thr_k'], 'loss_eV': ps.get('loss_eV'),
                        'angular': ps.get('angular', 'ISO'), 'sink': ps.get('sink', 'NONE'),
                        'kinetic_fraction': ps.get('kinetic_fraction', 1.0), 'orbitals': orbs,
                        'segments': [{'lo': s.get('lo_eff', s['lo']), 'hi': s['hi'], 'class': s['class'], 'tag': s['tag'], 'source': s['desc']} for s in segs],
                        'coverage': [{'lo': c[0], 'hi': c[1], 'class': c[2]} for c in cov],
                        'known_zero_below_eV': ps['_zb'], 'source_values_below_threshold_set_to_zero': ps.get('_cut', []),
                        'notes': ps.get('notes', ''), 'unresolved_bound': ps.get('unresolved_bound', {}).get('desc') if ps.get('unresolved_bound') else None}
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
                          g['bound_off'], g['bound_n'], g['rmq_off'], g['rmq_levels'], g['xi_off'], 0, 0, 0, 0)
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
        b += struct.pack('<11d', *([0.0] * 11))
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
        'processes': manifest_procs, 'regions': regions, 'seams': seams,
        'majorant_safety_factor': 1 + 2.0 ** -20,
    }
    json.dump(manifest, open(out_prefix + '-manifest.json', 'w'), indent=1)
    json.dump({'blob_sha256': blob_sha, 'checks': checks}, open(out_prefix + '-checks.json', 'w'))
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
