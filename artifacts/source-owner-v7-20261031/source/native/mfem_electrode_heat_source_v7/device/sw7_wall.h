// sw7_wall.h -- wall-interaction kernel of heat/source owner v7 (one wall hit -> outcome).
//
// Compiled by the host C++ compiler and by NVRTC from this one file (see sw7_math.h for the
// floating-point rules). One call evaluates one first-hit event of a macro-particle of weight
// w on a foil patch: every emitted, reflected or re-emitted particle is a macro-particle of the
// same weight w; stuck molecules are credited to the patch's adsorbed-water stock in units of
// w; the stimulated-desorption cross section of the event is returned for the patch hazard.
//
// Randomness: Philox4x32-10 (pfa::philox4x32_10 of pf_arrival/philox_core.h, read in place),
// surface stream layout v7 (keyed by event, independent of rank/partition/order):
//   key     (k0, k1) = seed XOR 0x5357375357375357 (low, high 32 bits)   -- distinct from the
//           arrival owner's key (seed) of the same run
//   counter c0 = block (bits 0..15) | purpose (bits 16..23) | step bits 32..39 (bits 24..31)
//           c1 = step bits 0..31
//           c2 = carrier id bits 0..31,  c3 = carrier id bits 32..63
//   domain  step < 2^40, block < 2^16, purpose 1..255; the map is injective on it.
//   words   block b yields R(2b) = (w1<<32 | w0), R(2b+1) = (w3<<32 | w2); uniforms are
//           pfa::u01_open(R) (open interval, exact conversion).
// Energies: eV in binary64 for kinematics; every ledger energy is an integer number of quanta
// of 2^-30 eV (int64), so per-event identities hold exactly and sums over events are exact and
// independent of order and partition. Charges and counts are integers.
#pragma once
#include "sw7_math.h"
#include "pf_arrival/philox_core.h"

#define SW7_MAXOUT 12
#define SW7_QUANTA_PER_EV 1073741824.0   /* 2^30 */
#define SW7_SIGMA_QUANTA_PER_M2 1267650600228229401496703205376.0   /* 2^100 */
#define SW7_RNG_WORD_CAP 4096u
#define SW7_SALT 0x5357375357375357ull
#define SW7_PURPOSE_WALL 0x21u
#define SW7_PURPOSE_PHOTO 0x22u

// Species indices (fixed; the host species table must follow this order).
#define SW7_E 0
#define SW7_H2O 1
#define SW7_H2OP 2
#define SW7_H3OP 3
#define SW7_OHP 4
#define SW7_OP 5
#define SW7_HP 6
#define SW7_H2P 7
#define SW7_OH 8
#define SW7_O 9
#define SW7_H 10
#define SW7_H2 11
#define SW7_NSPECIES 12

// Output channels.
#define SW7_CH_ELASTIC 1
#define SW7_CH_REDIFFUSED 2
#define SW7_CH_TRUE 3
#define SW7_CH_POTENTIAL 4
#define SW7_CH_KINETIC 5
#define SW7_CH_REFLECTED 6
#define SW7_CH_THERMAL 7
#define SW7_CH_SPECULAR 8
#define SW7_CH_FRAGMENT 9

// Status bits.
#define SW7_ST_ENERGY_CAP 1u      /* emitted count reduced to the energetically allowed maximum */
#define SW7_ST_RNG_BUDGET 2u      /* word budget exhausted (event refused; never resampled) */
#define SW7_ST_INVALID 4u         /* invalid input (event refused) */
#define SW7_ST_YIELD_CLAMP 8u     /* delta_e + delta_r >= 1 at this energy/angle (clamped; recorded) */
#define SW7_ST_OUTSIDE_DATA 16u   /* incident energy outside the cited data range of a channel (model extrapolation; recorded) */

struct SW7Species {
    int kind;                 // 0 electron, 1 positive ion (+1), 2 neutral
    int product;              // neutral species after neutralization (ions), itself for neutrals, -1 electron
    int fragment;             // second neutral product of neutralization (H3O+ -> H2O + H), else -1
    int water;                // 1 if the (neutralized) product is H2O (adsorbable into the stock)
    double mass_kg;
    double recombination_eV;  // ions: energy released by neutralization with a vacuum-level electron (product IE)
    double n_H, n_O;          // constituent atoms (kinetic emission additivity)
    int refl_set;             // Thomas-Janev-Smith constant set: 0 H-like (mu >= 20), 1 O-like (mu ~ 3)
    double refl_epsL_per_keV; // Lindhard reduced energy per keV of the reflecting atom on the target
    double refl_energy_share; // fraction of the projectile energy carried by the reflecting atom
};

struct SW7Params {
    // surface state of the patch (owner supplied at the step)
    double phi_eV, T_K, theta;          // work function, wall temperature, water coverage (monolayers)
    double adsorption_eV;               // binding energy released when a molecule sticks (adsorbed reference)
    // Furman-Pivi electron model (stainless steel table)
    double P1e_inf, P1e_hat, Ee_hat, W, p, sigma_e, e1, e2;
    double P1r_inf, Er, r, q, r1, r2;
    double dts_hat, Ets_hat, s, t1, t2, t3, t4;
    double pn[10], epsn[10];
    int M;
    double see_scale;                   // multiplies delta_ts (surface-condition scenario; 1 = fitted table)
    double fp_fit_E_max;                // upper end of the fitted data (eV); above: backscatter continuation
    double bs_blend_E;                  // energy (eV) at which the backscatter equals the measured eta(E)
    double bs_E[16], bs_eta[16]; int bs_n;   // measured backscatter coefficient eta(E) (Fe), eV
    // ion-induced electron emission
    double pe_a, pe_b, pe_c;            // gamma_PE = pe_a (pe_b E_rec - pe_c phi), when E_rec > 2 phi
    double ke_H_E[16], ke_H_g[16];      // kinetic emission per H atom vs energy per atom (eV), table (log-log interp.)
    double ke_O_E[16], ke_O_g[16];      // per O atom
    int ke_H_n, ke_O_n;
    double ke_scale;
    // ion reflection (Thomas-Janev-Smith 1991/1992 Eq. (4) with reduced-energy constants per mass-ratio group)
    double rn_A[6], re_A[6];            // R_N(eps), R_E(eps) constants, H-like set
    double rn_O_A[6], re_O_A[6];        // O-like set
    double reflect_eps_min, reflect_eps_max;   // fitted range of eps (outside: formula used, event flagged)
    double reflect_scale;
    // neutral-wall
    double stick0;                      // H2O initial sticking coefficient at theta -> 0
    double stick_theta_max;             // coverage at which sticking vanishes (Langmuir-like (1 - theta/theta_max))
    double accommodation;               // thermal accommodation (diffuse fraction) for reflected neutrals (molecules, O, OH)
    double accommodation_light;         // same for H and H2
    // stimulated desorption cross sections (per adsorbed molecule, m^2) vs incident energy (eV)
    double esd_E[16], esd_sigma[16]; int esd_n;    // electrons
    double isd_E[16], isd_sigma[16]; int isd_n;    // ions (per incident ion)
    double esd_scale, isd_scale;
    // channel switches
    int enable_see, enable_iiee, enable_reflection, enable_sticking, enable_esd, enable_isd;
};

struct SW7Hit {
    sw7_u64 id;          // carrier id (event key)
    sw7_u32 species;     // SW7_* index
    sw7_u32 patch;       // patch index (aggregation only)
    double weight;       // macro-particle weight (physical particles)
    double u[3];         // incident proper velocity (m/s), pointing into the wall
    double n[3];         // outward unit normal of the wall (into vacuum), axis aligned
};

struct SW7Out {
    sw7_u32 status, n_out, channel_main, words;
    sw7_u32 sp[SW7_MAXOUT], ch[SW7_MAXOUT];
    sw7_i64 kq[SW7_MAXOUT];          // kinetic energy quanta of each output
    double u[SW7_MAXOUT][3];         // output proper velocity (m/s), pointing into vacuum
    sw7_i64 k_in_q, e_in_q, e_out_q, wall_q, ads_q, phi_q;
    int charge_in, charge_out;
    sw7_u32 stuck;                   // adsorbed molecules in units of the hit weight (0/1)
    sw7_u32 electrons_out;           // emitted electrons
    sw7_i64 sigma_q;                 // stimulated desorption cross section, quanta of 2^-100 m^2
    double k_in_eV, cos_in;
};

namespace sw7 {

struct Stream {
    pfa::Key2 key;
    sw7_u32 c0hi, c1, c2, c3;
    sw7_u32 next;
    bool exhausted;
};
SW7_FN Stream stream(sw7_u64 seed, sw7_u64 step, sw7_u64 id, sw7_u32 purpose) {
    Stream s;
    const sw7_u64 k = seed ^ SW7_SALT;
    s.key.v[0] = (sw7_u32)k;
    s.key.v[1] = (sw7_u32)(k >> 32);
    s.c0hi = ((purpose & 0xFFu) << 16) | ((sw7_u32)((step >> 32) & 0xFFu) << 24);
    s.c1 = (sw7_u32)step;
    s.c2 = (sw7_u32)id;
    s.c3 = (sw7_u32)(id >> 32);
    s.next = 0;
    s.exhausted = false;
    return s;
}
SW7_FN sw7_u64 word(Stream& s) {
    if (s.next >= SW7_RNG_WORD_CAP) { s.exhausted = true; return 0x8000000000000000ull; }
    const sw7_u32 j = s.next++;
    pfa::Block c;
    c.v[0] = (j >> 1) | s.c0hi;
    c.v[1] = s.c1;
    c.v[2] = s.c2;
    c.v[3] = s.c3;
    const pfa::Block b = pfa::philox4x32_10(c, s.key);
    return (j & 1u) ? pfa::high_word(b) : pfa::low_word(b);
}
SW7_FN double u01(Stream& s) { return pfa::u01_open(word(s)); }

#define qe 1.602176634e-19
#define me 9.1093837139e-31
#define kb 1.380649e-23
#define c_light 299792458.0
#define amu 1.66053906660e-27

SW7_FN double dot3(const double* a, const double* b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
// Kinetic energy (J) of mass m with proper velocity u: m u^2 / (gamma + 1).
SW7_FN double kinetic_J(double m, const double* u) {
    const double u2 = dot3(u, u);
    return m * u2 / (sw7_sqrt(1.0 + u2 / (c_light * c_light)) + 1.0);
}
// Proper speed of kinetic energy K (J) and mass m.
SW7_FN double proper_speed(double m, double K) {
    return sw7_sqrt(K * (K + 2.0 * m * c_light * c_light)) / (m * c_light);
}
SW7_FN sw7_i64 quanta(double E_eV) {   // floor(E 2^30) for 0 <= E < 2^32 eV
    if (!(E_eV > 0)) return 0;
    return (sw7_i64)(E_eV * SW7_QUANTA_PER_EV);
}
SW7_FN double from_quanta(sw7_i64 q) { return (double)q / SW7_QUANTA_PER_EV; }
// Orthonormal frame (n, t1, t2) of an axis-aligned outward normal (same construction as the sealed kernel).
SW7_FN void frame(const double* n, double* t1, double* t2) {
    double ax[3] = {1, 0, 0};
    if (!(fabs_(n[0]) < 0.8)) { ax[0] = 0; ax[1] = 1; }
    t1[0] = n[1] * ax[2] - n[2] * ax[1];
    t1[1] = n[2] * ax[0] - n[0] * ax[2];
    t1[2] = n[0] * ax[1] - n[1] * ax[0];
    const double l = sw7_sqrt(dot3(t1, t1));
    t1[0] /= l; t1[1] /= l; t1[2] /= l;
    t2[0] = n[1] * t1[2] - n[2] * t1[1];
    t2[1] = n[2] * t1[0] - n[0] * t1[2];
    t2[2] = n[0] * t1[1] - n[1] * t1[0];
}
// Cosine-law (Lambertian flux) direction: cos(theta) = sqrt(u1), azimuth 2 pi u2.
SW7_FN void cosine_direction(Stream& s, const double* n, double* d) {
    double t1[3], t2[3];
    frame(n, t1, t2);
    const double a = u01(s), b = u01(s);
    const double ct = sw7_sqrt(a), st = sw7_sqrt(1.0 - a);
    double sn, cs;
    sincos_turn(b, &sn, &cs);
    for (int k = 0; k < 3; ++k) d[k] = ct * n[k] + st * (cs * t1[k] + sn * t2[k]);
}
SW7_FN void standard_normal_pair(Stream& s, double* z0, double* z1) {
    const double a = u01(s), b = u01(s);
    const double r = sw7_sqrt(-2.0 * log_(a));
    double sn, cs;
    sincos_turn(b, &sn, &cs);
    *z0 = r * cs;
    *z1 = r * sn;
}
// Gamma(a, 1) variate (Marsaglia-Tsang 2000; a < 1 by the a + 1 boost).
SW7_FN double gamma_variate(Stream& s, double a) {
    double boost = 1.0;
    if (a < 1.0) { boost = pow_(u01(s), 1.0 / a); a += 1.0; }
    const double d = a - 1.0 / 3.0, c = 1.0 / sw7_sqrt(9.0 * d);
    for (int it = 0; it < 256 && !s.exhausted; ++it) {
        double x, spare;
        standard_normal_pair(s, &x, &spare);
        double v = 1.0 + c * x;
        if (v <= 0) continue;
        v = v * v * v;
        const double u = u01(s);
        const double x2 = x * x;
        if (u < 1.0 - 0.0331 * x2 * x2) return d * v * boost;
        if (log_(u) < 0.5 * x2 + d * (1.0 - v + log_(v))) return d * v * boost;
    }
    s.exhausted = true;
    return 0.0;
}

SW7_FN double interp_loglog(const double* X, const double* Y, int n, double x, bool* outside) {
    if (n <= 0) return 0.0;
    if (x < X[0]) return 0.0;                                                                     // below the first point (threshold): 0
    if (x == X[0]) return Y[0];
    if (x >= X[n - 1]) { if (x > X[n - 1]) *outside = true; return Y[n - 1]; }                  // above: held (bounded, recorded)
    int i = 0;
    while (i + 1 < n && X[i + 1] < x) ++i;
    if (x == X[i + 1]) return Y[i + 1];                                                           // table node: exact
    if (!(Y[i] > 0) || !(Y[i + 1] > 0)) return Y[i] + (Y[i + 1] - Y[i]) * (x - X[i]) / (X[i + 1] - X[i]);
    const double t = (log_(x) - log_(X[i])) / (log_(X[i + 1]) - log_(X[i]));
    return exp_(log_(Y[i]) + t * (log_(Y[i + 1]) - log_(Y[i])));
}
// ---- Furman-Pivi (2002) yields
SW7_FN double fp_D(double x, double s_) { return s_ * x / (s_ - 1.0 + pow_(x, s_)); }
SW7_FN void fp_yields(const SW7Params& P, double E0, double cth, double* de, double* dr, double* dts) {
    const double c = fmin_(fmax_(cth, 0.0), 1.0);
    const double ae = fabs_(E0 - P.Ee_hat) / P.W;
    double e = P.P1e_inf + (P.P1e_hat - P.P1e_inf) * exp_(-pow_(ae, P.p) / P.p);
    e = e * (1.0 + P.e1 * (1.0 - pow_(c, P.e2)));
    double r = P.P1r_inf * (1.0 - exp_(-pow_(E0 / P.Er, P.r)));
    r = r * (1.0 + P.r1 * (1.0 - pow_(c, P.r2)));
    const double dhat = P.dts_hat * (1.0 + P.t1 * (1.0 - pow_(c, P.t2)));
    const double Ehat = P.Ets_hat * (1.0 + P.t3 * (1.0 - pow_(c, P.t4)));
    const double t = E0 > 0 ? dhat * fp_D(E0 / Ehat, P.s) * P.see_scale : 0.0;
    *de = e; *dr = r; *dts = t;
}
// Above the fitted range the FP backscatter plateau (P1e_inf + P1r_inf) is replaced by the measured
// backscatter coefficient eta(E): log-linear blend in E from the FP value at fp_fit_E_max to eta at
// bs_blend_E, eta(E) above (held beyond the last table point). The elastic/rediffused ratio is kept.
SW7_FN double bs_continuation(const SW7Params& P, double E0, double de, double dr, bool* outside) {
    if (!(E0 > P.fp_fit_E_max) || P.bs_n <= 0) return 1.0;
    *outside = true;
    const double fp = de + dr;
    if (!(fp > 0)) return 1.0;
    bool o2 = false;
    const double eta = interp_loglog(P.bs_E, P.bs_eta, P.bs_n, E0 < P.bs_E[0] ? P.bs_E[0] : E0, &o2);
    double b = log_(E0 / P.fp_fit_E_max) / log_(P.bs_blend_E / P.fp_fit_E_max);
    if (b > 1) b = 1;
    const double target = exp_((1.0 - b) * log_(fp) + b * log_(eta));
    return target / fp;
}

// Sample n Dirichlet(p, .., p) fractions times a total S ~ Gamma(n p, eps) truncated to [0, Amax]:
// exactly the FP joint density prod E_k^(p-1) exp(-E_k/eps) restricted to sum E_k <= Amax.
SW7_FN void fp_true_energies(Stream& s, double p, double eps, int n, double Amax, double* E) {
    const double a = n * p;
    const double xmax = Amax / eps;
    const double Pmax = gammainc_P(a, xmax);
    double S;
    if (!(Pmax > 1e-300)) S = Amax * u01(s);   // degenerate (Amax ~ 0): bounded, total <= Amax
    else S = eps * gammainc_P_inverse(a, u01(s) * Pmax, xmax);
    if (S > Amax) S = Amax;
    double g[SW7_MAXOUT], sum = 0;
    for (int k = 0; k < n; ++k) { g[k] = gamma_variate(s, p); sum += g[k]; }
    for (int k = 0; k < n; ++k) E[k] = sum > 0 ? S * (g[k] / sum) : S / n;
}

SW7_FN void emit(SW7Out* o, sw7_u32 species, sw7_u32 channel, sw7_i64 kq, double mass, const double* dir) {
    if (o->n_out >= SW7_MAXOUT) { o->status |= SW7_ST_INVALID; return; }
    const sw7_u32 i = o->n_out++;
    o->sp[i] = species;
    o->ch[i] = channel;
    o->kq[i] = kq;
    const double K = from_quanta(kq) * qe;
    const double us = K > 0 ? proper_speed(mass, K) : 0.0;
    for (int k = 0; k < 3; ++k) o->u[i][k] = us * dir[k];
}
// Maxwell flux at wall temperature: normal sqrt(-2kT/m ln a), tangential sqrt(-2kT/m ln b) at azimuth c
// (the v6.1 desorption packet law); returns the direction and kinetic energy (J, nonrelativistic).
SW7_FN double thermal_flux_direction(Stream& s, const double* n, double T, double m, double* dir) {
    double t1[3], t2[3];
    frame(n, t1, t2);
    const double a = u01(s), b = u01(s), c = u01(s);
    const double vn = sw7_sqrt(-2.0 * kb * T / m * log_(a)), vt = sw7_sqrt(-2.0 * kb * T / m * log_(b));
    double sn, cs;
    sincos_turn(c, &sn, &cs);
    double v[3];
    for (int k = 0; k < 3; ++k) v[k] = vn * n[k] + vt * (cs * t1[k] + sn * t2[k]);
    const double v2 = dot3(v, v), vm = sw7_sqrt(v2);
    for (int k = 0; k < 3; ++k) dir[k] = vm > 0 ? v[k] / vm : n[k];
    return 0.5 * m * v2;
}
// Thomas-Janev-Smith type particle/energy reflection coefficient: A0 ln(A1 eps + e) / (1 + A2 eps^A3 + A4 eps^A5).
SW7_FN double tjs(const double* A, double eps) {
    if (!(eps > 0)) return 0.0;
    return A[0] * log_(A[1] * eps + 2.718281828459045) / (1.0 + A[2] * pow_(eps, A[3]) + A[4] * pow_(eps, A[5]));
}

// ---- one wall event
SW7_FN void wall_event(const SW7Params& P, const SW7Species* table, const SW7Hit& h, sw7_u64 seed, sw7_u64 step, SW7Out* o) {
    o->status = 0; o->n_out = 0; o->channel_main = 0; o->words = 0;
    o->k_in_q = o->e_in_q = o->e_out_q = o->wall_q = o->ads_q = o->phi_q = 0;
    o->charge_in = o->charge_out = 0; o->stuck = 0; o->electrons_out = 0; o->sigma_q = 0;
    o->k_in_eV = 0; o->cos_in = 0;
    for (int i = 0; i < SW7_MAXOUT; ++i) { o->sp[i] = 0; o->ch[i] = 0; o->kq[i] = 0; o->u[i][0] = o->u[i][1] = o->u[i][2] = 0; }
    if (h.species >= SW7_NSPECIES || !(h.weight > 0) || !finite(h.weight) || step >= 0x10000000000ull) { o->status |= SW7_ST_INVALID; return; }
    const SW7Species sp = table[h.species];
    const double nn = dot3(h.n, h.n);
    const double un = dot3(h.u, h.n);
    if (!(fabs_(nn - 1.0) < 1e-12) || !finite(dot3(h.u, h.u)) || !(un < 0)) { o->status |= SW7_ST_INVALID; return; }
    Stream s = stream(seed, step, h.id, SW7_PURPOSE_WALL);
    const double K = kinetic_J(sp.mass_kg, h.u);
    const double E0 = K / qe;
    const double umag = sw7_sqrt(dot3(h.u, h.u));
    const double cth = -un / umag;
    o->k_in_eV = E0; o->cos_in = cth;
    const double phi = P.phi_eV;
    const sw7_i64 phi_q = quanta(phi);
    o->phi_q = phi_q;
    o->k_in_q = quanta(E0);
    sw7_i64 avail_q = 0;   // energy above the reference available to the outputs
    if (sp.kind == 0) {
        // ---------------- electron
        o->charge_in = -1;
        o->e_in_q = o->k_in_q + phi_q;
        if (P.enable_esd && P.esd_n > 0) {
            bool out = false;
            const double sg = interp_loglog(P.esd_E, P.esd_sigma, P.esd_n, E0, &out) * P.esd_scale;
            o->sigma_q = (sw7_i64)(sg * SW7_SIGMA_QUANTA_PER_M2);
            if (out) o->status |= SW7_ST_OUTSIDE_DATA;
        }
        if (P.enable_see) {
            double de, dr, dts;
            fp_yields(P, E0, cth, &de, &dr, &dts);
            bool outside = false;
            const double f_bs = bs_continuation(P, E0, de, dr, &outside);
            de *= f_bs; dr *= f_bs;
            if (outside) o->status |= SW7_ST_OUTSIDE_DATA;
            if (de + dr >= 1.0) { const double f = 0.999999 / (de + dr); de *= f; dr *= f; o->status |= SW7_ST_YIELD_CLAMP; }
            const double u0 = u01(s);
            double dir[3];
            if (u0 < de) {
                // elastic: Gaussian of width sigma_e below E0, truncated to [0, E0]
                double E = -1;
                for (int it = 0; it < 64 && !(E >= 0 && E <= E0); ++it) {
                    double z0, z1;
                    standard_normal_pair(s, &z0, &z1);
                    E = E0 - P.sigma_e * fabs_(z0);
                }
                if (!(E >= 0)) E = 0;
                cosine_direction(s, h.n, dir);
                sw7_i64 kq = quanta(E);
                if (kq > o->k_in_q) kq = o->k_in_q;
                emit(o, SW7_E, SW7_CH_ELASTIC, kq, me, dir);
                o->channel_main = SW7_CH_ELASTIC;
            } else if (u0 < de + dr) {
                const double E = E0 * pow_(u01(s), 1.0 / (P.q + 1.0));
                cosine_direction(s, h.n, dir);
                sw7_i64 kq = quanta(E);
                if (kq > o->k_in_q) kq = o->k_in_q;
                emit(o, SW7_E, SW7_CH_REDIFFUSED, kq, me, dir);
                o->channel_main = SW7_CH_REDIFFUSED;
            } else {
                double pts = dts / (1.0 - de - dr) / P.M;   // per-trial probability of the binomial (FP Eqs. 39, 46)
                if (pts > 1.0) { pts = 1.0; o->status |= SW7_ST_YIELD_CLAMP; }   // delta'_ts > M: clamped at M (recorded)
                int n = 0;
                for (int k = 0; k < P.M; ++k) n += (u01(s) < pts) ? 1 : 0;
                // energy: n electrons need n (phi + K_k) <= E0 + phi  ->  n_max = floor(E0/phi) + 1
                const int nmax = phi > 0 ? (int)floor_(E0 / phi) + 1 : n;
                if (n > nmax) { n = nmax; o->status |= SW7_ST_ENERGY_CAP; }
                if (n > 10) n = 10;
                while (n > 1 && (sw7_i64)(n - 1) * phi_q > o->k_in_q) { --n; o->status |= SW7_ST_ENERGY_CAP; }
                if (n > 0) {
                    const sw7_i64 Aq = o->k_in_q - (sw7_i64)(n - 1) * phi_q;
                    double E[SW7_MAXOUT];
                    fp_true_energies(s, P.pn[n - 1], P.epsn[n - 1], n, from_quanta(Aq), E);
                    sw7_i64 kq[SW7_MAXOUT], tot = 0;
                    int big = 0;
                    for (int k = 0; k < n; ++k) { kq[k] = quanta(E[k]); tot += kq[k]; if (kq[k] > kq[big]) big = k; }
                    if (tot > Aq) kq[big] -= (tot - Aq);   // exact cap: sum of quanta <= available quanta
                    for (int k = 0; k < n; ++k) { cosine_direction(s, h.n, dir); emit(o, SW7_E, SW7_CH_TRUE, kq[k], me, dir); }
                }
                o->channel_main = n > 0 ? SW7_CH_TRUE : 0;
            }
        }
    } else if (sp.kind == 1) {
        // ---------------- positive ion: neutralized at the surface by a wall electron
        o->charge_in = 1;
        const double Erec = sp.recombination_eV;
        o->e_in_q = o->k_in_q + quanta(Erec) - phi_q;   // K + (E_rec - phi)
        if (P.enable_isd && P.isd_n > 0 && sp.water) {
            bool out = false;
            const double sg = interp_loglog(P.isd_E, P.isd_sigma, P.isd_n, E0, &out) * P.isd_scale;
            o->sigma_q = (sw7_i64)(sg * SW7_SIGMA_QUANTA_PER_M2);
            if (out) o->status |= SW7_ST_OUTSIDE_DATA;
        }
        avail_q = o->e_in_q;
        double dir[3];
        if (P.enable_iiee) {
            // potential emission (Auger neutralization): requires E_rec > 2 phi; electrons up to E_rec - 2 phi
            double gpe = 0;
            if (Erec > 2.0 * phi) gpe = P.pe_a * (P.pe_b * Erec - P.pe_c * phi);
            if (gpe < 0) gpe = 0;
            // kinetic emission: additivity over constituent atoms at equal velocity
            bool out = false;
            const double mtot = sp.mass_kg / amu;
            const double EH = E0 * (1.00794 / mtot), EO = E0 * (15.9994 / mtot);
            double gke = 0;
            if (sp.n_H > 0) gke += sp.n_H * interp_loglog(P.ke_H_E, P.ke_H_g, P.ke_H_n, EH, &out);
            if (sp.n_O > 0) gke += sp.n_O * interp_loglog(P.ke_O_E, P.ke_O_g, P.ke_O_n, EO, &out);
            gke *= P.ke_scale;
            if (out) o->status |= SW7_ST_OUTSIDE_DATA;
            // counts: stochastic rounding of the means (exact means)
            const int npe = (int)floor_(gpe) + (u01(s) < gpe - floor_(gpe) ? 1 : 0);
            int nke = (int)floor_(gke) + (u01(s) < gke - floor_(gke) ? 1 : 0);
            if (nke > 8) nke = 8;
            for (int k = 0; k < npe && o->n_out < SW7_MAXOUT - 3; ++k) {
                const double Emax = Erec - 2.0 * phi;   // Auger upper edge
                const sw7_i64 kq = quanta(Emax * u01(s));
                if (kq + phi_q > avail_q) { o->status |= SW7_ST_ENERGY_CAP; break; }
                cosine_direction(s, h.n, dir);
                emit(o, SW7_E, SW7_CH_POTENTIAL, kq, me, dir);
                avail_q -= kq + phi_q;
            }
            for (int k = 0; k < nke && o->n_out < SW7_MAXOUT - 3; ++k) {
                // kinetic-emission electron spectrum: FP single true-secondary spectrum (p_1, eps_1), truncated
                double E1;
                const double Amax = from_quanta(avail_q - phi_q);
                if (!(Amax > 0)) { o->status |= SW7_ST_ENERGY_CAP; break; }
                fp_true_energies(s, P.pn[0], P.epsn[0], 1, Amax, &E1);
                sw7_i64 kq = quanta(E1);
                if (kq + phi_q > avail_q) kq = avail_q - phi_q;
                cosine_direction(s, h.n, dir);
                emit(o, SW7_E, SW7_CH_KINETIC, kq, me, dir);
                avail_q -= kq + phi_q;
            }
        }
        // reflection of the (neutralized) projectile
        const SW7Species prod = table[sp.product];
        bool reflected = false;
        if (P.enable_reflection) {
            const bool hlike = sp.refl_set == 0;
            const double eps = sp.refl_epsL_per_keV * (E0 * sp.refl_energy_share * 1e-3);
            if (eps < P.reflect_eps_min || eps > P.reflect_eps_max) o->status |= SW7_ST_OUTSIDE_DATA;
            const double RN = fmin_(fmax_(tjs(hlike ? P.rn_A : P.rn_O_A, eps) * P.reflect_scale, 0.0), 0.99);
            const double RE = fmin_(fmax_(tjs(hlike ? P.re_A : P.re_O_A, eps) * P.reflect_scale, 0.0), RN);
            if (u01(s) < RN) {
                const double Er = RN > 0 ? E0 * (RE / RN) : 0;
                sw7_i64 kq = quanta(Er);
                if (kq > avail_q) kq = avail_q > 0 ? avail_q : 0;
                cosine_direction(s, h.n, dir);
                emit(o, (sw7_u32)sp.product, SW7_CH_REFLECTED, kq, prod.mass_kg, dir);
                avail_q -= kq;
                reflected = true;
            }
        }
        if (!reflected) {
            // thermalized product: water sticks (stock) or is re-emitted thermally; other products re-emitted
            bool stuck = false;
            if (prod.water && P.enable_sticking) {
                const double sc = P.stick0 * fmax_(0.0, 1.0 - P.theta / P.stick_theta_max);
                if (u01(s) < sc) stuck = true;
            }
            if (stuck) {
                o->stuck = 1;
                o->ads_q = quanta(P.adsorption_eV);
            } else {
                const double Kt = thermal_flux_direction(s, h.n, P.T_K, prod.mass_kg, dir);
                const sw7_i64 kq = quanta(Kt / qe);
                emit(o, (sw7_u32)sp.product, SW7_CH_THERMAL, kq, prod.mass_kg, dir);
            }
        }
        o->channel_main = reflected ? SW7_CH_REFLECTED : (o->stuck ? 0 : SW7_CH_THERMAL);
        if (sp.fragment >= 0) {
            const SW7Species fr = table[sp.fragment];
            const double Kt = thermal_flux_direction(s, h.n, P.T_K, fr.mass_kg, dir);
            emit(o, (sw7_u32)sp.fragment, SW7_CH_FRAGMENT, quanta(Kt / qe), fr.mass_kg, dir);
        }
    } else {
        // ---------------- neutral
        o->e_in_q = o->k_in_q;
        bool stuck = false;
        if (sp.water && P.enable_sticking) {
            const double sc = P.stick0 * fmax_(0.0, 1.0 - P.theta / P.stick_theta_max);
            if (u01(s) < sc) stuck = true;
        }
        double dir[3];
        if (stuck) {
            o->stuck = 1;
            o->ads_q = quanta(P.adsorption_eV);
            o->channel_main = 0;
        } else if (u01(s) < ((h.species == SW7_H || h.species == SW7_H2) ? P.accommodation_light : P.accommodation)) {
            const double Kt = thermal_flux_direction(s, h.n, P.T_K, sp.mass_kg, dir);
            emit(o, h.species, SW7_CH_THERMAL, quanta(Kt / qe), sp.mass_kg, dir);
            o->channel_main = SW7_CH_THERMAL;
        } else {
            // specular: tangential velocity kept, normal reversed (kinetic energy unchanged)
            for (int k = 0; k < 3; ++k) dir[k] = (h.u[k] - 2.0 * un * h.n[k]) / umag;
            emit(o, h.species, SW7_CH_SPECULAR, o->k_in_q, sp.mass_kg, dir);
            o->channel_main = SW7_CH_SPECULAR;
        }
    }
    // ledger: outputs relative to the reference (vacuum electron: +phi; adsorbed molecule: -E_ads)
    sw7_i64 out_q = -o->ads_q;
    int qout = 0;
    sw7_u32 ne = 0;
    for (sw7_u32 i = 0; i < o->n_out; ++i) {
        out_q += o->kq[i];
        if (o->sp[i] == SW7_E) { out_q += phi_q; qout -= 1; ++ne; }
    }
    o->e_out_q = out_q;
    o->wall_q = o->e_in_q - out_q;
    o->charge_out = qout;
    o->electrons_out = ne;
    o->words = s.next;
    if (s.exhausted) o->status |= SW7_ST_RNG_BUDGET;
}

}  // namespace sw7
#undef qe
#undef me
#undef kb
#undef c_light
#undef amu
