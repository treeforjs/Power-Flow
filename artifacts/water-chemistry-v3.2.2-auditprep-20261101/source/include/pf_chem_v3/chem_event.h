// chem_event.h -- chemistry v3 per-carrier collision step (null-collision
// Monte Carlo), shared by the CPU reference and the NVRTC kernel.
//
// Model of one step (operator split; the carrier sits at its post-push
// position in one cell during the collision phase):
//  * The carrier's collision process is a Poisson process whose rate depends
//    only on its current state. Candidate times are drawn with the majorant
//    rate nu_max(state) = sum_k n_k M_k(state), M_k an upper bound of
//    sigma_tot(E') v' over every target-frame energy E' reachable with a target
//    velocity inside the truncated thermal sphere |u_t - drift| <= trunc*s_k.
//    Waiting time tau = -ln(U0)/nu_max (fresh exponential after every
//    candidate; exact by memorylessness, and the majorant is re-evaluated
//    whenever a real event changes the state).
//  * At a candidate: target slot k with probability n_k M_k / nu_max (U1); a
//    target velocity from the truncated Maxwellian (Box-Muller, U4..U19);
//    real/null and channel with one uniform U2: channel j is chosen iff
//    U2 M_k falls in [sum_{i<j} sigma_i v', sum_{i<=j} sigma_i v'); null if
//    U2 M_k >= sigma_tot v'.
//  * Every number is drawn from the chemistry stream (seed, step, carrier id,
//    candidate index, word), so results are independent of the order and
//    partition in which carriers are processed.
// The kernel writes one CarrierOut and up to max_events EventRec per carrier.
// Ledgers, target-inventory debits and neutral releases are formed on the
// host from these records (chem_host.hpp).
#pragma once
#include "chem_kinematics.h"
#include "chem_rng.h"
#include "chem_tables.h"

namespace pfc {

// Carrier status codes.
#define PFC_ST_OK 0u
#define PFC_ST_CONSUMED 1u          // projectile removed by a real event (still OK)
#define PFC_ST_CAND_CAP 2u          // candidate cap reached before dt was exhausted (refuse step)
#define PFC_ST_EVENT_CAP 3u         // real-event record cap reached (refuse step)
#define PFC_ST_KEY_DOMAIN 4u        // id/step/event outside the stream domain (refuse step)
#define PFC_ST_MAJORANT 5u          // sigma_tot v' exceeded the majorant (refuse step; never expected)
#define PFC_ST_COVER_CAP 6u         // more distinct unresolved regions than slots (refuse step)
#define PFC_ST_BAD_INPUT 7u         // nonfinite input, unknown species or cell (refuse step)

#define PFC_MAX_COVER_SLOTS 16

struct StepConfig {
    pfc_u64 seed, step;
    double dt_s;
    double trunc_sigmas;      // truncation radius of the target Maxwellian in standard deviations
    pfc_u32 purpose;
    pfc_u32 max_candidates;   // per carrier per step
    pfc_u32 max_events;       // real events per carrier per step (record slots)
    pfc_u32 pad;
};
struct CarrierIn {           // 48 bytes
    pfc_u64 id;
    double u[3];             // proper velocity, m/s
    double weight;           // physical particles represented
    pfc_u32 species, cell;
};
struct CellTargets {         // targets present in one cell
    double n[PFC_MAX_TARGETS];           // m^-3
    double drift[PFC_MAX_TARGETS][3];    // mean proper velocity, m/s
    double s[PFC_MAX_TARGETS];           // 1D thermal standard deviation sqrt(kT/m), m/s
    pfc_u32 species[PFC_MAX_TARGETS];
    pfc_u32 n_targets, pad;
};
struct EventRec {            // one real event (sizes asserted on the host)
    pfc_u64 carrier;
    pfc_u32 candidate, process, cell, target_species;
    pfc_u32 n_out, cover_class, internal_sink, flags;
    double weight;
    double time_s;           // event time within the step
    double u_in[3];          // projectile proper velocity before (m/s)
    double ut[3];            // sampled target proper velocity (m/s)
    double E_tf;             // projectile kinetic energy in the target frame (eV)
    double ke_in_p, ke_in_t; // canonical kinetic energies (eV, per physical particle)
    double q_eV, e_internal; // reference-energy change and internal sink (eV)
    pfc_u32 out_species[PFC_MAX_OUT];
    pfc_u32 out_dest[PFC_MAX_OUT];
    double out_u[PFC_MAX_OUT][3];
    double ke_out[PFC_MAX_OUT];
};
struct CarrierOut {
    pfc_u64 id;
    pfc_u32 status, n_events, n_candidates, species_final;
    double u_final[3];
    double nu_dt_max;            // max over the step of nu_max * dt (time-step guard quantity)
    pfc_u32 n_cover, pad;
    pfc_u32 cover_region[PFC_MAX_COVER_SLOTS];
    double cover_exposure[PFC_MAX_COVER_SLOTS];   // sum W n v dt (sigma-free), m^-2
    double cover_bound[PFC_MAX_COVER_SLOTS];      // sum W n v sigma_bound dt (expected events bound)
};

// ---------------------------------------------------------------- BEB / RBEB
// One orbital, Q = 1 (Kim & Rudd 1994; Kim, Santos & Parente 2000).
struct SdcsShape {
    double t, phi, cb, L, P;    // P in m^2
    double B;
};
PFC_FN double beta2_of(double x) { return x * (x + 2.0) / ((1.0 + x) * (1.0 + x)); }
PFC_FN SdcsShape sdcs_shape(double B, double U, double N, double T, double a0, double R, double mc2, int rel) {
    SdcsShape s;
    s.B = B;
    s.t = T / B;
    const double pi = 3.141592653589793;
    if (!rel) {
        s.phi = 1.0;
        s.cb = 0.0;
        s.L = pfc::log(s.t);
        const double u = U / B;
        const double S = 4.0 * pi * a0 * a0 * N * (R / B) * (R / B);
        s.P = S / (s.t + u + 1.0);
    } else {
        const double tp = T / mc2, bp = B / mc2, up = U / mc2;
        const double bt2 = beta2_of(tp), bb2 = beta2_of(bp), bu2 = beta2_of(up);
        const double alpha2 = 2.0 * R / mc2;
        s.P = 4.0 * pi * a0 * a0 * alpha2 * alpha2 * N / ((bt2 + bu2 + bb2) * 2.0 * bp);
        s.L = pfc::log(tp * (tp + 2.0)) - bt2 - pfc::log(2.0 * bp);
        const double h = 1.0 + tp / 2.0;
        s.phi = (1.0 + 2.0 * tp) / (h * h);
        s.cb = bp * bp / (h * h);
    }
    return s;
}
PFC_FN double sdcs_G(const SdcsShape& s, double w) {
    const double t = s.t, tw = t - w;
    return -s.phi / (t + 1.0) * (pfc::log1p(w) - pfc::log1p(-w / t)) + w / (w + 1.0) + w / (t * tw) + s.cb * w +
           s.L / 2.0 * (w * (w + 2.0) / ((w + 1.0) * (w + 1.0)) + w * (2.0 * t - w) / (t * t * tw * tw));
}
PFC_FN double sdcs_total(const SdcsShape& s) { return s.P * sdcs_G(s, (s.t - 1.0) / 2.0); }
// Reduced ejected energy w in [0, (t-1)/2] with G(w) = u G(wmax); 64 bisections.
PFC_FN double sdcs_sample_w(const SdcsShape& s, double u) {
    const double wm = (s.t - 1.0) / 2.0;
    const double target = u * sdcs_G(s, wm);
    double lo = 0.0, hi = wm;
    for (int i = 0; i < 64; ++i) {
        const double mid = lo + (hi - lo) / 2.0;
        if (!(mid > lo && mid < hi)) break;
        if (sdcs_G(s, mid) < target) lo = mid;
        else hi = mid;
    }
    return hi;
}

// ------------------------------------------------------------------ Rudd
// Singly differential cross section for ionization of one orbital (binding
// energy B, occupation N) by a heavy projectile of speed v_p (Rudd, Kim,
// Madison & Gay, Rev. Mod. Phys. 64, 441 (1992)):
//   dsigma/dW = (S/B) h(w) c(w),  w = W/B,  S = 4 pi a0^2 N (R/B)^2,
//   h(w) = (F1 + F2 w) / (1 + w)^3,  c(w) = 1 / (1 + exp(alpha (w - wc) / v)),
//   v = sqrt(T_e / B), T_e = m_e v_p^2 / 2 (electron with the projectile speed),
//   wc = 4 v^2 - 2 v - R / (4 B),
//   F1 = L1 + H1, F2 = L2 H2 / (L2 + H2),
//   L1 = C1 v^D1 / (1 + E1 v^(D1 + 4)), H1 = A1 ln(1 + v^2) / (v^2 + B1 / v^2),
//   L2 = C2 v^D2, H2 = A2 / v^2 + B2 / v^4.
// Sampling of w on [0, wmax] (composition-rejection over orbitals), with one of two envelopes per orbital,
// whichever has the smaller mass:
//  (P) power envelope g(w) = M h(w), M = min(1, exp(alpha wc / v)) >= c(w) for w >= 0, closed-form integral
//      H(w) = (F1 - F2)/2 (1 - (1+w)^-2) + F2 (1 - (1+w)^-1) and inverse (quadratic in y = 1/(1+w)); acceptance c(w)/M;
//  (E) exponential envelope (wc < 0 only, slow projectiles): c(w) <= exp(-a (w - wc)) and h(w) <= h_sup, so
//      g(w) = h_sup exp(a wc) exp(-a w), a = alpha / v, inverse in closed form; acceptance
//      (h(w)/h_sup) / (1 + exp(-a (w - wc))).
// The orbital is chosen by w_i S_i (envelope mass); the accepted (orbital, w) has density proportional to the summed
// SDCS (exact for either envelope).
struct RuddShape {
    double F1, F2, wc, a;   // a = alpha / v
    double S, M, wmax, H;   // power envelope: H = H(wmax), mass S M H
    double hsup, ea, qe;    // exponential envelope: h_sup, exp(a wc), 1 - exp(-a wmax)
    double mass;            // envelope mass of the chosen envelope divided by S
    int mode;               // 0 power, 1 exponential
};
PFC_FN double rudd_h(const RuddShape& r, double w) {
    const double q = 1.0 + w;
    return (r.F1 + r.F2 * w) / (q * q * q);
}
PFC_FN RuddShape rudd_shape(const double* P, double B, double N, double Te, double a0, double R, double wmax) {
    RuddShape r;
    const double pi = 3.141592653589793;
    const double v2 = Te / B;
    const double v = pfc_sqrt(v2);
    const double lv = pfc::log(v);
    const double L1 = P[2] * pfc::exp(P[3] * lv) / (1.0 + P[4] * pfc::exp((P[3] + 4.0) * lv));
    const double H1 = P[0] * pfc::log1p(v2) / (v2 + P[1] / v2);
    const double L2 = P[7] * pfc::exp(P[8] * lv);
    const double H2 = P[5] / v2 + P[6] / (v2 * v2);
    r.F1 = L1 + H1;
    r.F2 = L2 * H2 / (L2 + H2);
    r.wc = 4.0 * v2 - 2.0 * v - R / (4.0 * B);
    r.a = P[9] / v;
    r.S = 4.0 * pi * a0 * a0 * N * (R / B) * (R / B);
    r.M = r.wc < 0.0 ? pfc::exp(r.a * r.wc) : 1.0;
    r.wmax = wmax;
    const double w = wmax, q = 1.0 + wmax;
    r.H = (r.F1 - r.F2) / 2.0 * (w * (w + 2.0) / (q * q)) + r.F2 * (w / q);
    // h(w) has its maximum at w = 0 when F2 <= 3 F1, else at w* = (F2 - 3 F1) / (2 F2)
    if (r.F2 <= 3.0 * r.F1) r.hsup = r.F1;
    else r.hsup = rudd_h(r, (r.F2 - 3.0 * r.F1) / (2.0 * r.F2));
    r.ea = r.M;
    r.qe = 1.0 - pfc::exp(-r.a * wmax);
    r.mode = 0;
    r.mass = r.M * r.H;
    if (r.wc < 0.0 && r.a * wmax > 1.0) {
        const double me = r.hsup * r.ea * r.qe / r.a;
        if (me < r.mass) {
            r.mode = 1;
            r.mass = me;
        }
    }
    return r;
}
// Envelope inverse: w with G(w) = u G(wmax) for the chosen envelope.
PFC_FN double rudd_propose(const RuddShape& r, double u) {
    double w;
    if (r.mode == 1) {
        w = -pfc::log(1.0 - u * r.qe) / r.a;
    } else {
        const double A = (r.F1 - r.F2) / 2.0;
        const double c0 = A + r.F2 - u * r.H;   // > 0
        const double y = 2.0 * c0 / (r.F2 + pfc_sqrt(r.F2 * r.F2 + 4.0 * A * c0));
        w = 1.0 / y - 1.0;
    }
    if (!(w >= 0.0)) w = 0.0;
    if (w > r.wmax) w = r.wmax;
    return w;
}
// Acceptance probability of a proposal w.
PFC_FN double rudd_accept(const RuddShape& r, double w) {
    if (r.mode == 1) return rudd_h(r, w) / r.hsup / (1.0 + pfc::exp(-r.a * (w - r.wc)));
    if (r.wc < 0.0) return 1.0 / (pfc::exp(r.a * r.wc) + pfc::exp(r.a * w));
    return 1.0 / (1.0 + pfc::exp(r.a * (w - r.wc)));
}
// Universal inner-shell parameters A1 B1 C1 D1 E1 A2 B2 C2 D2 alpha (Rudd, Kim, Madison & Gay, Rev. Mod. Phys. 64,
// 441 (1992), Table V column "Inner shells"; data/extracted/ion-rudd1992-table5-sdcs-parameters-inner-shells.json).
#define PFC_RUDD_INNER_SHELL_PARAMS 1.25, 0.50, 1.00, 1.00, 3.0, 1.10, 1.30, 1.00, 0.0, 0.66
#define PFC_RUDD_WORD0 32u      // first word of the rejection trials (3 words per trial)
#define PFC_RUDD_TRIALS 160u    // words 32..511

// ------------------------------------------------------------- chain decay
// Decay a virtual particle of rest mass Mv (eV) with momentum Pv in frame X
// into the final species list f[0..n-1] (masses m[]), with total kinetic
// energy Kkin in its rest frame (flat sequential partition for n > 2,
// isotropic directions). Momenta in X are written to pout[]. Words are taken
// from word index *wp onward (3 per split). Returns 0 on kinematic failure.
PFC_FN int chain_decay(const ChemKey& key, pfc_u32 ev, pfc_u32* wp, double Mv, V3 Pv, const double* m, pfc_u32 n,
                       double Kkin, V3* pout) {
    if (n == 1) {
        pout[0] = Pv;
        return 1;
    }
    double rest = 0.0;
    for (pfc_u32 i = 0; i < n; ++i) rest += m[i];
    double Krem = Kkin;
    V3 P = Pv;
    double M = Mv;
    for (pfc_u32 i = 0; i + 1 < n; ++i) {
        rest -= m[i];   // rest mass of the remainder after removing fragment i
        double Ks;
        if (i + 2 == n) Ks = Krem;
        else Ks = Krem * chem_u01(key, ev, (*wp)++);
        if (!(Ks >= 0.0)) return 0;
        const double Mr = rest + (Krem - Ks);   // virtual mass of the remainder
        const double ps = pstar(Ks, m[i], Mr);
        const double ua = chem_u01(key, ev, (*wp)++), ub = chem_u01(key, ev, (*wp)++);
        const V3 d = iso_dir(ua, ub);
        const V3 w = scale(P, 1.0 / M);
        const double gw = pfc_sqrt(1.0 + norm2(w));
        const TwoBody tb = two_body_tf(w, gw, m[i], Mr, ps, d);
        pout[i] = tb.p3;
        P = tb.p4;
        M = Mr;
        Krem = Krem - Ks;
    }
    pout[n - 1] = P;
    return 1;
}

// ---------------------------------------------------------------- one event
// Inputs in the target frame (TF): projectile momentum pp (eV), its kinetic
// energy T, projectile mass mp, target mass mt, target proper velocity wt
// (dimensionless) with gamma gt (for the final boost to the lab). On success
// writes lab momenta plab[] and kinematic masses mk[] (rest energy plus the
// internal energy I_i the product carries) of the n_out final particles in
// record order, and the lab-frame internal-sink energy sum_i gamma_i I_i
// (the internal energy is released in each product's rest frame, so its lab
// energy is gamma_i I_i and the product keeps its velocity); returns 0 when
// the event is kinematically closed.
// Internal-energy carriers: EXCITE -> target (I = loss); DISSOC -> last
// fragment (I = (1 - f) (loss - Q_chain)); IONIZE and ION_IONIZE -> ion, first
// chain product (I = (1 - f)(B_eff - q)); CAPTURE -> first product
// (I = (1 - f)(K_cm - q)); CHARGE_X with fragments -> last fragment
// (I = (1 - f) loss); BREAKUP -> last projectile fragment (I = (1 - f) loss).
PFC_FN int event_kinematics(const TableView& t, const ProcessRec& pr, const GroupRec& gr, int gi, const ChemKey& key,
                            pfc_u32 ev, double mp, double mt, V3 pp, double T, V3 wt, double gt, V3* plab, double* mk,
                            double* e_int) {
    pfc_u32 wp = PFC_W_KIN;
    V3 ptf[PFC_MAX_OUT];
    double I[PFC_MAX_OUT];
    for (pfc_u32 i = 0; i < pr.n_out; ++i) {
        mk[i] = t.species[pr.out_species[i]].mass_eV;
        I[i] = 0.0;
    }
    *e_int = 0.0;
    const double Kcm = kcm_from_T(mp, mt, T);
    const double sqrt_s = mp + mt + Kcm;
    const V3 wcm = scale(pp, 1.0 / sqrt_s);
    const double gcm = pfc_sqrt(1.0 + norm2(wcm));
    const V3 din = unit_or_z(pp);
    if (pr.kind == PFC_P_ELASTIC || pr.kind == PFC_P_EXCITE || pr.kind == PFC_P_DISSOC || pr.kind == PFC_P_CHARGE_X ||
        pr.kind == PFC_P_BREAKUP) {
        // Scatter family: (out0-or-virtual, out1-or-virtual) in CM.
        double m3, m4, Kout;
        if (pr.kind == PFC_P_CHARGE_X) {
            // out[0] fast product; out[1..] the target-side system (one ion, or a virtual ion of
            // energy loss_eV above its fragments that decays into out[1..])
            m3 = mk[0];
            if (pr.n_out > 2) {
                m4 = pr.loss_eV;
                for (pfc_u32 i = 1; i < pr.n_out; ++i) m4 += mk[i];
                Kout = Kcm - pr.q_eV - pr.loss_eV;
            } else {
                m4 = mk[1];
                Kout = Kcm - pr.q_eV;
            }
        } else if (pr.kind == PFC_P_BREAKUP) {
            // virtual projectile system (energy loss_eV above its fragments out[0..n-2]) and the target out[n-1]
            m3 = pr.loss_eV;
            for (pfc_u32 i = 0; i + 1 < pr.n_out; ++i) m3 += mk[i];
            m4 = mk[pr.n_out - 1];
            Kout = Kcm - pr.q_eV - pr.loss_eV;
        } else {
            m3 = mp;
            m4 = mt + pr.loss_eV;
            Kout = Kcm - pr.loss_eV;
        }
        if (!(Kout >= 0.0)) return 0;
        const double ua = chem_u01(key, ev, wp++), ub = chem_u01(key, ev, wp++);
        V3 d;
        if (pr.angular == PFC_ANG_ISO) d = iso_dir(ua, ub);
        else if (pr.angular == PFC_ANG_FORWARD) d = din;
        else {
            const double xi = xi_eval(t, gr, gi, T);
            d = deflect(din, screened_mu(xi, ua), ub);
        }
        const double ps = pstar(Kout, m3, m4);
        const TwoBody tb = two_body_tf(wcm, gcm, m3, m4, ps, d);
        ptf[0] = tb.p3;
        if (pr.kind == PFC_P_BREAKUP) {
            // virtual projectile -> out[0..n-2]; internal energy on fragment n-2
            const double kin = pr.loss_eV * pr.kinetic_fraction;
            const pfc_u32 nf = pr.n_out - 1;
            I[nf - 1] = pr.loss_eV - kin;
            mk[nf - 1] += I[nf - 1];
            double Mv = kin;
            for (pfc_u32 i = 0; i < nf; ++i) Mv += mk[i];
            if (!chain_decay(key, ev, &wp, Mv, tb.p3, mk, nf, kin, ptf)) return 0;
            ptf[nf] = tb.p4;
        } else if (pr.kind == PFC_P_CHARGE_X && pr.n_out > 2) {
            // virtual target-side ion -> out[1..]; internal energy on the last fragment
            const double kin = pr.loss_eV * pr.kinetic_fraction;
            I[pr.n_out - 1] = pr.loss_eV - kin;
            mk[pr.n_out - 1] += I[pr.n_out - 1];
            double Mv = kin;
            for (pfc_u32 i = 1; i < pr.n_out; ++i) Mv += mk[i];
            if (!chain_decay(key, ev, &wp, Mv, tb.p4, mk + 1, pr.n_out - 1, kin, ptf + 1)) return 0;
        } else if (pr.kind == PFC_P_DISSOC) {
            // virtual excited target -> fragments out[1..]; internal on the last fragment
            const double avail = pr.loss_eV - pr.chain_q_eV;
            if (!(avail >= 0.0)) return 0;
            const double kin = avail * pr.kinetic_fraction;
            const pfc_u32 nf = pr.n_out - 1;
            I[pr.n_out - 1] = avail - kin;
            mk[pr.n_out - 1] += I[pr.n_out - 1];
            double Mv = kin;
            for (pfc_u32 i = 1; i < pr.n_out; ++i) Mv += mk[i];
            if (!chain_decay(key, ev, &wp, Mv, tb.p4, mk + 1, nf, kin, ptf + 1)) return 0;
        } else {
            ptf[1] = tb.p4;
            if (pr.kind == PFC_P_EXCITE) {
                I[1] = pr.loss_eV;
                mk[1] = m4;
            }
        }
    } else if (pr.kind == PFC_P_CAPTURE) {
        const double avail = Kcm - pr.q_eV;
        if (!(avail >= 0.0)) return 0;
        const double kin = avail * pr.kinetic_fraction;
        I[0] = avail - kin;
        mk[0] += I[0];
        double Mv = kin;
        for (pfc_u32 i = 0; i < pr.n_out; ++i) Mv += mk[i];
        // compound at rest in CM; decay there, then boost CM -> TF
        V3 pc[PFC_MAX_OUT];
        if (!chain_decay(key, ev, &wp, Mv, v3(0.0, 0.0, 0.0), mk, pr.n_out, kin, pc)) return 0;
        for (pfc_u32 i = 0; i < pr.n_out; ++i) ptf[i] = boost_out(wcm, gcm, pc[i], energy_of(mk[i], pc[i]));
    } else if (pr.kind == PFC_P_IONIZE) {
        // out[0] scattered electron (survivor), out[1] ejected electron (birth),
        // out[2..] ion-system products (chain; the ion is out[2]).
        const double me = mp;
        // orbital selection by w_i sigma_i(T) with B_eff = max(B_i, q)
        double terms[PFC_MAX_ORB];
        double total = 0.0;
        const int rel = T > pr.relativistic_from_eV ? 1 : 0;
        for (pfc_u32 i = 0; i < pr.n_orb; ++i) {
            terms[i] = 0.0;
            const double Be = pr.orb_B[i] > pr.q_eV ? pr.orb_B[i] : pr.q_eV;
            if (pr.orb_w[i] > 0.0 && T > Be) {
                const SdcsShape s = sdcs_shape(Be, pr.orb_U[i], pr.orb_N[i], T, pr.bohr_m, pr.rydberg_eV, me, rel);
                const double v = pr.orb_w[i] * sdcs_total(s);
                terms[i] = v > 0.0 ? v : 0.0;
            }
            total += terms[i];
        }
        if (!(total > 0.0)) return 0;
        const double uo = chem_u01(key, ev, wp++) * total;
        pfc_u32 io = 0;
        double cum = 0.0;
        for (pfc_u32 i = 0; i < pr.n_orb; ++i) {
            if (!(terms[i] > 0.0)) continue;
            io = i;
            cum += terms[i];
            if (uo < cum) break;
        }
        const double Be = pr.orb_B[io] > pr.q_eV ? pr.orb_B[io] : pr.q_eV;
        const SdcsShape s = sdcs_shape(Be, pr.orb_U[io], pr.orb_N[io], T, pr.bohr_m, pr.rydberg_eV, me, rel);
        const double W = Be * sdcs_sample_w(s, chem_u01(key, ev, wp++));
        // binary-encounter ejection angle about the incident direction
        double c2 = W * (T + 2.0 * me) / (T * (W + 2.0 * me));
        if (c2 > 1.0) c2 = 1.0;
        const double mu = pfc_sqrt(c2);
        const V3 dw = deflect(din, mu, chem_u01(key, ev, wp++));
        const V3 pw = scale(dw, pfc_sqrt(W * (W + 2.0 * me)));
        const double Ex = Be - pr.q_eV;               // ion-system excitation above ground products
        const double Sigma = mt + Be;                // m_scattered + M_x (M_x = mt - me + q + Ex)
        const double D = T - W - Be;
        const V3 pr_ = sub(pp, pw);
        const double x = D * (D + 2.0 * Sigma) - norm2(pr_);
        if (!(x >= 0.0)) return 0;
        const double Kr = x / (pfc_sqrt(Sigma * Sigma + x) + Sigma);
        const double sr = Sigma + Kr;
        const V3 wr = scale(pr_, 1.0 / sr);
        const double gr_ = pfc_sqrt(1.0 + norm2(wr));
        const double Mx = Sigma - me;                // ion-system mass (with Ex)
        const double ps = pstar(Kr, me, Mx);
        const TwoBody tb = two_body_tf(wr, gr_, me, Mx, ps, unit_or_z(pr_));
        ptf[0] = tb.p3;
        ptf[1] = pw;
        // ion system: products out[2..]; internal energy on the ion
        const double kin = Ex * pr.kinetic_fraction;
        I[2] = Ex - kin;
        mk[2] += I[2];
        double Mv = kin;
        for (pfc_u32 i = 2; i < pr.n_out; ++i) Mv += mk[i];
        if (!chain_decay(key, ev, &wp, Mv, tb.p4, mk + 2, pr.n_out - 2, kin, ptf + 2)) return 0;
    } else if (pr.kind == PFC_P_ION_IONIZE) {
        // out[0] scattered projectile (survivor), out[1] ejected electron (birth), out[2..] ion-system
        // products (chain; the ion is out[2]). Orbital and ejected energy W from the Rudd SDCS with
        // B_eff = max(B_i, q) for the energy balance, W <= K_cm - B_eff; isotropic ejection direction
        // in the target frame (low-energy secondaries dominate and are broadly distributed).
        // Shells with B > 2 B_min (B_min the least tightly bound listed orbital) use the universal inner-shell
        // parameters of Rudd et al. 1992 Table V; the others use the target set pr.rudd.
        const double me = mk[1];
        const double b2 = T * (T + 2.0 * mp) / ((T + mp) * (T + mp));
        const double Te = 0.5 * me * b2;
        const double inner[10] = {PFC_RUDD_INNER_SHELL_PARAMS};
        double Bmin = 1e300;
        for (pfc_u32 i = 0; i < pr.n_orb; ++i) Bmin = pr.orb_B[i] < Bmin ? pr.orb_B[i] : Bmin;
        RuddShape sh[PFC_MAX_ORB];
        double env[PFC_MAX_ORB];
        double total = 0.0;
        for (pfc_u32 i = 0; i < pr.n_orb; ++i) {
            env[i] = 0.0;
            const double Be = pr.orb_B[i] > pr.q_eV ? pr.orb_B[i] : pr.q_eV;
            const double Wm = Kcm - Be;
            if (pr.orb_w[i] > 0.0 && Wm > 0.0 && Te > 0.0) {
                const double* P = pr.orb_B[i] > 2.0 * Bmin ? inner : pr.rudd;
                sh[i] = rudd_shape(P, pr.orb_B[i], pr.orb_N[i], Te, pr.bohr_m, pr.rydberg_eV, Wm / pr.orb_B[i]);
                const double e = pr.orb_w[i] * sh[i].S * sh[i].mass;
                env[i] = e > 0.0 ? e : 0.0;
            }
            total += env[i];
        }
        if (!(total > 0.0)) return 0;
        int io = -1;
        double wred = 0.0;
        for (pfc_u32 k = 0; k < PFC_RUDD_TRIALS; ++k) {
            const pfc_u32 w0 = PFC_RUDD_WORD0 + 3u * k;
            const double uo = chem_u01(key, ev, w0) * total;
            pfc_u32 i_ = 0;
            double cum = 0.0;
            for (pfc_u32 i = 0; i < pr.n_orb; ++i) {
                if (!(env[i] > 0.0)) continue;
                i_ = i;
                cum += env[i];
                if (uo < cum) break;
            }
            const double w = rudd_propose(sh[i_], chem_u01(key, ev, w0 + 1u));
            if (chem_u01(key, ev, w0 + 2u) < rudd_accept(sh[i_], w)) {
                io = (int)i_;
                wred = w;
                break;
            }
        }
        if (io < 0) return 0;
        const double Be = pr.orb_B[io] > pr.q_eV ? pr.orb_B[io] : pr.q_eV;
        double W = pr.orb_B[io] * wred;
        if (W > Kcm - Be) W = Kcm - Be;
        const V3 dw = iso_dir(chem_u01(key, ev, wp), chem_u01(key, ev, wp + 1u));
        wp += 2u;
        const V3 pw = scale(dw, pfc_sqrt(W * (W + 2.0 * me)));
        const double Ex = Be - pr.q_eV;               // ion-system excitation above ground products
        const double Sigma = mp + mt - me + Be;      // m_scattered + M_x (M_x = mt - me + q + Ex)
        const double D = T - W - Be;
        const V3 pr_ = sub(pp, pw);
        const double x = D * (D + 2.0 * Sigma) - norm2(pr_);
        if (!(x >= 0.0)) return 0;
        const double Kr = x / (pfc_sqrt(Sigma * Sigma + x) + Sigma);
        const double sr = Sigma + Kr;
        const V3 wr = scale(pr_, 1.0 / sr);
        const double gr_ = pfc_sqrt(1.0 + norm2(wr));
        const double Mx = Sigma - mp;                // ion-system mass (with Ex)
        const double ps = pstar(Kr, mp, Mx);
        const TwoBody tb = two_body_tf(wr, gr_, mp, Mx, ps, unit_or_z(pr_));
        ptf[0] = tb.p3;
        ptf[1] = pw;
        const double kin = Ex * pr.kinetic_fraction;
        I[2] = Ex - kin;
        mk[2] += I[2];
        double Mv = kin;
        for (pfc_u32 i = 2; i < pr.n_out; ++i) Mv += mk[i];
        if (!chain_decay(key, ev, &wp, Mv, tb.p4, mk + 2, pr.n_out - 2, kin, ptf + 2)) return 0;
    } else {
        return 0;
    }
    double eint = 0.0;
    for (pfc_u32 i = 0; i < pr.n_out; ++i) {
        plab[i] = boost_out(wt, gt, ptf[i], energy_of(mk[i], ptf[i]));
        if (I[i] > 0.0) {
            const double g = pfc_sqrt(1.0 + norm2(plab[i]) / (mk[i] * mk[i]));
            eint += g * I[i];
        }
    }
    *e_int = eint;
    return 1;
}

// ------------------------------------------------------------- coverage add
PFC_FN int cover_add(CarrierOut* o, pfc_u32 region, double x, double b) {
    for (pfc_u32 i = 0; i < o->n_cover; ++i)
        if (o->cover_region[i] == region) {
            o->cover_exposure[i] += x;
            o->cover_bound[i] += b;
            return 1;
        }
    if (o->n_cover >= PFC_MAX_COVER_SLOTS) return 0;
    o->cover_region[o->n_cover] = region;
    o->cover_exposure[o->n_cover] = x;
    o->cover_bound[o->n_cover] = b;
    o->n_cover++;
    return 1;
}

// Coverage bins of a group are stored as rows of the bound table: the group's
// bin edges are its grid points; bin i (grid interval i) has region id
// bin_region[group.bin_off + i] (0xFFFFFFFF when every process is resolved
// there) and bound row bin_row[...] (0xFFFFFFFF when no bound exists). The
// two regions "below the grid" and "above the grid" have ids bin_low/bin_high.
struct CoverView {
    const pfc_u32* bin_region;   // per group: grid_n + 1 entries (below, intervals 0..n-2, above)
    const pfc_u32* bin_row;
    const pfc_u32* bin_off;      // per group offset
};

PFC_FN pfc_u32 bin_index(int i, pfc_u32 n) {
    if (i < 0) return 0u;
    if (i >= (int)n - 1) return n;   // at/above the last point
    return (pfc_u32)i + 1u;
}

// Conservative region-wise support bound. Actual MCC samples targets inside a
// truncated thermal sphere, so a drift-only exposure is not a bound. For every
// reachable coverage bin, bound sigma*v by max(endpoint sigma)*v(E_hi). Combine
// bins belonging to the same region with max, not sum: one target cannot occupy
// several energies at once. Outside-grid regions must be unbounded in policy.
// This deliberately overestimates thermal tail exposure rather than silently
// treating missing cross sections as zero. No new collision approximation is
// introduced; only the refusal/diagnostic ledger changes.
PFC_FN int cover_support(const TableView& t, const CoverView& cv, pfc_u32 group,
                         double lo, double hi, double Wn_dt, CarrierOut* out) {
    const GroupRec& gr = t.groups[group];
    const double* grid = t.grid + gr.grid_off;
    const double m = t.species[gr.projectile].mass_eV;
    pfc_u32 ids[PFC_MAX_COVER_SLOTS];
    double speed[PFC_MAX_COVER_SLOTS], rate[PFC_MAX_COVER_SLOTS];
    pfc_u32 nr = 0;
    const int first = grid_locate(grid, gr.grid_n, lo);
    const int last = grid_locate(grid, gr.grid_n, hi);
    for (int i = first; i <= last; ++i) {
        const pfc_u32 bin = bin_index(i, gr.grid_n);
        const pfc_u32 region = cv.bin_region[cv.bin_off[group] + bin];
        if (region == 0xFFFFFFFFu) continue;
        double a = lo, b = hi;
        if (i >= 0 && a < grid[i]) a = grid[i];
        if (i < (int)gr.grid_n - 1 && b > grid[i + 1]) b = grid[i + 1];
        if (b < a) continue;
        const double vmax = pfc_sqrt(b * (b + 2.0 * m)) / (m + b) * PFC_C_M_S;
        double sv = 0.0;
        const pfc_u32 row = cv.bin_row[cv.bin_off[group] + bin];
        if (row != 0xFFFFFFFFu && i >= 0 && i < (int)gr.grid_n - 1) {
            const double sa = bound_eval(t, gr, row, i, a);
            const double sb = bound_eval(t, gr, row, i, b);
            sv = (sa > sb ? sa : sb) * vmax * (1.0 + 1e-12);
        }
        pfc_u32 j = 0;
        while (j < nr && ids[j] != region) ++j;
        if (j == nr) {
            if (nr == PFC_MAX_COVER_SLOTS) return 0;
            ids[nr] = region; speed[nr] = 0.0; rate[nr] = 0.0; ++nr;
        }
        if (vmax > speed[j]) speed[j] = vmax;
        if (sv > rate[j]) rate[j] = sv;
    }
    for (pfc_u32 j = 0; j < nr; ++j)
        if (!cover_add(out, ids[j], Wn_dt * speed[j] * (1.0 + 1e-12), Wn_dt * rate[j])) return 0;
    return 1;
}

// --------------------------------------------------------------- the step
PFC_FN void mcc_carrier(const TableView& t, const CoverView& cv, const StepConfig& cfg, const CarrierIn& cin,
                        const CellTargets& cell, EventRec* evs, CarrierOut* out) {
    const double c = PFC_C_M_S;
    out->id = cin.id;
    out->status = PFC_ST_OK;
    out->n_events = 0;
    out->n_candidates = 0;
    out->species_final = cin.species;
    out->u_final[0] = cin.u[0];
    out->u_final[1] = cin.u[1];
    out->u_final[2] = cin.u[2];
    out->nu_dt_max = 0.0;
    out->n_cover = 0;
    out->pad = 0;
    for (int i = 0; i < PFC_MAX_COVER_SLOTS; ++i) {
        out->cover_region[i] = 0xFFFFFFFFu;
        out->cover_exposure[i] = 0.0;
        out->cover_bound[i] = 0.0;
    }
    ChemKey key;
    key.seed = cfg.seed;
    key.step = cfg.step;
    key.carrier = cin.id;
    key.purpose = cfg.purpose;
    if (!chem_key_valid(key, 0u) || cin.species >= t.n_species) {
        out->status = cin.species >= t.n_species ? PFC_ST_BAD_INPUT : PFC_ST_KEY_DOMAIN;
        return;
    }
    const double m = t.species[cin.species].mass_eV;
    V3 u = v3(cin.u[0], cin.u[1], cin.u[2]);
    if (!(norm2(u) < 1e300) || !(cin.weight > 0.0) || !(cin.weight < 1e300)) {
        out->status = PFC_ST_BAD_INPUT;
        return;
    }
    const double W = cin.weight;
    if (!(cfg.dt_s >= 0.0 && cfg.dt_s < 1e300) ||
        !(cfg.trunc_sigmas >= 0.0 && cfg.trunc_sigmas < 1e100) || cell.n_targets > PFC_MAX_TARGETS) {
        out->status = PFC_ST_BAD_INPUT;
        return;
    }
    for (pfc_u32 k = 0; k < cell.n_targets; ++k) {
        if (cell.species[k] >= t.n_species || !(cell.n[k] >= 0.0 && cell.n[k] < 1e300) ||
            !(cell.s[k] >= 0.0 && cell.s[k] < 1e150) ||
            !(norm2(v3(cell.drift[k][0], cell.drift[k][1], cell.drift[k][2])) < 1e300)) {
            out->status = PFC_ST_BAD_INPUT;
            return;
        }
    }
    double tacc = 0.0;
    pfc_u32 e = 0;
    for (;;) {
        const V3 p = p_from_u(m, u);
        const double pn = pfc_sqrt(norm2(p));
        const double K = ke_from_u(m, u);
        // majorant per target slot
        double nu[PFC_MAX_TARGETS], M[PFC_MAX_TARGETS];
        double cover_lo[PFC_MAX_TARGETS], cover_hi[PFC_MAX_TARGETS];
        pfc_u32 gidx[PFC_MAX_TARGETS];
        double numax = 0.0;
        for (pfc_u32 k = 0; k < PFC_MAX_TARGETS; ++k) {
            nu[k] = 0.0;
            M[k] = 0.0;
            gidx[k] = 0xFFFFFFFFu;
            cover_lo[k] = cover_hi[k] = 0.0;
        }
        for (pfc_u32 k = 0; k < cell.n_targets && k < PFC_MAX_TARGETS; ++k) {
            const pfc_u32 ts = cell.species[k];
            const pfc_u32 g = group_index(t, cin.species, ts);
            if (g == 0xFFFFFFFFu || !(cell.n[k] > 0.0)) continue;
            gidx[k] = g;
            const GroupRec& gr = t.groups[g];
            const double vcap = pfc_sqrt(norm2(v3(cell.drift[k][0], cell.drift[k][1], cell.drift[k][2]))) + cfg.trunc_sigmas * cell.s[k];
            const double wc = vcap / c;
            const double gc = pfc_sqrt(1.0 + wc * wc);
            const double base = gc * K + gamma_minus_one(wc * wc) * m;
            const double Ehi = (base + wc * pn) * (1.0 + 1e-12);
            if (!(Ehi >= 0.0 && Ehi < 1e300)) {
                out->status = PFC_ST_BAD_INPUT;
                return;
            }
            // projectile at rest in some admissible target frame?  (beta_cap >= beta_p)
            const double wp = pn / m;   // projectile proper speed (dimensionless)
            double Elo = (wc >= wp) ? 0.0 : (base - wc * pn) * (1.0 - 1e-12);
            if (Elo < 0.0) Elo = 0.0;
            M[k] = majorant(t, gr, Elo, Ehi);
            cover_lo[k] = Elo;
            cover_hi[k] = Ehi;
            // A cold drifting target has one exact target-frame energy. The
            // origin-centred speed envelope above is safe but unnecessarily
            // broad for coverage when the thermal radius is zero.
            if (cell.s[k] == 0.0) {
                const V3 wd = w_from_u(v3(cell.drift[k][0], cell.drift[k][1], cell.drift[k][2]));
                const V3 pd = boost_in(wd, pfc_sqrt(1.0 + norm2(wd)), p, m + K);
                cover_lo[k] = cover_hi[k] = kinetic_from_p2(m, norm2(pd));
            }
            nu[k] = cell.n[k] * M[k];
            numax += nu[k];
        }
        const double remaining = cfg.dt_s - tacc;
        // waiting time to the next candidate
        double tau;
        if (numax > 0.0) {
            const double g = numax * cfg.dt_s;
            if (g > out->nu_dt_max) out->nu_dt_max = g;
            if (e >= cfg.max_candidates) {
                out->status = PFC_ST_CAND_CAP;
                break;
            }
            if (!chem_key_valid(key, e)) {
                out->status = PFC_ST_KEY_DOMAIN;
                break;
            }
            tau = -pfc::log(chem_u01(key, e, PFC_W_TIME)) / numax;
        } else {
            tau = remaining;
        }
        const double seg = tau < remaining ? tau : remaining;
        // Region-wise upper bounds over the same thermal support used by MCC.
        for (pfc_u32 k = 0; k < cell.n_targets && k < PFC_MAX_TARGETS; ++k) {
            const pfc_u32 g = group_index(t, cin.species, cell.species[k]);
            if (g == 0xFFFFFFFFu || !(cell.n[k] > 0.0) || !(seg > 0.0)) continue;
            if (!cover_support(t, cv, g, cover_lo[k], cover_hi[k], W * cell.n[k] * seg, out)) {
                out->status = PFC_ST_COVER_CAP;
            }
        }
        if (out->status != PFC_ST_OK) break;
        if (!(numax > 0.0) || !(tau < remaining)) break;   // no further candidate in this step
        tacc += tau;
        // ---- candidate e
        double xs = chem_u01(key, e, PFC_W_TARGET) * numax;
        pfc_u32 k = 0xFFFFFFFFu, last = 0xFFFFFFFFu;
        double cum = 0.0;
        for (pfc_u32 j = 0; j < cell.n_targets && j < PFC_MAX_TARGETS; ++j) {
            if (!(nu[j] > 0.0)) continue;
            last = j;
            cum += nu[j];
            if (xs < cum) {
                k = j;
                break;
            }
        }
        if (k == 0xFFFFFFFFu) k = last;   // rounding at the top of the cumulative sum
        const GroupRec& gr = t.groups[gidx[k]];
        const pfc_u32 ts = cell.species[k];
        const double mt = t.species[ts].mass_eV;
        // truncated Maxwellian target velocity (up to 4 attempts, then the drift itself)
        V3 ut = v3(cell.drift[k][0], cell.drift[k][1], cell.drift[k][2]);
        {
            const double r2max = cfg.trunc_sigmas * cfg.trunc_sigmas;
            for (pfc_u32 a = 0; a < 4u; ++a) {
                const pfc_u32 w0 = PFC_W_TVEL + 4u * a;
                const double r1 = pfc_sqrt(-2.0 * pfc::log(chem_u01(key, e, w0)));
                const double r2 = pfc_sqrt(-2.0 * pfc::log(chem_u01(key, e, w0 + 2u)));
                double c1, s1, c2, s2;
                sincos_2pi(chem_u01(key, e, w0 + 1u), &c1, &s1);
                sincos_2pi(chem_u01(key, e, w0 + 3u), &c2, &s2);
                const double z0 = r1 * c1, z1 = r1 * s1, z2 = r2 * c2;
                if (z0 * z0 + z1 * z1 + z2 * z2 <= r2max) {
                    ut = add(ut, scale(v3(z0, z1, z2), cell.s[k]));
                    break;
                }
            }
        }
        const V3 wt = w_from_u(ut);
        const double gt = pfc_sqrt(1.0 + norm2(wt));
        const V3 pt = boost_in(wt, gt, p, m + K);
        const double T = kinetic_from_p2(m, norm2(pt));
        const double vt = pfc_sqrt(norm2(pt)) / (m + T) * c;
        const int gi = grid_locate(t.grid + gr.grid_off, gr.grid_n, T);
        const double y = chem_u01(key, e, PFC_W_CHANNEL) * M[k];
        pfc_u32 sel = 0xFFFFFFFFu;
        double F = 0.0;
        // below the grid: constant-rate continuation sigma v = sigma(E_0) v(E_0) for flagged processes
        double v0 = 0.0;
        if (gi < 0 && gr.low_rate) {
            const double E0 = t.grid[gr.grid_off];
            v0 = pfc_sqrt(E0 * (E0 + 2.0 * m)) / (E0 + m) * c;
        }
        for (pfc_u32 j = 0; j < gr.proc_n; ++j) {
            double sv;
            if (gi >= 0) sv = sigma_eval(t, gr, j, gi, T) * vt;
            else if (t.procs[gr.proc_off + j].flags & PFC_PFLAG_CONST_RATE_BELOW) sv = t.sigma[gr.sigma_off + (pfc_u64)j * gr.grid_n] * v0;
            else sv = 0.0;
            F += sv;
            if (sel == 0xFFFFFFFFu && y < F) sel = j;
        }
        if (F > M[k]) {
            out->status = PFC_ST_MAJORANT;
            break;
        }
        if (sel == 0xFFFFFFFFu) {   // null candidate
            ++e;
            continue;
        }
        if (out->n_events >= cfg.max_events) {
            out->status = PFC_ST_EVENT_CAP;
            break;
        }
        const pfc_u32 pidx = gr.proc_off + sel;
        const ProcessRec& pr = t.procs[pidx];
        V3 plab[PFC_MAX_OUT];
        double mk[PFC_MAX_OUT];
        double eint = 0.0;
        const int ok = event_kinematics(t, pr, gr, gi, key, e, m, mt, pt, T, wt, gt, plab, mk, &eint);
        if (!ok) {   // kinematically closed at this energy: counted as a null candidate
            ++e;
            continue;
        }
        EventRec& r = evs[out->n_events];
        r.carrier = cin.id;
        r.candidate = e;
        r.process = pidx;
        r.cell = cin.cell;
        r.target_species = ts;
        r.n_out = pr.n_out;
        {
            pfc_u32 reg;
            r.cover_class = coverage_class(pr, T, &reg);
        }
        r.internal_sink = pr.internal_sink;
        r.flags = 0;
        r.weight = W;
        r.time_s = tacc;
        r.u_in[0] = u.x;
        r.u_in[1] = u.y;
        r.u_in[2] = u.z;
        r.ut[0] = ut.x;
        r.ut[1] = ut.y;
        r.ut[2] = ut.z;
        r.E_tf = T;
        r.ke_in_p = K;
        r.ke_in_t = ke_from_u(mt, ut);
        r.q_eV = pr.q_eV;
        r.e_internal = eint;
        int survivor = -1;
        for (pfc_u32 i = 0; i < PFC_MAX_OUT; ++i) {
            if (i < pr.n_out) {
                const double mo = t.species[pr.out_species[i]].mass_eV;
                const V3 uo = u_from_p(mk[i], plab[i]);   // velocity of the (possibly excited) product
                r.out_species[i] = pr.out_species[i];
                r.out_dest[i] = pr.out_dest[i];
                r.out_u[i][0] = uo.x;
                r.out_u[i][1] = uo.y;
                r.out_u[i][2] = uo.z;
                r.ke_out[i] = ke_from_u(mo, uo);
                if (pr.out_dest[i] == PFC_DEST_SURVIVOR) survivor = (int)i;
            } else {
                r.out_species[i] = 0xFFFFFFFFu;
                r.out_dest[i] = 0u;
                r.out_u[i][0] = r.out_u[i][1] = r.out_u[i][2] = 0.0;
                r.ke_out[i] = 0.0;
            }
        }
        out->n_events++;
        ++e;
        if (survivor < 0) {
            out->status = PFC_ST_CONSUMED;
            out->species_final = 0xFFFFFFFFu;
            break;
        }
        u = v3(r.out_u[survivor][0], r.out_u[survivor][1], r.out_u[survivor][2]);
    }
    out->n_candidates = e;
    out->u_final[0] = u.x;
    out->u_final[1] = u.y;
    out->u_final[2] = u.z;
}

}  // namespace pfc
