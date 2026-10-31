#pragma once
// motion-work-gross-v1: tolerance of the per-step motion work identity |K - W| (Hall v3.2; specification
// MOTION-WORK-GROSS-V1.txt in the runtime directory, an amendment to native v7 API-007 1.1 / HALL-LINK 2c).
// The scales are exact sums in the native v7 superaccumulator form (pf_v7::ExactSum, products exact, one
// rounding), so they do not depend on order or partition and are bit-identical on every producer and on
// native; E^ is an exact maximum; the tolerance is evaluated in the fixed order of tolerance().
// Only the native v7 forms library (forms/include, read in place) is used; no Hall type appears here, so the
// CPU motion assembly and the device motion fold share this code.
#include "pf_v7/forms.hpp"
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

namespace motion_work_v1 {
inline constexpr char rule[] =
    "motion-work-gross-v1: |RN(K - W)| <= tol, tol = RN(RN(RN(1e-20 + RN(1e-10 * RN(G_K + G_W))) + RN(2^-47 * C_K)) + "
    "RN(RN(8e-18 * E^) * Q)); over the charged census carriers i (q_i != 0): G_K = RN(sum w_i |d_i|), C_K = RN(sum w_i K0_i + "
    "w_i RN(K0_i + d_i)), Q = RN(sum |q_i| w_i), d_i = RN(KE(m_i, u_end_i) - KE(m_i, u_start_i)), K0_i = KE(m_i, u_start_i); "
    "G_W = RN(dt * RN(sum_e |E_e| |J_e|)); E^ = max |E_e| over the edges e with an endpoint node on an edge carrying "
    "J != 0 (0 when none); exact sums (native v7 superaccumulator), exact maximum";

// Scales and tolerance of one motion (field for field the Hall receipt record WorkScales).
struct Scales {
    double kinetic_gross_J = 0, kinetic_content_J = 0, field_work_gross_J = 0, charge_abs_C = 0, field_bound_V_m = 0;
    double tolerance_J = 0;        // tolerance(*this)
    double v31_tolerance_J = 0;    // diagnostic: 1e-20 + 1e-10 (|K| + |W|) (the v3.1 / v6 net form)
    std::uint32_t status = 0;      // pf7 sum status bits of the four sums; E^ nonfinite -> pf7::sum_nonfinite_input
};

// Carrier side: one call per census carrier (any order, any partition: exact sums).
struct CarrierSums {
    pf_v7::ExactSum gk, ck, q;
    // weight w > 0, charge q, K0 = KE(m, u_start), d = RN(KE(m, u_end) - K0) (the kinetic change term of K)
    void add(double weight, double charge, double K0, double d) {
        if (charge == 0) return;   // neutral carriers exchange no work with the field
        const double K1 = K0 + d;  // RN(K0 + d)
        gk.add_product(weight, std::fabs(d));
        ck.add_product(weight, K0);
        ck.add_product(weight, K1);
        q.add_product(std::fabs(charge), weight);
    }
    void merge(const CarrierSums& o) { gk.merge(o.gk); ck.merge(o.ck); q.merge(o.q); }
};

// tol = RN(RN(RN(1e-20 + RN(1e-10 * RN(G_K + G_W))) + RN(2^-47 * C_K)) + RN(RN(8e-18 * E^) * Q)), binary64, no contraction.
inline double tolerance(const Scales& s) {
    const double gross = s.kinetic_gross_J + s.field_work_gross_J;
    const double relative = 1e-10 * gross;
    double t = 1e-20 + relative;
    const double rounding = 0x1p-47 * s.kinetic_content_J;   // 32 eps C_K (exact scaling)
    t = t + rounding;
    const double length = 8e-18 * s.field_bound_V_m;
    const double fixed_point = length * s.charge_abs_C;
    return t + fixed_point;
}

// Field side: G_W and E^ from the trial-field edge values E, the plasma current J (same length) and the edge
// endpoint nodes (edge_a[e], edge_b[e] < nn).
template <class IndexA, class IndexB>
void field_scales(Scales& s, double dt, std::span<const double> E, std::span<const double> J, const IndexA& edge_a, const IndexB& edge_b,
                  std::size_t nn) {
    pf_v7::ExactSum gw;
    std::vector<char> mark(nn, 0);
    for (std::size_t e = 0; e < J.size(); ++e) {
        gw.add_product(std::fabs(E[e]), std::fabs(J[e]));
        if (J[e] != 0) { mark[std::size_t(edge_a[e])] = 1; mark[std::size_t(edge_b[e])] = 1; }
    }
    double em = 0;
    for (std::size_t e = 0; e < J.size(); ++e)
        if (mark[std::size_t(edge_a[e])] || mark[std::size_t(edge_b[e])]) {
            const double a = std::fabs(E[e]);
            if (!(a <= em)) em = a;   // NaN propagates (refused below)
        }
    std::uint32_t st = 0;
    const double sum = gw.value(&st);
    s.status |= st;
    s.field_work_gross_J = dt * sum;
    if (!std::isfinite(s.field_work_gross_J)) s.status |= pf7::sum_overflow;
    s.field_bound_V_m = em;
    if (!std::isfinite(em)) s.status |= pf7::sum_nonfinite_input;
}

// Carrier sums into the scales; then the tolerance and the diagnostic v3.1 form.
inline void finish(Scales& s, const CarrierSums& c, double K, double W) {
    std::uint32_t st = 0;
    s.kinetic_gross_J = c.gk.value(&st); s.status |= st; st = 0;
    s.kinetic_content_J = c.ck.value(&st); s.status |= st; st = 0;
    s.charge_abs_C = c.q.value(&st); s.status |= st;
    s.tolerance_J = tolerance(s);
    s.v31_tolerance_J = 1e-20 + 1e-10 * (std::abs(K) + std::abs(W));
}
// The guard: the defect RN(K - W) within the tolerance, every sum finite.
inline bool admits(const Scales& s, double defect) {
    return s.status == 0 && std::isfinite(defect) && std::isfinite(s.tolerance_J) && std::abs(defect) <= s.tolerance_J;
}
}  // namespace motion_work_v1
