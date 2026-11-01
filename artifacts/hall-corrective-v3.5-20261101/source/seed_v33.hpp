#pragma once
// Hall corrective v3.3 DC seed: axial::dc_initial of Hall 8dddec0 (hall-finite-foil-em-20261029/axial_fields.hpp),
// statement for statement, with two changes in the magnetostatic part:
//  1. the Coulomb vector-potential solve keeps the sealed Jacobi CG (foil::solve arithmetic) but its independent
//     residual ||b - K x|| is evaluated in double-word arithmetic (Words, field_contract.hpp) and the correction
//     cycles use that residual (iterative refinement); target unchanged max(1e-30, 2e-13 ||b||); up to 8 cycles.
//     (A double-precision residual evaluation carries a rounding error of 2.3e-13 ||b|| on the 25 um mesh, above
//     the target, so the sealed check cannot pass there; resolution-series F1.)
//  2. B0 = curl(A) as double-word face values (exact products, compensated four-term sums): the seed's discrete
//     divergence is O(eps^2) instead of O(eps) (Hall v3.2 pre-run divB analysis: the refusing cell carried the
//     seed's rounding divergence).
// Everything else (scalar potentials, E, currents, background force from the rounded B) is the sealed code.
#include "field_contract.hpp"

namespace hall_corrective_v3 {
inline constexpr char seed_rule_v33[] =
    "dc-seed-v33: 8dddec0 dc_initial; Coulomb vector solve with the sealed Jacobi CG and a double-word independent residual "
    "(refinement cycles <= 8, target max(1e-30, 2e-13 ||b||)); B0 = double-word curl(A) (B, B_low); background force from "
    "the rounded B0";
struct SeedV33 {
    axial::Initial base;
    Vec B_low;                        // low words of B0 (base.B = high words)
    int coulomb_cycles = 0;
    double coulomb_residual = 0, coulomb_target = 0;
};
namespace seed_detail {
// Words curl of a double edge vector: per face sum_k Words(A[index_k]) * c_k in coefficient order.
inline WordVec curl_words(const axial::Mesh& g, const Vec& A) {
    WordVec B(std::size_t(g.nb));
    for (int b = 0; b < g.nb; ++b) { Words s; for (int k = 0; k < 4; ++k) s = s + Words(A[g.C[b].index[k]]) * g.C[b].coefficient[k]; B[b] = s; }
    return B;
}
// Double-word evaluation of the Coulomb-gauge vector-potential operator of dc_initial (same double coefficients).
inline WordVec coulomb_apply_words(const axial::Mesh& g, Vec A) {
    for (int e = 0; e < g.ne; ++e) if (g.fixed_A[e]) A[e] = 0;
    auto B = curl_words(g, A);
    for (int b = 0; b < g.nb; ++b) B[b] = B[b] * g.H[b];
    WordVec K(std::size_t(g.ne));
    for (int b = 0; b < g.nb; ++b) for (int k = 0; k < 4; ++k) { const int e = g.C[b].index[k]; K[e] = K[e] + B[b] * g.C[b].coefficient[k]; }
    WordVec p(std::size_t(g.nn));
    for (int e = 0; e < g.ne; ++e) { const auto q = (Words(A[e]) * g.VE[e]) / g.edges[e].length; p[g.edges[e].a] = p[g.edges[e].a] - q; p[g.edges[e].b] = p[g.edges[e].b] + q; }
    for (int n = 0; n < g.nn; ++n) p[n] = g.outer_node[n] ? Words() : p[n] / (axial::mu * g.VN[n]);
    for (int e = 0; e < g.ne; ++e) {
        if (g.fixed_A[e]) { K[e] = Words(); continue; }
        const auto grad = (p[g.edges[e].b] - p[g.edges[e].a]) / g.edges[e].length;
        K[e] = K[e] + grad * g.VE[e];
    }
    return K;
}
inline double norm_words(const WordVec& r) { long double s = 0; for (const auto& x : r) { const long double v = x.wide(); s += v * v; } return double(std::sqrt(s)); }
}  // namespace seed_detail

inline SeedV33 dc_initial_v33(const axial::Mesh& g, const Vec& sigma_cells, double current) {
    using namespace axial;
    hall::require(std::isfinite(current), "invalid supported DC start");
    SeedV33 seed; Initial& out = seed.base;
    out.current_A = current; out.feed = g.feed_shape(); auto S = g.conduction(sigma_cells);
    Vec impressed = out.feed; for (double& v : impressed) v *= current;
    auto rhs = g.gradient_transpose(impressed); std::vector<bool> fixed(g.nn);
    for (int n = 0; n < g.nn; ++n) fixed[n] = g.node_body[n] < 0 || g.short_node[n];
    double uniform_sigma = 0; bool uniform = true;
    for (std::size_t c = 0; c < g.cells.size(); ++c) if (g.cells[c].body >= 0) {
        hall::require(sigma_cells[c] > 0, "initial DC metal must conduct");
        if (uniform_sigma == 0) uniform_sigma = sigma_cells[c]; else if (sigma_cells[c] != uniform_sigma) uniform = false; }
    Vec phi(g.nn);
    if (uniform) {
        double gradient = current / (.005 * .0002 * uniform_sigma);
        for (int k = g.kshort; k <= g.kfeed; ++k) for (int j = 0; j <= g.ny; ++j) for (int i = 0; i <= g.nx; ++i) {
            int n = g.node(i, j, k); if (g.node_body[n] >= 0) phi[n] = (g.node_body[n] ? -1 : 1) * gradient * (k - g.kshort) * g.dz; }
        auto applied = g.gradient(phi); for (int e = 0; e < g.ne; ++e) applied[e] *= S[e]; auto residual = g.gradient_transpose(applied); Vec b = rhs;
        for (int n = 0; n < g.nn; ++n) { residual[n] = fixed[n] ? 0 : residual[n] - rhs[n]; if (fixed[n]) b[n] = 0; }
        hall::require(norm(residual) <= std::max(1e-30, 2e-13 * norm(b)), "analytic DC conductor independent residual failed");
    } else try { phi = scalar_dirichlet(g, S, rhs, fixed); } catch (const std::exception& e) { throw std::runtime_error(std::string("DC conductor scalar: ") + e.what()); }
    for (int n = 0; n < g.nn; ++n) fixed[n] = g.node_body[n] >= 0 || g.short_node[n];
    auto known = g.gradient(phi); for (int e = 0; e < g.ne; ++e) known[e] *= g.VE[e];
    rhs = g.gradient_transpose(known); for (double& v : rhs) v = -v;
    Vec vacuum; try { vacuum = scalar_dirichlet(g, g.VE, rhs, fixed); } catch (const std::exception& e) { throw std::runtime_error(std::string("DC vacuum scalar: ") + e.what()); }
    for (int n = 0; n < g.nn; ++n) phi[n] += vacuum[n]; out.potential = phi;
    out.E = g.gradient(phi); for (double& v : out.E) v = -v;
    out.terminal_V = -dot(out.E, out.feed);
    out.solid_current.resize(g.ne); Vec closed(g.ne);
    for (int e = 0; e < g.ne; ++e) closed[e] = (out.solid_current[e] = S[e] * out.E[e]) + impressed[e];
    auto residual = g.gradient_transpose(closed); out.short_current.resize(g.ne);
    if (uniform) {
        for (int j = 0; j < g.ny; ++j) for (int i = 0; i <= g.nx; ++i) out.short_current[g.ey(i, j, g.kshort)] = -impressed[g.ey(i, j, g.kfeed)];
        auto div = g.gradient_transpose(out.short_current); Vec b(g.nn), error(g.nn);
        for (int n = 0; n < g.nn; ++n) if (g.short_node[n]) { b[n] = residual[n]; error[n] = residual[n] + div[n]; }
        hall::require(norm(error) <= std::max(1e-30, 2e-13 * norm(b)), "analytic DC short independent residual failed");
    } else {
        Vec short_weight(g.ne); for (int e = 0; e < g.ne; ++e) if (g.fixed_E[e]) short_weight[e] = g.VE[e];
        bool anchored = false; for (int n = 0; n < g.nn; ++n) { fixed[n] = !g.short_node[n];
            if (!anchored && g.short_node[n]) { fixed[n] = true; anchored = true; } rhs[n] = -residual[n]; }
        Vec short_phi; try { short_phi = scalar_dirichlet(g, short_weight, rhs, fixed); } catch (const std::exception& e) { throw std::runtime_error(std::string("DC short scalar: ") + e.what()); }
        out.short_current = g.gradient(short_phi); for (int e = 0; e < g.ne; ++e) out.short_current[e] *= short_weight[e];
    }
    for (int e = 0; e < g.ne; ++e) closed[e] += out.short_current[e];
    out.dc_continuity_A = maximum(g.gradient_transpose(closed));
    hall::require(out.dc_continuity_A <= 2e-10 * std::max(1., std::abs(current)), "initial feed/solid/short KCL");
    // Magnetostatic vector potential (sealed operator and Jacobi CG; v3.3: double-word independent residual).
    Vec diagonal(g.ne);
    for (int b = 0; b < g.nb; ++b) for (int q = 0; q < 4; ++q) { int e = g.C[b].index[q]; diagonal[e] += g.H[b] * std::pow(g.C[b].coefficient[q], 2); }
    for (int e = 0; e < g.ne; ++e) { const auto& v = g.edges[e];
        diagonal[e] += g.VE[e] * g.VE[e] / (mu * v.length * v.length) * ((g.outer_node[v.a] ? 0 : 1 / g.VN[v.a]) + (g.outer_node[v.b] ? 0 : 1 / g.VN[v.b]));
        if (g.fixed_A[e]) { diagonal[e] = 1; closed[e] = 0; } }
    auto apply = [&](Vec A) { for (int e = 0; e < g.ne; ++e) if (g.fixed_A[e]) A[e] = 0;
        auto K = g.curl_adjoint(g.curl(A)); for (int e = 0; e < g.ne; ++e) A[e] *= g.VE[e];
        auto div = g.gradient_transpose(A); for (int n = 0; n < g.nn; ++n) div[n] = g.outer_node[n] ? 0 : div[n] / (mu * g.VN[n]); auto grad = g.gradient(div);
        for (int e = 0; e < g.ne; ++e) K[e] = g.fixed_A[e] ? 0 : K[e] + g.VE[e] * grad[e]; return K; };
    Vec A(closed.size()), r = closed;
    const double target = std::max(1e-30, 2e-13 * norm(closed));
    seed.coulomb_target = target;
    bool converged = false;
    for (int cycle = 0; cycle < 8 && !converged; ++cycle) {   // foil::solve iteration, refinement with the double-word residual
        Vec z(r.size()); for (std::size_t i = 0; i < r.size(); ++i) z[i] = r[i] / diagonal[i]; Vec p = z; double rz = dot(r, z);
        for (int iteration = 0; iteration < 5000 && norm(r) > target * .2; ++iteration) {
            auto Ap = apply(p); double denominator = dot(p, Ap);
            if (!(denominator > 0 && std::isfinite(denominator))) throw std::runtime_error("DC Coulomb vector: non-SPD field operator");
            double alpha = rz / denominator;
            for (std::size_t i = 0; i < r.size(); ++i) { A[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; z[i] = r[i] / diagonal[i]; }
            double nrz = dot(r, z), beta = nrz / rz; for (std::size_t i = 0; i < r.size(); ++i) p[i] = z[i] + beta * p[i]; rz = nrz; }
        const auto K = seed_detail::coulomb_apply_words(g, A);
        WordVec rw(K.size()); for (int e = 0; e < g.ne; ++e) rw[e] = Words(closed[e]) - K[e];
        seed.coulomb_residual = seed_detail::norm_words(rw); seed.coulomb_cycles = cycle + 1;
        for (int e = 0; e < g.ne; ++e) r[e] = rw[e].rounded();
        converged = seed.coulomb_residual <= target;
    }
    if (!converged || !finite(A)) throw std::runtime_error("DC Coulomb vector: field solve independent residual failed (double-word residual)");
    const auto Bw = seed_detail::curl_words(g, A);
    out.B.resize(g.nb); seed.B_low.resize(g.nb);
    for (int b = 0; b < g.nb; ++b) { out.B[b] = Bw[b].hi; seed.B_low[b] = Bw[b].lo; }
    Vec Br(g.nb); for (int b = 0; b < g.nb; ++b) Br[b] = Bw[b].rounded();
    auto ampere = g.curl_adjoint(Br); out.background_force.resize(g.ne);
    for (int e = 0; e < g.ne; ++e) {
        double defect = out.solid_current[e] + impressed[e] + out.short_current[e] - ampere[e];
        if (!g.fixed_A[e]) out.magnetic_residual_A_m = std::max(out.magnetic_residual_A_m, std::abs(defect));
        if (g.fixed_A[e]) out.background_force[e] = g.Y[e] * out.E[e] + defect;
    }
    hall::require(out.magnetic_residual_A_m <= 2e-9 * std::max(1., maximum(closed)), "initial interior Ampere residual");
    return seed;
}
}  // namespace hall_corrective_v3
