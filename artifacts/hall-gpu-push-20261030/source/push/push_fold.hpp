// push_fold.hpp -- host part of the motion record: refusal order of
// prepare_canonical_motion, the x87 long double folds (kinetic change in carrier-ID
// order, foil::dot for the field work), continuity with Mesh::gradient_transpose in
// edge order, and free_charge_guard. These are sequential O(N) / O(mesh) folds over
// downloaded arrays; extended precision is native here, so they run on the host.
#pragma once
#include "hall_gpu_push.hpp"

#include <cmath>
#include <limits>

namespace hall_gpu::push {

struct FoldInput {
    std::size_t n = 0;
    std::span<const std::uint64_t> ids;
    std::span<const double> weights, dke, E, current, initial, endpoint, gross, net_hi, net_lo, free_gross;
    const MeshTables* mt = nullptr;
    double dt = 0;
    std::uint32_t status[6] = {0xffffffffu, 0, 0xffffffffu, 0, 0xffffffffu, 0};   // first failing rank/code per stage
};

inline void fold_record(const FoldInput& in, Refusal& refusal, Refusal& charge_refusal, Budget& b) {
    constexpr std::uint32_t none = 0xffffffffu;
    const auto* st = in.status;
    auto carrier_refusal = [&](const char* stage, std::uint32_t rank, std::uint32_t code) {
        return Refusal{true, code_message(int(code)), stage, long(rank), in.ids[rank], int(code)};
    };
    if (st[0] != none) { refusal = carrier_refusal("census", st[0], st[1]); return; }
    if (!(std::isfinite(b.plasma_dt) && b.plasma_dt <= .1)) { refusal = {true, "motion full-census plasma guard", "plasma", -1, 0, 0}; return; }
    if (st[2] != none) { refusal = carrier_refusal("motion", st[2], st[3]); return; }
    long double kinetic = 0;
    for (std::size_t i = 0; i < in.n; ++i) kinetic += static_cast<long double>(in.weights[i]) * in.dke[i];
    b.kinetic_change_J = double(kinetic);
    const std::size_t ne = in.current.size(), nn = in.initial.size();
    long double dot = 0;
    for (std::size_t e = 0; e < ne; ++e) dot += (long double)in.E[e] * in.current[e];
    b.field_work_J = in.dt * double(dot);
    b.work_defect_J = b.kinetic_change_J - b.field_work_J;
    if (!(std::isfinite(b.work_defect_J) && std::abs(b.work_defect_J) <= 1e-20 + 1e-10 * (std::abs(b.kinetic_change_J) + std::abs(b.field_work_J)))) {
        refusal = {true, "canonical motion work mismatch", "work", -1, 0, 0};
        return;
    }
    std::vector<double> div(nn, 0.0);
    for (std::size_t e = 0; e < ne; ++e) {
        const double q = in.current[e] / in.mt->edge_length[e];
        div[std::size_t(in.mt->edge_a[e])] -= q;
        div[std::size_t(in.mt->edge_b[e])] += q;
    }
    b.continuity_C = 0;
    for (std::size_t k = 0; k < nn; ++k) b.continuity_C = std::max(b.continuity_C, std::abs(in.endpoint[k] - in.initial[k] - in.dt * div[k]));
    double gross_max = 0;
    for (double x : in.gross) gross_max = std::max(gross_max, std::abs(x));
    bool current_finite = true;
    for (double x : in.current) current_finite = current_finite && std::isfinite(x);
    if (!(current_finite && std::isfinite(b.continuity_C) && b.continuity_C <= 1e-12 * std::max(gross_max, 1e-30))) {
        refusal = {true, "canonical motion charge continuity mismatch", "continuity", -1, 0, 0};
        return;
    }
    if (st[4] != none) { charge_refusal = carrier_refusal("charge", st[4], st[5]); return; }
    for (std::size_t k = 0; k < nn; ++k) {
        const double hi = in.net_hi[k], lo = in.net_lo[k], gr = in.free_gross[k];
        const double s = hi + lo, v = s - hi, err = (hi - (s - v)) + (lo - v);
        const bool canonical = std::isfinite(hi) && std::isfinite(lo) && s == hi && err == lo;
        const long double wide = (long double)hi + lo;
        if (!(canonical && std::isfinite(gr) && gr >= 0 && std::abs(wide) <= gr * (1 + 32 * std::numeric_limits<double>::epsilon()))) {
            charge_refusal = {true, "free charge/gross operand invalid", "charge_guard", -1, 0, 0};
            return;
        }
    }
}

}  // namespace hall_gpu::push
