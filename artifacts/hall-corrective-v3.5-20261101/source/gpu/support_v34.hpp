// support_v34.hpp -- the deposit support of the endpoint continuity guard (census_charge.hpp
// DepositSupportBuilder, Hall v2/v3, unchanged rule endpoint-continuity-v2) and the x87 comparisons of the
// census admission, in a form compiled by NVRTC (device kernels of step_kernels_v34.cu) and by the host
// (CPU emulation for the qualification). Every long-double operation of the sealed builder is reproduced
// with x87_v34.hpp (64-bit significand, nearest even), every binary64 operation with the same binary64
// operation, in the sealed order:
//   per carrier p (ascending native ID; start = census position, end = endpoint, sink hit point, or the
//   newborn's endpoint/hit point):
//     Q  = |hi| + |lo| (x87) of Words(w) * q;  need(isfinite(double(Q)))      -> code charge_operand
//     Q == 0: no contribution
//     per axis d: a, b = cell coordinates of start/end (binary64);  need(finite) -> code nonfinite_coordinate
//       delta = eps (3 c_d + 8) (exact);  R += |b - a| + 2 delta (x87);  lo/hi = clamp(floor(min - delta),
//       floor(max + delta)) (x87)
//     need(R <= 1)                                                             -> code path_exceeds
//     t1 = Q (2C + 14 + (16C + 67) R), t2 = Q R (x87)
//   per node n, the carriers whose node box (lo..hi+1 per axis) holds n, in ascending native ID:
//     s1 += t1, s2 += t2, m += 1 (x87, sequential)
//   absolute_C(n) = RU(eps (s1 + 24 m s2) (1 + 1e-9L)) (x87; RU = nearest, then nextafter(+inf) if below)
// The sum per node is sequential in ascending carrier order (the sealed order); the device walks each
// node's carriers after a stable sort of (node, carrier) records, so the values are bit for bit those of
// the host builder (no identity change).
#pragma once
#include "x87_v34.hpp"
#include "push_core.hpp"

namespace sup34 {

enum Code : int { ok = 0, charge_operand = 1, nonfinite_coordinate = 2, path_exceeds = 3 };

// Geometry and constants of one mesh (host-filled, passed to the kernels by value).
struct Geo {
    int counts[3];
    double spacing[3], origin[3];
    double shift_x, shift_y;          // native -> field frame: x - width/2, y - gap/2
    x87::XBits delta[3], two_delta[3];
    x87::XBits c2p14, c16p67;         // 2C + 14, 16C + 67 (C = nx + ny + nz)
    x87::XBits k_margin;              // 1 + 1e-9L (x87 constant of the sealed finish)
    x87::XBits one;
};

// Per-carrier result.
struct Carrier {
    int code;
    int contributes;                  // Q != 0 and code == ok
    x87::X t1, t2;
    int lo[3], hi[3];
};

X87_FN bool finite_d(double x) { return ((x87_bits(x) >> 52) & 0x7ff) != 0x7ff; }

// Q = magnitude(Words(w) * q) (long double); returns false when double(Q) is not finite.
X87_FN bool charge_magnitude(double w, double q, x87::X& Q) {
    const hgp::Words ww = {w, 0};
    const hgp::Words p = hgp::words_mul(ww, q);
    if (!finite_d(p.hi) || !finite_d(p.lo)) return false;
    Q = x87::add(x87::from_double(HGP_FABS(p.hi)), x87::from_double(HGP_FABS(p.lo)));
    return finite_d(x87::to_double(Q, false));
}

X87_FN int clampi(x87_i64 v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : int(v)); }

// start_native/end_native: native-frame points (census position; endpoint or hit point).
X87_FN void carrier(const Geo& g, double w, double q, const double* start_native, const double* end_native, Carrier& out) {
    out.code = ok; out.contributes = 0; out.t1 = x87::zero(); out.t2 = x87::zero();
    for (int d = 0; d < 3; ++d) { out.lo[d] = 0; out.hi[d] = -1; }
    x87::X Q;
    if (!charge_magnitude(w, q, Q)) { out.code = charge_operand; return; }
    if (x87::is_zero(Q)) return;
    const double start[3] = {start_native[0] - g.shift_x, start_native[1] - g.shift_y, start_native[2]};
    const double end[3] = {end_native[0] - g.shift_x, end_native[1] - g.shift_y, end_native[2]};
    x87::X R = x87::zero();
    for (int d = 0; d < 3; ++d) {
        const double a = (start[d] - g.origin[d]) / g.spacing[d], b = (end[d] - g.origin[d]) / g.spacing[d];
        if (!(finite_d(a) && finite_d(b))) { out.code = nonfinite_coordinate; return; }
        const x87::X xa = x87::from_double(a), xb = x87::from_double(b);
        const x87::X delta = x87::unpack(g.delta[d]);
        R = x87::add(R, x87::add(x87::abs(x87::sub(xb, xa)), x87::unpack(g.two_delta[d])));
        const double mn = (b < a) ? b : a, mx = (a < b) ? b : a;   // std::min / std::max
        out.lo[d] = clampi(x87::floor_int(x87::sub(x87::from_double(mn), delta)), 0, g.counts[d] - 1);
        out.hi[d] = clampi(x87::floor_int(x87::add(x87::from_double(mx), delta)), 0, g.counts[d] - 1);
    }
    if (!x87::le(R, x87::unpack(g.one))) { out.code = path_exceeds; return; }
    out.t1 = x87::mul(Q, x87::add(x87::unpack(g.c2p14), x87::mul(x87::unpack(g.c16p67), R)));
    out.t2 = x87::mul(Q, R);
    out.contributes = 1;
}
// Number of nodes of the carrier's box.
X87_FN unsigned box_nodes(const Carrier& c) {
    if (!c.contributes) return 0;
    return unsigned(c.hi[0] - c.lo[0] + 2) * unsigned(c.hi[1] - c.lo[1] + 2) * unsigned(c.hi[2] - c.lo[2] + 2);
}
// Node index of slot k of the box, in the sealed loop order (k outermost, then j, then i).
X87_FN int box_node(const Carrier& c, unsigned slot, int nx, int ny) {
    const unsigned wi = unsigned(c.hi[0] - c.lo[0] + 2), wj = unsigned(c.hi[1] - c.lo[1] + 2);
    const unsigned i = slot % wi, j = (slot / wi) % wj, k = slot / (wi * wj);
    return (c.lo[2] + int(k)) * (nx + 1) * (ny + 1) + (c.lo[1] + int(j)) * (nx + 1) + (c.lo[0] + int(i));
}
// Finish of one node: RU(eps (s1 + 24 m s2) (1 + 1e-9L)); 0 for a node without carriers.
X87_FN double finish(const Geo& g, x87::X s1, x87::X s2, x87_u64 m) {
    if (m == 0) return 0.0;
    const x87::X inner = x87::add(s1, x87::mul(x87::from_u64(24 * m), s2));
    const x87::X bound = x87::mul(x87::scale2(inner, -52), x87::unpack(g.k_margin));
    return x87::to_double(bound, true);
}

// ---------------------------------------------------------------- census admission comparisons
// same_net_charge (census_charge.hpp) at one node:
//   |(a - b).wide()| > 256 eps ((long double)a.gross + b.gross)  ->  differs
X87_FN bool same_net_node(double ahi, double alo, double agross, double bhi, double blo, double bgross) {
    const hgp::Words a = {ahi, alo};
    const hgp::Words nb = hgp::two_sum_words(-bhi, -blo);   // operator-(Words)
    const hgp::Words d = hgp::words_add(a, nb);
    const x87::X wide = x87::add(x87::from_double(d.hi), x87::from_double(d.lo));
    const x87::X tol = x87::scale2(x87::add(x87::from_double(agross), x87::from_double(bgross)), -44);   // 256 eps = 2^-44
    return x87::le(x87::abs(wide), tol);
}
// free_charge_guard (field_contract.hpp) at one node: canonical net words, finite non-negative gross,
// |wide(net)| <= gross * (1 + 32 eps) (binary64 product, compared in long double).
X87_FN bool free_guard_node(double hi, double lo, double gross) {
    if (!(finite_d(hi) && finite_d(lo) && finite_d(gross) && gross >= 0)) return false;
    const hgp::Words n = hgp::two_sum_words(hi, lo);
    if (!(n.hi == hi && n.lo == lo)) return false;
    const x87::X wide = x87::add(x87::from_double(hi), x87::from_double(lo));
    const double limit = gross * (1 + 0x1p-47);   // 1 + 32 eps (exact)
    return x87::le(x87::abs(wide), x87::from_double(limit));
}

}  // namespace sup34
