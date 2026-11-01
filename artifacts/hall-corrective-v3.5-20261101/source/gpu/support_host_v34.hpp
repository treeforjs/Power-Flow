// support_host_v34.hpp -- host side of the device deposit support (support_v34.hpp): the constants of
// one mesh in the x87 interchange form (computed with the host's own long double, so they are the values
// of the sealed builder's expressions), and the CPU emulation of the whole builder (sequential per-node
// accumulation in carrier order; qualification against DepositSupportBuilder and against the device).
#pragma once
#include "support_v34.hpp"
#include <cstring>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace sup34 {
static_assert(std::numeric_limits<long double>::digits == 64, "x87 extended long double required");

inline x87::X from_long_double(long double v) {
    unsigned char b[sizeof(long double)] = {};
    std::memcpy(b, &v, sizeof(long double));
    x87_u64 mant = 0; std::memcpy(&mant, b, 8);
    unsigned short se = 0; std::memcpy(&se, b + 8, 2);
    if (mant == 0) return x87::zero();
    const int exponent = se & 0x7fff;
    if (exponent == 0 || exponent == 0x7fff || (mant >> 63) == 0) throw std::runtime_error("x87 operand outside the emulated range");
    x87::X r; r.m = mant; r.e = exponent - 16383 - 63; r.s = (se >> 15) & 1;
    return r;
}
inline long double to_long_double(x87::X a) {
    if (a.m == 0) return 0.0L;
    unsigned char b[sizeof(long double)] = {};
    std::memcpy(b, &a.m, 8);
    const int exponent = a.e + 63 + 16383;
    if (exponent <= 0 || exponent >= 0x7fff) throw std::runtime_error("x87 value outside the long double range");
    const unsigned short se = (unsigned short)(exponent | (a.s ? 0x8000 : 0));
    std::memcpy(b + 8, &se, 2);
    long double v; std::memcpy(&v, b, sizeof(long double));
    return v;
}

template <class Mesh, class Options>
Geo geo_of(const Mesh& g, const Options& o) {
    Geo out{};
    constexpr long double eps = std::numeric_limits<double>::epsilon();
    const int counts[3] = {g.nx, g.ny, g.nz};
    const double spacing[3] = {g.dx, g.dy, g.dz}, origin[3] = {-g.side / 2, -g.side / 2, g.zmin};
    for (int d = 0; d < 3; ++d) {
        out.counts[d] = counts[d]; out.spacing[d] = spacing[d]; out.origin[d] = origin[d];
        const long double delta = eps * (3 * counts[d] + 8);
        out.delta[d] = x87::pack(from_long_double(delta));
        out.two_delta[d] = x87::pack(from_long_double(2 * delta));
    }
    out.shift_x = o.width_m / 2; out.shift_y = o.gap_m / 2;
    const long double C = counts[0] + counts[1] + counts[2];
    out.c2p14 = x87::pack(from_long_double(2 * C + 14));
    out.c16p67 = x87::pack(from_long_double(16 * C + 67));
    out.k_margin = x87::pack(from_long_double(1 + 1e-9L));   // the constant of the sealed finish
    out.one = x87::pack(x87::from_u64(1));
    return out;
}

// CPU emulation: carriers in ascending native ID (the census order), start/end native points.
struct EmulatedSupport {
    std::vector<double> absolute_C;
    long first_failure = -1; int code = ok;
};
inline EmulatedSupport emulate(const Geo& g, int nn, int nx, int ny, std::span<const double> w, std::span<const double> q,
                               std::span<const double> start_xyz, std::span<const double> end_xyz) {
    EmulatedSupport out;
    const std::size_t n = w.size();
    std::vector<x87::X> s1(std::size_t(nn), x87::zero()), s2(std::size_t(nn), x87::zero());
    std::vector<x87_u64> m(std::size_t(nn), 0);
    for (std::size_t i = 0; i < n; ++i) {
        Carrier c;
        carrier(g, w[i], q[i], &start_xyz[3 * i], &end_xyz[3 * i], c);
        if (c.code != ok) { out.first_failure = long(i); out.code = c.code; return out; }
        const unsigned count = box_nodes(c);
        for (unsigned k = 0; k < count; ++k) {
            const int node = box_node(c, k, nx, ny);
            s1[std::size_t(node)] = x87::add(s1[std::size_t(node)], c.t1);
            s2[std::size_t(node)] = x87::add(s2[std::size_t(node)], c.t2);
            ++m[std::size_t(node)];
        }
    }
    out.absolute_C.assign(std::size_t(nn), 0.0);
    for (int k = 0; k < nn; ++k) out.absolute_C[std::size_t(k)] = finish(g, s1[std::size_t(k)], s2[std::size_t(k)], m[std::size_t(k)]);
    return out;
}
}  // namespace sup34
