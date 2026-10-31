// fixtures.hpp -- carrier sets and transaction policy for the merge/split tests.
// Electrons: Maxwellian 3 eV (92 %) + exponential keV tail (8 %, mean 2 keV, 0.3-30 keV,
// 70 % isotropic, 30 % in a forward cone cos(theta) > 0.9 along +z).
// H2O+: Maxwellian 0.1 eV (94 %) + exponential tail (6 %, mean 50 eV, 10-2000 eV, isotropic).
// Positions (native frame, 9 mm x 9 mm x 25 mm box, x,y in (0, 9 mm), z in (-12.5, 12.5) mm):
// 50 % uniform, 50 % Gaussian blob at the box centre (sigma 0.8, 0.8, 4 mm).
// Weights: 88 % w0 = 1e6, 6 % w0/2, 2 % w0/4, 2 % w0 (1 + k 2^-30), 2 % w0 U(0.5,1) with a full
// 53-bit mantissa (forces inexact sums, so groups stop early). Flags: 97 % merge+split
// eligible, 3 % ineligible. IDs 1..2n, species interleaved at random, exactly n each.
#pragma once
#include "merge_split.hpp"

#include <algorithm>
#include <cmath>
#include <random>

namespace fixtures {
using merge_split::Carriers;
inline constexpr double e = 1.602176634e-19, me = 9.1093837139e-31, c = 299792458., da = 1.66053906892e-27;
inline const double m_h2o_plus = 18.010564684 * da - me;
inline constexpr double side = .009, depth = .025;

struct Rng {
    std::mt19937_64 g;
    explicit Rng(std::uint64_t s) : g(s) {}
    double uniform() { return (double(g() >> 11) + .5) * 0x1p-53; }   // (0,1)
    double normal() {
        const double a = uniform(), b = uniform();
        return std::sqrt(-2 * std::log(a)) * std::cos(6.283185307179586 * b);
    }
    double exponential_trunc(double mean, double lo, double hi) {
        for (;;) {
            const double x = lo - mean * std::log(uniform());
            if (x <= hi) return x;
        }
    }
};

inline void speed_to_u(double m, double ke_J, double& umag) {
    const double gam = 1 + ke_J / (m * c * c);
    umag = c * std::sqrt(gam * gam - 1);
}
inline void isotropic(Rng& r, double mag, double u[3], double cos_min = -1) {
    const double ct = cos_min + (1 - cos_min) * r.uniform(), st = std::sqrt(std::max(0., 1 - ct * ct));
    const double ph = 6.283185307179586 * r.uniform();
    u[0] = mag * st * std::cos(ph);
    u[1] = mag * st * std::sin(ph);
    u[2] = mag * ct;
}
inline void velocity(Rng& r, unsigned species, double u[3]) {
    if (species == 0) {
        if (r.uniform() < .92) {
            const double s = std::sqrt(3 * e / me);
            for (int d = 0; d < 3; ++d) u[d] = s * r.normal();
        } else {
            double mag;
            speed_to_u(me, r.exponential_trunc(2000, 300, 30000) * e, mag);
            isotropic(r, mag, u, r.uniform() < .3 ? .9 : -1);
        }
    } else {
        if (r.uniform() < .94) {
            const double s = std::sqrt(.1 * e / m_h2o_plus);
            for (int d = 0; d < 3; ++d) u[d] = s * r.normal();
        } else {
            double mag;
            speed_to_u(m_h2o_plus, r.exponential_trunc(50, 10, 2000) * e, mag);
            isotropic(r, mag, u);
        }
    }
}
inline void position(Rng& r, double x[3]) {
    const double margin = 1e-7;
    if (r.uniform() < .5) {
        x[0] = margin + (side - 2 * margin) * r.uniform();
        x[1] = margin + (side - 2 * margin) * r.uniform();
        x[2] = -depth / 2 + margin + (depth - 2 * margin) * r.uniform();
        return;
    }
    for (;;) {
        x[0] = side / 2 + .0008 * r.normal();
        x[1] = side / 2 + .0008 * r.normal();
        x[2] = .004 * r.normal();
        if (x[0] > margin && x[0] < side - margin && x[1] > margin && x[1] < side - margin && std::abs(x[2]) < depth / 2 - margin)
            return;
    }
}
inline double weight(Rng& r) {
    const double w0 = 1e6, p = r.uniform();
    if (p < .88) return w0;
    if (p < .94) return w0 / 2;
    if (p < .96) return w0 / 4;
    if (p < .98) return w0 * (1 + double(1 + (r.g() % 1023)) * 0x1p-30);
    return w0 * (0.5 + 0.5 * double(r.g() >> 11) * 0x1p-53);
}

// n carriers per species (2 species), IDs 1..2n.
inline Carriers make_set(std::uint64_t seed, std::size_t n) {
    Rng r(seed);
    std::vector<std::uint32_t> sp(2 * n);
    for (std::size_t i = 0; i < 2 * n; ++i) sp[i] = i < n ? 0u : 1u;
    std::shuffle(sp.begin(), sp.end(), r.g);
    Carriers C;
    C.reserve(2 * n);
    for (std::size_t i = 0; i < 2 * n; ++i) {
        double x[3], u[3];
        position(r, x);
        velocity(r, sp[i], u);
        const double w = weight(r);
        const std::uint32_t f = r.uniform() < .97 ? (MS_FLAG_MERGE | MS_FLAG_SPLIT) : 0u;
        C.push(i + 1, sp[i], f, x, u, w);
    }
    return C;
}

// Policy used by the tests: energy bin edges 0.1 eV (electrons) / 0.01 eV (H2O+) times 2^k,
// k = 0..19; 6 direction faces; merge cells above 48 carriers of a species down to 32
// (groups of 3..16); split cells below 4 carriers toward 8 (two clones, weight >= 2).
inline MsParams policy(int n_transverse, std::uint64_t seed, std::uint64_t step, std::uint64_t base_id) {
    MsParams P = merge_split::hall_mesh_params(n_transverse);
    P.seed = seed;
    P.step = step;
    P.base_id = base_id;
    P.nsub = 1;
    merge_split::set_species(P, 0, me, -e);
    merge_split::set_species(P, 1, m_h2o_plus, e);
    for (unsigned s = 0; s < 2; ++s) {
        auto& sp = P.species[s];
        sp.n_thresh = 20;
        for (unsigned k = 0; k < 20; ++k) sp.e_thresh[k] = (s == 0 ? .1 : .01) * e * double(1u << k);
        sp.merge_above = 48;
        sp.merge_target = 32;
        sp.k_min = 3;
        sp.k_max = 16;
        sp.split_below = 4;
        sp.split_target = 8;
        sp.split_levels = 1;
        sp.split_min_weight = 2;
    }
    return P;
}

// Corner-case set: one merge group per cell (energy bins and direction bins off, so a
// cell's carriers form one group). Kinds by g % 5: 0 Maxwellian/tail mixture, 1 identical
// velocities (q^2 at rounding level, often negative: skipped), 2 symmetric +/- pairs
// (P = 0 exactly: isotropic axis), 3 relativistic electrons 50 keV - 2 MeV, 4 full-mantissa
// weights (inexact running sums: groups stop early, leftovers stay).
inline MsParams corner_policy(std::uint64_t seed) {
    MsParams P = policy(45, seed, 77, 1u << 30);
    P.nsub = 0;
    for (unsigned s = 0; s < 2; ++s) {
        auto& sp = P.species[s];
        sp.n_thresh = 0;
        sp.merge_above = 2;
        sp.merge_target = 2;
        sp.k_min = 3;
        sp.k_max = 64;
        sp.split_below = 0;
    }
    return P;
}
inline Carriers make_corner_set(std::uint64_t seed, unsigned groups) {
    Rng r(seed);
    Carriers C;
    const MsParams P = corner_policy(0);
    std::uint64_t id = 1;
    for (unsigned g = 0; g < groups && g < 4 * 45 * 45; ++g) {
        const unsigned kind = g % 5, species = kind == 3 ? 0u : (g / 5) & 1u;
        const unsigned cell = 45 * 45 + g, i = cell % 45, j = (cell / 45) % 45, k = cell / (45 * 45);
        const unsigned count = kind == 2 ? 4u : (kind == 4 ? 6u + unsigned(r.g() % 15) : 3u + unsigned(r.g() % 14));
        double same[3], other[3];
        velocity(r, species, same);
        velocity(r, species, other);
        for (unsigned m = 0; m < count; ++m) {
            double x[3] = {(i + .05 + .9 * r.uniform()) * P.h[0], (j + .05 + .9 * r.uniform()) * P.h[1],
                           P.origin[2] + (k + .05 + .9 * r.uniform()) * P.h[2]};
            double u[3];
            double w = 1e6;
            if (kind == 0) velocity(r, species, u);
            else if (kind == 1) for (int d = 0; d < 3; ++d) u[d] = same[d];
            else if (kind == 2) for (int d = 0; d < 3; ++d) u[d] = (m < 2 ? same[d] : other[d]) * (m % 2 ? -1 : 1);
            else if (kind == 3) {
                double mag;
                speed_to_u(me, (5e4 + 2e6 * r.uniform()) * e, mag);
                isotropic(r, mag, u, .5);
            } else {
                velocity(r, species, u);
                w = 1e6 * (0.5 + 0.5 * double(r.g() >> 11) * 0x1p-53);
            }
            C.push(id++, species, MS_FLAG_MERGE | MS_FLAG_SPLIT, x, u, w);
        }
    }
    return C;
}

}  // namespace fixtures
