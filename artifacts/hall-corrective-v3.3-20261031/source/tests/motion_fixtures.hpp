#pragma once
// Deterministic motion fixtures shared by the v2 tests and the sealed-provider
// reference program. Only native value types and the axial mesh are used, so
// the same header compiles against the sealed Hall headers and against v2.
#include <cmath>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace motion_fixtures {
namespace n = raikou::mfem_backend;
inline constexpr double me = 9.1093837139e-31, mn = 2.9915e-26, threshold_eV = 12.621;
inline double ion_mass() { return mn - me + threshold_eV * n::pic_e / (n::pic_c * n::pic_c); }
inline n::PicVector electron_u(double eV) {
    const double x = eV * n::pic_e / (me * n::pic_c * n::pic_c);
    return {n::pic_c * std::sqrt(x * (x + 2)), 0, 0};
}
template <class Mesh>
n::GapPicOptions geometry(const Mesh& g) {
    n::GapPicOptions o;
    o.nx = o.ny = g.nx; o.width_m = o.gap_m = g.side; o.depth_m = g.depth;
    o.finite_foils = o.finite_axial_extent = true; o.short_z_m = -g.depth / 2; o.feed_z_m = g.depth / 2;
    o.lower_foil = {g.side / 2 - .0025, g.side / 2 + .0025, g.side / 2 - .0027, g.side / 2 - .0025};
    o.upper_foil = {g.side / 2 - .0025, g.side / 2 + .0025, g.side / 2 + .0025, g.side / 2 + .0027};
    return o;
}
inline n::ArenaRateEntry row(std::uint64_t id, const std::string& packet, const std::string& state, n::PicVector x, n::PicVector u,
                             double charge, double mass, double weight) {
    n::ArenaRateEntry r;
    r.packet_id = packet; r.state = state; r.native.id = id; r.native.position_m = x;
    r.native.proper_velocity_m_s = u; r.native.charge_C = charge; r.native.mass_kg = mass; r.native.weight = weight;
    return r;
}
// Co-located electron/H2O+ pair at native position x (100 eV electron, thermal ion).
inline std::vector<n::ArenaRateEntry> pair(n::PicVector x, double weight, std::uint64_t first_id = 1, const std::string& tag = "p") {
    return {row(first_id, tag + "-e", "electron", x, electron_u(100), -n::pic_e, me, weight),
            row(first_id + 1, tag + "-i", "H2O+", x, {0, 372, 0}, n::pic_e, ion_mass(), weight)};
}
// N/2 co-located electron/ion pairs at random positions inside the vacuum gap
// region x,y in [3.5,5.5] mm (native frame), z in [-4,4] mm; electron energies
// 1..300 eV with random directions; thermal ions. Exactly neutral per node pair.
inline std::vector<n::ArenaRateEntry> pairs(std::size_t carriers, std::uint64_t seed, double weight) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> ux(.0035, .0055), uz(-.004, .004), ue(0, 1);
    std::vector<n::ArenaRateEntry> out;
    out.reserve(carriers);
    for (std::size_t k = 0; k < carriers / 2; ++k) {
        const n::PicVector x{ux(rng), ux(rng), uz(rng)};
        const double eV = 1 + 299 * ue(rng), ct = 2 * ue(rng) - 1, ph = 6.283185307179586 * ue(rng);
        const double st = std::sqrt(std::max(0., 1 - ct * ct));
        const double xe = eV * n::pic_e / (me * n::pic_c * n::pic_c), u = n::pic_c * std::sqrt(xe * (xe + 2));
        const std::string tag = "pair-" + std::to_string(k);
        out.push_back(row(2 * k + 1, tag + "-e", "electron", x, {u * st * std::cos(ph), u * st * std::sin(ph), u * ct}, -n::pic_e, me, weight));
        out.push_back(row(2 * k + 2, tag + "-i", "H2O+", x, {372 * (2 * ue(rng) - 1), 372 * (2 * ue(rng) - 1), 372 * (2 * ue(rng) - 1)},
                          n::pic_e, ion_mass(), weight));
    }
    return out;
}
// Positions on or next to node planes of the default mesh (x,y planes at
// multiples of 0.2 mm from -4.5 mm; z planes at -18.75 + 6.25 k mm), native frame.
inline std::vector<n::ArenaRateEntry> node_plane_pairs(std::size_t count, std::uint64_t seed, double weight) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> ix(18, 27), kz(2, 4), which(0, 6);
    std::uniform_real_distribution<double> ue(0, 1);
    const double offsets[] = {0., 1e-15, -1e-15, 1e-12, -1e-12, 1e-9, -1e-7};
    std::vector<n::ArenaRateEntry> out;
    for (std::size_t k = 0; k < count; ++k) {
        n::PicVector x{.0002 * ix(rng) + .0001 * ue(rng), .0002 * ix(rng) + .0001 * ue(rng), -.01875 + .00625 * kz(rng) + .001 * (ue(rng) - .5)};
        const int axis = int(k % 3);
        if (axis < 2) x[axis] = .0002 * ix(rng) + offsets[which(rng)];
        else x[2] = -.01875 + .00625 * kz(rng) + offsets[which(rng)];
        const double eV = 1 + 999 * ue(rng), ct = 2 * ue(rng) - 1, ph = 6.283185307179586 * ue(rng), st = std::sqrt(std::max(0., 1 - ct * ct));
        const double xe = eV * n::pic_e / (me * n::pic_c * n::pic_c), u = n::pic_c * std::sqrt(xe * (xe + 2));
        const std::string tag = "plane-" + std::to_string(k);
        out.push_back(row(2 * k + 1, tag + "-e", "electron", x, {u * st * std::cos(ph), u * st * std::sin(ph), u * ct}, -n::pic_e, me, weight));
        out.push_back(row(2 * k + 2, tag + "-i", "H2O+", x, {0, 372, 0}, n::pic_e, ion_mass(), weight));
    }
    return out;
}
}  // namespace motion_fixtures

namespace motion_fixtures {
struct Fixture { std::string name; std::vector<n::ArenaRateEntry> census; int field = 0; };  // field 0: DC seed midpoint, 1: E x100
template <class Mesh>
std::vector<Fixture> fixture_set(const Mesh& g) {
    std::vector<Fixture> f;
    const double cx = g.side / 2, cy = g.side / 2;
    f.push_back({"pair-z1mm-w1", pair({cx, cy, .001}, 1.)});
    f.push_back({"pair-z0-w1", pair({cx, cy, 0.}, 1.)});
    f.push_back({"pair-z0-w1e9", pair({cx, cy, 0.}, 1e9)});
    f.push_back({"audit-pair-w1e12", pair({.00021 + cx, -.00017 + cy, .001}, 1e12)});
    f.push_back({"pairs-1e3-w1e9", pairs(1000, 7, 1e9)});
    f.push_back({"pairs-1e4-w1", pairs(10000, 11, 1.)});
    f.push_back({"node-planes-500-w1", node_plane_pairs(250, 3, 1.)});
    f.push_back({"node-planes-500-w1e9-E100", node_plane_pairs(250, 5, 1e9), 1});
    f.push_back({"pairs-1e3-w1-E100", pairs(1000, 13, 1.), 1});
    {   auto c = pair({cx, cy, .001}, 1.);
        c.push_back(row(3, "n0", "H2O", {cx + 1e-4, cy, .002}, {300, -200, 100}, 0., mn, 1e3));
        f.push_back({"pair-plus-neutral", c}); }
    f.push_back({"carrier-in-lower-foil", {row(1, "m-e", "electron", {cx, cy - .0026, 0.}, electron_u(10), -n::pic_e, me, 1.)}});
    {   auto c = pair({cx, cy, .001}, 1.);
        c[1].packet_id = c[0].packet_id;  // same packet, different state/kinematics -> divergent Count2 packet
        f.push_back({"divergent-packet", c}); }
    return f;
}
}  // namespace motion_fixtures
