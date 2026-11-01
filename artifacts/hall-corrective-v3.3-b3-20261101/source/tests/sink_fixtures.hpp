#pragma once
// Kinetic sink fixtures (native frame, default 45-cell mesh, width = gap = 9 mm):
// lower foil x [2, 7] mm, y [1.8, 2.0] mm; upper foil y [7.0, 7.2] mm; target ends
// z = -12.5 mm (short plane; bridge |x-4.5| <= 2.5, 2.0 <= y <= 7.0 mm) and +12.5 mm.
// Shared by the CPU qualification and the GPU parity program.
#include "motion_fixtures.hpp"

namespace sink_fixtures {
namespace n = raikou::mfem_backend;
namespace mf = motion_fixtures;
inline n::PicVector u_dir(double eV, double mass, n::PicVector dir) {
    const double norm = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
    const double x = eV * n::pic_e / (mass * n::pic_c * n::pic_c), u = n::pic_c * std::sqrt(x * (x + 2));
    return {u * dir[0] / norm, u * dir[1] / norm, u * dir[2] / norm};
}
inline n::ArenaRateEntry electron(std::uint64_t id, n::PicVector x, double eV, n::PicVector dir, double w = 1.) {
    return mf::row(id, "e-" + std::to_string(id), "electron", x, u_dir(eV, mf::me, dir), -n::pic_e, mf::me, w);
}
inline n::ArenaRateEntry ion(std::uint64_t id, n::PicVector x, n::PicVector u, double w = 1.) {
    return mf::row(id, "i-" + std::to_string(id), "H2O+", x, u, n::pic_e, mf::ion_mass(), w);
}
// Exact kinetic face planes of the default geometry (same expressions as motion_fixtures::geometry).
inline const double lower_left = .009 / 2 - .0025, lower_right = .009 / 2 + .0025, lower_top = .009 / 2 - .0025, upper_bottom = .009 / 2 + .0025;
struct Fixture { std::string name; std::vector<n::ArenaRateEntry> census; int field = 0; std::vector<int> expected_sinks;
                 std::vector<std::pair<std::uint64_t, double>> newborns; };   // newborn native ID, birth time as a fraction of dt
// Absolute Gauss at bootstrap requires a neutral census at every node: every
// charged test carrier gets a co-located partner of opposite charge and equal
// weight that stays in the interior (slow H2O+ for electrons, slow electron for ions).
inline std::vector<n::ArenaRateEntry> with_partners(std::vector<n::ArenaRateEntry> c) {
    std::vector<n::ArenaRateEntry> out;
    std::uint64_t id = 1;
    for (auto& r : c) {
        r.native.id = id++; out.push_back(r);
        if (r.native.charge_C == 0) continue;
        n::ArenaRateEntry q = r.native.charge_C < 0 ? ion(id, r.native.position_m, {0, 0, 1}, r.native.weight)
                                                    : electron(id, r.native.position_m, 1e-3, {0, 0, 1}, r.native.weight);
        q.packet_id = r.packet_id + "-partner"; q.native.id = id++; out.push_back(q);
    }
    return out;
}
inline std::vector<int> partnered(std::vector<int> kinds) { std::vector<int> o; for (int k : kinds) { o.push_back(k); o.push_back(-1); } return o; }
// expected_sinks: KineticSink code per carrier in ID order (-1: survivor; -2: refusal expected).
inline std::vector<Fixture> fixture_set() {
    std::vector<Fixture> f;
    const double d = 3e-8;
    f.push_back({"lower-foil-top-face", {electron(1, {.0045, .002 + d, .001}, 100, {0, -1, 0})}, 0, {0}});
    f.push_back({"upper-foil-bottom-face", {electron(1, {.0045, .007 - d, .001}, 100, {0, 1, 0})}, 0, {1}});
    f.push_back({"upper-foil-ion-impact", {ion(1, {.0046, .007 - 2e-12, -.002}, {0, 400, 0}, 1e6)}, 0, {1}});
    f.push_back({"short-bridge", {electron(1, {.0045, .0045, -.0125 + d}, 100, {0, 0, -1})}, 0, {3}});
    f.push_back({"outer-xy-escape", {electron(1, {.009 - d, .0045, 0.}, 100, {1, 0, 0})}, 0, {2}});
    f.push_back({"axial-escape-feed", {electron(1, {.0045, .0045, .0125 - d}, 100, {0, 0, 1})}, 0, {4}});
    f.push_back({"axial-escape-short-plane-off-bridge", {electron(1, {.001, .0045, -.0125 + d}, 100, {0, 0, -1})}, 0, {4}});
    f.push_back({"lower-foil-side-face", {electron(1, {.002 - d, .0019, .001}, 100, {1, 0, 0})}, 0, {0}});
    f.push_back({"lower-foil-bottom-face", {electron(1, {.0045, .0018 - d, .001}, 100, {0, 1, 0})}, 0, {0}});
    f.push_back({"hit-on-x-and-z-node-planes", {electron(1, {.0046, .002 + d, 0.}, 100, {0, -1, 0})}, 0, {0}});
    f.push_back({"corner-aim", {electron(1, {.002 - 2e-8, .002 + 2e-8, .001}, 100, {1, -1, 0})}, 0, {0}});
    f.push_back({"first-of-two-surfaces", {electron(1, {.0045, .002 + 1e-8, .0125 - 2e-8}, 300, {0, -1, 1})}, 0, {0}});
    f.push_back({"grazing-parallel-to-face", {electron(1, {.0045, .002 + 1e-12, .001}, 100, {1, 0, 0})}, 0, {-1}});
    f.push_back({"near-miss-above-corner", {electron(1, {.002 - 2e-8, .002 + 1e-9, .001}, 100, {1, 0, 0})}, 0, {-1}});
    f.push_back({"pair-electron-hits-ion-stays", mf::pair({.0045, .002 + 2e-9, .001}, 1e9), 0, {}});
    f.back().census[0].native.proper_velocity_m_s = u_dir(100, mf::me, {0, -1, 0});
    f.back().expected_sinks = {0, -1};
    f.push_back({"lower-foil-E100", {electron(1, {.0045, .002 + d, .001}, 100, {0, -1, 0}, 1e9)}, 1, {0}});
    f.push_back({"neutral-first-hit-collected", {mf::row(1, "n0", "H2O", {.0045, .002 + 1e-12, .001}, {0, -300, 0}, 0., mf::mn, 1e3)}, 0, {0}});
    {   // newborns of the step (native v6 births): face starts leaving the foil, interior neutral, newborn sinks
        Fixture r{"newborns-faces-and-interior", {}, 0, {}, {}};
        auto pr = mf::pair({.0045, .0045, .001}, 1e6, 1, "res"); for (auto& x : pr) r.census.push_back(x);   // residents 1, 2
        r.census.push_back(electron(3, {.003, lower_top, .001}, 1000, {.1, 1, 0}, 1e6)); r.newborns.push_back({3, .25});
        r.census.push_back(electron(4, {.005, lower_top, -.003}, 2e6, {0, 1, .2}, 1e6)); r.newborns.push_back({4, .5});
        r.census.push_back(electron(5, {.004, upper_bottom, .002}, 500, {0, -1, 0}, 1e6)); r.newborns.push_back({5, 0.});
        r.census.push_back(electron(6, {lower_right, .0019, .0}, 1000, {1, 0, 0}, 1e6)); r.newborns.push_back({6, .9});   // right side face of the lower foil
        r.census.push_back(mf::row(7, "h2o-a", "H2O", {.0035, lower_top, .0}, {50, 650, 0}, 0., mf::mn, 1e9)); r.newborns.push_back({7, .1});
        r.census.push_back(mf::row(8, "h2o-b", "H2O", {.0055, upper_bottom, .004}, {0, -650, 30}, 0., mf::mn, 1e9)); r.newborns.push_back({8, .6});
        r.census.push_back(mf::row(9, "h2o-c", "H2O", {.0045, .004, .001}, {300, 0, 0}, 0., mf::mn, 1e9)); r.newborns.push_back({9, .3});   // interior neutral
        r.expected_sinks = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
        f.push_back(std::move(r));
    }
    {   Fixture r{"newborn-sinks-within-flight", {}, 0, {}, {}};
        auto pr = mf::pair({.0045, .0045, .001}, 1e6, 1, "res"); for (auto& x : pr) r.census.push_back(x);
        r.census.push_back(electron(3, {.0045, lower_top, -.0125 + 1e-8}, 1000, {0, 1, -3}, 1e6)); r.newborns.push_back({3, .2});   // face start, hits the short bridge
        r.census.push_back(mf::row(4, "h2o-escape", "H2O", {.0045, lower_top, .0125 - 2e-12}, {0, 100, 650}, 0., mf::mn, 1e9)); r.newborns.push_back({4, .1});   // axial escape
        r.census.push_back(mf::row(5, "h2o-across", "H2O", {.0045, upper_bottom, .0}, {0, -650, 0}, 0., mf::mn, 1e9)); r.newborns.push_back({5, 0.});   // survives
        r.expected_sinks = {-1, -1, 3, 4, -1};
        f.push_back(std::move(r));
    }
    {   // 400 electrons within 60 nm of the lower foil top face, random directions and energies.
        std::mt19937_64 rng(31);
        std::uniform_real_distribution<double> ux(.0021, .0069), uy(0, 6e-8), uz(-.01, .01), ue(0, 1);
        Fixture r{"random-near-lower-foil-400", {}, 0, {}};
        for (std::uint64_t k = 1; k <= 400; ++k) {
            const double ct = 2 * ue(rng) - 1, ph = 6.283185307179586 * ue(rng), st = std::sqrt(std::max(0., 1 - ct * ct));
            r.census.push_back(electron(k, {ux(rng), .002 + 1e-10 + uy(rng), uz(rng)}, 1 + 999 * ue(rng), {st * std::cos(ph), st * std::sin(ph), ct}, 1e3));
        }
        f.push_back(std::move(r));
    }
    {   // 300 electrons near every target surface (foil faces, target ends, outer xy) plus 100 interior pairs.
        std::mt19937_64 rng(37);
        std::uniform_real_distribution<double> ue(0, 1);
        Fixture r{"random-all-surfaces-500", {}, 1, {}};
        std::uint64_t id = 1;
        for (int k = 0; k < 300; ++k) {
            const int s = k % 6;
            n::PicVector x{.0021 + .0048 * ue(rng), .0045, -.01 + .02 * ue(rng)};
            const double off = 1e-10 + 6e-8 * ue(rng);
            if (s == 0) x[1] = .002 + off; else if (s == 1) x[1] = .007 - off; else if (s == 2) x[2] = -.0125 + off;
            else if (s == 3) x[2] = .0125 - off; else if (s == 4) x[0] = .009 - off; else { x[0] = .0005 + .0008 * ue(rng); x[1] = .0002 + .008 * ue(rng); x[2] = .0125 - off; }
            const double ct = 2 * ue(rng) - 1, ph = 6.283185307179586 * ue(rng), st = std::sqrt(std::max(0., 1 - ct * ct));
            r.census.push_back(electron(id++, x, 1 + 2999 * ue(rng), {st * std::cos(ph), st * std::sin(ph), ct}, 1e6));
        }
        for (auto& p : mf::pairs(200, 41, 1e6)) { p.native.id = id++; r.census.push_back(p); }
        f.push_back(std::move(r));
    }
    // Partner every charged carrier (expected classification: partner survives).
    for (auto& x : f) {
        if (x.name == "neutral-first-hit-refused" || x.name == "pair-electron-hits-ion-stays") continue;
        if (!x.newborns.empty()) continue;
        bool single = x.census.size() == 1 && x.census[0].native.charge_C != 0;
        x.census = with_partners(std::move(x.census));
        if (single) x.expected_sinks = partnered(x.expected_sinks);
    }
    return f;
}
}  // namespace sink_fixtures
