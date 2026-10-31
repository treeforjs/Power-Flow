// Sealed-provider reference: motion record digests of the shared fixture set,
// produced by the unchanged 8dddec0 Hall motion provider (native_arena_motion.cpp)
// and charge_of_endpoints. Compiled against the sealed Hall headers only.
//   sealed_motion_reference TRACE OUTFILE
#include "root_api.hpp"
#include "motion_fixtures.hpp"
#include <bit>
#include <chrono>
#include <fstream>
#include <iostream>

namespace hc = hall_corrective;
namespace hd = hall_driver;
namespace mf = motion_fixtures;
using J = hc::Json;

static std::string after_digest(const hc::FreeCharge& q) {
    // Same canonical encoding as the v2 test (tag, count, hi/lo pairs, count, gross).
    std::string bytes;
    auto u64 = [&](std::uint64_t v) { for (int i = 0; i < 8; ++i) bytes.push_back(char((v >> (8 * i)) & 255)); };
    auto f64 = [&](double v) { u64(std::bit_cast<std::uint64_t>(v)); };
    const std::string tag = "test-after-charge";
    u64(tag.size()); bytes += tag;
    u64(q.net_C.size()); for (const auto& w : q.net_C) { f64(w.hi); f64(w.lo); }
    u64(q.gross_C.size()); for (double v : q.gross_C) f64(v);
    return hc::sha256(bytes);
}
int main(int argc, char** argv) {
    int rc = 0;
    try {
        if (argc != 3) throw std::runtime_error("usage: sealed_motion_reference TRACE OUTFILE");
        auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, argv[1]);
        const auto& g = model->mesh();
        const auto o = mf::geometry(g);
        const auto mid = model->at(model->dc_seed());
        axial::Vec E100 = mid.E(); for (double& v : E100) v *= 100;
        J out = J::array();
        for (const auto& f : mf::fixture_set(g)) {
            const axial::Vec& E = f.field ? E100 : mid.E();
            hd::MotionBinding b;
            b.phase = hd::ParticlePhase::full_time_position_and_proper_velocity;
            b.accepted_step = 0; b.staged_arena_revision = 1; b.begin_s = model->drive().start(); b.dt_s = 1e-14; b.end_s = b.begin_s + b.dt_s;
            b.staged_arena_sha256 = hc::sha256("dump/staged"); b.trial_field_sha256 = hd::total_field_sha256(E, mid.B());
            b.mesh_sha256 = hd::motion_mesh_sha256(g, o); b.material_sha256 = hd::material_binding_sha256(); b.species_sha256 = hc::sha256("dump/species");
            b.common_clock_sha256 = hc::sha256("dump/clock"); b.kinetic_policy_sha256 = hd::kinetic_policy_sha256(hd::declared_kinetic_policy(o));
            J r{{"fixture", f.name}, {"carriers", f.census.size()}, {"entry", "sealed prepare_canonical_motion"}};
            const auto t = std::chrono::steady_clock::now();
            try {
                auto m = hd::prepare_canonical_motion(g, o, f.census, b, {}, E, mid.B());
                r["record_sha256"] = hc::sha256(hd::motion_record(m).dump()); r["seal_sha256"] = m.seal_sha256();
                try { r["after_sha256"] = after_digest(hc::charge_of_endpoints(g, o, m.endpoints())); }
                catch (const std::exception& e) { r["after_refusal"] = e.what(); }
            } catch (const std::exception& e) { r["refusal"] = e.what(); }
            r["wall_s"] = std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
            out.push_back(r);
        }
        std::ofstream(argv[2]) << out.dump(1) << '\n';
    } catch (const std::exception& e) { std::cerr << "fatal: " << e.what() << '\n'; rc = 2; }
    return rc;
}
