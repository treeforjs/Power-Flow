#pragma once
// v2/v3 physical-state parity scenarios. Included after `namespace hc = <owner
// namespace>;` by the v3 test and by a dump program compiled against the
// unchanged Hall v2 sources. The digest covers values only (no identity text):
// field words, clock, revision, the complete charge ledger, every heat slice
// (temperatures, inventory, cumulative ledger, watermarks, clock), and the
// motion record JSON text and arrays of each step.
#include "motion_fixtures.hpp"
#include <functional>

namespace parity {
namespace n = raikou::mfem_backend;
namespace mf = motion_fixtures;
using J = hc::Json;
inline std::string physical_sha256(const hc::Owner& o) {
    hc::Digest d("physical-state-v1");
    const auto& s = o.field();
    d.f64s(s.E); d.f64s(s.B); d.f64s(s.E_low); d.f64s(s.B_low); d.f64(s.time_s); d.u64(s.revision);
    const auto& q = o.charge();
    hc::encode(d, q.particles.net_C); d.f64s(q.particles.gross_C); hc::encode(d, q.conductor_C); hc::encode(d, q.exterior_C);
    for (const auto& h : o.thermal().states()) {
        d.f64s(h.T); d.u64(h.inventory.size()); for (const auto& c : h.inventory) { d.f64(c.hi); d.f64(c.lo); }
        d.u64(h.watermarks.size()); for (const auto& [k, v] : h.watermarks) { d.text(k); d.u64(v); }
        for (double x : h.cumulative) d.f64(x);
        d.f64(h.time); d.u64(h.revision);
    }
    return d.hex();
}
inline std::string arrays_sha256(const hall_driver::MotionReceipt& m) {
    hc::Digest d("motion-arrays-v1");
    d.f64s(m.current_A_m()); d.f64s(m.initial_charge_C()); d.f64s(m.endpoint_charge_C());
    d.u64(m.endpoints().size());
    for (const auto& p : m.endpoints()) { d.u64(p.id); for (double v : p.position_m) d.f64(v); for (double v : p.proper_velocity_m_s) d.f64(v); }
    const auto& b = m.budget();
    for (double x : {b.kinetic_change_J, b.field_work_J, b.work_defect_J, b.continuity_C, b.plasma_dt, b.maximum_gyro, b.maximum_displacement_cells}) d.f64(x);
    return d.hex();
}
inline hc::ChargeLedger declared(const hc::Model& m, const axial::State& s, hc::FreeCharge particles) {
    const auto& g = m.mesh(); hc::ChargeLedger q; q.particles = std::move(particles);
    q.conductor_C.assign(g.nn, hc::Words()); q.exterior_C.assign(g.nn, hc::Words());
    q.conductor_declaration = hc::sha256("parity: declared DC conductor/PEC induced charge");
    q.exterior_declaration = hc::sha256("parity: declared DC exterior termination charge");
    const auto field = hc::field_charge(g, s);
    for (int k = 0; k < g.nn; ++k) { if (g.node_body[k] >= 0 || g.short_node[k]) q.conductor_C[k] = field[k]; else if (g.outer_node[k]) q.exterior_C[k] = field[k]; }
    return q;
}
inline hall_driver::MotionBinding binding(const hc::Owner& o, int step) {
    hall_driver::MotionBinding b; b.phase = hall_driver::ParticlePhase::full_time_position_and_proper_velocity;
    b.accepted_step = o.revision(); b.staged_arena_revision = o.revision() + 1; b.begin_s = o.time_s(); b.dt_s = 1e-14; b.end_s = b.begin_s + b.dt_s;
    b.staged_arena_sha256 = hc::sha256("parity/staged/" + std::to_string(step)); b.mesh_sha256 = o.motion_mesh_sha256();
    b.material_sha256 = o.model().material_sha256(); b.species_sha256 = hc::sha256("parity/species");
    b.common_clock_sha256 = hc::sha256("parity/clock/" + std::to_string(step)); b.kinetic_policy_sha256 = o.motion_context().policy_sha256;
    return b;
}
struct Boundary { std::string field_sha256, heat_sha256, source_sha256, mesh_sha256, material_sha256, common_clock_sha256; double time_s = 0; std::uint64_t accepted_step = 0; };
// Runs the scenarios; returns one row per published step (values only).
inline J run(const std::string& trace) {
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace);
    const auto& g = model->mesh(); const auto o = mf::geometry(g); const auto seed = model->dc_seed();
    struct Scenario { std::string name; std::vector<n::ArenaRateEntry> census; int steps; };
    const double cx = g.side / 2, cy = g.side / 2;
    std::vector<Scenario> scenarios{{"pair-z0-w1", mf::pair({cx, cy, 0.}, 1.), 3}, {"pair-z1mm-w1", mf::pair({cx, cy, .001}, 1.), 2},
                                    {"pair-z0-w1e9", mf::pair({cx, cy, 0.}, 1e9), 2}, {"node-planes-500-w1", mf::node_plane_pairs(250, 21, 1.), 2},
                                    {"pairs-1e3-w1e9", mf::pairs(1000, 7, 1e9), 2}, {"vacuum", {}, 3}};
    J rows = J::array();
    for (auto& sc : scenarios) {
        hc::Owner owner(model, seed, declared(*model, seed, hc::charge_of_census(g, o, sc.census)), o, sc.census);
        auto census = sc.census;
        for (int step = 1; step <= sc.steps; ++step) {
            J r{{"scenario", sc.name}, {"step", step}};
            try {
                if (sc.census.empty()) {
                    auto f = owner.prepare_vacuum(1e-14);
                    Boundary bd{f.candidate_field_sha256(), f.candidate_heat_sha256(), "", "", "", "", f.end_s(), f.endpoint_revision()};
                    hc::need(owner.preflight_vacuum_final(f, bd), "vacuum preflight");
                    r["field_receipt"] = {f.field_receipt().field_delta_J, f.field_receipt().joule_J, f.field_receipt().balance_J};
                    owner.publish_vacuum_prevalidated(std::move(f));
                } else {
                    const std::string src = hc::sha256("parity/source");
                    auto p = hc::prepare_closed_proposal(owner, o, census, binding(owner, step), src);
                    const auto& b = p.motion().binding();
                    Boundary bd{p.candidate_field_sha256(), p.candidate_heat_sha256(), src, b.mesh_sha256, b.material_sha256, b.common_clock_sha256, p.end_s(), p.endpoint_revision()};
                    hc::need(owner.preflight_final(p, census, bd), "closed preflight");
                    J rec = hall_driver::motion_record(p.motion()); const auto text = rec.dump(); hc::drain(rec);
                    r["motion_record_text_sha256"] = hc::sha256(text); r["motion_seal_v2"] = p.motion().seal_sha256();
                    r["motion_record_v3_sha256"] = p.motion().record_v3_sha256(); r["motion_arrays_sha256"] = arrays_sha256(p.motion());
                    r["iterations"] = p.iterations(); r["combined_balance_J"] = p.combined_balance_J();
                    r["field_heat_balance_J"] = p.field_heat().field_heat_balance_J();
                    std::vector<n::ArenaRateEntry> next = census;
                    std::sort(next.begin(), next.end(), [](const auto& a, const auto& c) { return a.native.id < c.native.id; });
                    for (std::size_t i = 0; i < next.size(); ++i) next[i].native = p.motion().endpoints()[i];
                    owner.publish_prevalidated(std::move(p));
                    census = std::move(next);
                }
                r["physical_sha256"] = physical_sha256(owner);
            } catch (const std::exception& e) { r["refusal"] = e.what(); }
            rows.push_back(r);
            if (r.contains("refusal")) break;
        }
    }
    return rows;
}
}  // namespace parity
