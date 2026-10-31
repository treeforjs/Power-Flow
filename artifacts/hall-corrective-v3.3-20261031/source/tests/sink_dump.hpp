#pragma once
// Sink fixture dump shared by the CPU qualification (CpuMotionProvider) and the GPU
// program (GpuMotionProvider): per fixture the motion record digests, the sink
// record digest, the endpoint free-charge words, and the closed proposal digests
// (candidate field/heat/ledger) with the CPU field solve. CPU/GPU parity = equal rows.
#include "sink_fixtures.hpp"
#include "parity_scenarios.hpp"

namespace sink_dump {
namespace n = raikou::mfem_backend;
using J = hc::Json;
inline std::string free_sha(const hc::FreeCharge& q) { hc::Digest d("test-after-charge"); hc::encode(d, q.net_C); d.f64s(q.gross_C); return d.hex(); }
inline J run(const std::shared_ptr<const hc::Model>& model, hc::MotionProvider& provider) {
    const auto& g = model->mesh(); const auto o = motion_fixtures::geometry(g); const auto seed = model->dc_seed();
    J rows = J::array();
    for (const auto& f : sink_fixtures::fixture_set()) {
        J r{{"fixture", f.name}, {"carriers", f.census.size()}};
        try {
            std::vector<hc::Newborn> newborns;
            const auto t0 = model->drive().start();
            for (const auto& [id, frac] : f.newborns) newborns.push_back({id, t0 + frac * 1e-14});
            const auto residents = hc::residents(f.census, newborns);
            hc::Owner owner(model, seed, parity::declared(*model, seed, hc::charge_of_census(g, o, residents)), o, residents);
            const auto mid = model->at(seed);
            hc::Vec E = mid.E(); if (f.field) for (double& v : E) v *= 100;
            auto b = parity::binding(owner, 1); b.trial_field_sha256 = hall_driver::total_field_sha256(E, mid.B());
            provider.set_newborns(newborns); provider.begin(owner.motion_context(), f.census);
            try {
                hc::FreeCharge after;
                auto m = provider.motion(b, E, mid.B(), after);
                hc::MotionAssembly::seal_v2(m);
                J rec = hall_driver::motion_record(m); const auto text = rec.dump(); hc::drain(rec);
                r["record_text_sha256"] = hc::sha256(text); r["record_v3_sha256"] = m.record_v3_sha256(); r["seal"] = m.seal_sha256();
                r["arrays_sha256"] = parity::arrays_sha256(m); r["after_sha256"] = free_sha(after); r["sinks"] = m.sinks().size();
                r["sink_section_sha256"] = m.v2_identical() ? std::string() : hc::MotionAssembly::sink_birth_digest(m);
                { std::string nb; for (const auto& b : m.births()) nb += std::to_string(b.native_id) + ","; r["newborns"] = nb; }
                J kinds = J::array(); for (const auto& k : m.sinks()) kinds.push_back(k.boundary); r["sink_kinds"] = kinds;
            } catch (const std::exception& e) { r["motion_refusal"] = e.what(); }
            provider.end();
            if (f.field == 0 && !r.contains("motion_refusal")) {
                auto p = hc::prepare_closed_proposal(owner, o, f.census, parity::binding(owner, 1), hc::sha256("sink-dump/source"), provider,
                                                     hc::StepInputs{newborns, {}});
                r["candidate_field_sha256"] = p.candidate_field_sha256(); r["candidate_heat_sha256"] = p.candidate_heat_sha256();
                r["candidate_ledgers_sha256"] = p.field_heat().candidate_ledgers_sha256(); r["iterations"] = p.iterations();
                r["combined_balance_J"] = p.combined_balance_J(); r["energy_residual_J"] = p.energy().residual_J;
            }
        } catch (const std::exception& e) { r["refusal"] = e.what(); }
        rows.push_back(r);
    }
    return rows;
}
}  // namespace sink_dump
