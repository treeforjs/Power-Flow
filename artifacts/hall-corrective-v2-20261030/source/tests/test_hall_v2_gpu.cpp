// Hall corrective v2 GPU backends qualification (one device session, one stream).
//   test_hall_v2_gpu TRACE NVRTC UUID OUTDIR
// Writes OUTDIR/gpu.json and OUTDIR/gpu-motion-dump.json (same format as the CPU dump).
#include "root_api.hpp"
#include "gpu/gpu_backends.hpp"
#include "motion_fixtures.hpp"
#include <chrono>
#include <fstream>
#include <iostream>

namespace hc = hall_corrective_v2;
namespace hg = hall_corrective_v2::gpu;
namespace hd = hall_driver;
namespace n = raikou::mfem_backend;
namespace mf = motion_fixtures;
using J = hc::Json;
using steady = std::chrono::steady_clock;

namespace {
constexpr double dt = 1e-14;
double since(steady::time_point t) { return std::chrono::duration<double>(steady::now() - t).count(); }
struct Results {
    J cases = J::array(); int failed = 0;
    void record(const std::string& g, const std::string& name, bool ok, J ev = J::object()) {
        cases.push_back({{"group", g}, {"name", name}, {"passed", ok}, {"evidence", std::move(ev)}});
        if (!ok) ++failed;
        std::cout << g << " | " << name << ": " << (ok ? "pass" : "FAIL") << std::endl;
    }
};
hc::ChargeLedger declared_dc(const hc::Model& m, const axial::State& s, hc::FreeCharge particles) {
    const auto& g = m.mesh(); hc::ChargeLedger q; q.particles = std::move(particles);
    q.conductor_C.assign(g.nn, hc::Words()); q.exterior_C.assign(g.nn, hc::Words());
    q.conductor_declaration = hc::sha256("gpu test: declared DC conductor/PEC induced charge");
    q.exterior_declaration = hc::sha256("gpu test: declared DC exterior termination charge");
    const auto field = hc::field_charge(g, s);
    for (int k = 0; k < g.nn; ++k) { if (g.node_body[k] >= 0 || g.short_node[k]) q.conductor_C[k] = field[k]; else if (g.outer_node[k]) q.exterior_C[k] = field[k]; }
    return q;
}
hd::MotionBinding binding_for(const hc::Owner& owner) {
    hd::MotionBinding b; b.phase = hd::ParticlePhase::full_time_position_and_proper_velocity;
    b.accepted_step = owner.revision(); b.staged_arena_revision = 1; b.begin_s = owner.time_s(); b.dt_s = dt; b.end_s = b.begin_s + dt;
    b.staged_arena_sha256 = hc::sha256("gpu/staged"); b.mesh_sha256 = owner.motion_mesh_sha256(); b.material_sha256 = owner.model().material_sha256();
    b.species_sha256 = hc::sha256("gpu/species"); b.common_clock_sha256 = hc::sha256("gpu/clock"); b.kinetic_policy_sha256 = owner.motion_context().policy_sha256;
    return b;
}
std::string after_digest(const hc::FreeCharge& q) { hc::Digest d("test-after-charge"); hc::encode(d, q.net_C); d.f64s(q.gross_C); return d.hex(); }
J guard_receipt(const hc::JointProposal& p) {
    const auto& fh = p.field_heat(); const auto& r = fh.field_receipt(); const auto& c = fh.candidate_constraints(); const auto& s = fh.solve_record();
    return {{"backend", s.backend}, {"independent_residual", s.independent_residual}, {"residual_target", s.target},
            {"residual_over_target", s.target > 0 ? s.independent_residual / s.target : 0.}, {"backend_residual", s.backend_residual},
            {"iterations", s.iterations}, {"corrections", s.corrections}, {"max_gauss_C", c.max_gauss_C}, {"max_gauss_bound_C", c.max_gauss_bound_C},
            {"max_divB_T_m", c.max_divB_T_m}, {"max_divB_bound_T_m", c.max_divB_bound_T_m}, {"gauss_increment_C", r.gauss_increment_C},
            {"field_energy_balance_J", r.balance_J}, {"field_heat_balance_J", fh.field_heat_balance_J()}, {"combined_balance_J", p.combined_balance_J()},
            {"fixed_point_iterations", p.iterations()}, {"candidate_field_sha256", p.candidate_field_sha256()}, {"candidate_heat_sha256", p.candidate_heat_sha256()}};
}
J timings(const hc::ProposalTimings& t) {
    return {{"admit_census_s", t.admit_census_s}, {"motion_s", t.motion_s}, {"support_s", t.support_s}, {"field_prepare_s", t.field_prepare_s},
            {"solve_s", t.solve_s}, {"residual_check_s", t.residual_check_s}, {"seal_s", t.seal_s}, {"heat_s", t.heat_s}, {"total_s", t.total_s},
            {"iterations", t.iterations}, {"motion_provider", t.motion_provider}};
}
}  // namespace

int main(int argc, char** argv) {
    int rc = 0; J report = J::object(); Results R; std::filesystem::path out;
    try {
        if (argc != 5) throw std::runtime_error("usage: test_hall_v2_gpu TRACE NVRTC UUID OUTDIR");
        out = argv[4]; std::filesystem::create_directories(out);
        auto device = std::make_shared<hg::Device>(hg::DeviceOptions{argv[3], argv[2]});
        report["device"] = J::parse(device->device_json());
        auto cpu_model = std::make_shared<const hc::Model>(hc::MeshSpec{}, argv[1]);
        auto gpu_solver = hg::make_gpu_field_solver(device);
        auto gpu_model = std::make_shared<const hc::Model>(hc::MeshSpec{}, argv[1], gpu_solver);
        report["field_solver_compile"] = J::parse(hg::gpu_field_solver_compile_json(*gpu_solver));
        const auto& g = cpu_model->mesh(); const auto o = mf::geometry(g); const auto seed = cpu_model->dc_seed();
        R.record("identity", "GPU field backend is bound into the configuration identity",
                 gpu_model->configuration_sha256() != cpu_model->configuration_sha256(),
                 {{"cpu", cpu_model->configuration_sha256()}, {"gpu", gpu_model->configuration_sha256()},
                  {"gpu_backend_descriptor", gpu_solver->descriptor_sha256()}});
        // ---- motion parity on the shared fixture set (same digests as the CPU dump)
        {
            const auto ctx = hc::make_motion_context(g, o, cpu_model->mesh_json());
            hg::GpuMotionProvider gpu(device, ctx, 10000);
            hc::CpuMotionProvider cpu;
            report["push_compile"] = J::parse(gpu.compile_json());
            const auto mid = cpu_model->at(seed); hc::Vec E100 = mid.E(); for (double& v : E100) v *= 100;
            J dump = J::array(); bool same = true; J diffs = J::array();
            for (const auto& f : mf::fixture_set(g)) {
                const hc::Vec& E = f.field ? E100 : mid.E();
                hd::MotionBinding b; b.phase = hd::ParticlePhase::full_time_position_and_proper_velocity;
                b.accepted_step = 0; b.staged_arena_revision = 1; b.begin_s = cpu_model->drive().start(); b.dt_s = dt; b.end_s = b.begin_s + dt;
                b.staged_arena_sha256 = hc::sha256("dump/staged"); b.trial_field_sha256 = hd::total_field_sha256(E, mid.B());
                b.mesh_sha256 = ctx.mesh_sha256; b.material_sha256 = cpu_model->material_sha256(); b.species_sha256 = hc::sha256("dump/species");
                b.common_clock_sha256 = hc::sha256("dump/clock"); b.kinetic_policy_sha256 = ctx.policy_sha256;
                J pair = J::object();
                for (int which = 0; which < 2; ++which) {
                    hc::MotionProvider& p = which ? static_cast<hc::MotionProvider&>(gpu) : static_cast<hc::MotionProvider&>(cpu);
                    J r{{"fixture", f.name}, {"carriers", f.census.size()}, {"entry", p.name()}};
                    const auto t = steady::now();
                    try {
                        p.begin(ctx, f.census); hc::FreeCharge after;
                        auto m = p.motion(b, E, mid.B(), after); p.end();
                        r["record_v3_sha256"] = m.record_v3_sha256();
                        hc::MotionAssembly::seal_v2(m);
                        J rec = hd::motion_record(m); const auto text = rec.dump(); hc::drain(rec);
                        r["record_sha256"] = hc::sha256(text); r["seal_sha256"] = m.seal_sha256(); r["after_sha256"] = after_digest(after);
                    } catch (const std::exception& e) { p.end(); r["refusal"] = e.what(); }
                    r["wall_s"] = since(t);
                    if (which) dump.push_back(r);
                    pair[which ? "gpu" : "cpu"] = r;
                }
                bool eq = true;
                for (const char* k : {"record_sha256", "seal_sha256", "record_v3_sha256", "after_sha256", "refusal"}) eq = eq && pair["gpu"].value(k, J()) == pair["cpu"].value(k, J());
                same = same && eq;
                diffs.push_back({{"fixture", f.name}, {"identical", eq}, {"cpu_wall_s", pair["cpu"]["wall_s"]}, {"gpu_wall_s", pair["gpu"]["wall_s"]},
                                 {"refusal", pair["cpu"].value("refusal", "")}});
            }
            std::ofstream(out / "gpu-motion-dump.json") << dump.dump(1) << '\n';
            R.record("motion", "GPU provider record bytes, seal, endpoint free charge and refusal texts equal the CPU provider on every fixture", same, diffs);
        }
        // ---- above the sealed 1e5 cap: 1.5e5 carriers with a configured cap, GPU vs CPU provider
        {
            auto big = std::make_shared<const hc::Model>(hc::MeshSpec{}, argv[1], hc::make_cpu_pcg(), hc::CarrierLimits{2000000, 16ull << 30});
            const auto ctx = hc::make_motion_context(g, o, big->mesh_json(), big->limits().maximum_carriers);
            const auto mid = big->at(seed);
            J ev = J::object(); bool ok = true;
            for (std::size_t N : {150000ul, 1000000ul}) {
                const auto census = mf::pairs(N, 99, 1.);
                hd::MotionBinding b; b.phase = hd::ParticlePhase::full_time_position_and_proper_velocity;
                b.accepted_step = 0; b.staged_arena_revision = 1; b.begin_s = big->drive().start(); b.dt_s = dt; b.end_s = b.begin_s + dt;
                b.staged_arena_sha256 = hc::sha256("cap/staged"); b.trial_field_sha256 = hd::total_field_sha256(mid.E(), mid.B());
                b.mesh_sha256 = ctx.mesh_sha256; b.material_sha256 = big->material_sha256(); b.species_sha256 = hc::sha256("cap/species");
                b.common_clock_sha256 = hc::sha256("cap/clock"); b.kinetic_policy_sha256 = ctx.policy_sha256;
                hg::GpuMotionProvider gpu(device, ctx, N); hc::CpuMotionProvider cpu;
                J c{{"carriers", N}};
                auto t = steady::now(); hc::FreeCharge ga; gpu.begin(ctx, census); auto gm = gpu.motion(b, mid.E(), mid.B(), ga); gpu.end();
                c["gpu_motion_s"] = since(t); c["gpu_device_run_s"] = gpu.stats().device_run_s; c["gpu_census_admit_s"] = gpu.stats().admit_s;
                c["gpu_last_kernel_total_ms"] = gpu.stats().last_gpu_total_ms;
                t = steady::now(); hc::FreeCharge ca; cpu.begin(ctx, census); auto cm = cpu.motion(b, mid.E(), mid.B(), ca); cpu.end();
                c["cpu_motion_s"] = since(t);
                t = steady::now(); hc::MotionAssembly::seal_v2(gm); c["v2_json_seal_s"] = since(t);
                const bool same = gm.record_v3_sha256() == cm.record_v3_sha256() && after_digest(ga) == after_digest(ca) && gm.current_A_m() == cm.current_A_m();
                c["record_v3_identical"] = same; c["record_v3_sha256"] = gm.record_v3_sha256();
                ok = ok && same; ev[std::to_string(N)] = c;
            }
            R.record("carrier-cap", "above the sealed 1e5 cap (1.5e5 and 1e6 carriers, configured cap 2e6): GPU motion record v3, current and "
                                    "endpoint charge equal the CPU provider", ok, ev);
        }
        // ---- field backend: vacuum step and closed steps, accepted through the guard receipt
        {
            J ev = J::object(); bool ok = true;
            for (int rep = 0; rep < 3; ++rep) {
                hc::Owner owner(gpu_model, seed, declared_dc(*gpu_model, seed, {hc::WordVec(g.nn), hc::Vec(g.nn)}), o, {});
                auto f = owner.prepare_vacuum(dt);
                ev["vacuum_candidate_" + std::to_string(rep)] = f.candidate_field_sha256();
                ev["vacuum_residual_over_target"] = f.solve_record().independent_residual / f.solve_record().target;
            }
            ok = ev["vacuum_candidate_0"] == ev["vacuum_candidate_1"] && ev["vacuum_candidate_1"] == ev["vacuum_candidate_2"];
            R.record("field", "GPU vacuum step accepted by the guards and bit-identical over three in-process owners", ok, ev);
        }
        {
            const auto census = mf::pair({.0045, .0045, 0.}, 1.);
            const auto ctx_owner = [&](std::shared_ptr<const hc::Model> m) {
                return std::make_unique<hc::Owner>(m, seed, declared_dc(*m, seed, hc::charge_of_census(g, o, census)), o, census); };
            auto cpu_owner = ctx_owner(cpu_model);
            hc::CpuMotionProvider cpu;
            auto ref = hc::prepare_closed_proposal(*cpu_owner, o, census, binding_for(*cpu_owner), hc::sha256("gpu/source"), cpu);
            auto cpu_owner2 = ctx_owner(cpu_model);
            hg::GpuMotionProvider gpu_push(device, cpu_owner2->motion_context(), census.size());
            auto mixed = hc::prepare_closed_proposal(*cpu_owner2, o, census, binding_for(*cpu_owner2), hc::sha256("gpu/source"), gpu_push);
            R.record("closed-loop", "GPU push + CPU field: identical candidate field/heat, midpoint and iteration count to CPU push + CPU field",
                     mixed.candidate_field_sha256() == ref.candidate_field_sha256() && mixed.candidate_heat_sha256() == ref.candidate_heat_sha256() &&
                         mixed.midpoint().seal() == ref.midpoint().seal() && mixed.iterations() == ref.iterations() &&
                         mixed.motion().seal_sha256() == ref.motion().seal_sha256(),
                     {{"cpu", guard_receipt(ref)}, {"gpu_push", guard_receipt(mixed)}, {"cpu_timings", timings(ref.timings())}, {"gpu_push_timings", timings(mixed.timings())}});
            J runs = J::array(); std::vector<std::string> hashes;
            for (int rep = 0; rep < 3; ++rep) {
                auto owner = ctx_owner(gpu_model);
                hg::GpuMotionProvider push(device, owner->motion_context(), census.size());
                auto p = hc::prepare_closed_proposal(*owner, o, census, binding_for(*owner), hc::sha256("gpu/source"), push);
                hashes.push_back(p.candidate_field_sha256());
                double max_rel = 0; const auto& a = p.field_heat().candidate_field(); const auto& b = ref.field_heat().candidate_field();
                for (int e = 0; e < g.ne; ++e) { const double d = std::abs((a.E[e] + a.E_low[e]) - (b.E[e] + b.E_low[e])); max_rel = std::max(max_rel, d); }
                J r = guard_receipt(p); r["timings"] = timings(p.timings()); r["owner_preflight"] = owner->preflight(p, census);
                r["max_abs_E_difference_vs_cpu_V_m"] = max_rel;
                runs.push_back(r);
            }
            bool ok = hashes[0] == hashes[1] && hashes[1] == hashes[2];
            for (const auto& r : runs) ok = ok && r["owner_preflight"].get<bool>() && r["residual_over_target"].get<double>() <= 1;
            R.record("closed-loop", "GPU push + GPU field at z=0: accepted through the guard receipt (independent residual, Gauss, divB, energy), "
                                    "preflight passes, bit-identical over three runs", ok, runs);
        }
        report["solver_stats"] = {{"solves", hg::gpu_field_solver_stats(*gpu_solver).solves}, {"total_s", hg::gpu_field_solver_stats(*gpu_solver).total_s}};
        report["memory"] = J::parse(device->memory_json());
        gpu_solver.reset(); gpu_model.reset();
        report["teardown"] = J::parse(device->close_json());
    } catch (const std::exception& e) { report["error"] = e.what(); std::cerr << "fatal: " << e.what() << '\n'; rc = 2; }
    report["cases"] = R.cases; report["failed"] = R.failed;
    if (!out.empty()) std::ofstream(out / "gpu.json") << report.dump(1) << '\n';
    return rc ? rc : (R.failed ? 1 : 0);
}
