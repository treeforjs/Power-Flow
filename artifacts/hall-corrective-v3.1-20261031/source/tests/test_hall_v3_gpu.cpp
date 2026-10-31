// Hall corrective v3 GPU qualification (one device session, one stream).
//   test_hall_v3_gpu TRACE NVRTC UUID OUTDIR
// Writes OUTDIR/gpu.json, OUTDIR/sink-dump-gpu.json (same rows as the CPU sink dump)
// and OUTDIR/gpu-motion-dump.json (v2 fixture set, same rows as the CPU motion dump).
#include "root_api.hpp"
#include "gpu/gpu_backends.hpp"
namespace hc = hall_corrective_v3;
#include "parity_scenarios.hpp"
#include "sink_fixtures.hpp"
#include "sink_dump.hpp"
#include <chrono>
#include <fstream>
#include <iostream>

namespace hg = hall_corrective_v3::gpu;
namespace hd = hall_driver;
namespace n = raikou::mfem_backend;
namespace mf = motion_fixtures;
namespace sf = sink_fixtures;
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
std::string after_digest(const hc::FreeCharge& q) { hc::Digest d("test-after-charge"); hc::encode(d, q.net_C); d.f64s(q.gross_C); return d.hex(); }
}  // namespace

int main(int argc, char** argv) {
    int rc = 0; J report = J::object(); Results R; std::filesystem::path out;
    try {
        if (argc != 5) throw std::runtime_error("usage: test_hall_v3_gpu TRACE NVRTC UUID OUTDIR");
        out = argv[4]; std::filesystem::create_directories(out);
        auto device = std::make_shared<hg::Device>(hg::DeviceOptions{argv[3], argv[2]});
        report["device"] = J::parse(device->device_json());
        auto cpu_model = std::make_shared<const hc::Model>(hc::MeshSpec{}, argv[1]);
        const auto& g = cpu_model->mesh(); const auto o = mf::geometry(g); const auto seed = cpu_model->dc_seed();
        const auto ctx = hc::make_motion_context(g, o, cpu_model->mesh_json());
        // ---- 1. sink fixtures: GPU provider vs CPU provider (records, sink records, arrays, free charge, closed proposals)
        {
            hg::GpuMotionProvider gpu(device, ctx, 10000);
            report["push_compile"] = J::parse(gpu.compile_json());
            const auto t = steady::now();
            auto rows = sink_dump::run(cpu_model, gpu);
            report["sink_dump_gpu_s"] = since(t);
            std::ofstream(out / "sink-dump-gpu.json") << rows.dump(1) << '\n';
            hc::CpuMotionProvider cpu;
            auto ref = sink_dump::run(cpu_model, cpu);
            bool same = rows.size() == ref.size(); J diffs = J::array(); long sinks = 0;
            for (std::size_t i = 0; i < rows.size() && i < ref.size(); ++i) {
                const bool eq = rows[i] == ref[i];
                same = same && eq; sinks += rows[i].value("sinks", 0);
                diffs.push_back({{"fixture", rows[i]["fixture"]}, {"identical", eq}, {"sinks", rows[i].value("sinks", 0)},
                                 {"motion_refusal", rows[i].value("motion_refusal", "")}});
            }
            R.record("sink-parity", "GPU sink records, truncated paths, current, endpoint charge, free-charge words and closed proposals equal the CPU "
                                    "provider on every sink fixture (bitwise)", same && sinks > 0, {{"fixtures", diffs}, {"sinks_total", sinks}});
            auto again = sink_dump::run(cpu_model, gpu);
            R.record("sink-parity", "GPU sink dump is bit-identical when repeated in the same process", again == rows);
        }
        // ---- 2. v2 fixture set: unchanged records on the v3 device path
        {
            hg::GpuMotionProvider gpu(device, ctx, 10000); hc::CpuMotionProvider cpu;
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
                    J r{{"fixture", f.name}, {"carriers", f.census.size()}, {"entry", "MotionAssembly::cpu"}};
                    try {
                        p.begin(ctx, f.census); hc::FreeCharge after;
                        auto m = p.motion(b, E, mid.B(), after); p.end();
                        hc::MotionAssembly::seal_v2(m);
                        J rec = hd::motion_record(m); const auto text = rec.dump(); hc::drain(rec);
                        r["record_sha256"] = hc::sha256(text); r["seal_sha256"] = m.seal_sha256(); r["record_v3_sha256"] = m.record_v3_sha256();
                        r["after_sha256"] = after_digest(after);
                    } catch (const std::exception& e) { p.end(); r["refusal"] = e.what(); }
                    if (which) dump.push_back(r);
                    pair[which ? "gpu" : "cpu"] = r;
                }
                const bool eq = pair["gpu"] == pair["cpu"];
                same = same && eq; diffs.push_back({{"fixture", f.name}, {"identical", eq}});
            }
            std::ofstream(out / "gpu-motion-dump.json") << dump.dump(1) << '\n';
            R.record("v2-fixtures", "v2 fixture set: GPU v3 provider record bytes, seals, free charge and refusal texts equal the CPU v3 provider", same, diffs);
        }
        // ---- 3. GPU push + GPU field: a births/collection step and a short multi-step loop, closures checked
        {
            auto gpu_solver = hg::make_gpu_field_solver(device);
            auto gpu_model = std::make_shared<const hc::Model>(hc::MeshSpec{}, argv[1], gpu_solver);
            auto census = sf::with_partners({sf::electron(1, {.0045, .007 - 3e-8, .001}, 100, {0, 1, 0}, 1e6),
                                             sf::electron(2, {.0045, .002 + 3e-8, .002}, 100, {0, -1, 0}, 1e6)});
            hc::Owner owner(gpu_model, seed, parity::declared(*gpu_model, seed, hc::charge_of_census(g, o, census)), o, census);
            hg::GpuMotionProvider push(device, owner.motion_context(), 64);
            J steps = J::array(); bool ok = true;
            for (int s = 1; s <= 6; ++s) {
                std::vector<hc::Newborn> births;
                if (s <= 3) { auto e = sf::electron(100 + s, {.003 + .001 * s, sf::lower_top, 0.}, 2e6, {0, 1, 0}, 1e6);   // on the lower-foil face, leaving it
                              census.push_back(e); births.push_back({e.native.id, owner.time_s() + .3 * dt}); }
                J r{{"step", s}};
                try {
                    auto p = hc::prepare_closed_proposal(owner, o, census, parity::binding(owner, s), hc::sha256("gpu/source"), push, hc::StepInputs{births, {}});
                    const auto& c = p.field_heat().candidate_constraints();
                    r["sinks"] = p.sinks().size(); r["births"] = p.births().size(); r["energy_residual_J"] = p.energy().residual_J; r["energy_scale_J"] = p.energy().scale_J;
                    r["heat_roundoff_bound_J"] = p.field_heat().heat().roundoff_bound_J;
                    r["max_gauss_C"] = c.max_gauss_C; r["max_gauss_bound_C"] = c.max_gauss_bound_C;
                    r["residual_over_target"] = p.field_heat().solve_record().independent_residual / p.field_heat().solve_record().target;
                    const auto& b = p.motion().binding();
                    parity::Boundary bd{p.candidate_field_sha256(), p.candidate_heat_sha256(), hc::sha256("gpu/source"), b.mesh_sha256, b.material_sha256,
                                        b.common_clock_sha256, p.end_s(), p.endpoint_revision()};
                    const bool pre = owner.preflight_final(p, census, bd); r["preflight"] = pre;
                    auto next = hc::next_census(census, p.motion());
                    if (pre) owner.publish_prevalidated(std::move(p));
                    census = std::move(next);
                    ok = ok && pre && std::abs(r["energy_residual_J"].get<double>()) <= 4e-20 + 2e-10 * r["energy_scale_J"].get<double>() + r["heat_roundoff_bound_J"].get<double>();
                } catch (const std::exception& e) { r["refusal"] = e.what(); ok = false; }
                steps.push_back(r);
                if (!ok) break;
            }
            R.record("closed-loop", "GPU push + GPU field: births, collection at both foils, absolute Gauss, energy closure, preflight and publication over 6 steps",
                     ok && owner.kinetic_ledgers().sinks.count() >= 2, {{"steps", steps}, {"sink_ledger", owner.kinetic_ledgers().sinks.json()}});
            report["solver_stats"] = {{"solves", hg::gpu_field_solver_stats(*gpu_solver).solves}, {"total_s", hg::gpu_field_solver_stats(*gpu_solver).total_s}};
        }
        report["memory"] = J::parse(device->memory_json());
        report["teardown"] = J::parse(device->close_json());
    } catch (const std::exception& e) { report["error"] = e.what(); std::cerr << "fatal: " << e.what() << '\n'; rc = 2; }
    report["cases"] = R.cases; report["failed"] = R.failed;
    if (!out.empty()) std::ofstream(out / "gpu.json") << report.dump(1) << '\n';
    return rc ? rc : (R.failed ? 1 : 0);
}
