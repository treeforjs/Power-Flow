// Phase profile of Hall v3 (sealed sources, read in place) for the field-only
// (vacuum) step and the closed coupled step.
//   profile_v3_baseline vacuum  TRACE OUTDIR STEPS [gpu NVRTC UUID]
//   profile_v3_baseline coupled TRACE OUTDIR STEPS CARRIERS [gpu NVRTC UUID]
// Every step runs the owner protocol of record (prepare_vacuum or
// prepare_closed_proposal, preflight_*_final, publish). Every PROBE-th step the
// same inputs are also run through the phase-instrumented copy of
// Model::prepare (profile_prepare.hpp), whose candidate is compared word for
// word with the candidate of record, and the owner-level operations (admit,
// digests, conductivity, thermal prepare, preflight parts) are timed one by one.
#include "root_api.hpp"
#include "hall-finite-foil-em-20261029/axial_particles.hpp"
#ifdef PROFILE_GPU
#include "gpu/gpu_backends.hpp"
#endif
namespace hc = hall_corrective_v3;
#include "profile_prepare.hpp"
#include "motion_fixtures.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>

namespace n = raikou::mfem_backend;
namespace mf = motion_fixtures;
using J = hc::Json;
using steady = std::chrono::steady_clock;

namespace {
double ms_since(steady::time_point t) { return std::chrono::duration<double, std::milli>(steady::now() - t).count(); }
struct Stat {
    std::vector<double> v;
    void add(double x) { v.push_back(x); }
    J json() const {
        if (v.empty()) return J::object();
        auto s = v; std::sort(s.begin(), s.end());
        double sum = 0; for (double x : s) sum += x;
        return {{"n", s.size()}, {"median_ms", s[s.size() / 2]}, {"mean_ms", sum / s.size()}, {"min_ms", s.front()}, {"max_ms", s.back()}};
    }
};
struct Stats {
    std::map<std::string, Stat> m;
    void add(const std::string& k, double x) { m[k].add(x); }
    J json() const { J j = J::object(); for (const auto& [k, s] : m) j[k] = s.json(); return j; }
};
bool same(const hc::FieldCandidate& a, const hc::FieldCandidate& b) {
    const auto& x = a.field(); const auto& y = b.field();
    auto eqw = [](const hc::WordVec& p, const hc::WordVec& q) {
        if (p.size() != q.size()) return false;
        for (std::size_t i = 0; i < p.size(); ++i) if (p[i].hi != q[i].hi || p[i].lo != q[i].lo) return false;
        return true; };
    return x.E == y.E && x.B == y.B && x.E_low == y.E_low && x.B_low == y.B_low && x.time_s == y.time_s && x.revision == y.revision &&
           eqw(a.charge().conductor_C, b.charge().conductor_C) && eqw(a.charge().exterior_C, b.charge().exterior_C) &&
           eqw(a.charge().particles.net_C, b.charge().particles.net_C) && a.charge().particles.gross_C == b.charge().particles.gross_C &&
           a.receipt().cell_joule_W == b.receipt().cell_joule_W && a.midpoint().seal() == b.midpoint().seal();
}
struct Boundary { std::string field_sha256, heat_sha256, source_sha256, mesh_sha256, material_sha256, common_clock_sha256; double time_s = 0; std::uint64_t accepted_step = 0; };
hc::ChargeLedger declared(const hc::Model& m, const axial::State& s, hc::FreeCharge particles) {
    const auto& g = m.mesh(); hc::ChargeLedger q; q.particles = std::move(particles);
    q.conductor_C.assign(g.nn, hc::Words()); q.exterior_C.assign(g.nn, hc::Words());
    q.conductor_declaration = hc::sha256("profile: declared DC conductor/PEC induced charge");
    q.exterior_declaration = hc::sha256("profile: declared DC exterior termination charge");
    const auto field = hc::field_charge(g, s);
    for (int k = 0; k < g.nn; ++k) { if (g.node_body[k] >= 0 || g.short_node[k]) q.conductor_C[k] = field[k]; else if (g.outer_node[k]) q.exterior_C[k] = field[k]; }
    return q;
}
hall_driver::MotionBinding binding_for(const hc::Owner& o, long step) {
    hall_driver::MotionBinding b; b.phase = hall_driver::ParticlePhase::full_time_position_and_proper_velocity;
    b.accepted_step = o.revision(); b.staged_arena_revision = o.revision() + 1; b.begin_s = o.time_s(); b.dt_s = 1e-14; b.end_s = b.begin_s + b.dt_s;
    b.staged_arena_sha256 = hc::sha256("profile/staged/" + std::to_string(step)); b.mesh_sha256 = o.motion_mesh_sha256();
    b.material_sha256 = o.model().material_sha256(); b.species_sha256 = hc::sha256("profile/species");
    b.common_clock_sha256 = hc::sha256("profile/clock/" + std::to_string(step)); b.kinetic_policy_sha256 = o.motion_context().policy_sha256;
    return b;
}
// Owner-level operations of one field-only step, timed one by one on the
// accepted state and the candidate of record.
void probe_owner(const hc::Owner& owner, const hc::FieldHeatProposal& p, Stats& S) {
    const auto& model = owner.model();
    auto t = steady::now(); auto sigma = owner.thermal().conductivity(); S.add("owner.conductivity()", ms_since(t));
    t = steady::now(); (void)hc::Model::conductivity_sha256(sigma); S.add("owner.conductivity_sha256", ms_since(t));
    t = steady::now(); (void)model.admit(owner.field(), owner.charge()); S.add("owner.admit(one state)", ms_since(t));
    t = steady::now(); (void)model.state_sha256(owner.field(), owner.charge()); S.add("owner.state_sha256(one state)", ms_since(t));
    t = steady::now(); (void)hc::heat_sha256_of(model.material_sha256(), owner.thermal()); S.add("owner.heat_sha256(accepted)", ms_since(t));
    t = steady::now(); (void)hc::heat_sha256_of(model.material_sha256(), p.heat()); S.add("owner.heat_sha256(trial)", ms_since(t));
    t = steady::now(); (void)p.candidate_ledgers().sha256(); S.add("owner.ledgers_sha256", ms_since(t));
    t = steady::now();
    {
        const auto& f = p.candidate_field(); const auto& g = owner.mesh(); const auto& a = owner.field();
        hc::Vec E(g.ne), B(g.nb);
        for (int e = 0; e < g.ne; ++e) E[e] = ((hc::Words(a.E[e], a.E_low[e]) + hc::Words(f.E[e], f.E_low[e])) * .5).rounded();
        for (int b = 0; b < g.nb; ++b) B[b] = ((hc::Words(a.B[b], a.B_low[b]) + hc::Words(f.B[b], f.B_low[b])) * .5).rounded();
        S.add("preflight.midpoint_recompute", ms_since(t)); t = steady::now();
        (void)hc::Midpoint(std::move(E), std::move(B)).seal();
        S.add("preflight.midpoint_seal(hash)", ms_since(t));
    }
    t = steady::now(); (void)hall_driver::total_field_sha256(p.midpoint().E(), p.midpoint().B()); S.add("coupled.total_field_sha256(midpoint)", ms_since(t));
    t = steady::now(); (void)model.at(owner.field()); S.add("coupled.model.at(accepted)->Midpoint", ms_since(t));
    t = steady::now(); (void)owner.preflight_field_heat(p); S.add("preflight_field_heat(total)", ms_since(t));
}
// Thermal prepare on a separate ThermalV3 with the receipt of record.
void probe_thermal(hc::ThermalV3& th, const hc::FieldHeatProposal& p, Stats& S) {
    auto t = steady::now(); auto trial = th.prepare(p.field_receipt(), 1e-14, p.candidate_field_sha256());
    S.add("heat.thermal_prepare(8 slices)", ms_since(t));
    t = steady::now(); (void)th.preflight(trial); S.add("heat.thermal_preflight", ms_since(t));
}

std::shared_ptr<hc::FieldSolveBackend> backend(int argc, char** argv, int at, std::shared_ptr<void>& keep, J& meta) {
#ifdef PROFILE_GPU
    if (argc > at + 2 && std::string(argv[at]) == "gpu") {
        auto dev = std::make_shared<hc::gpu::Device>(hc::gpu::DeviceOptions{argv[at + 2], argv[at + 1]});
        keep = dev; meta["device"] = J::parse(dev->device_json());
        return hc::gpu::make_gpu_field_solver(dev);
    }
#else
    (void)argc; (void)argv; (void)at; (void)keep; (void)meta;
#endif
    return hc::make_cpu_pcg();
}

int vacuum(const std::string& trace, const std::filesystem::path& out, long steps, int argc, char** argv) {
    J meta = J::object(); std::shared_ptr<void> keep;
    auto solver = backend(argc, argv, 5, keep, meta);
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, solver);
    const auto& g = model->mesh(); const auto o = mf::geometry(g); const auto seed = model->dc_seed();
    std::vector<n::ArenaRateEntry> none;
    hc::Owner owner(model, seed, declared(*model, seed, hc::charge_of_census(g, o, none)), o, none);
    hc::ThermalV3 th(g, model->drive().start());
    Stats S; long probes = 0, probe_equal = 0;
    const long probe_every = std::max(1L, steps / 20);
    const auto t_all = steady::now();
    for (long step = 1; step <= steps; ++step) {
        const bool probe = step % probe_every == 0 || step == 1;
        auto t = steady::now();
        auto f = owner.prepare_vacuum(1e-14);
        S.add("step.prepare_vacuum", ms_since(t)); t = steady::now();
        Boundary bd{f.candidate_field_sha256(), f.candidate_heat_sha256(), "", "", "", "", f.end_s(), f.endpoint_revision()};
        const bool ok = owner.preflight_vacuum_final(f, bd);
        S.add("step.preflight_vacuum_final", ms_since(t));
        hc::need(ok, "vacuum preflight");
        if (probe) {
            ++probes;
            const auto sigma = owner.thermal().conductivity(); const auto sigma_sha = hc::Model::conductivity_sha256(sigma);
            hc::PhaseClock T;
            auto c = hc::ModelTestAccess::prepare(*model, T, owner.field(), owner.charge(), owner.field_sha256(), hc::Vec(g.ne), owner.charge().particles,
                                                  hc::DepositSupport::empty(g), sigma, sigma_sha, 1e-14);
            for (const auto& [k, v] : T.ms) S.add(k, v);
            auto t2 = steady::now(); (void)model->state_sha256(c.field(), c.charge()); S.add("heat_for.seal(state_sha256)", ms_since(t2));
            t2 = steady::now(); (void)model->admit(c.field(), c.charge()); S.add("heat_for.admit_candidate", ms_since(t2));
            probe_equal += (c.field().E == f.candidate_field().E && c.field().B == f.candidate_field().B && c.field().E_low == f.candidate_field().E_low &&
                            c.field().B_low == f.candidate_field().B_low && c.receipt().cell_joule_W == f.field_receipt().cell_joule_W &&
                            c.midpoint().seal() == f.midpoint().seal()) ? 1 : 0;
            probe_owner(owner, f, S);
            probe_thermal(th, f, S);
        }
        t = steady::now();
        owner.publish_vacuum_prevalidated(std::move(f));
        S.add("step.publish", ms_since(t));
        S.add("step.total(prepare+preflight+publish)", S.m["step.prepare_vacuum"].v.back() + S.m["step.preflight_vacuum_final"].v.back() + S.m["step.publish"].v.back());
    }
    const double wall = ms_since(t_all);
    J summary{{"schema", "hall-v3-vacuum-profile-v1"}, {"backend", solver->name()}, {"steps", steps}, {"wall_ms", wall},
              {"probes", probes}, {"probe_candidate_equal_to_record", probe_equal}, {"end_time_s", owner.time_s()},
              {"mesh", {{"ne", g.ne}, {"nb", g.nb}, {"nn", g.nn}, {"cells", g.cells.size()}}}, {"phases", S.json()}, {"meta", meta}};
    std::filesystem::create_directories(out);
    std::ofstream(out / "vacuum-profile.json") << summary.dump(1) << '\n';
    std::cout << summary.dump(1) << std::endl;
    return probe_equal == probes ? 0 : 3;
}

int coupled(const std::string& trace, const std::filesystem::path& out, long steps, std::size_t carriers, int argc, char** argv) {
    J meta = J::object(); std::shared_ptr<void> keep;
    auto solver = backend(argc, argv, 6, keep, meta);
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, solver, hc::CarrierLimits{std::max<std::size_t>(100000, carriers), 1ull << 34});
    const auto& g = model->mesh(); const auto o = mf::geometry(g); const auto seed = model->dc_seed();
    std::vector<n::ArenaRateEntry> census = carriers == 2 ? mf::pair({g.side / 2, g.side / 2, 0.}, 1.) : mf::pairs(carriers, 7, 1e3);
    hc::Owner owner(model, seed, declared(*model, seed, hc::charge_of_census(g, o, census, model->limits().maximum_carriers)), o, census);
    hc::ThermalV3 th(g, model->drive().start());
    std::unique_ptr<hc::MotionProvider> provider = std::make_unique<hc::CpuMotionProvider>();
#ifdef PROFILE_GPU
    if (keep) provider = std::make_unique<hc::gpu::GpuMotionProvider>(std::static_pointer_cast<hc::gpu::Device>(keep), owner.motion_context(), census.size());
#endif
    Stats S; long probes = 0, probe_equal = 0; const long probe_every = std::max(1L, steps / 10);
    const std::string src = hc::sha256("profile/source");
    for (long step = 1; step <= steps; ++step) {
        const auto h0 = hc::motion_host_timers();
        auto t = steady::now();
        auto p = hc::prepare_closed_proposal(owner, o, census, binding_for(owner, step), src, *provider, hc::StepInputs{});
        S.add("step.prepare_closed_proposal", ms_since(t));
        const auto& T = p.timings();
        S.add("proposal.admit_census", 1e3 * T.admit_census_s); S.add("proposal.motion(all iterates)", 1e3 * T.motion_s);
        S.add("proposal.deposit_support(all iterates)", 1e3 * T.support_s); S.add("proposal.field_prepare(all iterates)", 1e3 * T.field_prepare_s);
        S.add("proposal.solve(all iterates)", 1e3 * T.solve_s); S.add("proposal.residual_check(all iterates)", 1e3 * T.residual_check_s);
        S.add("proposal.motion_v2_json_seal", 1e3 * T.motion_v2_seal_s); S.add("proposal.field_seal", 1e3 * T.seal_s); S.add("proposal.heat_for", 1e3 * T.heat_s);
        S.add("proposal.iterations", T.iterations);
        const auto h1 = hc::motion_host_timers();
        S.add("motion_host.census_hash", 1e3 * (h1.census_hash_s - h0.census_hash_s)); S.add("motion_host.v3_digest", 1e3 * (h1.record_v3_s - h0.record_v3_s));
        const auto& b = p.motion().binding();
        Boundary bd{p.candidate_field_sha256(), p.candidate_heat_sha256(), src, b.mesh_sha256, b.material_sha256, b.common_clock_sha256, p.end_s(), p.endpoint_revision()};
        t = steady::now(); const bool ok = owner.preflight_final(p, census, bd); S.add("step.preflight_final", ms_since(t));
        hc::need(ok, "closed preflight");
        if (step % probe_every == 0 || step == 1) {
            ++probes;
            const auto sigma = owner.thermal().conductivity(); const auto sigma_sha = hc::Model::conductivity_sha256(sigma);
            const auto bank = hc::bank_transfer(g, o, p.motion(), model->exterior_support());
            auto t2 = steady::now();
            const auto support = hc::DepositSupportBuilder::build(g, o, census, p.motion());
            S.add("coupled.deposit_support(one iterate)", ms_since(t2));
            hc::PhaseClock PT;
            auto c = hc::ModelTestAccess::prepare(*model, PT, owner.field(), owner.charge(), owner.field_sha256(), p.motion().current_A_m(),
                                                  p.field_heat().candidate_charge().particles, support, sigma, sigma_sha, 1e-14, &bank);
            for (const auto& [k, v] : PT.ms) S.add(k, v);
            probe_equal += (c.field().E == p.field_heat().candidate_field().E && c.field().B == p.field_heat().candidate_field().B &&
                            c.field().E_low == p.field_heat().candidate_field().E_low && c.field().B_low == p.field_heat().candidate_field().B_low &&
                            c.midpoint().seal() == p.midpoint().seal()) ? 1 : 0;
            probe_owner(owner, p.field_heat(), S);
            probe_thermal(th, p.field_heat(), S);
            t2 = steady::now(); (void)hc::canonical_census_sha256(census, model->limits().maximum_carriers); S.add("coupled.canonical_census_sha256", ms_since(t2));
            t2 = steady::now(); (void)hc::MotionAssembly::record_v3(p.motion()); S.add("coupled.record_v3", ms_since(t2));
            t2 = steady::now(); (void)hc::charge_of_endpoints(g, o, p.motion().endpoints(), model->limits().maximum_carriers); S.add("coupled.charge_of_endpoints", ms_since(t2));
            t2 = steady::now(); (void)hc::charge_of_census(g, o, census, model->limits().maximum_carriers); S.add("coupled.charge_of_census(admit_census)", ms_since(t2));
        }
        auto next = hc::next_census(census, p.motion());
        t = steady::now(); owner.publish_prevalidated(std::move(p)); S.add("step.publish", ms_since(t));
        census = std::move(next);
    }
    J summary{{"schema", "hall-v3-coupled-profile-v1"}, {"backend", solver->name()}, {"motion_provider", provider->name()}, {"steps", steps},
              {"carriers", carriers}, {"probes", probes}, {"probe_candidate_equal_to_record", probe_equal}, {"phases", S.json()}, {"meta", meta}};
    std::filesystem::create_directories(out);
    std::ofstream(out / ("coupled-profile-" + std::to_string(carriers) + ".json")) << summary.dump(1) << '\n';
    std::cout << summary.dump(1) << std::endl;
    return probe_equal == probes ? 0 : 3;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc >= 5 && std::string(argv[1]) == "vacuum") return vacuum(argv[2], argv[3], std::stol(argv[4]), argc, argv);
        if (argc >= 6 && std::string(argv[1]) == "coupled") return coupled(argv[2], argv[3], std::stol(argv[4]), std::stoul(argv[5]), argc, argv);
        std::cerr << "usage: profile_v3_baseline vacuum TRACE OUT STEPS [gpu NVRTC UUID] | coupled TRACE OUT STEPS CARRIERS [gpu NVRTC UUID]\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }
}
