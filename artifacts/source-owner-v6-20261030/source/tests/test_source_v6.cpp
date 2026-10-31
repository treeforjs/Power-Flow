// Heat/source owner v6 bound to the Hall corrective v3 thermal owner (H5: one thermal owner).
//   test_source_v6 TRACE OUTDIR
// Groups: binding to the Hall slice view; Hall temperature consumed by emission and
// desorption; source loads applied exactly once by Hall (watermarks, closure);
// double-Joule detection (sealed source v5 private engine vs v6); typed refusals of
// the heat-engine extremes (1e200 stock, ~5e21 release); publication order and
// staleness; v6 component checkpoint bound to the Hall clock.
#include "root_api.hpp"
#include "hall-finite-foil-em-20261029/axial_surface_path_port.hpp"
#include "ss304_heat_material_v2.hpp"
#include "water_scenario.hpp"
#include "pic-electrode-heat-source-v5-20261030/native/mfem_electrode_heat_source_v5/source_owner.hpp"
#include "pic-electrode-heat-source-v6-20261030/native/mfem_electrode_heat_source_v6/source_owner.hpp"
#include "hall_source_binding.hpp"
#include "pf_arrival/arrival_owner.hpp"
#include "pf_arrival_v11/arrival_owner.hpp"
#include "pf_arrival_v11/source_v6_adapter.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <optional>

namespace hc = hall_corrective_v3;
namespace hb = hall_source_binding;
namespace s5 = electrode_heat_v5;
namespace s6 = electrode_heat_v6;
namespace eh = electrode_heat_v2;
namespace v3 = electrode_heat_v3;
namespace se = surface_emission;
namespace n = raikou::mfem_backend;
namespace pa1 = pf_arrival;
namespace pa = pf_arrival_v11;
namespace ps = pf_arrival_v11::source_v6;
using J = hc::Json;

namespace {
constexpr double dt = 1e-14;
struct Results {
    J cases = J::array(); int failed = 0;
    void record(const std::string& g, const std::string& name, bool ok, J ev = J::object()) {
        cases.push_back({{"group", g}, {"name", name}, {"passed", ok}, {"evidence", std::move(ev)}});
        if (!ok) ++failed;
        std::cout << g << " | " << name << ": " << (ok ? "pass" : "FAIL") << std::endl;
        if (!ok) std::cout << "  evidence: " << cases.back()["evidence"].dump().substr(0, 3000) << std::endl;
    }
};
template <class F> std::string refusal(F&& f) { try { f(); } catch (const std::exception& e) { return e.what(); } return ""; }
n::GapPicOptions geometry(const axial::Mesh& g) {
    n::GapPicOptions o; o.nx = o.ny = g.nx; o.width_m = o.gap_m = g.side; o.depth_m = g.depth;
    o.finite_foils = o.finite_axial_extent = true; o.short_z_m = -g.depth / 2; o.feed_z_m = g.depth / 2;
    o.lower_foil = {g.side / 2 - .0025, g.side / 2 + .0025, g.side / 2 - .0027, g.side / 2 - .0025};
    o.upper_foil = {g.side / 2 - .0025, g.side / 2 + .0025, g.side / 2 + .0025, g.side / 2 + .0027};
    return o;
}
hc::ChargeLedger declared(const hc::Model& m, const axial::State& s) {
    const auto& g = m.mesh(); hc::ChargeLedger q; q.particles = {hc::WordVec(g.nn), hc::Vec(g.nn)};
    q.conductor_C.assign(g.nn, hc::Words()); q.exterior_C.assign(g.nn, hc::Words());
    q.conductor_declaration = hc::sha256("source v6 test: declared DC conductor/PEC induced charge");
    q.exterior_declaration = hc::sha256("source v6 test: declared DC exterior termination charge");
    const auto field = hc::field_charge(g, s);
    for (int k = 0; k < g.nn; ++k) { if (g.node_body[k] >= 0 || g.short_node[k]) q.conductor_C[k] = field[k]; else if (g.outer_node[k]) q.exterior_C[k] = field[k]; }
    return q;
}
// Lower foil, inner face (field frame y = -2.5 mm: slice geometry y0 = -2.7 mm, HighY), source hash of the scenario.
const std::string source_hash = "source-v6-hall-slice-scenario";
eh::Patch metal_patch() { eh::Patch p; p.id = "lower-foil-inner-face"; p.face = eh::Face::HighY; p.source_hash = source_hash; return p; }
eh::Patch water(double N, double rate, double desorption_eV = .87) {
    auto p = metal_patch(); p.desorption_enabled = true; p.state_id = p.parent_id = "H2O"; p.initial_particles = N; p.initial_low_particles = 0;
    p.mass_kg = 18 * 1.66053906660e-27; p.attempt_s = rate; p.activation_eV = 0; p.desorption_eV = desorption_eV; p.internal_eV = .2; return p;
}
se::LayerPolicy layer_policy() { se::LayerPolicy p; p.local_linear_vacuum_extent_m = 25e-6; return p; }
hb::SliceSource slice_source(const axial::Mesh& g, eh::Patch p, int k) {
    hb::SliceSource s; s.body = 0; s.slice = k; s.patches = {p}; s.source_hash = source_hash;
    s.identities = {{"lower-foil-inner-surface-slice-" + std::to_string(k), "conditional-single-layer", .005 * g.dz}};
    s.settings = {{se::Config{4.3}, layer_policy(), "source-v6-test-sampling"}};
    return s;
}
v3::TrialFrame frame(const s6::Owner& o, std::uint64_t stage) {
    v3::TrialFrame f; f.role = v3::FieldRole::TotalMidpoint; f.accepted_revision = o.thermal_state().revision; f.trial_stage = stage;
    f.mesh_revision = 1; f.material_revision = 2; f.begin_s = o.thermal_state().time; f.dt_s = dt; f.end_s = f.begin_s + dt;
    f.field_stage_s = f.begin_s + .5 * dt; f.arena_hash = o.state_sha256(); f.field_hash = "hall-v3-authoritative-midpoint";
    f.mesh_hash = "hall-v3-axial-mesh"; f.material_hash = "global-material-v2-map"; f.species_hash = "actual-source-v1-electron";
    return f;
}
struct Arrays { foil::Vec E, B; };
Arrays arrays(const axial::Mesh& g, double F) { Arrays a{foil::Vec(g.ne), foil::Vec(g.nb)}; for (int e = 0; e < g.ne; ++e) if (g.edges[e].axis == 1) a.E[e] = F; return a; }
v3::PatchBinding binding(const s6::Owner& o, const axial::Mesh& g, int k) {
    v3::PatchBinding b; b.thermal_revision = o.thermal_state().revision; b.source_accepted_state = b.thermal_revision + 1; b.patch_id = 7;
    b.thermal_mesh_identity = hb::slice_geometry(g, 0, k).revision; b.thermal_material_hash = eh::material_headers_sha256;
    b.material_phase_admitted = b.planar_barrier_admitted = b.electron_supply_admitted = true; b.temperature = v3::HeldTemperature::AcceptedPatchMean;
    b.extrusion_center_z_m = hb::slice_center_z(g, k); b.extrusion_height_m = g.dz; return b;
}
v3::BirthSample sample(const s6::PreparedPatch& p, const s6::Owner& o, double w) {
    v3::BirthSample b; b.physical_weight = w; b.normal_uniform = .6; b.transverse_uniform = .35; b.azimuth_uniform = .1; b.time_uniform = .23; b.event_id = 123;
    b.canonical_particle_id = 42; b.count_owner_revision = o.thermal_state().revision; b.count_owner_base_hash = "external-accepted-count-base";
    b.count_owner_token = "external-staged-ticket"; b.exposure_key = p.exposure_key(); b.time_policy = v3::BirthTime::UniformHeldExposure;
    b.weight_provenance = v3::WeightProvenance::StagedCanonicalOwner; return b;
}
eh::SourceEventIdentity event(std::uint64_t seq) { return {"source-v6-test", "e-", "electrode-electron-reservoir", "native-layer", "canonical-electron-42", seq, seq}; }
std::vector<std::vector<s6::NeutralLeaf>> leaves(const std::vector<s6::ReleaseQuote>& q) {
    std::vector<std::vector<s6::NeutralLeaf>> out(q.size());
    for (std::size_t j = 0; j < q.size(); ++j) if (!q[j].released.zero()) for (double w : q[j].released.positive_carriers()) out[j].push_back({w, eh::Draw{.5, .5, .5, .5, .5, .5}});
    return out;
}
// Every accepted Active source layer must be advanced in each interval (source v5/v6 rule).
std::vector<s6::Candidate> advance_active(s6::Owner& src, const axial::Mesh& g, int k, const eh::Epoch& e, const std::shared_ptr<const s6::Field>& fld, std::uint64_t seq) {
    std::vector<s6::Candidate> out;
    if (src.active_layers().empty()) return out;
    const v3::Exposure ex{fld->frame().begin_s, dt, 1, "physical-quadrature"};
    auto pp = src.prepare_patch(e, 0, binding(src, g, k), fld, {.0001, hb::slice_geometry(g, 0, k).thickness + hb::slice_geometry(g, 0, k).y0, hb::slice_center_z(g, k)}, ex);
    for (const auto& [id, rec] : src.active_layers()) out.push_back(src.advance_active(pp, id, event(seq)));
    return out;
}
// Conditional oxide-water patch (heat-v2 water_scenario.hpp recipe: Temkin coverage,
// Eyring kT/h prefactor, 0.87 eV (1 - 0.5 theta) barrier and binding; 1.3e19 m^-2
// monolayer) on the lower-foil inner face of one Hall slice, initial coverage 0.33.
eh::Patch oxide_water(const eh::Geometry& gg) {
    auto p = eh::oxide_water_patch(gg, "lower-foil-inner-face", eh::Face::HighY, 0, 1, .33, 18 * 1.66053906660e-27, "H2O", "H2O", source_hash,
                                   eh::WaterRange::ConditionalAnalytic298To3000);
    return p;
}
struct Boundary { std::string field_sha256, heat_sha256, source_sha256, mesh_sha256, material_sha256, common_clock_sha256; double time_s = 0; std::uint64_t accepted_step = 0; };
// One Hall field-only step carrying loads; preflight; publish Hall then (optionally) the source trial.
bool publish(hc::Owner& hall, hc::FieldHeatProposal&& f, s6::Owner* src, std::optional<s6::Trial>* trial, const s6::FinalFieldView* view) {
    Boundary bd{f.candidate_field_sha256(), f.candidate_heat_sha256(), "", "", "", "", f.end_s(), f.endpoint_revision()};
    bool ok = hall.preflight_vacuum_final(f, bd);
    if (src && trial && *trial) ok = ok && src->preflight(**trial, *view);
    if (!ok) return false;
    hall.publish_vacuum_prevalidated(std::move(f));
    if (src && trial && *trial) { src->publish_prevalidated(std::move(**trial)); trial->reset(); }
    return true;
}
// Carrier growth of quantized neutral release (source v6 + arrival v1.1) over many
// steps at a fixed electrode temperature: one Hall-slice-sized oxide-water patch
// (lower-foil inner face, 5 mm x 6.25 mm, coverage 0.33), the one thermal owner a
// heat-v2 slice Coupler at T (desorption loads applied to it each step), dt = 10 fs.
J growth(double T, long steps, double w0) {
    const axial::Mesh g(45, 4, .009, 1);
    const int k = 2;
    const auto gg = hb::slice_geometry(g, 0, k);
    const auto patch = oxide_water(gg);
    eh::Coupler view(gg, hb::slice_material(), hb::heat_only({patch}), T, source_hash, 0.);
    pa::Configuration cfg; cfg.run_seed = 4242; cfg.streams = {{1, pa::Kind::neutral, 0, "H2O", w0, 4096}};
    auto P = std::make_shared<const pa::Policy>(cfg);
    pa1::Configuration cfg1; cfg1.run_seed = 4242; cfg1.streams = {{1, pa1::Kind::neutral, 0, "H2O", w0, 4096}};
    const pa1::Policy P1(cfg1);
    std::vector<s6::SourceSettings> settings{{se::Config{4.3}, layer_policy(), P->sampling_policy_identity(), w0}};
    s6::Owner src(gg, hb::slice_material(), {patch}, {{"lower-foil-inner-surface", "oxide-water-conditional", .005 * g.dz}}, settings, T, source_hash, 0., view);
    pa::Owner A(P, src.thermal_state().revision);
    const auto stock0 = src.inventory(0);
    mfem::ConstantCoefficient zero(0);
    const auto zeroE = arrays(g, 0);
    long double sumD = 0, sumx = 0, var = 0; std::uint64_t carriers = 0, carriers_v1_floor_rule = 0, steps_positive = 0;
    n::ExactParticleCount released;
    const auto t0 = std::chrono::steady_clock::now();
    J samples = J::array();
    for (long st = 0; st < steps; ++st) {
        const auto e = src.make_epoch(dt, "growth-isothermal-midpoint");
        const auto q = src.release_quotes(e);
        const double D = q[0].released.zero() ? 0 : q[0].released.rounded_for_rate();
        auto in = ps::neutral_input(*P, 1, q[0], src, e);
        auto plan = A.plan({src.thermal_state().revision, {}, {in}});
        auto staged = A.stage(plan, 1 + 100000 * std::uint64_t(st));
        const auto lv = ps::neutral_leaves(staged, 1);
        if (D > 0) { ++steps_positive; const double x = D / w0, f = x - std::floor(x); sumD += D; sumx += x; var += f * (1 - f);
            pa1::NeutralRelease r1; r1.stream_id = 1; r1.debit = q[0].released; r1.log_debit = q[0].log_released; r1.inventory = src.inventory(0);
            r1.begin_s = e.t; r1.end_s = e.t + e.dt; carriers_v1_floor_rule += pa1::decide(P1, std::uint64_t(st), r1).births; }
        carriers += lv[0].size();
        for (const auto& l : lv[0]) released = released.credited(l.physical_weight);
        v3::TrialFrame f; f.role = v3::FieldRole::TotalMidpoint; f.accepted_revision = src.thermal_state().revision; f.trial_stage = 1; f.mesh_revision = 1;
        f.material_revision = 2; f.begin_s = src.thermal_state().time; f.dt_s = dt; f.end_s = f.begin_s + dt; f.field_stage_s = f.begin_s + .5 * dt;
        f.arena_hash = "growth-arena"; f.field_hash = "growth-field"; f.mesh_hash = "growth-mesh"; f.material_hash = "global-material-v2-map"; f.species_hash = "actual-source-v1-electron";
        auto fld = std::make_shared<const s6::Field>(g, f, zeroE.E, zeroE.B);
        const s6::FinalFieldView fv{f, g, zeroE.E, zeroE.B};
        auto t = src.prepare(e, fld, fv, {}, {}, lv);
        staged.bind_source_times(1, ps::neutral_times(t, staged, 1, patch.id));
        auto heat = view.prepare(view.make_epoch(dt, "growth-isothermal-midpoint"), zero, t.loads(), std::vector<std::vector<eh::Draw>>(1));
        hc::need(src.preflight(t, fv) && view.preflight(heat), "growth step preflight");
        hc::need(view.commit(std::move(heat)), "growth thermal commit");
        src.publish_prevalidated(std::move(t));
        pa::NativeBirthReceipt r; r.base = staged.native_base(); std::uint64_t id = r.base;
        for (const auto& b : staged.births()) r.entries.push_back({b.packet_id(), "H2O", id++, b.weight});
        r.next_id_after = id;
        (void)A.commit(std::move(staged), r);
        if ((st + 1) % (steps / 10 > 0 ? steps / 10 : 1) == 0)
            samples.push_back({{"step", st + 1}, {"carriers", carriers}, {"carriers_v1_floor_rule", carriers_v1_floor_rule}, {"released_molecules", released.zero() ? 0. : released.rounded_for_rate()},
                               {"expected_molecules", double(sumD)}, {"surface_T_K", view.surface_temperatures().at(0)}});
    }
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    n::ExactParticleCount check = src.inventory(0); for (double w : released.positive_carriers()) check = check.credited(w);
    const double z = var > 0 ? double(((long double)carriers - sumx) / std::sqrt(var)) : 0.;
    return {{"T_K", T}, {"steps", steps}, {"simulated_ns", steps * dt * 1e9}, {"w0_molecules", w0}, {"initial_stock_molecules", stock0.rounded_for_rate()},
            {"expected_release_per_step_mean", double(sumD / std::max<long>(1, steps))}, {"expected_molecules", double(sumD)},
            {"released_molecules", released.zero() ? 0. : released.rounded_for_rate()}, {"carriers_created_v1_1", carriers},
            {"expected_carriers_sum_x", double(sumx)}, {"standardized_deviation", z}, {"carriers_v1_floor_rule_same_quotes", carriers_v1_floor_rule},
            {"steps_with_positive_quote", steps_positive}, {"carriers_per_ns_v1_1", carriers / (steps * dt * 1e9)},
            {"carriers_per_ns_v1_floor_rule", carriers_v1_floor_rule / (steps * dt * 1e9)}, {"stock_plus_released_equals_initial_exactly", check == stock0},
            {"released_equals_carriers_times_w0_exactly", carriers == 0 ? released.zero() : released == n::ExactParticleCount(double(carriers) * w0)},
            {"final_surface_T_K", view.surface_temperatures().at(0)}, {"wall_s", wall}, {"ms_per_step", 1e3 * wall / steps}, {"samples", samples}};
}
}  // namespace

int main(int argc, char** argv) {
    if (argc == 7 && std::string(argv[3]) == "growth") {
        try {
            const auto r = growth(std::stod(argv[4]), std::stol(argv[5]), std::stod(argv[6]));
            std::filesystem::create_directories(argv[2]);
            std::ofstream(std::filesystem::path(argv[2]) / ("growth-" + std::string(argv[4]) + "K.json")) << r.dump(1) << '\n';
            std::cout << r.dump(1) << std::endl;
            return 0;
        } catch (const std::exception& e) { std::cerr << "fatal: " << e.what() << '\n'; return 2; }
    }
    Results R; J report = J::object(); int rc = 0;
    try {
        if (argc != 3) throw std::runtime_error("usage: test_source_v6 TRACE OUTDIR [growth T_K STEPS W0]");
        const std::filesystem::path out = argv[2]; std::filesystem::create_directories(out);
        const int k = 2;   // lower foil, second foil slice
        auto base_model = std::make_shared<const hc::Model>(hc::MeshSpec{}, argv[1]);
        const auto& g0 = base_model->mesh();
        std::vector<hb::SliceSource> sources{slice_source(g0, water(1e15, 1e8), k)};
        auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, argv[1], hc::make_cpu_pcg(), hc::CarrierLimits{}, hb::declaration(sources));
        const auto& g = model->mesh(); const auto o = geometry(g); const auto seed = model->dc_seed();
        auto mesh_ptr = std::make_shared<const axial::Mesh>(g);
        report["hall_configuration_sha256"] = model->configuration_sha256(); report["hall_material_sha256"] = model->material_sha256();

        // ---- 1. binding
        hc::Owner hall(model, seed, declared(*model, seed), o, {});
        auto src = hb::bind(hall, sources[0]);
        {
            J ev;
            ev["wrong_slice_view"] = refusal([&] { auto s = sources[0]; s.slice = k + 1; (void)hb::bind(hall, s); });
            hc::Owner plain(base_model, seed, declared(*base_model, seed), o, {});
            ev["undeclared_hall_slice"] = refusal([&] { (void)hb::bind(plain, sources[0]); });
            ev["source_configuration_sha256"] = src->configuration_sha256();
            R.record("binding", "source v6 binds the declared heat-only Hall slice; another slice or an undeclared Hall slice is refused",
                     ev["wrong_slice_view"].get<std::string>().find("v6 thermal view") != std::string::npos &&
                         !ev["undeclared_hall_slice"].get<std::string>().empty(), ev);
        }
        // ---- 2. Hall heats the slice; the source reads the Hall temperature (desorption quote and emission exposure)
        {
            hc::SurfaceLoad heat; heat.body = 0; heat.slice = k; heat.channel = "test-preheat";
            auto& x = heat.transfer; x.producer = "test-preheat"; x.energy_owner = "test-preheat/energy"; x.sequence = x.owner_sequence = 1;
            x.state_id = x.parent_id = "electrode"; x.channel_id = "test-preheat"; x.event_time_s = hall.time_s() + .5 * dt; x.surface = true;
            x.x = .0; x.y = -.0025; x.incoming_K_J = 2e-3;
            std::vector<hc::SurfaceLoad> L{heat};
            auto f = hall.prepare_vacuum(dt, L);
            const bool pub = publish(hall, std::move(f), nullptr, nullptr, nullptr);
            const auto e = src->make_epoch(dt, "hall-v3-authoritative-midpoint");
            const auto q = src->release_quotes(e);
            const double hallT = hall.thermal().slice(0, k).surface_temperatures().at(0);
            auto fld = std::make_shared<const s6::Field>(g, frame(*src, 1), arrays(g, -3e9).E, arrays(g, -3e9).B);
            const s6::PatchIdentity id = src->patch_identity(0); (void)id;
            const v3::Exposure ex{fld->frame().begin_s, dt, 1, "physical-quadrature"};
            auto pp = src->prepare_patch(e, 0, binding(*src, g, k), fld, {.0001, hb::slice_geometry(g, 0, k).thickness + hb::slice_geometry(g, 0, k).y0, hb::slice_center_z(g, k)}, ex);
            J ev{{"hall_patch_mean_T_K", hallT}, {"source_quote_T_K", q[0].T_K}, {"source_exposure_T_K", pp.temperature_K()}, {"published", pub},
                 {"source_thermal_revision", src->thermal_state().revision}, {"hall_revision", hall.revision()}};
            R.record("temperature", "after Hall heats the slice, desorption and emission use the Hall temperature bit for bit (no private engine)",
                     pub && hallT > 300 && q[0].T_K == hallT && pp.temperature_K() == hallT && src->thermal_state().revision == hall.revision(), ev);
        }
        // ---- 3. desorption + emission loads applied exactly once by Hall
        {
            const auto e = src->make_epoch(dt, "hall-v3-authoritative-midpoint");
            auto fld = std::make_shared<const s6::Field>(g, frame(*src, 2), arrays(g, -3e9).E, arrays(g, -3e9).B);
            auto a = arrays(g, -3e9);
            const s6::FinalFieldView view{fld->frame(), g, a.E, a.B};
            const v3::Exposure ex{fld->frame().begin_s, dt, 1, "physical-quadrature"};
            auto pp = src->prepare_patch(e, 0, binding(*src, g, k), fld, {.0001, hb::slice_geometry(g, 0, k).thickness + hb::slice_geometry(g, 0, k).y0, hb::slice_center_z(g, k)}, ex);
            auto birth = src->begin_source(pp, sample(pp, *src, 1e9), event(1));
            std::vector<s6::Candidate> cands{birth};
            std::vector<s6::Handoff> hands;
            if (birth.state().native.status == se::LayerStatus::PICHandedOff) {
                s6::Particle canonical; canonical.id = 42; canonical.mass_kg = se::me; canonical.charge_C = -se::qe; canonical.weight = birth.state().native.particle.weight;
                hands.push_back(src->handoff(birth, view, canonical));
            } else if (birth.state().native.status == se::LayerStatus::Active) cands = {src->advance_birth(pp, birth, event(2))};
            const auto q = src->release_quotes(e);
            std::optional<s6::Trial> trial; trial.emplace(src->prepare(e, fld, view, cands, hands, leaves(q)));
            const auto L = hb::loads(*trial, sources[0]);
            double sum = 0; J channels = J::object();
            for (const auto& l : L) { sum += hc::transfer_energy_J(l.transfer); channels[l.channel] = channels.value(l.channel, 0.) + hc::transfer_energy_J(l.transfer); }
            auto f = hall.prepare_vacuum(dt, L);
            const auto& h = f.heat();
            J ev{{"loads", L.size()}, {"load_energy_J", sum}, {"channels", channels}, {"hall_thermal_delta_J", h.thermal_delta_J}, {"hall_joule_J", h.admitted_joule_J},
                 {"hall_external_J", h.external_J}, {"source_budget_joule_J", trial->heat_budget().joule_J}, {"source_budget_wall_net_J", trial->heat_budget().wall_net_J},
                 {"layer_status", int(birth.state().native.status)}, {"released_H2O", q[0].released.rounded_for_rate()}};
            const bool exact_once = std::abs(h.thermal_delta_J - h.admitted_joule_J - h.external_J) <= h.roundoff_bound_J + 2e-12 * (h.external_heating_J + h.external_cooling_J) &&
                                    std::abs(h.external_J - sum) <= 1e-15 * std::abs(sum);
            std::string audit = refusal([&] { std::vector<hb::OwnerReport> others{{"source-v6", trial->heat_budget().joule_J, trial->heat_budget().thermal_delta_J}};
                                              hb::audit_single_thermal_owner(h, others); });
            const auto L_copy = L;
            const bool pub = publish(hall, std::move(f), src.get(), &trial, &view);
            ev["published"] = pub; ev["audit"] = audit;
            ev["replayed_loads_refused"] = refusal([&] { (void)hall.prepare_vacuum(dt, L_copy); });
            ev["hall_ledger_loads"] = hall.kinetic_ledgers().loads;
            R.record("exactly-once", "emission layer and desorption loads of source v6 enter Hall once: thermal delta = Joule + loads, source reports no Joule, "
                                     "the single-owner audit passes, replaying the same loads next step is refused",
                     L.size() >= 2 && exact_once && audit.empty() && pub && !trial.has_value() &&
                         ev["replayed_loads_refused"].get<std::string>().find("duplicate/stale") != std::string::npos, ev);
        }
        // ---- 4. double-Joule detection: a sealed source v5 private engine driven by a Joule term is refused by the audit
        {
            eh::Geometry gg = hb::slice_geometry(g, 0, k);
            auto p5 = water(1e15, 1e8); p5.source_hash = source_hash;
            s5::Owner v5(gg, hb::slice_material(), {p5}, {{"lower-foil-inner-surface-slice-" + std::to_string(k), "conditional-single-layer", .005 * g.dz}},
                         {{se::Config{4.3}, layer_policy(), "source-v5-private-engine"}}, 300., source_hash, model->drive().start());
            v3::TrialFrame fr; fr.role = v3::FieldRole::TotalMidpoint; fr.accepted_revision = v5.thermal_state().revision; fr.trial_stage = 1; fr.mesh_revision = 1;
            fr.material_revision = 2; fr.begin_s = v5.thermal_state().time; fr.dt_s = dt; fr.end_s = fr.begin_s + dt; fr.field_stage_s = fr.begin_s + .5 * dt;
            fr.arena_hash = v5.state_sha256(); fr.field_hash = "f"; fr.mesh_hash = "m"; fr.material_hash = "global-material-v2-map"; fr.species_hash = "actual-source-v1-electron";
            auto a = arrays(g, 0); auto f5 = std::make_shared<const s5::Field>(g, fr, a.E, a.B); const s5::FinalFieldView view5{fr, g, a.E, a.B};
            mfem::ConstantCoefficient joule(1e15);   // W/m3: the Joule power Hall also integrates
            const auto e5 = v5.make_epoch(dt, "hall-midpoint");
            std::vector<std::vector<s5::NeutralLeaf>> l5(1);
            for (const auto& q : v5.release_quotes(e5)) if (!q.released.zero()) for (double w : q.released.positive_carriers()) l5[0].push_back({w, eh::Draw{.5, .5, .5, .5, .5, .5}});
            auto t5 = v5.prepare(e5, joule, f5, view5, {}, {}, l5);
            auto f = hall.prepare_vacuum(dt);
            J ev{{"v5_private_joule_J", t5.heat_budget().joule_J}, {"v5_private_thermal_delta_J", t5.heat_budget().thermal_delta_J}, {"hall_joule_J", f.heat().admitted_joule_J}};
            ev["audit_v5"] = refusal([&] { std::vector<hb::OwnerReport> others{{"source-v5-private-engine", t5.heat_budget().joule_J, t5.heat_budget().thermal_delta_J}};
                                           hb::audit_single_thermal_owner(f.heat(), others); });
            ev["audit_v6"] = refusal([&] { std::vector<hb::OwnerReport> others{{"source-v6", 0., 0.}}; hb::audit_single_thermal_owner(f.heat(), others); });
            ev["hall_external_joule_load"] = refusal([&] { hc::SurfaceLoad l; l.body = 0; l.slice = k; l.channel = "source-v5-joule"; auto& x = l.transfer;
                x.producer = "source-v5-private-heat"; x.energy_owner = "source-v5-private-heat/energy"; x.sequence = x.owner_sequence = 1; x.state_id = x.parent_id = "electrode";
                x.channel_id = "joule"; x.event_time_s = hall.time_s(); x.surface = true; x.x = 0; x.y = -.0025; x.authored_heat_J = t5.heat_budget().joule_J;
                std::vector<hc::SurfaceLoad> L{l}; (void)hall.prepare_vacuum(dt, L); });
            R.record("double-joule", "a source v5 private engine with its own Joule term is detected (audit refuses typed); v6 passes; an external Joule load is refused by Hall",
                     t5.heat_budget().joule_J > 0 && ev["audit_v5"].get<std::string>().find("[double-joule]") != std::string::npos && ev["audit_v6"] == "" &&
                         ev["hall_external_joule_load"].get<std::string>().find("[external-joule]") != std::string::npos, ev);
        }
        // ---- 5. heat-engine extremes are typed scoped refusals (no abort)
        {
            J ev = J::object(); bool ok = true;
            for (auto [label, N, rate] : std::vector<std::tuple<std::string, double, double>>{{"stock_1e200", 1e200, 1e300}, {"release_5e21", 5e21, 1e300}, {"release_1e20", 1e20, 1e300}}) {
                std::vector<hb::SliceSource> sx{slice_source(g0, water(N, rate), k)};
                auto mx = std::make_shared<const hc::Model>(hc::MeshSpec{}, argv[1], hc::make_cpu_pcg(), hc::CarrierLimits{}, hb::declaration(sx));
                hc::Owner hx(mx, seed, declared(*mx, seed), o, {});
                auto sxo = hb::bind(hx, sx[0]);
                J c;
                try {
                    const auto e = sxo->make_epoch(dt, "hall-v3-authoritative-midpoint");
                    auto fld = std::make_shared<const s6::Field>(g, frame(*sxo, 1), arrays(g, 0).E, arrays(g, 0).B);
                    auto a = arrays(g, 0); const s6::FinalFieldView view{fld->frame(), g, a.E, a.B};
                    const auto q = sxo->release_quotes(e);
                    c["released"] = q[0].released.rounded_for_rate();
                    auto t = sxo->prepare(e, fld, view, {}, {}, leaves(q));
                    const auto L = hb::loads(t, sx[0]);
                    double sum = 0; for (const auto& l : L) sum += hc::transfer_energy_J(l.transfer);
                    c["load_J"] = sum;
                    try { auto f = hx.prepare_vacuum(dt, L); c["admitted"] = true; c["thermal_delta_J"] = f.heat().thermal_delta_J; }
                    catch (const hc::ThermalRefusal& r) { c["typed"] = {{"kind", r.kind}, {"quantity", r.quantity}, {"identity", r.identity}, {"observed", r.observed}, {"bound", r.bound}}; }
                } catch (const s6::Refusal& r) { c["source_typed"] = {{"kind", s6::kind_name(r.kind)}, {"quantity", r.quantity}, {"log_scope", r.log_abs_scope}}; }
                catch (const std::exception& x) { c["untyped"] = x.what(); }
                c["hall_unchanged"] = hx.revision() == 0;
                ev[label] = c;
                if (label != "release_1e20") ok = ok && (c.contains("typed") || c.contains("source_typed")) && c["hall_unchanged"].get<bool>();
                else ok = ok && (c.contains("typed") || c.value("admitted", false));
            }
            // The sealed v5 engine with the same ~5e21 release refuses untyped (recorded; the 1e200 case aborts the v5 process and is not run here).
            {
                eh::Geometry gg = hb::slice_geometry(g, 0, k); auto p5 = water(5e21, 1e300);
                s5::Owner v5(gg, hb::slice_material(), {p5}, {{"s", "l", .005 * g.dz}}, {{se::Config{4.3}, layer_policy(), "v5"}}, 300., source_hash, model->drive().start());
                v3::TrialFrame fr; fr.role = v3::FieldRole::TotalMidpoint; fr.accepted_revision = 0; fr.trial_stage = 1; fr.mesh_revision = 1; fr.material_revision = 2;
                fr.begin_s = v5.thermal_state().time; fr.dt_s = dt; fr.end_s = fr.begin_s + dt; fr.field_stage_s = fr.begin_s + .5 * dt; fr.arena_hash = v5.state_sha256();
                fr.field_hash = "f"; fr.mesh_hash = "m"; fr.material_hash = "global-material-v2-map"; fr.species_hash = "actual-source-v1-electron";
                auto a = arrays(g, 0); auto f5 = std::make_shared<const s5::Field>(g, fr, a.E, a.B); const s5::FinalFieldView view5{fr, g, a.E, a.B};
                mfem::ConstantCoefficient zero(0); const auto e5 = v5.make_epoch(dt, "m");
                std::vector<std::vector<s5::NeutralLeaf>> l5(1);
                for (const auto& q : v5.release_quotes(e5)) for (double w : q.released.positive_carriers()) l5[0].push_back({w, eh::Draw{.5, .5, .5, .5, .5, .5}});
                ev["sealed_v5_release_5e21"] = refusal([&] { (void)v5.prepare(e5, zero, f5, view5, {}, {}, l5); });
            }
            R.record("typed-refusals", "1e200 stock and ~5e21 release: source v6 loads refused typed by Hall before any engine call (sealed v5: untyped engine error / abort)",
                     ok, ev);
        }
        // ---- 6. v6 component checkpoint (no heat section) restores against the same Hall clock only
        {
            const auto file = out / "source-v6.bin"; std::filesystem::remove(file); src->save(file.string());
            auto fresh = hb::bind(hall, sources[0]);
            const auto r1 = refusal([&] { fresh->load(file.string()); });
            J ev{{"restore", r1}, {"same_state", r1.empty() && fresh->state_sha256() == src->state_sha256()}};
            {   // one joint step: Hall with the source loads, the source trial advancing its Active layers
                const auto e = src->make_epoch(dt, "hall-v3-authoritative-midpoint");
                auto fld = std::make_shared<const s6::Field>(g, frame(*src, 4), arrays(g, 0).E, arrays(g, 0).B);
                auto a = arrays(g, 0); const s6::FinalFieldView view{fld->frame(), g, a.E, a.B};
                const auto active = advance_active(*src, g, k, e, fld, 20);
                std::optional<s6::Trial> t; t.emplace(src->prepare(e, fld, view, active, {}, leaves(src->release_quotes(e))));
                auto L = hb::loads(*t, sources[0]);
                auto f = hall.prepare_vacuum(dt, L);
                ev["joint_step_published"] = publish(hall, std::move(f), src.get(), &t, &view);
            }
            auto late = hb::bind(hall, sources[0]);
            ev["restore_after_hall_advanced"] = refusal([&] { late->load(file.string()); });
            R.record("checkpoint", "v6 checkpoint restores inventory/layers bound to the Hall clock; a different Hall clock is refused",
                     ev["same_state"].get<bool>() && ev["joint_step_published"].get<bool>() && ev["restore_after_hall_advanced"].get<std::string>().find("thermal clock") != std::string::npos, ev);
        }
        // ---- 7. staleness: a source trial prepared before a Hall publication is refused afterwards
        {
            const auto e = src->make_epoch(dt, "hall-v3-authoritative-midpoint");
            auto fld = std::make_shared<const s6::Field>(g, frame(*src, 3), arrays(g, 0).E, arrays(g, 0).B);
            auto a = arrays(g, 0); const s6::FinalFieldView view{fld->frame(), g, a.E, a.B};
            const auto active = advance_active(*src, g, k, e, fld, 30);
            auto t = src->prepare(e, fld, view, active, {}, leaves(src->release_quotes(e)));
            const bool before = src->preflight(t, view);
            auto f = hall.prepare_vacuum(dt);
            const bool pub = publish(hall, std::move(f), nullptr, nullptr, nullptr);
            const bool after = src->preflight(t, view);
            R.record("protocol", "source v6 trial preflights against the Hall thermal clock: valid before, stale after a Hall publication",
                     before && pub && !after, {{"before", before}, {"after", after}, {"active_layers", src->active_layers().size()}});
        }
        // ---- 8. quantized neutral release (source v6 mode + arrival policy v1.1)
        {
            const double w0 = 1e6;
            auto sq = slice_source(g0, water(1e15, 1e8), k); sq.settings[0].quantized_neutral_weight = w0;
            pa::Configuration cfg; cfg.run_seed = 20261030; cfg.streams = {{1, pa::Kind::neutral, 0, "H2O", w0, 4096}};
            auto P = std::make_shared<const pa::Policy>(cfg);
            sq.settings[0].sampling_policy_identity = P->sampling_policy_identity();
            std::vector<hb::SliceSource> sx{sq};
            auto mq = std::make_shared<const hc::Model>(hc::MeshSpec{}, argv[1], hc::make_cpu_pcg(), hc::CarrierLimits{}, hb::declaration(sx));
            hc::Owner hq(mq, seed, declared(*mq, seed), o, {});
            auto q6 = hb::bind(hq, sq);
            pa::Owner A(P, q6->thermal_state().revision);
            const auto stock0 = q6->inventory(0);
            n::ExactParticleCount released; J steps = J::array(); bool ok = true; double expected = 0; std::uint64_t leaves_total = 0;
            for (int st = 0; st < 6 && ok; ++st) {
                const auto e = q6->make_epoch(dt, "hall-v3-authoritative-midpoint");
                auto fld = std::make_shared<const s6::Field>(g, frame(*q6, st + 1), arrays(g, 0).E, arrays(g, 0).B);
                auto a = arrays(g, 0); const s6::FinalFieldView view{fld->frame(), g, a.E, a.B};
                const auto q = q6->release_quotes(e);
                auto plan = A.plan({q6->thermal_state().revision, {}, {ps::neutral_input(*P, 1, q[0], *q6, e)}});
                auto staged = A.stage(plan, 1000 + 10000 * std::uint64_t(st));
                const auto lv = ps::neutral_leaves(staged, 1);
                std::optional<s6::Trial> t; t.emplace(q6->prepare(e, fld, view, {}, {}, lv));
                staged.bind_source_times(1, ps::neutral_times(*t, staged, 1, sq.patches[0].id));
                n::ExactParticleCount debit; for (const auto& l : lv[0]) debit = debit.credited(l.physical_weight);
                n::ExactParticleCount pk; bool all_w0 = true; for (const auto& x : t->packets()) { pk = pk.credited(x.weight); all_w0 = all_w0 && x.weight == w0; }
                n::ExactParticleCount expect_inv = q6->inventory(0); for (double w : debit.positive_carriers()) expect_inv = expect_inv.debit(w);
                const bool exact = t->inventory(0) == expect_inv && pk == debit && all_w0 && lv[0].size() == staged.births().size();
                expected += q[0].released.rounded_for_rate(); leaves_total += lv[0].size();
                for (double w : debit.positive_carriers()) released = released.credited(w);
                steps.push_back({{"quote_D", q[0].released.rounded_for_rate()}, {"leaves", lv[0].size()}, {"debit_exact_carriers", hc::exact_json(debit)},
                                 {"exact_debit_equals_credit", exact}});
                ok = ok && exact;
                auto L = hb::loads(*t, sq);
                auto f = hq.prepare_vacuum(dt, L);
                ok = ok && publish(hq, std::move(f), q6.get(), &t, &view);
                pa::NativeBirthReceipt r; r.base = staged.native_base(); std::uint64_t id = r.base;
                for (const auto& b : staged.births()) r.entries.push_back({b.packet_id(), "H2O", id++, b.weight});
                r.next_id_after = id;
                (void)A.commit(std::move(staged), r);
            }
            n::ExactParticleCount stock_check = q6->inventory(0); for (double w : released.positive_carriers()) stock_check = stock_check.credited(w);
            J ev{{"w0", w0}, {"steps", steps}, {"leaves_total", leaves_total}, {"expected_total_D", expected},
                 {"stock_conserved_exactly", stock_check == stock0}};
            ev["non_w0_leaf_refused"] = refusal([&] { const auto e = q6->make_epoch(dt, "m"); auto fld = std::make_shared<const s6::Field>(g, frame(*q6, 50), arrays(g, 0).E, arrays(g, 0).B);
                auto a = arrays(g, 0); const s6::FinalFieldView view{fld->frame(), g, a.E, a.B};
                std::vector<std::vector<s6::NeutralLeaf>> bad{{{w0 * 1.5, eh::Draw{.5, .5, .5, .5, .5, .5}}}}; (void)q6->prepare(e, fld, view, {}, {}, bad); });
            R.record("quantized-release", "quantized release: leaves of exactly w0, exact source debit = arrival credit = packet sum, stock + released = initial "
                                          "stock exactly, a non-w0 leaf is refused", ok && ev["stock_conserved_exactly"].get<bool>() &&
                     ev["non_w0_leaf_refused"].get<std::string>().find("exactly w0") != std::string::npos, ev);
        }
        // ---- 9. unbiasedness of the keyed stochastic rounding of v1.1 (pure decisions)
        {
            pa::Configuration cfg; cfg.run_seed = 77; cfg.streams = {{3, pa::Kind::neutral, 0, "H2O", 1., 4096}};
            pa::Policy P(cfg);
            J rows = J::array(); bool ok = true;
            for (double x : {1e-3, .3, 1.7, 12.5}) {
                const long N = 200000; long double sum = 0; const double f = x - std::floor(x);
                n::ExactParticleCount inv(1e300);
                for (long st = 0; st < N; ++st) { pa::NeutralRelease r; r.stream_id = 3; r.debit = n::ExactParticleCount(x); r.log_debit = std::log(x);
                    r.inventory = inv; r.begin_s = 0; r.end_s = 1e-14; sum += pa::decide(P, std::uint64_t(st), r).births; }
                const double z = double((sum - (long double)N * x) / std::sqrt(double(N) * f * (1 - f)));
                rows.push_back({{"x", x}, {"decisions", N}, {"mean_n", double(sum / N)}, {"z", z}});
                ok = ok && std::abs(z) < 5;
            }
            R.record("unbiasedness", "v1.1 keyed stochastic rounding: mean n equals D/w0 within 5 standard deviations (2e5 decisions per x)", ok, rows);
        }
    } catch (const std::exception& e) { report["error"] = e.what(); std::cerr << "fatal: " << e.what() << '\n'; rc = 2; }
    report["cases"] = R.cases; report["failed"] = R.failed;
    std::ofstream(std::filesystem::path(argv[argc - 1]) / "source-v6.json") << report.dump(1) << '\n';
    return rc ? rc : (R.failed ? 1 : 0);
}
