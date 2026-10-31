// Hall corrective v2 component qualification (CPU backends).
//   test_hall_v2 components TRACE OUTDIR    identity, node-plane guard, audit A09/A10 ports,
//                                           backend negative controls, owner protocol
//   test_hall_v2 motion-dump TRACE OUTFILE  motion record digests of the shared fixture set
//                                           (compared with the sealed-provider reference)
#include "root_api.hpp"
#include "motion_fixtures.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <optional>

namespace hc = hall_corrective_v2;
namespace hd = hall_driver;
namespace n = raikou::mfem_backend;
namespace mf = motion_fixtures;
using J = hc::Json;
using steady = std::chrono::steady_clock;

namespace hall_corrective_v2 {
struct ModelTestAccess {
    static FieldCandidate prepare(const Model& m, const axial::State& old, const ChargeLedger& before, const Vec& plasma,
                                  const FreeCharge& after, const DepositSupport& support, double dt) {
        const auto sigma = m.initial_conductivity();
        return m.prepare(old, before, m.state_sha256(old, before), plasma, after, support, sigma, Model::conductivity_sha256(sigma), dt);
    }
};
}  // namespace hall_corrective_v2

namespace {
constexpr double dt = 1e-14;
double since(steady::time_point t) { return std::chrono::duration<double>(steady::now() - t).count(); }
struct Results {
    J cases = J::array();
    int failed = 0;
    void record(const std::string& group, const std::string& name, bool ok, J evidence = J::object()) {
        cases.push_back({{"group", group}, {"name", name}, {"passed", ok}, {"evidence", std::move(evidence)}});
        if (!ok) ++failed;
        std::cout << group << " | " << name << ": " << (ok ? "pass" : "FAIL") << std::endl;
    }
};
template <class F> std::string refusal(F&& f) {
    try { f(); } catch (const std::exception& e) { return e.what(); }
    return "";
}
hc::ChargeLedger declared_dc(const hc::Model& m, const axial::State& s, hc::FreeCharge particles, const std::string& tag) {
    const auto& g = m.mesh();
    hc::ChargeLedger q;
    q.particles = std::move(particles);
    q.conductor_C.assign(g.nn, hc::Words()); q.exterior_C.assign(g.nn, hc::Words());
    q.conductor_declaration = hc::sha256(tag + ": declared DC conductor/PEC induced charge");
    q.exterior_declaration = hc::sha256(tag + ": declared DC exterior termination charge");
    const auto field = hc::field_charge(g, s);
    for (int k = 0; k < g.nn; ++k) {
        if (g.node_body[k] >= 0 || g.short_node[k]) q.conductor_C[k] = field[k];
        else if (g.outer_node[k]) q.exterior_C[k] = field[k];
    }
    return q;
}
hd::MotionBinding binding_for(const hc::Owner& owner, const std::string& tag) {
    hd::MotionBinding b;
    b.phase = hd::ParticlePhase::full_time_position_and_proper_velocity;
    b.accepted_step = owner.revision(); b.staged_arena_revision = 1; b.begin_s = owner.time_s(); b.dt_s = dt; b.end_s = b.begin_s + dt;
    b.staged_arena_sha256 = hc::sha256(tag + "/staged"); b.mesh_sha256 = owner.motion_mesh_sha256();
    b.material_sha256 = owner.model().material_sha256(); b.species_sha256 = hc::sha256(tag + "/species");
    b.common_clock_sha256 = hc::sha256(tag + "/clock"); b.kinetic_policy_sha256 = owner.motion_context().policy_sha256;
    return b;
}
struct Boundary {
    std::string field_sha256, heat_sha256, source_sha256, mesh_sha256, material_sha256, common_clock_sha256;
    double time_s = 0; std::uint64_t accepted_step = 0;
};
Boundary boundary_of(const hc::JointProposal& p, const std::string& source) {
    const auto& b = p.motion().binding();
    return {p.candidate_field_sha256(), p.candidate_heat_sha256(), source, b.mesh_sha256, b.material_sha256, b.common_clock_sha256,
            p.end_s(), p.endpoint_revision()};
}
// Same layout as 8dddec0 Model::state_sha256 (JSON text), for the timing comparison only.
std::string sealed_style_json_hash(const std::string& cfg, const axial::State& s, const hc::ChargeLedger& q) {
    J j = J::object();
    j["schema"] = "hall-compensated-field-owner-v1"; j["configuration_sha256"] = cfg;
    j["field"] = hc::state_record(s); j["charge"] = hc::record(q);
    const auto text = j.dump(); hc::drain(j);
    return hc::sha256(text);
}
class ScaledBackend final : public hc::FieldSolveBackend {
    std::string name_, d_; double scale_; bool one_;
public:
    ScaledBackend(double scale, bool one_component) : name_(one_component ? "test-one-component" : "test-scaled"), scale_(scale), one_(one_component) {
        hc::Digest d("test-backend"); d.text(name_); d.f64(scale); d_ = d.hex();
    }
    const std::string& name() const noexcept override { return name_; }
    const std::string& descriptor_sha256() const noexcept override { return d_; }
    hc::Vec solve(const hc::Apply& a, const hc::Vec& b, int maxit, const hc::SolveContext& c, hc::SolveRecord& r) override {
        hc::CpuPcgBackend cpu; auto x = cpu.solve(a, b, maxit, c, r);
        if (one_) { std::size_t k = 0; for (std::size_t i = 0; i < x.size(); ++i) if (std::abs(x[i]) > std::abs(x[k])) k = i; x[k] *= scale_; }
        else for (double& v : x) v *= scale_;
        return x;
    }
};

// ------------------------------------------------------------------ components
int components(const std::string& trace, const std::filesystem::path& out) {
    Results R;
    J report = J::object();
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace);
    const auto& g = model->mesh();
    const auto o = mf::geometry(g);
    const auto seed = model->dc_seed();
    report["configuration_sha256"] = model->configuration_sha256();
    report["backend"] = {{"name", model->solver().name()}, {"descriptor_sha256", model->solver().descriptor_sha256()}};
    report["material_sha256"] = model->material_sha256();
    report["motion_mesh_sha256"] = model->motion_mesh_sha256(o);

    // ---- 1. binary identities
    {
        auto again = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace);
        auto other = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, std::make_shared<ScaledBackend>(1., false));
        R.record("identity", "configuration identity is reproducible and binds the field-solve backend",
                 again->configuration_sha256() == model->configuration_sha256() && other->configuration_sha256() != model->configuration_sha256(),
                 {{"cpu", model->configuration_sha256()}, {"cpu_again", again->configuration_sha256()}, {"other_backend", other->configuration_sha256()}});
        R.record("identity", "motion mesh identity equals the sealed function value (no per-call JSON)",
                 model->motion_mesh_sha256(o) == hd::motion_mesh_sha256(g, o), {{"cached", model->motion_mesh_sha256(o)}});
        const auto ledger = declared_dc(*model, seed, {hc::WordVec(g.nn), hc::Vec(g.nn)}, "identity");
        auto bumped = seed; bumped.E_low[100] = std::nextafter(0., 1.);
        const auto h0 = model->state_sha256(seed, ledger), h1 = model->state_sha256(seed, ledger), h2 = model->state_sha256(bumped, ledger);
        R.record("identity", "field-state digest is deterministic and sensitive to one subnormal low word", h0 == h1 && h0 != h2);
        // Timing: binary digest versus the 8dddec0 JSON-text layout of the same state.
        std::vector<double> tb, tj;
        for (int k = 0; k < 7; ++k) {
            auto t = steady::now(); (void)model->state_sha256(seed, ledger); tb.push_back(since(t));
            t = steady::now(); (void)sealed_style_json_hash(model->configuration_sha256(), seed, ledger); tj.push_back(since(t));
        }
        std::sort(tb.begin(), tb.end()); std::sort(tj.begin(), tj.end());
        report["field_state_hash_timing_s"] = {{"binary_median", tb[3]}, {"json_text_median", tj[3]}, {"speedup", tj[3] / tb[3]},
                                               {"mesh", {{"edges", g.ne}, {"faces", g.nb}, {"nodes", g.nn}}}};
        R.record("identity", "binary field-state digest is faster than the JSON-text digest of the same state", tb[3] < tj[3],
                 report["field_state_hash_timing_s"]);
    }

    // ---- 2. node-plane guard (defect 1 of Stage 1)
    J plane = J::array();
    for (double w : {1., 1e9})
        for (double z : {0., 1e-12, 1e-9, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3}) {
            const auto census = mf::pair({.0045, .0045, z}, w);
            J c{{"z_m", z}, {"weight", w}};
            try {
                hc::Owner owner(model, seed, declared_dc(*model, seed, hc::charge_of_census(g, o, census), "plane"), o, census);
                auto p = hc::prepare_closed_proposal(owner, o, census, binding_for(owner, "plane"), hc::sha256("plane/source"));
                c["admitted"] = true; c["iterations"] = p.iterations(); c["continuity_defect_over_tolerance"] = p.timings().continuity_ratio;
                c["owner_preflight"] = owner.preflight(p, census);
            } catch (const std::exception& e) { c["admitted"] = false; c["refusal"] = e.what(); }
            plane.push_back(c);
        }
    bool all = true; for (const auto& c : plane) all = all && c.at("admitted").get<bool>() && c.value("owner_preflight", false);
    R.record("node-plane", "charge-conserving pair motion on and next to the z=0 node plane is admitted (Stage 1 reproducer cases)", all, plane);
    for (double w : {1., 100., 1e6, -1.}) {
        // w = -1: seed 22 at weight 1, whose carrier 166 (electron 1e-9 m above a z node plane, 8.4e6 m/s) makes the
        // sealed axial::push work-path iteration exceed its 12 iterations already in the DC seed field.
        const int seed22 = w < 0;
        if (w < 0) w = 1;
        const auto census = mf::node_plane_pairs(250, seed22 ? 22 : w == 1 ? 21 : w == 100 ? 24 : 23, w);
        J c{{"carriers", census.size()}, {"weight", w}};
        try {
            hc::Owner owner(model, seed, declared_dc(*model, seed, hc::charge_of_census(g, o, census), "sweep"), o, census);
            auto p = hc::prepare_closed_proposal(owner, o, census, binding_for(owner, "sweep"), hc::sha256("sweep/source"));
            c["admitted"] = true; c["continuity_defect_over_tolerance"] = p.timings().continuity_ratio; c["owner_preflight"] = owner.preflight(p, census);
        } catch (const std::exception& e) { c["admitted"] = false; c["refusal"] = e.what(); }
        c["seed"] = seed22 ? 22 : w == 1 ? 21 : w == 100 ? 24 : 23;
        if (seed22)
            R.record("node-plane", "random node-plane pairs, seed 22: the refusal (if any) is the sealed push work-path limit, not the endpoint guard",
                     (c.at("admitted").get<bool>() && c.value("owner_preflight", false)) ||
                         c.value("refusal", "") == "axial native work path did not converge", c);
        else if (w < 1e6)
            R.record("node-plane", "random pairs on x, y and z node planes (offsets 0..1e-7 m, 1..1000 eV) are admitted",
                     c.at("admitted").get<bool>() && c.value("owner_preflight", false), c);
        else  // heavy carriers: the sealed self-consistent push may refuse; the endpoint guard must not be the reason
            R.record("node-plane", "heavy random node-plane pairs (weight 1e6): any refusal is not the endpoint charge/current guard",
                     (c.at("admitted").get<bool>() && c.value("owner_preflight", false)) ||
                         c.value("refusal", "") != "actual particle endpoint charge/current mismatch", c);
    }
    // Negative controls at z = 0 through the private preparation.
    {
        const auto census = mf::pair({.0045, .0045, 0.}, 1.);
        const auto before = declared_dc(*model, seed, hc::charge_of_census(g, o, census), "negative");
        hc::Owner owner(model, seed, before, o, census);
        auto b = binding_for(owner, "negative");
        const auto mid = model->at(seed);
        b.trial_field_sha256 = hd::total_field_sha256(mid.E(), mid.B());
        auto motion = hd::prepare_canonical_motion(g, o, census, b, {}, mid.E(), mid.B());
        auto after = hc::charge_of_endpoints(g, o, motion.endpoints());
        const auto support = hc::DepositSupportBuilder::build(g, o, census, motion.endpoints());
        J ev;
        ev["control_admitted"] = refusal([&] { (void)hc::ModelTestAccess::prepare(*model, seed, before, motion.current_A_m(), after, support, dt); }).empty();
        ev["without_absolute_term"] = refusal([&] { (void)hc::ModelTestAccess::prepare(*model, seed, before, motion.current_A_m(), after, hc::DepositSupport::empty(g), dt); });
        int worst = 0; for (int k = 0; k < g.nn; ++k) if (support.absolute_C()[k] > support.absolute_C()[worst]) worst = k;
        ev["absolute_term_max_C"] = support.absolute_C()[worst];
        auto invented = after; invented.net_C[worst] = invented.net_C[worst] + hc::Words(n::pic_e); invented.gross_C[worst] += n::pic_e;
        ev["one_elementary_charge_invented"] = refusal([&] { (void)hc::ModelTestAccess::prepare(*model, seed, before, motion.current_A_m(), invented, support, dt); });
        auto doubled = after; doubled.net_C[worst] = doubled.net_C[worst] + hc::Words(4 * support.absolute_C()[worst]);
        doubled.gross_C[worst] += 4 * support.absolute_C()[worst];
        ev["four_times_the_absolute_term_invented"] = refusal([&] { (void)hc::ModelTestAccess::prepare(*model, seed, before, motion.current_A_m(), doubled, support, dt); });
        std::vector<n::GapParticle> shifted(motion.endpoints().begin(), motion.endpoints().end());
        shifted[0].position_m[2] += 1e-9;
        const auto moved = hc::charge_of_endpoints(g, o, shifted);
        ev["endpoint_moved_1nm_without_current"] = refusal([&] { (void)hc::ModelTestAccess::prepare(*model, seed, before, motion.current_A_m(), moved, support, dt); });
        const std::string text = "actual particle endpoint charge/current mismatch";
        R.record("node-plane", "z=0 control admits; without the absolute term the 8dddec0 refusal reproduces; invented charge (1 e; 4x the term) and "
                               "a 1 nm endpoint shift without current refuse",
                 ev["control_admitted"].get<bool>() && ev["without_absolute_term"] == text && ev["one_elementary_charge_invented"] == text &&
                     ev["four_times_the_absolute_term_invented"] == text && ev["endpoint_moved_1nm_without_current"] == text, ev);
    }

    // ---- 3. audit A09/A10 counterexamples against v2
    {
        const axial::V3 audit_point{.00021, -.00017, .001};
        const double aw = 1e9;
        const std::vector<n::ArenaRateEntry> census{mf::row(1, "A09-e", "electron", {audit_point[0] + o.width_m / 2, audit_point[1] + o.gap_m / 2, audit_point[2]},
                                                            {0, 0, 0}, -n::pic_e, n::pic_me, aw)};
        const auto qc = hc::charge_of_census(g, o, census);
        const auto control = declared_dc(*model, seed, {hc::WordVec(g.nn), hc::Vec(g.nn)}, "A09");
        const auto honest = declared_dc(*model, seed, qc, "A09");
        auto absorbing = honest; {const auto Q = hc::field_charge(g, seed);
            for (int k = 0; k < g.nn; ++k) { if (g.node_body[k] >= 0 || g.short_node[k]) absorbing.conductor_C[k] = Q[k] - qc.net_C[k];
                else if (g.outer_node[k]) absorbing.exterior_C[k] = Q[k] - qc.net_C[k]; }}
        J ev;
        ev["control_admit"] = refusal([&] { (void)model->admit(seed, control); });
        ev["control_owner_empty_census"] = refusal([&] { hc::Owner x(model, seed, control, o, {}); });
        ev["admit_dc_declaration"] = refusal([&] { (void)model->admit(seed, honest); });
        ev["admit_absorbing_declaration"] = refusal([&] { (void)model->admit(seed, absorbing); });
        ev["prepare_from_inconsistent_state"] = refusal([&] { (void)hc::ModelTestAccess::prepare(*model, seed, honest, hc::Vec(g.ne), honest.particles, hc::DepositSupport::empty(g), dt); });
        ev["owner_matching_ledger"] = refusal([&] { hc::Owner x(model, seed, honest, o, census); });
        ev["owner_field_only_ledger"] = refusal([&] { hc::Owner x(model, seed, control, o, census); });
        ev["owner_empty_census_charged_ledger"] = refusal([&] { hc::Owner x(model, seed, honest, o, {}); });
        bool ok = ev["control_admit"] == "" && ev["control_owner_empty_census"] == "";
        for (const char* k : {"admit_dc_declaration", "admit_absorbing_declaration", "prepare_from_inconsistent_state", "owner_matching_ledger",
                              "owner_field_only_ledger", "owner_empty_census_charged_ledger"}) ok = ok && !ev[k].get<std::string>().empty();
        R.record("audit-A09", "manufactured -e x 1e9 inconsistent charge/field is rejected at admit, prepare and bootstrap; controls admit", ok, ev);
        const int face = g.bx(22, 22, g.kshort + 1);
        auto bad = seed; bad.B[face] += .01;
        J dv;
        dv["admit"] = refusal([&] { (void)model->admit(bad, control); });
        dv["prepare"] = refusal([&] { (void)hc::ModelTestAccess::prepare(*model, bad, control, hc::Vec(g.ne), control.particles, hc::DepositSupport::empty(g), dt); });
        dv["bootstrap"] = refusal([&] { hc::Owner x(model, bad, control, o, {}); });
        {
            const auto saved = out / "A09-divB-checkpoint", edited = out / "A09-divB-checkpoint-perturbed";
            std::filesystem::remove_all(saved); std::filesystem::remove_all(edited);
            hc::Owner live(model, seed, control, o, {}); live.save(saved);
            dv["restore_unmodified_control"] = refusal([&] { (void)hc::Owner::restore(model, o, {}, saved); });
            std::filesystem::copy(saved, edited, std::filesystem::copy_options::recursive);
            J j; {std::ifstream in(saved / "checkpoint.json"); in >> j;}
            j["field"]["B"][face] = j["field"]["B"][face].get<double>() + .01;
            j["field_sha256"] = model->state_sha256(hc::checkpoint_detail::field(j.at("field"), g), hc::checkpoint_detail::charge(j.at("charge"), g));
            {std::ofstream f(edited / "checkpoint.json", std::ios::binary); f << j.dump() << '\n';}
            dv["restore_resealed_perturbed"] = refusal([&] { (void)hc::Owner::restore(model, o, {}, edited); });
        }
        J sweep = J::array(); bool all_rejected = true;
        for (int p = 2; p <= 13; ++p) { auto s = seed; s.B[face] += std::pow(10., -p);
            const auto m = refusal([&] { (void)model->admit(s, control); }); all_rejected = all_rejected && !m.empty();
            sweep.push_back({{"delta_T", std::pow(10., -p)}, {"refusal", m}}); }
        dv["detection_sweep_1e-2_to_1e-13_T"] = sweep;
        R.record("audit-A09", "about 50 T/m divB (0.01 T face) is rejected at admit, prepare, bootstrap and resealed restore; 1e-13 T still rejected",
                 !dv["admit"].get<std::string>().empty() && !dv["prepare"].get<std::string>().empty() && !dv["bootstrap"].get<std::string>().empty() &&
                     dv["restore_unmodified_control"] == "" && !dv["restore_resealed_perturbed"].get<std::string>().empty() && all_rejected, dv);
        // A10: compensated midpoint is authoritative.
        J a10;
        {
            hc::Owner live(model, seed, control, o, {});
            auto f = live.prepare_vacuum(dt);
            const auto& s1 = f.candidate_field();
            std::size_t e_hw = 0, b_hw = 0, e_bad = 0, b_bad = 0;
            for (int e = 0; e < g.ne; ++e) {
                if (.5 * (seed.E[e] + s1.E[e]) != f.midpoint().E()[e]) ++e_hw;
                if (((hc::Words(seed.E[e], seed.E_low[e]) + hc::Words(s1.E[e], s1.E_low[e])) * .5).rounded() != f.midpoint().E()[e]) ++e_bad; }
            for (int b2 = 0; b2 < g.nb; ++b2) {
                if (.5 * (seed.B[b2] + s1.B[b2]) != f.midpoint().B()[b2]) ++b_hw;
                if (((hc::Words(seed.B[b2], seed.B_low[b2]) + hc::Words(s1.B[b2], s1.B_low[b2])) * .5).rounded() != f.midpoint().B()[b2]) ++b_bad; }
            a10["field_step"] = {{"E_edges_differing_from_high_word_average", e_hw}, {"B_faces_differing_from_high_word_average", b_hw},
                                 {"E_not_compensated_average", e_bad}, {"B_not_compensated_average", b_bad},
                                 {"receipt_mid_E_equals_authoritative", f.field_receipt().mid_E == f.midpoint().E()}};
            R.record("audit-A10", "field-only step: authoritative midpoint equals the compensated average and differs from the high-word average",
                     e_bad == 0 && b_bad == 0 && e_hw > 0 && f.field_receipt().mid_E == f.midpoint().E(), a10["field_step"]);
        }
        J closed = J::array(); bool closed_ok = true;
        for (double W : {1e9, 1e11, 1e12}) {
            const std::vector<n::ArenaRateEntry> rows{
                mf::row(1, "A10-e", "electron", {audit_point[0] + o.width_m / 2, audit_point[1] + o.gap_m / 2, audit_point[2]}, {2e6, -1e6, 5e5}, -n::pic_e, n::pic_me, W),
                mf::row(2, "A10-p", "proton", {audit_point[0] + o.width_m / 2, audit_point[1] + o.gap_m / 2, audit_point[2]}, {0, 0, 0}, n::pic_e, hall::proton_mass, W)};
            J c{{"weight", W}};
            try {
                hc::Owner owner(model, seed, declared_dc(*model, seed, hc::charge_of_census(g, o, rows), "A10"), o, rows);
                auto p = hc::prepare_closed_proposal(owner, o, rows, binding_for(owner, "A10"), hc::sha256("A10/source"));
                const auto md = hd::total_field_sha256(p.midpoint().E(), p.midpoint().B());
                c["iterations"] = p.iterations(); c["owner_preflight"] = owner.preflight(p, rows);
                c["motion_bound_to_midpoint_digest"] = p.motion().binding().trial_field_sha256 == md;
                c["digest_differs_from_role_seal"] = md != p.midpoint().seal();
                closed_ok = closed_ok && c["owner_preflight"].get<bool>() && c["motion_bound_to_midpoint_digest"].get<bool>() && c["digest_differs_from_role_seal"].get<bool>();
            } catch (const std::exception& e) { c["refusal"] = e.what(); closed_ok = false; }
            closed.push_back(c);
        }
        R.record("audit-A10", "closed proposals for e/p pairs (1e9..1e12) converge, preflight, and bind motion to the midpoint digest", closed_ok, closed);
    }

    // ---- 4. backend contract: the independent residual check is mandatory
    {
        J ev;
        for (auto [label, scale, one] : std::vector<std::tuple<std::string, double, bool>>{{"increment_scaled_1.001", 1.001, false},
                                                                                         {"one_component_scaled_1+1e-6", 1 + 1e-6, true}}) {
            auto bad = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, std::make_shared<ScaledBackend>(scale, one));
            const auto control = declared_dc(*bad, seed, {hc::WordVec(g.nn), hc::Vec(g.nn)}, "backend");
            hc::Owner owner(bad, seed, control, o, {});
            ev[label] = refusal([&] { (void)owner.prepare_vacuum(dt); });
        }
        {
            auto ok_model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, std::make_shared<ScaledBackend>(1., false));
            hc::Owner owner(ok_model, seed, declared_dc(*ok_model, seed, {hc::WordVec(g.nn), hc::Vec(g.nn)}, "backend"), o, {});
            ev["unscaled_test_backend_control"] = refusal([&] { (void)owner.prepare_vacuum(dt); });
        }
        R.record("backend", "a 0.1% increment error and a 1e-6 one-component error are refused by the independent FP64 residual check",
                 ev["unscaled_test_backend_control"] == "" &&
                     ev["increment_scaled_1.001"].get<std::string>().find("independent FP64 field residual check failed") == 0 &&
                     ev["one_component_scaled_1+1e-6"].get<std::string>().find("independent FP64 field residual check failed") == 0, ev);
    }

    // ---- 4b. configured Hall carrier cap with a memory budget
    {
        J ev;
        ev["cap_200000_budget_1GiB"] = refusal([&] { hc::Model m(hc::MeshSpec{}, trace, hc::make_cpu_pcg(), hc::CarrierLimits{200000, 1ull << 30}); });
        const auto census = mf::pairs(100002, 77, 1.);
        ev["default_cap_census_100002"] = refusal([&] {
            hc::Owner x(model, seed, declared_dc(*model, seed, {hc::WordVec(g.nn), hc::Vec(g.nn)}, "cap"), o, census); });
        auto big = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, hc::make_cpu_pcg(), hc::CarrierLimits{150000, 1ull << 30});
        ev["configuration_differs_from_default_cap"] = big->configuration_sha256() != model->configuration_sha256();
        const auto t = steady::now();
        ev["cap_150000_census_100002_closed_step"] = refusal([&] {
            hc::Owner x(big, seed, declared_dc(*big, seed, hc::charge_of_census(g, o, census, 150000), "cap"), o, census);
            auto p = hc::prepare_closed_proposal(x, o, census, binding_for(x, "cap"), hc::sha256("cap/source"));
            hc::need(x.preflight(p, census) && !p.motion().seal_sha256().empty() && !p.motion().record_v3_sha256().empty(), "cap step preflight");
            ev["cap_150000_timings"] = {{"motion_s", p.timings().motion_s}, {"motion_v2_seal_s", p.timings().motion_v2_seal_s},
                                        {"field_prepare_s", p.timings().field_prepare_s}, {"total_s", p.timings().total_s}};
        });
        ev["wall_s"] = since(t);
        R.record("carrier-cap", "Hall carrier cap is configured with a memory budget, bound into the configuration, and enforced "
                                "(1e5 default refuses 100002 carriers; a 1.5e5 cap admits and preflights them)",
                 ev["cap_200000_budget_1GiB"] == "Hall carrier limits exceed the memory budget" &&
                     ev["default_cap_census_100002"] == "motion census exceeds sealed native carrier cap" &&
                     ev["configuration_differs_from_default_cap"].get<bool>() && ev["cap_150000_census_100002_closed_step"] == "", ev);
    }

    // ---- 5. owner protocol: two consecutive published closed steps
    {
        auto census = mf::pair({.0045, .0045, 0.}, 1.);
        hc::Owner owner(model, seed, declared_dc(*model, seed, hc::charge_of_census(g, o, census), "protocol"), o, census);
        J ev = J::array(); bool ok = true;
        for (int step = 0; step < 2; ++step) {
            const std::string src = hc::sha256("protocol/source");
            auto p = hc::prepare_closed_proposal(owner, o, census, binding_for(owner, "protocol-" + std::to_string(step)), src);
            const auto bnd = boundary_of(p, src);
            const bool pre = owner.preflight_final(p, census, bnd);
            auto stale = bnd; stale.time_s = std::nextafter(stale.time_s, 1.);
            const bool stale_refused = !owner.preflight_final(p, census, stale);
            const bool pre2 = owner.preflight_final(p, census, bnd);
            const auto candidate = p.candidate_field_sha256();
            std::vector<n::ArenaRateEntry> next = census;
            for (std::size_t i = 0; i < next.size(); ++i) next[i].native = p.motion().endpoints()[i];
            if (pre2) owner.publish_prevalidated(std::move(p));
            const bool advanced = owner.revision() == std::uint64_t(step + 1) && owner.field_sha256() == candidate &&
                                  owner.field_sha256() == model->state_sha256(owner.field(), owner.charge());
            ev.push_back({{"step", step + 1}, {"preflight_final", pre}, {"stale_boundary_refused", stale_refused}, {"published", pre2}, {"advanced_with_cached_digest", advanced}});
            ok = ok && pre && stale_refused && pre2 && advanced;
            census = next;
        }
        R.record("protocol", "two consecutive closed steps preflight, refuse a stale boundary, publish, and the cached digest equals a recomputation", ok, ev);
    }
    report["cases"] = R.cases;
    report["failed"] = R.failed;
    std::ofstream(out / "components.json") << report.dump(1) << '\n';
    return R.failed ? 1 : 0;
}

// ------------------------------------------------------------------ motion dump
std::string after_digest(const hc::FreeCharge& q) {
    hc::Digest d("test-after-charge"); hc::encode(d, q.net_C); d.f64s(q.gross_C); return d.hex();
}
int motion_dump(const std::string& trace, const std::string& file) {
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace);
    const auto& g = model->mesh();
    const auto o = mf::geometry(g);
    const auto ctx = hc::make_motion_context(g, o, model->mesh_json());
    const auto mid = model->at(model->dc_seed());
    hc::Vec E100 = mid.E(); for (double& v : E100) v *= 100;
    J out = J::array();
    for (const auto& f : mf::fixture_set(g)) {
        const hc::Vec& E = f.field ? E100 : mid.E();
        hd::MotionBinding b;
        b.phase = hd::ParticlePhase::full_time_position_and_proper_velocity;
        b.accepted_step = 0; b.staged_arena_revision = 1; b.begin_s = model->drive().start(); b.dt_s = dt; b.end_s = b.begin_s + dt;
        b.staged_arena_sha256 = hc::sha256("dump/staged"); b.trial_field_sha256 = hd::total_field_sha256(E, mid.B());
        b.mesh_sha256 = ctx.mesh_sha256; b.material_sha256 = model->material_sha256(); b.species_sha256 = hc::sha256("dump/species");
        b.common_clock_sha256 = hc::sha256("dump/clock"); b.kinetic_policy_sha256 = ctx.policy_sha256;
        for (int api = 0; api < 2; ++api) {
            J r{{"fixture", f.name}, {"carriers", f.census.size()}, {"entry", api ? "prepare_canonical_motion" : "MotionAssembly::cpu"}};
            const auto t = steady::now();
            try {
                auto m = api ? hd::prepare_canonical_motion(g, o, f.census, b, {}, E, mid.B()) : hc::MotionAssembly::cpu(ctx, f.census, b, {}, E, mid.B());
                hc::MotionAssembly::seal_v2(m);
                J rec = hd::motion_record(m); const auto text = rec.dump(); hc::drain(rec);
                r["record_sha256"] = hc::sha256(text); r["seal_sha256"] = m.seal_sha256(); r["record_v3_sha256"] = m.record_v3_sha256();
                r["record_v3_recomputed_equal"] = hc::MotionAssembly::record_v3(m) == m.record_v3_sha256();
                try { r["after_sha256"] = after_digest(hc::charge_of_endpoints(g, o, m.endpoints())); }
                catch (const std::exception& e) { r["after_refusal"] = e.what(); }
            } catch (const std::exception& e) { r["refusal"] = e.what(); }
            r["wall_s"] = since(t);
            out.push_back(r);
        }
    }
    std::ofstream(file) << out.dump(1) << '\n';
    return 0;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 4) throw std::runtime_error("usage: test_hall_v2 components|motion-dump TRACE OUT");
        const std::string mode = argv[1];
        if (mode == "components") { std::filesystem::create_directories(argv[3]); return components(argv[2], argv[3]); }
        if (mode == "motion-dump") return motion_dump(argv[2], argv[3]);
        throw std::runtime_error("unknown mode");
    } catch (const std::exception& e) { std::cerr << "fatal: " << e.what() << '\n'; return 2; }
}
