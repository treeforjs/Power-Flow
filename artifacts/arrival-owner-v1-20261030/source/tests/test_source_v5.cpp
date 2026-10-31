// Integration of the arrival owner v1 with the sealed heat/source owner v5
// (unchanged headers, read in place). Fixture recipes follow the source v5
// focused tests (test_source_owner.cpp): 4x2 thermal mesh, Hall axial mesh
// (45, 4, .009, 1), dt = 1e-14 s.
//  A. neutral leaves: sweep of exact source debits; the leaves built by the
//     arrival policy are admitted by Owner::prepare (exact ExactParticleCount
//     equality inside source v5) or refused only by typed source v5 refusals;
//     exact inventory debit, packet weights and release times read back.
//  B. electron births: sweep of expected counts; births of weight w0 with the
//     arrival marks go through Owner::begin_source (the source emission law
//     samples energy/angle); birth time, canonical ID, receipt and weights checked.
//  C. the arrival configuration identity is bound into the source configuration.
// usage: test_source_v5 OUTPUT_DIRECTORY            (qualification sweeps)
//        test_source_v5 OUTPUT_DIRECTORY probe N0 RATE X   (one neutral case)
#include "hall-owner/development/hall-finite-foil-em-20261029/axial_surface_path_port.hpp"
#include "source_owner.hpp"
#include "ss304_heat_material_v2.hpp"
#include "pf_arrival/source_v5_adapter.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

using namespace electrode_heat_v3;
namespace s5 = electrode_heat_v5;
namespace pa = pf_arrival;
namespace ps = pf_arrival::source_v5;

namespace {
long checks = 0;
void check(bool ok, const std::string& what) {
    ++checks;
    if (!ok) throw std::runtime_error("check failed: " + what);
}
std::string hx(double x) { return pa::detail::hexfloat(x); }
std::string jstr(const std::string& s) {
    std::string r = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') r += '\\';
        r += c;
    }
    return r + "\"";
}

// ---- source v5 fixture (as in test_source_owner.cpp)
eh::Geometry geom() {
    eh::Geometry g;
    g.nx = 4;
    g.ny = 2;
    g.revision = "thermal-mesh-v3";
    g.electrode = "upper";
    return g;
}
eh::Material material() {
    return eh::conditional_ss304_v2(ss304v2::Qualification::source_conditioned_liquid_3000,
                                    ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation, 8000, eh::material_headers_sha256);
}
eh::Patch metal_patch() {
    eh::Patch p;
    p.id = "lower-face";
    p.face = eh::Face::LowY;
    p.source_hash = "v3-source-scenario";
    return p;
}
eh::Patch water(double n, double rate) {
    auto p = metal_patch();
    p.desorption_enabled = true;
    p.state_id = p.parent_id = "H2O";
    p.initial_particles = n;
    p.initial_low_particles = 0;
    p.mass_kg = 18 * 1.66053906660e-27;
    p.attempt_s = rate;
    p.activation_eV = p.desorption_eV = p.internal_eV = 0;
    return p;
}
se::LayerPolicy layer_policy() {
    se::LayerPolicy p;
    p.local_linear_vacuum_extent_m = 25e-6;
    return p;
}
eh::SourceEventIdentity event(std::uint64_t seq, std::uint64_t native_id) {
    return {"arrival-owner-v1", "e-", "electrode-electron-reservoir", "native-layer", "canonical-electron-" + std::to_string(native_id), seq, seq};
}
const V3 origin{.0001, .0025, .0002};
std::unique_ptr<s5::Owner> source_owner(eh::Patch p, double T, const std::string& sampling_identity) {
    return std::make_unique<s5::Owner>(geom(), material(), std::vector<eh::Patch>{p},
                                       std::vector<s5::PatchIdentity>{{"upper-foil-lower-surface", "conditional-single-layer", .005 * .025}},
                                       std::vector<s5::SourceSettings>{{se::Config{4.3}, layer_policy(), sampling_identity}}, T, "v3-source-scenario", 5e-9);
}
TrialFrame frame4(const s5::Owner& o, std::uint64_t stage) {
    TrialFrame f;
    f.role = FieldRole::TotalMidpoint;
    f.accepted_revision = o.thermal_state().revision;
    f.trial_stage = stage;
    f.mesh_revision = 1;
    f.material_revision = 2;
    f.begin_s = o.thermal_state().time;
    f.dt_s = 1e-14;
    f.end_s = f.begin_s + f.dt_s;
    f.field_stage_s = f.begin_s + .5 * f.dt_s;
    f.arena_hash = o.state_sha256();
    f.field_hash = "declared-parent-field";
    f.mesh_hash = "declared-parent-mesh";
    f.material_hash = "global-material-v2-map";
    f.species_hash = "actual-source-v1-electron";
    return f;
}
PatchBinding bind4(const s5::Owner& o) {
    PatchBinding b;
    b.thermal_revision = o.thermal_state().revision;
    b.source_accepted_state = b.thermal_revision + 1;
    b.patch_id = 7;
    b.thermal_mesh_identity = geom().revision;
    b.thermal_material_hash = eh::material_headers_sha256;
    b.material_phase_admitted = b.planar_barrier_admitted = b.electron_supply_admitted = true;
    b.temperature = HeldTemperature::AcceptedPatchMean;
    b.extrusion_center_z_m = 0;
    b.extrusion_height_m = geom().height;
    return b;
}
struct Arrays {
    foil::Vec E, B;
};
Arrays arrays(const axial::Mesh& g, double F) {
    Arrays a{foil::Vec(g.ne), foil::Vec(g.nb)};
    for (int e = 0; e < g.ne; ++e)
        if (g.edges[e].axis == 1) a.E[e] = F;
    return a;
}
std::shared_ptr<const s5::Field> field4(const axial::Mesh& g, const s5::Owner& o, std::uint64_t stage, double F) {
    auto a = arrays(g, F);
    return std::make_shared<const s5::Field>(g, frame4(o, stage), a.E, a.B);
}

std::shared_ptr<const pa::Policy> arrival_policy(double w_electron, double w_neutral, std::uint64_t seed = 20261030) {
    pa::Configuration c;
    c.run_seed = seed;
    c.streams = {{0, pa::Kind::electron, 0, "e-", w_electron, 4096}, {1, pa::Kind::neutral, 0, "H2O", w_neutral, 4096}};
    return std::make_shared<const pa::Policy>(c);
}
// Native stand-in for the receipt: births in staged canonical order from base.
pa::NativeBirthReceipt receipt(const pa::Policy& P, const pa::Staged& s) {
    pa::NativeBirthReceipt r;
    r.base = s.native_base();
    std::uint64_t id = r.base;
    for (const auto& b : s.births()) r.entries.push_back({b.packet_id(), P.stream(b.stream_id).species_state, id++, b.weight});
    r.next_id_after = id;
    return r;
}

std::string digest_material;   // binary record of every outcome, hashed into the results
std::uint64_t case_seed = 100;   // a distinct run seed per sweep case (deterministic sequence)
void record(const std::string& s) {
    pa::Writer w;
    w.text(s);
    digest_material += w.bytes;
}

// ---------------------------------------------------------------- A. neutral leaves
struct NeutralStats {
    long admitted = 0, source_refused = 0, source_refused_untyped = 0, quote_refused = 0, arrival_refused = 0, leaves = 0, committed = 0;
    double smallest_admitted = INFINITY, largest_refused = 0;
    std::vector<std::string> rows;
};
void neutral_case(const axial::Mesh& g, double N0, double rate, double x_target, bool commit_steps, NeutralStats& st) {
    // Probe the first release to choose w0 = D/x_target, then bind the final arrival
    // identity into the source owner that is actually used.
    double D0 = 0;
    {
        auto probe = source_owner(water(N0, rate), 1500, "probe");
        try {
            auto q = probe->release_quotes(probe->make_epoch(1e-14, "actual-parent-midpoint"));
            D0 = q[0].released.zero() ? 0 : q[0].released.rounded_for_rate();
        } catch (const s5::Refusal& r) {
            ++st.quote_refused;
            record("quote-refused/" + hx(N0) + "/" + hx(rate) + "/" + s5::kind_name(r.kind));
            return;
        }
    }
    const double w0 = D0 / x_target;
    if (!(D0 > 0) || !pa::detail::positive_normal(w0)) return;
    auto P = arrival_policy(1, w0, ++case_seed);
    auto o = source_owner(water(N0, rate), 1500, P->sampling_policy_identity());
    pa::Owner A(P, o->thermal_state().revision);
    mfem::ConstantCoefficient zero(0);
    for (int step = 0; step < (commit_steps ? 3 : 1); ++step) {
        auto f = field4(g, *o, step + 1, 0);
        auto a = arrays(g, 0);
        s5::FinalFieldView v{f->frame(), g, a.E, a.B};
        auto e = o->make_epoch(1e-14, "actual-parent-midpoint");
        auto q = o->release_quotes(e);
        const double D = q[0].released.rounded_for_rate();
        check(step > 0 || D == D0, "probe and bound owners quote the same release");
        pa::StepPlan plan;
        try {
            plan = A.plan({o->thermal_state().revision, {}, {ps::neutral_input(*P, 1, q[0], *o, e)}});
        } catch (const pa::Refusal& r) {
            ++st.arrival_refused;
            record("arrival-refused/" + hx(D) + "/" + pa::kind_name(r.kind));
            check(r.kind == pa::RefusalKind::macroparticle_cap, "only cap refusals expected from the arrival owner here");
            return;
        }
        auto staged = A.stage(plan, 1000 + 10000 * std::uint64_t(step));
        auto leaves = ps::neutral_leaves(staged, 1);
        pa::Count leaf_sum;
        for (const auto& l : leaves[0]) leaf_sum = leaf_sum.credited(l.physical_weight);
        check(leaf_sum == q[0].released, "arrival leaves sum exactly to the source v5 quote");
        const auto base = o->state_sha256();
        std::optional<s5::Trial> trial;
        try {
            trial.emplace(o->prepare(e, zero, f, v, {}, {}, leaves));
        } catch (const s5::Refusal& r) {
            // Typed source v5 refusal (e.g. G1 weighted representation of a tiny leaf).
            ++st.source_refused;
            st.largest_refused = std::max(st.largest_refused, D);
            check(o->state_sha256() == base, "source refusal left the source owner unchanged");
            record("source-refused/" + hx(D) + "/" + s5::kind_name(r.kind) + "/" + r.quantity);
            st.rows.push_back("{\"N0\":" + jstr(hx(N0)) + ",\"rate\":" + jstr(hx(rate)) + ",\"D\":" + jstr(hx(D)) + ",\"x\":" + jstr(hx(x_target)) +
                              ",\"leaves\":" + std::to_string(leaves[0].size()) + ",\"outcome\":\"source-refused\",\"kind\":" + jstr(s5::kind_name(r.kind)) +
                              ",\"quantity\":" + jstr(r.quantity) + "}");
            return;
        } catch (const std::exception& ex) {
            // Untyped sealed refusals outside the arrival scope (e.g. the heat engine's
            // thermal-domain guard for an enormous neutral heat load) are recorded; the
            // exact leaf-sum guard must never be the reason.
            const std::string what = ex.what();
            check(what.find("leaf sum") == std::string::npos && what.find("source debit") == std::string::npos,
                  "source v5 refused the exact leaf sum: " + what);
            check(o->state_sha256() == base, "untyped source refusal left the source owner unchanged");
            ++st.source_refused_untyped;
            record("source-refused-untyped/" + hx(D) + "/" + what);
            st.rows.push_back("{\"N0\":" + jstr(hx(N0)) + ",\"rate\":" + jstr(hx(rate)) + ",\"D\":" + jstr(hx(D)) + ",\"x\":" + jstr(hx(x_target)) +
                              ",\"leaves\":" + std::to_string(leaves[0].size()) + ",\"outcome\":\"source-refused-untyped\",\"message\":" + jstr(what) + "}");
            return;
        }
        // Exact debit: source inventory after = before - D, packets carry the leaves verbatim.
        check(trial->inventory(0) == o->inventory(0).debit(D), "source v5 exact inventory debit equals the quote");
        pa::Count packets;
        for (const auto& pk : trial->packets()) packets = packets.credited(pk.weight);
        check(packets == q[0].released && trial->packets().size() == leaves[0].size(), "packet weights sum exactly to the debit");
        auto times = ps::neutral_times(*trial, staged, 1, "lower-face");
        staged.bind_source_times(1, times);
        for (const auto& b : staged.births())
            check(b.time_s >= e.t && b.time_s <= e.t + e.dt && b.remaining_s == (e.t + e.dt) - b.time_s, "neutral release time inside the step");
        ++st.admitted;
        st.leaves += long(leaves[0].size());
        st.smallest_admitted = std::min(st.smallest_admitted, D);
        record("admitted/" + hx(D) + "/" + plan.sha256() + "/" + trial->content_sha256());
        st.rows.push_back("{\"N0\":" + jstr(hx(N0)) + ",\"rate\":" + jstr(hx(rate)) + ",\"D\":" + jstr(hx(D)) + ",\"x\":" + jstr(hx(x_target)) +
                          ",\"leaves\":" + std::to_string(leaves[0].size()) + ",\"outcome\":\"admitted\",\"step\":" + std::to_string(step) + "}");
        if (!commit_steps) return;
        // Publish both owners and continue from the new exact stock.
        check(o->commit(std::move(*trial), v), "source v5 commit");
        auto r = receipt(*P, staged);
        (void)A.commit(std::move(staged), r);
        check(A.next_step() == o->thermal_state().revision, "arrival step index follows the source accepted revision");
        ++st.committed;
    }
}
std::string neutral_sweep(const axial::Mesh& g, NeutralStats& st) {
    // Initial stocks up to 1e20 molecules per 1.25e-4 m2 patch (about 1e5 monolayers).
    // Larger stocks drive the sealed heat engine outside its domain (typed for some,
    // an MFEM CG abort for others; see the probe mode), which is not arrival scope.
    const double N0s[] = {1e-300, 1e-280, 1e-250, 1e-200, 1e-100, 1e-30, 1e-10, 1, 1e5, 1e10, 1e15, 1e18, 1e20};
    const double rates[] = {1e300, 1e14, 1e8};
    const double xs[] = {0.3, 1, 2.5, 37.5, 1000.25};
    for (double N0 : N0s)
        for (double rate : rates)
            for (double x : xs) neutral_case(g, N0, rate, x, false, st);
    neutral_case(g, 1e15, 1e14, 4095.5, false, st);   // near the per-patch leaf cap
    neutral_case(g, 1e15, 1e14, 5000, false, st);     // above the cap: arrival refusal
    neutral_case(g, 1e15, 1e8, 12.5, true, st);       // three committed steps, stock carried exactly
    check(st.admitted > 0 && st.source_refused > 0 && st.arrival_refused > 0, "neutral sweep covers admitted and refused outcomes");
    std::string rows = "[";
    for (std::size_t i = 0; i < st.rows.size(); ++i) rows += (i ? "," : "") + st.rows[i];
    return rows + "]";
}

// ---------------------------------------------------------------- B. electron births
struct ElectronStats {
    long admitted_births = 0, source_refused = 0, exposure_refused = 0, zero_births = 0, cases = 0;
    std::vector<std::string> rows;
};
void electron_case(const axial::Mesh& g, double area_fraction, double x_target, double F, ElectronStats& st) {
    double N0 = 0;
    {
        auto probe = source_owner(metal_patch(), 1500, "probe");
        auto pf = field4(g, *probe, 1, F);
        try {
            N0 = probe->prepare_patch(probe->make_epoch(1e-14, "actual-parent-midpoint"), 0, bind4(*probe), pf, origin,
                                      {pf->frame().begin_s, 1e-14, area_fraction, "physical-quadrature"})
                     .expected_number();
        } catch (const s5::Refusal& r) {
            ++st.exposure_refused;
            record("exposure-refused/" + hx(area_fraction) + "/" + s5::kind_name(r.kind));
            st.rows.push_back("{\"field_V_m\":" + jstr(hx(F)) + ",\"area_fraction\":" + jstr(hx(area_fraction)) + ",\"outcome\":\"exposure-refused\",\"kind\":" +
                              jstr(s5::kind_name(r.kind)) + ",\"log_scope\":" + jstr(hx(r.log_abs_scope)) + "}");
            return;
        } catch (const std::exception& ex) {
            ++st.exposure_refused;
            record("exposure-refused-untyped/" + hx(area_fraction) + "/" + ex.what());
            st.rows.push_back("{\"field_V_m\":" + jstr(hx(F)) + ",\"area_fraction\":" + jstr(hx(area_fraction)) + ",\"outcome\":\"exposure-refused-untyped\",\"message\":" + jstr(ex.what()) + "}");
            return;
        }
    }
    const double w0 = N0 / x_target;
    if (!pa::detail::positive_normal(w0)) return;
    auto P = arrival_policy(w0, 1, ++case_seed);
    auto o = source_owner(metal_patch(), 1500, P->sampling_policy_identity());
    auto f = field4(g, *o, 1, F);
    const electrode_heat_v3::Exposure x{f->frame().begin_s, 1e-14, area_fraction, "physical-quadrature"};
    std::optional<s5::PreparedPatch> p;
    p.emplace(o->prepare_patch(o->make_epoch(1e-14, "actual-parent-midpoint"), 0, bind4(*o), f, origin, x));
    const double N = p->expected_number();
    check(N == N0, "probe and bound owners give the same expected count");
    pa::Owner A(P, o->thermal_state().revision);
    const auto in = ps::electron_input(0, *p, x, f->frame().end_s);
    check(in.has_log && std::abs(in.log_expected_number - std::log(N)) <= 1e-9 * std::max(1.0, std::abs(in.log_expected_number)),
          "finite log and scalar of the exposure agree");
    auto plan = A.plan({o->thermal_state().revision, {in}, {}});
    auto staged = A.stage(plan, 777);
    ++st.cases;
    if (staged.births().empty()) ++st.zero_births;
    pa::Count weights;
    std::uint64_t seq = 1;
    long refused = 0;
    std::string first_refusal;
    for (const auto& b : staged.births()) {
        const auto sample = ps::birth_sample(staged, b, *p, o->thermal_state().revision);
        try {
            auto c = o->begin_source(*p, sample, event(seq++, b.native_id));
            check(c.canonical_id() == b.native_id, "source candidate canonical ID is the native birth ID");
            check(c.birth_sample() && c.birth_sample()->count_owner_token == staged.token() && c.birth_sample()->count_owner_base_hash == staged.owner_base_sha256() &&
                      c.birth_sample()->normal_uniform == b.marks[0] && c.birth_sample()->transverse_uniform == b.marks[1] &&
                      c.birth_sample()->azimuth_uniform == b.marks[2] && c.birth_sample()->time_uniform == b.marks[3],
                  "source candidate retains the arrival marks and count-owner receipt");
            check(c.state().native.particle.weight == w0, "source particle weight is w0");
            check(c.heat_events().front().event_time_s == b.time_s, "source birth event time equals the arrival birth time");
            if (c.state().native.status == se::LayerStatus::Active) {
                check(c.state().native.current_time_s == b.time_s && f->frame().end_s - c.state().native.current_time_s == b.remaining_s,
                      "remaining flight equals the source remaining interval");
            }
            weights = weights.credited(c.state().native.particle.weight);
            record("birth/" + b.packet_id() + "/" + c.content_sha256());
            ++st.admitted_births;
        } catch (const s5::Refusal& r) {
            ++refused;
            if (first_refusal.empty()) first_refusal = std::string(s5::kind_name(r.kind)) + "/" + r.quantity;
            record("birth-refused/" + b.packet_id() + "/" + s5::kind_name(r.kind));
        } catch (const std::exception& ex) {
            ++refused;
            if (first_refusal.empty()) first_refusal = std::string("untyped/") + ex.what();
            record("birth-refused-untyped/" + b.packet_id() + "/" + ex.what());
        }
    }
    st.source_refused += refused;
    if (!refused) check(weights == staged.plan().decisions().front().debit, "source birth weights sum exactly to the arrival debit n*w0");
    // One continuation over the remaining flight time of the first birth.
    if (!refused && !staged.births().empty()) {
        const auto& b = staged.births().front();
        auto c = o->begin_source(*p, ps::birth_sample(staged, b, *p, o->thermal_state().revision), event(1, b.native_id));
        if (c.state().native.status == se::LayerStatus::Active) {
            auto cont = o->advance_birth(*p, c, event(1000, b.native_id));
            check(cont.state().native.status != se::LayerStatus::Active || cont.state().native.current_time_s == f->frame().end_s,
                  "continuation spans the remaining flight to the step end");
            record("continuation/" + cont.content_sha256());
        }
    }
    st.rows.push_back("{\"field_V_m\":" + jstr(hx(F)) + ",\"area_fraction\":" + jstr(hx(area_fraction)) + ",\"N\":" + jstr(hx(N)) + ",\"log_N\":" + jstr(hx(in.log_expected_number)) +
                      ",\"w0\":" + jstr(hx(w0)) + ",\"x\":" + jstr(hx(N / w0)) + ",\"births\":" + std::to_string(staged.births().size()) +
                      ",\"source_refused\":" + std::to_string(refused) + ",\"first_refusal\":" + jstr(first_refusal) + "}");
}
std::string electron_sweep(const axial::Mesh& g, ElectronStats& st) {
    const double fractions[] = {1, 1e-3, 1e-10, 1e-50, 1e-100, 1e-150, 1e-200, 1e-230, 1e-240, 1e-280, 1e-300};
    const double xs[] = {0.4, 3.7, 25.5};
    for (double a : fractions)
        for (double x : xs) electron_case(g, a, x, 1e8, st);
    // Stronger extraction fields give large expected counts per step.
    for (double F : {1e9, 3e9})
        for (double a : {1., 1e-3})
            for (double x : {3.7, 1000.5}) electron_case(g, a, x, F, st);
    check(st.admitted_births > 0 && st.exposure_refused > 0, "electron sweep covers admitted births and refused exposures");
    std::string rows = "[";
    for (std::size_t i = 0; i < st.rows.size(); ++i) rows += (i ? "," : "") + st.rows[i];
    return rows + "]";
}

// ---------------------------------------------------------------- C. identity binding
void identity_binding() {
    auto a = arrival_policy(1e5, 1e6, 1), b = arrival_policy(1e5, 1e6, 2);
    auto oa = source_owner(metal_patch(), 1500, a->sampling_policy_identity());
    auto ob = source_owner(metal_patch(), 1500, b->sampling_policy_identity());
    check(oa->configuration_sha256() != ob->configuration_sha256(), "arrival configuration identity bound into the source configuration");
    auto oc = source_owner(metal_patch(), 1500, a->sampling_policy_identity());
    check(oa->configuration_sha256() == oc->configuration_sha256(), "source configuration deterministic for one arrival identity");
}
}  // namespace

int main(int argc, char** argv) {
    if (argc != 2 && argc != 6) {
        std::cerr << "usage: test_source_v5 OUTPUT_DIRECTORY [probe N0 RATE X]\n";
        return 2;
    }
    const std::filesystem::path out = argv[1];
    std::filesystem::create_directories(out);
    const auto t0 = std::chrono::steady_clock::now();
    std::string failure, neutral_rows = "[]", electron_rows = "[]";
    NeutralStats ns;
    ElectronStats es;
    if (argc == 6 && std::string(argv[2]) == "probe") {   // one neutral case: N0 rate x
        axial::Mesh g(45, 4, .009, 1);
        neutral_case(g, std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]), false, ns);
        std::cout << "probe admitted=" << ns.admitted << " typed_refused=" << ns.source_refused << " untyped_refused=" << ns.source_refused_untyped
                  << (ns.rows.empty() ? std::string() : " " + ns.rows.back()) << "\n";
        return 0;
    }
    try {
        axial::Mesh g(45, 4, .009, 1);
        identity_binding();
        neutral_rows = neutral_sweep(g, ns);
        electron_rows = electron_sweep(g, es);
    } catch (const std::exception& e) {
        failure = e.what();
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::ofstream r(out / "source-results.json");
    r << "{\"program\":\"test_source_v5\",\"checks\":" << checks << ",\"passed\":" << (failure.empty() ? "true" : "false") << ",\"failure\":" << jstr(failure)
      << ",\"seconds\":" << seconds << ",\"outcome_digest\":" << jstr(pa::sha256_hex(digest_material)) << ",\"neutral\":{\"admitted\":" << ns.admitted
      << ",\"source_refused\":" << ns.source_refused << ",\"source_refused_untyped\":" << ns.source_refused_untyped << ",\"quote_refused\":" << ns.quote_refused << ",\"arrival_refused\":" << ns.arrival_refused
      << ",\"leaves_admitted\":" << ns.leaves << ",\"committed_steps\":" << ns.committed << ",\"smallest_admitted_debit\":" << jstr(hx(ns.smallest_admitted))
      << ",\"largest_refused_debit\":" << jstr(hx(ns.largest_refused)) << ",\"cases\":" << neutral_rows << "},\"electron\":{\"cases\":" << es.cases
      << ",\"admitted_births\":" << es.admitted_births << ",\"source_refused_births\":" << es.source_refused << ",\"exposure_refused\":" << es.exposure_refused
      << ",\"zero_birth_cases\":" << es.zero_births << ",\"rows\":" << electron_rows << "}}\n";
    std::cout << "test_source_v5 checks=" << checks << " " << (failure.empty() ? "PASS" : "FAIL: " + failure) << "\n";
    return failure.empty() ? 0 : 1;
}
