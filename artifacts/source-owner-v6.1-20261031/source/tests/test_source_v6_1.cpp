// Heat/source owner v6.1 qualification against the sealed owner v6, in one program.
//   test_source_v6_1 OUTDIR
// Groups:
//   emission   exact fast spectrum: every Distribution field and node bit-identical to the
//              sealed surface_emission::evaluate over a deterministic sweep of temperatures,
//              signed fields and configurations; identical refusals outside the domain
//   a02        the A02 audit probes of the source v5 audit port (changed particle weight,
//              substituted spectrum, reused owner address, mutated mesh with old metadata,
//              mutation between preflight and publish), ported to the v6 API and run on the
//              sealed v6 owner and on v6.1 with the same inputs: same verdicts, v6.1 refuses
//              every mutation that v6 refuses
//   identity   configuration, state, prepared-patch, candidate, trial, field-content and mesh
//              digests and the component checkpoint bytes of v6.1 equal those of v6
//   timing     Field build/matches, mesh digest, prepare_patch: sealed v6 vs v6.1 (ms)
#include "hall-finite-foil-em-20261029/axial_surface_path_port.hpp"
#include "ss304_heat_material_v2.hpp"
#include "pic-electrode-heat-source-v6-20261030/native/mfem_electrode_heat_source_v6/source_owner.hpp"
#include "pic-electrode-heat-source-v6.1-20261031/native/mfem_electrode_heat_source_v6_1/source_owner.hpp"
#include "raikou/json.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <new>
#include <optional>
#include <random>

namespace eh = electrode_heat_v2;
namespace v3 = electrode_heat_v3;
namespace se = surface_emission;
namespace emi = electrode_heat_v6_1::emission;
using J = raikou::Json;
using V3 = se::Vector;

namespace {
struct Results {
  J cases = J::array();
  int failed = 0;
  void record(const std::string& g, const std::string& name, bool ok, J ev = J::object()) {
    cases.push_back({{"group", g}, {"name", name}, {"passed", ok}, {"evidence", std::move(ev)}});
    if (!ok) ++failed;
    std::cout << g << " | " << name << ": " << (ok ? "pass" : "FAIL") << std::endl;
    if (!ok) std::cout << "  evidence: " << cases.back()["evidence"].dump().substr(0, 3000) << std::endl;
  }
};
double since_ms(std::chrono::steady_clock::time_point a) { return 1e3 * std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count(); }
template <class F> double median_ms(F&& f, int reps) {
  std::vector<double> v;
  for (int r = 0; r < reps; ++r) { const auto t = std::chrono::steady_clock::now(); f(); v.push_back(since_ms(t)); }
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}
std::string hexs(double x) { std::ostringstream s; s << std::hexfloat << x; return s.str(); }

// ---------------------------------------------------------------- emission sweep
void emission_group(Results& R) {
  std::vector<se::Config> configs = {se::Config{4.3}, se::Config{2.0}, se::Config{6.0}, se::Config{4.5}, se::Config{3.1, 96, 48, 40}, se::Config{4.3, 16, 16, 32}};
  std::vector<double> temps = {0., 1., 77., 300., 300.16612345, 500., 1000., 1500., 2500., 3000.};
  std::vector<double> fields = {-5e9, -1e6, -2e4, -1., -0., 0., 1., 10., 256.0807, 1e3, 5503.48, 1e4, 1e5, 1e6, 3e6, 1e7, 3e7, 1e8, 3e8, 1e9, 3e9, 5e9, .5, 5e9 * 1.01};
  std::mt19937_64 rng(20261031);
  std::uniform_real_distribution<double> u(0, 1);
  for (int k = 0; k < 120; ++k) {
    temps.push_back(1 + 2999 * u(rng));
    fields.push_back((u(rng) < .5 ? -1 : 1) * std::exp(std::log(1.) + (std::log(5e9) - std::log(1.)) * u(rng)));
  }
  long compared = 0, identical = 0, refusals = 0, same_refusals = 0;
  double t_sealed = 0, t_fast = 0;
  J diffs = J::array();
  emi::reset_statistics();
  for (const auto& c : configs)
    for (std::size_t ti = 0; ti < temps.size(); ++ti)
      for (std::size_t fi = 0; fi < fields.size(); ++fi) {
        // Full grid for the fixed temperatures; random temperatures pair with random fields.
        if (ti >= 10 && fi >= 24 && ti - 10 != fi - 24) continue;
        if (ti >= 10 && fi < 24 && fi % 4) continue;
        const double T = temps[ti], F = fields[fi];
        std::string ws, wf;
        se::Distribution a, b;
        auto t0 = std::chrono::steady_clock::now();
        try { a = se::evaluate(T, F, c); } catch (const std::exception& e) { ws = e.what(); }
        t_sealed += since_ms(t0);
        t0 = std::chrono::steady_clock::now();
        try { b = emi::evaluate(T, F, c); } catch (const std::exception& e) { wf = e.what(); }
        t_fast += since_ms(t0);
        if (!ws.empty() || !wf.empty()) {
          ++refusals;
          if (ws == wf) ++same_refusals;
          else if (diffs.size() < 20) diffs.push_back({{"T", hexs(T)}, {"F", hexs(F)}, {"phi", c.work_function_eV}, {"sealed", ws}, {"fast", wf}});
          continue;
        }
        ++compared;
        if (emi::same_distribution(a, b)) ++identical;
        else if (diffs.size() < 20) diffs.push_back({{"T", hexs(T)}, {"F", hexs(F)}, {"phi", c.work_function_eV}, {"sealed_log_current", hexs(a.log_current_A_m2)},
                                                     {"fast_log_current", hexs(b.log_current_A_m2)}, {"nodes", {a.nodes.size(), b.nodes.size()}}});
      }
  // Repeated calls (retarding reuse) stay identical.
  bool repeat = true;
  for (int k = 0; k < 3; ++k) {
    const auto a = se::evaluate(300, -123.25, se::Config{4.3});
    const auto b = emi::evaluate(300, -123.25, se::Config{4.3});
    repeat = repeat && emi::same_distribution(a, b);
  }
  const auto st = emi::statistics();
  J ev = {{"inputs_compared", compared}, {"bit_identical", identical}, {"refused_inputs", refusals}, {"identical_refusals", same_refusals},
          {"differences", diffs}, {"sealed_total_ms", t_sealed}, {"fast_total_ms", t_fast},
          {"statistics", {{"extracting_fast", st.extracting_fast}, {"retarding_reused", st.retarding_reused}, {"retarding_sealed", st.retarding_sealed},
                          {"sealed_calls", st.sealed_calls}, {"nodes", st.nodes}, {"nodes_exact", st.nodes_exact}, {"nodes_bounded_zero", st.nodes_bounded_zero},
                          {"self_checks", st.self_checks}, {"sealed_fallback_configs", st.sealed_fallback_configs}}}};
  R.record("emission", "fast spectrum bit-identical to the sealed evaluate on every admitted input", compared > 1000 && identical == compared, ev);
  R.record("emission", "identical refusals (same message) on every refused input", refusals > 0 && same_refusals == refusals, {{"refused", refusals}});
  R.record("emission", "repeated retarding evaluations reuse the first and stay identical", repeat && st.retarding_reused > 0);
  R.record("emission", "no configuration fell back to the sealed path; nodes were dropped only as exact zeros", st.sealed_fallback_configs == 0 && st.nodes_bounded_zero > 0);
}

// ---------------------------------------------------------------- owner API traits
struct V6 {
  using Owner = electrode_heat_v6::Owner; using Field = electrode_heat_v6::Field; using View = electrode_heat_v6::FinalFieldView;
  using Candidate = electrode_heat_v6::Candidate; using Trial = electrode_heat_v6::Trial; using PreparedPatch = electrode_heat_v6::PreparedPatch;
  using Settings = electrode_heat_v6::SourceSettings; using Identity = electrode_heat_v6::PatchIdentity; using Particle = electrode_heat_v6::Particle;
  using Refusal = electrode_heat_v6::Refusal; using Leaf = electrode_heat_v6::NeutralLeaf;
  static constexpr const char* name = "v6";
};
struct V61 {
  using Owner = electrode_heat_v6_1::Owner; using Field = electrode_heat_v6_1::Field; using View = electrode_heat_v6_1::FinalFieldView;
  using Candidate = electrode_heat_v6_1::Candidate; using Trial = electrode_heat_v6_1::Trial; using PreparedPatch = electrode_heat_v6_1::PreparedPatch;
  using Settings = electrode_heat_v6_1::SourceSettings; using Identity = electrode_heat_v6_1::PatchIdentity; using Particle = electrode_heat_v6_1::Particle;
  using Refusal = electrode_heat_v6_1::Refusal; using Leaf = electrode_heat_v6_1::NeutralLeaf;
  static constexpr const char* name = "v6.1";
};

// ---- audit fixture (source v5 audit port inputs, unchanged)
eh::Geometry geom() { eh::Geometry g; g.nx = 4; g.ny = 2; g.revision = "thermal-mesh-v3"; g.electrode = "upper"; return g; }
eh::Material material() {
  return eh::conditional_ss304_v2(ss304v2::Qualification::source_conditioned_liquid_3000, ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation, 8000,
                                  eh::material_headers_sha256);
}
eh::Patch patch() { eh::Patch p; p.id = "lower-face"; p.face = eh::Face::LowY; p.source_hash = "v3-source-scenario"; return p; }
eh::Patch heat_only(eh::Patch p) { p.desorption_enabled = false; p.coverage_temkin = false; p.initial_particles = p.initial_low_particles = 0; return p; }
eh::Patch water_patch(double high, double low, double rate) {
  auto p = patch(); p.desorption_enabled = true; p.state_id = "H2O"; p.parent_id = "H2O"; p.initial_particles = high; p.initial_low_particles = low;
  p.mass_kg = 18 * 1.66053906660e-27; p.attempt_s = rate; p.activation_eV = p.desorption_eV = p.internal_eV = 0; return p;
}
eh::SourceEventIdentity event(std::uint64_t seq) { return {"source-v3", "e-", "electrode-electron-reservoir", "native-layer", "canonical-electron-42", seq, seq}; }
const V3 origin{.0001, .0025, .0002};
se::LayerPolicy policy() { se::LayerPolicy p; p.local_linear_vacuum_extent_m = 25e-6; return p; }
struct Arrays { foil::Vec E, B; };
Arrays arrays(const axial::Mesh& g, double F, double B, bool affine = false) {
  Arrays a{foil::Vec(g.ne), foil::Vec(g.nb)};
  std::vector<V3> nodes(g.nn);
  for (int k = 0; k <= g.nz; ++k) for (int j = 0; j <= g.ny; ++j) for (int i = 0; i <= g.nx; ++i) nodes[g.node(i, j, k)] = {i * g.dx - g.side / 2, j * g.dy - g.side / 2, g.zmin + k * g.dz};
  for (int e = 0; e < g.ne; ++e) {
    const auto& edge = g.edges[e];
    V3 x{};
    for (int d = 0; d < 3; ++d) x[d] = .5 * (nodes[edge.a][d] + nodes[edge.b][d]);
    if (edge.axis == 0) a.E[e] = affine ? 1e6 + 1e13 * (x[1] - origin[1]) : 0;
    if (edge.axis == 1) a.E[e] = F + (affine ? 1e13 * (x[0] - origin[0]) : 0);
    if (edge.axis == 2) a.E[e] = affine ? 2e6 + 1e13 * (x[1] - origin[1]) : 0;
  }
  for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i <= g.nx; ++i) a.B[g.bx(i, j, k)] = B + (affine ? 1e6 * (i * g.dx - g.side / 2 - origin[0]) : 0);
  for (int k = 0; k < g.nz; ++k) for (int j = 0; j <= g.ny; ++j) for (int i = 0; i < g.nx; ++i) a.B[g.by(i, j, k)] = affine ? -1e6 * (j * g.dy - g.side / 2 - origin[1]) : 0;
  return a;
}
// The heat-only thermal view of the patch configuration (the v6 one-thermal-owner binding).
std::unique_ptr<eh::Coupler> thermal_view(const eh::Patch& p, double T) { return std::make_unique<eh::Coupler>(geom(), material(), std::vector<eh::Patch>{heat_only(p)}, T, "v3-source-scenario", 5e-9); }
template <class A> typename A::Settings settings(se::Config c = se::Config{4.3}, double w0 = 0) { return {c, policy(), "audit-port-caller-authored-marks-no-arrival-policy", w0}; }
template <class A> std::unique_ptr<typename A::Owner> owner(const eh::Patch& p, double T, const eh::Coupler& view, se::Config c = se::Config{4.3}, double w0 = 0) {
  return std::make_unique<typename A::Owner>(geom(), material(), std::vector<eh::Patch>{p}, std::vector<typename A::Identity>{{"upper-foil-lower-surface", "conditional-single-layer", .005 * .025}},
                                             std::vector<typename A::Settings>{settings<A>(c, w0)}, T, "v3-source-scenario", 5e-9, view);
}
template <class O> v3::TrialFrame frame(const O& h, std::uint64_t stage = 1, double dt = 1e-14) {
  v3::TrialFrame f; f.role = v3::FieldRole::TotalMidpoint; f.accepted_revision = h.thermal_state().revision; f.trial_stage = stage; f.mesh_revision = 1; f.material_revision = 2;
  f.begin_s = h.thermal_state().time; f.end_s = f.begin_s + dt; f.dt_s = dt; f.field_stage_s = f.begin_s + .5 * dt;
  f.arena_hash = "accepted-arena-revision-" + std::to_string(h.thermal_state().revision); f.field_hash = "immutable-midpoint-trial-" + std::to_string(stage);
  f.mesh_hash = "Hall-axial-sealed-mesh"; f.material_hash = "global-material-v2-map"; f.species_hash = "actual-source-v1-electron"; return f;
}
template <class O> v3::PatchBinding binding(const O& h) {
  v3::PatchBinding b; b.thermal_revision = h.thermal_state().revision; b.source_accepted_state = h.thermal_state().revision + 1; b.patch_id = 7;
  b.thermal_mesh_identity = geom().revision; b.thermal_material_hash = eh::material_headers_sha256;
  b.material_phase_admitted = b.planar_barrier_admitted = b.electron_supply_admitted = true; b.temperature = v3::HeldTemperature::AcceptedPatchMean;
  b.extrusion_center_z_m = 0; b.extrusion_height_m = geom().height; return b;
}
template <class A> struct ParentField { std::shared_ptr<const typename A::Field> f; Arrays a; };
template <class A> ParentField<A> field(const axial::Mesh& g, const typename A::Owner& h, std::uint64_t stage, double F, double B = 0, bool affine = false) {
  auto a = arrays(g, F, B, affine);
  auto f = std::make_shared<const typename A::Field>(g, frame(h, stage), a.E, a.B);
  return {f, std::move(a)};
}
template <class A> typename A::View view(const ParentField<A>& p, const axial::Mesh& g) { return {p.f->frame(), g, p.a.E, p.a.B}; }
template <class A> typename A::PreparedPatch prepared(const typename A::Owner& h, const std::shared_ptr<const typename A::Field>& f, double area = 1) {
  return h.prepare_patch(h.make_epoch(f->frame().dt_s, "trial-midpoint-field"), 0, binding(h), f, origin, v3::Exposure{f->frame().begin_s, f->frame().dt_s, area, "native-surface-quadrature"});
}
template <class A> v3::BirthSample sample(const typename A::PreparedPatch& p, const typename A::Owner& h, double weight = 1e9, bool qualification = true) {
  v3::BirthSample b; b.physical_weight = weight; b.normal_uniform = .6; b.transverse_uniform = .35; b.azimuth_uniform = .1; b.time_uniform = .23; b.event_id = 123;
  b.canonical_particle_id = 42; b.count_owner_revision = h.thermal_state().revision; b.count_owner_base_hash = "unchanged-accepted-source-count-owner";
  b.count_owner_token = qualification ? "qualification-weight-not-an-event-rate" : "staged-canonical-fractional-count-receipt"; b.exposure_key = p.exposure_key();
  b.time_policy = v3::BirthTime::ExposureBeginLumped; b.weight_provenance = qualification ? v3::WeightProvenance::QualificationOnly : v3::WeightProvenance::StagedCanonicalOwner;
  return b;
}
template <class A> typename A::Particle canonical_electron(double w) { typename A::Particle c; c.id = 42; c.mass_kg = se::me; c.charge_C = -se::qe; c.weight = w; return c; }
struct Outcome { bool threw = false; std::string what; };
template <class F> Outcome attempt(F&& f) { Outcome o; try { f(); } catch (const std::exception& e) { o.threw = true; o.what = e.what(); } return o; }
J rec(const Outcome& o) { return {{"refused", o.threw}, {"message", o.what.substr(0, 300)}}; }

// ---- A02 probes (verdicts and recorded observations)
template <class A> J a02_weight(const axial::Mesh& g) {
  auto hv = thermal_view(patch(), 1500);
  auto o = owner<A>(patch(), 1500, *hv);
  auto pf = field<A>(g, *o, 1, 1e8);
  auto v = view<A>(pf, g);
  auto p = prepared<A>(*o, pf.f);
  auto b = o->begin_source(p, sample<A>(p, *o, 1e9, false), event(1));
  const bool genuine = o->preflight(b, pf.f, v);
  auto io = attempt([&] { (void)o->import(pf.f, v, std::span<const typename A::Candidate>(&b, 1), {42}); });
  auto cv = thermal_view(patch(), 300);
  auto cold = owner<A>(patch(), 300, *cv);
  auto cf = field<A>(g, *cold, 2, 0);
  auto cvw = view<A>(cf, g);
  auto cp = prepared<A>(*cold, cf.f);
  auto thermal = cold->begin_source(cp, sample<A>(cp, *cold, 1, false), event(1));
  const auto base = cold->state_sha256();
  auto ho = attempt([&] { (void)cold->handoff(thermal, cvw, canonical_electron<A>(2)); });
  const bool ok = genuine && !io.threw && ho.threw && cold->state_sha256() == base;
  return {{"genuine_candidate_preflight", genuine}, {"genuine_import", rec(io)}, {"handoff_with_doubled_canonical_weight", rec(ho)},
          {"owner_unchanged", cold->state_sha256() == base}, {"verdict", ok ? "runtime analogue refused" : "still-defective"}};
}
template <class A> J a02_spectrum(const axial::Mesh& g) {
  auto v15 = thermal_view(patch(), 1500), v16 = thermal_view(patch(), 1600), v15b = thermal_view(patch(), 1500);
  auto o15 = owner<A>(patch(), 1500, *v15), o16 = owner<A>(patch(), 1600, *v16);
  auto f15 = field<A>(g, *o15, 1, 1e8), f16 = field<A>(g, *o16, 1, 1e8);
  auto p15 = prepared<A>(*o15, f15.f), p16 = prepared<A>(*o16, f16.f);
  const auto b15 = o15->state_sha256(), b16 = o16->state_sha256();
  auto old_key = sample<A>(p15, *o15, 1e9, false);
  auto r1 = attempt([&] { (void)o16->begin_source(p16, old_key, event(1)); });
  auto r2 = attempt([&] { (void)o15->begin_source(p16, sample<A>(p16, *o16, 1e9, false), event(1)); });
  auto ctrl = attempt([&] { (void)o16->begin_source(p16, sample<A>(p16, *o16, 1e9, false), event(1)); });
  auto o43 = owner<A>(patch(), 1500, *v15b, se::Config{4.5});
  auto f43 = field<A>(g, *o43, 1, 1e8);
  auto p43 = prepared<A>(*o43, f43.f);
  auto r3 = attempt([&] { (void)o43->begin_source(p43, old_key, event(1)); });
  const bool unchanged = o15->state_sha256() == b15 && o16->state_sha256() == b16;
  const bool ok = p15.exposure_key() != p16.exposure_key() && p43.exposure_key() != p15.exposure_key() && r1.threw && r2.threw && r3.threw && !ctrl.threw && unchanged;
  return {{"exposure_keys_differ_1500_vs_1600", p15.exposure_key() != p16.exposure_key()}, {"1600K_patch_with_1500K_count_receipt", rec(r1)},
          {"1600K_patch_presented_to_1500K_owner", rec(r2)}, {"control_matching_1600K_receipt", rec(ctrl)}, {"4p5eV_patch_with_4p3eV_count_receipt", rec(r3)},
          {"owners_unchanged", unchanged}, {"verdict", ok ? "runtime analogues refused" : "still-defective"}};
}
template <class A> J a02_aba(const axial::Mesh& g) {
  auto hv = thermal_view(patch(), 1500);
  alignas(typename A::Owner) std::byte storage[sizeof(typename A::Owner)];
  auto make = [&] {
    return new (storage) typename A::Owner(geom(), material(), {patch()}, {{"upper-foil-lower-surface", "conditional-single-layer", .005 * .025}}, {settings<A>()}, 1500,
                                           "v3-source-scenario", 5e-9, *hv);
  };
  auto* old = make();
  auto fa = field<A>(g, *old, 3, 1e8);
  auto va = view<A>(fa, g);
  auto pa = prepared<A>(*old, fa.f);
  auto ca = old->begin_source(pa, sample<A>(pa, *old, 1, false), event(1));
  const bool before = old->preflight(ca, fa.f, va);
  const auto address = reinterpret_cast<std::uintptr_t>(old);
  std::destroy_at(old);
  auto* repl = make();
  const auto base = repl->state_sha256();
  const bool same_address = reinterpret_cast<std::uintptr_t>(repl) == address;
  const bool accepted = repl->preflight(ca, fa.f, va);
  auto r1 = attempt([&] { (void)repl->begin_source(pa, sample<A>(pa, *repl, 1, false), event(1)); });
  auto r2 = attempt([&] { (void)repl->prepare(repl->make_epoch(1e-14, "trial-midpoint-field"), fa.f, va, std::span<const typename A::Candidate>(&ca, 1), {}, {{}}); });
  auto r3 = attempt([&] { (void)repl->import(fa.f, va, std::span<const typename A::Candidate>(&ca, 1), {42}); });
  const bool unchanged = repl->state_sha256() == base;
  std::destroy_at(repl);
  const bool ok = before && !accepted && r1.threw && r2.threw && r3.threw && unchanged;
  return {{"same_address", same_address}, {"predecessor_preflight_before_destruction", before}, {"replacement_accepts_stale_candidate", accepted},
          {"replacement_begin_source_with_stale_patch", rec(r1)}, {"replacement_prepare_with_stale_candidate", rec(r2)},
          {"replacement_import_with_stale_candidate", rec(r3)}, {"replacement_owner_unchanged", unchanged}, {"verdict", ok ? "refused" : "still-defective"}};
}
template <class A> J a02_mesh() {
  auto hv = thermal_view(patch(), 1500);
  auto o = owner<A>(patch(), 1500, *hv);
  axial::Mesh gm(45, 4, .009, 1);
  auto fm = field<A>(gm, *o, 4, 1e8, 0, true);
  auto fu = field<A>(gm, *o, 5, 1e8);
  auto vm = view<A>(fm, gm);
  auto vu = view<A>(fu, gm);
  const auto before = fm.f->surface(origin, {0, -1, 0});
  auto pu = prepared<A>(*o, fu.f);
  auto cu = o->begin_source(pu, sample<A>(pu, *o, 1e9, false), event(1));
  auto e = o->make_epoch(1e-14, "trial-midpoint-field");
  std::optional<typename A::Trial> t;
  t.emplace(o->prepare(e, fm.f, vm, {}, {}, {{}}));
  J j = {{"pre_mutation_field_matches", fm.f->matches(vm)}, {"pre_mutation_candidate_preflight", o->preflight(cu, fu.f, vu)}, {"pre_mutation_trial_preflight", o->preflight(*t, vm)}};
  const auto base = o->state_sha256();
  const double dx = gm.dx;
  gm.dx *= 1.001;
  const auto after = fm.f->surface(origin, {0, -1, 0});
  j["sealed_field_gather_unchanged_by_alias_mutation"] = after.E == before.E && after.B == before.B;
  const bool m1 = fm.f->matches(vm), c1 = o->preflight(cu, fu.f, vu), t1 = o->preflight(*t, vm);
  j["mutated_mesh_field_matches"] = m1; j["mutated_mesh_candidate_preflight"] = c1; j["mutated_mesh_trial_preflight"] = t1;
  auto io = attempt([&] { (void)o->import(fu.f, vu, std::span<const typename A::Candidate>(&cu, 1), {42}); });
  auto ho = attempt([&] { (void)o->handoff(cu, vu, canonical_electron<A>(1e9)); });
  auto po = attempt([&] { (void)o->prepare(e, fm.f, vm, {}, {}, {{}}); });
  j["mutated_mesh_import"] = rec(io); j["mutated_mesh_handoff"] = rec(ho); j["mutated_mesh_new_prepare"] = rec(po);
  const bool committed_mutated = o->commit(std::move(*t), vm);
  j["mutated_mesh_commit"] = committed_mutated; j["owner_unchanged_after_refusals"] = o->state_sha256() == base;
  gm.dx = dx;
  const double zmin = gm.zmin; gm.zmin += 1e-9; const bool mz = fm.f->matches(vm); gm.zmin = zmin;
  const double ve = gm.VE[0]; gm.VE[0] *= 2; const bool mv = fm.f->matches(vm); gm.VE[0] = ve;
  const double len = gm.edges[0].length; gm.edges[0].length *= 1 + 1e-15; const bool ml = fm.f->matches(vm); gm.edges[0].length = len;
  // v6.1 additions: the last entries of the arrays and a bool member are bound as well.
  const double vn = gm.VN.back(); gm.VN.back() = std::nextafter(vn, INFINITY); const bool mn = fm.f->matches(vm); gm.VN.back() = vn;
  const bool on = gm.short_node.back(); gm.short_node.back() = !on; const bool mb = fm.f->matches(vm); gm.short_node.back() = on;
  j["mutated_zmin_matches"] = mz; j["mutated_VE0_matches"] = mv; j["mutated_edge0_length_matches"] = ml; j["mutated_VN_last_matches"] = mn; j["mutated_short_node_last_matches"] = mb;
  const bool restored = fm.f->matches(vm);
  const bool committed = o->commit(std::move(*t), vm);
  j["restored_mesh_matches"] = restored; j["restored_mesh_commit_positive_control"] = committed;
  j["predecessor_candidate_stale_after_commit"] = !o->preflight(cu, fu.f, vu);
  const bool ok = !m1 && !c1 && !t1 && io.threw && ho.threw && po.threw && !committed_mutated && !mz && !mv && !ml && !mn && !mb && restored && committed;
  j["verdict"] = ok ? "refused (mutation detected, owner unchanged; positive control commits)" : "still-defective";
  return j;
}
template <class A> J ext_publish_gap() {
  auto hv = thermal_view(patch(), 1500);
  auto o = owner<A>(patch(), 1500, *hv);
  axial::Mesh g(45, 4, .009, 1);
  auto pf = field<A>(g, *o, 1, 1e8);
  auto v = view<A>(pf, g);
  auto e = o->make_epoch(1e-14, "trial-midpoint-field");
  std::optional<typename A::Trial> t;
  t.emplace(o->prepare(e, pf.f, v, {}, {}, {{}}));
  const bool pre = o->preflight(*t, v);
  pf.a.E[0] += 1;
  const bool still = pf.f->matches(v);
  const auto base = o->state_sha256();
  o->publish_prevalidated(std::move(*t));
  const bool published = o->state_sha256() != base;
  return {{"preflight_before_mutation", pre}, {"field_matches_after_parent_array_mutation", still}, {"publish_prevalidated_published", published},
          {"verdict", published && !still ? "observation: undetected after preflight (root critical-section contract, not a component check)" : "observation: detected"}};
}
void a02_group(Results& R) {
  axial::Mesh g(45, 4, .009, 1);
  const std::vector<std::pair<std::string, std::pair<J, J>>> probes = {
      {"A02-C3-changed-particle-weight", {a02_weight<V6>(g), a02_weight<V61>(g)}},
      {"A02-C3-changed-spectrum", {a02_spectrum<V6>(g), a02_spectrum<V61>(g)}},
      {"A02-C3-reused-owner-address", {a02_aba<V6>(g), a02_aba<V61>(g)}},
      {"A02-C3-mutated-mesh-old-metadata", {a02_mesh<V6>(), a02_mesh<V61>()}},
      {"EXT-A02-mutation-between-preflight-and-publish", {ext_publish_gap<V6>(), ext_publish_gap<V61>()}}};
  for (const auto& [id, pr] : probes) {
    const auto& [a, b] = pr;
    const std::string va = a["verdict"], vb = b["verdict"];
    const bool refused = va.rfind("observation", 0) == 0 || vb.find("defective") == std::string::npos;
    // Same verdict and same boolean observations (messages may differ only in text that names the version).
    bool same_booleans = true;
    for (auto it = a.begin(); it != a.end(); ++it)
      if (it->is_boolean()) same_booleans = same_booleans && b.contains(it.key()) && b[it.key()] == *it;
      else if (it->is_object() && it->contains("refused")) same_booleans = same_booleans && b.contains(it.key()) && b[it.key()]["refused"] == (*it)["refused"];
    R.record("a02", id + ": v6.1 verdict and every observation equal v6", va == vb && refused && same_booleans, {{"v6", a}, {"v6.1", b}});
  }
}

// ---------------------------------------------------------------- digest identity v6 vs v6.1
std::string read_all(const std::filesystem::path& p) { std::ifstream f(p, std::ios::binary); return {std::istreambuf_iterator<char>(f), {}}; }
void identity_group(Results& R, const std::filesystem::path& out) {
  axial::Mesh g(45, 4, .009, 1);
  auto hv = thermal_view(patch(), 1500);
  auto a = owner<V6>(patch(), 1500, *hv);
  auto b = owner<V61>(patch(), 1500, *hv);
  J ev;
  ev["configuration"] = a->configuration_sha256() == b->configuration_sha256();
  ev["state"] = a->state_sha256() == b->state_sha256();
  auto fa = field<V6>(g, *a, 1, 1e8, 0, true);
  auto fb = field<V61>(g, *b, 1, 1e8, 0, true);
  ev["field_content"] = fa.f->content_sha256() == fb.f->content_sha256();
  ev["field_mesh"] = fa.f->mesh_sha256() == fb.f->mesh_sha256() && fb.f->mesh_sha256() == electrode_heat_v6::detail::mesh_hash(g);
  auto pa = prepared<V6>(*a, fa.f);
  auto pb = prepared<V61>(*b, fb.f);
  ev["prepared_patch_exposure_key"] = pa.exposure_key() == pb.exposure_key();
  ev["prepared_patch_expected_number"] = std::bit_cast<std::uint64_t>(pa.expected_number()) == std::bit_cast<std::uint64_t>(pb.expected_number());
  auto ca = a->begin_source(pa, sample<V6>(pa, *a, 1e9, false), event(1));
  auto cb = b->begin_source(pb, sample<V61>(pb, *b, 1e9, false), event(1));
  ev["candidate_content"] = ca.content_sha256() == cb.content_sha256();
  auto va = view<V6>(fa, g);
  auto vb = view<V61>(fb, g);
  auto e = a->make_epoch(1e-14, "trial-midpoint-field");
  // Continued births (Active at the common endpoint, or handed off / recollected) in the trial.
  auto da = a->advance_birth(pa, ca, event(2));
  auto db = b->advance_birth(pb, cb, event(2));
  ev["continued_candidate_content"] = da.content_sha256() == db.content_sha256();
  const bool active = da.state().native.status == se::LayerStatus::Active;
  auto ta = active ? a->prepare(e, fa.f, va, std::span<const V6::Candidate>(&da, 1), {}, {{}}) : a->prepare(e, fa.f, va, {}, {}, {{}});
  auto tb = active ? b->prepare(e, fb.f, vb, std::span<const V61::Candidate>(&db, 1), {}, {{}}) : b->prepare(e, fb.f, vb, {}, {}, {{}});
  ev["trial_content"] = ta.content_sha256() == tb.content_sha256();
  ev["trial_preflight"] = a->preflight(ta, va) && b->preflight(tb, vb);
  a->publish_prevalidated(std::move(ta));
  b->publish_prevalidated(std::move(tb));
  ev["state_after_publication"] = a->state_sha256() == b->state_sha256() && b->state_sha256() == b->state_sha256_uncached();
  a->save((out / "checkpoint-v6.bin").string());
  b->save((out / "checkpoint-v61.bin").string());
  ev["component_checkpoint_bytes"] = read_all(out / "checkpoint-v6.bin") == read_all(out / "checkpoint-v61.bin");
  // Quantized neutral release (desorbing patch) trial with leaves.
  auto wp = water_patch(1e15, 0, 1e13);
  auto wv = thermal_view(wp, 1500);
  auto wa = owner<V6>(wp, 1500, *wv, se::Config{4.3}, 2.);
  auto wb = owner<V61>(wp, 1500, *wv, se::Config{4.3}, 2.);
  auto ga = field<V6>(g, *wa, 1, 0);
  auto gb = field<V61>(g, *wb, 1, 0);
  auto ea = wa->make_epoch(1e-14, "trial-midpoint-field");
  std::vector<std::vector<V6::Leaf>> la(1);
  std::vector<std::vector<V61::Leaf>> lb(1);
  for (int k = 0; k < 5; ++k) { la[0].push_back({2., eh::Draw{.3, .4, .5, .6, .7, .1 + .1 * k}}); lb[0].push_back({2., eh::Draw{.3, .4, .5, .6, .7, .1 + .1 * k}}); }
  auto qa = wa->prepare(ea, ga.f, view<V6>(ga, g), {}, {}, la);
  auto qb = wb->prepare(ea, gb.f, view<V61>(gb, g), {}, {}, lb);
  ev["neutral_trial_content"] = qa.content_sha256() == qb.content_sha256();
  ev["release_quotes"] = wa->release_quotes(ea)[0].log_released == wb->release_quotes(ea)[0].log_released;
  bool all = true;
  for (auto it = ev.begin(); it != ev.end(); ++it) all = all && it->get<bool>();
  R.record("identity", "every v6.1 digest, trial content and checkpoint byte equals v6", all, ev);
}

// ---------------------------------------------------------------- timing
void timing_group(Results& R) {
  axial::Mesh g(45, 4, .009, 1);
  auto hv = thermal_view(patch(), 300);
  auto a = owner<V6>(patch(), 300, *hv);
  auto b = owner<V61>(patch(), 300, *hv);
  const auto arr = arrays(g, 2e6, 0, true);
  const auto f = frame(*a);
  std::shared_ptr<const V6::Field> fa;
  std::shared_ptr<const V61::Field> fb;
  J t;
  t["v6_field_build_ms"] = median_ms([&] { fa = std::make_shared<const V6::Field>(g, f, arr.E, arr.B); }, 5);
  t["v61_field_build_ms"] = median_ms([&] { fb = std::make_shared<const V61::Field>(g, f, arr.E, arr.B); }, 5);
  const V6::View va{fa->frame(), g, arr.E, arr.B};
  const V61::View vb{fb->frame(), g, arr.E, arr.B};
  bool ma = false, mb = false;
  t["v6_matches_ms"] = median_ms([&] { ma = fa->matches(va); }, 5);
  t["v61_matches_ms"] = median_ms([&] { mb = fb->matches(vb); }, 5);
  t["v6_mesh_hash_ms"] = median_ms([&] { (void)electrode_heat_v6::detail::mesh_hash(g); }, 5);
  t["v61_mesh_hash_ms"] = median_ms([&] { (void)electrode_heat_v6_1::detail::mesh_hash(g); }, 5);
  t["v6_prepare_patch_ms"] = median_ms([&] { (void)prepared<V6>(*a, fa); }, 5);
  t["v61_prepare_patch_ms"] = median_ms([&] { (void)prepared<V61>(*b, fb); }, 5);
  R.record("timing", "Field build, matches and prepare_patch timings (informational)", ma && mb, t);
}
}  // namespace

int main(int argc, char** argv) {
  Results R;
  int rc = 0;
  try {
    if (argc != 2) throw std::runtime_error("usage: test_source_v6_1 OUTDIR");
    const std::filesystem::path out = argv[1];
    std::filesystem::create_directories(out);
    emission_group(R);
    a02_group(R);
    identity_group(R, out);
    timing_group(R);
    std::ofstream(out / "source-v6.1.json") << J{{"schema", "source-v6.1-qualification-v1"}, {"cases", R.cases}, {"failed", R.failed}}.dump(1) << '\n';
    rc = R.failed == 0 ? 0 : 1;
  } catch (const std::exception& e) {
    std::cerr << "fatal: " << e.what() << '\n';
    rc = 2;
  }
  return rc;
}
