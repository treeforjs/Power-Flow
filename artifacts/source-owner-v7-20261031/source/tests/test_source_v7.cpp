// Heat/source owner v7 qualification against the sealed owner v6.1, in one program.
//   test_source_v7 OUTDIR
// Groups:
//   a02         the A02 audit probes (source v5 audit port rows, v6 API) on v6.1 and on v7 with every v7 option at
//               its v6.1 benchmark value: same verdicts and observations (ported from the v6.1 qualification)
//   identity    configuration, state, prepared-patch, candidate, trial, field-content and mesh digests, release
//               quotes and component checkpoint bytes of v7 (benchmark) equal those of v6.1
//   interval    the trial interval bound: default 1e-14 s (v6.1 decisions and texts), a configured 1e-13 s bound
//               (40 fs step published), refusals, and the error-target policy with its typed refusal and receipt
//   surface     work function from the exact stock, enhancement-class exposures, phase map, image-layer admission
//   births      matching-plane births (displacement, image tail, exact ledger) and sub-cell returns
//   wall        wall responses through the owner: exact ledgers, loads, stuck credits, ESD/ISD release via leaves of
//               w0, exact stock, coverage-dependent work function, staleness, checkpoint
#include "hall-finite-foil-em-20261029/axial_surface_path_port.hpp"
#include "ss304_heat_material_v2.hpp"
#include "pic-electrode-heat-source-v6.1-20261031/native/mfem_electrode_heat_source_v6_1/source_owner.hpp"
#include "pic-electrode-heat-source-v7-20261031/native/mfem_electrode_heat_source_v7/source_owner.hpp"
#include "event_set_v7.hpp"
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
using J = raikou::Json;
using V3 = se::Vector;
using Count = electrode_heat_v7::Count;

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
// ---------------------------------------------------------------- owner API traits
struct VA {
  using Owner = electrode_heat_v6_1::Owner; using Field = electrode_heat_v6_1::Field; using View = electrode_heat_v6_1::FinalFieldView;
  using Candidate = electrode_heat_v6_1::Candidate; using Trial = electrode_heat_v6_1::Trial; using PreparedPatch = electrode_heat_v6_1::PreparedPatch;
  using Settings = electrode_heat_v6_1::SourceSettings; using Identity = electrode_heat_v6_1::PatchIdentity; using Particle = electrode_heat_v6_1::Particle;
  using Refusal = electrode_heat_v6_1::Refusal; using Leaf = electrode_heat_v6_1::NeutralLeaf;
  static constexpr const char* name = "v6.1";
};
struct VB {
  using Owner = electrode_heat_v7::Owner; using Field = electrode_heat_v7::Field; using View = electrode_heat_v7::FinalFieldView;
  using Candidate = electrode_heat_v7::Candidate; using Trial = electrode_heat_v7::Trial; using PreparedPatch = electrode_heat_v7::PreparedPatch;
  using Settings = electrode_heat_v7::SourceSettings; using Identity = electrode_heat_v7::PatchIdentity; using Particle = electrode_heat_v7::Particle;
  using Refusal = electrode_heat_v7::Refusal; using Leaf = electrode_heat_v7::NeutralLeaf;
  static constexpr const char* name = "v7";
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
      {"A02-C3-changed-particle-weight", {a02_weight<VA>(g), a02_weight<VB>(g)}},
      {"A02-C3-changed-spectrum", {a02_spectrum<VA>(g), a02_spectrum<VB>(g)}},
      {"A02-C3-reused-owner-address", {a02_aba<VA>(g), a02_aba<VB>(g)}},
      {"A02-C3-mutated-mesh-old-metadata", {a02_mesh<VA>(), a02_mesh<VB>()}},
      {"EXT-A02-mutation-between-preflight-and-publish", {ext_publish_gap<VA>(), ext_publish_gap<VB>()}}};
  for (const auto& [id, pr] : probes) {
    const auto& [a, b] = pr;
    const std::string va = a["verdict"], vb = b["verdict"];
    const bool refused = va.rfind("observation", 0) == 0 || vb.find("defective") == std::string::npos;
    // Same verdict and same boolean observations (messages may differ only in text that names the version).
    bool same_booleans = true;
    for (auto it = a.begin(); it != a.end(); ++it)
      if (it->is_boolean()) same_booleans = same_booleans && b.contains(it.key()) && b[it.key()] == *it;
      else if (it->is_object() && it->contains("refused")) same_booleans = same_booleans && b.contains(it.key()) && b[it.key()]["refused"] == (*it)["refused"];
    R.record("a02", id + ": v7 (benchmark) verdict and every observation equal v6.1", va == vb && refused && same_booleans, {{"v6.1", a}, {"v7", b}});
  }
}

// ---------------------------------------------------------------- digest identity v6.1 vs v7 (benchmark)
std::string read_all(const std::filesystem::path& p) { std::ifstream f(p, std::ios::binary); return {std::istreambuf_iterator<char>(f), {}}; }
void identity_group(Results& R, const std::filesystem::path& out) {
  axial::Mesh g(45, 4, .009, 1);
  auto hv = thermal_view(patch(), 1500);
  auto a = owner<VA>(patch(), 1500, *hv);
  auto b = owner<VB>(patch(), 1500, *hv);
  J ev;
  ev["configuration"] = a->configuration_sha256() == b->configuration_sha256();
  ev["state"] = a->state_sha256() == b->state_sha256();
  auto fa = field<VA>(g, *a, 1, 1e8, 0, true);
  auto fb = field<VB>(g, *b, 1, 1e8, 0, true);
  ev["field_content"] = fa.f->content_sha256() == fb.f->content_sha256();
  ev["field_mesh"] = fa.f->mesh_sha256() == fb.f->mesh_sha256() && fb.f->mesh_sha256() == electrode_heat_v6_1::detail::mesh_hash(g);
  auto pa = prepared<VA>(*a, fa.f);
  auto pb = prepared<VB>(*b, fb.f);
  ev["prepared_patch_exposure_key"] = pa.exposure_key() == pb.exposure_key();
  ev["prepared_patch_expected_number"] = std::bit_cast<std::uint64_t>(pa.expected_number()) == std::bit_cast<std::uint64_t>(pb.expected_number());
  auto ca = a->begin_source(pa, sample<VA>(pa, *a, 1e9, false), event(1));
  auto cb = b->begin_source(pb, sample<VB>(pb, *b, 1e9, false), event(1));
  ev["candidate_content"] = ca.content_sha256() == cb.content_sha256();
  auto va = view<VA>(fa, g);
  auto vb = view<VB>(fb, g);
  auto e = a->make_epoch(1e-14, "trial-midpoint-field");
  // Continued births (Active at the common endpoint, or handed off / recollected) in the trial.
  auto da = a->advance_birth(pa, ca, event(2));
  auto db = b->advance_birth(pb, cb, event(2));
  ev["continued_candidate_content"] = da.content_sha256() == db.content_sha256();
  const bool active = da.state().native.status == se::LayerStatus::Active;
  auto ta = active ? a->prepare(e, fa.f, va, std::span<const VA::Candidate>(&da, 1), {}, {{}}) : a->prepare(e, fa.f, va, {}, {}, {{}});
  auto tb = active ? b->prepare(e, fb.f, vb, std::span<const VB::Candidate>(&db, 1), {}, {{}}) : b->prepare(e, fb.f, vb, {}, {}, {{}});
  ev["trial_content"] = ta.content_sha256() == tb.content_sha256();
  ev["trial_preflight"] = a->preflight(ta, va) && b->preflight(tb, vb);
  a->publish_prevalidated(std::move(ta));
  b->publish_prevalidated(std::move(tb));
  ev["state_after_publication"] = a->state_sha256() == b->state_sha256() && b->state_sha256() == b->state_sha256_uncached();
  a->save((out / "checkpoint-v61.bin").string());
  b->save((out / "checkpoint-v7.bin").string());
  ev["component_checkpoint_bytes"] = read_all(out / "checkpoint-v61.bin") == read_all(out / "checkpoint-v7.bin");
  // Quantized neutral release (desorbing patch) trial with leaves.
  auto wp = water_patch(1e15, 0, 1e13);
  auto wv = thermal_view(wp, 1500);
  auto wa = owner<VA>(wp, 1500, *wv, se::Config{4.3}, 2.);
  auto wb = owner<VB>(wp, 1500, *wv, se::Config{4.3}, 2.);
  auto ga = field<VA>(g, *wa, 1, 0);
  auto gb = field<VB>(g, *wb, 1, 0);
  auto ea = wa->make_epoch(1e-14, "trial-midpoint-field");
  std::vector<std::vector<VA::Leaf>> la(1);
  std::vector<std::vector<VB::Leaf>> lb(1);
  for (int k = 0; k < 5; ++k) { la[0].push_back({2., eh::Draw{.3, .4, .5, .6, .7, .1 + .1 * k}}); lb[0].push_back({2., eh::Draw{.3, .4, .5, .6, .7, .1 + .1 * k}}); }
  auto qa = wa->prepare(ea, ga.f, view<VA>(ga, g), {}, {}, la);
  auto qb = wb->prepare(ea, gb.f, view<VB>(gb, g), {}, {}, lb);
  ev["neutral_trial_content"] = qa.content_sha256() == qb.content_sha256();
  ev["release_quotes"] = wa->release_quotes(ea)[0].log_released == wb->release_quotes(ea)[0].log_released;
  bool all = true;
  for (auto it = ev.begin(); it != ev.end(); ++it) all = all && it->get<bool>();
  R.record("identity", "every v7 (benchmark) digest, trial content and checkpoint byte equals v6.1", all, ev);
}

// ---------------------------------------------------------------- v7: trial interval bound
namespace s7 = electrode_heat_v7;
namespace sv = surface_v7;
template <class F> std::string refusal_text(F&& f) { try { f(); } catch (const std::exception& e) { return e.what(); } return ""; }
void interval_group(Results& R) {
  axial::Mesh g(45, 4, .009, 1);
  auto hv = thermal_view(patch(), 1500);
  auto a = owner<VA>(patch(), 1500, *hv);
  auto b = owner<VB>(patch(), 1500, *hv);
  const auto arr = arrays(g, 1e8, 0);
  // (1) default bound: 1e-14 accepted, 2e-14 refused with the v6.1 text, by v6.1 and v7 alike
  const auto ta = refusal_text([&] { VA::Field f(g, frame(*a, 1, 2e-14), arr.E, arr.B); });
  const auto tb = refusal_text([&] { VB::Field f(g, frame(*b, 1, 2e-14), arr.E, arr.B); });
  const auto ok14 = refusal_text([&] { VB::Field f(g, frame(*b, 1, 1e-14), arr.E, arr.B); });
  R.record("interval", "default trial interval bound 1e-14 s: 2e-14 refused by v6.1 and v7 with the same text, 1e-14 accepted", !ta.empty() && ta == tb && ok14.empty(),
           {{"v61", ta}, {"v7", tb}});
  // (2) configured bound 1e-13: a 28.4 fs step (2^-45 s, so that t + dt is exact in binary64: the sealed v3
  //     advance_source compares end - t with dt_s and refuses when one rounding makes it larger) through
  //     prepare_patch, birth, continuation, trial, commit
  s7::OwnerOptions opt;
  opt.interval.cap_s = 1e-13;
  auto c = std::make_unique<s7::Owner>(geom(), material(), std::vector<eh::Patch>{patch()}, std::vector<s7::PatchIdentity>{{"upper-foil-lower-surface", "conditional-single-layer", .005 * .025}},
                                       std::vector<s7::SourceSettings>{settings<VB>()}, 1500, "v3-source-scenario", 5e-9, *hv, opt);
  const double dt = 0x1p-45;
  auto fc = std::make_shared<const s7::Field>(g, frame(*c, 1, dt), arr.E, arr.B, 1e-13);
  const VB::View vc{fc->frame(), g, arr.E, arr.B};
  J ev;
  ev["configuration_differs_from_default"] = c->configuration_sha256() != b->configuration_sha256();
  ev["default_owner_refuses_configured_field"] = !refusal_text([&] {
    auto e = b->make_epoch(dt, "trial-midpoint-field");
    (void)b->prepare_patch(e, 0, binding(*b), fc, origin, v3::Exposure{fc->frame().begin_s, dt, 1, "native-surface-quadrature"});
  }).empty();
  bool chain = false;
  std::string chain_error;
  try {
    auto e = c->make_epoch(dt, "trial-midpoint-field");
    auto p = c->prepare_patch(e, 0, binding(*c), fc, origin, v3::Exposure{fc->frame().begin_s, dt, 1, "native-surface-quadrature"});
    auto s = sample<VB>(p, *c, 1e9, false);
    auto birth = c->begin_source(p, s, event(1));
    auto cont = c->advance_birth(p, birth, event(2));
    const bool active = cont.state().native.status == se::LayerStatus::Active;
    auto t = active ? c->prepare(e, fc, vc, std::span<const VB::Candidate>(&cont, 1), {}, {{}}) : c->prepare(e, fc, vc, {}, {}, {{}});
    chain = c->commit(std::move(t), vc);
    ev["expected_count_28fs"] = p.expected_number();
  } catch (const std::exception& x) { chain_error = x.what(); }
  ev["chain_published"] = chain;
  ev["chain_error"] = chain_error;
  ev["cap_above_hard_limit_refused"] = !refusal_text([&] {
    s7::OwnerOptions o2; o2.interval.cap_s = 2e-12;
    s7::Owner x(geom(), material(), {patch()}, {{"upper-foil-lower-surface", "conditional-single-layer", .005 * .025}}, {settings<VB>()}, 1500, "v3-source-scenario", 5e-9, *hv, o2);
  }).empty();
  ev["field_above_configured_bound_refused"] = !refusal_text([&] { s7::Field f(g, frame(*c, 2, 2e-13), arr.E, arr.B, 1e-13); }).empty();
  R.record("interval", "configured bound 1e-13 s: 28.4 fs step published; bound enters the configuration; mismatched field, 2e-13 step and cap > 1e-12 refused",
           ev["configuration_differs_from_default"] && ev["default_owner_refuses_configured_field"] && chain && ev["cap_above_hard_limit_refused"] && ev["field_above_configured_bound_refused"], ev);
  // (3) error-target policy: evidence required; dt above the bound refused (typed), below published with its receipt
  s7::OwnerOptions op;
  op.interval.cap_s = 1e-13;
  op.interval.emission_relative_error = 1e-3;
  op.interval.count_floor = 0;
  auto d = std::make_unique<s7::Owner>(geom(), material(), std::vector<eh::Patch>{patch()}, std::vector<s7::PatchIdentity>{{"upper-foil-lower-surface", "conditional-single-layer", .005 * .025}},
                                       std::vector<s7::SourceSettings>{settings<VB>()}, 1500, "v3-source-scenario", 5e-9, *hv, op);
  J ev3;
  auto fd = std::make_shared<const s7::Field>(g, frame(*d, 1, 2e-14), arr.E, arr.B, 1e-13);
  const VB::View vd{fd->frame(), g, arr.E, arr.B};
  auto ed = d->make_epoch(2e-14, "trial-midpoint-field");
  ev3["no_evidence_refused"] = !refusal_text([&] { (void)d->prepare(ed, fd, vd, {}, {}, {{}}, {}, false, nullptr, nullptr); }).empty();
  const std::vector<V3> origins{origin};
  auto now = d->rate_sample(ed, fd, origins);
  auto prev = now;
  prev.t_s = now.t_s - 1e-14;
  for (auto& x : prev.log_emission) x -= 1e-14 * 0.9e11;   // d ln N/dt = 0.9e11 s^-1 -> bound 2 * 1e-3 * 0.5 / 0.9e11 = 1.11e-14 s
  const s7::IntervalEvidence evid{prev, now};
  const auto bound = d->interval_bound(evid);
  std::string typed;
  try { (void)d->prepare(ed, fd, vd, {}, {}, {{}}, {}, false, nullptr, &evid); } catch (const s7::Refusal& r) { typed = r.what(); } catch (const std::exception& r) { typed = std::string("untyped: ") + r.what(); }
  ev3["bound_s"] = bound.bound_s;
  ev3["limiting"] = bound.limiting;
  ev3["above_bound_typed_refusal"] = typed;
  auto fe = std::make_shared<const s7::Field>(g, frame(*d, 3, 1e-14), arr.E, arr.B, 1e-13);
  const VB::View ve{fe->frame(), g, arr.E, arr.B};
  auto ee = d->make_epoch(1e-14, "trial-midpoint-field");
  bool pub = false;
  std::string rec;
  try {
    auto now2 = d->rate_sample(ee, fe, origins);
    auto prev2 = now2; prev2.t_s -= 1e-14;
    for (auto& x : prev2.log_emission) x -= 1e-14 * 0.9e11;
    const s7::IntervalEvidence ev2{prev2, now2};
    auto t = d->prepare(ee, fe, ve, {}, {}, {{}}, {}, false, nullptr, &ev2);
    rec = t.interval_record();
    pub = d->commit(std::move(t), ve);
  } catch (const std::exception& x) { rec = x.what(); }
  ev3["within_bound_published"] = pub;
  ev3["receipt"] = rec;
  R.record("interval", "error-target policy: evidence required; dt = 2e-14 above the 1.11e-14 bound refused as a typed Refusal (kind interval-bound); dt = 1e-14 published with the bound receipt",
           ev3["no_evidence_refused"] && typed.find("source-owner interval bound") != std::string::npos && typed.find("kind=interval-bound") != std::string::npos && typed.rfind("untyped", 0) != 0 && std::abs(bound.bound_s / (1e-14 / 0.9) - 1) < 1e-6 && pub &&
               rec.find("bound=") == 0,
           ev3);
}

// ---------------------------------------------------------------- v7: surface state, classes, phase map
eh::Patch water_ss304(double molecules, double activation_eV) {
  auto p = patch(); p.desorption_enabled = true; p.state_id = "H2O"; p.parent_id = "H2O"; p.initial_particles = molecules; p.initial_low_particles = 0;
  p.mass_kg = 18.01528 * 1.66053906660e-27; p.attempt_s = 1e13; p.activation_eV = activation_eV; p.desorption_eV = 0.87; p.internal_eV = 0; p.monolayer_m2 = 1e19;
  return p;
}
std::unique_ptr<s7::Owner> model_owner(const eh::Patch& p, double T, const eh::Coupler& view, sv::Registry& reg, double w0 = 0, bool wall = true) {
  auto st = settings<VB>(se::Config{4.3}, w0);
  st.surface = sv::ss304_surface(reg, true, wall);
  return std::make_unique<s7::Owner>(geom(), material(), std::vector<eh::Patch>{p}, std::vector<s7::PatchIdentity>{{"upper-foil-lower-surface", "conditional-single-layer", .005 * .025}},
                                     std::vector<s7::SourceSettings>{st}, T, "v3-source-scenario", 5e-9, view);
}
void surface_group(Results& R) {
  axial::Mesh g(45, 4, .009, 1);
  sv::Registry reg;
  const double N = 0.33 * .005 * .025 * 1e19;   // 0.33 monolayer
  auto wp = water_ss304(N, 0.87);
  auto hv = thermal_view(wp, 300);
  auto o = model_owner(wp, 300, *hv, reg);
  auto bench = owner<VB>(wp, 300, *hv);
  const double theta = o->surface_coverage(0);
  const double phi_expected = sv::work_function_eV(sv::ss304_work_function(reg), N / (.005 * .025 * 1e19), 300);
  J ev;
  ev["theta"] = theta;
  ev["phi_eV"] = o->surface_work_function_eV(0);
  ev["phi_expected_eV"] = phi_expected;
  ev["configuration_differs_from_benchmark"] = o->configuration_sha256() != bench->configuration_sha256();
  ev["v7_identity"] = o->v7_identity();
  ev["benchmark_owner_v7_identity"] = bench->v7_identity();
  // class exposures: expected count = J(T, beta F, phi) a_k A dt / qe, bitwise from the same expressions
  const double F = 4e7;
  auto pf = field<VB>(g, *o, 1, F);
  auto e = o->make_epoch(1e-14, "trial-midpoint-field");
  J classes = J::array();
  bool exact = true;
  for (std::size_t k = 0; k < o->enhancement_classes(0); ++k) {
    J row;
    row["class"] = k;
    try {
      auto p = o->prepare_class(e, 0, k, binding(*o), pf.f, origin, v3::Exposure{e.t, 1e-14, 1, "native-surface-quadrature"});
      const auto cl = sv::ss304_enhancement(reg).classes[k];   // by value (the set is a temporary)
      const double Fm = p.macroscopic_field_V_m();   // the owner's one-sided gather of the field arrays
      const auto d = electrode_heat_v6_1::emission::evaluate(300, cl.beta == 1 ? Fm : Fm * cl.beta, se::Config{o->surface_work_function_eV(0)});
      const double area = hv->area(0) * cl.area_fraction;
      const double n = d.current_A_m2 * area * 1e-14 / se::qe;
      exact = exact && std::bit_cast<std::uint64_t>(n) == std::bit_cast<std::uint64_t>(p.expected_number()) && p.beta() == cl.beta && p.work_function_eV() == o->surface_work_function_eV(0);
      row["beta"] = cl.beta; row["area_fraction"] = cl.area_fraction; row["expected_number"] = p.expected_number(); row["log_expected"] = p.log_expected_number();
    } catch (const std::exception& x) { row["refusal"] = x.what(); exact = exact && k == 3; }   // beta 500 at 4e7: 2e10 V/m outside the kernel domain
    classes.push_back(row);
  }
  ev["classes"] = classes;
  // phase map regimes
  const auto pm = o->phase_map(0, 4e7, 300);
  const auto ph = o->phase_map(0, 1e5, 2500);
  J map = J::array();
  for (const auto& x : pm) map.push_back({{"beta", x.beta}, {"F_local", x.local_field_V_m}, {"regime", x.regime}, {"log_J", x.log_current_A_m2}, {"kT_over_d", x.kT_over_d}, {"refusal", x.refusal.substr(0, 80)}});
  for (const auto& x : ph) map.push_back({{"beta", x.beta}, {"F_local", x.local_field_V_m}, {"T", 2500}, {"regime", x.regime}, {"above_barrier", x.above_barrier_fraction}});
  ev["phase_map"] = map;
  const bool regimes = pm.size() == 4 && pm[2].regime == "field-emission" && pm[3].regime == "outside-kernel-domain" && !pm[3].refusal.empty() && ph[0].regime == "thermionic-schottky";
  // image layer refused for an enhanced class
  auto p2 = o->prepare_class(e, 0, 2, binding(*o), pf.f, origin, v3::Exposure{e.t, 1e-14, 1, "native-surface-quadrature"});
  ev["enhanced_class_begin_source_refused"] = !refusal_text([&] { (void)o->begin_source(p2, sample<VB>(p2, *o, 1e9, false), event(1)); }).empty();
  R.record("surface", "work function from the exact stock (theta, phi of the cited set); class exposures bit-equal to J(T, beta F, phi) a_k A dt/qe; phase-map regimes; enhanced class refused on the image layer",
           std::abs(theta - 0.33) < 1e-12 && ev["phi_eV"] == phi_expected && ev["configuration_differs_from_benchmark"] && !bench->v7_identity() && exact && regimes &&
               ev["enhanced_class_begin_source_refused"],
           ev);
}

// ---------------------------------------------------------------- v7: births at the macroscopic matching plane
void birth_group(Results& R) {
  axial::Mesh g(45, 4, .009, 1);
  sv::Registry reg;
  auto wp = water_ss304(0.33 * .005 * .025 * 1e19, 0.87);
  auto hv = thermal_view(wp, 300);
  auto o = model_owner(wp, 300, *hv, reg);
  J ev;
  bool ok = true;
  s7::BirthPolicy pol;
  pol.subcell_extent_m = 20e-6;   // 0.1 of a 200 um cell
  pol.return_extent_m = 20e-6;
  // extracting field, flat and enhanced classes
  for (std::size_t k : {std::size_t(0), std::size_t(2)}) {
    auto pf = field<VB>(g, *o, 1 + k, 3e9 / (k ? 100 : 1));
    auto e = o->make_epoch(1e-14, "trial-midpoint-field");
    auto p = o->prepare_class(e, 0, k, binding(*o), pf.f, origin, v3::Exposure{e.t, 1e-14, 1, "native-surface-quadrature"});
    const auto b = o->sample_birth(p, 1e6, .6, .35, .1, 7, pol);
    const double Fm = p.macroscopic_field_V_m(), L = b.packet.matching_distance_m, phi = p.work_function_eV();
    const double K = se::kinetic_energy_J(b.packet.particle.proper_velocity_m_s) / se::qe;
    const double expect = b.packet.particle.normal_metal_energy_eV + b.packet.particle.transverse_energy_eV - phi + Fm * L;
    const auto& l = b.packet.particle.ledger;
    const double scale = std::abs(l.wall_heat_J) + std::abs(l.outgoing_K_J) + std::abs(l.electron_reference_change_J) + std::abs(l.resolved_matching_field_work_J);
    const bool row_ok = b.outcome == s7::BirthOutcome::Handoff && L > 0 && std::abs(K / expect - 1) < 1e-9 && std::abs(l.residual_J()) <= 1e-12 * scale &&
                        std::abs(b.displacement_m[1] + L) <= 1e-18 && L >= (phi - b.packet.particle.normal_metal_energy_eV) / Fm * (1 - 1e-12);
    ok = ok && row_ok;
    ev[std::string("class_") + std::to_string(k)] = {{"macroscopic_F", Fm}, {"beta", p.beta()}, {"matching_distance_m", L}, {"K_eV", K}, {"expected_K_eV", expect},
                                                     {"ledger_residual_J", l.residual_J()}, {"field_work_J", l.resolved_matching_field_work_J}, {"ok", row_ok}};
  }
  // retarding field: sub-cell return analysis
  {
    auto pf = field<VB>(g, *o, 9, -1e6);
    auto hot = thermal_view(wp, 1500);
    auto oh = model_owner(wp, 1500, *hot, reg);
    auto pf2 = field<VB>(g, *oh, 9, -1e6);
    auto e = oh->make_epoch(1e-14, "trial-midpoint-field");
    auto p = oh->prepare_class(e, 0, 0, binding(*oh), pf2.f, origin, v3::Exposure{e.t, 1e-14, 1, "native-surface-quadrature"});
    const auto r = oh->sample_birth(p, 1e6, .6, .35, .1, 8, pol);
    s7::BirthPolicy tight = pol;
    tight.return_extent_m = 1e-12;
    const auto h = oh->sample_birth(p, 1e6, .6, .35, .1, 8, tight);
    const auto& c = r.cycle;
    const double scale = std::abs(r.packet.particle.ledger.wall_heat_J) + std::abs(r.packet.particle.ledger.outgoing_K_J);
    const bool rok = r.outcome == s7::BirthOutcome::ReturnedSubcell && r.turning_distance_m <= 20e-6 && std::abs(c.wall_heat_J) <= 1e-12 * scale &&
                     c.physical_electron_number == 0 && c.particle_charge_C == 0 && h.outcome == s7::BirthOutcome::Handoff;
    ok = ok && rok;
    ev["retarding"] = {{"turning_distance_m", r.turning_distance_m}, {"return_time_B0_s", r.return_time_s}, {"cycle_net_wall_heat_J", c.wall_heat_J},
                       {"cycle_electrons", c.physical_electron_number}, {"tight_extent_outcome_handoff", h.outcome == s7::BirthOutcome::Handoff}, {"ok", rok}};
  }
  R.record("births", "matching-plane births: displacement L n with K = E + Et - phi + F L (image tail restored), exact ledger closure; retarding sub-cell returns close with zero net heat and charge",
           ok, ev);
}

// ---------------------------------------------------------------- v7: wall responses through the owner
std::vector<s7::WallHit> wall_hits(int n, double w) {
  std::vector<s7::WallHit> h;
  v7test::SplitMix gsm{9};
  for (int i = 0; i < n; ++i) {
    s7::WallHit x;
    x.carrier_id = 1000 + std::uint64_t(i) * 7919;
    const int kind = i % 4;
    x.species = kind == 0 ? SW7_E : kind == 1 ? SW7_E : kind == 2 ? SW7_H2OP : SW7_H2O;
    x.patch = 0;
    x.weight = w;
    const double E = kind == 0 ? 300 : kind == 1 ? 30 + 3000 * gsm.u() : kind == 2 ? 500 : 0.04;
    const double m = x.species == SW7_E ? se::me : 18.01473 * 1.66053906660e-27;
    const double K = E * se::qe, c = 299792458.0, us = std::sqrt(K * (K + 2 * m * c * c)) / (m * c);
    const double ct = 0.2 + 0.8 * gsm.u();
    x.proper_velocity_m_s = {us * std::sqrt(1 - ct * ct), us * ct, 0};   // into the LowY face (outward normal -y)
    x.hit_point_m = {-.002 + .004 * gsm.u(), .0025, 0};
    x.hit_time_s = 5e-9 + 1e-15 * i;
    h.push_back(x);
  }
  return h;
}
void wall_group(Results& R, const std::filesystem::path& out) {
  axial::Mesh g(45, 4, .009, 1);
  sv::Registry reg;
  const double area = .005 * .025, N = 0.8 * area * 1e19;   // 0.8 monolayer (below saturation of g(theta))
  auto wp = water_ss304(N, 0.87);
  auto hv = thermal_view(wp, 300);
  const double w0 = 1e12;
  auto o = model_owner(wp, 300, *hv, reg, w0);
  auto pf = field<VB>(g, *o, 1, 1e6);
  const VB::View v = view<VB>(pf, g);
  auto e = o->make_epoch(1e-14, "trial-midpoint-field");
  const auto hits = wall_hits(1000, 1e12);
  J ev;
  auto wr = o->respond(e, hits, 4242, 77);
  // order and partition independence of the response
  auto rev = hits;
  std::reverse(rev.begin(), rev.end());
  auto wr2 = o->respond(e, rev, 4242, 77);
  ev["content_independent_of_input_order"] = wr.content_sha256() == wr2.content_sha256();
  // exact ledger and births
  const auto tot = wr.total();
  std::int64_t outs = 0;
  for (const auto& x : wr.outputs()) outs += x.n_out;
  ev["births"] = wr.births().size();
  ev["births_equal_outputs"] = std::int64_t(wr.births().size()) == outs;
  ev["ledger_closes"] = tot.e_in_q == tot.e_out_q + tot.wall_q && tot.events == 1000;
  ev["stuck_molecules"] = wr.stuck(0).rounded_for_rate();
  ev["stimulated_hazard"] = wr.hazard(0);
  ev["loads"] = wr.loads().size();
  double loadJ = 0;
  for (const auto& t : wr.loads()) loadJ += t.authored_heat_J;
  double expectJ = 0;
  for (std::size_t i = 0; i < wr.hits().size(); ++i) expectJ += wr.hits()[i].weight * (double(wr.outputs()[i].wall_q - wr.outputs()[i].k_in_q) / SW7_QUANTA_PER_EV) * se::qe;
  ev["wall_load_J"] = loadJ;
  ev["wall_load_matches_quanta"] = std::abs(loadJ - expectJ) <= 1e-12 * std::abs(expectJ);
  // stimulated release through the quantized route: n = floor(D / w0) leaves of exactly w0 (the arrival owner draws SR(D/w0))
  const auto qth = o->release_quotes(e)[0];
  const auto q = o->release_quotes(e, wr)[0];
  const double D = q.released.rounded_for_rate();
  const std::size_t nleaf = std::size_t(D / w0);
  std::vector<std::vector<VB::Leaf>> leaves(1);
  for (std::size_t k = 0; k < nleaf; ++k) leaves[0].push_back({w0, eh::Draw{.3, .4, .5, .6, .7, (k + .5) / double(nleaf)}});
  ev["thermal_release"] = qth.released.rounded_for_rate();
  ev["release_with_stimulated"] = D;
  ev["leaves"] = nleaf;
  const double phi_before = o->surface_work_function_eV(0), theta_before = o->surface_coverage(0);
  const auto stock_before = o->inventory(0);
  bool committed = false;
  std::string err;
  try {
    auto t = o->prepare(e, pf.f, v, {}, {}, leaves, {}, false, &wr, nullptr);
    ev["trial_carries_wall"] = t.wall() != nullptr;
    committed = o->commit(std::move(t), v);
  } catch (const std::exception& x) { err = x.what(); }
  // exact stock: before - n w0 + stuck
  Count expect = stock_before;
  for (std::size_t k = 0; k < nleaf; ++k) expect = expect.debit(w0);
  for (double w : wr.stuck(0).positive_carriers()) expect = expect.credited(w);
  ev["committed"] = committed;
  ev["error"] = err;
  ev["stock_exact"] = o->inventory(0) == expect;
  ev["phi_before"] = phi_before; ev["phi_after"] = o->surface_work_function_eV(0);
  ev["theta_before"] = theta_before; ev["theta_after"] = o->surface_coverage(0);
  // staleness: the response of the previous state is refused after publication
  ev["stale_response_refused"] = !refusal_text([&] { (void)o->release_quotes(o->make_epoch(1e-14, "trial-midpoint-field"), wr); }).empty();
  // a benchmark patch refuses wall hits (v6.1 behaviour: first hits stay with the Hall sink)
  auto bo = owner<VB>(wp, 300, *hv);
  ev["benchmark_refuses_wall_hits"] = !refusal_text([&] { (void)bo->respond(bo->make_epoch(1e-14, "trial-midpoint-field"), hits, 1, 1); }).empty();
  const bool ok = ev["content_independent_of_input_order"] && ev["births_equal_outputs"] && ev["ledger_closes"] && ev["wall_load_matches_quanta"] && wr.hazard(0) > 0 &&
                  D > qth.released.rounded_for_rate() && nleaf > 0 && committed && ev["stock_exact"] && ev["phi_after"].get<double>() > phi_before &&
                  ev["stale_response_refused"] && ev["benchmark_refuses_wall_hits"];
  R.record("wall", "owner wall responses: order-independent content, births = outputs, exact ledger, loads = quanta; ESD/ISD release through n w0 leaves; exact stock; phi follows the coverage; stale and benchmark refusals",
           ok, ev);
  // checkpoint of the model-mode owner
  o->save((out / "checkpoint-v7-model.bin").string());
  auto hv2 = thermal_view(wp, 300);
  auto o2 = model_owner(wp, 300, *hv2, reg, w0);
  bool loaded = false;
  try { o2->load((out / "checkpoint-v7-model.bin").string()); loaded = true; } catch (const std::exception&) {}
  R.record("wall", "model-mode component checkpoint restores the exact stock (thermal clock bound)", loaded && o2->state_sha256() == o->state_sha256() && o2->surface_work_function_eV(0) == o->surface_work_function_eV(0),
           {{"loaded", loaded}});
}

// ---------------------------------------------------------------- v7: phase-owner temperatures, 3000-5000 K, liquid
struct FakePhase : s7::PhaseView {
  std::uint64_t rev = 1;
  double t = 5e-9;
  std::string tag = "a";
  std::vector<s7::PatchPhase> p;
  std::uint64_t revision() const override { return rev; }
  double time_s() const override { return t; }
  std::string sha256() const override { return "phase-view-" + tag + "-" + std::to_string(rev); }
  std::size_t patches() const override { return p.size(); }
  s7::PatchPhase patch(std::size_t j) const override { return p.at(j); }
};
std::unique_ptr<s7::Owner> phase_owner(const eh::Patch& pt, const eh::Coupler& view, sv::Registry& reg, const s7::OwnerOptions& o, double w0 = 0) {
  auto st = settings<VB>(se::Config{4.3}, w0);
  st.surface = sv::ss304_surface(reg, true, true);
  return std::make_unique<s7::Owner>(geom(), material(), std::vector<eh::Patch>{pt}, std::vector<s7::PatchIdentity>{{"upper-foil-lower-surface", "conditional-single-layer", .005 * .025}},
                                     std::vector<s7::SourceSettings>{st}, 300, "v3-source-scenario", 5e-9, view, o);
}
void phase_group(Results& R) {
  axial::Mesh g(45, 4, .009, 1);
  sv::Registry reg;
  auto hv = thermal_view(patch(), 300);
  FakePhase fp;
  fp.p.resize(1);
  fp.p[0].mean_T = 4000; fp.p[0].max_T = 4100; fp.p[0].liquid_fraction = 1; fp.p[0].melt_depth_m = 5e-6; fp.p[0].evaporated_atoms = {1e10, 2e9, 1e9};
  s7::OwnerOptions o;
  o.emission_T_max_K = 5000; o.water_writeoff_T_K = 2900; o.phase = &fp;
  auto own = phase_owner(patch(), *hv, reg, o);
  J ev;
  auto pf = field<VB>(g, *own, 1, 1e8);
  auto e = own->make_epoch(1e-14, "trial-midpoint-field");
  const std::size_t K = own->enhancement_classes(0);
  ev["classes_with_liquid"] = K;
  ev["solid_class_zero_area_refused"] = !refusal_text([&] { (void)own->prepare_class(e, 0, 0, binding(*own), pf.f, origin, v3::Exposure{e.t, 1e-14, 1, "native-surface-quadrature"}); }).empty();
  bool liq_ok = false;
  try {
    auto p = own->prepare_class(e, 0, K - 1, binding(*own), pf.f, origin, v3::Exposure{e.t, 1e-14, 1, "native-surface-quadrature"});
    const double phl = std::nearbyint(sv::ss304_work_function(reg).phi_liquid_eV / 0x1p-12) * 0x1p-12;
    const auto d = electrode_heat_v7::emission_hot::evaluate(4000, p.macroscopic_field_V_m(), se::Config{phl}, 5000);
    const double n = d.current_A_m2 * hv->area(0) * 1e-14 / se::qe;
    liq_ok = p.temperature_K() == 4000 && p.work_function_eV() == phl && std::bit_cast<std::uint64_t>(n) == std::bit_cast<std::uint64_t>(p.expected_number());
    ev["liquid_expected_number"] = p.expected_number(); ev["liquid_phi"] = p.work_function_eV(); ev["T_used"] = p.temperature_K();
  } catch (const std::exception& x) { ev["liquid_error"] = x.what(); }
  const auto pm = own->phase_map(0, 1e8, 4000);
  ev["liquid_regime"] = pm.back().regime;
  // default T_max 3000 K: the sealed refusal above 3000 K
  s7::OwnerOptions o3;
  o3.phase = &fp;
  auto own3 = phase_owner(patch(), *hv, reg, o3);
  auto pf3 = field<VB>(g, *own3, 2, 1e8);
  const auto r3 = refusal_text([&] { (void)own3->prepare_class(e, 0, own3->enhancement_classes(0) - 1, binding(*own3), pf3.f, origin, v3::Exposure{e.t, 1e-14, 1, "native-surface-quadrature"}); });
  ev["default_cap_refusal"] = r3;
  R.record("phase", "phase-owner temperature drives emission; molten patch emits through the liquid class (beta 1, liquid phi) with the extended kernel at 4000 K; solid classes zero-area; default 3000 K cap keeps the sealed refusal",
           ev["classes_with_liquid"] == 5 && ev["solid_class_zero_area_refused"] && liq_ok && pm.back().regime == "thermionic-schottky" &&
               r3.find("1..3000K") != std::string::npos,
           ev);
  // water above the 3000 K law cap: typed refusal; residue write-off at T >= 2900 K; stale phase; metal vapor record
  J ew;
  {
    const double w0 = 1e12;
    FakePhase hot;
    hot.p.resize(1);
    hot.p[0].mean_T = 3500;
    s7::OwnerOptions oh;
    oh.emission_T_max_K = 5000; oh.phase = &hot;
    auto wp = water_ss304(2e11, 0.87);
    auto wv = thermal_view(wp, 300);
    auto ow = phase_owner(wp, *wv, reg, oh, w0);
    std::string typed;
    try { (void)ow->release_quotes(ow->make_epoch(1e-14, "trial-midpoint-field")); } catch (const s7::Refusal& r) { typed = r.what(); } catch (const std::exception& r) { typed = std::string("untyped ") + r.what(); }
    ew["water_above_3000K"] = typed;
    FakePhase warm;
    warm.p.resize(1);
    warm.p[0].mean_T = 2950; warm.p[0].liquid_fraction = 1; warm.p[0].evaporated_atoms = {3e8, 7e7, 4e7};
    s7::OwnerOptions ow2;
    ow2.emission_T_max_K = 5000; ow2.water_writeoff_T_K = 2900; ow2.phase = &warm;
    auto oz = phase_owner(wp, *wv, reg, ow2, w0);
    auto fz = field<VB>(g, *oz, 3, 1e6);
    const VB::View vz = view<VB>(fz, g);
    auto ez = oz->make_epoch(1e-14, "trial-midpoint-field");
    const auto q = oz->release_quotes(ez)[0];
    const double D = q.released.rounded_for_rate();
    std::vector<std::vector<VB::Leaf>> leaves(1);
    for (std::size_t k = 0; k < std::size_t(D / w0); ++k) leaves[0].push_back({w0, eh::Draw{.3, .4, .5, .6, .7, .5}});
    try {
      auto tr = oz->prepare(ez, fz.f, vz, {}, {}, leaves, {}, false, nullptr, nullptr);
      bool load = false, evap = false;
      for (const auto& x : tr.loads()) { load = load || x.channel_id == "surface-v7-water-writeoff"; evap = evap || x.channel_id.find("evapor") != std::string::npos; }
      ew["writeoff_load"] = load; ew["evaporation_load"] = evap;
      ew["writeoff_molecules"] = tr.water_writeoff().empty() ? 0. : tr.water_writeoff()[0].rounded_for_rate();
      ew["metal_vapor_atoms_reported"] = tr.phase_reports().empty() ? 0. : tr.phase_reports()[0].evaporated_atoms[0];
      warm.tag = "b";
      ew["stale_phase_preflight"] = oz->preflight(tr, vz);
      warm.tag = "a";
      ew["committed"] = oz->commit(std::move(tr), vz);
      ew["stock_after"] = oz->inventory(0).rounded_for_rate();
    } catch (const std::exception& x) { ew["error"] = x.what(); }
    ew["release_D"] = D;
  }
  R.record("phase", "water stock above 3000 K is a typed surface-domain refusal; sub-w0 residue written off explicitly at T >= 2900 K with its binding-energy load; stale phase report refused; metal vapor recorded, never a load",
           ew["water_above_3000K"].get<std::string>().find("surface-domain") != std::string::npos && ew.value("writeoff_load", false) && !ew.value("evaporation_load", true) &&
               ew.value("writeoff_molecules", 0.) > 0 && ew.value("metal_vapor_atoms_reported", 0.) == 3e8 && !ew.value("stale_phase_preflight", true) && ew.value("committed", false) &&
               ew.value("stock_after", 1.) == 0,
           ew);
}
}  // namespace

int main(int argc, char** argv) {
  Results R;
  int rc = 0;
  try {
    if (argc != 2) throw std::runtime_error("usage: test_source_v7 OUTDIR");
    const std::filesystem::path out = argv[1];
    std::filesystem::create_directories(out);
    a02_group(R);
    identity_group(R, out);
    interval_group(R);
    surface_group(R);
    birth_group(R);
    wall_group(R, out);
    phase_group(R);
    std::ofstream(out / "source-v7.json") << J{{"schema", "source-v7-qualification-v1"}, {"cases", R.cases}, {"failed", R.failed}}.dump(1) << '\n';
    rc = R.failed == 0 ? 0 : 1;
  } catch (const std::exception& e) {
    std::cerr << "fatal: " << e.what() << '\n';
    rc = 2;
  }
  return rc;
}
