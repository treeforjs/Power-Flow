// Sealed heat/source owner v5 harness for the native v6 source/arrival sweep.
// Fixture recipes follow the arrival owner v1 test (test_source_v5.cpp), which
// follows the source v5 focused tests: 4x2 thermal mesh, Hall axial mesh
// (45, 4, .009, 1), upper-electrode lower face, dt = 5e-15 s here.
// Compiled with the source v5 include paths only (no native headers).
#include "hall-owner/development/hall-finite-foil-em-20261029/axial_surface_path_port.hpp"
#include "source_owner.hpp"
#include "ss304_heat_material_v2.hpp"
#include "pf_arrival/source_v5_adapter.hpp"
#include "source_v5_harness.hpp"
#include <optional>

using namespace electrode_heat_v3;
namespace s5 = electrode_heat_v5;
namespace pa = pf_arrival;
namespace ps = pf_arrival::source_v5;

namespace source_harness {
namespace {
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
                                  ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,
                                  8000, eh::material_headers_sha256);
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
  return {"native-v6-source-sweep", "e-", "electrode-electron-reservoir", "native-layer",
          "canonical-electron-" + std::to_string(native_id), seq, seq};
}
const V3 origin{.0001, .0025, .0002};
std::unique_ptr<s5::Owner> source_owner(eh::Patch p, double T, const std::string &sampling,
                                        double start = 5e-9) {
  return std::make_unique<s5::Owner>(
      geom(), material(), std::vector<eh::Patch>{p},
      std::vector<s5::PatchIdentity>{{"upper-foil-lower-surface", "conditional-single-layer", .005 * .025}},
      std::vector<s5::SourceSettings>{{se::Config{4.3}, layer_policy(), sampling}}, T,
      "v3-source-scenario", start);
}
TrialFrame frame4(const s5::Owner &o, std::uint64_t stage, double dt) {
  TrialFrame f;
  f.role = FieldRole::TotalMidpoint;
  f.accepted_revision = o.thermal_state().revision;
  f.trial_stage = stage;
  f.mesh_revision = 1;
  f.material_revision = 2;
  f.begin_s = o.thermal_state().time;
  f.dt_s = dt;
  f.end_s = f.begin_s + f.dt_s;
  f.field_stage_s = f.begin_s + .5 * f.dt_s;
  f.arena_hash = o.state_sha256();
  f.field_hash = "declared-parent-field";
  f.mesh_hash = "declared-parent-mesh";
  f.material_hash = "global-material-v2-map";
  f.species_hash = "actual-source-v1-electron";
  return f;
}
PatchBinding bind4(const s5::Owner &o) {
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
const axial::Mesh &mesh() {
  static const axial::Mesh g(45, 4, .009, 1);
  return g;
}
struct Arrays {
  foil::Vec E, B;
};
Arrays arrays(double F) {
  Arrays a{foil::Vec(mesh().ne), foil::Vec(mesh().nb)};
  for (int e = 0; e < mesh().ne; ++e)
    if (mesh().edges[e].axis == 1)
      a.E[e] = F;
  return a;
}
std::shared_ptr<const s5::Field> field4(const s5::Owner &o, std::uint64_t stage, double F, double dt) {
  auto a = arrays(F);
  return std::make_shared<const s5::Field>(mesh(), frame4(o, stage, dt), a.E, a.B);
}
} // namespace

struct Neutral::Impl {
  std::unique_ptr<s5::Owner> owner;
  std::optional<eh::Epoch> epoch;
  std::optional<s5::ReleaseQuote> quote;
  std::optional<s5::Trial> trial;
  std::shared_ptr<const s5::Field> field;
  Arrays a = arrays(0);
  std::optional<TrialFrame> frame;
  double dt = 0;
  std::uint64_t stage = 0;
};
Neutral::Neutral(double N0, double rate, const std::string &sampling, double start)
    : impl_(std::make_unique<Impl>()) {
  impl_->owner = source_owner(water(N0, rate), 1500, sampling, start);
}
Neutral::~Neutral() = default;
std::uint64_t Neutral::revision() const { return impl_->owner->thermal_state().revision; }
std::string Neutral::state_sha256() const { return impl_->owner->state_sha256(); }
pa::NeutralRelease Neutral::quote(const pa::Policy &P, std::uint32_t stream, double dt) {
  auto &i = *impl_;
  i.dt = dt;
  i.epoch.emplace(i.owner->make_epoch(dt, "actual-parent-midpoint"));
  auto q = i.owner->release_quotes(*i.epoch);
  i.quote.emplace(q.at(0));
  return ps::neutral_input(P, stream, *i.quote, *i.owner, *i.epoch);
}
Step Neutral::prepare(pa::Staged &staged, std::uint32_t stream) {
  auto &i = *impl_;
  Step out;
  out.source_before = i.owner->state_sha256();
  out.debit = i.quote->released;
  out.inventory_before = i.owner->inventory(0);
  ++i.stage;
  i.field = field4(*i.owner, i.stage, 0, i.dt);
  i.frame.emplace(i.field->frame());
  mfem::ConstantCoefficient zero(0);
  s5::FinalFieldView v{*i.frame, mesh(), i.a.E, i.a.B};
  auto leaves = ps::neutral_leaves(staged, 1);
  try {
    i.trial.emplace(i.owner->prepare(*i.epoch, zero, i.field, v, {}, {}, leaves));
  } catch (const s5::Refusal &r) {
    out.outcome = "source-refused";
    out.refusal = std::string(s5::kind_name(r.kind)) + "/" + r.quantity;
    return out;
  } catch (const std::exception &e) {
    out.outcome = "source-refused-untyped";
    out.refusal = e.what();
    return out;
  }
  out.inventory_after = i.trial->inventory(0);
  staged.bind_source_times(stream, ps::neutral_times(*i.trial, staged, stream, "lower-face"));
  for (const auto &p : i.trial->packets()) {
    Birth b;
    b.weight = p.weight;
    b.time = p.time;
    b.x = p.x;
    b.u = p.proper_velocity_m_s;
    b.outward_normal = {0, -1, 0};
    b.event = p.event_id;
    out.births.push_back(b);
  }
  out.source_after = i.trial->content_sha256();
  out.admitted = true;
  out.outcome = "admitted";
  return out;
}
bool Neutral::commit() {
  auto &i = *impl_;
  s5::FinalFieldView v{*i.frame, mesh(), i.a.E, i.a.B};
  const bool ok = i.owner->commit(std::move(*i.trial), v);
  i.trial.reset();
  return ok;
}

struct Electron::Impl {
  std::unique_ptr<s5::Owner> owner;
  std::optional<s5::PreparedPatch> patch;
  std::shared_ptr<const s5::Field> field;
  double F = 0;
};
Electron::Electron(const std::string &sampling, double F, double start)
    : impl_(std::make_unique<Impl>()) {
  impl_->owner = source_owner(metal_patch(), 1500, sampling, start);
  impl_->F = F;
}
Electron::~Electron() = default;
std::uint64_t Electron::revision() const { return impl_->owner->thermal_state().revision; }
std::string Electron::state_sha256() const { return impl_->owner->state_sha256(); }
bool Electron::expose(std::uint32_t stream, double area_fraction, double dt, pa::ElectronExposure &out,
                      std::string &refusal) {
  auto &i = *impl_;
  i.field = field4(*i.owner, 1, i.F, dt);
  const electrode_heat_v3::Exposure x{i.field->frame().begin_s, dt, area_fraction, "physical-quadrature"};
  try {
    i.patch.emplace(i.owner->prepare_patch(i.owner->make_epoch(dt, "actual-parent-midpoint"), 0, bind4(*i.owner),
                                           i.field, origin, x));
  } catch (const s5::Refusal &r) {
    refusal = std::string(s5::kind_name(r.kind)) + "/" + r.quantity;
    return false;
  } catch (const std::exception &e) {
    refusal = std::string("untyped/") + e.what();
    return false;
  }
  out = ps::electron_input(stream, *i.patch, x, i.field->frame().end_s);
  return true;
}
Step Electron::births(const pa::Staged &staged) {
  auto &i = *impl_;
  Step out;
  out.source_before = i.owner->state_sha256();
  out.source_after = out.source_before; // candidates only; no source commit here
  pa::Count total;
  std::uint64_t seq = 1;
  for (const auto &b : staged.births()) {
    const auto sample = ps::birth_sample(staged, b, *i.patch, i.owner->thermal_state().revision);
    try {
      auto c = i.owner->begin_source(*i.patch, sample, event(seq++, b.native_id));
      const auto &layer = c.state().native;
      Birth x;
      x.weight = layer.particle.weight;
      x.time = layer.current_time_s;
      for (int d = 0; d < 3; ++d)
        x.x[d] = c.state().surface_origin_m[d] + layer.position_from_surface_m[d];
      x.u = layer.particle.proper_velocity_m_s;
      x.outward_normal = layer.outward_normal;
      x.event = c.content_sha256();
      x.active_layer = layer.status == se::LayerStatus::Active;
      total = total.credited(x.weight);
      out.births.push_back(x);
    } catch (const s5::Refusal &r) {
      out.outcome = "source-refused";
      out.refusal = std::string(s5::kind_name(r.kind)) + "/" + r.quantity;
      return out;
    } catch (const std::exception &e) {
      out.outcome = "source-refused-untyped";
      out.refusal = e.what();
      return out;
    }
  }
  out.debit = total;
  out.admitted = true;
  out.outcome = "admitted";
  return out;
}
double probe_neutral_debit(double N0, double rate, double dt) {
  auto o = source_owner(water(N0, rate), 1500, "probe");
  try {
    auto q = o->release_quotes(o->make_epoch(dt, "actual-parent-midpoint"));
    return q[0].released.zero() ? 0 : q[0].released.rounded_for_rate();
  } catch (const std::exception &) {
    return -1;
  }
}
double probe_electron_count(double F, double area, double dt) {
  auto o = source_owner(metal_patch(), 1500, "probe");
  auto f = field4(*o, 1, F, dt);
  try {
    return o->prepare_patch(o->make_epoch(dt, "actual-parent-midpoint"), 0, bind4(*o), f, origin,
                            {f->frame().begin_s, dt, area, "physical-quadrature"})
        .expected_number();
  } catch (const std::exception &) {
    return -1;
  }
}
} // namespace source_harness
