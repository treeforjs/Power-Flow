#pragma once
// Coordinator v6: identities derived from live owner objects (species catalog,
// native/Hall translation frame, configuration), the physical setup of the
// water-breakdown pilot (geometry, foil-slice source patches, arrival streams)
// and small shared helpers. Owners: native candidate006 (births, sinks, binary
// flows motion), Hall corrective v3 (sinks, banks, ThermalV3, newborns), heat/
// source owner v6 (one thermal owner, quantized neutral release), arrival owner
// v1.1, parent-coverage v2 chemistry. Sealed sources are included unchanged.
#include "raikou/mfem_pic_full_time.h"
#include "root_api.hpp"
#include "hall-finite-foil-em-20261029/axial_surface_path_port.hpp"
#include "ss304_heat_material_v2.hpp"
#include "water_scenario.hpp"
#include "pic-electrode-heat-source-v6-20261030/native/mfem_electrode_heat_source_v6/source_owner.hpp"
#include "hall_source_binding.hpp"
#include "pf_arrival_v11/arrival_owner.hpp"
#include "pf_arrival_v11/source_v6_adapter.hpp"
#include "pf_arrival_v11/native_v6_births.hpp"
#include "parent_coverage_v2_subflow.hpp"
#include "common_commit_v6.hpp"
#include "hall_v3_flows_adapter.hpp"
#include <bit>
#include <chrono>
#include <map>
#include <optional>
#include <span>

namespace raikou::native_water::v6 {
namespace h = ::hall_driver;
namespace hc = ::hall_corrective_v3;
namespace cv = ::raikou::native_water::coverage_v2;
namespace cc = ::raikou::native_water::common_commit;
namespace c6 = ::raikou::native_water::common_commit_v6;
namespace s6 = ::electrode_heat_v6;
namespace v3k = ::electrode_heat_v3;
namespace eh = ::electrode_heat_v2;
namespace se = ::surface_emission;
namespace hb = ::hall_source_binding;
namespace pa = ::pf_arrival_v11;
namespace pas = ::pf_arrival_v11::source_v6;
namespace av = ::raikou::mfem_backend::arrival_v11;
namespace hn = ::hall_v3_native_v6;

inline void need(bool ok, const std::string& what) {
  if (!ok) throw std::runtime_error("coordinator v6: " + what);
}

// Length-prefixed little-endian binary encoding (the v4/v5 coordinator convention).
struct Encoding {
  std::string bytes;
  void u64(std::uint64_t v) { for (int i = 0; i < 8; ++i) bytes.push_back(char((v >> (8 * i)) & 255)); }
  void i64(std::int64_t v) { u64(static_cast<std::uint64_t>(v)); }
  void f64(double v) { u64(std::bit_cast<std::uint64_t>(v)); }
  void text(const std::string& s) { u64(s.size()); bytes.append(s); }
  void count(const m::ArenaCount& c) { const auto l = c.positive_carriers(); u64(l.size()); for (double w : l) f64(w); }
  std::string finish() const { return sha(bytes); }
};

// ---------------------------------------------------------------- species and frame
inline std::map<std::string, const m::ArenaStateSpec*> state_index(const std::vector<m::ArenaStateSpec>& states) {
  std::map<std::string, const m::ArenaStateSpec*> out;
  for (const auto& s : states) need(!s.id.empty() && out.emplace(s.id, &s).second, "species catalog requires unique nonempty identities");
  need(!out.empty() && out.size() <= 256, "bounded nonempty species catalog");
  return out;
}
inline bool same_channel(const m::ArenaChannel& a, const m::ArenaChannel& b) {
  return a.id == b.id && a.primary_state == b.primary_state && a.target_state == b.target_state && a.product_states == b.product_states &&
         a.dataset_sha256 == b.dataset_sha256 && a.kernel_manifest_sha256 == b.kernel_manifest_sha256 &&
         a.kernel_binding_file_sha256 == b.kernel_binding_file_sha256 && a.model == b.model;
}
inline std::string species_sha256(const std::vector<m::ArenaStateSpec>& states, const std::vector<m::ArenaChannel>& channels) {
  const auto index = state_index(states);
  std::map<std::string, const m::ArenaChannel*> catalog;
  for (const auto& c : channels) {
    need(!c.id.empty() && catalog.emplace(c.id, &c).second, "channel catalog requires unique nonempty identities");
    need(index.count(c.primary_state) && index.count(c.target_state), "channel reactant species undeclared");
    for (const auto& p : c.product_states) need(index.count(p) > 0, "channel product species undeclared");
  }
  Encoding e;
  e.text("coordinator-v6-species-catalog-v1");
  e.u64(index.size());
  for (const auto& [id, s] : index) {
    e.text(id); e.i64(s->hydrogen); e.i64(s->oxygen); e.i64(s->charge_e);
    e.f64(s->mass_kg); e.f64(s->internal_J); e.f64(s->reference_J); e.f64(s->charge_e * m::pic_e);
  }
  e.u64(catalog.size());
  for (const auto& [id, c] : catalog) {
    e.text(id); e.text(c->primary_state); e.text(c->target_state); e.u64(c->product_states.size());
    for (const auto& p : c->product_states) e.text(p);
    e.text(c->dataset_sha256); e.text(c->kernel_manifest_sha256); e.text(c->kernel_binding_file_sha256); e.text(c->model);
  }
  return e.finish();
}
inline void require_census_species(std::span<const m::ArenaRateEntry> census, const std::vector<m::ArenaStateSpec>& states) {
  const auto index = state_index(states);
  for (const auto& row : census) {
    const auto it = index.find(row.state);
    need(it != index.end(), "census carrier has an undeclared species");
    need(row.native.charge_C == it->second->charge_e * m::pic_e && row.native.mass_kg == it->second->mass_kg,
         "census carrier charge/mass differs from its declared species");
  }
}
inline void require_arena_species(const m::PicArena& arena, const std::vector<m::ArenaStateSpec>& states,
                                  const std::vector<m::ArenaChannel>& channels, std::size_t cap) {
  const auto index = state_index(states);
  for (const auto& c : channels) need(same_channel(arena.channel_binding(c.id), c), "registered native channel differs from the declared catalog");
  for (const auto& [id, p] : arena.snapshot()->packets) { (void)id; need(index.count(p.state) > 0, "accepted native packet has an undeclared species"); }
  require_census_species(arena.rate_population(cap), states);
}
// Exactly the translation native admission compares (native -> Hall field).
inline m::PicVector native_to_field_shift(const m::GapPicOptions& o) { return {-.5 * o.width_m, -.5 * o.gap_m, 0.}; }
inline std::string frame_sha256(const m::GapPicOptions& o, const hc::Model& model) {
  Encoding e;
  e.text("native-to-hall-translation-frame-v4");
  for (double x : native_to_field_shift(o)) e.f64(x);
  e.text(model.motion_mesh_sha256(o));
  e.text(h::kinetic_policy_sha256(h::declared_kinetic_policy(o)));
  return e.finish();
}
// Hall corrective field box translated by (+w/2, +g/2, 0): z in (-depth/2, +depth/2).
// Foil planes: the binary64 values whose Hall field-frame cell coordinate
// ((x - w/2) + w/2)/dx is an exact node plane or rounds to the metal side, so a
// conductor-bank deposit at a foil face (face birth, collected sink) puts no
// rounding weight on a vacuum node. On the default 200 um mesh this moves the
// lower-foil outer face (side/2 - 2.7 mm) by 2 ulp to 1.8e-3 and the left x faces
// (side/2 - 2.5 mm) by 1 ulp to 2.0e-3; every other plane is unchanged. Any other
// mesh keeps the plain values.
inline m::GapPicOptions hall_geometry(const axial::Mesh& g) {
  m::GapPicOptions o;
  o.nx = o.ny = g.nx; o.width_m = o.gap_m = g.side; o.depth_m = g.depth;
  o.finite_foils = o.finite_axial_extent = true; o.short_z_m = -g.depth / 2; o.feed_z_m = g.depth / 2;
  double left = g.side / 2 - .0025, lower_bottom = g.side / 2 - .0027;
  if (g.nx == 45 && g.side == .009) { left = .002; lower_bottom = .0018; }
  o.lower_foil = {left, g.side / 2 + .0025, lower_bottom, g.side / 2 - .0025};
  o.upper_foil = {left, g.side / 2 + .0025, g.side / 2 + .0025, g.side / 2 + .0027};
  return o;
}

// ---------------------------------------------------------------- physical setup
// Species of the parent channel (V3 test_native_parent_subflow.cpp states).
inline constexpr double water_mass_kg = 2.9915e-26;
inline std::vector<m::ArenaStateSpec> pilot_states() {
  return {{"electron", 0, 0, -1, w::me, 0, 0}, {"H2O", 2, 1, 0, water_mass_kg, 0, 0},
          {"H2O+", 2, 1, 1, water_mass_kg - w::me + w::threshold_eV * pc::qe / (pc::c * pc::c), 0, w::threshold_eV * pc::qe}};
}
inline std::vector<m::ArenaChannel> pilot_channels(const cv::ArenaBinding& a) {
  return {{w::parent_channel, "electron", "H2O", {"electron", "electron", "H2O+"}, a.dataset_sha256, a.kernel_manifest_sha256,
           a.kernel_binding_file_sha256, a.model}};
}
inline const char* face_name(eh::Face f) {
  switch (f) { case eh::Face::LowY: return "low-y"; case eh::Face::HighY: return "high-y"; case eh::Face::LowX: return "low-x"; default: return "high-x"; }
}
// Run configuration (bound into the coordinator configuration identity).
struct PilotSettings {
  double work_function_eV = 4.3;     // clean-surface value of the source kernel (A15 limitation)
  double w0_electron = 1e6, w0_neutral = 1e6;
  double initial_coverage = .33;     // oxide-water submonolayer (water_scenario.hpp recipe)
  std::uint64_t arrival_seed = 20261031;
  int faces_per_slice = 4;           // 2: inner/outer faces only; 4: also the two x edges
  std::size_t maximum_carriers = 100000, maximum_native_records = 4096;
  std::uint64_t native_memory_budget = 2ull << 30, hall_memory_budget = 1ull << 30;
  std::string scenario = "power-flow-pilot-v6/oxide-water-scenario";
  std::string json_text() const {
    std::ostringstream s;
    s << std::hexfloat << "phi=" << work_function_eV << ";w0e=" << w0_electron << ";w0n=" << w0_neutral << ";theta0=" << initial_coverage
      << ";seed=" << arrival_seed << ";faces=" << faces_per_slice << ";carriers=" << maximum_carriers << ";records=" << maximum_native_records
      << ";native_budget=" << native_memory_budget << ";hall_budget=" << hall_memory_budget << ";scenario=" << scenario
      << ";emission=face-birth/kernel-phase-space-at-surface-step-or-SN-matching-plane";
    return s.str();
  }
};
// One source patch (global index g): Hall slice (body, k), local patch index j, face.
struct PatchRef {
  std::size_t slice = 0, local = 0;
  int body = 0, k = 0;
  eh::Face face = eh::Face::LowY;
  std::string id;
  std::uint32_t electron_stream = 0, neutral_stream = 0;
};
struct SourceSetup {
  std::vector<hb::SliceSource> slices;
  std::vector<PatchRef> patches;
  std::shared_ptr<const pa::Policy> arrival;
};
inline constexpr std::uint32_t electron_stream_base = 1000, neutral_stream_base = 2000;
inline SourceSetup make_source_setup(const axial::Mesh& g, const PilotSettings& S) {
  need(S.faces_per_slice == 2 || S.faces_per_slice == 4, "faces per slice must be 2 or 4");
  SourceSetup out;
  pa::Configuration ac;
  ac.run_seed = S.arrival_seed;
  const eh::Face faces[4] = {eh::Face::HighY, eh::Face::LowY, eh::Face::LowX, eh::Face::HighX};
  std::uint32_t gidx = 0;
  // Streams first (the sampling policy identity must be known to the source settings).
  for (int body = 0; body < 2; ++body)
    for (int k = g.kshort; k < g.kfeed; ++k)
      for (int f = 0; f < S.faces_per_slice; ++f) {
        ac.streams.push_back({electron_stream_base + gidx, pa::Kind::electron, gidx, "electron", S.w0_electron, 4096});
        ac.streams.push_back({neutral_stream_base + gidx, pa::Kind::neutral, gidx, "H2O", S.w0_neutral, 4096});
        ++gidx;
      }
  out.arrival = std::make_shared<const pa::Policy>(ac);
  gidx = 0;
  for (int body = 0; body < 2; ++body)
    for (int k = g.kshort; k < g.kfeed; ++k) {
      hb::SliceSource s;
      s.body = body; s.slice = k; s.source_hash = S.scenario;
      const auto gg = hb::slice_geometry(g, body, k);
      for (int f = 0; f < S.faces_per_slice; ++f) {
        PatchRef r;
        r.slice = out.slices.size(); r.local = std::size_t(f); r.body = body; r.k = k; r.face = faces[f];
        r.id = "b" + std::to_string(body) + "-k" + std::to_string(k) + "-" + face_name(faces[f]);
        r.electron_stream = electron_stream_base + gidx; r.neutral_stream = neutral_stream_base + gidx;
        auto p = eh::oxide_water_patch(gg, r.id, faces[f], 0, 1, S.initial_coverage, water_mass_kg, "H2O", "H2O", S.scenario,
                                       eh::WaterRange::ConditionalAnalytic298To3000);
        s.patches.push_back(p);
        s.settings.push_back({se::Config{S.work_function_eV}, se::LayerPolicy{}, out.arrival->sampling_policy_identity(), S.w0_neutral});
        s.settings.back().layer.local_linear_vacuum_extent_m = 25e-6;
        out.patches.push_back(r);
        ++gidx;
      }
      // Physical areas exactly as the slice Coupler computes them.
      eh::Coupler probe(gg, hb::slice_material(), s.patches, 300., s.source_hash, 0.);
      for (std::size_t j = 0; j < s.patches.size(); ++j)
        s.identities.push_back({"foil-" + std::to_string(body) + "/slice-" + std::to_string(k) + "/" + face_name(s.patches[j].face),
                                "oxide-water-conditional-temkin-eyring", probe.area(j)});
      out.slices.push_back(std::move(s));
    }
  return out;
}
// Global origin (field frame, global z) of a patch: the face centre.
inline v3k::V3 patch_origin(const axial::Mesh& g, const PatchRef& r) {
  const auto gg = hb::slice_geometry(g, r.body, r.k);
  const double zc = hb::slice_center_z(g, r.k);
  switch (r.face) {
    case eh::Face::LowY: return {gg.x0 + .5 * gg.width, gg.y0, zc};
    case eh::Face::HighY: return {gg.x0 + .5 * gg.width, gg.y0 + gg.thickness, zc};
    case eh::Face::LowX: return {gg.x0, gg.y0 + .5 * gg.thickness, zc};
    default: return {gg.x0 + gg.width, gg.y0 + .5 * gg.thickness, zc};
  }
}
// Native point on a foil face from a field-frame point with local z: tangential
// coordinates are translated, the face-normal coordinate is the native foil plane.
inline m::PicVector native_face_point(const m::GapPicOptions& o, const PatchRef& r, const std::array<double, 3>& x_field, double z_center) {
  m::PicVector p{x_field[0] + .5 * o.width_m, x_field[1] + .5 * o.gap_m, z_center + x_field[2]};
  const auto& f = r.body ? o.upper_foil : o.lower_foil;
  int axis = 1;
  double plane = 0;
  switch (r.face) {
    case eh::Face::LowY: axis = 1; plane = f.bottom_m; break;
    case eh::Face::HighY: axis = 1; plane = f.top_m; break;
    case eh::Face::LowX: axis = 0; plane = f.left_m; break;
    default: axis = 0; plane = f.right_m; break;
  }
  need(std::abs(p[std::size_t(axis)] - plane) <= 1e-12, "source face point is not on the native foil face");
  p[std::size_t(axis)] = plane;
  return p;
}
inline v3k::PatchBinding patch_binding(const eh::Coupler& heat, const axial::Mesh& g, const PatchRef& r) {
  v3k::PatchBinding b;
  b.thermal_revision = heat.state().revision;
  b.source_accepted_state = b.thermal_revision + 1;
  b.patch_id = 1 + r.local + 16 * (r.slice + 1);
  b.thermal_mesh_identity = hb::slice_geometry(g, r.body, r.k).revision;
  b.thermal_material_hash = eh::material_headers_sha256;
  b.material_phase_admitted = b.planar_barrier_admitted = b.electron_supply_admitted = true;
  b.temperature = v3k::HeldTemperature::AcceptedPatchMean;
  b.extrusion_center_z_m = hb::slice_center_z(g, r.k);
  b.extrusion_height_m = g.dz;
  return b;
}
// Electron position uniforms: words 4 and 5 of the electron marks stream of the
// birth (arrival RNG layout v1, purpose 2, blocks 0..1 are the four marks).
inline std::array<double, 2> electron_position_uniforms(const pa::Policy& P, const pa::Birth& b) {
  const auto key = P.key(b.step, b.stream_id, b.index, PFA_PURPOSE_ELECTRON_MARKS);
  const auto x = pfa::stream_block(key, 2u);
  return {pfa::u01_open(pfa::low_word(x)), pfa::u01_open(pfa::high_word(x))};
}
// Field-frame point (local z) on a patch for open uniforms (along, z), the same
// placement law as the source v6 neutral packets.
inline std::array<double, 3> face_point_local(const eh::Geometry& gg, const eh::Patch& p, double along_u, double z_u) {
  const double along = p.begin + (p.end - p.begin) * along_u;
  if (p.face == eh::Face::LowY || p.face == eh::Face::HighY)
    return {gg.x0 + along * gg.width, p.face == eh::Face::LowY ? gg.y0 : gg.y0 + gg.thickness, (z_u - .5) * gg.height};
  return {p.face == eh::Face::LowX ? gg.x0 : gg.x0 + gg.width, gg.y0 + along * gg.thickness, (z_u - .5) * gg.height};
}
inline double seconds_since(std::chrono::steady_clock::time_point a) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
}
inline double exact_value(const m::ArenaCount& c) { return c.zero() ? 0. : c.rounded_for_rate(); }
inline m::ArenaCount exact_add(m::ArenaCount a, const m::ArenaCount& b) {
  for (double w : b.positive_carriers()) a = a.credited(w);
  return a;
}
}  // namespace raikou::native_water::v6
