// Native v5 full-time canonical codecs, typed journal rows and bounded history
// storage. Included early in the arena TU (before the reaction builder).
// Identities are SHA-256 over canonical binary encodings (mfem_pic_canonical.h);
// JSON is produced only as human-readable export (audits, checkpoint export).
namespace full_time_detail {
namespace cn = canonical;
using cn::Reader;
using cn::Writer;
template <class IO> constexpr bool reading = std::is_same_v<IO, Reader>;
FullTimeResourceRefusal refusal(FullTimeResourceRefusal::Kind k,
                                const char *message) {
  return FullTimeResourceRefusal(k, message);
}
using RK = FullTimeResourceRefusal::Kind;
std::uint64_t checked_add(std::uint64_t a, std::uint64_t b, const char *what) {
  require(b <= std::numeric_limits<std::uint64_t>::max() - a, what);
  return a + b;
}
ArenaJson options_json(const GapPicOptions &o) {
  return {{"nx", o.nx},
          {"ny", o.ny},
          {"width_m", o.width_m},
          {"gap_m", o.gap_m},
          {"depth_m", o.depth_m},
          {"lower_potential_V", o.lower_potential_V},
          {"upper_potential_V", o.upper_potential_V},
          {"magnetic_field_T", o.magnetic_field_T},
          {"distributed_fields", o.distributed_fields},
          {"finite_foils", o.finite_foils},
          {"finite_axial_extent", o.finite_axial_extent},
          {"short_z_m", o.short_z_m},
          {"feed_z_m", o.feed_z_m},
          {"lower_foil",
           {o.lower_foil.left_m, o.lower_foil.right_m, o.lower_foil.bottom_m,
            o.lower_foil.top_m}},
          {"upper_foil",
           {o.upper_foil.left_m, o.upper_foil.right_m, o.upper_foil.bottom_m,
            o.upper_foil.top_m}}};
}
void options_guard(const GapPicOptions &o) {
  require(o.nx > 0 && o.nx <= 4096 && o.ny > 0 && o.ny <= 4096 &&
              std::size_t(o.nx) * std::size_t(o.ny) <= 262144,
          "bounded full-time immutable native mesh");
  require(!o.distributed_fields && std::isfinite(o.width_m) && o.width_m > 0 &&
              std::isfinite(o.gap_m) && o.gap_m > 0 &&
              std::isfinite(o.depth_m) && o.depth_m > 0 &&
              std::isfinite(o.lower_potential_V) &&
              std::isfinite(o.upper_potential_V) &&
              std::isfinite(o.short_z_m) && std::isfinite(o.feed_z_m),
          "finite full-time native options");
  for (double v : o.magnetic_field_T)
    require(std::isfinite(v), "finite prescribed-field option");
  for (const auto &f : {o.lower_foil, o.upper_foil})
    for (double v : {f.left_m, f.right_m, f.bottom_m, f.top_m})
      require(std::isfinite(v), "finite foil option");
  if (o.finite_foils) {
    for (const auto &f : {o.lower_foil, o.upper_foil}) {
      require(0 < f.left_m && f.left_m < f.right_m && f.right_m < o.width_m &&
                  0 < f.bottom_m && f.bottom_m < f.top_m && f.top_m < o.gap_m,
              "full-time native foil bounds");
      for (double x : {f.left_m, f.right_m})
        require(std::abs(x / o.width_m * o.nx -
                         std::round(x / o.width_m * o.nx)) <= 1e-9,
                "native foil x alignment");
      for (double y : {f.bottom_m, f.top_m})
        require(std::abs(y / o.gap_m * o.ny - std::round(y / o.gap_m * o.ny)) <=
                    1e-9,
                "native foil y alignment");
    }
    require(o.lower_foil.top_m < o.upper_foil.bottom_m, "native electrode gap");
  }
  if (o.finite_axial_extent)
    require(o.finite_foils && o.feed_z_m > o.short_z_m &&
                std::abs(o.feed_z_m - o.short_z_m - o.depth_m) <=
                    1e-12 * o.depth_m,
            "native axial geometry");
}
GapPicOptions read_options(const ArenaJson &j) {
  require(j.is_object() && j.size() == 15, "typed native option schema");
  GapPicOptions o;
  o.nx = static_cast<int>(json_size(j.at("nx"), 4096));
  o.ny = static_cast<int>(json_size(j.at("ny"), 4096));
  o.width_m = json_scalar(j.at("width_m"));
  o.gap_m = json_scalar(j.at("gap_m"));
  o.depth_m = json_scalar(j.at("depth_m"));
  o.lower_potential_V = json_scalar(j.at("lower_potential_V"));
  o.upper_potential_V = json_scalar(j.at("upper_potential_V"));
  o.magnetic_field_T = json_vector(j.at("magnetic_field_T"));
  o.distributed_fields = json_boolean(j.at("distributed_fields"));
  o.finite_foils = json_boolean(j.at("finite_foils"));
  o.finite_axial_extent = json_boolean(j.at("finite_axial_extent"));
  o.short_z_m = json_scalar(j.at("short_z_m"));
  o.feed_z_m = json_scalar(j.at("feed_z_m"));
  auto foil = [](const ArenaJson &a) {
    require(a.is_array() && a.size() == 4, "four-coordinate foil required");
    return GapFoil{json_scalar(a[0]), json_scalar(a[1]), json_scalar(a[2]),
                   json_scalar(a[3])};
  };
  o.lower_foil = foil(j.at("lower_foil"));
  o.upper_foil = foil(j.at("upper_foil"));
  options_guard(o);
  return o;
}
// ---------------------------------------------------------------------------
// Canonical codecs. One template serves both directions so every decoder is the
// exact inverse of its encoder.
template <class IO> void codec_digest(IO &io, std::string &s) {
  io.text(s, 64);
  if constexpr (reading<IO>)
    require(cn::is_digest(s), "canonical digest field");
}
template <class IO> void codec_options(IO &io, GapPicOptions &o) {
  io.tag("native-gap-options-v5");
  std::int64_t nx = o.nx, ny = o.ny;
  io.i64(nx);
  io.i64(ny);
  if constexpr (reading<IO>) {
    require(nx > 0 && nx <= 4096 && ny > 0 && ny <= 4096, "options codec");
    o.nx = int(nx);
    o.ny = int(ny);
  }
  for (double *x : {&o.width_m, &o.gap_m, &o.depth_m, &o.lower_potential_V,
                    &o.upper_potential_V})
    io.f64(*x);
  io.vec(o.magnetic_field_T);
  io.boolean(o.distributed_fields);
  io.boolean(o.finite_foils);
  io.boolean(o.finite_axial_extent);
  io.f64(o.short_z_m);
  io.f64(o.feed_z_m);
  for (auto *f : {&o.lower_foil, &o.upper_foil})
    for (double *x : {&f->left_m, &f->right_m, &f->bottom_m, &f->top_m})
      io.f64(*x);
}
template <class IO> void codec_count(IO &io, ArenaCount &c) {
  // Unit exponent -1074; canonical greedy positive-carrier decomposition.
  if constexpr (reading<IO>) {
    const auto n = io.count(ArenaCount::maximum_carriers);
    std::vector<double> w(n);
    for (auto &x : w)
      io.f64(x);
    c = ArenaCount::from_positive_carriers(w);
    require(w == c.positive_carriers(), "canonical count decomposition");
  } else {
    const auto w = c.positive_carriers();
    io.u64(w.size());
    for (double x : w)
      io.f64(x);
  }
}
template <class IO> void codec_common(IO &io, FullTimeStepBinding &b) {
  io.u64(b.accepted_step);
  io.u64(b.accepted_arena_revision);
  io.f64(b.begin_s);
  io.f64(b.end_s);
  io.f64(b.dt_s);
  for (auto *s : {&b.accepted_arena_sha256, &b.accepted_field_sha256,
                  &b.mesh_sha256, &b.material_sha256, &b.common_clock_sha256})
    io.text(*s, 1024);
}
template <class IO> void codec_exposure(IO &io, ExposureBinding &e) {
  codec_common(io, e.common);
  std::uint64_t slot = e.slot == FullTimeSlot::after_transport ? 1 : 0;
  io.u64(slot);
  if constexpr (reading<IO>) {
    require(slot <= 1, "canonical exposure slot");
    e.slot = slot ? FullTimeSlot::after_transport
                  : FullTimeSlot::before_transport;
  }
  io.f64(e.frozen_physical_stage_s);
  io.f64(e.collision_duration_s);
  io.u64(e.owner_stage_sequence);
  io.text(e.frozen_field_sha256, 1024);
}
template <class IO> void codec_completion(IO &io, ArenaStochasticCompletion &c) {
  for (auto *s : {&c.integration_owner_id, &c.before_rng_clock_sha256,
                  &c.after_rng_clock_sha256, &c.common_clock_sha256,
                  &c.channel_id, &c.dataset_sha256, &c.kernel_manifest_sha256,
                  &c.kernel_binding_file_sha256, &c.model})
    io.text(*s, 1024);
  io.u64(c.stage_sequence);
  io.f64(c.frozen_physical_stage_s);
  io.f64(c.common_interval_begin_s);
  io.f64(c.common_interval_end_s);
}
template <class IO>
void codec_boundary(IO &io, FullTimeAcceptedBoundary &b, bool arena) {
  io.u64(b.accepted_step);
  io.u64(b.accepted_arena_revision);
  io.f64(b.time_s);
  if (arena)
    io.text(b.accepted_arena_sha256, 1024);
  for (auto *s : {&b.field_sha256, &b.heat_sha256, &b.source_sha256,
                  &b.rng_sha256, &b.mesh_sha256, &b.material_sha256,
                  &b.frame_sha256, &b.common_clock_sha256})
    io.text(*s, 1024);
  io.boolean(b.authored_full_time_genesis);
}
template <class IO> void codec_genesis(IO &io, FullTimeGenesisBinding &g) {
  io.f64(g.time_s);
  for (auto *s : {&g.origin_provenance, &g.field_sha256, &g.mesh_sha256,
                  &g.material_sha256, &g.frame_sha256, &g.common_clock_sha256,
                  &g.heat_sha256, &g.source_sha256, &g.rng_sha256})
    io.text(*s, 1024);
}
template <class T, class F> std::string digest_of(const char *tag, T &value, F f) {
  Writer w;
  w.tag(tag);
  f(w, value);
  return cn::sha256_hex(w.bytes);
}
std::string completion_sha(const ArenaStochasticCompletion &c) {
  auto x = c;
  return digest_of("native-full-time-completion-v5", x,
                   [](Writer &w, auto &v) { codec_completion(w, v); });
}
std::string boundary_sha(const FullTimeAcceptedBoundary &b) {
  auto x = b;
  return digest_of("native-full-time-boundary-v5", x,
                   [](Writer &w, auto &v) { codec_boundary(w, v, true); });
}
std::string common_sha(const FullTimeStepBinding &b) {
  auto x = b;
  return digest_of("native-full-time-common-v5", x,
                   [](Writer &w, auto &v) { codec_common(w, v); });
}
bool same_common(const FullTimeStepBinding &a, const FullTimeStepBinding &b) {
  return a.accepted_step == b.accepted_step &&
         a.accepted_arena_revision == b.accepted_arena_revision &&
         a.begin_s == b.begin_s && a.end_s == b.end_s && a.dt_s == b.dt_s &&
         a.accepted_arena_sha256 == b.accepted_arena_sha256 &&
         a.accepted_field_sha256 == b.accepted_field_sha256 &&
         a.mesh_sha256 == b.mesh_sha256 &&
         a.material_sha256 == b.material_sha256 &&
         a.common_clock_sha256 == b.common_clock_sha256;
}
bool same_boundary_payload(const FullTimeAcceptedBoundary &a,
                           const FullTimeAcceptedBoundary &b) {
  return a.accepted_step == b.accepted_step &&
         a.accepted_arena_revision == b.accepted_arena_revision &&
         a.time_s == b.time_s && a.field_sha256 == b.field_sha256 &&
         a.heat_sha256 == b.heat_sha256 && a.source_sha256 == b.source_sha256 &&
         a.rng_sha256 == b.rng_sha256 && a.mesh_sha256 == b.mesh_sha256 &&
         a.material_sha256 == b.material_sha256 &&
         a.frame_sha256 == b.frame_sha256 &&
         a.common_clock_sha256 == b.common_clock_sha256 &&
         a.authored_full_time_genesis == b.authored_full_time_genesis;
}
// ---- JSON export/typed import of the same structures (human-readable) ----
ArenaJson common_json(const FullTimeStepBinding &b) {
  return {{"accepted_step", b.accepted_step},
          {"accepted_arena_revision", b.accepted_arena_revision},
          {"begin_s", b.begin_s},
          {"end_s", b.end_s},
          {"dt_s", b.dt_s},
          {"accepted_arena_sha256", b.accepted_arena_sha256},
          {"accepted_field_sha256", b.accepted_field_sha256},
          {"mesh_sha256", b.mesh_sha256},
          {"material_sha256", b.material_sha256},
          {"common_clock_sha256", b.common_clock_sha256}};
}
FullTimeStepBinding read_common(const ArenaJson &j) {
  require(j.is_object() && j.size() == 10, "full-time common schema");
  FullTimeStepBinding b;
  b.accepted_step = json_u64(j.at("accepted_step"));
  b.accepted_arena_revision = json_u64(j.at("accepted_arena_revision"));
  b.begin_s = json_scalar(j.at("begin_s"));
  b.end_s = json_scalar(j.at("end_s"));
  b.dt_s = json_scalar(j.at("dt_s"));
  b.accepted_arena_sha256 = json_text(j.at("accepted_arena_sha256"));
  b.accepted_field_sha256 = json_text(j.at("accepted_field_sha256"));
  b.mesh_sha256 = json_text(j.at("mesh_sha256"));
  b.material_sha256 = json_text(j.at("material_sha256"));
  b.common_clock_sha256 = json_text(j.at("common_clock_sha256"));
  return b;
}
ArenaJson exposure_json(const ExposureBinding &b) {
  return {{"common", common_json(b.common)},
          {"slot", b.slot == FullTimeSlot::before_transport
                       ? "before_transport"
                       : "after_transport"},
          {"frozen_physical_stage_s", b.frozen_physical_stage_s},
          {"collision_duration_s", b.collision_duration_s},
          {"owner_stage_sequence", b.owner_stage_sequence},
          {"frozen_field_sha256", b.frozen_field_sha256}};
}
ArenaJson boundary_json(const FullTimeAcceptedBoundary &b,
                        bool include_arena = true) {
  ArenaJson j = {{"accepted_step", b.accepted_step},
                 {"accepted_arena_revision", b.accepted_arena_revision},
                 {"time_s", b.time_s},
                 {"field_sha256", b.field_sha256},
                 {"heat_sha256", b.heat_sha256},
                 {"source_sha256", b.source_sha256},
                 {"rng_sha256", b.rng_sha256},
                 {"mesh_sha256", b.mesh_sha256},
                 {"material_sha256", b.material_sha256},
                 {"frame_sha256", b.frame_sha256},
                 {"common_clock_sha256", b.common_clock_sha256},
                 {"authored_full_time_genesis", b.authored_full_time_genesis}};
  if (include_arena)
    j["accepted_arena_sha256"] = b.accepted_arena_sha256;
  return j;
}
FullTimeAcceptedBoundary read_boundary(const ArenaJson &j) {
  require(j.is_object() && j.size() == 13, "full-time boundary schema");
  FullTimeAcceptedBoundary b;
  b.accepted_step = json_u64(j.at("accepted_step"));
  b.accepted_arena_revision = json_u64(j.at("accepted_arena_revision"));
  b.time_s = json_scalar(j.at("time_s"));
  b.accepted_arena_sha256 = json_text(j.at("accepted_arena_sha256"));
  b.field_sha256 = json_text(j.at("field_sha256"));
  b.heat_sha256 = json_text(j.at("heat_sha256"));
  b.source_sha256 = json_text(j.at("source_sha256"));
  b.rng_sha256 = json_text(j.at("rng_sha256"));
  b.mesh_sha256 = json_text(j.at("mesh_sha256"));
  b.material_sha256 = json_text(j.at("material_sha256"));
  b.frame_sha256 = json_text(j.at("frame_sha256"));
  b.common_clock_sha256 = json_text(j.at("common_clock_sha256"));
  b.authored_full_time_genesis =
      json_boolean(j.at("authored_full_time_genesis"));
  return b;
}
FullTimeAcceptedBoundary boundary_payload(const ArenaJson &j) {
  require(j.is_object() && j.size() == 12 &&
              !j.contains("accepted_arena_sha256"),
          "full-time boundary payload schema");
  auto full = j;
  full["accepted_arena_sha256"] = std::string(64, '0');
  auto b = read_boundary(full);
  b.accepted_arena_sha256.clear();
  return b;
}
ArenaJson genesis_json(const FullTimeGenesisBinding &g) {
  return {{"schema", "authored-full-time-genesis-v5"},
          {"time_s", g.time_s},
          {"origin_provenance", g.origin_provenance},
          {"field_sha256", g.field_sha256},
          {"mesh_sha256", g.mesh_sha256},
          {"material_sha256", g.material_sha256},
          {"frame_sha256", g.frame_sha256},
          {"common_clock_sha256", g.common_clock_sha256},
          {"heat_sha256", g.heat_sha256},
          {"source_sha256", g.source_sha256},
          {"rng_sha256", g.rng_sha256}};
}
FullTimeGenesisBinding read_full_genesis(const ArenaJson &j) {
  require(j.is_object() && j.size() == 11 &&
              j.at("schema") == "authored-full-time-genesis-v5",
          "full-time genesis schema");
  FullTimeGenesisBinding g;
  g.time_s = json_scalar(j.at("time_s"));
  g.origin_provenance = json_text(j.at("origin_provenance"));
  g.field_sha256 = json_text(j.at("field_sha256"));
  g.mesh_sha256 = json_text(j.at("mesh_sha256"));
  g.material_sha256 = json_text(j.at("material_sha256"));
  g.frame_sha256 = json_text(j.at("frame_sha256"));
  g.common_clock_sha256 = json_text(j.at("common_clock_sha256"));
  g.heat_sha256 = json_text(j.at("heat_sha256"));
  g.source_sha256 = json_text(j.at("source_sha256"));
  g.rng_sha256 = json_text(j.at("rng_sha256"));
  return g;
}
// ---------------------------------------------------------------------------
// Typed journal rows. Memory retains each row's canonical bytes and digest
// only; the motion row carries digests of the per-carrier endpoint and path
// records (block Merkle roots), never the records themselves.
enum class RowKind : std::uint64_t {
  outer_begin = 1,
  reaction = 2,
  exposure = 3,
  motion = 4,
  boundary = 5,
  birth = 6, // v6: birth credit (between first exposure and motion)
  sink = 7   // v6: first-hit removals (directly after a flows motion row)
};
struct ProductRow {
  std::string id, state;
  PicVector u{};
  bool primary = false;
};
struct ReactionRow {
  std::string event_id, channel, primary, target, untouched_primary_id;
  double event_weight = 0;
  ArenaCount primary_before, target_before, primary_remainder,
      target_remainder;
  PicVector primary_u{}, target_u{}, primary_position_m{}, target_position_m{};
  std::vector<ProductRow> products;
  double physical_stage_s = 0, previous_dt_s = 0, next_dt_s = 0,
         collision_elapsed_s = 0;
  std::string state_prefix, mesh_prefix, clock_prefix, phase, provenance_kind;
  std::uint64_t accepted_step = 0;
  std::string dataset_sha256, kernel_manifest_sha256,
      kernel_binding_file_sha256, conditional_model, selection_method,
      selection_state_prefix, reference_convention;
  double energy_defect_J = 0, energy_roundoff_bound_J = 0,
         external_particle_source_debit_J = 0;
  PicVector momentum_defect_kg_m_s{};
  bool mcc_heating_and_em_current_qualified = false;
};
template <class IO> void codec_reaction(IO &io, ReactionRow &r) {
  for (auto *s : {&r.event_id, &r.channel, &r.primary, &r.target,
                  &r.untouched_primary_id})
    io.text(*s, 1024);
  io.f64(r.event_weight);
  for (auto *c : {&r.primary_before, &r.target_before, &r.primary_remainder,
                  &r.target_remainder})
    codec_count(io, *c);
  for (auto *v : {&r.primary_u, &r.target_u, &r.primary_position_m,
                  &r.target_position_m})
    io.vec(*v);
  std::uint64_t n = r.products.size();
  if constexpr (reading<IO>) {
    n = io.count(8);
    r.products.resize(n);
  } else
    io.u64(n);
  for (auto &p : r.products) {
    io.text(p.id, 1024);
    io.text(p.state, 1024);
    io.vec(p.u);
    io.boolean(p.primary);
  }
  io.f64(r.physical_stage_s);
  io.f64(r.previous_dt_s);
  io.f64(r.next_dt_s);
  io.f64(r.collision_elapsed_s);
  for (auto *s : {&r.state_prefix, &r.mesh_prefix, &r.clock_prefix, &r.phase,
                  &r.provenance_kind})
    io.text(*s, 1024);
  io.u64(r.accepted_step);
  for (auto *s : {&r.dataset_sha256, &r.kernel_manifest_sha256,
                  &r.kernel_binding_file_sha256, &r.conditional_model,
                  &r.selection_method, &r.selection_state_prefix,
                  &r.reference_convention})
    io.text(*s, 1024);
  io.f64(r.energy_defect_J);
  io.f64(r.energy_roundoff_bound_J);
  io.f64(r.external_particle_source_debit_J);
  io.vec(r.momentum_defect_kg_m_s);
  io.boolean(r.mcc_heating_and_em_current_qualified);
}
struct MotionRow {
  std::string record_schema, kinetic_policy_schema;
  std::string transport_trial_sha256, census_sha256, endpoints_root,
      paths_root, species_sha256, trial_field_sha256, kinetic_policy_sha256,
      current_A_m_sha256, initial_charge_C_sha256, endpoint_charge_C_sha256,
      staged_arena_sha256;
  std::string common_clock_sha256, mesh_sha256, material_sha256,
      field_before_sha256, field_after_sha256, heat_after_sha256,
      source_after_sha256;
  std::uint64_t accepted_step = 0, staged_arena_revision = 0, trial_stage = 0,
                carriers = 0;
  double begin_s = 0, end_s = 0, dt_s = 0;
  std::array<double, 7> budget{};
  PicVector native_to_field_m{};
};
template <class IO> void codec_motion(IO &io, MotionRow &r) {
  io.text(r.record_schema, 256);
  io.text(r.kinetic_policy_schema, 256);
  for (auto *s :
       {&r.transport_trial_sha256, &r.census_sha256, &r.endpoints_root,
        &r.paths_root, &r.species_sha256, &r.trial_field_sha256,
        &r.kinetic_policy_sha256, &r.current_A_m_sha256,
        &r.initial_charge_C_sha256, &r.endpoint_charge_C_sha256,
        &r.staged_arena_sha256, &r.common_clock_sha256, &r.mesh_sha256,
        &r.material_sha256, &r.field_before_sha256, &r.field_after_sha256,
        &r.heat_after_sha256, &r.source_after_sha256})
    io.text(*s, 1024);
  io.u64(r.accepted_step);
  io.u64(r.staged_arena_revision);
  io.u64(r.trial_stage);
  io.u64(r.carriers);
  io.f64(r.begin_s);
  io.f64(r.end_s);
  io.f64(r.dt_s);
  for (auto &x : r.budget)
    io.f64(x);
  io.vec(r.native_to_field_m);
}
// v6 birth credit row: binds the arrival staged token, plan digest, B_s and
// n_s, the source transition and the exact per-stream credit = debit.
struct BirthStreamRow {
  std::uint64_t stream_id = 0, births = 0;
  std::string state;
  ArenaCount credit;
};
struct BirthRow {
  std::uint64_t accepted_step = 0, base = 0, count = 0, next_id_after = 0;
  double begin_s = 0, end_s = 0;
  std::string plan_sha256, staged_token_sha256, source_before_sha256,
      source_after_sha256, births_root;
  std::vector<BirthStreamRow> streams;
};
template <class IO> void codec_birth(IO &io, BirthRow &r) {
  io.u64(r.accepted_step);
  io.u64(r.base);
  io.u64(r.count);
  io.u64(r.next_id_after);
  io.f64(r.begin_s);
  io.f64(r.end_s);
  for (auto *s : {&r.plan_sha256, &r.staged_token_sha256,
                  &r.source_before_sha256, &r.source_after_sha256,
                  &r.births_root})
    codec_digest(io, *s);
  std::uint64_t n = r.streams.size();
  if constexpr (reading<IO>) {
    n = io.count(8192);
    r.streams.resize(n);
  } else
    io.u64(n);
  for (auto &x : r.streams) {
    io.u64(x.stream_id);
    io.u64(x.births);
    io.text(x.state, 1024);
    codec_count(io, x.credit);
  }
}
// v6 sink row: removals of one motion (same transport identity).
struct SinkRow {
  std::uint64_t accepted_step = 0, census = 0, survivors = 0, carriers = 0,
                packets = 0;
  std::string transport_trial_sha256, records_root, ledger_sha256,
      removed_packets_sha256;
};
template <class IO> void codec_sink(IO &io, SinkRow &r) {
  io.u64(r.accepted_step);
  io.u64(r.census);
  io.u64(r.survivors);
  io.u64(r.carriers);
  io.u64(r.packets);
  for (auto *s : {&r.transport_trial_sha256, &r.records_root, &r.ledger_sha256,
                  &r.removed_packets_sha256})
    codec_digest(io, *s);
}
struct FullRow {
  RowKind kind = RowKind::outer_begin;
  std::string previous_chain;
  FullTimeStepBinding common;
  std::string outer_sha;
  ReactionRow reaction;
  ExposureBinding exposure;
  ArenaStochasticCompletion completion;
  std::string stage_before_chain, stage_after_reaction_chain;
  MotionRow motion;
  std::string transport_trial_sha256, first_prefix, second_prefix;
  FullTimeAcceptedBoundary boundary;
  BirthRow birth;
  SinkRow sink;
};
template <class IO> void codec_row(IO &io, FullRow &r) {
  io.tag("native-full-time-journal-row-v5");
  std::uint64_t kind = static_cast<std::uint64_t>(r.kind);
  io.u64(kind);
  if constexpr (reading<IO>) {
    require(kind >= 1 && kind <= 7, "canonical row kind");
    r.kind = static_cast<RowKind>(kind);
  }
  codec_digest(io, r.previous_chain);
  switch (r.kind) {
  case RowKind::outer_begin:
    codec_common(io, r.common);
    io.text(r.outer_sha, 1024);
    break;
  case RowKind::reaction:
    codec_reaction(io, r.reaction);
    break;
  case RowKind::exposure:
    codec_exposure(io, r.exposure);
    codec_completion(io, r.completion);
    io.text(r.stage_before_chain, 1024);
    io.text(r.stage_after_reaction_chain, 1024);
    break;
  case RowKind::motion:
    codec_motion(io, r.motion);
    break;
  case RowKind::boundary:
    codec_common(io, r.common);
    for (auto *s : {&r.outer_sha, &r.transport_trial_sha256, &r.first_prefix,
                    &r.second_prefix})
      io.text(*s, 1024);
    codec_boundary(io, r.boundary, false);
    break;
  case RowKind::birth:
    codec_birth(io, r.birth);
    break;
  case RowKind::sink:
    codec_sink(io, r.sink);
    break;
  }
}
std::string encode_row(const FullRow &r) {
  Writer w;
  codec_row(w, const_cast<FullRow &>(r));
  return std::move(w.bytes);
}
FullRow decode_row(const std::string &bytes) {
  Reader in(bytes);
  FullRow r;
  codec_row(in, r);
  in.finish();
  return r;
}
const char *row_kind_name(RowKind k) {
  switch (k) {
  case RowKind::outer_begin:
    return "full_time_outer_begin";
  case RowKind::reaction:
    return "reaction";
  case RowKind::exposure:
    return "full_time_exposure";
  case RowKind::motion:
    return "full_time_motion";
  case RowKind::boundary:
    return "full_time_boundary";
  case RowKind::birth:
    return "full_time_birth";
  case RowKind::sink:
    return "full_time_sink";
  }
  return "unknown";
}
} // namespace full_time_detail
// ---------------------------------------------------------------------------
// Bounded history: immutable shared pieces. A snapshot copy shares all of it.
struct FullTimeStoredRow {
  std::string bytes, sha256;
  full_time_detail::RowKind kind{};
};
struct FullTimeWindowNode {
  FullTimeStoredRow row;
  std::shared_ptr<const FullTimeWindowNode> previous;
  std::uint64_t rows = 0, bytes = 0; // window totals up to and including row
  FullTimeWindowNode(FullTimeStoredRow r,
                     std::shared_ptr<const FullTimeWindowNode> p)
      : row(std::move(r)), previous(std::move(p)) {
    rows = (previous ? previous->rows : 0) + 1;
    bytes = (previous ? previous->bytes : 0) + row.bytes.size();
  }
  // Iterative release: a long window never recurses through destructors.
  ~FullTimeWindowNode() noexcept {
    auto p = std::move(const_cast<FullTimeWindowNode *>(this)->previous);
    while (p && p.use_count() == 1) {
      auto next = std::move(const_cast<FullTimeWindowNode &>(*p).previous);
      p = std::move(next);
    }
  }
};
struct FullTimeAnchorRecord {
  std::string bytes, sha256;
  FullTimeAnchorFields fields;
};
struct FullTimeAnchorChain {
  std::vector<std::shared_ptr<const FullTimeAnchorRecord>> anchors;
  std::uint64_t bytes = 0;
};
struct FullTimeHistoryBase {
  std::string chain;
  std::uint64_t revision = 0, bytes = 0, steps = 0;
  double previous_dt_s = 0;
  FullTimeAcceptedBoundary boundary; // includes accepted arena sha
  ArenaJson completions = ArenaJson::object();
  SharedSet<std::string> event_ids;
  ~FullTimeHistoryBase() noexcept { drain_json(completions); }
};
struct FullTimeHistoryPolicyState {
  FullTimeHistoryPolicy policy;
  std::string digest, root, genesis_chain;
  FullTimeAcceptedBoundary genesis_boundary;
};
struct FullTimeHistory {
  std::shared_ptr<const FullTimeHistoryPolicyState> policy;
  std::shared_ptr<const FullTimeAnchorChain> anchors;
  std::shared_ptr<const FullTimeHistoryBase> base;
  std::shared_ptr<const FullTimeWindowNode> tail;
  std::uint64_t window_rows() const noexcept { return tail ? tail->rows : 0; }
  std::uint64_t window_bytes() const noexcept {
    return tail ? tail->bytes : 0;
  }
  const std::string &head() const {
    return anchors->anchors.empty() ? policy->root
                                    : anchors->anchors.back()->sha256;
  }
};
// v6 cumulative flows of an accepted snapshot (immutable, shared).
struct FullTimeFlows {
  bool has_last_birth = false;
  std::uint64_t last_birth_step = 0, last_birth_base = 0, last_birth_count = 0;
  std::string last_birth_plan_sha256, last_birth_token_sha256;
  std::uint64_t births = 0, birth_rows = 0, sunk_carriers = 0, sink_rows = 0;
  std::array<FullTimeSinkSurfaceTotals, 5> surfaces{};
  mutable std::string digest;
  FullTimeFlows() {
    for (int k = 0; k < 5; ++k) {
      surfaces[k].surface = static_cast<FullTimeSinkSurface>(k);
      surfaces[k].kind = (k == 0 || k == 1 || k == 3)
                             ? FullTimeSinkKind::collected_at_conductor
                             : FullTimeSinkKind::escaped_exterior;
    }
  }
};
namespace full_time_detail {
// ---------------------------------------------------------------------------
// History operations.
const FullTimeHistory &history(const ArenaSnapshot &s) {
  require(bool(s.full_time_history), "full-time snapshot without history");
  return *s.full_time_history;
}
// Appends one typed row to a fresh (not yet shared) snapshot. Explicit window
// caps refuse before any change; nothing is ever dropped.
std::string full_time_append(ArenaSnapshot &s, FullRow row) {
  const auto &h = history(s);
  row.previous_chain = s.chain;
  FullTimeStoredRow stored{encode_row(row), {}, row.kind};
  stored.sha256 = cn::sha256_hex(stored.bytes);
  const auto &p = h.policy->policy;
  if (h.window_rows() >= p.maximum_window_rows)
    throw refusal(RK::window_rows, "full-time history window row cap");
  if (stored.bytes.size() > p.maximum_window_bytes - h.window_bytes())
    throw refusal(RK::window_bytes, "full-time history window byte cap");
  const auto revision = checked_add(s.revision, 1, "full-time revision exhaustion");
  const auto bytes =
      checked_add(s.journal_bytes, stored.bytes.size(), "full-time byte total");
  auto next = std::make_shared<FullTimeHistory>(h);
  const auto sha = stored.sha256;
  next->tail = std::make_shared<const FullTimeWindowNode>(std::move(stored), h.tail);
  s.full_time_history = std::move(next);
  s.chain = sha;
  s.revision = revision;
  s.journal_bytes = bytes;
  s.digest_cache.clear();
  return sha;
}
std::vector<const FullTimeStoredRow *> window_rows(const FullTimeHistory &h) {
  std::vector<const FullTimeStoredRow *> rows;
  rows.reserve(h.window_rows());
  for (auto n = h.tail.get(); n; n = n->previous.get())
    rows.push_back(&n->row);
  std::reverse(rows.begin(), rows.end());
  return rows;
}
// Rows appended after `from_revision` (all inside the unanchored window).
std::vector<FullRow> rows_since(const ArenaSnapshot &s,
                                std::uint64_t from_revision) {
  require(from_revision <= s.revision, "row range");
  const auto n = s.revision - from_revision;
  const auto &h = history(s);
  require(n <= h.window_rows(), "requested rows precede the anchored window");
  std::vector<FullRow> rows;
  rows.reserve(n);
  auto node = h.tail.get();
  for (std::uint64_t i = 0; i < n; ++i, node = node->previous.get())
    rows.push_back(decode_row(node->row.bytes));
  std::reverse(rows.begin(), rows.end());
  return rows;
}
} // namespace full_time_detail
