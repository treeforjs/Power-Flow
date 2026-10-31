// Private total-field frames. Included in the arena TU to reuse the exact
// reaction builder and count/history codecs without a second native owner.
class FullTimeExposureKey {};
struct FullTimePhaseKey {};
namespace full_time_detail {
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
struct FullFrame {
  std::shared_ptr<const ArenaSnapshot> snapshot;
  PicArena::NativeState native;
};
void swap_frame(FullFrame &a, FullFrame &b) noexcept {
  a.snapshot.swap(b.snapshot);
  std::swap(a.native.options, b.native.options);
  a.native.particles.swap(b.native.particles);
  a.native.horizontal.swap(b.native.horizontal);
  a.native.vertical.swap(b.native.vertical);
  a.native.potential.Swap(b.native.potential);
  std::swap(a.native.ledger, b.native.ledger);
  std::swap(a.native.next_id, b.native.next_id);
  std::swap(a.native.fields_current, b.native.fields_current);
}
void release_native(PicArena::NativeState &n) noexcept {
  std::vector<GapParticle>().swap(n.particles);
  std::vector<double>().swap(n.horizontal);
  std::vector<double>().swap(n.vertical);
  mfem::Vector empty;
  n.potential.Swap(empty);
}
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
  require(j.size() == 10, "full-time common schema");
  FullTimeStepBinding b;
  b.accepted_step = json_u64(j.at("accepted_step"));
  b.accepted_arena_revision = json_u64(j.at("accepted_arena_revision"));
  b.begin_s = json_scalar(j.at("begin_s"));
  b.end_s = json_scalar(j.at("end_s"));
  b.dt_s = json_scalar(j.at("dt_s"));
  b.accepted_arena_sha256 = j.at("accepted_arena_sha256");
  b.accepted_field_sha256 = j.at("accepted_field_sha256");
  b.mesh_sha256 = j.at("mesh_sha256");
  b.material_sha256 = j.at("material_sha256");
  b.common_clock_sha256 = j.at("common_clock_sha256");
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
ExposureBinding read_exposure(const ArenaJson &j) {
  require(j.size() == 6, "full-time exposure schema");
  ExposureBinding b;
  b.common = read_common(j.at("common"));
  const auto slot = j.at("slot").get<std::string>();
  require(slot == "before_transport" || slot == "after_transport",
          "full-time slot");
  b.slot = slot == "before_transport" ? FullTimeSlot::before_transport
                                      : FullTimeSlot::after_transport;
  b.frozen_physical_stage_s = json_scalar(j.at("frozen_physical_stage_s"));
  b.collision_duration_s = json_scalar(j.at("collision_duration_s"));
  b.owner_stage_sequence = json_u64(j.at("owner_stage_sequence"));
  b.frozen_field_sha256 = j.at("frozen_field_sha256");
  return b;
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
  require(j.size() == 13, "full-time boundary schema");
  FullTimeAcceptedBoundary b;
  b.accepted_step = json_u64(j.at("accepted_step"));
  b.accepted_arena_revision = json_u64(j.at("accepted_arena_revision"));
  b.time_s = json_scalar(j.at("time_s"));
  b.accepted_arena_sha256 = j.at("accepted_arena_sha256");
  b.field_sha256 = j.at("field_sha256");
  b.heat_sha256 = j.at("heat_sha256");
  b.source_sha256 = j.at("source_sha256");
  b.rng_sha256 = j.at("rng_sha256");
  b.mesh_sha256 = j.at("mesh_sha256");
  b.material_sha256 = j.at("material_sha256");
  b.frame_sha256 = j.at("frame_sha256");
  b.common_clock_sha256 = j.at("common_clock_sha256");
  b.authored_full_time_genesis =
      json_boolean(j.at("authored_full_time_genesis"));
  return b;
}
void common_guard(const FullTimeStepBinding &b,
                  const FullTimeAcceptedBoundary &a) {
  require(std::isfinite(b.begin_s) && b.begin_s >= 0 &&
              std::isfinite(b.end_s) && std::isfinite(b.dt_s) && b.dt_s > 0 &&
              b.end_s > b.begin_s && b.end_s == b.begin_s + b.dt_s,
          "full-time authoritative common interval");
  require(b.accepted_step == a.accepted_step &&
              b.accepted_arena_revision == a.accepted_arena_revision &&
              b.begin_s == a.time_s &&
              b.accepted_arena_sha256 == a.accepted_arena_sha256 &&
              b.accepted_field_sha256 == a.field_sha256 &&
              b.mesh_sha256 == a.mesh_sha256 &&
              b.material_sha256 == a.material_sha256 &&
              digest(b.common_clock_sha256),
          "full-time accepted common binding");
}
void closed_position(const GapPicOptions &o, const GapParticle &p) {
  native_particle_guard(p);
  require(p.position_m[0] > 0 && p.position_m[0] < o.width_m &&
              p.position_m[1] > 0 && p.position_m[1] < o.gap_m,
          "closed full-time frame requires strict interior endpoints");
  if (o.finite_axial_extent)
    require(p.position_m[2] > o.short_z_m && p.position_m[2] < o.feed_z_m,
            "closed full-time axial endpoint");
  if (o.finite_foils)
    for (const auto &f : {o.lower_foil, o.upper_foil})
      require(!(p.position_m[0] >= f.left_m && p.position_m[0] <= f.right_m &&
                p.position_m[1] >= f.bottom_m && p.position_m[1] <= f.top_m),
              "closed full-time conductor endpoint");
}
std::string frame_prefix(const std::string &initial,
                         const std::string &registry, const FullFrame &f) {
  return sha_bytes("native-full-time-exact-config-v2" + initial + registry +
                   options_json(f.native.options).dump() +
                   metadata(*f.snapshot).dump() +
                   native_rows_json(f.native.particles).dump() +
                   std::to_string(f.native.next_id) +
                   ArenaJson(f.native.ledger.time_s).dump());
}
std::vector<ArenaRateEntry>
frame_population(const FullFrame &f,
                 const std::map<std::string, ArenaStateSpec> &states,
                 std::size_t cap) {
  require(f.native.particles.size() <= cap, "full-time native census cap");
  std::size_t n = 0;
  for (const auto &[id, p] : f.snapshot->packets) {
    (void)id;
    require(p.native_ids.size() <= cap - n, "full-time canonical census cap");
    n += p.native_ids.size();
  }
  auto index = validate_rows(*f.snapshot, f.native.particles, states);
  std::vector<ArenaRateEntry> out;
  out.reserve(n);
  for (const auto &[id, p] : f.snapshot->packets)
    for (auto native : p.native_ids)
      out.push_back({id, p.state, *index.at(native)});
  return out;
}
void exposure_guard(const ExposureBinding &e, const FullTimeStepBinding &common,
                    FullTimeSlot slot, const std::string &field) {
  require(common_json(e.common) == common_json(common) && e.slot == slot &&
              e.collision_duration_s == 0.5 * common.dt_s &&
              std::isfinite(e.collision_duration_s) &&
              e.collision_duration_s > 0 && e.frozen_field_sha256 == field &&
              e.frozen_physical_stage_s ==
                  (slot == FullTimeSlot::before_transport ? common.begin_s
                                                          : common.end_s),
          "full-time explicit half-exposure binding");
}
void full_ssa(const ArenaJson &journal, std::size_t start,
              const ArenaStochasticCompletion &c, const ExposureBinding &e,
              std::size_t end = std::numeric_limits<std::size_t>::max()) {
  require(c.stage_sequence == e.owner_stage_sequence &&
              c.common_clock_sha256 == e.common.common_clock_sha256 &&
              c.common_interval_begin_s == e.common.begin_s &&
              c.common_interval_end_s == e.common.end_s &&
              c.frozen_physical_stage_s == e.frozen_physical_stage_s,
          "full-time completion envelope");
  if (end == std::numeric_limits<std::size_t>::max())
    end = journal.size();
  require(journal.is_array() && start <= end && end <= journal.size(),
          "full-time SSA journal bounds");
  double previous = 0;
  for (std::size_t i = start; i < end; ++i) {
    const auto &a = journal[i].at("payload");
    require(a.at("kind") == "reaction", "full-time collision children only");
    require(a.at("phase") == "full_time_endpoint" &&
                a.at("physical_stage_s") == e.frozen_physical_stage_s &&
                a.at("clock_prefix") == e.common.common_clock_sha256 &&
                a.at("mesh_prefix") == e.common.mesh_sha256,
            "full-time reaction phase binding");
    if (a.at("selection_method") == "forced_kernel_integration_only")
      continue;
    require(a.at("selection_method") == "ssa_weighted_pair_clock_v1" &&
                a.at("channel") == c.channel_id &&
                a.at("dataset_sha256") == c.dataset_sha256 &&
                a.at("kernel_manifest_sha256") == c.kernel_manifest_sha256 &&
                a.at("kernel_binding_file_sha256") ==
                    c.kernel_binding_file_sha256 &&
                a.at("conditional_model") == c.model,
            "full-time SSA child channel binding");
    double elapsed = a.at("collision_elapsed_s");
    require(std::isfinite(elapsed) && elapsed > previous &&
                elapsed <= e.collision_duration_s,
            "full-time SSA ordered explicit half-duration");
    previous = elapsed;
  }
}
std::optional<ArenaStochasticCompletion>
frame_completion(const FullFrame &f, const std::string &owner) {
  auto i = f.snapshot->stochastic_completions.find(owner);
  if (i == f.snapshot->stochastic_completions.end())
    return {};
  auto c = completion_read(*i);
  require(c.integration_owner_id == owner, "completion owner key");
  return c;
}
FullTimeAcceptedBoundary boundary_payload(const ArenaJson &j) {
  require(j.is_object() && j.size() == 12 &&
              !j.contains("accepted_arena_sha256"),
          "full-time boundary payload schema");
  auto full = j;
  full["accepted_arena_sha256"] = std::string(64, '0');
  return read_boundary(full);
}
void full_journal_validate(
    const ArenaSnapshot &s, const ArenaSnapshot &genesis,
    const std::map<std::string, ArenaChannel> &channels) {
  require(s.journal.is_array() && s.journal.size() <= maximum_journal_entries &&
              s.revision == s.journal.size() && s.sinks == ArenaJson::array(),
          "bounded closed full-time journal dimensions");
  auto boundary = boundary_payload(genesis.full_time_context.at("boundary"));
  require(boundary.accepted_step == 0 &&
              boundary.accepted_arena_revision == 0 &&
              boundary.authored_full_time_genesis,
          "full-time journal genesis marker");
  std::string chain = genesis.chain, outer_sha, range_before_chain, motion_sha;
  std::set<std::string> events;
  std::size_t bytes = 0, start = 0;
  unsigned state = 0;
  std::uint64_t steps = 0;
  double previous_dt = 0;
  FullTimeStepBinding common;
  std::optional<ArenaStochasticCompletion> first, second;
  ArenaJson completions = ArenaJson::object();
  std::string field = boundary.field_sha256, heat = boundary.heat_sha256,
              source = boundary.source_sha256;
  for (std::size_t i = 0; i < s.journal.size(); ++i) {
    const auto &row = s.journal[i];
    require(row.is_object() && row.size() == 2, "full-time journal row schema");
    const auto &a = row.at("payload");
    const auto prior = chain;
    const auto n = row.dump().size();
    require(n <= maximum_journal_bytes - bytes,
            "bounded full-time journal bytes");
    bytes += n;
    require(a.is_object() && a.at("previous_chain") == chain &&
                row.at("sha256") == sha_bytes(a.dump()),
            "full-time journal hash chain");
    chain = json_text(row.at("sha256"));
    const auto kind = json_text(a.at("kind"));
    if (state == 0) {
      require(kind == "full_time_outer_begin" && a.size() == 4,
              "full-time journal requires outer begin");
      common = read_common(a.at("common"));
      require(std::isfinite(common.dt_s) && common.dt_s > 0 &&
                  common.end_s > common.begin_s &&
                  common.end_s == common.begin_s + common.dt_s &&
                  common.begin_s == boundary.time_s &&
                  common.accepted_step == steps &&
                  common.accepted_arena_revision == i &&
                  common.accepted_field_sha256 == field &&
                  common.mesh_sha256 == boundary.mesh_sha256 &&
                  common.material_sha256 == boundary.material_sha256 &&
                  digest(common.accepted_arena_sha256) &&
                  digest(common.common_clock_sha256),
              "full-time journal common continuity");
      outer_sha = json_text(a.at("outer_step_sha256"));
      require(digest(outer_sha), "full-time outer commitment");
      state = 1;
      start = i + 1;
      range_before_chain = chain;
      first.reset();
      second.reset();
      continue;
    }
    if ((state == 1 || state == 3) && kind == "reaction") {
      const auto event = json_text(a.at("event_id"));
      nonempty(event);
      require(events.insert(event).second, "duplicate full-time event");
      require(a.at("phase") == "full_time_endpoint" &&
                  json_u64(a.at("accepted_step")) == steps &&
                  json_scalar(a.at("physical_stage_s")) ==
                      (state == 1 ? common.begin_s : common.end_s) &&
                  json_scalar(a.at("previous_dt_s")) == previous_dt &&
                  json_scalar(a.at("next_dt_s")) == common.dt_s,
              "full-time journal reaction phase");
      require(a.at("products").is_array() && !a.at("products").empty() &&
                  a.at("products").size() <= 8,
              "full-time journal product cap");
      for (const auto &p : a.at("products")) {
        require(p.is_object() && p.size() == 4, "full-time product schema");
        (void)json_vector(p.at("u"));
      }
      require(digest(json_text(a.at("state_prefix"))) &&
                  a.at("selection_state_prefix") == a.at("state_prefix"),
              "full-time selection prefix");
      (void)json_scalar(a.at("event_weight"));
      (void)json_scalar(a.at("collision_elapsed_s"));
      continue;
    }
    if (state == 1 || state == 3) {
      require(kind == "full_time_exposure" && a.size() == 6,
              "full-time journal requires scoped exposure");
      auto e = read_exposure(a.at("exposure"));
      exposure_guard(e, common,
                     state == 1 ? FullTimeSlot::before_transport
                                : FullTimeSlot::after_transport,
                     field);
      auto c = completion_read(a.at("completion"));
      std::optional<ArenaStochasticCompletion> previous;
      auto found = completions.find(c.integration_owner_id);
      if (found != completions.end())
        previous = completion_read(*found);
      completion_validate(c, channels.at(c.channel_id), previous,
                          e.frozen_physical_stage_s);
      completion_chain_guard(e, c, boundary, state == 1 ? std::nullopt : first);
      full_ssa(s.journal, start, c, e, i);
      require(a.at("stage_before_chain") == range_before_chain &&
                  a.at("stage_after_reaction_chain") == prior,
              "full-time journal scoped reaction range");
      if (state == 1) {
        require(c.before_rng_clock_sha256 == boundary.rng_sha256,
                "full-time journal accepted chemistry RNG binding");
        first = c;
        state = 2;
      } else {
        require(first &&
                    c.integration_owner_id == first->integration_owner_id &&
                    c.channel_id == first->channel_id &&
                    first->stage_sequence <
                        std::numeric_limits<std::uint64_t>::max() &&
                    c.stage_sequence == first->stage_sequence + 1 &&
                    c.before_rng_clock_sha256 == first->after_rng_clock_sha256,
                "full-time journal second chemistry chain");
        second = c;
        state = 4;
      }
      completions[c.integration_owner_id] = completion_json(c);
      continue;
    }
    if (state == 2) {
      require(kind == "full_time_motion" && a.size() == 7,
              "full-time journal requires one motion");
      const auto &r = a.at("receipt");
      require(r.is_object() && r.size() == 10 &&
                  r.at("schema") == "arena_total_field_motion_receipt_v1" &&
                  json_scalar(r.at("begin_s")) == common.begin_s &&
                  json_scalar(r.at("end_s")) == common.end_s &&
                  json_scalar(r.at("dt_s")) == common.dt_s &&
                  r.at("common_clock_sha256") == common.common_clock_sha256 &&
                  r.at("mesh_sha256") == common.mesh_sha256 &&
                  r.at("field_before_sha256") == field,
              "full-time journal motion envelope");
      motion_sha = json_text(r.at("transport_trial_sha256"));
      require(digest(motion_sha), "full-time journal motion identity");
      field = json_text(a.at("field_after_sha256"));
      heat = json_text(a.at("heat_after_sha256"));
      require(digest(field) && digest(heat) &&
                  a.at("source_after_sha256") == source &&
                  r.at("field_after_sha256") == field,
              "full-time journal closed owner seals");
      require(a.at("endpoints").is_array() &&
                  a.at("endpoints").size() <= maximum_native_carriers,
              "full-time motion endpoint cap");
      state = 3;
      start = i + 1;
      range_before_chain = chain;
      continue;
    }
    require(state == 4 && kind == "full_time_boundary" && a.size() == 8 &&
                common_json(read_common(a.at("common"))) ==
                    common_json(common) &&
                a.at("outer_step_sha256") == outer_sha &&
                a.at("transport_trial_sha256") == motion_sha,
            "full-time journal terminal boundary grammar");
    auto next = boundary_payload(a.at("boundary"));
    require(steps < std::numeric_limits<std::uint64_t>::max(),
            "full-time journal step exhaustion");
    ++steps;
    require(first && second && next.accepted_step == steps &&
                next.accepted_arena_revision == i + 1 &&
                next.time_s == common.end_s &&
                !next.authored_full_time_genesis &&
                next.field_sha256 == field && next.heat_sha256 == heat &&
                next.source_sha256 == source &&
                next.rng_sha256 == second->after_rng_clock_sha256 &&
                next.mesh_sha256 == boundary.mesh_sha256 &&
                next.material_sha256 == boundary.material_sha256 &&
                next.frame_sha256 == boundary.frame_sha256 &&
                next.common_clock_sha256 == common.common_clock_sha256 &&
                digest(json_text(a.at("first_completed_state_prefix"))) &&
                digest(json_text(a.at("second_completed_state_prefix"))),
            "full-time journal final owner/clock/count seals");
    boundary = next;
    previous_dt = common.dt_s;
    state = 0;
  }
  require(state == 0, "full-time journal ends inside private step");
  require(bytes == s.journal_bytes && chain == s.chain &&
              events == s.event_ids && steps == s.accepted_steps &&
              previous_dt == s.previous_dt_s &&
              completions.dump() == s.stochastic_completions.dump(),
          "full-time journal counters/completion replay mismatch");
}
struct FullOuter {
  PicArena *owner = nullptr;
  std::weak_ptr<void> lifetime;
  FullTimeStepBinding common;
  FullFrame before, frame;
  std::shared_ptr<const FullTimeAcceptedBoundary> boundary;
  std::size_t cap = 0, records = 0;
  unsigned phase = 0;
  bool first_reserved = false, second_reserved = false, builder_active = false;
  std::size_t transport_forks = 0, slot_reaction_budget = 0;
  bool sealed = false;
  std::string outer_sha, motion_sha;
  FullTimeExposureReceipt first, second;
};
struct FullCollision {
  std::shared_ptr<FullOuter> outer;
  ExposureBinding exposure;
  std::shared_ptr<const FullTimeExposureKey> key =
      std::make_shared<FullTimeExposureKey>();
  std::shared_ptr<const FullTimePhaseKey> phase_key =
      std::make_shared<FullTimePhaseKey>();
  std::shared_ptr<void> token_identity = std::make_shared<int>(0);
  std::weak_ptr<void> outstanding_token;
  std::string completed_prefix;
  std::uint64_t generation = 0;
  bool building = false, closed = false;
  std::size_t limit = 0, base_revision = 0;
};
struct FullTrial {
  std::shared_ptr<FullOuter> outer;
  FullFrame base;
  std::string prefix;
  bool consumed = false;
};
// Hall record uses length-prefixed little-endian encoding, not bare JSON SHA.
struct FullEncoding {
  std::string bytes;
  void integer(std::uint64_t n) {
    for (int i = 0; i < 8; ++i)
      bytes.push_back(char((n >> (i * 8)) & 255));
  }
  void scalar(double x) {
    std::uint64_t u;
    std::memcpy(&u, &x, 8);
    integer(u);
  }
  void text(const std::string &s) {
    integer(s.size());
    bytes += s;
  }
};
std::string kinetic_policy_hash(const GapPicOptions &o) {
  FullEncoding e;
  e.bytes = "native-finite-target-kinetic-policy-v1";
  e.bytes.push_back('\0');
  for (double x :
       {o.width_m, o.gap_m, o.short_z_m, o.feed_z_m, o.lower_foil.left_m,
        o.lower_foil.right_m, o.lower_foil.bottom_m, o.lower_foil.top_m,
        o.upper_foil.left_m, o.upper_foil.right_m, o.upper_foil.bottom_m,
        o.upper_foil.top_m})
    e.scalar(x);
  return sha_bytes(e.bytes);
}
void segment_guard(const GapPicOptions &o, const PicVector &a,
                   const PicVector &b) {
  if (!o.finite_foils)
    return;
  for (const auto &f : {o.lower_foil, o.upper_foil}) {
    double enter = 0, leave = 1;
    bool possible = true;
    for (int d = 0; d < 2; ++d) {
      const double low = d ? f.bottom_m : f.left_m,
                   high = d ? f.top_m : f.right_m, D = b[d] - a[d];
      if (D == 0) {
        if (a[d] < low || a[d] > high)
          possible = false;
        continue;
      }
      double l = (low - a[d]) / D, h = (high - a[d]) / D;
      if (l > h)
        std::swap(l, h);
      enter = std::max(enter, l);
      leave = std::min(leave, h);
    }
    require(!possible || enter > leave,
            "closed motion path intersects native foil");
  }
}
std::string census_sha(const std::vector<ArenaRateEntry> &census) {
  std::vector<const ArenaRateEntry *> rows;
  for (const auto &r : census)
    rows.push_back(&r);
  std::sort(rows.begin(), rows.end(),
            [](auto a, auto b) { return a->native.id < b->native.id; });
  FullEncoding e;
  e.text("canonical-motion-census-v1");
  e.integer(rows.size());
  for (auto r : rows) {
    e.text(r->packet_id);
    e.text(r->state);
    const auto &p = r->native;
    e.integer(p.id);
    for (double x : p.position_m)
      e.scalar(x);
    for (double x : p.proper_velocity_m_s)
      e.scalar(x);
    e.scalar(p.charge_C);
    e.scalar(p.mass_kg);
    e.scalar(p.weight);
  }
  return sha_bytes(e.bytes);
}
} // namespace full_time_detail
using namespace full_time_detail;
struct CollisionScope::Impl {
  std::shared_ptr<FullCollision> scope;
};
struct PrivateFullTimeStep::Impl {
  std::shared_ptr<FullOuter> outer;
};
struct PrivateTransportTrial::Impl {
  std::shared_ptr<FullTrial> trial;
};
struct PreparedPrivateMotion::Impl {
  ~Impl() noexcept { drain_json(audit); }
  std::shared_ptr<FullTrial> trial;
  FullFrame after;
  ArenaJson audit;
  FullTimeMotionBinding binding;
  bool consumed = false;
};
class FullTimeAccess {
public:
  static void live(const FullOuter &o) {
    require(!o.lifetime.expired() && o.owner &&
                o.owner->full_time_boundary_ == o.boundary &&
                o.owner->accepted_ == o.before.snapshot &&
                PicArena::matches(*o.owner->pic_, o.before.native) && !o.sealed,
            "stale full-time outer owner");
  }
  static void live(const FullCollision &s) {
    live(*s.outer);
    require(
        !s.closed &&
            s.outer->phase ==
                (s.exposure.slot == FullTimeSlot::before_transport ? 0u : 2u),
        "closed or wrong collision scope");
  }
  static std::string prefix(const FullOuter &o, const FullFrame &f) {
    return frame_prefix(o.owner->initial_prefix_, o.owner->registry_sha_, f);
  }
  static std::vector<ArenaRateEntry>
  population(const FullOuter &o, const FullFrame &f, std::size_t cap) {
    return frame_population(f, o.owner->states_, std::min(cap, o.cap));
  }
  static PreparedArena reaction(FullCollision &s, const ArenaEvent &e,
                                int fault) {
    live(s);
    require(s.building &&
                s.outer->frame.snapshot->revision - s.base_revision < s.limit,
            "reaction outside bounded private compound");
    auto phase = phase_for(s);
    require(e.phase.kind == ArenaPhase::Kind::FullTimeEndpoint &&
                e.phase.full_time_key == s.phase_key &&
                e.phase.accepted_step == phase.accepted_step &&
                e.phase.previous_dt_s == phase.previous_dt_s &&
                e.phase.next_dt_s == phase.next_dt_s &&
                e.phase.mesh_prefix == phase.mesh_prefix &&
                e.phase.clock_prefix == phase.clock_prefix,
            "foreign full-time collision phase");
    auto &o = *s.outer;
    require(s.outstanding_token.expired(),
            "only one unconsumed private token may retain a particle frame");
    require(e.primary != e.target && !e.products.empty() &&
                e.products.size() <= 8,
            "bounded distinct reactants/products before resource arithmetic");
    const auto &primary = o.frame.snapshot->packets.at(e.primary);
    const auto &target = o.frame.snapshot->packets.at(e.target);
    auto pr = PicArena::debit(primary.count, e.event_weight),
         tr = PicArena::debit(target.count, e.event_weight);
    std::size_t n = o.frame.native.particles.size() -
                    primary.native_ids.size() - target.native_ids.size();
    auto add = [&](std::size_t amount) {
      require(amount <= o.cap - n,
              "outer reaction result carrier budget before allocation");
      n += amount;
    };
    add(PicArena::positive_carriers(pr).size());
    add(PicArena::positive_carriers(tr).size());
    add(e.products.size());
    auto built = o.owner->build_reaction(e, *o.frame.snapshot, o.frame.native,
                                         prefix(o, o.frame), true);
    require(!fault, "manufactured private reaction fault");
    auto token = std::make_unique<PreparedArena::Impl>();
    token->private_identity = s.token_identity;
    token->private_generation = s.generation;
    token->owner = o.owner;
    token->before = o.frame.snapshot;
    token->after = built.snapshot;
    token->native_before = o.frame.native;
    token->native_after = std::move(built.native);
    token->work = token->native_after;
    token->audit = std::move(built.audit);
    token->candidate_prefix = prefix(o, {token->after, token->native_after});
    token->private_memory_lease = std::make_shared<int>(0);
    s.outstanding_token = token->private_memory_lease;
    return PreparedArena(std::move(token));
  }
  static ArenaPhase phase_for(const FullCollision &s) {
    live(s);
    ArenaPhase p;
    p.kind = ArenaPhase::Kind::FullTimeEndpoint;
    p.accepted_step = s.exposure.common.accepted_step;
    p.physical_stage_s = s.exposure.frozen_physical_stage_s;
    p.previous_dt_s = s.outer->frame.snapshot->previous_dt_s;
    p.next_dt_s = s.exposure.common.dt_s;
    p.state_prefix = prefix(*s.outer, s.outer->frame);
    p.mesh_prefix = s.exposure.common.mesh_sha256;
    p.clock_prefix = s.exposure.common.common_clock_sha256;
    p.full_time_key = s.phase_key;
    return p;
  }
  static bool publishable(const FullCollision &s, const PreparedArena &t) {
    live(s);
    return t.impl_ && t.impl_->private_identity == s.token_identity &&
           t.impl_->private_generation == s.generation &&
           t.impl_->before == s.outer->frame.snapshot &&
           t.impl_->status == PreparedArena::Impl::Status::Staged;
  }
  static bool publish(FullCollision &s, PreparedArena &t) {
    if (!publishable(s, t))
      return false;
    auto &o = *s.outer;
    const bool final =
        t.impl_->audit.at("kind").get_ref<const std::string &>() ==
        "full_time_exposure";
    if ((!s.building && !final) || (s.building && final))
      return false;
    // All validation, hashing and allocation precede the first private swap.
    std::optional<FullTimeExposureReceipt> receipt;
    std::string completed;
    if (final) {
      receipt.emplace(FullTimeExposureReceipt{
          s.exposure, completion_read(t.impl_->audit.at("completion")),
          t.impl_->candidate_prefix,
          sha_bytes(t.impl_->audit.at("completion").dump()), s.key});
      completed = t.impl_->candidate_prefix;
    }
    FullFrame prepared{t.impl_->after, t.impl_->native_after};
    static_assert(std::is_nothrow_move_assignable_v<FullTimeExposureReceipt>);
    swap_frame(o.frame, prepared);
    t.impl_->status = PreparedArena::Impl::Status::Published;
    if (final) {
      s.closed = true;
      s.completed_prefix.swap(completed);
      if (s.exposure.slot == FullTimeSlot::before_transport) {
        o.first = std::move(*receipt);
        o.slot_reaction_budget = s.limit;
        o.phase = 1;
      } else {
        o.second = std::move(*receipt);
        o.phase = 3;
      }
    }
    // Consumed private tokens retain audit/prefix only, never another
    // population.
    release_native(t.impl_->native_before);
    release_native(t.impl_->native_after);
    release_native(t.impl_->work);
    t.impl_->before.reset();
    t.impl_->after.reset();
    t.impl_->private_memory_lease.reset();
    return true;
  }
  static PreparedArena compound(
      FullCollision &s, std::size_t maximum,
      const std::function<ArenaStochasticCompletion(CollisionScope &)> &build,
      CollisionScope &facade) {
    live(s);
    require(!s.building && !s.outer->builder_active && maximum > 0 &&
                maximum <= 4096 && bool(build),
            "private compound cap or reentry");
    auto &o = *s.outer;
    const std::size_t used =
        o.frame.snapshot->journal.size() - o.before.snapshot->journal.size();
    const bool first_slot = s.exposure.slot == FullTimeSlot::before_transport;
    if (first_slot) {
      require(used <= o.records && o.records - used >= 4 &&
                  maximum <= (o.records - used - 4) / 2,
              "reserve both slot reaction budgets before callbacks");
      require(o.frame.snapshot->journal.size() <= maximum_journal_entries - 4 &&
                  maximum <= (maximum_journal_entries -
                              o.frame.snapshot->journal.size() - 4) /
                                 2,
              "reserve both slot global history budgets before callbacks");
    } else
      require(maximum <= o.slot_reaction_budget && used <= o.records &&
                  o.records - used >= 2 && maximum <= o.records - used - 2,
              "second slot exceeds reserved reaction/history budget");
    require(s.outstanding_token.expired() &&
                s.generation <= std::numeric_limits<std::uint64_t>::max() - 2,
            "private token/generation resource gate");
    auto base = o.frame;
    s.building = true;
    o.builder_active = true;
    o.owner->full_time_builder_active_ = true;
    s.limit = maximum;
    s.base_revision = base.snapshot->revision;
    ++s.generation;
    try {
      auto c = build(facade);
      live(o);
      completion_chain_guard(
          s.exposure, c, *o.boundary,
          s.exposure.slot == FullTimeSlot::before_transport
              ? std::nullopt
              : std::optional<ArenaStochasticCompletion>(o.first.completion));
      completion_validate(c, o.owner->channels_.at(c.channel_id),
                          frame_completion(base, c.integration_owner_id),
                          o.frame.native.ledger.time_s);
      full_ssa(o.frame.snapshot->journal, base.snapshot->journal.size(), c,
               s.exposure);
      require(s.outstanding_token.expired(),
              "unconsumed private event token at completion");
      auto next = std::make_shared<ArenaSnapshot>(*o.frame.snapshot);
      next->stochastic_completions[c.integration_owner_id] = completion_json(c);
      ArenaJson a = {{"kind", "full_time_exposure"},
                     {"exposure", exposure_json(s.exposure)},
                     {"completion", completion_json(c)},
                     {"stage_before_chain", base.snapshot->chain},
                     {"stage_after_reaction_chain", o.frame.snapshot->chain}};
      journal_append(*next, a);
      require(next->journal.size() - o.before.snapshot->journal.size() <=
                  o.records,
              "outer record resource cap");
      auto t = std::make_unique<PreparedArena::Impl>();
      t->owner = o.owner;
      t->private_identity = s.token_identity;
      t->private_generation = s.generation + 1;
      t->before = base.snapshot;
      t->after = next;
      t->native_before = base.native;
      t->native_after = o.frame.native;
      t->work = t->native_after;
      t->audit = a;
      t->candidate_prefix = prefix(o, {next, t->native_after});
      t->private_memory_lease = std::make_shared<int>(0);
      swap_frame(o.frame, base);
      s.building = false;
      o.builder_active = false;
      if (!o.lifetime.expired())
        o.owner->full_time_builder_active_ = false;
      ++s.generation;
      s.outstanding_token = t->private_memory_lease;
      return PreparedArena(std::move(t));
    } catch (...) {
      swap_frame(o.frame, base);
      s.building = false;
      o.builder_active = false;
      if (!o.lifetime.expired())
        o.owner->full_time_builder_active_ = false;
      ++s.generation;
      throw;
    }
  }
  static std::unique_ptr<PreparedPrivateMotion::Impl>
  motion(const std::shared_ptr<FullTrial> &t, const FullTimeMotionBinding &b,
         const std::vector<GapParticle> &endpoints, const ArenaJson &receipt,
         const ArenaJson &sources, const ArenaJson &boundaries) {
    auto &o = *t->outer;
    live(o);
    require(o.phase == 1 && !t->consumed && sources == ArenaJson::array() &&
                boundaries == ArenaJson::array(),
            "closed motion source/boundary ownership unsupported");
    require(common_json(b.common) == common_json(o.common) &&
                b.field_before_sha256 == o.common.accepted_field_sha256 &&
                digest(b.transport_trial_sha256) &&
                digest(b.field_after_sha256) && digest(b.heat_after_sha256) &&
                b.source_after_sha256 == o.boundary->source_sha256,
            "motion common/owner seals");
    require(
        receipt.is_object() && receipt.size() == 10 &&
            receipt.at("schema") == "arena_total_field_motion_receipt_v1" &&
            receipt.at("transport_trial_sha256") == b.transport_trial_sha256 &&
            receipt.at("common_clock_sha256") == o.common.common_clock_sha256 &&
            receipt.at("mesh_sha256") == o.common.mesh_sha256 &&
            receipt.at("field_before_sha256") == b.field_before_sha256 &&
            receipt.at("field_after_sha256") == b.field_after_sha256 &&
            receipt.at("begin_s") == o.common.begin_s &&
            receipt.at("end_s") == o.common.end_s &&
            receipt.at("dt_s") == o.common.dt_s,
        "motion envelope binding");
    auto audit = receipt.at("owner_audit");
    const auto ownerseal = audit.at("seal_sha256").get<std::string>();
    audit.erase("seal_sha256");
    FullEncoding enc;
    enc.text("motion-receipt-record-v1");
    enc.text(audit.dump());
    require(ownerseal == b.transport_trial_sha256 &&
                sha_bytes(enc.bytes) == ownerseal,
            "motion genuine owner record seal");
    require(audit.is_object() && audit.size() == 15 &&
                audit.at("schema") ==
                    "canonical-full-time-total-field-motion-v2",
            "recognized Hall v2 record schema required");
    const auto &bind = audit.at("binding");
    require(bind.is_object() && bind.size() == 14, "Hall v2 binding schema");
    (void)json_u64(bind.at("accepted_step"));
    (void)json_u64(bind.at("staged_arena_revision"));
    (void)json_u64(bind.at("trial_stage"));
    for (const auto *key : {"begin_s", "end_s", "dt_s"})
      (void)json_scalar(bind.at(key));
    require(digest(json_text(bind.at("species_sha256"))) &&
                bind.at("kinetic_policy_sha256") ==
                    kinetic_policy_hash(t->base.native.options),
            "Hall native species/census and kinetic-policy identity");
    const auto &opt = t->base.native.options;
    require(opt.finite_foils && opt.finite_axial_extent,
            "Hall v2 finite kinetic geometry required");
    const auto &policy = audit.at("kinetic_policy");
    ArenaJson expected_policy = {
        {"schema", "native-finite-target-kinetic-policy-v1"},
        {"width_m", opt.width_m},
        {"gap_m", opt.gap_m},
        {"short_z_m", opt.short_z_m},
        {"feed_z_m", opt.feed_z_m},
        {"lower_foil",
         {opt.lower_foil.left_m, opt.lower_foil.right_m,
          opt.lower_foil.bottom_m, opt.lower_foil.top_m}},
        {"upper_foil",
         {opt.upper_foil.left_m, opt.upper_foil.right_m,
          opt.upper_foil.bottom_m, opt.upper_foil.top_m}},
        {"short_bridge", "lower foil x span and inner y gap only"},
        {"xy_policy", "absorbing escape"},
        {"feed_policy", "axial escape"},
        {"tie_order", "x y lower upper axial; later exact tie wins"}};
    require(policy.dump() == expected_policy.dump(),
            "Hall complete native kinetic policy operands");
    const PicVector shift = {-.5 * opt.width_m, -.5 * opt.gap_m, 0};
    require(json_vector(audit.at("native_to_field_m")) == shift,
            "Hall translation-only frame");
    for (const auto *key : {"current_A_m_sha256", "initial_charge_C_sha256",
                            "endpoint_charge_C_sha256"})
      require(digest(json_text(audit.at(key))),
              "Hall charge/current array identity");
    require(
        bind.at("phase") == "full_time_position_and_proper_velocity" &&
            bind.at("accepted_step") == o.common.accepted_step &&
            bind.at("staged_arena_revision") == t->base.snapshot->revision &&
            bind.at("staged_arena_sha256") == t->prefix &&
            bind.at("begin_s") == o.common.begin_s &&
            bind.at("end_s") == o.common.end_s &&
            bind.at("dt_s") == o.common.dt_s &&
            bind.at("mesh_sha256") == o.common.mesh_sha256 &&
            bind.at("material_sha256") == o.common.material_sha256 &&
            bind.at("common_clock_sha256") == o.common.common_clock_sha256 &&
            digest(bind.at("trial_field_sha256")),
        "motion owner staged identity");
    require(audit.at("census_sha256") ==
                    census_sha(population(o, t->base, o.cap)) &&
                audit.at("coverage") ==
                    "all_native_carriers_no_sources_no_boundaries" &&
                audit.at("neutral_motion_owner") == "same_full_time_q0_drift" &&
                !json_boolean(audit.at("native_Poisson_called")) &&
                !json_boolean(audit.at("native_DKD_called")),
            "motion coverage/census/mode");
    require(endpoints.size() == t->base.native.particles.size() &&
                endpoints.size() <= o.cap &&
                audit.at("endpoints").size() == endpoints.size(),
            "motion endpoint dimensions");
    std::map<std::uint64_t, GapParticle> previous, actual;
    for (const auto &p : t->base.native.particles)
      previous.emplace(p.id, p);
    for (const auto &p : endpoints) {
      closed_position(t->base.native.options, p);
      require(actual.emplace(p.id, p).second, "duplicate motion endpoint");
      const auto &old = previous.at(p.id);
      require(p.charge_C == old.charge_C && p.mass_kg == old.mass_kg &&
                  p.weight == old.weight,
              "motion altered native species/count");
    }
    require(audit.at("endpoints").is_array() && audit.at("paths").is_array() &&
                audit.at("paths").size() == endpoints.size(),
            "Hall complete endpoint/path census");
    std::set<std::uint64_t> audited_ids, path_ids;
    for (const auto &p : audit.at("endpoints")) {
      require(p.is_object() && p.size() == 6, "Hall endpoint row schema");
      const auto id = json_u64(p.at("native_id"));
      require(audited_ids.insert(id).second, "duplicate Hall endpoint row");
      const auto &n = actual.at(id);
      (void)json_vector(p.at("position_m"));
      (void)json_vector(p.at("proper_velocity_m_s"));
      for (const auto *k : {"charge_C", "mass_kg", "weight"})
        (void)json_scalar(p.at(k));
      require(p.at("position_m") == ArenaJson(n.position_m) &&
                  p.at("proper_velocity_m_s") ==
                      ArenaJson(n.proper_velocity_m_s) &&
                  p.at("charge_C") == n.charge_C &&
                  p.at("mass_kg") == n.mass_kg && p.at("weight") == n.weight,
              "motion endpoint disagrees with sealed owner");
    }
    for (const auto &path : audit.at("paths")) {
      require(path.is_object() && path.size() == 5 &&
                  path.at("motion_owner") == "ordinary_total_field",
              "Hall path row schema/owner");
      const auto id = json_u64(path.at("native_id"));
      require(path_ids.insert(id).second, "duplicate Hall path native ID");
      const auto &old = previous.at(id);
      const auto &endpoint = actual.at(id);
      const auto start = json_vector(path.at("field_start_m")),
                 D = json_vector(path.at("field_displacement_m"));
      require(json_scalar(path.at("flight_s")) == o.common.dt_s,
              "Hall complete common flight");
      for (int d = 0; d < 3; ++d)
        require(start[d] == old.position_m[d] + shift[d] &&
                    D[d] == (endpoint.position_m[d] + shift[d]) - start[d],
                "Hall canonical path start/endpoint operands");
      segment_guard(opt, old.position_m, endpoint.position_m);
    }
    require(audited_ids.size() == actual.size() && path_ids == audited_ids,
            "Hall endpoint/path ID bijection");
    const auto &budget = audit.at("budget");
    require(budget.is_object() && budget.size() == 7,
            "Hall complete budget schema");
    for (const auto *key :
         {"kinetic_change_J", "field_work_J", "work_defect_J", "continuity_C",
          "plasma_dt", "maximum_gyro", "maximum_displacement_cells"})
      (void)json_scalar(budget.at(key));
    long double kinetic_change = 0;
    for (const auto &[id, old] : previous)
      kinetic_change +=
          static_cast<long double>(old.weight) *
          (pic_kinetic_energy_J(actual.at(id)) - pic_kinetic_energy_J(old));
    const double K = json_scalar(budget.at("kinetic_change_J")),
                 W = json_scalar(budget.at("field_work_J")),
                 defect = json_scalar(budget.at("work_defect_J"));
    require(K == double(kinetic_change) && defect == K - W &&
                std::abs(defect) <= 1e-20 + 1e-10 * (std::abs(K) + std::abs(W)),
            "Hall retained exact kinetic/work defect and unchanged guard");
    require(json_scalar(budget.at("continuity_C")) >= 0 &&
                json_scalar(budget.at("plasma_dt")) >= 0 &&
                json_scalar(budget.at("plasma_dt")) <= .1 &&
                json_scalar(budget.at("maximum_gyro")) >= 0 &&
                json_scalar(budget.at("maximum_gyro")) <= .05 &&
                json_scalar(budget.at("maximum_displacement_cells")) >= 0 &&
                json_scalar(budget.at("maximum_displacement_cells")) <= .2,
            "Hall retained phase/plasma/gyro/displacement guards");
    auto result = std::make_unique<PreparedPrivateMotion::Impl>();
    result->trial = t;
    result->binding = b;
    result->after = t->base;
    result->after.native.particles.clear();
    for (const auto &old : t->base.native.particles)
      result->after.native.particles.push_back(actual.at(old.id));
    result->after.native.ledger.time_s = o.common.end_s;
    result->after.native.fields_current = false;
    auto next = std::make_shared<ArenaSnapshot>(*t->base.snapshot);
    next->full_time_context["field_sha256"] = b.field_after_sha256;
    next->full_time_context["heat_sha256"] = b.heat_after_sha256;
    next->full_time_context["source_sha256"] = b.source_after_sha256;
    result->audit = {
        {"kind", "full_time_motion"},
        {"receipt", receipt},
        {"endpoints", native_rows_json(result->after.native.particles)},
        {"field_after_sha256", b.field_after_sha256},
        {"heat_after_sha256", b.heat_after_sha256},
        {"source_after_sha256", b.source_after_sha256}};
    journal_append(*next, result->audit);
    result->after.snapshot = next;
    validate_rows(*next, result->after.native.particles, o.owner->states_);
    return result;
  }
};
FullTimeMotionAdmission::~FullTimeMotionAdmission() noexcept {
  drain_json(record_);
}
FullTimeMotionAdmission &
FullTimeMotionAdmission::operator=(FullTimeMotionAdmission &&other) noexcept {
  if (this != &other) {
    drain_json(record_);
    record_.swap(other.record_);
  }
  return *this;
}
CollisionScope::CollisionScope(std::unique_ptr<Impl> p) : impl_(std::move(p)) {}
CollisionScope::CollisionScope(CollisionScope &&) noexcept = default;
CollisionScope &CollisionScope::operator=(CollisionScope &&) noexcept = default;
CollisionScope::~CollisionScope() = default;
const ExposureBinding &CollisionScope::exposure() const {
  FullTimeAccess::live(*impl_->scope);
  return impl_->scope->exposure;
}
std::shared_ptr<const FullTimeExposureKey>
CollisionScope::exposure_key() const {
  FullTimeAccess::live(*impl_->scope);
  return impl_->scope->key;
}
ArenaPhase CollisionScope::phase() const {
  return FullTimeAccess::phase_for(*impl_->scope);
}
const GapPicOptions &CollisionScope::native_options() const noexcept {
  return impl_->scope->outer->before.native.options;
}
std::shared_ptr<const ArenaSnapshot> CollisionScope::snapshot() const {
  FullTimeAccess::live(*impl_->scope);
  return impl_->scope->outer->frame.snapshot;
}
std::vector<ArenaRateEntry>
CollisionScope::rate_population(std::size_t cap) const {
  FullTimeAccess::live(*impl_->scope);
  return FullTimeAccess::population(*impl_->scope->outer,
                                    impl_->scope->outer->frame, cap);
}
const ArenaChannel &
CollisionScope::channel_binding(const std::string &id) const {
  FullTimeAccess::live(*impl_->scope);
  return impl_->scope->outer->owner->channel_binding(id);
}
std::optional<ArenaStochasticCompletion>
CollisionScope::stochastic_completion(const std::string &id) const {
  FullTimeAccess::live(*impl_->scope);
  return frame_completion(impl_->scope->outer->frame, id);
}
std::string CollisionScope::state_prefix() const {
  if (impl_->scope->closed)
    return impl_->scope->completed_prefix;
  FullTimeAccess::live(*impl_->scope);
  return FullTimeAccess::prefix(*impl_->scope->outer,
                                impl_->scope->outer->frame);
}
PreparedArena CollisionScope::prepare(const ArenaEvent &e, int fault) {
  return FullTimeAccess::reaction(*impl_->scope, e, fault);
}
PreparedArena CollisionScope::prepare_compound_with_completion(
    std::size_t cap,
    const std::function<ArenaStochasticCompletion(CollisionScope &)> &f) {
  auto scope = impl_->scope;
  return FullTimeAccess::compound(*scope, cap, f, *this);
}
bool CollisionScope::can_publish(const PreparedArena &t) const noexcept {
  if (!impl_ || !impl_->scope)
    return false;
  try {
    return FullTimeAccess::publishable(*impl_->scope, t);
  } catch (...) {
    return false;
  }
}
bool CollisionScope::try_publish(PreparedArena &t) noexcept {
  if (!impl_ || !impl_->scope)
    return false;
  try {
    return FullTimeAccess::publish(*impl_->scope, t);
  } catch (...) {
    return false;
  }
}
PreparedPrivateMotion::PreparedPrivateMotion(std::unique_ptr<Impl> p)
    : impl_(std::move(p)) {}
PreparedPrivateMotion::PreparedPrivateMotion(
    PreparedPrivateMotion &&) noexcept = default;
PreparedPrivateMotion &
PreparedPrivateMotion::operator=(PreparedPrivateMotion &&) noexcept = default;
PreparedPrivateMotion::~PreparedPrivateMotion() = default;
const ArenaJson &PreparedPrivateMotion::audit() const {
  require(bool(impl_), "empty motion token");
  return impl_->audit;
}
PrivateTransportTrial::PrivateTransportTrial(std::unique_ptr<Impl> p)
    : impl_(std::move(p)) {}
PrivateTransportTrial::PrivateTransportTrial(
    PrivateTransportTrial &&) noexcept = default;
PrivateTransportTrial &
PrivateTransportTrial::operator=(PrivateTransportTrial &&) noexcept = default;
PrivateTransportTrial::~PrivateTransportTrial() = default;
std::vector<ArenaRateEntry>
PrivateTransportTrial::native_census(std::size_t cap) const {
  FullTimeAccess::live(*impl_->trial->outer);
  require(!impl_->trial->consumed && impl_->trial->outer->phase == 1,
          "stale transport census");
  return FullTimeAccess::population(*impl_->trial->outer, impl_->trial->base,
                                    cap);
}
std::shared_ptr<const ArenaSnapshot> PrivateTransportTrial::snapshot() const {
  (void)native_census(impl_->trial->outer->cap);
  return impl_->trial->base.snapshot;
}
std::uint64_t PrivateTransportTrial::staged_arena_revision() const {
  return snapshot()->revision;
}
std::string PrivateTransportTrial::state_prefix() const {
  (void)snapshot();
  return impl_->trial->prefix;
}
PreparedPrivateMotion PrivateTransportTrial::prepare_total_field_motion(
    const FullTimeMotionBinding &b, FullTimeMotionAdmission admission) {
  FullTimeAccess::live(*impl_->trial->outer);
  const auto &record = admission.record_;
  require(record.is_object() && record.contains("endpoints") &&
              record.at("endpoints").is_array() &&
              record.at("endpoints").size() <= impl_->trial->outer->cap,
          "bounded typed Hall receipt before copying");
  std::vector<GapParticle> endpoints;
  endpoints.reserve(record.at("endpoints").size());
  for (const auto &row : record.at("endpoints")) {
    GapParticle p;
    p.id = json_u64(row.at("native_id"));
    p.position_m = json_vector(row.at("position_m"));
    p.proper_velocity_m_s = json_vector(row.at("proper_velocity_m_s"));
    p.charge_C = json_scalar(row.at("charge_C"));
    p.mass_kg = json_scalar(row.at("mass_kg"));
    p.weight = json_scalar(row.at("weight"));
    endpoints.push_back(p);
  }
  ArenaJson envelope = {{"schema", "arena_total_field_motion_receipt_v1"},
                        {"transport_trial_sha256", b.transport_trial_sha256},
                        {"common_clock_sha256", b.common.common_clock_sha256},
                        {"mesh_sha256", b.common.mesh_sha256},
                        {"field_before_sha256", b.field_before_sha256},
                        {"field_after_sha256", b.field_after_sha256},
                        {"begin_s", b.common.begin_s},
                        {"end_s", b.common.end_s},
                        {"dt_s", b.common.dt_s},
                        {"owner_audit", std::move(admission.record_)}};
  return PreparedPrivateMotion(
      FullTimeAccess::motion(impl_->trial, b, endpoints, envelope,
                             ArenaJson::array(), ArenaJson::array()));
}
PreparedPrivateMotion PrivateTransportTrial::prepare_total_field_motion(
    const FullTimeMotionBinding &, const std::vector<GapParticle> &,
    const ArenaJson &, const ArenaJson &, const ArenaJson &) {
  throw std::invalid_argument("live JSON motion admission forbidden; use "
                              "genuine typed Hall MotionReceipt adapter");
}
PrivateFullTimeStep::PrivateFullTimeStep(std::unique_ptr<Impl> p)
    : impl_(std::move(p)) {}
PrivateFullTimeStep::PrivateFullTimeStep(PrivateFullTimeStep &&) noexcept =
    default;
PrivateFullTimeStep &
PrivateFullTimeStep::operator=(PrivateFullTimeStep &&) noexcept = default;
PrivateFullTimeStep::~PrivateFullTimeStep() = default;
CollisionScope PrivateFullTimeStep::first_exposure(const ExposureBinding &e) {
  auto o = impl_->outer;
  FullTimeAccess::live(*o);
  require(o->phase == 0 && !o->first_reserved && !o->builder_active,
          "first exposure already reserved or reentered");
  exposure_guard(e, o->common, FullTimeSlot::before_transport,
                 o->common.accepted_field_sha256);
  auto s = std::make_shared<FullCollision>();
  s->outer = o;
  s->exposure = e;
  auto p = std::make_unique<CollisionScope::Impl>();
  p->scope = std::move(s);
  o->first_reserved = true;
  return CollisionScope(std::move(p));
}
PrivateTransportTrial PrivateFullTimeStep::fork_transport_trial() {
  auto o = impl_->outer;
  FullTimeAccess::live(*o);
  require(o->phase == 1 && !o->builder_active && o->transport_forks < 4,
          "transport requires completed first exposure and bounded forks");
  auto p = std::make_unique<PrivateTransportTrial::Impl>();
  p->trial = std::make_shared<FullTrial>();
  p->trial->outer = o;
  p->trial->base = o->frame;
  p->trial->prefix = FullTimeAccess::prefix(*o, o->frame);
  ++o->transport_forks;
  return PrivateTransportTrial(std::move(p));
}
void PrivateFullTimeStep::accept_transport_trial(PreparedPrivateMotion &&t) {
  auto &o = *impl_->outer;
  FullTimeAccess::live(o);
  require(t.impl_ && t.impl_->trial->outer == impl_->outer &&
              !t.impl_->consumed && !t.impl_->trial->consumed && o.phase == 1 &&
              o.frame.snapshot == t.impl_->trial->base.snapshot,
          "foreign/stale/double motion token");
  require(t.impl_->after.snapshot->journal.size() -
                  o.before.snapshot->journal.size() <=
              o.records,
          "outer motion record cap");
  swap_frame(o.frame, t.impl_->after);
  o.motion_sha.swap(t.impl_->binding.transport_trial_sha256);
  o.phase = 2;
  t.impl_->consumed = true;
  t.impl_->trial->consumed = true;
}
CollisionScope PrivateFullTimeStep::second_exposure(const ExposureBinding &e) {
  auto o = impl_->outer;
  FullTimeAccess::live(*o);
  require(o->phase == 2 && !o->second_reserved && !o->builder_active,
          "second exposure already reserved or reentered");
  exposure_guard(e, o->common, FullTimeSlot::after_transport,
                 o->frame.snapshot->full_time_context.at("field_sha256"));
  auto s = std::make_shared<FullCollision>();
  s->outer = o;
  s->exposure = e;
  auto p = std::make_unique<CollisionScope::Impl>();
  p->scope = std::move(s);
  o->second_reserved = true;
  return CollisionScope(std::move(p));
}
PreparedArena PrivateFullTimeStep::seal(const FullTimeCommitSeals &seals) {
  auto &o = *impl_->outer;
  FullTimeAccess::live(o);
  require(o.phase == 3 && !o.builder_active,
          "outer seal needs both exact collision exposures");
  auto &owner = *o.owner;
  final_rng_guard(seals, o.second.completion);
  require(digest(o.motion_sha),
          "outer candidate lacks genuine consumed motion identity");
  require(seals.field_sha256 ==
                  o.frame.snapshot->full_time_context.at("field_sha256")
                      .get<std::string>() &&
              seals.heat_sha256 ==
                  o.frame.snapshot->full_time_context.at("heat_sha256")
                      .get<std::string>() &&
              seals.source_sha256 ==
                  o.frame.snapshot->full_time_context.at("source_sha256")
                      .get<std::string>() &&
              seals.rng_sha256 == o.second.completion.after_rng_clock_sha256,
          "outer final field/heat/source/RNG seals");
  auto boundary = std::make_shared<FullTimeAcceptedBoundary>(*o.boundary);
  require(boundary->accepted_step < std::numeric_limits<std::uint64_t>::max(),
          "outer step exhaustion");
  ++boundary->accepted_step;
  boundary->accepted_arena_revision = o.frame.snapshot->revision + 1;
  boundary->time_s = o.common.end_s;
  boundary->field_sha256 = seals.field_sha256;
  boundary->heat_sha256 = seals.heat_sha256;
  boundary->source_sha256 = seals.source_sha256;
  boundary->rng_sha256 = seals.rng_sha256;
  boundary->common_clock_sha256 = o.common.common_clock_sha256;
  boundary->authored_full_time_genesis = false;
  auto next = std::make_shared<ArenaSnapshot>(*o.frame.snapshot);
  next->accepted_steps = boundary->accepted_step;
  next->previous_dt_s = o.common.dt_s;
  next->full_time_context = {{"mode", "accepted_total_field_v1"},
                             {"boundary", boundary_json(*boundary, false)}};
  ArenaJson audit = {
      {"kind", "full_time_boundary"},
      {"common", common_json(o.common)},
      {"outer_step_sha256", o.outer_sha},
      {"transport_trial_sha256", o.motion_sha},
      {"first_completed_state_prefix", o.first.completed_state_prefix},
      {"second_completed_state_prefix", o.second.completed_state_prefix},
      {"boundary", boundary_json(*boundary, false)}};
  journal_append(*next, audit);
  require(next->journal.size() - o.before.snapshot->journal.size() <= o.records,
          "outer final record cap");
  auto t = std::make_unique<PreparedArena::Impl>();
  t->owner = &owner;
  t->before = o.before.snapshot;
  t->after = next;
  t->native_before = o.before.native;
  t->native_after = o.frame.native;
  t->work = t->native_after;
  t->audit = audit;
  t->boundary_before = o.boundary;
  t->boundary_after = boundary;
  t->full_outer_identity = impl_->outer;
  t->candidate_prefix = FullTimeAccess::prefix(o, {next, t->native_after});
  boundary->accepted_arena_sha256 = t->candidate_prefix;
  auto exposures = std::make_shared<FullTimeCandidateExposures>();
  exposures->common = o.common;
  exposures->outer_step_sha256 = o.outer_sha;
  exposures->transport_trial_sha256 = o.motion_sha;
  exposures->first = o.first;
  exposures->second = o.second;
  t->exposures = exposures;
  o.sealed = true;
  return PreparedArena(std::move(t));
}
std::optional<ArenaStochasticCompletion>
PreparedArena::candidate_stochastic_completion(const std::string &owner) const {
  require(bool(impl_), "empty prepared token");
  require(bool(impl_->after),
          "consumed private token retains audit/prefix only");
  return frame_completion({impl_->after, impl_->native_after}, owner);
}
std::shared_ptr<const FullTimeAcceptedBoundary>
PreparedArena::candidate_full_time_boundary() const {
  require(bool(impl_), "empty prepared token");
  return impl_->boundary_after;
}
std::shared_ptr<const FullTimeAcceptedBoundary>
PicArena::full_time_boundary() const {
  require(bool(full_time_boundary_), "arena is not full-time total-field mode");
  return full_time_boundary_;
}
PrivateFullTimeStep PicArena::begin_full_time_step(const FullTimeStepBinding &b,
                                                   std::size_t carriers,
                                                   std::size_t records) {
  require(full_time_boundary_ && !active_scope_ &&
              active_full_time_outer_.expired(),
          "full-time outer requires accepted full-time boundary");
  common_guard(b, *full_time_boundary_);
  require(carriers > 0 && carriers <= maximum_native_carriers && records >= 7 &&
              records <= maximum_journal_entries,
          "outer bounded resources");
  (void)rate_population(carriers);
  auto o = std::make_shared<FullOuter>();
  o->owner = this;
  o->lifetime = full_time_lifetime_;
  o->common = b;
  o->boundary = full_time_boundary_;
  o->before = {accepted_, capture(*pic_)};
  o->frame = o->before;
  o->cap = carriers;
  o->records = records;
  o->outer_sha = sha_bytes("arena-full-time-outer-v1" + state_prefix() +
                           common_json(b).dump());
  auto next = std::make_shared<ArenaSnapshot>(*accepted_);
  next->full_time_context = {
      {"mode", "private_total_field_outer_v1"},
      {"common", common_json(b)},
      {"field_sha256", b.accepted_field_sha256},
      {"heat_sha256", full_time_boundary_->heat_sha256},
      {"source_sha256", full_time_boundary_->source_sha256}};
  journal_append(*next, {{"kind", "full_time_outer_begin"},
                         {"common", common_json(b)},
                         {"outer_step_sha256", o->outer_sha}});
  o->frame.snapshot = next;
  auto p = std::make_unique<PrivateFullTimeStep::Impl>();
  p->outer = o;
  active_full_time_outer_ = o;
  return PrivateFullTimeStep(std::move(p));
}
std::unique_ptr<PicArena>
PicArena::create_full_time(GapPicOptions options,
                           std::vector<ArenaStateSpec> states,
                           std::vector<ArenaChannel> channels,
                           std::vector<FullTimeAuthoredPacket> packets,
                           const FullTimeGenesisBinding &g) {
  require(!options.distributed_fields && std::isfinite(g.time_s) &&
              g.time_s >= 0,
          "full-time fresh serial origin");
  nonempty(g.origin_provenance);
  for (const auto *s : {&g.field_sha256, &g.mesh_sha256, &g.material_sha256,
                        &g.frame_sha256, &g.common_clock_sha256, &g.heat_sha256,
                        &g.source_sha256, &g.rng_sha256})
    require(digest(*s), "full-time genesis requires sealed owners/frame");
  options_guard(options);
  require(states.size() <= 256 && channels.size() <= 128 &&
              packets.size() <= maximum_native_carriers,
          "bounded full-time genesis catalog/population");
  auto pic = std::make_unique<GapPic>(options);
  std::map<std::string, ArenaStateSpec> registry;
  for (const auto &s : states)
    require(registry.emplace(s.id, s).second, "duplicate genesis species");
  std::vector<GapParticle> native;
  std::vector<ArenaPacket> logical;
  require(packets.size() <= maximum_native_carriers, "genesis packet cap");
  std::uint64_t id = 1;
  for (const auto &p : packets) {
    nonempty(p.origin);
    nonempty(p.id);
    const auto &spec = registry.at(p.state);
    require(!p.count.zero(), "zero authored full-time packet");
    ArenaPacket r{p.id, p.state, p.origin, "", p.count, {}};
    auto weights = positive_carriers(p.count);
    require(weights.size() <= maximum_native_carriers - native.size(),
            "genesis native carrier cap");
    for (double w : weights) {
      GapParticle n;
      n.id = id++;
      n.position_m = p.position_m;
      n.proper_velocity_m_s = p.proper_velocity_m_s;
      n.charge_C = spec.charge_e * pic_e;
      n.mass_kg = spec.mass_kg;
      n.weight = w;
      closed_position(options, n);
      native.push_back(n);
      r.native_ids.push_back(n.id);
    }
    logical.push_back(std::move(r));
  }
  pic->particles_ = native;
  pic->next_particle_id_ = id;
  pic->ledger_.time_s = g.time_s;
  for (const auto &p : native) {
    pic->ledger_.authored_charge_C += p.weight * p.charge_C;
    pic->ledger_.authored_kinetic_energy_J +=
        p.weight * pic_kinetic_energy_J(p);
    require(std::isfinite(pic->ledger_.authored_charge_C) &&
                std::isfinite(pic->ledger_.authored_kinetic_energy_J),
            "authored aggregate ledger overflow");
  }
  pic->fields_current_ = false;
  pic->horizontal_transfers_.clear();
  pic->vertical_transfers_.clear();
  pic->potential_.SetSize(0);
  ArenaJson origin = {{"schema", "authored-full-time-genesis-v1"},
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
  auto owner = std::unique_ptr<PicArena>(
      new PicArena(*pic, std::move(states), std::move(channels),
                   std::move(logical), sha_bytes(origin.dump()), true));
  owner->owned_pic_ = std::move(pic);
  auto boundary = std::make_shared<FullTimeAcceptedBoundary>();
  boundary->time_s = g.time_s;
  boundary->field_sha256 = g.field_sha256;
  boundary->mesh_sha256 = g.mesh_sha256;
  boundary->material_sha256 = g.material_sha256;
  boundary->frame_sha256 = g.frame_sha256;
  boundary->common_clock_sha256 = g.common_clock_sha256;
  boundary->heat_sha256 = g.heat_sha256;
  boundary->source_sha256 = g.source_sha256;
  boundary->rng_sha256 = g.rng_sha256;
  boundary->authored_full_time_genesis = true;
  auto snapshot = std::make_shared<ArenaSnapshot>(*owner->accepted_);
  snapshot->full_time_context = {{"mode", "accepted_total_field_v1"},
                                 {"boundary", boundary_json(*boundary, false)},
                                 {"genesis", origin}};
  owner->accepted_ = snapshot;
  owner->genesis_ = snapshot;
  boundary->accepted_arena_sha256 =
      frame_prefix(owner->initial_prefix_, owner->registry_sha_,
                   {snapshot, capture(*owner->pic_)});
  owner->full_time_boundary_ = boundary;
  return owner;
}
namespace {
ArenaJson clean_payload(ArenaJson a) {
  a.erase("previous_chain");
  return a;
}
ArenaEvent replay_full_event(const ArenaJson &a, CollisionScope &scope) {
  ArenaEvent e;
  e.id = a.at("event_id");
  e.channel = a.at("channel");
  e.primary = a.at("primary");
  e.target = a.at("target");
  e.untouched_primary_id = a.at("untouched_primary_id");
  e.event_weight = json_scalar(a.at("event_weight"));
  e.phase = scope.phase();
  require(a.at("phase") == "full_time_endpoint" &&
              a.at("state_prefix") == e.phase.state_prefix &&
              a.at("selection_state_prefix") == e.phase.state_prefix &&
              a.at("physical_stage_s") == e.phase.physical_stage_s &&
              a.at("accepted_step") == e.phase.accepted_step &&
              a.at("previous_dt_s") == e.phase.previous_dt_s &&
              a.at("next_dt_s") == e.phase.next_dt_s &&
              a.at("mesh_prefix") == e.phase.mesh_prefix &&
              a.at("clock_prefix") == e.phase.clock_prefix,
          "full-time replay exact phase/prefix");
  const auto kind = a.at("provenance_kind").get<std::string>();
  require(kind == "manufactured_forced" || kind == "conditional_physical_data",
          "full-time replay provenance kind");
  e.provenance.kind = kind == "conditional_physical_data"
                          ? ArenaProvenance::Kind::ConditionalPhysicalData
                          : ArenaProvenance::Kind::ManufacturedForced;
  e.provenance.dataset_sha256 = a.at("dataset_sha256");
  e.provenance.kernel_manifest_sha256 = a.at("kernel_manifest_sha256");
  e.provenance.kernel_binding_file_sha256 = a.at("kernel_binding_file_sha256");
  e.provenance.model = a.at("conditional_model");
  e.provenance.selection_method = a.at("selection_method");
  e.provenance.selection_state_prefix = a.at("selection_state_prefix");
  e.provenance.collision_elapsed_s = json_scalar(a.at("collision_elapsed_s"));
  for (const auto &p : a.at("products")) {
    ArenaProduct product;
    product.id = p.at("id");
    product.state = p.at("state");
    product.proper_velocity_m_s = json_vector(p.at("u"));
    const auto at = p.at("location").get<std::string>();
    require(at == "primary" || at == "target",
            "full-time replay product location");
    product.location =
        at == "primary" ? ArenaProduct::At::Primary : ArenaProduct::At::Target;
    if (product.id == e.primary)
      e.scattered_primary_id = e.primary;
    e.products.push_back(product);
  }
  return e;
}
FullTimeGenesisBinding read_full_genesis(const ArenaJson &j) {
  require(j.size() == 11 && j.at("schema") == "authored-full-time-genesis-v1",
          "full-time genesis schema");
  FullTimeGenesisBinding g;
  g.time_s = json_scalar(j.at("time_s"));
  g.origin_provenance = j.at("origin_provenance");
  g.field_sha256 = j.at("field_sha256");
  g.mesh_sha256 = j.at("mesh_sha256");
  g.material_sha256 = j.at("material_sha256");
  g.frame_sha256 = j.at("frame_sha256");
  g.common_clock_sha256 = j.at("common_clock_sha256");
  g.heat_sha256 = j.at("heat_sha256");
  g.source_sha256 = j.at("source_sha256");
  g.rng_sha256 = j.at("rng_sha256");
  return g;
}
} // namespace
void PicArena::save_full_time(const std::string &path) const {
  require(full_time_boundary_ && !active_scope_ && owned_pic_ &&
              accepted_->full_time_context.at("mode") ==
                  "accepted_total_field_v1",
          "only accepted complete full-time boundary may be saved");
  validate_population(*accepted_);
  validate_full_time_history();
  ArenaJson body = {{"schema", "canonical_native_arena_full_time_exact_v2"},
                    {"source_sha256", source_identity()},
                    {"registry", registry()},
                    {"registry_sha256", registry_sha_},
                    {"initial_coupling_prefix", initial_prefix_},
                    {"options", options_json(pic_->options_)},
                    {"ABI", {80, 168, 200, sizeof(GapPic)}},
                    {"genesis", metadata(*genesis_)},
                    {"genesis_native", native_rows_json(genesis_native_)},
                    {"genesis_next_id", genesis_next_id_},
                    {"metadata", metadata(*accepted_)},
                    {"native", native_rows_json(pic_->particles_)},
                    {"next_id", pic_->next_particle_id_},
                    {"boundary", boundary_json(*full_time_boundary_)}};
  require(body.dump().size() <= 64 * 1024 * 1024 - 1024,
          "bounded full-time output checkpoint bytes");
  atomic_json({{"payload", body}, {"sha256", sha_bytes(body.dump())}}, path);
}
std::unique_ptr<PicArena>
PicArena::restore_full_time(const std::string &path, GapPicOptions options,
                            std::vector<ArenaStateSpec> states,
                            std::vector<ArenaChannel> channels,
                            const FullTimeAcceptedBoundary &expected) {
  const auto doc = read_full_json(path);
  require(doc.is_object() && doc.size() == 2, "full-time restart envelope");
  const auto &b = doc.at("payload");
  require(b.size() == 14 &&
              b.at("schema") == "canonical_native_arena_full_time_exact_v2" &&
              b.at("source_sha256") == source_identity() &&
              doc.at("sha256") == sha_bytes(b.dump()),
          "full-time restart identity/integrity");
  require(b.at("ABI") == ArenaJson::array({80, 168, 200, sizeof(GapPic)}),
          "full-time restart ABI");
  const auto saved = read_options(b.at("options"));
  options_guard(options);
  require(options_equal(saved, options),
          "full-time restart kinetic geometry/options");
  require(boundary_json(read_boundary(b.at("boundary"))) ==
              boundary_json(expected),
          "full-time expected accepted owner seals");
  auto genesis = read_metadata(b.at("genesis"));
  auto rows = native_rows_read(b.at("genesis_native"));
  auto g = read_full_genesis(genesis.full_time_context.at("genesis"));
  std::map<std::uint64_t, GapParticle> index;
  for (const auto &p : rows)
    require(index.emplace(p.id, p).second,
            "full-time restart genesis native ID");
  std::vector<const ArenaPacket *> order;
  for (const auto &[id, p] : genesis.packets) {
    (void)id;
    require(!p.native_ids.empty(), "empty full-time genesis carriers");
    order.push_back(&p);
  }
  std::sort(order.begin(), order.end(), [](auto a, auto c) {
    return a->native_ids.front() < c->native_ids.front();
  });
  std::vector<FullTimeAuthoredPacket> authored;
  for (const auto *p : order) {
    const auto &native = index.at(p->native_ids.front());
    authored.push_back({p->id, p->state, p->origin, p->count, native.position_m,
                        native.proper_velocity_m_s});
  }
  auto owner = create_full_time(options, std::move(states), std::move(channels),
                                std::move(authored), g);
  require(
      owner->registry().dump() == b.at("registry").dump() &&
          owner->registry_sha_ == b.at("registry_sha256").get<std::string>() &&
          owner->initial_prefix_ ==
              b.at("initial_coupling_prefix").get<std::string>() &&
          metadata(*owner->genesis_) == b.at("genesis") &&
          native_rows_json(owner->genesis_native_) == b.at("genesis_native") &&
          owner->genesis_next_id_ == json_u64(b.at("genesis_next_id")),
      "full-time restart immutable genesis/catalog");
  const auto final = read_metadata(b.at("metadata"));
  full_journal_validate(final, *owner->genesis_, owner->channels_);
  std::size_t i = 0;
  auto payload = [&](std::size_t k) -> const ArenaJson & {
    require(k < final.journal.size(), "incomplete full-time step journal");
    return final.journal[k].at("payload");
  };
  while (i < final.journal.size()) {
    require(payload(i).at("kind") == "full_time_outer_begin",
            "full-time restart begins at outer boundary");
    const auto common = read_common(payload(i).at("common"));
    const auto outer_start = i;
    auto outer = owner->begin_full_time_step(common, maximum_native_carriers,
                                             maximum_journal_entries);
    require(clean_payload(payload(i)) ==
                ArenaJson({{"kind", "full_time_outer_begin"},
                           {"common", common_json(common)},
                           {"outer_step_sha256",
                            sha_bytes("arena-full-time-outer-v1" +
                                      owner->state_prefix() +
                                      common_json(common).dump())}}),
            "full-time restart outer identity");
    ++i;
    std::size_t replay_limit = 1, batch = 0;
    for (std::size_t look = i; look < final.journal.size(); ++look) {
      const auto k = payload(look).at("kind").get<std::string>();
      if (k == "full_time_boundary")
        break;
      if (k == "reaction")
        ++batch;
      else if (k == "full_time_exposure") {
        replay_limit = std::max(replay_limit, batch);
        batch = 0;
      }
    }
    auto replay_exposure = [&](FullTimeSlot slot) {
      std::size_t end = i;
      while (end < final.journal.size() &&
             payload(end).at("kind") == "reaction")
        ++end;
      require(payload(end).at("kind") == "full_time_exposure",
              "full-time restart missing exposure completion");
      auto e = read_exposure(payload(end).at("exposure"));
      require(e.slot == slot, "full-time restart slot ordering");
      auto scope = slot == FullTimeSlot::before_transport
                       ? outer.first_exposure(e)
                       : outer.second_exposure(e);
      auto t = scope.prepare_compound_with_completion(
          replay_limit, [&](CollisionScope &s) {
            for (; i < end; ++i) {
              auto ev = replay_full_event(payload(i), s);
              auto event = s.prepare(ev);
              require(event.audit() == clean_payload(payload(i)),
                      "full-time replay reaction ledgers");
              require(s.try_publish(event),
                      "full-time replay private reaction publication");
            }
            return completion_read(payload(end).at("completion"));
          });
      require(t.audit() == clean_payload(payload(end)),
              "full-time replay exposure receipt");
      require(scope.try_publish(t),
              "full-time replay private exposure publication");
      i = end + 1;
    };
    replay_exposure(FullTimeSlot::before_transport);
    const auto &motion = payload(i);
    require(motion.at("kind") == "full_time_motion",
            "full-time restart missing motion");
    auto trial = outer.fork_transport_trial();
    const auto &wrapper = motion.at("receipt");
    FullTimeMotionBinding mb;
    mb.common = common;
    mb.transport_trial_sha256 = wrapper.at("transport_trial_sha256");
    mb.field_before_sha256 = wrapper.at("field_before_sha256");
    mb.field_after_sha256 = motion.at("field_after_sha256");
    mb.heat_after_sha256 = motion.at("heat_after_sha256");
    mb.source_after_sha256 = motion.at("source_after_sha256");
    auto moved = PreparedPrivateMotion(FullTimeAccess::motion(
        trial.impl_->trial, mb, native_rows_read(motion.at("endpoints")),
        wrapper, ArenaJson::array(), ArenaJson::array()));
    require(moved.audit() == clean_payload(motion),
            "full-time replay motion receipt");
    outer.accept_transport_trial(std::move(moved));
    ++i;
    replay_exposure(FullTimeSlot::after_transport);
    const auto &boundary = payload(i);
    require(boundary.at("kind") == "full_time_boundary",
            "full-time restart missing terminal boundary");
    const auto &seals = boundary.at("boundary");
    auto t = outer.seal({seals.at("field_sha256"), seals.at("heat_sha256"),
                         seals.at("source_sha256"), seals.at("rng_sha256")});
    require(t.audit() == clean_payload(boundary),
            "full-time replay exact slot prefixes/final seals");
    require(owner->try_publish(t), "full-time replay outer publication");
    ++i;
    require(i - outer_start <= maximum_journal_entries,
            "full-time replay record bound");
  }
  require(metadata(*owner->accepted_) == b.at("metadata") &&
              native_rows_json(owner->pic_->particles_) == b.at("native") &&
              owner->pic_->next_particle_id_ == json_u64(b.at("next_id")) &&
              boundary_json(*owner->full_time_boundary_) == b.at("boundary") &&
              boundary_json(*owner->full_time_boundary_) ==
                  boundary_json(expected),
          "full-time restart final native/count/history/owner seals");
  return owner;
}
void PicArena::validate_full_time_history() const {
  require(full_time_boundary_ && owned_pic_ && genesis_ &&
              genesis_->journal.empty() && genesis_->revision == 0 &&
              genesis_->accepted_steps == 0 && genesis_->previous_dt_s == 0,
          "full-time immutable genesis");
  full_journal_validate(*accepted_, *genesis_, channels_);
  validate_rows(*genesis_, genesis_native_, states_);
  validate_population(*accepted_);
  require(full_time_boundary_->accepted_step == accepted_->accepted_steps &&
              full_time_boundary_->accepted_arena_revision ==
                  accepted_->revision &&
              full_time_boundary_->time_s == pic_->ledger_.time_s &&
              full_time_boundary_->accepted_arena_sha256 ==
                  frame_prefix(initial_prefix_, registry_sha_,
                               {accepted_, capture(*pic_)}) &&
              !pic_->fields_current_,
          "full-time accepted boundary/current native fingerprint");
  if (accepted_->journal.empty())
    require(full_time_boundary_->authored_full_time_genesis,
            "full-time empty journal genesis");
  else
    require(accepted_->journal.back().at("payload").at("kind") ==
                "full_time_boundary",
            "full-time only complete outer checkpoint");
}
