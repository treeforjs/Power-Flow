// Private total-field frames (native v5 + v6 births/sinks/binary motion v3).
// Included at the end of the arena TU to reuse the exact reaction builder and
// count codecs without a second owner.
class FullTimeExposureKey {};
struct FullTimePhaseKey {};
#include "mfem_pic_full_time_flows.inc"
namespace full_time_detail {
// ---------------------------------------------------------------------------
// Cached digests of immutable pieces.
const std::string &packets_digest(const SharedMap<std::string, ArenaPacket> &m) {
  if (!m.cached_digest().empty())
    return m.cached_digest();
  cn::Hasher h;
  h.w.tag("native-full-time-packets-v5");
  h.w.u64(m.size());
  for (const auto &[id, p] : m) {
    require(id == p.id, "packet map key/identity");
    h.w.text(p.id);
    h.w.text(p.state);
    h.w.text(p.origin);
    h.w.text(p.last_event);
    codec_count(h.w, const_cast<ArenaCount &>(p.count));
    h.w.u64(p.native_ids.size());
    for (auto n : p.native_ids)
      h.w.u64(n);
    h.flush_if_large();
  }
  m.cache_digest(cn::hex(h.finish()));
  return m.cached_digest();
}
const std::string &id_set_digest(const SharedSet<std::string> &s) {
  if (!s.cached_digest().empty())
    return s.cached_digest();
  cn::Hasher h;
  h.w.tag("native-v5/string-set");
  h.w.u64(s.size());
  for (const auto &x : s) {
    h.w.text(x);
    h.flush_if_large();
  }
  s.cache_digest(cn::hex(h.finish()));
  return s.cached_digest();
}
std::string completions_digest(const ArenaJson &j) {
  require(j.is_object() && j.size() <= 128, "bounded completion owners");
  Writer w;
  w.tag("native-full-time-completions-v5");
  w.u64(j.size());
  for (auto it = j.begin(); it != j.end(); ++it) {
    auto c = completion_read(it.value());
    require(c.integration_owner_id == it.key(), "completion owner key");
    w.text(it.key());
    codec_completion(w, c);
  }
  return cn::sha256_hex(w.bytes);
}
std::string context_digest(const ArenaJson &j) {
  Writer w;
  w.tag("native-full-time-context-v5");
  if (j.is_null()) {
    w.u8(0);
  } else {
    require(j.is_object() && j.contains("mode"), "full-time context schema");
    const auto mode = json_text(j.at("mode"));
    if (mode == "accepted_total_field_v1") {
      require(j.size() == 2 || (j.size() == 3 && j.contains("genesis")),
              "accepted context schema");
      w.u8(1);
      auto b = boundary_payload(j.at("boundary"));
      codec_boundary(w, b, false);
      w.boolean(j.contains("genesis"));
      if (j.contains("genesis")) {
        auto g = read_full_genesis(j.at("genesis"));
        codec_genesis(w, g);
      }
    } else {
      require(mode == "private_total_field_outer_v1" && j.size() == 5,
              "private outer context schema");
      w.u8(2);
      auto c = read_common(j.at("common"));
      codec_common(w, c);
      for (const auto *k : {"field_sha256", "heat_sha256", "source_sha256"})
        w.text(json_text(j.at(k)));
    }
  }
  return cn::sha256_hex(w.bytes);
}
std::string ledger_digest(const GapPicDiagnostics &d) {
  static_assert(sizeof(GapPicDiagnostics) == 25 * sizeof(double) &&
                std::is_standard_layout_v<GapPicDiagnostics>);
  double v[25];
  std::memcpy(v, &d, sizeof v);
  Writer w;
  w.tag("native-gap-ledger-v5");
  for (double x : v)
    w.f64(x);
  return cn::sha256_hex(w.bytes);
}
std::string options_digest(const GapPicOptions &o) {
  auto x = o;
  Writer w;
  codec_options(w, x);
  return cn::sha256_hex(w.bytes);
}
// Logical snapshot identity. Excludes full_time_history (storage only).
const std::string &snapshot_digest(const ArenaSnapshot &s) {
  if (!s.digest_cache.empty())
    return s.digest_cache;
  require(s.journal.is_array() && s.journal.empty() && s.sinks.is_array() &&
              s.sinks.empty(),
          "full-time snapshot has no legacy JSON journal or sinks");
  Writer w;
  w.tag("native-full-time-snapshot-v5");
  w.text(packets_digest(s.packets));
  w.text(id_set_digest(s.ever_packet_ids));
  w.text(id_set_digest(s.event_ids));
  w.text(completions_digest(s.stochastic_completions));
  w.text(context_digest(s.full_time_context));
  w.text(s.chain);
  w.u64(s.revision);
  w.u64(s.journal_bytes);
  w.u64(s.accepted_steps);
  w.f64(s.previous_dt_s);
  // v6: births/sinks extend the identity only when they ever occurred, so a
  // history without them keeps every v5 snapshot digest bit for bit.
  if (s.full_time_flows) {
    w.text("native-full-time-flows-v6");
    w.text(flows_digest(*s.full_time_flows));
  }
  s.digest_cache = cn::sha256_hex(w.bytes);
  return s.digest_cache;
}
// ---------------------------------------------------------------------------
using NativeDigest = std::shared_ptr<const cn::BlockDigest>;
struct FullFrame {
  std::shared_ptr<const ArenaSnapshot> snapshot;
  PicArena::NativeState native;
  mutable NativeDigest digest;
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
  a.digest.swap(b.digest);
}
void release_native(PicArena::NativeState &n) noexcept {
  std::vector<GapParticle>().swap(n.particles);
  std::vector<double>().swap(n.horizontal);
  std::vector<double>().swap(n.vertical);
  mfem::Vector empty;
  n.potential.Swap(empty);
}
const cn::BlockDigest &native_digest(const PicArena::NativeState &n,
                                     NativeDigest &cache,
                                     const PicArena::NativeState *parent = nullptr,
                                     const NativeDigest *parent_digest = nullptr) {
  if (!cache)
    cache = cn::carrier_rows_digest(
        n.particles, parent ? &parent->particles : nullptr,
        parent_digest && *parent_digest ? parent_digest->get() : nullptr);
  return *cache;
}
const cn::BlockDigest &frame_digest(const FullFrame &f,
                                    const FullFrame *parent = nullptr) {
  return native_digest(f.native, f.digest, parent ? &parent->native : nullptr,
                       parent ? &parent->digest : nullptr);
}
// Native v5 accepted-state identity (replaces the v2 JSON frame prefix).
std::string state_identity(const std::string &initial,
                           const std::string &registry,
                           const std::string &options,
                           const ArenaSnapshot &s,
                           const PicArena::NativeState &n,
                           const cn::BlockDigest &rows) {
  Writer w;
  w.tag("native-full-time-exact-config-v5");
  w.text(initial);
  w.text(registry);
  w.text(options);
  w.text(snapshot_digest(s));
  w.text(rows.root);
  w.u64(n.next_id);
  w.text(ledger_digest(n.ledger));
  w.boolean(n.fields_current);
  w.u64(n.horizontal.size());
  for (double x : n.horizontal)
    w.f64(x);
  w.u64(n.vertical.size());
  for (double x : n.vertical)
    w.f64(x);
  w.u64(std::uint64_t(n.potential.Size()));
  for (int i = 0; i < n.potential.Size(); ++i)
    w.f64(n.potential[i]);
  return cn::sha256_hex(w.bytes);
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
void exposure_guard(const ExposureBinding &e, const FullTimeStepBinding &common,
                    FullTimeSlot slot, const std::string &field) {
  require(same_common(e.common, common) && e.slot == slot &&
              e.collision_duration_s == 0.5 * common.dt_s &&
              std::isfinite(e.collision_duration_s) &&
              e.collision_duration_s > 0 && e.frozen_field_sha256 == field &&
              e.frozen_physical_stage_s ==
                  (slot == FullTimeSlot::before_transport ? common.begin_s
                                                          : common.end_s),
          "full-time explicit half-exposure binding");
}
// Typed port of the candidate003 SSA child check over decoded reaction rows.
void full_ssa(const std::vector<FullRow> &children,
              const ArenaStochasticCompletion &c, const ExposureBinding &e) {
  require(c.stage_sequence == e.owner_stage_sequence &&
              c.common_clock_sha256 == e.common.common_clock_sha256 &&
              c.common_interval_begin_s == e.common.begin_s &&
              c.common_interval_end_s == e.common.end_s &&
              c.frozen_physical_stage_s == e.frozen_physical_stage_s,
          "full-time completion envelope");
  double previous = 0;
  for (const auto &row : children) {
    require(row.kind == RowKind::reaction, "full-time collision children only");
    const auto &a = row.reaction;
    require(a.phase == "full_time_endpoint" &&
                a.physical_stage_s == e.frozen_physical_stage_s &&
                a.clock_prefix == e.common.common_clock_sha256 &&
                a.mesh_prefix == e.common.mesh_sha256,
            "full-time reaction phase binding");
    if (a.selection_method == "forced_kernel_integration_only")
      continue;
    require(a.selection_method == "ssa_weighted_pair_clock_v1" &&
                a.channel == c.channel_id &&
                a.dataset_sha256 == c.dataset_sha256 &&
                a.kernel_manifest_sha256 == c.kernel_manifest_sha256 &&
                a.kernel_binding_file_sha256 == c.kernel_binding_file_sha256 &&
                a.conditional_model == c.model,
            "full-time SSA child channel binding");
    const double elapsed = a.collision_elapsed_s;
    require(std::isfinite(elapsed) && elapsed > previous &&
                elapsed <= e.collision_duration_s,
            "full-time SSA ordered explicit half-duration");
    previous = elapsed;
  }
}
std::optional<ArenaStochasticCompletion>
frame_completion(const ArenaSnapshot &s, const std::string &owner) {
  auto i = s.stochastic_completions.find(owner);
  if (i == s.stochastic_completions.end())
    return {};
  auto c = completion_read(*i);
  require(c.integration_owner_id == owner, "completion owner key");
  return c;
}
struct WindowSummary {
  std::uint64_t reactions = 0, motions = 0, steps = 0;
  std::uint64_t birth_rows = 0, births = 0, sink_rows = 0, sunk = 0;
  FullTimeAcceptedBoundary last_boundary;
  bool any_boundary = false;
};
// Dedicated full-time grammar over the unanchored window, starting from the
// exact accepted state committed by the last anchor (or genesis).
WindowSummary window_validate(const ArenaSnapshot &s,
                              const std::map<std::string, ArenaChannel> &channels,
                              std::size_t carrier_cap) {
  const auto &h = history(s);
  const auto &base = *h.base;
  auto boundary = base.boundary;
  require(boundary.accepted_step == base.steps &&
              boundary.accepted_arena_revision == base.revision,
          "full-time history base boundary/counters");
  std::string chain = base.chain, outer_sha, range_before_chain, motion_sha;
  std::set<std::string> events;
  std::uint64_t bytes = base.bytes, revision = base.revision, steps = base.steps;
  unsigned state = 0;
  double previous_dt = base.previous_dt_s;
  FullTimeStepBinding common;
  std::optional<ArenaStochasticCompletion> first, second;
  ArenaJson completions = base.completions;
  std::string field = boundary.field_sha256, heat = boundary.heat_sha256,
              source = boundary.source_sha256;
  std::vector<FullRow> children;
  WindowSummary summary;
  // v6: at most one birth row per step (after the first exposure, before the
  // motion, the only place the source seal may change) and at most one sink
  // row directly after a flows motion row.
  bool step_birth = false, after_motion = false;
  std::string motion_schema;
  std::uint64_t motion_carriers = 0;
  const auto rows = window_rows(h);
  for (const auto *stored : rows) {
    require(cn::sha256_hex(stored->bytes) == stored->sha256,
            "full-time journal row digest");
    const auto r = decode_row(stored->bytes);
    require(r.kind == stored->kind && r.previous_chain == chain,
            "full-time journal hash chain");
    const auto prior = chain;
    chain = stored->sha256;
    const auto index = revision;
    revision = checked_add(revision, 1, "full-time revision exhaustion");
    bytes = checked_add(bytes, stored->bytes.size(), "full-time byte total");
    if (state == 0) {
      require(r.kind == RowKind::outer_begin,
              "full-time journal requires outer begin");
      common = r.common;
      require(std::isfinite(common.dt_s) && common.dt_s > 0 &&
                  common.end_s > common.begin_s &&
                  common.end_s == common.begin_s + common.dt_s &&
                  common.begin_s == boundary.time_s &&
                  common.accepted_step == steps &&
                  common.accepted_arena_revision == index &&
                  common.accepted_field_sha256 == field &&
                  common.mesh_sha256 == boundary.mesh_sha256 &&
                  common.material_sha256 == boundary.material_sha256 &&
                  digest(common.accepted_arena_sha256) &&
                  digest(common.common_clock_sha256),
              "full-time journal common continuity");
      outer_sha = r.outer_sha;
      require(digest(outer_sha), "full-time outer commitment");
      step_birth = after_motion = false;
      state = 1;
      range_before_chain = chain;
      children.clear();
      first.reset();
      second.reset();
      continue;
    }
    if (state == 2 && r.kind == RowKind::birth) {
      const auto &bb = r.birth;
      std::uint64_t births = 0, previous_stream = 0;
      for (std::size_t i = 0; i < bb.streams.size(); ++i) {
        require(i == 0 || bb.streams[i].stream_id > previous_stream,
                "full-time journal birth streams order");
        previous_stream = bb.streams[i].stream_id;
        births = checked_add(births, bb.streams[i].births, "birth count");
      }
      require(!step_birth && bb.accepted_step == common.accepted_step &&
                  bb.begin_s == common.begin_s && bb.end_s == common.end_s &&
                  bb.source_before_sha256 == source && bb.base > 0 &&
                  births == bb.count &&
                  bb.count <= std::numeric_limits<std::uint64_t>::max() - bb.base &&
                  bb.next_id_after == bb.base + bb.count,
              "full-time journal birth binding");
      source = bb.source_after_sha256;
      step_birth = true;
      ++summary.birth_rows;
      summary.births += bb.count;
      continue;
    }
    if (state == 3 && r.kind == RowKind::sink) {
      const auto &k = r.sink;
      require(after_motion && motion_schema == motion_v3_flows_schema &&
                  k.transport_trial_sha256 == motion_sha &&
                  k.accepted_step == common.accepted_step &&
                  k.census == motion_carriers && k.carriers > 0 &&
                  k.packets > 0 && k.survivors + k.carriers == k.census,
              "full-time journal sink binding");
      after_motion = false;
      range_before_chain = chain; // second-slot range starts after the sink row
      ++summary.sink_rows;
      summary.sunk += k.carriers;
      continue;
    }
    after_motion = false;
    if ((state == 1 || state == 3) && r.kind == RowKind::reaction) {
      const auto &a = r.reaction;
      nonempty(a.event_id);
      require(events.insert(a.event_id).second &&
                  !base.event_ids.count(a.event_id),
              "duplicate full-time event");
      require(a.phase == "full_time_endpoint" && a.accepted_step == steps &&
                  a.physical_stage_s ==
                      (state == 1 ? common.begin_s : common.end_s) &&
                  a.previous_dt_s == previous_dt && a.next_dt_s == common.dt_s,
              "full-time journal reaction phase");
      require(!a.products.empty() && a.products.size() <= 8,
              "full-time journal product cap");
      require(digest(a.state_prefix) &&
                  a.selection_state_prefix == a.state_prefix,
              "full-time selection prefix");
      require(std::isfinite(a.event_weight) &&
                  std::isfinite(a.collision_elapsed_s),
              "full-time reaction operands");
      children.push_back(r);
      ++summary.reactions;
      continue;
    }
    if (state == 1 || state == 3) {
      require(r.kind == RowKind::exposure,
              "full-time journal requires scoped exposure");
      const auto &e = r.exposure;
      exposure_guard(e, common,
                     state == 1 ? FullTimeSlot::before_transport
                                : FullTimeSlot::after_transport,
                     field);
      const auto &c = r.completion;
      std::optional<ArenaStochasticCompletion> previous;
      auto found = completions.find(c.integration_owner_id);
      if (found != completions.end())
        previous = completion_read(*found);
      completion_validate(c, channels.at(c.channel_id), previous,
                          e.frozen_physical_stage_s);
      completion_chain_guard(e, c, boundary, state == 1 ? std::nullopt : first);
      full_ssa(children, c, e);
      require(r.stage_before_chain == range_before_chain &&
                  r.stage_after_reaction_chain == prior,
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
      require(r.kind == RowKind::motion, "full-time journal requires one motion");
      const auto &m = r.motion;
      require((m.record_schema == "canonical-full-time-total-field-motion-v2" ||
                   m.record_schema == motion_v3_schema ||
                   m.record_schema == motion_v3_flows_schema) &&
                  (!step_birth || m.record_schema == motion_v3_flows_schema) &&
                  m.kinetic_policy_schema ==
                      "native-finite-target-kinetic-policy-v1" &&
                  m.begin_s == common.begin_s && m.end_s == common.end_s &&
                  m.dt_s == common.dt_s &&
                  m.common_clock_sha256 == common.common_clock_sha256 &&
                  m.mesh_sha256 == common.mesh_sha256 &&
                  m.material_sha256 == common.material_sha256 &&
                  m.accepted_step == common.accepted_step &&
                  m.field_before_sha256 == field,
              "full-time journal motion envelope");
      require(m.carriers <= carrier_cap, "full-time motion endpoint cap");
      for (const auto *d :
           {&m.transport_trial_sha256, &m.census_sha256, &m.endpoints_root,
            &m.paths_root, &m.species_sha256, &m.trial_field_sha256,
            &m.kinetic_policy_sha256, &m.current_A_m_sha256,
            &m.initial_charge_C_sha256, &m.endpoint_charge_C_sha256,
            &m.staged_arena_sha256})
        require(digest(*d), "full-time journal motion digests");
      motion_sha = m.transport_trial_sha256;
      motion_schema = m.record_schema;
      motion_carriers = m.carriers;
      after_motion = true;
      field = m.field_after_sha256;
      heat = m.heat_after_sha256;
      require(digest(field) && digest(heat) && m.source_after_sha256 == source,
              "full-time journal closed owner seals");
      state = 3;
      range_before_chain = chain;
      children.clear();
      ++summary.motions;
      continue;
    }
    require(state == 4 && r.kind == RowKind::boundary &&
                same_common(r.common, common) && r.outer_sha == outer_sha &&
                r.transport_trial_sha256 == motion_sha,
            "full-time journal terminal boundary grammar");
    const auto &next = r.boundary;
    steps = checked_add(steps, 1, "full-time journal step exhaustion");
    require(first && second && next.accepted_step == steps &&
                next.accepted_arena_revision == revision &&
                next.time_s == common.end_s &&
                !next.authored_full_time_genesis &&
                next.field_sha256 == field && next.heat_sha256 == heat &&
                next.source_sha256 == source &&
                next.rng_sha256 == second->after_rng_clock_sha256 &&
                next.mesh_sha256 == boundary.mesh_sha256 &&
                next.material_sha256 == boundary.material_sha256 &&
                next.frame_sha256 == boundary.frame_sha256 &&
                next.common_clock_sha256 == common.common_clock_sha256 &&
                digest(r.first_prefix) && digest(r.second_prefix),
            "full-time journal final owner/clock/count seals");
    boundary = next;
    previous_dt = common.dt_s;
    summary.last_boundary = next;
    summary.any_boundary = true;
    ++summary.steps;
    state = 0;
  }
  require(state == 0, "full-time journal ends inside private step");
  require(bytes == s.journal_bytes && chain == s.chain &&
              revision == s.revision && steps == s.accepted_steps &&
              previous_dt == s.previous_dt_s &&
              completions == s.stochastic_completions &&
              s.event_ids.size() == base.event_ids.size() + events.size(),
          "full-time journal counters/completion replay mismatch");
  for (const auto &e : events)
    require(s.event_ids.count(e), "full-time window event identity");
  if (summary.any_boundary)
    require(same_boundary_payload(
                summary.last_boundary,
                boundary_payload(s.full_time_context.at("boundary"))),
            "full-time window terminal boundary differs from accepted context");
  return summary;
}
// ---------------------------------------------------------------------------
// Anchors.
template <class IO> void codec_anchor(IO &io, FullTimeAnchorFields &a) {
  io.tag("native-full-time-history-anchor-v6");
  codec_digest(io, a.source_sha256);
  codec_digest(io, a.history_policy_sha256);
  io.u64(a.anchor_interval_steps);
  io.u64(a.sequence);
  codec_digest(io, a.previous_anchor_sha256);
  io.u64(a.first_step);
  io.u64(a.last_step);
  io.u64(a.base_revision);
  io.u64(a.tail_revision);
  codec_digest(io, a.base_history_chain);
  codec_digest(io, a.tail_history_chain);
  io.u64(a.window_rows);
  io.u64(a.window_bytes);
  io.u64(a.cumulative_bytes);
  io.u64(a.window_reactions);
  io.u64(a.window_motions);
  io.u64(a.window_birth_rows);
  io.u64(a.window_births);
  io.u64(a.window_sink_rows);
  io.u64(a.window_sunk_carriers);
  codec_boundary(io, a.boundary, true);
  for (auto *s : {&a.snapshot_sha256, &a.native_rows_sha256,
                  &a.native_ledger_sha256, &a.completions_sha256,
                  &a.event_ids_sha256, &a.packet_ids_sha256})
    codec_digest(io, *s);
  io.u64(a.next_native_id);
  io.f64(a.previous_dt_s);
  io.text(a.flows_sha256, 64);
  if constexpr (reading<IO>) {
    require(a.flows_sha256.empty() || cn::is_digest(a.flows_sha256),
            "anchor flows digest");
    a.accepted_arena_sha256 = a.boundary.accepted_arena_sha256;
  }
}
std::string policy_digest(const FullTimeHistoryPolicy &p) {
  Writer w;
  w.tag("native-full-time-history-policy-v5");
  w.u64(p.anchor_interval_steps);
  w.u64(p.maximum_anchors);
  w.u64(p.maximum_anchor_bytes);
  w.u64(p.maximum_window_rows);
  w.u64(p.maximum_window_bytes);
  w.u64(cn::carrier_block);
  return cn::sha256_hex(w.bytes);
}
std::string anchor_root(const std::string &policy,
                        const FullTimeAcceptedBoundary &genesis,
                        const std::string &genesis_chain) {
  auto b = genesis;
  Writer w;
  w.tag("native-full-time-anchor-root-v6");
  w.text(PicArena::source_identity());
  w.text(policy);
  codec_boundary(w, b, true);
  w.text(genesis_chain);
  return cn::sha256_hex(w.bytes);
}
constexpr std::uint64_t anchor_byte_bound = 4096;
} // namespace full_time_detail
FullTimeAnchorFields decode_full_time_anchor(const std::string &canonical) {
  using namespace full_time_detail;
  require(canonical.size() <= anchor_byte_bound, "bounded anchor record");
  Reader in(canonical);
  FullTimeAnchorFields a;
  codec_anchor(in, a);
  in.finish();
  a.sha256 = cn::sha256_hex(canonical);
  return a;
}
std::vector<FullTimeAnchorFields>
verify_full_time_anchor_chain(const std::string &serialized,
                              const std::string &expected_root,
                              const std::string &expected_head) {
  using namespace full_time_detail;
  Reader in(serialized);
  in.tag("native-full-time-anchor-chain-v6");
  std::string policy, genesis_chain;
  FullTimeHistoryPolicy p;
  codec_digest(in, policy);
  in.u64(p.anchor_interval_steps);
  std::uint64_t count = 0;
  in.u64(count);
  FullTimeAcceptedBoundary genesis;
  codec_boundary(in, genesis, true);
  codec_digest(in, genesis_chain);
  require(count <= (1ull << 32), "bounded anchor count");
  require(anchor_root(policy, genesis, genesis_chain) == expected_root,
          "anchor chain root differs from the independently held root");
  require(genesis.accepted_step == 0 && genesis.accepted_arena_revision == 0 &&
              genesis.authored_full_time_genesis,
          "anchor chain genesis boundary");
  std::vector<FullTimeAnchorFields> out;
  out.reserve(std::size_t(count));
  std::string previous = expected_root;
  FullTimeAnchorFields prior;
  prior.tail_history_chain = genesis_chain;
  prior.boundary = genesis;
  for (std::uint64_t i = 0; i < count; ++i) {
    std::string bytes;
    in.text(bytes, anchor_byte_bound);
    auto a = decode_full_time_anchor(bytes);
    const auto K = a.anchor_interval_steps;
    require(a.source_sha256 == PicArena::source_identity() &&
                a.history_policy_sha256 == policy &&
                K == p.anchor_interval_steps && K > 0 && a.sequence == i &&
                a.previous_anchor_sha256 == previous,
            "anchor identity/policy/sequence/previous link");
    require(a.first_step == prior.boundary.accepted_step &&
                a.last_step == a.first_step + K &&
                a.base_revision == prior.boundary.accepted_arena_revision &&
                a.tail_revision == a.base_revision + a.window_rows &&
                a.base_history_chain == prior.tail_history_chain &&
                a.cumulative_bytes ==
                    (i ? prior.cumulative_bytes : 0) + a.window_bytes &&
                a.window_motions == K && a.window_birth_rows <= K &&
                a.window_sink_rows <= K &&
                (a.window_birth_rows > 0 || a.window_births == 0) &&
                (a.window_sink_rows > 0 || a.window_sunk_carriers == 0),
            "anchor step/revision/chain/byte continuity");
    const auto &b = a.boundary;
    require(b.accepted_step == a.last_step &&
                b.accepted_arena_revision == a.tail_revision &&
                std::isfinite(b.time_s) && b.time_s > prior.boundary.time_s &&
                !b.authored_full_time_genesis &&
                b.mesh_sha256 == genesis.mesh_sha256 &&
                b.material_sha256 == genesis.material_sha256 &&
                b.frame_sha256 == genesis.frame_sha256 &&
                digest(b.source_sha256) &&
                (a.window_birth_rows > 0 ||
                 b.source_sha256 == prior.boundary.source_sha256) &&
                digest(b.accepted_arena_sha256) && digest(b.rng_sha256) &&
                digest(b.field_sha256) && digest(b.heat_sha256) &&
                digest(b.common_clock_sha256) && std::isfinite(a.previous_dt_s) &&
                a.previous_dt_s > 0,
            "anchor accepted boundary continuity");
    previous = a.sha256;
    prior = a;
    out.push_back(std::move(a));
  }
  in.finish();
  require(previous == expected_head,
          "anchor chain head differs from the independently held head");
  return out;
}
namespace full_time_detail {
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
  // v6: birth scope reservation (B_s) and the accepted motion's flows.
  bool births_reserved = false;
  std::uint64_t birth_base = 0;
  std::shared_ptr<const FullTimeStepFlows> flows;
};
// v6 births attached to a transport trial.
struct BirthInfo {
  std::shared_ptr<const FullTimeBirthReceipt> receipt;
  std::vector<FullTimeBirthFlight> flights; // ascending native ID
  std::string source_after;
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
  std::size_t limit = 0;
  std::uint64_t base_revision = 0;
};
struct FullTrial {
  std::shared_ptr<FullOuter> outer;
  FullFrame base;
  std::string prefix;
  bool consumed = false;
  std::shared_ptr<const BirthInfo> births; // v6 (null: no birth credit)
  // Outer frame the trial was forked from (= base.snapshot without births).
  std::shared_ptr<const ArenaSnapshot> fork_frame;
  // v6 per-trial cache (computed at most once): census digest, the logical
  // packet of every base row, and multi-carrier packets (row indices).
  mutable bool census_cached = false;
  mutable std::string census;
  mutable std::vector<const ArenaPacket *> owner_of_row;
  mutable std::vector<std::vector<std::size_t>> multi;
};
// Hall v2 census encoding (identical to hall_driver::canonical_census_sha256).
std::string census_sha(const std::vector<ArenaRateEntry> &census) {
  std::vector<const ArenaRateEntry *> rows;
  rows.reserve(census.size());
  for (const auto &r : census)
    rows.push_back(&r);
  std::sort(rows.begin(), rows.end(),
            [](auto a, auto b) { return a->native.id < b->native.id; });
  cn::Hasher h;
  h.w.tag("canonical-motion-census-v1");
  h.w.u64(rows.size());
  for (auto r : rows) {
    h.w.text(r->packet_id);
    h.w.text(r->state);
    const auto &p = r->native;
    h.w.u64(p.id);
    h.w.vec(p.position_m);
    h.w.vec(p.proper_velocity_m_s);
    h.w.f64(p.charge_C);
    h.w.f64(p.mass_kg);
    h.w.f64(p.weight);
    h.flush_if_large();
  }
  return cn::hex(h.finish());
}
std::string kinetic_policy_hash(const GapPicOptions &o) {
  std::string bytes = "native-finite-target-kinetic-policy-v1";
  bytes.push_back('\0');
  Writer e;
  for (double x :
       {o.width_m, o.gap_m, o.short_z_m, o.feed_z_m, o.lower_foil.left_m,
        o.lower_foil.right_m, o.lower_foil.bottom_m, o.lower_foil.top_m,
        o.upper_foil.left_m, o.upper_foil.right_m, o.upper_foil.bottom_m,
        o.upper_foil.top_m})
    e.f64(x);
  return cn::sha256_hex(bytes + e.bytes);
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
struct PathRecord {
  std::uint64_t id = 0;
  PicVector start{}, displacement{};
  double flight_s = 0;
};
// Block Merkle root over Hall path records [u64 id, vec start, vec D, f64 flight].
std::string path_rows_root(const std::vector<PathRecord> &rows) {
  std::vector<cn::BlockLeaf> leaves;
  std::size_t i = 0;
  while (i < rows.size()) {
    require(rows[i].id > 0, "path record ID");
    const auto block = (rows[i].id - 1) / cn::carrier_block;
    std::size_t j = i;
    Writer w;
    w.tag("native-v5/path-block");
    w.u64(block);
    while (j < rows.size() && (rows[j].id - 1) / cn::carrier_block == block)
      ++j;
    w.u64(j - i);
    for (std::size_t k = i; k < j; ++k) {
      require(k == i || rows[k].id > rows[k - 1].id, "path record order");
      w.u64(rows[k].id);
      w.vec(rows[k].start);
      w.vec(rows[k].displacement);
      w.f64(rows[k].flight_s);
    }
    leaves.push_back({block, i, j - i, cn::sha256_raw(w.bytes)});
    i = j;
  }
  return cn::merkle_root("native-v5/path-root", leaves, rows.size());
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
struct BirthScope::Impl {
  std::shared_ptr<FullOuter> outer;
  std::shared_ptr<const ArenaSnapshot> frame;
  std::uint64_t base = 0;
  std::string source_before;
};
struct PreparedBirths::Impl {
  std::shared_ptr<FullOuter> outer;
  std::shared_ptr<const ArenaSnapshot> base;
  FullFrame after;
  std::shared_ptr<const BirthInfo> info;
};
struct PreparedPrivateMotion::Impl {
  ~Impl() noexcept {
    drain_json(envelope);
    drain_json(audit);
  }
  std::shared_ptr<FullTrial> trial;
  FullFrame after;
  // v6 binary admission: record summary and sink ledger of this motion.
  std::shared_ptr<const FullTimeSinkLedger> sinks;
  std::string v3_schema, v3_sha256;
  std::uint64_t v3_endpoints = 0, v3_paths = 0;
  // The complete Hall record for this step stays with the caller's token; the
  // journal retains only its digests.
  ArenaJson envelope;
  mutable ArenaJson audit;
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
  static std::string prefix(const PicArena &a, const ArenaSnapshot &s,
                            const PicArena::NativeState &n,
                            const cn::BlockDigest &rows) {
    return state_identity(a.initial_prefix_, a.registry_v5_sha_,
                          a.options_v5_sha_, s, n, rows);
  }
  static std::string prefix(const FullOuter &o, const FullFrame &f) {
    return prefix(*o.owner, *f.snapshot, f.native, frame_digest(f));
  }
  static std::vector<ArenaRateEntry>
  population(const FullOuter &o, const FullFrame &f, std::size_t cap) {
    return o.owner->full_time_census(*f.snapshot, f.native.particles,
                                     std::min(cap, o.cap), true);
  }
  static std::uint64_t window_room(const ArenaSnapshot &s) {
    const auto &h = history(s);
    return h.policy->policy.maximum_window_rows - h.window_rows();
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
    const auto &limits = o.owner->identity_cap_;
    if (o.frame.snapshot->ever_packet_ids.size() + e.products.size() + 1 >
            limits ||
        o.frame.snapshot->event_ids.size() + 1 > limits)
      throw refusal(RK::identity_count,
                    "full-time logical identity cap before allocation");
    auto built = o.owner->build_reaction(e, *o.frame.snapshot, o.frame.native,
                                         prefix(o, o.frame), true);
    require(!fault, "manufactured private reaction fault");
    auto token = std::make_unique<PreparedArena::Impl>();
    token->private_identity = s.token_identity;
    token->private_generation = s.generation;
    token->owner = o.owner;
    token->before = o.frame.snapshot;
    token->after = built.snapshot;
    token->native_after = std::move(built.native);
    token->native_after_digest = cn::carrier_rows_digest(
        token->native_after.particles, &o.frame.native.particles,
        frame_digest(o.frame).root.empty() ? nullptr : o.frame.digest.get());
    token->audit = std::move(built.audit);
    token->candidate_prefix = prefix(*o.owner, *token->after,
                                     token->native_after,
                                     *token->native_after_digest);
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
    const bool final = t.impl_->full_time_exposure;
    if ((!s.building && !final) || (s.building && final))
      return false;
    // All validation, hashing and allocation precede the first private swap.
    std::optional<FullTimeExposureReceipt> receipt;
    std::string completed;
    if (final) {
      receipt.emplace(FullTimeExposureReceipt{
          s.exposure, t.impl_->completion, t.impl_->candidate_prefix,
          completion_sha(t.impl_->completion), s.key});
      completed = t.impl_->candidate_prefix;
    }
    FullFrame prepared{t.impl_->after, {}, t.impl_->native_after_digest};
    static_assert(std::is_nothrow_move_assignable_v<FullTimeExposureReceipt>);
    std::swap(prepared.native.options, t.impl_->native_after.options);
    prepared.native.particles.swap(t.impl_->native_after.particles);
    prepared.native.horizontal.swap(t.impl_->native_after.horizontal);
    prepared.native.vertical.swap(t.impl_->native_after.vertical);
    prepared.native.potential.Swap(t.impl_->native_after.potential);
    std::swap(prepared.native.ledger, t.impl_->native_after.ledger);
    std::swap(prepared.native.next_id, t.impl_->native_after.next_id);
    std::swap(prepared.native.fields_current,
              t.impl_->native_after.fields_current);
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
    // Consumed private tokens retain audit/prefix only, never a population.
    release_native(prepared.native);
    release_native(t.impl_->native_before);
    release_native(t.impl_->native_after);
    release_native(t.impl_->work);
    prepared.snapshot.reset();
    prepared.digest.reset();
    t.impl_->native_after_digest.reset();
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
    const std::uint64_t used =
        o.frame.snapshot->revision - o.before.snapshot->revision;
    const bool first_slot = s.exposure.slot == FullTimeSlot::before_transport;
    const auto room = window_room(*o.frame.snapshot);
    if (first_slot) {
      require(used <= o.records && o.records - used >= 4 &&
                  maximum <= (o.records - used - 4) / 2,
              "reserve both slot reaction budgets before callbacks");
      if (room < 4 || maximum > (room - 4) / 2)
        throw refusal(RK::window_rows,
                      "reserve both slot history window budgets before "
                      "callbacks");
    } else {
      require(maximum <= o.slot_reaction_budget && used <= o.records &&
                  o.records - used >= 2 && maximum <= o.records - used - 2,
              "second slot exceeds reserved reaction/history budget");
      if (room < 2 || maximum > room - 2)
        throw refusal(RK::window_rows,
                      "second slot exceeds reserved history window budget");
    }
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
                          frame_completion(*base.snapshot, c.integration_owner_id),
                          o.frame.native.ledger.time_s);
      full_ssa(rows_since(*o.frame.snapshot, base.snapshot->revision), c,
               s.exposure);
      require(s.outstanding_token.expired(),
              "unconsumed private event token at completion");
      auto next = std::make_shared<ArenaSnapshot>(*o.frame.snapshot);
      next->stochastic_completions[c.integration_owner_id] = completion_json(c);
      FullRow row;
      row.kind = RowKind::exposure;
      row.exposure = s.exposure;
      row.completion = c;
      row.stage_before_chain = base.snapshot->chain;
      row.stage_after_reaction_chain = o.frame.snapshot->chain;
      full_time_append(*next, row);
      require(next->revision - o.before.snapshot->revision <= o.records,
              "outer record resource cap");
      ArenaJson a = {{"kind", "full_time_exposure"},
                     {"exposure", exposure_json(s.exposure)},
                     {"completion", completion_json(c)},
                     {"stage_before_chain", base.snapshot->chain},
                     {"stage_after_reaction_chain", o.frame.snapshot->chain},
                     {"row_sha256", next->chain}};
      auto t = std::make_unique<PreparedArena::Impl>();
      t->owner = o.owner;
      t->private_identity = s.token_identity;
      t->private_generation = s.generation + 1;
      t->before = base.snapshot;
      t->after = next;
      t->native_after = o.frame.native;
      t->native_after_digest = o.frame.digest
                                   ? o.frame.digest
                                   : cn::carrier_rows_digest(
                                         o.frame.native.particles,
                                         &base.native.particles,
                                         base.digest.get());
      t->full_time_exposure = true;
      t->completion = c;
      t->audit = std::move(a);
      t->candidate_prefix =
          prefix(*o.owner, *next, t->native_after, *t->native_after_digest);
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

  // ---- v6 helpers --------------------------------------------------------
  static const std::string &trial_source(const FullTrial &t) {
    return t.births ? t.births->source_after : t.outer->boundary->source_sha256;
  }
  // Census digest, packet of every base row and multi-carrier packets of a
  // trial, computed at most once (no per-carrier allocation).
  static void census_cache(const FullTrial &t) {
    if (t.census_cached)
      return;
    const auto &rows = t.base.native.particles;
    const auto &s = *t.base.snapshot;
    std::vector<const ArenaPacket *> owner(rows.size(), nullptr);
    std::vector<std::vector<std::size_t>> multi;
    for (const auto &[id, p] : s.packets) {
      (void)id;
      if (p.native_ids.size() > 1)
        multi.emplace_back();
      for (auto native_id : p.native_ids) {
        const auto it = std::lower_bound(
            rows.begin(), rows.end(), native_id,
            [](const GapParticle &a, std::uint64_t b) { return a.id < b; });
        require(it != rows.end() && it->id == native_id,
                "native canonical carrier missing");
        const auto k = std::size_t(it - rows.begin());
        require(!owner[k], "native carrier owned twice");
        owner[k] = &p;
        if (p.native_ids.size() > 1)
          multi.back().push_back(k);
      }
    }
    cn::Hasher h;
    h.w.tag("canonical-motion-census-v1");
    h.w.u64(rows.size());
    for (std::size_t k = 0; k < rows.size(); ++k) {
      require(owner[k] != nullptr, "unowned native particle");
      h.w.text(owner[k]->id);
      h.w.text(owner[k]->state);
      h.w.raw(&rows[k], sizeof(GapParticle));
      h.flush_if_large();
    }
    t.census = cn::hex(h.finish());
    t.owner_of_row.swap(owner);
    t.multi.swap(multi);
    t.census_cached = true;
  }
  static std::size_t rows_used(const FullOuter &o, const ArenaSnapshot &s) {
    return std::size_t(s.revision - o.before.snapshot->revision);
  }
  // Births: validates the credit against contract (a)-(d) and builds the
  // post-birth frame (pure; the outer frame is not changed).
  static std::unique_ptr<PreparedBirths::Impl>
  births(const BirthScope::Impl &scope, const FullTimeBirthCredit &c) {
    auto &o = *scope.outer;
    live(o);
    flow_need(o.phase == 1 && !o.builder_active &&
                  o.frame.snapshot == scope.frame &&
                  o.frame.native.next_id == scope.base,
              FK::birth_scope,
              "birth scope is stale (outer frame changed or wrong phase)");
    const auto &base = o.frame;
    const auto &opt = base.native.options;
    const std::uint64_t B = scope.base;
    flow_need(c.arrival_step == o.common.accepted_step, FK::birth_binding,
              "arrival step differs from the outer accepted step");
    flow_need(c.base_native_id == B, FK::birth_binding,
              "birth base differs from native B_s of this step");
    flow_need(digest(c.plan_sha256) && digest(c.staged_token_sha256) &&
                  digest(c.source_before_sha256) &&
                  digest(c.source_after_sha256),
              FK::birth_binding, "birth plan/token/source digests required");
    flow_need(c.source_before_sha256 == scope.source_before,
              FK::birth_binding,
              "birth source_before differs from the accepted source seal");
    const std::size_t n = c.births.size();
    flow_need(n <= 8192 && c.streams.size() <= 8192, FK::birth_order,
              "bounded births/streams per step");
    birth_id_block_guard(B, n);
    if (n > o.cap - base.native.particles.size())
      throw refusal(RK::carrier_cap,
                    "births exceed the outer carrier budget before allocation");
    // Rows still needed after the birth row: motion, optional sink row,
    // second-slot reactions + exposure, boundary.
    const auto used = rows_used(o, *base.snapshot);
    const std::size_t need = 4 + o.slot_reaction_budget + 1;
    if (used > o.records || o.records - used < need)
      throw refusal(RK::step_records,
                    "birth credit exceeds the declared outer record budget");
    if (window_room(*base.snapshot) < need)
      throw refusal(RK::window_rows,
                    "birth credit exceeds the history window budget");
    // Streams: ascending, registered species.
    std::map<std::uint32_t, std::size_t> stream_index;
    for (std::size_t i = 0; i < c.streams.size(); ++i) {
      const auto &x = c.streams[i];
      flow_need(i == 0 || x.stream_id > c.streams[i - 1].stream_id,
                FK::birth_order, "birth streams must ascend by stream id");
      flow_need(o.owner->states_.count(x.state) == 1, FK::birth_binding,
                "birth stream species is not registered");
      stream_index.emplace(x.stream_id, i);
    }
    std::vector<ExactSum> credit(c.streams.size());
    std::vector<std::uint64_t> per_stream(c.streams.size(), 0);
    auto info = std::make_shared<BirthInfo>();
    auto receipt = std::make_shared<FullTimeBirthReceipt>();
    info->flights.reserve(n);
    receipt->entries.reserve(n);
    std::vector<GapParticle> born;
    born.reserve(n);
    for (std::size_t r = 0; r < n; ++r) {
      const auto &b = c.births[r];
      if (r) {
        const auto &a = c.births[r - 1];
        flow_need(b.stream_id > a.stream_id ||
                      (b.stream_id == a.stream_id && b.index > a.index),
                  FK::birth_order,
                  "births not in canonical (stream id, index) order");
      }
      const auto si = stream_index.find(b.stream_id);
      flow_need(si != stream_index.end(), FK::birth_binding,
                "birth of an undeclared stream");
      flow_need(b.index == per_stream[si->second], FK::birth_order,
                "birth indices of a stream must be 0..n-1 in order");
      ++per_stream[si->second];
      flow_need(b.step == c.arrival_step, FK::birth_binding,
                "birth step differs from the arrival step");
      flow_need(b.state == c.streams[si->second].state, FK::birth_binding,
                "birth species differs from its stream");
      flow_need(b.packet_id ==
                    full_time_birth_packet_id(b.step, b.stream_id, b.index),
                FK::birth_identity,
                "birth packet ID is not the canonical arrival-v1 ID of its key");
      flow_need(!base.snapshot->packets.count(b.packet_id) &&
                    !base.snapshot->ever_packet_ids.count(b.packet_id),
                FK::birth_identity, "duplicate birth packet ID");
      flow_need(b.native_id == B + r, FK::birth_identity,
                "birth native ID differs from B_s + canonical rank");
      flow_need(!b.origin.empty() && b.origin.size() < 1024,
                FK::birth_identity, "birth origin (source event) required");
      flow_need(std::isfinite(b.weight) && b.weight > 0, FK::birth_count,
                "birth weight must be one positive finite carrier");
      const auto &spec = o.owner->states_.at(b.state);
      GapParticle p;
      p.id = b.native_id;
      p.position_m = b.position_m;
      p.proper_velocity_m_s = b.proper_velocity_m_s;
      p.charge_C = spec.charge_e * pic_e;
      p.mass_kg = spec.mass_kg;
      p.weight = b.weight;
      // Strictly interior vacuum point, or on an emitting foil face with an
      // outward proper velocity (source v5 neutral packets).
      bool interior = true;
      std::string why;
      try {
        closed_position(opt, p);
      } catch (const std::exception &e) {
        interior = false;
        why = e.what();
      }
      if (!interior && !surface_birth(opt, p))
        flow_refusal(FK::birth_kinematics,
                     "birth position/velocity: " + why +
                         " (and not a foil-face birth leaving into vacuum)");
      flow_need(std::isfinite(b.birth_time_s) &&
                    b.birth_time_s >= o.common.begin_s &&
                    b.birth_time_s <= o.common.end_s &&
                    b.remaining_flight_s == o.common.end_s - b.birth_time_s &&
                    b.remaining_flight_s >= 0,
                FK::birth_kinematics,
                "birth time outside the step or remaining flight != end - tau");
      flow_need(std::isfinite(b.weight * pic_kinetic_energy_J(p)),
                FK::birth_kinematics, "birth kinetic energy operand");
      credit[si->second].add(b.weight);
      born.push_back(p);
      info->flights.push_back({p.id, b.birth_time_s, b.remaining_flight_s});
      receipt->entries.push_back({b.packet_id, b.state, p.id, b.weight});
    }
    receipt->streams = c.streams;
    for (std::size_t i = 0; i < c.streams.size(); ++i)
      flow_need(credit[i].count() == c.streams[i].source_debit,
                FK::birth_count,
                "native credit differs from the exact source debit of a stream");
    // Post-birth frame: births are the largest native IDs (B_s ... B_s+n-1).
    auto after = std::make_unique<PreparedBirths::Impl>();
    after->outer = scope.outer;
    after->base = base.snapshot;
    auto &native = after->after.native;
    native = base.native;
    native.particles.reserve(native.particles.size() + n);
    for (const auto &p : born)
      native.particles.push_back(p);
    native.next_id = B + n;
    auto next = std::make_shared<ArenaSnapshot>(*base.snapshot);
    for (std::size_t r = 0; r < n; ++r) {
      const auto &b = c.births[r];
      ArenaPacket packet{b.packet_id, b.state, b.origin, "",
                         ArenaCount(b.weight), {b.native_id}};
      require(next->packets.emplace(b.packet_id, packet).second,
              "birth logical insertion collision");
    }
    auto flows = next->full_time_flows
                     ? std::make_shared<FullTimeFlows>(*next->full_time_flows)
                     : std::make_shared<FullTimeFlows>();
    flows->digest.clear();
    // Last nonempty birth block (the arrival owner's anchor semantics).
    if (n) {
      flows->has_last_birth = true;
      flows->last_birth_step = c.arrival_step;
      flows->last_birth_base = B;
      flows->last_birth_count = n;
      flows->last_birth_plan_sha256 = c.plan_sha256;
      flows->last_birth_token_sha256 = c.staged_token_sha256;
    }
    flows->births = checked_add(flows->births, n, "birth total");
    flows->birth_rows = checked_add(flows->birth_rows, 1, "birth rows");
    next->full_time_flows = flows;
    next->full_time_context["source_sha256"] = c.source_after_sha256;
    std::vector<FullTimeBirth> sorted(c.births.begin(), c.births.end());
    FullRow row;
    row.kind = RowKind::birth;
    auto &br = row.birth;
    br.accepted_step = c.arrival_step;
    br.base = B;
    br.count = n;
    br.next_id_after = B + n;
    br.begin_s = o.common.begin_s;
    br.end_s = o.common.end_s;
    br.plan_sha256 = c.plan_sha256;
    br.staged_token_sha256 = c.staged_token_sha256;
    br.source_before_sha256 = c.source_before_sha256;
    br.source_after_sha256 = c.source_after_sha256;
    br.births_root = id_records_root("native-v6/birth-block",
                                     "native-v6/birth-root", sorted,
                                     encode_birth_record);
    for (std::size_t i = 0; i < c.streams.size(); ++i)
      br.streams.push_back({c.streams[i].stream_id, per_stream[i],
                            c.streams[i].state, credit[i].count()});
    receipt->accepted_step = c.arrival_step;
    receipt->base = B;
    receipt->next_id_after = B + n;
    receipt->plan_sha256 = c.plan_sha256;
    receipt->staged_token_sha256 = c.staged_token_sha256;
    receipt->source_before_sha256 = c.source_before_sha256;
    receipt->source_after_sha256 = c.source_after_sha256;
    receipt->births_root_sha256 = br.births_root;
    receipt->row_sha256 = full_time_append(*next, std::move(row));
    (void)o.owner->full_time_census(*next, native.particles, o.cap, false);
    after->after.snapshot = next;
    (void)native_digest(native, after->after.digest, &base.native, &base.digest);
    info->receipt = receipt;
    info->source_after = c.source_after_sha256;
    after->info = info;
    return after;
  }
  // Binary motion v3 admission (plain Hall corrective v2 record or the v6
  // flows record with births and first-hit sinks). O(carriers) over the
  // binary arrays; one allocation for the endpoint rows.
  static std::unique_ptr<PreparedPrivateMotion::Impl>
  motion_v3(const std::shared_ptr<FullTrial> &t, const FullTimeMotionBinding &b,
            const FullTimeMotionRecordV3 &r) {
    auto &o = *t->outer;
    live(o);
    require(o.phase == 1 && !t->consumed,
            "closed motion source/boundary ownership unsupported");
    require(same_common(b.common, o.common) &&
                b.field_before_sha256 == o.common.accepted_field_sha256 &&
                digest(b.transport_trial_sha256) &&
                digest(b.field_after_sha256) && digest(b.heat_after_sha256) &&
                b.source_after_sha256 == trial_source(*t),
            "motion common/owner seals");
    const bool flows = r.schema == motion_v3_flows_schema;
    flow_need(flows || r.schema == motion_v3_schema, FK::motion_record,
              "unrecognized binary motion record schema");
    flow_need(r.coverage == (flows ? flows_coverage : plain_coverage) &&
                  r.neutral_motion_owner == "same_full_time_q0_drift" &&
                  !r.native_DKD_called && !r.native_Poisson_called,
              FK::motion_record, "binary motion coverage/mode");
    flow_need(flows || (r.sinks.empty() && r.births.empty()),
              FK::motion_record, "plain v3 record carries no sinks/births");
    flow_need(!t->births || flows, FK::motion_record,
              "births require the binary flows motion record");
    const auto &opt = t->base.native.options;
    require(opt.finite_foils && opt.finite_axial_extent,
            "Hall v2 finite kinetic geometry required");
    const auto &base = t->base.native.particles;
    const std::size_t N = base.size();
    flow_need(N <= o.cap && r.paths.size() == N &&
                  r.endpoints.size() <= N &&
                  r.sinks.size() == N - r.endpoints.size(),
              FK::motion_record,
              "binary motion endpoint/sink/path dimensions");
    const auto record_sha = full_time_motion_record_v3_sha256(r);
    flow_need(record_sha == r.record_sha256 &&
                  record_sha == b.transport_trial_sha256,
              FK::motion_record,
              "binary motion record digest differs from its seal/binding");
    flow_need(r.accepted_step == o.common.accepted_step &&
                  r.staged_arena_revision == t->base.snapshot->revision &&
                  r.staged_arena_sha256 == t->prefix &&
                  r.begin_s == o.common.begin_s && r.end_s == o.common.end_s &&
                  r.dt_s == o.common.dt_s &&
                  r.mesh_sha256 == o.common.mesh_sha256 &&
                  r.material_sha256 == o.common.material_sha256 &&
                  r.common_clock_sha256 == o.common.common_clock_sha256 &&
                  digest(r.trial_field_sha256) && digest(r.species_sha256) &&
                  r.kinetic_policy_sha256 == kinetic_policy_hash(opt),
              FK::motion_record, "motion owner staged identity");
    const std::array<double, 12> policy = {
        opt.width_m,           opt.gap_m,             opt.short_z_m,
        opt.feed_z_m,          opt.lower_foil.left_m, opt.lower_foil.right_m,
        opt.lower_foil.bottom_m, opt.lower_foil.top_m, opt.upper_foil.left_m,
        opt.upper_foil.right_m, opt.upper_foil.bottom_m, opt.upper_foil.top_m};
    flow_need(std::memcmp(policy.data(), r.kinetic_policy.data(),
                          sizeof policy) == 0,
              FK::motion_record, "Hall complete native kinetic policy operands");
    const PicVector shift = {-.5 * opt.width_m, -.5 * opt.gap_m, 0};
    flow_need(r.native_to_field_m == shift, FK::motion_record,
              "Hall translation-only frame");
    flow_need(digest(r.current_A_m_sha256) &&
                  digest(r.initial_charge_C_sha256) &&
                  digest(r.endpoint_charge_C_sha256),
              FK::motion_record, "Hall charge/current array identity");
    census_cache(*t);
    flow_need(r.census_sha256 == t->census, FK::motion_record,
              "binary motion census digest differs from the trial census");
    // Births stated by the record must equal the trial's newborns.
    static const std::vector<FullTimeBirthFlight> no_births;
    const auto &newborn = t->births ? t->births->flights : no_births;
    flow_need(r.births.size() == newborn.size(), FK::motion_record,
              "binary motion births differ from the trial births");
    for (std::size_t k = 0; k < newborn.size(); ++k)
      flow_need(r.births[k].native_id == newborn[k].native_id &&
                    r.births[k].birth_time_s == newborn[k].birth_time_s &&
                    r.births[k].flight_s == newborn[k].flight_s,
                FK::motion_record,
                "binary motion births differ from the trial births");
    auto result = std::make_unique<PreparedPrivateMotion::Impl>();
    auto &after = result->after.native;
    after.particles.reserve(r.endpoints.size());
    std::size_t j = 0, s = 0, nb = 0;
    long double kinetic_change = 0;
    std::vector<std::size_t> sunk_rows;
    sunk_rows.reserve(r.sinks.size());
    for (std::size_t i = 0; i < N; ++i) {
      const auto &old = base[i];
      const auto &path = r.paths[i];
      flow_need(path.native_id == old.id, FK::motion_record,
                "Hall path rows must cover the census once in ascending ID");
      double F = o.common.dt_s;
      const bool is_newborn =
          nb < newborn.size() && newborn[nb].native_id == old.id;
      if (is_newborn)
        F = newborn[nb++].flight_s;
      for (int d = 0; d < 3; ++d)
        flow_need(path.field_start_m[d] == old.position_m[d] + shift[d],
                  FK::motion_record, "Hall canonical path start operand");
      if (s < r.sinks.size() && r.sinks[s].native_id == old.id) {
        const auto &k = r.sinks[s++];
        sink_check(opt, old, path, k, F, shift);
        sunk_rows.push_back(i);
        GapParticle hit = old;
        hit.proper_velocity_m_s = k.incident_proper_velocity_m_s;
        kinetic_change += static_cast<long double>(old.weight) *
                          (pic_kinetic_energy_J(hit) - pic_kinetic_energy_J(old));
        continue;
      }
      flow_need(j < r.endpoints.size() && r.endpoints[j].id == old.id,
                FK::sink_duplicate,
                "census carrier without exactly one endpoint or sink record");
      const auto &e = r.endpoints[j++];
      closed_position(opt, e);
      require(e.charge_C == old.charge_C && e.mass_kg == old.mass_kg &&
                  e.weight == old.weight,
              "motion altered native species/count");
      flow_need(path.flight_s == F, FK::motion_record,
                F == o.common.dt_s ? "Hall complete common flight"
                                   : "newborn flight differs from end - tau");
      for (int d = 0; d < 3; ++d)
        flow_need(path.field_displacement_m[d] ==
                      (e.position_m[d] + shift[d]) - path.field_start_m[d],
                  FK::motion_record, "Hall canonical path endpoint operand");
      if (!is_newborn)
        segment_guard(opt, old.position_m, e.position_m);
      else {
        // Newborn (possibly leaving its emitting foil face): the straight
        // segment to the endpoint must not reach any kinetic boundary.
        PicVector D{};
        for (int d = 0; d < 3; ++d)
          D[d] = e.position_m[d] - old.position_m[d];
        flow_need(native_first_hit(opt, old.position_m, D).sink < 0,
                  FK::sink_duplicate,
                  "newborn reaches a kinetic boundary without a sink record");
      }
      kinetic_change += static_cast<long double>(old.weight) *
                        (pic_kinetic_energy_J(e) - pic_kinetic_energy_J(old));
      after.particles.push_back(e);
    }
    flow_need(j == r.endpoints.size() && s == r.sinks.size(),
              FK::sink_duplicate,
              "sink/endpoint records repeat, overlap or name a carrier "
              "outside the census");
    const auto &budget = r.budget;
    const double K = budget[0], W = budget[1], defect = budget[2];
    flow_need(K == double(kinetic_change) && defect == K - W &&
                  std::abs(defect) <=
                      1e-20 + 1e-10 * (std::abs(K) + std::abs(W)),
              FK::motion_record,
              "Hall retained exact kinetic/work defect and unchanged guard");
    flow_need(budget[3] >= 0 && budget[4] >= 0 && budget[4] <= .1 &&
                  budget[5] >= 0 && budget[5] <= .05 && budget[6] >= 0 &&
                  budget[6] <= .2,
              FK::motion_record,
              "Hall retained phase/plasma/gyro/displacement guards");
    // Logical packets: a sink removes every carrier of its packet; carriers
    // of a surviving multi-carrier packet stay co-located.
    std::vector<const ArenaPacket *> removed;
    for (auto i : sunk_rows) {
      const auto *p = t->owner_of_row[i];
      if (removed.empty() || removed.back() != p) {
        bool seen = false;
        for (const auto *q : removed)
          seen = seen || q == p;
        if (!seen)
          removed.push_back(p);
      }
    }
    std::size_t removed_carriers = 0;
    for (const auto *p : removed)
      removed_carriers += p->native_ids.size();
    flow_need(removed_carriers == sunk_rows.size(), FK::sink_partial_packet,
              "sink records remove only some carriers of a logical packet");
    for (const auto &rows : t->multi) {
      const GapParticle *first = nullptr;
      for (auto k : rows) {
        const auto it = std::lower_bound(
            after.particles.begin(), after.particles.end(), base[k].id,
            [](const GapParticle &a, std::uint64_t id) { return a.id < id; });
        if (it == after.particles.end() || it->id != base[k].id)
          continue; // removed packet (checked above)
        if (first)
          require(first->position_m == it->position_m &&
                      first->proper_velocity_m_s == it->proper_velocity_m_s,
                  "divergent carrier positions/velocities");
        else
          first = &*it;
      }
    }
    result->trial = t;
    result->binding = b;
    after.options = t->base.native.options;
    after.horizontal = t->base.native.horizontal;
    after.vertical = t->base.native.vertical;
    after.potential = t->base.native.potential;
    after.ledger = t->base.native.ledger;
    after.next_id = t->base.native.next_id;
    after.ledger.time_s = o.common.end_s;
    after.fields_current = false;
    const auto &rows = native_digest(after, result->after.digest,
                                     &t->base.native, &t->base.digest);
    auto next = std::make_shared<ArenaSnapshot>(*t->base.snapshot);
    next->full_time_context["field_sha256"] = b.field_after_sha256;
    next->full_time_context["heat_sha256"] = b.heat_after_sha256;
    next->full_time_context["source_sha256"] = b.source_after_sha256;
    MotionRow row;
    row.record_schema = r.schema;
    row.kinetic_policy_schema = "native-finite-target-kinetic-policy-v1";
    row.accepted_step = r.accepted_step;
    row.staged_arena_revision = r.staged_arena_revision;
    row.trial_stage = r.trial_stage;
    row.species_sha256 = r.species_sha256;
    row.kinetic_policy_sha256 = r.kinetic_policy_sha256;
    row.native_to_field_m = r.native_to_field_m;
    row.current_A_m_sha256 = r.current_A_m_sha256;
    row.initial_charge_C_sha256 = r.initial_charge_C_sha256;
    row.endpoint_charge_C_sha256 = r.endpoint_charge_C_sha256;
    row.trial_field_sha256 = r.trial_field_sha256;
    row.staged_arena_sha256 = r.staged_arena_sha256;
    row.census_sha256 = r.census_sha256;
    for (int k = 0; k < 7; ++k)
      row.budget[k] = budget[k];
    row.transport_trial_sha256 = b.transport_trial_sha256;
    row.endpoints_root = rows.root;
    {
      // Block Merkle root over the path records straight from the span.
      std::vector<cn::BlockLeaf> leaves;
      std::size_t i = 0;
      while (i < N) {
        const auto block = (r.paths[i].native_id - 1) / cn::carrier_block;
        std::size_t k = i;
        while (k < N && (r.paths[k].native_id - 1) / cn::carrier_block == block)
          ++k;
        cn::Hasher h;
        h.w.tag("native-v5/path-block");
        h.w.u64(block);
        h.w.u64(k - i);
        h.flush_if_large(0);
        for (std::size_t m = i; m < k; ++m) {
          h.w.raw(&r.paths[m], sizeof(FullTimeMotionPath));
          h.flush_if_large();
        }
        leaves.push_back({block, i, k - i, h.finish()});
        i = k;
      }
      row.paths_root = cn::merkle_root("native-v5/path-root", leaves, N);
    }
    row.common_clock_sha256 = o.common.common_clock_sha256;
    row.mesh_sha256 = o.common.mesh_sha256;
    row.material_sha256 = o.common.material_sha256;
    row.field_before_sha256 = b.field_before_sha256;
    row.field_after_sha256 = b.field_after_sha256;
    row.heat_after_sha256 = b.heat_after_sha256;
    row.source_after_sha256 = b.source_after_sha256;
    row.carriers = N;
    row.begin_s = o.common.begin_s;
    row.end_s = o.common.end_s;
    row.dt_s = o.common.dt_s;
    FullRow full;
    full.kind = RowKind::motion;
    full.motion = std::move(row);
    full_time_append(*next, std::move(full));
    if (!r.sinks.empty()) {
      const auto used = rows_used(o, *next);
      const std::size_t need = 1 + o.slot_reaction_budget + 2;
      if (used > o.records || o.records - used < need)
        throw refusal(RK::step_records,
                      "sink row exceeds the declared outer record budget");
      if (window_room(*next) < need)
        throw refusal(RK::window_rows,
                      "sink row exceeds the history window budget");
      auto ledger = std::make_shared<FullTimeSinkLedger>();
      ledger->accepted_step = o.common.accepted_step;
      ledger->transport_trial_sha256 = b.transport_trial_sha256;
      ledger->carriers = r.sinks.size();
      ledger->records.assign(r.sinks.begin(), r.sinks.end());
      std::map<std::pair<int, std::string>, std::size_t> tally_index;
      std::vector<std::array<ExactSum, 2>> tally_sums;
      std::array<std::array<ExactSum, 4>, 5> surface_sums;
      std::array<std::uint64_t, 5> surface_carriers{};
      for (std::size_t q = 0; q < sunk_rows.size(); ++q) {
        const auto &old = base[sunk_rows[q]];
        const auto &k = r.sinks[q];
        const auto *packet = t->owner_of_row[sunk_rows[q]];
        const auto &spec = o.owner->states_.at(packet->state);
        GapParticle hit = old;
        hit.proper_velocity_m_s = k.incident_proper_velocity_m_s;
        const double energy = old.weight * pic_kinetic_energy_J(hit);
        flow_need(std::isfinite(energy) && energy >= 0, FK::sink_record,
                  "sink incident kinetic energy operand");
        const int surface = int(k.surface);
        auto key = std::make_pair(surface, packet->state);
        auto [where, fresh] = tally_index.emplace(key, ledger->tallies.size());
        if (fresh) {
          FullTimeSinkTally tally;
          tally.surface = k.surface;
          tally.kind = k.kind;
          tally.state = packet->state;
          tally.charge_e = spec.charge_e;
          ledger->tallies.push_back(tally);
          tally_sums.emplace_back();
        }
        auto &tally = ledger->tallies[where->second];
        ++tally.carriers;
        tally_sums[where->second][0].add(old.weight);
        tally_sums[where->second][1].add(energy);
        ++surface_carriers[std::size_t(surface)];
        auto &ss = surface_sums[std::size_t(surface)];
        ss[0].add(old.weight);
        ss[1].add(energy);
        for (int e = 0; e < std::abs(spec.charge_e); ++e)
          ss[spec.charge_e > 0 ? 2 : 3].add(old.weight);
      }
      for (std::size_t q = 0; q < ledger->tallies.size(); ++q) {
        ledger->tallies[q].weight = tally_sums[q][0].count();
        ledger->tallies[q].kinetic_energy_J = tally_sums[q][1].count();
      }
      std::sort(ledger->tallies.begin(), ledger->tallies.end(),
                [](const auto &a, const auto &c) {
                  return std::make_pair(int(a.surface), a.state) <
                         std::make_pair(int(c.surface), c.state);
                });
      auto flows_state =
          next->full_time_flows
              ? std::make_shared<FullTimeFlows>(*next->full_time_flows)
              : std::make_shared<FullTimeFlows>();
      flows_state->digest.clear();
      for (int k = 0; k < 5; ++k) {
        FullTimeSinkSurfaceTotals totals = flows_state->surfaces[std::size_t(k)];
        totals.carriers = surface_carriers[std::size_t(k)];
        totals.weight = surface_sums[std::size_t(k)][0].count();
        totals.kinetic_energy_J = surface_sums[std::size_t(k)][1].count();
        totals.positive_charge_e = surface_sums[std::size_t(k)][2].count();
        totals.negative_charge_e = surface_sums[std::size_t(k)][3].count();
        ledger->surfaces.push_back(totals);
        auto &cum = flows_state->surfaces[std::size_t(k)];
        cum.carriers = checked_add(cum.carriers, totals.carriers, "sink total");
        for (auto [a, c] : {std::pair{&cum.weight, &totals.weight},
                            std::pair{&cum.kinetic_energy_J, &totals.kinetic_energy_J},
                            std::pair{&cum.positive_charge_e, &totals.positive_charge_e},
                            std::pair{&cum.negative_charge_e, &totals.negative_charge_e}}) {
          ExactSum sum;
          sum.add(*a);
          sum.add(*c);
          *a = sum.count();
        }
      }
      flows_state->sunk_carriers =
          checked_add(flows_state->sunk_carriers, r.sinks.size(), "sunk total");
      flows_state->sink_rows = checked_add(flows_state->sink_rows, 1, "sink rows");
      for (const auto *p : removed)
        ledger->removed_packets.push_back(p->id);
      std::sort(ledger->removed_packets.begin(), ledger->removed_packets.end());
      for (const auto &id : ledger->removed_packets)
        next->packets.erase(id);
      next->full_time_flows = flows_state;
      ledger->records_root_sha256 =
          id_records_root("native-v6/sink-block", "native-v6/sink-root",
                          ledger->records, encode_sink_record);
      ledger->ledger_sha256 = sink_ledger_digest(*ledger);
      FullRow sink;
      sink.kind = RowKind::sink;
      sink.sink.accepted_step = o.common.accepted_step;
      sink.sink.census = N;
      sink.sink.survivors = r.endpoints.size();
      sink.sink.carriers = r.sinks.size();
      sink.sink.packets = removed.size();
      sink.sink.transport_trial_sha256 = b.transport_trial_sha256;
      sink.sink.records_root = ledger->records_root_sha256;
      sink.sink.ledger_sha256 = ledger->ledger_sha256;
      sink.sink.removed_packets_sha256 = text_list_digest(
          "native-v6/removed-packets", ledger->removed_packets);
      ledger->row_sha256 = full_time_append(*next, std::move(sink));
      result->sinks = ledger;
    }
    result->after.snapshot = next;
    result->v3_schema = r.schema;
    result->v3_sha256 = record_sha;
    result->v3_endpoints = r.endpoints.size();
    result->v3_paths = r.paths.size();
    return result;
  }
  // One first-hit sink record against the native kinetic policy port.
  static void sink_check(const GapPicOptions &opt, const GapParticle &old,
                         const FullTimeMotionPath &path,
                         const FullTimeSinkRecord &k, double F,
                         const PicVector &shift) {
    const int surface = int(k.surface);
    flow_need(surface >= 0 && surface <= 4 &&
                  (k.kind == FullTimeSinkKind::collected_at_conductor ||
                   k.kind == FullTimeSinkKind::escaped_exterior) &&
                  k.kind == sink_kind_of(surface) && k.axis <= 2,
              FK::sink_record, "sink kind/surface/axis inconsistent");
    for (double v : {k.plane_m, k.segment_m[0], k.segment_m[1], k.segment_m[2],
                     k.hit_point_m[0], k.hit_point_m[1], k.hit_point_m[2],
                     k.incident_proper_velocity_m_s[0],
                     k.incident_proper_velocity_m_s[1],
                     k.incident_proper_velocity_m_s[2]})
      flow_need(std::isfinite(v), FK::sink_record, "nonfinite sink record");
    flow_need(std::isfinite(k.fraction) && k.fraction >= 0 && k.fraction <= 1 &&
                  std::isfinite(k.flight_s) && k.flight_s == F * k.fraction &&
                  k.flight_s >= 0 && k.flight_s <= F &&
                  k.discarded_flight_s == F - k.flight_s,
              FK::sink_outside_step,
              "sink hit outside the carrier flight of this step");
    const auto hit = native_first_hit(opt, old.position_m, k.segment_m);
    flow_need(hit.sink >= 0, FK::sink_not_crossing,
              "sink record for a carrier whose segment crosses no kinetic "
              "boundary");
    flow_need(hit.sink == surface && hit.fraction == k.fraction &&
                  hit.axis == int(k.axis) && hit.plane == k.plane_m,
              FK::sink_not_crossing,
              "sink record differs from the native first kinetic hit");
    PicVector point{};
    for (int d = 0; d < 3; ++d)
      point[d] = old.position_m[d] + k.fraction * k.segment_m[d];
    point[k.axis] = k.plane_m;
    flow_need(point == k.hit_point_m, FK::sink_record,
              "sink hit point differs from start + fraction * segment");
    if (old.charge_C == 0) {
      double u2 = 0;
      for (double v : old.proper_velocity_m_s)
        u2 += v * v;
      const double gamma = std::sqrt(1 + u2 / (pic_c * pic_c));
      for (int d = 0; d < 3; ++d)
        flow_need(k.segment_m[d] == F * old.proper_velocity_m_s[d] / gamma,
                  FK::sink_record,
                  "neutral sink segment differs from its straight drift");
      flow_need(k.incident_proper_velocity_m_s == old.proper_velocity_m_s,
                FK::sink_record, "neutral incident velocity changed");
    }
    flow_need(path.flight_s == k.flight_s, FK::sink_record,
              "sink path flight differs from the hit flight");
    for (int d = 0; d < 3; ++d)
      flow_need(path.field_displacement_m[d] ==
                    (k.hit_point_m[d] + shift[d]) - path.field_start_m[d],
                FK::sink_record, "sink path displacement differs from the hit");
  }
  static std::unique_ptr<PreparedPrivateMotion::Impl>
  motion(const std::shared_ptr<FullTrial> &t, const FullTimeMotionBinding &b,
         ArenaJson envelope) {
    auto &o = *t->outer;
    live(o);
    require(o.phase == 1 && !t->consumed,
            "closed motion source/boundary ownership unsupported");
    require(same_common(b.common, o.common) &&
                b.field_before_sha256 == o.common.accepted_field_sha256 &&
                digest(b.transport_trial_sha256) &&
                digest(b.field_after_sha256) && digest(b.heat_after_sha256) &&
                b.source_after_sha256 == trial_source(*t),
            "motion common/owner seals");
    // v6: the v2 JSON record has one common flight per carrier and no sink
    // records; births need the binary flows record.
    flow_need(!t->births, FK::motion_record,
              "v2 JSON motion record cannot carry births (partial flights)");
    const auto &receipt = envelope;
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
    auto &audit = envelope.at("owner_audit");
    require(audit.is_object() && audit.contains("seal_sha256"),
            "motion genuine owner record seal");
    const auto ownerseal = json_text(audit.at("seal_sha256"));
    audit.erase("seal_sha256");
    {
      // Hall v2 seals its record as text("motion-receipt-record-v1") +
      // text(JSON dump); verified once here, then only the digest is kept.
      Writer enc;
      enc.tag("motion-receipt-record-v1");
      enc.text(audit.dump());
      const bool sealed = ownerseal == b.transport_trial_sha256 &&
                          cn::sha256_hex(enc.bytes) == ownerseal;
      audit["seal_sha256"] = ownerseal;
      require(sealed, "motion genuine owner record seal");
    }
    require(audit.size() == 16 &&
                audit.at("schema") ==
                    "canonical-full-time-total-field-motion-v2",
            "recognized Hall v2 record schema required");
    const auto &bind = audit.at("binding");
    require(bind.is_object() && bind.size() == 14, "Hall v2 binding schema");
    MotionRow row;
    row.record_schema = "canonical-full-time-total-field-motion-v2";
    row.kinetic_policy_schema = "native-finite-target-kinetic-policy-v1";
    row.accepted_step = json_u64(bind.at("accepted_step"));
    row.staged_arena_revision = json_u64(bind.at("staged_arena_revision"));
    row.trial_stage = json_u64(bind.at("trial_stage"));
    for (const auto *key : {"begin_s", "end_s", "dt_s"})
      (void)json_scalar(bind.at(key));
    row.species_sha256 = json_text(bind.at("species_sha256"));
    row.kinetic_policy_sha256 = json_text(bind.at("kinetic_policy_sha256"));
    require(digest(row.species_sha256) &&
                row.kinetic_policy_sha256 ==
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
    require(policy == expected_policy,
            "Hall complete native kinetic policy operands");
    const PicVector shift = {-.5 * opt.width_m, -.5 * opt.gap_m, 0};
    row.native_to_field_m = json_vector(audit.at("native_to_field_m"));
    require(row.native_to_field_m == shift, "Hall translation-only frame");
    row.current_A_m_sha256 = json_text(audit.at("current_A_m_sha256"));
    row.initial_charge_C_sha256 = json_text(audit.at("initial_charge_C_sha256"));
    row.endpoint_charge_C_sha256 =
        json_text(audit.at("endpoint_charge_C_sha256"));
    for (const auto *d : {&row.current_A_m_sha256, &row.initial_charge_C_sha256,
                          &row.endpoint_charge_C_sha256})
      require(digest(*d), "Hall charge/current array identity");
    row.trial_field_sha256 = json_text(bind.at("trial_field_sha256"));
    row.staged_arena_sha256 = json_text(bind.at("staged_arena_sha256"));
    require(
        bind.at("phase") == "full_time_position_and_proper_velocity" &&
            row.accepted_step == o.common.accepted_step &&
            row.staged_arena_revision == t->base.snapshot->revision &&
            row.staged_arena_sha256 == t->prefix &&
            bind.at("begin_s") == o.common.begin_s &&
            bind.at("end_s") == o.common.end_s &&
            bind.at("dt_s") == o.common.dt_s &&
            bind.at("mesh_sha256") == o.common.mesh_sha256 &&
            bind.at("material_sha256") == o.common.material_sha256 &&
            bind.at("common_clock_sha256") == o.common.common_clock_sha256 &&
            digest(row.trial_field_sha256),
        "motion owner staged identity");
    row.census_sha256 = json_text(audit.at("census_sha256"));
    require(row.census_sha256 == census_sha(population(o, t->base, o.cap)) &&
                audit.at("coverage") ==
                    "all_native_carriers_no_sources_no_boundaries" &&
                audit.at("neutral_motion_owner") == "same_full_time_q0_drift" &&
                !json_boolean(audit.at("native_Poisson_called")) &&
                !json_boolean(audit.at("native_DKD_called")),
            "motion coverage/census/mode");
    const auto &base = t->base.native.particles;
    const auto &endpoints = audit.at("endpoints");
    const auto &paths = audit.at("paths");
    require(endpoints.is_array() && paths.is_array() &&
                endpoints.size() == base.size() && base.size() <= o.cap &&
                paths.size() == base.size(),
            "motion endpoint dimensions");
    std::vector<GapParticle> moved(base.size());
    std::vector<PathRecord> path_rows(base.size());
    long double kinetic_change = 0;
    for (std::size_t i = 0; i < base.size(); ++i) {
      const auto &old = base[i];
      const auto &e = endpoints[i];
      require(e.is_object() && e.size() == 6, "Hall endpoint row schema");
      GapParticle p;
      p.id = json_u64(e.at("native_id"));
      require(p.id == old.id,
              "Hall endpoint rows must cover the census once in canonical "
              "ascending native-ID order");
      p.position_m = json_vector(e.at("position_m"));
      p.proper_velocity_m_s = json_vector(e.at("proper_velocity_m_s"));
      p.charge_C = json_scalar(e.at("charge_C"));
      p.mass_kg = json_scalar(e.at("mass_kg"));
      p.weight = json_scalar(e.at("weight"));
      closed_position(opt, p);
      require(p.charge_C == old.charge_C && p.mass_kg == old.mass_kg &&
                  p.weight == old.weight,
              "motion altered native species/count");
      const auto &path = paths[i];
      require(path.is_object() && path.size() == 5 &&
                  path.at("motion_owner") == "ordinary_total_field" &&
                  json_u64(path.at("native_id")) == old.id,
              "Hall path row schema/owner/order");
      auto &pr = path_rows[i];
      pr.id = old.id;
      pr.start = json_vector(path.at("field_start_m"));
      pr.displacement = json_vector(path.at("field_displacement_m"));
      pr.flight_s = json_scalar(path.at("flight_s"));
      require(pr.flight_s == o.common.dt_s, "Hall complete common flight");
      for (int d = 0; d < 3; ++d)
        require(pr.start[d] == old.position_m[d] + shift[d] &&
                    pr.displacement[d] ==
                        (p.position_m[d] + shift[d]) - pr.start[d],
                "Hall canonical path start/endpoint operands");
      segment_guard(opt, old.position_m, p.position_m);
      kinetic_change += static_cast<long double>(old.weight) *
                        (pic_kinetic_energy_J(p) - pic_kinetic_energy_J(old));
      moved[i] = p;
    }
    const auto &budget = audit.at("budget");
    require(budget.is_object() && budget.size() == 7,
            "Hall complete budget schema");
    const char *keys[7] = {"kinetic_change_J", "field_work_J",
                           "work_defect_J",    "continuity_C",
                           "plasma_dt",        "maximum_gyro",
                           "maximum_displacement_cells"};
    for (int k = 0; k < 7; ++k)
      row.budget[k] = json_scalar(budget.at(keys[k]));
    const double K = row.budget[0], W = row.budget[1], defect = row.budget[2];
    require(K == double(kinetic_change) && defect == K - W &&
                std::abs(defect) <= 1e-20 + 1e-10 * (std::abs(K) + std::abs(W)),
            "Hall retained exact kinetic/work defect and unchanged guard");
    require(row.budget[3] >= 0 && row.budget[4] >= 0 && row.budget[4] <= .1 &&
                row.budget[5] >= 0 && row.budget[5] <= .05 &&
                row.budget[6] >= 0 && row.budget[6] <= .2,
            "Hall retained phase/plasma/gyro/displacement guards");
    auto result = std::make_unique<PreparedPrivateMotion::Impl>();
    result->trial = t;
    result->binding = b;
    auto &after = result->after.native;
    after.options = t->base.native.options;
    after.horizontal = t->base.native.horizontal;
    after.vertical = t->base.native.vertical;
    after.potential = t->base.native.potential;
    after.ledger = t->base.native.ledger;
    after.next_id = t->base.native.next_id;
    after.particles = std::move(moved);
    after.ledger.time_s = o.common.end_s;
    after.fields_current = false;
    const auto &rows = native_digest(after, result->after.digest,
                                     &t->base.native, &t->base.digest);
    auto next = std::make_shared<ArenaSnapshot>(*t->base.snapshot);
    next->full_time_context["field_sha256"] = b.field_after_sha256;
    next->full_time_context["heat_sha256"] = b.heat_after_sha256;
    next->full_time_context["source_sha256"] = b.source_after_sha256;
    row.transport_trial_sha256 = b.transport_trial_sha256;
    row.endpoints_root = rows.root;
    row.paths_root = path_rows_root(path_rows);
    row.common_clock_sha256 = o.common.common_clock_sha256;
    row.mesh_sha256 = o.common.mesh_sha256;
    row.material_sha256 = o.common.material_sha256;
    row.field_before_sha256 = b.field_before_sha256;
    row.field_after_sha256 = b.field_after_sha256;
    row.heat_after_sha256 = b.heat_after_sha256;
    row.source_after_sha256 = b.source_after_sha256;
    row.carriers = base.size();
    row.begin_s = o.common.begin_s;
    row.end_s = o.common.end_s;
    row.dt_s = o.common.dt_s;
    FullRow full;
    full.kind = RowKind::motion;
    full.motion = std::move(row);
    full_time_append(*next, std::move(full));
    result->after.snapshot = next;
    (void)o.owner->full_time_census(*next, after.particles, o.owner->carrier_cap_,
                                    false);
    result->envelope = std::move(envelope);
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
  return frame_completion(*impl_->scope->outer->frame.snapshot, id);
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
// Human-readable export of this step's complete motion (Hall record, native
// endpoints and owner seals), built on first request. Not journal storage.
const ArenaJson &PreparedPrivateMotion::audit() const {
  require(bool(impl_), "empty motion token");
  auto &i = *impl_;
  if (i.audit.is_null() && !i.v3_schema.empty()) {
    // Binary admission: the record arrays stay with the caller; the audit
    // names the record and the native result (endpoints = accepted rows).
    ArenaJson sinks = ArenaJson::array();
    if (i.sinks)
      for (const auto &k : i.sinks->records)
        sinks.push_back({{"native_id", k.native_id},
                         {"kind", int(k.kind)},
                         {"surface", int(k.surface)},
                         {"fraction", k.fraction},
                         {"hit_point_m", k.hit_point_m},
                         {"flight_s", k.flight_s},
                         {"discarded_flight_s", k.discarded_flight_s}});
    i.audit = {{"kind", "full_time_motion"},
               {"record_schema", i.v3_schema},
               {"record_sha256", i.v3_sha256},
               {"endpoints", i.v3_endpoints},
               {"paths", i.v3_paths},
               {"sinks", sinks},
               {"native_rows_sha256", i.after.digest ? i.after.digest->root : ""},
               {"field_after_sha256", i.binding.field_after_sha256},
               {"heat_after_sha256", i.binding.heat_after_sha256},
               {"source_after_sha256", i.binding.source_after_sha256}};
  }
  if (i.audit.is_null()) {
    ArenaJson endpoints = ArenaJson::array();
    for (const auto &e : i.envelope.at("owner_audit").at("endpoints"))
      endpoints.push_back({{"id", e.at("native_id")},
                           {"position_m", e.at("position_m")},
                           {"u", e.at("proper_velocity_m_s")},
                           {"q", e.at("charge_C")},
                           {"m", e.at("mass_kg")},
                           {"weight", e.at("weight")}});
    i.audit = {{"kind", "full_time_motion"},
               {"receipt", i.envelope},
               {"endpoints", std::move(endpoints)},
               {"field_after_sha256", i.binding.field_after_sha256},
               {"heat_after_sha256", i.binding.heat_after_sha256},
               {"source_after_sha256", i.binding.source_after_sha256}};
  }
  return i.audit;
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
  FullTimeAccess::live(*impl_->trial->outer);
  require(!impl_->trial->consumed && impl_->trial->outer->phase == 1,
          "stale transport snapshot");
  return impl_->trial->base.snapshot;
}
std::uint64_t PrivateTransportTrial::staged_arena_revision() const {
  return snapshot()->revision;
}
std::string PrivateTransportTrial::state_prefix() const {
  (void)snapshot();
  return impl_->trial->prefix;
}
const std::string &PrivateTransportTrial::census_sha256() const {
  (void)snapshot();
  FullTimeAccess::census_cache(*impl_->trial);
  return impl_->trial->census;
}
std::vector<FullTimeBirthFlight> PrivateTransportTrial::birth_flights() const {
  (void)snapshot();
  return impl_->trial->births ? impl_->trial->births->flights
                              : std::vector<FullTimeBirthFlight>{};
}
std::string PrivateTransportTrial::source_after_sha256() const {
  (void)snapshot();
  return FullTimeAccess::trial_source(*impl_->trial);
}
PreparedBirths::PreparedBirths(std::unique_ptr<Impl> p) : impl_(std::move(p)) {}
PreparedBirths::PreparedBirths(PreparedBirths &&) noexcept = default;
PreparedBirths &PreparedBirths::operator=(PreparedBirths &&) noexcept = default;
PreparedBirths::~PreparedBirths() = default;
const FullTimeBirthReceipt &PreparedBirths::receipt() const {
  require(bool(impl_), "empty prepared births");
  return *impl_->info->receipt;
}
const std::vector<FullTimeBirthFlight> &PreparedBirths::flights() const {
  require(bool(impl_), "empty prepared births");
  return impl_->info->flights;
}
BirthScope::BirthScope(std::unique_ptr<Impl> p) : impl_(std::move(p)) {}
BirthScope::BirthScope(BirthScope &&) noexcept = default;
BirthScope &BirthScope::operator=(BirthScope &&) noexcept = default;
BirthScope::~BirthScope() = default;
std::uint64_t BirthScope::base_native_id() const {
  require(bool(impl_), "empty birth scope");
  return impl_->base;
}
std::uint64_t BirthScope::accepted_step() const {
  require(bool(impl_), "empty birth scope");
  return impl_->outer->common.accepted_step;
}
const std::string &BirthScope::source_before_sha256() const {
  require(bool(impl_), "empty birth scope");
  return impl_->source_before;
}
PreparedBirths BirthScope::prepare(const FullTimeBirthCredit &c) const {
  require(bool(impl_), "empty birth scope");
  return PreparedBirths(FullTimeAccess::births(*impl_, c));
}
BirthScope PrivateFullTimeStep::birth_scope() {
  auto o = impl_->outer;
  FullTimeAccess::live(*o);
  if (o->phase != 1 || o->builder_active || o->births_reserved)
    throw FullTimeFlowRefusal(
        FullTimeFlowRefusal::Kind::birth_scope,
        "birth scope needs the completed first exposure, before transport, "
        "once per outer step");
  auto p = std::make_unique<BirthScope::Impl>();
  p->outer = o;
  p->frame = o->frame.snapshot;
  p->base = o->frame.native.next_id;
  p->source_before = o->boundary->source_sha256;
  require(p->base > 0, "native next ID");
  o->births_reserved = true;
  o->birth_base = p->base;
  return BirthScope(std::move(p));
}
PrivateTransportTrial
PrivateFullTimeStep::fork_transport_trial(const PreparedBirths &births) {
  auto o = impl_->outer;
  FullTimeAccess::live(*o);
  require(o->phase == 1 && !o->builder_active && o->transport_forks < 4,
          "transport requires completed first exposure and bounded forks");
  if (!births.impl_ || births.impl_->outer != o ||
      births.impl_->base != o->frame.snapshot)
    throw FullTimeFlowRefusal(FullTimeFlowRefusal::Kind::birth_scope,
                              "foreign or stale prepared births");
  auto p = std::make_unique<PrivateTransportTrial::Impl>();
  p->trial = std::make_shared<FullTrial>();
  p->trial->outer = o;
  p->trial->base = births.impl_->after;
  p->trial->base.digest = births.impl_->after.digest;
  p->trial->births = births.impl_->info;
  p->trial->fork_frame = o->frame.snapshot;
  p->trial->prefix = FullTimeAccess::prefix(*o, p->trial->base);
  ++o->transport_forks;
  return PrivateTransportTrial(std::move(p));
}
PreparedPrivateMotion PrivateTransportTrial::prepare_total_field_motion(
    const FullTimeMotionBinding &b, FullTimeMotionAdmission admission) {
  FullTimeAccess::live(*impl_->trial->outer);
  if (admission.v3_) {
    // Binary motion v3: typed arrays, no JSON record, no seal re-dump.
    auto keep = std::move(admission.v3_storage_);
    auto record = std::move(admission.v3_);
    return PreparedPrivateMotion(
        FullTimeAccess::motion_v3(impl_->trial, b, *record));
  }
  const auto &record = admission.record_;
  require(record.is_object() && record.contains("endpoints") &&
              record.at("endpoints").is_array() &&
              record.at("endpoints").size() <= impl_->trial->outer->cap,
          "bounded typed Hall receipt before copying");
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
      FullTimeAccess::motion(impl_->trial, b, std::move(envelope)));
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
  p->trial->base.digest = o->frame.digest;
  p->trial->fork_frame = o->frame.snapshot;
  ++o->transport_forks;
  return PrivateTransportTrial(std::move(p));
}
void PrivateFullTimeStep::accept_transport_trial(PreparedPrivateMotion &&t) {
  auto &o = *impl_->outer;
  FullTimeAccess::live(o);
  require(t.impl_ && t.impl_->trial->outer == impl_->outer &&
              !t.impl_->consumed && !t.impl_->trial->consumed && o.phase == 1 &&
              o.frame.snapshot == t.impl_->trial->fork_frame,
          "foreign/stale/double motion token");
  require(t.impl_->after.snapshot->revision - o.before.snapshot->revision <=
              o.records,
          "outer motion record cap");
  std::shared_ptr<const FullTimeStepFlows> flows;
  if (t.impl_->trial->births || t.impl_->sinks) {
    auto f = std::make_shared<FullTimeStepFlows>();
    if (t.impl_->trial->births)
      f->births = t.impl_->trial->births->receipt;
    f->sinks = t.impl_->sinks;
    flows = std::move(f);
  }
  swap_frame(o.frame, t.impl_->after);
  o.flows.swap(flows);
  o.motion_sha = t.impl_->binding.transport_trial_sha256;
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
                 json_text(o->frame.snapshot->full_time_context.at("field_sha256")));
  auto s = std::make_shared<FullCollision>();
  s->outer = o;
  s->exposure = e;
  auto p = std::make_unique<CollisionScope::Impl>();
  p->scope = std::move(s);
  o->second_reserved = true;
  return CollisionScope(std::move(p));
}
// Full-time population validation over canonical ascending-ID native rows.
// Same checks as validate_rows; the exact BIGNUM count sums are memoized per
// (immutable packets payload, SHA-256 of the (ID, weight) sequence).
std::vector<ArenaRateEntry>
PicArena::full_time_census(const ArenaSnapshot &s,
                           const std::vector<GapParticle> &native,
                           std::size_t cap, bool entries) const {
  require(native.size() <= cap && cap <= carrier_cap_,
          "full-time native census cap");
  std::size_t n = 0;
  for (const auto &[id, p] : s.packets) {
    (void)id;
    require(p.native_ids.size() <= cap - n, "full-time canonical census cap");
    n += p.native_ids.size();
  }
  cn::Hasher weights;
  weights.w.tag("native-v5/count-memo");
  for (std::size_t i = 0; i < native.size(); ++i) {
    require(native[i].id > 0 && (i == 0 || native[i].id > native[i - 1].id),
            "duplicate native ID");
    weights.w.u64(native[i].id);
    weights.w.f64(native[i].weight);
    weights.flush_if_large();
  }
  const auto key = s.packets.identity();
  const auto weight_digest = cn::hex(weights.finish());
  bool memo = false;
  for (const auto &[packets, digest] : count_memo_)
    memo = memo || (packets == key && digest == weight_digest);
  std::vector<std::uint8_t> owned(native.size(), 0);
  std::vector<ArenaRateEntry> out;
  if (entries)
    out.reserve(n);
  for (const auto &[id, p] : s.packets) {
    // v6: arrival-v1 birth packets are injective by construction (birth
    // step = accepted step, canonical order) and are not kept in the
    // ever-ID set; every other packet must be a known identity.
    require(id == p.id &&
                (s.ever_packet_ids.count(id) ||
                 (full_time_reserved_packet_id(id) && canonical_birth_id(id) &&
                  p.native_ids.size() == 1)) &&
                states_.count(p.state) && !p.native_ids.empty(),
            "logical packet identity/state");
    if (!memo)
      count_guard(p.count);
    require(!zero(p.count), "zero resident logical count");
    const auto &spec = states_.at(p.state);
    auto total = memo ? Big(nullptr, BN_free) : big();
    const GapParticle *first = nullptr;
    for (auto native_id : p.native_ids) {
      const auto it = std::lower_bound(
          native.begin(), native.end(), native_id,
          [](const GapParticle &a, std::uint64_t b) { return a.id < b; });
      require(it != native.end() && it->id == native_id,
              "native canonical carrier missing");
      auto &mark = owned[std::size_t(it - native.begin())];
      require(!mark, "native carrier owned twice");
      mark = 1;
      const auto &row = *it;
      native_particle_guard(row);
      require(row.charge_C == spec.charge_e * pic_e &&
                  row.mass_kg == spec.mass_kg,
              "native species binding");
      if (!memo) {
        auto weight = exact(row.weight);
        require(BN_add(total.get(), total.get(), weight.get()),
                "carrier count accumulation");
      }
      if (first)
        require(first->position_m == row.position_m &&
                    first->proper_velocity_m_s == row.proper_velocity_m_s,
                "divergent carrier positions/velocities");
      else
        first = &row;
      if (entries)
        out.push_back({id, p.state, row});
    }
    if (!memo) {
      auto expected = exact(p.count);
      require(BN_cmp(total.get(), expected.get()) == 0,
              "native carriers do not exactly represent Count");
    }
  }
  require(n == native.size(), "unowned native particle");
  if (!memo) {
    if (count_memo_.size() >= 4)
      count_memo_.erase(count_memo_.begin());
    count_memo_.emplace_back(key, weight_digest);
  }
  return out;
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
  const auto &context = o.frame.snapshot->full_time_context;
  require(seals.field_sha256 == json_text(context.at("field_sha256")) &&
              seals.heat_sha256 == json_text(context.at("heat_sha256")) &&
              seals.source_sha256 == json_text(context.at("source_sha256")) &&
              seals.rng_sha256 == o.second.completion.after_rng_clock_sha256,
          "outer final field/heat/source/RNG seals");
  auto boundary = std::make_shared<FullTimeAcceptedBoundary>(*o.boundary);
  require(boundary->accepted_step < std::numeric_limits<std::uint64_t>::max(),
          "outer step exhaustion");
  ++boundary->accepted_step;
  boundary->accepted_arena_revision =
      checked_add(o.frame.snapshot->revision, 1, "outer revision exhaustion");
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
  FullRow row;
  row.kind = RowKind::boundary;
  row.common = o.common;
  row.outer_sha = o.outer_sha;
  row.transport_trial_sha256 = o.motion_sha;
  row.first_prefix = o.first.completed_state_prefix;
  row.second_prefix = o.second.completed_state_prefix;
  row.boundary = *boundary;
  row.boundary.accepted_arena_sha256.clear();
  full_time_append(*next, row);
  require(next->revision - o.before.snapshot->revision <= o.records,
          "outer final record cap");
  require(next->revision == boundary->accepted_arena_revision,
          "outer terminal revision");
  auto t = std::make_unique<PreparedArena::Impl>();
  t->owner = &owner;
  t->before = o.before.snapshot;
  t->native_before = o.before.native;
  t->native_after = o.frame.native;
  t->work = t->native_after;
  t->native_before_digest = o.before.digest;
  t->native_after_digest = o.frame.digest;
  if (!t->native_after_digest)
    t->native_after_digest = cn::carrier_rows_digest(t->native_after.particles);
  t->boundary_before = o.boundary;
  t->boundary_after = boundary;
  t->full_outer_identity = impl_->outer;
  t->candidate_prefix = FullTimeAccess::prefix(owner, *next, t->native_after,
                                               *t->native_after_digest);
  boundary->accepted_arena_sha256 = t->candidate_prefix;
  ArenaJson audit = {
      {"kind", "full_time_boundary"},
      {"common", common_json(o.common)},
      {"outer_step_sha256", o.outer_sha},
      {"transport_trial_sha256", o.motion_sha},
      {"first_completed_state_prefix", o.first.completed_state_prefix},
      {"second_completed_state_prefix", o.second.completed_state_prefix},
      {"boundary", boundary_json(*boundary, false)},
      {"row_sha256", next->chain}};
  // Anchor every K accepted steps: validate the complete unanchored window
  // from the previous anchor, then compact it into one chained commitment.
  // Storage only; the logical identity computed above is unchanged.
  const auto &h = history(*next);
  const auto &policy = h.policy->policy;
  if (boundary->accepted_step % policy.anchor_interval_steps == 0) {
    const auto summary = window_validate(*next, owner.channels_, owner.carrier_cap_);
    require(summary.steps == policy.anchor_interval_steps &&
                summary.motions == policy.anchor_interval_steps,
            "anchor window must contain exactly K complete steps");
    FullTimeAnchorFields a;
    a.source_sha256 = PicArena::source_identity();
    a.history_policy_sha256 = h.policy->digest;
    a.anchor_interval_steps = policy.anchor_interval_steps;
    a.sequence = h.anchors->anchors.size();
    a.previous_anchor_sha256 = h.head();
    a.first_step = h.base->steps;
    a.last_step = boundary->accepted_step;
    a.base_revision = h.base->revision;
    a.tail_revision = next->revision;
    a.base_history_chain = h.base->chain;
    a.tail_history_chain = next->chain;
    a.window_rows = h.window_rows();
    a.window_bytes = h.window_bytes();
    a.cumulative_bytes = next->journal_bytes;
    a.window_reactions = summary.reactions;
    a.window_motions = summary.motions;
    a.window_birth_rows = summary.birth_rows;
    a.window_births = summary.births;
    a.window_sink_rows = summary.sink_rows;
    a.window_sunk_carriers = summary.sunk;
    a.flows_sha256 =
        next->full_time_flows ? flows_digest(*next->full_time_flows) : "";
    a.boundary = *boundary;
    a.accepted_arena_sha256 = boundary->accepted_arena_sha256;
    a.snapshot_sha256 = snapshot_digest(*next);
    a.native_rows_sha256 = t->native_after_digest->root;
    a.native_ledger_sha256 = ledger_digest(t->native_after.ledger);
    a.completions_sha256 = completions_digest(next->stochastic_completions);
    a.event_ids_sha256 = id_set_digest(next->event_ids);
    a.packet_ids_sha256 = id_set_digest(next->ever_packet_ids);
    a.next_native_id = t->native_after.next_id;
    a.previous_dt_s = next->previous_dt_s;
    auto record = std::make_shared<FullTimeAnchorRecord>();
    Writer w;
    codec_anchor(w, a);
    record->bytes = std::move(w.bytes);
    record->sha256 = cn::sha256_hex(record->bytes);
    a.sha256 = record->sha256;
    record->fields = a;
    if (h.anchors->anchors.size() >= policy.maximum_anchors)
      throw refusal(RK::anchor_count, "full-time anchor count cap");
    if (record->bytes.size() > policy.maximum_anchor_bytes - h.anchors->bytes)
      throw refusal(RK::anchor_bytes, "full-time anchor byte cap");
    auto chain = std::make_shared<FullTimeAnchorChain>(*h.anchors);
    chain->anchors.push_back(record);
    chain->bytes += record->bytes.size();
    auto base = std::make_shared<FullTimeHistoryBase>();
    base->chain = next->chain;
    base->revision = next->revision;
    base->bytes = next->journal_bytes;
    base->steps = next->accepted_steps;
    base->previous_dt_s = next->previous_dt_s;
    base->boundary = *boundary;
    base->completions = next->stochastic_completions;
    base->event_ids = next->event_ids;
    auto compacted = std::make_shared<FullTimeHistory>(h);
    compacted->anchors = std::move(chain);
    compacted->base = std::move(base);
    compacted->tail.reset();
    // Storage replacement only; snapshot_digest excludes history by design.
    next->full_time_history = std::move(compacted);
    audit["anchor"] = {{"sequence", a.sequence},
                       {"sha256", record->sha256},
                       {"first_step", a.first_step},
                       {"last_step", a.last_step},
                       {"window_rows", a.window_rows},
                       {"window_bytes", a.window_bytes}};
  }
  t->after = next;
  t->audit = std::move(audit);
  auto exposures = std::make_shared<FullTimeCandidateExposures>();
  exposures->common = o.common;
  exposures->outer_step_sha256 = o.outer_sha;
  exposures->transport_trial_sha256 = o.motion_sha;
  exposures->first = o.first;
  exposures->second = o.second;
  t->exposures = exposures;
  t->flows = o.flows;
  o.sealed = true;
  return PreparedArena(std::move(t));
}
std::optional<ArenaStochasticCompletion>
PreparedArena::candidate_stochastic_completion(const std::string &owner) const {
  require(bool(impl_), "empty prepared token");
  require(bool(impl_->after),
          "consumed private token retains audit/prefix only");
  return frame_completion(*impl_->after, owner);
}
std::shared_ptr<const FullTimeAcceptedBoundary>
PreparedArena::candidate_full_time_boundary() const {
  require(bool(impl_), "empty prepared token");
  return impl_->boundary_after;
}
std::shared_ptr<const FullTimeStepFlows>
PreparedArena::candidate_full_time_flows() const {
  require(bool(impl_), "empty prepared token");
  return impl_->flows;
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
  require(carriers > 0 && carriers <= carrier_cap_ && records >= 7 &&
              records <= 4096,
          "outer bounded resources");
  const auto &h = history(*accepted_);
  const auto &policy = h.policy->policy;
  if (FullTimeAccess::window_room(*accepted_) < 5)
    throw refusal(RK::window_rows,
                  "full-time history window cannot hold one more step");
  if ((full_time_boundary_->accepted_step + 1) % policy.anchor_interval_steps ==
      0) {
    if (h.anchors->anchors.size() >= policy.maximum_anchors)
      throw refusal(RK::anchor_count,
                    "full-time anchor count cap reached before step");
    if (anchor_byte_bound > policy.maximum_anchor_bytes - h.anchors->bytes)
      throw refusal(RK::anchor_bytes,
                    "full-time anchor byte cap reached before step");
  }
  require(pic_->particles_.size() <= carriers,
          "native carrier census resource cap");
  (void)full_time_census(*accepted_, pic_->particles_, carriers, false);
  auto o = std::make_shared<FullOuter>();
  o->owner = this;
  o->lifetime = full_time_lifetime_;
  o->common = b;
  o->boundary = full_time_boundary_;
  o->before.snapshot = accepted_;
  o->before.native = capture(*pic_);
  o->before.digest = accepted_native_digest_;
  o->cap = carriers;
  o->records = records;
  {
    Writer w;
    w.tag("arena-full-time-outer-v5");
    w.text(state_prefix());
    auto common = b;
    codec_common(w, common);
    o->outer_sha = cn::sha256_hex(w.bytes);
  }
  auto next = std::make_shared<ArenaSnapshot>(*accepted_);
  next->full_time_context = {
      {"mode", "private_total_field_outer_v1"},
      {"common", common_json(b)},
      {"field_sha256", b.accepted_field_sha256},
      {"heat_sha256", full_time_boundary_->heat_sha256},
      {"source_sha256", full_time_boundary_->source_sha256}};
  FullRow row;
  row.kind = RowKind::outer_begin;
  row.common = b;
  row.outer_sha = o->outer_sha;
  full_time_append(*next, row);
  o->frame.snapshot = next;
  o->frame.native = o->before.native;
  o->frame.digest = o->before.digest;
  auto p = std::make_unique<PrivateFullTimeStep::Impl>();
  p->outer = o;
  active_full_time_outer_ = o;
  return PrivateFullTimeStep(std::move(p));
}
namespace full_time_detail {
void configuration_guard(const FullTimeConfiguration &c) {
  const auto &l = c.limits;
  if (l.maximum_carriers == 0 ||
      l.maximum_carriers > full_time_absolute_carrier_cap)
    throw refusal(RK::carrier_cap,
                  "full-time carrier cap outside [1, absolute cap]");
  if (l.maximum_identity_entries == 0 || l.maximum_identity_entries > (1u << 30))
    throw refusal(RK::configuration, "full-time identity cap");
  const long double required =
      static_cast<long double>(l.maximum_carriers) * full_time_bytes_per_carrier +
      full_time_fixed_bytes;
  if (required > static_cast<long double>(l.memory_budget_bytes))
    throw refusal(RK::memory_budget,
                  "full-time carrier cap exceeds the stated memory budget");
  const auto &h = c.history;
  if (h.anchor_interval_steps == 0 || h.anchor_interval_steps > (1ull << 32) ||
      h.maximum_anchors == 0 || h.maximum_anchor_bytes < anchor_byte_bound ||
      h.maximum_window_rows / 5 < h.anchor_interval_steps ||
      h.maximum_window_bytes < (1u << 20))
    throw refusal(RK::configuration,
                  "history policy cannot hold one complete anchor window");
}
std::string configuration_digest(const FullTimeNativeLimits &l) {
  Writer w;
  w.tag("native-full-time-configuration-v5");
  w.u64(l.maximum_carriers);
  w.u64(l.memory_budget_bytes);
  w.u64(l.maximum_identity_entries);
  w.u64(full_time_bytes_per_carrier);
  w.u64(full_time_fixed_bytes);
  w.u64(cn::carrier_block);
  return cn::sha256_hex(w.bytes);
}
std::string registry_digest(const std::map<std::string, ArenaStateSpec> &states,
                            const std::map<std::string, ArenaChannel> &channels) {
  Writer w;
  w.tag("native-full-time-registry-v5");
  w.u64(states.size());
  for (const auto &[id, s] : states) {
    w.text(id);
    w.i64(s.hydrogen);
    w.i64(s.oxygen);
    w.i64(s.charge_e);
    w.f64(s.mass_kg);
    w.f64(s.internal_J);
    w.f64(s.reference_J);
  }
  w.u64(channels.size());
  for (const auto &[id, c] : channels) {
    w.text(id);
    w.text(c.primary_state);
    w.text(c.target_state);
    w.u64(c.product_states.size());
    for (const auto &p : c.product_states)
      w.text(p);
    w.text(c.dataset_sha256);
    w.text(c.kernel_manifest_sha256);
    w.text(c.kernel_binding_file_sha256);
    w.text(c.model);
  }
  return cn::sha256_hex(w.bytes);
}
} // namespace full_time_detail
std::unique_ptr<PicArena>
PicArena::create_full_time(GapPicOptions options,
                           std::vector<ArenaStateSpec> states,
                           std::vector<ArenaChannel> channels,
                           std::vector<FullTimeAuthoredPacket> packets,
                           const FullTimeGenesisBinding &g) {
  return create_full_time(std::move(options), std::move(states),
                          std::move(channels), std::move(packets), g,
                          FullTimeConfiguration{});
}
std::unique_ptr<PicArena>
PicArena::create_full_time(GapPicOptions options,
                           std::vector<ArenaStateSpec> states,
                           std::vector<ArenaChannel> channels,
                           std::vector<FullTimeAuthoredPacket> packets,
                           const FullTimeGenesisBinding &g,
                           const FullTimeConfiguration &configuration) {
  configuration_guard(configuration);
  const auto cap = configuration.limits.maximum_carriers;
  require(!options.distributed_fields && std::isfinite(g.time_s) &&
              g.time_s >= 0,
          "full-time fresh serial origin");
  nonempty(g.origin_provenance);
  for (const auto *s : {&g.field_sha256, &g.mesh_sha256, &g.material_sha256,
                        &g.frame_sha256, &g.common_clock_sha256, &g.heat_sha256,
                        &g.source_sha256, &g.rng_sha256})
    require(digest(*s), "full-time genesis requires sealed owners/frame");
  options_guard(options);
  require(states.size() <= 256 && channels.size() <= 128,
          "bounded full-time genesis catalog/population");
  if (packets.size() > cap)
    throw refusal(RK::carrier_cap, "genesis packet count exceeds carrier cap");
  if (packets.size() > configuration.limits.maximum_identity_entries)
    throw refusal(RK::identity_count, "genesis packet identities exceed cap");
  auto pic = std::make_unique<GapPic>(options);
  std::map<std::string, ArenaStateSpec> registry;
  for (const auto &s : states)
    require(registry.emplace(s.id, s).second, "duplicate genesis species");
  std::vector<GapParticle> native;
  std::vector<ArenaPacket> logical;
  native.reserve(packets.size());
  logical.reserve(packets.size());
  std::uint64_t id = 1;
  for (const auto &p : packets) {
    nonempty(p.origin);
    nonempty(p.id);
    if (full_time_reserved_packet_id(p.id))
      throw FullTimeFlowRefusal(
          FullTimeFlowRefusal::Kind::birth_namespace,
          "genesis packet ID in the reserved arrival-v1/ namespace");
    const auto &spec = registry.at(p.state);
    require(!p.count.zero(), "zero authored full-time packet");
    ArenaPacket r{p.id, p.state, p.origin, "", p.count, {}};
    auto weights = positive_carriers(p.count);
    if (weights.size() > cap - native.size())
      throw refusal(RK::carrier_cap, "genesis native carrier cap");
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
  pic->particles_ = std::move(native);
  pic->next_particle_id_ = id;
  pic->ledger_.time_s = g.time_s;
  for (const auto &p : pic->particles_) {
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
  std::string initial;
  {
    auto genesis = g;
    Writer w;
    w.tag("authored-full-time-genesis-v5");
    codec_genesis(w, genesis);
    w.text(configuration_digest(configuration.limits));
    initial = cn::sha256_hex(w.bytes);
  }
  auto owner = std::unique_ptr<PicArena>(
      new PicArena(*pic, std::move(states), std::move(channels),
                   std::move(logical), initial, true, cap));
  owner->owned_pic_ = std::move(pic);
  owner->memory_budget_bytes_ = configuration.limits.memory_budget_bytes;
  owner->identity_cap_ = configuration.limits.maximum_identity_entries;
  owner->configuration_ = std::make_shared<const FullTimeConfiguration>(configuration);
  owner->registry_v5_sha_ = registry_digest(owner->states_, owner->channels_);
  owner->options_v5_sha_ = options_digest(owner->pic_->options_);
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
  {
    Writer w;
    w.tag("native-full-time-journal-genesis-v5");
    w.text(owner->initial_prefix_);
    w.text(owner->registry_v5_sha_);
    snapshot->chain = cn::sha256_hex(w.bytes);
  }
  snapshot->full_time_context = {{"mode", "accepted_total_field_v1"},
                                 {"boundary", boundary_json(*boundary, false)},
                                 {"genesis", genesis_json(g)}};
  snapshot->digest_cache.clear();
  owner->accepted_native_digest_ =
      cn::carrier_rows_digest(owner->pic_->particles_);
  boundary->accepted_arena_sha256 =
      FullTimeAccess::prefix(*owner, *snapshot, capture(*owner->pic_),
                             *owner->accepted_native_digest_);
  auto policy = std::make_shared<FullTimeHistoryPolicyState>();
  policy->policy = configuration.history;
  policy->digest = policy_digest(configuration.history);
  policy->genesis_chain = snapshot->chain;
  policy->genesis_boundary = *boundary;
  policy->root = anchor_root(policy->digest, *boundary, snapshot->chain);
  auto base = std::make_shared<FullTimeHistoryBase>();
  base->chain = snapshot->chain;
  base->boundary = *boundary;
  base->event_ids = snapshot->event_ids;
  auto history = std::make_shared<FullTimeHistory>();
  history->policy = std::move(policy);
  history->anchors = std::make_shared<const FullTimeAnchorChain>();
  history->base = std::move(base);
  snapshot->full_time_history = std::move(history);
  owner->accepted_ = snapshot;
  owner->genesis_ = snapshot;
  owner->full_time_boundary_ = boundary;
  return owner;
}
namespace full_time_detail {
FullTimeFlowTotals flow_totals(const FullTimeFlows *f, std::uint64_t step,
                               std::uint64_t next_id) {
  FullTimeFlowTotals t;
  t.accepted_step = step;
  t.next_native_id = next_id;
  if (!f) {
    FullTimeFlows empty;
    t.surfaces.assign(empty.surfaces.begin(), empty.surfaces.end());
    return t;
  }
  t.births = f->births;
  t.birth_rows = f->birth_rows;
  t.sunk_carriers = f->sunk_carriers;
  t.sink_rows = f->sink_rows;
  t.has_last_birth = f->has_last_birth;
  t.last_birth_step = f->last_birth_step;
  t.last_birth_base = f->last_birth_base;
  t.last_birth_count = f->last_birth_count;
  t.last_birth_plan_sha256 = f->last_birth_plan_sha256;
  t.last_birth_token_sha256 = f->last_birth_token_sha256;
  t.surfaces.assign(f->surfaces.begin(), f->surfaces.end());
  t.flows_sha256 = flows_digest(*f);
  return t;
}
} // namespace full_time_detail
FullTimeFlowTotals PicArena::full_time_flow_totals() const {
  require(bool(full_time_boundary_), "arena is not full-time total-field mode");
  return flow_totals(accepted_->full_time_flows.get(),
                     full_time_boundary_->accepted_step,
                     pic_->next_particle_id_);
}
FullTimeFlowTotals PicArena::read_full_time_flow_totals(const std::string &path) {
  const auto doc = read_full_json(path);
  require(doc.is_object() && doc.size() == 2, "full-time restart envelope");
  const auto &b = doc.at("payload");
  require(b.is_object() && b.size() == 16 &&
              b.at("schema") == "canonical_native_arena_full_time_exact_v5" &&
              b.at("source_sha256") == source_identity() &&
              doc.at("sha256") == sha_bytes(b.dump()),
          "full-time checkpoint identity/integrity");
  const auto boundary = read_boundary(b.at("boundary"));
  const auto &meta = b.at("metadata");
  require(meta.is_object() && json_u64(meta.at("accepted_steps")) ==
                                  boundary.accepted_step,
          "full-time checkpoint step");
  std::shared_ptr<FullTimeFlows> flows;
  if (meta.contains("flows"))
    flows = read_flows_json(meta.at("flows"));
  return flow_totals(flows.get(), boundary.accepted_step,
                     json_u64(b.at("next_id")));
}
FullTimeHistoryStatus PicArena::full_time_history_status() const {
  require(bool(full_time_boundary_), "arena is not full-time total-field mode");
  const auto &h = history(*accepted_);
  FullTimeHistoryStatus s;
  s.anchor_interval_steps = h.policy->policy.anchor_interval_steps;
  s.anchors = h.anchors->anchors.size();
  s.anchor_bytes = h.anchors->bytes;
  s.window_rows = h.window_rows();
  s.window_bytes = h.window_bytes();
  s.cumulative_rows = accepted_->revision;
  s.cumulative_bytes = accepted_->journal_bytes;
  s.base_step = h.base->steps;
  s.accepted_step = accepted_->accepted_steps;
  s.history_policy_sha256 = h.policy->digest;
  s.anchor_root_sha256 = h.policy->root;
  s.anchor_head_sha256 = h.head();
  return s;
}
std::string PicArena::export_full_time_anchor_chain() const {
  require(bool(full_time_boundary_), "arena is not full-time total-field mode");
  const auto &h = history(*accepted_);
  Writer w;
  w.tag("native-full-time-anchor-chain-v6");
  w.text(h.policy->digest);
  w.u64(h.policy->policy.anchor_interval_steps);
  w.u64(h.anchors->anchors.size());
  auto genesis = h.policy->genesis_boundary;
  codec_boundary(w, genesis, true);
  w.text(h.policy->genesis_chain);
  for (const auto &a : h.anchors->anchors)
    w.text(a->bytes);
  return std::move(w.bytes);
}
void PicArena::validate_full_time_history() const {
  require(full_time_boundary_ && owned_pic_ && genesis_ &&
              genesis_->journal.empty() && genesis_->revision == 0 &&
              genesis_->accepted_steps == 0 && genesis_->previous_dt_s == 0,
          "full-time immutable genesis");
  const auto &h = history(*accepted_);
  const auto anchors = verify_full_time_anchor_chain(
      export_full_time_anchor_chain(), h.policy->root, h.head());
  const auto &base = *h.base;
  if (anchors.empty())
    require(base.chain == h.policy->genesis_chain && base.revision == 0 &&
                base.bytes == 0 && base.steps == 0 &&
                base.previous_dt_s == 0 &&
                boundary_sha(base.boundary) ==
                    boundary_sha(h.policy->genesis_boundary) &&
                base.completions == ArenaJson::object() &&
                base.event_ids.empty(),
            "full-time history base equals genesis");
  else {
    const auto &a = anchors.back();
    require(base.chain == a.tail_history_chain &&
                base.revision == a.tail_revision &&
                base.bytes == a.cumulative_bytes && base.steps == a.last_step &&
                base.previous_dt_s == a.previous_dt_s &&
                boundary_sha(base.boundary) == boundary_sha(a.boundary) &&
                completions_digest(base.completions) == a.completions_sha256 &&
                id_set_digest(base.event_ids) == a.event_ids_sha256,
            "full-time history base equals the last anchor commitment");
  }
  (void)window_validate(*accepted_, channels_, carrier_cap_);
  (void)full_time_census(*genesis_, genesis_native_, carrier_cap_, false);
  validate_population(*accepted_);
  require(accepted_native_digest_ &&
              accepted_native_digest_->root ==
                  cn::carrier_rows_digest(pic_->particles_)->root,
          "accepted native row digest cache");
  require(full_time_boundary_->accepted_step == accepted_->accepted_steps &&
              full_time_boundary_->accepted_arena_revision ==
                  accepted_->revision &&
              full_time_boundary_->time_s == pic_->ledger_.time_s &&
              full_time_boundary_->accepted_arena_sha256 ==
                  FullTimeAccess::prefix(*this, *accepted_, capture(*pic_),
                                         *accepted_native_digest_) &&
              !pic_->fields_current_,
          "full-time accepted boundary/current native fingerprint");
  if (accepted_->accepted_steps == 0)
    require(full_time_boundary_->authored_full_time_genesis &&
                h.window_rows() == 0,
            "full-time empty journal genesis");
  else
    require(h.window_rows() == 0 ||
                h.tail->row.kind == RowKind::boundary,
            "full-time only complete outer checkpoint");
}
namespace {
ArenaJson metadata_v5(const ArenaSnapshot &s) {
  ArenaJson p = ArenaJson::array();
  for (const auto &[id, x] : s.packets) {
    (void)id;
    p.push_back(packet_json(x));
  }
  ArenaJson out = {{"packets", p},
                   {"ever_packet_ids", s.ever_packet_ids.get()},
                   {"event_ids", s.event_ids.get()},
                   {"journal", ArenaJson::array()},
                   {"sinks", ArenaJson::array()},
                   {"chain", s.chain},
                   {"revision", s.revision},
                   {"accepted_steps", s.accepted_steps},
                   {"previous_dt_s", s.previous_dt_s},
                   {"stochastic_completions", s.stochastic_completions},
                   {"journal_bytes", s.journal_bytes},
                   {"full_time_context", s.full_time_context}};
  if (s.full_time_flows)
    out["flows"] = flows_json(*s.full_time_flows);
  if (s.full_time_history) {
    const auto &h = *s.full_time_history;
    ArenaJson window = ArenaJson::array();
    for (const auto *r : window_rows(h))
      window.push_back({{"kind", row_kind_name(r->kind)},
                        {"sha256", r->sha256},
                        {"canonical_hex", hex_bytes(r->bytes)}});
    ArenaJson anchors = ArenaJson::array();
    for (const auto &a : h.anchors->anchors)
      anchors.push_back(
          {{"sha256", a->sha256}, {"canonical_hex", hex_bytes(a->bytes)}});
    const auto &b = *h.base;
    out["history"] = {
        {"policy_sha256", h.policy->digest},
        {"anchor_root_sha256", h.policy->root},
        {"anchors", anchors},
        {"base",
         {{"chain", b.chain},
          {"revision", b.revision},
          {"bytes", b.bytes},
          {"steps", b.steps},
          {"previous_dt_s", b.previous_dt_s},
          {"boundary", boundary_json(b.boundary)},
          {"completions", b.completions},
          {"event_ids", b.event_ids.get()}}},
        {"window", window}};
  }
  return out;
}
} // namespace
void PicArena::save_full_time(const std::string &path) const {
  require(full_time_boundary_ && !active_scope_ && owned_pic_ &&
              accepted_->full_time_context.at("mode") ==
                  "accepted_total_field_v1",
          "only accepted complete full-time boundary may be saved");
  validate_full_time_history();
  const auto &c = *configuration_;
  ArenaJson body = {
      {"schema", "canonical_native_arena_full_time_exact_v5"},
      {"source_sha256", source_identity()},
      {"registry", registry()},
      {"registry_sha256", registry_v5_sha_},
      {"initial_coupling_prefix", initial_prefix_},
      {"options", options_json(pic_->options_)},
      {"ABI", {80, 168, 200, sizeof(GapPic)}},
      {"configuration",
       {{"maximum_carriers", c.limits.maximum_carriers},
        {"memory_budget_bytes", c.limits.memory_budget_bytes},
        {"maximum_identity_entries", c.limits.maximum_identity_entries}}},
      {"history_policy",
       {{"anchor_interval_steps", c.history.anchor_interval_steps},
        {"maximum_anchors", c.history.maximum_anchors},
        {"maximum_anchor_bytes", c.history.maximum_anchor_bytes},
        {"maximum_window_rows", c.history.maximum_window_rows},
        {"maximum_window_bytes", c.history.maximum_window_bytes}}},
      {"genesis", metadata_v5(*genesis_)},
      {"genesis_native", native_rows_json(genesis_native_)},
      {"genesis_next_id", genesis_next_id_},
      {"metadata", metadata_v5(*accepted_)},
      {"native", native_rows_json(pic_->particles_)},
      {"next_id", pic_->next_particle_id_},
      {"boundary", boundary_json(*full_time_boundary_)}};
  require(body.dump().size() <= 64 * 1024 * 1024 - 1024,
          "bounded full-time output checkpoint bytes");
  atomic_json({{"payload", body}, {"sha256", sha_bytes(body.dump())}}, path);
}
namespace {
FullTimeConfiguration read_configuration(const ArenaJson &limits,
                                         const ArenaJson &history) {
  require(limits.is_object() && limits.size() == 3 && history.is_object() &&
              history.size() == 5,
          "typed configuration schema");
  FullTimeConfiguration c;
  c.limits.maximum_carriers =
      json_size(limits.at("maximum_carriers"), full_time_absolute_carrier_cap);
  c.limits.memory_budget_bytes = json_u64(limits.at("memory_budget_bytes"));
  c.limits.maximum_identity_entries =
      json_size(limits.at("maximum_identity_entries"), 1u << 30);
  c.history.anchor_interval_steps =
      json_u64(history.at("anchor_interval_steps"));
  c.history.maximum_anchors = json_size(history.at("maximum_anchors"), 1u << 30);
  c.history.maximum_anchor_bytes =
      json_size(history.at("maximum_anchor_bytes"), std::size_t(1) << 40);
  c.history.maximum_window_rows =
      json_size(history.at("maximum_window_rows"), std::size_t(1) << 40);
  c.history.maximum_window_bytes =
      json_size(history.at("maximum_window_bytes"), std::size_t(1) << 40);
  configuration_guard(c);
  return c;
}
ArenaSnapshot read_metadata_v5(const ArenaJson &j, std::size_t cap,
                               std::size_t identities) {
  json_depth_guard(j);
  require(j.is_object() &&
              (j.size() == 13 || (j.size() == 14 && j.contains("flows"))) &&
              j.contains("history"),
          "v5 metadata schema");
  require(j.at("packets").is_array() && j.at("packets").size() <= cap &&
              j.at("journal") == ArenaJson::array() &&
              j.at("sinks") == ArenaJson::array(),
          "bounded v5 metadata collections");
  require(j.at("stochastic_completions").is_object() &&
              j.at("stochastic_completions").size() <= 128,
          "bounded stochastic owner collection");
  std::size_t total_native_ids = 0;
  for (const auto &packet : j.at("packets")) {
    require(packet.is_object() && packet.contains("native_ids") &&
                packet.at("native_ids").is_array() &&
                packet.at("native_ids").size() <= cap - total_native_ids,
            "canonical native ID cap before copy");
    total_native_ids += packet.at("native_ids").size();
  }
  ArenaSnapshot s;
  s.full_time_context = j.at("full_time_context");
  s.ever_packet_ids = json_id_set(j.at("ever_packet_ids"), identities);
  s.event_ids = json_id_set(j.at("event_ids"), identities);
  s.chain = json_text(j.at("chain"));
  require(digest(s.chain), "v5 metadata chain digest");
  s.revision = json_u64(j.at("revision"));
  s.accepted_steps = json_u64(j.at("accepted_steps"));
  s.previous_dt_s = json_scalar(j.at("previous_dt_s"));
  s.stochastic_completions = j.at("stochastic_completions");
  s.journal_bytes = json_u64(j.at("journal_bytes"));
  if (j.contains("flows"))
    s.full_time_flows = read_flows_json(j.at("flows"));
  for (const auto &x : j.at("packets")) {
    require(x.size() == 6, "packet schema");
    require(!x.at("native_ids").empty() &&
                x.at("native_ids").size() <= ArenaCount::maximum_carriers,
            "bounded canonical native ID array");
    std::vector<std::uint64_t> ids;
    for (const auto &id : x.at("native_ids"))
      ids.push_back(json_u64(id));
    ArenaPacket p{json_text(x.at("id")),     json_text(x.at("state")),
                  json_text(x.at("origin")), json_text(x.at("last_event")),
                  read_count(x.at("count")), std::move(ids)};
    require(s.packets.emplace(p.id, p).second, "duplicate logical packet");
  }
  return s;
}
} // namespace
std::unique_ptr<PicArena>
PicArena::restore_full_time(const std::string &path, GapPicOptions options,
                            std::vector<ArenaStateSpec> states,
                            std::vector<ArenaChannel> channels,
                            const FullTimeAcceptedBoundary &expected) {
  const auto doc = read_full_json(path);
  require(doc.is_object() && doc.size() == 2, "full-time restart envelope");
  const auto &b = doc.at("payload");
  require(b.is_object() && b.size() == 16 &&
              b.at("schema") == "canonical_native_arena_full_time_exact_v5" &&
              b.at("source_sha256") == source_identity() &&
              doc.at("sha256") == sha_bytes(b.dump()),
          "full-time restart identity/integrity");
  require(b.at("ABI") == ArenaJson::array({80, 168, 200, sizeof(GapPic)}),
          "full-time restart ABI");
  const auto saved = read_options(b.at("options"));
  options_guard(options);
  require(options_equal(saved, options),
          "full-time restart kinetic geometry/options");
  const auto configuration =
      read_configuration(b.at("configuration"), b.at("history_policy"));
  const auto cap = configuration.limits.maximum_carriers;
  const auto identities = configuration.limits.maximum_identity_entries;
  require(boundary_json(read_boundary(b.at("boundary"))) ==
              boundary_json(expected),
          "full-time expected accepted owner seals");
  auto genesis = read_metadata_v5(b.at("genesis"), cap, identities);
  auto rows = native_rows_read(b.at("genesis_native"));
  require(rows.size() <= cap, "full-time restart genesis carrier cap");
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
  const auto final = read_metadata_v5(b.at("metadata"), cap, identities);
  (void)json_u64(b.at("next_id"));
  (void)json_u64(b.at("genesis_next_id"));
  auto owner = create_full_time(options, std::move(states), std::move(channels),
                                std::move(authored), g, configuration);
  require(
      owner->registry() == b.at("registry") &&
          owner->registry_v5_sha_ == json_text(b.at("registry_sha256")) &&
          owner->initial_prefix_ ==
              json_text(b.at("initial_coupling_prefix")) &&
          metadata_v5(*owner->genesis_) == b.at("genesis") &&
          native_rows_json(owner->genesis_native_) == b.at("genesis_native") &&
          owner->genesis_next_id_ == json_u64(b.at("genesis_next_id")),
      "full-time restart immutable genesis/catalog");
  if (final.revision == 0 && final.accepted_steps == 0) {
    require(metadata_v5(*owner->accepted_) == b.at("metadata") &&
                native_rows_json(owner->pic_->particles_) == b.at("native") &&
                owner->pic_->next_particle_id_ == json_u64(b.at("next_id")) &&
                boundary_json(*owner->full_time_boundary_) == b.at("boundary") &&
                boundary_json(*owner->full_time_boundary_) ==
                    boundary_json(expected),
            "full-time restart final native/count/history/owner seals");
    return owner;
  }
  // Nonempty checkpoint: decode and validate the exported anchor chain and
  // unanchored window (canonical bytes, digests, grammar), then refuse.
  // Restoring a nonempty full-time owner from disk is deferred (A05/A11).
  const auto &h = b.at("metadata").at("history");
  require(h.is_object() && h.size() == 5 && h.at("anchors").is_array() &&
              h.at("window").is_array() &&
              h.at("anchors").size() <= configuration.history.maximum_anchors &&
              h.at("window").size() <= configuration.history.maximum_window_rows,
          "bounded exported history");
  const auto &live = history(*owner->accepted_);
  require(json_text(h.at("policy_sha256")) == live.policy->digest &&
              json_text(h.at("anchor_root_sha256")) == live.policy->root,
          "exported history policy/root");
  Writer chain;
  chain.tag("native-full-time-anchor-chain-v6");
  chain.text(live.policy->digest);
  chain.u64(configuration.history.anchor_interval_steps);
  chain.u64(h.at("anchors").size());
  auto genesis_boundary = live.policy->genesis_boundary;
  codec_boundary(chain, genesis_boundary, true);
  chain.text(live.policy->genesis_chain);
  std::string head = live.policy->root;
  for (const auto &a : h.at("anchors")) {
    require(a.is_object() && a.size() == 2, "exported anchor schema");
    const auto bytes = unhex(json_text_any(a.at("canonical_hex")));
    require(cn::sha256_hex(bytes) == json_text(a.at("sha256")),
            "exported anchor digest");
    chain.text(bytes);
    head = json_text(a.at("sha256"));
  }
  const auto anchors =
      verify_full_time_anchor_chain(chain.bytes, live.policy->root, head);
  const auto &eb = h.at("base");
  require(eb.is_object() && eb.size() == 8, "exported history base schema");
  auto base = std::make_shared<FullTimeHistoryBase>();
  base->chain = json_text(eb.at("chain"));
  base->revision = json_u64(eb.at("revision"));
  base->bytes = json_u64(eb.at("bytes"));
  base->steps = json_u64(eb.at("steps"));
  base->previous_dt_s = json_scalar(eb.at("previous_dt_s"));
  base->boundary = read_boundary(eb.at("boundary"));
  base->completions = eb.at("completions");
  require(base->completions.is_object() && base->completions.size() <= 128,
          "exported base completions");
  base->event_ids = json_id_set(eb.at("event_ids"), identities);
  if (anchors.empty())
    require(base->chain == live.policy->genesis_chain && base->revision == 0 &&
                base->steps == 0 && base->bytes == 0,
            "exported base equals genesis");
  else
    require(base->chain == anchors.back().tail_history_chain &&
                base->revision == anchors.back().tail_revision &&
                base->steps == anchors.back().last_step &&
                boundary_sha(base->boundary) ==
                    boundary_sha(anchors.back().boundary),
            "exported base equals last anchor");
  auto rebuilt = std::make_shared<FullTimeHistory>(live);
  rebuilt->base = base;
  for (const auto &r : h.at("window")) {
    require(r.is_object() && r.size() == 3, "exported row schema");
    FullTimeStoredRow stored;
    stored.bytes = unhex(json_text_any(r.at("canonical_hex")));
    stored.sha256 = json_text(r.at("sha256"));
    stored.kind = decode_row(stored.bytes).kind;
    require(json_text(r.at("kind")) == row_kind_name(stored.kind),
            "exported row kind");
    rebuilt->tail =
        std::make_shared<const FullTimeWindowNode>(std::move(stored), rebuilt->tail);
  }
  auto probe = std::make_shared<ArenaSnapshot>(final);
  probe->full_time_history = rebuilt;
  (void)window_validate(*probe, owner->channels_, cap);
  throw refusal(RK::restore_deferred,
                "nonempty full-time disk restore is deferred (A05/A11); "
                "exported anchor chain and window validated");
}
