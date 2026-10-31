// Native v6 births and first-hit sinks; native v7 motion record (API-007.txt).
// Included at the top of mfem_pic_full_time.inc (same TU as the arena).
// Every identity below is SHA-256 over the canonical binary encoding of
// mfem_pic_canonical.h; no JSON is hashed.
std::string full_time_birth_packet_id(std::uint64_t step,
                                      std::uint32_t stream_id,
                                      std::uint32_t index) {
  return std::string(full_time_birth_namespace) + "s" + std::to_string(step) +
         "/p" + std::to_string(stream_id) + "/k" + std::to_string(index);
}
bool full_time_reserved_packet_id(const std::string &id) noexcept {
  return id.rfind(full_time_birth_namespace, 0) == 0;
}
namespace full_time_detail {
using FK = FullTimeFlowRefusal::Kind;
[[noreturn]] void flow_refusal(FK kind, const std::string &message) {
  throw FullTimeFlowRefusal(kind, message);
}
void flow_need(bool ok, FK kind, const char *message) {
  if (!ok)
    flow_refusal(kind, message);
}
// Parses "arrival-v1/s<step>/p<stream>/k<index>" (decimal, no leading zeros,
// stream and index < 2^32) and requires the canonical spelling.
bool canonical_birth_id(const std::string &id, std::uint64_t *step = nullptr,
                        std::uint32_t *stream = nullptr,
                        std::uint32_t *index = nullptr) noexcept {
  if (!full_time_reserved_packet_id(id) || id.size() > 128)
    return false;
  std::size_t at = sizeof(full_time_birth_namespace) - 1;
  std::uint64_t parts[3] = {0, 0, 0};
  const char keys[3] = {'s', 'p', 'k'};
  for (int k = 0; k < 3; ++k) {
    if (k && (at >= id.size() || id[at++] != '/'))
      return false;
    if (at >= id.size() || id[at++] != keys[k])
      return false;
    const std::size_t first = at;
    std::uint64_t v = 0;
    while (at < id.size() && id[at] >= '0' && id[at] <= '9') {
      const std::uint64_t d = std::uint64_t(id[at] - '0');
      if (v > (std::numeric_limits<std::uint64_t>::max() - d) / 10)
        return false;
      v = 10 * v + d;
      ++at;
    }
    if (at == first || (at - first > 1 && id[first] == '0'))
      return false;
    parts[k] = v;
  }
  if (at != id.size() || parts[1] > 0xffffffffull || parts[2] > 0xffffffffull)
    return false;
  if (step)
    *step = parts[0];
  if (stream)
    *stream = std::uint32_t(parts[1]);
  if (index)
    *index = std::uint32_t(parts[2]);
  return true;
}
// ---------------------------------------------------------------------------
// Exact nonnegative dyadic sums (units 2^-1074), converted to ArenaCount
// through its canonical positive-carrier decomposition.
class ExactSum {
public:
  ExactSum() : n_(big()) {}
  void add(double x) {
    require(std::isfinite(x) && x >= 0, "exact sum operand");
    if (x == 0)
      return;
    auto w = exact(x);
    require(BN_add(n_.get(), n_.get(), w.get()) == 1, "exact sum");
  }
  void add(const ArenaCount &c) {
    auto w = c.copy_units();
    require(BN_add(n_.get(), n_.get(), w.get()) == 1, "exact sum");
  }
  ArenaCount count() const {
    auto remaining = big();
    require(BN_copy(remaining.get(), n_.get()) != nullptr, "exact sum copy");
    std::vector<double> out;
    while (!BN_is_zero(remaining.get())) {
      require(out.size() < ArenaCount::maximum_carriers,
              "bounded exact sum representation");
      const int shift = std::max(0, BN_num_bits(remaining.get()) - 53);
      auto top = big();
      require(BN_rshift(top.get(), remaining.get(), shift) == 1,
              "exact sum shift");
      const double weight =
          std::ldexp(double(BN_get_word(top.get())), shift - 1074);
      require(std::isfinite(weight) && weight > 0,
              "exact sum carrier unrepresentable");
      auto word = exact(weight);
      require(BN_sub(remaining.get(), remaining.get(), word.get()) == 1,
              "exact sum decomposition");
      out.push_back(weight);
    }
    return ArenaCount::from_positive_carriers(out);
  }

private:
  Big n_;
};
// ---------------------------------------------------------------------------
// Cumulative flows digest (part of the logical snapshot when present).
void codec_surface_totals(Writer &w, const FullTimeSinkSurfaceTotals &t) {
  w.u64(std::uint64_t(t.surface));
  w.u64(std::uint64_t(t.kind));
  w.u64(t.carriers);
  for (const auto *c : {&t.weight, &t.kinetic_energy_J, &t.positive_charge_e,
                        &t.negative_charge_e})
    codec_count(w, const_cast<ArenaCount &>(*c));
}
const std::string &flows_digest(const FullTimeFlows &f) {
  if (!f.digest.empty())
    return f.digest;
  Writer w;
  w.tag("native-full-time-flows-v6");
  w.boolean(f.has_last_birth);
  w.u64(f.last_birth_step);
  w.u64(f.last_birth_base);
  w.u64(f.last_birth_count);
  w.text(f.last_birth_plan_sha256);
  w.text(f.last_birth_token_sha256);
  w.u64(f.births);
  w.u64(f.birth_rows);
  w.u64(f.sunk_carriers);
  w.u64(f.sink_rows);
  for (const auto &t : f.surfaces)
    codec_surface_totals(w, t);
  f.digest = cn::sha256_hex(w.bytes);
  return f.digest;
}
ArenaJson count_export(const ArenaCount &c) { return count_json(c); }
ArenaJson flows_json(const FullTimeFlows &f) {
  ArenaJson surfaces = ArenaJson::array();
  for (const auto &t : f.surfaces)
    surfaces.push_back({{"surface", int(t.surface)},
                        {"kind", int(t.kind)},
                        {"carriers", t.carriers},
                        {"weight", count_export(t.weight)},
                        {"kinetic_energy_J", count_export(t.kinetic_energy_J)},
                        {"positive_charge_e", count_export(t.positive_charge_e)},
                        {"negative_charge_e", count_export(t.negative_charge_e)}});
  return {{"schema", "native-full-time-flows-v6"},
          {"has_last_birth", f.has_last_birth},
          {"last_birth_step", f.last_birth_step},
          {"last_birth_base", f.last_birth_base},
          {"last_birth_count", f.last_birth_count},
          {"last_birth_plan_sha256", f.last_birth_plan_sha256},
          {"last_birth_token_sha256", f.last_birth_token_sha256},
          {"births", f.births},
          {"birth_rows", f.birth_rows},
          {"sunk_carriers", f.sunk_carriers},
          {"sink_rows", f.sink_rows},
          {"surfaces", surfaces},
          {"sha256", flows_digest(f)}};
}
std::shared_ptr<FullTimeFlows> read_flows_json(const ArenaJson &j) {
  require(j.is_object() && j.size() == 13 &&
              j.at("schema") == "native-full-time-flows-v6" &&
              j.at("surfaces").is_array() && j.at("surfaces").size() == 5,
          "flows checkpoint schema");
  auto f = std::make_shared<FullTimeFlows>();
  f->has_last_birth = json_boolean(j.at("has_last_birth"));
  f->last_birth_step = json_u64(j.at("last_birth_step"));
  f->last_birth_base = json_u64(j.at("last_birth_base"));
  f->last_birth_count = json_u64(j.at("last_birth_count"));
  f->last_birth_plan_sha256 = json_text_any(j.at("last_birth_plan_sha256"));
  f->last_birth_token_sha256 = json_text_any(j.at("last_birth_token_sha256"));
  f->births = json_u64(j.at("births"));
  f->birth_rows = json_u64(j.at("birth_rows"));
  f->sunk_carriers = json_u64(j.at("sunk_carriers"));
  f->sink_rows = json_u64(j.at("sink_rows"));
  for (int k = 0; k < 5; ++k) {
    const auto &x = j.at("surfaces")[std::size_t(k)];
    require(x.is_object() && x.size() == 7 && json_u64(x.at("surface")) == std::uint64_t(k),
            "flows surface schema");
    auto &t = f->surfaces[std::size_t(k)];
    require(json_u64(x.at("kind")) == std::uint64_t(t.kind), "flows surface kind");
    t.carriers = json_u64(x.at("carriers"));
    t.weight = read_count(x.at("weight"));
    t.kinetic_energy_J = read_count(x.at("kinetic_energy_J"));
    t.positive_charge_e = read_count(x.at("positive_charge_e"));
    t.negative_charge_e = read_count(x.at("negative_charge_e"));
  }
  require(json_text(j.at("sha256")) == flows_digest(*f), "flows checkpoint digest");
  return f;
}
// ---------------------------------------------------------------------------
// Block Merkle root over ID-sorted records (same blocking as carriers).
template <class Row, class Encode>
std::string id_records_root(const char *leaf_tag, const char *root_tag,
                            const std::vector<Row> &rows, Encode encode) {
  std::vector<cn::BlockLeaf> leaves;
  std::size_t i = 0;
  while (i < rows.size()) {
    require(rows[i].native_id > 0, "record native ID");
    const auto block = (rows[i].native_id - 1) / cn::carrier_block;
    std::size_t j = i;
    while (j < rows.size() && (rows[j].native_id - 1) / cn::carrier_block == block)
      ++j;
    Writer w;
    w.tag(leaf_tag);
    w.u64(block);
    w.u64(j - i);
    for (std::size_t k = i; k < j; ++k) {
      require(k == i || rows[k].native_id > rows[k - 1].native_id,
              "record order");
      encode(w, rows[k]);
    }
    leaves.push_back({block, i, j - i, cn::sha256_raw(w.bytes)});
    i = j;
  }
  return cn::merkle_root(root_tag, leaves, rows.size());
}
void encode_birth_record(Writer &w, const FullTimeBirth &b) {
  w.u64(b.native_id);
  w.u64(b.step);
  w.u64(b.stream_id);
  w.u64(b.index);
  w.text(b.packet_id);
  w.text(b.state);
  w.text(b.origin);
  w.f64(b.weight);
  w.vec(b.position_m);
  w.vec(b.proper_velocity_m_s);
  w.f64(b.birth_time_s);
  w.f64(b.remaining_flight_s);
}
void encode_sink_record(Writer &w, const FullTimeSinkRecord &r) {
  w.u64(r.native_id);
  w.u64(std::uint64_t(r.kind));
  w.u64(std::uint64_t(r.surface));
  w.u64(r.axis);
  w.f64(r.plane_m);
  w.f64(r.fraction);
  w.vec(r.segment_m);
  w.vec(r.hit_point_m);
  w.vec(r.incident_proper_velocity_m_s);
  w.f64(r.flight_s);
  w.f64(r.discarded_flight_s);
}
std::string sink_ledger_digest(const FullTimeSinkLedger &l) {
  Writer w;
  w.tag("native-full-time-sink-ledger-v6");
  w.u64(l.accepted_step);
  w.text(l.transport_trial_sha256);
  w.u64(l.carriers);
  w.text(l.records_root_sha256);
  w.u64(l.removed_packets.size());
  for (const auto &p : l.removed_packets)
    w.text(p);
  w.u64(l.tallies.size());
  for (const auto &t : l.tallies) {
    w.u64(std::uint64_t(t.surface));
    w.u64(std::uint64_t(t.kind));
    w.text(t.state);
    w.i64(t.charge_e);
    w.u64(t.carriers);
    codec_count(w, const_cast<ArenaCount &>(t.weight));
    codec_count(w, const_cast<ArenaCount &>(t.kinetic_energy_J));
  }
  w.u64(l.surfaces.size());
  for (const auto &t : l.surfaces)
    codec_surface_totals(w, t);
  return cn::sha256_hex(w.bytes);
}
std::string text_list_digest(const char *tag, const std::vector<std::string> &v) {
  Writer w;
  w.tag(tag);
  w.u64(v.size());
  for (const auto &x : v)
    w.text(x);
  return cn::sha256_hex(w.bytes);
}
// ---------------------------------------------------------------------------
// Native port of the sealed Hall kinetic first-hit classification
// (hall-canonical-driver-20261029/kinetic_boundary.cpp first_kinetic_hit,
// native-finite-target-kinetic-policy-v1): same operands, same operation
// order, same tie rule (later exact tie wins; fraction 1 is a hit).
struct NativeHit {
  int sink = -1; // hall_driver::KineticSink value, -1 none
  double fraction = 1;
  int axis = -1;
  double plane = 0;
};
NativeHit native_first_hit(const GapPicOptions &o, const PicVector &x,
                           const PicVector &D) {
  for (double v : D)
    require(std::isfinite(v), "nonfinite kinetic displacement");
  NativeHit hit;
  auto take = [&](double f, int sink, int axis, double plane) {
    if (f < 0 || f > hit.fraction)
      return;
    hit = {sink, f, axis, plane};
  };
  for (int axis = 0; axis < 2; ++axis)
    if (D[axis] != 0) {
      const double plane = D[axis] > 0 ? (axis ? o.gap_m : o.width_m) : 0;
      take((plane - x[axis]) / D[axis], 2, axis, plane);
    }
  int body = 0;
  for (const auto &foil : {o.lower_foil, o.upper_foil}) {
    double enter = -std::numeric_limits<double>::infinity(),
           leave = std::numeric_limits<double>::infinity(), plane = 0;
    int enter_axis = -1;
    bool possible = true;
    for (int axis = 0; axis < 2; ++axis) {
      const double low = axis ? foil.bottom_m : foil.left_m,
                   high = axis ? foil.top_m : foil.right_m;
      if (D[axis] == 0) {
        if (x[axis] < low || x[axis] > high)
          possible = false;
        continue;
      }
      double a = (low - x[axis]) / D[axis], b = (high - x[axis]) / D[axis],
             entry_plane = low;
      if (a > b) {
        std::swap(a, b);
        entry_plane = high;
      }
      if (a >= enter) {
        enter = a;
        enter_axis = axis;
        plane = entry_plane;
      }
      leave = std::min(leave, b);
    }
    if (possible && leave > 0 && enter >= 0 && enter <= leave &&
        enter_axis >= 0)
      take(enter, body ? 1 : 0, enter_axis, plane);
    ++body;
  }
  if (D[2] != 0) {
    const double plane = D[2] > 0 ? o.feed_z_m : o.short_z_m,
                 f = (plane - x[2]) / D[2];
    const double xx = x[0] + f * D[0], yy = x[1] + f * D[1];
    const bool bridge = D[2] < 0 && xx >= o.lower_foil.left_m &&
                        xx <= o.lower_foil.right_m && yy >= o.lower_foil.top_m &&
                        yy <= o.upper_foil.bottom_m;
    take(f, bridge ? 3 : 4, 2, plane);
  }
  return hit;
}
// Birth on an emitting foil face: on exactly one face of one foil (not an
// edge or corner), strictly inside the kinetic box, proper velocity pointing
// out of the metal into vacuum. Source v5 places neutral packets there.
bool surface_birth(const GapPicOptions &o, const GapParticle &p) {
  native_particle_guard(p);
  const auto &x = p.position_m;
  const auto &u = p.proper_velocity_m_s;
  if (!o.finite_foils || !(x[0] > 0 && x[0] < o.width_m && x[1] > 0 &&
                           x[1] < o.gap_m))
    return false;
  if (o.finite_axial_extent && !(x[2] > o.short_z_m && x[2] < o.feed_z_m))
    return false;
  int faces = 0;
  bool outward = false;
  for (const auto &f : {o.lower_foil, o.upper_foil}) {
    const bool in_x = x[0] > f.left_m && x[0] < f.right_m,
               in_y = x[1] > f.bottom_m && x[1] < f.top_m;
    const bool on_x = x[0] == f.left_m || x[0] == f.right_m,
               on_y = x[1] == f.bottom_m || x[1] == f.top_m;
    if (in_x && in_y)
      return false; // inside metal
    if (in_x && on_y) {
      ++faces;
      outward = x[1] == f.bottom_m ? u[1] < 0 : u[1] > 0;
    } else if (in_y && on_x) {
      ++faces;
      outward = x[0] == f.left_m ? u[0] < 0 : u[0] > 0;
    } else if (on_x && on_y &&
               (x[0] == f.left_m || x[0] == f.right_m) &&
               (x[1] == f.bottom_m || x[1] == f.top_m))
      return false; // foil corner: source/boundary owner
  }
  return faces == 1 && outward;
}
FullTimeSinkKind sink_kind_of(int surface) {
  return (surface == 0 || surface == 1 || surface == 3)
             ? FullTimeSinkKind::collected_at_conductor
             : FullTimeSinkKind::escaped_exterior;
}
// ---------------------------------------------------------------------------
// Streaming canonical encoder into SHA-256 (bounded buffer, no per-record
// allocation after the first flush).
class StreamDigest {
public:
  void raw(const void *p, std::size_t n) {
    h_.w.raw(p, n);
    h_.flush_if_large(1u << 16);
  }
  void u64(std::uint64_t v) { h_.w.u64(v); }
  void f64(double v) { h_.w.f64(v); }
  void text(std::string_view s) {
    h_.w.u64(s.size());
    h_.w.raw(s.data(), s.size());
  }
  void vec(const PicVector &v) { h_.w.vec(v); }
  std::string hex() { return cn::hex(h_.finish()); }

private:
  cn::Hasher h_;
};
static_assert(sizeof(FullTimeMotionPath) == 64 &&
              std::is_trivially_copyable_v<FullTimeMotionPath> &&
              offsetof(FullTimeMotionPath, field_start_m) == 8 &&
              offsetof(FullTimeMotionPath, field_displacement_m) == 32 &&
              offsetof(FullTimeMotionPath, flight_s) == 56);
} // namespace full_time_detail
FullTimeKineticHit full_time_first_kinetic_hit(const GapPicOptions &o,
                                               const PicVector &x,
                                               const PicVector &D) {
  const auto h = full_time_detail::native_first_hit(o, x, D);
  return {h.sink, h.fraction, h.axis, h.plane};
}
namespace full_time_detail {
constexpr const char *flows_coverage =
    "all_native_carriers_births_and_first_hit_sinks";
} // namespace full_time_detail
namespace full_time_detail {
static_assert(sizeof(GapParticle) == sizeof(pf_v7::Carrier) &&
              sizeof(FullTimeMotionPath) == sizeof(pf_v7::Path) &&
              sizeof(FullTimeBirthFlight) == sizeof(pf_v7::Birth) &&
              offsetof(FullTimeBirthFlight, birth_time_s) == 8 &&
              offsetof(FullTimeBirthFlight, flight_s) == 16 &&
              std::is_trivially_copyable_v<FullTimeBirthFlight>);
pf_v7::MotionHeaderV7 v7_header(const FullTimeMotionRecordV7 &r) {
  pf_v7::MotionHeaderV7 h;
  h.accepted_step = r.accepted_step;
  h.staged_arena_revision = r.staged_arena_revision;
  h.trial_stage = r.trial_stage;
  h.begin_s = r.begin_s;
  h.end_s = r.end_s;
  h.dt_s = r.dt_s;
  h.staged_arena_sha256 = r.staged_arena_sha256;
  h.trial_field_sha256 = r.trial_field_sha256;
  h.mesh_sha256 = r.mesh_sha256;
  h.material_sha256 = r.material_sha256;
  h.species_sha256 = r.species_sha256;
  h.common_clock_sha256 = r.common_clock_sha256;
  h.kinetic_policy_sha256 = r.kinetic_policy_sha256;
  h.budget = r.budget;
  h.census_sha256 = r.census_sha256;
  h.coverage = r.coverage;
  h.current_A_m_sha256 = r.current_A_m_sha256;
  h.endpoint_charge_C_sha256 = r.endpoint_charge_C_sha256;
  h.initial_charge_C_sha256 = r.initial_charge_C_sha256;
  h.kinetic_policy = r.kinetic_policy;
  h.native_to_field_m = r.native_to_field_m;
  h.neutral_motion_owner = r.neutral_motion_owner;
  h.native_DKD_called = r.native_DKD_called;
  h.native_Poisson_called = r.native_Poisson_called;
  h.neutral_stream_dt_max_s = r.neutral_stream_dt_max_s;
  return h;
}
// Canonical 136-byte sink records (enums widened to u64).
std::vector<pf_v7::Sink> v7_sinks(std::span<const FullTimeSinkRecord> sinks) {
  std::vector<pf_v7::Sink> out(sinks.size());
  for (std::size_t i = 0; i < sinks.size(); ++i) {
    const auto &k = sinks[i];
    auto &c = out[i];
    c.id = k.native_id;
    c.kind = std::uint64_t(k.kind);
    c.surface = std::uint64_t(k.surface);
    c.axis = k.axis;
    c.plane_m = k.plane_m;
    c.fraction = k.fraction;
    for (int d = 0; d < 3; ++d) {
      c.segment_m[d] = k.segment_m[d];
      c.hit_point_m[d] = k.hit_point_m[d];
      c.incident_proper_velocity_m_s[d] = k.incident_proper_velocity_m_s[d];
    }
    c.flight_s = k.flight_s;
    c.discarded_flight_s = k.discarded_flight_s;
  }
  return out;
}
pf_v7::MotionSectionsV7 v7_sections(const FullTimeMotionRecordV7 &r) {
  const auto sinks = v7_sinks(r.sinks);
  return pf_v7::motion_sections_v7(
      {reinterpret_cast<const pf_v7::Carrier *>(r.endpoints.data()), r.endpoints.size()},
      {reinterpret_cast<const pf_v7::Path *>(r.paths.data()), r.paths.size()}, sinks,
      {reinterpret_cast<const pf_v7::Birth *>(r.births.data()), r.births.size()});
}
} // namespace full_time_detail
// Canonical digest of a motion record v7 (API-007.txt 1.2): header with the
// four section block-list roots.
std::string full_time_motion_record_v7_sha256(const FullTimeMotionRecordV7 &r) {
  using namespace full_time_detail;
  require(r.schema == motion_v7_schema, "motion record v7 schema");
  return pf_v7::motion_record_v7_sha256(v7_header(r), v7_sections(r));
}
FullTimeMotionAdmission
bind_total_field_motion_v7(const FullTimeMotionRecordV7 &record,
                           std::shared_ptr<const void> storage) {
  FullTimeMotionAdmission a;
  a.v7_ = std::make_shared<const FullTimeMotionRecordV7>(record);
  a.v7_storage_ = std::move(storage);
  return a;
}
