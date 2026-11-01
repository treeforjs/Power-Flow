// Native v8 chemistry v3 integration (API-008.txt 3): configuration and
// registration, the chemistry stage of an exposure (FullTimeAccess::chemistry)
// and the accepted chemistry status. Included in mfem_pic_full_time.inc after
// the FullTimeAccess class (same TU as the arena).
namespace full_time_detail {
// v8 restart state of the published sinks (source v7 B1/D8): bit-exact JSON.
std::string bits_hex(double x) {
  char b[20];
  std::snprintf(b, sizeof b, "%016llx", (unsigned long long)std::bit_cast<std::uint64_t>(x));
  return b;
}
double hex_bits(const ArenaJson &j) {
  const auto t = json_text(j);
  require(t.size() == 16 && t.find_first_not_of("0123456789abcdef") == std::string::npos, "binary64 bit text");
  return std::bit_cast<double>(std::uint64_t(std::stoull(t, nullptr, 16)));
}
ArenaJson step_sinks_json(const FullTimeSinkLedger &l) {
  ArenaJson ex = ArenaJson::array();
  for (const auto &e : l.exposed) {
    const auto &r = e.record;
    ArenaJson v = ArenaJson::array();
    for (double x : {r.plane_m, r.fraction, r.segment_m[0], r.segment_m[1], r.segment_m[2], r.hit_point_m[0],
                     r.hit_point_m[1], r.hit_point_m[2], r.incident_proper_velocity_m_s[0],
                     r.incident_proper_velocity_m_s[1], r.incident_proper_velocity_m_s[2], r.flight_s,
                     r.discarded_flight_s, e.weight, e.charge_C, e.mass_kg, e.start_position_m[0],
                     e.start_position_m[1], e.start_position_m[2], e.internal_energy_J})
      v.push_back(bits_hex(x));
    ex.push_back({{"native_id", r.native_id}, {"kind", int(r.kind)}, {"surface", int(r.surface)}, {"axis", r.axis},
                  {"packet_id", e.packet_id}, {"state", e.state}, {"values", v}});
  }
  return {{"accepted_step", l.accepted_step}, {"records_root_sha256", l.records_root_sha256},
          {"transport_trial_sha256", l.transport_trial_sha256}, {"ledger_sha256", l.ledger_sha256},
          {"row_sha256", l.row_sha256}, {"exposed", ex}};
}
std::shared_ptr<FullTimeSinkLedger> read_step_sinks(const ArenaJson &j) {
  require(j.is_object() && j.size() == 6 && j.at("exposed").is_array(), "published sinks schema");
  auto l = std::make_shared<FullTimeSinkLedger>();
  l->accepted_step = json_u64(j.at("accepted_step"));
  l->records_root_sha256 = json_text(j.at("records_root_sha256"));
  l->transport_trial_sha256 = json_text(j.at("transport_trial_sha256"));
  l->ledger_sha256 = json_text(j.at("ledger_sha256"));
  l->row_sha256 = json_text(j.at("row_sha256"));
  for (const auto &x : j.at("exposed")) {
    require(x.is_object() && x.size() == 7 && x.at("values").is_array() && x.at("values").size() == 20,
            "published sink schema");
    FullTimeSinkExposure e;
    auto &r = e.record;
    r.native_id = json_u64(x.at("native_id"));
    r.kind = FullTimeSinkKind(std::uint8_t(json_u64(x.at("kind"))));
    r.surface = FullTimeSinkSurface(std::uint8_t(json_u64(x.at("surface"))));
    r.axis = std::uint32_t(json_u64(x.at("axis")));
    std::array<double, 20> v;
    for (int k = 0; k < 20; ++k)
      v[k] = hex_bits(x.at("values")[k]);
    r.plane_m = v[0];
    r.fraction = v[1];
    r.segment_m = {v[2], v[3], v[4]};
    r.hit_point_m = {v[5], v[6], v[7]};
    r.incident_proper_velocity_m_s = {v[8], v[9], v[10]};
    r.flight_s = v[11];
    r.discarded_flight_s = v[12];
    e.weight = v[13];
    e.charge_C = v[14];
    e.mass_kg = v[15];
    e.start_position_m = {v[16], v[17], v[18]};
    e.internal_energy_J = v[19];
    e.packet_id = json_text(x.at("packet_id"));
    e.state = json_text(x.at("state"));
    l->records.push_back(r);
    l->exposed.push_back(e);
  }
  l->carriers = l->records.size();
  require(id_records_root("native-v6/sink-block", "native-v6/sink-root", l->records, encode_sink_record) ==
              l->records_root_sha256,
          "published sink records differ from their records root");
  return l;
}
ArenaJson chemistry_configuration_json(const FullTimeConfiguration &c) {
  ArenaJson j = {{"identity", "v8.2-auditfix2"}};
  if (c.merge_split.enabled) {
    const auto &P = *c.merge_split.params;
    j["merge_split"] = {{"cadence_steps", c.merge_split.cadence_steps},
                        {"species_states", c.merge_split.species_states},
                        {"params_hex", hex_bytes(std::string(reinterpret_cast<const char *>(&P), sizeof P))},
                        {"policy_sha256", merge_split::policy_sha256(P)}};
  }
  const auto &k = c.chemistry;
  if (!k.enabled) {
    j["chemistry"] = nullptr;
    return j;
  }
  j["chemistry"] = {{"diagnostic_field_target_approximation", k.diagnostic_field_target_approximation},
                    {"blob_path", k.blob_path},
                    {"blob_sha256", k.blob_sha256},
                    {"run_seed", k.run_seed},
                    {"release_weight", k.release_weight},
                    {"trunc_sigmas", k.trunc_sigmas},
                    {"max_candidates", k.max_candidates},
                    {"max_events", k.max_events},
                    {"grid_counts", k.grid_counts},
                    {"preset", k.preset},
                    {"policy",
                     {{"event_limit", k.policy.event_limit},
                      {"event_fraction", k.policy.event_fraction},
                      {"admitted_regions", k.policy.admitted_regions},
                      {"policy_file", k.policy.policy_file},
                      {"policy_file_sha256", k.policy.policy_file_sha256}}}};
  return j;
}
void read_chemistry_configuration_json(const ArenaJson &j, FullTimeConfiguration &c) {
  const bool ms = j.is_object() && j.contains("merge_split");
  require(j.is_object() && j.size() == (ms ? 3u : 2u) && j.at("identity") == "v8.2-auditfix2", "v8 configuration schema");
  c.identity = FullTimeIdentitySchema::v8;
  if (ms) {
    const auto &m = j.at("merge_split");
    require(m.is_object() && m.size() == 4 && m.at("species_states").is_array(), "merge/split configuration schema");
    const auto bytes = unhex(json_text_any(m.at("params_hex")));
    require(bytes.size() == sizeof(MsParams), "merge/split policy size");
    auto P = std::make_shared<MsParams>();
    std::memcpy(P.get(), bytes.data(), sizeof(MsParams));
    require(merge_split::policy_sha256(*P) == json_text(m.at("policy_sha256")), "merge/split policy digest");
    c.merge_split.enabled = true;
    c.merge_split.cadence_steps = json_u64(m.at("cadence_steps"));
    for (const auto &st : m.at("species_states"))
      c.merge_split.species_states.push_back(json_text(st));
    c.merge_split.params = std::move(P);
  }
  const auto &k = j.at("chemistry");
  if (k.is_null())
    return;
  require(k.is_object() && k.size() == 11, "chemistry configuration schema");
  auto &x = c.chemistry;
  x.enabled = true;
  require(k.at("diagnostic_field_target_approximation").is_boolean(), "target energy mode boolean");
  x.diagnostic_field_target_approximation = k.at("diagnostic_field_target_approximation").get<bool>();
  x.blob_path = json_text(k.at("blob_path"));
  x.blob_sha256 = json_text(k.at("blob_sha256"));
  x.run_seed = json_u64(k.at("run_seed"));
  x.release_weight = json_scalar(k.at("release_weight"));
  x.trunc_sigmas = json_scalar(k.at("trunc_sigmas"));
  x.max_candidates = std::uint32_t(json_u64(k.at("max_candidates")));
  x.max_events = std::uint32_t(json_u64(k.at("max_events")));
  x.preset = json_text_any(k.at("preset"));
  require(k.at("grid_counts").is_array() && k.at("grid_counts").size() == 3, "chemistry grid counts");
  for (int d = 0; d < 3; ++d)
    x.grid_counts[d] = std::uint32_t(json_u64(k.at("grid_counts")[d]));
  const auto &p = k.at("policy");
  require(p.is_object() && p.size() == 5 && p.at("admitted_regions").is_array(), "chemistry policy schema");
  x.policy.policy_file = json_text_any(p.at("policy_file"));
  x.policy.policy_file_sha256 = json_text_any(p.at("policy_file_sha256"));
  x.policy.event_limit = json_scalar(p.at("event_limit"));
  x.policy.event_fraction = json_scalar(p.at("event_fraction"));
  for (const auto &r : p.at("admitted_regions"))
    // Region labels come from the SHA-bound table and may exceed packet-ID
    // lengths. make_setup validates exact membership in that registry.
    x.policy.admitted_regions.push_back(json_text_any(r));
}
// Identity extension of the configuration digest: empty in schema v5 (v7.1
// identities unchanged); in schema v8 the chemistry configuration identity.
std::string v8_configuration_extension(const FullTimeConfiguration &c, const GapPicOptions &o) {
  if (c.identity != FullTimeIdentitySchema::v8)
    return "";
  Writer w;
  w.tag("native-v8.2-auditfix2/identity-schema-v8");
  if (c.chemistry.enabled)
    w.text(chem8::make_setup(c.chemistry, o)->configuration_sha256);
  else
    w.text("no-chemistry");
  // v8.1: merge/split policy (absent when disabled: v8 identities unchanged)
  if (c.merge_split.enabled) {
    w.text("native-v8.1/merge-split");
    w.text(merge_split::policy_sha256(*c.merge_split.params));
    w.u64(c.merge_split.cadence_steps);
    w.u64(c.merge_split.species_states.size());
    for (const auto &st : c.merge_split.species_states)
      w.text(st);
  }
  return cn::sha256_hex(w.bytes);
}
void chemistry_registration_guard(const chem8::Setup &s, const std::map<std::string, ArenaStateSpec> &states,
                                  const std::map<std::string, ArenaChannel> &channels) {
  using CK = FullTimeChemistryRefusal::Kind;
  const auto &r = *s.reg;
  for (std::size_t i = 0; i < r.states.size(); ++i) {
    auto it = states.find(r.states[i]);
    if (it == states.end() || it->second.charge_e != r.charge[i] ||
        std::bit_cast<std::uint64_t>(it->second.mass_kg) != std::bit_cast<std::uint64_t>(r.mass_kg[i]) ||
        it->second.hydrogen != r.nuclei[i][0] || it->second.oxygen != r.nuclei[i][1])
      throw FullTimeChemistryRefusal(CK::configuration,
                                     "chemistry_v3/configuration: species " + r.states[i] +
                                         " not registered with the blob mass/charge/nuclei");
  }
  auto ch = channels.find(chem8::channel_id);
  const auto &w = s.channel;
  if (ch == channels.end() || ch->second.primary_state != w.primary_state ||
      ch->second.target_state != w.target_state || ch->second.product_states != w.product_states ||
      ch->second.dataset_sha256 != w.dataset_sha256 || ch->second.kernel_manifest_sha256 != w.kernel_manifest_sha256 ||
      ch->second.kernel_binding_file_sha256 != w.kernel_binding_file_sha256 || ch->second.model != w.model)
    throw FullTimeChemistryRefusal(CK::configuration,
                                   "chemistry_v3/configuration: channel pfc-chem-v3 not registered with the "
                                   "blob/registry/configuration binding");
}
std::string births_root(const std::vector<std::tuple<std::uint64_t, std::string, double>> &born) {
  Writer w;
  w.tag("native-v8/chemistry-births");
  w.u64(born.size());
  for (const auto &[id, packet, weight] : born) {
    w.u64(id);
    w.text(packet);
    w.f64(weight);
  }
  return cn::sha256_hex(w.bytes);
}
} // namespace full_time_detail

struct FullTimeAccess::ChemOutcome {
  std::shared_ptr<ArenaSnapshot> next;
  PicArena::NativeState native;
  std::shared_ptr<const cn::BlockDigest> digest; // device root (null: host rows)
  std::shared_ptr<FullTimeChemistryStageReceipt> receipt;
  std::string inputs_sha256;
  std::uint64_t purpose = 0;
};
namespace full_time_detail {
// Census state table of an arena (registered states, ascending = state index)
// and the registry species of each state index.
std::vector<std::uint32_t> state_species_of(const std::map<std::string, ArenaStateSpec> &states,
                                            const chem8::Setup &setup) {
  std::vector<std::uint32_t> out;
  for (const auto &[id, spec] : states) {
    (void)spec;
    auto it = std::find(setup.reg->states.begin(), setup.reg->states.end(), id);
    out.push_back(it == setup.reg->states.end() ? chem8::none : std::uint32_t(it - setup.reg->states.begin()));
  }
  return out;
}
// Thrown by the device outcome for a stage it does not cover (multi-carrier
// target packets, non-binary64 remainders, refusing statuses): the host path runs.
struct ChemDeviceFallback {};
} // namespace full_time_detail
std::shared_ptr<FullTimeAccess::ChemOutcome>
FullTimeAccess::chemistry_host(FullCollision &s, const chem8::Setup &setup_ref, const chem8::StageRequest &req) {
  auto &o = *s.outer;
  const chem8::Setup *setup = &setup_ref;
  const auto &base = o.frame;
  const auto &snap = *base.snapshot;
  require(snap.full_time_identity == 1 && bool(snap.full_time_chemistry), "frame without chemistry state");
  const auto &before = *snap.full_time_chemistry;
  // ---- census view of the private frame (rows ascending native ID)
  const auto &rows = base.native.particles.host();
  const std::size_t n = rows.size();
  std::vector<std::uint32_t> species(n, chem8::none), first(n, chem8::none);
  std::vector<std::uint8_t> single(n, 0);
  std::vector<const ArenaPacket *> owner_of(n, nullptr);
  std::map<std::string, std::uint32_t> index;
  for (std::uint32_t i = 0; i < setup->reg->states.size(); ++i)
    index.emplace(setup->reg->states[i], i);
  auto row_of = [&](std::uint64_t id) {
    const auto it = std::lower_bound(rows.begin(), rows.end(), id,
                                     [](const GapParticle &a, std::uint64_t b) { return a.id < b; });
    require(it != rows.end() && it->id == id, "native canonical carrier missing");
    return std::uint32_t(it - rows.begin());
  };
  for (const auto &[id, p] : snap.packets) {
    (void)id;
    const auto it = index.find(p.state);
    const auto sp = it == index.end() ? chem8::none : it->second;
    std::uint32_t f = chem8::none;
    for (auto nid : p.native_ids) {
      const auto r = row_of(nid);
      if (f == chem8::none)
        f = r;
      species[r] = sp;
      single[r] = p.native_ids.size() == 1;
      first[r] = f;
      owner_of[r] = &p;
    }
  }
  chem8::RowView view{rows.data(), n, species.data(), single.data()};
  auto in = chem8::build_inputs(*setup, view, req);
  auto rec = chem8::run_host(*setup, in, setup->config.host_partitions);
  const auto active = chem8::active_of(*setup, in, rec, view);
  auto count_of = [&](std::uint32_t r) {
    chem8::ESum e;
    for (double w : PicArena::positive_carriers(owner_of[r]->count))
      e.add(w);
    return e;
  };
  chem8::DebitAccess da;
  da.inventory = [&](std::uint32_t cell, std::uint32_t sp) {
    const auto it = in.inventory.find(pfc::InventoryKey{cell, sp});
    return it == in.inventory.end() ? chem8::ESum() : it->second;
  };
  da.packets = [&](std::uint32_t atlas, std::uint32_t sp) {
    std::vector<std::pair<std::uint32_t, chem8::ESum>> out;
    for (std::uint32_t r = 0; r < n; ++r)
      if (species[r] == sp && in.row_cell[r] == atlas && first[r] == r)
        out.emplace_back(r, count_of(r));
    return out;
  };
  auto res = chem8::resolve(*setup, active, req, before, da);
  const auto &plan = res.plan;
  // ---- apply the plan to a new private frame
  auto next = std::make_shared<ArenaSnapshot>(snap);
  next->full_time_chemistry = res.state;
  next->digest_cache.clear();
  PicArena::NativeState native = base.native;
  std::uint64_t next_id = native.next_id;
  const std::uint64_t base_id = next_id;
  const std::string stage_tag = std::string(chem8::packet_namespace) + "s" + std::to_string(req.step) + "/g" +
                                std::to_string(req.stage);
  std::vector<std::uint8_t> drop(n, 0);
  std::vector<std::int64_t> vel(n, -1);
  for (std::size_t k = 0; k < plan.velocity.size(); ++k)
    vel[plan.velocity[k].first] = std::int64_t(k);
  for (auto r : plan.consumed_rows) {
    drop[r] = 1;
    require(next->packets.erase(owner_of[r]->id) == 1 && next->full_time_packet_index.erase(rows[r].id) == 1,
            "chemistry consumed packet");
  }
  struct Remainder {
    std::uint32_t row;
    std::vector<double> weights;
    ArenaCount count;
  };
  std::vector<Remainder> remainders;
  for (const auto &d : plan.debits) {
    const auto &p = *owner_of[d.row];
    for (auto nid : p.native_ids)
      drop[row_of(nid)] = 1;
    if (d.whole) {
      require(next->packets.erase(p.id) == 1 &&
                  next->full_time_packet_index.erase(*std::min_element(p.native_ids.begin(), p.native_ids.end())) == 1,
              "chemistry debited packet");
      continue;
    }
    // exact remainder count = count - amount (credits first, then debits)
    ArenaCount c = p.count;
    std::vector<double> pieces;
    chem8::ESum left = d.amount;
    while (!left.is_zero()) {
      const double h = left.approx();
      pieces.push_back(h);
      chem8::ESum t;
      t.add(h);
      left.sub(t);
    }
    for (double x : pieces)
      if (x < 0)
        c = c.credited(-x);
    for (double x : pieces)
      if (x > 0)
        c = c.debit(x);
    require(!c.zero(), "chemistry partial debit leaves a positive count");
    remainders.push_back({d.row, PicArena::positive_carriers(c), c});
  }
  const auto &reg = *setup->reg;
  auto carrier_of = [&](std::uint32_t sp, const PicVector &x, const double *u, double w) {
    GapParticle g;
    g.id = 0;
    g.position_m = x;
    g.proper_velocity_m_s = {u[0], u[1], u[2]};
    g.charge_C = reg.charge[sp] * pic_e;
    g.mass_kg = reg.mass_kg[sp];
    g.weight = w;
    return g;
  };
  auto issue = [&]() {
    require(next_id < std::numeric_limits<std::uint64_t>::max(), "native ID exhaustion");
    return next_id++;
  };
  std::vector<GapParticle> born;
  std::vector<std::tuple<std::uint64_t, std::string, double>> born_list;
  for (const auto &b : plan.births) {
    auto g = carrier_of(b.species, rows[b.row].position_m, b.u, b.weight);
    g.id = issue();
    closed_position(native.options, g);
    const auto pid = chem8::packet_id_birth(req.step, req.stage, b);
    require(next->packets.emplace(pid, ArenaPacket{pid, reg.states[b.species],
                                                   std::string(chem8::packet_namespace) + "p" + std::to_string(b.process),
                                                   "", ArenaCount(b.weight), {g.id}})
                .second,
            "chemistry birth packet collision");
    require(next->full_time_packet_index.emplace(g.id, pid).second, "chemistry birth key collision");
    born.push_back(g);
    born_list.emplace_back(g.id, pid, b.weight);
  }
  std::uint64_t remainder_carriers = 0;
  for (auto &rm : remainders) {
    const auto &p = *owner_of[rm.row];
    auto &pp = next->packets.modify(p.id);
    pp.count = rm.count;
    pp.last_event = stage_tag;
    pp.native_ids = {p.native_ids.front()};
    for (std::size_t k = 1; k < rm.weights.size(); ++k) {
      GapParticle g = rows[rm.row];
      if (vel[rm.row] >= 0)
        g.proper_velocity_m_s = plan.velocity[vel[rm.row]].second;
      g.weight = rm.weights[k];
      g.id = issue();
      pp.native_ids.push_back(g.id);
      born.push_back(g);
      ++remainder_carriers;
    }
  }
  for (const auto &rl : plan.releases) {
    auto g = carrier_of(rl.species, PicVector{rl.x[0], rl.x[1], rl.x[2]}, rl.u, rl.weight);
    g.id = issue();
    closed_position(native.options, g);
    const auto pid = chem8::packet_id_release(req.step, rl);
    require(next->packets.emplace(pid, ArenaPacket{pid, reg.states[rl.species],
                                                   std::string(chem8::packet_namespace) + "release", "",
                                                   ArenaCount(rl.weight), {g.id}})
                .second,
            "chemistry release packet collision");
    require(next->full_time_packet_index.emplace(g.id, pid).second, "chemistry release key collision");
    born.push_back(g);
    born_list.emplace_back(g.id, pid, rl.weight);
  }
  // remainder first rows keep their native ID
  std::vector<std::int64_t> remainder_of(n, -1);
  for (std::size_t k = 0; k < remainders.size(); ++k)
    remainder_of[remainders[k].row] = std::int64_t(k);
  std::vector<GapParticle> out;
  out.reserve(n - plan.consumed_rows.size() + born.size());
  for (std::size_t i = 0; i < n; ++i) {
    if (drop[i] && remainder_of[i] < 0)
      continue;
    GapParticle g = rows[i];
    if (vel[i] >= 0)
      g.proper_velocity_m_s = plan.velocity[vel[i]].second;
    if (remainder_of[i] >= 0)
      g.weight = remainders[remainder_of[i]].weights.front();
    native_particle_guard(g);
    out.push_back(g);
  }
  for (const auto &g : born) {
    native_particle_guard(g);
    out.push_back(g);
  }
  if (out.size() > o.cap)
    throw refusal(RK::carrier_cap, "chemistry stage result exceeds the outer carrier budget");
  native.particles = NativeRows(std::move(out));
  native.next_id = next_id;
  native.fields_current = false;
  auto oc = std::make_shared<ChemOutcome>();
  oc->receipt = std::make_shared<FullTimeChemistryStageReceipt>();
  auto &receipt = oc->receipt;
  chem8::fill_receipt(*receipt, *setup, in, active, res, req);
  receipt->records_root = chem8::records_root(rec);
  receipt->records_sha256 = chem8::records_stream_sha256(rec);
  receipt->state_before_sha256 = chem8::state_digest(before);
  receipt->base_native_id = base_id;
  receipt->next_native_id = next_id;
  receipt->remainder_carriers = remainder_carriers;
  receipt->births_root = births_root(born_list);
  oc->next = next;
  oc->native = std::move(native);
  oc->inputs_sha256 = in.sha256;
  oc->purpose = in.cfg.purpose;
  return oc;
}
// Device outcome (API-008.txt 3.7): the stage on the device rows and census
// columns; the host receives the events of carriers with events, exact
// totals and digests, edits the packet table (O(changed packets)) and sends
// the row edits back; the new rows, columns and v5 carrier root stay on the
// device. Stages the device does not cover throw ChemDeviceFallback.
namespace full_time_detail {
// Development profile of the device outcome (PF_V8_PROFILE set): thread CPU ms per phase on stderr.
struct ChemProfile {
  bool on = std::getenv("PF_V8_PROFILE") != nullptr;
  double t0 = now();
  static double now() {
    timespec t;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
    return 1e3 * double(t.tv_sec) + 1e-6 * double(t.tv_nsec);
  }
  void mark(const char *what) {
    if (!on)
      return;
    const double t = now();
    std::fprintf(stderr, "chem-profile %s %.3f\n", what, t - t0);
    t0 = t;
  }
};
} // namespace full_time_detail
std::shared_ptr<FullTimeAccess::ChemOutcome>
FullTimeAccess::chemistry_device(FullCollision &s, const chem8::Setup &setup, const chem8::StageRequest &req,
                                 FullTimeLane &lane) {
  ChemProfile prof;
  auto &o = *s.outer;
  auto &owner = *o.owner;
  auto &dev = *lane.device;
  const auto &base = o.frame;
  const auto &snap = *base.snapshot;
  const auto &before = *snap.full_time_chemistry;
  const auto &reg = *setup.reg;
  const auto state_species = state_species_of(owner.states_, setup);
  std::vector<std::uint64_t> species_state(reg.states.size(), ~std::uint64_t(0));
  for (std::uint64_t k = 0; k < state_species.size(); ++k)
    if (state_species[k] != chem8::none)
      species_state[state_species[k]] = k;
  // frame rows on the device with census columns (one upload when absent)
  const auto &rows_ref = base.native.particles;
  if (!(rows_ref.device() && rows_ref.device()->columns)) {
    const auto &rows = rows_ref.host();
    std::vector<std::uint64_t> keys(rows.size()), states(rows.size());
    std::map<std::string, std::uint64_t> idx;
    std::uint64_t k = 0;
    for (const auto &[id, spec] : owner.states_) {
      (void)spec;
      idx[id] = k++;
    }
    for (const auto &[id, p] : snap.packets) {
      (void)id;
      const auto key = *std::min_element(p.native_ids.begin(), p.native_ids.end());
      for (auto nid : p.native_ids) {
        const auto it = std::lower_bound(rows.begin(), rows.end(), nid,
                                         [](const GapParticle &a, std::uint64_t b) { return a.id < b; });
        require(it != rows.end() && it->id == nid, "native canonical carrier missing");
        keys[std::size_t(it - rows.begin())] = key;
        states[std::size_t(it - rows.begin())] = idx.at(p.state);
      }
    }
    rows_ref.attach_device(device_upload_rows(dev, rows, &keys, &states));
  }
  const auto &drows = *rows_ref.device();
  ChemDeviceRequest rq;
  rq.rows = &drows;
  rq.setup = &setup;
  rq.state_species = state_species;
  rq.partitions = setup.config.device_partitions;
  auto &cfg = rq.cfg;
  std::memset(&cfg, 0, sizeof cfg);
  cfg.seed = setup.config.run_seed;
  cfg.step = req.step;
  cfg.dt_s = req.dt_s;
  cfg.trunc_sigmas = setup.config.trunc_sigmas;
  cfg.purpose = req.stage == 1 ? PFC_PURPOSE_MCC : PFC_PURPOSE_MCC_SECOND;
  cfg.max_candidates = setup.config.max_candidates;
  cfg.max_events = setup.config.max_events;
  if (!(req.step < PFA_STEP_LIMIT && std::isfinite(req.dt_s) && req.dt_s > 0))
    throw ChemDeviceFallback{};
  prof.mark("setup");
  auto ds = device_chem_stage(dev, rq);
  prof.mark("device_stage");
  if (ds.error || ds.active.totals.bad_carrier != ~0ull)
    throw ChemDeviceFallback{}; // the host path reports the refusal text
  // debit access: inventories and segment rows from the device
  std::map<std::uint32_t, ChemSegmentRow> seg_rows; // row -> segment row of debit candidates
  chem8::DebitAccess da;
  da.inventory = [&](std::uint32_t cell, std::uint32_t sp) { return device_chem_inventory(dev, ds, cell, sp); };
  da.packets = [&](std::uint32_t atlas, std::uint32_t sp) {
    const auto &ca = ds.active.cell_atlas;
    const auto it = std::lower_bound(ca.begin(), ca.end(), atlas);
    require(it != ca.end() && *it == atlas, "chemistry device segment cell");
    std::vector<std::pair<std::uint32_t, chem8::ESum>> out;
    for (const auto &r : device_chem_segment(dev, ds, std::uint32_t(it - ca.begin()), sp)) {
      if (r.key != r.id)
        continue;
      if (r.flags & 4u)
        throw ChemDeviceFallback{}; // multi-carrier target packet
      chem8::ESum c;
      c.add(r.weight);
      out.emplace_back(r.row, c);
      seg_rows[r.row] = r;
    }
    return out;
  };
  auto res = chem8::resolve(setup, ds.active, req, before, da);
  prof.mark("resolve");
  const auto &plan = res.plan;
  std::map<std::uint32_t, const chem8::ActiveCarrier *> active_of;
  for (const auto &c : ds.active.carriers)
    active_of[c.row] = &c;
  auto next = std::make_shared<ArenaSnapshot>(snap);
  next->full_time_chemistry = res.state;
  next->digest_cache.clear();
  PicArena::NativeState native = base.native;
  std::uint64_t next_id = native.next_id;
  const std::uint64_t base_id = next_id;
  const std::string stage_tag = std::string(chem8::packet_namespace) + "s" + std::to_string(req.step) + "/g" +
                                std::to_string(req.stage);
  ChemDeviceEdits ed;
  auto packet_of_key = [&](std::uint64_t key) -> const ArenaPacket & {
    const auto it = snap.full_time_packet_index.find(key);
    require(it != snap.full_time_packet_index.end(), "v8 packet key index entry");
    return snap.packets.at(it->second);
  };
  for (auto r : plan.consumed_rows) {
    const auto *c = active_of.at(r);
    const auto &p = packet_of_key(c->id);
    require(p.native_ids.size() == 1 && next->packets.erase(p.id) == 1 &&
                next->full_time_packet_index.erase(c->id) == 1,
            "chemistry consumed packet");
    ed.drop.push_back(r);
  }
  std::uint64_t remainder_carriers = 0;
  for (const auto &d : plan.debits) {
    const auto &sr = seg_rows.at(d.row);
    const auto &p = packet_of_key(sr.key);
    if (p.native_ids.size() != 1)
      throw ChemDeviceFallback{};
    if (d.whole) {
      require(next->packets.erase(p.id) == 1 && next->full_time_packet_index.erase(sr.key) == 1,
              "chemistry debited packet");
      ed.drop.push_back(d.row);
      continue;
    }
    ArenaCount c = p.count;
    std::vector<double> pieces;
    chem8::ESum left = d.amount;
    while (!left.is_zero()) {
      const double h = left.approx();
      pieces.push_back(h);
      chem8::ESum t;
      t.add(h);
      left.sub(t);
    }
    for (double x : pieces)
      if (x < 0)
        c = c.credited(-x);
    for (double x : pieces)
      if (x > 0)
        c = c.debit(x);
    const auto w = PicArena::positive_carriers(c);
    if (w.size() != 1)
      throw ChemDeviceFallback{}; // remainder is not one binary64 carrier
    auto &pp = next->packets.modify(p.id);
    pp.count = c;
    pp.last_event = stage_tag;
    ed.weight.push_back({d.row, w.front()});
  }
  for (const auto &v : plan.velocity)
    ed.velocity.push_back(v);
  auto carrier_of = [&](std::uint32_t sp, const PicVector &x, const double *u, double w) {
    GapParticle g;
    g.id = 0;
    g.position_m = x;
    g.proper_velocity_m_s = {u[0], u[1], u[2]};
    g.charge_C = reg.charge[sp] * pic_e;
    g.mass_kg = reg.mass_kg[sp];
    g.weight = w;
    return g;
  };
  auto issue = [&]() {
    require(next_id < std::numeric_limits<std::uint64_t>::max(), "native ID exhaustion");
    return next_id++;
  };
  std::vector<std::tuple<std::uint64_t, std::string, double>> born_list;
  auto append = [&](GapParticle g, std::uint32_t sp, const std::string &pid, const std::string &origin, double w) {
    g.id = issue();
    closed_position(native.options, g);
    native_particle_guard(g);
    require(next->packets.emplace(pid, ArenaPacket{pid, reg.states[sp], origin, "", ArenaCount(w), {g.id}}).second &&
                next->full_time_packet_index.emplace(g.id, pid).second,
            "chemistry new packet collision");
    ed.append.push_back(g);
    ed.append_keys.push_back(g.id);
    ed.append_states.push_back(species_state.at(sp));
    born_list.emplace_back(g.id, pid, w);
  };
  for (const auto &b : plan.births) {
    const auto *c = active_of.at(b.row);
    append(carrier_of(b.species, PicVector{c->x[0], c->x[1], c->x[2]}, b.u, b.weight), b.species,
           chem8::packet_id_birth(req.step, req.stage, b),
           std::string(chem8::packet_namespace) + "p" + std::to_string(b.process), b.weight);
  }
  for (const auto &rl : plan.releases)
    append(carrier_of(rl.species, PicVector{rl.x[0], rl.x[1], rl.x[2]}, rl.u, rl.weight), rl.species,
           chem8::packet_id_release(req.step, rl), std::string(chem8::packet_namespace) + "release", rl.weight);
  for (auto sp : species_state)
    if (sp == ~std::uint64_t(0))
      throw ChemDeviceFallback{}; // a chemistry species without a census state
  const std::size_t n_out = drows.n - ed.drop.size() + ed.append.size();
  if (n_out > o.cap)
    throw refusal(RK::carrier_cap, "chemistry stage result exceeds the outer carrier budget");
  prof.mark("packets");
  auto nrows = device_chem_apply(dev, drows, ed);
  prof.mark("apply");
  native.particles = NativeRows::from_device(nrows, nrows->n);
  native.next_id = next_id;
  native.fields_current = false;
  auto oc = std::make_shared<ChemOutcome>();
  {
    auto dg = std::make_shared<cn::BlockDigest>();
    dg->root = nrows->carrier_root;
    dg->records = nrows->n;
    oc->digest = dg;
  }
  oc->receipt = std::make_shared<FullTimeChemistryStageReceipt>();
  auto &receipt = *oc->receipt;
  chem8::StageInputs in;
  in.cfg = cfg;
  in.sha256 = ds.inputs_sha256;
  chem8::fill_receipt(receipt, setup, in, ds.active, res, req);
  receipt.records_root = ds.records_root;
  receipt.state_before_sha256 = chem8::state_digest(before);
  receipt.base_native_id = base_id;
  receipt.next_native_id = next_id;
  receipt.remainder_carriers = remainder_carriers;
  receipt.births_root = births_root(born_list);
  receipt.device = true;
  receipt.device_ms = ds.device_wall_ms;
  prof.mark("receipt");
  oc->next = next;
  oc->native = std::move(native);
  oc->inputs_sha256 = ds.inputs_sha256;
  oc->purpose = cfg.purpose;
  return oc;
}
PreparedArena FullTimeAccess::chemistry(FullCollision &s) {
  live(s);
  auto &o = *s.outer;
  auto &owner = *o.owner;
  const auto setup = std::static_pointer_cast<const chem8::Setup>(owner.chemistry_);
  require(bool(setup), "chemistry v3 is not configured for this arena");
  require(!s.building && !o.builder_active, "private compound cap or reentry");
  require(s.outstanding_token.expired() &&
              s.generation <= std::numeric_limits<std::uint64_t>::max() - 2,
          "private token/generation resource gate");
  const bool first_slot = s.exposure.slot == FullTimeSlot::before_transport;
  const std::uint64_t used = o.frame.snapshot->revision - o.before.snapshot->revision;
  const auto room = window_room(*o.frame.snapshot);
  if (first_slot) {
    require(used <= o.records && o.records - used >= 4, "reserve both chemistry stage record budgets");
    if (room < 4)
      throw refusal(RK::window_rows, "reserve both chemistry stage history window rows");
  } else {
    require(used <= o.records && o.records - used >= 2, "second chemistry stage record budget");
    if (room < 2)
      throw refusal(RK::window_rows, "second chemistry stage history window rows");
  }
  const auto &base = o.frame;
  const auto &snap = *base.snapshot;
  require(snap.full_time_identity == 1 && bool(snap.full_time_chemistry), "frame without chemistry state");
  const auto &before = *snap.full_time_chemistry;
  const chem8::StageRequest req{o.common.accepted_step, first_slot ? 1u : 2u, s.exposure.collision_duration_s,
                              std::min<std::size_t>(o.cap, 1000000)};
  const bool base_valid = owner.full_time_validated(snap, base.native.particles, o.cap);
  auto *lane = lane_of(o);
  const auto exec = setup->config.execution;
  std::shared_ptr<ChemOutcome> oc;
  if (exec != FullTimeChemistryExecution::host && lane && lane->device) {
    try {
      oc = chemistry_device(s, *setup, req, *lane);
    } catch (const ChemDeviceFallback &) {
      oc.reset();
    }
    if (oc && exec == FullTimeChemistryExecution::device_verified) {
      auto h = chemistry_host(s, *setup, req);
      const auto &a = *oc->receipt, &b = *h->receipt;
      const auto hroot = cn::carrier_rows_digest(h->native.particles)->root;
      if (a.records_root != b.records_root || a.inputs_sha256 != b.inputs_sha256 || a.ledger_sha256 != b.ledger_sha256 ||
          a.state_after_sha256 != b.state_after_sha256 || a.births_root != b.births_root ||
          a.next_native_id != b.next_native_id || oc->digest->root != hroot ||
          snapshot_digest(*oc->next) != snapshot_digest(*h->next) || oc->native.particles.size() != h->native.particles.size())
        throw std::logic_error("native v8 chemistry device stage and host path differ");
      ++device_counters(*lane->device).verified_comparisons;
      (void)owner.full_time_census(*h->next, h->native.particles.host(), o.cap, false);
    }
  }
  if (!oc)
    oc = chemistry_host(s, *setup, req);
  ChemProfile prof;
  auto next = oc->next;
  auto &native = oc->native;
  auto receipt = oc->receipt;
  const std::uint64_t base_id = receipt->base_native_id, next_id = receipt->next_native_id;
  const auto remainder_carriers = receipt->remainder_carriers;
  struct InShim {
    std::string sha256;
    struct {
      std::uint64_t purpose;
    } cfg;
  } in{oc->inputs_sha256, {oc->purpose}};
  ArenaStochasticCompletion c;
  c.integration_owner_id = chem8::owner_id;
  c.before_rng_clock_sha256 = receipt->state_before_sha256;
  c.after_rng_clock_sha256 = receipt->state_after_sha256;
  c.common_clock_sha256 = s.exposure.common.common_clock_sha256;
  c.channel_id = setup->channel.id;
  c.dataset_sha256 = setup->channel.dataset_sha256;
  c.kernel_manifest_sha256 = setup->channel.kernel_manifest_sha256;
  c.kernel_binding_file_sha256 = setup->channel.kernel_binding_file_sha256;
  c.model = setup->channel.model;
  c.stage_sequence = before.stages;
  c.frozen_physical_stage_s = s.exposure.frozen_physical_stage_s;
  c.common_interval_begin_s = s.exposure.common.begin_s;
  c.common_interval_end_s = s.exposure.common.end_s;
  completion_chain_guard(s.exposure, c, *o.boundary,
                         first_slot ? std::nullopt
                                    : std::optional<ArenaStochasticCompletion>(o.first.completion));
  completion_validate(c, owner.channels_.at(c.channel_id), frame_completion(snap, c.integration_owner_id),
                      base.native.ledger.time_s);
  FullRow row;
  row.kind = RowKind::chemistry;
  auto &k = row.chemistry;
  k.accepted_step = req.step;
  k.stage = req.stage;
  k.purpose = in.cfg.purpose;
  k.dt_s = req.dt_s;
  k.base_native_id = base_id;
  k.next_native_id = next_id;
  k.projectiles = receipt->projectiles;
  k.events = receipt->events;
  k.admitted_events = receipt->admitted_events;
  k.refused_debits = receipt->refused_debits;
  k.births = receipt->births;
  k.consumed = receipt->consumed;
  k.debited_packets = receipt->debited_packets;
  k.removed_packets = receipt->removed_packets;
  k.released_packets = receipt->released_packets;
  k.remainder_carriers = remainder_carriers;
  k.configuration_sha256 = setup->configuration_sha256;
  k.inputs_sha256 = in.sha256;
  k.records_root = receipt->records_root;
  k.ledger_sha256 = receipt->ledger_sha256;
  k.state_before_sha256 = receipt->state_before_sha256;
  k.state_after_sha256 = receipt->state_after_sha256;
  k.births_root = receipt->births_root;
  full_time_append(*next, row);
  const auto stage_chain = next->chain;
  next->stochastic_completions[c.integration_owner_id] = completion_json(c);
  FullRow er;
  er.kind = RowKind::exposure;
  er.exposure = s.exposure;
  er.completion = c;
  er.stage_before_chain = snap.chain;
  er.stage_after_reaction_chain = stage_chain;
  full_time_append(*next, er);
  require(next->revision - o.before.snapshot->revision <= o.records, "outer record resource cap");
  ++s.generation;
  auto t = std::make_unique<PreparedArena::Impl>();
  t->owner = o.owner;
  t->private_identity = s.token_identity;
  t->private_generation = s.generation + 1;
  t->before = base.snapshot;
  t->after = next;
  t->native_after = native;
  t->native_after_digest = oc->digest ? oc->digest
                                       : cn::carrier_rows_digest(t->native_after.particles,
                                                                 host_rows_if_present(&base.native), base.digest.get());
  if (base_valid)
    owner.full_time_mark_validated(*next, t->native_after.particles, o.cap);
  t->full_time_exposure = true;
  t->completion = c;
  t->chemistry = receipt;
  t->audit = {{"kind", "full_time_exposure"},
              {"exposure", exposure_json(s.exposure)},
              {"completion", completion_json(c)},
              {"chemistry_stage",
               {{"stage", req.stage},
                {"records_root", receipt->records_root},
                {"inputs_sha256", receipt->inputs_sha256},
                {"ledger_sha256", receipt->ledger_sha256},
                {"projectiles", receipt->projectiles},
                {"events", receipt->events},
                {"admitted_events", receipt->admitted_events},
                {"births", receipt->births},
                {"released_packets", receipt->released_packets}}},
              {"row_sha256", next->chain}};
  prof.mark("token_rows");
  t->candidate_prefix = prefix(owner, *next, t->native_after, *t->native_after_digest);
  prof.mark("token_prefix");
  t->private_memory_lease = std::make_shared<int>(0);
  ++s.generation;
  s.outstanding_token = t->private_memory_lease;
  return PreparedArena(std::move(t));
}
std::shared_ptr<const FullTimeChemistryStageReceipt> PreparedArena::candidate_chemistry_stage() const {
  require(bool(impl_), "empty prepared token");
  return impl_->chemistry;
}
std::vector<ArenaStateSpec> full_time_metal_vapour_states(bool with_fe_cr_ni) {
  constexpr double u = 1.66053906660e-27; // kg (CODATA 2018)
  std::vector<ArenaStateSpec> out;
  if (with_fe_cr_ni)
    for (const auto &[name, m] : {std::pair{"metal-v1/Fe", 55.845}, {"metal-v1/Cr", 51.9961}, {"metal-v1/Ni", 58.6934}})
      out.push_back({name, 0, 0, 0, m * u, 0, 0});
  for (const auto &[name, m] : {std::pair{"metal-v1/Mn", 54.938043}, {"metal-v1/Si", 28.085}, {"metal-v1/C", 12.011}})
    out.push_back({name, 0, 0, 0, m * u, 0, 0});
  return out;
}
FullTimeFlowsV81Status PicArena::full_time_flows_v81() const {
  FullTimeFlowsV81Status out;
  require(bool(accepted_), "accepted state");
  const auto &fp = accepted_->full_time_v81;
  if (!fp)
    return out;
  const auto &f = *fp;
  out.present = true;
  out.wall_births = f.wall_births;
  out.returned_births = f.returned_births;
  out.bridges = f.bridges;
  out.carrying_packets = f.internal.size();
  out.internal_born_J = chem8::exact_value(f.internal_born_J);
  out.internal_sunk_J = chem8::exact_value(f.internal_sunk_J);
  out.internal_released_J = chem8::exact_value(f.internal_released_J);
  chem8::ESum flight;
  for (const auto &[id, x] : f.internal) {
    (void)id;
    flight.add(x.energy_J);
  }
  out.internal_in_flight_J = chem8::exact_value(flight);
  out.extra_work_J = chem8::exact_value(f.extra_work_J);
  for (const auto &[st, e] : f.wall_weight)
    out.wall_weight.emplace_back(st, chem8::exact_value(e));
  for (const auto &[st, e] : f.returned_weight)
    out.returned_weight.emplace_back(st, chem8::exact_value(e));
  out.merge_transactions = f.merge_transactions;
  out.merge_groups = f.merge_groups;
  out.merge_inputs = f.merge_inputs;
  out.merge_outputs = f.merge_outputs;
  out.merge_relocations = f.merge_relocations;
  out.merge_formation_refused = f.merge_formation_refused;
  out.merge_centroid_fallback = f.merge_centroid_fallback;
  out.merge_status = f.merge_status;
  for (const auto &[st, r] : f.merge_residual)
    out.merge_residual.emplace_back(st, std::array<FullTimeExactValue, 4>{chem8::exact_value(r[0]), chem8::exact_value(r[1]),
                                                                         chem8::exact_value(r[2]), chem8::exact_value(r[3])});
  out.has_last_merge = f.has_last_merge;
  out.last_merge_step = f.last_merge_step;
  out.last_merge_base = f.last_merge_base;
  out.last_merge_outputs = f.last_merge_outputs;
  out.last_merge_sha256 = f.last_merge_sha256;
  out.merge_policy_sha256 = f.merge_policy_sha256;
  out.sha256 = chem8::flows_v81_digest(f);
  return out;
}
FullTimeChemistryStatus PicArena::full_time_chemistry() const {
  FullTimeChemistryStatus out;
  const auto setup = std::static_pointer_cast<const chem8::Setup>(chemistry_);
  if (!setup || !accepted_->full_time_chemistry)
    return out;
  const auto &st = *accepted_->full_time_chemistry;
  out.enabled = true;
  out.diagnostic_field_target_approximation = setup->config.diagnostic_field_target_approximation;
  out.full_system_energy_closure = false;
  out.state_sha256 = chem8::state_digest(st);
  out.configuration_sha256 = setup->configuration_sha256;
  out.registry_sha256 = setup->reg->sha256;
  out.blob_sha256 = setup->reg->blob->sha256;
  out.stages = st.stages;
  out.reservoirs = st.reservoirs.size();
  out.events = st.cumulative.admitted_events;
  out.cumulative_energy = chem8::energy_values(st.cumulative);
  out.cumulative_counts = chem8::count_values(*setup, st.cumulative);
  out.cumulative_ledger_sha256 = chem8::ledger_digest(st.cumulative);
  std::map<std::uint32_t, chem8::ESum> stock;
  for (const auto &[k, r] : st.reservoirs)
    stock[k.species].add(r.count);
  for (const auto &[sp, v] : stock)
    out.reservoir_counts.push_back({setup->reg->states[sp], chem8::exact_value(v)});
  out.coverage_json = chem8::coverage_json(*setup, st);
  if (const auto &ml = accepted_->full_time_melt) {
    out.melt.credits = ml->credits;
    for (int k = 0; k < 6; ++k)
      out.melt.atoms[k] = chem8::exact_value(ml->atoms[k]);
    out.melt.energy_J = chem8::exact_value(ml->energy_J);
    out.melt.kinetic_J = chem8::exact_value(ml->kinetic_J);
    out.melt.last_record_sha256 = ml->last_record_sha256;
    out.melt.sha256 = chem8::melt_digest(*ml);
  }
  return out;
}
FullTimeChemistryRegistry full_time_chemistry_registry(const FullTimeChemistryConfiguration &c,
                                                       const GapPicOptions &o) {
  const auto s = chem8::make_setup(c, o);
  FullTimeChemistryRegistry r;
  const auto &reg = *s->reg;
  const auto &b = *reg.blob;
  r.blob_sha256 = b.sha256;
  r.registry_sha256 = reg.sha256;
  r.configuration_sha256 = s->configuration_sha256;
  r.policy_sha256 = s->policy_sha256;
  r.grid_sha256 = s->grid.sha256;
  r.initial_rng_sha256 = chem8::state_digest(*chem8::initial_state(*s));
  for (std::size_t i = 0; i < reg.states.size(); ++i) {
    FullTimeChemistrySpecies sp;
    sp.state = reg.states[i];
    sp.mass_kg = reg.mass_kg[i];
    sp.mass_eV = b.view.species[i].mass_eV;
    sp.reference_eV = b.view.species[i].ref_eV;
    sp.charge_e = reg.charge[i];
    sp.nuclei = reg.nuclei[i];
    sp.projectile = reg.projectile[i];
    sp.target = reg.target[i];
    r.species.push_back(sp);
  }
  r.channel = s->channel;
  r.region_names = b.region_names;
  for (auto x : reg.region_bounded)
    r.region_bounded.push_back(x != 0);
  return r;
}
std::vector<ArenaStateSpec> full_time_chemistry_states(const FullTimeChemistryRegistry &r) {
  std::vector<ArenaStateSpec> out;
  for (const auto &s : r.species) {
    ArenaStateSpec a;
    a.id = s.state;
    a.hydrogen = s.nuclei[0];
    a.oxygen = s.nuclei[1];
    a.charge_e = s.charge_e;
    a.mass_kg = s.mass_kg;
    a.internal_J = 0;
    a.reference_J = s.reference_eV * 1.602176634e-19;
    out.push_back(a);
  }
  return out;
}
std::vector<std::string> full_time_chemistry_recommended_admitted(const FullTimeChemistryRegistry &r) {
  (void)r;
  // No unbounded cross-section region is recommended for production.
  // Explicit diagnostic policies may list them, with exposure recorded.
  return {};
}
FullTimeChemistryConfiguration full_time_chemistry_preset(const std::string &name, const GapPicOptions &o,
                                                          const std::string &blob_path, std::uint64_t run_seed,
                                                          double release_weight,
                                                          const std::array<std::uint32_t, 3> &grid_counts) {
  using CK = FullTimeChemistryRefusal::Kind;
  if (name != "production-v1")
    throw FullTimeChemistryRefusal(CK::configuration, "chemistry_v3/configuration: unknown preset " + name);
  FullTimeChemistryConfiguration c;
  c.enabled = true;
  const std::string dir = full_time_chemistry_v31_tables;
  c.blob_path = blob_path.empty() ? dir + "chemistry-v3.1.blob" : blob_path;
  c.blob_sha256 = full_time_chemistry_production_v1_blob_sha256;
  c.run_seed = run_seed;
  c.release_weight = release_weight;
  c.trunc_sigmas = 7;
  c.max_candidates = 4096;
  c.max_events = 64;
  c.grid_counts = grid_counts;
  c.policy.event_limit = 1.0;   // checked equal to the policy file's values
  c.policy.event_fraction = 0.01;
  c.policy.policy_file = dir + "chemistry-v3.1-coverage-policy.json";
  c.policy.policy_file_sha256 = full_time_chemistry_production_v1_policy_sha256;
  c.preset = name;
  (void)full_time_chemistry_registry(c, o); // validates the blob and the policy file now
  return c;
}
std::shared_ptr<const FullTimeStepFlows> PicArena::full_time_published_flows() const {
  return accepted_->full_time_step_flows;
}
void PrivateFullTimeStep::credit_ejected_melt(const FullTimeEjectedMelt &m) {
  auto &o = *impl_->outer;
  FullTimeAccess::live(o);
  require(o.frame.snapshot->full_time_identity == 1, "ejected-melt credit requires identity schema v8");
  require(o.phase == 3 && !o.builder_active, "ejected-melt credit after the second exposure, before seal");
  require(o.frame.snapshot->full_time_melt_step != o.common.accepted_step + 1, "ejected-melt credit twice in one step");
  for (double a : m.atoms)
    require(std::isfinite(a) && a >= 0, "ejected-melt atoms");
  require(std::isfinite(m.energy_J) && std::isfinite(m.kinetic_J) && m.kinetic_J >= 0 &&
              digest(m.phase_record_sha256),
          "ejected-melt energies and record digest");
  auto melt = o.frame.snapshot->full_time_melt ? std::make_shared<FullTimeMeltState>(*o.frame.snapshot->full_time_melt)
                                               : std::make_shared<FullTimeMeltState>();
  melt->digest.clear();
  ++melt->credits;
  for (int k = 0; k < 6; ++k)
    melt->atoms[k].add(m.atoms[k]);
  melt->energy_J.add(m.energy_J);
  melt->kinetic_J.add(m.kinetic_J);
  melt->last_record_sha256 = m.phase_record_sha256;
  auto next = std::make_shared<ArenaSnapshot>(*o.frame.snapshot);
  next->full_time_melt = melt;
  next->full_time_melt_step = o.common.accepted_step + 1;
  next->digest_cache.clear();
  o.frame.snapshot = next; // rows (and their digest) unchanged
}
