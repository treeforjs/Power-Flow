// Canonical diagnostic facts for a future prepared gas-owner transaction.
// No gas owner is bound here and the strict eventful-stage refusal is unchanged.
namespace full_time_detail {
void export_count(Writer &w, const ArenaCount &c) {
  auto value = c;
  codec_count(w, value);
}
void export_packet(Writer &w, const FullTimeChemistryPacketSnapshot &p) {
  w.text(p.packet_id); w.text(p.state); export_count(w, p.count);
  w.f64(p.carried_internal_energy_J); w.u64(p.carriers.size());
  for (const auto &r : p.carriers) {
    w.u64(r.id);
    for (double v : r.position_m) w.f64(v);
    for (double v : r.proper_velocity_m_s) w.f64(v);
    w.f64(r.charge_C); w.f64(r.mass_kg); w.f64(r.weight);
  }
}
std::string export_admitted_root(const FullTimeChemistryStageExport &x) {
  Writer w; w.tag("native-v8.3/admitted-prefixes-v1");
  w.u64(x.admissions.size());
  for (const auto &p : x.admissions) {
    w.u64(p.carrier); w.u64(p.raw_begin); w.u64(p.raw_count); w.u64(p.admitted_count);
    w.u64(p.compact_cell); w.u64(p.atlas_cell);
  }
  w.u64(x.events.size());
  for (const auto &e : x.events) {
    w.u64(e.raw_index); w.u64(e.atlas_cell); w.u64(e.target_species);
    w.boolean(e.consumes_tracked_target); w.f64(e.weight); w.text(e.event_bytes);
  }
  return cn::sha256_hex(w.bytes);
}
std::string export_debits_root(const FullTimeChemistryStageExport &x) {
  Writer w; w.tag("native-v8.3/actual-target-debits-v1"); w.u64(x.target_debits.size());
  for (const auto &d : x.target_debits) {
    w.u64(d.atlas_cell); w.u64(d.species); export_count(w, d.amount);
    w.boolean(d.also_projectile); w.boolean(d.internal_energy_deferred);
    export_packet(w, d.before); export_packet(w, d.remainder);
  }
  return cn::sha256_hex(w.bytes);
}
std::string export_releases_root(const FullTimeChemistryStageExport &x) {
  Writer w; w.tag("native-v8.3/actual-parcel-releases-v1"); w.u64(x.releases.size());
  for (const auto &r : x.releases) export_packet(w, r);
  return cn::sha256_hex(w.bytes);
}
std::string export_digest(const FullTimeChemistryStageExport &x) {
  Writer w; w.tag(x.schema);
  for (auto v : {x.accepted_step, x.stage, x.purpose, x.base_native_id, x.next_native_id, x.owner_stage_sequence}) w.u64(v);
  w.f64(x.begin_s); w.f64(x.end_s); w.f64(x.step_duration_s); w.f64(x.collision_duration_s);
  w.f64(x.frozen_physical_stage_s);
  for (const auto *p : {&x.accepted_native_sha256, &x.common_clock_sha256, &x.frame_snapshot_sha256,
                       &x.frame_history_sha256, &x.configuration_sha256, &x.registry_sha256,
                       &x.inputs_sha256, &x.chemistry_before_sha256, &x.chemistry_after_sha256,
                       &x.frozen_field_sha256, &x.ledger_sha256, &x.births_root, &x.raw_records_root,
                       &x.admitted_root, &x.target_debits_root, &x.releases_root, &x.host_target_feedback_sha256})
    w.text(*p);
  return cn::sha256_hex(w.bytes);
}
ArenaCount export_sum(const ArenaCount &a, const ArenaCount &b) {
  auto c = a;
  for (double v : b.positive_carriers()) c = c.credited(v);
  return c;
}
ArenaCount export_difference(const ArenaCount &a, const ArenaCount &b) {
  auto c = a;
  for (double v : b.positive_carriers()) c = c.debit(v);
  return c;
}
template<class RowAt>
FullTimeChemistryPacketSnapshot export_snapshot(const ArenaSnapshot &s, const ArenaPacket &p, RowAt row_at) {
  FullTimeChemistryPacketSnapshot out;
  out.packet_id = p.id; out.state = p.state; out.count = p.count;
  for (auto id : p.native_ids) out.carriers.push_back(row_at(id));
  std::sort(out.carriers.begin(), out.carriers.end(), [](const auto &a, const auto &b) { return a.id < b.id; });
  if (s.full_time_v81) {
    const auto it = s.full_time_v81->internal.find(p.id);
    if (it != s.full_time_v81->internal.end()) out.carried_internal_energy_J = it->second.energy_J;
  }
  return out;
}
std::shared_ptr<FullTimeChemistryStageExport> export_begin(const chem8::Resolution &res) {
  auto x = std::make_shared<FullTimeChemistryStageExport>();
  x->admissions = res.admissions; x->events = res.admitted_records;
  return x;
}
void export_finish(FullTimeChemistryStageExport &x, const FullCollision &s, const chem8::Setup &setup,
                   const FullTimeChemistryStageReceipt &r) {
  const auto &o = *s.outer;
  x.accepted_step = r.accepted_step; x.stage = r.stage; x.purpose = r.purpose;
  x.base_native_id = r.base_native_id; x.next_native_id = r.next_native_id;
  x.begin_s = o.common.begin_s; x.end_s = o.common.end_s;
  x.step_duration_s = o.common.dt_s;
  x.collision_duration_s = s.exposure.collision_duration_s;
  x.owner_stage_sequence = s.exposure.owner_stage_sequence;
  x.frozen_physical_stage_s = s.exposure.frozen_physical_stage_s;
  x.frozen_field_sha256 = s.exposure.frozen_field_sha256;
  x.accepted_native_sha256 = o.common.accepted_arena_sha256;
  x.common_clock_sha256 = o.common.common_clock_sha256;
  x.frame_snapshot_sha256 = snapshot_digest(*o.frame.snapshot);
  x.frame_history_sha256 = o.frame.snapshot->chain;
  x.configuration_sha256 = setup.configuration_sha256; x.registry_sha256 = setup.reg->sha256;
  x.inputs_sha256 = r.inputs_sha256; x.chemistry_before_sha256 = r.state_before_sha256;
  x.chemistry_after_sha256 = r.state_after_sha256; x.ledger_sha256 = r.ledger_sha256; x.births_root = r.births_root;
  x.raw_records_root = r.records_root;
  if(r.host_target_feedback)x.host_target_feedback_sha256=r.host_target_feedback->sha256;
  std::uint64_t raw_count = 0;
  for (const auto &p : x.admissions) raw_count = checked_add(raw_count, p.raw_count, "export raw events");
  require(raw_count == r.events && x.events.size() == r.admitted_events &&
          x.target_debits.size() == r.debited_packets && x.releases.size() == r.released_packets,
          "export counts differ from the applied chemistry receipt");
  x.admitted_root = export_admitted_root(x); x.target_debits_root = export_debits_root(x);
  x.releases_root = export_releases_root(x); x.sha256 = export_digest(x);
  validate_full_time_chemistry_stage_export(x, false);
}
} // namespace full_time_detail

void validate_full_time_chemistry_stage_export(const FullTimeChemistryStageExport &x, bool unambiguous) {
  using namespace full_time_detail;
  using CK = FullTimeChemistryRefusal::Kind;
  auto check = [&](bool ok, const char *message) {
    if (!ok) throw FullTimeChemistryRefusal(CK::event, std::string("admitted_export: ") + message);
  };
  check(x.schema == "native-v8.5/admitted-target-export-v2", "schema");
  check((x.stage == 1 || x.stage == 2) && std::isfinite(x.begin_s) && std::isfinite(x.end_s) &&
        x.end_s > x.begin_s && std::isfinite(x.step_duration_s) && x.step_duration_s > 0 &&
        x.end_s == x.begin_s + x.step_duration_s && x.collision_duration_s == .5 * x.step_duration_s &&
        x.purpose == (x.stage == 1 ? PFC_PURPOSE_MCC : PFC_PURPOSE_MCC_SECOND) &&
        x.frozen_physical_stage_s == (x.stage == 1 ? x.begin_s : x.end_s) &&
        x.base_native_id > 0 && x.next_native_id >= x.base_native_id, "stage binding");
  for (const auto *p : {&x.accepted_native_sha256, &x.common_clock_sha256, &x.frame_snapshot_sha256,
                       &x.frame_history_sha256, &x.configuration_sha256, &x.registry_sha256,
                       &x.inputs_sha256, &x.chemistry_before_sha256, &x.chemistry_after_sha256,
                       &x.frozen_field_sha256, &x.ledger_sha256, &x.births_root, &x.raw_records_root})
    check(digest(*p), "identity digest");
  check(x.host_target_feedback_sha256.empty() || digest(x.host_target_feedback_sha256),"host target feedback digest");
  using Key = std::pair<std::uint32_t, std::uint32_t>;
  std::map<Key, ArenaCount> admitted, debited;
  std::size_t ei = 0;
  std::uint64_t raw_end = 0, previous = 0;
  for (const auto &p : x.admissions) {
    check(p.carrier > previous && p.raw_begin == raw_end && p.raw_count > 0 &&
          p.admitted_count <= p.raw_count && p.raw_count <= ~std::uint64_t(0) - raw_end, "prefix ordering");
    previous = p.carrier; raw_end += p.raw_count;
    for (std::uint64_t j = 0; j < p.admitted_count; ++j) {
      check(ei < x.events.size(), "missing admitted record");
      const auto &r = x.events[ei++];
      check(r.raw_index == p.raw_begin + j && r.atlas_cell == p.atlas_cell &&
            r.event_bytes.size() == sizeof(pfc::EventRec) && std::isfinite(r.weight) && r.weight > 0,
            "admitted record layout/order");
      pfc::EventRec e{}; std::memcpy(&e, r.event_bytes.data(), sizeof e);
      check(e.carrier == p.carrier && e.cell == p.compact_cell && e.target_species == r.target_species &&
            std::bit_cast<std::uint64_t>(e.weight) == std::bit_cast<std::uint64_t>(r.weight) &&
            e.n_out > 0 && e.n_out <= PFC_MAX_OUT, "admitted record binding");
      bool gas_return = false;
      for (pfc_u32 k = 0; k < e.n_out; ++k) gas_return |= e.out_dest[k] == PFC_DEST_GAS;
      check(!gas_return || !r.consumes_tracked_target, "GAS-return event debits inventory");
      if (r.consumes_tracked_target) {
        auto &c = admitted[{r.atlas_cell, r.target_species}]; c = c.credited(r.weight);
      }
    }
  }
  check(ei == x.events.size(), "unselected admitted records");
  std::set<std::uint64_t> before_ids, after_ids;
  std::set<std::string> packet_ids;
  auto packet = [&](const FullTimeChemistryPacketSnapshot &p, bool before) {
    check(std::isfinite(p.carried_internal_energy_J) && p.carried_internal_energy_J >= 0,
          "packet internal energy");
    if (p.packet_id.empty()) {
      check(p.state.empty() && p.count.zero() && p.carriers.empty() && p.carried_internal_energy_J == 0,
            "absent remainder"); return;
    }
    check(!p.state.empty() && !p.count.zero() && !p.carriers.empty(), "resident packet");
    ArenaCount sum; std::uint64_t prior = 0;
    const auto &first = p.carriers.front();
    for (const auto &r : p.carriers) {
      check(r.id > prior && (before ? before_ids : after_ids).insert(r.id).second,
            "duplicate carrier ownership"); prior = r.id;
      native_particle_guard(r);
      check(r.position_m == first.position_m && r.proper_velocity_m_s == first.proper_velocity_m_s &&
            r.mass_kg == first.mass_kg && r.charge_C == first.charge_C, "divergent packet phase space");
      check(before ? r.id < x.base_native_id : r.id < x.next_native_id, "carrier ID range");
      sum = sum.credited(r.weight);
    }
    check(sum == p.count, "packet exact count");
  };
  for (const auto &d : x.target_debits) {
    check(!d.before.packet_id.empty() && packet_ids.insert(d.before.packet_id).second && !d.amount.zero(),
          "duplicate or empty target debit");
    packet(d.before, true); packet(d.remainder, false);
    check(d.remainder.packet_id.empty() || (d.remainder.packet_id == d.before.packet_id &&
          d.remainder.state == d.before.state), "remainder identity/species");
    check(export_sum(d.amount, d.remainder.count) == d.before.count, "debit does not telescope exactly");
    check(d.internal_energy_deferred == (d.before.carried_internal_energy_J != 0), "internal ownership label");
    auto &c = debited[{d.atlas_cell, d.species}]; c = export_sum(c, d.amount);
    if (unambiguous && (d.also_projectile || d.internal_energy_deferred))
      throw FullTimeChemistryRefusal(CK::target_energy_owner, "admitted_export: target debit attribution unresolved");
  }
  check(admitted == debited, "admitted consuming count differs from actual target debit");
  for (const auto &p : x.releases) {
    check(x.stage == 2 && !p.packet_id.empty() && packet_ids.insert(p.packet_id).second, "release identity/stage");
    packet(p, false);
    for (const auto &r : p.carriers) check(r.id >= x.base_native_id, "release native ID");
  }
  check(x.admitted_root == export_admitted_root(x) && x.target_debits_root == export_debits_root(x) &&
        x.releases_root == export_releases_root(x) && x.sha256 == export_digest(x), "commitment mismatch");
}
