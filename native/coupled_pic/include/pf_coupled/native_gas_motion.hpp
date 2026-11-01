#pragma once
#include "raikou/mfem_pic_full_time.h"
#include "pf_gas_v1/gas_owner.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>

// Pure checked records conversion. These inputs and outputs are NOT live
// publication capabilities. The common bridge must obtain the view from its
// exact native prepared-motion token and revalidate that token before mutation.
namespace native_gas_motion {
namespace n = raikou::mfem_backend;
namespace g = pf_gas_v1;
inline void need(bool x, const char *why) {
  if (!x) throw std::invalid_argument(std::string("native_gas_motion: ") + why);
}
inline bool same(double x, double y) {
  return std::bit_cast<std::uint64_t>(x) == std::bit_cast<std::uint64_t>(y);
}
struct View {
  std::uint64_t accepted_step = 0;
  std::span<const n::ArenaRateEntry> before, start, after;
  std::span<const n::FullTimeStoredPacketEnergy> before_internal, start_internal, after_internal;
  const n::FullTimeBirthReceipt *births = nullptr;
  const n::FullTimeSinkLedger *sinks = nullptr;
  bool has_merge = false;
};
struct Records {
  std::vector<g::FlowRecord> births, sinks;
  std::vector<g::MoverRecord> movers;
  // Gas owner's historical wall/* aggregate combines all population exits.
  // Preserve the native kind/surface so the common ledger can distinguish
  // conductor heating from feed/exterior escape; never assign both as heat.
  std::vector<n::FullTimeSinkRecord> sink_destinations;
};
namespace detail {
using Rows = std::map<std::uint64_t, const n::ArenaRateEntry *>;
using Energy = std::map<std::string, const n::FullTimeStoredPacketEnergy *>;
using Species = std::map<std::string, std::uint32_t>;
inline std::string native_digits(const g::Exact &x) {
  // Native chem8::exact_digits uses every normalised forms digit, high first.
  // Gas digits_text instead omits zero spans and inserts separators.
  const auto d = x.normalised(); std::string out; out.reserve(d.size() * 16);
  for (auto i = d.rbegin(); i != d.rend(); ++i) {
    char s[17]; std::snprintf(s, sizeof s, "%016llx", static_cast<unsigned long long>(*i)); out += s;
  }
  return out;
}
inline bool same_sink(const n::FullTimeSinkRecord &a, const n::FullTimeSinkRecord &b) {
  if (a.native_id != b.native_id || a.kind != b.kind || a.surface != b.surface || a.axis != b.axis ||
      !same(a.plane_m,b.plane_m) || !same(a.fraction,b.fraction) || !same(a.flight_s,b.flight_s) ||
      !same(a.discarded_flight_s,b.discarded_flight_s)) return false;
  for (int d=0;d<3;++d)
    if (!same(a.segment_m[d],b.segment_m[d]) || !same(a.hit_point_m[d],b.hit_point_m[d]) ||
        !same(a.incident_proper_velocity_m_s[d],b.incident_proper_velocity_m_s[d])) return false;
  return true;
}
inline Species species(const g::Configuration &cfg) {
  Species out;
  for (std::uint32_t k = 0; k < cfg.species.size(); ++k) {
    const auto &s = cfg.species[k];
    need(!s.state.empty() && std::isfinite(s.mass_kg) && s.mass_kg > 0 &&
             out.emplace(s.state, k).second, "invalid species catalogue");
  }
  return out;
}
inline Rows rows(const g::Configuration &cfg, const Species &sp,
                 std::span<const n::ArenaRateEntry> input, std::size_t cap) {
  need(input.size() <= cap, "carrier capacity");
  Rows out;
  for (const auto &r : input) {
    auto it = sp.find(r.state);
    need(it != sp.end() && !r.packet_id.empty(), "unknown species or missing packet");
    const auto &s = cfg.species[it->second];
    const auto &p = r.native;
    need(p.id > 0 && p.weight > 0 && std::isfinite(p.weight) &&
             same(p.mass_kg, s.mass_kg) && same(p.charge_C, s.charge * n::pic_e),
         "carrier mass/charge/weight");
    for (int d = 0; d < 3; ++d)
      need(std::isfinite(p.position_m[d]) && std::isfinite(p.proper_velocity_m_s[d]),
           "nonfinite phase space");
    need(cfg.grid.locate(p.position_m.data()) != pfg::kNone, "carrier outside gas grid");
    need(out.emplace(p.id, &r).second, "duplicate native ID");
  }
  return out;
}
inline bool identity(const n::ArenaRateEntry &a, const n::ArenaRateEntry &b) {
  return a.packet_id == b.packet_id && a.state == b.state && a.native.id == b.native.id &&
         same(a.native.weight, b.native.weight) && same(a.native.mass_kg, b.native.mass_kg) &&
         same(a.native.charge_C, b.native.charge_C);
}
inline bool phase(const n::ArenaRateEntry &a, const n::ArenaRateEntry &b) {
  if (!identity(a, b)) return false;
  for (int d = 0; d < 3; ++d)
    if (!same(a.native.position_m[d], b.native.position_m[d]) ||
        !same(a.native.proper_velocity_m_s[d], b.native.proper_velocity_m_s[d])) return false;
  return true;
}
inline Energy energies(std::span<const n::FullTimeStoredPacketEnergy> input, const Rows &rows,
                       std::size_t cap) {
  need(input.size() <= cap, "internal-stock capacity");
  std::map<std::string, std::size_t> counts;
  std::map<std::string, const n::ArenaRateEntry *> packet;
  for (const auto &[id, r] : rows) { (void)id; ++counts[r->packet_id]; packet[r->packet_id] = r; }
  Energy out;
  for (const auto &e : input) {
    need(std::isfinite(e.energy_J) && e.energy_J >= 0 && std::isfinite(e.original_birth_weight) &&
             e.original_birth_weight > 0 && counts.count(e.packet_id) && counts.at(e.packet_id) == 1,
         "packet internal energy requires one present carrier");
    need(same(packet.at(e.packet_id)->native.weight, e.original_birth_weight) &&
             packet.at(e.packet_id)->native.charge_C == 0,
         "packet internal weight/neutral mismatch");
    need(out.emplace(e.packet_id, &e).second, "duplicate internal-energy packet");
  }
  return out;
}
inline double internal(const Energy &e, const std::string &packet) {
  auto it = e.find(packet); return it == e.end() ? 0. : it->second->energy_J;
}
inline std::vector<g::CensusRow> census_rows(const Rows &rows, const Species &sp) {
  std::vector<g::CensusRow> out; out.reserve(rows.size());
  for (const auto &[id, r] : rows) {
    g::CensusRow x{}; x.native_id = id; x.species = sp.at(r->state); x.weight = r->native.weight;
    for (int d = 0; d < 3; ++d) {
      x.position_m[d] = r->native.position_m[d];
      x.proper_velocity_m_s[d] = r->native.proper_velocity_m_s[d];
    }
    out.push_back(x);
  }
  return out;
}
inline g::FlowRecord flow(const g::Configuration &cfg, const Species &sp, const n::ArenaRateEntry &r,
                          std::uint32_t kind, double internal) {
  g::FlowRecord f{}; f.kind = kind; f.negative = kind == pfg::kFlowSink;
  f.grid_cell = cfg.grid.locate(r.native.position_m.data()); f.species = sp.at(r.state);
  f.native_id = r.native.id; f.weight = r.native.weight; f.internal_J = internal;
  for (int d = 0; d < 3; ++d) f.u[d] = r.native.proper_velocity_m_s[d];
  return f;
}
} // namespace detail

inline Records prepare(const g::Configuration &cfg, const View &v, std::size_t cap) {
  need(cap > 0 && cap <= 1000000 && cfg.transport == g::Transport::carrier_shares,
       "bounded carrier-shares transport required");
  need(!v.has_merge, "merge/split gas transport not yet qualified");
  const auto sp = detail::species(cfg);
  const auto before = detail::rows(cfg, sp, v.before, cap);
  const auto start = detail::rows(cfg, sp, v.start, cap);
  const auto after = detail::rows(cfg, sp, v.after, cap);
  const auto e0 = detail::energies(v.before_internal, before, cap);
  const auto e1 = detail::energies(v.start_internal, start, cap);
  const auto e2 = detail::energies(v.after_internal, after, cap);
  Records out;
  std::set<std::uint64_t> new_ids;
  if (v.births) {
    need(v.births->accepted_step == v.accepted_step && v.births->entries.size() <= cap &&
             v.births->wall_entries.empty(), "birth step/cap or unsupported wall-response births");
    g::Exact internal_sum;
    for (const auto &b : v.births->entries) {
      need(!b.returned, "returned subcell birth requires separate flow policy");
      auto it = start.find(b.native_id);
      need(!before.count(b.native_id) && it != start.end() && new_ids.insert(b.native_id).second,
           "birth native ID partition");
      const auto &r = *it->second;
      need(b.packet_id == r.packet_id && b.state == r.state && same(b.weight, r.native.weight),
           "native birth receipt differs from actual start carrier");
      const double U = detail::internal(e1, r.packet_id);
      internal_sum.add_product(U, 1);
      out.births.push_back(detail::flow(cfg, sp, r, pfg::kFlowBirth, U));
    }
    need(same(internal_sum.value(), v.births->internal_energy_J), "birth internal-energy sum");
    if (!v.births->internal_energy_J_digits.empty()) {
      need(detail::native_digits(internal_sum) == v.births->internal_energy_J_digits, "birth exact internal-energy sum");
    } else need(internal_sum.zero(), "missing exact birth internal-energy stock");
  }
  need(start.size() == before.size() + new_ids.size(), "unreported source birth/debit");
  for (const auto &[id, r] : before) {
    auto it = start.find(id);
    need(it != start.end() && detail::phase(*r, *it->second), "pre-motion carrier changed outside births");
    need(same(detail::internal(e0, r->packet_id), detail::internal(e1, r->packet_id)),
         "pre-motion carried internal energy changed");
  }
  std::vector<std::uint64_t> sunk;
  std::set<std::uint64_t> sinks;
  if (v.sinks) {
    need(v.sinks->accepted_step == v.accepted_step && v.sinks->records.size() <= cap &&
             v.sinks->records.size() == v.sinks->exposed.size() &&
             v.sinks->carriers == v.sinks->records.size(), "sink step/cardinality");
    for (std::size_t i = 0; i < v.sinks->exposed.size(); ++i) {
      const auto &x = v.sinks->exposed[i]; const auto &k = x.record;
      need(detail::same_sink(k,v.sinks->records[i]) && sinks.insert(k.native_id).second &&
               !after.count(k.native_id), "sink native ID partition");
      auto it = start.find(k.native_id); need(it != start.end(), "sink missing from actual start");
      const auto &r = *it->second;
      need(r.packet_id == x.packet_id && r.state == x.state && same(r.native.weight, x.weight) &&
               same(r.native.mass_kg, x.mass_kg) && same(r.native.charge_C, x.charge_C),
           "sink exposure species/weight");
      need(same(detail::internal(e1, r.packet_id), x.internal_energy_J), "sink internal-energy stock");
      need(k.kind == n::FullTimeSinkKind::collected_at_conductor ||
               k.kind == n::FullTimeSinkKind::escaped_exterior, "unknown sink destination");
      const bool conductor = k.surface == n::FullTimeSinkSurface::lower_foil ||
                             k.surface == n::FullTimeSinkSurface::upper_foil ||
                             k.surface == n::FullTimeSinkSurface::short_bridge;
      const bool exterior = k.surface == n::FullTimeSinkSurface::outer_xy ||
                            k.surface == n::FullTimeSinkSurface::axial;
      need((conductor || exterior) && conductor == (k.kind == n::FullTimeSinkKind::collected_at_conductor) &&
               k.axis < 3 && std::isfinite(k.plane_m) && std::isfinite(k.fraction) &&
               k.fraction >= 0 && k.fraction <= 1 && std::isfinite(k.flight_s) && k.flight_s >= 0 &&
               std::isfinite(k.discarded_flight_s) && k.discarded_flight_s >= 0,
           "sink geometry/destination metadata");
      for (int d = 0; d < 3; ++d)
        need(same(r.native.position_m[d], x.start_position_m[d]) &&
                 std::isfinite(k.incident_proper_velocity_m_s[d]) && std::isfinite(k.segment_m[d]) &&
                 std::isfinite(k.hit_point_m[d]), "sink start/incident phase space");
      auto f = detail::flow(cfg, sp, r, pfg::kFlowSink, x.internal_energy_J);
      // Sink K/P are evaluated at the actual incident velocity, not at the
      // beginning of the field push. Grid attribution uses the origin cell.
      for (int d = 0; d < 3; ++d) f.u[d] = k.incident_proper_velocity_m_s[d];
      out.sinks.push_back(f); out.sink_destinations.push_back(k); sunk.push_back(k.native_id);
    }
  }
  need(after.size() + sinks.size() == start.size(), "unreported transport birth/sink");
  for (const auto &[id, r] : after) {
    auto it = start.find(id);
    need(it != start.end() && detail::identity(*r, *it->second), "transport changed species/weight/packet");
    need(same(detail::internal(e1, r->packet_id), detail::internal(e2, r->packet_id)),
         "transport changed carried internal energy");
  }
  const auto from = detail::census_rows(start, sp), to = detail::census_rows(after, sp);
  out.movers = g::movers(cfg, from, to, sunk);
  std::sort(out.births.begin(), out.births.end(), [](const auto &a, const auto &b) { return a.native_id < b.native_id; });
  // Native sink arrays are canonical. Refuse reordered arrays rather than
  // separating the one-to-one destination and flow correspondence.
  need(std::is_sorted(out.sinks.begin(), out.sinks.end(), [](const auto &a, const auto &b) { return a.native_id < b.native_id; }),
       "native sink order");
  return out;
}
} // namespace native_gas_motion
