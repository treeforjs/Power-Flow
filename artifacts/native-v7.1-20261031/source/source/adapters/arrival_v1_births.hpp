#pragma once
// Native v6 <-> arrival owner v1 (runtimes/pic-arrival-owner-v1-20261030,
// API-v1.txt section 8) adapter. Header only; include after
// raikou/mfem_pic_full_time.h and pf_arrival/arrival_owner.hpp (built with
// PF_ARRIVAL_EXACT_COUNT_HEADER="raikou/exact_particle_count.h", so
// pf_arrival::Count and native ArenaCount are the same type).
//  * credit():   Staged step (canonical births, IDs B_s + rank, bound times)
//                + per-birth source kinematics + per-stream exact source
//                debits -> FullTimeBirthCredit for BirthScope::prepare.
//  * receipt():  native FullTimeBirthReceipt -> pf_arrival::NativeBirthReceipt
//                for Owner::commit.
//  * require_restart_alignment(): native checkpoint values (next ID, last
//                nonempty birth block, cumulative births) against a restored
//                arrival owner (contract (f)).
#include "raikou/mfem_pic_full_time.h"
#include "pf_arrival/arrival_owner.hpp"
#include <stdexcept>
#include <string>
#include <vector>

namespace raikou::mfem_backend::arrival_v1 {
static_assert(std::is_same_v<pf_arrival::Count, ArenaCount>,
              "arrival and native exact counts must be the same type");

// Source-owner kinematics of one birth, in the staged canonical order.
struct BirthKinematics {
  PicVector position_m{}, proper_velocity_m_s{};
  std::string origin; // the source event (e.g. source v5 event id / patch)
};
// Exact source debit of one arrival stream (neutral: ReleaseQuote released =
// source v5 Trial inventory decrease; electron: sum of candidate weights).
struct StreamDebit {
  std::uint32_t stream_id = 0;
  ArenaCount debit;
};

inline FullTimeBirthCredit credit(const pf_arrival::Policy &P,
                                  const pf_arrival::Staged &s,
                                  const std::vector<BirthKinematics> &kinematics,
                                  const std::vector<StreamDebit> &debits,
                                  const std::string &source_before_sha256,
                                  const std::string &source_after_sha256) {
  if (kinematics.size() != s.births().size())
    throw std::invalid_argument("one source kinematics row per staged birth");
  FullTimeBirthCredit c;
  c.arrival_step = s.plan().step();
  c.base_native_id = s.native_base();
  c.plan_sha256 = s.plan().sha256();
  c.staged_token_sha256 = s.token();
  c.source_before_sha256 = source_before_sha256;
  c.source_after_sha256 = source_after_sha256;
  for (const auto &d : debits)
    c.streams.push_back({d.stream_id, P.stream(d.stream_id).species_state, d.debit});
  for (std::size_t k = 0; k < s.births().size(); ++k) {
    const auto &b = s.births()[k];
    FullTimeBirth x;
    x.step = b.step;
    x.stream_id = b.stream_id;
    x.index = b.index;
    x.packet_id = b.packet_id();
    x.state = P.stream(b.stream_id).species_state;
    x.origin = kinematics[k].origin;
    x.weight = b.weight;
    x.position_m = kinematics[k].position_m;
    x.proper_velocity_m_s = kinematics[k].proper_velocity_m_s;
    x.birth_time_s = b.time_s;
    x.remaining_flight_s = b.remaining_s;
    x.native_id = b.native_id;
    c.births.push_back(std::move(x));
  }
  return c;
}
// Every decided stream of the plan with its exact debit (Decision::debit is
// the exact total of the births' weights = the source debit the arrival owner
// partitioned); callers that hold a source v5 Trial pass the source-side
// values instead (test_v6_source_arrival does both and requires equality).
inline std::vector<StreamDebit> plan_debits(const pf_arrival::StepPlan &plan) {
  std::vector<StreamDebit> out;
  for (const auto &d : plan.decisions())
    out.push_back({d.stream_id, d.debit});
  return out;
}
inline pf_arrival::NativeBirthReceipt receipt(const FullTimeBirthReceipt &r) {
  pf_arrival::NativeBirthReceipt out;
  out.base = r.base;
  out.next_id_after = r.next_id_after;
  for (const auto &e : r.entries)
    out.entries.push_back({e.packet_id, e.state, e.native_id, e.weight});
  return out;
}
// Contract (f): after restoring the arrival owner from its checkpoint next to
// a native checkpoint, both must describe the same accepted step, the same
// cumulative birth ordinal and the same last nonempty birth block; the next
// birth base B_s is then read from native (next_native_id).
inline void require_restart_alignment(const FullTimeFlowTotals &native,
                                      const pf_arrival::Owner &arrival,
                                      std::uint64_t arrival_first_step = 0,
                                      std::uint64_t native_first_step = 0) {
  auto need = [](bool ok, const char *what) {
    if (!ok)
      throw std::invalid_argument(std::string("arrival/native restart: ") + what);
  };
  need(arrival.next_step() - arrival_first_step ==
           native.accepted_step - native_first_step,
       "accepted step differs");
  need(arrival.next_ordinal() == native.births, "birth ordinal differs");
  const auto a = arrival.last_anchor();
  need(arrival.has_anchor() == native.has_last_birth, "anchor presence differs");
  if (arrival.has_anchor())
    need(a[0] == native.last_birth_step && a[1] == native.last_birth_base &&
             a[2] == native.last_birth_count,
         "last birth block differs");
  need(!native.has_last_birth ||
           native.next_native_id >= native.last_birth_base + native.last_birth_count,
       "native next ID precedes the last birth block");
}
} // namespace raikou::mfem_backend::arrival_v1
