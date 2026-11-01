// Native v8.1 (stage B) controls (API-008.txt 8): wall-response births (source owner v7 B2),
// matching-plane births with bridge paths (B3), born-and-collected births (B4), internal energy
// carried by metal-vapour macros (phase owner v1.1 3.2: published with the sink, released when
// chemistry or a transaction removes the packet), the merge/split transaction (merge-split
// INTEGRATION 2: native recomputation, outputs, relocation records in the v8.1 motion record,
// exact residuals) and the restore of all of it. With NVRTC and UUID the same sequence also runs in
// device_verified and device modes and must give the same accepted identities every step.
// usage: test_v8_stage_b OUT_DIR [NVRTC UUID]
#include "v6_fixture.hpp"
#include "arrival_v1_births.hpp"
#include "forms_gpu.hpp"
#include "hall_gpu_common.hpp"
#include "merge_split.hpp"
#include "raikou/mfem_pic_full_time_device.h"
#include <filesystem>
#include <fstream>
#include <iostream>

using namespace v6;
namespace pa = pf_arrival;
namespace av = m::arrival_v1;
using FK = m::FullTimeFlowRefusal::Kind;
namespace {
int n_checks = 0, n_fail = 0;
J report = J::object();
void need(bool ok, const std::string &what, const J &d = J::object()) {
  ++n_checks;
  if (!ok) {
    ++n_fail;
    std::cout << "FAIL " << what << " " << d.dump() << "\n";
  }
}
template <class F> std::pair<int, std::string> refusal_of(F f) {
  try {
    f();
  } catch (const m::FullTimeFlowRefusal &e) {
    return {int(e.kind()), e.what()};
  } catch (const std::exception &e) {
    return {-2, e.what()};
  }
  return {-1, "none"};
}
std::shared_ptr<const pa::Policy> arrival_policy() {
  pa::Configuration c;
  c.run_seed = 20261101;
  c.streams = {{0, pa::Kind::electron, 0, "e", 4, 4096}, {1, pa::Kind::neutral, 0, "n", 2, 4096}};
  return std::make_shared<const pa::Policy>(c);
}
void bind_neutral_times(pa::Staged &s) {
  for (const auto &d : s.plan().decisions())
    if (d.kind == pa::Kind::neutral && d.births) {
      std::vector<double> t;
      for (const auto &b : s.births())
        if (b.stream_id == d.stream_id)
          t.push_back(d.begin_s + b.marks[5] * (d.end_s - d.begin_s));
      s.bind_source_times(d.stream_id, t);
    }
}
// Placement: 0 emitting face (upper foil lower face y = .007, into the gap); 1 sink-prone
// (neutrals just above the lower foil top face moving -y); 2 interior.
std::vector<av::BirthKinematics> placement(const pa::Staged &s, int kind) {
  std::vector<av::BirthKinematics> out;
  for (const auto &b : s.births()) {
    av::BirthKinematics k;
    const double dx = 1e-5 * double(b.index % 100), dz = 1e-5 * double(b.index / 100);
    const bool e = b.kind == pa::Kind::electron;
    if (kind == 0) {
      k.position_m = {.0030 + dx, .007, -.001 + dz + (e ? 0 : .002)};
      k.proper_velocity_m_s = {0, e ? -2e6 : -650, 0};
    } else if (kind == 1) {
      k.position_m = e ? m::PicVector{.004 + dx, .0045, .001 + dz} : m::PicVector{.004 + dx, .002 + 1e-14, dz};
      k.proper_velocity_m_s = e ? m::PicVector{1e6, 0, 0} : m::PicVector{0, -650, 0};
    } else {
      k.position_m = {.004 + dx, .0045, .001 + dz + (e ? 0 : .001)};
      k.proper_velocity_m_s = {e ? 1e6 : 300., 0, 0};
    }
    k.origin = std::string(e ? "source/electron-emission/event-" : "source/neutral-release/event-") +
               std::to_string(b.ordinal + 1);
    out.push_back(k);
  }
  return out;
}
// Policy: species 0 = "e", 1 = "n" on the 45 x 45 Hall mesh (4 axial cells), one energy bin, no
// direction bins; merge above 4 per (cell, species), split below 2.
std::shared_ptr<const MsParams> merge_policy() {
  auto P = std::make_shared<MsParams>(merge_split::hall_mesh_params(45, 4, .009, 1, .025));
  P->seed = 20261101;
  P->nsub = 0;
  merge_split::set_species(*P, 0, m::pic_me, -m::pic_e);
  merge_split::set_species(*P, 1, 2.9915e-26, 0.0);
  for (unsigned s = 0; s < 2; ++s) {
    auto &sp = P->species[s];
    sp.n_thresh = 0;
    sp.merge_above = 4;
    sp.merge_target = 2;
    sp.split_below = 2;
    sp.split_target = 2;
    sp.k_min = 3;
    sp.k_max = 8;
    sp.split_levels = 1;
    sp.split_min_weight = 1.0;
  }
  return P;
}
std::vector<m::FullTimeAuthoredPacket> population() {
  std::vector<m::FullTimeAuthoredPacket> p = {
      // sinks onto the lower foil's top face in the first step (parent of wall-response births)
      {"e-sink", "e", "manufactured/source", {4, 0}, {.004, .0020001, 0}, {0, -3e7, 0}},
      {"neutral", "n", "manufactured/source", {8, 0}, {.00403, .004, 0}, {0, 0, 0}}};
  // ten electrons in one Hall cell (merge candidates), different velocities
  for (int k = 0; k < 10; ++k)
    p.push_back({"m-" + std::to_string(k), "e", "manufactured/source", {2, 0},
                 {.00451 + 1e-6 * k, .00451 + 2e-6 * k, .00101 + 1e-5 * k},
                 {1e5 * (k % 3 - 1), 2e5 * (k % 2), 3e4 * (k - 4)}});
  // a lone electron in another cell (split candidate)
  p.push_back({"lone", "e", "manufactured/source", {8, 0}, {.0060, .0050, -.004}, {1e4, 0, 0}});
  return p;
}
std::unique_ptr<m::PicArena> fresh(bool v8 = true, std::uint64_t cadence = 3) {
  m::FullTimeConfiguration c;
  c.identity = v8 ? m::FullTimeIdentitySchema::v8 : m::FullTimeIdentitySchema::v5;
  c.limits.maximum_carriers = 1024;
  c.limits.memory_budget_bytes = 1024 * m::full_time_bytes_per_carrier + m::full_time_fixed_bytes;
  c.history.anchor_interval_steps = 4;
  if (v8) {
    c.merge_split.enabled = true;
    c.merge_split.cadence_steps = cadence;
    c.merge_split.params = merge_policy();
    c.merge_split.species_states = {"e", "n"};
  }
  return m::PicArena::create_full_time(options(), states(), channels(), population(), genesis(), c);
}
// v8.1 motion record producer: the reference flows record with the extras of the step.
struct Extras {
  std::shared_ptr<ReferenceMotion> ref;
  m::FullTimeMotionRecordV7 record;
  std::vector<m::FullTimeMotionRelocation> reloc;
  std::vector<m::FullTimeBirthBridge> bridges;
  std::vector<double> xJ;
};
using Mutate = std::function<void(Extras &)>;
m::PreparedPrivateMotion v81_motion(m::PrivateTransportTrial &trial, m::FullTimeMotionBinding &binding,
                                    const std::vector<m::FullTimeMotionRelocation> &reloc,
                                    const std::vector<m::FullTimeBirthBridge> &bridges, bool v81, const Mutate &mutate) {
  auto h = std::make_shared<Extras>();
  h->ref = reference_motion(trial, binding.common, options(), 1024, trial.census_sha256(), true);
  h->record = h->ref->record;
  if (v81) {
    h->record.schema = m::motion_v81_schema;
    h->reloc = reloc;
    h->bridges = bridges;
    if (!reloc.empty() || !bridges.empty()) {
      h->xJ.assign(zero_field().J.size(), 0.0);
      h->record.extra_current_A_m_sha256 = pf_v7::array_sha256(h->xJ);
      h->record.extra_work_J = 0;
    }
  }
  if (mutate)
    mutate(*h);
  h->record.relocations = h->reloc;
  h->record.bridges = h->bridges;
  h->record.extra_current_A_m = h->xJ;
  h->record.record_sha256 =
      v81 ? m::full_time_motion_record_v81_sha256(h->record) : m::full_time_motion_record_v7_sha256(h->record);
  binding.transport_trial_sha256 = h->record.record_sha256;
  return trial.prepare_total_field_motion(binding, m::bind_total_field_motion_v7(h->record, h));
}
// One outer step with optional merge/split, arrival + wall-response births and the v8.1 record.
struct Plan {
  int place = 2;
  double electrons = 12, neutrals = 3;
  bool merge = false, v81 = true;
  std::string expected_digest;
  std::function<void(m::FullTimeBirthCredit &, const m::PicArena &, std::uint64_t base)> customize;
  Mutate mutate;
};
struct Out {
  std::shared_ptr<const m::FullTimeStepFlows> flows;
  std::optional<m::FullTimeMergeSplitReceipt> merge;
  std::optional<m::FullTimeBirthReceipt> births;
};
struct Run {
  std::shared_ptr<const pa::Policy> P = arrival_policy();
  std::unique_ptr<pa::Owner> owner = std::make_unique<pa::Owner>(P, 0);
  std::uint64_t source_index = 0;
  StepState st;
  Out step(m::PicArena &a, const Plan &plan) {
    const auto before = a.full_time_boundary();
    const std::uint64_t k = before->accepted_step;
    auto b = binding_of(a, 5e-15, sha(step_text("v81-clock", k)));
    auto outer = a.begin_full_time_step(b, 1024, 32);
    m::ExposureBinding e1{b, m::FullTimeSlot::before_transport, b.begin_s, .5 * b.dt_s, st.sequence,
                          b.accepted_field_sha256};
    const auto rng_a = sha(step_text("v81-rng-first", k)), rng_b = sha(step_text("v81-rng-second", k));
    auto s1 = outer.first_exposure(e1);
    auto t1 = s1.prepare_compound_with_completion(1, [&](m::CollisionScope &) { return completion(e1, st.rng, rng_a); });
    if (!s1.try_publish(t1))
      throw std::runtime_error("first exposure publication");
    Out out;
    std::optional<m::PreparedMergeSplit> ms;
    if (plan.merge) {
      auto scope = outer.merge_scope();
      ms.emplace(scope.prepare(plan.expected_digest));
      out.merge = ms->receipt();
    }
    auto scope = ms ? outer.birth_scope(*ms) : outer.birth_scope();
    pa::StepInput in;
    in.step = k;
    in.electrons = {{0, plan.electrons, std::log(plan.electrons), true, b.begin_s, b.dt_s, b.end_s}};
    in.neutrals = {{1, pa::Count(plan.neutrals), std::log(plan.neutrals), pa::Count(1e20), b.begin_s, b.end_s}};
    auto staged = owner->stage(owner->plan(in), scope.base_native_id());
    bind_neutral_times(staged);
    auto credit = av::credit(*P, staged, placement(staged, plan.place), av::plan_debits(staged.plan()),
                             scope.source_before_sha256(), sha("source/" + std::to_string(source_index + 1)));
    if (plan.customize)
      plan.customize(credit, a, scope.base_native_id());
    auto births = scope.prepare(credit);
    out.births = births.receipt();
    auto trial = outer.fork_transport_trial(births);
    const auto field_after = sha(step_text("v81-field-after", k)), heat_after = sha(step_text("v81-heat-after", k));
    m::FullTimeMotionBinding mb{b, "", b.accepted_field_sha256, field_after, heat_after, credit.source_after_sha256};
    static const std::vector<m::FullTimeMotionRelocation> none;
    auto moved = v81_motion(trial, mb, ms ? ms->receipt().relocation_records : none, births.receipt().bridges,
                            plan.v81, plan.mutate);
    outer.accept_transport_trial(std::move(moved));
    m::ExposureBinding e2 = e1;
    e2.slot = m::FullTimeSlot::after_transport;
    e2.frozen_physical_stage_s = b.end_s;
    e2.owner_stage_sequence = st.sequence + 1;
    e2.frozen_field_sha256 = field_after;
    auto s2 = outer.second_exposure(e2);
    auto t2 = s2.prepare_compound_with_completion(1, [&](m::CollisionScope &) { return completion(e2, rng_a, rng_b); });
    if (!s2.try_publish(t2))
      throw std::runtime_error("second exposure publication");
    auto token = outer.seal({field_after, heat_after, credit.source_after_sha256, rng_b});
    out.flows = token.candidate_full_time_flows();
    if (!a.can_publish(token) || !a.try_publish(token))
      throw std::runtime_error("outer publication refused");
    (void)owner->commit(std::move(staged), av::receipt(*out.flows->births));
    ++source_index;
    st.sequence += 2;
    st.rng = rng_b;
    return out;
  }
};
// Wall-response births from the previous step's published sinks: an electron (outward velocity)
// and a desorbed neutral per collected parent on a foil face.
void add_wall_births(m::FullTimeBirthCredit &c, const m::PicArena &a, std::uint64_t base) {
  const auto f = a.full_time_published_flows();
  if (!f || !f->sinks)
    return;
  std::uint64_t id = base + c.births.size();
  for (const auto &e : f->sinks->exposed) {
    if (e.record.kind != m::FullTimeSinkKind::collected_at_conductor || e.record.surface == m::FullTimeSinkSurface::short_bridge)
      continue;
    const double out_sign = e.record.surface == m::FullTimeSinkSurface::lower_foil ? 1.0 : -1.0;
    for (std::uint32_t k = 0; k < 2; ++k) {
      m::FullTimeWallBirth w;
      w.wall_response_sha256 = sha("wall-response/" + std::to_string(e.record.native_id));
      w.parent_native_id = e.record.native_id;
      w.ordinal = k;
      w.packet_id = m::full_time_wall_birth_packet_id(a.full_time_boundary()->accepted_step, w.parent_native_id, k);
      w.state = k == 0 ? "e" : "n";
      w.origin = k == 0 ? "surface-v7-wall-response/true-secondary" : "surface-v7-wall-response/stuck-release";
      w.weight = e.weight;
      w.position_m = e.record.hit_point_m;
      w.proper_velocity_m_s = {0, 0, 0};
      w.proper_velocity_m_s[e.record.axis] = out_sign * (k == 0 ? 1e6 : 500.0);
      w.parent_hit_time_s = a.full_time_boundary()->time_s - 5e-15 + e.record.flight_s;
      w.native_id = id++;
      c.wall.push_back(w);
    }
  }
}
// The first electron birth becomes a matching-plane birth (face point + 1 um along the outward
// normal), the second a returned (born-and-collected) birth; the first neutral carries internal energy.
void mark_v81(m::FullTimeBirthCredit &c, const m::PicArena &, std::uint64_t) {
  int electrons = 0;
  bool neutral = false;
  for (auto &b : c.births) {
    if (b.state == "e") {
      if (electrons == 0) {
        b.bridge = true;
        b.bridge_start_m = b.position_m;
        b.position_m[1] -= 1e-6;
      } else if (electrons == 1)
        b.returned = true;
      ++electrons;
    } else if (!neutral) {
      b.internal_energy_J = 3.2e-20;
      neutral = true;
    }
  }
}
std::string trace_of(const m::PicArena &a) {
  return a.full_time_boundary()->accepted_arena_sha256 + "/" + a.full_time_flow_totals().flows_sha256 + "/" +
         a.full_time_flows_v81().sha256;
}
// The main sequence: 0 sink step; 1 wall-response births + bridge/returned/internal births on the
// emitting face; 2 sink-prone births carrying internal energy (sunk in their first flight); 3
// merge/split transaction (+ wall births of step 2's sinks); 4 interior births.
std::vector<std::string> sequence(m::PicArena &a, Run &R, J *log) {
  std::vector<std::string> trace;
  auto both = [](m::FullTimeBirthCredit &c, const m::PicArena &x, std::uint64_t base) {
    mark_v81(c, x, base);
    add_wall_births(c, x, base);
  };
  Plan p0;
  p0.place = 2;
  auto o0 = R.step(a, p0);
  trace.push_back(trace_of(a));
  Plan p1;
  p1.place = 0;
  p1.customize = both;
  auto o1 = R.step(a, p1);
  trace.push_back(trace_of(a));
  Plan p2;
  p2.place = 1;
  p2.customize = [](m::FullTimeBirthCredit &c, const m::PicArena &x, std::uint64_t base) {
    for (auto &b : c.births)
      if (b.state == "n")
        b.internal_energy_J = 1.6e-20;
    add_wall_births(c, x, base);
  };
  auto o2 = R.step(a, p2);
  trace.push_back(trace_of(a));
  Plan p3;
  p3.place = 2;
  p3.merge = true;
  p3.customize = add_wall_births;
  auto o3 = R.step(a, p3);
  trace.push_back(trace_of(a));
  Plan p4;
  p4.place = 2;
  auto o4 = R.step(a, p4);
  trace.push_back(trace_of(a));
  (void)o0;
  (void)o4;
  if (log) {
    const auto &b1 = *o1.births;
    (*log)["step1"] = {{"wall_entries", b1.wall_entries.size()}, {"bridges", b1.bridges.size()},
                       {"returned", std::count_if(b1.entries.begin(), b1.entries.end(), [](auto &e) { return e.returned; })},
                       {"internal_J", b1.internal_energy_J}, {"next_id_after", b1.next_id_after},
                       {"arrival_next_id_after", b1.arrival_next_id_after}};
    J sunk = J::array();
    if (o2.flows && o2.flows->sinks)
      for (const auto &e : o2.flows->sinks->exposed)
        if (e.internal_energy_J != 0)
          sunk.push_back({{"native_id", e.record.native_id}, {"state", e.state}, {"internal_J", e.internal_energy_J}});
    (*log)["step2_sunk_internal"] = sunk;
    const auto &mr = *o3.merge;
    (*log)["step3_merge"] = {{"outputs", mr.outputs}, {"groups", mr.groups}, {"accepted_merge", mr.accepted_merge},
                             {"accepted_split", mr.accepted_split}, {"inputs", mr.inputs_consumed},
                             {"relocations", mr.relocation_records.size()}, {"transaction", mr.transaction_sha256},
                             {"base", mr.base_native_id}, {"birth_base", o3.births->base}};
    J res = J::array();
    for (const auto &[st, r] : mr.residuals)
      res.push_back({{"state", st}, {"energy_J", r[3].value}, {"energy_units", r[3].digits},
                     {"momentum_x", r[0].value}});
    (*log)["step3_residuals"] = res;
  }
  return trace;
}
} // namespace

int main(int argc, char **argv) {
  try {
    if (argc != 2 && argc != 4)
      throw std::runtime_error("usage: test_v8_stage_b OUT_DIR [NVRTC UUID]");
    const std::filesystem::path out = argv[1];
    std::filesystem::create_directories(out);
    // ---- main sequence (host)
    auto a = fresh();
    Run R;
    J log = J::object();
    const auto host_trace = sequence(*a, R, &log);
    report["sequence"] = log;
    const auto st = a->full_time_flows_v81();
    need(log["step1"]["wall_entries"] == 2 && log["step1"]["bridges"] == 1 && log["step1"]["returned"] == 1 &&
             log["step1"]["internal_J"] == 3.2e-20,
         "births: two wall-response births of the sunk parent, one bridge, one returned, internal energy", log["step1"]);
    need(log["step1"]["arrival_next_id_after"].get<std::uint64_t>() + 2 == log["step1"]["next_id_after"].get<std::uint64_t>(),
         "wall-response births follow the arrival births (arrival receipt block unchanged)", log["step1"]);
    need(!log["step2_sunk_internal"].empty(), "sunk macros publish their internal energy", log["step2_sunk_internal"]);
    need(log["step3_merge"]["accepted_merge"].get<std::uint64_t>() >= 1 && log["step3_merge"]["relocations"].get<std::uint64_t>() >= 3 &&
             log["step3_merge"]["accepted_split"].get<std::uint64_t>() >= 1 &&
             log["step3_merge"]["birth_base"].get<std::uint64_t>() ==
                 log["step3_merge"]["base"].get<std::uint64_t>() + log["step3_merge"]["outputs"].get<std::uint64_t>(),
         "merge/split: merge and split groups, relocation records, B_s = M_s + n_out", log["step3_merge"]);
    // exact internal-energy closure: born = sunk + released + in flight
    {
      const double b = st.internal_born_J.value, s = st.internal_sunk_J.value, r = st.internal_released_J.value,
                   f = st.internal_in_flight_J.value;
      need(st.present && b > 0 && s > 0 && std::abs(b - (s + r + f)) <= 1e-12 * b,
           "internal energy: born = sunk + released + in flight",
           {{"born", b}, {"sunk", s}, {"released", r}, {"in_flight", f}});
    }
    need(st.wall_births >= 2 && st.returned_births == 1 && st.bridges == 1 && st.merge_transactions == 1 &&
             st.merge_outputs == log["step3_merge"]["outputs"].get<std::uint64_t>(),
         "v8.1 cumulative flows (wall, returned, bridges, merge)");
    a->validate_full_time_history();
    // ---- restore
    const auto saved = *a->full_time_boundary();
    a->save((out / "stage-b.json").string());
    auto r = m::PicArena::restore_full_time((out / "stage-b.json").string(), options(), states(), channels(), saved);
    need(r->full_time_boundary()->accepted_arena_sha256 == saved.accepted_arena_sha256 &&
             r->full_time_flows_v81().sha256 == st.sha256,
         "restore: accepted identity and v8.1 flows state");
    // ---- refusals (each on a fresh arena at the state after step 0 or 2; accepted state unchanged)
    auto at = [&](int steps, Run &q) {
      auto x = fresh();
      Plan p0;
      p0.place = 2;
      (void)q.step(*x, p0);
      if (steps > 1) {
        Plan p1;
        p1.place = 0;
        p1.customize = [](m::FullTimeBirthCredit &c, const m::PicArena &y, std::uint64_t base) {
          mark_v81(c, y, base);
          add_wall_births(c, y, base);
        };
        (void)q.step(*x, p1);
        Plan p2;
        p2.place = 1;
        p2.customize = add_wall_births;
        (void)q.step(*x, p2);
      }
      return x;
    };
    struct Case {
      const char *name;
      int steps;
      int kind;
      std::function<void(Plan &)> plan;
    };
    const std::vector<Case> cases = {
        {"wall birth weight differs from the parent", 1, int(FK::birth_count),
         [](Plan &p) {
           p.customize = [](m::FullTimeBirthCredit &c, const m::PicArena &x, std::uint64_t base) {
             add_wall_births(c, x, base);
             if (!c.wall.empty())
               c.wall[0].weight = std::nextafter(c.wall[0].weight, 1e9);
           };
         }},
        {"wall birth position differs from the hit point", 1, int(FK::birth_kinematics),
         [](Plan &p) {
           p.customize = [](m::FullTimeBirthCredit &c, const m::PicArena &x, std::uint64_t base) {
             add_wall_births(c, x, base);
             if (!c.wall.empty())
               c.wall[0].position_m[0] += 1e-9;
           };
         }},
        {"wall birth of an unpublished parent", 1, int(FK::birth_binding),
         [](Plan &p) {
           p.customize = [](m::FullTimeBirthCredit &c, const m::PicArena &x, std::uint64_t base) {
             add_wall_births(c, x, base);
             if (!c.wall.empty()) {
               c.wall.resize(1);
               c.wall[0].parent_native_id = 2; // a resident carrier, not a sink of the previous step
               c.wall[0].packet_id = m::full_time_wall_birth_packet_id(x.full_time_boundary()->accepted_step, 2, 0);
             }
           };
         }},
        {"wall birth ordinals not 0..k-1", 1, int(FK::birth_order),
         [](Plan &p) {
           p.customize = [](m::FullTimeBirthCredit &c, const m::PicArena &x, std::uint64_t base) {
             add_wall_births(c, x, base);
             if (c.wall.size() > 1)
               c.wall.erase(c.wall.begin());
           };
         }},
        {"bridge not along the outward face normal", 1, int(FK::birth_kinematics),
         [](Plan &p) {
           p.place = 0;
           p.customize = [](m::FullTimeBirthCredit &c, const m::PicArena &x, std::uint64_t base) {
             mark_v81(c, x, base);
             for (auto &b : c.births)
               if (b.bridge)
                 b.position_m[0] += 1e-7;
           };
         }},
        {"returned birth off the emitting face", 1, int(FK::birth_kinematics),
         [](Plan &p) {
           p.place = 2;
           p.customize = [](m::FullTimeBirthCredit &c, const m::PicArena &, std::uint64_t) {
             for (auto &b : c.births)
               if (b.state == "e") {
                 b.returned = true;
                 break;
               }
           };
         }},
        {"bridge births with a v7 motion record", 1, int(FK::motion_record),
         [](Plan &p) {
           p.place = 0;
           p.customize = mark_v81;
           p.v81 = false;
         }},
        {"merge/split on a non-transaction step", 1, int(FK::merge_scope), [](Plan &p) { p.merge = true; }},
        {"merge/split producer digest differs", 3, int(FK::merge_binding),
         [](Plan &p) {
           p.merge = true;
           p.expected_digest = sha("foreign transaction");
         }},
        {"merge/split relocations differ in the motion record", 3, int(FK::motion_record),
         [](Plan &p) {
           p.merge = true;
           p.mutate = [](Extras &x) {
             if (!x.reloc.empty())
               x.reloc[0].field_displacement_m[0] = std::nextafter(x.reloc[0].field_displacement_m[0], 1.0);
           };
         }},
        {"merge/split step with a v7 motion record", 3, int(FK::motion_record),
         [](Plan &p) {
           p.merge = true;
           p.v81 = false;
         }}};
    J refusals = J::array();
    for (const auto &c : cases) {
      Run q;
      auto x = at(c.steps, q);
      const auto before = x->state_prefix();
      Plan p;
      c.plan(p);
      const auto [kind, text] = refusal_of([&] { (void)q.step(*x, p); });
      refusals.push_back({{"case", c.name}, {"kind", kind}, {"text", text}});
      need(kind == c.kind && x->state_prefix() == before, std::string("refusal: ") + c.name, refusals.back());
    }
    report["refusals"] = refusals;
    {
      auto v5 = fresh(false);
      Run q;
      Plan p;
      p.place = 0;
      p.customize = mark_v81;
      const auto [kind, text] = refusal_of([&] { (void)q.step(*v5, p); });
      need(kind == int(FK::birth_kinematics), "schema v5 refuses v8.1 birth fields", {{"text", text}});
    }
    // ---- device modes
    if (argc == 4) {
      hall_gpu::SessionOptions so;
      so.require_uuid = argv[3];
      so.nvrtc_library = argv[2];
      auto session = std::make_unique<hall_gpu::Session>(so);
      auto forms = std::make_unique<pf_v7::gpu::Forms>(*session);
      auto dev = std::make_shared<m::FullTimeDevice>(*session, *forms);
      for (auto mode : {m::FullTimeExecution::device_verified, m::FullTimeExecution::device}) {
        auto d = fresh();
        d->attach_full_time_device(dev, mode);
        Run q;
        const auto c0 = dev->counters();
        const auto tr = sequence(*d, q, nullptr);
        const auto c1 = dev->counters();
        const char *name = mode == m::FullTimeExecution::device ? "device" : "device_verified";
        report[std::string("counters_") + name] = {{"device_admissions", c1.device_admissions - c0.device_admissions},
                                                   {"host_admissions", c1.host_admissions - c0.host_admissions},
                                                   {"verified", c1.verified_comparisons - c0.verified_comparisons}};
        need(tr == host_trace, std::string(name) + ": every step equals the host run (identity, flows, v8.1 flows)");
        need(c1.device_admissions - c0.device_admissions == host_trace.size() && c1.host_admissions == c0.host_admissions,
             std::string(name) + ": every admission on the device lane (merge, wall, bridge, returned)",
             report[std::string("counters_") + name]);
      }
      dev.reset();
      forms.reset();
      (void)session->close();
    }
    report["checks"] = n_checks;
    report["failures"] = n_fail;
    std::ofstream(out / "v8-stage-b.json") << report.dump(2) << "\n";
    std::cout << (n_fail ? "FAILED " : "") << n_checks - n_fail << " native v8.1 stage B controls passed\n";
    return n_fail ? 1 : 0;
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 2;
  }
}
