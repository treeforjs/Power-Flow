// Native v8 device admission of steps with births and first-hit sinks (API-008.txt 4.1):
//   CPU/GPU identity   host, device_verified and device runs of the same arrival-birth
//                      sequence (births every step: emitting-face, sink-prone and interior
//                      placements; newborns that sink in their first flight; co-located
//                      multi-carrier packets sunk whole) end every step in the same accepted
//                      identity, flows digest and sink ledger
//   device residency   in device mode every admission runs on the device lane (no host
//                      admission), births are appended to the device rows
//   refusal parity     resealed records with a shifted sink fraction, a changed newborn flight
//                      and a changed incident velocity refuse with the same typed refusal and
//                      text in device mode as on the host (the host path supplies the text)
// usage: test_v8_flows_gpu OUT_DIR NVRTC UUID [STEPS]
#include "v6_fixture.hpp"
#include "arrival_v1_births.hpp"
#include "forms_gpu.hpp"
#include "hall_gpu_common.hpp"
#include "raikou/mfem_pic_full_time_device.h"
#include <filesystem>
#include <fstream>
#include <iostream>

using namespace v6;
namespace pa = pf_arrival;
namespace av = m::arrival_v1;
namespace {
int n_checks = 0, n_fail = 0;
void need(bool ok, const std::string &what, const J &d = J::object()) {
  ++n_checks;
  if (!ok) {
    ++n_fail;
    std::cout << "FAIL " << what << " " << d.dump() << "\n";
  }
}
std::shared_ptr<const pa::Policy> arrival_policy() {
  pa::Configuration c;
  c.run_seed = 20261031;
  c.streams = {{0, pa::Kind::electron, 0, "e", 4, 4096}, {1, pa::Kind::neutral, 0, "n", 2, 4096}};
  return std::make_shared<const pa::Policy>(c);
}
pa::StepInput arrival_input(std::uint64_t step, double begin, double dt, double electrons, double neutral_debit) {
  pa::StepInput in;
  in.step = step;
  in.electrons = {{0, electrons, std::log(electrons), true, begin, dt, begin + dt}};
  in.neutrals = {{1, pa::Count(neutral_debit), std::log(neutral_debit), pa::Count(1e20), begin, begin + dt}};
  return in;
}
// Placement 0: emitting-face births (upper foil lower face, into the gap); 1: sink-prone births
// (electrons at the outer x face moving +x, sinking in their first flight; neutrals just above the
// lower foil moving -y); 2: interior births.
std::vector<av::BirthKinematics> placement(const pa::Staged &s, int kind) {
  std::vector<av::BirthKinematics> out;
  for (const auto &b : s.births()) {
    av::BirthKinematics k;
    const double dx = 1e-6 * double(b.index % 1000), dz = 1e-6 * double(b.index / 1000);
    const bool e = b.kind == pa::Kind::electron;
    if (kind == 0) {
      k.position_m = {.0030 + dx, .007, -.001 + dz + (e ? 0 : .002)};
      k.proper_velocity_m_s = {0, e ? -2e6 : -650, 0};
    } else if (kind == 1) {
      if (e) {
        k.position_m = {.009 - 1e-7, .0045 + dx, dz};
        k.proper_velocity_m_s = {1e8, 0, 0};
      } else {
        k.position_m = {.004 + dx, .002 + 1e-12, dz};
        k.proper_velocity_m_s = {0, -650, 0};
      }
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
std::vector<m::FullTimeAuthoredPacket> population() {
  auto p = base_population();
  // co-located two-carrier packets: one leaving through the outer x face (sunk whole), one staying
  p.push_back({"pair-out", "e", "manufactured/source", {1, 1e-20}, {.009 - 2e-7, .0050, .002}, {1e8, 0, 0}});
  p.push_back({"pair-in", "n", "manufactured/source", {2, 1e-19}, {.0050, .0050, -.002}, {0, 0, 300}});
  for (int k = 0; k < 24; ++k)
    p.push_back({"bulk-" + std::to_string(k), k % 2 ? "e" : "i", "manufactured/source", {1.5 + k, 0},
                 {.0025 + 1e-4 * k, .0045 + 2e-5 * k, -.003 + 2e-4 * k},
                 {k % 3 ? 2e5 : -1e5, 1e4 * (k % 5), 3e4}});
  return p;
}
struct Run {
  std::vector<std::string> trace;
  std::uint64_t births = 0, sunk = 0, sunk_newborn_steps = 0;
  m::FullTimeDeviceCounters before, after;
  std::string final_prefix;
};
struct Harness {
  std::shared_ptr<const pa::Policy> P = arrival_policy();
  std::unique_ptr<pa::Owner> owner = std::make_unique<pa::Owner>(P, 0);
  std::optional<pa::Staged> staged;
  std::uint64_t source_index = 0;
  int place = 0;
  double begin = 0, dt = 0;
  Hooks hooks(Hooks::Motion motion) {
    Hooks h;
    h.births = [this](m::BirthScope &scope) -> std::optional<m::FullTimeBirthCredit> {
      auto plan = owner->plan(arrival_input(scope.accepted_step(), begin, dt, 6, 1.5));
      staged.emplace(owner->stage(std::move(plan), scope.base_native_id()));
      bind_neutral_times(*staged);
      return av::credit(*P, *staged, placement(*staged, place), av::plan_debits(staged->plan()),
                        scope.source_before_sha256(), sha("source/" + std::to_string(++source_index)));
    };
    h.motion = std::move(motion);
    return h;
  }
  StepResult step(m::PicArena &a, StepState &st, Hooks::Motion motion = reference_hook(512)) {
    StepPlan plan;
    plan.dt_s = 5e-15;
    plan.maximum_carriers = 512;
    begin = a.full_time_boundary()->time_s;
    dt = plan.dt_s;
    auto r = full_step(a, plan, st, hooks(std::move(motion)));
    (void)owner->commit(std::move(*staged), av::receipt(*r.flows->births));
    staged.reset();
    return r;
  }
};
std::unique_ptr<m::PicArena> fresh_v8() {
  m::FullTimeConfiguration c;
  c.identity = m::FullTimeIdentitySchema::v8;
  c.limits.maximum_carriers = 512;
  c.limits.memory_budget_bytes = 512 * m::full_time_bytes_per_carrier + m::full_time_fixed_bytes;
  c.history.anchor_interval_steps = 16;
  return m::PicArena::create_full_time(options(), states(), channels(), population(), genesis(), c);
}
std::string sink_digest(const std::shared_ptr<const m::FullTimeStepFlows> &f) {
  if (!f || !f->sinks)
    return "-";
  return f->sinks->ledger_sha256 + "/" + f->sinks->records_root_sha256 + "/" + std::to_string(f->sinks->exposed.size());
}
} // namespace

int main(int argc, char **argv) {
  try {
    if (argc < 4)
      throw std::runtime_error("usage: test_v8_flows_gpu OUT_DIR NVRTC UUID [STEPS]");
    const std::filesystem::path out = argv[1];
    std::filesystem::create_directories(out);
    const int steps = argc > 4 ? std::atoi(argv[4]) : 48;
    hall_gpu::SessionOptions so;
    so.require_uuid = argv[3];
    so.nvrtc_library = argv[2];
    auto session = std::make_unique<hall_gpu::Session>(so);
    auto forms = std::make_unique<pf_v7::gpu::Forms>(*session);
    auto dev = std::make_shared<m::FullTimeDevice>(*session, *forms);
    J report = J::object();
    auto run = [&](m::FullTimeExecution mode) {
      Run R;
      auto a = fresh_v8();
      if (mode != m::FullTimeExecution::host)
        a->attach_full_time_device(dev, mode);
      Harness H;
      StepState st;
      R.before = dev->counters();
      for (int k = 0; k < steps; ++k) {
        H.place = k % 3;
        const auto r = H.step(*a, st);
        R.births += r.flows->births->entries.size();
        if (r.flows->sinks) {
          R.sunk += r.flows->sinks->carriers;
          bool newborn_sunk = false;
          for (const auto &x : r.flows->sinks->records)
            newborn_sunk = newborn_sunk || x.native_id >= r.birth_base;
          R.sunk_newborn_steps += newborn_sunk;
        }
        const auto b = a->full_time_boundary();
        R.trace.push_back(b->accepted_arena_sha256 + "/" + a->full_time_flow_totals().flows_sha256 + "/" +
                          sink_digest(r.flows));
      }
      a->validate_full_time_history();
      R.after = dev->counters();
      R.final_prefix = a->state_prefix();
      return R;
    };
    const auto H = run(m::FullTimeExecution::host);
    const auto V = run(m::FullTimeExecution::device_verified);
    const auto D = run(m::FullTimeExecution::device);
    std::size_t first_diff_v = steps, first_diff_d = steps;
    for (int k = 0; k < steps; ++k) {
      if (first_diff_v == std::size_t(steps) && H.trace[k] != V.trace[k])
        first_diff_v = k;
      if (first_diff_d == std::size_t(steps) && H.trace[k] != D.trace[k])
        first_diff_d = k;
    }
    need(first_diff_v == std::size_t(steps) && H.final_prefix == V.final_prefix,
         "device_verified run equals the host run every step (accepted identity, flows, sink ledger)",
         {{"first_difference", first_diff_v}});
    need(first_diff_d == std::size_t(steps) && H.final_prefix == D.final_prefix,
         "device run equals the host run every step (accepted identity, flows, sink ledger)",
         {{"first_difference", first_diff_d}});
    need(2 * H.births > std::uint64_t(steps) && H.sunk >= std::uint64_t(steps) / 2 && H.sunk_newborn_steps > 0,
         "the sequence has births every step, sinks and newborns sunk in their first flight",
         {{"births", H.births}, {"sunk", H.sunk}, {"newborn_sink_steps", H.sunk_newborn_steps}});
    const auto dv = [](const Run &R) {
      return J{{"device_admissions", R.after.device_admissions - R.before.device_admissions},
               {"host_admissions", R.after.host_admissions - R.before.host_admissions},
               {"verified_comparisons", R.after.verified_comparisons - R.before.verified_comparisons},
               {"host_refusal_texts", R.after.host_refusal_texts - R.before.host_refusal_texts},
               {"row_uploads", R.after.row_uploads - R.before.row_uploads}};
    };
    const auto cv = dv(V), cd = dv(D);
    need(cv["verified_comparisons"] == std::uint64_t(steps) && cv["host_admissions"] == 0,
         "device_verified: every admission verified against the host path", cv);
    need(cd["device_admissions"] == std::uint64_t(steps) && cd["host_admissions"] == 0 && cd["host_refusal_texts"] == 0,
         "device: every admission (births and sinks included) on the device lane", cd);
    // Refusal parity (device mode vs host) on the state after two steps.
    struct Mut {
      const char *name;
      std::function<void(ReferenceMotion &)> f;
    };
    const std::vector<Mut> muts = {
        {"sink fraction shifted", [](ReferenceMotion &x) {
           if (!x.sinks.empty())
             x.sinks[0].fraction = std::nextafter(x.sinks[0].fraction, 2.0);
         }},
        {"newborn flight changed", [](ReferenceMotion &x) {
           if (!x.births.empty())
             x.births[0].flight_s *= 0.5;
         }},
        {"incident velocity changed", [](ReferenceMotion &x) {
           if (!x.sinks.empty())
             x.sinks[0].incident_proper_velocity_m_s[0] *= 1.5;
         }}};
    J parity = J::array();
    for (const auto &mu : muts) {
      std::string text[2];
      int kind[2] = {-1, -1};
      bool had[2] = {false, false};
      for (int side = 0; side < 2; ++side) {
        auto a = fresh_v8();
        if (side)
          a->attach_full_time_device(dev, m::FullTimeExecution::device);
        Harness h;
        StepState st;
        h.place = 1;
        (void)h.step(*a, st);
        (void)h.step(*a, st);
        bool applicable = false;
        auto probe = [&](ReferenceMotion &x) {
          applicable = !x.sinks.empty() && !x.births.empty();
          mu.f(x);
        };
        const auto before = a->state_prefix();
        try {
          (void)h.step(*a, st, reference_hook(512, probe, true));
        } catch (const m::FullTimeFlowRefusal &e) {
          text[side] = e.what();
          kind[side] = int(e.kind());
        } catch (const std::exception &e) {
          text[side] = std::string("untyped: ") + e.what();
        }
        had[side] = applicable && a->state_prefix() == before;
      }
      parity.push_back({{"mutation", mu.name}, {"host", text[0]}, {"device", text[1]}, {"kind", kind[0]}});
      need(had[0] && had[1] && kind[0] >= 0 && kind[0] == kind[1] && text[0] == text[1],
           std::string("refusal parity: ") + mu.name, parity.back());
    }
    report["steps"] = steps;
    report["host"] = {{"births", H.births}, {"sunk", H.sunk}, {"newborn_sink_steps", H.sunk_newborn_steps}};
    report["device_verified_counters"] = cv;
    report["device_counters"] = cd;
    report["refusal_parity"] = parity;
    report["final_prefix"] = H.final_prefix;
    report["checks"] = n_checks;
    report["failures"] = n_fail;
    std::ofstream(out / "v8-flows-gpu.json") << report.dump(2) << "\n";
    std::cout << (n_fail ? "FAILED " : "") << n_checks - n_fail << " native v8 device flow controls passed\n";
    dev.reset();
    forms.reset();
    (void)session->close();
    return n_fail ? 1 : 0;
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 2;
  }
}
