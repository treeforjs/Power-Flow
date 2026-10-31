// Native v6 births from the sealed heat/source owner v5 through arrival owner
// v1 plans: one exact ExactParticleCount equality per stream (native credit =
// source debit) across a sweep of neutral releases and electron exposures,
// with genuine source v5 kinematics (neutral packets on the emitting face,
// electron layer states at the birth time), native B_s staging, motion of the
// newborns over [tau, end], publication and arrival commit with the native
// receipt (and source commit for the multi-step neutral cases).
// Source frame -> native frame: the source/field frame is centered; native =
// source + (width/2, gap/2, 0), with points on the source patch face
// (y_s = y0) mapped exactly onto the upper foil lower face.
// usage: test_v6_source_arrival OUTPUT_DIRECTORY [reduced]
#include "v6_fixture.hpp"
#include "arrival_v1_births.hpp"
#include "source_v5_harness.hpp"

using namespace v6;
namespace pa = pf_arrival;
namespace av = m::arrival_v1;
namespace sh = source_harness;
int forbidden_calls = 0;
extern "C" void
__wrap__ZN6raikou12mfem_backend6GapPic12solve_fieldsEv(m::GapPic *) {
  ++forbidden_calls;
  throw std::runtime_error("forbidden native solver trap");
}
extern "C" void __wrap__ZN6raikou12mfem_backend6GapPic7advanceEd(m::GapPic *,
                                                                 double) {
  ++forbidden_calls;
  throw std::runtime_error("forbidden native advance trap");
}
constexpr double dt = 5e-15;
bool reduced = false;
std::string hx(double x) {
  std::ostringstream s;
  s << std::hexfloat << x;
  return s.str();
}
m::PicVector native_position(const sh::V3d &x) {
  const auto o = options();
  m::PicVector p{x[0] + .5 * o.width_m, x[1] + .5 * o.gap_m, x[2]};
  if (x[1] == sh::patch_face_y_m)
    p[1] = o.upper_foil.bottom_m;
  return p;
}
std::unique_ptr<m::PicArena> arena_for(const std::string &source_seal) {
  auto g = genesis();
  g.source_sha256 = source_seal;
  return m::PicArena::create_full_time(options(), states(), channels(), base_population(), g);
}
struct Stats {
  long cases = 0, admitted_steps = 0, births = 0, exact_equalities = 0, source_refused = 0,
       arrival_refused = 0, exposure_refused = 0, surface_births = 0, interior_births = 0,
       committed_multi = 0, active_layers = 0;
  J rows = J::array();
};
// One native step that credits the births of `staged` with the source debit.
struct NativeStep {
  m::PicArena &arena;
  StepState &st;
};
m::FullTimeStepBinding open(m::PicArena &a, StepState &st, std::optional<m::PrivateFullTimeStep> &outer,
                            std::optional<m::CollisionScope> &first) {
  const auto k = a.full_time_boundary()->accepted_step;
  auto b = binding_of(a, dt, sha(step_text("v6-common-clock", k)));
  outer.emplace(a.begin_full_time_step(b, 16384, 32));
  m::ExposureBinding e1{b, m::FullTimeSlot::before_transport, b.begin_s, .5 * b.dt_s, st.sequence,
                        b.accepted_field_sha256};
  first.emplace(outer->first_exposure(e1));
  auto t = first->prepare_compound_with_completion(1, [&](m::CollisionScope &) {
    return completion(e1, st.rng, sha(step_text("v6-rng-first", k)));
  });
  if (!first->try_publish(t))
    throw std::runtime_error("first exposure");
  return b;
}
std::shared_ptr<const m::FullTimeStepFlows>
finish(m::PicArena &a, StepState &st, m::PrivateFullTimeStep &outer, const m::FullTimeStepBinding &b,
       const m::PreparedBirths &births, const std::string &source_after) {
  auto trial = outer.fork_transport_trial(births);
  auto mb = m::FullTimeMotionBinding{b, "", b.accepted_field_sha256,
                                     sha(step_text("v6-field-after", b.accepted_step)),
                                     sha(step_text("v6-heat-after", b.accepted_step)), source_after};
  outer.accept_transport_trial(reference_hook(16384)(trial, mb));
  const auto k = b.accepted_step;
  m::ExposureBinding e2{b, m::FullTimeSlot::after_transport, b.end_s, .5 * b.dt_s, st.sequence + 1,
                        sha(step_text("v6-field-after", k))};
  auto s2 = outer.second_exposure(e2);
  auto t2 = s2.prepare_compound_with_completion(1, [&](m::CollisionScope &) {
    return completion(e2, sha(step_text("v6-rng-first", k)), sha(step_text("v6-rng-second", k)));
  });
  if (!s2.try_publish(t2))
    throw std::runtime_error("second exposure");
  auto token = outer.seal({sha(step_text("v6-field-after", k)), sha(step_text("v6-heat-after", k)),
                           source_after, sha(step_text("v6-rng-second", k))});
  auto flows = token.candidate_full_time_flows();
  if (!a.try_publish(token))
    throw std::runtime_error("publication");
  st.sequence += 2;
  st.rng = sha(step_text("v6-rng-second", k));
  return flows;
}
std::vector<av::BirthKinematics> kinematics(const sh::Step &s, Stats &stats) {
  std::vector<av::BirthKinematics> out;
  for (const auto &b : s.births) {
    av::BirthKinematics k;
    k.position_m = native_position(b.x);
    k.proper_velocity_m_s = b.u;
    k.origin = "source-v5/" + b.event;
    (b.x[1] == sh::patch_face_y_m ? stats.surface_births : stats.interior_births) += 1;
    stats.active_layers += b.active_layer;
    out.push_back(k);
  }
  return out;
}
void neutral_case(double N0, double rate, double x_target, int steps, Stats &stats) {
  const double D0 = sh::probe_neutral_debit(N0, rate, dt);
  if (!(D0 > 0))
    return;
  const double w0 = D0 / x_target;
  if (!pa::detail::positive_normal(w0))
    return;
  ++stats.cases;
  pa::Configuration c;
  c.run_seed = 4000 + std::uint64_t(stats.cases);
  c.streams = {{0, pa::Kind::electron, 0, "e", 1, 4096}, {1, pa::Kind::neutral, 0, "n", w0, 4096}};
  auto P = std::make_shared<const pa::Policy>(c);
  sh::Neutral src(N0, rate, P->sampling_policy_identity(), genesis().time_s);
  auto arena = arena_for(src.state_sha256());
  pa::Owner A(P, 0);
  StepState st;
  for (int s = 0; s < steps; ++s) {
    J row = {{"kind", "neutral"}, {"N0", hx(N0)}, {"rate", hx(rate)}, {"x", hx(x_target)}, {"step", s}};
    const auto before = arena->state_prefix();
    std::optional<m::PrivateFullTimeStep> outer;
    std::optional<m::CollisionScope> first;
    const auto b = open(*arena, st, outer, first);
    auto scope = outer->birth_scope();
    check(src.revision() == scope.accepted_step() && A.next_step() == scope.accepted_step(),
          "source revision, arrival step and native accepted step are the same index");
    const auto release = src.quote(*P, 1, dt);
    check(release.begin_s == b.begin_s && release.end_s == b.end_s,
          "source epoch equals the native common interval");
    std::optional<pa::Staged> staged;
    try {
      staged.emplace(A.stage(A.plan({scope.accepted_step(), {}, {release}}), scope.base_native_id()));
    } catch (const pa::Refusal &r) {
      ++stats.arrival_refused;
      row["outcome"] = std::string("arrival-refused/") + pa::kind_name(r.kind);
      stats.rows.push_back(row);
      return;
    }
    auto step = src.prepare(*staged, 1);
    row["leaves"] = staged->births().size();
    if (!step.admitted) {
      ++stats.source_refused;
      row["outcome"] = step.outcome;
      row["refusal"] = step.refusal;
      outer.reset();
      first.reset();
      check(arena->state_prefix() == before, "source refusal: dropped native step left native unchanged");
      stats.rows.push_back(row);
      return;
    }
    // Exact source-side identities.
    const double D = release.debit.positive_carriers().front();
    pa::Count leaves;
    for (const auto &x : staged->births())
      leaves = leaves.credited(x.weight);
    bool bits = step.births.size() == staged->births().size();
    for (std::size_t k = 0; bits && k < step.births.size(); ++k)
      bits = std::bit_cast<std::uint64_t>(step.births[k].weight) ==
             std::bit_cast<std::uint64_t>(staged->births()[k].weight);
    check(step.inventory_after == step.inventory_before.debit(D) && leaves == step.debit && bits,
          "source v5: inventory after = before - D exactly; arrival leaves sum to D; packets carry the leaves");
    const auto kin = kinematics(step, stats);
    auto credit = av::credit(*P, *staged, kin, {{1, step.debit}}, scope.source_before_sha256(),
                             step.source_after);
    auto births = scope.prepare(credit); // native exact equality per stream inside
    ++stats.exact_equalities;
    check(births.receipt().streams.size() == 1 && births.receipt().streams[0].source_debit == step.debit &&
              births.receipt().entries.size() == staged->births().size(),
          "native credit = source v5 debit (exact) for the neutral stream");
    const auto flows = finish(*arena, st, *outer, b, births, step.source_after);
    check(flows && flows->births && flows->births->staged_token_sha256 == staged->token(),
          "native step with source v5 neutral births publishes");
    (void)A.commit(std::move(*staged), av::receipt(*flows->births));
    check(src.commit(), "source v5 commits its trial");
    ++stats.admitted_steps;
    stats.births += long(flows->births->entries.size());
    if (steps > 1 && s + 1 == steps)
      ++stats.committed_multi;
    row["outcome"] = "admitted";
    row["debit"] = hx(D);
    row["native_base"] = flows->births->base;
    stats.rows.push_back(row);
  }
  arena->validate_full_time_history();
  check(A.next_ordinal() == arena->full_time_flow_totals().births,
        "arrival ordinal equals native cumulative births (neutral)");
}
void electron_case(double F, double area, double x_target, Stats &stats) {
  const double N0 = sh::probe_electron_count(F, area, dt);
  J row = {{"kind", "electron"}, {"field_V_m", hx(F)}, {"area_fraction", hx(area)}, {"x", hx(x_target)}};
  if (!(N0 >= 0)) {
    ++stats.exposure_refused;
    row["outcome"] = "exposure-refused";
    stats.rows.push_back(row);
    return;
  }
  const double w0 = N0 / x_target;
  if (!pa::detail::positive_normal(w0))
    return;
  ++stats.cases;
  pa::Configuration c;
  c.run_seed = 9000 + std::uint64_t(stats.cases);
  c.streams = {{0, pa::Kind::electron, 0, "e", w0, 4096}, {1, pa::Kind::neutral, 0, "n", 1, 4096}};
  auto P = std::make_shared<const pa::Policy>(c);
  sh::Electron src(P->sampling_policy_identity(), F, genesis().time_s);
  auto arena = arena_for(src.state_sha256());
  pa::Owner A(P, 0);
  StepState st;
  std::optional<m::PrivateFullTimeStep> outer;
  std::optional<m::CollisionScope> first;
  const auto b = open(*arena, st, outer, first);
  auto scope = outer->birth_scope();
  pa::ElectronExposure input;
  std::string refusal;
  if (!src.expose(0, area, dt, input, refusal)) {
    ++stats.exposure_refused;
    row["outcome"] = "exposure-refused";
    row["refusal"] = refusal;
    stats.rows.push_back(row);
    return;
  }
  check(input.begin_s == b.begin_s && input.end_s == b.end_s, "source exposure inside the native step");
  auto staged = A.stage(A.plan({0, {input}, {}}), scope.base_native_id());
  const auto step = src.births(staged);
  row["births"] = staged.births().size();
  row["expected"] = hx(input.expected_number);
  if (!step.admitted) {
    ++stats.source_refused;
    row["outcome"] = step.outcome;
    row["refusal"] = step.refusal;
    stats.rows.push_back(row);
    return;
  }
  check(step.debit == staged.plan().decisions().front().debit,
        "source v5 candidate weights sum exactly to the arrival debit n * w0");
  bool times = true;
  for (std::size_t k = 0; k < step.births.size(); ++k)
    times = times && step.births[k].time == staged.births()[k].time_s;
  check(times, "source layer time equals the arrival birth time tau for every birth");
  const auto kin = kinematics(step, stats);
  auto births = scope.prepare(av::credit(*P, staged, kin, {{0, step.debit}}, scope.source_before_sha256(),
                                         scope.source_before_sha256()));
  ++stats.exact_equalities;
  check(births.receipt().streams[0].source_debit == step.debit,
        "native credit = source v5 debit (exact) for the electron stream");
  const auto flows = finish(*arena, st, *outer, b, births, scope.source_before_sha256());
  (void)A.commit(std::move(staged), av::receipt(*flows->births));
  ++stats.admitted_steps;
  stats.births += long(flows->births->entries.size());
  row["outcome"] = "admitted";
  stats.rows.push_back(row);
}
int main(int argc, char **argv) {
  try {
    if (argc != 2 && !(argc == 3 && std::string(argv[2]) == "reduced"))
      throw std::runtime_error("usage: test_v6_source_arrival OUTPUT_DIRECTORY [reduced]");
    const std::string dir = argv[1];
    reduced = argc == 3;
    const auto t0 = std::chrono::steady_clock::now();
    Stats stats;
    std::vector<double> N0s = {1e-300, 1e-250, 1, 1e5, 1e10, 1e15, 1e18, 1e20};
    std::vector<double> xs = {0.3, 2.5, 37.5, 1000.25};
    if (reduced) {
      N0s = {1e-300, 1e10, 1e18};
      xs = {0.3, 37.5};
    }
    for (double N0 : N0s)
      for (double rate : {1e14, 1e8})
        for (double x : xs)
          neutral_case(N0, rate, x, 1, stats);
    neutral_case(1e15, 1e14, 12.5, 3, stats);
    neutral_case(1e18, 1e8, 3.7, reduced ? 2 : 3, stats);
    for (double F : {1e8, 1e9, 3e9})
      for (double area : {1., 1e-3, 1e-10})
        for (double x : {0.4, 3.7, 25.5}) {
          if (reduced && x == 3.7)
            continue;
          electron_case(F, area, x, stats);
        }
    check(stats.exact_equalities > 10 && stats.admitted_steps == stats.exact_equalities &&
              stats.surface_births > 0 && stats.births > 100 && stats.committed_multi >= 2,
          "sweep: every admitted step closed with exact per-stream equalities; multi-step stocks carried");
    check(forbidden_calls == 0, "no native solver/advance call");
    J out = {{"status", "passed"},
             {"check_count", checks},
             {"checks", results},
             {"cases", stats.cases},
             {"admitted_steps", stats.admitted_steps},
             {"births", stats.births},
             {"exact_stream_equalities", stats.exact_equalities},
             {"source_refused", stats.source_refused},
             {"arrival_refused", stats.arrival_refused},
             {"exposure_refused", stats.exposure_refused},
             {"surface_births", stats.surface_births},
             {"interior_births", stats.interior_births},
             {"active_layer_births", stats.active_layers},
             {"multi_step_committed_cases", stats.committed_multi},
             {"rows", stats.rows},
             {"reduced", reduced},
             {"seconds", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()}};
    std::ofstream(dir + "/v6-source-arrival-tests.json") << out.dump(2) << '\n';
    std::cout << checks << " native v6 source/arrival controls passed (" << stats.exact_equalities
              << " exact stream equalities, " << stats.births << " births)\n";
    return 0;
  } catch (const std::exception &x) {
    std::cerr << "FAILED after " << checks << " checks: " << x.what() << '\n';
    if (argc >= 2)
      std::ofstream(std::string(argv[1]) + "/v6-source-arrival-tests.json")
          << J({{"status", "failed"}, {"failure", x.what()}, {"check_count", checks}, {"checks", results}}).dump(2)
          << '\n';
    return 1;
  }
}
