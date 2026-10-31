#include "parent_coverage_v2_subflow.hpp"
#include <fstream>
#include <iostream>

// Coverage-scope checks for the water parent-ionization channel. Conditional
// parent-only operator at a frozen particle stage; no physical trajectory time
// advances and no Hall coupling is constructed.
namespace nw = raikou::native_water;
namespace cv = raikou::native_water::coverage_v2;
namespace m = raikou::mfem_backend;
namespace pc = raikou::pic_collisions;
namespace w = raikou::pic_water;
using J = nw::Json;
J checks = J::array();
void check(bool ok, const std::string& name) {
  checks.push_back({{"name", name}, {"passed", ok}});
  if (!ok) throw std::runtime_error("check failed: " + name);
}
template<class F> void reject(F f, const std::string& s) {
  bool threw = false;try { f(); } catch (const std::exception&) { threw = true; }check(threw, s);
}
template<class F> cv::CoverageRefusal refusal(F f, const std::string& s) {
  try { f(); } catch (const cv::CoverageRefusal& r) { check(true, s);return r; } catch (const std::exception& e) {
    check(false, s + " (other error: " + e.what() + ")");
  }
  check(false, s + " (no refusal)");throw std::logic_error("unreachable");
}
constexpr double mn = 2.9915e-26;
const double inf = std::numeric_limits<double>::infinity();
m::GapPicOptions options() {
  m::GapPicOptions o;o.nx = o.ny = 45;o.width_m = o.gap_m = .009;o.depth_m = .025;
  o.finite_foils = true;o.finite_axial_extent = true;o.short_z_m = 0;o.feed_z_m = .025;
  o.lower_foil = {.002, .007, .0018, .002};o.upper_foil = {.002, .007, .007, .0072};return o;
}
m::PicVector electron_u(double E) { double x = E * pc::qe / (w::me * pc::c * pc::c);return {pc::c * std::sqrt(x * (x + 2)), 0, 0}; }
std::vector<m::ArenaStateSpec> states() {
  return {{"electron", 0, 0, -1, w::me, 0, 0}, {"H2O", 2, 1, 0, mn, 0, 0},
          {"H2O+", 2, 1, 1, mn - w::me + w::threshold_eV * pc::qe / (pc::c * pc::c), 0, w::threshold_eV * pc::qe}};
}
std::vector<m::ArenaChannel> channels(const nw::Binding& b) {
  return {{w::parent_channel, "electron", "H2O", {"electron", "electron", "H2O+"}, b.dataset_sha256, b.kernel_sha256,
           b.kernel_binding_file_sha256, w::parent_model}};
}
const m::PicVector here{.004, .004, .0125}, water_at{.00403, .004, .0125}, elsewhere{.005, .004, .0125};
struct Electron { double weight, energy_eV;m::PicVector position; };
// Carriers are added electrons-first, then water, then extra electrons, so
// native IDs agree between a reference and a variant that only moves extras.
struct Fixture {
  m::GapPic pic{options()};std::unique_ptr<m::PicArena> arena;
  Fixture(const nw::Binding& b, std::vector<Electron> lead, double water_weight, std::vector<Electron> extra = {}) {
    std::vector<m::GapParticle> batch;
    for (const auto& x : lead) { m::GapParticle e;e.weight = x.weight;e.mass_kg = w::me;e.position_m = x.position;e.proper_velocity_m_s = electron_u(x.energy_eV);batch.push_back(e); }
    m::GapParticle n;n.weight = water_weight;n.mass_kg = mn;n.charge_C = 0;n.position_m = water_at;batch.push_back(n);
    for (const auto& x : extra) { m::GapParticle e;e.weight = x.weight;e.mass_kg = w::me;e.position_m = x.position;e.proper_velocity_m_s = electron_u(x.energy_eV);batch.push_back(e); }
    auto ids = pic.add_authored_batch(batch);pic.solve_fields();
    std::vector<m::ArenaPacket> packets;std::size_t k = 0;
    for (std::size_t i = 0; i < lead.size(); ++i, ++k)
      packets.push_back({i ? "e" + std::to_string(i) : "e", "electron", "authored_fixture", "", {lead[i].weight, 0}, {ids[k]}});
    packets.push_back({"n", "H2O", "authored_fixture", "", {water_weight, 0}, {ids[k++]}});
    for (std::size_t i = 0; i < extra.size(); ++i, ++k)
      packets.push_back({"x" + std::to_string(i), "electron", "authored_fixture", "", {extra[i].weight, 0}, {ids[k]}});
    arena = std::make_unique<m::PicArena>(pic, states(), channels(b), packets, "native-water-subflow-fixture/v1");
  }
};
template<class Flow>
nw::Interval interval(m::PicArena& arena, const Flow& flow, double duration, std::uint64_t sequence = 0) {
  nw::Interval i;i.sequence = sequence;i.end_s = duration;i.common_clock_sha256 = nw::sha("frozen-native-parent-fixture/" + std::to_string(sequence));
  i.phase.state_prefix = arena.state_prefix();i.phase.mesh_prefix = flow.grid_identity();i.phase.clock_prefix = i.common_clock_sha256;return i;
}
long double weight(const m::PicArena& a, const std::string& state) {
  long double total = 0;
  for (const auto& [_, p] : a.snapshot()->packets) if (p.state == state) for (double part : m::PicArena::positive_carriers(p.count)) total += part;
  return total;
}
J events_without_native(const J& audit) {
  J out = J::array();for (auto e : audit.at("events")) { e.erase("native");out.push_back(e); }return out;
}
J rng_fields(const J& cp) { return {cp.at("scheduler"), cp.at("product_rng"), cp.at("next_event"), cp.at("next_interval")}; }
// Every carrier except the named moved packets, with exact identity/kinematics.
J census(const m::PicArena& a, const std::set<std::string>& moved = {}) {
  J out = J::array();
  for (const auto& r : a.rate_population(100000)) {
    J row{{"packet", r.packet_id}, {"state", r.state}, {"id", r.native.id}, {"weight", r.native.weight},
          {"mass", r.native.mass_kg}, {"u", r.native.proper_velocity_m_s}};
    if (!moved.count(r.packet_id)) row["x"] = r.native.position_m;
    out.push_back(row);
  }
  return out;
}
double value(const J& exact) { return exact.at("value").get<double>(); }
cv::CoveragePolicy record_only() { cv::CoveragePolicy p;p.max_exposure_unresolved_events = inf;p.max_cumulative_unresolved_events = inf;p.max_unresolved_reaction_fraction = 1;p.admit_unbounded = true;return p; }

int main(int argc, char** argv) {
  try {
    if (argc != 4) throw std::runtime_error("dataset binding output-directory");
    std::ifstream f(argv[2]);J binding;f >> binding;
    nw::Binding b{binding.at("dataset_sha256"), binding.at("kernel_sha256"), binding.at("kernel_binding_file_sha256")};
    const std::string dir = argv[3];
    auto data = std::make_shared<const w::Data>(w::Data::read(argv[1]));
    w::ParentBridge bridge(data, w::parent_model, b.dataset_sha256, b.kernel_sha256, mn);
    nw::Limits limits;nw::Grid grid(options(), {45, 45, 5}, limits.cells);
    const double V = grid.cells[0].volume_m3;
    J receipts = J::array(), coverage_map = J::array();

    // ---------------- exact accumulator ----------------
    {
      cv::ExactSum a, b2;
      for (int i = 0; i < 10; ++i) a.add(.1);
      check(a.value() == 1.0 && a.mantissa_hex() != "1", "ten binary64 0.1 values: exact sum kept, correctly rounded to 1");
      std::vector<double> xs{1e300, 3.5, 1e-300, 2.5e-310, 7.0, 1e16, 1.0, 0x1p-1074, 123456.789};
      for (double x : xs) b2.add(x);
      cv::ExactSum c2;for (auto it = xs.rbegin(); it != xs.rend(); ++it) c2.add(*it);
      check(b2 == c2 && b2.mantissa_hex() == c2.mantissa_hex(), "exact sum independent of addition order");
      cv::ExactSum big;big.add(0x1p53);big.add(1.0);
      check(big.value() == 0x1p53 && big.mantissa_hex() == "20000000000001" && big.exponent2() == 0, "2^53+1 held exactly, rounds half to even");
      big.add(1.0);check(big.value() == 0x1p53 + 2, "exact carry into representable sum");
      cv::ExactSum sub;sub.add(0x1p-1074);sub.add(0x1p-1074);
      check(sub.value() == 0x1p-1073 && sub.mantissa_hex() == "1" && sub.exponent2() == -1073, "subnormal addends exact");
      cv::ExactSum huge;huge.add(std::numeric_limits<double>::max());huge.add(std::numeric_limits<double>::max());
      check(std::isinf(huge.value()) && cv::ExactSum::from_canonical(huge.mantissa_hex(), huge.exponent2()) == huge,
            "beyond-binary64 exact sum retained with infinite rounded report");
      check(cv::ExactSum::from_canonical(b2.mantissa_hex(), b2.exponent2()) == b2 && cv::ExactSum::from_canonical("0", 0).zero(),
            "canonical exact round trip");
      reject([&] { (void)cv::ExactSum::from_canonical("2", 0); }, "even canonical mantissa refused");
      reject([&] { (void)cv::ExactSum::from_canonical("0f", 0); }, "leading-zero canonical mantissa refused");
      reject([&] { (void)cv::ExactSum::from_canonical("F", 0); }, "uppercase canonical mantissa refused");
      reject([&] { cv::ExactSum s;s.add(-1.0); }, "negative addend refused");
      reject([&] { cv::ExactSum s;s.add(std::nan("")); }, "NaN addend refused");
      cv::ExactSum tie;tie.add(0x1p53);tie.add(3.0);check(tie.value() == 0x1p53 + 4, "round half to even upward");
    }

    // ---------------- scope classification and coverage map ----------------
    const cv::CoverageScope scope(data, b.dataset_sha256, bridge);
    const cv::CoverageScope open(data, b.dataset_sha256, bridge, record_only());
    check(scope.identity() != open.identity() && pc::digest(scope.identity()), "policy is bound into the coverage scope identity");
    {
      const double s_inel = 2.800285205390781e-21;
      struct Row { double E;cv::Kind kind;cv::Sigma sigma;double upper; };
      const std::vector<Row> rows{
        {0, cv::Kind::known_zero, cv::Sigma::none, 0}, {12.621, cv::Kind::known_zero, cv::Sigma::none, 0},
        {std::nextafter(12.621, 13.), cv::Kind::onset_gap, cv::Sigma::bounded, 1.6660000000000001 * s_inel},
        {13., cv::Kind::onset_gap, cv::Sigma::bounded, 1.6660000000000001 * s_inel},
        {std::nextafter(13.5, 0.), cv::Kind::onset_gap, cv::Sigma::bounded, 1.6660000000000001 * s_inel},
        {13.5, cv::Kind::supported, cv::Sigma::tabulated, 0}, {38.3, cv::Kind::supported, cv::Sigma::tabulated, 0},
        {1000., cv::Kind::supported, cv::Sigma::tabulated, 0},
        {std::nextafter(1000., 2000.), cv::Kind::above_table, cv::Sigma::bounded, 2.9828200000000002 * s_inel},
        {8000., cv::Kind::above_table, cv::Sigma::bounded, 0.62408 * s_inel},
        {10000., cv::Kind::above_table, cv::Sigma::bounded, 0.49404000000000003 * s_inel},
        {std::nextafter(10000., 2e4), cv::Kind::above_table, cv::Sigma::unbounded, 0},
        {1e7, cv::Kind::above_table, cv::Sigma::unbounded, 0}};
      for (const auto& r : rows) {
        const auto c = scope.classify(r.E);
        const bool sigma_ok = r.kind == cv::Kind::supported ? c.sigma_m2 == data->observable_value_SI(w::parent_channel, r.E) : true;
        check(c.kind == r.kind && c.sigma == r.sigma && c.upper_m2 == r.upper && sigma_ok &&
              scope.excluded(r.E) == cv::unresolved(r.kind), std::string("coverage class at ") + std::to_string(r.E) + " eV");
        bool coverage_error = false;
        try { (void)bridge.channel().sigma_m2(r.E); } catch (const raikou::kinetics::CollisionCoverageError&) { coverage_error = true; }
        check(coverage_error == cv::unresolved(r.kind), std::string("sealed bridge refuses exactly the unresolved energy ") + std::to_string(r.E));
      }
      // Product energies of the equal-pair kernel land in the onset gap only
      // for incident energies between 2*12.621+12.621 and 2*13.5+12.621 eV.
      check(scope.classify((38.3 - w::threshold_eV) / 2).kind == cv::Kind::onset_gap &&
            scope.classify((37.8 - w::threshold_eV) / 2).kind == cv::Kind::known_zero &&
            scope.classify((39.7 - w::threshold_eV) / 2).kind == cv::Kind::supported, "equal-pair products of 37.86-39.62 eV incidents land in the onset gap");
      const auto& d = scope.description();
      check(d.at("table_support_eV") == J::array({13.5, 1000.}) && d.at("kernel_support_eV") == J::array({13.5, 1000.}) &&
            d.at("known_zero_at_or_below_eV") == 12.621 && d.at("upper_bound").at("support_eV") == J::array({7.5, 10000.}),
            "scope description records parent table, kernel and bound supports");
      coverage_map = d.at("bins");
      reject([&] { cv::CoveragePolicy p;p.max_unresolved_reaction_fraction = 1.5;cv::CoverageScope bad(data, b.dataset_sha256, bridge, p); }, "fraction limit above one refused");
      reject([&] { cv::CoveragePolicy p;p.max_exposure_unresolved_events = std::nan("");cv::CoverageScope bad(data, b.dataset_sha256, bridge, p); }, "NaN event limit refused");
      reject([&] { cv::CoverageScope bad(data, nw::sha("other"), bridge); }, "scope bound to another dataset digest refused");
      reject([&] { cv::CoverageScope bad(data, b.dataset_sha256, bridge, {}, "H2O:munoz2007:recommended_total"); }, "non-inelastic bound observable refused");
    }

    // ---------------- pure limit and ledger transitions ----------------
    {
      cv::Measures u;cv::Applied a;
      auto d = cv::evaluate(scope.policy(), u, a, .01, "exposure");
      check(!d.refused && d.fraction_upper == 0, "empty exposure passes the limit");
      cv::Measures::Terms t;t.we = 1;t.wn = 1;t.phys = 1;t.ev = 1;t.up_phys = .05;t.up_ev = .05;t.sigma = cv::Sigma::bounded;u.add(t);
      d = cv::evaluate(scope.policy(), u, a, .01, "exposure");
      check(d.refused && d.absolute_exceeded && d.fraction_exceeded && d.criterion == "events_and_fraction_upper_bounds", "both upper bounds exceeded refuses");
      a.physical.add(100);d = cv::evaluate(scope.policy(), u, a, .01, "exposure");
      check(!d.refused && d.absolute_exceeded && !d.fraction_exceeded, "small unresolved fraction admits despite event bound");
      cv::Measures small;t.up_phys = t.up_ev = .005;small.add(t);d = cv::evaluate(scope.policy(), small, cv::Applied{}, .01, "exposure");
      check(!d.refused && !d.absolute_exceeded && d.fraction_exceeded && d.fraction_upper == 1, "small event bound admits despite unit fraction");
      cv::Measures ub;t.sigma = cv::Sigma::unbounded;t.up_phys = t.up_ev = 0;ub.add(t);
      d = cv::evaluate(scope.policy(), ub, a, .01, "exposure");
      check(d.refused && d.criterion == "unbounded_unresolved_sigma", "unbounded unresolved sigma refused by default");
      check(!cv::evaluate(open.policy(), ub, a, inf, "exposure").refused, "record-only policy admits unbounded exposure");
      // Full-time two-exposure transition, pure.
      cv::ExposureLedger l1;l1.slot = cv::Slot::before_transport;l1.sequence = 4;l1.unresolved = small;
      l1.bins[cv::BinKey{3, cv::Kind::onset_gap, cv::Sigma::bounded}] = small;
      cv::Ledger base;
      auto staged = scope.advance(base, l1);
      check(staged.pending && staged.last_step.empty() && staged.cumulative.accepted_exposures == 0, "first exposure ledger is staged, not accepted");
      reject([&] { (void)scope.advance(staged, l1); }, "second first-exposure while one is pending refused");
      cv::ExposureLedger l2;l2.slot = cv::Slot::after_transport;l2.sequence = 5;
      auto closed = scope.advance(staged, l2);
      check(!closed.pending && closed.last_step.size() == 2 && closed.cumulative.accepted_steps == 1 &&
            closed.cumulative.accepted_exposures == 2 && closed.cumulative.exposures_with_exclusions == 1 &&
            closed.cumulative.steps_with_exclusions == 1 && closed.cumulative.unresolved.json() == small.json(),
            "second exposure closes the step and accumulates both exposures");
      cv::ExposureLedger orphan = l2;reject([&] { (void)scope.advance(base, orphan); }, "second exposure without staged first refused");
      cv::ExposureLedger wrong = l2;wrong.sequence = 9;reject([&] { (void)scope.advance(staged, wrong); }, "second exposure with foreign sequence refused");
      const auto again = scope.ledger_from(scope.ledger_json(closed));
      check(scope.ledger_json(again).dump() == scope.ledger_json(closed).dump(), "ledger JSON round trip is canonical");
      // Cumulative refusal at the second exposure.
      cv::CoveragePolicy tight;tight.max_cumulative_unresolved_events = .008;
      const cv::CoverageScope cum(data, b.dataset_sha256, bridge, tight);
      auto staged2 = cum.advance(base, l1);
      cv::ExposureLedger l2b = l2;l2b.unresolved = small;
      const auto r = refusal([&] { (void)cum.advance(staged2, l2b); }, "cumulative limit refuses the closing second exposure");
      check(r.criterion == "events_and_fraction_upper_bounds" && r.exposure_ledger.at("cumulative_limit").at("refused") == true &&
            r.exposure_ledger.at("exposure_limit").at("refused") == false, "refusal names the cumulative criterion");
    }

    // ---------------- each unsupported range, single pair ----------------
    struct Range { const char* label;double E;cv::Kind kind;cv::Sigma sigma; };
    for (const Range& r : {Range{"13 eV onset gap", 13., cv::Kind::onset_gap, cv::Sigma::bounded},
                           Range{"1000.5 eV above parent table", 1000.5, cv::Kind::above_table, cv::Sigma::bounded},
                           Range{"8000 eV above parent table", 8000., cv::Kind::above_table, cv::Sigma::bounded},
                           Range{"20000 eV above every bounding table", 20000., cv::Kind::above_table, cv::Sigma::unbounded}}) {
      const std::string L = r.label;const double duration = 1e3;
      {
        Fixture old(b, {{1, r.E, here}}, 16);nw::ParentSubflow flow(*old.arena, grid.counts, b, limits, bridge, 1, 2);
        const auto cp = flow.checkpoint();
        reject([&] { flow.prepare(interval(*old.arena, flow, duration)); }, L + ": sealed v3 subflow aborts the exposure");
        check(flow.checkpoint() == cp, L + ": v3 abort leaves v3 owners unchanged");
      }
      {
        Fixture c(b, {{1, r.E, here}}, 16);cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, scope, 1, 2);
        const auto cp = flow.checkpoint();const auto prefix = c.arena->state_prefix();const auto native = m::PicArena::fingerprint(c.pic);
        const auto x = refusal([&] { flow.prepare(interval(*c.arena, flow, duration)); }, L + ": default limit refuses with typed scoped refusal");
        check(x.criterion == (r.sigma == cv::Sigma::unbounded ? "unbounded_unresolved_sigma" : "events_and_fraction_upper_bounds") &&
              x.scope_identity == scope.identity() && x.exposure_ledger.at("unresolved").at("pair_evaluations") == 1,
              L + ": refusal carries criterion, scope and the unresolved ledger");
        check(flow.checkpoint() == cp && c.arena->state_prefix() == prefix && m::PicArena::fingerprint(c.pic) == native,
              L + ": refusal leaves native, RNG and ledger unchanged");
      }
      Fixture c(b, {{1, r.E, here}}, 16);cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, open, 1, 2);
      auto trial = flow.prepare(interval(*c.arena, flow, duration));
      const auto& cov = trial.audit().at("coverage");const auto& u = cov.at("unresolved");
      check(trial.audit().at("events").empty() && trial.audit().at("pairs") == 1, L + ": excluded pair yields no event and the step proceeds");
      check(u.at("pair_evaluations") == 1 && cov.at("distinct_excluded_pairs") == 1 && cov.at("distinct_excluded_electron_carriers") == 1 &&
            cov.at("distinct_excluded_water_carriers") == 1 && value(cov.at("distinct_electron_weight_sum_at_first_exclusion")) == 1 &&
            value(cov.at("distinct_water_weight_sum_at_first_exclusion")) == 16 && value(u.at("electron_weight_sum")) == 1 &&
            value(u.at("water_weight_sum")) == 16 && value(u.at("event_weight_sum")) == 1, L + ": exact count and weight sums");
      pc::Packet e{"e", "electron", 0, 1, w::me, electron_u(r.E)}, n{"n", "H2O", 0, 16, mn, {}};
      const auto kin = pc::relative(e, n);const auto cls = scope.classify(kin.target_frame_energy_eV);
      const long double dt = duration;
      const long double phys = 1.L * 16 * kin.moller_flux_m_s / V * dt, ev = 16.L * kin.moller_flux_m_s / V * dt;
      check(value(u.at("physical_exposure_inv_m2")) == static_cast<double>(phys) && value(u.at("event_exposure_inv_m2")) == static_cast<double>(ev) &&
            std::abs(value(u.at("physical_exposure_inv_m2")) - 16 * kin.moller_flux_m_s * duration / V) <= 1e-14 * value(u.at("physical_exposure_inv_m2")),
            L + ": integrated sigma-free exposure equals the independent pair integral");
      if (r.sigma == cv::Sigma::bounded)
        check(u.at("upper_bounded_sigma_pair_evaluations") == 1 && value(u.at("upper_bound_physical_reactions")) == static_cast<double>(phys * cls.upper_m2) &&
              value(u.at("tabulated_expected_physical_reactions")) == 0 && u.at("lower_bound_physical_reactions") == 0,
              L + ": sigma unknown, recorded with zero lower bound and inelastic upper bound");
      else
        check(u.at("unbounded_sigma_pair_evaluations") == 1 && value(u.at("upper_bound_physical_reactions")) == 0 &&
              value(u.at("unbounded_physical_exposure_inv_m2")) == static_cast<double>(phys) && cov.at("exposure_limit").at("unbounded_admitted") == true,
              L + ": no data upper bound; exposure recorded as unbounded and admitted only by explicit policy");
      const auto& bin = cov.at("bins").at(0);
      check(cov.at("bins").size() == 1 && bin.at("coverage") == cv::name(r.kind) && bin.at("sigma") == cv::name(r.sigma) &&
            bin.at("lower_eV_exclusive").get<double>() < kin.target_frame_energy_eV &&
            (bin.at("upper_eV_inclusive").is_string() || kin.target_frame_energy_eV <= bin.at("upper_eV_inclusive").get<double>()) &&
            bin.at("energy_min_eV") == kin.target_frame_energy_eV, L + ": energy bin and class recorded");
      check(flow.can_publish(trial) && flow.try_publish(trial) && weight(*c.arena, "electron") == 1 && weight(*c.arena, "H2O") == 16 &&
            weight(*c.arena, "H2O+") == 0, L + ": publication changes no species inventory");
      const auto acc = flow.accepted_coverage();
      check(acc.at("cumulative").at("accepted_steps") == 1 && acc.at("cumulative").at("steps_with_exclusions") == 1 &&
            acc.at("last_step").at("exposures").size() == 1 && acc.at("last_step").at("exposures").at(0) == cov,
            L + ": ledger is part of the accepted state");
      receipts.push_back({{"case", L}, {"coverage", cov}});
    }

    // ---------------- 38.3 eV: supported event, products in the gap ----------------
    {
      {
        Fixture old(b, {{1, 38.3, here}}, 16);nw::ParentSubflow flow(*old.arena, grid.counts, b, limits, bridge, 123, 456);
        reject([&] { flow.prepare(interval(*old.arena, flow, 1e3)); }, "38.3 eV: sealed v3 subflow aborts after its first event");
      }
      {
        Fixture c(b, {{1, 38.3, here}}, 16);cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, scope, 123, 456);
        const auto cp = flow.checkpoint();const auto native = m::PicArena::fingerprint(c.pic);
        refusal([&] { flow.prepare(interval(*c.arena, flow, 1e3)); }, "38.3 eV: default limit refuses gap products");
        check(flow.checkpoint() == cp && m::PicArena::fingerprint(c.pic) == native, "38.3 eV: late refusal restores the reacted native population");
      }
      Fixture c(b, {{1, 38.3, here}}, 16);cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, open, 123, 456);
      const double duration = 1e3;
      auto trial = flow.prepare(interval(*c.arena, flow, duration));
      const auto& ev = trial.audit().at("events");const auto& cov = trial.audit().at("coverage");const auto& u = cov.at("unresolved");
      check(ev.size() == 1 && cov.at("segments") == 2 && cov.at("applied").at("candidate_evaluations") == 1,
            "38.3 eV: supported first event, then a terminal segment with only excluded pairs");
      check(flow.try_publish(trial) && weight(*c.arena, "electron") == 2 && weight(*c.arena, "H2O+") == 1 && weight(*c.arena, "H2O") == 15,
            "38.3 eV: supported event published, gap products kept unreacted");
      const double t1 = ev.at(0).at("collision_elapsed_s");
      long double phys = 0;std::size_t gap = 0;
      for (const auto& r : c.arena->rate_population(100)) if (r.state == "electron") {
        pc::Packet e{"e", "electron", 0, r.native.weight, w::me, r.native.proper_velocity_m_s}, n{"n", "H2O", 0, 15, mn, {}};
        const auto k = pc::relative(e, n);gap += scope.classify(k.target_frame_energy_eV).kind == cv::Kind::onset_gap;
        phys += static_cast<long double>(r.native.weight) * 15 * k.moller_flux_m_s / V * (static_cast<long double>(duration) - t1);
      }
      check(gap == 2 && u.at("pair_evaluations") == 2 && cov.at("distinct_excluded_electron_carriers") == 2 &&
            value(cov.at("distinct_electron_weight_sum_at_first_exclusion")) == 2 && value(u.at("water_weight_sum")) == 30,
            "38.3 eV: both product electrons recorded exactly in the onset gap");
      check(std::abs(value(u.at("physical_exposure_inv_m2")) - static_cast<double>(phys)) <= 1e-13 * static_cast<double>(phys),
            "38.3 eV: excluded exposure equals the independent post-event integral");
      check(value(cov.at("applied").at("expected_physical_reactions")) > 0, "38.3 eV: applied expected reactions recorded as the limit denominator");
      receipts.push_back({{"case", "38.3 eV post-reaction onset gap"}, {"coverage", cov}});
    }

    // ---------------- mixed populations: supported events identical to v3 ----------------
    for (const double Ex : {13., 8000.}) {
      const std::string L = "mixed 100 eV + " + std::to_string(static_cast<int>(Ex)) + " eV";
      const double duration = 1e8;
      Fixture ref(b, {{1, 100, here}}, 16, {{.5, Ex, elsewhere}});
      nw::ParentSubflow old(*ref.arena, grid.counts, b, limits, bridge, 123, 456);
      auto old_trial = old.prepare(interval(*ref.arena, old, duration));
      check(old.try_publish(old_trial), L + ": v3 reference with the extra electron in a disjoint cell");
      {
        Fixture same(b, {{1, 100, here}}, 16, {{.5, Ex, here}});
        nw::ParentSubflow v3(*same.arena, grid.counts, b, limits, bridge, 123, 456);
        reject([&] { v3.prepare(interval(*same.arena, v3, duration)); }, L + ": v3 aborts the whole mixed exposure");
      }
      Fixture mix(b, {{1, 100, here}}, 16, {{.5, Ex, here}});
      cv::ParentSubflow flow(*mix.arena, grid.counts, b, limits, bridge, open, 123, 456);
      auto trial = flow.prepare(interval(*mix.arena, flow, duration));
      check(trial.audit().at("events").size() >= 2 && events_without_native(trial.audit()) == events_without_native(old_trial.audit()),
            L + ": supported event physics identical to v3");
      check(flow.try_publish(trial), L + ": coverage subflow publishes");
      check(rng_fields(flow.checkpoint()) == rng_fields(old.checkpoint()), L + ": scheduler and product RNG streams identical to v3");
      check(census(*mix.arena, {"x0"}) == census(*ref.arena, {"x0"}), L + ": every carrier identity, weight and proper velocity identical to v3");
      const auto& cov = trial.audit().at("coverage");const auto& u = cov.at("unresolved");
      const auto segments = cov.at("segments").get<std::uint64_t>();
      check(segments == trial.audit().at("events").size() + 1 && u.at("pair_evaluations") == segments &&
            cov.at("distinct_excluded_electron_carriers") == 1 && value(cov.at("distinct_electron_weight_sum_at_first_exclusion")) == .5 &&
            value(u.at("electron_weight_sum")) == .5 * segments, L + ": excluded electron recorded exactly in every segment");
      // Independent integral: water weight per segment from the event weights.
      pc::Packet e{"e", "electron", 0, .5, w::me, electron_u(Ex)}, n{"n", "H2O", 0, 1, mn, {}};
      const auto k = pc::relative(e, n);
      long double phys = 0, water = 16, last = 0;
      for (const auto& x : trial.audit().at("events")) {
        const double t = x.at("collision_elapsed_s");phys += .5L * water * k.moller_flux_m_s / V * (t - last);
        water -= x.at("weight").get<double>();last = t;
      }
      phys += .5L * water * k.moller_flux_m_s / V * (static_cast<long double>(duration) - last);
      check(std::abs(value(u.at("physical_exposure_inv_m2")) - static_cast<double>(phys)) <= 1e-12 * static_cast<double>(phys),
            L + ": excluded exposure equals the independent segment integral");
      receipts.push_back({{"case", L}, {"events", events_without_native(trial.audit())}, {"coverage", cov}});
    }

    // ---------------- supported-only parity, including the compensated reservoir ----------------
    for (const auto& [we, nwgt, duration] : std::vector<std::tuple<double, double, double>>{{1, 16, 1e8}, {.25, 1e20, 1e-5}}) {
      const std::string L = "supported-only parity w_e=" + std::to_string(we);
      Fixture a(b, {{we, 100, here}}, nwgt), c(b, {{we, 100, here}}, nwgt);
      nw::ParentSubflow old(*a.arena, grid.counts, b, limits, bridge, 123, 456);
      cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, scope, 123, 456);
      auto ot = old.prepare(interval(*a.arena, old, duration));auto nt = flow.prepare(interval(*c.arena, flow, duration));
      check(events_without_native(ot.audit()) == events_without_native(nt.audit()) && nt.audit().at("events").size() >= 2,
            L + ": event physics identical under the default policy");
      check(old.try_publish(ot) && flow.try_publish(nt) && rng_fields(old.checkpoint()) == rng_fields(flow.checkpoint()) &&
            census(*a.arena) == census(*c.arena), L + ": RNG streams and native carriers identical");
      const auto& cov = nt.audit().at("coverage");
      check(cov.at("unresolved").at("pair_evaluations") == 0 && cov.at("bins").empty() &&
            value(cov.at("applied").at("expected_events")) > 0 && cov.at("exposure_limit").at("refused") == false,
            L + ": empty unresolved ledger with applied measure recorded");
    }

    // ---------------- hard limit branches with the default policy ----------------
    {
      // Unresolved weight comparable to supported weight: both bounds exceeded.
      Fixture c(b, {{1, 100, here}}, 1e20, {{1, 13, here}});
      cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, scope, 123, 456);
      const auto cp = flow.checkpoint();const auto native = m::PicArena::fingerprint(c.pic);const auto prefix = c.arena->state_prefix();
      const auto r = refusal([&] { flow.prepare(interval(*c.arena, flow, 3e-17)); }, "hard limit: comparable gap weight refuses");
      const auto& d = r.exposure_ledger.at("exposure_limit");
      check(d.at("absolute_exceeded") == true && d.at("fraction_exceeded") == true && d.at("fraction_upper").get<double>() > .01 &&
            d.at("events_upper").get<double>() > .01, "hard limit: refusal records both exceeded upper bounds");
      check(flow.checkpoint() == cp && m::PicArena::fingerprint(c.pic) == native && c.arena->state_prefix() == prefix,
            "hard limit: refused exposure leaves native, RNG and ledger unchanged");
      receipts.push_back({{"case", "hard-limit refusal"}, {"refusal", {{"what", r.what()}, {"criterion", r.criterion}, {"ledger", r.exposure_ledger}}}});
    }
    {
      // Small unresolved weight: event bound exceeded but fraction below 1%.
      Fixture c(b, {{1, 100, here}}, 1e20, {{1e-3, 13, here}});
      cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, scope, 123, 456);
      auto t = flow.prepare(interval(*c.arena, flow, 3e-17));const auto& d = t.audit().at("coverage").at("exposure_limit");
      check(d.at("refused") == false && d.at("absolute_exceeded") == true && d.at("fraction_exceeded") == false, "hard limit: fraction branch admits");
      check(flow.try_publish(t), "hard limit: fraction-admitted exposure publishes");
    }
    {
      // Short exposure: fraction is one but the expected excluded events are tiny.
      Fixture c(b, {{1, 13, here}}, 16);cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, scope, 1, 2);
      auto t = flow.prepare(interval(*c.arena, flow, 1.));const auto& d = t.audit().at("coverage").at("exposure_limit");
      check(d.at("refused") == false && d.at("absolute_exceeded") == false && d.at("fraction_upper") == 1.0, "hard limit: event-bound branch admits");
      check(flow.try_publish(t), "hard limit: event-bound admitted exposure publishes");
    }
    {
      // Cumulative: each exposure admitted, the history bound is not.
      cv::CoveragePolicy p;p.max_cumulative_unresolved_events = .002;
      const cv::CoverageScope cum(data, b.dataset_sha256, bridge, p);
      Fixture c(b, {{1, 13, here}}, 16);cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, cum, 1, 2);
      std::uint64_t k = 0;bool refused = false;J last;
      for (; k < 10 && !refused; ++k) {
        try { auto t = flow.prepare(interval(*c.arena, flow, 1., k));check(flow.try_publish(t), "cumulative: admitted exposure publishes"); }
        catch (const cv::CoverageRefusal& r) { refused = true;last = r.exposure_ledger; }
      }
      check(refused && k == 3 && last.at("exposure_limit").at("refused") == false && last.at("cumulative_limit").at("refused") == true,
            "hard limit: cumulative bound refuses the third exposure");
      check(flow.accepted_coverage().at("cumulative").at("accepted_exposures") == 2, "cumulative: refused exposure not accumulated");
    }

    // ---------------- restart and rollback through the chemistry checkpoint ----------------
    {
      Fixture c(b, {{1, 100, here}}, 1e20, {{1e-3, 13, here}});
      cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, open, 77, 78);
      auto t0 = flow.prepare(interval(*c.arena, flow, 5e-17, 0));check(flow.try_publish(t0), "restart: first exposure with exclusions");
      c.arena->save(dir + "/native-k1");const auto k1 = flow.checkpoint();
      std::ofstream(dir + "/subflow-k1.json") << k1.dump(2);
      check(k1.at("coverage").at("cumulative").at("unresolved").at("pair_evaluations").get<std::uint64_t>() > 0, "restart: checkpoint carries the ledger");
      auto t1 = flow.prepare(interval(*c.arena, flow, 5e-17, 1));
      {
        const auto unchanged = flow.checkpoint();
        { auto discarded = flow.prepare(interval(*c.arena, flow, 5e-17, 1));(void)discarded; }
        check(flow.checkpoint() == unchanged, "rollback: unpublished preparation leaves the ledger unchanged");
      }
      check(flow.try_publish(t1), "restart: second exposure");
      const auto k2 = flow.checkpoint();
      check(k2.at("coverage") != k1.at("coverage") && k2.at("coverage").at("cumulative").at("accepted_steps") == 2, "restart: ledger advanced");
      // Restore the k1 native snapshot and k1 chemistry checkpoint (rollback to k1).
      m::GapPic back{options()};
      auto a1 = m::PicArena::restore(back, dir + "/native-k1", states(), channels(b), "native-water-subflow-fixture/v1");
      cv::ParentSubflow r1(*a1, grid.counts, b, limits, bridge, open, 9, 9);
      reject([&] { r1.restore(k2); }, "rollback: newer ledger refused against the older native history");
      r1.restore(k1);check(r1.checkpoint() == k1, "rollback: older checkpoint and ledger restored exactly");
      auto tampered = k1;tampered["coverage"]["cumulative"]["accepted_exposures"] = 7;
      reject([&] { r1.restore(tampered); }, "restart: tampered ledger counter refused");
      tampered = k1;tampered["coverage"]["cumulative"]["unresolved"]["electron_weight_sum"]["mantissa_hex"] = "3";
      reject([&] { r1.restore(tampered); }, "restart: tampered exact weight refused");
      cv::ParentSubflow r2(*a1, grid.counts, b, limits, bridge, scope, 9, 9);
      reject([&] { r2.restore(k1); }, "restart: ledger of another coverage policy refused");
      // Replaying the second exposure from the restored state reproduces it.
      auto replay = r1.prepare(interval(*a1, r1, 5e-17, 1));
      check(replay.audit().at("coverage") == t1.audit().at("coverage") && events_without_native(replay.audit()) == events_without_native(t1.audit()),
            "restart: next exposure and its ledger byte-equivalent after restart");
      check(r1.try_publish(replay) && r1.checkpoint().at("coverage") == k2.at("coverage") && rng_fields(r1.checkpoint()) == rng_fields(k2),
            "restart: accepted ledger and RNG identical after replay");
      // A refused exposure after restart leaves the restored state intact.
      Fixture d(b, {{1, 13, here}}, 1e20);cv::ParentSubflow rf(*d.arena, grid.counts, b, limits, bridge, scope, 1, 1);
      const auto cp = rf.checkpoint();refusal([&] { rf.prepare(interval(*d.arena, rf, 1e-3)); }, "rollback: refusal raised");
      check(rf.checkpoint() == cp, "rollback: refused exposure leaves the checkpoint unchanged");
    }

    // ---------------- full-time first exposure ----------------
    {
      const auto ft_rng = cv::ParentSubflow::initial_full_time_rng_sha256(options(), grid.counts, b, limits, scope, 5, 6);
      const auto v3_rng = nw::ParentSubflow::initial_full_time_rng_sha256(options(), grid.counts, b, limits, 5, 6);
      check(ft_rng != v3_rng && ft_rng != cv::ParentSubflow::initial_full_time_rng_sha256(options(), grid.counts, b, limits, open, 5, 6),
            "full-time genesis RNG identity binds the coverage scope and policy");
      auto make = [&](double water, const cv::CoverageScope& s) {
        m::FullTimeGenesisBinding g;g.time_s = 7.4978805549199997e-10;g.origin_provenance = "authored/coverage-fixture-no-simulation";
        g.field_sha256 = nw::sha("coverage-fixture/field");g.mesh_sha256 = nw::sha("coverage-fixture/mesh");
        g.material_sha256 = nw::sha("coverage-fixture/material");g.frame_sha256 = nw::sha("coverage-fixture/frame");
        g.common_clock_sha256 = nw::sha("coverage-fixture/clock");g.heat_sha256 = nw::sha("coverage-fixture/heat");
        g.source_sha256 = nw::sha("coverage-fixture/source");
        g.rng_sha256 = cv::ParentSubflow::initial_full_time_rng_sha256(options(), grid.counts, b, limits, s, 5, 6);
        return m::PicArena::create_full_time(options(), states(), channels(b),
          {{"e", "electron", "authored/fixture", {1e10, 0}, here, electron_u(13)},
           {"n", "H2O", "authored/fixture", {water, 0}, water_at, {0, 0, 0}}}, g);
      };
      auto step = [&](m::PicArena& a) {
        const auto bd = a.full_time_boundary();
        m::FullTimeStepBinding c{bd->accepted_step, bd->accepted_arena_revision, bd->time_s, bd->time_s + 5e-15, 5e-15,
                                 bd->accepted_arena_sha256, bd->field_sha256, bd->mesh_sha256, bd->material_sha256, nw::sha("coverage-common")};
        return c;
      };
      for (const auto& [water, admitted] : std::vector<std::pair<double, bool>>{{1e16, true}, {1e21, false}}) {
        auto a = make(water, scope);
        cv::ParentSubflow flow(*a, grid.counts, b, limits, bridge, scope, 5, 6, true);
        const auto cp = flow.checkpoint();const auto rng = flow.accepted_rng_sha256();
        {
          const auto common = step(*a);auto outer = a->begin_full_time_step(common, 100000, 4096);
          m::ExposureBinding e{common, m::FullTimeSlot::before_transport, common.begin_s, .5 * common.dt_s, 0, common.accepted_field_sha256};
          auto first_scope = outer.first_exposure(e);
          if (admitted) {
            auto stage = flow.prepare_on(first_scope);
            const auto& cov = stage.audit().at("coverage");
            check(cov.at("slot") == "before_transport" && cov.at("unresolved").at("pair_evaluations") == 1 &&
                  cov.at("exposure_limit").at("refused") == false && stage.audit().at("events").empty(),
                  "full-time first exposure: excluded pair recorded, exposure proceeds");
            check(flow.staged_rng_sha256(stage) != rng && pc::digest(flow.staged_rng_sha256(stage)),
                  "full-time first exposure: staged ledger bound into the staged RNG/clock identity");
            receipts.push_back({{"case", "full-time first exposure"}, {"coverage", cov}});
          } else {
            const auto r = refusal([&] { (void)flow.prepare_on(first_scope); }, "full-time first exposure: hard limit refuses");
            check(r.exposure_ledger.at("slot") == "before_transport", "full-time refusal identifies the exposure slot");
          }
        }
        check(flow.checkpoint() == cp && flow.accepted_rng_sha256() == rng, "full-time: abandoned private step leaves the accepted ledger unchanged");
      }
      auto a = make(1e16, scope);
      reject([&] { cv::ParentSubflow bad(*a, grid.counts, b, limits, bridge, open, 5, 6, true);(void)bad.checkpoint(); },
             "full-time genesis bound to another coverage policy refused");
    }

    J result{{"status", "passed"}, {"checks", checks}, {"check_count", checks.size()}, {"coverage_scope", scope.description()},
             {"coverage_scope_sha256", scope.identity()}, {"record_only_scope_sha256", open.identity()}, {"coverage_bins", coverage_map},
             {"receipts", receipts}, {"conditional_parent_only", true}, {"physical_trajectory_advanced", false}, {"full_coupled_run", false}};
    std::ofstream(dir + "/result.json") << result.dump(2) << '\n';
    std::cout << "Passed " << checks.size() << " checks\n";
  } catch (const std::exception& e) { std::cerr << e.what() << '\n';return 1; }
}
