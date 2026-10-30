#pragma once
#include "exact_sum.hpp"
#include "native_parent_subflow.hpp"
#include <algorithm>
#include <map>
#include <optional>
#include <set>
#include <tuple>
#include <unordered_map>

// Versioned coverage scope for the conditional water parent-ionization channel
// (H2O:ion_yield:H2O+, parent_ground_equal_pair_CM_v1).
//
// A (carrier, H2O) pair whose target-rest energy lies where the parent cross
// section is positive or unknown but no qualified event kernel exists is
// EXCLUDED from the reaction exposure and recorded in an exact unresolved-
// coverage ledger. Other pairs, rates and branches are untouched: nothing is
// renormalized, no rate is invented and no unresolved pair is set to zero
// without a ledger entry. A configurable hard limit refuses the exposure with a
// typed CoverageRefusal when the unresolved part is too large.
namespace raikou::native_water::coverage_v1 {
namespace m = ::raikou::mfem_backend;
namespace pc = ::raikou::pic_collisions;
namespace w = ::raikou::pic_water;
using Json = m::ArenaJson;
using ::raikou::native_water::sha;

inline constexpr const char* scope_schema = "parent-channel-coverage-scope-v1";
inline constexpr const char* ledger_schema = "parent-channel-unresolved-coverage-ledger-v1";
// Tabulated aggregate inelastic cross section of the SAME bound dataset. Parent
// ionization is one inelastic process, so sigma_parent <= sigma_inelastic. It is
// used only as an upper bound for recording, never as a rate.
inline const std::string default_upper_bound_observable = "H2O:munoz2007:recommended_inelastic";

inline Json number_or_inf(double x) {
  if (std::isinf(x) && x > 0) return "inf";
  pc::require(std::isfinite(x), "nonfinite coverage number");
  return x;
}
inline Json exact_json(const ExactSum& s) {
  const double v = s.value();
  return Json{{"mantissa_hex", s.mantissa_hex()}, {"exponent2", s.exponent2()},
              {"value", std::isfinite(v) ? Json(v) : Json(nullptr)}};
}
inline ExactSum exact_from(const Json& j) {
  pc::require(j.is_object() && j.size() == 3 && j.at("mantissa_hex").is_string() &&
                  j.at("exponent2").is_number_integer(),
              "malformed exact ledger value");
  auto s = ExactSum::from_canonical(j.at("mantissa_hex").get<std::string>(), j.at("exponent2").get<long>());
  const double v = s.value();
  pc::require(std::isfinite(v) ? (j.at("value").is_number() && j.at("value").get<double>() == v)
                               : j.at("value").is_null(),
              "exact ledger value differs from its rounded record");
  return s;
}

// Hard-limit policy. Bound into the scope identity, and through it into the
// subflow configuration and every RNG/clock state hash.
struct CoveragePolicy {
  // Upper bound on expected excluded weighted-pair (scheduler) events in one
  // collision exposure. By Poisson thinning this also bounds the probability
  // that the excluded hazard would have changed the exposure's sampled event
  // sequence.
  double max_exposure_unresolved_events = 0.01;
  // Same bound accumulated over all accepted exposures plus the candidate.
  double max_cumulative_unresolved_events = 1.0;
  // Upper-bound share of expected physical parent ionizations left unresolved:
  // U/(A+U) with A applied (qualified) and U the unresolved upper bound.
  double max_unresolved_reaction_fraction = 0.01;
  // Pairs with no data-supported sigma upper bound (above every tabulated
  // bounding observable) make U unbounded. Refused unless explicitly admitted.
  bool admit_unbounded = false;
  Json json() const {
    return {{"max_exposure_unresolved_events", number_or_inf(max_exposure_unresolved_events)},
            {"max_cumulative_unresolved_events", number_or_inf(max_cumulative_unresolved_events)},
            {"max_unresolved_reaction_fraction", max_unresolved_reaction_fraction},
            {"admit_unbounded", admit_unbounded},
            {"refusal_rule", "refuse iff unbounded-and-not-admitted, or (events_upper > event limit AND "
                             "physical_upper*(1-f) > f*applied_physical); evaluated per exposure and cumulatively"}};
  }
  void validate() const {
    for (double x : {max_exposure_unresolved_events, max_cumulative_unresolved_events})
      pc::require(!std::isnan(x) && x >= 0, "coverage event limit must be a nonnegative number or +inf");
    pc::require(std::isfinite(max_unresolved_reaction_fraction) && max_unresolved_reaction_fraction >= 0 &&
                    max_unresolved_reaction_fraction <= 1,
                "coverage reaction fraction limit must lie in [0,1]");
  }
};

enum class Kind { known_zero, supported, onset_gap, below_table, above_table, kernel_domain };
enum class Sigma { none, tabulated, bounded, unbounded };
inline const char* name(Kind k) {
  switch (k) {
  case Kind::known_zero: return "known_zero";
  case Kind::supported: return "supported";
  case Kind::onset_gap: return "onset_gap";
  case Kind::below_table: return "below_table_without_threshold";
  case Kind::above_table: return "above_table";
  case Kind::kernel_domain: return "tabulated_outside_kernel_domain";
  }
  return "";
}
inline const char* name(Sigma s) {
  switch (s) {
  case Sigma::none: return "none";
  case Sigma::tabulated: return "tabulated";
  case Sigma::bounded: return "unknown_upper_bounded";
  case Sigma::unbounded: return "unknown_unbounded";
  }
  return "";
}
inline bool unresolved(Kind k) {
  return k == Kind::onset_gap || k == Kind::below_table || k == Kind::above_table || k == Kind::kernel_domain;
}
struct Classified {
  Kind kind = Kind::known_zero;
  Sigma sigma = Sigma::none;
  double sigma_m2 = 0, upper_m2 = 0;
  std::size_t bin = 0;
};
// Bins are (lower, upper] energy intervals. An energy exactly on an edge can
// carry a different class than the interval interior (closed support ends), so
// a bin is identified by interval, coverage class and sigma status together.
using BinKey = std::tuple<std::size_t, Kind, Sigma>;
inline Kind kind_named(const std::string& s) {
  for (Kind k : {Kind::onset_gap, Kind::below_table, Kind::above_table, Kind::kernel_domain})
    if (s == name(k)) return k;
  throw pc::Error("unknown unresolved coverage class");
}
inline Sigma sigma_named(const std::string& s) {
  for (Sigma k : {Sigma::tabulated, Sigma::bounded, Sigma::unbounded})
    if (s == name(k)) return k;
  throw pc::Error("unknown unresolved sigma status");
}

// Additive unresolved measures. "Physical" quantities use w_e*w_n (physical
// reactions); "event" quantities use max(w_e,w_n) (scheduler events with event
// weight min(w_e,w_n)), exactly as the sealed propensity. Exposures have unit
// m^-2: multiplied by a cross section they give an expected count.
struct Measures {
  std::uint64_t pair_evaluations = 0, tabulated = 0, bounded = 0, unbounded = 0;
  ExactSum electron_weight, water_weight, event_weight;
  ExactSum physical_exposure, event_exposure;
  ExactSum tabulated_physical, tabulated_events;
  ExactSum upper_physical, upper_events;
  ExactSum unbounded_physical_exposure, unbounded_event_exposure;
  double energy_min_eV = 0, energy_max_eV = 0;
  struct Terms {
    double we = 0, wn = 0, phys = 0, ev = 0, tab_phys = 0, tab_ev = 0, up_phys = 0, up_ev = 0, energy = 0;
    Sigma sigma = Sigma::none;
  };
  void add(const Terms& t) {
    pc::require(pair_evaluations < UINT64_MAX, "coverage pair counter exhausted");
    if (!pair_evaluations) energy_min_eV = energy_max_eV = t.energy;
    else { energy_min_eV = std::min(energy_min_eV, t.energy); energy_max_eV = std::max(energy_max_eV, t.energy); }
    ++pair_evaluations;
    electron_weight.add(t.we);water_weight.add(t.wn);event_weight.add(std::min(t.we, t.wn));
    physical_exposure.add(t.phys);event_exposure.add(t.ev);
    if (t.sigma == Sigma::tabulated) {
      ++tabulated;tabulated_physical.add(t.tab_phys);tabulated_events.add(t.tab_ev);
      upper_physical.add(t.up_phys);upper_events.add(t.up_ev);
    } else if (t.sigma == Sigma::bounded) {
      ++bounded;upper_physical.add(t.up_phys);upper_events.add(t.up_ev);
    } else {
      pc::require(t.sigma == Sigma::unbounded, "unresolved pair without sigma status");
      ++unbounded;unbounded_physical_exposure.add(t.phys);unbounded_event_exposure.add(t.ev);
    }
  }
  void merge(const Measures& o) {
    if (!o.pair_evaluations) return;
    pc::require(pair_evaluations <= UINT64_MAX - o.pair_evaluations, "coverage pair counter exhausted");
    if (!pair_evaluations) { energy_min_eV = o.energy_min_eV; energy_max_eV = o.energy_max_eV; }
    else { energy_min_eV = std::min(energy_min_eV, o.energy_min_eV); energy_max_eV = std::max(energy_max_eV, o.energy_max_eV); }
    pair_evaluations += o.pair_evaluations; tabulated += o.tabulated; bounded += o.bounded; unbounded += o.unbounded;
    electron_weight.add(o.electron_weight);water_weight.add(o.water_weight);event_weight.add(o.event_weight);
    physical_exposure.add(o.physical_exposure);event_exposure.add(o.event_exposure);
    tabulated_physical.add(o.tabulated_physical);tabulated_events.add(o.tabulated_events);
    upper_physical.add(o.upper_physical);upper_events.add(o.upper_events);
    unbounded_physical_exposure.add(o.unbounded_physical_exposure);unbounded_event_exposure.add(o.unbounded_event_exposure);
  }
  Json json() const {
    Json j{{"pair_evaluations", pair_evaluations}, {"tabulated_sigma_pair_evaluations", tabulated},
           {"upper_bounded_sigma_pair_evaluations", bounded}, {"unbounded_sigma_pair_evaluations", unbounded},
           {"electron_weight_sum", exact_json(electron_weight)}, {"water_weight_sum", exact_json(water_weight)},
           {"event_weight_sum", exact_json(event_weight)},
           {"physical_exposure_inv_m2", exact_json(physical_exposure)},
           {"event_exposure_inv_m2", exact_json(event_exposure)},
           {"tabulated_expected_physical_reactions", exact_json(tabulated_physical)},
           {"tabulated_expected_events", exact_json(tabulated_events)},
           {"upper_bound_physical_reactions", exact_json(upper_physical)},
           {"upper_bound_events", exact_json(upper_events)},
           {"unbounded_physical_exposure_inv_m2", exact_json(unbounded_physical_exposure)},
           {"unbounded_event_exposure_inv_m2", exact_json(unbounded_event_exposure)},
           {"lower_bound_physical_reactions", 0}};
    if (pair_evaluations) { j["energy_min_eV"] = energy_min_eV; j["energy_max_eV"] = energy_max_eV; }
    else { j["energy_min_eV"] = nullptr; j["energy_max_eV"] = nullptr; }
    return j;
  }
  static Measures from(const Json& j) {
    Measures s;
    s.pair_evaluations = j.at("pair_evaluations").get<std::uint64_t>();
    s.tabulated = j.at("tabulated_sigma_pair_evaluations").get<std::uint64_t>();
    s.bounded = j.at("upper_bounded_sigma_pair_evaluations").get<std::uint64_t>();
    s.unbounded = j.at("unbounded_sigma_pair_evaluations").get<std::uint64_t>();
    pc::require(s.tabulated + s.bounded + s.unbounded == s.pair_evaluations && s.tabulated <= s.pair_evaluations &&
                    s.bounded <= s.pair_evaluations && s.unbounded <= s.pair_evaluations,
                "coverage pair partition inconsistent");
    s.electron_weight = exact_from(j.at("electron_weight_sum"));
    s.water_weight = exact_from(j.at("water_weight_sum"));
    s.event_weight = exact_from(j.at("event_weight_sum"));
    s.physical_exposure = exact_from(j.at("physical_exposure_inv_m2"));
    s.event_exposure = exact_from(j.at("event_exposure_inv_m2"));
    s.tabulated_physical = exact_from(j.at("tabulated_expected_physical_reactions"));
    s.tabulated_events = exact_from(j.at("tabulated_expected_events"));
    s.upper_physical = exact_from(j.at("upper_bound_physical_reactions"));
    s.upper_events = exact_from(j.at("upper_bound_events"));
    s.unbounded_physical_exposure = exact_from(j.at("unbounded_physical_exposure_inv_m2"));
    s.unbounded_event_exposure = exact_from(j.at("unbounded_event_exposure_inv_m2"));
    if (s.pair_evaluations) {
      s.energy_min_eV = j.at("energy_min_eV").get<double>();
      s.energy_max_eV = j.at("energy_max_eV").get<double>();
    }
    return s;
  }
};
// Expected reactions actually applied by the qualified kernel (the denominator
// part A of the limit), from the sealed candidate propensities.
struct Applied {
  std::uint64_t candidate_evaluations = 0;
  ExactSum physical, events;
  void merge(const Applied& o) {
    pc::require(candidate_evaluations <= UINT64_MAX - o.candidate_evaluations, "applied counter exhausted");
    candidate_evaluations += o.candidate_evaluations; physical.add(o.physical); events.add(o.events);
  }
  Json json() const {
    return {{"candidate_evaluations", candidate_evaluations}, {"expected_physical_reactions", exact_json(physical)},
            {"expected_events", exact_json(events)}};
  }
  static Applied from(const Json& j) {
    Applied a;a.candidate_evaluations = j.at("candidate_evaluations").get<std::uint64_t>();
    a.physical = exact_from(j.at("expected_physical_reactions"));a.events = exact_from(j.at("expected_events"));
    return a;
  }
};
struct Decision {
  std::string scope;
  double limit_events = 0, limit_fraction = 0, events_upper = 0, physical_upper = 0, applied_physical = 0,
         fraction_upper = 0;
  std::uint64_t unbounded_pairs = 0;
  bool admit_unbounded = false, absolute_exceeded = false, fraction_exceeded = false, refused = false;
  std::string criterion = "none";
  Json json() const {
    return {{"evaluated_over", scope}, {"event_limit", number_or_inf(limit_events)}, {"fraction_limit", limit_fraction},
            {"events_upper", number_or_inf(events_upper)}, {"physical_upper", number_or_inf(physical_upper)},
            {"applied_physical", number_or_inf(applied_physical)}, {"fraction_upper", fraction_upper},
            {"unbounded_sigma_pair_evaluations", unbounded_pairs}, {"unbounded_admitted", admit_unbounded},
            {"absolute_exceeded", absolute_exceeded}, {"fraction_exceeded", fraction_exceeded},
            {"refused", refused}, {"criterion", criterion}};
  }
};
inline Decision evaluate(const CoveragePolicy& p, const Measures& u, const Applied& a, double event_limit,
                         const char* over) {
  Decision d;d.scope = over;d.limit_events = event_limit;d.limit_fraction = p.max_unresolved_reaction_fraction;
  d.admit_unbounded = p.admit_unbounded;d.unbounded_pairs = u.unbounded;
  d.events_upper = u.upper_events.value();d.physical_upper = u.upper_physical.value();
  d.applied_physical = a.physical.value();
  const long double U = d.physical_upper, A = d.applied_physical, f = d.limit_fraction;
  d.fraction_upper = U == 0 ? 0. : (std::isinf(U) ? 1. : static_cast<double>(U / (A + U)));
  d.absolute_exceeded = static_cast<long double>(d.events_upper) > static_cast<long double>(event_limit);
  d.fraction_exceeded = U > 0 && (std::isinf(U) || U * (1 - f) > f * A);
  if (u.unbounded && !p.admit_unbounded) { d.refused = true; d.criterion = "unbounded_unresolved_sigma"; }
  else if (d.absolute_exceeded && d.fraction_exceeded) { d.refused = true; d.criterion = "events_and_fraction_upper_bounds"; }
  return d;
}

enum class Slot { legacy, before_transport, after_transport };
inline const char* name(Slot s) {
  return s == Slot::legacy ? "legacy_frozen_stage" : (s == Slot::before_transport ? "before_transport" : "after_transport");
}

struct ExposureLedger {
  Slot slot = Slot::legacy;
  std::uint64_t sequence = 0, segments = 0;
  double duration_s = 0;
  Measures unresolved;
  Applied applied;
  std::map<BinKey, Measures> bins;
  std::uint64_t distinct_pairs = 0, distinct_electrons = 0, distinct_waters = 0;
  ExactSum distinct_electron_weight, distinct_water_weight;
  Decision exposure_decision, cumulative_decision;
};
struct Cumulative {
  std::uint64_t accepted_steps = 0, steps_with_exclusions = 0, accepted_exposures = 0, exposures_with_exclusions = 0;
  Measures unresolved;
  Applied applied;
  std::map<BinKey, Measures> bins;
};
struct Ledger {
  std::optional<ExposureLedger> pending; // staged first full-time exposure only; never accepted
  std::vector<ExposureLedger> last_step;
  Cumulative cumulative;
};

struct CoverageRefusal : pc::Error {
  std::string criterion, scope_identity;
  Json exposure_ledger;
  CoverageRefusal(const std::string& what, std::string c, std::string s, Json l)
      : pc::Error(what), criterion(std::move(c)), scope_identity(std::move(s)), exposure_ledger(std::move(l)) {}
};

class CoverageScope {
  std::shared_ptr<const w::Data> data_;
  std::string dataset_sha_, kernel_sha_, bound_id_, identity_;
  CoveragePolicy policy_;
  bool zero_present_ = false;
  double zero_below_ = 0, table_lo_ = 0, table_hi_ = 0, kernel_lo_ = 0, kernel_hi_ = 0, bound_lo_ = 0, bound_hi_ = 0,
         bound_scale_ = 0;
  std::vector<double> bound_e_, bound_y_, edges_;
  std::vector<Classified> bin_class_;
  std::vector<double> probes_;
  Json description_;

  static double num(const Json& j) {
    pc::require(j.is_number(), "coverage metadata number");
    const double x = j.get<double>();pc::require(std::isfinite(x), "coverage metadata finite");return x;
  }
  double envelope(double E) const {
    if (E == bound_e_.back()) return bound_y_.back() * bound_scale_;
    const auto k = static_cast<std::size_t>(std::upper_bound(bound_e_.begin(), bound_e_.end(), E) - bound_e_.begin() - 1);
    return std::max(bound_y_[k], bound_y_[k + 1]) * bound_scale_;
  }
  Classified classify_unbinned(double E) const {
    Classified c;
    pc::require(std::isfinite(E) && E >= 0, "coverage classification needs a finite nonnegative energy");
    auto unknown = [&](Kind k) {
      c.kind = k;
      if (E >= bound_lo_ && E <= bound_hi_) { c.sigma = Sigma::bounded; c.upper_m2 = envelope(E); }
      else c.sigma = Sigma::unbounded;
    };
    // Same order as StrictCollisionDataset::observable_value_SI and
    // require_MCC_event_kernel; the bridge probe checks the correspondence.
    if (E > table_hi_) unknown(Kind::above_table);
    else if (zero_present_ && E <= zero_below_) c.kind = Kind::known_zero;
    else if (E < table_lo_) unknown(zero_present_ ? Kind::onset_gap : Kind::below_table);
    else {
      const double s = data_->observable_value_SI(w::parent_channel, E);
      if (!(s > 0)) c.kind = Kind::known_zero;
      else if (E >= kernel_lo_ && E <= kernel_hi_) { c.kind = Kind::supported; c.sigma = Sigma::tabulated; c.sigma_m2 = s; }
      else { c.kind = Kind::kernel_domain; c.sigma = Sigma::tabulated; c.sigma_m2 = s; c.upper_m2 = s; }
    }
    return c;
  }

public:
  CoverageScope(std::shared_ptr<const w::Data> data, std::string dataset_sha256, const w::ParentBridge& bridge,
                CoveragePolicy policy = {}, std::string bound_observable = default_upper_bound_observable)
      : data_(std::move(data)), dataset_sha_(std::move(dataset_sha256)), bound_id_(std::move(bound_observable)),
        policy_(policy) {
    pc::require(bool(data_) && pc::digest(dataset_sha_), "coverage scope needs bound dataset bytes");
    policy_.validate();
    const auto channel = bridge.channel();
    pc::require(channel.id == w::parent_channel && channel.data_sha256 == dataset_sha_, "coverage scope and bridge bind different data");
    kernel_sha_ = channel.kernel_sha256;
    const auto& p = data_->channel_metadata(w::parent_channel);
    pc::require(p.at("model") == "linear_table" && p.at("unit_SI") == "m2" && p.at("outside_support") == "error",
                "parent coverage requires a strict tabulated area");
    table_lo_ = num(p.at("support_eV")[0]);table_hi_ = num(p.at("support_eV")[1]);
    zero_present_ = !p.at("zero_below_eV").is_null();
    if (zero_present_) zero_below_ = num(p.at("zero_below_eV"));
    const auto status = data_->kernel_status(w::parent_channel, w::parent_model);
    pc::require(status.at("status") == "conditional_available" && status.at("model").at("products").is_object(),
                "parent kernel must be explicitly available with products");
    kernel_lo_ = num(status.at("model").at("data_support_eV")[0]);kernel_hi_ = num(status.at("model").at("data_support_eV")[1]);
    pc::require(table_lo_ < table_hi_ && kernel_lo_ <= kernel_hi_ && (!zero_present_ || zero_below_ <= table_lo_),
                "parent coverage support ordering");
    const auto& b = data_->channel_metadata(bound_id_);
    pc::require(b.at("model") == "linear_table" && b.at("unit_SI") == "m2" && b.at("observable") == "aggregate_inelastic",
                "upper-bound observable must be a tabulated aggregate inelastic area");
    bound_lo_ = num(b.at("support_eV")[0]);bound_hi_ = num(b.at("support_eV")[1]);bound_scale_ = num(b.at("ordinate_scale_to_SI"));
    for (const auto& x : b.at("energy_eV")) bound_e_.push_back(num(x));
    for (const auto& x : b.at("values")) bound_y_.push_back(num(x));
    pc::require(bound_e_.size() == bound_y_.size() && bound_e_.size() >= 2 && bound_e_.front() == bound_lo_ &&
                    bound_e_.back() == bound_hi_ && bound_scale_ > 0,
                "upper-bound table shape");
    // Plausibility of the bound on its overlap with the qualified parent table.
    for (const auto& x : p.at("energy_eV")) {
      const double E = num(x);
      if (E >= bound_lo_ && E <= bound_hi_)
        pc::require(data_->observable_value_SI(w::parent_channel, E) <= envelope(E),
                    "upper-bound observable does not dominate the parent table");
    }
    // Versioned energy bins (eV). Class boundaries are always edges, so every
    // bin has exactly one coverage class and sigma status.
    std::set<double> e{0., 12.75, 13., 13.25, 1500., 2000., 3000., 5000., 8000., 20000., 50000., 1e5, 1e6,
                       table_lo_, table_hi_, kernel_lo_, kernel_hi_, bound_lo_, bound_hi_};
    if (zero_present_) e.insert(zero_below_);
    edges_.assign(e.begin(), e.end());
    edges_.push_back(std::numeric_limits<double>::infinity());
    for (std::size_t i = 0; i + 1 < edges_.size(); ++i) {
      const double lo = edges_[i], hi = edges_[i + 1];
      const double mid = std::isinf(hi) ? 2 * lo : (i == 0 ? hi / 2 : lo + (hi - lo) / 2);
      auto c = classify_unbinned(mid);c.bin = i;c.sigma_m2 = c.upper_m2 = 0;bin_class_.push_back(c);
    }
    // Bind the scope to the sealed bridge behavior: every probe classified
    // known-zero returns 0, supported returns the same positive sigma, and
    // unresolved raises the strict coverage error.
    std::set<double> probes{0., 1., 10., 12.7, 13.4999, 37.9, 38.3, 39.6, 100., 999.999, 1000.0001, 8000., 20000.};
    for (std::size_t i = 0; i + 1 < edges_.size(); ++i) {
      const double lo = edges_[i], hi = edges_[i + 1];
      probes.insert(lo);
      if (lo > 0) probes.insert(std::nextafter(lo, 0.));
      probes.insert(std::nextafter(lo, std::numeric_limits<double>::infinity()));
      if (!std::isinf(hi)) probes.insert(lo + (hi - lo) / 2);
    }
    for (double E : probes) {
      const auto c = classify_unbinned(E);
      bool coverage_error = false;double s = -1;
      try { s = channel.sigma_m2(E); } catch (const ::raikou::kinetics::CollisionCoverageError&) { coverage_error = true; }
      if (c.kind == Kind::known_zero) pc::require(!coverage_error && s == 0, "coverage probe: known-zero energy differs from bridge");
      else if (c.kind == Kind::supported)
        pc::require(!coverage_error && s > 0 && s == c.sigma_m2, "coverage probe: supported energy differs from bridge");
      else pc::require(coverage_error, "coverage probe: unresolved energy is accepted by the bridge");
      probes_.push_back(E);
    }
    Json bins = Json::array();
    for (std::size_t i = 0; i + 1 < edges_.size(); ++i)
      bins.push_back({{"bin", i}, {"lower_eV_exclusive", edges_[i]}, {"upper_eV_inclusive", number_or_inf(edges_[i + 1])},
                      {"coverage", name(bin_class_[i].kind)}, {"sigma", name(bin_class_[i].sigma)}});
    description_ = {{"schema", scope_schema}, {"channel", w::parent_channel}, {"model", w::parent_model},
      {"dataset_sha256", dataset_sha_}, {"kernel_sha256", kernel_sha_},
      {"table_support_eV", {table_lo_, table_hi_}}, {"known_zero_at_or_below_eV", zero_present_ ? Json(zero_below_) : Json(nullptr)},
      {"kernel_support_eV", {kernel_lo_, kernel_hi_}}, {"kernel_products", status.at("model").at("products")},
      {"upper_bound", {{"observable", bound_id_}, {"source", b.at("source")}, {"support_eV", {bound_lo_, bound_hi_}},
                       {"energy_eV", bound_e_}, {"values", bound_y_}, {"ordinate_scale_to_SI", bound_scale_},
                       {"envelope", "max of the two bracketing tabulated values; sigma_parent <= sigma_inelastic"}}},
      {"bins", bins}, {"probe_energies_eV", probes_}, {"policy", policy_.json()},
      {"exclusion", "pair excluded from the exposure and recorded; no renormalization, no rate invention"}};
    identity_ = sha(description_.dump());
  }
  const std::string& identity() const noexcept { return identity_; }
  const Json& description() const noexcept { return description_; }
  const CoveragePolicy& policy() const noexcept { return policy_; }
  const std::string& dataset_sha256() const noexcept { return dataset_sha_; }
  const std::string& kernel_sha256() const noexcept { return kernel_sha_; }
  std::size_t bin_count() const noexcept { return bin_class_.size(); }
  std::size_t bin_of(double E) const {
    auto i = static_cast<std::size_t>(std::lower_bound(edges_.begin(), edges_.end(), E) - edges_.begin());
    return i == 0 ? 0 : i - 1;
  }
  Classified classify(double E) const {
    auto c = classify_unbinned(E);c.bin = bin_of(E);
    // A pair can only fall in an interval whose closed hull contains it, and an
    // unresolved class never occurs at a supported interior.
    pc::require(c.bin + 1 < edges_.size() && (c.bin == 0 ? E >= edges_[0] : E > edges_[c.bin]) && E <= edges_[c.bin + 1],
                "coverage bin lookup");
    return c;
  }
  // Fast exclusion predicate for the candidate query. Invalid energies are
  // never excluded here, so the sealed bridge still refuses them.
  bool excluded(double E) const {
    if (!(std::isfinite(E) && E >= 0)) return false;
    if (E > table_hi_) return true;
    if (zero_present_ && E <= zero_below_) return false;
    if (E < table_lo_) return true;
    if (E >= kernel_lo_ && E <= kernel_hi_) return false;
    return data_->observable_value_SI(w::parent_channel, E) > 0;
  }
  // The scope must bind the same bridge identity the subflow uses.
  void require_bridge(const w::ParentBridge& bridge) const {
    const auto channel = bridge.channel();
    pc::require(channel.data_sha256 == dataset_sha_ && channel.kernel_sha256 == kernel_sha_, "coverage scope bound to another bridge");
    for (double E : probes_) {
      const auto c = classify_unbinned(E);bool coverage_error = false;double s = -1;
      try { s = channel.sigma_m2(E); } catch (const ::raikou::kinetics::CollisionCoverageError&) { coverage_error = true; }
      pc::require(c.kind == Kind::known_zero ? (!coverage_error && s == 0)
                  : c.kind == Kind::supported ? (!coverage_error && s == c.sigma_m2) : coverage_error,
                  "coverage scope classification differs from the subflow bridge");
    }
  }

  // ---- ledger serialization (bin metadata comes from this scope) ----
  Json bins_json(const std::map<BinKey, Measures>& bins) const {
    Json out = Json::array();
    for (const auto& [key, s] : bins) {
      const auto& [i, kind, sigma] = key;
      pc::require(i + 1 < edges_.size() && unresolved(kind) && sigma != Sigma::none && s.pair_evaluations > 0,
                  "coverage bin record");
      auto j = s.json();j["bin"] = i;j["lower_eV_exclusive"] = edges_[i];j["upper_eV_inclusive"] = number_or_inf(edges_[i + 1]);
      j["coverage"] = name(kind);j["sigma"] = name(sigma);
      out.push_back(std::move(j));
    }
    return out;
  }
  std::map<BinKey, Measures> bins_from(const Json& j) const {
    std::map<BinKey, Measures> out;
    pc::require(j.is_array() && j.size() <= 4 * bin_class_.size(), "coverage bins");
    for (const auto& b : j) {
      const BinKey key{b.at("bin").get<std::size_t>(), kind_named(b.at("coverage").get<std::string>()),
                       sigma_named(b.at("sigma").get<std::string>())};
      pc::require(std::get<0>(key) + 1 < edges_.size() && !out.count(key), "coverage bin index");
      out.emplace(key, Measures::from(b));
    }
    return out;
  }
  Json exposure_json(const ExposureLedger& l) const {
    return {{"slot", name(l.slot)}, {"sequence", l.sequence}, {"segments", l.segments}, {"duration_s", l.duration_s},
            {"unresolved", l.unresolved.json()}, {"applied", l.applied.json()}, {"bins", bins_json(l.bins)},
            {"distinct_excluded_pairs", l.distinct_pairs}, {"distinct_excluded_electron_carriers", l.distinct_electrons},
            {"distinct_excluded_water_carriers", l.distinct_waters},
            {"distinct_electron_weight_sum_at_first_exclusion", exact_json(l.distinct_electron_weight)},
            {"distinct_water_weight_sum_at_first_exclusion", exact_json(l.distinct_water_weight)},
            {"exposure_limit", l.exposure_decision.json()}, {"cumulative_limit", l.cumulative_decision.json()}};
  }
  ExposureLedger exposure_from(const Json& j) const {
    ExposureLedger l;
    const auto slot = j.at("slot").get<std::string>();
    l.slot = slot == name(Slot::legacy) ? Slot::legacy : slot == name(Slot::before_transport) ? Slot::before_transport : Slot::after_transport;
    pc::require(slot == name(l.slot), "coverage exposure slot");
    l.sequence = j.at("sequence").get<std::uint64_t>();l.segments = j.at("segments").get<std::uint64_t>();
    l.duration_s = j.at("duration_s").get<double>();
    l.unresolved = Measures::from(j.at("unresolved"));l.applied = Applied::from(j.at("applied"));l.bins = bins_from(j.at("bins"));
    l.distinct_pairs = j.at("distinct_excluded_pairs").get<std::uint64_t>();
    l.distinct_electrons = j.at("distinct_excluded_electron_carriers").get<std::uint64_t>();
    l.distinct_waters = j.at("distinct_excluded_water_carriers").get<std::uint64_t>();
    l.distinct_electron_weight = exact_from(j.at("distinct_electron_weight_sum_at_first_exclusion"));
    l.distinct_water_weight = exact_from(j.at("distinct_water_weight_sum_at_first_exclusion"));
    // Decisions are recomputed from the measures and this scope's policy.
    const auto exposure = evaluate(policy_, l.unresolved, l.applied, policy_.max_exposure_unresolved_events, "exposure");
    pc::require(exposure.json() == j.at("exposure_limit") && !exposure.refused, "coverage exposure decision differs");
    l.exposure_decision = exposure;
    l.cumulative_decision = decision_from(j.at("cumulative_limit"));
    Measures merged;for (const auto& [_, s] : l.bins) merged.merge(s);
    pc::require(exposure_json_measures_equal(merged, l.unresolved), "coverage bins do not sum to exposure totals");
    return l;
  }
  static bool exposure_json_measures_equal(const Measures& a, const Measures& b) { return a.json() == b.json(); }
  Decision decision_from(const Json& j) const {
    Decision d;d.scope = j.at("evaluated_over").get<std::string>();
    d.limit_events = policy_.max_cumulative_unresolved_events;d.limit_fraction = policy_.max_unresolved_reaction_fraction;
    d.events_upper = j.at("events_upper").is_string() ? std::numeric_limits<double>::infinity() : j.at("events_upper").get<double>();
    d.physical_upper = j.at("physical_upper").is_string() ? std::numeric_limits<double>::infinity() : j.at("physical_upper").get<double>();
    d.applied_physical = j.at("applied_physical").is_string() ? std::numeric_limits<double>::infinity() : j.at("applied_physical").get<double>();
    d.fraction_upper = j.at("fraction_upper").get<double>();d.unbounded_pairs = j.at("unbounded_sigma_pair_evaluations").get<std::uint64_t>();
    d.admit_unbounded = j.at("unbounded_admitted").get<bool>();d.absolute_exceeded = j.at("absolute_exceeded").get<bool>();
    d.fraction_exceeded = j.at("fraction_exceeded").get<bool>();d.refused = j.at("refused").get<bool>();
    d.criterion = j.at("criterion").get<std::string>();
    pc::require(d.scope == "cumulative" && !d.refused && d.json() == j, "coverage cumulative decision record");
    return d;
  }
  Json totals_json(const std::vector<ExposureLedger>& exposures) const {
    Measures u;Applied a;std::map<BinKey, Measures> bins;
    for (const auto& l : exposures) {
      u.merge(l.unresolved);a.merge(l.applied);
      for (const auto& [k, s] : l.bins) bins[k].merge(s);
    }
    return {{"unresolved", u.json()}, {"applied", a.json()}, {"bins", bins_json(bins)}};
  }
  Json ledger_json(const Ledger& l) const {
    Json last = Json::array();
    for (const auto& x : l.last_step) last.push_back(exposure_json(x));
    const auto& c = l.cumulative;
    return {{"schema", ledger_schema}, {"scope", identity_},
            {"pending_first_exposure", l.pending ? exposure_json(*l.pending) : Json(nullptr)},
            {"last_step", {{"exposures", last}, {"totals", totals_json(l.last_step)}}},
            {"cumulative", {{"accepted_steps", c.accepted_steps}, {"steps_with_exclusions", c.steps_with_exclusions},
                            {"accepted_exposures", c.accepted_exposures}, {"exposures_with_exclusions", c.exposures_with_exclusions},
                            {"unresolved", c.unresolved.json()}, {"applied", c.applied.json()}, {"bins", bins_json(c.bins)}}}};
  }
  Ledger ledger_from(const Json& j) const {
    pc::require(j.is_object() && j.at("schema") == ledger_schema && j.at("scope") == identity_, "coverage ledger belongs to another scope");
    Ledger l;
    if (!j.at("pending_first_exposure").is_null()) l.pending = exposure_from(j.at("pending_first_exposure"));
    for (const auto& x : j.at("last_step").at("exposures")) l.last_step.push_back(exposure_from(x));
    pc::require(l.last_step.size() <= 2, "coverage step has at most two exposures");
    const auto& c = j.at("cumulative");
    l.cumulative.accepted_steps = c.at("accepted_steps").get<std::uint64_t>();
    l.cumulative.steps_with_exclusions = c.at("steps_with_exclusions").get<std::uint64_t>();
    l.cumulative.accepted_exposures = c.at("accepted_exposures").get<std::uint64_t>();
    l.cumulative.exposures_with_exclusions = c.at("exposures_with_exclusions").get<std::uint64_t>();
    l.cumulative.unresolved = Measures::from(c.at("unresolved"));l.cumulative.applied = Applied::from(c.at("applied"));
    l.cumulative.bins = bins_from(c.at("bins"));
    Measures merged;for (const auto& [_, s] : l.cumulative.bins) merged.merge(s);
    pc::require(merged.json() == l.cumulative.unresolved.json(), "cumulative coverage bins do not sum to totals");
    pc::require(l.cumulative.steps_with_exclusions <= l.cumulative.accepted_steps &&
                l.cumulative.exposures_with_exclusions <= l.cumulative.accepted_exposures, "coverage counters");
    // Canonical: a restored ledger must serialize to exactly the same bytes.
    pc::require(ledger_json(l).dump() == j.dump(), "coverage ledger is not canonical");
    return l;
  }

  // Pure ledger transition for one completed exposure. Throws CoverageRefusal
  // when the hard limit is exceeded; the caller's transaction is discarded.
  Ledger advance(const Ledger& before, ExposureLedger l) const {
    Measures cu = before.cumulative.unresolved;Applied ca = before.cumulative.applied;
    if (l.slot == Slot::after_transport) {
      pc::require(before.pending && before.pending->slot == Slot::before_transport && before.pending->sequence + 1 == l.sequence,
                  "second exposure coverage requires the staged first exposure ledger");
      cu.merge(before.pending->unresolved);ca.merge(before.pending->applied);
    } else pc::require(!before.pending, "unpublished first exposure coverage ledger is pending");
    cu.merge(l.unresolved);ca.merge(l.applied);
    l.exposure_decision = evaluate(policy_, l.unresolved, l.applied, policy_.max_exposure_unresolved_events, "exposure");
    l.cumulative_decision = evaluate(policy_, cu, ca, policy_.max_cumulative_unresolved_events, "cumulative");
    for (const auto* d : {&l.exposure_decision, &l.cumulative_decision})
      if (d->refused)
        throw CoverageRefusal(std::string("parent-channel coverage limit refused the ") + d->scope + " (" + d->criterion +
                                  "; exposure sequence " + std::to_string(l.sequence) + ", slot " + name(l.slot) + ")",
                              d->criterion, identity_, exposure_json(l));
    Ledger next = before;
    auto absorb = [&](const ExposureLedger& x) {
      auto& c = next.cumulative;
      c.unresolved.merge(x.unresolved);c.applied.merge(x.applied);
      for (const auto& [k, s] : x.bins) c.bins[k].merge(s);
      pc::require(c.accepted_exposures < UINT64_MAX, "coverage exposure counter exhausted");
      ++c.accepted_exposures;if (x.unresolved.pair_evaluations) ++c.exposures_with_exclusions;
    };
    auto close_step = [&](std::vector<ExposureLedger> step) {
      bool excluded = false;
      for (const auto& x : step) { absorb(x);excluded = excluded || x.unresolved.pair_evaluations > 0; }
      pc::require(next.cumulative.accepted_steps < UINT64_MAX, "coverage step counter exhausted");
      ++next.cumulative.accepted_steps;if (excluded) ++next.cumulative.steps_with_exclusions;
      next.last_step = std::move(step);next.pending.reset();
    };
    if (l.slot == Slot::legacy) close_step({std::move(l)});
    else if (l.slot == Slot::before_transport) next.pending = std::move(l);
    else close_step({*before.pending, std::move(l)});
    return next;
  }
};

// Accumulates one exposure. A segment is one constant-population interval of
// the SSA loop; its exposure time is known only when the scheduler closes it.
class ExposureAccumulator {
  const CoverageScope& scope_;
  ExposureLedger l_;
  std::set<std::pair<std::uint64_t, std::uint64_t>> pairs_;
  std::map<std::uint64_t, double> electrons_, waters_;

public:
  struct Row {
    std::uint64_t electron_id, water_id;
    double we, wn, energy_eV, flux_m_s, volume_m3;
    Classified c;
  };
  struct Segment {
    std::vector<Row> rows;
    std::vector<std::pair<double, double>> applied; // (rate_s, event_weight)
    std::size_t pair_visits = 0;
  };
  ExposureAccumulator(const CoverageScope& scope, Slot slot, std::uint64_t sequence, double duration)
      : scope_(scope) { l_.slot = slot;l_.sequence = sequence;l_.duration_s = duration; }

  // Mirrors the sealed candidate pairing: projectiles in input order, targets
  // of the same cell in input order, and the sealed relative kinematics.
  Segment segment(const std::vector<pc::Packet>& electrons, const std::vector<std::uint64_t>& electron_ids,
                  const std::vector<pc::Packet>& water, const std::vector<std::uint64_t>& water_ids,
                  const std::unordered_map<std::uint64_t, double>& volumes,
                  const std::vector<pc::Candidate>& candidates) const {
    Segment s;
    std::unordered_map<std::uint64_t, std::vector<std::size_t>> local;
    for (std::size_t j = 0; j < water.size(); ++j) local[water[j].cell].push_back(j);
    for (std::size_t i = 0; i < electrons.size(); ++i) {
      const auto it = local.find(electrons[i].cell);
      if (it == local.end()) continue;
      const double volume = volumes.at(electrons[i].cell);
      for (const auto j : it->second) {
        ++s.pair_visits;
        const auto kin = pc::relative(electrons[i], water[j]);
        const auto c = scope_.classify(kin.target_frame_energy_eV);
        if (!unresolved(c.kind)) continue;
        s.rows.push_back({electron_ids[i], water_ids[j], electrons[i].weight, water[j].weight,
                          kin.target_frame_energy_eV, kin.moller_flux_m_s, volume, c});
      }
    }
    s.applied.reserve(candidates.size());
    for (const auto& c : candidates) s.applied.emplace_back(c.rate_s, c.event_weight);
    return s;
  }
  static double rounded(long double x) {
    const double v = static_cast<double>(x);
    pc::require(std::isfinite(v) && v >= 0, "coverage measure overflow");
    return v;
  }
  void close(const Segment& s, double begin_s, double end_s) {
    const long double dt = static_cast<long double>(end_s) - static_cast<long double>(begin_s);
    pc::require(dt >= 0 && l_.segments < UINT64_MAX, "coverage segment clock");
    ++l_.segments;
    for (const auto& r : s.rows) {
      // Each term is evaluated once in extended precision and rounded once to
      // binary64; the sums of the rounded terms are exact.
      const long double phys = static_cast<long double>(r.we) * r.wn * r.flux_m_s / r.volume_m3 * dt;
      const long double ev = static_cast<long double>(std::max(r.we, r.wn)) * r.flux_m_s / r.volume_m3 * dt;
      Measures::Terms t;t.we = r.we;t.wn = r.wn;t.energy = r.energy_eV;t.sigma = r.c.sigma;
      t.phys = rounded(phys);t.ev = rounded(ev);
      if (r.c.sigma == Sigma::tabulated) {
        t.tab_phys = t.up_phys = rounded(phys * r.c.sigma_m2);t.tab_ev = t.up_ev = rounded(ev * r.c.sigma_m2);
      } else if (r.c.sigma == Sigma::bounded) {
        t.up_phys = rounded(phys * r.c.upper_m2);t.up_ev = rounded(ev * r.c.upper_m2);
      }
      l_.unresolved.add(t);l_.bins[BinKey{r.c.bin, r.c.kind, r.c.sigma}].add(t);
      pairs_.emplace(r.electron_id, r.water_id);
      electrons_.emplace(r.electron_id, r.we);waters_.emplace(r.water_id, r.wn);
    }
    for (const auto& [rate, weight] : s.applied) {
      pc::require(l_.applied.candidate_evaluations < UINT64_MAX, "applied counter exhausted");
      ++l_.applied.candidate_evaluations;
      l_.applied.events.add(rounded(static_cast<long double>(rate) * dt));
      l_.applied.physical.add(rounded(static_cast<long double>(rate) * weight * dt));
    }
  }
  ExposureLedger finish() {
    l_.distinct_pairs = pairs_.size();l_.distinct_electrons = electrons_.size();l_.distinct_waters = waters_.size();
    for (const auto& [_, x] : electrons_) l_.distinct_electron_weight.add(x);
    for (const auto& [_, x] : waters_) l_.distinct_water_weight.add(x);
    return l_;
  }
};
} // namespace raikou::native_water::coverage_v1
