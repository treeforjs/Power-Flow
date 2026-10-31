#pragma once
#include "exact_sum.hpp"
#include "beb_sdcs_parent_kernel.hpp"
#include <algorithm>
#include <map>
#include <optional>
#include <set>
#include <tuple>
#include <unordered_map>

// Versioned coverage scope v2 for the water parent-ionization channel.
//
// Bands (target-rest electron energy E):
//   E <= 12.621 eV          known zero (v1)
//   12.621 < E < 13.5 eV    unresolved onset gap (v1, unchanged)
//   13.5 <= E <= 1000 eV    dataset: Song 2021 Table 12 + parent_ground_equal_pair_CM_v1 (v1, unchanged)
//   1000 < E <= 10000 eV    dataset: outer3_BEB_Tan_parent_v1 rate + BEB SDCS kernel
//   10000 < E <= E_max      PROVISIONAL model-derived: RBEB x Tan parent weights + RBEB SDCS kernel
//   E > E_max               unresolved, no bound (refused unless admit_unbounded)
// Without HighEnergyBands the scope is the v1 coverage map exactly (every
// energy above 1 keV unresolved) and binds the v1 channel registration.
//
// All three supported bands share one registered channel and one stochastic
// owner (the sealed native arena admits one channel binding per owner and one
// owner per full-time exposure, so a single SSA over the combined hazard is
// kept). The banded channel binding is a composite identity; each event's band
// is recorded in its id prefix, its receipt and the exact per-category ledger.
namespace raikou::native_water::coverage_v2 {
using ExactSum = ::raikou::native_water::coverage_v1::ExactSum;

inline constexpr const char* scope_schema = "parent-channel-coverage-scope-v2";
inline constexpr const char* ledger_schema = "parent-channel-coverage-ledger-v2";
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
  pc::require(j.is_object() && j.size() == 3 && j.at("mantissa_hex").is_string() && j.at("exponent2").is_number_integer(),
              "malformed exact ledger value");
  auto s = ExactSum::from_canonical(j.at("mantissa_hex").get<std::string>(), j.at("exponent2").get<long>());
  const double v = s.value();
  pc::require(std::isfinite(v) ? (j.at("value").is_number() && j.at("value").get<double>() == v) : j.at("value").is_null(),
              "exact ledger value differs from its rounded record");
  return s;
}

struct CoveragePolicy {
  double max_exposure_unresolved_events = 0.01;
  double max_cumulative_unresolved_events = 1.0;
  double max_unresolved_reaction_fraction = 0.01;
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

enum class Kind { known_zero, supported, onset_gap, below_table, above_table, kernel_domain, beb_supported, provisional_supported };
enum class Sigma { none, tabulated, bounded, unbounded, dataset_formula, provisional_model };
inline const char* name(Kind k) {
  switch (k) {
  case Kind::known_zero: return "known_zero";
  case Kind::supported: return "supported";
  case Kind::onset_gap: return "onset_gap";
  case Kind::below_table: return "below_table_without_threshold";
  case Kind::above_table: return "above_table";
  case Kind::kernel_domain: return "tabulated_outside_kernel_domain";
  case Kind::beb_supported: return "supported_dataset_outer3_BEB";
  case Kind::provisional_supported: return "supported_provisional_RBEB_model";
  }
  return "";
}
inline const char* name(Sigma s) {
  switch (s) {
  case Sigma::none: return "none";
  case Sigma::tabulated: return "tabulated";
  case Sigma::bounded: return "unknown_upper_bounded";
  case Sigma::unbounded: return "unknown_unbounded";
  case Sigma::dataset_formula: return "dataset_formula";
  case Sigma::provisional_model: return "provisional_model";
  }
  return "";
}
inline bool unresolved(Kind k) {
  return k == Kind::onset_gap || k == Kind::below_table || k == Kind::above_table || k == Kind::kernel_domain;
}
// Supported categories of executed/applied reactions.
enum class Category { song_table = 0, outer3_beb = 1, provisional_rbeb = 2 };
inline constexpr std::size_t category_count = 3;
inline const char* name(Category c) {
  switch (c) {
  case Category::song_table: return "dataset_song2021_parent_table";
  case Category::outer3_beb: return "dataset_outer3_BEB_Tan_parent";
  case Category::provisional_rbeb: return "provisional_RBEB_Tan_parent_model";
  }
  return "";
}
inline const char* basis(Category c) { return c == Category::provisional_rbeb ? "model_derived" : "dataset"; }
inline Category category_of(Kind k) {
  if (k == Kind::supported) return Category::song_table;
  if (k == Kind::beb_supported) return Category::outer3_beb;
  pc::require(k == Kind::provisional_supported, "category of an unsupported coverage class");
  return Category::provisional_rbeb;
}
inline const char* event_prefix(Category c) {
  switch (c) {
  case Category::song_table: return "conditional/native-parent/";
  case Category::outer3_beb: return "conditional/native-parent-outer3-beb/";
  case Category::provisional_rbeb: return "model-derived/native-parent-rbeb/";
  }
  return "";
}
struct Classified {
  Kind kind = Kind::known_zero;
  Sigma sigma = Sigma::none;
  double sigma_m2 = 0, upper_m2 = 0;
  std::size_t bin = 0;
};
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

// Unresolved measures (v1 definitions, unchanged).
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
    else { energy_min_eV = std::min(energy_min_eV, t.energy);energy_max_eV = std::max(energy_max_eV, t.energy); }
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
    if (!pair_evaluations) { energy_min_eV = o.energy_min_eV;energy_max_eV = o.energy_max_eV; }
    else { energy_min_eV = std::min(energy_min_eV, o.energy_min_eV);energy_max_eV = std::max(energy_max_eV, o.energy_max_eV); }
    pair_evaluations += o.pair_evaluations;tabulated += o.tabulated;bounded += o.bounded;unbounded += o.unbounded;
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
    if (pair_evaluations) { j["energy_min_eV"] = energy_min_eV;j["energy_max_eV"] = energy_max_eV; }
    else { j["energy_min_eV"] = nullptr;j["energy_max_eV"] = nullptr; }
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
    s.electron_weight = exact_from(j.at("electron_weight_sum"));s.water_weight = exact_from(j.at("water_weight_sum"));
    s.event_weight = exact_from(j.at("event_weight_sum"));
    s.physical_exposure = exact_from(j.at("physical_exposure_inv_m2"));s.event_exposure = exact_from(j.at("event_exposure_inv_m2"));
    s.tabulated_physical = exact_from(j.at("tabulated_expected_physical_reactions"));
    s.tabulated_events = exact_from(j.at("tabulated_expected_events"));
    s.upper_physical = exact_from(j.at("upper_bound_physical_reactions"));s.upper_events = exact_from(j.at("upper_bound_events"));
    s.unbounded_physical_exposure = exact_from(j.at("unbounded_physical_exposure_inv_m2"));
    s.unbounded_event_exposure = exact_from(j.at("unbounded_event_exposure_inv_m2"));
    if (s.pair_evaluations) { s.energy_min_eV = j.at("energy_min_eV").get<double>();s.energy_max_eV = j.at("energy_max_eV").get<double>(); }
    return s;
  }
};
// Applied (qualified) hazard totals, the limit denominator (v1 definition).
struct Applied {
  std::uint64_t candidate_evaluations = 0;
  ExactSum physical, events;
  void merge(const Applied& o) {
    pc::require(candidate_evaluations <= UINT64_MAX - o.candidate_evaluations, "applied counter exhausted");
    candidate_evaluations += o.candidate_evaluations;physical.add(o.physical);events.add(o.events);
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
// Supported (applied and executed) measures of one category. With the
// sigma-free exposures a replacement cross section sigma'(E) gives the
// reprocessed expectation per energy bin; the executed events are counted
// and weighted exactly and listed individually in the exposure receipts.
struct Supported {
  Applied applied;
  ExactSum sigma_free_physical, sigma_free_events;
  std::uint64_t events = 0;
  ExactSum executed_event_weight;
  double energy_min_eV = 0, energy_max_eV = 0; // over applied candidate rows
  struct Terms { double phys_sigma = 0, ev_sigma = 0, phys = 0, ev = 0, energy = 0; };
  void add(const Terms& t) {
    pc::require(applied.candidate_evaluations < UINT64_MAX, "applied counter exhausted");
    if (!applied.candidate_evaluations) energy_min_eV = energy_max_eV = t.energy;
    else { energy_min_eV = std::min(energy_min_eV, t.energy);energy_max_eV = std::max(energy_max_eV, t.energy); }
    ++applied.candidate_evaluations;applied.physical.add(t.phys_sigma);applied.events.add(t.ev_sigma);
    sigma_free_physical.add(t.phys);sigma_free_events.add(t.ev);
  }
  void executed(double weight) {
    pc::require(events < UINT64_MAX, "executed event counter exhausted");
    ++events;executed_event_weight.add(weight);
  }
  bool empty() const { return !applied.candidate_evaluations && !events; }
  void merge(const Supported& o) {
    if (o.applied.candidate_evaluations) {
      if (!applied.candidate_evaluations) { energy_min_eV = o.energy_min_eV;energy_max_eV = o.energy_max_eV; }
      else { energy_min_eV = std::min(energy_min_eV, o.energy_min_eV);energy_max_eV = std::max(energy_max_eV, o.energy_max_eV); }
    }
    applied.merge(o.applied);sigma_free_physical.add(o.sigma_free_physical);sigma_free_events.add(o.sigma_free_events);
    pc::require(events <= UINT64_MAX - o.events, "executed event counter exhausted");
    events += o.events;executed_event_weight.add(o.executed_event_weight);
  }
  Json json() const {
    Json j{{"applied", applied.json()}, {"sigma_free_physical_exposure_inv_m2", exact_json(sigma_free_physical)},
           {"sigma_free_event_exposure_inv_m2", exact_json(sigma_free_events)}, {"executed_events", events},
           {"executed_event_weight_sum", exact_json(executed_event_weight)}};
    if (applied.candidate_evaluations) { j["energy_min_eV"] = energy_min_eV;j["energy_max_eV"] = energy_max_eV; }
    else { j["energy_min_eV"] = nullptr;j["energy_max_eV"] = nullptr; }
    return j;
  }
  static Supported from(const Json& j) {
    Supported s;s.applied = Applied::from(j.at("applied"));
    s.sigma_free_physical = exact_from(j.at("sigma_free_physical_exposure_inv_m2"));
    s.sigma_free_events = exact_from(j.at("sigma_free_event_exposure_inv_m2"));
    s.events = j.at("executed_events").get<std::uint64_t>();
    s.executed_event_weight = exact_from(j.at("executed_event_weight_sum"));
    if (s.applied.candidate_evaluations) { s.energy_min_eV = j.at("energy_min_eV").get<double>();s.energy_max_eV = j.at("energy_max_eV").get<double>(); }
    else pc::require(j.at("energy_min_eV").is_null() && j.at("energy_max_eV").is_null(), "empty category energy range");
    return s;
  }
};
struct CategoryLedger {
  Supported total;
  std::map<std::size_t, Supported> bins;
  void merge(const CategoryLedger& o) { total.merge(o.total);for (const auto& [k, s] : o.bins) bins[k].merge(s); }
};
using Categories = std::array<CategoryLedger, category_count>;

struct Decision {
  std::string scope;
  double limit_events = 0, limit_fraction = 0, events_upper = 0, physical_upper = 0, applied_physical = 0, fraction_upper = 0;
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
inline Decision evaluate(const CoveragePolicy& p, const Measures& u, const Applied& a, double event_limit, const char* over) {
  Decision d;d.scope = over;d.limit_events = event_limit;d.limit_fraction = p.max_unresolved_reaction_fraction;
  d.admit_unbounded = p.admit_unbounded;d.unbounded_pairs = u.unbounded;
  d.events_upper = u.upper_events.value();d.physical_upper = u.upper_physical.value();
  d.applied_physical = a.physical.value();
  const long double U = d.physical_upper, A = d.applied_physical, f = d.limit_fraction;
  d.fraction_upper = U == 0 ? 0. : (std::isinf(U) ? 1. : static_cast<double>(U / (A + U)));
  d.absolute_exceeded = static_cast<long double>(d.events_upper) > static_cast<long double>(event_limit);
  d.fraction_exceeded = U > 0 && (std::isinf(U) || U * (1 - f) > f * A);
  if (u.unbounded && !p.admit_unbounded) { d.refused = true;d.criterion = "unbounded_unresolved_sigma"; }
  else if (d.absolute_exceeded && d.fraction_exceeded) { d.refused = true;d.criterion = "events_and_fraction_upper_bounds"; }
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
  Categories categories;
  std::map<BinKey, Measures> bins;
  std::uint64_t distinct_pairs = 0, distinct_electrons = 0, distinct_waters = 0;
  ExactSum distinct_electron_weight, distinct_water_weight;
  Decision exposure_decision, cumulative_decision;
};
struct Cumulative {
  std::uint64_t accepted_steps = 0, steps_with_exclusions = 0, accepted_exposures = 0, exposures_with_exclusions = 0;
  Measures unresolved;
  Applied applied;
  Categories categories;
  std::map<BinKey, Measures> bins;
};
struct Ledger {
  std::optional<ExposureLedger> pending;
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
  std::shared_ptr<const HighEnergyBands> he_;
  std::string dataset_sha_, kernel_sha_, bound_id_, identity_;
  ArenaBinding arena_;
  CoveragePolicy policy_;
  bool zero_present_ = false;
  double zero_below_ = 0, table_lo_ = 0, table_hi_ = 0, kernel_lo_ = 0, kernel_hi_ = 0, bound_lo_ = 0, bound_hi_ = 0, bound_scale_ = 0;
  double beb_hi_ = 0, provisional_hi_ = 0;
  bool provisional_ = false;
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
  // High-energy band of E > table_hi_, if any.
  std::optional<Kind> band_kind(double E) const {
    if (!he_ || !(E > table_hi_)) return std::nullopt;
    if (E <= beb_hi_) return Kind::beb_supported;
    if (provisional_ && E <= provisional_hi_) return Kind::provisional_supported;
    return std::nullopt;
  }
  Classified classify_unbinned(double E) const {
    Classified c;
    pc::require(std::isfinite(E) && E >= 0, "coverage classification needs a finite nonnegative energy");
    auto unknown = [&](Kind k) {
      c.kind = k;
      if (E >= bound_lo_ && E <= bound_hi_) { c.sigma = Sigma::bounded;c.upper_m2 = envelope(E); }
      else c.sigma = Sigma::unbounded;
    };
    if (E > table_hi_) {
      if (const auto band = band_kind(E)) {
        c.kind = *band;
        if (*band == Kind::beb_supported) { c.sigma = Sigma::dataset_formula;c.sigma_m2 = he_->outer3_parent_m2(E); }
        else { c.sigma = Sigma::provisional_model;c.sigma_m2 = he_->provisional_parent_m2(E); }
        pc::require(c.sigma_m2 > 0, "high-energy band with a nonpositive parent rate");
      } else unknown(Kind::above_table);
    } else if (zero_present_ && E <= zero_below_) c.kind = Kind::known_zero;
    else if (E < table_lo_) unknown(zero_present_ ? Kind::onset_gap : Kind::below_table);
    else {
      const double s = data_->observable_value_SI(w::parent_channel, E);
      if (!(s > 0)) c.kind = Kind::known_zero;
      else if (E >= kernel_lo_ && E <= kernel_hi_) { c.kind = Kind::supported;c.sigma = Sigma::tabulated;c.sigma_m2 = s; }
      else { c.kind = Kind::kernel_domain;c.sigma = Sigma::tabulated;c.sigma_m2 = s;c.upper_m2 = s; }
    }
    return c;
  }
  void construct(const w::ParentBridge& bridge, const std::string& bound_observable) {
    bound_id_ = bound_observable;
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
    for (const auto& x : p.at("energy_eV")) {
      const double E = num(x);
      if (E >= bound_lo_ && E <= bound_hi_)
        pc::require(data_->observable_value_SI(w::parent_channel, E) <= envelope(E), "upper-bound observable does not dominate the parent table");
    }
    std::set<double> e{0., 12.75, 13., 13.25, 1500., 2000., 3000., 5000., 8000., 20000., 50000., 1e5, 1e6,
                       table_lo_, table_hi_, kernel_lo_, kernel_hi_, bound_lo_, bound_hi_};
    if (zero_present_) e.insert(zero_below_);
    if (he_) {
      // The high-energy bands start exactly where the Song table and the v1
      // kernel end, and the Song table of the outer3 freeze is the bound one.
      pc::require(he_->song_dataset_sha256() == dataset_sha_ && he_->song_kernel_sha256() == kernel_sha_ &&
                      he_->beb_low_eV() == table_hi_ && kernel_hi_ == table_hi_,
                  "high-energy bands bound to another Song binding or not adjacent to the Song table");
      beb_hi_ = he_->beb_high_eV();provisional_ = he_->provisional().admit;provisional_hi_ = he_->provisional().maximum_eV;
      for (double x : {2e5, 5e5, 2e6, 5e6, beb_hi_}) e.insert(x);
      if (provisional_) e.insert(provisional_hi_);
      arena_ = he_->arena_binding();
    } else
      arena_ = {dataset_sha_, kernel_sha_, "", w::parent_model}; // binding file supplied by the subflow
    edges_.assign(e.begin(), e.end());
    edges_.push_back(std::numeric_limits<double>::infinity());
    for (std::size_t i = 0; i + 1 < edges_.size(); ++i) {
      const double lo = edges_[i], hi = edges_[i + 1];
      const double mid = std::isinf(hi) ? 2 * lo : (i == 0 ? hi / 2 : lo + (hi - lo) / 2);
      auto c = classify_unbinned(mid);c.bin = i;c.sigma_m2 = c.upper_m2 = 0;bin_class_.push_back(c);
    }
    std::set<double> probes{0., 1., 10., 12.7, 13.4999, 37.9, 38.3, 39.6, 100., 999.999, 1000.0001, 8000., 20000.};
    if (he_) for (double x : {1000.5, 9999., 10000.5, 1e5, 1e6, 1.5e7}) probes.insert(x);
    for (std::size_t i = 0; i + 1 < edges_.size(); ++i) {
      const double lo = edges_[i], hi = edges_[i + 1];
      probes.insert(lo);
      if (lo > 0) probes.insert(std::nextafter(lo, 0.));
      probes.insert(std::nextafter(lo, std::numeric_limits<double>::infinity()));
      if (!std::isinf(hi)) probes.insert(lo + (hi - lo) / 2);
    }
    probes_.assign(probes.begin(), probes.end());
    require_bridge(bridge);
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
      {"exclusion", "pair excluded from the exposure and recorded; no renormalization, no rate invention"},
      {"high_energy_bands", he_ ? he_->description() : Json(nullptr)},
      {"arena_channel_binding", he_ ? arena_.json() : Json("v1 registration (dataset, kernel, binding file, parent_ground_equal_pair_CM_v1)")},
      {"categories", {{name(Category::song_table), {{"basis", "dataset"}, {"band_eV", {table_lo_, table_hi_}}}},
                      {name(Category::outer3_beb), {{"basis", "dataset"}, {"admitted", bool(he_)},
                                                    {"band_eV", he_ ? Json{{"lower_exclusive", table_hi_}, {"upper_inclusive", beb_hi_}} : Json(nullptr)}}},
                      {name(Category::provisional_rbeb), {{"basis", "model_derived"}, {"admitted", bool(he_) && provisional_},
                                                          {"band_eV", he_ && provisional_ ? Json{{"lower_exclusive", beb_hi_}, {"upper_inclusive", provisional_hi_}} : Json(nullptr)}}}}},
      {"seams", seams()}};
    identity_ = sha(description_.dump());
  }

public:
  // v1-equivalent scope: the v1 coverage map, every energy above 1 keV unresolved.
  CoverageScope(std::shared_ptr<const w::Data> data, std::string dataset_sha256, const w::ParentBridge& bridge,
                CoveragePolicy policy = {}, std::string bound_observable = default_upper_bound_observable)
      : data_(std::move(data)), dataset_sha_(std::move(dataset_sha256)), policy_(policy) {
    construct(bridge, bound_observable);
  }
  // Coverage v2 with the dataset BEB band and (optionally) the provisional band.
  CoverageScope(std::shared_ptr<const w::Data> data, std::string dataset_sha256, const w::ParentBridge& bridge,
                std::shared_ptr<const HighEnergyBands> bands, CoveragePolicy policy = {},
                std::string bound_observable = default_upper_bound_observable)
      : data_(std::move(data)), he_(std::move(bands)), dataset_sha_(std::move(dataset_sha256)), policy_(policy) {
    pc::require(bool(he_), "coverage v2 high-energy bands missing");
    construct(bridge, bound_observable);
  }
  const std::string& identity() const noexcept { return identity_; }
  const Json& description() const noexcept { return description_; }
  const CoveragePolicy& policy() const noexcept { return policy_; }
  const std::string& dataset_sha256() const noexcept { return dataset_sha_; }
  const std::string& kernel_sha256() const noexcept { return kernel_sha_; }
  bool high_energy() const noexcept { return bool(he_); }
  const HighEnergyBands* bands() const noexcept { return he_.get(); }
  // Registered native channel binding. Without bands it is the v1
  // registration; the binding file digest then comes from the subflow binding.
  ArenaBinding arena_binding(const std::string& v1_binding_file) const {
    if (he_) {
      pc::require(v1_binding_file == he_->song_kernel_binding_file_sha256(), "Song-band binding file differs from the bound one");
      return arena_;
    }
    auto a = arena_;a.kernel_binding_file_sha256 = v1_binding_file;return a;
  }
  std::size_t bin_count() const noexcept { return bin_class_.size(); }
  const std::vector<double>& edges() const noexcept { return edges_; }
  std::size_t bin_of(double E) const {
    auto i = static_cast<std::size_t>(std::lower_bound(edges_.begin(), edges_.end(), E) - edges_.begin());
    return i == 0 ? 0 : i - 1;
  }
  Classified classify(double E) const {
    auto c = classify_unbinned(E);c.bin = bin_of(E);
    pc::require(c.bin + 1 < edges_.size() && (c.bin == 0 ? E >= edges_[0] : E > edges_[c.bin]) && E <= edges_[c.bin + 1], "coverage bin lookup");
    return c;
  }
  // Supported category of a candidate energy (positive rate only).
  Category category(double E) const {
    if (E > table_hi_) {
      const auto band = band_kind(E);
      pc::require(bool(band), "candidate above the supported bands");
      return category_of(*band);
    }
    return Category::song_table;
  }
  bool excluded(double E) const {
    if (!(std::isfinite(E) && E >= 0)) return false;
    if (E > table_hi_) return !band_kind(E);
    if (zero_present_ && E <= zero_below_) return false;
    if (E < table_lo_) return true;
    if (E >= kernel_lo_ && E <= kernel_hi_) return false;
    return data_->observable_value_SI(w::parent_channel, E) > 0;
  }
  // Composite channel: the sealed Song bridge for E <= 1 keV, the high-energy
  // band rates above; unresolved energies are counted and excluded.
  pc::Channel channel(const w::ParentBridge& bridge, std::uint64_t* excluded_count, const std::string& v1_binding_file) const {
    auto ch = bridge.channel();
    const auto arena = arena_binding(v1_binding_file);
    ch.data_sha256 = arena.dataset_sha256;ch.kernel_sha256 = arena.kernel_manifest_sha256;
    auto inner = ch.sigma_m2;
    const CoverageScope* scope = this;
    ch.sigma_m2 = [inner, scope, excluded_count](double E) {
      if (scope->excluded(E)) { if (excluded_count) ++*excluded_count;return 0.; }
      if (scope->he_ && E > scope->table_hi_) {
        return *scope->band_kind(E) == Kind::beb_supported ? scope->he_->outer3_parent_m2(E) : scope->he_->provisional_parent_m2(E);
      }
      return inner(E);
    };
    return ch;
  }
  void require_bridge(const w::ParentBridge& bridge) const {
    const auto channel = bridge.channel();
    pc::require(channel.data_sha256 == dataset_sha_ && channel.kernel_sha256 == kernel_sha_, "coverage scope bound to another bridge");
    for (double E : probes_) {
      const auto c = classify_unbinned(E);bool coverage_error = false;double s = -1;
      try { s = channel.sigma_m2(E); } catch (const ::raikou::kinetics::CollisionCoverageError&) { coverage_error = true; }
      const bool band = c.kind == Kind::beb_supported || c.kind == Kind::provisional_supported;
      pc::require(c.kind == Kind::known_zero ? (!coverage_error && s == 0)
                  : c.kind == Kind::supported ? (!coverage_error && s == c.sigma_m2)
                  : band ? (coverage_error && c.sigma_m2 > 0) : coverage_error,
                  "coverage scope classification differs from the subflow bridge");
      pc::require(excluded(E) == unresolved(c.kind), "coverage exclusion predicate differs from the classification");
    }
  }
  Json seams() const {
    if (!he_) return nullptr;
    const double above1 = std::nextafter(table_hi_, std::numeric_limits<double>::infinity());
    const double song = data_->observable_value_SI(w::parent_channel, table_hi_);
    const double o_at = he_->outer3_parent_m2(table_hi_), o_above = he_->outer3_parent_m2(above1);
    Json s{{"song_to_outer3_BEB", {{"energy_eV", table_hi_}, {"song_parent_m2_at_edge", song},
                                   {"outer3_parent_m2_at_edge", o_at}, {"outer3_parent_m2_next_binary64_above", o_above},
                                   {"relative_jump", (o_above - song) / song}, {"handling", "recorded, not smoothed"}}}};
    const double o10 = he_->outer3_parent_m2(beb_hi_);
    const double above10 = std::nextafter(beb_hi_, std::numeric_limits<double>::infinity());
    const double r10 = static_cast<double>(he_->rbeb_parent_any(above10)), r10_at = static_cast<double>(he_->rbeb_parent_any(beb_hi_));
    s["outer3_BEB_to_provisional_RBEB"] = {{"energy_eV", beb_hi_}, {"outer3_parent_m2_at_edge", o10},
                                           {"provisional_parent_m2_next_binary64_above", r10},
                                           {"provisional_parent_m2_at_edge_diagnostic", r10_at},
                                           {"relative_jump", (r10 - o10) / o10}, {"admitted", provisional_},
                                           {"handling", "recorded, not smoothed"}};
    s["provisional_upper_end"] = {{"energy_eV", provisional_hi_}, {"admitted", provisional_},
                                  {"above", "unresolved, no data bound"}};
    return s;
  }

  // ---- ledger serialization ----
  Json bins_json(const std::map<BinKey, Measures>& bins) const {
    Json out = Json::array();
    for (const auto& [key, s] : bins) {
      const auto& [i, kind, sigma] = key;
      pc::require(i + 1 < edges_.size() && unresolved(kind) && (sigma == Sigma::tabulated || sigma == Sigma::bounded || sigma == Sigma::unbounded) &&
                      s.pair_evaluations > 0, "coverage bin record");
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
      const BinKey key{b.at("bin").get<std::size_t>(), kind_named(b.at("coverage").get<std::string>()), sigma_named(b.at("sigma").get<std::string>())};
      pc::require(std::get<0>(key) + 1 < edges_.size() && !out.count(key), "coverage bin index");
      out.emplace(key, Measures::from(b));
    }
    return out;
  }
  Json categories_json(const Categories& cats) const {
    Json out = Json::array();
    for (std::size_t c = 0; c < category_count; ++c) {
      const auto cat = static_cast<Category>(c);
      Json bins = Json::array();
      for (const auto& [i, s] : cats[c].bins) {
        pc::require(i + 1 < edges_.size() && !s.empty(), "category bin record");
        auto j = s.json();j["bin"] = i;j["lower_eV_exclusive"] = edges_[i];j["upper_eV_inclusive"] = number_or_inf(edges_[i + 1]);
        bins.push_back(std::move(j));
      }
      auto j = cats[c].total.json();
      j["category"] = name(cat);j["basis"] = basis(cat);j["model_derived"] = cat == Category::provisional_rbeb;
      j["bins"] = std::move(bins);
      out.push_back(std::move(j));
    }
    return out;
  }
  Categories categories_from(const Json& j) const {
    pc::require(j.is_array() && j.size() == category_count, "coverage categories");
    Categories out;
    for (std::size_t c = 0; c < category_count; ++c) {
      const auto cat = static_cast<Category>(c);const auto& r = j[c];
      pc::require(r.at("category") == name(cat) && r.at("basis") == basis(cat) && r.at("model_derived") == (cat == Category::provisional_rbeb),
                  "coverage category identity");
      out[c].total = Supported::from(r);
      pc::require(r.at("bins").is_array() && r.at("bins").size() <= bin_class_.size(), "category bins");
      Supported merged;
      for (const auto& b : r.at("bins")) {
        const auto i = b.at("bin").get<std::size_t>();
        pc::require(i + 1 < edges_.size() && !out[c].bins.count(i), "category bin index");
        auto s = Supported::from(b);pc::require(!s.empty(), "empty category bin");
        merged.merge(s);out[c].bins.emplace(i, s);
      }
      pc::require(merged.json() == out[c].total.json(), "category bins do not sum to the category total");
      if (!he_) pc::require(cat == Category::song_table || out[c].total.empty(), "high-energy category without high-energy bands");
      if (cat == Category::provisional_rbeb && !provisional_) pc::require(out[c].total.empty(), "provisional category not admitted");
    }
    return out;
  }
  static Applied applied_of(const Categories& c) {
    Applied a;for (const auto& x : c) a.merge(x.total.applied);return a;
  }
  Json exposure_json(const ExposureLedger& l) const {
    return {{"slot", name(l.slot)}, {"sequence", l.sequence}, {"segments", l.segments}, {"duration_s", l.duration_s},
            {"unresolved", l.unresolved.json()}, {"applied", l.applied.json()}, {"bins", bins_json(l.bins)},
            {"categories", categories_json(l.categories)},
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
    l.categories = categories_from(j.at("categories"));
    pc::require(applied_of(l.categories).json() == l.applied.json(), "coverage categories do not sum to the applied total");
    l.distinct_pairs = j.at("distinct_excluded_pairs").get<std::uint64_t>();
    l.distinct_electrons = j.at("distinct_excluded_electron_carriers").get<std::uint64_t>();
    l.distinct_waters = j.at("distinct_excluded_water_carriers").get<std::uint64_t>();
    l.distinct_electron_weight = exact_from(j.at("distinct_electron_weight_sum_at_first_exclusion"));
    l.distinct_water_weight = exact_from(j.at("distinct_water_weight_sum_at_first_exclusion"));
    const auto exposure = evaluate(policy_, l.unresolved, l.applied, policy_.max_exposure_unresolved_events, "exposure");
    pc::require(exposure.json() == j.at("exposure_limit") && !exposure.refused, "coverage exposure decision differs");
    l.exposure_decision = exposure;
    l.cumulative_decision = decision_from(j.at("cumulative_limit"));
    Measures merged;for (const auto& [_, s] : l.bins) merged.merge(s);
    pc::require(merged.json() == l.unresolved.json(), "coverage bins do not sum to exposure totals");
    return l;
  }
  Decision decision_from(const Json& j) const {
    Decision d;d.scope = j.at("evaluated_over").get<std::string>();
    d.limit_events = policy_.max_cumulative_unresolved_events;d.limit_fraction = policy_.max_unresolved_reaction_fraction;
    auto value = [](const Json& x) { return x.is_string() ? std::numeric_limits<double>::infinity() : x.get<double>(); };
    d.events_upper = value(j.at("events_upper"));d.physical_upper = value(j.at("physical_upper"));d.applied_physical = value(j.at("applied_physical"));
    d.fraction_upper = j.at("fraction_upper").get<double>();d.unbounded_pairs = j.at("unbounded_sigma_pair_evaluations").get<std::uint64_t>();
    d.admit_unbounded = j.at("unbounded_admitted").get<bool>();d.absolute_exceeded = j.at("absolute_exceeded").get<bool>();
    d.fraction_exceeded = j.at("fraction_exceeded").get<bool>();d.refused = j.at("refused").get<bool>();
    d.criterion = j.at("criterion").get<std::string>();
    pc::require(d.scope == "cumulative" && !d.refused && d.json() == j, "coverage cumulative decision record");
    return d;
  }
  Json totals_json(const std::vector<ExposureLedger>& exposures) const {
    Measures u;Applied a;std::map<BinKey, Measures> bins;Categories cats;
    for (const auto& l : exposures) {
      u.merge(l.unresolved);a.merge(l.applied);
      for (const auto& [k, s] : l.bins) bins[k].merge(s);
      for (std::size_t c = 0; c < category_count; ++c) cats[c].merge(l.categories[c]);
    }
    return {{"unresolved", u.json()}, {"applied", a.json()}, {"bins", bins_json(bins)}, {"categories", categories_json(cats)}};
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
                            {"unresolved", c.unresolved.json()}, {"applied", c.applied.json()}, {"bins", bins_json(c.bins)},
                            {"categories", categories_json(c.categories)}}}};
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
    l.cumulative.bins = bins_from(c.at("bins"));l.cumulative.categories = categories_from(c.at("categories"));
    pc::require(applied_of(l.cumulative.categories).json() == l.cumulative.applied.json(), "cumulative categories do not sum to the applied total");
    Measures merged;for (const auto& [_, s] : l.cumulative.bins) merged.merge(s);
    pc::require(merged.json() == l.cumulative.unresolved.json(), "cumulative coverage bins do not sum to totals");
    pc::require(l.cumulative.steps_with_exclusions <= l.cumulative.accepted_steps &&
                l.cumulative.exposures_with_exclusions <= l.cumulative.accepted_exposures, "coverage counters");
    pc::require(ledger_json(l).dump() == j.dump(), "coverage ledger is not canonical");
    return l;
  }

  // Pure ledger transition for one completed exposure.
  Ledger advance(const Ledger& before, ExposureLedger l) const {
    Measures cu = before.cumulative.unresolved;Applied ca = before.cumulative.applied;
    if (l.slot == Slot::after_transport) {
      pc::require(before.pending && before.pending->slot == Slot::before_transport && before.pending->sequence + 1 == l.sequence,
                  "second exposure coverage requires the staged first exposure ledger");
      cu.merge(before.pending->unresolved);ca.merge(before.pending->applied);
    } else pc::require(!before.pending, "unpublished first exposure coverage ledger is pending");
    pc::require(applied_of(l.categories).json() == l.applied.json(), "exposure categories do not sum to the applied total");
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
      for (std::size_t i = 0; i < category_count; ++i) c.categories[i].merge(x.categories[i]);
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

// Accumulates one exposure (v1 accumulator plus per-category applied and
// executed measures).
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
  struct AppliedRow {
    double rate_s, event_weight, we, wn, flux_m_s, volume_m3, energy_eV;
    Category category;
    std::size_t bin;
  };
  struct Segment {
    std::vector<Row> rows;
    std::vector<AppliedRow> applied;
    std::size_t pair_visits = 0;
  };
  ExposureAccumulator(const CoverageScope& scope, Slot slot, std::uint64_t sequence, double duration) : scope_(scope) {
    l_.slot = slot;l_.sequence = sequence;l_.duration_s = duration;
  }
  Segment segment(const std::vector<pc::Packet>& electrons, const std::vector<std::uint64_t>& electron_ids,
                  const std::vector<pc::Packet>& water, const std::vector<std::uint64_t>& water_ids,
                  const std::unordered_map<std::uint64_t, double>& volumes, const std::vector<pc::Candidate>& candidates) const {
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
        s.rows.push_back({electron_ids[i], water_ids[j], electrons[i].weight, water[j].weight, kin.target_frame_energy_eV,
                          kin.moller_flux_m_s, volume, c});
      }
    }
    s.applied.reserve(candidates.size());
    for (const auto& c : candidates)
      s.applied.push_back({c.rate_s, c.event_weight, c.projectile_before, c.target_before, c.flux_m_s, c.volume_m3, c.energy_eV,
                           scope_.category(c.energy_eV), scope_.bin_of(c.energy_eV)});
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
    for (const auto& a : s.applied) {
      pc::require(l_.applied.candidate_evaluations < UINT64_MAX, "applied counter exhausted");
      // v1 applied terms (unchanged expressions).
      const double events = rounded(static_cast<long double>(a.rate_s) * dt);
      const double physical = rounded(static_cast<long double>(a.rate_s) * a.event_weight * dt);
      ++l_.applied.candidate_evaluations;l_.applied.events.add(events);l_.applied.physical.add(physical);
      Supported::Terms t;t.ev_sigma = events;t.phys_sigma = physical;t.energy = a.energy_eV;
      t.phys = rounded(static_cast<long double>(a.we) * a.wn * a.flux_m_s / a.volume_m3 * dt);
      t.ev = rounded(static_cast<long double>(std::max(a.we, a.wn)) * a.flux_m_s / a.volume_m3 * dt);
      auto& cat = l_.categories[static_cast<std::size_t>(a.category)];
      cat.total.add(t);cat.bins[a.bin].add(t);
    }
  }
  // One executed event of a supported category.
  void executed(double energy_eV, double event_weight) {
    const auto c = scope_.category(energy_eV);
    auto& cat = l_.categories[static_cast<std::size_t>(c)];
    cat.total.executed(event_weight);cat.bins[scope_.bin_of(energy_eV)].executed(event_weight);
  }
  ExposureLedger finish() {
    l_.distinct_pairs = pairs_.size();l_.distinct_electrons = electrons_.size();l_.distinct_waters = waters_.size();
    for (const auto& [_, x] : electrons_) l_.distinct_electron_weight.add(x);
    for (const auto& [_, x] : waters_) l_.distinct_water_weight.add(x);
    return l_;
  }
};
} // namespace raikou::native_water::coverage_v2
