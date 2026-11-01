// Arrival / count / RNG / ID owner v1.1 (versioned from pic-arrival-owner-v1-20261030,
// unchanged and read in place for canonical.hpp and philox_core.h). Changes: the
// neutral stream rule (quantized release, see neutral_release_rule and API-v1.1.txt),
// the policy/configuration/checkpoint identities and the sampling identity prefix.
// arrival_owner.hpp -- arrival / count / RNG / ID owner v1 for source births.
// Normative contract: API-v1.txt. Summary:
//  * Policy (pf-arrival-policy-v1): per stream a target macroparticle weight w0.
//    Electrons: n = SR(N / w0) births of weight exactly w0 (debit = n*w0 exactly).
//    Neutrals:  the exact source debit D (one binary64 value) is split into
//               n = min(max(1, SR(D / w0)), M) leaves whose weights are integer
//               multiples of D's own quantum 2^E (D = M * 2^E) and sum to D exactly.
//    SR = stochastic rounding with the exact Bernoulli of philox_core.h, keyed by
//    (run seed, step, stream id, candidate index, purpose tag).
//  * Birth marks (open uniforms) are keyed per birth; the source emission law of
//    source owner v5 turns them into energy/angle/position (not reimplemented here).
//  * Electron birth time: tau = begin + u_time * duration (uniform held exposure,
//    the same binary64 expression as source v3); remaining flight = end - tau.
//    Neutral release times are source-v5 owned and bound back after its prepare.
//  * IDs: birth key (step, stream id, index) -> logical packet id
//    "arrival-v1/s<step>/p<stream>/k<index>", global ordinal = next_ordinal + rank,
//    native id = native_birth_base + rank, rank = position in canonical order.
//  * Ledger: one canonical binary record per accepted step, chained by SHA-256.
// Every function is pure with respect to the accepted Owner state; Owner changes
// only in commit() and restore().
#pragma once
#include "pf_arrival/canonical.hpp"
#include "pf_arrival/philox_core.h"
// The existing sealed ExactParticleCount header (SHA256 46e01fe2...), included as is.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include PF_ARRIVAL_EXACT_COUNT_HEADER
#pragma GCC diagnostic pop
#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <span>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace pf_arrival_v11 {
using ::pf_arrival::Writer;using ::pf_arrival::Reader;using ::pf_arrival::EncodingError;using ::pf_arrival::sha256_hex;using ::pf_arrival::hex_bytes;

using Count = raikou::mfem_backend::ExactParticleCount;

inline constexpr const char* policy_version = "pf-arrival-policy-v1.1";
inline constexpr const char* rng_version = "philox4x32-10/stream-layout-v1/u01-open-52/bernoulli-exact-17x64/stochastic-round-v1";
inline constexpr const char* electron_time_law = "uniform-held-exposure-v1";
inline constexpr const char* neutral_time_law = "source-v5-conditional-release-time";
// v1.1: quantized neutral release (no floor leaf): n = SR(D / w0) leaves of weight
// exactly w0, the source debits exactly n * w0 from its exact stock and the
// remainder stays in the stock (E[n w0] = D); n is capped at the largest m with
// m w0 <= stock (recorded as capped).
inline constexpr const char* neutral_release_rule = "quantized-release-v1.1: n=SR(D/w0), leaves of exactly w0, debit n*w0, remainder kept, n<=floor(stock/w0)";
inline constexpr const char* id_namespace = "arrival-v1";
inline constexpr std::uint32_t source_cap = 4096;   // source v5 candidate/leaf caps
inline constexpr std::uint32_t patch_cap = 256;     // source v5 patch cap

// ---------------------------------------------------------------- refusals
enum class RefusalKind : std::uint8_t {
    configuration = 1,          // invalid policy configuration
    step_index = 2,             // step outside 2^40 or not the next accepted step
    stream_input = 3,           // unknown/duplicate stream or kind mismatch
    expected_count = 4,         // nonfinite/negative expected count, finite log with zero scalar
    macroparticle_cap = 5,      // N/w0 or D/w0 above the per-stream cap
    step_cap = 6,               // per-step electron/neutral total above the cap
    debit_representation = 7,   // debit is not one positive binary64 carrier
    inventory = 8,              // debit exceeds the exact inventory
    clock = 9,                  // invalid exposure/step clock, or a bound time outside it
    id_space = 10,              // native ID or ordinal exhaustion
    stale = 11,                 // staged plan does not bind the accepted owner state
    native_receipt = 12,        // native birth receipt differs from the staged mapping
    checkpoint = 13,            // malformed or foreign checkpoint
    incomplete = 14             // births missing/duplicated, or times not bound at commit
};
inline const char* kind_name(RefusalKind k) {
    switch (k) {
        case RefusalKind::configuration: return "configuration";
        case RefusalKind::step_index: return "step-index";
        case RefusalKind::stream_input: return "stream-input";
        case RefusalKind::expected_count: return "expected-count";
        case RefusalKind::macroparticle_cap: return "macroparticle-cap";
        case RefusalKind::step_cap: return "step-cap";
        case RefusalKind::debit_representation: return "debit-representation";
        case RefusalKind::inventory: return "inventory";
        case RefusalKind::clock: return "clock";
        case RefusalKind::id_space: return "id-space";
        case RefusalKind::stale: return "stale";
        case RefusalKind::native_receipt: return "native-receipt";
        case RefusalKind::checkpoint: return "checkpoint";
        case RefusalKind::incomplete: return "incomplete";
    }
    return "unknown";
}

// Typed refusal. log_abs_scope = log|quantity| where meaningful (NaN otherwise);
// identity names step/stream/birth; observed/reference are the values compared.
struct Refusal : std::invalid_argument {
    RefusalKind kind;
    std::string quantity, identity, detail;
    double log_abs_scope, observed, reference;
    Refusal(RefusalKind k, std::string q, std::string id, std::string more = {},
            double scope = std::numeric_limits<double>::quiet_NaN(),
            double obs = std::numeric_limits<double>::quiet_NaN(),
            double ref = std::numeric_limits<double>::quiet_NaN())
        : std::invalid_argument(message(k, q, id, more)), kind(k), quantity(std::move(q)), identity(std::move(id)),
          detail(std::move(more)), log_abs_scope(scope), observed(obs), reference(ref) {}
    static std::string message(RefusalKind k, const std::string& q, const std::string& id, const std::string& more) {
        std::string s = "arrival owner v1.1 refusal: kind=";
        s += kind_name(k);
        s += "; quantity=" + q + "; identity=" + id;
        if (!more.empty()) s += "; " + more;
        return s + "; whole step proposal refused";
    }
};

namespace detail {
inline std::string hexfloat(double x) {
    std::ostringstream s;
    s << std::hexfloat << x;
    return s.str();
}
inline void count(Writer& w, const Count& n) { w.f64s(n.positive_carriers()); }
inline Count read_count(Reader& r) {
    const auto w = r.f64s(Count::maximum_carriers);
    Count n;
    try {
        n = Count::from_positive_carriers(w);
    } catch (const std::exception& e) {
        throw EncodingError(std::string("exact count decode: ") + e.what());
    }
    if (n.positive_carriers() != w) throw EncodingError("noncanonical exact count carriers");
    return n;
}
inline Count credit_n(double w, std::uint64_t n) {
    Count c;
    for (std::uint64_t i = 0; i < n; ++i) c = c.credited(w);
    return c;
}
inline bool positive_normal(double x) { return std::isfinite(x) && x >= DBL_MIN; }
}  // namespace detail

// ---------------------------------------------------------------- configuration
enum class Kind : std::uint8_t { electron = 1, neutral = 2 };
inline const char* kind_name(Kind k) { return k == Kind::electron ? "electron" : "neutral"; }

struct StreamConfig {
    std::uint32_t stream_id = 0;      // RNG "source patch id", < 2^16, unique
    Kind kind = Kind::electron;
    std::uint32_t patch_index = 0;    // source v5 patch index (< 256)
    std::string species_state;        // native species state of the births
    double target_weight = 0;         // w0: positive normal binary64
    std::uint32_t maximum_births = 0; // per step and stream, 1..4096
};
struct Configuration {
    std::uint64_t run_seed = 0;
    std::vector<StreamConfig> streams;
    std::uint32_t maximum_electron_births_per_step = source_cap;
    std::uint32_t maximum_neutral_leaves_per_step = source_cap;
};

class Policy {
public:
    explicit Policy(Configuration c) : c_(std::move(c)) {
        auto bad = [](const std::string& why) { throw Refusal(RefusalKind::configuration, "arrival configuration", policy_version, why); };
        if (c_.streams.empty() || c_.streams.size() > PFA_STREAM_LIMIT) bad("1..65536 streams required");
        std::sort(c_.streams.begin(), c_.streams.end(), [](const auto& a, const auto& b) { return a.stream_id < b.stream_id; });
        std::set<std::uint32_t> neutral_patches;
        for (std::size_t i = 0; i < c_.streams.size(); ++i) {
            const auto& s = c_.streams[i];
            if (s.stream_id >= PFA_STREAM_LIMIT) bad("stream id must be < 2^16");
            if (i && c_.streams[i - 1].stream_id == s.stream_id) bad("duplicate stream id");
            if (s.kind != Kind::electron && s.kind != Kind::neutral) bad("stream kind");
            if (s.patch_index >= patch_cap) bad("patch index must be < 256");
            if (s.species_state.empty() || s.species_state.size() > 256) bad("species state id required");
            if (!detail::positive_normal(s.target_weight)) bad("target weight must be a positive normal binary64");
            if (s.maximum_births < 1 || s.maximum_births > source_cap) bad("per-stream maximum births must be 1..4096");
            if (s.kind == Kind::neutral && !neutral_patches.insert(s.patch_index).second)
                bad("one neutral stream per source patch (source v5 releases one exact debit per patch)");
            index_.emplace(s.stream_id, i);
        }
        if (c_.maximum_electron_births_per_step < 1 || c_.maximum_electron_births_per_step > source_cap ||
            c_.maximum_neutral_leaves_per_step < 1 || c_.maximum_neutral_leaves_per_step > source_cap)
            bad("per-step caps must be 1..4096");
        Writer w;
        w.tag("pf-arrival-configuration-v1.1");
        w.text(policy_version);
        w.text(rng_version);
        w.text(electron_time_law);
        w.text(neutral_time_law);
        w.text(neutral_release_rule);
        w.text(id_namespace);
        w.u64(c_.run_seed);
        w.u64(c_.maximum_electron_births_per_step);
        w.u64(c_.maximum_neutral_leaves_per_step);
        w.u64(c_.streams.size());
        for (const auto& s : c_.streams) {
            w.u64(s.stream_id);
            w.u8(std::uint8_t(s.kind));
            w.u64(s.patch_index);
            w.text(s.species_state);
            w.f64(s.target_weight);
            w.u64(s.maximum_births);
        }
        sha_ = w.sha256();
    }
    const Configuration& configuration() const noexcept { return c_; }
    const std::string& configuration_sha256() const noexcept { return sha_; }
    // Bind this string as source v5 SourceSettings::sampling_policy_identity.
    std::string sampling_policy_identity() const { return std::string("pf-arrival-owner-v1.1/") + sha_; }
    bool has_stream(std::uint32_t id) const { return index_.contains(id); }
    const StreamConfig& stream(std::uint32_t id) const {
        auto it = index_.find(id);
        if (it == index_.end()) throw Refusal(RefusalKind::stream_input, "stream id", std::to_string(id), "not configured");
        return c_.streams[it->second];
    }
    pfa::StreamKey key(std::uint64_t step, std::uint32_t stream, std::uint32_t candidate, std::uint32_t purpose) const {
        if (step >= PFA_STEP_LIMIT) throw Refusal(RefusalKind::step_index, "accepted step index", std::to_string(step), "must be < 2^40");
        return {c_.run_seed, step, stream, candidate, purpose};
    }

private:
    Configuration c_;
    std::string sha_;
    std::map<std::uint32_t, std::size_t> index_;
};

// ---------------------------------------------------------------- per-step inputs
// One prepared electron exposure (source v5 PreparedPatch + its v3 Exposure).
struct ElectronExposure {
    std::uint32_t stream_id = 0;
    double expected_number = 0;                 // PreparedPatch::expected_number()
    double log_expected_number = -INFINITY;     // PreparedPatch::log_expected_number()
    bool has_log = false;                       // finite log count defined
    double begin_s = 0, duration_s = 0;         // exposure
    double end_s = 0;                           // common step end (frame end_s)
};
// One neutral release (source v5 ReleaseQuote) with the exact accepted inventory.
struct NeutralRelease {
    std::uint32_t stream_id = 0;
    Count debit;                                // ReleaseQuote::released (exact)
    double log_debit = -INFINITY;               // ReleaseQuote::log_released
    Count inventory;                            // Owner::inventory(patch) (exact)
    double begin_s = 0, end_s = 0;              // epoch t and t + dt
};
struct StepInput {
    std::uint64_t step = 0;
    std::vector<ElectronExposure> electrons;
    std::vector<NeutralRelease> neutrals;
};

// ---------------------------------------------------------------- decisions and births
struct Decision {
    std::uint64_t step = 0;
    std::uint32_t stream_id = 0, patch_index = 0;
    Kind kind = Kind::electron;
    double expected = 0;              // electron: N; neutral: d (binary64 value of the debit)
    double log_expected = -INFINITY;  // electron: L; neutral: supplied log debit
    bool has_log = false;
    double target_weight = 0;         // w0
    double ratio = 0;                 // x = expected / w0 (one IEEE division)
    std::uint64_t rounded = 0;        // SR(x)
    std::uint64_t target = 0;         // electron: rounded; neutral: max(1, rounded) when d > 0
    std::uint64_t births = 0;         // after the merge rule
    bool merged = false, floored = false;
    std::uint64_t quanta = 0;         // neutral: M with d = M * 2^exponent
    std::int64_t exponent = 0;
    Count debit;                      // exact total of the births' weights
    Count inventory_before;           // neutral: exact inventory before the debit
    double begin_s = 0, duration_s = 0, end_s = 0;
};

struct Birth {
    std::uint64_t step = 0;
    std::uint32_t stream_id = 0, index = 0, patch_index = 0;
    Kind kind = Kind::electron;
    double weight = 0;
    std::array<double, 6> marks{};    // electron: normal, transverse, azimuth, time (4 used); neutral: 6
    double time_s = std::numeric_limits<double>::quiet_NaN();
    double remaining_s = std::numeric_limits<double>::quiet_NaN();
    std::uint64_t ordinal = 0, native_id = 0;   // assigned at staging
    std::string packet_id() const {
        return std::string(id_namespace) + "/s" + std::to_string(step) + "/p" + std::to_string(stream_id) + "/k" + std::to_string(index);
    }
    bool before(const Birth& o) const { return stream_id != o.stream_id ? stream_id < o.stream_id : index < o.index; }
};

namespace detail {
inline std::string stream_identity(std::uint64_t step, std::uint32_t stream) {
    return "step=" + std::to_string(step) + ";stream=" + std::to_string(stream);
}
inline void need_clock(bool ok, std::uint64_t step, std::uint32_t stream, const char* why) {
    if (!ok) throw Refusal(RefusalKind::clock, "birth clock", stream_identity(step, stream), why);
}
}  // namespace detail

inline Decision decide(const Policy& P, std::uint64_t step, const ElectronExposure& e) {
    const auto& sc = P.stream(e.stream_id);
    const auto id = detail::stream_identity(step, e.stream_id);
    if (sc.kind != Kind::electron) throw Refusal(RefusalKind::stream_input, "stream kind", id, "electron exposure on a neutral stream");
    (void)P.key(step, e.stream_id, 0, PFA_PURPOSE_COUNT);   // step range
    detail::need_clock(std::isfinite(e.begin_s) && std::isfinite(e.duration_s) && std::isfinite(e.end_s) && e.duration_s > 0 &&
                           e.begin_s + e.duration_s <= e.end_s && e.begin_s < e.end_s,
                       step, e.stream_id, "finite exposure inside the common step required");
    Decision d;
    d.step = step;
    d.stream_id = e.stream_id;
    d.patch_index = sc.patch_index;
    d.kind = Kind::electron;
    d.expected = e.expected_number;
    d.log_expected = e.has_log ? e.log_expected_number : -INFINITY;
    d.has_log = e.has_log;
    d.target_weight = sc.target_weight;
    d.begin_s = e.begin_s;
    d.duration_s = e.duration_s;
    d.end_s = e.end_s;
    if (!(std::isfinite(e.expected_number) && e.expected_number >= 0))
        throw Refusal(RefusalKind::expected_count, "electron expected exposure count", id, "nonfinite or negative scalar",
                      e.has_log ? e.log_expected_number : std::numeric_limits<double>::quiet_NaN(), e.expected_number);
    if (e.has_log && !std::isfinite(e.log_expected_number))
        throw Refusal(RefusalKind::expected_count, "electron expected exposure count", id, "log count flagged but not finite",
                      std::numeric_limits<double>::quiet_NaN(), e.expected_number);
    if (e.has_log && e.expected_number == 0)
        throw Refusal(RefusalKind::expected_count, "electron expected exposure count", id,
                      "finite log count with a zero scalar (unrepresented positive count)", e.log_expected_number, 0.);
    if (e.expected_number == 0) return d;   // no flux defined: exact zero, no draw
    d.ratio = e.expected_number / sc.target_weight;
    if (!(d.ratio <= double(sc.maximum_births)))
        throw Refusal(RefusalKind::macroparticle_cap, "electron macroparticles N/w0", id,
                      "per-stream cap " + std::to_string(sc.maximum_births),
                      e.has_log ? e.log_expected_number - std::log(sc.target_weight) : std::log(d.ratio), d.ratio,
                      double(sc.maximum_births));
    d.rounded = d.target = d.births = pfa::stochastic_round(P.key(step, e.stream_id, 0, PFA_PURPOSE_COUNT), d.ratio);
    d.debit = detail::credit_n(sc.target_weight, d.births);
    return d;
}

inline Decision decide(const Policy& P, std::uint64_t step, const NeutralRelease& r) {
    const auto& sc = P.stream(r.stream_id);
    const auto id = detail::stream_identity(step, r.stream_id);
    if (sc.kind != Kind::neutral) throw Refusal(RefusalKind::stream_input, "stream kind", id, "neutral release on an electron stream");
    (void)P.key(step, r.stream_id, 0, PFA_PURPOSE_COUNT);
    detail::need_clock(std::isfinite(r.begin_s) && std::isfinite(r.end_s) && r.begin_s < r.end_s, step, r.stream_id,
                       "finite step clock required");
    Decision d;
    d.step = step;
    d.stream_id = r.stream_id;
    d.patch_index = sc.patch_index;
    d.kind = Kind::neutral;
    d.log_expected = r.log_debit;
    d.has_log = std::isfinite(r.log_debit);
    d.target_weight = sc.target_weight;
    d.begin_s = r.begin_s;
    d.end_s = r.end_s;
    d.duration_s = r.end_s - r.begin_s;
    d.inventory_before = r.inventory;
    const auto carriers = r.debit.positive_carriers();
    if (carriers.size() > 1)
        throw Refusal(RefusalKind::debit_representation, "neutral release expectation", id,
                      "quoted release has " + std::to_string(carriers.size()) + " positive carriers; v1.1 reads one binary64 value",
                      std::log(r.debit.rounded_for_rate()), r.debit.rounded_for_rate());
    if (carriers.empty()) return d;   // zero quote: no leaves, no draw
    const double v = carriers.front();
    d.expected = v;
    d.ratio = v / sc.target_weight;
    if (!(d.ratio <= double(sc.maximum_births)))
        throw Refusal(RefusalKind::macroparticle_cap, "neutral leaves D/w0", id, "per-stream cap " + std::to_string(sc.maximum_births),
                      std::log(v) - std::log(sc.target_weight), d.ratio, double(sc.maximum_births));
    d.rounded = pfa::stochastic_round(P.key(step, r.stream_id, 0, PFA_PURPOSE_COUNT), d.ratio);
    d.target = d.rounded;
    // cap at the exact stock: largest m <= target with m * w0 <= inventory (exact count)
    std::uint64_t m = d.target;
    Count debit = detail::credit_n(sc.target_weight, m);
    auto covers = [&](const Count& c) { Count rest = r.inventory;
        try { for (double w : c.positive_carriers()) rest = rest.debit(w); } catch (const std::exception&) { return false; } return true; };
    if (m > 0 && !covers(debit)) {
        const double fit = std::floor(r.inventory.rounded_for_rate() / sc.target_weight);
        m = std::min<std::uint64_t>(m, fit > 0 ? std::uint64_t(fit) : 0);
        debit = detail::credit_n(sc.target_weight, m);
        while (m > 0 && !covers(debit)) { --m; debit = detail::credit_n(sc.target_weight, m); }
    }
    d.merged = m < d.target;   // v1.1: "capped at the exact stock"
    d.floored = false;
    d.births = m;
    d.debit = debit;
    d.quanta = 0;
    d.exponent = 0;
    return d;
}

// Birth k of a decision: a pure function of (policy, decision, k).
inline Birth evaluate_birth(const Policy& P, const Decision& d, std::uint32_t k) {
    if (k >= d.births) throw Refusal(RefusalKind::incomplete, "birth index", detail::stream_identity(d.step, d.stream_id), "index beyond births");
    Birth b;
    b.step = d.step;
    b.stream_id = d.stream_id;
    b.index = k;
    b.patch_index = d.patch_index;
    b.kind = d.kind;
    if (d.kind == Kind::electron) {
        b.weight = d.target_weight;
        pfa::electron_marks(P.key(d.step, d.stream_id, k, PFA_PURPOSE_ELECTRON_MARKS), b.marks.data());
        // Same binary64 expression as source v3 BirthTime::UniformHeldExposure.
        b.time_s = d.begin_s + b.marks[3] * d.duration_s;
        b.remaining_s = d.end_s - b.time_s;
        detail::need_clock(b.time_s >= d.begin_s && b.time_s <= d.begin_s + d.duration_s && b.remaining_s >= 0, d.step, d.stream_id,
                           "electron birth time outside its exposure");
    } else {
        b.weight = d.target_weight;   // v1.1: every neutral leaf is exactly w0
        pfa::neutral_marks(P.key(d.step, d.stream_id, k, PFA_PURPOSE_NEUTRAL_MARKS), b.marks.data());
    }
    return b;
}

// ---------------------------------------------------------------- step plan
struct WorkItem {
    std::uint32_t stream_id = 0, index = 0;
};

namespace detail {
inline void encode_decision(Writer& w, const Decision& d) {
    w.u64(d.stream_id);
    w.u8(std::uint8_t(d.kind));
    w.u64(d.patch_index);
    w.f64(d.expected);
    w.f64(d.log_expected);
    w.boolean(d.has_log);
    w.f64(d.target_weight);
    w.f64(d.ratio);
    w.u64(d.rounded);
    w.u64(d.target);
    w.u64(d.births);
    w.boolean(d.merged);
    w.boolean(d.floored);
    w.u64(d.quanta);
    w.i64(d.exponent);
    count(w, d.debit);
    count(w, d.inventory_before);
    w.f64(d.begin_s);
    w.f64(d.duration_s);
    w.f64(d.end_s);
}
inline void encode_birth(Writer& w, const Birth& b, bool assigned) {
    w.u64(b.step);
    w.u64(b.stream_id);
    w.u64(b.index);
    w.u8(std::uint8_t(b.kind));
    w.f64(b.weight);
    for (double m : b.marks) w.f64(m);
    w.f64(b.time_s);
    w.f64(b.remaining_s);
    if (assigned) {
        w.text(b.packet_id());
        w.u64(b.ordinal);
        w.u64(b.native_id);
    }
}
}  // namespace detail

class StepPlan {
public:
    std::uint64_t step() const noexcept { return step_; }
    const std::string& configuration_sha256() const noexcept { return config_; }
    const std::vector<Decision>& decisions() const noexcept { return decisions_; }
    const std::vector<Birth>& births() const noexcept { return births_; }
    std::uint64_t electron_births() const noexcept { return electrons_; }
    std::uint64_t neutral_leaves() const noexcept { return neutrals_; }
    const std::string& sha256() const noexcept { return sha_; }
    std::string encode() const {
        Writer w;
        w.tag("pf-arrival-plan-v1");
        w.text(config_);
        w.u64(step_);
        w.u64(decisions_.size());
        for (const auto& d : decisions_) detail::encode_decision(w, d);
        w.u64(births_.size());
        for (const auto& b : births_) detail::encode_birth(w, b, false);
        return w.bytes;
    }

private:
    friend StepPlan assemble(const Policy&, std::uint64_t, std::vector<Decision>, std::vector<Birth>);
    friend class Staged;
    friend class Owner;
    std::uint64_t step_ = 0;
    std::string config_;
    std::vector<Decision> decisions_;
    std::vector<Birth> births_;
    std::uint64_t electrons_ = 0, neutrals_ = 0;
    std::string sha_;
};

// All stream decisions of a step, validated and sorted by stream id.
inline std::vector<Decision> decide_all(const Policy& P, const StepInput& in) {
    (void)P.key(in.step, 0, 0, PFA_PURPOSE_COUNT);
    std::vector<Decision> out;
    std::set<std::uint32_t> seen;
    auto unique = [&](std::uint32_t id) {
        if (!seen.insert(id).second)
            throw Refusal(RefusalKind::stream_input, "stream id", detail::stream_identity(in.step, id), "duplicate stream input in one step");
    };
    for (const auto& e : in.electrons) {
        unique(e.stream_id);
        out.push_back(decide(P, in.step, e));
    }
    for (const auto& r : in.neutrals) {
        unique(r.stream_id);
        out.push_back(decide(P, in.step, r));
    }
    std::sort(out.begin(), out.end(), [](const auto& a, const auto& b) { return a.stream_id < b.stream_id; });
    return out;
}
// Every (stream, index) birth of the decisions, in canonical order.
inline std::vector<WorkItem> work_items(const std::vector<Decision>& ds) {
    std::vector<WorkItem> w;
    for (const auto& d : ds)
        for (std::uint64_t k = 0; k < d.births; ++k) w.push_back({d.stream_id, std::uint32_t(k)});
    return w;
}

// Assembles births evaluated in any order/partition into the canonical plan and
// re-verifies coverage, caps and the exact debit equality of every stream.
inline StepPlan assemble(const Policy& P, std::uint64_t step, std::vector<Decision> ds, std::vector<Birth> births) {
    std::sort(ds.begin(), ds.end(), [](const auto& a, const auto& b) { return a.stream_id < b.stream_id; });
    std::sort(births.begin(), births.end(), [](const Birth& a, const Birth& b) { return a.before(b); });
    StepPlan p;
    p.step_ = step;
    p.config_ = P.configuration_sha256();
    std::size_t at = 0;
    for (std::size_t i = 0; i < ds.size(); ++i) {
        const auto& d = ds[i];
        const auto id = detail::stream_identity(step, d.stream_id);
        if (d.step != step || (i && ds[i - 1].stream_id == d.stream_id))
            throw Refusal(RefusalKind::stream_input, "stream decision", id, "foreign step or duplicate stream decision");
        Count total;
        for (std::uint64_t k = 0; k < d.births; ++k, ++at) {
            if (at >= births.size() || births[at].stream_id != d.stream_id || births[at].index != k || births[at].step != step)
                throw Refusal(RefusalKind::incomplete, "births", id, "missing, duplicated or foreign birth at index " + std::to_string(k));
            const auto& b = births[at];
            if (!(std::isfinite(b.weight) && b.weight > 0))
                throw Refusal(RefusalKind::incomplete, "birth weight", id, "nonpositive weight");
            total = total.credited(b.weight);
            (d.kind == Kind::electron ? p.electrons_ : p.neutrals_) += 1;
        }
        if (total != d.debit)
            throw Refusal(RefusalKind::debit_representation, "exact leaf sum", id, "leaves do not sum exactly to the debit");
    }
    if (at != births.size()) throw Refusal(RefusalKind::incomplete, "births", "step=" + std::to_string(step), "births of no decided stream");
    const auto& c = P.configuration();
    if (p.electrons_ > c.maximum_electron_births_per_step)
        throw Refusal(RefusalKind::step_cap, "electron births per step", "step=" + std::to_string(step), "", std::log(double(p.electrons_)),
                      double(p.electrons_), double(c.maximum_electron_births_per_step));
    if (p.neutrals_ > c.maximum_neutral_leaves_per_step)
        throw Refusal(RefusalKind::step_cap, "neutral leaves per step", "step=" + std::to_string(step), "", std::log(double(p.neutrals_)),
                      double(p.neutrals_), double(c.maximum_neutral_leaves_per_step));
    p.decisions_ = std::move(ds);
    p.births_ = std::move(births);
    p.sha_ = sha256_hex(p.encode());
    return p;
}

// Serial reference evaluation.
inline StepPlan plan_step(const Policy& P, const StepInput& in) {
    auto ds = decide_all(P, in);
    std::vector<Birth> births;
    for (const auto& d : ds)
        for (std::uint64_t k = 0; k < d.births; ++k) births.push_back(evaluate_birth(P, d, std::uint32_t(k)));
    return assemble(P, in.step, std::move(ds), std::move(births));
}

// ---------------------------------------------------------------- native receipt
// What the native birth credit (H2) returns for one step: the birth base it used,
// its next native ID afterwards and, in ascending native ID order, each birth's
// logical packet ID, species state, native ID and carrier weight.
struct NativeBirthEntry {
    std::string packet_id, state;
    std::uint64_t native_id = 0;
    double weight = 0;
};
struct NativeBirthReceipt {
    std::uint64_t base = 0, next_id_after = 0;
    std::vector<NativeBirthEntry> entries;
};

// ---------------------------------------------------------------- staged step
class Owner;
class Staged {
public:
    Staged(Staged&&) = default;
    Staged& operator=(Staged&&) = default;
    Staged(const Staged&) = delete;
    Staged& operator=(const Staged&) = delete;
    const StepPlan& plan() const noexcept { return plan_; }
    const std::vector<Birth>& births() const noexcept { return plan_.births_; }
    std::uint64_t native_base() const noexcept { return base_; }
    std::uint64_t first_ordinal() const noexcept { return ordinal_; }
    const std::string& owner_base_sha256() const noexcept { return owner_sha_; }
    // Opaque staged count-owner token (source v3 BirthSample::count_owner_token).
    const std::string& token() const noexcept { return token_; }
    // Binds source-owned birth times (neutral release times from source v5 packets),
    // in index order for one stream; remaining flight = end - time.
    void bind_source_times(std::uint32_t stream_id, std::span<const double> times) {
        const Decision* d = nullptr;
        for (const auto& x : plan_.decisions_)
            if (x.stream_id == stream_id) d = &x;
        const auto id = detail::stream_identity(plan_.step_, stream_id);
        if (!d || d->kind != Kind::neutral) throw Refusal(RefusalKind::stream_input, "source time binding", id, "no neutral stream decision");
        if (times.size() != d->births) throw Refusal(RefusalKind::incomplete, "source birth times", id, "one time per leaf required");
        std::size_t k = 0;
        for (auto& b : plan_.births_) {
            if (b.stream_id != stream_id) continue;
            const double t = times[k++];
            detail::need_clock(std::isfinite(t) && t >= d->begin_s && t <= d->end_s, plan_.step_, stream_id, "source birth time outside the step");
            b.time_s = t;
            b.remaining_s = d->end_s - t;
        }
    }

private:
    friend class Owner;
    Staged() = default;
    StepPlan plan_;
    std::uint64_t base_ = 0, ordinal_ = 0, generation_ = 0;
    std::string owner_sha_, token_;
};

// ---------------------------------------------------------------- ledger record
struct LedgerRecord {
    std::string bytes, record_sha256, chain_head;
};
struct StreamTotals {
    std::uint32_t stream_id = 0;
    Count expected_sum, sampled_sum;     // exact sums of expected (N or D) and of sampled debits
    std::uint64_t births = 0, steps = 0, merges = 0, floors = 0;
    double variance_sum = 0;             // sum of f(1-f), f = frac(N/w0), electrons (binary64, step order)
};

// ---------------------------------------------------------------- owner (accepted state)
class Owner {
public:
    Owner(std::shared_ptr<const Policy> p, std::uint64_t first_step) : p_(std::move(p)) {
        if (!p_) throw Refusal(RefusalKind::configuration, "policy", "owner", "policy required");
        if (first_step >= PFA_STEP_LIMIT) throw Refusal(RefusalKind::step_index, "first step", std::to_string(first_step), "must be < 2^40");
        next_step_ = first_step;
        Writer w;
        w.tag("pf-arrival-ledger-genesis-v1");
        w.text(p_->configuration_sha256());
        w.u64(first_step);
        head_ = w.sha256();
        for (const auto& s : p_->configuration().streams) {
            StreamTotals t;
            t.stream_id = s.stream_id;
            totals_.push_back(std::move(t));
        }
    }
    Owner(const Owner&) = delete;
    Owner& operator=(const Owner&) = delete;

    const Policy& policy() const noexcept { return *p_; }
    std::uint64_t next_step() const noexcept { return next_step_; }
    std::uint64_t next_ordinal() const noexcept { return next_ordinal_; }
    const std::string& ledger_head() const noexcept { return head_; }
    const std::vector<StreamTotals>& totals() const noexcept { return totals_; }
    bool has_anchor() const noexcept { return has_anchor_; }
    std::array<std::uint64_t, 4> last_anchor() const noexcept { return {anchor_step_, anchor_base_, anchor_count_, anchor_ordinal_}; }

    std::string state_sha256() const { return sha256_hex(encode_state()); }

    // Plan of the next accepted step (pure). Refuses any other step index.
    StepPlan plan(const StepInput& in) const {
        require_step(in.step);
        return plan_step(*p_, in);
    }
    // Validates a plan assembled from partitioned evaluation against this owner.
    void require_plan(const StepPlan& s) const {
        require_step(s.step());
        if (s.configuration_sha256() != p_->configuration_sha256() || sha256_hex(s.encode()) != s.sha256())
            throw Refusal(RefusalKind::stale, "step plan", "step=" + std::to_string(s.step()), "foreign configuration or altered plan");
    }

    // Assigns ordinals and native IDs: rank r in canonical order (stream id, index)
    // gets ordinal next_ordinal + r and native ID native_birth_base + r.
    Staged stage(StepPlan plan, std::uint64_t native_birth_base) const {
        require_plan(plan);
        const std::uint64_t n = plan.births().size();
        const auto id = "step=" + std::to_string(plan.step());
        if (native_birth_base == 0 || native_birth_base > std::numeric_limits<std::uint64_t>::max() - n)
            throw Refusal(RefusalKind::id_space, "native birth ID block", id, "base must be >= 1 and base + n <= 2^64 - 1",
                          std::numeric_limits<double>::quiet_NaN(), double(native_birth_base), double(n));
        if (next_ordinal_ > std::numeric_limits<std::uint64_t>::max() - n)
            throw Refusal(RefusalKind::id_space, "birth ordinal", id, "ordinal exhaustion");
        Staged s;
        s.plan_ = std::move(plan);
        s.base_ = native_birth_base;
        s.ordinal_ = next_ordinal_;
        s.generation_ = generation_;
        s.owner_sha_ = state_sha256();
        for (std::uint64_t r = 0; r < n; ++r) {
            auto& b = s.plan_.births_[r];
            b.ordinal = next_ordinal_ + r;
            b.native_id = native_birth_base + r;
        }
        Writer w;
        w.tag("pf-arrival-staged-v1");
        w.text(s.owner_sha_);
        w.text(s.plan_.sha256());
        w.u64(native_birth_base);
        w.u64(next_ordinal_);
        s.token_ = w.sha256();
        return s;
    }

    // Verifies the native receipt against the staged mapping, then publishes the
    // step: ledger record, chain head, totals, ordinal counter, anchor, next step.
    LedgerRecord commit(Staged&& s, const NativeBirthReceipt& r) {
        const auto& plan = s.plan_;
        const auto id = "step=" + std::to_string(plan.step());
        if (s.generation_ != generation_ || s.owner_sha_ != state_sha256() || plan.step() != next_step_)
            throw Refusal(RefusalKind::stale, "staged step", id, "staged against another accepted owner state");
        const std::uint64_t n = plan.births().size();
        for (const auto& b : plan.births())
            if (!(std::isfinite(b.time_s) && std::isfinite(b.remaining_s)))
                throw Refusal(RefusalKind::incomplete, "birth time", b.packet_id(), "time not bound before commit");
        auto bad = [&](const std::string& why) { throw Refusal(RefusalKind::native_receipt, "native birth receipt", id, why); };
        if (r.base != s.base_) bad("native base differs from the staged base");
        if (r.next_id_after != s.base_ + n) bad("native next ID after births differs from base + births");
        if (r.entries.size() != n) bad("native birth count differs");
        for (std::uint64_t k = 0; k < n; ++k) {
            const auto& b = plan.births()[k];
            const auto& e = r.entries[k];
            if (e.native_id != b.native_id || e.packet_id != b.packet_id() || e.state != p_->stream(b.stream_id).species_state ||
                std::bit_cast<std::uint64_t>(e.weight) != std::bit_cast<std::uint64_t>(b.weight))
                bad("entry " + std::to_string(k) + " differs from the staged canonical mapping");
        }
        // Build the record and the next state without touching accepted state.
        Writer w;
        w.tag("pf-arrival-ledger-step-v1");
        w.text(p_->configuration_sha256());
        w.text(head_);
        w.u64(plan.step());
        w.text(plan.sha256());
        w.text(s.token_);
        w.u64(s.base_);
        w.u64(n);
        w.u64(s.ordinal_);
        w.u64(plan.electron_births());
        w.u64(plan.neutral_leaves());
        w.u64(plan.decisions().size());
        for (const auto& d : plan.decisions()) detail::encode_decision(w, d);
        Writer rows;
        rows.tag("pf-arrival-birth-rows-v1");
        rows.u64(n);
        for (const auto& b : plan.births()) detail::encode_birth(rows, b, true);
        w.text(rows.sha256());
        LedgerRecord rec;
        rec.bytes = std::move(w.bytes);
        rec.record_sha256 = sha256_hex(rec.bytes);
        Writer chain;
        chain.tag("pf-arrival-ledger-chain-v1");
        chain.text(head_);
        chain.text(rec.record_sha256);
        rec.chain_head = chain.sha256();
        auto totals = totals_;
        for (const auto& d : plan.decisions()) {
            auto it = std::lower_bound(totals.begin(), totals.end(), d.stream_id, [](const auto& t, std::uint32_t v) { return t.stream_id < v; });
            if (it == totals.end() || it->stream_id != d.stream_id) bad("decision of an unknown stream");
            it->expected_sum = it->expected_sum.credited(d.expected);
            for (double c : d.debit.positive_carriers()) it->sampled_sum = it->sampled_sum.credited(c);
            it->births += d.births;
            it->steps += 1;
            it->merges += d.merged ? 1 : 0;
            it->floors += d.floored ? 1 : 0;
            if (d.kind == Kind::electron && d.expected > 0) {
                const double f = d.ratio - double(std::uint64_t(d.ratio));
                it->variance_sum += f * (1 - f);
            }
        }
        // Publication (no failure after this point other than allocation inside swaps of
        // already built values).
        totals_.swap(totals);
        head_ = rec.chain_head;
        if (n) {
            has_anchor_ = true;
            anchor_step_ = plan.step();
            anchor_base_ = s.base_;
            anchor_count_ = n;
            anchor_ordinal_ = s.ordinal_;
        }
        next_ordinal_ += n;
        next_step_ += 1;
        ++generation_;
        return rec;
    }

    // Checkpoint: small canonical state with a SHA-256 trailer.
    std::string save() const {
        Writer w;
        w.bytes = encode_state();
        Writer t;
        t.tag("pf-arrival-owner-checkpoint-v1.1");
        t.text(w.bytes);
        t.text(sha256_hex(t.bytes));
        return t.bytes;
    }
    static std::unique_ptr<Owner> restore(std::shared_ptr<const Policy> p, std::string_view bytes) {
        auto bad = [](const std::string& why) -> void { throw Refusal(RefusalKind::checkpoint, "arrival owner checkpoint", "restore", why); };
        if (!p) bad("policy required");
        try {
            Reader outer(bytes);
            outer.expect_tag("pf-arrival-owner-checkpoint-v1.1");
            const auto state = outer.text(64u << 20);
            const auto body_end = outer.position();
            const auto sum = outer.text(64);
            outer.end();
            if (sum != sha256_hex(bytes.substr(0, body_end))) bad("checksum mismatch");
            Reader r(state);
            r.expect_tag("pf-arrival-owner-state-v1.1");
            if (r.text(64) != p->configuration_sha256()) bad("configuration identity differs");
            auto o = std::make_unique<Owner>(p, 0);
            o->next_step_ = r.u64();
            o->next_ordinal_ = r.u64();
            o->has_anchor_ = r.boolean();
            o->anchor_step_ = r.u64();
            o->anchor_base_ = r.u64();
            o->anchor_count_ = r.u64();
            o->anchor_ordinal_ = r.u64();
            o->head_ = r.text(64);
            const auto n = r.u64();
            if (n != p->configuration().streams.size()) bad("stream totals dimension");
            for (std::uint64_t i = 0; i < n; ++i) {
                auto& t = o->totals_[i];
                if (r.u64() != t.stream_id) bad("stream totals order/identity");
                t.expected_sum = detail::read_count(r);
                t.sampled_sum = detail::read_count(r);
                t.births = r.u64();
                t.steps = r.u64();
                t.merges = r.u64();
                t.floors = r.u64();
                t.variance_sum = r.f64();
                if (!(std::isfinite(t.variance_sum) && t.variance_sum >= 0)) bad("variance sum");
            }
            r.end();
            if (o->next_step_ >= PFA_STEP_LIMIT || o->head_.size() != 64) bad("state ranges");
            if (!o->has_anchor_ && (o->anchor_step_ || o->anchor_base_ || o->anchor_count_ || o->anchor_ordinal_)) bad("anchor fields");
            if (o->has_anchor_ && (o->anchor_count_ == 0 || o->anchor_ordinal_ + o->anchor_count_ != o->next_ordinal_ ||
                                   o->anchor_step_ >= o->next_step_))
                bad("anchor inconsistent with the ordinal counter");
            if (o->encode_state() != state) bad("noncanonical state encoding");
            return o;
        } catch (const Refusal&) {
            throw;
        } catch (const std::exception& e) {
            throw Refusal(RefusalKind::checkpoint, "arrival owner checkpoint", "restore", e.what());
        }
    }

private:
    void require_step(std::uint64_t step) const {
        if (step != next_step_)
            throw Refusal(RefusalKind::step_index, "accepted step index", "step=" + std::to_string(step),
                          "the next accepted step is " + std::to_string(next_step_));
    }
    std::string encode_state() const {
        Writer w;
        w.tag("pf-arrival-owner-state-v1.1");
        w.text(p_->configuration_sha256());
        w.u64(next_step_);
        w.u64(next_ordinal_);
        w.boolean(has_anchor_);
        w.u64(anchor_step_);
        w.u64(anchor_base_);
        w.u64(anchor_count_);
        w.u64(anchor_ordinal_);
        w.text(head_);
        w.u64(totals_.size());
        for (const auto& t : totals_) {
            w.u64(t.stream_id);
            detail::count(w, t.expected_sum);
            detail::count(w, t.sampled_sum);
            w.u64(t.births);
            w.u64(t.steps);
            w.u64(t.merges);
            w.u64(t.floors);
            w.f64(t.variance_sum);
        }
        return w.bytes;
    }

    std::shared_ptr<const Policy> p_;
    std::uint64_t next_step_ = 0, next_ordinal_ = 0, generation_ = 0;
    bool has_anchor_ = false;
    std::uint64_t anchor_step_ = 0, anchor_base_ = 0, anchor_count_ = 0, anchor_ordinal_ = 0;
    std::string head_;
    std::vector<StreamTotals> totals_;
};

// Recomputes a chain head from the genesis head and step records (audit).
inline std::string chain(std::string head, const std::vector<LedgerRecord>& records) {
    for (const auto& r : records) {
        if (sha256_hex(r.bytes) != r.record_sha256) throw EncodingError("ledger record digest mismatch");
        Reader rd(r.bytes);
        rd.expect_tag("pf-arrival-ledger-step-v1");
        (void)rd.text(64);
        if (rd.text(64) != head) throw EncodingError("ledger record does not link to the previous head");
        Writer c;
        c.tag("pf-arrival-ledger-chain-v1");
        c.text(head);
        c.text(r.record_sha256);
        head = c.sha256();
        if (head != r.chain_head) throw EncodingError("ledger chain head mismatch");
    }
    return head;
}

// ---------------------------------------------------------------- refusal log (diagnostic)
// Attempt-level log of refused step proposals; never part of accepted state.
class RefusalLog {
public:
    explicit RefusalLog(const Policy& p) {
        Writer w;
        w.tag("pf-arrival-refusal-log-genesis-v1");
        w.text(p.configuration_sha256());
        head_ = w.sha256();
    }
    std::string append(std::uint64_t step, const Refusal& r) {
        Writer w;
        w.tag("pf-arrival-refusal-v1");
        w.text(head_);
        w.u64(step);
        w.u8(std::uint8_t(r.kind));
        w.text(r.quantity);
        w.text(r.identity);
        w.text(r.detail);
        w.f64(r.log_abs_scope);
        w.f64(r.observed);
        w.f64(r.reference);
        head_ = w.sha256();
        ++count_;
        return w.bytes;
    }
    const std::string& head() const noexcept { return head_; }
    std::uint64_t count() const noexcept { return count_; }

private:
    std::string head_;
    std::uint64_t count_ = 0;
};

}  // namespace pf_arrival_v11
