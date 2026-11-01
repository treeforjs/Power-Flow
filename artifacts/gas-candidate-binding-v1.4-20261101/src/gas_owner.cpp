// gas_owner.cpp -- gas target owner v1, host implementation (CPU reference,
// state, transactions, digests, codec). The device lane (gas_device.cpp)
// mirrors every mutation with the same gas_core.h functions.
#include "pf_gas_v1/gas_owner.hpp"
#ifdef PF_GAS_WITH_DEVICE
#include "pf_gas_v1/gas_device.hpp"
#else
namespace pf_gas_v1 {
class DeviceLane {};   // host-only build
}
#endif

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <map>
#include <set>
#include <unordered_map>

namespace pf_gas_v1 {
using pfg::i64;
using pfg::kDigits;
using pfg::u32;
using pfg::u64;

namespace {
constexpr std::size_t D = kDigits;
constexpr double kClosureTol = 0x1p-40;

[[noreturn]] void refuse(Refusal::Kind k, const std::string &what) { throw Refusal(k, what); }
void need(bool ok, Refusal::Kind k, const std::string &what) {
  if (!ok)
    refuse(k, what);
}

void validate_movers(const Configuration &cfg, std::span<const MoverRecord> movs) {
  std::set<u64> ids;
  struct Counts { double denominator = 0; Exact moved; };
  std::map<std::pair<u32, u32>, Counts> counts;
  const u64 nc = u64(cfg.grid.n[0]) * cfg.grid.n[1] * cfg.grid.n[2];
  for (const auto &m : movs) {
    need(ids.insert(m.native_id).second && m.species < cfg.species.size() &&
             std::isfinite(m.weight) && m.weight > 0 && std::isfinite(m.cell_count) &&
             m.cell_count >= m.weight && m.from_cell < nc &&
             (m.to_cell == pfg::kNone || m.to_cell < nc),
         Refusal::Kind::input, "input: duplicate or invalid mover");
    if (cfg.species[m.species].charge != 0) continue;
    auto &c = counts[{m.from_cell, m.species}];
    need(c.denominator == 0 || c.denominator == m.cell_count, Refusal::Kind::input,
         "input: inconsistent source-cell count");
    c.denominator = m.cell_count;
    c.moved.add_product(m.weight, 1.0);
  }
  for (const auto &[key, c] : counts) {
    (void)key;
    // Both the specified denominator and this sum are RN(exact count).
    // This validates inventory bounds; it does not alter rounded share transport.
    need(c.moved.value() <= c.denominator, Refusal::Kind::input,
         "input: movers exceed source-cell count");
  }
}

// Host SHA-256 on the forms core (same bytes as the device leaves).
struct Hash {
  pf7::Sha s;
  Hash() { pf7::sha_init(s); }
  void text(const std::string &t) { pf7::sha_text(s, t.data(), (pf7_u32)t.size()); }
  void u64v(u64 v) { pf7::sha_u64(s, v); }
  void f64(double x) { pf7::sha_f64(s, x); }
  void i64v(i64 v) { pf7::sha_u64(s, (u64)v); }
  void words(const pf7_u32 *w) { pf7::sha_digest_words(s, w); }
  void bytes(const void *p, std::size_t n) {
    const auto *b = static_cast<const unsigned char *>(p);
    for (std::size_t i = 0; i < n; ++i)
      pf7::sha_byte(s, b[i]);
  }
  std::array<pf7_u32, 8> finish() {
    std::array<pf7_u32, 8> o{};
    pf7::sha_finish(s, o.data());
    return o;
  }
  std::string hex() {
    const auto o = finish();
    static const char *hx = "0123456789abcdef";
    std::string r;
    for (auto w : o)
      for (int k = 7; k >= 0; --k)
        r += hx[(w >> (4 * k)) & 15];
    return r;
  }
};

std::vector<std::string> input_digests(std::span<const RecoilRecord> recs, std::span<const FlowRecord> flows) {
  std::vector<std::string> v; v.reserve(recs.size() + flows.size());
  for (const auto &r : recs) {
    Hash h; h.text("gas-v1.1/recoil-input"); h.u64v(r.accepted_step); h.u64v(r.stage); h.u64v(r.grid_cell);
    const auto &e = r.ev;
    h.u64v(e.carrier); h.u64v(e.candidate); h.u64v(e.process); h.u64v(e.cell); h.u64v(e.target_species);
    h.u64v(e.n_out); h.u64v(e.cover_class); h.u64v(e.internal_sink); h.u64v(e.flags);
    h.f64(e.weight); h.f64(e.time_s); for (auto x : e.u_in) h.f64(x); for (auto x : e.ut) h.f64(x);
    h.f64(e.E_tf); h.f64(e.ke_in_p); h.f64(e.ke_in_t); h.f64(e.q_eV); h.f64(e.e_internal);
    for (auto x : e.out_species) h.u64v(x);
    for (auto x : e.out_dest) h.u64v(x);
    for (const auto &u : e.out_u) for (auto x : u) h.f64(x);
    for (auto x : e.ke_out) h.f64(x);
    v.push_back(h.hex());
  }
  for (const auto &f : flows) {
    Hash h; h.text("gas-v1.1/flow-input"); h.u64v(f.kind); h.u64v(f.negative); h.u64v(f.grid_cell);
    h.u64v(f.species); h.u64v(f.native_id); h.f64(f.weight); for (auto x : f.u) h.f64(x); h.f64(f.internal_J);
    v.push_back(h.hex());
  }
  std::sort(v.begin(), v.end()); return v;
}
std::string key_table_digest(const std::vector<u64> &keys) {
  Hash h; h.text("gas-v1.1/key-table"); h.u64v(keys.size());
  for (auto k : keys) h.u64v(k);
  return h.hex();
}
std::string partial_digest(const Partial &p) {
  Hash h; h.text("gas-v1.1/partial"); h.text(p.configuration_sha256); h.text(p.key_table_sha256);
  h.u64v(p.accepted_step); h.u64v(p.stage); h.u64v(p.keys); h.text(p.candidate_sha256);
  h.u64v(p.input_sha256.size()); for (const auto &x : p.input_sha256) h.text(x);
  h.u64v(p.entry.size()); for (auto x : p.entry) h.i64v(x);
  h.u64v(p.global.size()); for (auto x : p.global) h.i64v(x);
  for (auto x : p.counter) h.u64v(x);
  return h.hex();
}

const char *fixed_name(u32 g) {
  static const char *names[] = {"W", "W_gas", "W_consume", "E_in_projectile", "E_in_target", "E_out_survivor",
                                "E_out_birth", "E_out_reservoir", "E_out_gas", "Q", "I_gas", "I_un_none",
                                "I_un_rot", "I_un_vib", "I_un_elec", "I_un_ion_internal", "I_un_fragment",
                                "defect_E", "P_in_projectile_x", "P_in_projectile_y", "P_in_projectile_z",
                                "P_in_target_x", "P_in_target_y", "P_in_target_z", "P_out_survivor_x",
                                "P_out_survivor_y", "P_out_survivor_z", "P_out_birth_x", "P_out_birth_y",
                                "P_out_birth_z", "P_out_reservoir_x", "P_out_reservoir_y", "P_out_reservoir_z",
                                "P_out_gas_x", "P_out_gas_y", "P_out_gas_z", "P_int_gas_x", "P_int_gas_y",
                                "P_int_gas_z", "P_int_un_x", "P_int_un_y", "P_int_un_z", "defect_P_x",
                                "defect_P_y", "defect_P_z", "reservoir_E", "reservoir_P_x", "reservoir_P_y",
                                "reservoir_P_z", "moved_E", "moved_P_x", "moved_P_y", "moved_P_z", "moved_U"};
  static_assert(sizeof(names) / sizeof(names[0]) == pfg::G_FIXED);
  return names[g];
}
}  // namespace

// ------------------------------------------------------------------ state
struct Owner::State {
  u64 step = 0;
  std::vector<u64> keys;
  std::vector<i64> entry;      // keys * 5 * D (normalised after every operation)
  std::vector<double> value;   // keys * 5
  std::vector<i64> global;     // G * D
  std::array<u64, pfg::C_COUNTERS> counter{};
};

namespace {
struct HostSink {
  i64 *e;
  i64 *g;
  u64 *c;
  void entry(i64 idx, u32 q, const pfg::Chunks &x, int negate) {
    pfg::add_chunks(e + ((u64)idx * pfg::kEntryQ + q) * D, x, negate);
  }
  void global(u32 idx, const pfg::Chunks &x, int negate) { pfg::add_chunks(g + (u64)idx * D, x, negate); }
  void status(u64 bits) { c[pfg::C_STATUS] |= bits; }
  void closure(double relE, double relP) {
    if (!(relE <= kClosureTol))
      ++c[pfg::C_CLOSURE_E];
    if (!(relP <= kClosureTol))
      ++c[pfg::C_CLOSURE_P];
    const u64 be = std::bit_cast<u64>(relE), bp = std::bit_cast<u64>(relP);
    if (be > c[pfg::C_MAXREL_E])
      c[pfg::C_MAXREL_E] = be;
    if (bp > c[pfg::C_MAXREL_P])
      c[pfg::C_MAXREL_P] = bp;
  }
};
void merge_counters(std::array<u64, pfg::C_COUNTERS> &a, const std::array<u64, pfg::C_COUNTERS> &b) {
  for (u32 k = 0; k < pfg::C_COUNTERS; ++k) {
    if (k == pfg::C_STATUS)
      a[k] |= b[k];
    else if (k == pfg::C_MAXREL_E || k == pfg::C_MAXREL_P)
      a[k] = std::max(a[k], b[k]);
    else
      a[k] += b[k];
  }
}
}  // namespace

// ------------------------------------------------------------------ text helpers
std::string digits_text(const std::int64_t *d) {
  int hi = (int)D - 1;
  while (hi > 0 && d[hi] == 0)
    --hi;
  int lo = 0;
  while (lo < hi && d[lo] == 0)
    ++lo;
  std::string s = "@" + std::to_string(lo) + ":";
  char buf[32];
  for (int k = hi; k >= lo; --k) {
    if (k == (int)D - 1 || d[k] < 0)
      std::snprintf(buf, sizeof buf, "%lld", (long long)d[k]);
    else
      std::snprintf(buf, sizeof buf, "%08llx", (unsigned long long)d[k]);
    s += buf;
    if (k > lo)
      s += ".";
  }
  return s;
}
ExactValue exact_value(const std::int64_t *digits) {
  std::array<i64, D> x;
  std::memcpy(x.data(), digits, sizeof(i64) * D);
  pf7::normalise(x.data());
  ExactValue v;
  v.digits = digits_text(x.data());
  pf7_u32 st = 0;
  auto y = x;
  v.value = pf7::round_sum(y.data(), &st);
  return v;
}

// ------------------------------------------------------------------ grid / configuration
std::uint32_t Grid::locate(const double *x) const noexcept {
  u64 idx[3];
  for (int d = 0; d < 3; ++d) {
    if (!(x[d] >= low[d] && x[d] <= high[d]))
      return pfg::kNone;
    const double f = (x[d] - low[d]) / (high[d] - low[d]);
    u64 i = u64(f * double(n[d]));
    if (i >= n[d])
      i = n[d] - 1;
    idx[d] = i;
  }
  return u32(idx[0] + u64(n[0]) * (idx[1] + u64(n[1]) * idx[2]));
}

Configuration preset(const std::string &name, Configuration base) {
  if (name == "production-gas-feedback") {
    base.feedback = true;
    base.transport = Transport::carrier_shares;
  } else if (name == "production-gas-ledger") {
    base.feedback = false;
    base.transport = Transport::carrier_shares;
  } else
    refuse(Refusal::Kind::configuration, "configuration: unknown preset " + name);
  base.preset = name;
  return base;
}

std::string configuration_sha256(const Configuration &c) {
  Hash h;
  h.text("gas-target-owner-v1.3-publication/configuration");
  h.text(c.preset);
  h.u64v(c.feedback ? 1 : 0);
  h.u64v(u64(c.transport));
  h.u64v(c.species.size());
  for (const auto &s : c.species) {
    h.text(s.state);
    h.f64(s.mass_kg);
    h.i64v(s.charge);
    h.u64v(s.chemistry_target ? 1 : 0);
  }
  h.u64v(c.processes.size());
  for (const auto &p : c.processes) {
    h.u64v(p.projectile);
    h.u64v(p.n_out);
    h.u64v(p.internal_carrier);
    h.u64v(p.kind);
  }
  h.text(c.chemistry_blob_sha256);
  h.text(c.chemistry_registry_sha256);
  for (int d = 0; d < 3; ++d) {
    h.f64(c.grid.low[d]);
    h.f64(c.grid.high[d]);
    h.u64v(c.grid.n[d]);
  }
  h.f64(c.grid.volume);
  h.u64v(c.max_entries);
  return h.hex();
}

// ------------------------------------------------------------------ owner
namespace {
std::string state_digest(const std::string &cfg_sha, u64 step, const std::array<u64, pfg::C_COUNTERS> &counter,
                         const std::vector<std::array<pf7_u32, 8>> &entry_leaves,
                         const std::vector<std::array<pf7_u32, 8>> &global_leaves) {
  Hash e;
  e.text("gas-target-owner-v1/entries");
  e.u64v(entry_leaves.size());
  for (const auto &l : entry_leaves)
    e.words(l.data());
  const auto er = e.finish();
  Hash g;
  g.text("gas-target-owner-v1/globals");
  g.u64v(global_leaves.size());
  for (const auto &l : global_leaves)
    g.words(l.data());
  const auto gr = g.finish();
  Hash h;
  h.text("gas-target-owner-v1/state-v1");
  h.text(cfg_sha);
  h.u64v(step);
  for (auto c : counter)
    h.u64v(c);
  h.words(er.data());
  h.words(gr.data());
  return h.hex();
}
void host_leaves(const Owner::State &s, u32 G, std::vector<std::array<pf7_u32, 8>> &el,
                 std::vector<std::array<pf7_u32, 8>> &gl) {
  el.clear();
  gl.clear();
  for (u64 i = 0; i < s.keys.size(); ++i) {
    const i64 *d = s.entry.data() + i * pfg::kEntryQ * D;
    if (pfg::digits_zero(d, int(pfg::kEntryQ * D)))
      continue;
    std::array<pf7_u32, 8> w;
    pfg::entry_leaf(s.keys[i], d, w.data());
    el.push_back(w);
  }
  for (u32 g = 0; g < G; ++g) {
    const i64 *d = s.global.data() + u64(g) * D;
    if (pfg::digits_zero(d, int(D)))
      continue;
    std::array<pf7_u32, 8> w;
    pfg::global_leaf(g, d, w.data());
    gl.push_back(w);
  }
}
void normalise_state(Owner::State &s, u32 G) {
  const u64 ne = s.keys.size() * pfg::kEntryQ;
  s.value.assign(ne, 0.0);
  for (u64 a = 0; a < ne; ++a) {
    i64 *d = s.entry.data() + a * D;
    pf7::normalise(d);
    std::array<i64, D> x;
    std::memcpy(x.data(), d, sizeof x);
    pf7_u32 st = 0;
    s.value[a] = pf7::round_sum(x.data(), &st);
    if (st)
      s.counter[pfg::C_STATUS] |= pfg::kStOverflow;
  }
  for (u32 g = 0; g < G; ++g)
    pf7::normalise(s.global.data() + u64(g) * D);
}
}  // namespace

Owner::Owner(const Configuration &c) : cfg_(c) {
  need(!c.species.empty() && c.species.size() < (1u << 20), Refusal::Kind::configuration, "configuration: species");
  need(c.processes.size() < (1u << 20), Refusal::Kind::configuration, "configuration: processes");
  for (int d = 0; d < 3; ++d)
    need(c.grid.n[d] > 0 && c.grid.high[d] > c.grid.low[d] && std::isfinite(c.grid.low[d]) &&
             std::isfinite(c.grid.high[d]),
         Refusal::Kind::configuration, "configuration: grid");
  need(std::isfinite(c.grid.volume) && c.grid.volume > 0, Refusal::Kind::configuration, "configuration: grid volume");
  u64 cell_count = 1;
  for (auto n : c.grid.n) {
    need(cell_count <= u64(pfg::kNone) / n, Refusal::Kind::configuration, "configuration: grid cell ID overflow");
    cell_count *= n;
  }
  const double expected_volume = ((c.grid.high[0] - c.grid.low[0]) / c.grid.n[0]) *
                                 ((c.grid.high[1] - c.grid.low[1]) / c.grid.n[1]) *
                                 ((c.grid.high[2] - c.grid.low[2]) / c.grid.n[2]);
  need(c.grid.volume == expected_volume, Refusal::Kind::configuration,
       "configuration: volume differs from native grid formula");
  need(c.max_entries > 0, Refusal::Kind::configuration, "configuration: max_entries");
  std::set<std::string> names;
  for (const auto &s : c.species) {
    need(std::isfinite(s.mass_kg) && s.mass_kg > 0 && names.insert(s.state).second, Refusal::Kind::configuration,
         "configuration: species " + s.state);
    pfg::SpeciesInfo x;
    x.mass_kg = s.mass_kg;
    x.charge = s.charge;
    x.flags = (s.charge == 0 ? pfg::kSpGas : 0u) | (s.chemistry_target ? pfg::kSpTarget : 0u) |
              ((s.charge == 0 && s.chemistry_target) ? pfg::kSpFeedback : 0u);
    species_.push_back(x);
  }
  for (const auto &p : c.processes) {
    need(p.projectile < c.species.size() && p.n_out >= 1 && p.n_out <= pfg::kMaxOut &&
             (p.internal_carrier == pfg::kNone || p.internal_carrier < p.n_out),
         Refusal::Kind::configuration, "configuration: process table");
    procs_.push_back({p.projectile, p.n_out, p.internal_carrier, p.kind});
  }
  if (procs_.empty())
    procs_.push_back({0, 1, pfg::kNone, 0});   // placeholder; a record naming it is refused (n_proc check)
  L_.S = u32(species_.size());
  L_.NP = u32(c.processes.size());
  cfg_sha_ = pf_gas_v1::configuration_sha256(c);
  acc_ = std::make_unique<State>();
  if (host())
    acc_->global.assign(u64(L_.total()) * D, 0);
#ifdef PF_GAS_WITH_DEVICE
  if (c.execution != Execution::host)
    dev_ = std::make_unique<DeviceLane>(c, species_, procs_, L_);
#else
  need(c.execution == Execution::host, Refusal::Kind::device, "device: built without the device lane");
#endif
}
Owner::~Owner() = default;

bool Owner::host() const noexcept { return cfg_.execution != Execution::device; }
std::uint64_t Owner::accepted_step() const noexcept { return acc_->step; }

void Owner::ensure_open(std::uint64_t step, const char *what) const {
  need(open_, Refusal::Kind::transaction, std::string("transaction: ") + what + " without an open step");
  need(!publication_ready_, Refusal::Kind::transaction, "transaction: candidate sealed for publication");
  need(!failed_, Refusal::Kind::transaction, "transaction: failed candidate requires abandon");
  need(step == open_step_, Refusal::Kind::transaction, std::string("transaction: ") + what + " for another step");
}

void Owner::invalidate_candidate_bindings() {
  need(candidate_revision_ != ~u64(0), Refusal::Kind::transaction,
       "transaction: candidate mutation revision exhausted");
  ++candidate_revision_;
}
Owner::CandidateBinding Owner::candidate_binding() const {
  ensure_open(open_step_, "candidate binding");
  need(bool(cand_), Refusal::Kind::transaction, "transaction: no candidate state");
  CandidateBinding token;
  token.lifetime_ = lifetime_; token.generation_ = generation_;
  token.revision_ = candidate_revision_; token.step_ = open_step_;
  return token;
}
bool Owner::matches_candidate(const CandidateBinding &token) const noexcept {
  return open_ && bool(cand_) && !publication_ready_ && !failed_ &&
         token.lifetime_ == lifetime_ && token.generation_ == generation_ &&
         token.revision_ == candidate_revision_ && token.step_ == open_step_;
}
bool Owner::matches_candidate(const CandidateBinding &token, std::uint64_t expected_step) const noexcept {
  return expected_step == open_step_ && matches_candidate(token);
}

void Owner::begin_step(std::uint64_t accepted_step) {
  need(!open_, Refusal::Kind::transaction, "transaction: step already open");
  need(generation_ != ~u64(0), Refusal::Kind::transaction, "transaction: lifetime generation exhausted");
  const bool genesis = acc_->counter[pfg::C_STAGES] == 0;
  need(genesis ? accepted_step <= 1 : acc_->step != ~u64(0) && accepted_step == acc_->step + 1,
       Refusal::Kind::transaction, "transaction: step must follow the accepted step (genesis permits 0 or 1)");
  cand_ = std::make_unique<State>(*acc_);
  ++generation_;
  candidate_revision_ = 0;
  publication_ready_ = false;
  cand_->step = accepted_step;
  open_ = true;
  open_step_ = accepted_step;
  next_stage_ = 1;
  partial_batch_ready_ = false;
  partial_batch_inputs_.clear();
  failed_ = false;
#ifdef PF_GAS_WITH_DEVICE
  if (dev_) {
    try { dev_->begin_step(accepted_step); }
    catch (...) { failed_ = true; throw; }
  }
#endif
}

void Owner::insert_keys(State &s, std::vector<u64> add) {
  std::sort(add.begin(), add.end());
  add.erase(std::unique(add.begin(), add.end()), add.end());
  std::vector<u64> fresh;
  for (u64 k : add)
    if (!std::binary_search(s.keys.begin(), s.keys.end(), k))
      fresh.push_back(k);
  if (fresh.empty())
    return;
  std::vector<u64> keys;
  keys.reserve(s.keys.size() + fresh.size());
  std::merge(s.keys.begin(), s.keys.end(), fresh.begin(), fresh.end(), std::back_inserter(keys));
  need(keys.size() <= cfg_.max_entries, Refusal::Kind::capacity, "capacity: more than max_entries gas entries");
  if (host()) {
    std::vector<i64> entry(keys.size() * pfg::kEntryQ * D, 0);
    std::vector<double> value(keys.size() * pfg::kEntryQ, 0.0);
    u64 j = 0;
    for (u64 i = 0; i < keys.size(); ++i)
      if (j < s.keys.size() && s.keys[j] == keys[i]) {
        std::memcpy(entry.data() + i * pfg::kEntryQ * D, s.entry.data() + j * pfg::kEntryQ * D,
                    sizeof(i64) * pfg::kEntryQ * D);
        std::memcpy(value.data() + i * pfg::kEntryQ, s.value.data() + j * pfg::kEntryQ, sizeof(double) * pfg::kEntryQ);
        ++j;
      }
    s.entry = std::move(entry);
    s.value = std::move(value);
  }
  s.keys = std::move(keys);
#ifdef PF_GAS_WITH_DEVICE
  if (dev_)
    dev_->ensure_keys(s.keys);
#endif
}

void Owner::prepare_keys(std::span<const RecoilRecord> recs, std::span<const FlowRecord> flows,
                         std::span<const MoverRecord> movs) {
  ensure_open(open_step_, "keys");
  invalidate_candidate_bindings();
  prepare_keys_impl(recs, flows, movs, true);
}
void Owner::prepare_keys_impl(std::span<const RecoilRecord> recs, std::span<const FlowRecord> flows,
                              std::span<const MoverRecord> movs, bool bind) {
  ensure_open(open_step_, "keys");
  validate_movers(cfg_, movs);
  std::vector<u64> add;
  const u32 S = L_.S;
  const u64 cells = u64(cfg_.grid.n[0]) * cfg_.grid.n[1] * cfg_.grid.n[2];
  auto gas = [&](u32 s) { return s < S && (species_[s].flags & pfg::kSpGas); };
  for (const auto &r : recs) {
    need(r.grid_cell < cells, Refusal::Kind::input, "input: recoil cell outside grid");
    const auto &ev = r.ev;
    if (gas(ev.target_species))
      add.push_back(pfg::entry_key(r.grid_cell, ev.target_species));
    for (u32 i = 0; i < ev.n_out && i < pfg::kMaxOut; ++i)
      if ((ev.out_dest[i] == pfg::kDestGas || ev.out_dest[i] == pfg::kDestReservoir) && gas(ev.out_species[i]))
        add.push_back(pfg::entry_key(r.grid_cell, ev.out_species[i]));
  }
  for (const auto &f : flows)
    if (f.kind == pfg::kFlowDebit && gas(f.species)) {
      need(f.grid_cell < cells, Refusal::Kind::input, "input: debit cell outside grid");
      add.push_back(pfg::entry_key(f.grid_cell, f.species));
    }
  for (const auto &m : movs) {
    if (gas(m.species))
      need(m.from_cell < cells && (m.to_cell == pfg::kNone || m.to_cell < cells), Refusal::Kind::input,
           "input: mover cell outside grid");
    if (gas(m.species) && m.to_cell != pfg::kNone &&
        std::binary_search(cand_->keys.begin(), cand_->keys.end(), pfg::entry_key(m.from_cell, m.species)))
      add.push_back(pfg::entry_key(m.to_cell, m.species));
  }
  failed_ = true; // Allocation/device failure may leave a partial key preparation.
  insert_keys(*cand_, std::move(add));
  partial_batch_ready_ = bind;
  partial_batch_inputs_ = bind ? input_digests(recs, flows) : std::vector<std::string>{};
  failed_ = false;
}

Partial Owner::accumulate(std::span<const RecoilRecord> recs, std::span<const FlowRecord> flows) const {
  return accumulate_impl(recs, flows, true);
}
Partial Owner::accumulate_impl(std::span<const RecoilRecord> recs, std::span<const FlowRecord> flows, bool bind) const {
  ensure_open(open_step_, "accumulate");
  if (bind) need(partial_batch_ready_, Refusal::Kind::transaction, "transaction: prepare complete input batch before partial accumulation");
  need(host(), Refusal::Kind::configuration, "configuration: host accumulation in device execution");
  const State &s = *cand_;
  Partial p;
  p.keys = s.keys.size();
  p.accepted_step = open_step_;
  p.stage = next_stage_;
  if (bind) p.key_table_sha256 = key_table_digest(s.keys);
  p.configuration_sha256 = cfg_sha_;
  if (bind) {
    p.candidate_sha256 = digest(s, true);
    p.input_sha256 = input_digests(recs, flows);
    need(std::includes(partial_batch_inputs_.begin(), partial_batch_inputs_.end(), p.input_sha256.begin(), p.input_sha256.end()),
         Refusal::Kind::transaction, "transaction: partial input outside prepared batch");
  }
  for (const auto &r : recs)
    need(r.accepted_step == open_step_ && r.stage == next_stage_, Refusal::Kind::input,
         "input: partial record of another step/stage");
  for (const auto &f : flows)
    if (!(f.kind == pfg::kFlowDebit || (next_stage_ == 2 && f.kind == pfg::kFlowRelease)))
      p.stage = 0; // External population flows cannot masquerade as a chemistry partial.

  p.entry.assign(s.keys.size() * pfg::kEntryQ * D, 0);
  p.global.assign(u64(L_.total()) * D, 0);
  HostSink sink{p.entry.data(), p.global.data(), p.counter.data()};
  pfg::Tables T{species_.data(), L_.S, procs_.data(), L_.NP, s.keys.data(), s.keys.size(), L_};
  for (const auto &r : recs)
    pfg::apply_record(sink, T, r);
  for (const auto &f : flows)
    pfg::apply_flow(sink, T, f);
  p.counter[pfg::C_RECORDS] += recs.size();
  p.counter[pfg::C_FLOWS] += flows.size();
  if (bind) p.payload_sha256 = partial_digest(p);
  return p;
}

// Device application (device and device_verified), compared with the host
// counters when both run; returns the counters of the operation.
std::array<u64, pfg::C_COUNTERS> Owner::device_apply(std::span<const RecoilRecord> recs,
                                                    std::span<const FlowRecord> flows,
                                                    std::span<const MoverRecord> movs,
                                                    const std::array<u64, pfg::C_COUNTERS> *host_counters) {
  std::array<u64, pfg::C_COUNTERS> c{};
#ifdef PF_GAS_WITH_DEVICE
  if (dev_) {
    c = dev_->apply(recs, flows, movs);
    c[pfg::C_RECORDS] = recs.size();
    c[pfg::C_FLOWS] = flows.size();
    if (host_counters) {
      auto h = *host_counters;
      h[pfg::C_STATUS] = 0;
      auto d = c;
      d[pfg::C_STATUS] = 0;
      need(h == d, Refusal::Kind::device, "device: operation counters differ from the host reference");
    }
  }
#else
  (void)recs;
  (void)flows;
  (void)movs;
  (void)host_counters;
#endif
  return c;
}

std::string Owner::finish(State &s, const char *what) {
  if (host()) {
    normalise_state(s, L_.total());
    need(s.counter[pfg::C_STATUS] == 0, Refusal::Kind::exactness,
         std::string("exactness: ") + what + " status " + std::to_string(s.counter[pfg::C_STATUS]));
  } else
    need(s.counter[pfg::C_STATUS] == 0, Refusal::Kind::exactness,
         std::string("exactness: ") + what + " status " + std::to_string(s.counter[pfg::C_STATUS]));
  const std::string d = digest(s, true);
  verify_device(what, d);
  return d;
}

std::string Owner::merge(std::uint64_t accepted_step, std::uint64_t stage, std::span<const Partial> parts) {
  ensure_open(accepted_step, "merge");
  invalidate_candidate_bindings();
  need((stage == 1 || stage == 2) && stage == next_stage_, Refusal::Kind::transaction,
       "transaction: chemistry stages must occur once in order");
  need(host() && !dev_, Refusal::Kind::configuration, "configuration: partial merge is the host reference path");
  State &s = *cand_;
  std::array<u64, pfg::C_COUNTERS> c{};
  need(partial_batch_ready_, Refusal::Kind::transaction, "transaction: no prepared partial batch");
  std::vector<std::string> inputs;
  const auto state_sha = digest(s, true);
  for (const auto &p : parts) {
    need(p.accepted_step == accepted_step && p.stage == stage && p.configuration_sha256 == cfg_sha_ &&
             p.key_table_sha256 == key_table_digest(s.keys) && p.payload_sha256 == partial_digest(p) &&
             p.candidate_sha256 == state_sha,
         Refusal::Kind::transaction, "transaction: stale or modified partial binding");
    need(p.keys == s.keys.size() && p.entry.size() == s.entry.size() && p.global.size() == s.global.size(),
         Refusal::Kind::transaction, "transaction: partial built on another key table");
    merge_counters(c, p.counter);
    inputs.insert(inputs.end(), p.input_sha256.begin(), p.input_sha256.end());
  }
  std::sort(inputs.begin(), inputs.end());
  need(inputs == partial_batch_inputs_, Refusal::Kind::transaction, "transaction: incomplete or duplicated partial inputs");
  need(c[pfg::C_STATUS] == 0, Refusal::Kind::input,
       "input: stage " + std::to_string(stage) + " status " + std::to_string(c[pfg::C_STATUS]));
  failed_ = true;
  for (const auto &p : parts) {
    for (u64 i = 0; i < s.entry.size(); ++i)
      s.entry[i] += p.entry[i];
    for (u64 i = 0; i < s.global.size(); ++i)
      s.global[i] += p.global[i];
  }
  merge_counters(s.counter, c);
  if (stage)
    ++s.counter[pfg::C_STAGES];
  const auto result = finish(s, "merge");
  ++next_stage_;
  partial_batch_ready_ = false;
  failed_ = false;
  return result;
}

std::string Owner::run(std::uint64_t stage, std::span<const RecoilRecord> recs, std::span<const FlowRecord> flows,
                       const char *what) {
  State &s = *cand_;
  Partial p;
  if (host()) {
    p = accumulate_impl(recs, flows, false);
    need(p.counter[pfg::C_STATUS] == 0, Refusal::Kind::input,
         std::string("input: ") + what + " status " + std::to_string(p.counter[pfg::C_STATUS]));
  }
  failed_ = true;
  const auto dc = device_apply(recs, flows, {}, host() ? &p.counter : nullptr);
  if (host()) {
    for (u64 i = 0; i < s.entry.size(); ++i)
      s.entry[i] += p.entry[i];
    for (u64 i = 0; i < s.global.size(); ++i)
      s.global[i] += p.global[i];
    merge_counters(s.counter, p.counter);
  } else
    merge_counters(s.counter, dc);
  if (stage)
    ++s.counter[pfg::C_STAGES];
  const auto result = finish(s, what);
  failed_ = false;
  return result;
}

std::string Owner::take(std::uint64_t accepted_step, std::uint64_t stage, std::span<const RecoilRecord> recs,
                        std::span<const FlowRecord> flows) {
  ensure_open(accepted_step, "take");
  invalidate_candidate_bindings();
  need(stage == 1 || stage == 2, Refusal::Kind::input, "input: stage");
  need(stage == next_stage_, Refusal::Kind::transaction, "transaction: chemistry stages must occur once in order");
  for (const auto &r : recs)
    need(r.accepted_step == accepted_step && r.stage == stage, Refusal::Kind::input,
         "input: record of another step/stage");
  for (const auto &f : flows)
    need(f.kind == pfg::kFlowDebit || (stage == 2 && f.kind == pfg::kFlowRelease), Refusal::Kind::input,
         "input: stage flows are debits (and stage-2 releases)");
  prepare_keys_impl(recs, flows, {}, false);
  const auto result = run(stage, recs, flows, "take");
  ++next_stage_;
  return result;
}

std::string Owner::apply_flows(std::uint64_t accepted_step, std::span<const FlowRecord> flows) {
  ensure_open(accepted_step, "flows");
  invalidate_candidate_bindings();
  for (const auto &f : flows)
    need(f.kind == pfg::kFlowBirth || f.kind == pfg::kFlowSink || f.kind == pfg::kFlowTransaction,
         Refusal::Kind::input, "input: chemistry flows require take");
  prepare_keys_impl({}, flows, {}, false);
  return run(0, {}, flows, "flows");
}

std::string Owner::transport(std::uint64_t accepted_step, std::span<const MoverRecord> movs) {
  ensure_open(accepted_step, "transport");
  invalidate_candidate_bindings();
  need(cfg_.transport == Transport::carrier_shares, Refusal::Kind::configuration,
       "configuration: transport is none in this configuration");
  prepare_keys_impl({}, {}, movs, false);
  State &s = *cand_;
  std::vector<i64> e, g;
  std::array<u64, pfg::C_COUNTERS> c{};
  if (host()) {
    e.assign(s.entry.size(), 0);
    g.assign(s.global.size(), 0);
    HostSink sink{e.data(), g.data(), c.data()};
    pfg::Tables T{species_.data(), L_.S, procs_.data(), L_.NP, s.keys.data(), s.keys.size(), L_};
    for (const auto &m : movs)
      pfg::apply_mover(sink, T, m, s.value.data());
    need(c[pfg::C_STATUS] == 0, Refusal::Kind::input, "input: transport status " + std::to_string(c[pfg::C_STATUS]));
  }
  failed_ = true;
  auto dc = device_apply({}, {}, movs, host() ? &c : nullptr);
  if (host()) {
    for (u64 i = 0; i < e.size(); ++i)
      s.entry[i] += e[i];
    for (u64 i = 0; i < g.size(); ++i)
      s.global[i] += g[i];
  } else
    c = dc;
  c[pfg::C_MOVERS] = movs.size();
  c[pfg::C_TRANSPORTS] = 1;
  merge_counters(s.counter, c);
  const auto result = finish(s, "transport");
  failed_ = false;
  return result;
}

FeedbackReceipt Owner::feedback(std::uint64_t accepted_step, std::uint64_t stage, std::span<const std::uint32_t> atlas,
                                std::span<CellTargetsAbi> cells) const {
  need(stage == 1 || stage == 2, Refusal::Kind::input, "input: feedback stage");
  if (open_) {
    ensure_open(accepted_step, "feedback");
    need(stage == next_stage_, Refusal::Kind::transaction, "transaction: feedback stage order");
  } else
    need(accepted_step == acc_->step, Refusal::Kind::transaction, "transaction: feedback for another accepted state");
  need(atlas.size() == cells.size(), Refusal::Kind::input, "input: feedback atlas size");
  if (open_)
    need(accepted_step == open_step_, Refusal::Kind::transaction, "transaction: feedback for another step");
  std::vector<CellTargetsAbi> candidate(cells.begin(), cells.end());
  auto cells_sha = [&](std::span<const CellTargetsAbi> data) {
    Hash h;
    h.text("gas-target-owner-v1/cell-targets");
    h.u64v(cells.size());
    h.bytes(atlas.data(), atlas.size() * 4);
    h.bytes(data.data(), data.size() * sizeof(CellTargetsAbi));
    return h.hex();
  };
  FeedbackReceipt r;
  r.cells_before_sha256 = cells_sha(cells);
  std::set<u32> occupied;
  const u64 total_cells = u64(cfg_.grid.n[0]) * cfg_.grid.n[1] * cfg_.grid.n[2];
  for (std::size_t i = 0; i < cells.size(); ++i) {
    const auto &c = cells[i];
    need(atlas[i] < total_cells && occupied.insert(atlas[i]).second, Refusal::Kind::input,
         "input: duplicate or out-of-grid feedback atlas");
    need(c.n_targets <= pfg::kMaxTargets, Refusal::Kind::input, "input: cell with more than 12 targets");
    std::set<u32> slots;
    for (u32 k = 0; k < c.n_targets; ++k)
      need(c.species[k] < L_.S && slots.insert(c.species[k]).second && std::isfinite(c.n[k]) && c.n[k] >= 0 &&
               pfg::finite3(c.drift[k]) && std::isfinite(c.s[k]) && c.s[k] >= 0,
           Refusal::Kind::input, "input: invalid or duplicate feedback target");
    r.slots += c.n_targets;
  }
  if (cfg_.feedback) {
#ifdef PF_GAS_WITH_DEVICE
    std::vector<CellTargetsAbi> dev_cells;
    FeedbackReceipt dr;
    if (dev_) {
      dev_cells.assign(cells.begin(), cells.end());
      dr = dev_->feedback(atlas, dev_cells, cfg_.grid.volume, open_);
    }
#endif
    if (host()) {
      const State &s = open_ ? *cand_ : *acc_;
      for (std::size_t i = 0; i < cells.size(); ++i) {
        auto &c = candidate[i];
        for (u32 k = 0; k < c.n_targets; ++k) {
          const u32 sp = c.species[k];
          if (sp >= L_.S || !(species_[sp].flags & pfg::kSpFeedback))
            continue;
          ++r.eligible;
          const i64 e = pfg::find_key(s.keys.data(), s.keys.size(), pfg::entry_key(atlas[i], sp));
          if (e < 0)
            continue;
          int clamped = 0;
          const int result = pfg::feedback_slot(c.n[k], cfg_.grid.volume, species_[sp].mass_kg,
                                                s.value.data() + e * pfg::kEntryQ, c.drift[k], &c.s[k], &clamped);
          need(result >= 0, Refusal::Kind::exactness, "exactness: nonrealizable gas feedback");
          if (result > 0) {
            ++r.changed;
            r.clamped += u64(clamped);
          }
        }
      }
    }
#ifdef PF_GAS_WITH_DEVICE
    if (dev_) {
      if (!host()) {
        std::copy(dev_cells.begin(), dev_cells.end(), candidate.begin());
        r.eligible = dr.eligible;
        r.changed = dr.changed;
        r.clamped = dr.clamped;
      } else
        need(std::memcmp(dev_cells.data(), candidate.data(), cells.size() * sizeof(CellTargetsAbi)) == 0 &&
                 dr.changed == r.changed && dr.clamped == r.clamped && dr.eligible == r.eligible,
             Refusal::Kind::device, "device: feedback differs from the host reference");
    }
#endif
  }
  r.cells_after_sha256 = cells_sha(candidate);
  std::copy(candidate.begin(), candidate.end(), cells.begin());
  return r;
}

Owner::PreparedCommit Owner::prepare_commit(const std::string &candidate_sha256) {
  need(open_, Refusal::Kind::transaction, "transaction: commit without an open step");
  need(!publication_ready_, Refusal::Kind::transaction, "transaction: publication already prepared");
  invalidate_candidate_bindings();
  need(!failed_ && next_stage_ == 3 && cand_->counter[pfg::C_STATUS] == 0,
       Refusal::Kind::transaction, "transaction: commit requires two successful ordered stages");
  const std::string d = digest(*cand_, true);
  need(d == candidate_sha256, Refusal::Kind::transaction, "transaction: commit digest differs from the candidate");
  State &s = *cand_;
  if (host()) {   // compaction: zero entries are not part of the state
    std::vector<u64> keys;
    std::vector<i64> entry;
    std::vector<double> value;
    for (u64 i = 0; i < s.keys.size(); ++i) {
      const i64 *dd = s.entry.data() + i * pfg::kEntryQ * D;
      if (pfg::digits_zero(dd, int(pfg::kEntryQ * D)))
        continue;
      keys.push_back(s.keys[i]);
      entry.insert(entry.end(), dd, dd + pfg::kEntryQ * D);
      value.insert(value.end(), s.value.begin() + i64(i * pfg::kEntryQ),
                   s.value.begin() + i64((i + 1) * pfg::kEntryQ));
    }
    s.keys = std::move(keys);
    s.entry = std::move(entry);
    s.value = std::move(value);
  }
  PreparedCommit token;
  token.lifetime_ = lifetime_;
  token.generation_ = generation_;
  token.digest_ = d;
  publication_ready_ = true;
  return token;
}

bool Owner::can_publish(const PreparedCommit &token) const noexcept {
  return open_ && publication_ready_ && !failed_ && token.lifetime_ == lifetime_ &&
         token.generation_ == generation_;
}

bool Owner::try_publish(PreparedCommit &token) noexcept {
  if (!can_publish(token)) return false;
  // Publication performs only scalar/shared-token changes, buffer swaps and
  // destruction of the old host state. It allocates nothing and launches no work.
#ifdef PF_GAS_WITH_DEVICE
  if (dev_)
    dev_->commit();
#endif
  acc_ = std::move(cand_);
  open_ = false;
  publication_ready_ = false;
  token.lifetime_.reset();
  return true;
}

void Owner::commit(const std::string &candidate_sha256) {
  auto token = prepare_commit(candidate_sha256);
  if (!try_publish(token))
    refuse(Refusal::Kind::transaction, "transaction: prepared publication failed");
}

void Owner::abandon() {
  cand_.reset();
  open_ = false;
  publication_ready_ = false;
  failed_ = false;
  next_stage_ = 1;
  partial_batch_ready_ = false;
  partial_batch_inputs_.clear();
#ifdef PF_GAS_WITH_DEVICE
  if (dev_)
    dev_->abandon();
#endif
}

// State digest: host leaves (host execution) or device leaves.
std::string Owner::digest(const State &s, bool candidate) const {
  std::vector<std::array<pf7_u32, 8>> el, gl;
#ifdef PF_GAS_WITH_DEVICE
  if (!host()) {
    dev_->leaves(candidate, el, gl);
    return state_digest(cfg_sha_, s.step, s.counter, el, gl);
  }
#endif
  (void)candidate;
  host_leaves(s, L_.total(), el, gl);
  return state_digest(cfg_sha_, s.step, s.counter, el, gl);
}
std::string Owner::accepted_sha256() const { return digest(*acc_, false); }
std::string Owner::candidate_sha256() const {
  need(open_, Refusal::Kind::transaction, "transaction: no candidate");
  return digest(*cand_, true);
}

void Owner::verify_device(const char *what, const std::string &host_digest) const {
#ifdef PF_GAS_WITH_DEVICE
  if (dev_ && cfg_.execution == Execution::device_verified) {
    std::vector<std::array<pf7_u32, 8>> el, gl;
    dev_->leaves(true, el, gl);
    const std::string d = state_digest(cfg_sha_, cand_->step, cand_->counter, el, gl);
    need(host_digest == d, Refusal::Kind::device,
         std::string("device: ") + what + " state differs from the host reference");
  }
#else
  (void)what;
  (void)host_digest;
#endif
}

// Host view of a state (device execution: downloaded).
const Owner::State &Owner::view(bool candidate, State &tmp) const {
  const State &s = candidate ? *cand_ : *acc_;
  if (host())
    return s;
#ifdef PF_GAS_WITH_DEVICE
  tmp.step = s.step;
  tmp.counter = s.counter;
  dev_->download(candidate, tmp.keys, tmp.entry, tmp.global);
  normalise_state(tmp, L_.total());
#endif
  return tmp;
}

std::array<std::int64_t, pfg::kDigits> Owner::global_digits(std::uint32_t g, bool candidate) const {
  need(g < L_.total(), Refusal::Kind::input, "input: global index");
  State tmp;
  const State &s = view(candidate, tmp);
  std::array<i64, D> x;
  std::memcpy(x.data(), s.global.data() + u64(g) * D, sizeof x);
  pf7::normalise(x.data());
  return x;
}

std::vector<EntryValue> Owner::entries(bool candidate) const {
  State tmp;
  const State &s = view(candidate, tmp);
  std::vector<EntryValue> v;
  for (u64 i = 0; i < s.keys.size(); ++i) {
    const i64 *dd = s.entry.data() + i * pfg::kEntryQ * D;
    if (pfg::digits_zero(dd, int(pfg::kEntryQ * D)))
      continue;
    EntryValue e;
    e.cell = u32(s.keys[i] >> 32);
    e.species = u32(s.keys[i] & 0xffffffffu);
    for (u32 q = 0; q < pfg::kEntryQ; ++q)
      e.value[q] = s.value[i * pfg::kEntryQ + q];
    v.push_back(e);
  }
  return v;
}

Exact Owner::entry_total(std::uint32_t q, bool candidate) const {
  State tmp;
  const State &s = view(candidate, tmp);
  Exact t;
  for (u64 i = 0; i < s.keys.size(); ++i)
    t.add_digits(s.entry.data() + (i * pfg::kEntryQ + q) * D);
  return t;
}

Ledgers Owner::ledgers(bool candidate) const {
  State tmp;
  const State &s = view(candidate, tmp);
  Ledgers l;
  l.state_sha256 = candidate ? candidate_sha256() : accepted_sha256();
  l.counters = s.counter;
  auto ev = [&](u32 g) { return exact_value(s.global.data() + u64(g) * D); };
  for (u32 g = 0; g < pfg::G_FIXED; ++g)
    l.fixed.push_back({fixed_name(g), ev(g)});
  static const char *fq[] = {"count", "E", "Px", "Py", "Pz", "internal"};
  static const char *kinds[] = {"debit", "release", "birth", "sink", "transaction"};
  static const char *cq[] = {"E", "Px", "Py", "Pz"};
  static const char *wq[] = {"E", "Px", "Py", "Pz", "U"};
  static const char *pq[] = {"W", "E_gas", "Px_gas", "Py_gas", "Pz_gas"};
  for (u32 sp = 0; sp < L_.S; ++sp) {
    const std::string &n = cfg_.species[sp].state;
    for (u32 k = 0; k < pfg::kFlowKinds; ++k)
      for (u32 q = 0; q < 6; ++q) {
        const u32 g = L_.flow(k, sp, q);
        if (!pfg::digits_zero(s.global.data() + u64(g) * D, int(D)))
          l.per_species.push_back({std::string("flow/") + kinds[k] + "/" + n + "/" + fq[q], ev(g)});
      }
    for (u32 q = 0; q < 4; ++q) {
      const u32 g = L_.chg(sp, q);
      if (!pfg::digits_zero(s.global.data() + u64(g) * D, int(D)))
        l.per_species.push_back({std::string("charged_substitution/") + n + "/" + cq[q], ev(g)});
    }
    for (u32 q = 0; q < 5; ++q) {
      const u32 g = L_.wall(sp, q);
      if (!pfg::digits_zero(s.global.data() + u64(g) * D, int(D)))
        l.per_species.push_back({std::string("wall_share/") + n + "/" + wq[q], ev(g)});
    }
  }
  for (u32 p = 0; p < L_.NP; ++p)
    for (u32 q = 0; q < 5; ++q) {
      const u32 g = L_.proc(p, q);
      if (!pfg::digits_zero(s.global.data() + u64(g) * D, int(D)))
        l.per_process.push_back({"process/" + std::to_string(p) + "/" + pq[q], ev(g)});
    }
  for (u64 i = 0; i < s.keys.size(); ++i) {
    const i64 *dd = s.entry.data() + i * pfg::kEntryQ * D;
    if (pfg::digits_zero(dd, int(pfg::kEntryQ * D)))
      continue;
    EntryValue e;
    e.cell = u32(s.keys[i] >> 32);
    e.species = u32(s.keys[i] & 0xffffffffu);
    for (u32 q = 0; q < pfg::kEntryQ; ++q)
      e.value[q] = s.value[i * pfg::kEntryQ + q];
    l.entries.push_back(e);
  }
  for (u32 q = 0; q < pfg::kEntryQ; ++q) {
    Exact t;
    for (u64 i = 0; i < s.keys.size(); ++i)
      t.add_digits(s.entry.data() + (i * pfg::kEntryQ + q) * D);
    const auto x = t.normalised();
    l.entry_totals[q] = exact_value(x.data());
  }
  return l;
}

std::string Owner::device_json() const {
#ifdef PF_GAS_WITH_DEVICE
  if (dev_)
    return dev_->session_json();
#endif
  return "null";
}
double Owner::device_kernel_ms() const {
#ifdef PF_GAS_WITH_DEVICE
  if (dev_)
    return dev_->kernel_ms_total();
#endif
  return 0;
}

// ------------------------------------------------------------------ codec
namespace {
const char kMagic[8] = {'P', 'F', 'G', 'A', 'S', 'V', '1', 0};
struct Writer {
  std::vector<std::uint8_t> b;
  void raw(const void *p, std::size_t n) {
    const auto *c = static_cast<const std::uint8_t *>(p);
    b.insert(b.end(), c, c + n);
  }
  void u64v(u64 v) { raw(&v, 8); }
  void text(const std::string &s) {
    u64v(s.size());
    raw(s.data(), s.size());
  }
};
struct Reader {
  std::span<const std::uint8_t> b;
  std::size_t p = 0;
  void raw(void *o, std::size_t n) {
    need(p + n <= b.size() && p + n >= p, Refusal::Kind::restore, "restore: truncated checkpoint");
    std::memcpy(o, b.data() + p, n);
    p += n;
  }
  u64 u64v() {
    u64 v;
    raw(&v, 8);
    return v;
  }
  std::string text() {
    const u64 n = u64v();
    need(n <= 4096, Refusal::Kind::restore, "restore: text field");
    std::string s(n, '\0');
    raw(s.data(), n);
    return s;
  }
};
bool normalised_digits(const i64 *d) {
  for (std::size_t k = 0; k + 1 < D; ++k)
    if (d[k] < 0 || d[k] > 0xffffffffll)
      return false;
  return true;
}
}  // namespace

std::vector<std::uint8_t> Owner::encode() const {
  need(!open_, Refusal::Kind::transaction, "transaction: checkpoint with an open step");
  State tmp;
  const State &s = view(false, tmp);
  Writer w;
  w.raw(kMagic, 8);
  w.u64v(1);   // codec version
  w.text(cfg_sha_);
  w.u64v(s.step);
  for (auto c : s.counter)
    w.u64v(c);
  std::vector<u64> nz;
  for (u64 i = 0; i < s.keys.size(); ++i)
    if (!pfg::digits_zero(s.entry.data() + i * pfg::kEntryQ * D, int(pfg::kEntryQ * D)))
      nz.push_back(i);
  w.u64v(nz.size());
  for (u64 i : nz) {
    w.u64v(s.keys[i]);
    w.raw(s.entry.data() + i * pfg::kEntryQ * D, sizeof(i64) * pfg::kEntryQ * D);
  }
  std::vector<u32> gz;
  for (u32 g = 0; g < L_.total(); ++g)
    if (!pfg::digits_zero(s.global.data() + u64(g) * D, int(D)))
      gz.push_back(g);
  w.u64v(gz.size());
  for (u32 g : gz) {
    w.u64v(g);
    w.raw(s.global.data() + u64(g) * D, sizeof(i64) * D);
  }
  w.text(accepted_sha256());
  return w.b;
}

std::unique_ptr<Owner> Owner::restore(std::span<const std::uint8_t> bytes, const Configuration &cfg,
                                      const std::string &expected) {
  auto o = std::make_unique<Owner>(cfg);
  Reader r{bytes};
  char magic[8];
  r.raw(magic, 8);
  need(std::memcmp(magic, kMagic, 8) == 0, Refusal::Kind::restore, "restore: not a gas owner v1 checkpoint");
  need(r.u64v() == 1, Refusal::Kind::restore, "restore: codec version");
  need(r.text() == o->cfg_sha_, Refusal::Kind::restore, "restore: configuration identity differs");
  State &s = *o->acc_;
  s.global.assign(u64(o->L_.total()) * D, 0);
  s.step = r.u64v();
  for (auto &c : s.counter)
    c = r.u64v();
  need(s.counter[pfg::C_STATUS] == 0, Refusal::Kind::restore, "restore: status bits set");
  const u64 ne = r.u64v();
  need(ne <= cfg.max_entries, Refusal::Kind::restore, "restore: entry count");
  s.keys.resize(ne);
  s.entry.assign(ne * pfg::kEntryQ * D, 0);
  for (u64 i = 0; i < ne; ++i) {
    s.keys[i] = r.u64v();
    const u32 sp = u32(s.keys[i] & 0xffffffffu);
    need(sp < o->L_.S && (o->species_[sp].flags & pfg::kSpGas), Refusal::Kind::restore, "restore: entry species");
    need(i == 0 || s.keys[i] > s.keys[i - 1], Refusal::Kind::restore, "restore: entry order");
    r.raw(s.entry.data() + i * pfg::kEntryQ * D, sizeof(i64) * pfg::kEntryQ * D);
    for (u32 q = 0; q < pfg::kEntryQ; ++q)
      need(normalised_digits(s.entry.data() + (i * pfg::kEntryQ + q) * D), Refusal::Kind::restore,
           "restore: entry digits not normalised");
    need(!pfg::digits_zero(s.entry.data() + i * pfg::kEntryQ * D, int(pfg::kEntryQ * D)), Refusal::Kind::restore,
         "restore: zero entry");
  }
  const u64 ng = r.u64v();
  need(ng <= o->L_.total(), Refusal::Kind::restore, "restore: global count");
  i64 last = -1;
  for (u64 i = 0; i < ng; ++i) {
    const u64 g = r.u64v();
    need(g < o->L_.total() && i64(g) > last, Refusal::Kind::restore, "restore: global index");
    last = i64(g);
    r.raw(s.global.data() + g * D, sizeof(i64) * D);
    need(normalised_digits(s.global.data() + g * D) && !pfg::digits_zero(s.global.data() + g * D, int(D)),
         Refusal::Kind::restore, "restore: global digits");
  }
  const std::string stored = r.text();
  need(r.p == bytes.size(), Refusal::Kind::restore, "restore: trailing bytes");
  normalise_state(s, o->L_.total());
#ifdef PF_GAS_WITH_DEVICE
  if (o->dev_) {
    o->dev_->ensure_keys(s.keys);
    o->dev_->load(s.keys, s.entry, s.global);
  }
#endif
  if (!o->host()) {
    s.entry.clear();
    s.value.clear();
    s.global.clear();
  }
  const std::string d = o->accepted_sha256();
  need(d == stored, Refusal::Kind::restore, "restore: state digest differs from the stored digest (tampered)");
  need(d == expected, Refusal::Kind::restore, "restore: state digest differs from the expected boundary");
  if (o->dev_ && o->host())
    need(o->digest(s, false) == d, Refusal::Kind::device, "device: restored state differs from the host reference");
  return o;
}

// ------------------------------------------------------------------ census and movers
namespace {
void census_configuration(const Configuration &cfg) {
  need(!cfg.species.empty(), Refusal::Kind::configuration, "configuration: census species");
  u64 nc = 1;
  for (int a = 0; a < 3; ++a) {
    need(cfg.grid.n[a] && std::isfinite(cfg.grid.low[a]) && std::isfinite(cfg.grid.high[a]) &&
             cfg.grid.high[a] > cfg.grid.low[a] && nc <= u64(pfg::kNone) / cfg.grid.n[a],
         Refusal::Kind::configuration, "configuration: census grid");
    nc *= cfg.grid.n[a];
  }
  const double volume = ((cfg.grid.high[0]-cfg.grid.low[0])/cfg.grid.n[0]) *
                        ((cfg.grid.high[1]-cfg.grid.low[1])/cfg.grid.n[1]) *
                        ((cfg.grid.high[2]-cfg.grid.low[2])/cfg.grid.n[2]);
  need(std::isfinite(volume) && volume > 0 && cfg.grid.volume == volume,
       Refusal::Kind::configuration, "configuration: census volume");
  for (const auto &s : cfg.species)
    need(std::isfinite(s.mass_kg) && s.mass_kg > 0,
         Refusal::Kind::configuration, "configuration: census mass");
}
std::map<u64, const CensusRow *> census_inventory(const Configuration &cfg, std::span<const CensusRow> rows) {
  std::map<u64, const CensusRow *> inventory;
  for (const auto &r : rows) {
    need(r.species < cfg.species.size() && inventory.emplace(r.native_id, &r).second &&
             std::isfinite(r.weight) && r.weight > 0 && pfg::finite3(r.position_m) &&
             pfg::finite3(r.proper_velocity_m_s),
         Refusal::Kind::input, "input: invalid or duplicate census row");
    need(cfg.grid.locate(r.position_m) != pfg::kNone,
         Refusal::Kind::input, "input: census carrier outside the collision grid");
  }
  return inventory;
}
} // namespace
const CensusCell *Census::find(std::uint32_t cell, std::uint32_t species) const {
  auto it = std::lower_bound(cells.begin(), cells.end(), std::make_pair(cell, species),
                             [](const CensusCell &c, const std::pair<u32, u32> &k) {
                               return c.cell != k.first ? c.cell < k.first : c.species < k.second;
                             });
  return (it != cells.end() && it->cell == cell && it->species == species) ? &*it : nullptr;
}

Census census(const Configuration &cfg, std::span<const CensusRow> rows) {
  census_configuration(cfg);
  census_inventory(cfg, rows);
  Census c;
  c.species_totals.assign(cfg.species.size(), {});
  std::map<std::pair<u32, u32>, CensusCell> m;
  for (const auto &r : rows) {
    need(r.species < cfg.species.size(), Refusal::Kind::input, "input: census species");
    if (cfg.species[r.species].charge != 0)
      continue;
    const u32 cell = cfg.grid.locate(r.position_m);
    need(cell != pfg::kNone, Refusal::Kind::input, "input: census carrier outside the collision grid");
    auto &x = m[{cell, r.species}];
    x.cell = cell;
    x.species = r.species;
    const double mass = cfg.species[r.species].mass_kg;
    const double K = pfg::kinetic(mass, r.proper_velocity_m_s);
    need(std::isfinite(K), Refusal::Kind::exactness, "exactness: nonfinite census kinetic energy");
    x.count.add_product(r.weight, 1.0);
    x.kinetic.add_product(r.weight, K);
    auto &t = c.species_totals[r.species];
    t[0].add_product(r.weight, 1.0);
    t[1].add_product(r.weight, K);
    for (int a = 0; a < 3; ++a) {
      x.momentum[a].add_triple(r.weight, mass, r.proper_velocity_m_s[a]);
      t[2 + a].add_triple(r.weight, mass, r.proper_velocity_m_s[a]);
    }
  }
  for (auto &[k, v] : m)
    c.cells.push_back(std::move(v));
  return c;
}

std::vector<MoverRecord> movers(const Configuration &cfg, std::span<const CensusRow> before,
                                std::span<const CensusRow> after, std::span<const std::uint64_t> sunk) {
  const Census cb = census(cfg, before);
  const auto before_inventory = census_inventory(cfg, before);
  const auto after_inventory = census_inventory(cfg, after);
  const std::set<u64> sunk_set(sunk.begin(), sunk.end());
  need(sunk_set.size() == sunk.size(), Refusal::Kind::input, "input: duplicate sink ID");
  for (u64 id : sunk_set)
    need(before_inventory.count(id) && !after_inventory.count(id), Refusal::Kind::input,
         "input: sink absent from before inventory or also survives");
  for (const auto &[id, row] : after_inventory) {
    const auto b = before_inventory.find(id);
    need(b != before_inventory.end() && row->species == b->second->species &&
             std::bit_cast<u64>(row->weight) == std::bit_cast<u64>(b->second->weight),
         Refusal::Kind::input, "input: motion changes particle inventory, species or weight");
  }
  need(before_inventory.size() == after_inventory.size() + sunk_set.size(),
       Refusal::Kind::input, "input: incomplete motion inventory");
  std::unordered_map<u64, u32> after_cell;
  after_cell.reserve(after.size());
  for (const auto &r : after)
    if (r.species < cfg.species.size() && cfg.species[r.species].charge == 0)
      after_cell[r.native_id] = cfg.grid.locate(r.position_m);
  std::vector<MoverRecord> out;
  std::vector<const CensusRow *> order;
  for (const auto &r : before)
    if (r.species < cfg.species.size() && cfg.species[r.species].charge == 0)
      order.push_back(&r);
  std::sort(order.begin(), order.end(), [](const CensusRow *a, const CensusRow *b) { return a->native_id < b->native_id; });
  for (const auto *r : order) {
    const u32 from = cfg.grid.locate(r->position_m);
    u32 to;
    if (sunk_set.count(r->native_id))
      to = pfg::kNone;
    else {
      auto it = after_cell.find(r->native_id);
      need(it != after_cell.end(), Refusal::Kind::input,
           "input: carrier " + std::to_string(r->native_id) + " neither moved nor sunk");
      to = it->second;
      need(to != pfg::kNone, Refusal::Kind::input, "input: carrier left the grid without a sink");
      if (to == from)
        continue;
    }
    const CensusCell *cc = cb.find(from, r->species);
    MoverRecord m{};
    m.native_id = r->native_id;
    m.species = r->species;
    m.from_cell = from;
    m.to_cell = to;
    m.weight = r->weight;
    m.cell_count = cc->count.value();
    out.push_back(m);
  }
  return out;
}

}  // namespace pf_gas_v1
