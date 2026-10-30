#include "raikou/mfem_pic_arena.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <openssl/bn.h>
#include <openssl/sha.h>
#include <sstream>
#include <system_error>
#include <unistd.h>
#ifndef PIC_REACTION_NATIVE_SOURCE_SHA256
#error "native source identity must be supplied"
#endif
namespace raikou::mfem_backend {
namespace {
std::string sha_bytes(const std::string &bytes) {
  unsigned char d[SHA256_DIGEST_LENGTH];
  SHA256(reinterpret_cast<const unsigned char *>(bytes.data()), bytes.size(),
         d);
  std::ostringstream o;
  for (auto x : d)
    o << std::hex << std::setfill('0') << std::setw(2) << unsigned(x);
  return o.str();
}
void require(bool condition, const char *message) {
  if (!condition)
    throw std::invalid_argument(message);
}
bool particle_equal(const GapParticle &a, const GapParticle &b) noexcept {
  return a.id == b.id && a.position_m == b.position_m &&
         a.proper_velocity_m_s == b.proper_velocity_m_s &&
         a.charge_C == b.charge_C && a.mass_kg == b.mass_kg &&
         a.weight == b.weight;
}
void native_particle_guard(const GapParticle &p) {
  double u2 = 0;
  for (double x : p.position_m)
    require(std::isfinite(x), "nonfinite native product position");
  for (double u : p.proper_velocity_m_s) {
    require(std::isfinite(u), "nonfinite native product velocity");
    u2 += u * u;
  }
  require(std::isfinite(u2) && std::isfinite(p.charge_C) &&
              std::isfinite(p.mass_kg) && p.mass_kg > 0 &&
              std::isfinite(p.weight) && p.weight > 0 &&
              std::isfinite(p.charge_C * p.weight) &&
              std::isfinite(p.charge_C / p.mass_kg) &&
              std::isfinite(pic_kinetic_energy_J(p) * p.weight),
          "native product guard");
}
bool options_equal(const GapPicOptions &a, const GapPicOptions &b) noexcept {
  const auto foil = [](const GapFoil &x, const GapFoil &y) {
    return x.left_m == y.left_m && x.right_m == y.right_m &&
           x.bottom_m == y.bottom_m && x.top_m == y.top_m;
  };
  return a.nx == b.nx && a.ny == b.ny && a.width_m == b.width_m &&
         a.gap_m == b.gap_m && a.depth_m == b.depth_m &&
         a.lower_potential_V == b.lower_potential_V &&
         a.upper_potential_V == b.upper_potential_V &&
         a.magnetic_field_T == b.magnetic_field_T &&
         a.distributed_fields == b.distributed_fields &&
         a.finite_foils == b.finite_foils &&
         a.finite_axial_extent == b.finite_axial_extent &&
         a.short_z_m == b.short_z_m && a.feed_z_m == b.feed_z_m &&
         foil(a.lower_foil, b.lower_foil) && foil(a.upper_foil, b.upper_foil);
}
struct TemporaryNativeFile {
  std::string path;
  explicit TemporaryNativeFile(const std::string &near_path) {
    std::string pattern = near_path + ".native.XXXXXX";
    std::vector<char> name(pattern.begin(), pattern.end());
    name.push_back('\0');
    int fd = mkstemp(name.data());
    if (fd < 0)
      throw std::system_error(errno, std::generic_category(),
                              "native temporary create");
    ::close(fd);
    path = name.data();
  }
  ~TemporaryNativeFile() {
    unlink(path.c_str());
    unlink((path + ".tmp").c_str());
  }
};
std::string hex_bytes(const std::string &bytes) {
  static constexpr char digits[] = "0123456789abcdef";
  std::string out;
  out.reserve(2 * bytes.size());
  for (unsigned char x : bytes) {
    out.push_back(digits[x >> 4]);
    out.push_back(digits[x & 15]);
  }
  return out;
}
std::string unhex(const std::string &text) {
  require(text.size() % 2 == 0 && text.size() <= 32 * 1024 * 1024,
          "native bundle byte bound");
  std::string bytes;
  bytes.reserve(text.size() / 2);
  auto digit = [](char c) -> unsigned {
    if (c >= '0' && c <= '9')
      return c - '0';
    if (c >= 'a' && c <= 'f')
      return c - 'a' + 10;
    throw std::invalid_argument("native hex encoding");
  };
  for (std::size_t i = 0; i < text.size(); i += 2)
    bytes.push_back(char(16 * digit(text[i]) + digit(text[i + 1])));
  return bytes;
}
void atomic_json(const ArenaJson &document, const std::string &path) {
  const auto bytes = document.dump() + "\n";
  require(bytes.size() <= 64 * 1024 * 1024, "native bundle size bound");
  std::string pattern = path + ".tmp.XXXXXX";
  std::vector<char> name(pattern.begin(), pattern.end());
  name.push_back('\0');
  int fd = mkstemp(name.data());
  if (fd < 0)
    throw std::system_error(errno, std::generic_category(),
                            "bundle temporary create");
  try {
    std::size_t offset = 0;
    while (offset < bytes.size()) {
      auto n = write(fd, bytes.data() + offset, bytes.size() - offset);
      if (n < 0 && errno == EINTR)
        continue;
      if (n <= 0)
        throw std::system_error(errno, std::generic_category(), "bundle write");
      offset += std::size_t(n);
    }
    if (fsync(fd))
      throw std::system_error(errno, std::generic_category(), "bundle fsync");
    ::close(fd);
    fd = -1;
    if (rename(name.data(), path.c_str()))
      throw std::system_error(errno, std::generic_category(), "bundle rename");
    auto directory = std::filesystem::path(path).parent_path();
    if (directory.empty())
      directory = ".";
    int d = open(directory.c_str(), O_RDONLY | O_DIRECTORY);
    if (d < 0)
      throw std::system_error(errno, std::generic_category(),
                              "bundle directory open");
    int status = fsync(d), saved = errno;
    ::close(d);
    if (status)
      throw std::system_error(saved, std::generic_category(),
                              "bundle directory fsync");
  } catch (...) {
    if (fd >= 0)
      ::close(fd);
    unlink(name.data());
    throw;
  }
}

} // namespace
struct PicArena::NativeState {
  GapPicOptions options;
  std::vector<GapParticle> particles;
  std::vector<double> horizontal, vertical;
  mfem::Vector potential;
  GapPicDiagnostics ledger;
  std::uint64_t next_id = 0;
  bool fields_current = false;
};
PicArena::NativeState PicArena::capture(const GapPic &p) {
  NativeState s;
  s.options = p.options_;
  s.particles = p.particles_;
  s.horizontal = p.horizontal_transfers_;
  s.vertical = p.vertical_transfers_;
  s.potential = p.potential_;
  s.ledger = p.ledger_;
  s.next_id = p.next_particle_id_;
  s.fields_current = p.fields_current_;
  return s;
}
void PicArena::swap_native(GapPic &p, NativeState &s) noexcept {
  p.particles_.swap(s.particles);
  p.horizontal_transfers_.swap(s.horizontal);
  p.vertical_transfers_.swap(s.vertical);
  p.potential_.Swap(s.potential);
  std::swap(p.ledger_, s.ledger);
  std::swap(p.next_particle_id_, s.next_id);
  std::swap(p.fields_current_, s.fields_current);
}
bool PicArena::matches(const GapPic &p, const NativeState &s) noexcept {
  static_assert(sizeof(GapPicDiagnostics) == 200 && sizeof(GapParticle) == 80 &&
                sizeof(GapPicOptions) == 168);
  if (!options_equal(p.options_, s.options) ||
      p.particles_.size() != s.particles.size() ||
      p.potential_.Size() != s.potential.Size() ||
      p.next_particle_id_ != s.next_id ||
      p.fields_current_ != s.fields_current ||
      p.horizontal_transfers_ != s.horizontal ||
      p.vertical_transfers_ != s.vertical ||
      std::memcmp(&p.ledger_, &s.ledger, sizeof(s.ledger)))
    return false;
  for (std::size_t i = 0; i < s.particles.size(); ++i)
    if (!particle_equal(p.particles_[i], s.particles[i]))
      return false;
  for (int i = 0; i < s.potential.Size(); ++i)
    if (p.potential_[i] != s.potential[i])
      return false;
  return true;
}
std::string PicArena::fingerprint(const GapPic &p) {
  std::string data;
  auto add = [&](const auto &x) {
    data.append(reinterpret_cast<const char *>(&x), sizeof(x));
  };
  add(p.options_.nx);
  add(p.options_.ny);
  for (double x : {p.options_.width_m, p.options_.gap_m, p.options_.depth_m,
                   p.options_.lower_potential_V, p.options_.upper_potential_V,
                   p.options_.short_z_m, p.options_.feed_z_m})
    add(x);
  for (double x : p.options_.magnetic_field_T)
    add(x);
  for (bool x : {p.options_.distributed_fields, p.options_.finite_foils,
                 p.options_.finite_axial_extent})
    add(x);
  for (auto f : {p.options_.lower_foil, p.options_.upper_foil})
    for (double x : {f.left_m, f.right_m, f.bottom_m, f.top_m})
      add(x);
  add(p.ledger_);
  add(p.next_particle_id_);
  add(p.fields_current_);
  std::uint64_t count = p.particles_.size();
  add(count);
  for (const auto &particle : p.particles_)
    add(particle);
  for (int i = 0; i < p.potential_.Size(); ++i)
    add(p.potential_[i]);
  for (double x : p.horizontal_transfers_)
    add(x);
  for (double x : p.vertical_transfers_)
    add(x);
  return sha_bytes(data);
}
const char *PicArena::source_identity() noexcept {
  return PIC_REACTION_NATIVE_SOURCE_SHA256;
}
namespace {
using Big = std::unique_ptr<BIGNUM, decltype(&BN_free)>;
Big big() {
  Big n(BN_new(), BN_free);
  if (!n)
    throw std::bad_alloc();
  return n;
}
Big exact(double x) {
  require(std::isfinite(x), "nonfinite exact count");
  std::uint64_t bits = 0;
  std::memcpy(&bits, &x, 8);
  unsigned e = (bits >> 52) & 2047;
  std::uint64_t m = bits & ((1ULL << 52) - 1);
  if (e)
    m |= 1ULL << 52;
  auto n = big();
  static_assert(sizeof(BN_ULONG) >= 8);
  require(BN_set_word(n.get(), m), "count allocation");
  if (e)
    require(BN_lshift(n.get(), n.get(), int(e) - 1), "count shift");
  BN_set_negative(n.get(), (bits >> 63) && m);
  return n;
}
Big exact(ArenaCount n) {
  auto a = exact(n.hi), b = exact(n.lo);
  require(BN_add(a.get(), a.get(), b.get()), "count sum");
  return a;
}
void count_guard(ArenaCount n) {
  require(electrode_heat::canonical(n), "noncanonical heat-owner Count");
  auto a = exact(n);
  require(!BN_is_negative(a.get()), "negative exact Count");
}
bool zero(ArenaCount n) {
  auto a = exact(n);
  return BN_is_zero(a.get());
}
ArenaJson count_json(ArenaCount n) { return ArenaJson::array({n.hi, n.lo}); }
ArenaCount read_count(const ArenaJson &j) {
  require(j.is_array() && j.size() == 2, "count schema");
  ArenaCount n{j.at(0).get<double>(), j.at(1).get<double>()};
  count_guard(n);
  return n;
}
bool digest(const std::string &s) {
  return s.size() == 64 && std::all_of(s.begin(), s.end(), [](char c) {
           return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
         });
}
void nonempty(const std::string &s) {
  require(!s.empty() && s.size() < 1024, "bounded explicit identity required");
}
long double kinetic(const GapParticle &p) {
  long double u2 = 0;
  for (double x : p.proper_velocity_m_s)
    u2 += static_cast<long double>(x) * x;
  constexpr long double c = 299792458.L;
  return p.mass_kg * u2 / (std::sqrt(1 + u2 / (c * c)) + 1);
}
ArenaJson packet_json(const ArenaPacket &p) {
  return {{"id", p.id},
          {"state", p.state},
          {"origin", p.origin},
          {"last_event", p.last_event},
          {"count", count_json(p.count)},
          {"native_ids", p.native_ids}};
}
ArenaJson metadata(const ArenaSnapshot &s) {
  ArenaJson p = ArenaJson::array();
  for (const auto &[id, x] : s.packets) {
    (void)id;
    p.push_back(packet_json(x));
  }
  return {{"packets", p},
          {"ever_packet_ids", s.ever_packet_ids},
          {"event_ids", s.event_ids},
          {"journal", s.journal},
          {"sinks", s.sinks},
          {"chain", s.chain},
          {"revision", s.revision},
          {"accepted_steps", s.accepted_steps},
          {"previous_dt_s", s.previous_dt_s}};
}
ArenaSnapshot read_metadata(const ArenaJson &j) {
  require(j.size() == 9, "metadata schema");
  ArenaSnapshot s;
  s.ever_packet_ids = j.at("ever_packet_ids").get<std::set<std::string>>();
  s.event_ids = j.at("event_ids").get<std::set<std::string>>();
  s.journal = j.at("journal");
  s.sinks = j.at("sinks");
  s.chain = j.at("chain");
  s.revision = j.at("revision");
  s.accepted_steps = j.at("accepted_steps");
  s.previous_dt_s = j.at("previous_dt_s");
  for (const auto &x : j.at("packets")) {
    require(x.size() == 6, "packet schema");
    ArenaPacket p{x.at("id"),
                  x.at("state"),
                  x.at("origin"),
                  x.at("last_event"),
                  read_count(x.at("count")),
                  x.at("native_ids").get<std::vector<std::uint64_t>>()};
    require(s.packets.emplace(p.id, p).second, "duplicate logical packet");
  }
  return s;
}
void journal_append(ArenaSnapshot &s, ArenaJson audit) {
  audit["previous_chain"] = s.chain;
  s.chain = sha_bytes(audit.dump());
  s.journal.push_back({{"payload", audit}, {"sha256", s.chain}});
  ++s.revision;
}
void journal_validate(const ArenaSnapshot &s, const std::string &genesis) {
  require(s.journal.is_array() && s.sinks.is_array() &&
              s.revision == s.journal.size(),
          "journal dimensions");
  std::string chain = genesis;
  std::set<std::string> events;
  std::uint64_t steps = 0;
  double previous = 0;
  for (const auto &x : s.journal) {
    require(x.size() == 2 && x.at("payload").at("previous_chain") == chain &&
                x.at("sha256") == sha_bytes(x.at("payload").dump()),
            "journal chain");
    chain = x.at("sha256");
    auto type = x.at("payload").at("kind").get<std::string>();
    if (type == "reaction") {
      auto id = x.at("payload").at("event_id").get<std::string>();
      require(events.insert(id).second, "duplicate accepted event");
    } else if (type == "native_transport") {
      ++steps;
      previous = x.at("payload").at("dt_s");
    } else
      throw std::invalid_argument("journal operation");
  }
  require(chain == s.chain && events == s.event_ids &&
              steps == s.accepted_steps && previous == s.previous_dt_s,
          "restart prefix mismatch");
}
} // namespace
std::string PicArena::sha256(const std::string &s) { return sha_bytes(s); }
bool PicArena::exact_count_equal(ArenaCount a, ArenaCount b) {
  auto x = exact(a), y = exact(b);
  return BN_cmp(x.get(), y.get()) == 0;
}
ArenaCount PicArena::debit(ArenaCount n, double w) {
  count_guard(n);
  require(std::isfinite(w) && w > 0, "positive event weight");
  auto before = exact(n), amount = exact(w);
  require(BN_cmp(before.get(), amount.get()) >= 0,
          "insufficient exact inventory");
  require(BN_sub(before.get(), before.get(), amount.get()),
          "count subtraction");
  auto out = electrode_heat::add(n, -w);
  count_guard(out);
  auto represented = exact(out);
  require(BN_cmp(before.get(), represented.get()) == 0,
          "exact remainder exceeds two-word Count; unchanged-state refusal");
  return out;
}
std::vector<double> PicArena::positive_carriers(ArenaCount n) {
  count_guard(n);
  auto remaining = exact(n);
  std::vector<double> out;
  while (!BN_is_zero(remaining.get())) {
    require(out.size() < 42, "native count carrier bound");
    int shift = std::max(0, BN_num_bits(remaining.get()) - 53);
    auto top = big();
    require(BN_rshift(top.get(), remaining.get(), shift), "carrier shift");
    double w = std::ldexp(double(BN_get_word(top.get())), shift - 1074);
    require(std::isfinite(w) && w > 0,
            "positive native carrier unrepresentable");
    auto word = exact(w);
    require(BN_cmp(word.get(), remaining.get()) <= 0 &&
                BN_sub(remaining.get(), remaining.get(), word.get()),
            "native carrier exact decomposition");
    out.push_back(w);
  }
  return out;
}
struct PreparedArena::Impl {
  PicArena *owner = nullptr;
  std::uint64_t scope = 0;
  std::shared_ptr<const ArenaSnapshot> before, after;
  PicArena::NativeState native_before, native_after, work;
  ArenaJson audit;
  enum class Status { Staged, Published, RolledBack };
  Status status = Status::Staged;
};
PreparedArena::PreparedArena(std::unique_ptr<Impl> p) : impl_(std::move(p)) {}
PreparedArena::PreparedArena(PreparedArena &&) noexcept = default;
PreparedArena &PreparedArena::operator=(PreparedArena &&) noexcept = default;
PreparedArena::~PreparedArena() = default;
const ArenaJson &PreparedArena::audit() const { return impl_->audit; }
ArenaJson PicArena::registry() const {
  ArenaJson states = ArenaJson::array(), channels = ArenaJson::array();
  for (const auto &[id, s] : states_)
    states.push_back({{"id", id},
                      {"H", s.hydrogen},
                      {"O", s.oxygen},
                      {"q_e", s.charge_e},
                      {"mass_kg", s.mass_kg},
                      {"internal_J", s.internal_J},
                      {"reference_J", s.reference_J}});
  for (const auto &[id, c] : channels_)
    channels.push_back({{"id", id},
                        {"primary", c.primary_state},
                        {"target", c.target_state},
                        {"products", c.product_states},
                        {"dataset", c.dataset_sha256},
                        {"kernel_manifest", c.kernel_manifest_sha256},
                        {"kernel_binding_file", c.kernel_binding_file_sha256},
                        {"model", c.model}});
  return {{"states", states}, {"channels", channels}};
}
PicArena::PicArena(GapPic &p, std::vector<ArenaStateSpec> states,
                   std::vector<ArenaChannel> channels,
                   std::vector<ArenaPacket> packets, std::string prefix)
    : pic_(&p), initial_prefix_(std::move(prefix)) {
  nonempty(initial_prefix_);
  require(!p.options_.distributed_fields && p.fields_current_,
          "arena currently requires solved serial CPU fields");
  for (auto s : states) {
    nonempty(s.id);
    require(s.hydrogen >= 0 && s.oxygen >= 0 && std::abs(s.charge_e) <= 8 &&
                std::isfinite(s.mass_kg) && s.mass_kg > 0 &&
                std::isfinite(s.internal_J) && s.internal_J >= 0 &&
                std::isfinite(s.reference_J),
            "state specification");
    require(states_.emplace(s.id, s).second, "duplicate state");
  }
  for (auto c : channels) {
    nonempty(c.id);
    require(states_.count(c.primary_state) && states_.count(c.target_state) &&
                states_.at(c.target_state).charge_e == 0 &&
                !c.product_states.empty(),
            "channel reactants");
    int h = -states_.at(c.primary_state).hydrogen -
            states_.at(c.target_state).hydrogen,
        o = -states_.at(c.primary_state).oxygen -
            states_.at(c.target_state).oxygen,
        q = -states_.at(c.primary_state).charge_e;
    for (auto id : c.product_states) {
      require(states_.count(id), "unknown product state");
      h += states_.at(id).hydrogen;
      o += states_.at(id).oxygen;
      q += states_.at(id).charge_e;
    }
    require(h == 0 && o == 0 && q == 0, "channel H/O/charge stoichiometry");
    if (c.model.empty())
      require(c.id.rfind("manufactured/", 0) == 0 && c.dataset_sha256.empty() &&
                  c.kernel_manifest_sha256.empty() &&
                  c.kernel_binding_file_sha256.empty(),
              "manufactured channel tagging");
    else
      require(c.id.rfind("manufactured/", 0) != 0 && digest(c.dataset_sha256) &&
                  digest(c.kernel_manifest_sha256) &&
                  digest(c.kernel_binding_file_sha256),
              "conditional channel immutable data/kernel binding");
    require(channels_.emplace(c.id, c).second, "duplicate channel");
  }
  registry_sha_ = sha_bytes(registry().dump());
  auto s = std::make_shared<ArenaSnapshot>();
  for (auto packet : packets) {
    nonempty(packet.id);
    nonempty(packet.origin);
    require(s->packets.emplace(packet.id, packet).second &&
                s->ever_packet_ids.insert(packet.id).second,
            "initial packet IDs");
  }
  s->chain = sha_bytes(initial_prefix_ + registry_sha_);
  accepted_ = s;
  validate_population(*s);
}
GapParticle PicArena::view(const std::string &id) const {
  const auto &packet = accepted_->packets.at(id);
  const GapParticle *first = nullptr;
  for (auto native : packet.native_ids) {
    auto it = std::find_if(pic_->particles_.begin(), pic_->particles_.end(),
                           [&](const auto &p) { return p.id == native; });
    require(it != pic_->particles_.end(), "native canonical carrier missing");
    if (first)
      require(first->position_m == it->position_m &&
                  first->proper_velocity_m_s == it->proper_velocity_m_s,
              "logical carrier trajectories differ; no averaging");
    else
      first = &*it;
  }
  require(first, "empty logical packet");
  auto p = *first;
  p.weight = packet.count.rounded_for_rate();
  return p;
}
std::vector<GapParticle> PicArena::rate_carriers(const std::string &id) const {
  validate_population(*accepted_);
  std::vector<GapParticle> out;
  for (auto native : accepted_->packets.at(id).native_ids) {
    auto it = std::find_if(pic_->particles_.begin(), pic_->particles_.end(),
                           [&](const auto &p) { return p.id == native; });
    require(it != pic_->particles_.end(), "rate carrier missing");
    out.push_back(*it);
  }
  return out;
}
void PicArena::validate_population(const ArenaSnapshot &s) const {
  std::set<std::uint64_t> all;
  for (const auto &[id, p] : s.packets) {
    require(id == p.id && s.ever_packet_ids.count(id) &&
                states_.count(p.state) && !p.native_ids.empty(),
            "logical packet identity/state");
    count_guard(p.count);
    require(!zero(p.count), "zero resident logical count");
    auto total = big();
    const GapParticle *first = nullptr;
    const auto &spec = states_.at(p.state);
    for (auto native : p.native_ids) {
      require(all.insert(native).second, "native carrier owned twice");
      auto it = std::find_if(pic_->particles_.begin(), pic_->particles_.end(),
                             [&](const auto &x) { return x.id == native; });
      require(it != pic_->particles_.end(),
              "missing canonical native particle");
      native_particle_guard(*it);
      require(it->charge_C == spec.charge_e * pic_e &&
                  it->mass_kg == spec.mass_kg,
              "native species binding");
      auto weight = exact(it->weight);
      require(BN_add(total.get(), total.get(), weight.get()),
              "carrier count accumulation");
      if (first)
        require(first->position_m == it->position_m &&
                    first->proper_velocity_m_s == it->proper_velocity_m_s,
                "divergent carrier positions/velocities");
      else
        first = &*it;
    }
    auto expected = exact(p.count);
    require(BN_cmp(total.get(), expected.get()) == 0,
            "native carriers do not exactly represent Count");
  }
  require(all.size() == pic_->particles_.size(), "unowned native particle");
}
PreparedArena PicArena::stage(std::shared_ptr<ArenaSnapshot> next,
                              NativeState &candidate, ArenaJson audit,
                              int fault) {
  auto token = std::make_unique<PreparedArena::Impl>();
  token->owner = this;
  token->scope = active_scope_;
  token->before = accepted_;
  token->native_before = capture(*pic_);
  bool swapped = false;
  try {
    swap_native(*pic_, candidate);
    swapped = true;
    if (fault == 1)
      throw std::runtime_error("manufactured fault after native stage");
    if (!pic_->fields_current_)
      pic_->solve_fields();
    validate_population(*next);
    if (fault == 2)
      throw std::runtime_error("manufactured fault after solve");
    token->native_after = capture(*pic_);
    token->work = token->native_after;
    swap_native(*pic_, candidate);
    swapped = false;
    require(matches(*pic_, token->native_before),
            "prepare failed to restore original native state");
    token->after = std::move(next);
    token->audit = std::move(audit);
    if (fault == 3)
      throw std::runtime_error("manufactured fault before return");
    return PreparedArena(std::move(token));
  } catch (...) {
    if (swapped)
      swap_native(*pic_, candidate);
    throw;
  }
}
bool PicArena::can_publish(const PreparedArena &t) const noexcept {
  return t.impl_ && t.impl_->owner == this && t.impl_->scope == active_scope_ &&
         t.impl_->status == PreparedArena::Impl::Status::Staged &&
         accepted_ == t.impl_->before && matches(*pic_, t.impl_->native_before);
}
bool PicArena::try_publish(PreparedArena &t) noexcept {
  if (!can_publish(t))
    return false;
  swap_native(*pic_, t.impl_->work);
  accepted_ = t.impl_->after;
  t.impl_->status = PreparedArena::Impl::Status::Published;
  return true;
}
bool PicArena::try_rollback(PreparedArena &t) noexcept {
  if (!t.impl_ || t.impl_->owner != this || t.impl_->scope != active_scope_ ||
      t.impl_->status != PreparedArena::Impl::Status::Published ||
      accepted_ != t.impl_->after || !matches(*pic_, t.impl_->native_after))
    return false;
  swap_native(*pic_, t.impl_->work);
  accepted_ = t.impl_->before;
  t.impl_->status = PreparedArena::Impl::Status::RolledBack;
  return true;
}

std::string PicArena::state_prefix() const {
  return sha_bytes(initial_prefix_ + registry_sha_ +
                   metadata(*accepted_).dump() + fingerprint(*pic_));
}
PreparedArena PicArena::prepare(const ArenaEvent &e, int fault) {
  require(!active_scope_ ||
              accepted_->revision - compound_base_revision_ < compound_limit_,
          "compound reaction resource cap reached");
  validate_population(*accepted_);
  require(pic_->fields_current_ && !pic_->options_.distributed_fields,
          "reaction requires solved serial synchronized fields");
  nonempty(e.id);
  require(!accepted_->event_ids.count(e.id) && e.primary != e.target &&
              accepted_->packets.count(e.primary) &&
              accepted_->packets.count(e.target),
          "duplicate event or invalid reactant identity");
  const auto prefix = state_prefix();
  const auto &phase = e.phase;
  require(phase.state_prefix == prefix &&
              phase.physical_stage_s == pic_->ledger_.time_s,
          "stale physical state/clock");
  nonempty(phase.mesh_prefix);
  nonempty(phase.clock_prefix);
  if (phase.kind == ArenaPhase::Kind::InitialFieldProof)
    require(phase.accepted_step == 0 && accepted_->accepted_steps == 0 &&
                phase.physical_stage_s == 0 && phase.previous_dt_s == 0 &&
                phase.next_dt_s == 0,
            "initial proof phase");
  else
    require(phase.accepted_step == accepted_->accepted_steps &&
                phase.accepted_step > 0 &&
                phase.previous_dt_s == accepted_->previous_dt_s &&
                phase.previous_dt_s > 0 && std::isfinite(phase.next_dt_s) &&
                phase.next_dt_s == phase.previous_dt_s,
            "unqualified variable-dt/non-mid-x placement");
  const auto &ch = channels_.at(e.channel);
  const auto &primary = accepted_->packets.at(e.primary);
  const auto &target = accepted_->packets.at(e.target);
  require(primary.state == ch.primary_state && target.state == ch.target_state,
          "channel/state mismatch");
  const auto &prov = e.provenance;
  require(std::isfinite(prov.collision_elapsed_s) &&
              prov.collision_elapsed_s >= 0 &&
              prov.selection_state_prefix == prefix,
          "stale selection clock/prefix");
  if (prov.kind == ArenaProvenance::Kind::ManufacturedForced)
    require(ch.model.empty() && prov.model.rfind("manufactured/", 0) == 0 &&
                prov.selection_method == "manufactured_forced" &&
                prov.dataset_sha256.empty() &&
                prov.kernel_manifest_sha256.empty() &&
                prov.kernel_binding_file_sha256.empty(),
            "manufactured provenance");
  else
    require(!ch.model.empty() && prov.dataset_sha256 == ch.dataset_sha256 &&
                prov.kernel_manifest_sha256 == ch.kernel_manifest_sha256 &&
                prov.kernel_binding_file_sha256 ==
                    ch.kernel_binding_file_sha256 &&
                prov.model == ch.model &&
                (prov.selection_method == "ssa_weighted_pair_clock_v1" ||
                 prov.selection_method == "forced_kernel_integration_only"),
            "conditional physical-data provenance not bound to channel");
  auto primary_remaining = debit(primary.count, e.event_weight),
       target_remaining = debit(target.count, e.event_weight);
  const auto pview = view(e.primary), tview = view(e.target);
  require(!e.products.empty() && e.products.size() <= 8,
          "bounded product list");
  std::vector<std::string> actual_states;
  std::set<std::string> newids;
  bool scattered = false;
  long double before = kinetic(pview) + states_.at(primary.state).internal_J +
                       states_.at(primary.state).reference_J + kinetic(tview) +
                       states_.at(target.state).internal_J +
                       states_.at(target.state).reference_J,
              after = 0;
  std::array<long double, 3> dp{}, momentum_activity{};
  for (int d = 0; d < 3; ++d) {
    dp[d] =
        -static_cast<long double>(pview.mass_kg) *
            pview.proper_velocity_m_s[d] -
        static_cast<long double>(tview.mass_kg) * tview.proper_velocity_m_s[d];
    momentum_activity[d] = std::abs(static_cast<long double>(pview.mass_kg) *
                                    pview.proper_velocity_m_s[d]) +
                           std::abs(static_cast<long double>(tview.mass_kg) *
                                    tview.proper_velocity_m_s[d]);
  }
  int primary_q = -states_.at(primary.state).charge_e, target_q = 0;
  for (const auto &product : e.products) {
    nonempty(product.id);
    require(states_.count(product.state) && newids.insert(product.id).second,
            "unknown or duplicate product");
    const bool is_primary = product.id == e.scattered_primary_id;
    if (is_primary) {
      require(product.id == e.primary &&
                  product.location == ArenaProduct::At::Primary && !scattered,
              "scattered primary identity/placement");
      scattered = true;
    } else
      require(!accepted_->ever_packet_ids.count(product.id),
              "reused immutable logical product ID");
    const auto &state = states_.at(product.state);
    actual_states.push_back(product.state);
    GapParticle p;
    p.mass_kg = state.mass_kg;
    p.charge_C = state.charge_e * pic_e;
    p.weight = e.event_weight;
    p.proper_velocity_m_s = product.proper_velocity_m_s;
    p.position_m = product.location == ArenaProduct::At::Primary
                       ? pview.position_m
                       : tview.position_m;
    native_particle_guard(p);
    after += kinetic(p) + state.internal_J + state.reference_J;
    for (int d = 0; d < 3; ++d) {
      long double v =
          static_cast<long double>(p.mass_kg) * p.proper_velocity_m_s[d];
      dp[d] += v;
      momentum_activity[d] += std::abs(v);
    }
    if (product.location == ArenaProduct::At::Primary)
      primary_q += state.charge_e;
    else
      target_q += state.charge_e;
  }
  require(e.scattered_primary_id.empty() || scattered,
          "missing scattered primary");
  auto expected = ch.product_states;
  std::sort(expected.begin(), expected.end());
  std::sort(actual_states.begin(), actual_states.end());
  require(actual_states == expected, "channel products mismatch");
  require(
      (primary_q == 0 && target_q == 0) || pview.position_m == tview.position_m,
      "reaction creates nonlocal charge without qualified current/work path");
  const long double eps = 256 * std::numeric_limits<double>::epsilon();
  const long double e_bound = eps * (std::abs(before) + std::abs(after));
  require(std::isfinite(before) && std::isfinite(after) && e_bound > 0 &&
              std::abs(after - before) <= e_bound,
          "event kinetic/internal/reference energy closure");
  long double pscale = 0;
  for (auto v : momentum_activity)
    pscale += v;
  for (auto v : dp)
    require(std::abs(v) <= eps * pscale, "event momentum closure");
  auto next = std::make_shared<ArenaSnapshot>(*accepted_);
  auto candidate = capture(*pic_);
  std::set<std::uint64_t> remove;
  remove.insert(primary.native_ids.begin(), primary.native_ids.end());
  remove.insert(target.native_ids.begin(), target.native_ids.end());
  candidate.particles.erase(
      std::remove_if(candidate.particles.begin(), candidate.particles.end(),
                     [&](const auto &p) { return remove.count(p.id); }),
      candidate.particles.end());
  next->packets.erase(e.primary);
  next->packets.erase(e.target);
  auto insert = [&](std::string id, const std::string &state, ArenaCount count,
                    const GapParticle &pose, std::uint64_t retained,
                    const std::string &origin) {
    if (zero(count))
      return;
    ArenaPacket record{id, state, origin, e.id, count, {}};
    const auto &spec = states_.at(state);
    auto weights = positive_carriers(count);
    for (std::size_t i = 0; i < weights.size(); ++i) {
      GapParticle p = pose;
      p.weight = weights[i];
      p.mass_kg = spec.mass_kg;
      p.charge_C = spec.charge_e * pic_e;
      if (i == 0 && retained)
        p.id = retained;
      else {
        require(candidate.next_id < std::numeric_limits<std::uint64_t>::max(),
                "native ID exhaustion");
        p.id = candidate.next_id++;
      }
      native_particle_guard(p);
      record.native_ids.push_back(p.id);
      candidate.particles.push_back(p);
    }
    require(next->packets.emplace(id, record).second,
            "logical insertion collision");
    next->ever_packet_ids.insert(id);
  };
  if (!zero(target_remaining))
    insert(e.target, target.state, target_remaining, tview,
           target.native_ids.front(), target.origin);
  if (!zero(primary_remaining)) {
    nonempty(e.untouched_primary_id);
    require(e.untouched_primary_id != e.primary &&
                !next->ever_packet_ids.count(e.untouched_primary_id) &&
                !newids.count(e.untouched_primary_id),
            "partial-primary untouched remainder ID");
    insert(e.untouched_primary_id, primary.state, primary_remaining, pview, 0,
           primary.origin);
  } else
    require(e.untouched_primary_id.empty(), "unexpected untouched primary ID");
  for (const auto &product : e.products) {
    auto pose = product.location == ArenaProduct::At::Primary ? pview : tview;
    pose.proper_velocity_m_s = product.proper_velocity_m_s;
    insert(product.id, product.state, {e.event_weight, 0}, pose,
           product.id == e.primary ? primary.native_ids.front() : 0,
           "reaction/" + e.id);
  }
  candidate.fields_current = false;
  next->event_ids.insert(e.id);
  ArenaJson products = ArenaJson::array();
  for (const auto &product : e.products)
    products.push_back(
        {{"id", product.id},
         {"state", product.state},
         {"u", product.proper_velocity_m_s},
         {"location", product.location == ArenaProduct::At::Primary
                          ? "primary"
                          : "target"}});
  ArenaJson audit = {
      {"kind", "reaction"},
      {"event_id", e.id},
      {"channel", e.channel},
      {"primary", e.primary},
      {"target", e.target},
      {"event_weight", e.event_weight},
      {"primary_before", count_json(primary.count)},
      {"primary_u", pview.proper_velocity_m_s},
      {"target_u", tview.proper_velocity_m_s},
      {"primary_position_m", pview.position_m},
      {"target_position_m", tview.position_m},
      {"target_before", count_json(target.count)},
      {"primary_remainder", count_json(primary_remaining)},
      {"target_remainder", count_json(target_remaining)},
      {"untouched_primary_id", e.untouched_primary_id},
      {"products", products},
      {"physical_stage_s", phase.physical_stage_s},
      {"state_prefix", prefix},
      {"mesh_prefix", phase.mesh_prefix},
      {"clock_prefix", phase.clock_prefix},
      {"phase", phase.kind == ArenaPhase::Kind::InitialFieldProof
                    ? "initial_field_proof"
                    : "provisional_equal_dt_mid_x"},
      {"provenance_kind", prov.kind == ArenaProvenance::Kind::ManufacturedForced
                              ? "manufactured_forced"
                              : "conditional_physical_data"},
      {"dataset_sha256", prov.dataset_sha256},
      {"kernel_manifest_sha256", prov.kernel_manifest_sha256},
      {"kernel_binding_file_sha256", prov.kernel_binding_file_sha256},
      {"conditional_model", prov.model},
      {"selection_method", prov.selection_method},
      {"selection_state_prefix", prov.selection_state_prefix},
      {"collision_elapsed_s", prov.collision_elapsed_s},
      {"energy_defect_J", double((after - before) * e.event_weight)},
      {"energy_roundoff_bound_J", double(e_bound * e.event_weight)},
      {"reference_convention",
       "kinetic+state_internal+state_reference; no extra threshold heat"},
      {"external_particle_source_debit_J", 0},
      {"MCC_heating_and_EM_current_qualified", false}};
  for (int d = 0; d < 3; ++d)
    audit["momentum_defect_kg_m_s"][d] = double(dp[d] * e.event_weight);
  require(std::isfinite(audit.at("energy_defect_J").get<double>()) &&
              audit.at("energy_roundoff_bound_J").get<double>() > 0,
          "weighted event receipt unrepresentable");
  journal_append(*next, audit);
  return stage(std::move(next), candidate, std::move(audit), fault);
}
PreparedArena PicArena::prepare_transport(double dt, int fault) {
  require(!active_scope_,
          "transport forbidden during frozen collision compound");
  validate_population(*accepted_);
  require(pic_->fields_current_ && !pic_->options_.distributed_fields,
          "transport requires solved serial fields");
  require(std::isfinite(dt) && dt > 0, "transport interval");
  const auto prefix = state_prefix();
  auto before = capture(*pic_);
  auto next = std::make_shared<ArenaSnapshot>(*accepted_);
  struct Sink {
    GapParticle particle;
    int kind;
    double time;
  };
  std::map<std::uint64_t, Sink> sinks;
  bool restored = false;
  try {
    pic_->advance_observed(dt, [&](const GapParticle &p, int hit, double time) {
      require(sinks.emplace(p.id, Sink{p, hit, time}).second,
              "duplicate native sink");
    });
    std::set<std::uint64_t> resident;
    for (const auto &p : pic_->particles_)
      resident.insert(p.id);
    ArenaJson new_sinks = ArenaJson::array();
    for (const auto &[id, packet] : accepted_->packets) {
      std::size_t lost = 0;
      for (auto native : packet.native_ids)
        if (sinks.count(native))
          ++lost;
      require(lost == 0 || lost == packet.native_ids.size(),
              "partial native carrier loss; atomic refusal");
      if (!lost) {
        for (auto native : packet.native_ids)
          require(resident.count(native), "unreported native removal");
        continue;
      }
      const auto &first = sinks.at(packet.native_ids.front());
      auto total = big();
      ArenaJson carriers = ArenaJson::array();
      long double K = 0;
      for (auto native : packet.native_ids) {
        const auto &row = sinks.at(native);
        require(row.kind == first.kind && row.time == first.time &&
                    row.particle.position_m == first.particle.position_m &&
                    row.particle.proper_velocity_m_s ==
                        first.particle.proper_velocity_m_s,
                "inconsistent logical sink");
        auto weight = exact(row.particle.weight);
        require(BN_add(total.get(), total.get(), weight.get()),
                "sink number accumulation");
        K += kinetic(row.particle) * row.particle.weight;
        carriers.push_back(
            {{"native_id", native},
             {"weight", row.particle.weight},
             {"incident_K_J",
              pic_kinetic_energy_J(row.particle) * row.particle.weight}});
      }
      auto expected = exact(packet.count);
      require(BN_cmp(total.get(), expected.get()) == 0, "sink Count mismatch");
      const auto &spec = states_.at(packet.state);
      std::string kind = first.kind == 0   ? "lower_foil"
                         : first.kind == 1 ? "upper_foil"
                         : first.kind == 3 ? "short"
                         : first.kind == 4 ? "axial_escape"
                                           : "outer_xy_escape";
      ArenaJson receipt = {
          {"receipt_id",
           "native_sink/" + std::to_string(next->revision + 1) + "/" + id},
          {"packet_id", id},
          {"state", packet.state},
          {"origin", packet.origin},
          {"last_event", packet.last_event},
          {"count", count_json(packet.count)},
          {"sink", kind},
          {"time_s", first.time},
          {"position_m", first.particle.position_m},
          {"proper_velocity_m_s", first.particle.proper_velocity_m_s},
          {"incident_K_J", double(K)},
          {"internal_per_particle_J", spec.internal_J},
          {"reference_per_particle_J", spec.reference_J},
          {"native_carrier_receipts", carriers},
          {"wall_heat_assumed", false}};
      new_sinks.push_back(receipt);
      next->sinks.push_back(receipt);
      next->packets.erase(id);
    }
    validate_population(*next);
    ++next->accepted_steps;
    next->previous_dt_s = dt;
    auto candidate = capture(*pic_);
    swap_native(*pic_, before);
    restored = true;
    ArenaJson audit = {
        {"kind", "native_transport"},
        {"state_prefix", prefix},
        {"dt_s", dt},
        {"start_time_s", pic_->ledger_.time_s},
        {"end_time_s", candidate.ledger.time_s},
        {"accepted_step", next->accepted_steps},
        {"species_sink_receipts", new_sinks},
        {"neutral_motion_owner", "actual native GapParticle drift/Boris(q=0)"},
        {"internal_reference_owner",
         "canonical Count/state receipt; no wall heat assumed"}};
    journal_append(*next, audit);
    return stage(std::move(next), candidate, std::move(audit), fault);
  } catch (...) {
    if (!restored)
      swap_native(*pic_, before);
    throw;
  }
}
void PicArena::save(const std::string &path) const {
  require(!active_scope_,
          "restart publication forbidden during private compound");
  validate_population(*accepted_);
  require(pic_->fields_current_, "restart requires current solved fields");
  journal_validate(*accepted_, sha_bytes(initial_prefix_ + registry_sha_));
  TemporaryNativeFile temporary(path);
  pic_->save_checkpoint(temporary.path);
  std::ifstream input(temporary.path, std::ios::binary);
  std::string bytes((std::istreambuf_iterator<char>(input)), {});
  require(input.eof() || bool(input), "native checkpoint read");
  std::vector<double> potential;
  for (int i = 0; i < pic_->potential_.Size(); ++i)
    potential.push_back(pic_->potential_[i]);
  ArenaJson body = {{"schema", "canonical_native_arena_v3"},
                    {"source_sha256", source_identity()},
                    {"registry_sha256", registry_sha_},
                    {"registry", registry()},
                    {"initial_coupling_prefix", initial_prefix_},
                    {"ABI", {80, 168, 200, sizeof(GapPic)}},
                    {"native_bytes_hex", hex_bytes(bytes)},
                    {"native_bytes_sha256", sha_bytes(bytes)},
                    {"metadata", metadata(*accepted_)},
                    {"saved_potential_V", potential}};
  atomic_json({{"payload", body}, {"sha256", sha_bytes(body.dump())}}, path);
}
std::unique_ptr<PicArena>
PicArena::restore(GapPic &p, const std::string &path,
                  const std::vector<ArenaStateSpec> &states,
                  const std::vector<ArenaChannel> &channels,
                  const std::string &prefix) {
  require(p.particles_.empty(), "restore requires fresh native population");
  auto before = capture(p);
  try {
    std::ifstream input(path);
    ArenaJson doc;
    input >> doc;
    require(bool(input) && doc.size() == 2, "arena restart envelope");
    input >> std::ws;
    require(input.eof(), "trailing arena restart data");
    const auto &body = doc.at("payload");
    require(doc.at("sha256") == sha_bytes(body.dump()) && body.size() == 10 &&
                body.at("schema") == "canonical_native_arena_v3" &&
                body.at("source_sha256") == source_identity() &&
                body.at("initial_coupling_prefix") == prefix,
            "arena restart identity/integrity");
    require(body.at("ABI") == ArenaJson::array({80, 168, 200, sizeof(GapPic)}),
            "arena restart ABI");
    const auto bytes = unhex(body.at("native_bytes_hex"));
    require(body.at("native_bytes_sha256") == sha_bytes(bytes),
            "arena native checkpoint digest");
    auto state =
        std::make_shared<ArenaSnapshot>(read_metadata(body.at("metadata")));
    TemporaryNativeFile temporary(path);
    {
      std::ofstream output(temporary.path, std::ios::binary);
      output.write(bytes.data(), bytes.size());
      require(bool(output), "native restore temporary write");
    }
    p.load_checkpoint(temporary.path);
    std::vector<ArenaPacket> packets;
    for (const auto &[id, x] : state->packets) {
      (void)id;
      packets.push_back(x);
    }
    auto arena =
        std::make_unique<PicArena>(p, states, channels, packets, prefix);
    require(body.at("registry_sha256") == arena->registry_sha_ &&
                arena->registry() == body.at("registry"),
            "arena immutable registry");
    journal_validate(*state, sha_bytes(prefix + arena->registry_sha_));
    for (const auto &entry : state->journal) {
      const auto &a = entry.at("payload");
      if (a.at("kind") == "reaction") {
        const auto &ch =
            arena->channels_.at(a.at("channel").get<std::string>());
        if (ch.model.empty())
          require(a.at("provenance_kind") == "manufactured_forced" &&
                      a.at("conditional_model")
                              .get<std::string>()
                              .rfind("manufactured/", 0) == 0,
                  "restart manufactured provenance");
        else
          require(a.at("provenance_kind") == "conditional_physical_data" &&
                      a.at("dataset_sha256") == ch.dataset_sha256 &&
                      a.at("kernel_manifest_sha256") ==
                          ch.kernel_manifest_sha256 &&
                      a.at("kernel_binding_file_sha256") ==
                          ch.kernel_binding_file_sha256 &&
                      a.at("conditional_model") == ch.model,
                  "restart conditional data/model provenance");
      }
    }

    arena->accepted_ = state;
    arena->validate_population(*state);
    const auto values = body.at("saved_potential_V").get<std::vector<double>>();
    require(values.size() == std::size_t(p.potential_.Size()),
            "restart potential dimensions");
    long double delta = 0, scale = 0;
    for (std::size_t i = 0; i < values.size(); ++i) {
      require(std::isfinite(values[i]), "restart potential finite");
      long double d = p.potential_[i] - values[i];
      delta += d * d;
      scale += static_cast<long double>(values[i]) * values[i];
    }
    require(std::sqrt(delta) <= 1e-10L * std::max(std::sqrt(scale), 1e-30L),
            "restart cold native solve differs from saved-state field budget");
    return arena;
  } catch (...) {
    swap_native(p, before);
    throw;
  }
}
PreparedArena
PicArena::prepare_compound(std::size_t maximum_reactions,
                           const std::function<void(PicArena &)> &builder) {
  require(!active_scope_, "nested compound preparation forbidden");
  require(bool(builder) && maximum_reactions > 0 && maximum_reactions <= 4096,
          "explicit bounded compound resource cap required");
  require(scope_counter_ < std::numeric_limits<std::uint64_t>::max(),
          "compound scope exhaustion");
  validate_population(*accepted_);
  auto native_before = capture(*pic_);
  auto before = accepted_;
  const auto prefix = state_prefix();
  active_scope_ = ++scope_counter_;
  compound_base_revision_ = before->revision;
  compound_limit_ = maximum_reactions;
  try {
    builder(*this);
    require(accepted_->revision >= before->revision &&
                accepted_->revision - before->revision <= maximum_reactions,
            "compound revision bound");
    require(accepted_->accepted_steps == before->accepted_steps &&
                accepted_->previous_dt_s == before->previous_dt_s &&
                accepted_->sinks == before->sinks &&
                pic_->ledger_.time_s == native_before.ledger.time_s &&
                options_equal(pic_->options_, native_before.options),
            "compound changed transport/source stage");
    const auto &a = pic_->ledger_;
    const auto &b = native_before.ledger;
    require(a.authored_charge_C == b.authored_charge_C &&
                a.authored_kinetic_energy_J == b.authored_kinetic_energy_J &&
                a.lower_injected_charge_C == b.lower_injected_charge_C &&
                a.upper_injected_charge_C == b.upper_injected_charge_C &&
                a.injected_energy_J == b.injected_energy_J &&
                a.lower_absorbed_charge_C == b.lower_absorbed_charge_C &&
                a.upper_absorbed_charge_C == b.upper_absorbed_charge_C &&
                a.lower_absorbed_energy_J == b.lower_absorbed_energy_J &&
                a.upper_absorbed_energy_J == b.upper_absorbed_energy_J &&
                a.short_absorbed_charge_C == b.short_absorbed_charge_C &&
                a.short_absorbed_energy_J == b.short_absorbed_energy_J &&
                a.escaped_charge_C == b.escaped_charge_C &&
                a.escaped_kinetic_energy_J == b.escaped_kinetic_energy_J &&
                a.axial_escaped_charge_C == b.axial_escaped_charge_C &&
                pic_->horizontal_transfers_ == native_before.horizontal &&
                pic_->vertical_transfers_ == native_before.vertical,
            "compound changed external source/sink/current prefix");
    require(accepted_->journal.size() >= before->journal.size(),
            "compound journal truncated");
    for (std::size_t i = 0; i < before->journal.size(); ++i)
      require(accepted_->journal[i] == before->journal[i],
              "compound rewrote accepted prefix");
    for (std::size_t i = before->journal.size(); i < accepted_->journal.size();
         ++i)
      require(accepted_->journal[i].at("payload").at("kind") == "reaction",
              "nonreaction subevent in collision compound");
    validate_population(*accepted_);
    auto final = std::make_shared<ArenaSnapshot>(*accepted_);
    auto candidate = capture(*pic_);
    const auto subevents = accepted_->revision - before->revision;
    swap_native(*pic_, native_before);
    accepted_ = before;
    active_scope_ = 0;
    compound_limit_ = 0;
    ArenaJson audit = {{"kind", "private_reaction_compound"},
                       {"state_prefix", prefix},
                       {"subevents", subevents},
                       {"maximum_reactions", maximum_reactions},
                       {"before_chain", before->chain},
                       {"after_chain", final->chain},
                       {"native_trajectory_time_unchanged", true},
                       {"external_source_sink_current_prefix_unchanged", true},
                       {"escaped_child_tokens_inert", true}};
    return stage(std::move(final), candidate, std::move(audit), 0);
  } catch (...) {
    // native_before contains the original state until the explicit successful
    // restore above; stage() itself restores the base on any later failure.
    if (active_scope_) {
      swap_native(*pic_, native_before);
      accepted_ = before;
      active_scope_ = 0;
      compound_limit_ = 0;
    }
    throw;
  }
}
} // namespace raikou::mfem_backend
