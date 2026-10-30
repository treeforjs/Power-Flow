#include "raikou/mfem_pic_arena.h"
#include "raikou/mfem_pic_full_time.h"
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
constexpr std::size_t maximum_native_carriers = 100000;
constexpr std::size_t maximum_journal_entries = 4096;
constexpr std::size_t maximum_journal_bytes = 16 * 1024 * 1024;
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
  require(p.charge_C == 0 || p.charge_C * p.weight != 0,
          "charged native carrier underflow: unchanged-state refusal");
  require(pic_kinetic_energy_J(p) == 0 ||
              pic_kinetic_energy_J(p) * p.weight != 0,
          "native kinetic carrier underflow: unchanged-state refusal");
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
Big exact(ArenaCount n) { return n.copy_units(); }
void count_guard(ArenaCount n) {
  auto a = n.copy_units();
  require(!BN_is_negative(a.get()) &&
              BN_num_bits(a.get()) <= ArenaCount::maximum_bits,
          "bounded exact count guard");
}
bool zero(ArenaCount n) { return n.zero(); }
ArenaJson count_json(ArenaCount n) {
  return {{"unit_exponent", -1074},
          {"positive_carriers", n.positive_carriers()}};
}
ArenaCount read_count(const ArenaJson &j) {
  require(j.is_object() && j.size() == 2 && j.at("unit_exponent") == -1074,
          "exact count schema");
  require(j.at("positive_carriers").is_array() &&
              j.at("positive_carriers").size() <= ArenaCount::maximum_carriers,
          "bounded exact count codec");
  auto weights = j.at("positive_carriers").get<std::vector<double>>();
  auto count = ArenaCount::from_positive_carriers(weights);
  require(weights == count.positive_carriers(),
          "exact count canonical positive representation");
  return count;
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
ArenaJson native_rows_json(const std::vector<GapParticle> &rows) {
  ArenaJson j = ArenaJson::array();
  for (const auto &p : rows)
    j.push_back({{"id", p.id},
                 {"position_m", p.position_m},
                 {"u", p.proper_velocity_m_s},
                 {"q", p.charge_C},
                 {"m", p.mass_kg},
                 {"weight", p.weight}});
  return j;
}
std::vector<GapParticle> native_rows_read(const ArenaJson &j) {
  require(j.is_array() && j.size() <= maximum_native_carriers,
          "native history row dimensions");
  std::vector<GapParticle> out;
  for (const auto &x : j) {
    require(x.size() == 6, "native history row schema");
    GapParticle p;
    p.id = x.at("id");
    p.position_m = x.at("position_m").get<PicVector>();
    p.proper_velocity_m_s = x.at("u").get<PicVector>();
    p.charge_C = x.at("q");
    p.mass_kg = x.at("m");
    p.weight = x.at("weight");
    native_particle_guard(p);
    out.push_back(p);
  }
  return out;
}
ArenaJson metadata(const ArenaSnapshot &s) {
  ArenaJson p = ArenaJson::array();
  for (const auto &[id, x] : s.packets) {
    (void)id;
    p.push_back(packet_json(x));
  }
  ArenaJson result = {{"packets", p},
                      {"ever_packet_ids", s.ever_packet_ids},
                      {"event_ids", s.event_ids},
                      {"journal", s.journal},
                      {"sinks", s.sinks},
                      {"chain", s.chain},
                      {"revision", s.revision},
                      {"accepted_steps", s.accepted_steps},
                      {"previous_dt_s", s.previous_dt_s},
                      {"stochastic_completions", s.stochastic_completions},
                      {"journal_bytes", s.journal_bytes}};
  if (!s.full_time_context.is_null())
    result["full_time_context"] = s.full_time_context;
  return result;
}
ArenaSnapshot read_metadata(const ArenaJson &j) {
  require(j.size() == 11 || (j.size() == 12 && j.contains("full_time_context")),
          "metadata schema");
  ArenaSnapshot s;
  if (j.contains("full_time_context"))
    s.full_time_context = j.at("full_time_context");
  s.ever_packet_ids = j.at("ever_packet_ids").get<std::set<std::string>>();
  s.event_ids = j.at("event_ids").get<std::set<std::string>>();
  s.journal = j.at("journal");
  s.sinks = j.at("sinks");
  s.chain = j.at("chain");
  s.revision = j.at("revision");
  s.accepted_steps = j.at("accepted_steps");
  s.previous_dt_s = j.at("previous_dt_s");
  s.stochastic_completions = j.at("stochastic_completions");
  s.journal_bytes = j.at("journal_bytes");
  require(s.stochastic_completions.is_object() &&
              s.journal_bytes <= maximum_journal_bytes,
          "completion/history byte schema");
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
  require(s.journal.size() < maximum_journal_entries,
          "bounded prototype journal entry ceiling; validated anchors required "
          "for longer runs");
  audit["previous_chain"] = s.chain;
  const auto chain = sha_bytes(audit.dump());
  ArenaJson row = {{"payload", audit}, {"sha256", chain}};
  const auto bytes = row.dump().size();
  require(s.journal_bytes <= maximum_journal_bytes &&
              bytes <= maximum_journal_bytes - s.journal_bytes,
          "bounded prototype receipt byte ceiling; history is never truncated");
  s.journal.push_back(std::move(row));
  s.journal_bytes += bytes;
  s.chain = chain;
  ++s.revision;
}
void journal_validate(const ArenaSnapshot &s, const std::string &genesis) {
  require(s.journal.is_array() && s.sinks.is_array() &&
              s.revision == s.journal.size(),
          "journal dimensions");
  std::string chain = genesis;
  std::set<std::string> events;
  std::size_t bytes = 0;
  std::uint64_t steps = 0;
  double previous = 0;
  require(s.journal.size() <= maximum_journal_entries,
          "bounded history entries");
  for (const auto &x : s.journal) {
    const auto n = x.dump().size();
    require(bytes <= maximum_journal_bytes &&
                n <= maximum_journal_bytes - bytes,
            "bounded receipt bytes");
    bytes += n;
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
    } else if (type != "stochastic_exposure")
      throw std::invalid_argument("journal operation");
  }
  require(bytes == s.journal_bytes && chain == s.chain &&
              events == s.event_ids && steps == s.accepted_steps &&
              previous == s.previous_dt_s,
          "restart prefix mismatch");
}
} // namespace
std::string PicArena::sha256(const std::string &s) { return sha_bytes(s); }
bool PicArena::exact_count_equal(ArenaCount a, ArenaCount b) { return a == b; }
ArenaCount PicArena::debit(ArenaCount n, double w) { return n.debit(w); }
std::vector<double> PicArena::positive_carriers(ArenaCount n) {
  return n.positive_carriers();
}
struct PreparedArena::Impl {
  PicArena *owner = nullptr;
  std::uint64_t scope = 0;
  std::shared_ptr<const ArenaSnapshot> before, after;
  PicArena::NativeState native_before, native_after, work;
  ArenaJson audit;
  std::string candidate_prefix;
  std::shared_ptr<const FullTimeCandidateExposures> exposures;
  std::shared_ptr<const FullTimeAcceptedBoundary> boundary_before,
      boundary_after;
  std::shared_ptr<void> private_identity;
  std::uint64_t private_generation = 0;
  enum class Status { Staged, Published, RolledBack };
  Status status = Status::Staged;
};
PreparedArena::PreparedArena(std::unique_ptr<Impl> p) : impl_(std::move(p)) {}
PreparedArena::PreparedArena(PreparedArena &&) noexcept = default;
PreparedArena &PreparedArena::operator=(PreparedArena &&) noexcept = default;
PreparedArena::~PreparedArena() = default;
const ArenaJson &PreparedArena::audit() const { return impl_->audit; }
std::string PreparedArena::candidate_state_prefix() const {
  require(bool(impl_), "empty prepared token");
  return impl_->candidate_prefix;
}
std::shared_ptr<const FullTimeCandidateExposures>
PreparedArena::candidate_full_time_exposures() const {
  require(bool(impl_), "empty prepared token");
  return impl_->exposures;
}

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
    : PicArena(p, std::move(states), std::move(channels), std::move(packets),
               std::move(prefix), false) {}
PicArena::PicArena(GapPic &p, std::vector<ArenaStateSpec> states,
                   std::vector<ArenaChannel> channels,
                   std::vector<ArenaPacket> packets, std::string prefix,
                   bool total_field)
    : pic_(&p), initial_prefix_(std::move(prefix)) {
  nonempty(initial_prefix_);
  require(!p.options_.distributed_fields && (p.fields_current_ || total_field),
          "arena currently requires solved serial CPU fields");
  for (auto s : states) {
    nonempty(s.id);
    require(s.hydrogen >= 0 && s.oxygen >= 0 && s.charge_e >= -8 &&
                s.charge_e <= 8 && std::isfinite(s.mass_kg) && s.mass_kg > 0 &&
                std::isfinite(s.internal_J) && s.internal_J >= 0 &&
                std::isfinite(s.reference_J),
            "state specification");
    require(states_.emplace(s.id, s).second, "duplicate state");
  }
  for (auto c : channels) {
    nonempty(c.id);
    require(states_.count(c.primary_state) && states_.count(c.target_state) &&
                states_.at(c.target_state).charge_e == 0 &&
                !c.product_states.empty() && c.product_states.size() <= 8,
            "channel reactants");
    static_assert(std::numeric_limits<int>::digits <= 31,
                  "bounded catalog integer representation");
    auto checked_add = [](std::int64_t a, std::int64_t b) {
      require((b <= 0 || a <= std::numeric_limits<std::int64_t>::max() - b) &&
                  (b >= 0 || a >= std::numeric_limits<std::int64_t>::min() - b),
              "stoichiometry overflow");
      return a + b;
    };
    std::int64_t h =
        checked_add(-std::int64_t(states_.at(c.primary_state).hydrogen),
                    -std::int64_t(states_.at(c.target_state).hydrogen));
    std::int64_t o =
        checked_add(-std::int64_t(states_.at(c.primary_state).oxygen),
                    -std::int64_t(states_.at(c.target_state).oxygen));
    std::int64_t q = -std::int64_t(states_.at(c.primary_state).charge_e);
    for (const auto &id : c.product_states) {
      require(states_.count(id), "unknown product state");
      h = checked_add(h, states_.at(id).hydrogen);
      o = checked_add(o, states_.at(id).oxygen);
      q = checked_add(q, states_.at(id).charge_e);
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
  genesis_ = s;
  genesis_native_ = p.particles_;
  genesis_next_id_ = p.next_particle_id_;
  genesis_time_s_ = p.ledger_.time_s;
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
namespace {
using NativeIndex = std::map<std::uint64_t, const GapParticle *>;
NativeIndex validate_rows(const ArenaSnapshot &s,
                          const std::vector<GapParticle> &native,
                          const std::map<std::string, ArenaStateSpec> &states) {
  require(native.size() <= maximum_native_carriers,
          "bounded prototype native carrier ceiling");
  NativeIndex index;
  for (const auto &p : native)
    require(p.id > 0 && index.emplace(p.id, &p).second, "duplicate native ID");
  std::set<std::uint64_t> all;
  for (const auto &[id, p] : s.packets) {
    require(id == p.id && s.ever_packet_ids.count(id) &&
                states.count(p.state) && !p.native_ids.empty(),
            "logical packet identity/state");
    count_guard(p.count);
    require(!zero(p.count), "zero resident logical count");
    auto total = big();
    const GapParticle *first = nullptr;
    const auto &spec = states.at(p.state);
    for (auto native_id : p.native_ids) {
      require(all.insert(native_id).second, "native carrier owned twice");
      const auto *row = index.at(native_id);
      native_particle_guard(*row);
      require(row->charge_C == spec.charge_e * pic_e &&
                  row->mass_kg == spec.mass_kg,
              "native species binding");
      auto weight = exact(row->weight);
      require(BN_add(total.get(), total.get(), weight.get()),
              "carrier count accumulation");
      if (first)
        require(first->position_m == row->position_m &&
                    first->proper_velocity_m_s == row->proper_velocity_m_s,
                "divergent carrier positions/velocities");
      else
        first = row;
    }
    auto expected = exact(p.count);
    require(BN_cmp(total.get(), expected.get()) == 0,
            "native carriers do not exactly represent Count");
  }
  require(all.size() == native.size(), "unowned native particle");
  return index;
}
} // namespace
void PicArena::validate_population(const ArenaSnapshot &s) const {
  (void)validate_rows(s, pic_->particles_, states_);
}
std::vector<ArenaRateEntry>
PicArena::rate_population(std::size_t maximum_entries) const {
  require(pic_->particles_.size() <= maximum_entries,
          "native carrier census resource cap");
  std::size_t total_ids = 0;
  for (const auto &[id, p] : accepted_->packets) {
    (void)id;
    require(p.native_ids.size() <= maximum_entries - total_ids,
            "canonical carrier census resource cap");
    total_ids += p.native_ids.size();
  }
  auto index = validate_rows(*accepted_, pic_->particles_, states_);
  std::vector<ArenaRateEntry> out;
  out.reserve(pic_->particles_.size());
  for (const auto &[id, packet] : accepted_->packets)
    for (auto native : packet.native_ids)
      out.push_back({id, packet.state, *index.at(native)});
  return out;
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
    token->candidate_prefix =
        sha_bytes(initial_prefix_ + registry_sha_ + metadata(*next).dump() +
                  fingerprint(*pic_));
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
  return t.impl_ && !t.impl_->private_identity && t.impl_->owner == this &&
         t.impl_->scope == active_scope_ &&
         t.impl_->status == PreparedArena::Impl::Status::Staged &&
         accepted_ == t.impl_->before && matches(*pic_, t.impl_->native_before);
}
bool PicArena::try_publish(PreparedArena &t) noexcept {
  if (!can_publish(t))
    return false;
  swap_native(*pic_, t.impl_->work);
  accepted_ = t.impl_->after;
  if (t.impl_->boundary_after)
    full_time_boundary_ = t.impl_->boundary_after;
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
  if (t.impl_->boundary_after)
    full_time_boundary_ = t.impl_->boundary_before;
  t.impl_->status = PreparedArena::Impl::Status::RolledBack;
  return true;
}

std::string PicArena::state_prefix() const {
  if (full_time_boundary_)
    return full_time_boundary_->accepted_arena_sha256;
  return sha_bytes(initial_prefix_ + registry_sha_ +
                   metadata(*accepted_).dump() + fingerprint(*pic_));
}
struct PicArena::ReactionBuild {
  std::shared_ptr<ArenaSnapshot> snapshot;
  NativeState native;
  ArenaJson audit;
};
PreparedArena PicArena::prepare(const ArenaEvent &e, int fault) {
  require(!full_time_boundary_,
          "full-time reactions require a private collision scope");
  require(!active_scope_ ||
              accepted_->revision - compound_base_revision_ < compound_limit_,
          "compound reaction resource cap reached");
  require(pic_->fields_current_ && !pic_->options_.distributed_fields,
          "reaction requires solved serial synchronized fields");
  auto built =
      build_reaction(e, *accepted_, capture(*pic_), state_prefix(), false);
  return stage(std::move(built.snapshot), built.native, std::move(built.audit),
               fault);
}
PicArena::ReactionBuild PicArena::build_reaction(const ArenaEvent &e,
                                                 const ArenaSnapshot &base,
                                                 const NativeState &native,
                                                 const std::string &prefix,
                                                 bool full_time) const {
  require(base.journal.size() < maximum_journal_entries,
          "bounded prototype history is full");
  const ArenaSnapshot *accepted_ = &base;
  validate_rows(base, native.particles, states_);
  nonempty(e.id);
  require(!accepted_->event_ids.count(e.id) && e.primary != e.target &&
              accepted_->packets.count(e.primary) &&
              accepted_->packets.count(e.target),
          "duplicate event or invalid reactant identity");
  const auto &phase = e.phase;
  require(phase.state_prefix == prefix &&
              phase.physical_stage_s == native.ledger.time_s,
          "stale physical state/clock");
  nonempty(phase.mesh_prefix);
  nonempty(phase.clock_prefix);
  if (phase.kind == ArenaPhase::Kind::InitialFieldProof)
    require(phase.accepted_step == 0 && accepted_->accepted_steps == 0 &&
                phase.physical_stage_s == 0 && phase.previous_dt_s == 0 &&
                phase.next_dt_s == 0,
            "initial proof phase");
  else if (phase.kind == ArenaPhase::Kind::FullTimeEndpoint) {
    require(full_time && bool(phase.full_time_key),
            "unauthenticated full-time phase");
  } else
    require(!full_time && phase.accepted_step == accepted_->accepted_steps &&
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
  auto pose_of = [&](const ArenaPacket &p) {
    auto i = std::find_if(
        native.particles.begin(), native.particles.end(),
        [&](const auto &n) { return n.id == p.native_ids.front(); });
    require(i != native.particles.end(), "missing reactant carrier");
    return *i;
  };
  const auto pview = pose_of(primary), tview = pose_of(target);
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
  auto next = std::make_shared<ArenaSnapshot>(base);
  auto candidate = native;
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
      {"accepted_step", phase.accepted_step},
      {"previous_dt_s", phase.previous_dt_s},
      {"next_dt_s", phase.next_dt_s},
      {"phase",
       phase.kind == ArenaPhase::Kind::InitialFieldProof
           ? "initial_field_proof"
           : (full_time ? "full_time_endpoint" : "provisional_equal_dt_mid_x")},
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
  return {std::move(next), std::move(candidate), std::move(audit)};
}
PreparedArena PicArena::prepare_transport(double dt, int fault) {
  require(!full_time_boundary_,
          "full-time transport requires genuine total-field motion port");
  require(!active_scope_,
          "transport forbidden during frozen collision compound");
  require(accepted_->journal.size() < maximum_journal_entries,
          "bounded prototype history is full");
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
    ArenaJson endpoint_kinematics = ArenaJson::array();
    for (const auto &[id, packet] : next->packets) {
      const auto it = std::find_if(
          pic_->particles_.begin(), pic_->particles_.end(),
          [&](const auto &p) { return p.id == packet.native_ids.front(); });
      require(it != pic_->particles_.end(), "transport endpoint missing");
      endpoint_kinematics.push_back({{"packet_id", id},
                                     {"position_m", it->position_m},
                                     {"u", it->proper_velocity_m_s}});
    }
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
        {"endpoint_kinematics", endpoint_kinematics},
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
  if (full_time_boundary_) {
    save_full_time(path);
    return;
  }
  require(!active_scope_,
          "restart publication forbidden during private compound");
  validate_population(*accepted_);
  require(pic_->fields_current_, "restart requires current solved fields");
  journal_validate(*accepted_, sha_bytes(initial_prefix_ + registry_sha_));
  validate_history();
  TemporaryNativeFile temporary(path);
  pic_->save_checkpoint(temporary.path);
  std::ifstream input(temporary.path, std::ios::binary);
  std::string bytes((std::istreambuf_iterator<char>(input)), {});
  require(input.eof() || bool(input), "native checkpoint read");
  std::vector<double> potential;
  for (int i = 0; i < pic_->potential_.Size(); ++i)
    potential.push_back(pic_->potential_[i]);
  ArenaJson body = {{"schema", "canonical_native_arena_exact_v1"},
                    {"genesis", metadata(*genesis_)},
                    {"genesis_native", native_rows_json(genesis_native_)},
                    {"genesis_next_id", genesis_next_id_},
                    {"genesis_time_s", genesis_time_s_},
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
    require(std::filesystem::file_size(path) <= 64 * 1024 * 1024,
            "bounded arena restart size");
    std::ifstream input(path);
    ArenaJson doc;
    input >> doc;
    require(bool(input) && doc.size() == 2, "arena restart envelope");
    input >> std::ws;
    require(input.eof(), "trailing arena restart data");
    const auto &body = doc.at("payload");
    require(doc.at("sha256") == sha_bytes(body.dump()) && body.size() == 14 &&
                body.at("schema") == "canonical_native_arena_exact_v1" &&
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

    arena->genesis_ =
        std::make_shared<ArenaSnapshot>(read_metadata(body.at("genesis")));
    arena->genesis_native_ = native_rows_read(body.at("genesis_native"));
    arena->genesis_next_id_ = body.at("genesis_next_id");
    arena->genesis_time_s_ = body.at("genesis_time_s");
    arena->accepted_ = state;
    arena->validate_population(*state);
    arena->validate_history();
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
  return compound_impl(maximum_reactions, builder, false);
}
PreparedArena PicArena::prepare_compound_with_completion(
    std::size_t maximum_reactions,
    const std::function<ArenaStochasticCompletion(PicArena &)> &builder) {
  require(bool(builder), "completion builder required");
  return compound_impl(
      maximum_reactions,
      [&](PicArena &a) {
        auto c = builder(a);
        a.complete_stochastic(c);
      },
      true);
}
PreparedArena
PicArena::compound_impl(std::size_t maximum_reactions,
                        const std::function<void(PicArena &)> &builder,
                        bool with_completion) {
  require(!active_scope_, "nested compound preparation forbidden");
  require(bool(builder) && maximum_reactions > 0 && maximum_reactions <= 4096,
          "explicit bounded compound resource cap required");
  require(scope_counter_ < std::numeric_limits<std::uint64_t>::max(),
          "compound scope exhaustion");
  require(accepted_->journal.size() < maximum_journal_entries,
          "bounded prototype history is full");
  require(accepted_->journal.size() + maximum_reactions +
                  (with_completion ? 1 : 0) <=
              maximum_journal_entries,
          "compound declared resource budget exceeds bounded history");
  validate_population(*accepted_);
  auto native_before = capture(*pic_);
  auto before = accepted_;
  const auto prefix = state_prefix();
  active_scope_ = ++scope_counter_;
  compound_base_revision_ = before->revision;
  compound_limit_ = maximum_reactions;
  completion_compound_ = with_completion;
  compound_base_chain_ = before->chain;
  try {
    builder(*this);
    require(accepted_->revision >= before->revision &&
                accepted_->revision - before->revision <=
                    maximum_reactions + (with_completion ? 1 : 0),
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
      require(accepted_->journal[i].at("payload").at("kind") == "reaction" ||
                  (with_completion && i + 1 == accepted_->journal.size() &&
                   accepted_->journal[i].at("payload").at("kind") ==
                       "stochastic_exposure"),
              "nonreaction subevent in collision compound");
    require(!with_completion || accepted_->journal.back().at("payload").at(
                                    "kind") == "stochastic_exposure",
            "missing terminal stochastic completion");
    validate_population(*accepted_);
    auto final = std::make_shared<ArenaSnapshot>(*accepted_);
    auto candidate = capture(*pic_);
    const auto subevents =
        accepted_->revision - before->revision - (with_completion ? 1 : 0);
    swap_native(*pic_, native_before);
    accepted_ = before;
    active_scope_ = 0;
    compound_limit_ = 0;
    completion_compound_ = false;
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
      completion_compound_ = false;
    }
    throw;
  }
}

namespace {
ArenaJson completion_json(const ArenaStochasticCompletion &c) {
  return {{"integration_owner_id", c.integration_owner_id},
          {"before_rng_clock_sha256", c.before_rng_clock_sha256},
          {"after_rng_clock_sha256", c.after_rng_clock_sha256},
          {"common_clock_sha256", c.common_clock_sha256},
          {"channel_id", c.channel_id},
          {"dataset_sha256", c.dataset_sha256},
          {"kernel_manifest_sha256", c.kernel_manifest_sha256},
          {"kernel_binding_file_sha256", c.kernel_binding_file_sha256},
          {"model", c.model},
          {"stage_sequence", c.stage_sequence},
          {"frozen_physical_stage_s", c.frozen_physical_stage_s},
          {"common_interval_begin_s", c.common_interval_begin_s},
          {"common_interval_end_s", c.common_interval_end_s}};
}
ArenaStochasticCompletion completion_read(const ArenaJson &j) {
  require(j.is_object() && j.size() == 13, "stochastic completion schema");
  ArenaStochasticCompletion c;
  c.integration_owner_id = j.at("integration_owner_id");
  c.before_rng_clock_sha256 = j.at("before_rng_clock_sha256");
  c.after_rng_clock_sha256 = j.at("after_rng_clock_sha256");
  c.common_clock_sha256 = j.at("common_clock_sha256");
  c.channel_id = j.at("channel_id");
  c.dataset_sha256 = j.at("dataset_sha256");
  c.kernel_manifest_sha256 = j.at("kernel_manifest_sha256");
  c.kernel_binding_file_sha256 = j.at("kernel_binding_file_sha256");
  c.model = j.at("model");
  c.stage_sequence = j.at("stage_sequence");
  c.frozen_physical_stage_s = j.at("frozen_physical_stage_s");
  c.common_interval_begin_s = j.at("common_interval_begin_s");
  c.common_interval_end_s = j.at("common_interval_end_s");
  return c;
}
void completion_validate(
    const ArenaStochasticCompletion &c, const ArenaChannel &channel,
    const std::optional<ArenaStochasticCompletion> &previous, double time) {
  nonempty(c.integration_owner_id);
  require(digest(c.before_rng_clock_sha256) &&
              digest(c.after_rng_clock_sha256) && digest(c.common_clock_sha256),
          "stochastic RNG/clock hashes must bind canonical owner state");
  require(c.channel_id == channel.id &&
              c.dataset_sha256 == channel.dataset_sha256 &&
              c.kernel_manifest_sha256 == channel.kernel_manifest_sha256 &&
              c.kernel_binding_file_sha256 ==
                  channel.kernel_binding_file_sha256 &&
              c.model == channel.model && !channel.model.empty(),
          "stochastic exposure channel/data/kernel/model binding");
  require(c.frozen_physical_stage_s == time &&
              std::isfinite(c.common_interval_begin_s) &&
              std::isfinite(c.common_interval_end_s) &&
              c.common_interval_begin_s >= 0 &&
              c.common_interval_end_s > c.common_interval_begin_s &&
              c.frozen_physical_stage_s >= c.common_interval_begin_s &&
              c.frozen_physical_stage_s <= c.common_interval_end_s,
          "stochastic frozen-stage/common-interval clock");
  if (previous)
    require(previous->stage_sequence <
                    std::numeric_limits<std::uint64_t>::max() &&
                c.stage_sequence == previous->stage_sequence + 1 &&
                c.before_rng_clock_sha256 == previous->after_rng_clock_sha256 &&
                c.channel_id == previous->channel_id,
            "stale stochastic owner sequence/RNG-clock prefix");
  else
    require(c.stage_sequence == 0,
            "initial stochastic exposure sequence must be zero");
}
} // namespace
std::optional<ArenaStochasticCompletion>
PicArena::stochastic_completion(const std::string &owner) const {
  auto it = accepted_->stochastic_completions.find(owner);
  if (it == accepted_->stochastic_completions.end())
    return std::nullopt;
  auto c = completion_read(*it);
  require(c.integration_owner_id == owner, "stochastic completion owner key");
  return c;
}

namespace {
void validate_scoped_ssa(const ArenaJson &journal, std::size_t start,
                         const ArenaStochasticCompletion &c) {
  require(journal.is_array() && start <= journal.size(),
          "stochastic scoped journal range");
  const double duration = c.common_interval_end_s - c.common_interval_begin_s;
  require(std::isfinite(duration) && duration > 0,
          "stochastic scoped duration");
  double previous_elapsed = 0;
  for (std::size_t i = start; i < journal.size(); ++i) {
    const auto &a = journal[i].at("payload");
    require(a.at("kind") == "reaction",
            "stochastic compound overlaps transport/exposure history");
    const auto method = a.at("selection_method").get<std::string>();
    // Forced conditional-kernel comparisons retain their explicit non-SSA tag.
    // Their declared test offset is not certified as an actual waiting clock.
    if (method == "forced_kernel_integration_only")
      continue;
    require(method == "ssa_weighted_pair_clock_v1",
            "unsupported stochastic scoped selection");
    require(a.at("physical_stage_s") == c.frozen_physical_stage_s &&
                a.at("channel") == c.channel_id &&
                a.at("dataset_sha256") == c.dataset_sha256 &&
                a.at("kernel_manifest_sha256") == c.kernel_manifest_sha256 &&
                a.at("kernel_binding_file_sha256") ==
                    c.kernel_binding_file_sha256 &&
                a.at("conditional_model") == c.model &&
                a.at("clock_prefix") == c.common_clock_sha256,
            "SSA child differs from bound common "
            "stage/channel/data/kernel/model/clock");
    const double elapsed = a.at("collision_elapsed_s");
    require(std::isfinite(elapsed) && elapsed > previous_elapsed &&
                elapsed <= duration,
            "SSA child elapsed offset outside ordered common duration");
    previous_elapsed = elapsed;
  }
}
} // namespace
void PicArena::complete_stochastic(const ArenaStochasticCompletion &c) {
  require(active_scope_ && completion_compound_,
          "stochastic completion only inside its private compound");
  completion_validate(c, channels_.at(c.channel_id),
                      stochastic_completion(c.integration_owner_id),
                      pic_->ledger_.time_s);
  validate_scoped_ssa(accepted_->journal, compound_base_revision_, c);
  auto next = std::make_shared<ArenaSnapshot>(*accepted_);
  next->stochastic_completions[c.integration_owner_id] = completion_json(c);
  ArenaJson audit = {{"kind", "stochastic_exposure"},
                     {"completion", completion_json(c)},
                     {"stage_before_chain", compound_base_chain_},
                     {"stage_after_reaction_chain", accepted_->chain}};
  journal_append(*next, audit);
  accepted_ = std::move(next);
}
void PicArena::validate_history() const {
  require(!full_time_boundary_, "full-time history requires typed replay");
  require(bool(genesis_) && genesis_->journal.empty() &&
              genesis_->sinks.empty() && genesis_->event_ids.empty() &&
              genesis_->revision == 0 && genesis_->accepted_steps == 0 &&
              genesis_->previous_dt_s == 0 &&
              genesis_->chain == sha_bytes(initial_prefix_ + registry_sha_),
          "invalid arena genesis prefix");
  auto replay = *genesis_;
  auto rows = genesis_native_;
  auto next_id = genesis_next_id_;
  double time = genesis_time_s_;
  require(std::isfinite(time) && time >= 0 && next_id > 0,
          "genesis native clock/ID");
  std::set<std::string> initial_ids;
  for (const auto &[id, p] : replay.packets) {
    (void)p;
    initial_ids.insert(id);
  }
  require(initial_ids == replay.ever_packet_ids, "genesis ever-ID set");
  validate_rows(replay, rows, states_);
  for (const auto &p : rows)
    require(p.id < next_id, "genesis native ID prefix");
  auto lookup = [&](std::uint64_t id) -> GapParticle & {
    auto it = std::find_if(rows.begin(), rows.end(),
                           [&](const auto &p) { return p.id == id; });
    require(it != rows.end(), "history native ID missing");
    return *it;
  };
  for (const auto &entry : accepted_->journal) {
    const auto &a = entry.at("payload");
    require(a.at("previous_chain") == replay.chain &&
                entry.at("sha256") == sha_bytes(a.dump()),
            "history hash chain");
    if (a.at("kind") == "reaction") {
      auto id = a.at("event_id").get<std::string>(),
           primary_id = a.at("primary").get<std::string>(),
           target_id = a.at("target").get<std::string>();
      nonempty(id);
      require(!replay.event_ids.count(id) && primary_id != target_id,
              "history replayed event identity");
      const auto primary = replay.packets.at(primary_id),
                 target = replay.packets.at(target_id);
      const auto &ch = channels_.at(a.at("channel").get<std::string>());
      require(primary.state == ch.primary_state &&
                  target.state == ch.target_state,
              "history channel/reactant states");
      const auto pview = lookup(primary.native_ids.front()),
                 tview = lookup(target.native_ids.front());
      require(a.at("physical_stage_s") == time &&
                  a.at("primary_u") == ArenaJson(pview.proper_velocity_m_s) &&
                  a.at("target_u") == ArenaJson(tview.proper_velocity_m_s) &&
                  a.at("primary_position_m") == ArenaJson(pview.position_m) &&
                  a.at("target_position_m") == ArenaJson(tview.position_m),
              "history reactant kinematics/clock continuity");
      const auto phase = a.at("phase").get<std::string>();
      if (phase == "initial_field_proof")
        require(time == 0 && replay.accepted_steps == 0 &&
                    a.at("accepted_step") == 0 && a.at("previous_dt_s") == 0 &&
                    a.at("next_dt_s") == 0,
                "history initial collision phase");
      else
        require(phase == "provisional_equal_dt_mid_x" &&
                    replay.accepted_steps > 0 &&
                    a.at("accepted_step") == replay.accepted_steps &&
                    a.at("previous_dt_s") == replay.previous_dt_s &&
                    a.at("next_dt_s") == replay.previous_dt_s,
                "history collision split phase");
      const auto elapsed = a.at("collision_elapsed_s").get<double>();
      require(std::isfinite(elapsed) && elapsed >= 0 &&
                  a.at("selection_state_prefix") == a.at("state_prefix") &&
                  digest(a.at("state_prefix").get<std::string>()),
              "history selection prefix/elapsed clock");
      nonempty(a.at("mesh_prefix"));
      nonempty(a.at("clock_prefix"));
      if (ch.model.empty())
        require(a.at("provenance_kind") == "manufactured_forced" &&
                    a.at("selection_method") == "manufactured_forced" &&
                    a.at("conditional_model")
                            .get<std::string>()
                            .rfind("manufactured/", 0) == 0 &&
                    a.at("dataset_sha256") == "" &&
                    a.at("kernel_manifest_sha256") == "" &&
                    a.at("kernel_binding_file_sha256") == "",
                "history manufactured provenance");
      else
        require(
            a.at("provenance_kind") == "conditional_physical_data" &&
                a.at("dataset_sha256") == ch.dataset_sha256 &&
                a.at("kernel_manifest_sha256") == ch.kernel_manifest_sha256 &&
                a.at("kernel_binding_file_sha256") ==
                    ch.kernel_binding_file_sha256 &&
                a.at("conditional_model") == ch.model &&
                (a.at("selection_method") == "ssa_weighted_pair_clock_v1" ||
                 a.at("selection_method") == "forced_kernel_integration_only"),
            "history conditional provenance/selection");
      const double weight = a.at("event_weight");
      require(exact_count_equal(read_count(a.at("primary_before")),
                                primary.count) &&
                  exact_count_equal(read_count(a.at("target_before")),
                                    target.count),
              "history reactant count continuity");
      auto pr = debit(primary.count, weight), tr = debit(target.count, weight);
      require(exact_count_equal(pr, read_count(a.at("primary_remainder"))) &&
                  exact_count_equal(tr, read_count(a.at("target_remainder"))),
              "history one-debit/remainder ledger");
      std::vector<std::string> product_states;
      std::set<std::string> product_ids;
      const auto &products = a.at("products");
      require(products.is_array() && !products.empty() && products.size() <= 8,
              "history product dimensions");
      long double before =
                      kinetic(pview) + states_.at(primary.state).internal_J +
                      states_.at(primary.state).reference_J + kinetic(tview) +
                      states_.at(target.state).internal_J +
                      states_.at(target.state).reference_J,
                  after = 0;
      std::array<long double, 3> dp{}, activity{};
      for (int d = 0; d < 3; ++d) {
        dp[d] = -static_cast<long double>(pview.mass_kg) *
                    pview.proper_velocity_m_s[d] -
                static_cast<long double>(tview.mass_kg) *
                    tview.proper_velocity_m_s[d];
        activity[d] = std::abs(static_cast<long double>(pview.mass_kg) *
                               pview.proper_velocity_m_s[d]) +
                      std::abs(static_cast<long double>(tview.mass_kg) *
                               tview.proper_velocity_m_s[d]);
      }
      int primary_q = -states_.at(primary.state).charge_e, target_q = 0;
      for (const auto &x : products) {
        require(x.size() == 4, "history product schema");
        auto pid = x.at("id").get<std::string>(),
             state = x.at("state").get<std::string>(),
             location = x.at("location").get<std::string>();
        nonempty(pid);
        require(product_ids.insert(pid).second && states_.count(state) &&
                    (location == "primary" || location == "target"),
                "history product ID/state/location");
        if (pid == primary_id)
          require(location == "primary", "history scattered primary location");
        else
          require(!replay.ever_packet_ids.count(pid),
                  "history reused product ID");
        product_states.push_back(state);
        const auto &spec = states_.at(state);
        GapParticle p = location == "primary" ? pview : tview;
        p.mass_kg = spec.mass_kg;
        p.charge_C = spec.charge_e * pic_e;
        p.weight = weight;
        p.proper_velocity_m_s = x.at("u").get<PicVector>();
        native_particle_guard(p);
        after += kinetic(p) + spec.internal_J + spec.reference_J;
        for (int d = 0; d < 3; ++d) {
          long double v =
              static_cast<long double>(p.mass_kg) * p.proper_velocity_m_s[d];
          dp[d] += v;
          activity[d] += std::abs(v);
        }
        if (location == "primary")
          primary_q += spec.charge_e;
        else
          target_q += spec.charge_e;
      }
      auto expected = ch.product_states;
      std::sort(expected.begin(), expected.end());
      std::sort(product_states.begin(), product_states.end());
      require(product_states == expected &&
                  ((primary_q == 0 && target_q == 0) ||
                   pview.position_m == tview.position_m),
              "history stoichiometry/local charge ledger");
      const long double factor = 256 * std::numeric_limits<double>::epsilon(),
                        bound = factor * (std::abs(before) + std::abs(after));
      require(bound > 0 && std::abs(after - before) <= bound &&
                  a.at("energy_defect_J") ==
                      double((after - before) * weight) &&
                  a.at("energy_roundoff_bound_J") == double(bound * weight) &&
                  a.at("external_particle_source_debit_J") == 0 &&
                  a.at("reference_convention") ==
                      "kinetic+state_internal+state_reference; no extra "
                      "threshold heat",
              "history local total-energy receipt");
      long double pscale = 0;
      for (auto v : activity)
        pscale += v;
      for (int d = 0; d < 3; ++d)
        require(std::abs(dp[d]) <= factor * pscale &&
                    a.at("momentum_defect_kg_m_s").at(d) ==
                        double(dp[d] * weight),
                "history momentum receipt");
      std::set<std::uint64_t> removed(primary.native_ids.begin(),
                                      primary.native_ids.end());
      removed.insert(target.native_ids.begin(), target.native_ids.end());
      rows.erase(
          std::remove_if(rows.begin(), rows.end(),
                         [&](const auto &p) { return removed.count(p.id); }),
          rows.end());
      replay.packets.erase(primary_id);
      replay.packets.erase(target_id);
      auto insert = [&](const std::string &pid, const std::string &state,
                        ArenaCount count, GapParticle pose,
                        std::uint64_t retained, const std::string &origin) {
        if (zero(count))
          return;
        ArenaPacket p{pid, state, origin, id, count, {}};
        const auto &spec = states_.at(state);
        const auto weights = positive_carriers(count);
        for (std::size_t i = 0; i < weights.size(); ++i) {
          pose.weight = weights[i];
          pose.mass_kg = spec.mass_kg;
          pose.charge_C = spec.charge_e * pic_e;
          if (i == 0 && retained)
            pose.id = retained;
          else {
            require(next_id < std::numeric_limits<std::uint64_t>::max(),
                    "history native ID exhaustion");
            pose.id = next_id++;
          }
          native_particle_guard(pose);
          p.native_ids.push_back(pose.id);
          rows.push_back(pose);
        }
        require(replay.packets.emplace(pid, p).second,
                "history logical insertion");
        replay.ever_packet_ids.insert(pid);
      };
      if (!zero(tr))
        insert(target_id, target.state, tr, tview, target.native_ids.front(),
               target.origin);
      const auto untouched = a.at("untouched_primary_id").get<std::string>();
      if (!zero(pr)) {
        nonempty(untouched);
        require(!replay.ever_packet_ids.count(untouched) &&
                    !product_ids.count(untouched),
                "history partial-primary remainder ID");
        insert(untouched, primary.state, pr, pview, 0, primary.origin);
      } else
        require(untouched.empty(), "history unexpected remainder ID");
      for (const auto &x : products) {
        auto pid = x.at("id").get<std::string>();
        auto pose = x.at("location") == "primary" ? pview : tview;
        pose.proper_velocity_m_s = x.at("u").get<PicVector>();
        insert(pid, x.at("state"), {weight, 0}, pose,
               pid == primary_id ? primary.native_ids.front() : 0,
               "reaction/" + id);
      }
      replay.event_ids.insert(id);
    } else if (a.at("kind") == "native_transport") {
      const double dt = a.at("dt_s");
      require(std::isfinite(dt) && dt > 0 && a.at("start_time_s") == time &&
                  a.at("end_time_s") == time + dt &&
                  a.at("accepted_step") == replay.accepted_steps + 1,
              "history native transport clock");
      const double upper_time = std::max(time + dt, (time + .5 * dt) + .5 * dt);
      std::set<std::string> lost;
      std::set<std::uint64_t> removed;
      for (const auto &receipt : a.at("species_sink_receipts")) {
        require(receipt.is_object() && receipt.size() == 15,
                "history sink schema");
        const auto id = receipt.at("packet_id").get<std::string>();
        require(lost.insert(id).second, "duplicate history sink");
        const auto packet = replay.packets.at(id);
        const auto &spec = states_.at(packet.state);
        require(receipt.at("receipt_id") ==
                        "native_sink/" + std::to_string(replay.revision + 1) +
                            "/" + id &&
                    receipt.at("state") == packet.state &&
                    receipt.at("origin") == packet.origin &&
                    receipt.at("last_event") == packet.last_event &&
                    exact_count_equal(read_count(receipt.at("count")),
                                      packet.count) &&
                    receipt.at("internal_per_particle_J") == spec.internal_J &&
                    receipt.at("reference_per_particle_J") ==
                        spec.reference_J &&
                    receipt.at("wall_heat_assumed") == false,
                "history sink identity/species/count/reference");
        const auto sink = receipt.at("sink").get<std::string>();
        require(sink == "lower_foil" || sink == "upper_foil" ||
                    sink == "short" || sink == "axial_escape" ||
                    sink == "outer_xy_escape",
                "history sink classification");
        double t = receipt.at("time_s");
        require(std::isfinite(t) && t >= time && t <= upper_time,
                "history sink clock");
        const auto x = receipt.at("position_m").get<PicVector>(),
                   u = receipt.at("proper_velocity_m_s").get<PicVector>();
        const auto &carriers = receipt.at("native_carrier_receipts");
        require(carriers.size() == packet.native_ids.size(),
                "history sink carriers");
        long double K = 0;
        for (std::size_t i = 0; i < packet.native_ids.size(); ++i) {
          auto p = lookup(packet.native_ids[i]);
          p.position_m = x;
          p.proper_velocity_m_s = u;
          native_particle_guard(p);
          require(carriers[i].size() == 3 &&
                      carriers[i].at("native_id") == p.id &&
                      carriers[i].at("weight") == p.weight &&
                      carriers[i].at("incident_K_J") ==
                          pic_kinetic_energy_J(p) * p.weight,
                  "history sink native carrier receipt");
          K += kinetic(p) * p.weight;
          removed.insert(p.id);
        }
        require(receipt.at("incident_K_J") == double(K),
                "history sink total kinetic receipt");
        replay.sinks.push_back(receipt);
        replay.packets.erase(id);
      }
      rows.erase(
          std::remove_if(rows.begin(), rows.end(),
                         [&](const auto &p) { return removed.count(p.id); }),
          rows.end());
      std::set<std::string> seen;
      for (const auto &endpoint : a.at("endpoint_kinematics")) {
        require(endpoint.size() == 3, "history endpoint schema");
        const auto id = endpoint.at("packet_id").get<std::string>();
        require(seen.insert(id).second, "duplicate endpoint");
        const auto &packet = replay.packets.at(id);
        for (auto native : packet.native_ids) {
          auto &p = lookup(native);
          p.position_m = endpoint.at("position_m").get<PicVector>();
          p.proper_velocity_m_s = endpoint.at("u").get<PicVector>();
          native_particle_guard(p);
        }
      }
      require(seen.size() == replay.packets.size(),
              "missing survivor endpoint");
      time += dt;
      ++replay.accepted_steps;
      replay.previous_dt_s = dt;
    } else if (a.at("kind") == "stochastic_exposure") {
      require(a.size() == 5 &&
                  a.at("stage_after_reaction_chain") == replay.chain,
              "stochastic receipt reaction-chain binding");
      auto c = completion_read(a.at("completion"));
      std::optional<ArenaStochasticCompletion> previous;
      auto prior = replay.stochastic_completions.find(c.integration_owner_id);
      if (prior != replay.stochastic_completions.end())
        previous = completion_read(*prior);
      completion_validate(c, channels_.at(c.channel_id), previous, time);
      const auto base = a.at("stage_before_chain").get<std::string>();
      std::size_t start = 0;
      if (base != genesis_->chain) {
        bool found = false;
        for (std::size_t i = 0; i < replay.journal.size(); ++i)
          if (replay.journal[i].at("sha256") == base) {
            start = i + 1;
            found = true;
            break;
          }
        require(found, "unknown stochastic base chain");
      }
      validate_scoped_ssa(replay.journal, start, c);
      replay.stochastic_completions[c.integration_owner_id] =
          completion_json(c);
    } else
      throw std::invalid_argument("unknown semantic history operation");
    journal_append(replay, a);
    require(replay.chain == entry.at("sha256").get<std::string>(),
            "history receipt hash");
    validate_rows(replay, rows, states_);
  }
  require(metadata(replay) == metadata(*accepted_) &&
              next_id == pic_->next_particle_id_ &&
              time == pic_->ledger_.time_s &&
              rows.size() == pic_->particles_.size(),
          "restart logical/neutral/sink/history continuity");
  for (std::size_t i = 0; i < rows.size(); ++i)
    require(particle_equal(rows[i], pic_->particles_[i]),
            "restart final native history mismatch");
}
#include "mfem_pic_full_time.inc"
} // namespace raikou::mfem_backend
