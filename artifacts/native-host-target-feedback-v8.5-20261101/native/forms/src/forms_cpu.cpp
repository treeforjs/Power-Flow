// forms_cpu.cpp -- CPU reference of the native v7 producer forms.
// SHA-256 through OpenSSL EVP (the portable core in forms_core.h is the
// device implementation; both are FIPS 180-4 and are compared in the tests).
#include "pf_v7/forms.hpp"
#include <algorithm>
#include <cstring>
#include <map>
#include <memory>
#include <numeric>
#include <openssl/evp.h>
#include <stdexcept>

namespace pf_v7 {
namespace {
static_assert(std::endian::native == std::endian::little);
static_assert(offsetof(Carrier, position_m) == 8 && offsetof(Carrier, proper_velocity_m_s) == 32 &&
              offsetof(Carrier, charge_C) == 56 && offsetof(Carrier, mass_kg) == 64 &&
              offsetof(Carrier, weight) == 72);
static_assert(offsetof(Path, field_start_m) == 8 && offsetof(Path, field_displacement_m) == 32 &&
              offsetof(Path, flight_s) == 56);
static_assert(offsetof(Sink, plane_m) == 32 && offsetof(Sink, segment_m) == 48 &&
              offsetof(Sink, hit_point_m) == 72 && offsetof(Sink, incident_proper_velocity_m_s) == 96 &&
              offsetof(Sink, flight_s) == 120 && offsetof(Sink, discarded_flight_s) == 128);
static_assert(offsetof(CensusRecord, packet_key) == 80 && offsetof(CensusRecord, state_index) == 88);
static_assert(pf7::c_squared == 299792458.0 * 299792458.0);

// Canonical byte writer (v5 primitives) feeding one reusable EVP context.
class Stream {
public:
  Stream() : ctx_(EVP_MD_CTX_new()) {
    if (!ctx_)
      throw std::bad_alloc();
    start();
  }
  ~Stream() { EVP_MD_CTX_free(ctx_); }
  Stream(const Stream &) = delete;
  Stream &operator=(const Stream &) = delete;
  void start() {
    if (EVP_DigestInit_ex(ctx_, EVP_sha256(), nullptr) != 1)
      throw std::runtime_error("SHA-256 init");
  }
  void raw(const void *p, std::size_t n) {
    if (n && EVP_DigestUpdate(ctx_, p, n) != 1)
      throw std::runtime_error("SHA-256 update");
  }
  void u64(std::uint64_t v) {
    unsigned char b[8];
    for (int i = 0; i < 8; ++i)
      b[i] = static_cast<unsigned char>(v >> (8 * i));
    raw(b, 8);
  }
  void f64(double x) { u64(std::bit_cast<std::uint64_t>(x)); }
  void text(std::string_view s) {
    u64(s.size());
    raw(s.data(), s.size());
  }
  Digest finish() {
    Digest d{};
    unsigned n = 0;
    if (EVP_DigestFinal_ex(ctx_, d.data(), &n) != 1 || n != 32)
      throw std::runtime_error("SHA-256 final");
    start();
    return d;
  }

private:
  EVP_MD_CTX *ctx_;
};
std::string_view tag(std::uint32_t section, int level) {
  pf7_u32 n = 0;
  const char *t = pf7::section_tag(section, level, n);
  return {t, n};
}
Digest node(Stream &s, const Digest &a, const Digest &b) {
  s.text("native-v5/merkle-node");
  s.raw(a.data(), 32);
  s.raw(b.data(), 32);
  return s.finish();
}
Digest tree_of(Stream &s, std::span<const Digest> leaves) {
  if (leaves.empty()) {
    s.text("native-v5/merkle-empty");
    return s.finish();
  }
  std::vector<Digest> level(leaves.begin(), leaves.end());
  while (level.size() > 1) {
    std::vector<Digest> up;
    up.reserve((level.size() + 1) / 2);
    for (std::size_t i = 0; i + 1 < level.size(); i += 2)
      up.push_back(node(s, level[i], level[i + 1]));
    if (level.size() % 2)
      up.push_back(level.back());
    level.swap(up);
  }
  return level.front();
}
Digest root_of(Stream &s, std::uint32_t section, std::uint64_t records, std::span<const Digest> leaves) {
  const auto tree = tree_of(s, leaves);
  s.text(tag(section, pf7::level_root));
  s.u64(pf7::leaf_records);
  s.u64(records);
  s.u64(leaves.size());
  s.raw(tree.data(), 32);
  return s.finish();
}
// Block list over a record source: fill(j, m, buffer) writes chunk j's m
// canonical records into buffer (or returns a pointer to them).
template <class Source>
BlockList block_list_of(std::uint32_t section, std::uint64_t n, Source source) {
  const std::uint64_t R = pf7::section_record_bytes(section);
  const auto chunk_tag = tag(section, pf7::level_chunk), block_tag = tag(section, pf7::level_block);
  Stream s;
  BlockList out;
  out.records = n;
  const std::uint64_t chunks = (n + pf7::chunk_records - 1) / pf7::chunk_records;
  const std::uint64_t leaves = (n + pf7::leaf_records - 1) / pf7::leaf_records;
  out.leaves.reserve(leaves);
  std::vector<Digest> chunk_digests;
  chunk_digests.reserve(pf7::leaf_chunks);
  for (std::uint64_t k = 0; k < leaves; ++k) {
    chunk_digests.clear();
    const std::uint64_t first_chunk = k * pf7::leaf_chunks;
    const std::uint64_t last_chunk = std::min(chunks, first_chunk + pf7::leaf_chunks);
    for (std::uint64_t j = first_chunk; j < last_chunk; ++j) {
      const std::uint64_t m = std::min<std::uint64_t>(pf7::chunk_records, n - j * pf7::chunk_records);
      s.text(chunk_tag);
      s.u64(j);
      s.u64(m);
      s.raw(source(j, m), m * R);
      chunk_digests.push_back(s.finish());
    }
    const std::uint64_t nk = std::min<std::uint64_t>(pf7::leaf_records, n - k * pf7::leaf_records);
    s.text(block_tag);
    s.u64(k);
    s.u64(nk);
    s.u64(chunk_digests.size());
    for (const auto &d : chunk_digests)
      s.raw(d.data(), 32);
    out.leaves.push_back(s.finish());
  }
  out.root_digest = root_of(s, section, n, out.leaves);
  out.root = hex(out.root_digest);
  return out;
}
} // namespace

std::string hex(const Digest &d) {
  static constexpr char digits[] = "0123456789abcdef";
  std::string out(64, '0');
  for (std::size_t i = 0; i < 32; ++i) {
    out[2 * i] = digits[d[i] >> 4];
    out[2 * i + 1] = digits[d[i] & 15];
  }
  return out;
}
Digest digest_of_words(const std::uint32_t *w) {
  Digest d{};
  for (int i = 0; i < 8; ++i)
    for (int b = 0; b < 4; ++b)
      d[4 * i + b] = static_cast<unsigned char>(w[i] >> (24 - 8 * b));
  return d;
}
void words_of_digest(const Digest &d, std::uint32_t *w) {
  for (int i = 0; i < 8; ++i)
    w[i] = (std::uint32_t(d[4 * i]) << 24) | (std::uint32_t(d[4 * i + 1]) << 16) |
           (std::uint32_t(d[4 * i + 2]) << 8) | std::uint32_t(d[4 * i + 3]);
}

// ---------------------------------------------------------------- (a)
void ExactSum::add_product(double a, double b) {
  const auto t = pf7::product_term(a, b);
  if (t.bad)
    status_ |= pf7::sum_nonfinite_input;
  pf7::add_term(digits_.data(), t);
  ++terms_;
  if (++pending_ >= pf7::acc_terms_per_normalisation) {
    pf7::normalise(digits_.data());
    pending_ = 0;
  }
}
void ExactSum::merge(const ExactSum &o) {
  auto a = digits(), b = o.digits();
  for (int k = 0; k < pf7::acc_digits; ++k)
    a[k] += b[k];
  pf7::normalise(a.data());
  digits_ = a;
  pending_ = 0;
  terms_ += o.terms_;
  status_ |= o.status_;
}
double ExactSum::value(std::uint32_t *status) const {
  auto d = digits_;
  std::uint32_t st = status_;
  const double v = pf7::round_sum(d.data(), &st);
  if (status)
    *status = st;
  return v;
}
ExactSum::Digits ExactSum::digits() const {
  auto d = digits_;
  pf7::normalise(d.data());
  return d;
}
ExactSum ExactSum::from_digits(const Digits &d, std::uint32_t status, std::uint64_t terms) {
  ExactSum s;
  s.digits_ = d;
  pf7::normalise(s.digits_.data());
  s.status_ = status;
  s.terms_ = terms;
  return s;
}
double kinetic_energy_J(double mass_kg, const double *u) {
  return pf7::kinetic_energy(mass_kg, u[0], u[1], u[2]);
}
namespace {
SumValue finish(const ExactSum &s, std::uint32_t extra) {
  SumValue v;
  v.value = v.sum = s.value(&v.status);
  v.status |= extra;
  return v;
}
} // namespace
SumValue kinetic_change_J(std::span<const Carrier> census, std::span<const Carrier> endpoints,
                          std::span<const Sink> sinks) {
  ExactSum s;
  std::uint32_t extra = 0;
  std::size_t j = 0, k = 0;
  for (const auto &c : census) {
    const double *u_end = nullptr;
    if (k < sinks.size() && sinks[k].id == c.id)
      u_end = sinks[k++].incident_proper_velocity_m_s;
    else if (j < endpoints.size() && endpoints[j].id == c.id)
      u_end = endpoints[j++].proper_velocity_m_s;
    else {
      extra |= pf7::sum_census_missing;
      continue;
    }
    s.add_product(c.weight, pf7::kinetic_change(c.mass_kg, c.proper_velocity_m_s, u_end));
  }
  if (j != endpoints.size() || k != sinks.size())
    extra |= pf7::sum_census_missing;
  return finish(s, extra);
}
SumValue kinetic_change_J_aligned(std::span<const Carrier> census, std::span<const double> u_end) {
  if (u_end.size() != 3 * census.size())
    throw std::invalid_argument("aligned kinetic change: 3 end components per census carrier");
  ExactSum s;
  for (std::size_t i = 0; i < census.size(); ++i)
    s.add_product(census[i].weight,
                  pf7::kinetic_change(census[i].mass_kg, census[i].proper_velocity_m_s, &u_end[3 * i]));
  return finish(s, 0);
}
SumValue exact_dot(std::span<const double> a, std::span<const double> b) {
  if (a.size() != b.size())
    throw std::invalid_argument("exact dot: equal lengths required");
  ExactSum s;
  for (std::size_t i = 0; i < a.size(); ++i)
    s.add_product(a[i], b[i]);
  return finish(s, 0);
}
SumValue field_work_J(double dt_s, std::span<const double> E, std::span<const double> J) {
  auto v = exact_dot(E, J);
  v.value = dt_s * v.sum;
  if (!std::isfinite(v.value))
    v.status |= pf7::sum_overflow;
  return v;
}

// ---------------------------------------------------------------- (b)
BlockList block_list(std::uint32_t section, const void *records, std::uint64_t n) {
  if (section >= pf7::section_count)
    throw std::invalid_argument("unknown block-list section");
  const auto *base = static_cast<const unsigned char *>(records);
  const std::uint64_t R = pf7::section_record_bytes(section);
  return block_list_of(section, n, [&](std::uint64_t j, std::uint64_t) {
    return base + j * pf7::chunk_records * R;
  });
}
std::string block_root(std::uint32_t section, std::uint64_t records, std::span<const Digest> leaves) {
  if (leaves.size() != (records + pf7::leaf_records - 1) / pf7::leaf_records)
    throw std::invalid_argument("block list: leaf count differs from the record count");
  Stream s;
  return hex(root_of(s, section, records, leaves));
}
MotionSectionsV7 motion_sections_v7(std::span<const Carrier> endpoints, std::span<const Path> paths,
                                    std::span<const Sink> sinks, std::span<const Birth> births) {
  MotionSectionsV7 m;
  m.endpoints = endpoints.size();
  m.paths = paths.size();
  m.sinks = sinks.size();
  m.births = births.size();
  m.endpoints_root = block_list(pf7::section_endpoints, endpoints.data(), endpoints.size()).root;
  m.paths_root = block_list(pf7::section_paths, paths.data(), paths.size()).root;
  m.sinks_root = block_list(pf7::section_sinks, sinks.data(), sinks.size()).root;
  m.births_root = block_list(pf7::section_births, births.data(), births.size()).root;
  return m;
}
std::string motion_record_v7_sha256(const MotionHeaderV7 &r, const MotionSectionsV7 &m) {
  Stream e;
  e.text(motion_v7_schema);
  e.text("binding");
  e.text("full_time_position_and_proper_velocity");
  e.u64(r.accepted_step);
  e.u64(r.staged_arena_revision);
  e.u64(r.trial_stage);
  e.f64(r.begin_s);
  e.f64(r.end_s);
  e.f64(r.dt_s);
  for (const auto *x : {&r.staged_arena_sha256, &r.trial_field_sha256, &r.mesh_sha256, &r.material_sha256,
                        &r.species_sha256, &r.common_clock_sha256, &r.kinetic_policy_sha256})
    e.text(*x);
  e.text("budget");
  for (double x : r.budget)
    e.f64(x);
  e.text("census_sha256");
  e.text(r.census_sha256);
  e.text("coverage");
  e.text(r.coverage);
  e.text("current_A_m_sha256");
  e.text(r.current_A_m_sha256);
  e.text("endpoint_charge_C_sha256");
  e.text(r.endpoint_charge_C_sha256);
  e.text("initial_charge_C_sha256");
  e.text(r.initial_charge_C_sha256);
  e.text("kinetic_policy");
  e.text("native-finite-target-kinetic-policy-v1");
  for (double x : r.kinetic_policy)
    e.f64(x);
  for (const char *x : {"lower foil x span and inner y gap only", "absorbing escape", "axial escape",
                        "x y lower upper axial; later exact tie wins"})
    e.text(x);
  e.text("native_DKD_called");
  e.u64(r.native_DKD_called ? 1 : 0);
  e.text("native_Poisson_called");
  e.u64(r.native_Poisson_called ? 1 : 0);
  e.text("native_to_field_m");
  for (double x : r.native_to_field_m)
    e.f64(x);
  e.text("neutral_motion_owner");
  e.text(r.neutral_motion_owner);
  e.text("neutral_stream_dt_max_s");
  e.f64(r.neutral_stream_dt_max_s);
  e.text("motion_owner");
  e.text("ordinary_total_field");
  const std::pair<const char *, std::pair<std::uint64_t, const std::string *>> sections[4] = {
      {"endpoints", {m.endpoints, &m.endpoints_root}},
      {"paths", {m.paths, &m.paths_root}},
      {"sinks", {m.sinks, &m.sinks_root}},
      {"births", {m.births, &m.births_root}}};
  for (const auto &[name, v] : sections) {
    e.text(name);
    e.u64(v.first);
    e.text(*v.second);
  }
  return hex(e.finish());
}
std::string array_sha256(std::span<const double> v) {
  Stream e;
  e.text("motion-array-binary64-v1");
  e.u64(v.size());
  e.raw(v.data(), v.size() * sizeof(double));
  return hex(e.finish());
}

// ---------------------------------------------------------------- (c)
std::string census_states_sha256(std::span<const std::string> states) {
  for (std::size_t i = 1; i < states.size(); ++i)
    if (!(states[i - 1] < states[i]))
      throw std::invalid_argument("census states must be unique and in ascending byte order");
  Stream e;
  e.text("canonical-motion-census-v2/states");
  e.u64(states.size());
  for (const auto &s : states)
    e.text(s);
  return hex(e.finish());
}
CensusCheck check_census(std::span<const CensusRecord> c, std::uint64_t states) {
  const auto *w = reinterpret_cast<const pf7_u64 *>(c.data());
  for (std::uint64_t i = 0; i < c.size(); ++i)
    if (const auto code = pf7::census_check(w, c.size(), i, states))
      return {code, i};
  return {};
}
std::string census_v2_sha256(std::uint64_t carriers, const std::string &states_sha256,
                             const std::string &carriers_root) {
  Stream e;
  e.text(census_v2_tag);
  e.u64(carriers);
  e.text(states_sha256);
  e.text(carriers_root);
  return hex(e.finish());
}
std::string census_v2_sha256(std::span<const CensusRecord> c, std::span<const std::string> states) {
  const auto check = check_census(c, states.size());
  if (!check.ok())
    throw std::invalid_argument("census record " + std::to_string(check.index) + " fails canonical check " +
                                std::to_string(check.code));
  return census_v2_sha256(c.size(), census_states_sha256(states),
                          block_list(pf7::section_census, c.data(), c.size()).root);
}
std::vector<std::uint64_t> census_order(std::span<const CensusRecord> c) {
  std::vector<std::uint64_t> p(c.size());
  std::iota(p.begin(), p.end(), 0);
  std::stable_sort(p.begin(), p.end(),
                   [&](std::uint64_t a, std::uint64_t b) { return c[a].carrier.id < c[b].carrier.id; });
  return p;
}
std::vector<CensusRecord> census_records(std::span<const HostCensusEntry> entries,
                                         std::span<const std::string> states) {
  for (std::size_t i = 1; i < states.size(); ++i)
    if (!(states[i - 1] < states[i]))
      throw std::invalid_argument("census states must be unique and in ascending byte order");
  std::map<std::string_view, std::uint64_t> key;
  for (const auto &e : entries) {
    auto [it, inserted] = key.emplace(e.packet_id, e.carrier->id);
    if (!inserted)
      it->second = std::min(it->second, e.carrier->id);
  }
  std::vector<CensusRecord> out;
  out.reserve(entries.size());
  for (const auto &e : entries) {
    const auto s = std::lower_bound(states.begin(), states.end(), e.state);
    if (s == states.end() || *s != e.state)
      throw std::invalid_argument("census state not in the state table");
    out.push_back({*e.carrier, key.at(e.packet_id), std::uint64_t(s - states.begin())});
  }
  std::stable_sort(out.begin(), out.end(),
                   [](const CensusRecord &a, const CensusRecord &b) { return a.carrier.id < b.carrier.id; });
  return out;
}
BlockList census_block_list(std::span<const Carrier> rows, std::span<const std::uint64_t> keys,
                            std::span<const std::uint64_t> states) {
  if (keys.size() != rows.size() || states.size() != rows.size())
    throw std::invalid_argument("census block list: one key and one state per row");
  CensusRecord buffer[pf7::chunk_records];
  return block_list_of(pf7::section_census, rows.size(), [&](std::uint64_t j, std::uint64_t m) {
    const std::uint64_t first = j * pf7::chunk_records;
    for (std::uint64_t i = 0; i < m; ++i)
      buffer[i] = {rows[first + i], keys[first + i], states[first + i]};
    return reinterpret_cast<const unsigned char *>(buffer);
  });
}

// ---------------------------------------------------------------- work tolerance
void HallAxialMesh::edge_nodes(std::uint64_t e, std::uint64_t &a, std::uint64_t &b) const noexcept {
  pf7_u64 x = 0, y = 0;
  pf7::hall_edge_nodes(n, nz, e, x, y);
  a = x;
  b = y;
}
bool hall_axial_mesh_of(std::uint64_t edges, std::uint64_t faces, HallAxialMesh &out) {
  int found = 0;
  for (std::uint64_t n = 1; n <= 4096; ++n) {
    const std::uint64_t base = 2 * n * (n + 1), per = 2 * n * (n + 1) + (n + 1) * (n + 1);
    if (edges < base + per || (edges - base) % per)
      continue;
    const HallAxialMesh m{n, (edges - base) / per};
    if (m.nz >= 1 && m.ne() == edges && m.nb() == faces) {
      out = m;
      ++found;
    }
  }
  return found == 1;
}
std::string total_field_sha256(std::span<const double> E, std::span<const double> B) {
  Stream e;
  e.text("total-midpoint-E-B-binary64-v1");
  e.u64(E.size());
  e.raw(E.data(), E.size() * sizeof(double));
  e.u64(B.size());
  e.raw(B.data(), B.size() * sizeof(double));
  return hex(e.finish());
}
std::string trial_field_blocks_sha256(std::span<const double> E, std::span<const double> B) {
  Stream e;
  e.text("total-midpoint-E-B-blocks-v1");
  e.text(block_list(pf7::section_field_E, E.data(), E.size()).root);
  e.text(block_list(pf7::section_field_B, B.data(), B.size()).root);
  return hex(e.finish());
}
std::string trial_field_blocks_sha256(const std::string &E_root, const std::string &B_root) {
  Stream e;
  e.text("total-midpoint-E-B-blocks-v1");
  e.text(E_root);
  e.text(B_root);
  return hex(e.finish());
}
std::string current_blocks_sha256(std::span<const double> J) {
  return current_blocks_sha256(block_list(pf7::section_current_J, J.data(), J.size()).root);
}
std::string current_blocks_sha256(const std::string &J_root) {
  Stream e;
  e.text("motion-array-blocks-v1");
  e.text(J_root);
  return hex(e.finish());
}
void WorkCarrierSums::add(double weight, double charge, double K0, double d) {
  if (charge == 0)
    return; // neutral carriers exchange no work with the field
  const double K1 = K0 + d;
  gk.add_product(weight, std::fabs(d));
  ck.add_product(weight, K0);
  ck.add_product(weight, K1);
  q.add_product(std::fabs(charge), weight);
}
void WorkCarrierSums::merge(const WorkCarrierSums &o) {
  gk.merge(o.gk);
  ck.merge(o.ck);
  q.merge(o.q);
}
void work_field_scales(WorkScales &s, double dt, std::span<const double> E, std::span<const double> J,
                       const HallAxialMesh &mesh) {
  if (E.size() != J.size() || E.size() != mesh.ne())
    throw std::invalid_argument("work scales: E and J must hold one value per Hall axial mesh edge");
  ExactSum gw;
  std::vector<char> mark(mesh.nn(), 0);
  std::uint64_t a = 0, b = 0;
  for (std::size_t e = 0; e < J.size(); ++e) {
    gw.add_product(std::fabs(E[e]), std::fabs(J[e]));
    if (J[e] != 0) {
      mesh.edge_nodes(e, a, b);
      mark[a] = mark[b] = 1;
    }
  }
  double em = 0;
  for (std::size_t e = 0; e < J.size(); ++e) {
    mesh.edge_nodes(e, a, b);
    if (mark[a] || mark[b]) {
      const double v = std::fabs(E[e]);
      if (!(v <= em))
        em = v; // NaN propagates (refused)
    }
  }
  std::uint32_t st = 0;
  const double sum = gw.value(&st);
  s.status |= st;
  s.field_work_gross_J = dt * sum;
  if (!std::isfinite(s.field_work_gross_J))
    s.status |= pf7::sum_overflow;
  s.field_bound_V_m = em;
  if (!std::isfinite(em))
    s.status |= pf7::sum_nonfinite_input;
}
double work_tolerance(const WorkScales &s) {
  const double gross = s.kinetic_gross_J + s.field_work_gross_J;
  const double relative = 1e-10 * gross;
  double t = 1e-20 + relative;
  const double rounding = 0x1p-47 * s.kinetic_content_J;
  t = t + rounding;
  const double length = 8e-18 * s.field_bound_V_m;
  const double fixed_point = length * s.charge_abs_C;
  return t + fixed_point;
}
void work_finish(WorkScales &s, const WorkCarrierSums &c) {
  std::uint32_t st = 0;
  s.kinetic_gross_J = c.gk.value(&st);
  s.status |= st;
  st = 0;
  s.kinetic_content_J = c.ck.value(&st);
  s.status |= st;
  st = 0;
  s.charge_abs_C = c.q.value(&st);
  s.status |= st;
  s.tolerance_J = work_tolerance(s);
}
bool work_admits(const WorkScales &s, double defect) {
  return s.status == 0 && std::isfinite(defect) && std::isfinite(s.tolerance_J) && std::abs(defect) <= s.tolerance_J;
}

// ---------------------------------------------------------------- neutral stream
void require_kinetic_policy(const KineticPolicy &p) {
  auto need = [](bool ok, const char *m) {
    if (!ok)
      throw std::invalid_argument(m);
  };
  for (double v : {p.width_m, p.gap_m, p.short_z_m, p.feed_z_m})
    need(std::isfinite(v), "nonfinite kinetic boundary");
  need(p.width_m > 0 && p.gap_m > 0 && p.feed_z_m > p.short_z_m, "invalid kinetic extents");
  for (const double *f : {p.lower, p.upper})
    need(std::isfinite(f[0]) && std::isfinite(f[1]) && std::isfinite(f[2]) && std::isfinite(f[3]) && f[0] >= 0 &&
             f[1] <= p.width_m && f[0] < f[1] && f[2] >= 0 && f[3] <= p.gap_m && f[2] < f[3],
         "invalid kinetic foil");
  need(p.lower[3] < p.upper[2] && p.lower[0] == p.upper[0] && p.lower[1] == p.upper[1],
       "kinetic short bridge geometry");
}
NeutralStreamResult neutral_stream(const KineticPolicy &k, const Carrier &c, double begin_s, double end_s,
                                   double dt_s, double dt_max_s) {
  if (!(std::isfinite(dt_max_s) && dt_max_s > 0))
    throw std::invalid_argument("neutral stream interval bound must be positive and finite");
  static const char *const text[] = {"",
                                     "neutral authoritative interval",
                                     "neutral stream identity/species",
                                     "neutral momentum finite",
                                     "neutral momentum outside resolved range",
                                     "particle requires native boundary owner",
                                     "particle in/on native metal",
                                     "nonfinite kinetic displacement",
                                     "neutral incident kinetic operand"};
  const pf7::Policy p{k.width_m, k.gap_m, k.short_z_m, k.feed_z_m,
                      {k.lower[0], k.lower[1], k.lower[2], k.lower[3]},
                      {k.upper[0], k.upper[1], k.upper[2], k.upper[3]}};
  pf7::NeutralStream s{};
  // The sealed order: interval, identity, momentum, gamma, then the policy
  // guard and interior of first_kinetic_hit.
  auto code = pf7::neutral_stream(p, c.id, c.position_m, c.proper_velocity_m_s, c.charge_C, c.mass_kg, c.weight,
                                  begin_s, end_s, dt_s, dt_max_s, s);
  if (code == pf7::neutral_ok || code >= pf7::neutral_boundary_owner)
    require_kinetic_policy(k);
  if (code != pf7::neutral_ok)
    throw std::invalid_argument(text[code]);
  NeutralStreamResult r;
  r.endpoint = c;
  for (int d = 0; d < 3; ++d) {
    r.endpoint.position_m[d] = s.endpoint[d];
    r.displacement_m[d] = s.D[d];
  }
  r.hit = s.hit;
  r.begin_s = begin_s;
  r.end_s = end_s;
  r.dt_s = dt_s;
  r.flight_s = s.flight_s;
  r.remaining_s = s.remaining_s;
  r.incident_kinetic_J = s.incident_kinetic_J;
  return r;
}
} // namespace pf_v7
