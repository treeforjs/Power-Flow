// Native v5 controls: canonical binary encodings (round trip, ambiguity,
// refusals), block Merkle digests, complete closed steps with the genuine Hall
// v2 motion provider (zero total field, no field solve), bounded history
// anchors (chain verification, tamper detection, caps and typed refusals,
// K-independence of accepted states), configurable carrier cap/budget, and
// resident populations of 1e5 and 1e6 carriers.
#include "resident_fixture.hpp"
#include "raikou/mfem_pic_canonical.h"
#include <random>
namespace cn = m::canonical;
using Kind = m::FullTimeResourceRefusal::Kind;
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
template <class F> void refuses_kind(F f, Kind kind, const char *name) {
  bool ok = false;
  try {
    f();
  } catch (const m::FullTimeResourceRefusal &e) {
    ok = e.kind() == kind;
  } catch (const std::exception &) {
  }
  check(ok, name);
}
template <class F> std::string message_of(F f) {
  try {
    f();
  } catch (const std::exception &e) {
    return e.what();
  }
  return {};
}
J measurements = J::object();
// Sanitized runs use reduced sizes (same checks): an 8000-step unanchored
// window and the 1e5-carrier population only (1e6 runs in the normal build).
bool reduced = false;
double seconds_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
      .count();
}
std::unique_ptr<m::PicArena> resident_arena(const m::FullTimeConfiguration &c) {
  return m::PicArena::create_full_time(resident::hall_options(), states(),
                                       channels(), resident::population(),
                                       resident::hall_genesis(), c);
}
m::FullTimeConfiguration with_k(std::uint64_t k) {
  m::FullTimeConfiguration c;
  c.history.anchor_interval_steps = k;
  return c;
}
struct Trace {
  std::vector<std::string> prefixes;
  std::vector<std::uint64_t> revisions;
  std::vector<std::string> rng;
};
resident::StepOptions mixed_events(std::uint64_t k) {
  resident::StepOptions o;
  o.first_events = k % 3 == 1 ? 1 : 0;
  o.second_events = k % 4 == 2 ? 2 : 0;
  return o;
}
Trace run_steps(m::PicArena &a, std::uint64_t steps, bool events) {
  Trace t;
  resident::StepState st;
  for (std::uint64_t k = 0; k < steps; ++k) {
    resident::full_step(a, events ? mixed_events(k) : resident::StepOptions{}, st);
    const auto b = a.full_time_boundary();
    t.prefixes.push_back(b->accepted_arena_sha256);
    t.revisions.push_back(b->accepted_arena_revision);
    t.rng.push_back(b->rng_sha256);
  }
  return t;
}
std::string flip(std::string s, std::size_t at) {
  s[at] = char(s[at] ^ 1);
  return s;
}
void encoding_controls() {
  {
    cn::Writer w;
    w.u8(7);
    w.u64(0x0102030405060708ull);
    w.i64(-5);
    w.f64(-0.0);
    w.f64(std::numeric_limits<double>::denorm_min());
    w.boolean(true);
    w.text("abc");
    w.vec({1.5, -2.25, 1e-300});
    cn::Reader r(w.bytes);
    std::uint8_t a;
    std::uint64_t b;
    std::int64_t c;
    double d, e;
    bool f;
    std::string g;
    m::PicVector h;
    r.u8(a);
    r.u64(b);
    r.i64(c);
    r.f64(d);
    r.f64(e);
    r.boolean(f);
    r.text(g);
    r.vec(h);
    r.finish();
    check(a == 7 && b == 0x0102030405060708ull && c == -5 && d == 0 &&
              std::signbit(d) && e == std::numeric_limits<double>::denorm_min() &&
              f && g == "abc" && h == m::PicVector{1.5, -2.25, 1e-300},
          "canonical primitive round trip (u8/u64/i64/f64 incl -0 and "
          "subnormal/bool/text/vector)");
    check(w.bytes.substr(1, 8) == std::string("\x08\x07\x06\x05\x04\x03\x02\x01", 8),
          "u64 is little-endian");
  }
  {
    cn::Writer x, y;
    x.text("ab");
    x.text("c");
    y.text("a");
    y.text("bc");
    check(x.bytes != y.bytes && std::string("ab") + "c" == std::string("a") + "bc",
          "length prefixes separate (ab,c) from (a,bc); bare concatenation "
          "would collide");
    cn::Writer p, n;
    p.f64(0.0);
    n.f64(-0.0);
    check(p.bytes != n.bytes, "+0.0 and -0.0 have distinct exact encodings");
    cn::Writer t1, t2;
    t1.tag("native-full-time-boundary-v5");
    t1.u64(1);
    t2.tag("native-full-time-common-v5");
    t2.u64(1);
    check(cn::sha256_hex(t1.bytes) != cn::sha256_hex(t2.bytes),
          "record version tags separate equal payloads");
    rejects(
        [&] {
          cn::Reader r(t1.bytes);
          r.tag("native-full-time-common-v5");
        },
        "decoder refuses a foreign record tag");
  }
  {
    cn::Writer w;
    w.u64(5);
    const std::string truncated = w.bytes.substr(0, 7), extended = w.bytes + "x",
                      two(1, '\x02');
    rejects([&] { cn::Reader r(truncated); std::uint64_t v; r.u64(v); },
            "truncated fixed-width field refuses");
    rejects([&] { cn::Reader r(extended); std::uint64_t v; r.u64(v); r.finish(); },
            "trailing bytes after a record refuse");
    rejects([&] { cn::Reader r(two); bool b; r.boolean(b); },
            "boolean byte other than 0/1 refuses");
    check(!std::is_constructible_v<cn::Reader, std::string &&>,
          "decoder cannot be constructed over a temporary byte string");
    cn::Writer l;
    l.u64(1u << 30);
    rejects([&] { cn::Reader r(l.bytes); std::string s; r.text(s, 1024); },
            "declared text length beyond bound refuses before allocation");
    rejects([&] { cn::Reader r(l.bytes); (void)r.count(42); },
            "declared sequence count beyond bound refuses");
  }
}
std::vector<m::GapParticle> rows(std::size_t n, std::uint64_t first = 1) {
  std::vector<m::GapParticle> v(n);
  for (std::size_t i = 0; i < n; ++i) {
    v[i].id = first + i;
    v[i].position_m = {1e-3 + 1e-9 * double(i), 2e-3, 3e-3};
    v[i].proper_velocity_m_s = {double(i), 0, 1e6};
    v[i].charge_C = -m::pic_e;
    v[i].mass_kg = m::pic_me;
    v[i].weight = 4;
  }
  return v;
}
void merkle_controls() {
  auto a = rows(3000);
  auto d = cn::carrier_rows_digest(a);
  check(d->leaves.size() == 3 && d->leaves[2].count == 3000 - 2048 &&
            d->records == 3000,
        "carrier blocks are ID ranges of 1024 (odd leaf count)");
  auto b = a;
  b[1500].position_m[0] = std::nextafter(b[1500].position_m[0], 1.0);
  std::uint64_t rehashed = 0;
  auto inc = cn::carrier_rows_digest(b, &a, d.get(), &rehashed);
  auto full = cn::carrier_rows_digest(b);
  check(inc->root == full->root && rehashed == 1 && inc->root != d->root,
        "incremental digest rehashes only the changed block and equals a full "
        "recompute; one-ulp change alters the root");
  auto e = a;
  e.erase(e.begin() + 10);
  check(cn::carrier_rows_digest(e)->root != d->root,
        "removing one carrier changes the root");
  auto swapped = a;
  std::swap(swapped[3], swapped[4]);
  rejects([&] { cn::carrier_rows_digest(swapped); },
          "unsorted native rows refuse");
  auto dup = a;
  dup[5].id = dup[4].id;
  rejects([&] { cn::carrier_rows_digest(dup); }, "duplicate native ID refuses");
  check(cn::carrier_rows_digest({})->root != cn::carrier_rows_digest(rows(1))->root,
        "empty population has its own root");
  // Host verification of externally computed leaves (device-side design).
  check(cn::verify_carrier_leaves(a, d->leaves, d->root, {0, 2}),
        "host accepts correct leaves, sampled leaves and root");
  auto bad = d->leaves;
  bad[1].hash[0] ^= 1;
  check(!cn::verify_carrier_leaves(a, bad, d->root, {0}),
        "tampered unsampled leaf changes the recomputed root");
  check(!cn::verify_carrier_leaves(b, d->leaves, d->root, {1}),
        "sampled leaf recomputation detects changed rows");
  auto reordered = d->leaves;
  std::swap(reordered[0], reordered[1]);
  check(!cn::verify_carrier_leaves(a, reordered, d->root, {}),
        "reordered leaves refuse");
  // A two-leaf tree cannot be reproduced by a single leaf over the children:
  // leaf and node hashes use different tags.
  auto two = rows(2048);
  auto td = cn::carrier_rows_digest(two);
  std::vector<cn::BlockLeaf> single{{0, 0, 2048, td->leaves[0].hash}};
  check(cn::merkle_root("native-v5/carrier-root", single, 2048) != td->root,
        "leaf/node domain separation");
}

struct ChainParts {
  std::string header;               // tag .. genesis chain, without count
  std::string policy;
  std::uint64_t k = 0;
  std::string genesis;              // boundary + chain bytes
  std::vector<std::string> anchors; // canonical anchor records
};
ChainParts parse_chain(const std::string &chain) {
  ChainParts p;
  cn::Reader in(chain);
  in.tag("native-full-time-anchor-chain-v6");
  in.text(p.policy, 64);
  in.u64(p.k);
  std::uint64_t n = in.count(1u << 20);
  cn::Writer g;
  std::uint64_t u;
  double t;
  std::string txt;
  bool flag;
  for (int i = 0; i < 2; ++i) {
    in.u64(u);
    g.u64(u);
  }
  in.f64(t);
  g.f64(t);
  for (int i = 0; i < 9; ++i) {
    in.text(txt);
    g.text(txt);
  }
  in.boolean(flag);
  g.boolean(flag);
  in.text(txt);
  g.text(txt);
  p.genesis = g.bytes;
  for (std::uint64_t i = 0; i < n; ++i) {
    in.text(txt, 8192);
    p.anchors.push_back(txt);
  }
  in.finish();
  return p;
}
std::string serialize_chain(const ChainParts &p) {
  cn::Writer w;
  w.tag("native-full-time-anchor-chain-v6");
  w.text(p.policy);
  w.u64(p.k);
  w.u64(p.anchors.size());
  w.raw(p.genesis.data(), p.genesis.size());
  for (const auto &a : p.anchors)
    w.text(a);
  return w.bytes;
}
void anchor_controls(const std::string &dir) {
  const std::uint64_t S = 23;
  Trace reference;
  J by_k = J::object();
  std::string chain_k3, chain_k2, root_k3, head_k3, root_k2, head_k2;
  for (std::uint64_t K : {1ull, 2ull, 3ull, 5ull, 8ull, 64ull, 1024ull}) {
    auto a = resident_arena(with_k(K));
    const auto t = run_steps(*a, S, true);
    const auto s = a->full_time_history_status();
    if (K == 1)
      reference = t;
    check(t.prefixes == reference.prefixes &&
              t.revisions == reference.revisions && t.rng == reference.rng,
          "accepted states identical for every anchor interval K");
    check(s.anchors == S / K && s.base_step == (S / K) * K &&
              s.anchor_interval_steps == K && s.cumulative_rows == t.revisions.back(),
          "anchor count, base step and cumulative rows follow K");
    a->validate_full_time_history();
    const auto chain = a->export_full_time_anchor_chain();
    const auto fields = m::verify_full_time_anchor_chain(
        chain, s.anchor_root_sha256, s.anchor_head_sha256);
    bool linked = fields.size() == s.anchors;
    for (std::size_t i = 0; i < fields.size(); ++i)
      linked = linked && fields[i].sequence == i &&
               fields[i].last_step == (i + 1) * K &&
               fields[i].boundary.accepted_arena_sha256 ==
                   t.prefixes[(i + 1) * K - 1] &&
               fields[i].boundary.accepted_arena_revision ==
                   t.revisions[(i + 1) * K - 1] &&
               fields[i].boundary.rng_sha256 == t.rng[(i + 1) * K - 1];
    check(linked, "anchor chain verifies and each anchor commits the exact "
                  "accepted state at its step");
    a->save(dir + "/anchored-K" + std::to_string(K) + ".json");
    refuses_kind(
        [&] {
          m::PicArena::restore_full_time(dir + "/anchored-K" + std::to_string(K) + ".json",
                                         resident::hall_options(), states(),
                                         channels(), *a->full_time_boundary());
        },
        Kind::restore_deferred,
        "nonempty restore validates exported anchors/window then refuses "
        "(deferred)");
    by_k[std::to_string(K)] = {{"anchors", s.anchors},
                               {"anchor_bytes", s.anchor_bytes},
                               {"window_rows", s.window_rows},
                               {"window_bytes", s.window_bytes},
                               {"cumulative_rows", s.cumulative_rows},
                               {"cumulative_bytes", s.cumulative_bytes},
                               {"history_policy_sha256", s.history_policy_sha256}};
    if (K == 3) {
      chain_k3 = chain;
      root_k3 = s.anchor_root_sha256;
      head_k3 = s.anchor_head_sha256;
      // Export file tamper: an altered exported anchor refuses differently.
      auto doc = J::parse(std::ifstream(dir + "/anchored-K3.json"));
      auto &hex = doc["payload"]["metadata"]["history"]["anchors"][1]["canonical_hex"];
      auto text = hex.get<std::string>();
      text[100] = text[100] == '0' ? '1' : '0';
      hex = text;
      doc["sha256"] = sha(doc["payload"].dump());
      std::ofstream(dir + "/anchored-K3-tampered.json") << doc.dump() << '\n';
      const auto msg = message_of([&] {
        m::PicArena::restore_full_time(dir + "/anchored-K3-tampered.json",
                                       resident::hall_options(), states(),
                                       channels(), *a->full_time_boundary());
      });
      check(msg.find("exported anchor digest") != std::string::npos,
            "tampered exported anchor is detected before the deferral refusal");
    }
    if (K == 2) {
      chain_k2 = chain;
      root_k2 = s.anchor_root_sha256;
      head_k2 = s.anchor_head_sha256;
    }
  }
  measurements["anchors_by_K_23_steps_with_events"] = by_k;
  auto refuses = [&](const std::string &x, const std::string &r,
                     const std::string &h) {
    return !message_of([&] { m::verify_full_time_anchor_chain(x, r, h); }).empty();
  };
  check(!refuses(chain_k3, root_k3, head_k3), "untampered chain verifies");
  std::size_t refused = 0, positions = 0;
  for (std::size_t at = 0; at < chain_k3.size(); at += 37, ++positions)
    refused += refuses(flip(chain_k3, at), root_k3, head_k3);
  measurements["single_bit_flip_positions_tested"] = positions;
  check(positions > 50 && refused == positions,
        "every sampled single-bit flip in the serialized chain refuses");
  check(refuses(chain_k3.substr(0, chain_k3.size() - 1), root_k3, head_k3) &&
            refuses(chain_k3 + "x", root_k3, head_k3),
        "truncated or extended chain refuses");
  auto parts = parse_chain(chain_k3);
  check(serialize_chain(parts) == chain_k3 && parts.anchors.size() == 7,
        "chain parse/serialize round trip");
  check(refuses(chain_k3, sha("foreign-root"), head_k3) &&
            refuses(chain_k3, root_k3,
                    m::decode_full_time_anchor(parts.anchors[5]).sha256),
        "foreign root or stale head refuses");
  {
    auto p = parts;
    p.anchors.erase(p.anchors.begin() + 2);
    check(refuses(serialize_chain(p), root_k3, head_k3), "dropped anchor refuses");
    p = parts;
    std::swap(p.anchors[2], p.anchors[3]);
    check(refuses(serialize_chain(p), root_k3, head_k3),
          "reordered anchors refuse");
    p = parts;
    p.anchors.pop_back();
    check(refuses(serialize_chain(p), root_k3,
                  m::decode_full_time_anchor(p.anchors.back()).sha256) == false,
          "a verified prefix of the chain verifies against its own head");
    p = parts;
    p.anchors[3] = flip(p.anchors[3], p.anchors[3].size() - 20);
    check(refuses(serialize_chain(p), root_k3, head_k3),
          "altered anchor field refuses even though its own digest is "
          "recomputed by the verifier");
    auto k2 = parse_chain(chain_k2);
    p = parts;
    p.anchors[6] = k2.anchors[6];
    check(refuses(serialize_chain(p), root_k3, head_k3),
          "anchor spliced from a run with another K refuses");
    check(refuses(chain_k2, root_k3, head_k2),
          "chain of another history policy refuses against this root");
  }
  {
    const auto a = m::decode_full_time_anchor(parts.anchors[4]);
    check(a.sequence == 4 && a.first_step == 12 && a.last_step == 15 &&
              a.anchor_interval_steps == 3 && a.window_motions == 3 &&
              a.sha256 == cn::sha256_hex(parts.anchors[4]) &&
              cn::is_digest(a.native_rows_sha256) &&
              cn::is_digest(a.native_ledger_sha256) &&
              a.boundary.accepted_step == 15,
          "anchor decode round trip exposes step range, digests and boundary");
    rejects([&] { m::decode_full_time_anchor(parts.anchors[4] + "x"); },
            "anchor with trailing bytes refuses");
    rejects([&] { m::decode_full_time_anchor(parts.anchors[4].substr(0, 50)); },
            "truncated anchor refuses");
  }
}
void cap_controls(const std::string &dir) {
  (void)dir;
  {
    auto c = with_k(1);
    c.history.maximum_anchors = 2;
    auto a = resident_arena(c);
    resident::StepState st;
    resident::full_step(*a, {}, st);
    resident::full_step(*a, {}, st);
    const auto prefix = a->state_prefix();
    refuses_kind([&] { resident::full_step(*a, {}, st); }, Kind::anchor_count,
                 "anchor count cap refuses the next step before callbacks");
    check(a->state_prefix() == prefix && a->full_time_history_status().anchors == 2,
          "anchor cap refusal leaves accepted state and chain unchanged");
  }
  {
    auto c = with_k(1);
    c.history.maximum_anchor_bytes = 4096;
    auto a = resident_arena(c);
    resident::StepState st;
    resident::full_step(*a, {}, st);
    const auto prefix = a->state_prefix();
    refuses_kind([&] { resident::full_step(*a, {}, st); }, Kind::anchor_bytes,
                 "anchor byte cap refuses before the step");
    check(a->state_prefix() == prefix, "anchor byte refusal is state-preserving");
  }
  {
    auto c = with_k(4);
    c.history.maximum_window_rows = 20;
    auto a = resident_arena(c);
    resident::StepState st;
    resident::StepOptions o;
    o.first_events = o.second_events = 2;
    std::string prefix;
    Kind kind = Kind::configuration;
    std::uint64_t steps = 0;
    try {
      for (; steps < 8; ++steps) {
        prefix = a->state_prefix();
        resident::full_step(*a, o, st);
      }
    } catch (const m::FullTimeResourceRefusal &e) {
      kind = e.kind();
    }
    const auto s = a->full_time_history_status();
    check(kind == Kind::window_rows && a->state_prefix() == prefix &&
              s.window_rows <= 20 && s.cumulative_rows == a->full_time_boundary()->accepted_arena_revision,
          "window row cap gives a typed refusal; nothing is truncated");
    a->validate_full_time_history();
    check(true, "history remains valid after a window refusal");
  }
  {
    auto c = with_k(100);
    c.history.maximum_window_rows = 499;
    refuses_kind([&] { resident_arena(c); }, Kind::configuration,
                 "window cap below one anchor interval refuses at creation");
    c = with_k(0);
    refuses_kind([&] { resident_arena(c); }, Kind::configuration,
                 "zero anchor interval refuses");
    c = with_k(1);
    c.history.maximum_anchor_bytes = 100;
    refuses_kind([&] { resident_arena(c); }, Kind::configuration,
                 "anchor byte cap below one record refuses");
  }
  {
    auto c = with_k(1024);
    c.limits.maximum_identity_entries = 5;
    auto a = resident_arena(c);
    resident::StepState st;
    resident::StepOptions o;
    o.first_events = 1;
    resident::full_step(*a, o, st);
    const auto prefix = a->state_prefix();
    refuses_kind([&] { resident::full_step(*a, o, st); }, Kind::identity_count,
                 "logical identity cap gives a typed refusal");
    check(a->state_prefix() == prefix, "identity refusal is state-preserving");
  }
  {
    // Iterative release of a long unanchored window (no recursive teardown).
    auto c = with_k(1ull << 20);
    c.history.maximum_window_rows = 6u << 20;
    auto a = resident_arena(c);
    resident::StepState st;
    const auto t0 = std::chrono::steady_clock::now();
    const std::uint64_t steps = reduced ? 8000 : 40000;
    for (std::uint64_t k = 0; k < steps; ++k)
      resident::full_step(*a, {}, st);
    const auto run_s = seconds_since(t0);
    const auto rows = a->full_time_history_status().window_rows;
    const auto t1 = std::chrono::steady_clock::now();
    a.reset();
    measurements["long_window_release"] = {{"steps", steps},
                                           {"window_rows", rows},
                                           {"run_s", run_s},
                                           {"release_s", seconds_since(t1)}};
    check(rows == 5 * steps, "long unanchored window released iteratively");
  }
}
void carrier_controls(const std::string &dir) {
  (void)dir;
  {
    m::FullTimeConfiguration c;
    c.limits.maximum_carriers = 1;
    refuses_kind([&] { resident_arena(c); }, Kind::carrier_cap,
                 "genesis above the configured carrier cap refuses");
    c.limits.maximum_carriers = 1000000;
    c.limits.memory_budget_bytes = 1ull << 30;
    refuses_kind([&] { resident_arena(c); }, Kind::memory_budget,
                 "carrier cap beyond the stated memory budget refuses");
    c.limits.maximum_carriers = m::full_time_absolute_carrier_cap + 1;
    c.limits.memory_budget_bytes = ~0ull;
    refuses_kind([&] { resident_arena(c); }, Kind::carrier_cap,
                 "carrier cap above the absolute cap refuses");
    m::FullTimeConfiguration small, big;
    big.limits.maximum_carriers = 200000;
    big.limits.memory_budget_bytes = 1ull << 32;
    auto a = resident_arena(small), b = resident_arena(big), c2 = resident_arena(small);
    check(a->state_prefix() == c2->state_prefix() &&
              a->state_prefix() != b->state_prefix() &&
              b->full_time_carrier_cap() == 200000,
          "carrier cap and budget are bound into the configuration identity");
    auto k1 = resident_arena(with_k(1)), k9 = resident_arena(with_k(9));
    check(k1->state_prefix() == k9->state_prefix() &&
              k1->full_time_history_status().history_policy_sha256 !=
                  k9->full_time_history_status().history_policy_sha256,
          "anchor interval is recorded in the policy digest but not in the "
          "accepted-state identity");
    const auto bnd = a->full_time_boundary();
    m::FullTimeStepBinding sb{bnd->accepted_step, bnd->accepted_arena_revision,
                              bnd->time_s, bnd->time_s + 5e-15, 5e-15,
                              bnd->accepted_arena_sha256, bnd->field_sha256,
                              bnd->mesh_sha256, bnd->material_sha256,
                              sha("cap-clock")};
    rejects([&] { a->begin_full_time_step(sb, 100001, 32); },
            "outer carrier budget above the configured cap refuses");
  }
  J scaling = J::array();
  for (std::size_t n : {std::size_t(100000), std::size_t(1000000)}) {
    if (reduced && n > 100000) {
      // Not executed in the reduced sanitized run (recorded, not counted):
      // the four 1e6-carrier controls run in the normal build only.
      scaling.push_back({{"carriers", n}, {"skipped_in_reduced_sanitized_run", true}});
      continue;
    }
    m::FullTimeConfiguration c;
    c.limits.maximum_carriers = n;
    c.limits.memory_budget_bytes =
        n * m::full_time_bytes_per_carrier + m::full_time_fixed_bytes;
    c.limits.maximum_identity_entries = 2 * n;
    const auto t0 = std::chrono::steady_clock::now();
    auto a = m::PicArena::create_full_time(resident::hall_options(), states(),
                                           channels(), resident::lattice(n),
                                           resident::hall_genesis(), c);
    const double create_s = seconds_since(t0);
    const long rss_created = resident::status_kib("VmRSS:");
    J row = {{"carriers", n}, {"create_s", create_s},
             {"VmRSS_after_create_KiB", rss_created}};
    resident::StepState st;
    resident::StepOptions o;
    o.maximum_carriers = n;
    if (n <= 100000) {
      const auto t1 = std::chrono::steady_clock::now();
      resident::full_step(*a, o, st);
      row["full_step_s"] = seconds_since(t1);
      const auto t2 = std::chrono::steady_clock::now();
      resident::full_step(*a, o, st);
      row["second_full_step_s"] = seconds_since(t2);
      check(a->full_time_boundary()->accepted_step == 2,
            "1e5 resident carriers: complete genuine steps publish");
    } else {
      // Hall motion v2 refuses a census above its own 1e5 carrier bound; the
      // native owner admits the configured 1e6 cap through census/sealing.
      const auto bnd = a->full_time_boundary();
      m::FullTimeStepBinding b{bnd->accepted_step, bnd->accepted_arena_revision,
                               bnd->time_s, bnd->time_s + 5e-15, 5e-15,
                               bnd->accepted_arena_sha256, bnd->field_sha256,
                               bnd->mesh_sha256, bnd->material_sha256,
                               sha("million-clock")};
      const auto t1 = std::chrono::steady_clock::now();
      auto outer = a->begin_full_time_step(b, n, 32);
      row["begin_s"] = seconds_since(t1);
      m::ExposureBinding e1{b, m::FullTimeSlot::before_transport, b.begin_s,
                            .5 * b.dt_s, 0, b.accepted_field_sha256};
      auto s1 = outer.first_exposure(e1);
      const auto t2 = std::chrono::steady_clock::now();
      auto t = s1.prepare_compound_with_completion(1, [&](auto &) { return completion(e1); });
      check(s1.try_publish(t), "1e6 resident carriers: first exposure publishes");
      row["first_exposure_s"] = seconds_since(t2);
      const auto t3 = std::chrono::steady_clock::now();
      auto trial = outer.fork_transport_trial();
      const auto census = trial.native_census(n);
      row["fork_and_census_s"] = seconds_since(t3);
      check(census.size() == n, "1e6 resident carriers: complete native census");
      const std::vector<double> E(resident::mesh().ne, 0.0), B(resident::mesh().nb, 0.0);
      hall_driver::MotionBinding mb;
      mb.phase = hall_driver::ParticlePhase::full_time_position_and_proper_velocity;
      mb.accepted_step = 0;
      mb.staged_arena_revision = trial.staged_arena_revision();
      mb.begin_s = b.begin_s; mb.end_s = b.end_s; mb.dt_s = b.dt_s;
      mb.staged_arena_sha256 = trial.state_prefix();
      mb.trial_field_sha256 = hall_driver::total_field_sha256(E, B);
      mb.mesh_sha256 = b.mesh_sha256; mb.material_sha256 = b.material_sha256;
      mb.species_sha256 = sha("species"); mb.common_clock_sha256 = b.common_clock_sha256;
      mb.kinetic_policy_sha256 = hall_driver::kinetic_policy_sha256(
          hall_driver::declared_kinetic_policy(a->native_options()));
      const auto msg = message_of([&] {
        (void)hall_driver::prepare_canonical_motion(resident::mesh(), a->native_options(),
                                                    census, mb, {}, E, B);
      });
      row["hall_v2_motion_refusal"] = msg;
      check(msg.find("sealed native carrier cap") != std::string::npos,
            "Hall v2 motion refuses a 1e6 census (provider cap, not native)");
    }
    a->validate_full_time_history();
    row["VmHWM_KiB"] = resident::status_kib("VmHWM:");
    row["budget_bytes"] = c.limits.memory_budget_bytes;
    scaling.push_back(row);
    check(true, "resident carrier population validates at the configured cap");
  }
  measurements["resident_carrier_scaling"] = scaling;
}
int main(int argc, char **argv) {
  try {
    if (argc != 2 && !(argc == 3 && std::string(argv[2]) == "reduced"))
      throw std::runtime_error("output directory [reduced]");
    const std::string dir = argv[1];
    reduced = argc == 3;
    encoding_controls();
    merkle_controls();
    // --- genuine closed step ---------------------------------------------
    {
      auto a = resident_arena({});
      const auto before = *a->full_time_boundary();
      resident::StepState st;
      resident::full_step(*a, {}, st);
      const auto after = *a->full_time_boundary();
      check(after.accepted_step == 1 && after.accepted_arena_revision == 5 &&
                after.time_s == before.time_s + 5e-15 &&
                after.field_sha256 != before.field_sha256 &&
                after.heat_sha256 != before.heat_sha256 &&
                after.rng_sha256 != before.rng_sha256 &&
                after.source_sha256 == before.source_sha256 &&
                after.accepted_arena_sha256 != before.accepted_arena_sha256 &&
                !after.authored_full_time_genesis,
            "complete closed step with genuine Hall v2 motion publishes");
      a->validate_full_time_history();
      check(true, "accepted history revalidates after a genuine step");
      const auto s = a->full_time_history_status();
      check(s.window_rows == 5 && s.cumulative_rows == 5 && s.anchors == 0,
            "five canonical rows per zero-event step");
    }
    // Motion rows hold digests: size independent of carrier count; the full
    // record stays with the caller's motion token for the step.
    {
      auto measure = [&](std::size_t extra) {
        auto pop = resident::population();
        for (std::size_t i = 0; i < extra; ++i)
          pop.push_back({"extra/" + std::to_string(i),
                         "n",
                         "manufactured/source",
                         {2, 0},
                         {.0030 + 1e-5 * double(i), .005, 0},
                         {0, 0, 0}});
        auto a = m::PicArena::create_full_time(resident::hall_options(), states(),
                                               channels(), pop,
                                               resident::hall_genesis());
        resident::StepState st;
        resident::StepOptions o;
        o.maximum_carriers = 64;
        resident::full_step(*a, o, st);
        a->save(dir + "/motion-" + std::to_string(extra) + ".json");
        const auto doc = J::parse(std::ifstream(dir + "/motion-" + std::to_string(extra) + ".json"));
        for (const auto &r : doc["payload"]["metadata"]["history"]["window"])
          if (r["kind"] == "full_time_motion")
            return r["canonical_hex"].get<std::string>().size() / 2;
        return std::size_t(0);
      };
      const auto small = measure(0), large = measure(40);
      measurements["motion_row_bytes"] = {{"carriers_2", small}, {"carriers_42", large}};
      check(small > 0 && small == large,
            "journal motion row stores digests only (size independent of "
            "carrier count)");
      auto a = resident_arena({});
      const auto b0 = a->full_time_boundary();
      m::FullTimeStepBinding b{b0->accepted_step, b0->accepted_arena_revision,
                               b0->time_s,        b0->time_s + 5e-15,
                               5e-15,             b0->accepted_arena_sha256,
                               b0->field_sha256,  b0->mesh_sha256,
                               b0->material_sha256, sha("audit-clock")};
      auto outer = a->begin_full_time_step(b, 8, 32);
      m::ExposureBinding e1{b, m::FullTimeSlot::before_transport, b.begin_s,
                            .5 * b.dt_s, 0, b.accepted_field_sha256};
      auto s1 = outer.first_exposure(e1);
      auto t1 = s1.prepare_compound_with_completion(1, [&](auto &) {
        auto c = completion(e1);
        return c;
      });
      check(s1.try_publish(t1), "audit fixture first exposure");
      auto trial = outer.fork_transport_trial();
      const auto census = trial.native_census(8);
      namespace h = hall_driver;
      const std::vector<double> E(resident::mesh().ne, 0.0), B(resident::mesh().nb, 0.0);
      h::MotionBinding mb;
      mb.phase = h::ParticlePhase::full_time_position_and_proper_velocity;
      mb.accepted_step = b.accepted_step;
      mb.staged_arena_revision = trial.staged_arena_revision();
      mb.begin_s = b.begin_s;
      mb.end_s = b.end_s;
      mb.dt_s = b.dt_s;
      mb.staged_arena_sha256 = trial.state_prefix();
      mb.trial_field_sha256 = h::total_field_sha256(E, B);
      mb.mesh_sha256 = b.mesh_sha256;
      mb.material_sha256 = b.material_sha256;
      mb.species_sha256 = sha("species");
      mb.common_clock_sha256 = b.common_clock_sha256;
      mb.kinetic_policy_sha256 = h::kinetic_policy_sha256(h::declared_kinetic_policy(a->native_options()));
      auto receipt = h::prepare_canonical_motion(resident::mesh(), a->native_options(), census, mb, {}, E, B);
      check(receipt.census_sha256() == h::canonical_census_sha256(census),
            "Hall canonical census digest of the native census");
      auto moved = trial.prepare_total_field_motion(
          {b, receipt.seal_sha256(), b.accepted_field_sha256, sha("f1"), sha("h1"), b0->source_sha256},
          m::bind_hall_total_field_motion(receipt));
      const auto &audit = moved.audit();
      check(audit.at("receipt").at("owner_audit").at("endpoints").size() == 2 &&
                audit.at("endpoints").size() == 2 &&
                audit.at("receipt").at("owner_audit").at("seal_sha256") == receipt.seal_sha256(),
            "complete motion record remains available to the caller for the step");
      // Native census digest admission equals Hall's: a foreign census refuses.
      auto wrong = mb;
      wrong.staged_arena_sha256 = sha("other-prefix");
      auto other = h::prepare_canonical_motion(resident::mesh(), a->native_options(), census, wrong, {}, E, B);
      auto trial2 = outer.fork_transport_trial();
      rejects([&] {
        trial2.prepare_total_field_motion(
            {b, other.seal_sha256(), b.accepted_field_sha256, sha("f1"), sha("h1"), b0->source_sha256},
            m::bind_hall_total_field_motion(other));
      }, "genuine Hall receipt bound to another staged prefix refuses");
      outer.accept_transport_trial(std::move(moved));
      check(true, "genuine motion accepted after audit export");
    }
    anchor_controls(dir);
    cap_controls(dir);
    carrier_controls(dir);
    check(forbidden_calls == 0,
          "no native solver/advance call in any v5 control");
    J output = {{"status", "passed"},
                {"check_count", checks},
                {"checks", results},
                {"measurements", measurements},
                {"forbidden_native_solver_or_advance_calls", forbidden_calls},
                {"GPU", false},
                {"motion_provider",
                 "genuine Hall v2 prepare_canonical_motion at identically zero "
                 "total field; no field solve"},
                {"disk_restore_of_nonempty_history", "deferred (typed refusal)"},
                {"reduced_sanitizer_sizes", reduced}};
    std::ofstream(dir + "/v5-history-tests.json") << output.dump(2) << '\n';
    std::cout << checks << " native v5 history/encoding controls passed\n";
    return 0;
  } catch (const std::exception &x) {
    std::cerr << x.what() << '\n';
    return 1;
  }
}
