#include "electrode_heat_feedback.hpp"
#include "heat-v2-count.hpp"
#include "raikou/mfem_pic_arena.h"
#include "water_parent_bridge.hpp"
#include <climits>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <openssl/sha.h>
#include <sstream>
#include <type_traits>
namespace m = raikou::mfem_backend;
namespace w = raikou::pic_water;
namespace pc = raikou::pic_collisions;
using J = m::ArenaJson;
J checks = J::array(), proofs = J::array();
void check(bool ok, const std::string &name) {
  checks.push_back({{"name", name}, {"passed", ok}});
  if (!ok)
    throw std::runtime_error(name);
}
template <class F> void reject(F fn, const std::string &name) {
  bool caught = false;
  try {
    fn();
  } catch (const std::exception &) {
    caught = true;
  }
  check(caught, name);
}
std::string sha(const std::string &bytes) {
  unsigned char d[SHA256_DIGEST_LENGTH];
  SHA256(reinterpret_cast<const unsigned char *>(bytes.data()), bytes.size(),
         d);
  std::ostringstream o;
  for (auto x : d)
    o << std::hex << std::setfill('0') << std::setw(2) << unsigned(x);
  return o.str();
}
constexpr double mn = 2.9915e-26;
const std::string datahash =
    "840663cdbcf4c39f1ab6f5407f9ecacb232a07c3c48002a4a182c8187792e4bf";
const std::string manifest =
    "0394f9fd49be9753d2258652e6e5bf5e489d39dda43b5ccfa95fd7b34ec67b60";
const std::string bindingfile =
    "f9de77cd4c52dbf85069994d17c3808196ee30333116472ae98a45540d1d453f";
m::GapPicOptions options() {
  m::GapPicOptions o;
  o.nx = o.ny = 45;
  o.width_m = o.gap_m = .009;
  o.depth_m = .025;
  o.finite_foils = true;
  o.finite_axial_extent = true;
  o.short_z_m = 0;
  o.feed_z_m = .025;
  o.lower_foil = {.002, .007, .0018, .002};
  o.upper_foil = {.002, .007, .007, .0072};
  return o;
}
std::vector<m::ArenaStateSpec> states() {
  return {{"electron", 0, 0, -1, w::me, 0, 0},
          {"H2O", 2, 1, 0, mn, 0, 0},
          {"H2O+", 2, 1, 1,
           mn - w::me + w::threshold_eV * pc::qe / (pc::c * pc::c), 0,
           w::threshold_eV * pc::qe}};
}
std::vector<m::ArenaChannel> channels() {
  return {{w::parent_channel,
           "electron",
           "H2O",
           {"electron", "electron", "H2O+"},
           datahash,
           manifest,
           bindingfile,
           w::parent_model}};
}
m::PicVector eu(double E) {
  double x = E * pc::qe / (w::me * pc::c * pc::c);
  return {pc::c * std::sqrt(x * (x + 2)), 0, 0};
}
struct Case {
  m::GapPic pic{options()};
  std::unique_ptr<m::PicArena> arena;
  Case(double ew, double nw) {
    m::GapParticle e, n;
    e.weight = ew;
    e.mass_kg = w::me;
    e.position_m = {.004, .004, .0125};
    e.proper_velocity_m_s = eu(100);
    n.weight = nw;
    n.mass_kg = mn;
    n.charge_C = 0;
    n.position_m = {.00403, .004, .0125};
    auto ids = pic.add_authored_batch({e, n});
    pic.solve_fields();
    arena = std::make_unique<m::PicArena>(
        pic, states(), channels(),
        std::vector<m::ArenaPacket>{
            {"e", "electron", "authored_fixture", "", {ew, 0}, {ids[0]}},
            {"n", "H2O", "authored_fixture", "", {nw, 0}, {ids[1]}}},
        "test/initial-joint-state");
  }
  m::ArenaEvent event(const w::ParentBridge &bridge) {
    auto e = arena->view("e"), n = arena->view("n");
    pc::Packet ep{"e",      "electron", 920,
                  e.weight, e.mass_kg,  e.proper_velocity_m_s},
        np{"n", "H2O", 920, n.weight, n.mass_kg, n.proper_velocity_m_s};
    const auto candidates =
        pc::candidates({ep}, {np}, {{920, (.009 / 45) * (.009 / 45) * .025}},
                       {bridge.channel()});
    check(candidates.size() == 1, "actual strict-support100eV candidate");
    auto out = bridge.sample(candidates[0], ep, np, {0, 1, 0});
    m::ArenaEvent event;
    event.id = "conditional/parent/1";
    event.channel = w::parent_channel;
    event.primary = "e";
    event.target = "n";
    event.scattered_primary_id = "e";
    event.event_weight = candidates[0].event_weight;
    if (ew_greater(candidates[0]))
      event.untouched_primary_id = "e/remainder";
    event.phase.state_prefix = arena->state_prefix();
    event.phase.mesh_prefix = "native/45x45-fixed";
    event.phase.clock_prefix = "joint/initial";
    event.provenance = {m::ArenaProvenance::Kind::ConditionalPhysicalData,
                        datahash,
                        manifest,
                        bindingfile,
                        w::parent_model,
                        "forced_kernel_integration_only",
                        event.phase.state_prefix,
                        3e-14};
    for (int i = 0; i < 3; ++i)
      event.products.push_back({i == 0   ? "e"
                                : i == 1 ? "secondary"
                                         : "water_ion",
                                out.particles[i].state,
                                out.particles[i].proper_velocity_m_s,
                                i == 0 ? m::ArenaProduct::At::Primary
                                       : m::ArenaProduct::At::Target});
    return event;
  }
  bool ew_greater(const pc::Candidate &c) {
    return c.projectile_before > c.event_weight;
  }
};
int main(int argc, char **argv) {
  try {
    if (argc != 3)
      throw std::runtime_error("inputdataset and reportdir");
    const std::string dir = argv[2];
    auto data = std::make_shared<const w::Data>(w::Data::read(argv[1]));
    w::ParentBridge bridge(data, w::parent_model, datahash, manifest, mn);
    for (double wgt : {.25, .1, 1e-20, 1e-300}) {
      auto count = m::PicArena::debit({1, 0}, wgt);
      auto carriers = m::PicArena::positive_carriers(count);
      check(!carriers.empty() && carriers.size() <= 42,
            "positive native limbs for exact compensated remainder");
      check(m::PicArena::exact_count_equal(count.credited(wgt), {1, 0}),
            "count exact debit/add restoration");
    }
    {
      auto exact = m::PicArena::debit({1, -1e-20}, 1e-100);
      check(exact.credited(1e-100) == m::ArenaCount(1, -1e-20),
            "bounded exact count retains third residual");
      reject([&] { exact.require_two_limb_values(); },
             "lossy export back to Count2 explicitly refused");
    }
    reject([&] { m::PicArena::debit({1, -1e-20}, 1.); },
           "rounded rate weight above exact inventory refused");
    for (const auto weights :
         {std::pair<double, double>{.25, 64}, {4, .25}, {.25, 1e20}, {.1, 1}}) {
      Case s(weights.first, weights.second);
      auto e = s.event(bridge);
      auto old = m::PicArena::fingerprint(s.pic);
      auto prefix = s.arena->state_prefix();
      auto source = s.pic.diagnostics();
      for (int fault : {1, 2, 3}) {
        reject([&] { s.arena->prepare(e, fault); },
               "fault-injected actual native event prepare");
        check(old == m::PicArena::fingerprint(s.pic) &&
                  prefix == s.arena->state_prefix(),
              "event failure leaves native/arena unchanged");
      }
      auto bad = e;
      bad.provenance.kind = m::ArenaProvenance::Kind::ManufacturedForced;
      bad.provenance.model = "manufactured/fake";
      reject([&] { s.arena->prepare(bad); },
             "physical data cannot masquerade as manufactured");
      bad = e;
      bad.provenance.dataset_sha256 = std::string(64, 'a');
      reject([&] { s.arena->prepare(bad); },
             "foreign physical data digest refused");
      bad = e;
      bad.phase.physical_stage_s = bad.provenance.collision_elapsed_s;
      reject([&] { s.arena->prepare(bad); },
             "collision subflow clock cannot become native trajectory clock");
      bad = e;
      bad.products[1].proper_velocity_m_s[0] += 1000;
      reject([&] { s.arena->prepare(bad); },
             "bad energy or momentum rejected before publication");
      bad = e;
      bad.event_weight = weights.second * 2;
      reject([&] { s.arena->prepare(bad); },
             "insufficient target or primary refused");
      auto t = s.arena->prepare(e), other = s.arena->prepare(e);
      check(s.arena->try_publish(t),
            "actual parent products atomically published");
      check(!s.arena->try_publish(other), "stale prepared event refused");
      reject([&] { s.arena->prepare(e); }, "accepted event replay refused");
      check(s.pic.diagnostics().authored_charge_C == source.authored_charge_C &&
                s.pic.diagnostics().authored_kinetic_energy_J ==
                    source.authored_kinetic_energy_J,
            "reaction preserves external authored C/K ledgers");
      check(s.pic.diagnostics().time_s == 0,
            "native frozen physical clock unchanged by SSA elapsed time");
      check(s.arena->snapshot()->packets.at("water_ion").count ==
                m::ArenaCount(e.event_weight, 0),
            "matched represented product weights");
      if (weights.first > weights.second) {
        auto remainder = s.arena->view("e/remainder");
        check(remainder.proper_velocity_m_s == eu(100) &&
                  remainder.weight == 3.75,
              "untouched primary remainder preserves exact pre-collision "
              "u/count");
      }
      if (weights.second == 1e20) {
        const auto &n = s.arena->snapshot()->packets.at("n");
        check(n.count.require_two_limb_values()[0] == 1e20 &&
                  n.count.require_two_limb_values()[1] == -.25 &&
                  n.native_ids.size() > 1,
              "subULP neutral debit retained in canonical Count AND native "
              "carriers");
      }
      proofs.push_back(t.audit());
      const std::string bundle =
          dir + "/case-" + std::to_string(proofs.size()) + ".json";
      s.arena->save(bundle);
      m::GapPic restored(options());
      auto a = m::PicArena::restore(restored, bundle, states(), channels(),
                                    "test/initial-joint-state");
      check(a->snapshot()->journal == s.arena->snapshot()->journal &&
                a->snapshot()->event_ids == s.arena->snapshot()->event_ids,
            "typed event provenance/replay prefix persists restart");
      check(a->view("water_ion").proper_velocity_m_s ==
                s.arena->view("water_ion").proper_velocity_m_s,
            "native product state persists restart");
      reject([&] { a->prepare(e); }, "accepted event refused after restart");
      for (int corruption = 0; corruption < 12; ++corruption) {
        std::ifstream in(bundle);
        J doc;
        in >> doc;
        auto &body = doc["payload"];
        if (corruption == 0)
          body["source_sha256"] = std::string(64, 'a');
        if (corruption == 1)
          body["ABI"][0] = 81;
        if (corruption == 2)
          body["registry"]["states"][0]["mass_kg"] = 1.;
        if (corruption == 3)
          body["metadata"]["packets"][0]["count"]["positive_carriers"] =
              J::array({2.});
        if (corruption == 4)
          body["metadata"]["packets"][0]["native_ids"][0] = 99999;
        if (corruption == 5)
          body["native_bytes_sha256"] = std::string(64, 'b');
        if (corruption == 6) {
          auto &row = body["metadata"]["journal"][0];
          row["payload"]["dataset_sha256"] = std::string(64, 'a');
          row["sha256"] = sha(row["payload"].dump());
          body["metadata"]["chain"] = row["sha256"];
        }
        if (corruption == 7)
          body["metadata"]["event_ids"] = J::array();
        if (corruption >= 8) {
          auto &row = body["metadata"]["journal"][0];
          if (corruption == 8)
            row["payload"]["event_weight"] = 123.;
          if (corruption == 9)
            row["payload"]["target_remainder"] = {-1., 0.};
          if (corruption == 10)
            row["payload"]["selection_method"] = "unsupported_method";
          if (corruption == 11)
            row["payload"]["collision_elapsed_s"] = -1.;
          row["sha256"] = sha(row["payload"].dump());
          body["metadata"]["chain"] = row["sha256"];
        }
        doc["sha256"] = sha(body.dump());
        const auto file = dir + "/corrupt.json";
        {
          std::ofstream o(file);
          o << doc.dump();
        }
        m::GapPic untouched(options());
        auto original = m::PicArena::fingerprint(untouched);
        reject(
            [&] {
              m::PicArena::restore(untouched, file, states(), channels(),
                                   "test/initial-joint-state");
            },
            "rehash corruption refused " + std::to_string(corruption));
        check(original == m::PicArena::fingerprint(untouched),
              "restart refusal leaves fresh native state exact");
      }

      check(s.arena->try_rollback(t) && s.arena->state_prefix() == prefix &&
                m::PicArena::fingerprint(s.pic) == old,
            "joint rollback restores native and canonical inventory");
      check(!s.arena->try_publish(t), "rolled-back token cannot be reused");
    }
    {
      Case s(.25, 64);
      const double step = 1e-14;
      auto transport = s.arena->prepare_transport(step);
      check(s.arena->try_publish(transport),
            "preceding native step establishes synchronized phase");
      auto e = s.event(bridge);
      e.phase.kind = m::ArenaPhase::Kind::ProvisionalEqualDtMidX;
      e.phase.accepted_step = 1;
      e.phase.physical_stage_s = s.pic.diagnostics().time_s;
      e.phase.previous_dt_s = step;
      e.phase.next_dt_s = step;
      auto bad = e;
      bad.phase.next_dt_s = 2 * step;
      reject([&] { s.arena->prepare(bad); },
             "variable dt collision placement explicitly unqualified");
      auto t = s.arena->prepare(e);
      check(
          s.arena->try_publish(t),
          "actual kernel event accepted at explicit equal-dt native boundary");
      proofs.push_back(t.audit());
    }
    {
      Case s(.25, 64);
      auto old = m::PicArena::fingerprint(s.pic),
           prefix = s.arena->state_prefix();
      int callbacks = 0;
      reject(
          [&] { s.arena->prepare_compound(0, [&](auto &) { ++callbacks; }); },
          "zero compound cap refused before builder");
      reject(
          [&] {
            s.arena->prepare_compound(4097, [&](auto &) { ++callbacks; });
          },
          "unbounded compound cap refused before builder");
      check(callbacks == 0, "resource admission precedes compound callback");
      auto make = [&](int n) {
        auto e = s.event(bridge);
        e.id = "conditional/compound/" + std::to_string(n);
        e.products[1].id = "secondary/" + std::to_string(n);
        e.products[2].id = "water_ion/" + std::to_string(n);
        return e;
      };
      for (int failure : {0, 1, 2, 3}) {
        reject(
            [&] {
              s.arena->prepare_compound(1, [&](m::PicArena &a) {
                auto t = a.prepare(make(1));
                check(a.try_publish(t), "private first reaction accepted");
                if (failure == 0)
                  throw std::runtime_error("manufactured builder failure");
                if (failure == 1)
                  a.prepare_compound(1, [](auto &) {});
                if (failure == 2)
                  a.prepare_transport(1e-14);
                if (failure == 3)
                  a.prepare(make(2));
              });
            },
            "compound exception/nesting/transport/cap refusal");
        check(m::PicArena::fingerprint(s.pic) == old &&
                  s.arena->state_prefix() == prefix &&
                  s.arena->snapshot()->event_ids.empty(),
              "compound failure restores complete native and journal/count "
              "prefix");
      }
      std::unique_ptr<m::PreparedArena> escaped_prepared, escaped_published;
      auto compound = s.arena->prepare_compound(2, [&](m::PicArena &a) {
        auto first = make(1);
        escaped_prepared = std::make_unique<m::PreparedArena>(a.prepare(first));
        escaped_published =
            std::make_unique<m::PreparedArena>(a.prepare(first));
        check(a.try_publish(*escaped_published), "private subevent1 accepted");
        auto t = a.prepare(make(2));
        check(a.try_publish(t),
              "private subevent2 accepted with recomputed native input");
        reject([&] { a.save(dir + "/forbidden-private-restart.json"); },
               "private compound cannot publish checkpoint");
      });
      check(s.arena->state_prefix() == prefix &&
                m::PicArena::fingerprint(s.pic) == old &&
                s.arena->snapshot()->event_ids.empty(),
            "successful compound prepare restores base before return");
      check(!s.arena->try_publish(*escaped_prepared) &&
                !s.arena->try_rollback(*escaped_published),
            "escaped child prepared/published tokens are inert outside unique "
            "scope");
      check(s.arena->try_publish(compound) &&
                s.arena->snapshot()->event_ids.size() == 2,
            "single compound publication accepts both reactions");
      proofs.push_back(compound.audit());
      s.arena->save(dir + "/compound.json");
      m::GapPic fresh(options());
      auto restored =
          m::PicArena::restore(fresh, dir + "/compound.json", states(),
                               channels(), "test/initial-joint-state");
      check(restored->snapshot()->event_ids.size() == 2 &&
                restored->snapshot()->journal.size() == 2,
            "compound accepted event prefix persists restart");
      check(s.arena->try_rollback(compound) &&
                s.arena->state_prefix() == prefix &&
                m::PicArena::fingerprint(s.pic) == old,
            "single joint rollback restores entire compound");
      check(!s.arena->try_publish(compound),
            "rolled-back compound cannot replay");
      auto empty = s.arena->prepare_compound(1, [](auto &) {});
      check(s.arena->try_publish(empty) &&
                s.arena->snapshot()->event_ids.empty(),
            "zero selected events can publish an unchanged collision subflow");
      check(s.arena->native_options().nx == 45 &&
                s.arena->native_options().short_z_m == 0,
            "collision geometry is read from owning native arena");
    }
    // Exact native neutral transport, all three positive limbs share real
    // drift.
    m::GapPic pic(options()), reference(options());
    std::vector<m::GapParticle> particles;
    std::vector<m::ArenaPacket> records;
    auto append = [&](std::string id, m::ArenaCount count, m::PicVector x,
                      m::PicVector u) {
      m::ArenaPacket record{id, "H2O", "authored_fixture", "", count, {}};
      for (double weight : m::PicArena::positive_carriers(count)) {
        m::GapParticle p;
        p.position_m = x;
        p.proper_velocity_m_s = u;
        p.weight = weight;
        p.mass_kg = mn;
        p.charge_C = 0;
        particles.push_back(p);
      }
      records.push_back(record);
    };
    append("wall", {1, -1e-20}, {.004, .0021, .0125}, {0, -1e6, 0});
    append("escape", {2, 0}, {.008, .008, .02499}, {0, 0, 1e6});
    append("resident", {3, 0}, {.001, .005, .0125}, {1e5, 0, 0});
    auto ids = pic.add_authored_batch(particles);
    auto refids = reference.add_authored_batch(particles);
    check(ids == refids, "neutral reference native identity match");
    std::size_t ix = 0;
    for (auto &record : records) {
      auto count = m::PicArena::positive_carriers(record.count).size();
      for (std::size_t k = 0; k < count; ++k)
        record.native_ids.push_back(ids[ix++]);
    }
    pic.solve_fields();
    reference.solve_fields();
    m::PicArena arena(pic, states(), channels(), records, "neutral/initial");
    const double dt = 3e-10;
    auto old = m::PicArena::fingerprint(pic), prefix = arena.state_prefix();
    for (int fault : {1, 2, 3}) {
      reject([&] { arena.prepare_transport(dt, fault); },
             "native transport prepare fault");
      check(m::PicArena::fingerprint(pic) == old &&
                arena.state_prefix() == prefix,
            "transport failure rollback exact");
    }
    auto transport = arena.prepare_transport(dt);
    reference.advance(dt);
    check(m::PicArena::fingerprint(pic) == old,
          "native transport prepare unpublished");
    check(arena.try_publish(transport), "native transport published");
    check(m::PicArena::fingerprint(pic) == m::PicArena::fingerprint(reference),
          "stack observer leaves native drift/force/deposit/Poisson arithmetic "
          "exact");
    check(arena.snapshot()->packets.size() == 1 &&
              arena.snapshot()->packets.count("resident"),
          "actual native q0 collection/escape depletes canonical species");
    check(arena.snapshot()->sinks.size() == 2, "species-resolved sinks once");
    check(m::ArenaCount::from_positive_carriers(
              arena.snapshot()
                  ->sinks[0]["count"]["positive_carriers"]
                  .get<std::vector<double>>()) == m::ArenaCount(2, 0) ||
              m::ArenaCount::from_positive_carriers(
                  arena.snapshot()
                      ->sinks[1]["count"]["positive_carriers"]
                      .get<std::vector<double>>()) == m::ArenaCount(2, 0),
          "exact escaped neutral number receipt");
    bool wall = false;
    for (auto receipt : arena.snapshot()->sinks)
      if (receipt["packet_id"] == "wall") {
        wall = true;
        check(m::ArenaCount::from_positive_carriers(
                  receipt["count"]["positive_carriers"]
                      .get<std::vector<double>>()) ==
                      m::ArenaCount(1, -1e-20) &&
                  receipt["sink"] == "lower_foil" &&
                  receipt["wall_heat_assumed"] == false,
              "neutral wall receipt preserves two-word count and distinguishes "
              "incident K from wall heat");
      }
    check(wall, "actual wall sink received");
    proofs.push_back(transport.audit());
    arena.save(dir + "/neutral.json");
    m::GapPic restored(options());
    auto restored_arena =
        m::PicArena::restore(restored, dir + "/neutral.json", states(),
                             channels(), "neutral/initial");
    check(restored_arena->snapshot()->sinks == arena.snapshot()->sinks,
          "species/reference sink provenance survives restart");
    auto continuation = restored_arena->prepare_transport(1e-12);
    check(restored_arena->try_publish(continuation),
          "native neutral continuation after restart");
    check(arena.try_rollback(transport) && arena.state_prefix() == prefix,
          "transport joint rollback exact");
    static_assert(!std::is_copy_constructible_v<m::PicArena> &&
                  !std::is_move_constructible_v<m::PicArena>);
    check(!std::is_copy_assignable_v<m::PicArena> &&
              !std::is_move_assignable_v<m::PicArena>,
          "exclusive native owner cannot be copied/moved");
    for (int mutation = 0; mutation < 4; ++mutation) {
      std::ifstream in(dir + "/neutral.json");
      J doc;
      in >> doc;
      auto &body = doc["payload"];
      if (mutation == 0)
        body["metadata"]["sinks"][0]["count"]["positive_carriers"] =
            J::array({1e30});
      if (mutation == 1)
        body["metadata"]["sinks"][0]["state"] = "electron";
      if (mutation == 2)
        body["metadata"]["sinks"][0]["reference_per_particle_J"] = 1e20;
      if (mutation == 3)
        body["metadata"]["sinks"].push_back(
            {{"fake", "not a species receipt"}});
      doc["sha256"] = sha(body.dump());
      const auto file = dir + "/corrupt-sink.json";
      {
        std::ofstream out(file);
        out << doc.dump();
      }
      m::GapPic fresh(options());
      const auto before = m::PicArena::fingerprint(fresh);
      reject(
          [&] {
            m::PicArena::restore(fresh, file, states(), channels(),
                                 "neutral/initial");
          },
          "inconsistent rehashed sink/catalog/number/reference refused");
      check(m::PicArena::fingerprint(fresh) == before,
            "failed sink history restore leaves native unchanged");
    }
    {
      m::GapPic p(options());
      m::GapParticle n;
      n.position_m = {.004, .004, .0125};
      n.mass_kg = mn;
      n.charge_C = 0;
      auto ids = p.add_authored_batch({n, n});
      p.solve_fields();
      reject(
          [&] {
            m::PicArena bad(
                p, {{"A", INT_MAX, 0, 0, mn, 0, 0}, {"B", 2, 0, 0, mn, 0, 0}},
                {{"manufactured/overflow",
                  "A",
                  "A",
                  {"A", "A", "A", "A", "B"},
                  "",
                  "",
                  "",
                  ""}},
                {{"a", "A", "authored", "", {1, 0}, {ids[0]}},
                 {"b", "A", "authored", "", {1, 0}, {ids[1]}}},
                "integer/initial");
          },
          "wide stoichiometry rejects2^32 atom mismatch");
      reject(
          [&] {
            m::PicArena bad(p, {{"A", 0, 0, INT_MIN, mn, 0, 0}}, {}, {},
                            "integer/initial");
          },
          "INT_MIN charge refused without abs undefined behavior");
      auto ch = channels();
      ch[0].product_states = std::vector<std::string>(9, "electron");
      reject([&] { m::PicArena bad(p, states(), ch, {}, "integer/initial"); },
             "catalog product list capped at event bound8");
    }
    {
      m::GapPic p(options());
      m::GapParticle n;
      n.position_m = {.004, .004, .0125};
      n.mass_kg = w::me;
      n.charge_C = -m::pic_e;
      n.weight = std::numeric_limits<double>::denorm_min();
      n.proper_velocity_m_s = eu(100);
      auto ids = p.add_authored_batch({n});
      p.solve_fields();
      reject(
          [&] {
            m::PicArena a(
                p, states(), channels(),
                {{"tiny", "electron", "authored", "", {n.weight, 0}, ids}},
                "tiny/initial");
          },
          "charged carrier with zero represented native charge/energy refused "
          "unchanged");
    }
    {
      Case c(.25, 64);
      reject([&] { c.arena->rate_population(1); },
             "bulk census rejects cap before allocation/truncation");
      auto bulk = c.arena->rate_population(100);
      std::vector<std::uint64_t> expected;
      for (const auto &[id, p] : c.arena->snapshot()->packets) {
        for (auto native : p.native_ids)
          expected.push_back(native);
      }
      std::vector<std::uint64_t> actual;
      for (const auto &row : bulk) {
        actual.push_back(row.native.id);
        check(c.arena->snapshot()->packets.at(row.packet_id).state == row.state,
              "bulk rate entry has canonical state");
      }
      check(actual == expected && bulk.size() == c.pic.particles().size(),
            "one validated bulk rate population preserves canonical/native "
            "carrier order");
    }

    {
      Case c(.25, 64);
      const auto native = m::PicArena::fingerprint(c.pic),
                 initial = c.arena->state_prefix();
      m::ArenaStochasticCompletion record;
      record.integration_owner_id = "waterflow/conditional-parent";
      record.before_rng_clock_sha256 = std::string(64, 'a');
      record.after_rng_clock_sha256 = std::string(64, 'b');
      record.common_clock_sha256 = std::string(64, '9');
      record.channel_id = w::parent_channel;
      record.dataset_sha256 = datahash;
      record.kernel_manifest_sha256 = manifest;
      record.kernel_binding_file_sha256 = bindingfile;
      record.model = w::parent_model;
      record.stage_sequence = 0;
      record.frozen_physical_stage_s = 0;
      record.common_interval_begin_s = 0;
      record.common_interval_end_s = 1e-12;
      const auto &bound = c.arena->channel_binding(w::parent_channel);
      check(bound.dataset_sha256 == record.dataset_sha256 &&
                bound.kernel_manifest_sha256 == record.kernel_manifest_sha256 &&
                bound.model == record.model,
            "no-event constructor can verify exact immutable channel binding");
      for (int invalid = 0; invalid < 3; ++invalid) {
        auto bad = record;
        if (invalid == 0)
          bad.common_interval_begin_s = -1e-12;
        if (invalid == 1) {
          bad.common_interval_begin_s = 1e-12;
          bad.common_interval_end_s = 2e-12;
        }
        if (invalid == 2)
          bad.common_clock_sha256 = "";
        reject(
            [&] {
              c.arena->prepare_compound_with_completion(
                  1, [&](auto &) { return bad; });
            },
            "standalone completion rejects negative/outside/unbound common "
            "clock");
        check(c.arena->state_prefix() == initial &&
                  m::PicArena::fingerprint(c.pic) == native,
              "invalid common clock completion leaves all state unchanged");
      }
      int calls = 0;
      reject(
          [&] {
            c.arena->prepare_compound_with_completion(4096, [&](auto &) {
              ++calls;
              return record;
            });
          },
          "completion reserves bounded terminal-history entry before callback");
      check(calls == 0, "history budget fails before stochastic builder");
      auto zero = c.arena->prepare_compound_with_completion(1, [&](auto &) {
        ++calls;
        return record;
      });
      check(c.arena->state_prefix() == initial &&
                !c.arena->stochastic_completion(record.integration_owner_id),
            "zero-event completion prepare restores native/owner base");
      check(c.arena->try_publish(zero),
            "zero-event stochastic exposure atomically accepted");
      auto prefix0 = c.arena->state_prefix();
      check(prefix0 != initial && m::PicArena::fingerprint(c.pic) == native &&
                c.arena->snapshot()->event_ids.empty(),
            "accepted zero-event exposure changes arena prefix without native "
            "motion");
      auto accepted =
          c.arena->stochastic_completion(record.integration_owner_id);
      check(accepted && accepted->stage_sequence == 0 &&
                accepted->after_rng_clock_sha256 == std::string(64, 'b'),
            "typed completion binds accepted RNG+clock payload and sequence");
      check(accepted->after_rng_clock_sha256 != record.before_rng_clock_sha256,
            "mixed-age component payload would fail root receipt binding");
      c.arena->save(dir + "/completion-zero.json");
      m::GapPic fresh(options());
      auto restored =
          m::PicArena::restore(fresh, dir + "/completion-zero.json", states(),
                               channels(), "test/initial-joint-state");
      auto receipt =
          restored->stochastic_completion(record.integration_owner_id);
      check(receipt &&
                receipt->after_rng_clock_sha256 ==
                    record.after_rng_clock_sha256 &&
                receipt->stage_sequence == 0,
            "zero-event completion survives semantic restart");
      auto bad = record;
      bad.stage_sequence = 1;
      bad.before_rng_clock_sha256 = std::string(64, 'f');
      bad.after_rng_clock_sha256 = std::string(64, 'c');
      reject(
          [&] {
            c.arena->prepare_compound_with_completion(
                1, [&](auto &) { return bad; });
          },
          "wrong prior RNG-clock hash refused evenzeroevents");
      bad = record;
      bad.stage_sequence = 1;
      bad.before_rng_clock_sha256 = std::string(64, 'b');
      bad.after_rng_clock_sha256 = std::string(64, 'c');
      bad.dataset_sha256 = std::string(64, 'e');
      reject(
          [&] {
            c.arena->prepare_compound_with_completion(
                1, [&](auto &) { return bad; });
          },
          "wrong physical-data binding refused evenzeroevents");
      check(c.arena->state_prefix() == prefix0 &&
                m::PicArena::fingerprint(c.pic) == native,
            "failed completion leaves accepted state exact");
      record.stage_sequence = 1;
      record.before_rng_clock_sha256 = std::string(64, 'b');
      record.after_rng_clock_sha256 = std::string(64, 'c');
      auto second = c.arena->prepare_compound_with_completion(
          1, [&](auto &) { return record; });
      check(c.arena->try_publish(second) && c.arena->state_prefix() != prefix0,
            "second no-event exposure advances sequence/prefix atomically");
      const auto prefix1 = c.arena->state_prefix();
      check(!c.arena->try_rollback(zero),
            "older completion cannot roll back a newer accepted exposure");
      auto next = record;
      next.stage_sequence = 2;
      next.before_rng_clock_sha256 = std::string(64, 'c');
      next.after_rng_clock_sha256 = std::string(64, 'd');
      auto invalid = next;
      invalid.model = "foreign_model";
      reject(
          [&] {
            c.arena->prepare_compound_with_completion(1, [&](auto &a) {
              auto child = a.prepare(c.event(bridge));
              check(a.try_publish(child),
                    "reaction privately accepted before bad completion");
              return invalid;
            });
          },
          "bad terminal completion rolls back all private reactions");
      check(c.arena->state_prefix() == prefix1 &&
                m::PicArena::fingerprint(c.pic) == native &&
                c.arena->snapshot()->event_ids.empty(),
            "invalid completion restores full private event/RNG-prefix base");
      auto reaction =
          c.arena->prepare_compound_with_completion(1, [&](auto &a) {
            auto child = a.prepare(c.event(bridge));
            check(a.try_publish(child),
                  "private actual kernel event before valid completion");
            return next;
          });
      check(c.arena->try_publish(reaction) &&
                c.arena->snapshot()->event_ids.size() == 1,
            "reaction plus terminal RNG-clock completion publish together");
      c.arena->save(dir + "/completion-reaction.json");
      m::GapPic fresh2(options());
      auto restored2 = m::PicArena::restore(
          fresh2, dir + "/completion-reaction.json", states(), channels(),
          "test/initial-joint-state");
      check(restored2->stochastic_completion(record.integration_owner_id)
                        ->stage_sequence == 2 &&
                restored2->snapshot()->event_ids.size() == 1,
            "reaction+completion semantic history roundtrip");
      for (int tamper = 0; tamper < 4; ++tamper) {
        std::ifstream in(dir + "/completion-reaction.json");
        J doc;
        in >> doc;
        auto &b = doc["payload"];
        auto &rows = b["metadata"]["journal"];
        auto &completion = rows.back()["payload"]["completion"];
        if (tamper == 0)
          completion["before_rng_clock_sha256"] = std::string(64, 'f');
        if (tamper == 1)
          completion["stage_sequence"] = 99;
        if (tamper == 2)
          completion["kernel_manifest_sha256"] = std::string(64, 'e');
        if (tamper == 3)
          completion["common_interval_end_s"] = -1.;
        std::string chain = rows[0]["payload"]["previous_chain"];
        std::size_t bytes = 0;
        for (auto &row : rows) {
          row["payload"]["previous_chain"] = chain;
          row["sha256"] = sha(row["payload"].dump());
          chain = row["sha256"];
          bytes += row.dump().size();
        }
        b["metadata"]["chain"] = chain;
        b["metadata"]["journal_bytes"] = bytes;
        b["metadata"]["stochastic_completions"][record.integration_owner_id] =
            completion;
        doc["sha256"] = sha(b.dump());
        const auto file = dir + "/corrupt-completion.json";
        {
          std::ofstream out(file);
          out << doc.dump();
        }
        m::GapPic untouched(options());
        const auto before = m::PicArena::fingerprint(untouched);
        reject(
            [&] {
              m::PicArena::restore(untouched, file, states(), channels(),
                                   "test/initial-joint-state");
            },
            "rehash RNG/sequence/model/clock completion corruption refused");
        check(before == m::PicArena::fingerprint(untouched),
              "bad completion restart leaves native unchanged");
      }
      check(c.arena->try_rollback(reaction) &&
                c.arena->state_prefix() == prefix1 &&
                m::PicArena::fingerprint(c.pic) == native,
            "common rollback restores events and stochastic receipt together");
      proofs.push_back({{"kind", "zero_event_completion"},
                        {"native_unchanged", true},
                        {"prefix_changed", prefix0 != initial},
                        {"integration_owner", record.integration_owner_id},
                        {"sequence_checked", true}});
    }

    {
      auto stamp = []() {
        m::ArenaStochasticCompletion c;
        c.integration_owner_id = "waterflow/scoped-ssa";
        c.before_rng_clock_sha256 = std::string(64, 'a');
        c.after_rng_clock_sha256 = std::string(64, 'b');
        c.common_clock_sha256 = std::string(64, '1');
        c.channel_id = w::parent_channel;
        c.dataset_sha256 = datahash;
        c.kernel_manifest_sha256 = manifest;
        c.kernel_binding_file_sha256 = bindingfile;
        c.model = w::parent_model;
        c.common_interval_end_s = 1e-12;
        return c;
      };
      for (int failure = 0; failure < 5; ++failure) {
        Case c(.25, 64);
        const auto initial = c.arena->state_prefix(),
                   native = m::PicArena::fingerprint(c.pic);
        auto completion = stamp();
        reject(
            [&] {
              c.arena->prepare_compound_with_completion(2, [&](auto &a) {
                for (int i = 0; i < 2; ++i) {
                  auto event = c.event(bridge);
                  event.id = "scoped/" + std::to_string(i);
                  event.products[1].id =
                      "scoped-secondary/" + std::to_string(i);
                  event.products[2].id = "scoped-ion/" + std::to_string(i);
                  event.provenance.selection_method =
                      "ssa_weighted_pair_clock_v1";
                  event.phase.clock_prefix = completion.common_clock_sha256;
                  event.provenance.collision_elapsed_s =
                      (i == 0 ? .25 : .75) * 1e-12;
                  if (failure == 0)
                    event.provenance.collision_elapsed_s = 2e-12;
                  if (failure == 1)
                    event.provenance.collision_elapsed_s =
                        (i == 0 ? .75 : .25) * 1e-12;
                  if (failure == 2)
                    event.phase.clock_prefix = std::string(64, '2');
                  if (failure == 3)
                    event.provenance.collision_elapsed_s = .25e-12;
                  if (failure == 4)
                    event.provenance.collision_elapsed_s = 0;
                  auto t = a.prepare(event);
                  check(a.try_publish(t),
                        "private hostile SSA-clock fixture event accepted "
                        "before terminal validation");
                }
                return completion;
              });
            },
            "SSA scoped elapsed range/order/common clock refused live");
        check(c.arena->state_prefix() == initial &&
                  m::PicArena::fingerprint(c.pic) == native &&
                  c.arena->snapshot()->event_ids.empty(),
              "bad SSA completion rolls back every private event and receipt");
      }
      Case c(.25, 64);
      auto completion = stamp();
      auto token = c.arena->prepare_compound_with_completion(2, [&](auto &a) {
        for (int i = 0; i < 2; ++i) {
          auto event = c.event(bridge);
          event.id = "valid-scoped/" + std::to_string(i);
          event.products[1].id = "valid-secondary/" + std::to_string(i);
          event.products[2].id = "valid-ion/" + std::to_string(i);
          event.provenance.selection_method = "ssa_weighted_pair_clock_v1";
          event.phase.clock_prefix = completion.common_clock_sha256;
          event.provenance.collision_elapsed_s = (i == 0 ? .25 : .75) * 1e-12;
          auto t = a.prepare(event);
          check(a.try_publish(t), "private ordered SSA-schema fixture event");
        }
        return completion;
      });
      check(c.arena->try_publish(token),
            "ordered in-range same-context SSA receipts accepted");
      c.arena->save(dir + "/scope-ssa.json");
      m::GapPic fresh(options());
      auto restored =
          m::PicArena::restore(fresh, dir + "/scope-ssa.json", states(),
                               channels(), "test/initial-joint-state");
      check(restored->snapshot()->event_ids.size() == 2,
            "valid scoped SSA history restores");
      for (int bad = 0; bad < 3; ++bad) {
        std::ifstream in(dir + "/scope-ssa.json");
        J doc;
        in >> doc;
        auto &b = doc["payload"];
        auto &rows = b["metadata"]["journal"];
        if (bad == 0)
          rows[0]["payload"]["collision_elapsed_s"] = 2e-12;
        if (bad == 1) {
          rows[0]["payload"]["collision_elapsed_s"] = .75e-12;
          rows[1]["payload"]["collision_elapsed_s"] = .25e-12;
        }
        if (bad == 2)
          rows[0]["payload"]["clock_prefix"] = std::string(64, '2');
        std::string chain = rows[0]["payload"]["previous_chain"];
        std::size_t bytes = 0;
        for (auto &row : rows) {
          row["payload"]["previous_chain"] = chain;
          if (row["payload"]["kind"] == "stochastic_exposure")
            row["payload"]["stage_after_reaction_chain"] = chain;
          row["sha256"] = sha(row["payload"].dump());
          chain = row["sha256"];
          bytes += row.dump().size();
        }
        b["metadata"]["chain"] = chain;
        b["metadata"]["journal_bytes"] = bytes;
        doc["sha256"] = sha(b.dump());
        const auto file = dir + "/corrupt-scoped-ssa.json";
        {
          std::ofstream out(file);
          out << doc.dump();
        }
        m::GapPic untouched(options());
        const auto before = m::PicArena::fingerprint(untouched);
        reject(
            [&] {
              m::PicArena::restore(untouched, file, states(), channels(),
                                   "test/initial-joint-state");
            },
            "rehash out-of-range/decreasing/foreign-clock SSA history refused");
        check(m::PicArena::fingerprint(untouched) == before,
              "invalid scoped history leaves native unchanged");
      }
      {
        Case forced(.25, 64);
        auto record = stamp();
        record.integration_owner_id = "explicit-forced-kernel-comparison";
        auto t =
            forced.arena->prepare_compound_with_completion(1, [&](auto &a) {
              auto e = forced.event(bridge);
              e.provenance.collision_elapsed_s = 2e-12;
              e.phase.clock_prefix = "forced/uncertified-waiting-clock";
              auto child = a.prepare(e);
              check(a.try_publish(child),
                    "explicit forced-kernel comparison privately accepted");
              return record;
            });
        check(forced.arena->try_publish(t),
              "forced-kernel comparison kept distinct from SSA waiting-clock "
              "certification");
        forced.arena->save(dir + "/scope-forced-comparison.json");
        m::GapPic copy(options());
        auto saved = m::PicArena::restore(
            copy, dir + "/scope-forced-comparison.json", states(), channels(),
            "test/initial-joint-state");
        check(saved->snapshot()->event_ids.size() == 1,
              "explicit forced provenance remains distinct after restore");
      }
      {
        Case alias(.25, 64);
        std::vector<m::ArenaPacket> packets;
        for (const auto &[id, p] : alias.arena->snapshot()->packets) {
          (void)id;
          packets.push_back(p);
        }
        alias.arena.reset();
        auto catalog = channels();
        auto other = catalog[0];
        other.id = "conditional/other-registered-channel";
        catalog.push_back(other);
        alias.arena = std::make_unique<m::PicArena>(
            alias.pic, states(), catalog, packets, "test/initial-joint-state");
        auto record = stamp();
        record.channel_id = other.id;
        const auto before = alias.arena->state_prefix();
        reject(
            [&] {
              alias.arena->prepare_compound_with_completion(1, [&](auto &a) {
                auto e = alias.event(bridge);
                e.provenance.selection_method = "ssa_weighted_pair_clock_v1";
                e.provenance.collision_elapsed_s = .5e-12;
                e.phase.clock_prefix = record.common_clock_sha256;
                auto t = a.prepare(e);
                check(a.try_publish(t), "private parent before foreign "
                                        "registered completion channel");
                return record;
              });
            },
            "SSA event must match completion's registered channel and "
            "data/model context");
        check(alias.arena->state_prefix() == before,
              "foreign completion channel restores base");
      }
    }
    {
      electrode_heat::Count old{1, -1e-20};
      electrode_heat_v2::Count heat{1, -1e-20};
      auto a = m::ArenaCount::from_two_word(old),
           b = m::ArenaCount::from_two_word(heat);
      check(a == b,
            "heat-v1 and heat-v2 Count2 inputs import identical exact fields");
      auto out = b.require_two_word<electrode_heat_v2::Count>();
      check(out.hi == heat.hi && out.lo == heat.lo,
            "representable heat-v2 Count2 export is exact");
      auto value = m::ArenaCount(1e20, 0);
      const double tiny = 9.7827104392643903e-59;
      for (int i = 0; i < 100; ++i)
        value = value.debit(tiny);
      check(value.credited(tiny).debit(tiny) == value,
            "repeated actual cold source weight exact count remains stable");
      reject(
          [&] { value.require_two_limb_values(); },
          "heterogeneous/multiple cold debits cannot be collapsed to Count2");
      auto leaves = value.positive_carriers();
      check(leaves.size() <= m::ArenaCount::maximum_carriers &&
                m::ArenaCount::from_positive_carriers(leaves) == value,
            "bounded positive carriers exactly preserve extended gas count");
      auto owned = value.copy_units();
      BN_zero(owned.get());
      check(!value.zero(),
            "mutable integer copy cannot mutate immutable count value");
      auto scalar = m::ArenaCount(std::numeric_limits<double>::denorm_min(), 0);
      check(scalar.positive_carriers() ==
                std::vector<double>{std::numeric_limits<double>::denorm_min()},
            "binary64 minimum positive count encodes without a source floor");
      reject(
          [&] {
            m::ArenaCount(std::numeric_limits<double>::max(), 0)
                .credited(std::numeric_limits<double>::max());
          },
          "count range remains explicitly bounded");
    }
    std::ofstream out(dir + "/arena-tests.json");
    out << J({{"status", "passed"},
              {"check_count", checks.size()},
              {"checks", checks},
              {"proofs", proofs},
              {"actual_native_CPU_fields", true},
              {"GPU", false},
              {"physical_parent_rate_used_for_selection", false},
              {"conditional_kernel_support_eV", {13.5, 1000.}},
              {"full_water_mechanism", false},
              {"MCC_heating_or_EM_current_qualified", false},
              {"ABI",
               {sizeof(m::GapParticle), sizeof(m::GapPicOptions),
                sizeof(m::GapPicDiagnostics), sizeof(m::GapPic)}}})
               .dump(2)
        << '\n';
    check(bool(out), "test result write");
    std::cout << checks.size() << " arena checks passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
