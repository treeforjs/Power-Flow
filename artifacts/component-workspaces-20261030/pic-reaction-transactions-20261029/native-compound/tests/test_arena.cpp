#include "raikou/mfem_pic_arena.h"
#include "water_parent_bridge.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <openssl/sha.h>
#include <sstream>
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
      check(m::PicArena::exact_count_equal(electrode_heat::add(count, wgt),
                                           {1, 0}),
            "count exact debit/add restoration");
    }
    reject([&] { m::PicArena::debit({1, -1e-20}, 1e-100); },
           "third residual beyond two words refused");
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
      check(s.arena->snapshot()->packets.at("water_ion").count.hi ==
                e.event_weight,
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
        check(n.count.hi == 1e20 && n.count.lo == -.25 &&
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
      for (int corruption = 0; corruption < 8; ++corruption) {
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
          body["metadata"]["packets"][0]["count"][0] = 2.;
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
    check(arena.snapshot()->sinks[0]["count"] == J::array({2., 0.}) ||
              arena.snapshot()->sinks[1]["count"] == J::array({2., 0.}),
          "exact escaped neutral number receipt");
    bool wall = false;
    for (auto receipt : arena.snapshot()->sinks)
      if (receipt["packet_id"] == "wall") {
        wall = true;
        check(receipt["count"] == J::array({1., -1e-20}) &&
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
