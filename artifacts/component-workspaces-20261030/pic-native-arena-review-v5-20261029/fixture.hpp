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
