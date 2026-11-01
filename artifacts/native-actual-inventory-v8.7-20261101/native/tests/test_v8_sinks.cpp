// Native v8 published first-hit sinks (source owner v7 INTEGRATION B1/D8; API-008.txt 5):
// an electron crossing the lower foil's top face in one step is a collected sink; after
// publication PicArena::full_time_published_flows() exposes the record with packet ID, state,
// weight, charge, mass and start position; the checkpoint carries it and a nonempty restore
// gives the same exposure (checked against the journal sink row); a tampered exposure refuses.
// usage: test_v8_sinks OUT_DIR
#include "v6_fixture.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <openssl/sha.h>

namespace m = raikou::mfem_backend;
using JJ = nlohmann::json;
namespace {
int fails = 0, n_checks = 0;
void need(bool ok, const std::string &what) {
  ++n_checks;
  if (!ok) {
    ++fails;
    std::cout << "FAIL " << what << "\n";
  }
}
std::string sha_hex(const std::string &s) {
  unsigned char d[32];
  SHA256(reinterpret_cast<const unsigned char *>(s.data()), s.size(), d);
  char hx[65];
  for (int i = 0; i < 32; ++i)
    std::snprintf(hx + 2 * i, 3, "%02x", d[i]);
  return std::string(hx, 64);
}
} // namespace
int main(int argc, char **argv) {
  try {
    const std::filesystem::path out = argc > 1 ? argv[1] : ".";
    std::filesystem::create_directories(out);
    for (const unsigned anchor_interval : {1u, 2u}) {
    m::FullTimeConfiguration c;
    c.identity = m::FullTimeIdentitySchema::v8;
    c.limits.maximum_carriers = 64;
    c.limits.memory_budget_bytes = 64 * m::full_time_bytes_per_carrier + m::full_time_fixed_bytes;
    c.history.anchor_interval_steps = anchor_interval;
    // electron 0.1 um above the lower foil's top face (y = 2 mm), moving down (u = 3e7 m/s)
    std::vector<m::FullTimeAuthoredPacket> pop = {
        {"electron", "e", "manufactured/source", {4, 0}, {.004, .0020001, 0}, {0, -3e7, 0}},
        {"neutral", "n", "manufactured/source", {8, 0}, {.00403, .004, 0}, {0, 0, 0}}};
    auto arena = m::PicArena::create_full_time(v6::options(), v6::states(), v6::channels(), pop, v6::genesis(), c);
    v6::StepPlan plan;
    plan.dt_s = 1e-14;
    v6::StepState st;
    v6::Hooks hooks;
    hooks.motion = v6::reference_hook(64);
    (void)v6::full_step(*arena, plan, st, hooks);
    const auto f = arena->full_time_published_flows();
    need(f && f->sinks && f->sinks->exposed.size() == 1, "published flows expose the step's sink");
    if (f && f->sinks && !f->sinks->exposed.empty()) {
      const auto &e = f->sinks->exposed.front();
      need(e.packet_id == "electron" && e.state == "e" && e.weight == 4 && e.charge_C == -m::pic_e &&
               e.mass_kg == m::pic_me && e.record.surface == m::FullTimeSinkSurface::lower_foil &&
               e.record.kind == m::FullTimeSinkKind::collected_at_conductor,
           "exposure fields (packet, state, weight, charge, mass, surface, kind)");
    }
    // a second step without sinks replaces the published flows
    const auto saved_boundary = *arena->full_time_boundary();
    arena->save((out / "sinks.json").string());
    auto restored = m::PicArena::restore_full_time((out / "sinks.json").string(), v6::options(), v6::states(),
                                                   v6::channels(), saved_boundary);
    const auto g = restored->full_time_published_flows();
    need(g && g->sinks && f && f->sinks && g->sinks->records_root_sha256 == f->sinks->records_root_sha256 &&
             g->sinks->exposed.size() == f->sinks->exposed.size() &&
             g->sinks->exposed.front().packet_id == f->sinks->exposed.front().packet_id &&
             g->sinks->exposed.front().weight == f->sinks->exposed.front().weight &&
             g->sinks->ledger_sha256 == f->sinks->ledger_sha256 &&
             g->sinks->row_sha256 == f->sinks->row_sha256 &&
             g->sinks->tallies.size() == f->sinks->tallies.size() &&
             g->sinks->surfaces.size() == f->sinks->surfaces.size(),
         "restore: published sinks equal (records root, exposure)");
    need(restored->full_time_boundary()->accepted_arena_sha256 == saved_boundary.accepted_arena_sha256,
         "restore: accepted identity");
    // Both unanchored and compacted checkpoints bind complete exposures.
    for (const std::string mutation : {"weight", "start", "missing", "internal"}) {
      std::ifstream in(out / "sinks.json");
      JJ doc;
      in >> doc;
      auto p = doc.at("payload");
      if (mutation == "missing") p["metadata"].erase("step_sinks");
      else {
        const int index = mutation == "weight" ? 13 : mutation == "start" ? 16 : 19;
        p["metadata"]["step_sinks"]["exposed"][0]["values"][index] = "4014000000000000";
      }
      const auto bad = JJ{{"payload", p}, {"sha256", sha_hex(p.dump())}}.dump();
      std::ofstream(out / "tampered.json") << bad;
      std::string what = "none";
      try {
        (void)m::PicArena::restore_full_time((out / "tampered.json").string(), v6::options(), v6::states(),
                                             v6::channels(), saved_boundary);
      } catch (const std::exception &e) {
        what = e.what();
      }
      need(what != "none", "restore: " + mutation + " sink exposure refused at anchor interval " + std::to_string(anchor_interval) + " (" + what + ")");
    }
    (void)v6::full_step(*arena, plan, st, hooks);
    const auto h = arena->full_time_published_flows();
    need(!h || !h->sinks, "published flows replaced by the next step (no sinks)");
    }
    std::cout << (fails ? "FAILED " : "") << n_checks - fails << " native v8 published-sink controls passed\n";
    return fails ? 1 : 0;
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 2;
  }
}
