#pragma once
#include "raikou/mfem_pic_full_time.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <openssl/sha.h>
#include <sstream>
namespace m = raikou::mfem_backend;
using J = m::ArenaJson;
int checks = 0;
J results = J::array();
void check(bool x, const char *name) {
  ++checks;
  results.push_back({{"name", name}, {"passed", x}});
  if (!x)
    throw std::runtime_error(name);
}
template <class F> void rejects(F f, const char *name) {
  bool yes = false;
  try {
    f();
  } catch (const std::exception &) {
    yes = true;
  }
  check(yes, name);
}
std::string sha(std::string s) {
  unsigned char b[32];
  SHA256(reinterpret_cast<const unsigned char *>(s.data()), s.size(), b);
  std::ostringstream out;
  for (auto c : b)
    out << std::hex << std::setfill('0') << std::setw(2) << unsigned(c);
  return out.str();
}
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
  return {{"e", 0, 0, -1, m::pic_me, 0, 0}, {"n", 2, 1, 0, 2.9915e-26, 0, 0}};
}
std::vector<m::ArenaChannel> channels() {
  return {{"conditional-storage-control",
           "e",
           "n",
           {"e", "n"},
           sha("data"),
           sha("kernel"),
           sha("binding"),
           "manufactured-storage-control-only"}};
}
m::FullTimeGenesisBinding genesis() {
  m::FullTimeGenesisBinding g;
  g.time_s = 26.31715e-9;
  g.origin_provenance = "manufactured/full-time-storage-control-no-simulation";
  g.field_sha256 = sha("field");
  g.mesh_sha256 = sha("mesh");
  g.material_sha256 = sha("material");
  g.frame_sha256 = sha("frame");
  g.common_clock_sha256 = sha("initial-clock");
  g.heat_sha256 = sha("heat");
  g.source_sha256 = sha("source");
  g.rng_sha256 = sha("chemistry-RNG-and-config");
  return g;
}
auto fresh() {
  return m::PicArena::create_full_time(options(), states(), channels(),
                                       {{"electron",
                                         "e",
                                         "manufactured/source",
                                         {4, 0},
                                         {.004, .004, .0125},
                                         {1e6, 0, 0}},
                                        {"neutral",
                                         "n",
                                         "manufactured/source",
                                         {8, 0},
                                         {.00403, .004, .0125},
                                         {0, 0, 0}}},
                                       genesis());
}
m::FullTimeStepBinding common(m::PicArena &a) {
  const auto b = a.full_time_boundary();
  return {b->accepted_step,
          b->accepted_arena_revision,
          b->time_s,
          b->time_s + 5e-15,
          5e-15,
          b->accepted_arena_sha256,
          b->field_sha256,
          b->mesh_sha256,
          b->material_sha256,
          sha("new-common")};
}
m::ExposureBinding exposure(const m::FullTimeStepBinding &b) {
  return {b,         m::FullTimeSlot::before_transport,
          b.begin_s, .5 * b.dt_s,
          0,         b.accepted_field_sha256};
}
m::ArenaStochasticCompletion completion(const m::ExposureBinding &b) {
  const auto c = channels()[0];
  m::ArenaStochasticCompletion x;
  x.integration_owner_id = "storage-control-owner";
  x.before_rng_clock_sha256 = genesis().rng_sha256;
  x.after_rng_clock_sha256 = sha("after-chemistry");
  x.common_clock_sha256 = b.common.common_clock_sha256;
  x.channel_id = c.id;
  x.dataset_sha256 = c.dataset_sha256;
  x.kernel_manifest_sha256 = c.kernel_manifest_sha256;
  x.kernel_binding_file_sha256 = c.kernel_binding_file_sha256;
  x.model = c.model;
  x.stage_sequence = b.owner_stage_sequence;
  x.frozen_physical_stage_s = b.frozen_physical_stage_s;
  x.common_interval_begin_s = b.common.begin_s;
  x.common_interval_end_s = b.common.end_s;
  return x;
}
