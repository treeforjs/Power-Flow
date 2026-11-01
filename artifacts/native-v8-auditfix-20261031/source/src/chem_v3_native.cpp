// chem_v3_native.cpp -- native v8 chemistry v3 owner core (API-008.txt 3).
#include "chem_v3_native.hpp"
#include "raikou/json.hpp"
#include <fstream>
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <limits>
#include <openssl/sha.h>
#include <set>
#include <sstream>

namespace raikou::mfem_backend {
namespace chem8 {
namespace {
using RK = FullTimeChemistryRefusal::Kind;
[[noreturn]] void refuse(RK k, const std::string &m) { throw FullTimeChemistryRefusal(k, "chemistry_v3/" + m); }
void need(bool ok, RK k, const std::string &m) {
  if (!ok)
    refuse(k, m);
}
struct H {
  pfc::Sha256 s;
  void text(const std::string &t) { s.text(t); }
  void u64(std::uint64_t v) { s.pod(v); }
  void u32(std::uint32_t v) { s.pod(v); }
  void f64(double v) { s.pod(std::bit_cast<std::uint64_t>(v)); }
  void bytes(const void *p, std::size_t n) { s.update(p, n); }
  void sum(const ESum &e) {
    if (e.is_zero()) {
      s.pod(std::uint8_t(0));
      return;
    }
    s.pod(std::uint8_t(1));
    s.pod(e.limbs());
  }
  std::string hex() { return s.hex(); }
};
const double kElementary = 1.602176634e-19, kC = 299792458.0;
} // namespace

std::string sha_hex(const std::string &bytes) { return pfc::sha256_hex(bytes.data(), bytes.size()); }

// ------------------------------------------------------------- registry
std::shared_ptr<const Registry> load_registry(const std::string &path, const std::string &expected) {
  need(!path.empty() && expected.size() == 64, RK::configuration, "configuration: blob path and sha256 required");
  auto reg = std::make_shared<Registry>();
  try {
    reg->blob = std::make_shared<const pfc::Blob>(pfc::load_blob(path, expected));
  } catch (const std::exception &e) {
    refuse(RK::configuration, std::string("configuration: table blob: ") + e.what());
  }
  const auto &b = *reg->blob;
  need(b.n_species == 22, RK::configuration, "configuration: the blob must register 22 species");
  const double c2 = kC * kC;
  for (pfc_u32 i = 0; i < b.n_species; ++i) {
    const auto &sp = b.view.species[i];
    reg->states.push_back(std::string(state_prefix) + b.species_names[i]);
    reg->mass_kg.push_back((sp.mass_eV * kElementary) / c2);
    reg->charge.push_back(sp.charge);
    std::array<std::int32_t, nuclei_slots> n{};
    for (std::uint32_t k = nucleus_H; k <= nucleus_Ni; ++k)
      n[k] = std::int32_t(sp.nuclei[k]); // v3 slots 0 H, 1 O, 2 Fe, 3 Cr, 4 Ni = ABI v4 slots
    need(sp.nuclei[5] == 0, RK::configuration, "configuration: v3 spare nuclei slot is not zero");
    reg->nuclei.push_back(n);
  }
  reg->projectile.assign(b.n_species, 0);
  reg->target.assign(b.n_species, 0);
  for (pfc_u32 g = 0; g < b.n_groups; ++g) {
    reg->projectile[b.view.groups[g].projectile] = 1;
    reg->target[b.view.groups[g].target] = 1;
  }
  std::vector<int> bounded(b.n_regions, 0), unbounded(b.n_regions, 0);
  for (pfc_u32 g = 0; g < b.n_groups; ++g)
    for (pfc_u32 i = 0; i <= b.view.groups[g].grid_n; ++i) {
      const auto r = b.cover.bin_region[b.cover.bin_off[g] + i];
      if (r == none)
        continue;
      (b.cover.bin_row[b.cover.bin_off[g] + i] != none ? bounded : unbounded)[r]++;
    }
  for (pfc_u32 r = 0; r < b.n_regions; ++r)
    reg->region_bounded.push_back(bounded[r] > 0 && unbounded[r] == 0);
  H h;
  h.text("native-v8/chemistry-v3-registry-v1");
  h.text(b.sha256);
  h.text("nuclei-abi-v4");
  h.u64(b.n_species);
  for (pfc_u32 i = 0; i < b.n_species; ++i) {
    h.text(reg->states[i]);
    h.f64(reg->mass_kg[i]);
    h.u64(std::uint64_t(std::int64_t(reg->charge[i])));
    for (auto v : reg->nuclei[i])
      h.u32(std::uint32_t(v));
  }
  reg->sha256 = h.hex();
  return reg;
}
// ------------------------------------------------------------- grid
std::uint32_t Grid::locate(const double *x) const noexcept {
  std::uint64_t idx[3];
  for (int d = 0; d < 3; ++d) {
    if (!(x[d] >= low[d] && x[d] <= high[d]))
      return none;
    const double f = (x[d] - low[d]) / (high[d] - low[d]);
    std::uint64_t i = std::uint64_t(f * double(n[d]));
    if (i >= n[d])
      i = n[d] - 1;
    idx[d] = i;
  }
  const std::uint64_t id = idx[0] + std::uint64_t(n[0]) * (idx[1] + std::uint64_t(n[1]) * idx[2]);
  return vacuum[id] ? std::uint32_t(id) : none;
}
Grid make_grid(const GapPicOptions &o, const std::array<std::uint32_t, 3> &counts) {
  need(o.finite_axial_extent, RK::configuration, "configuration: collision grid needs the finite axial extent");
  Grid g;
  g.low = {0.0, 0.0, o.short_z_m};
  g.high = {o.width_m, o.gap_m, o.feed_z_m};
  g.n = counts;
  std::uint64_t total = 1;
  for (int d = 0; d < 3; ++d) {
    need(counts[d] > 0 && counts[d] <= 1u << 20 && std::isfinite(g.low[d]) && std::isfinite(g.high[d]) &&
             g.high[d] > g.low[d],
         RK::configuration, "configuration: collision grid extent/counts");
    total *= counts[d];
  }
  need(total <= (1ull << 28), RK::configuration, "configuration: collision grid cell cap");
  g.volume = ((g.high[0] - g.low[0]) / counts[0]) * ((g.high[1] - g.low[1]) / counts[1]) *
             ((g.high[2] - g.low[2]) / counts[2]);
  need(std::isfinite(g.volume) && g.volume > 0, RK::configuration, "configuration: collision cell volume");
  g.vacuum.assign(total, 1);
  std::vector<std::array<std::uint64_t, 4>> solid;
  auto aligned = [&](double x, int d) {
    const long double q = (static_cast<long double>(x) - g.low[d]) * counts[d] /
                          (static_cast<long double>(g.high[d]) - g.low[d]);
    const long double r = std::round(q);
    need(r >= 0 && r <= counts[d] &&
             std::abs(q - r) <= 64 * std::numeric_limits<double>::epsilon() * std::max(1.0L, std::abs(q)),
         RK::configuration, "configuration: collision grid cuts a foil");
    return std::uint64_t(r);
  };
  if (o.finite_foils)
    for (const auto &f : {o.lower_foil, o.upper_foil})
      solid.push_back({aligned(f.left_m, 0), aligned(f.right_m, 0), aligned(f.bottom_m, 1), aligned(f.top_m, 1)});
  for (std::uint64_t k = 0; k < counts[2]; ++k)
    for (std::uint64_t j = 0; j < counts[1]; ++j)
      for (std::uint64_t i = 0; i < counts[0]; ++i)
        for (const auto &r : solid)
          if (i >= r[0] && i < r[1] && j >= r[2] && j < r[3])
            g.vacuum[i + counts[0] * (j + std::uint64_t(counts[1]) * k)] = 0;
  H h;
  h.text("native-v8/collision-grid-v1");
  for (int d = 0; d < 3; ++d) {
    h.f64(g.low[d]);
    h.f64(g.high[d]);
    h.u32(g.n[d]);
  }
  h.f64(g.volume);
  h.u64(solid.size());
  for (const auto &r : solid)
    for (auto v : r)
      h.u64(v);
  g.sha256 = h.hex();
  return g;
}
// ------------------------------------------------------------- setup
std::shared_ptr<const Setup> make_setup(const FullTimeChemistryConfiguration &c, const GapPicOptions &o) {
  need(c.enabled, RK::configuration, "configuration: chemistry not enabled");
  need(std::isfinite(c.release_weight) && c.release_weight > 0 && std::isnormal(c.release_weight), RK::configuration,
       "configuration: release weight w0");
  need(std::isfinite(c.trunc_sigmas) && c.trunc_sigmas > 0 && c.trunc_sigmas <= 64, RK::configuration,
       "configuration: truncation radius");
  need(c.max_events >= 1 && c.max_events <= 4096 && c.max_candidates >= 1 && c.max_candidates < (1u << 24),
       RK::configuration, "configuration: candidate/event caps");
  need(c.host_partitions >= 1 && c.host_partitions <= 4096, RK::configuration, "configuration: host partitions");
  need(std::isfinite(c.policy.event_limit) && c.policy.event_limit >= 0 && std::isfinite(c.policy.event_fraction) &&
           c.policy.event_fraction >= 0,
       RK::configuration, "configuration: refusal policy");
  if (c.preset == "production-v1")
    need(!c.diagnostic_field_target_approximation && !c.policy.policy_file.empty(), RK::configuration,
         "configuration: production preset requires strict target-energy and file-backed coverage policy");
  auto s = std::make_shared<Setup>();
  s->config = c;
  s->reg = load_registry(c.blob_path, c.blob_sha256);
  s->grid = make_grid(o, c.grid_counts);
  const auto &b = *s->reg->blob;
  s->admitted.assign(b.n_regions, 0);
  std::set<std::string> names(c.policy.admitted_regions.begin(), c.policy.admitted_regions.end());
  for (const auto &n : names) {
    auto it = std::find(b.region_names.begin(), b.region_names.end(), n);
    need(it != b.region_names.end(), RK::configuration, "configuration: unknown admitted region " + n);
    s->admitted[std::size_t(it - b.region_names.begin())] = 1;
  }
  if (!c.policy.policy_file.empty()) {
    // chemistry v3.1 coverage-policy file (pfc-chem-v3.1-coverage-policy-1)
    need(names.empty(), RK::configuration, "configuration: admitted regions given both inline and by the policy file");
    std::ifstream f(c.policy.policy_file, std::ios::binary);
    need(bool(f), RK::configuration, "configuration: cannot read the coverage-policy file " + c.policy.policy_file);
    const std::string text((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    need(sha_hex(text) == c.policy.policy_file_sha256, RK::configuration,
         "configuration: coverage-policy file sha256 differs");
    nlohmann::json j;
    try {
      j = nlohmann::json::parse(text);
    } catch (const std::exception &e) {
      refuse(RK::configuration, std::string("configuration: coverage-policy file: ") + e.what());
    }
    need(j.value("schema", "") == "pfc-chem-v3.1-coverage-policy-1" && j.value("blob_sha256", "") == b.sha256 &&
             j.value("region_count", 0u) == b.n_regions &&
             j.value("reachability_schema", "") == "pfc-region-reachable-support-v1",
         RK::configuration, "configuration: coverage-policy file schema/blob/region count");
    if (c.preset == "production-v1")
      need(j.value("model_mode", "") == "strict-production" &&
               j.at("admitted_unbounded_regions").is_array() && j.at("admitted_unbounded_regions").empty(),
           RK::configuration, "configuration: production preset cannot admit unknown cross sections");
    const auto &pol = j.at("policy");
    need(std::bit_cast<std::uint64_t>(pol.at("event_limit_physical_events").get<double>()) ==
                 std::bit_cast<std::uint64_t>(c.policy.event_limit) &&
             std::bit_cast<std::uint64_t>(pol.at("fraction").get<double>()) ==
                 std::bit_cast<std::uint64_t>(c.policy.event_fraction),
         RK::configuration, "configuration: event_limit/event_fraction differ from the coverage-policy file");
    auto region = [&](const nlohmann::json &e) {
      const auto i = e.at("index").get<std::uint32_t>();
      need(i < b.n_regions && e.at("name").get<std::string>() == b.region_names[i], RK::configuration,
           "configuration: coverage-policy region index/name differs from the blob");
      return i;
    };
    for (const auto &e : j.at("admitted_unbounded_regions")) {
      const auto i = region(e);
      need(!s->reg->region_bounded[i], RK::configuration, "configuration: admitted region is bounded");
      s->admitted[i] = 1;
      names.insert(b.region_names[i]);
    }
    for (const auto &e : j.at("refusing_unbounded_regions")) {
      const auto i = region(e);
      need(!s->reg->region_bounded[i] && !s->admitted[i], RK::configuration,
           "configuration: refusing region is bounded or admitted");
    }
    std::vector<std::uint8_t> bounded(b.n_regions, 0);
    for (const auto &e : j.at("bounded_regions"))
      bounded[region(e)] = 1;
    for (pfc_u32 r = 0; r < b.n_regions; ++r)
      need((bounded[r] != 0) == (s->reg->region_bounded[r] != 0), RK::configuration,
           "configuration: bounded region set differs from the blob (region " + std::to_string(r) + ")");
  }
  {
    H h;
    h.text("pfc-chem-v3.1.1-coverage-policy-reachable-support-v1");
    h.text("A = admitted physical events (sum of W)");
    h.f64(c.policy.event_limit);
    h.f64(c.policy.event_fraction);
    h.u64(names.size());
    for (const auto &n : names)
      h.text(n);
    if (!c.policy.policy_file.empty()) {
      h.text("policy-file");
      h.text(c.policy.policy_file_sha256);
    }
    s->policy_sha256 = h.hex();
  }
  {
    H h;
    h.text("native-v8.2/chemistry-v3-configuration-v2");
    h.u32(c.diagnostic_field_target_approximation ? 1 : 0);
    h.text(b.sha256);
    h.text(s->reg->sha256);
    h.text(s->grid.sha256);
    h.text(s->policy_sha256);
    h.u64(c.run_seed);
    h.f64(c.release_weight);
    h.f64(c.trunc_sigmas);
    h.u32(c.max_candidates);
    h.u32(c.max_events);
    h.text("pfc-chem-stream-v3;release=pfc-neutral-parcels-v1-no-rng");
    for (std::uint32_t p : {PFC_PURPOSE_MCC, PFC_PURPOSE_MCC_SECOND, PFC_PURPOSE_RELEASE_COUNT, PFC_PURPOSE_RELEASE_PICK})
      h.u32(p);
    h.text("stage-split dt/2 x 2, step = accepted step of the base state");
    if (!c.preset.empty()) {
      h.text("preset");
      h.text(c.preset);
    }
    s->configuration_sha256 = h.hex();
  }
  auto &ch = s->channel;
  ch.id = channel_id;
  ch.primary_state = std::string(state_prefix) + "e";
  ch.target_state = std::string(state_prefix) + "H2O";
  ch.product_states = {std::string(state_prefix) + "e", std::string(state_prefix) + "H2O"};
  ch.dataset_sha256 = b.sha256;
  ch.kernel_manifest_sha256 = s->reg->sha256;
  ch.kernel_binding_file_sha256 = s->configuration_sha256;
  ch.model = model_id;
  return s;
}
// ------------------------------------------------------------- ledgers
void Ledger::merge(const Ledger &o) {
  auto &a = c3;
  const auto &b = o.c3;
  for (auto [x, y] : {std::pair{&a.ke_in_projectile, &b.ke_in_projectile}, {&a.ke_in_target, &b.ke_in_target},
                      {&a.ke_out_survivor, &b.ke_out_survivor}, {&a.ke_out_birth, &b.ke_out_birth},
                      {&a.ke_out_reservoir, &b.ke_out_reservoir}, {&a.ke_out_gas, &b.ke_out_gas},
                      {&a.q_reference, &b.q_reference}, {&a.defect, &b.defect},
                      {&a.charge_created, &b.charge_created}, {&a.charge_destroyed, &b.charge_destroyed}})
    x->add(*y);
  for (std::size_t k = 0; k < a.internal.size(); ++k)
    a.internal[k].add(b.internal[k]);
  for (std::size_t k = 0; k < a.created.size(); ++k) {
    a.created[k].add(b.created[k]);
    a.destroyed[k].add(b.destroyed[k]);
  }
  for (int k = 0; k < 6; ++k) {
    a.nuclei_created[k].add(b.nuclei_created[k]);
    a.nuclei_destroyed[k].add(b.nuclei_destroyed[k]);
  }
  for (std::size_t k = 0; k < a.process.size(); ++k) {
    a.process[k].events += b.process[k].events;
    a.process[k].weight.add(b.process[k].weight);
  }
  for (int k = 0; k < 4; ++k)
    a.events_by_cover[k] += b.events_by_cover[k];
  for (std::size_t k = 0; k < a.region_exposure.size(); ++k) {
    a.region_exposure[k].add(b.region_exposure[k]);
    a.region_bound[k].add(b.region_bound[k]);
    need(b.region_support[k] <= UINT64_MAX - a.region_support[k], RK::configuration,
         "configuration: chemistry region support count overflow");
    a.region_support[k] += b.region_support[k];
  }
  a.events += b.events;
  a.closure_violations += b.closure_violations;
  a.max_rel_defect = std::max(a.max_rel_defect, b.max_rel_defect);
  a.carriers += b.carriers;
  a.candidates += b.candidates;
  for (const auto &[k, v] : b.status_counts)
    a.status_counts[k] += v;
  a.nu_dt_max = std::max(a.nu_dt_max, b.nu_dt_max);
  for (int k = 0; k < 16; ++k) {
    nuclei_created[k].add(o.nuclei_created[k]);
    nuclei_destroyed[k].add(o.nuclei_destroyed[k]);
  }
  for (const auto &[k, v] : o.gas_heating)
    gas_heating[k].add(v);
  admitted_events += o.admitted_events;
  refused_debits += o.refused_debits;
  births += o.births;
  consumed += o.consumed;
  released_packets += o.released_packets;
  debited_packets += o.debited_packets;
  removed_packets += o.removed_packets;
  remainder_carriers += o.remainder_carriers;
  stages += o.stages;
  for (auto [x, y] : {std::pair{&released_weight, &o.released_weight}, {&birth_weight, &o.birth_weight},
                      {&consumed_weight, &o.consumed_weight}, {&debited_weight, &o.debited_weight}})
    x->add(*y);
}
std::string ledger_digest(const Ledger &L) {
  H h;
  h.text("native-v8.2/chemistry-ledger-v2-reachable-support");
  const auto &a = L.c3;
  for (const auto *x : {&a.ke_in_projectile, &a.ke_in_target, &a.ke_out_survivor, &a.ke_out_birth,
                        &a.ke_out_reservoir, &a.ke_out_gas, &a.q_reference, &a.defect, &a.charge_created,
                        &a.charge_destroyed})
    h.sum(*x);
  for (const auto &x : a.internal)
    h.sum(x);
  h.u64(a.created.size());
  for (std::size_t k = 0; k < a.created.size(); ++k) {
    h.sum(a.created[k]);
    h.sum(a.destroyed[k]);
  }
  h.u64(a.process.size());
  for (const auto &p : a.process) {
    h.u64(p.events);
    h.sum(p.weight);
  }
  for (auto v : a.events_by_cover)
    h.u64(v);
  h.u64(a.region_exposure.size());
  for (std::size_t k = 0; k < a.region_exposure.size(); ++k) {
    h.sum(a.region_exposure[k]);
    h.sum(a.region_bound[k]);
    h.u64(a.region_support[k]);
  }
  for (auto v : {a.events, a.closure_violations, a.carriers, a.candidates})
    h.u64(v);
  h.f64(a.max_rel_defect);
  h.f64(a.nu_dt_max);
  h.u64(a.status_counts.size());
  for (const auto &[k, v] : a.status_counts) {
    h.u32(k);
    h.u64(v);
  }
  for (int k = 0; k < 16; ++k) {
    h.sum(L.nuclei_created[k]);
    h.sum(L.nuclei_destroyed[k]);
  }
  h.u64(L.gas_heating.size());
  for (const auto &[k, v] : L.gas_heating) {
    h.u32(k.cell);
    h.u32(k.species);
    h.sum(v);
  }
  for (auto v : {L.admitted_events, L.refused_debits, L.births, L.consumed, L.released_packets, L.debited_packets,
                 L.removed_packets, L.remainder_carriers, L.stages})
    h.u64(v);
  for (const auto *x : {&L.released_weight, &L.birth_weight, &L.consumed_weight, &L.debited_weight})
    h.sum(*x);
  return h.hex();
}
std::string state_digest(const FullTimeChemistryState &st) {
  if (!st.digest.empty())
    return st.digest;
  H h;
  h.text("native-v8.2/chemistry-state-v3-parcels-support");
  h.text(st.configuration_sha256);
  h.u64(st.stages);
  h.u32(st.pending ? 1 : 0);
  h.u64(st.reservoirs.size());
  for (const auto &[k, r] : st.reservoirs) {
    h.u32(k.cell);
    h.u32(k.species);
    h.sum(r.count);
    for (const auto &m : r.momentum)
      h.sum(m);
    h.sum(r.kinetic);
    h.u64(r.parcels.size());
    for (const auto &e : r.parcels) {
      h.f64(e.weight);
      for (double v : e.u)
        h.f64(v);
      for (double v : e.x)
        h.f64(v);
    }
  }
  h.text(ledger_digest(st.cumulative));
  h.text(ledger_digest(st.step));
  h.text(ledger_digest(st.last_step));
  st.digest = h.hex();
  return st.digest;
}
std::shared_ptr<FullTimeChemistryState> initial_state(const Setup &s) {
  auto st = std::make_shared<FullTimeChemistryState>();
  st->configuration_sha256 = s.configuration_sha256;
  st->cumulative.init(*s.reg->blob);
  st->step.init(*s.reg->blob);
  st->last_step.init(*s.reg->blob);
  return st;
}
// ------------------------------------------------------------- inputs
StageInputs build_inputs(const Setup &s, const RowView &v, const StageRequest &q) {
  const auto &reg = *s.reg;
  StageInputs in;
  need(q.stage == 1 || q.stage == 2, RK::configuration, "configuration: stage");
  need(q.step < PFA_STEP_LIMIT, RK::carrier_domain, "carrier_domain: accepted step outside the stream domain");
  need(std::isfinite(q.dt_s) && q.dt_s > 0, RK::configuration, "configuration: stage duration");
  in.cfg.seed = s.config.run_seed;
  in.cfg.step = q.step;
  in.cfg.dt_s = q.dt_s;
  in.cfg.trunc_sigmas = s.config.trunc_sigmas;
  in.cfg.purpose = q.stage == 1 ? PFC_PURPOSE_MCC : PFC_PURPOSE_MCC_SECOND;
  in.cfg.max_candidates = s.config.max_candidates;
  in.cfg.max_events = s.config.max_events;
  in.cfg.pad = 0;
  in.row_cell.assign(v.n, none);
  std::vector<std::uint32_t> occupied;
  for (std::size_t i = 0; i < v.n; ++i) {
    const auto sp = v.species[i];
    if (sp == none)
      continue;
    const auto c = s.grid.locate(v.rows[i].position_m.data());
    need(c != none, RK::carrier_domain, "carrier_domain: carrier " + std::to_string(v.rows[i].id) +
                                            " outside the collision grid or in a solid cell");
    in.row_cell[i] = c;
    if (reg.projectile[sp]) {
      need(v.rows[i].id < PFC_CARRIER_LIMIT, RK::carrier_domain, "carrier_domain: native ID >= 2^48");
      need(v.single[i], RK::projectile_packet,
           "projectile_packet: projectile carrier " + std::to_string(v.rows[i].id) + " of a multi-carrier packet");
      occupied.push_back(c);
    }
  }
  std::sort(occupied.begin(), occupied.end());
  occupied.erase(std::unique(occupied.begin(), occupied.end()), occupied.end());
  in.cell_atlas = occupied;
  auto compact = [&](std::uint32_t atlas) -> std::uint32_t {
    auto it = std::lower_bound(occupied.begin(), occupied.end(), atlas);
    return (it != occupied.end() && *it == atlas) ? std::uint32_t(it - occupied.begin()) : none;
  };
  // targets present in occupied cells: (compact, species, row), canonical order
  struct T {
    std::uint32_t cell, species, row;
  };
  std::vector<T> targets;
  for (std::size_t i = 0; i < v.n; ++i) {
    const auto sp = v.species[i];
    if (sp == none || !reg.target[sp])
      continue;
    const auto c = compact(in.row_cell[i]);
    if (c != none)
      targets.push_back({c, sp, std::uint32_t(i)});
  }
  std::sort(targets.begin(), targets.end(), [](const T &a, const T &b) {
    return a.cell != b.cell ? a.cell < b.cell : (a.species != b.species ? a.species < b.species : a.row < b.row);
  });
  in.cells.assign(occupied.size(), pfc::CellTargets{});
  for (auto &c : in.cells)
    std::memset(&c, 0, sizeof c);
  std::size_t i = 0;
  while (i < targets.size()) {
    std::size_t j = i;
    while (j < targets.size() && targets[j].cell == targets[i].cell && targets[j].species == targets[i].species)
      ++j;
    pf_v7::ExactSum sw, su[3], s2;
    ESum inv;
    for (std::size_t k = i; k < j; ++k) {
      const auto &r = v.rows[targets[k].row];
      sw.add(r.weight);
      inv.add(r.weight);
      for (int a = 0; a < 3; ++a)
        su[a].add_product(r.weight, r.proper_velocity_m_s[a]);
    }
    const double W = sw.value();
    double drift[3];
    for (int a = 0; a < 3; ++a)
      drift[a] = su[a].value() / W;
    for (std::size_t k = i; k < j; ++k) {
      const auto &r = v.rows[targets[k].row];
      const double t0 = r.proper_velocity_m_s[0] - drift[0], t1 = r.proper_velocity_m_s[1] - drift[1],
                   t2 = r.proper_velocity_m_s[2] - drift[2];
      const double q2 = t0 * t0 + t1 * t1 + t2 * t2;
      s2.add_product(r.weight, q2);
    }
    auto &cell = in.cells[targets[i].cell];
    need(cell.n_targets < PFC_MAX_TARGETS, RK::cell_targets,
         "cell_targets: more than 12 target species in collision cell " + std::to_string(occupied[targets[i].cell]));
    const auto slot = cell.n_targets++;
    cell.species[slot] = targets[i].species;
    cell.n[slot] = W / s.grid.volume;
    for (int a = 0; a < 3; ++a)
      cell.drift[slot][a] = drift[a];
    cell.s[slot] = std::sqrt((s2.value() / W) / 3.0);
    in.inventory[pfc::InventoryKey{targets[i].cell, targets[i].species}] = inv;
    i = j;
  }
  for (std::size_t r = 0; r < v.n; ++r) {
    const auto sp = v.species[r];
    if (sp == none || !reg.projectile[sp])
      continue;
    pfc::CarrierIn c;
    std::memset(&c, 0, sizeof c);
    c.id = v.rows[r].id;
    for (int a = 0; a < 3; ++a)
      c.u[a] = v.rows[r].proper_velocity_m_s[a];
    c.weight = v.rows[r].weight;
    c.species = sp;
    c.cell = compact(in.row_cell[r]);
    in.carriers.push_back(c);
    in.carrier_row.push_back(std::uint32_t(r));
  }
  std::vector<std::uint64_t> atlas(in.cell_atlas.begin(), in.cell_atlas.end());
  in.sha256 = inputs_digest(in.cfg, block_root("native-v8/chemistry-carriers", in.carriers.data(), sizeof(pfc::CarrierIn),
                                               in.carriers.size(), 1024),
                            in.carriers.size(),
                            block_root("native-v8/chemistry-cells", in.cells.data(), sizeof(pfc::CellTargets),
                                       in.cells.size(), 64),
                            in.cells.size(), block_root("native-v8/chemistry-atlas", atlas.data(), 8, atlas.size(), 4096));
  return in;
}
std::string block_root_of(const std::string &tag, std::uint64_t n, const std::vector<std::array<unsigned char, 32>> &leaves) {
  pfc::Sha256 h;
  h.text(tag + "/root");
  h.pod(n);
  h.pod(std::uint64_t(leaves.size()));
  for (const auto &l : leaves)
    h.update(l.data(), 32);
  return h.hex();
}
std::string block_root(const std::string &tag, const void *data, std::size_t record_bytes, std::uint64_t n,
                       std::uint64_t per_leaf) {
  std::vector<std::array<unsigned char, 32>> leaves;
  const auto *bytes = static_cast<const unsigned char *>(data);
  const std::string leaf_tag = tag + "/leaf";
  for (std::uint64_t l = 0; l * per_leaf < n; ++l) {
    const std::uint64_t lo = l * per_leaf, hi = std::min<std::uint64_t>(n, lo + per_leaf), tl = leaf_tag.size(),
                        cnt = hi - lo;
    std::string b;
    b.append(reinterpret_cast<const char *>(&tl), 8);
    b += leaf_tag;
    b.append(reinterpret_cast<const char *>(&l), 8);
    b.append(reinterpret_cast<const char *>(&cnt), 8);
    b.append(reinterpret_cast<const char *>(bytes + lo * record_bytes), cnt * record_bytes);
    std::array<unsigned char, 32> d;
    SHA256(reinterpret_cast<const unsigned char *>(b.data()), b.size(), d.data());
    leaves.push_back(d);
  }
  return block_root_of(tag, n, leaves);
}
std::string inputs_digest(const pfc::StepConfig &cfg, const std::string &carriers_root, std::uint64_t n_carriers,
                          const std::string &cells_root, std::uint64_t n_cells, const std::string &atlas_root) {
  H h;
  h.text("native-v8/chemistry-inputs-v2");
  h.bytes(&cfg, sizeof cfg);
  h.u64(n_carriers);
  h.text(carriers_root);
  h.u64(n_cells);
  h.text(cells_root);
  h.text(atlas_root);
  return h.hex();
}
// ------------------------------------------------------------- host stage
StageRecords run_host(const Setup &s, const StageInputs &in, unsigned partitions) {
  const auto &b = *s.reg->blob;
  const std::size_t n = in.carriers.size();
  StageRecords out;
  out.outs.assign(n, pfc::CarrierOut{});
  const unsigned P = std::max(1u, std::min<unsigned>(partitions, unsigned(std::max<std::size_t>(n, 1))));
  std::vector<std::vector<pfc::EventRec>> part(P);
  std::vector<std::vector<std::uint32_t>> counts(P);
  std::vector<pfc::EventRec> scratch(in.cfg.max_events);
  // Partitions are contiguous carrier ranges processed in reverse order; the
  // records depend only on (seed, step, carrier, candidate), so any partition
  // and order give the same bytes.
  for (unsigned p = P; p-- > 0;) {
    const std::size_t lo = n * p / P, hi = n * (p + 1) / P;
    for (std::size_t i = lo; i < hi; ++i) {
      const auto &c = in.carriers[i];
      std::memset(scratch.data(), 0, scratch.size() * sizeof(pfc::EventRec));
      if (c.cell >= in.cells.size()) {
        auto &o = out.outs[i];
        std::memset(&o, 0, sizeof o);
        o.id = c.id;
        o.status = PFC_ST_BAD_INPUT;
        continue;
      }
      std::memset(&out.outs[i], 0, sizeof(pfc::CarrierOut));
      pfc::mcc_carrier(b.view, b.cover, in.cfg, c, in.cells[c.cell], scratch.data(), &out.outs[i]);
      counts[p].push_back(out.outs[i].n_events);
      part[p].insert(part[p].end(), scratch.begin(), scratch.begin() + out.outs[i].n_events);
    }
  }
  out.event_begin.assign(n + 1, 0);
  for (std::size_t i = 0; i < n; ++i)
    out.event_begin[i + 1] = out.event_begin[i] + out.outs[i].n_events;
  out.events.reserve(out.event_begin[n]);
  for (unsigned p = 0; p < P; ++p)
    out.events.insert(out.events.end(), part[p].begin(), part[p].end());
  if (out.events.size() != out.event_begin[n])
    throw std::logic_error("chemistry host stage event count");
  return out;
}
std::string records_stream_sha256(const StageRecords &r) {
  pfc::Sha256 h;
  h.text("pfc-chem-v3-step-records-v1");
  for (std::size_t i = 0; i < r.outs.size(); ++i) {
    h.pod(r.outs[i]);
    for (auto k = r.event_begin[i]; k < r.event_begin[i + 1]; ++k)
      h.pod(r.events[k]);
  }
  return h.hex();
}
std::string records_root(const StageRecords &r) {
  const std::size_t n = r.outs.size(), leaf = 256;
  pfc::Sha256 root;
  root.text("native-v8/chemistry-records-root");
  root.pod(std::uint64_t(n));
  root.pod(std::uint64_t((n + leaf - 1) / leaf));
  for (std::size_t k = 0; k * leaf < n; ++k) {
    const std::size_t lo = k * leaf, hi = std::min(n, lo + leaf);
    std::string bytes;
    auto put = [&](const void *p, std::size_t m) { bytes.append(static_cast<const char *>(p), m); };
    const std::string tag = "native-v8/chemistry-records-leaf";
    const std::uint64_t tl = tag.size(), kk = k, cnt = hi - lo;
    put(&tl, 8);
    put(tag.data(), tag.size());
    put(&kk, 8);
    put(&cnt, 8);
    for (std::size_t i = lo; i < hi; ++i) {
      put(&r.outs[i], sizeof(pfc::CarrierOut));
      for (auto e = r.event_begin[i]; e < r.event_begin[i + 1]; ++e)
        put(&r.events[e], sizeof(pfc::EventRec));
    }
    unsigned char d[32];
    SHA256(reinterpret_cast<const unsigned char *>(bytes.data()), bytes.size(), d);
    root.update(d, 32);
  }
  return root.hex();
}
// ------------------------------------------------------------- resolution
std::string status_name(std::uint32_t s) {
  switch (s) {
  case PFC_ST_OK: return "OK";
  case PFC_ST_CONSUMED: return "CONSUMED";
  case PFC_ST_CAND_CAP: return "CAND_CAP";
  case PFC_ST_EVENT_CAP: return "EVENT_CAP";
  case PFC_ST_KEY_DOMAIN: return "KEY_DOMAIN";
  case PFC_ST_MAJORANT: return "MAJORANT";
  case PFC_ST_COVER_CAP: return "COVER_CAP";
  case PFC_ST_BAD_INPUT: return "BAD_INPUT";
  }
  return "UNKNOWN";
}
namespace {
// An event consumes its target iff the target does not reappear as a GAS
// product (INTEGRATION.txt 4.3: elastic, excitation and detachment partners
// return to the gas; ionization, dissociation, attachment, ion-molecule
// reactions and recombination consume it).
bool consumes_target(const pfc::EventRec &e) {
  for (pfc_u32 p = 0; p < e.n_out; ++p)
    if (e.out_dest[p] == PFC_DEST_GAS)
      return false;
  return true;
}
void add_nuclei(std::array<ESum, 16> &acc, const std::array<std::int32_t, 16> &n, double w) {
  for (int k = 0; k < 16; ++k)
    if (n[k])
      acc[k].add_product(w, double(n[k]));
}
} // namespace

StageActive active_of(const Setup &s, const StageInputs &in, const StageRecords &rec, const RowView &v) {
  StageActive a;
  a.cell_atlas = in.cell_atlas;
  a.totals.exposure.assign(s.reg->blob->n_regions, ESum{});
  a.totals.bound.assign(s.reg->blob->n_regions, ESum{});
  pfc::Ledger tmp;
  tmp.init(*s.reg->blob);
  for (std::size_t i = 0; i < in.carriers.size(); ++i) {
    const auto &o = rec.outs[i];
    pfc::ledger_add_carrier(tmp, o);
    need(o.id == in.carriers[i].id && o.n_events == rec.event_begin[i + 1] - rec.event_begin[i], RK::event,
         "event: carrier record identity");
    if (o.status != PFC_ST_OK && o.status != PFC_ST_CONSUMED && a.totals.bad_carrier == ~0ull) {
      a.totals.bad_carrier = o.id;
      a.totals.bad_status = o.status;
    }
    if (!o.n_events)
      continue;
    ActiveCarrier c;
    c.row = in.carrier_row[i];
    c.cell = in.carriers[i].cell;
    c.species = in.carriers[i].species;
    c.id = o.id;
    c.status = o.status;
    c.begin = a.events.size();
    a.events.insert(a.events.end(), rec.events.begin() + std::ptrdiff_t(rec.event_begin[i]),
                    rec.events.begin() + std::ptrdiff_t(rec.event_begin[i + 1]));
    c.end = a.events.size();
    for (int d = 0; d < 3; ++d)
      c.x[d] = v.rows[c.row].position_m[d];
    c.weight = v.rows[c.row].weight;
    a.carriers.push_back(c);
  }
  a.totals.carriers = tmp.carriers;
  a.totals.candidates = tmp.candidates;
  a.totals.status_counts = tmp.status_counts;
  a.totals.nu_dt_max = tmp.nu_dt_max;
  a.totals.exposure = tmp.region_exposure;
  a.totals.bound = tmp.region_bound;
  a.totals.support = tmp.region_support;
  return a;
}
Resolution resolve(const Setup &s, const StageActive &act, const StageRequest &q, const FullTimeChemistryState &before,
                   const DebitAccess &da) {
  const auto &reg = *s.reg;
  const auto &b = *reg.blob;
  Resolution out;
  out.stage.init(b);
  // 1. statuses and record checks
  if (act.totals.bad_carrier != ~0ull)
    refuse(RK::status, "status " + status_name(act.totals.bad_status) + " (carrier " +
                           std::to_string(act.totals.bad_carrier) + ")");
  for (const auto &c : act.carriers) {
    need(c.end > c.begin && c.end <= act.events.size(), RK::event, "event: active carrier range");
    for (auto k = c.begin; k < c.end; ++k) {
      const auto &e = act.events[k];
      need(e.carrier == c.id && e.process < b.n_proc && e.cell == c.cell && e.target_species < b.n_species &&
               e.n_out >= 1 && e.n_out <= PFC_MAX_OUT,
           RK::event, "event: record fields");
      need(pfc::event_conserves(b, e), RK::event, "event: charge/nuclei conservation");
      const auto proj = b.view.groups[b.view.procs[e.process].group].projectile;
      need(proj == c.species, RK::event, "event: projectile species");
      for (pfc_u32 p = 0; p < e.n_out; ++p) {
        const auto dest = e.out_dest[p];
        need(dest >= PFC_DEST_SURVIVOR && dest <= PFC_DEST_GAS && e.out_species[p] < b.n_species, RK::event,
             "event: product destination");
        if (dest == PFC_DEST_SURVIVOR)
          need(e.out_species[p] == proj, RK::event, "event: survivor species change");
        if (dest == PFC_DEST_RESERVOIR)
          need(reg.charge[e.out_species[p]] == 0, RK::event, "event: charged reservoir product");
      }
    }
  }
  // Event-local closure is not closure of actual tracked targets: target
  // velocities are sampled from cell moments, whereas inventory is debited
  // from physical native packets; returned-gas heat also has no receiving
  // energy/momentum owner. Do not silently admit that approximation.
  if (!act.events.empty() && !s.config.diagnostic_field_target_approximation)
    refuse(RK::target_energy_owner,
           "target_energy_owner: no target energy/momentum owner; only explicit diagnostic field-target mode is available");
  // 2. tracked-target debits (C3 order: cell, carrier ID, candidate; prefix rule;
  // target-consuming events only)
  std::vector<std::size_t> order(act.carriers.size());
  for (std::size_t i = 0; i < order.size(); ++i)
    order[i] = i;
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t c) {
    return act.carriers[a].cell != act.carriers[c].cell ? act.carriers[a].cell < act.carriers[c].cell
                                                        : act.carriers[a].id < act.carriers[c].id;
  });
  std::map<pfc::InventoryKey, ESum> left;
  std::vector<std::uint32_t> admitted(act.carriers.size(), 0);
  auto &D = out.debits;
  for (auto i : order) {
    const auto &c = act.carriers[i];
    std::uint32_t k = 0;
    for (; k < c.end - c.begin; ++k) {
      const auto &e = act.events[c.begin + k];
      if (!reg.target[e.target_species] || !consumes_target(e))
        continue;
      const pfc::InventoryKey key{e.cell, e.target_species};
      auto it = left.find(key);
      if (it == left.end())
        it = left.emplace(key, da.inventory(e.cell, e.target_species)).first;
      ESum rem = it->second, w;
      w.add(e.weight);
      rem.sub(w);
      if (rem.negative()) {
        D.refused++;
        break;
      }
      it->second = rem;
      D.debited[key].add(e.weight);
    }
    admitted[i] = k;
  }
  auto st = std::make_shared<FullTimeChemistryState>(before);
  st->digest.clear();
  Ledger &L = out.stage;
  L.stages = 1;
  L.refused_debits = D.refused;
  {
    auto &a = L.c3;
    a.carriers = act.totals.carriers;
    a.candidates = act.totals.candidates;
    a.status_counts = act.totals.status_counts;
    a.nu_dt_max = act.totals.nu_dt_max;
    a.region_exposure = act.totals.exposure;
    a.region_bound = act.totals.bound;
    need(act.totals.support.size() == b.n_regions, RK::event, "event: region support count size");
    a.region_support = act.totals.support;
  }
  std::set<std::uint32_t> consumed_rows;
  // 3. per carrier outcome (ascending carrier ID), ledgers, births, deposits
  for (std::size_t i = 0; i < act.carriers.size(); ++i) {
    const auto &c = act.carriers[i];
    const auto k = admitted[i];
    if (!k)
      continue;
    const auto atlas = act.cell_atlas[c.cell];
    const auto &last = act.events[c.begin + k - 1];
    int survivor = -1;
    for (pfc_u32 p = 0; p < last.n_out; ++p)
      if (last.out_dest[p] == PFC_DEST_SURVIVOR)
        survivor = int(p);
    if (survivor >= 0)
      out.plan.velocity.push_back(
          {c.row, {last.out_u[survivor][0], last.out_u[survivor][1], last.out_u[survivor][2]}});
    else {
      out.plan.consumed_rows.push_back(c.row);
      consumed_rows.insert(c.row);
      L.consumed++;
      L.consumed_weight.add(c.weight);
    }
    for (pfc_u32 j = 0; j < k; ++j) {
      const auto &e = act.events[c.begin + j];
      pfc::ledger_add_event(b, L.c3, e);
      L.admitted_events++;
      const auto proj = b.view.groups[b.view.procs[e.process].group].projectile;
      add_nuclei(L.nuclei_destroyed, reg.nuclei[proj], e.weight);
      add_nuclei(L.nuclei_destroyed, reg.nuclei[e.target_species], e.weight);
      for (pfc_u32 p = 0; p < e.n_out; ++p) {
        add_nuclei(L.nuclei_created, reg.nuclei[e.out_species[p]], e.weight);
        switch (e.out_dest[p]) {
        case PFC_DEST_BIRTH: {
          Birth bb;
          bb.row = c.row;
          bb.carrier = e.carrier;
          bb.candidate = e.candidate;
          bb.product = p;
          bb.species = e.out_species[p];
          bb.process = e.process;
          bb.weight = e.weight;
          for (int a = 0; a < 3; ++a)
            bb.u[a] = e.out_u[p][a];
          out.plan.births.push_back(bb);
          L.births++;
          L.birth_weight.add(e.weight);
          break;
        }
        case PFC_DEST_RESERVOIR: {
          auto &R = st->reservoirs[Key{atlas, e.out_species[p]}];
          const double m = b.view.species[e.out_species[p]].mass_eV;
          const pfc::V3 u = pfc::v3(e.out_u[p][0], e.out_u[p][1], e.out_u[p][2]);
          const pfc::V3 pp = pfc::p_from_u(m, u);
          R.count.add(e.weight);
          R.momentum[0].add_product(e.weight, pp.x);
          R.momentum[1].add_product(e.weight, pp.y);
          R.momentum[2].add_product(e.weight, pp.z);
          const double phase_ke = pfc::ke_from_u(m, u);
          need(e.ke_out[p] == phase_ke, RK::event, "event: reservoir product energy differs from phase space");
          R.kinetic.add_product(e.weight, phase_ke);
          R.parcels.push_back(pfc::ReservoirEntry{e.weight, {u.x, u.y, u.z}, {c.x[0], c.x[1], c.x[2]}});
          break;
        }
        case PFC_DEST_GAS: {
          auto &g = L.gas_heating[Key{atlas, e.out_species[p]}];
          g.add_product(e.weight, e.ke_out[p]);
          g.add_product(-e.weight, e.ke_in_t);
          break;
        }
        default:
          break;
        }
      }
    }
  }
  pfc::ledger_finish(L.c3);
  // 4. target debit allocation (ascending native ID among the packets the stage
  // did not consume as projectiles)
  for (const auto &[k, d] : D.debited) {
    const auto atlas = act.cell_atlas[k.cell];
    ESum left_d = d;
    for (const auto &[r, c] : da.packets(atlas, k.species)) {
      if (left_d.is_zero())
        break;
      if (consumed_rows.count(r))
        continue;
      ESum diff = left_d;
      diff.sub(c);
      TargetDebit td;
      td.row = r;
      if (!diff.negative()) {
        td.amount = c;
        td.whole = true;
        left_d = diff;
        L.removed_packets++;
      } else {
        td.amount = left_d;
        left_d = ESum();
      }
      L.debited_packets++;
      L.debited_weight.add(td.amount);
      out.plan.debits.push_back(td);
    }
    need(left_d.is_zero(), RK::debit_conflict,
         "debit_conflict: admitted debits exceed the target carriers left in cell " + std::to_string(atlas));
  }
  std::sort(out.plan.debits.begin(), out.plan.debits.end(),
            [](const TargetDebit &a, const TargetDebit &c) { return a.row < c.row; });
  for (std::size_t i = 1; i < out.plan.debits.size(); ++i)
    need(out.plan.debits[i].row != out.plan.debits[i - 1].row, RK::debit_conflict, "debit_conflict: packet debited twice");
  // 5. ledgers, coverage policy
  if (q.stage == 1) {
    need(!before.pending, RK::configuration, "configuration: stage 1 with a pending stage");
    st->step = L;
    st->pending = true;
  } else {
    need(before.pending, RK::configuration, "configuration: stage 2 without stage 1");
    st->step.merge(L);
  }
  const auto &step = st->step;
  {
    ESum physical; // A = applied physical events of the step: exact sum of the admitted events' weights
    for (const auto &pt : step.c3.process)
      physical.add(pt.weight);
    const double A = physical.approx();
    const double limit = std::max(s.config.policy.event_limit, s.config.policy.event_fraction * A);
    for (pfc_u32 r = 0; r < b.n_regions; ++r) {
      const double B = step.c3.region_bound[r].approx();
      if (reg.region_bounded[r]) {
        if (B > limit) {
          out.policy_refused = true;
          out.policy_message = "coverage_policy: bounded region " + std::to_string(r) + " expected unresolved events " +
                               std::to_string(B) + " > " + std::to_string(limit) + " (" + b.region_names[r] + ")";
          break;
        }
      } else if (step.c3.region_support[r] > 0 && !s.admitted[r]) {
        out.policy_refused = true;
        out.policy_message = "coverage_policy: unbounded region " + std::to_string(r) + " reachable and not admitted (" +
                             b.region_names[r] + ")";
        break;
      }
    }
  }
  if (out.policy_refused)
    refuse(RK::coverage_policy, out.policy_message);
  // 6. quantized release (stage 2): C3 release_reservoirs with atlas keys
  if (q.stage == 2) {
    Ledger &S = st->step;
    std::map<pfc::InventoryKey, pfc::Reservoir> stock;
    for (const auto &[key, R] : st->reservoirs)
      stock.emplace(pfc::InventoryKey{key.cell, key.species}, R);
    std::vector<pfc::ReleasedPacket> released;
    try {
      released = pfc::release_reservoirs(b, stock, s.config.release_weight, s.config.run_seed, q.step, q.max_release_packets);
    } catch (const std::exception &e) {
      refuse(RK::event, std::string("event: conservative neutral release: ") + e.what());
    }
    for (auto &[key, R] : stock)
      st->reservoirs[Key{key.cell, key.species}] = std::move(R);
    for (const auto &r : released) {
      Release rl;
      rl.cell = r.cell; rl.species = r.species; rl.index = r.index; rl.weight = r.weight;
      for (int a = 0; a < 3; ++a) { rl.u[a] = r.u[a]; rl.x[a] = r.x[a]; }
      out.plan.releases.push_back(rl);
      S.released_packets++;
      S.released_weight.add(r.weight);
    }
    out.stage.released_packets = S.released_packets;
    out.stage.released_weight = S.released_weight;
    st->last_step = st->step;
    st->cumulative.merge(st->step);
    st->step = Ledger();
    st->step.init(b);
    st->pending = false;
  }
  st->stages = before.stages + 1;
  st->digest.clear();
  out.state = std::move(st);
  return out;
}
// ------------------------------------------------------------- ids, values
std::string packet_id_birth(std::uint64_t step, unsigned stage, const Birth &b) {
  return std::string(packet_namespace) + "s" + std::to_string(step) + "/g" + std::to_string(stage) + "/c" +
         std::to_string(b.carrier) + "/e" + std::to_string(b.candidate) + "/k" + std::to_string(b.product);
}
std::string packet_id_release(std::uint64_t step, const Release &r) {
  return std::string(packet_namespace) + "s" + std::to_string(step) + "/release/c" + std::to_string(r.cell) + "/p" +
         std::to_string(r.species) + "/k" + std::to_string(r.index);
}
bool reserved_packet_id(const std::string &id) noexcept { return id.rfind(packet_namespace, 0) == 0; }
namespace {
bool number(const std::string &s, std::size_t &i, char lead) {
  if (i >= s.size() || s[i] != lead)
    return false;
  ++i;
  const std::size_t j = i;
  while (i < s.size() && s[i] >= '0' && s[i] <= '9')
    ++i;
  return i > j && i - j <= 20 && (s[j] != '0' || i - j == 1);
}
} // namespace
bool canonical_packet_id(const std::string &id) noexcept {
  if (!reserved_packet_id(id))
    return false;
  std::size_t i = std::strlen(packet_namespace);
  if (!number(id, i, 's') || i >= id.size() || id[i] != '/')
    return false;
  ++i;
  if (id.compare(i, 8, "release/") == 0) {
    i += 8;
    return number(id, i, 'c') && i < id.size() && id[i++] == '/' && number(id, i, 'p') && i < id.size() &&
           id[i++] == '/' && number(id, i, 'k') && i == id.size();
  }
  return number(id, i, 'g') && i < id.size() && id[i++] == '/' && number(id, i, 'c') && i < id.size() &&
         id[i++] == '/' && number(id, i, 'e') && i < id.size() && id[i++] == '/' && number(id, i, 'k') &&
         i == id.size();
}
pf_v7::ExactSum to_forms(const ESum &e) {
  // value = sum_k limb_k 2^(64k - 1074) (two's complement); exact as a sum of
  // 32-bit pieces times powers of two (each an exact binary64 product).
  ESum a = e;
  const bool neg = a.negative();
  if (neg)
    a.negate();
  pf_v7::ExactSum out;
  const auto &l = a.limbs();
  for (int k = 0; k < ESum::kLimbs; ++k)
    for (int half = 0; half < 2; ++half) {
      const std::uint64_t piece = half ? (l[k] >> 32) : (l[k] & 0xFFFFFFFFull);
      if (!piece)
        continue;
      const int ex = 64 * k + 32 * half - 1074; // piece * 2^ex
      const double p = double(piece) * (neg ? -1.0 : 1.0);
      // split 2^ex into two representable factors
      const int e1 = ex / 2, e2 = ex - e1;
      out.add_product(std::ldexp(p, e1), std::ldexp(1.0, e2));
    }
  return out;
}
std::string exact_digits(const ESum &e) {
  const auto d = to_forms(e).digits();
  std::ostringstream o;
  o << std::hex;
  for (int k = int(d.size()) - 1; k >= 0; --k) {
    char buf[20];
    std::snprintf(buf, sizeof buf, "%016llx", (unsigned long long)d[k]);
    o << buf;
  }
  return o.str();
}
FullTimeExactValue exact_value(const ESum &e) { return {to_forms(e).value(), exact_digits(e)}; }
std::vector<std::pair<std::string, FullTimeExactValue>> energy_values(const Ledger &L) {
  const auto &a = L.c3;
  std::vector<std::pair<std::string, FullTimeExactValue>> v = {
      {"ke_in_projectile_eV", exact_value(a.ke_in_projectile)}, {"ke_in_target_eV", exact_value(a.ke_in_target)},
      {"ke_out_survivor_eV", exact_value(a.ke_out_survivor)},   {"ke_out_birth_eV", exact_value(a.ke_out_birth)},
      {"ke_out_reservoir_eV", exact_value(a.ke_out_reservoir)}, {"ke_out_gas_eV", exact_value(a.ke_out_gas)},
      {"q_reference_eV", exact_value(a.q_reference)}};
  const char *sinks[] = {"none", "rot", "vib", "elec", "ion_internal", "fragment"};
  for (int k = 0; k < 6; ++k)
    v.push_back({std::string("internal_") + sinks[k] + "_eV", exact_value(a.internal[k])});
  v.push_back({"defect_eV", exact_value(a.defect)});
  ESum heat;
  for (const auto &[k, g] : L.gas_heating)
    heat.add(g);
  v.push_back({"gas_heating_eV", exact_value(heat)});
  return v;
}
std::vector<std::pair<std::string, FullTimeExactValue>> count_values(const Setup &s, const Ledger &L) {
  std::vector<std::pair<std::string, FullTimeExactValue>> v;
  const auto &a = L.c3;
  for (std::size_t k = 0; k < a.created.size(); ++k) {
    if (!a.created[k].is_zero())
      v.push_back({"created/" + s.reg->states[k], exact_value(a.created[k])});
    if (!a.destroyed[k].is_zero())
      v.push_back({"destroyed/" + s.reg->states[k], exact_value(a.destroyed[k])});
  }
  v.push_back({"charge_created_e", exact_value(a.charge_created)});
  v.push_back({"charge_destroyed_e", exact_value(a.charge_destroyed)});
  const char *nuc[] = {"H", "O", "Fe", "Cr", "Ni", "Mn", "Si", "C", "N"};
  for (int k = 0; k < 9; ++k) {
    v.push_back({std::string("nuclei_created/") + nuc[k], exact_value(L.nuclei_created[k])});
    v.push_back({std::string("nuclei_destroyed/") + nuc[k], exact_value(L.nuclei_destroyed[k])});
  }
  v.push_back({"released_weight", exact_value(L.released_weight)});
  v.push_back({"birth_weight", exact_value(L.birth_weight)});
  v.push_back({"consumed_weight", exact_value(L.consumed_weight)});
  v.push_back({"debited_weight", exact_value(L.debited_weight)});
  return v;
}
std::string coverage_json(const Setup &s, const FullTimeChemistryState &st) {
  const auto &b = *s.reg->blob;
  std::ostringstream o;
  o.precision(17);
  auto ledger = [&](const Ledger &L) {
    o << "[";
    bool first = true;
    ESum physical;
    for (const auto &pt : L.c3.process)
      physical.add(pt.weight);
    const double threshold = std::max(s.config.policy.event_limit, s.config.policy.event_fraction * physical.approx());
    for (pfc_u32 r = 0; r < b.n_regions; ++r) {
      if (L.c3.region_support[r] == 0 && L.c3.region_exposure[r].is_zero() && L.c3.region_bound[r].is_zero())
        continue;
      const double B = to_forms(L.c3.region_bound[r]).value();
      o << (first ? "" : ",") << "{\"region\":" << r << ",\"name\":\"" << b.region_names[r] << "\""
        << ",\"threshold\":" << threshold << ",\"bound_over_threshold\":" << (threshold > 0 ? B / threshold : 0.0)
        << ",\"bounded\":" << (s.reg->region_bounded[r] ? "true" : "false")
        << ",\"admitted\":" << (s.admitted[r] ? "true" : "false")
        << ",\"reachable_support\":" << L.c3.region_support[r]
        << ",\"exposure\":" << to_forms(L.c3.region_exposure[r]).value()
        << ",\"bound\":" << to_forms(L.c3.region_bound[r]).value() << ",\"exposure_digits\":\""
        << exact_digits(L.c3.region_exposure[r]) << "\",\"bound_digits\":\"" << exact_digits(L.c3.region_bound[r])
        << "\"}";
      first = false;
    }
    o << "]";
  };
  o << "{\"schema\":\"pfc-chem-v3.1.1-coverage-ledger-reachable-support-1\",\"policy_sha256\":\"" << s.policy_sha256
    << "\",\"last_step\":";
  ledger(st.last_step);
  o << ",\"cumulative\":";
  ledger(st.cumulative);
  o << "}";
  return o.str();
}
void fill_receipt(FullTimeChemistryStageReceipt &r, const Setup &s, const StageInputs &in, const StageActive &act,
                  const Resolution &res, const StageRequest &q) {
  r.accepted_step = q.step;
  r.diagnostic_field_target_approximation = s.config.diagnostic_field_target_approximation;
  r.full_system_energy_closure = false;
  r.stage = q.stage;
  r.purpose = in.cfg.purpose;
  r.dt_s = q.dt_s;
  r.projectiles = act.totals.carriers;
  r.cells = act.cell_atlas.size();
  r.candidates = act.totals.candidates;
  r.nu_dt_max = act.totals.nu_dt_max;
  r.events = act.events.size();
  r.admitted_events = res.stage.admitted_events;
  r.refused_debits = res.stage.refused_debits;
  r.births = res.stage.births;
  r.consumed = res.stage.consumed;
  r.debited_packets = res.stage.debited_packets;
  r.removed_packets = res.stage.removed_packets;
  r.released_packets = res.plan.releases.size();
  r.inputs_sha256 = in.sha256;
  r.ledger_sha256 = ledger_digest(res.stage);
  r.state_after_sha256 = state_digest(*res.state);
  r.energy = energy_values(res.stage);
  r.counts = count_values(s, res.stage);
  for (std::size_t k = 0; k < res.stage.c3.region_exposure.size(); ++k)
    r.coverage_regions_exposed += res.stage.c3.region_support[k] > 0;
  if (s.config.retain_stage_inputs) {
    r.step_config_bytes.assign(reinterpret_cast<const char *>(&in.cfg), sizeof in.cfg);
    r.carrier_bytes.assign(reinterpret_cast<const char *>(in.carriers.data()),
                           in.carriers.size() * sizeof(pfc::CarrierIn));
    r.cell_bytes.assign(reinterpret_cast<const char *>(in.cells.data()), in.cells.size() * sizeof(pfc::CellTargets));
  }
}
} // namespace chem8
} // namespace raikou::mfem_backend
// ------------------------------------------------------------- state codec
namespace raikou::mfem_backend::chem8 {
namespace {
struct Out {
  std::string b;
  void u64(std::uint64_t v) { b.append(reinterpret_cast<const char *>(&v), 8); }
  void f64(double v) { u64(std::bit_cast<std::uint64_t>(v)); }
  void text(const std::string &t) {
    u64(t.size());
    b += t;
  }
  void sum(const ESum &e) {
    for (auto l : e.limbs())
      u64(l);
  }
};
struct In {
  const std::string &b;
  std::size_t at = 0;
  std::uint64_t u64() {
    need(at + 8 <= b.size(), RK::configuration, "configuration: chemistry state codec length");
    std::uint64_t v;
    std::memcpy(&v, b.data() + at, 8);
    at += 8;
    return v;
  }
  double f64() { return std::bit_cast<double>(u64()); }
  std::uint64_t count(std::uint64_t cap) {
    const auto n = u64();
    need(n <= cap, RK::configuration, "configuration: chemistry state codec count");
    return n;
  }
  std::string text() {
    const auto n = count(1 << 20);
    need(at + n <= b.size(), RK::configuration, "configuration: chemistry state codec text");
    std::string t = b.substr(at, n);
    at += n;
    return t;
  }
  ESum sum() {
    std::array<std::uint64_t, ESum::kLimbs> l;
    for (auto &x : l)
      x = u64();
    const bool neg = l[ESum::kLimbs - 1] >> 63;
    if (neg) { // magnitude = ~l + 1
      std::uint64_t carry = 1;
      for (auto &x : l) {
        x = ~x;
        const std::uint64_t s = x + carry;
        carry = (carry && s == 0) ? 1 : 0;
        x = s;
      }
    }
    ESum e;
    for (int k = 0; k < ESum::kLimbs; ++k)
      for (int half = 0; half < 2; ++half) {
        const std::uint64_t piece = half ? (l[k] >> 32) : (l[k] & 0xFFFFFFFFull);
        if (!piece)
          continue;
        const int ex = 64 * k + 32 * half - 1074;
        need(ex <= 960, RK::configuration, "configuration: chemistry state value range");
        e.add(std::ldexp(double(piece), ex));
      }
    if (neg)
      e.negate();
    return e;
  }
};
void put_ledger(Out &o, const Ledger &L) {
  const auto &a = L.c3;
  for (const auto *x : {&a.ke_in_projectile, &a.ke_in_target, &a.ke_out_survivor, &a.ke_out_birth,
                        &a.ke_out_reservoir, &a.ke_out_gas, &a.q_reference, &a.defect, &a.charge_created,
                        &a.charge_destroyed})
    o.sum(*x);
  for (const auto &x : a.internal)
    o.sum(x);
  o.u64(a.created.size());
  for (std::size_t k = 0; k < a.created.size(); ++k) {
    o.sum(a.created[k]);
    o.sum(a.destroyed[k]);
  }
  for (int k = 0; k < 6; ++k) {
    o.sum(a.nuclei_created[k]);
    o.sum(a.nuclei_destroyed[k]);
  }
  o.u64(a.process.size());
  for (const auto &p : a.process) {
    o.u64(p.events);
    o.sum(p.weight);
  }
  for (auto v : a.events_by_cover)
    o.u64(v);
  o.u64(a.region_exposure.size());
  for (std::size_t k = 0; k < a.region_exposure.size(); ++k) {
    o.sum(a.region_exposure[k]);
    o.sum(a.region_bound[k]);
    o.u64(a.region_support[k]);
  }
  for (auto v : {a.events, a.closure_violations, a.carriers, a.candidates})
    o.u64(v);
  o.f64(a.max_rel_defect);
  o.f64(a.nu_dt_max);
  o.u64(a.status_counts.size());
  for (const auto &[k, v] : a.status_counts) {
    o.u64(k);
    o.u64(v);
  }
  for (int k = 0; k < 16; ++k) {
    o.sum(L.nuclei_created[k]);
    o.sum(L.nuclei_destroyed[k]);
  }
  o.u64(L.gas_heating.size());
  for (const auto &[k, v] : L.gas_heating) {
    o.u64(k.cell);
    o.u64(k.species);
    o.sum(v);
  }
  for (auto v : {L.admitted_events, L.refused_debits, L.births, L.consumed, L.released_packets, L.debited_packets,
                 L.removed_packets, L.remainder_carriers, L.stages})
    o.u64(v);
  for (const auto *x : {&L.released_weight, &L.birth_weight, &L.consumed_weight, &L.debited_weight})
    o.sum(*x);
}
void get_ledger(In &in, Ledger &L, const pfc::Blob &b) {
  L = Ledger();
  L.init(b);
  auto &a = L.c3;
  for (auto *x : {&a.ke_in_projectile, &a.ke_in_target, &a.ke_out_survivor, &a.ke_out_birth, &a.ke_out_reservoir,
                  &a.ke_out_gas, &a.q_reference, &a.defect, &a.charge_created, &a.charge_destroyed})
    *x = in.sum();
  for (auto &x : a.internal)
    x = in.sum();
  need(in.u64() == a.created.size(), RK::configuration, "configuration: chemistry state species count");
  for (std::size_t k = 0; k < a.created.size(); ++k) {
    a.created[k] = in.sum();
    a.destroyed[k] = in.sum();
  }
  for (int k = 0; k < 6; ++k) {
    a.nuclei_created[k] = in.sum();
    a.nuclei_destroyed[k] = in.sum();
  }
  need(in.u64() == a.process.size(), RK::configuration, "configuration: chemistry state process count");
  for (auto &p : a.process) {
    p.events = in.u64();
    p.weight = in.sum();
  }
  for (auto &v : a.events_by_cover)
    v = in.u64();
  need(in.u64() == a.region_exposure.size(), RK::configuration, "configuration: chemistry state region count");
  for (std::size_t k = 0; k < a.region_exposure.size(); ++k) {
    a.region_exposure[k] = in.sum();
    a.region_bound[k] = in.sum();
    a.region_support[k] = in.u64();
  }
  for (auto *v : {&a.events, &a.closure_violations, &a.carriers, &a.candidates})
    *v = in.u64();
  a.max_rel_defect = in.f64();
  a.nu_dt_max = in.f64();
  const auto ns = in.count(64);
  for (std::uint64_t k = 0; k < ns; ++k) {
    const auto key = std::uint32_t(in.u64());
    a.status_counts[key] = in.u64();
  }
  for (int k = 0; k < 16; ++k) {
    L.nuclei_created[k] = in.sum();
    L.nuclei_destroyed[k] = in.sum();
  }
  const auto ng = in.count(1ull << 32);
  for (std::uint64_t k = 0; k < ng; ++k) {
    Key key{std::uint32_t(in.u64()), std::uint32_t(in.u64())};
    L.gas_heating[key] = in.sum();
  }
  for (auto *v : {&L.admitted_events, &L.refused_debits, &L.births, &L.consumed, &L.released_packets,
                  &L.debited_packets, &L.removed_packets, &L.remainder_carriers, &L.stages})
    *v = in.u64();
  for (auto *x : {&L.released_weight, &L.birth_weight, &L.consumed_weight, &L.debited_weight})
    *x = in.sum();
}
} // namespace
std::string encode_state(const FullTimeChemistryState &st) {
  Out o;
  o.text("native-v8.2/chemistry-state-codec-v3-parcels-support");
  o.text(st.configuration_sha256);
  o.u64(st.stages);
  o.u64(st.pending ? 1 : 0);
  o.u64(st.reservoirs.size());
  for (const auto &[k, r] : st.reservoirs) {
    o.u64(k.cell);
    o.u64(k.species);
    o.sum(r.count);
    for (const auto &m : r.momentum)
      o.sum(m);
    o.sum(r.kinetic);
    o.u64(r.parcels.size());
    for (const auto &e : r.parcels) {
      o.f64(e.weight);
      for (double v : e.u)
        o.f64(v);
      for (double v : e.x)
        o.f64(v);
    }
  }
  put_ledger(o, st.cumulative);
  put_ledger(o, st.step);
  put_ledger(o, st.last_step);
  return o.b;
}
std::shared_ptr<FullTimeChemistryState> decode_state(const Setup &s, const std::string &bytes) {
  In in{bytes};
  need(in.text() == "native-v8.2/chemistry-state-codec-v3-parcels-support", RK::configuration, "configuration: chemistry state codec tag");
  auto st = std::make_shared<FullTimeChemistryState>();
  st->configuration_sha256 = in.text();
  need(st->configuration_sha256 == s.configuration_sha256, RK::configuration,
       "configuration: checkpoint chemistry configuration (blob, registry, policy, w0) differs");
  st->stages = in.u64();
  const auto pend = in.u64();
  need(pend <= 1, RK::configuration, "configuration: chemistry state flag");
  st->pending = pend == 1;
  const auto nr = in.count(1ull << 32);
  for (std::uint64_t i = 0; i < nr; ++i) {
    Key k{std::uint32_t(in.u64()), std::uint32_t(in.u64())};
    need(k.species < s.reg->blob->n_species && k.cell < s.grid.cells(), RK::configuration,
         "configuration: chemistry reservoir key");
    auto &r = st->reservoirs[k];
    r.count = in.sum();
    for (auto &m : r.momentum)
      m = in.sum();
    r.kinetic = in.sum();
    const auto ne = in.count(1ull << 32);
    for (std::uint64_t e = 0; e < ne; ++e) {
      pfc::ReservoirEntry en;
      en.weight = in.f64();
      for (double &v : en.u)
        v = in.f64();
      for (double &v : en.x)
        v = in.f64();
      r.parcels.push_back(en);
    }
    pfc::validate_reservoir(*s.reg->blob, pfc::InventoryKey{k.cell, k.species}, r);
  }
  get_ledger(in, st->cumulative, *s.reg->blob);
  get_ledger(in, st->step, *s.reg->blob);
  get_ledger(in, st->last_step, *s.reg->blob);
  need(in.at == bytes.size(), RK::configuration, "configuration: chemistry state codec trailing bytes");
  return st;
}
} // namespace raikou::mfem_backend::chem8
namespace raikou::mfem_backend::chem8 {
std::string melt_digest(const FullTimeMeltState &m) {
  if (!m.digest.empty())
    return m.digest;
  H h;
  h.text("native-v8/ejected-melt-v1");
  h.u64(m.credits);
  for (const auto &a : m.atoms)
    h.sum(a);
  h.sum(m.energy_J);
  h.sum(m.kinetic_J);
  h.text(m.last_record_sha256);
  m.digest = h.hex();
  return m.digest;
}
std::string encode_melt(const FullTimeMeltState &m) {
  Out o;
  o.text("native-v8/ejected-melt-codec-v1");
  o.u64(m.credits);
  for (const auto &a : m.atoms)
    o.sum(a);
  o.sum(m.energy_J);
  o.sum(m.kinetic_J);
  o.text(m.last_record_sha256);
  return o.b;
}
std::shared_ptr<FullTimeMeltState> decode_melt(const std::string &bytes) {
  In in{bytes};
  need(in.text() == "native-v8/ejected-melt-codec-v1", RK::configuration, "configuration: melt codec tag");
  auto m = std::make_shared<FullTimeMeltState>();
  m->credits = in.u64();
  for (auto &a : m->atoms)
    a = in.sum();
  m->energy_J = in.sum();
  m->kinetic_J = in.sum();
  m->last_record_sha256 = in.text();
  need(in.at == bytes.size(), RK::configuration, "configuration: melt codec trailing bytes");
  return m;
}
} // namespace raikou::mfem_backend::chem8

namespace raikou::mfem_backend::chem8 {
// v8.1 flows state: digest (tag native-v8.1/flows-v1) and canonical codec
// (native-v8.1/flows-codec-v1); maps in key order.
namespace {
template <class W> void put_flows_v81(W &o, const FullTimeFlowsV81 &f) {
  o.u64(f.internal.size());
  for (const auto &[id, x] : f.internal) {
    o.text(id);
    o.f64(x.weight);
    o.f64(x.energy_J);
  }
  o.sum(f.internal_born_J);
  o.sum(f.internal_sunk_J);
  o.sum(f.internal_released_J);
  o.u64(f.wall_births);
  o.u64(f.returned_births);
  o.u64(f.bridges);
  for (const auto *m : {&f.wall_weight, &f.returned_weight}) {
    o.u64(m->size());
    for (const auto &[st, e] : *m) {
      o.text(st);
      o.sum(e);
    }
  }
  o.sum(f.extra_work_J);
  for (auto v : {f.merge_transactions, f.merge_groups, f.merge_inputs, f.merge_outputs, f.merge_relocations,
                 f.merge_formation_refused, f.merge_centroid_fallback})
    o.u64(v);
  for (auto v : f.merge_status)
    o.u64(v);
  o.u64(f.merge_residual.size());
  for (const auto &[st, r] : f.merge_residual) {
    o.text(st);
    for (const auto &e : r)
      o.sum(e);
  }
  o.u64(f.has_last_merge ? 1 : 0);
  o.u64(f.last_merge_step);
  o.u64(f.last_merge_base);
  o.u64(f.last_merge_outputs);
  o.text(f.last_merge_sha256);
  o.text(f.merge_policy_sha256);
  o.u64(f.last_sink_step);
  o.u64(f.last_sink_internal.size());
  for (const auto &[id, e] : f.last_sink_internal) {
    o.u64(id);
    o.f64(e);
  }
}
} // namespace
std::string flows_v81_digest(const FullTimeFlowsV81 &f) {
  if (!f.digest.empty())
    return f.digest;
  H h;
  h.text("native-v8.1/flows-v1");
  put_flows_v81(h, f);
  f.digest = h.hex();
  return f.digest;
}
std::string encode_flows_v81(const FullTimeFlowsV81 &f) {
  Out o;
  o.text("native-v8.1/flows-codec-v1");
  put_flows_v81(o, f);
  return o.b;
}
std::shared_ptr<FullTimeFlowsV81> decode_flows_v81(const std::string &bytes) {
  In in{bytes};
  need(in.text() == "native-v8.1/flows-codec-v1", RK::configuration, "configuration: v8.1 flows codec tag");
  auto f = std::make_shared<FullTimeFlowsV81>();
  const auto ni = in.count(1u << 27);
  for (std::uint64_t k = 0; k < ni; ++k) {
    auto id = in.text();
    FullTimeFlowsV81::Internal x;
    x.weight = in.f64();
    x.energy_J = in.f64();
    need(f->internal.emplace(std::move(id), x).second, RK::configuration, "configuration: v8.1 flows internal key");
  }
  f->internal_born_J = in.sum();
  f->internal_sunk_J = in.sum();
  f->internal_released_J = in.sum();
  f->wall_births = in.u64();
  f->returned_births = in.u64();
  f->bridges = in.u64();
  for (auto *m : {&f->wall_weight, &f->returned_weight}) {
    const auto n = in.count(1u << 16);
    for (std::uint64_t k = 0; k < n; ++k) {
      auto st = in.text();
      need(m->emplace(std::move(st), in.sum()).second, RK::configuration, "configuration: v8.1 flows state key");
    }
  }
  f->extra_work_J = in.sum();
  for (auto *v : {&f->merge_transactions, &f->merge_groups, &f->merge_inputs, &f->merge_outputs,
                  &f->merge_relocations, &f->merge_formation_refused, &f->merge_centroid_fallback})
    *v = in.u64();
  for (auto &v : f->merge_status)
    v = in.u64();
  const auto nr = in.count(1u << 16);
  for (std::uint64_t k = 0; k < nr; ++k) {
    auto st = in.text();
    std::array<ESum, 4> r;
    for (auto &e : r)
      e = in.sum();
    need(f->merge_residual.emplace(std::move(st), r).second, RK::configuration, "configuration: v8.1 merge key");
  }
  const auto has = in.u64();
  need(has <= 1, RK::configuration, "configuration: v8.1 flows flag");
  f->has_last_merge = has == 1;
  f->last_merge_step = in.u64();
  f->last_merge_base = in.u64();
  f->last_merge_outputs = in.u64();
  f->last_merge_sha256 = in.text();
  f->merge_policy_sha256 = in.text();
  f->last_sink_step = in.u64();
  const auto ns = in.count(1u << 27);
  for (std::uint64_t k = 0; k < ns; ++k) {
    const auto id = in.u64();
    f->last_sink_internal.emplace_back(id, in.f64());
  }
  need(in.at == bytes.size(), RK::configuration, "configuration: v8.1 flows codec trailing bytes");
  return f;
}
} // namespace raikou::mfem_backend::chem8
