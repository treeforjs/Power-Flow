// Native v8 chemistry v3 controls (API-008.txt 3, 6):
//   stage parity   native census -> CarrierIn/CellTargets -> records equal to a
//                  C3 replay (cpu_step, step_digest) from the saved inputs; C3
//                  ledger of the replay equals the native stage ledger
//   partitions     host partitions 1/2/5/8 give identical records, ledgers,
//                  states and accepted identities
//   ledgers        per-event closure, exact charge and nuclei conservation,
//                  population balance (carriers + reservoirs) per step
//   mutations      blob byte, registry mass, purpose, step index, packet-table
//                  identity
//   refusals       coverage policy, status (event cap), multi-carrier projectile;
//                  the outer step is unchanged after a refusal
// usage: test_v8_chemistry OUT_DIR [reduced|full|long=EVENTS]
#include "chem_v8_fixture.hpp"
#include "chem_host.hpp"
#include "raikou/json.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

namespace m = raikou::mfem_backend;
using J = nlohmann::json;
namespace {
int checks = 0, failures = 0;
J check_log = J::array();
void check(bool ok, const std::string &what, const J &detail = J::object()) {
  ++checks;
  if (!ok) {
    ++failures;
    std::cout << "FAIL " << what << " " << detail.dump() << "\n";
  }
  check_log.push_back({{"check", what}, {"ok", ok}, {"detail", detail}});
}
template <class F> std::string refusal_kind(F f) {
  try {
    f();
  } catch (const m::FullTimeChemistryRefusal &e) {
    return std::string("chemistry:") + e.what();
  } catch (const std::exception &e) {
    return std::string("other:") + e.what();
  }
  return "none";
}
pf_v7::ExactSum to_forms(const pfc::ExactSum &e) {
  pfc::ExactSum a = e;
  const bool neg = a.negative();
  if (neg)
    a.negate();
  pf_v7::ExactSum out;
  const auto &l = a.limbs();
  for (int k = 0; k < pfc::ExactSum::kLimbs; ++k)
    for (int half = 0; half < 2; ++half) {
      const std::uint64_t piece = half ? (l[k] >> 32) : (l[k] & 0xFFFFFFFFull);
      if (!piece)
        continue;
      const int ex = 64 * k + 32 * half - 1074, e1 = ex / 2, e2 = ex - e1;
      out.add_product(std::ldexp(double(piece) * (neg ? -1.0 : 1.0), e1), std::ldexp(1.0, e2));
    }
  return out;
}
std::string digits_hex(const pf_v7::ExactSum &x) {
  const auto d = x.digits();
  std::string s;
  for (int k = int(d.size()) - 1; k >= 0; --k) {
    char buf[20];
    std::snprintf(buf, sizeof buf, "%016llx", (unsigned long long)d[k]);
    s += buf;
  }
  return s;
}
const m::FullTimeExactValue *find(const std::vector<std::pair<std::string, m::FullTimeExactValue>> &v,
                                  const std::string &k) {
  for (const auto &[n, x] : v)
    if (n == k)
      return &x;
  return nullptr;
}
// Exact per-species carrier weight of an arena's accepted population.
std::map<std::string, pf_v7::ExactSum> population_weights(const m::PicArena &a) {
  std::map<std::string, pf_v7::ExactSum> out;
  for (const auto &[id, p] : a.snapshot()->packets) {
    (void)id;
    for (double w : m::PicArena::positive_carriers(p.count))
      out[p.state].add(w);
  }
  return out;
}
chemfix::Population dense(std::size_t n, std::uint64_t seed) {
  chemfix::Population P;
  P.carriers = n;
  P.seed = seed;
  P.water_weight = 2e13;
  P.clustered = true;
  P.electron_eV_min = 5;
  P.electron_eV_max = 300;
  return P;
}
} // namespace

int main(int argc, char **argv) {
  try {
    if (argc < 2)
      throw std::runtime_error("usage: test_v8_chemistry OUT_DIR [reduced|full|long=EVENTS]");
    const std::filesystem::path out = argv[1];
    std::filesystem::create_directories(out);
    const std::string mode = argc > 2 ? argv[2] : "full";
    const auto o = chemfix::options(45);
    const double dt = 1e-12;
    J results = J::object();
    // ------------------------------------------------------------ 1. stage parity with a C3 replay
    {
      auto chem = chemfix::chemistry_configuration(o);
      chem.retain_stage_inputs = true;
      auto A = chemfix::make_arena(o, chem, chemfix::population(dense(3000, 11)), 20000);
      const auto blob = pfc::load_blob(chemfix::blob_path(), chemfix::blob_sha());
      J stages = J::array();
      for (int k = 0; k < 3; ++k) {
        const auto r = chemfix::step(A, dt, k);
        for (const auto &st : {r.stage1, r.stage2}) {
          pfc::StepConfig cfg;
          check(st->step_config_bytes.size() == sizeof cfg, "parity: saved StepConfig");
          std::memcpy(&cfg, st->step_config_bytes.data(), sizeof cfg);
          std::vector<pfc::CarrierIn> carriers(st->carrier_bytes.size() / sizeof(pfc::CarrierIn));
          std::memcpy(carriers.data(), st->carrier_bytes.data(), st->carrier_bytes.size());
          std::vector<pfc::CellTargets> cells(st->cell_bytes.size() / sizeof(pfc::CellTargets));
          std::memcpy(cells.data(), st->cell_bytes.data(), st->cell_bytes.size());
          pfc::StepResult ref;
          pfc::cpu_step(blob, cfg, carriers, cells, ref);
          const auto d = pfc::step_digest(ref.outs, ref.events, cfg.max_events, carriers);
          check(d == st->records_sha256, "parity: C3 replay step_digest equals the native records",
                {{"step", k}, {"stage", st->stage}, {"c3", d}, {"native", st->records_sha256}});
          check(cfg.purpose == (st->stage == 1 ? 0x30u : 0x34u) && cfg.step == st->accepted_step &&
                    cfg.dt_s == st->dt_s && cfg.seed == chem.run_seed,
                "parity: stage configuration (purpose, step, dt/2, seed)");
          // C3 ledger of the replay (every event admitted when no debit was refused)
          pfc::Ledger L;
          L.init(blob);
          std::uint64_t events = 0;
          for (std::size_t i = 0; i < carriers.size(); ++i) {
            pfc::ledger_add_carrier(L, ref.outs[i]);
            for (pfc_u32 e = 0; e < ref.outs[i].n_events; ++e) {
              pfc::ledger_add_event(blob, L, ref.events[i * cfg.max_events + e]);
              ++events;
            }
          }
          pfc::ledger_finish(L);
          if (st->refused_debits == 0) {
            const auto *kin = find(st->energy, "ke_in_projectile_eV");
            const auto *q = find(st->energy, "q_reference_eV");
            const auto *def = find(st->energy, "defect_eV");
            check(kin && q && def && kin->digits == digits_hex(to_forms(L.ke_in_projectile)) &&
                      q->digits == digits_hex(to_forms(L.q_reference)) &&
                      def->digits == digits_hex(to_forms(L.defect)) && events == st->admitted_events,
                  "parity: C3 ledger of the replay equals the native stage ledger (exact digits)",
                  {{"events", events}, {"admitted", st->admitted_events}});
          }
          // mutation: purpose and step index change the records
          auto cfg2 = cfg;
          cfg2.purpose = cfg.purpose == 0x30u ? 0x34u : 0x30u;
          pfc::StepResult r2;
          pfc::cpu_step(blob, cfg2, carriers, cells, r2);
          auto cfg3 = cfg;
          cfg3.step += 1;
          pfc::StepResult r3;
          pfc::cpu_step(blob, cfg3, carriers, cells, r3);
          if (st->events > 0) {
            check(pfc::step_digest(r2.outs, r2.events, cfg.max_events, carriers) != d,
                  "mutation: another purpose gives other records");
            check(pfc::step_digest(r3.outs, r3.events, cfg.max_events, carriers) != d,
                  "mutation: another step index gives other records");
          }
          stages.push_back({{"step", k}, {"stage", st->stage}, {"projectiles", st->projectiles},
                            {"events", st->events}, {"admitted", st->admitted_events},
                            {"births", st->births}, {"consumed", st->consumed}, {"released", st->released_packets},
                            {"refused_debits", st->refused_debits}, {"nu_dt_max", st->nu_dt_max},
                            {"records_sha256", st->records_sha256}, {"records_root", st->records_root}});
        }
      }
      results["stage_parity"] = stages;
    }
    // ------------------------------------------------------------ 2. partitions 1/2/5/8
    {
      std::vector<std::string> ids, ledgers;
      for (unsigned P : {1u, 2u, 5u, 8u}) {
        auto chem = chemfix::chemistry_configuration(o, P);
        auto A = chemfix::make_arena(o, chem, chemfix::population(dense(2500, 21)), 20000);
        std::string trace;
        for (int k = 0; k < 3; ++k) {
          const auto r = chemfix::step(A, dt, k);
          trace += r.stage1->records_root + r.stage2->records_root + r.boundary_sha256;
        }
        ids.push_back(chemfix::sha(trace));
        ledgers.push_back(A.arena->full_time_chemistry().cumulative_ledger_sha256);
      }
      check(std::all_of(ids.begin(), ids.end(), [&](auto &x) { return x == ids[0]; }) &&
                std::all_of(ledgers.begin(), ledgers.end(), [&](auto &x) { return x == ledgers[0]; }),
            "partitions 1/2/5/8: identical records, ledgers and accepted identities", {{"traces", ids}});
      results["partitions"] = {{"traces", ids}, {"ledgers", ledgers}};
    }
    // ------------------------------------------------------------ 3. ledgers: conservation and balance
    {
      auto chem = chemfix::chemistry_configuration(o);
      auto A = chemfix::make_arena(o, chem, chemfix::population(dense(4000, 31)), 40000);
      J steps = J::array();
      const auto reg = A.reg;
      for (int k = 0; k < 4; ++k) {
        const auto before = population_weights(*A.arena);
        const auto res_before = A.arena->full_time_chemistry();
        const auto r = chemfix::step(A, dt, k);
        const auto after = population_weights(*A.arena);
        const auto st = A.arena->full_time_chemistry();
        // population balance per species: (carriers after - before) + (reservoir after - before)
        // == created - destroyed of the two stages (exact)
        bool ok = true;
        std::map<std::string, pf_v7::ExactSum> res0, res1;
        for (const auto &[n, v] : res_before.reservoir_counts)
          res0[n].add(v.value); // reservoir stocks are counts of integer-valued weights here
        for (const auto &[n, v] : st.reservoir_counts)
          res1[n].add(v.value);
        J balance = J::object();
        for (const auto &s : reg.species) {
          pf_v7::ExactSum d;
          auto add = [&](const std::map<std::string, pf_v7::ExactSum> &mm, double sign) {
            auto it = mm.find(s.state);
            if (it != mm.end())
              d.add_product(sign, it->second.value());
          };
          add(after, 1);
          add(before, -1);
          add(res1, 1);
          add(res0, -1);
          for (const auto &stg : {r.stage1, r.stage2}) {
            if (const auto *c = find(stg->counts, "created/" + s.state))
              d.add_product(-1.0, c->value);
            if (const auto *x = find(stg->counts, "destroyed/" + s.state))
              d.add_product(1.0, x->value);
          }
          const double v = d.value();
          if (v != 0.0)
            balance[s.state] = v;
          ok = ok && v == 0.0;
        }
        check(ok, "ledgers: population balance per species (carriers + reservoirs = created - destroyed)",
              {{"step", k}, {"residual", balance}});
        ok = true;
        // charge and nuclei: created == destroyed (exact digits) per stage
        for (const auto &stg : {r.stage1, r.stage2}) {
          const auto *qc = find(stg->counts, "charge_created_e"), *qd = find(stg->counts, "charge_destroyed_e");
          ok = ok && qc && qd && qc->digits == qd->digits;
          for (const char *n : {"H", "O", "Fe", "Cr", "Ni", "Mn", "Si", "C", "N"}) {
            const auto *a = find(stg->counts, std::string("nuclei_created/") + n),
                       *b = find(stg->counts, std::string("nuclei_destroyed/") + n);
            ok = ok && a && b && a->digits == b->digits;
          }
        }
        check(ok, "ledgers: exact charge and nuclei (ABI v4, 9 slots) conservation per stage", {{"step", k}});
        steps.push_back({{"step", k}, {"events", r.stage1->admitted_events + r.stage2->admitted_events},
                         {"births", r.stage1->births + r.stage2->births},
                         {"released", r.stage2->released_packets}, {"reservoirs", st.reservoirs},
                         {"census", r.census}});
      }
      const auto st = A.arena->full_time_chemistry();
      const auto *def = find(st.cumulative_energy, "defect_eV");
      check(def != nullptr, "ledgers: cumulative defect present");
      results["ledgers"] = {{"steps", steps}, {"cumulative_defect_eV", def ? def->value : 0.0},
                            {"cumulative_energy", [&] {
                               J e = J::object();
                               for (const auto &[n, v] : st.cumulative_energy)
                                 e[n] = v.value;
                               return e;
                             }()}};
    }
    // ------------------------------------------------------------ 3b. conservative parcel release (w <= 1.5)
    {
      auto chem = chemfix::chemistry_configuration(o);
      chem.release_weight = 1.5;
      auto A = chemfix::make_arena(o, chem, chemfix::population(dense(2000, 33)), 40000);
      std::uint64_t released = 0;
      bool weights = true;
      for (int k = 0; k < 3; ++k) {
        const auto r = chemfix::step(A, dt, k);
        released += r.stage2->released_packets;
        check(r.stage1->released_packets == 0, "release: none in stage 1");
      }
      std::uint64_t carriers = 0;
      for (const auto &[id, p] : A.arena->snapshot()->packets)
        if (id.find("/release/") != std::string::npos) {
          ++carriers;
          const auto w = p.count.positive_carriers();
          weights = weights && w.size() == 1 && w.front() > 0 && w.front() <= 1.5 && p.native_ids.size() == 1;
        }
      check(released > 0 && carriers == released && weights,
            "release: positive parcel packets bounded by w0 become native carriers", {{"released", released}});
      results["release"] = {{"released_packets", released}};
    }
    // ------------------------------------------------------------ 4. mutations
    {
      auto chem = chemfix::chemistry_configuration(o);
      const auto pop = chemfix::population(dense(500, 41));
      // 4a. registry mass changed by one ulp
      {
        auto reg = m::full_time_chemistry_registry(chem, o);
        auto states = m::full_time_chemistry_states(reg);
        for (auto &s : states)
          if (s.id == "pfc-v3/H2O")
            s.mass_kg = std::nextafter(s.mass_kg, 1.0);
        m::FullTimeConfiguration c;
        c.identity = m::FullTimeIdentitySchema::v8;
        c.chemistry = chem;
        m::FullTimeGenesisBinding g;
        auto A = chemfix::make_arena(o, chem, pop, 2000);
        g = A.g;
        const auto k = refusal_kind([&] { (void)m::PicArena::create_full_time(o, states, A.channels, pop, g, A.cfg); });
        check(k.find("chemistry:chemistry_v3/configuration") == 0, "mutation: registry mass (one ulp) refused",
              {{"refusal", k}});
      }
      // 4b. blob byte changed
      {
        const auto copy = out / "mutated-blob.bin";
        std::filesystem::copy_file(chemfix::blob_path(), copy, std::filesystem::copy_options::overwrite_existing);
        {
          std::fstream f(copy, std::ios::in | std::ios::out | std::ios::binary);
          f.seekg(200000);
          char c0;
          f.read(&c0, 1);
          c0 ^= 0x01;
          f.seekp(200000);
          f.write(&c0, 1);
        }
        auto bad = chem;
        bad.blob_path = copy.string();
        const auto k = refusal_kind([&] { (void)m::full_time_chemistry_registry(bad, o); });
        check(k.find("chemistry:chemistry_v3/configuration") == 0, "mutation: changed blob byte refused (sha256)",
              {{"refusal", k}});
        std::filesystem::remove(copy);
      }
      // 4c. genesis RNG seal differs from the chemistry state digest
      {
        auto A = chemfix::make_arena(o, chem, pop, 2000);
        auto g = A.g;
        g.rng_sha256 = chemfix::sha("other");
        const auto k = refusal_kind([&] { (void)m::PicArena::create_full_time(o, A.states, A.channels, pop, g, A.cfg); });
        check(k.find("chemistry:chemistry_v3/configuration") == 0, "mutation: genesis RNG seal refused", {{"refusal", k}});
      }
      // 4d. packet-table identity: equal content -> equal identity; one packet's origin -> other identity;
      // schema v5 vs v8 identities differ (identity change)
      {
        auto A = chemfix::make_arena(o, chem, pop, 2000);
        auto B = chemfix::make_arena(o, chem, pop, 2000);
        auto pop2 = pop;
        pop2[17].origin = "manufactured/other-origin";
        auto C = chemfix::make_arena(o, chem, pop2, 2000);
        const auto a = A.arena->full_time_boundary()->accepted_arena_sha256,
                   b = B.arena->full_time_boundary()->accepted_arena_sha256,
                   c = C.arena->full_time_boundary()->accepted_arena_sha256;
        check(a == b && a != c, "mutation: packet-table identity (equal content equal; one origin changed differs)");
        results["identity"] = {{"v8_accepted_arena_sha256", a}, {"v8_one_origin_changed", c}};
      }
      // 4e. step index of the exposure binding (stale accepted step) refused by native
      {
        auto A = chemfix::make_arena(o, chem, pop, 2000);
        const auto bnd = A.arena->full_time_boundary();
        m::FullTimeStepBinding b{bnd->accepted_step + 1, bnd->accepted_arena_revision, bnd->time_s, bnd->time_s + dt, dt,
                                 bnd->accepted_arena_sha256, bnd->field_sha256, bnd->mesh_sha256, bnd->material_sha256,
                                 chemfix::sha("clock")};
        const auto k = refusal_kind([&] { (void)A.arena->begin_full_time_step(b, 2000, 64); });
        check(k.find("other:") == 0, "mutation: another step index refused at the step binding", {{"refusal", k}});
      }
    }
    // ------------------------------------------------------------ 5. refusals (typed; outer unchanged)
    {
      // 5a. coverage policy: nothing admitted
      {
        auto chem = chemfix::chemistry_configuration(o);
        chem.policy.admitted_regions.clear();
        auto A = chemfix::make_arena(o, chem, chemfix::population(dense(800, 51)), 4000);
        const auto before = A.arena->full_time_boundary()->accepted_arena_sha256;
        const auto k = refusal_kind([&] { (void)chemfix::step(A, dt, 0); });
        check(k.find("chemistry:chemistry_v3/coverage_policy") == 0, "refusal: coverage policy (unadmitted region)",
              {{"refusal", k}});
        check(A.arena->full_time_boundary()->accepted_arena_sha256 == before, "refusal: accepted state unchanged");
      }
      // 5b. status: event cap
      {
        auto chem = chemfix::chemistry_configuration(o);
        chem.max_events = 1;
        auto P = dense(800, 52);
        P.water_weight = 1e15;
        auto A = chemfix::make_arena(o, chem, chemfix::population(P), 4000);
        const auto k = refusal_kind([&] { (void)chemfix::step(A, 4e-12, 0); });
        check(k.find("chemistry:chemistry_v3/status EVENT_CAP") == 0, "refusal: status EVENT_CAP", {{"refusal", k}});
      }
      // 5c. projectile of a multi-carrier packet
      {
        auto chem = chemfix::chemistry_configuration(o);
        auto pop = chemfix::population(dense(300, 53));
        for (auto &p : pop)
          if (p.state == "pfc-v3/e") {
            p.count = m::ArenaCount(3.0, 0x1p-60);
            break;
          }
        auto A = chemfix::make_arena(o, chem, pop, 4000);
        const auto k = refusal_kind([&] { (void)chemfix::step(A, dt, 0); });
        check(k.find("chemistry:chemistry_v3/projectile_packet") == 0, "refusal: multi-carrier projectile packet",
              {{"refusal", k}});
      }
    }
    // ------------------------------------------------------------ 5d. production-v1 preset (chemistry v3.1)
    {
      std::string what = "none";
      m::FullTimeChemistryConfiguration pre;
      try {
        pre = m::full_time_chemistry_preset("production-v1", o, "", 0x5EED5EED20261031ull, 1e6, {45, 45, 5});
      } catch (const std::exception &e) {
        what = e.what();
      }
      check(what == "none" && pre.policy.policy_file_sha256 ==
                                  std::string(m::full_time_chemistry_production_v1_policy_sha256),
            "preset production-v1: v3.1 blob and coverage-policy file validated", {{"refusal", what}});
      if (what == "none") {
        check(!pre.diagnostic_field_target_approximation, "production default requires target energy owner");
        pre.preset.clear(); // explicit component diagnostic, not a production preset
        pre.diagnostic_field_target_approximation = true; // isolate the coverage refusal in this test
        // ions 2 keV relative to the water: an unbounded region refused by the production list; the
        // refusal names the region
        auto P = dense(400, 81);
        P.electron_fraction = 0;
        P.ion_fraction = .9;
        auto pop = chemfix::population(P);
        for (auto &x : pop)
          if (x.state == "pfc-v3/H2O+")
            x.proper_velocity_m_s = {0, 0, 1.5e5}; // about 2.1 keV for H2O+
        auto A = chemfix::make_arena(o, pre, pop, 4000);
        const auto k = refusal_kind([&] { (void)chemfix::step(A, 1e-13, 0); });
        check(k.find("chemistry:chemistry_v3/coverage_policy: unbounded region") == 0 &&
                  k.find("|decade 1e3 eV") != std::string::npos,
              "preset production-v1: keV ion region refused with its name", {{"refusal", k}});
        results["production_v1_refusal"] = k;
      }
    }
    // ------------------------------------------------------------ 6. long run (exact ledger closure over many events)
    if (mode.rfind("long=", 0) == 0) {
      const std::uint64_t want = std::stoull(mode.substr(5));
      auto chem = chemfix::chemistry_configuration(o);
      auto P = dense(100000, 61);
      P.water_weight = 2e11;
      auto A = chemfix::make_arena(o, chem, chemfix::population(P), 400000);
      std::uint64_t events = 0;
      int k = 0;
      J steps = J::array();
      while (events < want && k < 2000) {
        const auto r = chemfix::step(A, dt, k);
        events += r.stage1->admitted_events + r.stage2->admitted_events;
        steps.push_back({r.stage1->admitted_events, r.stage2->admitted_events, r.census, r.stage1->nu_dt_max});
        ++k;
      }
      const auto st = A.arena->full_time_chemistry();
      const auto *qc = find(st.cumulative_counts, "charge_created_e"), *qd = find(st.cumulative_counts, "charge_destroyed_e");
      bool nuclei = true;
      for (const char *n : {"H", "O", "Fe", "Cr", "Ni"}) {
        const auto *a = find(st.cumulative_counts, std::string("nuclei_created/") + n),
                   *b = find(st.cumulative_counts, std::string("nuclei_destroyed/") + n);
        nuclei = nuclei && a && b && a->digits == b->digits;
      }
      check(st.events >= want && qc && qd && qc->digits == qd->digits && nuclei,
            "long run: exact charge/nuclei closure over the cumulative ledger", {{"events", st.events}, {"steps", k}});
      J e = J::object();
      for (const auto &[n, v] : st.cumulative_energy)
        e[n] = {{"value", v.value}};
      results["long_run"] = {{"events", st.events}, {"steps", k}, {"energy", e}, {"per_step", steps},
                             {"coverage", J::parse(st.coverage_json)}};
    }
    results["checks"] = checks;
    results["failures"] = failures;
    results["log"] = check_log;
    results["source_identity"] = m::PicArena::source_identity();
    std::ofstream(out / "v8-chemistry-tests.json") << results.dump(1) << "\n";
    std::cout << (failures ? "FAILED " : "") << checks - failures << " native v8 chemistry controls passed"
              << (failures ? " (" + std::to_string(failures) + " failed)" : std::string()) << "\n";
    return failures ? 1 : 0;
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 2;
  }
}
