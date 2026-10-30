#include "../source/src/full_time_contract_validation.hpp"
#include "full_time_fixture.hpp"
#include <cstdlib>
#include <limits>
#include <new>
#include <type_traits>
static long fail_allocation = -1;
void *operator new(std::size_t n) {
  if (fail_allocation == 0) {
    fail_allocation = -1;
    throw std::bad_alloc();
  }
  if (fail_allocation > 0)
    --fail_allocation;
  if (void *p = std::malloc(n ? n : 1))
    return p;
  throw std::bad_alloc();
}
void *operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void *p) noexcept { std::free(p); }
void operator delete[](void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }
void operator delete[](void *p, std::size_t) noexcept { std::free(p); }
int forbidden_solve_calls = 0, forbidden_advance_calls = 0;
extern "C" void
__wrap__ZN6raikou12mfem_backend6GapPic12solve_fieldsEv(m::GapPic *) {
  ++forbidden_solve_calls;
  throw std::runtime_error("forbidden full-time native solver trap");
}
extern "C" void __wrap__ZN6raikou12mfem_backend6GapPic7advanceEd(m::GapPic *,
                                                                 double) {
  ++forbidden_advance_calls;
  throw std::runtime_error("forbidden full-time native advance trap");
}
extern "C" void
__wrap__ZN6raikou12mfem_backend6GapPic13advance_emptyEd(m::GapPic *, double) {
  ++forbidden_advance_calls;
  throw std::runtime_error("forbidden full-time native empty advance trap");
}
J read(const std::string &p) {
  std::ifstream in(p);
  J j;
  in >> j;
  if (!in)
    throw std::runtime_error("fixture read");
  return j;
}
void write(const std::string &p, J j) {
  j["sha256"] = sha(j.at("payload").dump());
  std::ofstream out(p);
  out << j.dump(2) << '\n';
  if (!out)
    throw std::runtime_error("fixture write");
}
m::ArenaEvent elastic(m::CollisionScope &s, const std::string &id) {
  m::ArenaEvent e;
  e.id = id;
  e.channel = channels()[0].id;
  e.primary = "electron";
  e.target = "neutral";
  e.event_weight = 4;
  e.scattered_primary_id = e.primary;
  e.phase = s.phase();
  const auto c = channels()[0];
  e.provenance = {m::ArenaProvenance::Kind::ConditionalPhysicalData,
                  c.dataset_sha256,
                  c.kernel_manifest_sha256,
                  c.kernel_binding_file_sha256,
                  c.model,
                  "forced_kernel_integration_only",
                  s.state_prefix(),
                  0};
  e.products = {{e.primary, "e", {1e6, 0, 0}, m::ArenaProduct::At::Primary},
                {id + "/neutral", "n", {0, 0, 0}, m::ArenaProduct::At::Target}};
  return e;
}
int main(int argc, char **argv) {
  try {
    if (argc != 2)
      throw std::runtime_error("output directory");
    const std::string dir = argv[1];
    auto a = fresh();
    const auto base = a->state_prefix();
    const auto b = common(*a);
    const auto e = exposure(b);
    const auto c = completion(e);
    check(!std::is_default_constructible_v<m::FullTimeMotionAdmission> &&
              !std::is_constructible_v<m::FullTimeMotionAdmission, J> &&
              !std::is_copy_constructible_v<m::FullTimeMotionAdmission>,
          "live motion admission cannot be default/JSON constructed or copied");
    int callbacks = 0;
    rejects([&] { a->prepare_compound(1, [&](auto &) { ++callbacks; }); },
            "legacy empty compound refuses full-time");
    rejects(
        [&] {
          a->prepare_compound_with_completion(1, [&](auto &) {
            ++callbacks;
            return c;
          });
        },
        "legacy completion compound refuses full-time");
    rejects([&] { a->prepare_transport(b.dt_s); },
            "legacy transport refuses full-time");
    rejects(
        [&] {
          m::ArenaEvent x;
          a->prepare(x);
        },
        "legacy event refuses full-time");
    check(callbacks == 0 && forbidden_solve_calls == 0 &&
              forbidden_advance_calls == 0,
          "legacy routes call no callbacks/Poisson/advance");
    check(a->state_prefix() == base,
          "legacy refusals preserve accepted prefix");
    a->save(dir + "/genesis.json");
    auto restored = m::PicArena::restore_full_time(
        dir + "/genesis.json", options(), states(), channels(),
        *a->full_time_boundary());
    check(restored->state_prefix() == base, "typed v2 genesis restart exact");
    auto doc = read(dir + "/genesis.json");
    check(doc["payload"].contains("options") &&
              !doc["payload"].contains("options_bytes_hex"),
          "restart uses typed option codec only");
    auto mutate = [&](const char *name, const std::function<void(J &)> &fn) {
      auto changed = doc;
      fn(changed["payload"]);
      const auto path = dir + "/" + name + ".json";
      write(path, changed);
      rejects(
          [&] {
            m::PicArena::restore_full_time(path, options(), states(),
                                           channels(),
                                           *a->full_time_boundary());
          },
          name);
      check(a->state_prefix() == base,
            "corrupt restore does not mutate accepted owner");
    };
    mutate("option_boolean_integer_refused",
           [](J &x) { x["options"]["finite_foils"] = 1; });
    mutate("option_nx_fraction_refused",
           [](J &x) { x["options"]["nx"] = 45.5; });
    mutate("option_nx_negative_refused", [](J &x) { x["options"]["nx"] = -1; });
    mutate("option_nx_overflow_refused", [](J &x) {
      x["options"]["nx"] = std::numeric_limits<std::uint64_t>::max();
    });
    mutate("metadata_revision_fraction_refused",
           [](J &x) { x["metadata"]["revision"] = 0.5; });
    mutate("metadata_revision_negative_refused",
           [](J &x) { x["metadata"]["revision"] = -1; });
    mutate("metadata_journal_bytes_fraction_refused",
           [](J &x) { x["metadata"]["journal_bytes"] = 0.5; });
    mutate("genesis_id_fraction_refused",
           [](J &x) { x["genesis_native"][0]["id"] = 1.5; });
    mutate("canonical_native_id_negative_refused",
           [](J &x) { x["genesis"]["packets"][0]["native_ids"][0] = -1; });
    mutate("next_id_fraction_refused", [](J &x) { x["next_id"] = 3.5; });
    mutate("genesis_next_id_fraction_refused",
           [](J &x) { x["genesis_next_id"] = 3.5; });
    mutate("boundary_revision_fraction_refused",
           [](J &x) { x["boundary"]["accepted_arena_revision"] = 0.5; });
    mutate("boundary_genesis_boolean_integer_refused",
           [](J &x) { x["boundary"]["authored_full_time_genesis"] = 1; });
    mutate("legacy_full_time_schema_refused", [](J &x) {
      x["schema"] = "canonical_native_arena_full_time_exact_v1";
    });
    {
      const auto path = dir + "/duplicate-key.json";
      std::ofstream(path) << "{\"sha256\":\"duplicate\","
                          << doc.dump().substr(1);
      rejects(
          [&] {
            m::PicArena::restore_full_time(path, options(), states(),
                                           channels(),
                                           *a->full_time_boundary());
          },
          "duplicate restart keys refuse");
    }
    {
      auto deep = doc;
      J x = 0;
      for (int i = 0; i < 40; ++i)
        x = J::array({x});
      deep["payload"]["deep"] = x;
      write(dir + "/deep.json", deep);
      rejects(
          [&] {
            m::PicArena::restore_full_time(dir + "/deep.json", options(),
                                           states(), channels(),
                                           *a->full_time_boundary());
          },
          "restart nesting bound refuses before owner construction");
    }
    {
      auto opt = options();
      opt.nx = 90;
      auto different = m::PicArena::create_full_time(opt, states(), channels(),
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
      check(different->state_prefix() != base,
            "actual native options independently bind prefix despite same "
            "declared seals");
    }
    {
      std::vector<m::FullTimeAuthoredPacket> huge;
      for (int i = 0; i < 4; ++i)
        huge.push_back({"huge" + std::to_string(i),
                        "e",
                        "manufactured/overflow-control",
                        {2e179, 0},
                        {.004, .004, .0125},
                        {1e150, 0, 0}});
      rejects(
          [&] {
            m::PicArena::create_full_time(options(), states(), channels(), huge,
                                          genesis());
          },
          "fresh aggregate authored energy overflow refuses");
    }
    {
      auto opt = options();
      opt.nx = opt.ny = 4096;
      rejects(
          [&] {
            m::PicArena::create_full_time(opt, states(), channels(), {},
                                          genesis());
          },
          "immutable geometry allocation bound refuses");
    }
    {
      auto outer = a->begin_full_time_step(b, 8, 32);
      rejects([&] { a->begin_full_time_step(b, 8, 32); },
              "one live private outer per owner");
      auto scope = outer.first_exposure(e);
      const auto before = scope.state_prefix();
      rejects([&] { outer.first_exposure(e); },
              "first exposure reservation rejects sibling handle");
      rejects(
          [&] {
            scope.prepare_compound_with_completion(
                1, [&](auto &) -> m::ArenaStochasticCompletion {
                  rejects([&] { outer.first_exposure(e); },
                          "outer first handle reentry refuses");
                  rejects([&] { outer.fork_transport_trial(); },
                          "transport during collision builder refuses");
                  rejects(
                      [&] {
                        scope.prepare_compound_with_completion(1, [&](auto &) {
                          ++callbacks;
                          return c;
                        });
                      },
                      "nested scoped builder refuses before callback");
                  throw std::runtime_error("manufactured callback failure");
                });
          },
          "builder exception restores private frame without allocation");
      check(scope.state_prefix() == before && a->state_prefix() == base,
            "reentry/exception preserve both states");
      auto bad = c;
      bad.before_rng_clock_sha256 = sha("foreign accepted random");
      rejects(
          [&] {
            scope.prepare_compound_with_completion(1,
                                                   [&](auto &) { return bad; });
          },
          "first completion binds accepted chemistry RNG");
      bad = c;
      bad.common_clock_sha256 = sha("foreign clock");
      rejects(
          [&] {
            scope.prepare_compound_with_completion(1,
                                                   [&](auto &) { return bad; });
          },
          "full-time common clock refuses");
      auto token =
          scope.prepare_compound_with_completion(1, [&](auto &) { return c; });
      check(!token.candidate_full_time_boundary(),
            "private child has no accepted final-boundary authority");
      rejects(
          [&] {
            scope.prepare_compound_with_completion(1, [&](auto &) {
              ++callbacks;
              return c;
            });
          },
          "one unconsumed private population token bounds retained frames");
      check(!a->can_publish(token), "child cannot publish accepted owner");
      auto foreign = fresh();
      auto other = foreign->begin_full_time_step(common(*foreign), 8, 32);
      auto fs = other.first_exposure(exposure(common(*foreign)));
      check(!fs.can_publish(token),
            "exact foreign scope token identity refuses");
      const auto prefix = token.candidate_state_prefix();
      check(scope.try_publish(token),
            "zero-event exposure publishes privately");
      check(scope.state_prefix() == prefix,
            "completed exposure prefix equals exact candidate");
      check(!scope.try_publish(token), "double private publish refuses");
      rejects(
          [&] {
            token.candidate_stochastic_completion(c.integration_owner_id);
          },
          "consumed child retains no duplicate population snapshot");
      check(token.audit().at("kind") == "full_time_exposure",
            "consumed child audit retained");
      auto trial = outer.fork_transport_trial();
      check(trial.state_prefix() == prefix &&
                trial.staged_arena_revision() == 2 &&
                trial.native_census(8).size() == 2,
            "transport starts exact first-exposure native frame");
      for (int i = 0; i < 3; ++i) {
        auto extra = outer.fork_transport_trial();
        check(extra.state_prefix() == prefix,
              "bounded trial retry starts original first frame");
      }
      rejects([&] { outer.fork_transport_trial(); },
              "fifth transport fork refuses before frame allocation");
      m::FullTimeMotionBinding mb{b,
                                  sha("not-provider"),
                                  b.accepted_field_sha256,
                                  sha("newfield"),
                                  sha("newheat"),
                                  genesis().source_sha256};
      std::vector<m::GapParticle> rows;
      for (const auto &row : trial.native_census(8))
        rows.push_back(row.native);
      rejects([&] { trial.prepare_total_field_motion(mb, rows, J::object()); },
              "public JSON cannot authorize live motion");
      check(a->state_prefix() == base,
            "all private motion refusals preserve accepted state");
    }
    check(a->state_prefix() == base,
          "dropping entire private transaction preserves accepted state");
    {
      auto outer = a->begin_full_time_step(b, 2, 32);
      auto scope = outer.first_exposure(e);
      auto token = scope.prepare_compound_with_completion(1, [&](auto &s) {
        rejects([&] { s.prepare(elastic(s, "manufactured/cap")); },
                "reaction output cap rejects before particle copy");
        return c;
      });
      check(scope.try_publish(token),
            "resource refusal leaves zero-event private exposure usable");
    }
    {
      auto outer = a->begin_full_time_step(b, 8, 32);
      auto scope = outer.first_exposure(e);
      const auto before = scope.state_prefix();
      std::optional<m::PreparedArena> escaped;
      rejects(
          [&] {
            scope.prepare_compound_with_completion(2, [&](auto &s) {
              escaped.emplace(s.prepare(elastic(s, "manufactured/escaped")));
              rejects(
                  [&] { s.prepare(elastic(s, "manufactured/second-pending")); },
                  "second pending event token refuses");
              return c;
            });
          },
          "unconsumed event cannot escape compound completion");
      check(!scope.can_publish(*escaped) && scope.state_prefix() == before,
            "escaped child is inert and private base restored");
      escaped.reset();
      auto token = scope.prepare_compound_with_completion(2, [&](auto &s) {
        auto child = s.prepare(elastic(s, "manufactured/elastic"));
        check(s.try_publish(child),
              "manufactured conservative event publishes only privately");
        check(child.audit().at("event_id") == "manufactured/elastic",
              "published child keeps event receipt");
        check(s.rate_population(8).size() == 3,
              "exact target remainder and product carriers preserved");
        return c;
      });
      check(scope.try_publish(token) && a->state_prefix() == base,
            "manufactured event compound leaves accepted owner unchanged");
    }
    {
      auto owner = fresh();
      auto outer = owner->begin_full_time_step(common(*owner), 8, 32);
      auto scope = outer.first_exposure(exposure(common(*owner)));
      auto cc = completion(exposure(common(*owner)));
      auto token =
          scope.prepare_compound_with_completion(1, [&](auto &) { return cc; });
      owner.reset();
      check(!scope.can_publish(token),
            "destroyed owner scope lifetime refuses safely");
    }
    namespace v = m::full_time_detail;
    {
      auto first = c;
      auto second = c;
      auto endpoint = e;
      endpoint.slot = m::FullTimeSlot::after_transport;
      endpoint.frozen_physical_stage_s = b.end_s;
      endpoint.owner_stage_sequence = 1;
      second.stage_sequence = 1;
      second.frozen_physical_stage_s = b.end_s;
      second.before_rng_clock_sha256 = first.after_rng_clock_sha256;
      second.after_rng_clock_sha256 = sha("second-after");
      v::completion_chain_guard(endpoint, second, *a->full_time_boundary(),
                                first);
      check(true, "pure two-slot chemistry chain control");
      auto bad = second;
      bad.integration_owner_id = "different";
      rejects(
          [&] {
            v::completion_chain_guard(endpoint, bad, *a->full_time_boundary(),
                                      first);
          },
          "second chemistry owner mismatch refuses");
      bad = second;
      bad.channel_id = "different";
      rejects(
          [&] {
            v::completion_chain_guard(endpoint, bad, *a->full_time_boundary(),
                                      first);
          },
          "second chemistry channel mismatch refuses");
      bad = second;
      bad.dataset_sha256 = sha("different");
      rejects(
          [&] {
            v::completion_chain_guard(endpoint, bad, *a->full_time_boundary(),
                                      first);
          },
          "second chemistry dataset mismatch refuses");
      bad = second;
      bad.before_rng_clock_sha256 = genesis().rng_sha256;
      rejects(
          [&] {
            v::completion_chain_guard(endpoint, bad, *a->full_time_boundary(),
                                      first);
          },
          "second cannot reuse accepted RNG");
      m::FullTimeCommitSeals seals{sha("field"), sha("heat"), sha("source"),
                                   second.after_rng_clock_sha256};
      v::final_rng_guard(seals, second);
      check(true, "pure final RNG seal control");
      seals.rng_sha256 = first.after_rng_clock_sha256;
      rejects([&] { v::final_rng_guard(seals, second); },
              "final RNG must equal actual second completion");
    }
    {
      auto changed = doc;
      J common_json = {{"accepted_step", b.accepted_step},
                       {"accepted_arena_revision", b.accepted_arena_revision},
                       {"begin_s", b.begin_s},
                       {"end_s", b.end_s},
                       {"dt_s", b.dt_s},
                       {"accepted_arena_sha256", b.accepted_arena_sha256},
                       {"accepted_field_sha256", b.accepted_field_sha256},
                       {"mesh_sha256", b.mesh_sha256},
                       {"material_sha256", b.material_sha256},
                       {"common_clock_sha256", b.common_clock_sha256}};
      auto &meta = changed["payload"]["metadata"];
      J payload = {{"kind", "full_time_outer_begin"},
                   {"common", common_json},
                   {"outer_step_sha256", sha("arena-full-time-outer-v1" + base +
                                             common_json.dump())},
                   {"previous_chain", meta.at("chain")}};
      J row = {{"payload", payload}, {"sha256", sha(payload.dump())}};
      meta["journal"] = J::array({row});
      meta["journal_bytes"] = row.dump().size();
      meta["revision"] = 1;
      meta["chain"] = row.at("sha256");
      write(dir + "/incomplete-new-grammar.json", changed);
      std::string message;
      try {
        m::PicArena::restore_full_time(dir + "/incomplete-new-grammar.json",
                                       options(), states(), channels(),
                                       *a->full_time_boundary());
      } catch (const std::exception &error) {
        message = error.what();
      }
      check(message.find("full-time journal ends inside private step") !=
                std::string::npos,
            "full-time begin recognized by dedicated grammar; incomplete step "
            "refuses");
    }
    {
      auto owner = fresh();
      const auto binding = common(*owner);
      auto outer = owner->begin_full_time_step(binding, 8, 32);
      auto scope = outer.first_exposure(exposure(binding));
      auto cc = completion(exposure(binding));
      rejects(
          [&] {
            scope.prepare_compound_with_completion(1, [&](auto &) {
              owner.reset();
              return cc;
            });
          },
          "owner destroyed during callback refuses without dangling cleanup");
    }
    int allocation_failures = 0, allocation_successes = 0;
    for (long index = 0; index < 48; ++index) {
      auto outer = a->begin_full_time_step(b, 8, 32);
      auto scope = outer.first_exposure(e);
      auto token =
          scope.prepare_compound_with_completion(1, [&](auto &) { return c; });
      const auto before = scope.state_prefix();
      fail_allocation = index;
      const bool ok = scope.try_publish(token);
      fail_allocation = -1;
      if (ok) {
        ++allocation_successes;
      } else {
        ++allocation_failures;
        check(scope.can_publish(token) && scope.state_prefix() == before &&
                  a->state_prefix() == base,
              "allocation failure before private publish leaves token/base "
              "reusable");
        check(scope.try_publish(token),
              "retry after private publication allocation failure succeeds");
      }
    }
    check(allocation_failures > 0 && allocation_successes > 0,
          "allocation failure sweep crosses complete prepublication path");
    check(forbidden_solve_calls == 0 && forbidden_advance_calls == 0 &&
              callbacks == 0,
          "all component controls retain zero native solver/advance/callback "
          "traps");
    check(a->state_prefix() == base,
          "all component controls preserve accepted genesis");
    J output = {{"status", "passed"},
                {"check_count", checks},
                {"checks", results},
                {"allocation_failure_sweep",
                 {{"cases", 48},
                  {"refusals", allocation_failures},
                  {"successes", allocation_successes}}},
                {"forbidden_native_solver_calls", forbidden_solve_calls},
                {"forbidden_native_advance_calls", forbidden_advance_calls},
                {"GPU", false},
                {"full_joint_steps", 0},
                {"genuine_completed_step_restart_qualified", false},
                {"scope",
                 "Manufactured storage and pure contract controls only; no "
                 "fabricated positive Hall JSON or coupled trajectory."}};
    std::ofstream(dir + "/candidate003-tests.json") << output.dump(2) << '\n';
    std::cout << checks << " candidate003 component controls passed\n";
    return 0;
  } catch (const std::exception &x) {
    fail_allocation = -1;
    std::cerr << x.what() << '\n';
    return 1;
  }
}
