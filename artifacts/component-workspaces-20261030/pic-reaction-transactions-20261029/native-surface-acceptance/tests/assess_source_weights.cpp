#include "raikou/mfem_pic_arena.h"
#include "surface_emission.hpp"
#include <fstream>
#include <iostream>
#include <sched.h>
#include <sys/resource.h>
namespace m = raikou::mfem_backend;
namespace se = surface_emission;
using J = m::ArenaJson;
int main(int argc, char **argv) {
  try {
    if (argc != 2)
      throw std::runtime_error("report path required");
    cpu_set_t mask;
    CPU_ZERO(&mask);
    if (sched_getaffinity(0, sizeof(mask), &mask) || CPU_COUNT(&mask) != 1 ||
        !CPU_ISSET(7, &mask) || getpriority(PRIO_PROCESS, 0) != 15)
      throw std::runtime_error("CPU7/nice15 required");
    J rows = J::array();
    for (double T : {300., 700., 1500.})
      for (double F : {0., 1e7, 1e8, 1e9})
        for (double B : {0., 30.}) {
          se::Stage stage{1,         1,     1,       1,    1,    1,
                          0,         2e-12, 0,       0,    0,    {0, -F, 0},
                          {0, 1, 0}, T,     1.25e-4, true, true, true};
          stage.total_B_T = {B, 0, 0};
          se::Config config{4.3};
          se::LayerPolicy policy;
          policy.local_linear_vacuum_extent_m = 25e-6;
          J row = {
              {"T_K", T},
              {"F_V_m", F},
              {"B_T", B},
              {"area_m2", stage.area_m2},
              {"duration_s", stage.dt},
              {"barrier_guard",
               policy.maximum_quantum_magnetic_energy_fraction},
              {"source_macro_packets_per_stage", 1},
              {"debit_is_forced_domain_stress_not_sampled_collision", true}};
          try {
            const auto flux = se::prepare(stage, config);
            const auto source = se::begin_layer(
                flux, policy, flux.expected_physical_electrons, .6, .35, .1, 1);
            row["admitted_source"] = true;
            row["source_weight"] = source.state.particle.weight;
            row["source_current_A_m2"] = flux.spectrum.current_A_m2;
            row["source_log_current_A_m2"] = flux.spectrum.log_current_A_m2;
            row["source_status"] = int(source.state.status);
            row["source_exchange_residual_J"] = source.exchange.residual_J();
            row["count_stresses"] = J::array();
            for (const std::string sequence :
                 {"constant_accepted_source_weight",
                  "heterogeneous_accepted_source_weights"}) {
              m::ArenaCount count{1e20, 0};
              J receipts = J::array();
              bool stopped = false;
              std::string reason;
              for (int step = 0; step < 16; ++step) {
                try {
                  auto varied = stage;
                  if (sequence[0] == 'h')
                    varied.T_K = T * (1 + 1e-4 * step);
                  auto f = se::prepare(varied, config);
                  auto particle =
                      se::begin_layer(f, policy, f.expected_physical_electrons,
                                      .6, .35, .1, step + 2);
                  const double weight = particle.state.particle.weight;
                  const auto before = count;
                  auto after = m::PicArena::debit(count, weight);
                  auto carriers = m::PicArena::positive_carriers(after);
                  receipts.push_back({{"step", step},
                                      {"source_T_K", varied.T_K},
                                      {"event_weight", weight},
                                      {"before", {before.hi, before.lo}},
                                      {"after", {after.hi, after.lo}},
                                      {"positive_native_carriers", carriers}});
                  count = after;
                } catch (const std::exception &e) {
                  stopped = true;
                  reason = e.what();
                  break;
                }
              }
              row["count_stresses"].push_back(
                  {{"sequence", sequence},
                   {"initial_neutral_macro_count", 1e20},
                   {"maximum_attempts", 16},
                   {"accepted_debits", receipts.size()},
                   {"stopped", stopped},
                   {"reason", reason},
                   {"final_count", {count.hi, count.lo}},
                   {"receipts", receipts}});
            }
          } catch (const std::exception &e) {
            row["admitted_source"] = false;
            row["reason"] = e.what();
          }
          rows.push_back(row);
        }
    J report = {
        {"status", "completed"},
        {"cases", rows},
        {"source_kernel", "sealed usable-api-v1 called by heat-v2"},
        {"native_arena_base", "sealedv6.1 Count2"},
        {"physics_scope",
         "conditional component controls; actual prepare/begin_layer weights, "
         "not actual C19 fields or a coupled trajectory pilot"},
        {"selection_scope", "forced debit stress only, no collision "
                            "probability or cross section sampled"},
        {"no_floor_clamp_rounding_fallback", true},
        {"GPU", false},
        {"PIC_evolution", false}};
    std::ofstream out(argv[1]);
    out << report.dump(2) << '\n';
    if (!out)
      throw std::runtime_error("write");
    std::cout << rows.size() << " admitted-source weight controls assessed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
