#include "raikou/exact_particle_count.h"
#include "raikou/json.hpp"
#include "surface_emission.hpp"
#include <fstream>
#include <iostream>
using C = raikou::mfem_backend::ExactParticleCount;
using J = nlohmann::json;
namespace se = surface_emission;
int main(int argc, char **argv) {
  try {
    if (argc != 3)
      throw std::runtime_error("Count2 cases and output required");
    J input;
    {
      std::ifstream in(argv[1]);
      in >> input;
      if (!in)
        throw std::runtime_error("input");
    }
    J cases = J::array(), controls = J::array();
    for (const auto &item : input) {
      J row = {{"id", item["id"]}, {"hi", item["hi"]}, {"lo", item["lo"]}};
      try {
        C count(item["hi"].get<double>(), item["lo"].get<double>());
        row["accepted"] = true;
        row["carriers"] = count.positive_carriers();
        row["rounded"] = count.rounded_for_rate();
        row["export_two_limbs"] = count.require_two_limb_values();
        if (item.contains("debit")) {
          row["debit"] = item["debit"];
          try {
            auto after = count.debit(item["debit"]);
            row["debit_accepted"] = true;
            row["after_carriers"] = after.positive_carriers();
            row["round_after"] = after.rounded_for_rate();
            try {
              row["after_two_limbs"] = after.require_two_limb_values();
            } catch (const std::exception &e) {
              row["after_export_refused"] = e.what();
            }
          } catch (const std::exception &e) {
            row["debit_accepted"] = false;
            row["debit_refusal"] = e.what();
          }
        }
      } catch (const std::exception &e) {
        row["accepted"] = false;
        row["refusal"] = e.what();
      }
      cases.push_back(row);
    }
    const double volume = 25e-6 * 25e-6 * .025, patch_area = 25e-6 * .025;
    const std::vector<std::pair<std::string, double>> inventories = {
        {"background_cell_n1e18_m-3", 1e18 * volume},
        {"background_cell_n1e21_m-3", 1e21 * volume},
        {"background_cell_n1e24_m-3", 1e24 * volume},
        {"source_patch_theta006", .06 * 1.3e19 * patch_area},
        {"source_patch_theta033", .33 * 1.3e19 * patch_area},
        {"source_patch_theta100", 1.3e19 * patch_area},
        {"whole_gap_theta033", .33 * 1.3e19 * .005 * .025},
        {"stress_only_1e20", 1e20}};
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
          J control = {{"T_K", T},
                       {"F_V_m", F},
                       {"B_T", B},
                       {"source_area_m2", stage.area_m2},
                       {"interval_s", stage.dt},
                       {"native_cell_volume_m3", volume},
                       {"source_patch_area_m2", patch_area},
                       {"nominal_barrier_guard",
                        policy.maximum_quantum_magnetic_energy_fraction}};
          try {
            std::vector<double> constant, varied;
            const auto f = se::prepare(stage, config);
            const auto p = se::begin_layer(
                f, policy, f.expected_physical_electrons, .6, .35, .1, 1);
            for (int step = 0; step < 16; ++step) {
              constant.push_back(p.state.particle.weight);
              auto s = stage;
              s.T_K = T * (1 + 1e-4 * step);
              const auto flux = se::prepare(s, config);
              varied.push_back(se::begin_layer(flux, policy,
                                               flux.expected_physical_electrons,
                                               .6, .35, .1, step + 2)
                                   .state.particle.weight);
            }
            control["admitted_source"] = true;
            control["source_weight"] = p.state.particle.weight;
            control["sequences"] = J::array();
            for (int which = 0; which < 2; ++which)
              for (const auto &[name, initial] : inventories) {
                C count(initial, 0);
                J records = J::array();
                const auto &weights = which ? varied : constant;
                for (std::size_t i = 0; i < weights.size(); ++i) {
                  auto before = count;
                  count = count.debit(weights[i]);
                  J receipt = {{"step", i},
                               {"weight", weights[i]},
                               {"before", before.positive_carriers()},
                               {"after", count.positive_carriers()}};
                  try {
                    receipt["two_limb_export"] =
                        count.require_two_limb_values();
                  } catch (const std::exception &) {
                    receipt["two_limb_export_refused"] = true;
                  }
                  records.push_back(receipt);
                }
                control["sequences"].push_back(
                    {{"inventory", name},
                     {"initial", initial},
                     {"sequence", which ? "heterogeneous_admitted_source"
                                        : "constant_admitted_source"},
                     {"records", records},
                     {"final", count.positive_carriers()}});
              }
          } catch (const std::exception &e) {
            control["admitted_source"] = false;
            control["reason"] = e.what();
          }
          controls.push_back(control);
        }
    J report = {
        {"status", "completed"},
        {"count2_cases", cases},
        {"source_controls", controls},
        {"capacity_bits", C::maximum_bits},
        {"maximum_positive_carriers", C::maximum_carriers},
        {"inventory_scope",
         "illustrative background densities and heat-v2 ideal-oxide coverage "
         "assumptions, not measured shot inventories;1e20 stress-only"},
        {"source_scope",
         "actual sealed source-kernel admitted weights at declared controls, "
         "forced count-domain debits only; no collision selection/trajectory"},
        {"no_sampling_policy_change", true}};
    std::ofstream out(argv[2]);
    out << report.dump(2) << '\n';
    if (!out)
      throw std::runtime_error("output");
    std::cout << cases.size() << " Count2 inputs and " << controls.size()
              << " source controls completed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
