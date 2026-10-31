// Reference values of the outer3 parent rate from the dataset's own
// implementation (h2o-collision-data-20261029 include tree: collision_data.h
// with orbital weights, h2o_outer_valence.h). Built separately from the
// coverage-v2 sources; run_coverage_v2_checks.py compares its output with the
// coverage-v2 evaluation bit for bit.
#include "raikou/h2o_outer_valence.h"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

std::vector<double> grid() {
  std::vector<double> e{1000., std::nextafter(1000., 2000.), 1000.5, 1500., 2000., 3000., 5000., 8000., 8084.470604111418,
                        9999.999, std::nextafter(10000., 0.), 10000.};
  for (int i = 0; i <= 200; ++i) e.push_back(1000. * std::pow(10., i / 200.));
  return e;
}
std::string hex(double x) { char b[64];std::snprintf(b, sizeof b, "%a", x);return b; }

int main(int argc, char** argv) {
  try {
    if (argc != 3) throw std::runtime_error("outer3-dataset output-json");
    const auto data = raikou::kinetics::StrictCollisionDataset::read(argv[1]);
    raikou::Json out = raikou::Json::array();
    for (double E : grid()) {
      if (E > 10000.) E = 10000.; // the last grid point rounds above 10 keV
      const auto r = raikou::kinetics::conditional_outer3_ionization_rates(data, E);
      out.push_back({{"energy", hex(E)}, {"parent", hex(r.parent_m2)}, {"valence", hex(r.four_valence_total_m2)},
                     {"unmapped_2a1", hex(r.unmapped_2a1_m2)}, {"OHp_H", hex(r.OHplus_H_m2)}, {"Hp_OH", hex(r.Hplus_OH_m2)}});
    }
    std::ofstream(argv[2]) << out.dump(1) << '\n';
    std::cout << "wrote " << out.size() << " reference rows\n";
  } catch (const std::exception& e) { std::cerr << e.what() << '\n';return 1; }
}
