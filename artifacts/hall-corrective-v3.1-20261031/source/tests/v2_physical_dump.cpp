// Physical-state parity reference: the parity scenarios run on the UNCHANGED Hall
// corrective v2 sources (hall-corrective-v2-20261030/source, read in place).
//   v2_physical_dump TRACE OUTFILE
#include "root_api.hpp"
namespace hc = hall_corrective_v2;
#include "parity_scenarios.hpp"
#include <fstream>
#include <iostream>

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::runtime_error("usage: v2_physical_dump TRACE OUTFILE");
        auto rows = parity::run(argv[1]);
        std::ofstream(argv[2]) << rows.dump(1) << '\n';
        return 0;
    } catch (const std::exception& e) { std::cerr << "fatal: " << e.what() << '\n'; return 2; }
}
