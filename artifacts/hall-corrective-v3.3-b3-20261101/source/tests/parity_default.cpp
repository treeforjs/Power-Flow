#include "root_api.hpp"
namespace hc=hall_corrective_v3;
#include "parity_scenarios.hpp"
#include <fstream>
#include <iostream>
int main(int argc,char**argv){try{if(argc!=3)throw std::runtime_error("TRACE OUT required");auto rows=parity::run(argv[1]);std::ofstream(argv[2])<<rows.dump(1)<<'\n';return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
