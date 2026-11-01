#include "fast_emission.hpp"
#include "distribution_bytes.hpp"
#include <iostream>
int main(){surface_emission::Config c;c.work_function_eV=4.3;c.energy_panels=32;c.action_nodes=32;for(double T:{1.,300.,1000.,3000.})for(double F:{-3.,-0.,0.,1.,300.,2e6,3e7}){const auto b=distribution_bytes(electrode_heat_v6_1::emission::evaluate(T,F,c));std::cout.write(b.data(),b.size());}surface_emission::Config defaults;defaults.work_function_eV=4.3;for(double F:{-3.,300.,2e6}){const auto b=distribution_bytes(electrode_heat_v6_1::emission::evaluate(1000.,F,defaults));std::cout.write(b.data(),b.size());}return !std::cout;}
