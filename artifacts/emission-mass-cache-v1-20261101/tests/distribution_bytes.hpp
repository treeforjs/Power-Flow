#pragma once
#include "surface_emission.hpp"
#include <bit>
#include <cstdint>
#include <string>
inline std::string distribution_bytes(const surface_emission::Distribution&d){
 std::string b;auto word=[&](std::uint64_t x){for(unsigned i=0;i<8;++i)b.push_back(char(x>>(8*i)));};
 auto number=[&](double x){word(std::bit_cast<std::uint64_t>(x));};const auto&c=d.config;
 number(c.work_function_eV);word(c.energy_panels);word(c.action_nodes);number(c.tail_widths);number(c.electron_mass_kg);word(c.electron_mass_identity.size());b+=c.electron_mass_identity;
 for(double x:{d.T_K,d.extraction_field_V_m,d.current_A_m2,d.log_current_A_m2,d.mean_metal_energy_eV,d.wall_heat_W_m2,d.barrier_lowering_eV,d.fermi_decay_width_eV})number(x);
 word(d.has_log_current);word(d.current_underflow);word(d.nodes.size());for(const auto&n:d.nodes)for(double x:{n.normal_energy_eV,n.probability,n.cdf,n.mean_transverse_eV})number(x);return b;
}
