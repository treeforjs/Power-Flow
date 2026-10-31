#pragma once
// Heat/source owner v7: sealed thermal-field emission evaluate with the temperature domain extended to
// T_max (3000 <= T_max <= 5000 K); see emission_hot_v7.cpp. Returns the sealed Distribution type.
#include "surface_emission.hpp"
namespace electrode_heat_v7::emission_hot {
surface_emission::Distribution evaluate(double temperature_K, double signed_extraction_field_V_m, const surface_emission::Config&, double T_max_K);
}
