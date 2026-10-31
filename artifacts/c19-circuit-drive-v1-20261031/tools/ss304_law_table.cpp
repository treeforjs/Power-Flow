// Tabulate the Hall SS304 law (heat-v2 material-v2 ss304v2::Model, evaluated_solid_1500, fixed reference volume,
// rho_ref 8000 kg/m3: the arguments of Hall corrective v3/v3.1 Model) for the reduced Hall vacuum load.
// Output CSV: T_K, sigma_S_m, capacity_J_m3_K, energy_J_m3 (read-only include of the sealed header).
#include "ss304_thermal_v2.hpp"
#include <cstdio>
int main(){ss304v2::Model m(ss304v2::Qualification::evaluated_solid_1500,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,8000.);
    std::printf("T_K,sigma_S_m,capacity_J_m3_K,energy_J_m3\n");
    for(int i=0;i<=1200;++i){double T=300.+i;auto v=m.evaluate(T);std::printf("%.17g,%.17g,%.17g,%.17g\n",T,v.sigma_s_m,v.capacity_j_m3_reference_k,v.energy_j_m3_reference);}
    return 0;}
