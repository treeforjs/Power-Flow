#include "reviewed-collision_scheduler.hpp"
#include <iostream>
#include <iomanip>
using namespace raikou::pic_collisions;
int main(){const std::string H(64,'a'),K(64,'b');Packet p{"p","e",1,1e200,9.1093837139e-31,{3.5774131702061213e+208,0,0}},n{"n","H2O",1,1e200,2.9915e-26,{2.99792458e+208,0,0}};auto k=relative(p,n);Channel ch{"c","e","H2O",H,K,Observable::ExclusiveEvent,[](double E){require(E>7000&&E<9000,"8keV support");return 1e-20;}};auto rows=candidates({p},{n},{Cell{1,1.}},{ch});std::cout<<std::setprecision(17)<<"native_energy_eV="<<k.target_frame_energy_eV<<" native_flux="<<k.moller_flux_m_s<<" candidates="<<rows.size()<<'\n';}
