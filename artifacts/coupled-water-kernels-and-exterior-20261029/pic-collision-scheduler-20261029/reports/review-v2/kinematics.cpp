#include "reviewed-collision_scheduler.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace raikou::pic_collisions;
int main(int argc,char**argv){if(argc!=2)return 2;std::ifstream f(argv[1]);std::string line;std::cout<<std::setprecision(17)<<"label,status,energy_eV,moller_m_s\n";while(std::getline(f,line)){std::istringstream s(line);std::string label;Packet p,n;std::getline(s,label,',');p.mass_kg=9.1093837139e-31;n.mass_kg=2.9915e-26;for(int i=0;i<6;++i){std::string x;std::getline(s,x,',');(i<3?p.proper_velocity_m_s[i]:n.proper_velocity_m_s[i-3])=std::stod(x);}try{auto k=relative(p,n);std::cout<<label<<",ok,"<<k.target_frame_energy_eV<<','<<k.moller_flux_m_s<<'\n';}catch(const Error&){std::cout<<label<<",rejected,0,0\n";}}}
