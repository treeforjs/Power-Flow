#include "source/water_parent_bridge.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace raikou::pic_water;
namespace pc=raikou::pic_collisions;
void vec(std::ostream& o,const Vec& v){o<<'['<<v[0]<<','<<v[1]<<','<<v[2]<<']';}
int main(int argc,char** argv){
  if(argc!=5)return 2;
  const auto data=std::make_shared<const Data>(Data::read(argv[1]));
  ParentBridge bridge(data,parent_model,argv[2],argv[3],2.9915e-26);
  std::ofstream o(argv[4]);o<<std::setprecision(17);const double mn=2.9915e-26;
  std::size_t index=0;
  for(double E:{14.,100.,999.})for(double boost:{0.,1e-6,1e-4,1e-3,.8,1e5})
  for(int shape=0;shape<4;++shape){
    const double q=E*pc::qe/(me*pc::c*pc::c),u=pc::c*std::sqrt(q*(q+2));
    Vec direction{},U{};
    if(shape==0){direction={1,0,0};U={0,boost*pc::c,0};}
    if(shape==1){direction={1,0,0};U={boost*pc::c,0,0};}
    if(shape==2){direction={.6,.8,0};U={.3*boost*pc::c,-.4*boost*pc::c,std::sqrt(.75)*boost*pc::c};}
    if(shape==3){direction={0,1,0};U={.6*boost*pc::c,.8*boost*pc::c,0};}
    Vec rest{};for(int d=0;d<3;++d)rest[d]=u*direction[d];
    pc::Packet e{"e","electron",1,.25,me,to_lab(rest,U)},n{"n","H2O",1,64,mn,U};
    o<<"{\"index\":"<<index++<<",\"input_eV\":"<<E<<",\"boost_u_over_c\":"<<boost<<",\"shape\":"<<shape<<",\"e_u\":";vec(o,e.proper_velocity_m_s);o<<",\"n_u\":";vec(o,U);
    try{
      const auto kin=pc::relative(e,n);o<<",\"energy_eV\":"<<kin.target_frame_energy_eV<<",\"flux_m_s\":"<<kin.moller_flux_m_s;
      const auto dir=incident_direction(e,n,kin.target_frame_energy_eV);o<<",\"direction\":";vec(o,dir);
      const auto rows=pc::candidates({e},{n},{{1,1e-12}},{bridge.channel()},1,1,1);
      if(rows.empty())throw pc::Error("empty rows");
      const auto out=bridge.sample(rows[0],e,n,{0,0,1});
      o<<",\"accepted\":true,\"energy_defect_J\":"<<out.energy_defect_J<<",\"energy_bound_J\":"<<out.energy_roundoff_limit_J<<",\"momentum_bound\":"<<out.momentum_roundoff_limit_kg_m_s<<",\"products\":[";
      for(int i=0;i<3;++i){if(i)o<<',';const auto& p=out.particles[i];o<<"{\"mass\":"<<p.mass_kg<<",\"weight\":"<<p.weight<<",\"reference_J\":"<<p.reference_J<<",\"u\":";vec(o,p.proper_velocity_m_s);o<<'}';}o<<']';
    }catch(const std::exception& ex){o<<",\"accepted\":false,\"exception\":\""<<ex.what()<<'"';}
    o<<"}\n";
  }
  pc::Packet e{"e","electron",1,.25,me,{pc::c*std::sqrt((100*pc::qe/(me*pc::c*pc::c))*(2+100*pc::qe/(me*pc::c*pc::c))),0,0}},n{"n","H2O",1,64,mn,{0,0,0}};
  const auto row=pc::candidates({e},{n},{{1,1e-12}},{bridge.channel()})[0];
  auto mutate=[&](const char* label,auto e2,auto n2,auto event){o<<"{\"mutation\":\""<<label<<'"';try{const auto p=bridge.sample(event,e2,n2,{0,0,1});o<<",\"accepted\":true,\"energy_bound_J\":"<<p.energy_roundoff_limit_J<<",\"energy_defect_J\":"<<p.energy_defect_J<<",\"momentum_bound\":"<<p.momentum_roundoff_limit_kg_m_s<<",\"momentum_defect\":";vec(o,p.momentum_defect_kg_m_s);}catch(const std::exception& ex){o<<",\"accepted\":false,\"exception\":\""<<ex.what()<<'"';}o<<"}\n";};
  auto e2=e;std::swap(e2.proper_velocity_m_s[0],e2.proper_velocity_m_s[1]);mutate("rotated_projectile_same_rate",e2,n,row);
  auto n2=n;n2.mass_kg*=2;mutate("changed_target_mass_same_rate",e,n2,row);
  auto event=row;event.volume_m3*=2;mutate("changed_volume_without_rate",e,n,event);
  event=row;event.volume_m3*=2;event.rate_s/=2;mutate("changed_volume_with_rate",e,n,event);
  e2=e;n2=n;e2.weight=1e-300;n2.weight=64e-300;
  event=pc::candidates({e2},{n2},{{1,1e-12}},{bridge.channel()})[0];mutate("tiny_weight",e2,n2,event);
  auto axis_probe=[&](const char* label,const pc::Candidate& a,const pc::Packet& ep,const pc::Packet& np){o<<"{\"mutation\":\""<<label<<'"';try{const auto out=bridge.sample(a,ep,np,{0,0,1+4e-13});o<<",\"accepted\":true,\"energy_bound_J\":"<<out.energy_roundoff_limit_J<<",\"energy_defect_J\":"<<out.energy_defect_J;}catch(const std::exception& ex){o<<",\"accepted\":false,\"exception\":\""<<ex.what()<<'"';}o<<"}\n";};
  axis_probe("near_unit_axis_ordinary_weight",row,e,n);
  axis_probe("near_unit_axis_tiny_weight",event,e2,n2);
  for(double weight:{1e-292,1e-289,1e-284,1e-280}){e2=e;n2=n;e2.weight=weight;n2.weight=64*weight;event=pc::candidates({e2},{n2},{{1,1e-12}},{bridge.channel()})[0];std::string name="valid_weight_"+std::to_string(static_cast<int>(std::log10(weight)));mutate(name.c_str(),e2,n2,event);}

}
