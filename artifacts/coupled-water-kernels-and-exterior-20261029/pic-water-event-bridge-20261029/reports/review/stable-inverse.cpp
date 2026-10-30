#include "source/water_parent_bridge.hpp"
#include <fstream>
#include <iomanip>
using namespace raikou::pic_water;
namespace pc=raikou::pic_collisions;
void vec(std::ostream& o,const Vec& v){o<<'['<<v[0]<<','<<v[1]<<','<<v[2]<<']';}
int main(int argc,char** argv){
 if(argc!=2)return 2;
 std::ofstream o(argv[1]);o<<std::setprecision(17);
 for(double E:{14.,100.,999.})for(double boost:{0.,.8,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e10,1e12})for(int shape=0;shape<4;++shape){
  const double x=E*pc::qe/(me*pc::c*pc::c),speed=pc::c*std::sqrt(x*(x+2));Vec direction{},U{};
  if(shape==0){direction={1,0,0};U={0,boost*pc::c,0};}
  if(shape==1){direction={1,0,0};U={boost*pc::c,0,0};}
  if(shape==2){direction={.6,.8,0};U={.3*boost*pc::c,-.4*boost*pc::c,std::sqrt(.75)*boost*pc::c};}
  if(shape==3){direction={0,1,0};U={.6*boost*pc::c,.8*boost*pc::c,0};}
  Vec rest{};for(int i=0;i<3;++i)rest[i]=speed*direction[i];
  pc::Packet e{"e","electron",1,.25,me,to_lab(rest,U)},n{"n","H2O",1,64,2.9915e-26,U};
  const double Eactual=pc::relative(e,n).target_frame_energy_eV;
  const long double grel=1+static_cast<long double>(Eactual)*pc::qe/(static_cast<long double>(me)*pc::c*pc::c);
  const long double k=(gamma(e.proper_velocity_m_s)+grel)/(gamma(U)+1);
  std::array<long double,3> r{};long double norm2=0;
  for(int i=0;i<3;++i){r[i]=std::fma(-k,static_cast<long double>(U[i]),static_cast<long double>(e.proper_velocity_m_s[i]));norm2+=r[i]*r[i];}
  Vec dir{};for(int i=0;i<3;++i)dir[i]=static_cast<double>(r[i]/std::sqrt(norm2));
  o<<"{\"input_eV\":"<<E<<",\"boost_u_over_c\":"<<boost<<",\"shape\":"<<shape<<",\"e_u\":";vec(o,e.proper_velocity_m_s);o<<",\"n_u\":";vec(o,U);o<<",\"direction\":";vec(o,dir);o<<"}\n";
 }
}
