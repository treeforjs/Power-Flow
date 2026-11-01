#include "atomic_elastic.hpp"
#include <iomanip>
#include <iostream>
#include <type_traits>
using namespace atomic_elastic;
constexpr auto diagnostic=Mode::diagnostic_reported_total_normalized_shape;
unsigned controls=0;
void check(bool b,const char*s){++controls;if(!b)throw std::runtime_error(s);}
template<class F>void refuses(F f){bool bad=false;try{f();}catch(const std::exception&){bad=true;}check(bad,"missing refusal");}
void vec(pfc::V3 v){std::cout<<'['<<v.x<<','<<v.y<<','<<v.z<<']';}
void emit(const Event&r,const char*type,int Z,double mass,unsigned i){
 std::cout<<"{\"type\":\""<<type<<"\",\"Z\":"<<Z<<",\"M\":"<<mass<<",\"i\":"<<i<<",\"E\":"<<r.E_tf<<",\"mu\":"<<r.mu_tf<<",\"sigma_m2\":"<<r.sigma_m2<<",\"recoil\":"<<r.recoil_tf_eV<<",\"incoming\":[";
 vec(r.incoming_u[0]);std::cout<<',';vec(r.incoming_u[1]);std::cout<<"],\"outgoing\":[";vec(r.outgoing_u[0]);std::cout<<',';vec(r.outgoing_u[1]);std::cout<<"],\"K0\":["<<r.incoming_ke_eV[0]<<','<<r.incoming_ke_eV[1]<<"],\"K1\":["<<r.outgoing_ke_eV[0]<<','<<r.outgoing_ke_eV[1]<<"]}\n";
}
int main(int argc,char**argv){
 check(argc==5,"arguments");std::cout<<std::setprecision(17);
 static_assert(!std::is_default_constructible_v<Model>);
 const auto m=Model::load(argv[1],argv[2],argv[3]);
 if(std::string(argv[4])=="load"){std::cout<<m.identity()<<'\n';return 0;}
 const pfc::V3 z{0,0,1},zero{};
 const pfc::ChemKey key{123,4,5,PFC_PURPOSE_FIXTURE_TEST};
 const double nan=std::numeric_limits<double>::quiet_NaN(),inf=std::numeric_limits<double>::infinity();
 std::cout<<"{\"type\":\"identity\",\"model_sha256\":\""<<m.identity()<<"\",\"registry_sha256\":\""<<m.registry_blob_sha256()<<"\",\"electron_mass_eV\":"<<m.electron_mass_eV()<<"}\n";
 for(double E:{0.,-1.,999.,std::nextafter(1000.,0.),std::nextafter(10000.,inf),10001.,nan,inf})refuses([&]{m.cross_section(Target::H,E,diagnostic);});
 refuses([&]{m.cross_section(Target::H,1000,Mode(0));});refuses([&]{m.cross_section(Target(2),1000,diagnostic);});
 for(double q:{-.001,1.001,nan,inf})refuses([&]{m.quantile(Target::H,1000,q,diagnostic);});
 for(double mu:{-1.001,1.001,nan,inf})refuses([&]{m.cdf(Target::H,1000,mu,diagnostic);});
 refuses([&]{m.controlled_event(Target::H,1000,zero,zero,0,0,0,diagnostic);});
 refuses([&]{m.controlled_event(Target::H,1000,{0,0,2},zero,0,0,0,diagnostic);});
 refuses([&]{m.controlled_event(Target::H,1000,z,{inf,0,0},0,0,0,diagnostic);});
 refuses([&]{m.controlled_event(Target::H,1000,z,zero,1,0,0,diagnostic);});
 refuses([&]{m.controlled_event(Target::H,1000,z,zero,nan,0,0,diagnostic);});
 refuses([&]{m.controlled_event(Target::H,1000,z,zero,0,1.01,0,diagnostic);});
 refuses([&]{m.controlled_event(Target::H,1000,z,zero,0,0,1,diagnostic);});
 auto bad=key;bad.purpose=PFC_PURPOSE_MCC;refuses([&]{m.sample(Target::H,1000,z,zero,0,bad,0,diagnostic);});
 bad=key;bad.carrier=PFC_CARRIER_LIMIT;refuses([&]{m.sample(Target::H,1000,z,zero,0,bad,0,diagnostic);});
 bad=key;bad.step=PFA_STEP_LIMIT;refuses([&]{m.sample(Target::H,1000,z,zero,0,bad,0,diagnostic);});
 refuses([&]{m.sample(Target::H,1000,z,zero,0,key,PFC_EVENT_LIMIT,diagnostic);});
 for(const auto&s:m.summaries())std::cout<<"{\"type\":\"summary\",\"Z\":"<<(s.target==Target::H?1:8)<<",\"E\":"<<s.E<<",\"reported\":"<<s.reported<<",\"integrated\":"<<s.integrated<<",\"transport\":"<<s.transport<<",\"difference\":"<<s.relative_integral_difference<<"}\n";
 unsigned events=0;
 for(Target t:{Target::H,Target::O}){
  const int Z=t==Target::H?1:8;
  for(unsigned k=0;k<19;++k){const double E=1000+500*k;
   std::cout<<"{\"type\":\"cross_section\",\"Z\":"<<Z<<",\"E\":"<<E<<",\"sigma_m2\":"<<m.cross_section(t,E,diagnostic)<<"}\n";
   for(unsigned j=0;j<=1000;++j){const double q=j/1000.;const double mu=m.quantile(t,E,q,diagnostic),got=m.cdf(t,E,mu,diagnostic);check(std::abs(got-q)<1e-11,"CDF inverse");}
   check(m.quantile(t,E,0,diagnostic)==-1&&m.quantile(t,E,1,diagnostic)==1,"angular endpoints");
   for(double mu:{-1.,-.9,0.,.5,1.}){const auto r=m.controlled_event(t,E,z,zero,0,mu,.375,diagnostic);emit(r,"analytic",Z,m.target_mass_eV(t),events++);}
  }
  for(unsigned i=0;i<4096;++i){const double E=1000+500*(i%19);
   const pfc::V3 d=i%2?pfc::V3{.6,0,.8}:z;
   const pfc::V3 u=i%3?pfc::V3{12345,-67890,100000}:zero;
   const pfc::ChemKey k{0xabcdef1234567890ULL,17,1+i,PFC_PURPOSE_FIXTURE_TEST};
   const auto r=m.sample(t,E,d,u,0,k,0,diagnostic),repeat=m.sample(t,E,d,u,0,k,0,diagnostic);
   check(r.mu_tf==repeat.mu_tf&&r.outgoing_ke_eV==repeat.outgoing_ke_eV,"deterministic replay");
   const double q=pfc::chem_u01(k,0,PFC_W_KIN);check(std::abs(m.cdf(t,E,r.mu_tf,diagnostic)-q)<1e-11,"sample word and CDF");
   emit(r,"event",Z,m.target_mass_eV(t),events++);
  }
  // Independent fixed-energy population for distribution sanity, not source accuracy.
  for(unsigned i=0;i<8192;++i){const pfc::ChemKey k{7654321,9,1+i,PFC_PURPOSE_FIXTURE_TEST};
   std::cout<<"{\"type\":\"angle\",\"Z\":"<<Z<<",\"mu\":"<<m.quantile(t,5000,pfc::chem_u01(k,0,PFC_W_KIN),diagnostic)<<"}\n";}
 }
 std::cout<<"{\"type\":\"controls\",\"assertions\":"<<controls<<",\"events\":"<<events<<"}\n";
}
