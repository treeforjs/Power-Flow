#pragma once
#include "chem_host.hpp"
#include "source_pins.hpp"
#include <array>
#include <bit>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>

namespace atomic_elastic {
inline void need(bool b,const char*s){if(!b)throw std::invalid_argument(std::string("atomic elastic: ")+s);}
enum class Target {H=0,O=1};
enum class Mode {diagnostic_reported_total_normalized_shape=1};
inline constexpr double a0_squared_m2=2.8002852e-21;
inline constexpr char model_contract[]=
 "atomic-elastic-nist64-v1;calculated-free-neutral-ground-atoms;source-reported-total;linear-energy;"
 "181-point-DCS-linear-mu;normalized-angular-shape;full-minus1-plus1;no-extrapolation;"
 "target-rest-polar-angle;finite-relativistic-two-body-recoil;uniform-azimuth;"
 "RNG-layout2-fixture-purpose-words32,33;no-internal-energy;diagnostic-only;incomplete-network";
struct Summary {Target target;double E,reported,integrated,transport,relative_integral_difference;};
struct Event {
 std::array<pfc::V3,2> incoming_u{},outgoing_u{};
 std::array<double,2> incoming_ke_eV{},outgoing_ke_eV{};
 double E_tf=0,mu_tf=0,sigma_m2=0,recoil_tf_eV=0,energy_residual_eV=0;
 pfc::V3 momentum_residual_eV_c{};
 std::string model_sha256;
};
class Model {
 struct Table {double total=0;std::array<double,181>dcs{};};
 struct Shape {double reported=0,integrated=0,transport=0;std::array<double,181>dcs{},cumulative{};};
 std::array<std::array<Table,10>,2>tables_{};
 std::array<double,181>mu_{};
 double electron_mass_=0;std::array<double,2>target_mass_{};
 std::string identity_,registry_blob_sha_;
 Model()=default;
 static unsigned species(Target t){need(t==Target::H||t==Target::O,"unsupported target state");return unsigned(t);}
 static void mode(Mode m){need(m==Mode::diagnostic_reported_total_normalized_shape,"explicit diagnostic model required");}
 static bool finite(pfc::V3 v){return std::isfinite(v.x)&&std::isfinite(v.y)&&std::isfinite(v.z)&&std::isfinite(pfc::norm2(v));}
 static std::string read(const std::filesystem::path&p){std::ifstream f(p,std::ios::binary);need(bool(f),"missing source file");std::array<char,32769>b{};f.read(b.data(),b.size());const auto n=f.gcount();need(n<=32768&&!f.bad(),"source read bound");return std::string(b.data(),std::size_t(n));}
 static double number(std::string_view v){while(!v.empty()&&v.front()==' ')v.remove_prefix(1);while(!v.empty()&&v.back()==' ')v.remove_suffix(1);double x=0;const auto r=std::from_chars(v.data(),v.data()+v.size(),x);need(r.ec==std::errc{}&&r.ptr==v.data()+v.size()&&std::isfinite(x),"nonfinite/invalid source operand");return x;}
 static Table parse(const std::string&raw,const SourcePin&pin){
  need(pfc::sha256_hex(raw.data(),raw.size())==pin.sha256,"source SHA differs");
  std::istringstream in(raw);std::vector<std::string>line;std::string s;while(std::getline(in,s)){if(!s.empty()&&s.back()=='\r')s.pop_back();line.push_back(s);}
  need(line.size()==191,"source line cardinality");
  need(line[0]=="RELATIVISTIC"&&line[1]=="ELASTIC SCATTERING CROSS SECTION"&&line[2]=="d sigma/d omega vs theta","source process/angle definition");
  need(line[4]=="Atomic number: "+std::to_string(pin.Z)&&line[5]=="Energy: "+std::to_string(pin.energy)+" eV","source target/energy");
  need(line[7]=="(a0**2 = 2.8002852E-21 m**2)"&&line[9]=="Theta (deg), d sigma/d omega (a0**2/sr)","source units");
  constexpr std::string_view prefix="Total cross section: ",suffix=" a0**2";std::string_view total=line[6];
  need(total.starts_with(prefix)&&total.ends_with(suffix),"reported total field");Table t;t.total=number(total.substr(prefix.size(),total.size()-prefix.size()-suffix.size()));need(t.total>0,"positive reported total");
  for(unsigned i=0;i<181;++i){const std::string_view row=line[i+10];const auto k=row.find(',');need(k!=row.npos&&row.find(',',k+1)==row.npos,"DCS columns");need(number(row.substr(0,k))==i,"complete primary angular grid");const double d=number(row.substr(k+1));need(d>0,"positive DCS");t.dcs[180-i]=d;}
  return t;
 }
 Shape shape(Target t,double E,Mode m)const {
  mode(m);const auto k=species(t);need(std::isfinite(E)&&E>=1000&&E<=10000,"unsupported incident energy");
  const unsigned lo=std::min(9u,unsigned((E-1000)/1000));const unsigned hi=std::min(9u,lo+1);const double a=(E-(1000+1000*lo))/1000;
  Shape s;s.reported=(1-a)*tables_[k][lo].total+a*tables_[k][hi].total;
  for(unsigned i=0;i<181;++i)s.dcs[i]=(1-a)*tables_[k][lo].dcs[i]+a*tables_[k][hi].dcs[i];
  for(unsigned i=0;i<180;++i){const double h=mu_[i+1]-mu_[i],x=s.dcs[i],y=s.dcs[i+1];s.cumulative[i+1]=s.cumulative[i]+.5*h*(x+y);
   s.transport+=h/6*((1-mu_[i])*x+4*(1-.5*(mu_[i]+mu_[i+1]))*.5*(x+y)+(1-mu_[i+1])*y);}
  s.integrated=s.cumulative.back();need(s.integrated>0&&std::isfinite(s.integrated),"invalid angular normalization");return s;
 }
 double inverse(const Shape&s,double q)const {
  need(std::isfinite(q)&&q>=0&&q<=1,"CDF quantile domain");if(q==0)return -1;if(q==1)return 1;
  const double y=q*s.integrated;auto it=std::upper_bound(s.cumulative.begin(),s.cumulative.end(),y);
  if(it==s.cumulative.end())return 1;
  const unsigned i=unsigned(it-s.cumulative.begin()-1);const double h=mu_[i+1]-mu_[i],target=y-s.cumulative[i];
  const double slope=(s.dcs[i+1]-s.dcs[i])/h;double lo=0,hi=h;
  for(unsigned k=0;k<64;++k){const double x=.5*(lo+hi);if(x==lo||x==hi)break;const double area=x*(s.dcs[i]+.5*slope*x);if(area<target)lo=x;else hi=x;}
  return mu_[i]+.5*(lo+hi);
 }
 Event event(Target t,double E,pfc::V3 d,pfc::V3 u,double internal,double mu,double az,Mode m)const {
  mode(m);const auto k=species(t);need(internal==0,"excited/internal target state unsupported");
  need(finite(d)&&std::abs(pfc::norm2(d)-1)<=64*std::numeric_limits<double>::epsilon(),"unit incident direction required");
  need(finite(u)&&std::isfinite(mu)&&mu>=-1&&mu<=1&&std::isfinite(az)&&az>=0&&az<1,"finite target/angular operands");
  Event r;r.E_tf=E;r.mu_tf=mu;r.sigma_m2=cross_section(t,E,m);r.model_sha256=identity_;
  const double me=electron_mass_,M=target_mass_[k],p=pfc_sqrt(E*(E+2*me)),Et=me+E,C=M+Et;
  // Exact relativistic elastic solution at a specified TARGET-REST LAB angle.
  // sqrt(A^2-me^2 D)=p sqrt(M^2-me^2(1-mu^2)); this avoids subtracting rest-energy squares.
  const double D=C*C-p*p*mu*mu;
  const double po=p*((M*Et+me*me)*mu+C*pfc_sqrt(M*M-me*me*(1-mu*mu)))/D;
  need(po>=0&&std::isfinite(po),"finite-recoil root");
  const auto pi=pfc::scale(d,p),pe=pfc::scale(pfc::deflect(d,mu,az),po),pt=pfc::sub(pi,pe);
  const auto w=pfc::w_from_u(u);const double gamma=pfc_sqrt(1+pfc::norm2(w));
  const auto il=pfc::boost_out(w,gamma,pi,me+E),el=pfc::boost_out(w,gamma,pe,pfc::energy_of(me,pe)),tl=pfc::boost_out(w,gamma,pt,pfc::energy_of(M,pt));
  r.incoming_u={pfc::u_from_p(me,il),u};r.outgoing_u={pfc::u_from_p(me,el),pfc::u_from_p(M,tl)};
  const std::array<double,2>mass{me,M};pfc::V3 P0{},P1{};double K0=0,K1=0;
  for(unsigned j=0;j<2;++j){r.incoming_ke_eV[j]=pfc::ke_from_u(mass[j],r.incoming_u[j]);r.outgoing_ke_eV[j]=pfc::ke_from_u(mass[j],r.outgoing_u[j]);K0+=r.incoming_ke_eV[j];K1+=r.outgoing_ke_eV[j];P0=pfc::add(P0,pfc::p_from_u(mass[j],r.incoming_u[j]));P1=pfc::add(P1,pfc::p_from_u(mass[j],r.outgoing_u[j]));}
  r.energy_residual_eV=K1-K0;r.momentum_residual_eV_c=pfc::sub(P1,P0);r.recoil_tf_eV=pfc::kinetic_from_p2(M,pfc::norm2(pt));
  const double eps=2048*std::numeric_limits<double>::epsilon();
  need(std::isfinite(K0)&&std::isfinite(K1)&&std::abs(r.energy_residual_eV)<=eps*std::max(1.,K0+K1),"phase-space energy closure");
  need(pfc_sqrt(pfc::norm2(r.momentum_residual_eV_c))<=eps*std::max(1.,pfc_sqrt(pfc::norm2(P0))+pfc_sqrt(pfc::norm2(P1))),"phase-space momentum closure");return r;
 }
public:
 Model(const Model&)=default;Model(Model&&)=default;Model&operator=(const Model&)=default;Model&operator=(Model&&)=default;
 static Model load(const std::filesystem::path&sources,const std::string&registry_path,const std::string&registry_sha){
  need(registry_sha=="8435e873b7c68146fb4c429d96c3e8051334c22525fd508729299e8a2d65ea0f"||registry_sha=="70fcc2adc2fdaa62abed7ed02d68ffc09b1b18554eda8fa99c5451120ae454b8","unqualified registry identity");
  const auto blob=pfc::load_blob(registry_path,registry_sha);Model m;m.registry_blob_sha_=registry_sha;
  auto sid=[&](const char*name){const auto it=std::find(blob.species_names.begin(),blob.species_names.end(),name);need(it!=blob.species_names.end(),"registry species missing");return std::size_t(it-blob.species_names.begin());};
  m.electron_mass_=blob.view.species[sid("e")].mass_eV;m.target_mass_={blob.view.species[sid("H")].mass_eV,blob.view.species[sid("O")].mass_eV};
  need(blob.view.species[sid("e")].charge==-1&&blob.view.species[sid("H")].charge==0&&blob.view.species[sid("O")].charge==0,"registry charge operands");
  for(unsigned i=0;i<181;++i){double c,s;pfc::sincos_2pi(double(180-i)/360,&c,&s);m.mu_[i]=c;}
  need(m.mu_.front()==-1&&m.mu_.back()==1,"full angular endpoints");for(unsigned i=0;i<180;++i)need(m.mu_[i]<m.mu_[i+1],"strict angular grid");
  pfc::Sha256 id;id.text(model_contract);id.text(registry_sha);id.pod(a0_squared_m2);id.pod(m.electron_mass_);for(double v:m.target_mass_)id.pod(v);
  for(const auto&pin:source_pins){auto t=parse(read(sources/std::string(pin.filename)),pin);const unsigned target=pin.Z==1?0:1;m.tables_[target][pin.energy/1000-1]=t;id.text(std::string(pin.filename));id.text(std::string(pin.sha256));}
  m.identity_=id.hex();return m;
 }
 const std::string&identity()const noexcept{return identity_;}
 const std::string&registry_blob_sha256()const noexcept{return registry_blob_sha_;}
 double electron_mass_eV()const noexcept{return electron_mass_;}
 double target_mass_eV(Target t)const{return target_mass_[species(t)];}
 double cross_section(Target t,double E,Mode m)const{return shape(t,E,m).reported*a0_squared_m2;}
 double quantile(Target t,double E,double q,Mode m)const{return inverse(shape(t,E,m),q);}
 double cdf(Target t,double E,double mu,Mode m)const {need(std::isfinite(mu)&&mu>=-1&&mu<=1,"CDF angle domain");const auto s=shape(t,E,m);if(mu==-1)return 0;if(mu==1)return 1;const auto it=std::upper_bound(mu_.begin(),mu_.end(),mu);const unsigned i=unsigned(it-mu_.begin()-1);const double x=mu-mu_[i],h=mu_[i+1]-mu_[i];return (s.cumulative[i]+x*(s.dcs[i]+.5*(s.dcs[i+1]-s.dcs[i])/h*x))/s.integrated;}
 std::vector<Summary>summaries()const {std::vector<Summary>v;for(Target t:{Target::H,Target::O})for(unsigned k=0;k<10;++k){const double E=1000+1000*k;const auto s=shape(t,E,Mode::diagnostic_reported_total_normalized_shape);const double factor=2*3.141592653589793238462643383279502884;v.push_back({t,E,s.reported,s.integrated*factor,s.transport*factor,s.integrated*factor/s.reported-1});}return v;}
 Event controlled_event(Target t,double E,pfc::V3 d,pfc::V3 u,double internal,double mu,double az,Mode m)const{return event(t,E,d,u,internal,mu,az,m);}
 Event sample(Target t,double E,pfc::V3 d,pfc::V3 u,double internal,const pfc::ChemKey&key,pfc_u32 candidate,Mode m)const{
  need(pfc::chem_key_valid(key,candidate)&&key.purpose==PFC_PURPOSE_FIXTURE_TEST,"diagnostic random key domain");static_assert(PFC_W_KIN==32&&PFC_W_KIN+1<PFC_W_KIN_END&&PFC_W_KIN_END<=PFC_W_H2_PARTITION);
  return event(t,E,d,u,internal,quantile(t,E,pfc::chem_u01(key,candidate,PFC_W_KIN),m),pfc::chem_u01(key,candidate,PFC_W_KIN+1),m);
 }
};
} // namespace atomic_elastic
