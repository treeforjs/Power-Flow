#include "chem_host.hpp"
#include <iostream>
#include <iomanip>
#include <array>
using namespace pfc;
struct Stats {unsigned n=0;long double x=0,y=0,xx=0,yy=0,xy=0;void add(double a,double b){++n;x+=a;y+=b;xx+=a*a;yy+=b*b;xy+=a*b;}double mean()const{return double(y/n);}double corr()const{return double((xy-x*y/n)/std::sqrt((xx-x*x/n)*(yy-y*y/n)));}};
int main(){
 unsigned interval_controls=0;
#ifdef PFC_RNG_LAYOUT_VERSION
 static_assert(PFC_W_TVEL_END<=PFC_W_KIN && PFC_W_KIN_END<=PFC_W_H2_PARTITION && PFC_W_H2_PARTITION_END<=PFC_EVENT_WORDS);
 require(kinematic_words_fit(PFC_W_KIN,PFC_W_KIN_END-PFC_W_KIN),"complete ordinary word span");
 require(!kinematic_words_fit(PFC_W_KIN-1,1) && !kinematic_words_fit(PFC_W_KIN_END,1) && !kinematic_words_fit(PFC_W_KIN,0xffffffffu),"kinematic range refusal");interval_controls=4;
 for(unsigned a=0;a<PFC_TVEL_ATTEMPTS;++a)for(unsigned j=0;j<PFC_TVEL_WORDS_PER_ATTEMPT;++j){const unsigned w=PFC_W_TVEL+a*PFC_TVEL_WORDS_PER_ATTEMPT+j;require(w<PFC_W_KIN,"target/ordinary alias");++interval_controls;}
 for(unsigned a=0;a<PFC_H2_PARTITION_ATTEMPTS;++a)for(unsigned j=0;j<2;++j){const unsigned w=PFC_W_H2_PARTITION+2*a+j;require(w>=PFC_W_KIN_END && w<PFC_EVENT_WORDS,"H2 word alias/domain");++interval_controls;}
#else
 require(PFC_W_TVEL+4==PFC_W_KIN,"baseline overlap not reproduced");
#endif
 SpeciesRec species[2]{};species[0].mass_eV=510998.95069;species[0].charge=-1;species[0].kind=PFC_KIND_ELECTRON;species[1].mass_eV=1.677e10;species[1].kind=PFC_KIND_NEUTRAL;
 ProcessRec pr{};pr.kind=PFC_P_ELASTIC;pr.angular=PFC_ANG_ISO;pr.tag=4;pr.n_out=2;pr.out_species[0]=0;pr.out_species[1]=1;pr.out_dest[0]=PFC_DEST_SURVIVOR;pr.out_dest[1]=PFC_DEST_GAS;pr.n_cover=1;pr.cover_lo[0]=0;pr.cover_hi[0]=1e300;pr.cover_class[0]=PFC_COVER_MODEL;pr.cover_region[0]=0xffffffffu;
 GroupRec gr{};gr.projectile=0;gr.target=1;gr.grid_n=2;gr.proc_n=1;gr.rmq_levels=1;gr.xi_off=0xffffffffu;
 double grid[2]={90,110},sigma[2]={1e-20,1e-20},rmq[1]={1e-20*PFC_C_M_S*std::sqrt(110*(110+2*species[0].mass_eV))/(110+species[0].mass_eV)*(1+0x1p-20)};
 pfc_u32 go[4]={0xffffffffu,0,0xffffffffu,0xffffffffu},bins[3]={0xffffffffu,0xffffffffu,0xffffffffu},off[1]={0};
 TableView table{};table.n_species=2;table.n_groups=1;table.n_proc=1;table.species=species;table.groups=&gr;table.procs=&pr;table.grid=grid;table.sigma=sigma;table.rmq=rmq;table.group_of=go;
 CoverView cover{bins,bins,off};CellTargets cell{};cell.n_targets=1;cell.species[0]=1;cell.n[0]=1e22;cell.s[0]=1e5;
 StepConfig cfg{};cfg.seed=0x123456789abcd;cfg.step=9;cfg.dt_s=5e-10;cfg.trunc_sigmas=1;cfg.purpose=PFC_PURPOSE_FIXTURE_TEST;cfg.max_candidates=64;cfg.max_events=16;
 std::array<Stats,5> stats{};Sha256 first_identity;unsigned populated=0,replays=0;double max_angle_error=0,max_energy_error=0;
 for(unsigned id=1;id<=40000;++id){
  CarrierIn c{};c.id=id;c.species=0;c.weight=1;c.u[2]=PFC_C_M_S*std::sqrt(100*(100+2*species[0].mass_eV))/species[0].mass_eV;
  EventRec ev[16]{};CarrierOut out{};mcc_carrier(table,cover,cfg,c,cell,ev,&out);require(out.status==PFC_ST_OK,"manufactured MCC guard failure");if(!out.n_events)continue;++populated;
  const auto& e=ev[0];const ChemKey key{cfg.seed,cfg.step,c.id,cfg.purpose};int accepted=4;V3 expected=v3(0,0,0);
  for(unsigned a=0;a<4;++a){const unsigned w=PFC_W_TVEL+4*a;const double r1=std::sqrt(-2*pfc::log(chem_u01(key,e.candidate,w))),r2=std::sqrt(-2*pfc::log(chem_u01(key,e.candidate,w+2)));double c1,s1,c2,s2;sincos_2pi(chem_u01(key,e.candidate,w+1),&c1,&s1);sincos_2pi(chem_u01(key,e.candidate,w+3),&c2,&s2);const V3 z=v3(r1*c1,r1*s1,r2*c2);if(norm2(z)<=1){accepted=int(a);expected=scale(z,cell.s[0]);break;}}
  require(e.ut[0]==expected.x && e.ut[1]==expected.y && e.ut[2]==expected.z,"actual target rejection replay differs");
  first_identity.pod(e.carrier);first_identity.pod(e.candidate);for(double u:e.ut)first_identity.pod(u);
  const double m=species[0].mass_eV,M=species[1].mass_eV;const V3 wt=w_from_u(expected);const double gt=std::sqrt(1+norm2(wt));
  const V3 pin=p_from_u(m,v3(e.u_in[0],e.u_in[1],e.u_in[2])),ptf=boost_in(wt,gt,pin,m+e.ke_in_p);const double W=m+M+kcm_from_T(m,M,e.E_tf);const V3 wc=scale(ptf,1/W);const double gc=std::sqrt(1+norm2(wc));
  const V3 po=p_from_u(m,v3(e.out_u[0][0],e.out_u[0][1],e.out_u[0][2]));const V3 pot=boost_in(wt,gt,po,m+e.ke_out[0]);const V3 poc=boost_in(wc,gc,pot,m+kinetic_from_p2(m,norm2(pot)));const double mu=poc.z/std::sqrt(norm2(poc));
  const double expected_mu=2*chem_u01(key,e.candidate,PFC_W_KIN)-1;max_angle_error=std::max(max_angle_error,std::fabs(mu-expected_mu));require(std::fabs(mu-expected_mu)<1e-12,"actual CM angle does not follow reserved word");
  const double scaleE=e.ke_in_p+e.ke_in_t+e.ke_out[0]+e.ke_out[1];const double de=std::fabs(e.ke_in_p+e.ke_in_t-e.ke_out[0]-e.ke_out[1])/scaleE;max_energy_error=std::max(max_energy_error,de);require(de<0x1p-40,"elastic energy closure");
  stats[accepted].add(chem_u01(key,e.candidate,PFC_W_TVEL+4),mu);
  if(accepted>=1 && accepted<4 && replays<32){EventRec again[16]{};CarrierOut second{};mcc_carrier(table,cover,cfg,c,cell,again,&second);require(std::memcmp(&out,&second,sizeof out)==0 && std::memcmp(ev,again,sizeof ev)==0,"late-rejection actual replay differs");++replays;}
 }
 require(stats[1].n>1000 && replays==32,"insufficient late-rejection controls");
#ifdef PFC_RNG_LAYOUT_VERSION
 require(std::fabs(stats[1].mean())<0.06 && std::fabs(stats[1].corr())<0.06,"corrected angular/target independence sanity failure");
 for(unsigned a=0;a<5;++a){require(stats[a].n>500 && std::fabs(stats[a].mean())<0.08,"conditional angular mean sanity");}
 const ChemKey key{1,1,1,PFC_PURPOSE_FIXTURE_TEST};pfc_u32 wp=PFC_W_KIN_END-1;double mass[2]={1,1};V3 pout[2]{};require(!chain_decay(key,0,&wp,3,v3(0,0,0),mass,2,1,pout) && wp==PFC_W_KIN_END-1,"chain word-budget refusal");++interval_controls;
#else
 require(stats[1].mean()>0.4 && stats[1].corr()>0.999999,"old overlap bias not reproduced");
#endif
 std::cout<<std::setprecision(17)<<"{\"ordinary_first_word\":"<<PFC_W_KIN<<",\"interval_controls\":"<<interval_controls<<",\"populated_actual_MCC\":"<<populated<<",\"late_rejection_replays\":"<<replays<<",\"first_event_target_identity\":\""<<first_identity.hex()<<"\",\"max_angle_reconstruction_error\":"<<max_angle_error<<",\"max_relative_energy_defect\":"<<max_energy_error<<",\"conditional\":[";
 for(unsigned a=0;a<5;++a){if(a)std::cout<<',';std::cout<<"{\"attempt\":"<<a<<",\"count\":"<<stats[a].n<<",\"mean_CM_cosine\":"<<stats[a].mean()<<",\"correlation_with_attempt1_first_word\":"<<stats[a].corr()<<'}';}
 std::cout<<"],\"passed\":true,\"physical_model_validated\":false,\"production_qualified\":false}\n";
}
