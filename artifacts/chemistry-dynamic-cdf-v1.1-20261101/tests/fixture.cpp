#include "chem_host.hpp"
#include "component.hpp"
#include <iostream>
#include <iomanip>
using namespace pfc;
int main(int argc,char**argv) {
 require(argc==3,"blob/expected SHA required");auto b=load_blob(argv[1],argv[2]);const auto& t=b.view;
 if(t.blob_version!=2) {std::cout<<"legacy blob loaded\n";return 0;}
 require(t.n_cdf==1 && t.n_cdf_points==100,"descriptor counts");
 const auto& d=t.cdf_desc[0];const auto& g=t.groups[group_index(t,d.projectile,d.target)];const auto& pr=t.procs[g.proc_off];
 require(pr.kind==PFC_P_DISSOC_CDF && sizeof(ProcessRec)==512 && sizeof(CdfRec)==256,"explicit CDF process ABI");
 unsigned direct=0,low=0,scope=0,populated=0,selected_refusal=0;
 double worst=0;
 for(unsigned i=0;i<1024;++i) {
  const double E=i%4==0?7.2:i%4==1?10.:i%4==2?20.:150.;const double m=t.species[d.projectile].mass_eV,M=t.species[d.target].mass_eV;
  const V3 pp=v3(0,0,std::sqrt(E*(E+2*m))),ut=i%5==0?v3(12345,-6789,4321):v3(0,0,0),wt=w_from_u(ut);const double gt=std::sqrt(1+norm2(wt));
  const ChemKey k{0x1234abcdefULL,17,1+i,PFC_PURPOSE_FIXTURE_TEST};V3 out[PFC_MAX_OUT]{};double mk[PFC_MAX_OUT],internal;
  require(event_kinematics(t,pr,g,0,k,0,m,M,pp,E,wt,gt,out,mk,&internal),"supported CDF event refused");
  const auto reference=h2diag::sample(E,h2diag::EnergyFrame::target_rest,v3(0,0,1),ut,k,0,h2diag::Mode::diagnostic_primary_cdf_iso);
  require(std::memcmp(out,reference.outgoing,3*sizeof(V3))==0 && internal==reference.internal_eV,"descriptor path differs from qualified component");
  double ki=reference.kinetic_in[0]+reference.kinetic_in[1],ko=0;for(unsigned j=0;j<3;++j)ko+=kinetic_from_p2(mk[j],norm2(out[j]));
  const double defect=std::fabs(ki-ko-pr.q_eV)/(ki+ko+pr.q_eV);require(defect<=0x1p-40,"energy closure");worst=std::max(worst,defect);++direct;
 }
 const double m=t.species[d.projectile].mass_eV,M=t.species[d.target].mass_eV;
 for(double E:{4.5,5.,7.,d.event_threshold_eV-1e-8,151.,1e4}) {
  const ChemKey k{12,4,7,PFC_PURPOSE_FIXTURE_TEST};V3 out[PFC_MAX_OUT]{};double mk[PFC_MAX_OUT],internal;
  require(!event_kinematics(t,pr,g,0,k,0,m,M,v3(0,0,std::sqrt(E*(E+2*m))),E,v3(0,0,0),1,out,mk,&internal),"unsupported CDF energy admitted");
  unsigned region;require(coverage_class(pr,E,&region)==PFC_COVER_UNRESOLVED,"unsupported interval mislabeled known zero");++low;
 }
 // Actual MCC: every positive encounter retains missing-process support, even
 // zero-speed records. A selected unsupported CDF event must refuse, never null.
 for(double E:{0.,7.,20.,150.,1000.})for(unsigned id=1;id<=512;++id) {
  CarrierIn c{};c.id=id;c.species=d.projectile;c.weight=1.25;c.u[2]=PFC_C_M_S*std::sqrt(E*(E+2*m))/m;const auto saved=c;
  CellTargets cell{};cell.n_targets=1;cell.species[0]=d.target;cell.n[0]=(E==7.?1e30:1e22);
  StepConfig cfg{};cfg.seed=12345;cfg.step=91;cfg.dt_s=1e-10;cfg.trunc_sigmas=4;cfg.purpose=PFC_PURPOSE_FIXTURE_TEST;cfg.max_candidates=128;cfg.max_events=8;
  EventRec events[8]{};CarrierOut out{};mcc_carrier(t,b.cover,cfg,c,cell,events,&out);
  require(out.n_cover>0,"missing-channel support dropped");for(unsigned j=0;j<out.n_cover;++j)require(out.cover_bound[j]==0,"unbounded omission acquired finite bound");
  require(std::memcmp(&c,&saved,sizeof c)==0,"input mutation");
  if(out.status==PFC_ST_CDF_SUPPORT) {if(E==7.){require(out.n_events==0,"unsupported low-energy event emitted");++selected_refusal;}}
  else {require(out.status==PFC_ST_OK,"MCC execution guard");populated+=out.n_events;}
  EventRec ev2[8]{};CarrierOut out2{};mcc_carrier(t,b.cover,cfg,c,cell,ev2,&out2);
  require(std::memcmp(&out,&out2,sizeof out)==0 && std::memcmp(events,ev2,sizeof events)==0,"MCC replay");
  cell.n[0]=0;mcc_carrier(t,b.cover,cfg,c,cell,ev2,&out2);require(out2.status==PFC_ST_OK && out2.n_cover==0 && out2.n_events==0,"zero density encounter");++scope;
 }
 require(populated>0 && selected_refusal>100,"insufficient actual MCC population/refusals");
 std::cout<<std::setprecision(17)<<"{\"direct_CDF_events\":"<<direct<<",\"unsupported_energy_controls\":"<<low<<",\"positive_encounter_controls\":"<<scope<<",\"MCC_events\":"<<populated<<",\"selected_CDF_refusals\":"<<selected_refusal<<",\"maximum_energy_defect\":"<<worst<<",\"passed\":true,\"production_qualified\":false,\"network_complete\":false,\"native_adopted\":false}\n";
}
