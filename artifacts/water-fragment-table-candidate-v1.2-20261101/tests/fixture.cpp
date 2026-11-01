#include "chem_host.hpp"
#include <iostream>
#include <iomanip>
#include <limits>
using namespace pfc;
unsigned species(const Blob& b,const std::string& n){for(unsigned i=0;i<b.species_names.size();++i)if(b.species_names[i]==n)return i;throw Error("species missing");}
int main(int argc,char**argv){
 require(argc==3,"blob and SHA required");const auto b=load_blob(argv[1],argv[2]);const auto& t=b.view;
 const unsigned electron=species(b,"e"),h=species(b,"H"),o=species(b,"O"),oh=species(b,"OH"),h2=species(b,"H2");
 require(b.n_species==22,"species ABI changed");unsigned controls=0;
 require(t.group_of[electron*t.n_species+h2]==0xffffffffu,"H2 rate-only bank became a native encounter group");++controls;
 for(unsigned target:{oh}){
  const auto gidx=t.group_of[electron*t.n_species+target];require(gidx!=0xffffffffu,"new encounter group absent");
  for(double E:{0.,4.48,4.5,7.3605,8.5,10.,150.,151.,10000.,1e8}){
   const auto u=u_from_p(t.species[electron].mass_eV,v3(0,0,std::sqrt(E*(E+2*t.species[electron].mass_eV))));
   CarrierIn c{};c.id=42;c.species=electron;c.cell=0;c.weight=1.25;c.u[0]=u.x;c.u[1]=u.y;c.u[2]=u.z;const auto saved=c;
   StepConfig cfg{};cfg.seed=123;cfg.step=4;cfg.dt_s=1e-16;cfg.trunc_sigmas=4;cfg.purpose=PFC_PURPOSE_MCC;cfg.max_candidates=32;cfg.max_events=8;
   CellTargets cell{};cell.n_targets=1;cell.species[0]=target;cell.n[0]=1e20;
   EventRec events[8]{};CarrierOut out{};mcc_carrier(t,b.cover,cfg,c,cell,events,&out);
   require(out.status==PFC_ST_OK && out.n_cover>0,"positive encounter lacked missing-channel support");
   require(std::memcmp(&c,&saved,sizeof c)==0,"input carrier mutated");
   for(unsigned i=0;i<out.n_cover;++i){require(out.cover_region[i]<b.n_regions,"bad region");require(out.cover_bound[i]==0,"missing region received finite bound");}
   if(E==0)for(unsigned i=0;i<out.n_cover;++i)require(out.cover_exposure[i]==0,"cold rest exposure is not zero");
   cell.n[0]=0;mcc_carrier(t,b.cover,cfg,c,cell,events,&out);require(out.status==PFC_ST_OK && out.n_cover==0 && out.n_events==0,"absent target generated encounter");
   controls+=2;
  }
 }
 const auto gi=t.group_of[electron*t.n_species+oh];const auto& g=t.groups[gi];const auto& pr=t.procs[g.proc_off];
 require(pr.kind==PFC_P_DISSOC && pr.n_out==3 && pr.out_species[0]==electron && pr.out_species[1]==h && pr.out_species[2]==o,"OH process mapping");
 require(pr.loss_eV==7.36 && pr.kinetic_fraction==1 && pr.internal_sink==PFC_SINK_NONE,"OH model changed");
 for(double E:{std::nextafter(7.3605,0.),std::nextafter(10.,11.),1e4}){const int at=grid_locate(t.grid+g.grid_off,g.grid_n,E);require(sigma_eval(t,g,0,at,E)==0,"OH data extrapolation");++controls;}
 for(double E:{7.3605,8.5,10.}){const int at=grid_locate(t.grid+g.grid_off,g.grid_n,E);require(sigma_eval(t,g,0,at,E)>0,"OH supported knot missing");++controls;}
 double maximum=0;unsigned n=0;
 for(unsigned i=0;i<256;++i){
  const double E=(i%3==0?7.3605:i%3==1?8.5:10.);const double m=t.species[electron].mass_eV,M=t.species[oh].mass_eV;
  const V3 pp=v3(0,0,std::sqrt(E*(E+2*m))),ut=(i%2?v3(12345,-6789,4321):v3(0,0,0)),wt=w_from_u(ut);const double gt=std::sqrt(1+norm2(wt));
  const ChemKey key{12345,8,i+1,PFC_PURPOSE_FIXTURE_TEST};V3 out[PFC_MAX_OUT];double masses[PFC_MAX_OUT],internal;
  require(event_kinematics(t,pr,g,grid_locate(t.grid+g.grid_off,g.grid_n,E),key,0,m,M,pp,E,wt,gt,out,masses,&internal),"OH event refused");
  const V3 pin=boost_out(wt,gt,pp,m+E),pt=p_from_u(M,ut);double ki=kinetic_from_p2(m,norm2(pin))+ke_from_u(M,ut),ko=0;V3 po=v3(0,0,0);
  for(unsigned j=0;j<3;++j){require(masses[j]==t.species[pr.out_species[j]].mass_eV,"unspecified product internal energy");ko+=kinetic_from_p2(masses[j],norm2(out[j]));po=add(po,out[j]);}
  const double scale=ki+ko+std::fabs(pr.q_eV);const double defect=std::fabs(ki-ko-pr.q_eV)/scale;maximum=std::max(maximum,defect);
  require(internal==0 && defect<=0x1p-40,"OH energy closure");const auto dp=sub(add(pin,pt),po);require(std::sqrt(norm2(dp))<=0x1p-40*(std::sqrt(norm2(pin))+std::sqrt(norm2(pt))+std::sqrt(norm2(po))),"OH momentum closure");++n;
 }
 std::cout<<std::setprecision(17)<<"{\"scope_controls\":"<<controls<<",\"OH_actual_blob_events\":"<<n<<",\"maximum_relative_energy_defect\":"<<maximum<<",\"H2_native_event_enabled\":false,\"native_transaction_tested\":false,\"production_qualified\":false,\"network_complete\":false,\"passed\":true}\n";
}
