#include "thermoelastic_v3.h"
#include "ss304_heat_material_v2.hpp"
#include <fstream>
#include <iostream>
#include <numeric>

namespace {
int checks=0;
void test(bool ok,const std::string& name){++checks;if(!ok)throw std::runtime_error(name);}
void call(int code){if(code)throw std::runtime_error(tfv3_last_error());}
void near(double a,double b,double abs_tol,const std::string& name){test(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=abs_tol,name);}
template<class F>void rejects(F fn,const std::string& name){bool caught=false;try{fn();}catch(const std::exception&){caught=true;}test(caught,name);}
struct Case {
 tfv3_handle handle=nullptr;tfv3_config config{};std::vector<double>initial;
 std::string source=std::string(64,'1'),mesh="thermal-reference-grid-v3",electrode="upper",material=electrode_heat::ss304_v2_source_seal_sha256();
 explicit Case(int nx=2,int ny=2,double T=300.){
  config.abi_version=3;config.nx=nx;config.ny=ny;config.x0_m=-.0025;config.y0_m=.0025;config.width_m=.005;config.thickness_m=.0002;config.height_m=.025;
  config.reference_density_kg_m3=7897.403;config.initial_temperature_K=300.;config.initial_time_s=0;
  initial.assign((nx+1)*(ny+1),T);
 }
 void create(){config.mesh_revision=mesh.c_str();config.material_sha256=material.c_str();config.source_sha256=source.c_str();config.electrode_id=electrode.c_str();config.initial_nodal_temperature_K=initial.data();config.initial_node_count=initial.size();call(tfv3_create(&config,&handle));}
 ~Case(){tfv3_destroy(handle);}
 Case(const Case&)=delete;Case&operator=(const Case&)=delete;
 tfv3_state_info info(){tfv3_state_info i{};call(tfv3_get_state_info(handle,&i));return i;}
 std::vector<double>T(int trial=0){std::vector<double>v(initial.size());call(tfv3_copy_nodal_temperature(handle,trial,v.data(),v.size()));return v;}
 std::vector<double>energy(int trial=0){std::vector<double>v(config.nx*config.ny);call(tfv3_copy_cell_energy(handle,trial,v.data(),v.size()));return v;}
};
struct Request {
 tfv3_epoch epoch{};std::string state,id,em=std::string(64,'2');
 Request(Case& c,double dt){auto i=c.info();state=i.state_sha256;id=state+"/mechanical/"+std::to_string(i.revision+1);epoch.base_revision=i.revision;epoch.mechanical_sequence=i.revision+1;epoch.start_time_s=i.time_s;epoch.dt_s=dt;epoch.stage_time_s=i.time_s+dt;epoch.mesh_revision=c.mesh.c_str();epoch.material_sha256=c.material.c_str();epoch.source_sha256=c.source.c_str();epoch.em_state_sha256=em.c_str();epoch.base_state_sha256=state.c_str();epoch.transaction_id=id.c_str();}
 tfv3_budget prepare(Case& c,const std::vector<double>& energy){tfv3_budget b{};call(tfv3_prepare(c.handle,&epoch,energy.data(),energy.size(),&b));return b;}
 void discard(Case& c){call(tfv3_discard(c.handle,id.c_str()));}
 void commit(Case& c){call(tfv3_preflight(c.handle,id.c_str(),state.c_str()));call(tfv3_commit(c.handle,id.c_str(),state.c_str()));}
};

void material_contract(){
 using namespace electrode_heat;
 test(thermal_v3_detail::sha256("abc")=="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad","SHA256 abc vector");
 test(thermal_v3_detail::sha256("")=="e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855","SHA256 empty vector");
 auto m=conditional_ss304_v2(ss304v2::Qualification::nominal_1800,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,7897.403,ss304_v2_source_seal_sha256());
 test(m.energy_reference(298.2)<0,"negative cold-reference energy allowed");near(m.energy_reference(300.),0,0,"caloric reference300K");
 for(auto pair:std::vector<std::pair<double,double>>{{299.,300.},{400.,300.},{1750.,1600.},{1800.,298.2}}){double a=pair.first,b=pair.second,target=m.energy_reference(a)-m.energy_reference(b);near(m.delta_h(a,b),target,32*std::numeric_limits<double>::epsilon()*std::abs(target)+1e-7,"stable delta h versus owner endpoints");near(m.delta_h(a,b),-m.delta_h(b,a),0,"signed enthalpy antisymmetry");}
 double tiny=300.+1e-9;near(m.delta_h(tiny,300.)/(tiny-300.),m.cv(300.),1e-3,"tiny enthalpy increment capacity");
 tiny=std::nextafter(300.,INFINITY);near(m.delta_h(tiny,300.)/(tiny-300.),m.cv(300.),1e-7,"nextafter enthalpy increment capacity");
 for(double boundary:{ss304::Model::t_solidus,ss304::Model::t_liquidus}){double a=boundary-2,b=boundary+2,target=m.energy_reference(b)-m.energy_reference(a);near(m.delta_h(b,a),target,1e-5,"phase-boundary exact capacity integral");near(m.delta_h(b,a),-m.delta_h(a,b),0,"phase-boundary reverse signed integral");}
 auto half=conditional_ss304_v2(ss304v2::Qualification::nominal_1800,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,7897.403/2,ss304_v2_source_seal_sha256());near(half.delta_h(1750.,299.),m.delta_h(1750.,299.)/2,1e-6,"reference density applied once");
 rejects([&]{m.cv(298.199);},"cold support rejected");
 rejects([&]{conditional_ss304_v2(ss304v2::Qualification::source_conditioned_liquid_3000,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,7897.403,ss304_v2_source_seal_sha256());},"wrong material qualification rejected");
 rejects([&]{conditional_ss304_v2(ss304v2::Qualification::nominal_1800,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,7897.403,std::string(64,'0'));},"wrong material seal rejected");
}

void signed_cell_port(){
 Case c;c.create();const auto original=c.info();const auto initial=c.T();Request r(c,1e-9);
 std::vector<double>input{1e-5,-2e-6,3e-6,-1e-6};auto b=r.prepare(c,input);
 near(b.mechanical_exchange_J,std::accumulate(input.begin(),input.end(),0.),1e-20,"mechanical budget signed sum");
 test(b.joule_J==0&&b.authored_heat_J==0&&b.wall_net_J==0,"mechanical exchange is not external/Joule/wall heat");
 near(b.thermal_delta_J,b.mechanical_exchange_J,b.roundoff_bound_J,"native enthalpy balances signed cell load");
 test(c.info().revision==original.revision&&c.T()==initial,"prepare does not mutate accepted state");
 test(c.T(1)!=initial,"native thermal trial actually evolved");
 auto e=c.energy(1);near(std::accumulate(e.begin(),e.end(),0.),b.trial_reference_energy_J,1e-15,"percell enthalpy sums to native energy");
 std::vector<double>d(initial.size()),m(d.size()),q(d.size()),res(d.size());call(tfv3_copy_trial_nodal_ledger(c.handle,d.data(),m.data(),q.data(),res.data(),d.size()));
 near(std::accumulate(q.begin(),q.end(),0.),0.,1e-15,"insulated conduction globally redistributes only");
 near(std::accumulate(m.begin(),m.end(),0.),b.mechanical_exchange_J,1e-19,"cell to nodal mechanical conservation");
 for(std::size_t i=0;i<d.size();++i)near(d[i]-m[i]+q[i],res[i],1e-16,"native nodewise first law");
 // Bottom-left corner receives one quarter of cell0; central node receives all4 quarters.
 near(m[0],input[0]/4,1e-18,"mechanical Q1 corner load");near(m[4],b.mechanical_exchange_J/4,1e-18,"mechanical Q1 shared node load");
 std::vector<double>copied(input.size());call(tfv3_copy_trial_cell_mechanical(c.handle,copied.data(),copied.size()));test(copied==input,"original cell debit identity retained");
 r.discard(c);test(c.info().revision==0&&c.T()==initial,"discard preserves accepted state");
 r.prepare(c,input);r.commit(c);test(c.info().revision==1,"commit advances state once");near(c.info().cumulative_mechanical_exchange_J,b.mechanical_exchange_J,1e-20,"cumulative mechanical ledger");
 tfv3_budget scratch{};test(tfv3_prepare(c.handle,&r.epoch,input.data(),input.size(),&scratch)!=0,"committed replay rejected");
 Request cool(c,1e-9);std::vector<double>negative(input.size(),-1e-6);auto cb=cool.prepare(c,negative);test(cb.mechanical_exchange_J<0&&cb.thermal_delta_J<0,"negative reversible mechanical exchange cools");cool.commit(c);
}

void identities_and_failure_atomicity(){
 Case c;c.create();Request r(c,1e-9);std::vector<double>zeros(4);tfv3_budget b{};auto base=c.info();
 auto reject=[&](tfv3_epoch epoch,const std::string& name){test(tfv3_prepare(c.handle,&epoch,zeros.data(),zeros.size(),&b)!=0,name);test(c.info().revision==base.revision&&std::string(c.info().state_sha256)==base.state_sha256,"failed prepare leaves accepted identity");};
 auto e=r.epoch;e.stage_time_s=e.start_time_s;reject(e,"wrong thermal stage rejected");
 e=r.epoch;e.base_revision=3;reject(e,"wrong base revision rejected");
 e=r.epoch;e.mechanical_sequence=2;reject(e,"wrong deterministic sequence rejected");
 e=r.epoch;e.mesh_revision="other-mesh";reject(e,"wrong mesh rejected");
 e=r.epoch;e.transaction_id="invented-id";reject(e,"wrong deterministic transaction rejected");
 const std::string wrong(64,'9');e=r.epoch;e.base_state_sha256=wrong.c_str();reject(e,"wrong state digest rejected");
 e=r.epoch;e.material_sha256=wrong.c_str();reject(e,"wrong material rejected");
 e=r.epoch;e.source_sha256=wrong.c_str();reject(e,"wrong source rejected");
 zeros[0]=std::numeric_limits<double>::quiet_NaN();test(tfv3_prepare(c.handle,&r.epoch,zeros.data(),zeros.size(),&b)!=0,"nonfinite signed cell load rejected");zeros[0]=-1e4;test(tfv3_prepare(c.handle,&r.epoch,zeros.data(),zeros.size(),&b)!=0,"below material floor cooling rejected");test(c.info().revision==0,"invalid cooling never commits");
}

void conduction_and_restart(const std::string& output){
 Case c(4,2);for(int j=0;j<=2;++j)for(int i=0;i<=4;++i)c.initial[j*5+i]=300+20.*i;c.create();
 const auto before=c.T();const double E0=c.info().stored_reference_energy_J;Request r(c,1e-5);std::vector<double>zero(8);auto b=r.prepare(c,zero);
 auto after=c.T(1);test(after.front()>before.front()&&after.back()<before.back(),"actual native MFEM diffusion changes nonuniform state");test(b.newton_iterations>0,"native nonlinear solve ran");near(b.thermal_delta_J,0,b.roundoff_bound_J,"zero external/internal heat conduction conservation");r.commit(c);
 near(c.info().stored_reference_energy_J,E0,2*b.roundoff_bound_J,"accepted reference energy conserved by conduction");
 const auto path=output+"/native-v3.chk";call(tfv3_save(c.handle,path.c_str()));
 Case resumed(4,2);resumed.initial=c.initial;resumed.create();call(tfv3_load(resumed.handle,path.c_str()));test(resumed.T()==c.T(),"native checkpoint temperatures exact");test(std::string(resumed.info().state_sha256)==c.info().state_sha256,"native restart state digest exact");
 Request a(c,1e-6),z(resumed,1e-6);std::vector<double>mechanical(8,1e-6);a.prepare(c,mechanical);z.prepare(resumed,mechanical);test(c.T(1)==resumed.T(1),"restart trial deterministic");a.commit(c);z.commit(resumed);test(std::string(c.info().state_sha256)==resumed.info().state_sha256,"restart committed replay identity exact");
 Case wrong(4,2);wrong.initial=c.initial;wrong.source=std::string(64,'8');wrong.create();test(tfv3_load(wrong.handle,path.c_str())!=0,"restart wrong source rejected");
 std::ifstream f(path,std::ios::binary);std::string bytes((std::istreambuf_iterator<char>(f)),std::istreambuf_iterator<char>());bytes[bytes.size()/2]^=1;std::ofstream corrupt(output+"/corrupt.chk",std::ios::binary);corrupt.write(bytes.data(),bytes.size());corrupt.close();Case damaged(4,2);damaged.initial=c.initial;damaged.create();test(tfv3_load(damaged.handle,(output+"/corrupt.chk").c_str())!=0,"checkpoint checksum rejects corrupted valid-looking state");
}
}

int main(int argc,char** argv){try{const std::string output=argc>1?argv[1]:".";std::filesystem::create_directories(output);material_contract();signed_cell_port();identities_and_failure_atomicity();conduction_and_restart(output);std::cout<<"{\"status\":\"passed\",\"checks\":"<<checks<<",\"native_MFEM\":\"4.10.0\",\"internal_mechanical_port\":true,\"physical_production\":false}"<<std::endl;return 0;}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
