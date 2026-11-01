#include "thermoelastic_v4.h"
#include "surface_source_binding.hpp"
#include "sha256.hpp"
#include <array>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sched.h>
#include <sys/resource.h>
namespace {
int checks=0;void check(bool b,const std::string& msg){++checks;if(!b)throw std::runtime_error(msg);}
void call(int code){check(code==0,tfv4_last_error());}
void close(double a,double b,double bound,const char* msg){check(std::abs(a-b)<=bound,msg);}
const std::string material="1a1845e04f1f6fda33338d9bd284f5e1ead9cbc4b7747d0c6355e32cb5661c2d",source(64,'a'),field(64,'b');
struct Fixture {
 tfv4_handle h=nullptr;std::vector<double>T=std::vector<double>(15,700);
 tfv4_patch patch{"oxide","H2O(g,0)","adsorbed-H2O",source.c_str(),1,0,1,1,0,.33,2.9915076246e-26};
 tfv4_config config{4,4,2,-.0025,.0025,.005,.0002,.025,8000,700,0,"test-native-v4",material.c_str(),source.c_str(),"upper",T.data(),T.size(),&patch,1,"e-","electrode-electrons"};
 Fixture(){call(tfv4_create(&config,&h));}~Fixture(){tfv4_destroy(h);}
 tfv4_state_info info(int w=0){tfv4_state_info i{};call(tfv4_get_state_info(h,w,&i));return i;}
};
struct Epoch {
 std::string base,tx;tfv4_epoch e{};
 Epoch(Fixture& f,double dt){auto i=f.info();base=i.state_sha256;tx=base+"/mechanical/"+std::to_string(i.revision+1);e.base_revision=i.revision;e.mechanical_sequence=i.revision+1;e.start_time_s=i.time_s;e.dt_s=dt;e.stage_time_s=i.time_s+dt;e.mesh_revision=f.config.mesh_revision;e.material_sha256=material.c_str();e.source_sha256=source.c_str();e.em_state_sha256=field.c_str();e.base_state_sha256=base.c_str();e.transaction_id=tx.c_str();e.common_state=i.revision+1;e.field_state=e.common_state;e.source_mesh_revision=1;e.source_material_revision=1;e.field_time_s=i.time_s;e.E_V_m[1]=-1e8;e.B_T[0]=30;e.material_phase_admitted=e.planar_barrier_admitted=e.electron_supply_admitted=1;}
};
tfv4_source_action begin(std::uint64_t id=77){tfv4_source_action a{};a.kind=1;a.layer_id=id;a.patch_index=0;a.patch_id=1;a.weight=1e9;a.normal_uniform=.6;a.transverse_uniform=.35;a.azimuth_uniform=.1;a.origin_m[0]=0;a.origin_m[1]=.0027;a.work_function_eV=4.3;a.energy_panels=128;a.action_nodes=64;a.tail_widths=48;a.local_linear_vacuum_extent_m=25e-6;a.image_force_fraction=.03;a.image_energy_fraction=.03;a.maximum_quantum_magnetic_energy_fraction=.3;a.maximum_gyroangle=.02;a.maximum_relative_normal_step=.02;a.maximum_substeps=4096;a.own_return_unit_accommodation=1;return a;}
tfv4_water_draw draw(){return {0,.4,.6,.3,.5,.5,.5};}
tfv4_budget prepare(Fixture& f,Epoch& e,const std::vector<tfv4_source_action>& actions={},double mechanical=1e-8){std::vector<double>cells(8,mechanical/8);auto d=draw();tfv4_budget b{};call(tfv4_prepare(f.h,&e.e,cells.data(),cells.size(),actions.data(),actions.size(),&d,1,&b));return b;}
void ledger(Fixture& f,const tfv4_budget& b){
 std::array<std::vector<double>,7> a;for(auto& v:a)v.resize(15);call(tfv4_copy_trial_nodal_ledger(f.h,a[0].data(),a[1].data(),a[2].data(),a[3].data(),a[4].data(),a[5].data(),a[6].data(),15));
 for(int i=0;i<15;++i)close(a[0][i]-a[1][i]-a[4][i]+a[2][i],a[3][i],1e-15,"local nodal law");
 close(std::accumulate(a[5].begin(),a[5].end(),0.),b.source_wall_J,1e-21,"electron load lost");close(std::accumulate(a[6].begin(),a[6].end(),0.),-b.neutral_K_J-b.neutral_internal_J-b.adsorbed_reference_change_J,1e-21,"water load lost");
 close(b.thermal_delta_J-b.mechanical_exchange_J-b.wall_net_J,b.residual_J,1e-21,"global external/mechanical partition");check(b.joule_J==0&&b.authored_heat_J==0&&b.gas_reference_J==0,"internal energy mislabeled");
 close(b.source_wall_J+b.source_particle_K_change_J+b.source_electron_reference_change_J+b.source_image_energy_change_J,b.source_field_work_J,1e-10*(std::abs(b.source_wall_J)+std::abs(b.source_particle_K_change_J)+std::abs(b.source_electron_reference_change_J)+std::abs(b.source_image_energy_change_J)+std::abs(b.source_field_work_J)),"certified source closure");
}
void water_and_identity(const std::filesystem::path& out){
 Fixture f;auto original=f.info();Epoch e(f,1e-12);auto b=prepare(f,e);ledger(f,b);check(b.neutral_K_J>0&&b.adsorbed_reference_change_J>0,"actual water inactive");check(f.info().revision==0,"prepare mutated accepted state");auto candidate=f.info(1);check(candidate.gas_packet_count==1&&std::strcmp(candidate.state_sha256,original.state_sha256)!=0,"candidate joint state absent");
 call(tfv4_discard(f.h,e.tx.c_str()));check(std::strcmp(f.info().state_sha256,original.state_sha256)==0,"discard mutated accepted");auto repeat=prepare(f,e);check(std::strcmp(f.info(1).state_sha256,candidate.state_sha256)==0,"repeat not deterministic");call(tfv4_preflight(f.h,e.tx.c_str(),e.base.c_str()));call(tfv4_commit(f.h,e.tx.c_str(),e.base.c_str()));check(std::strcmp(f.info().state_sha256,candidate.state_sha256)==0,"committed digest differs candidate");check(tfv4_commit(f.h,e.tx.c_str(),e.base.c_str())!=0,"committed replay accepted");
 auto path=(out/"water-joint.bin").string();call(tfv4_save(f.h,path.c_str()));Fixture r;call(tfv4_load(r.h,path.c_str()));check(std::strcmp(r.info().state_sha256,f.info().state_sha256)==0,"joint restart differs");
 std::ifstream saved(path,std::ios::binary);std::string corrupted((std::istreambuf_iterator<char>(saved)),{});corrupted.resize(corrupted.size()-64);double bad_K=1;std::memcpy(corrupted.data()+corrupted.size()-48,&bad_K,sizeof(bad_K));auto bad_path=(out/"rehashed-bad-gas.bin").string();std::ofstream bad_file(bad_path,std::ios::binary);bad_file<<corrupted<<thermal_v4_detail::sha256(corrupted);bad_file.close();Fixture rejected;auto rejected_before=rejected.info();check(tfv4_load(rejected.h,bad_path.c_str())!=0,"rehashed inconsistent gas K admitted");check(std::strcmp(rejected_before.state_sha256,rejected.info().state_sha256)==0,"bad restart changed accepted state");
 Epoch next(f,1e-12),next_r(r,1e-12);auto n=prepare(f,next,{},-1e-8),nr=prepare(r,next_r,{},-1e-8);check(n.mechanical_exchange_J<0&&nr.mechanical_exchange_J==n.mechanical_exchange_J,"signed cooling absent");check(std::strcmp(f.info(1).state_sha256,r.info(1).state_sha256)==0,"restarted continuation differs");call(tfv4_discard(f.h,next.tx.c_str()));call(tfv4_discard(r.h,next_r.tx.c_str()));
 // Overlap rejects before acquiring a handle; adjacent partitions remain legal.
 std::array<tfv4_patch,2> overlap{f.patch,f.patch};overlap[1].id="duplicate-water";auto config=f.config;config.patches=overlap.data();config.patch_count=2;tfv4_handle bad=nullptr;check(tfv4_create(&config,&bad)!=0&&!bad,"overlap duplicated physical reservoir");
}
void source_layers(const std::filesystem::path& out){
 Fixture f;Epoch e(f,1e-15);auto birth=begin();auto fly=birth;fly.kind=2;fly.duration_s=e.e.dt_s;auto b=prepare(f,e,{birth,fly});ledger(f,b);auto trial=f.info(1);check(trial.layer_count==1&&trial.path_count>1&&b.source_wall_J!=0,"actual sealed begin/advance inactive");std::array<tfv4_layer_info,1> layers{};call(tfv4_copy_layers(f.h,1,layers.data(),1));check(layers[0].status==0&&layers[0].time_s==e.e.stage_time_s,"short step did not retain active layer");
 call(tfv4_commit(f.h,e.tx.c_str(),e.base.c_str()));auto path=(out/"active-layer-joint.bin").string();call(tfv4_save(f.h,path.c_str()));Fixture restored;call(tfv4_load(restored.h,path.c_str()));check(std::strcmp(f.info().state_sha256,restored.info().state_sha256)==0,"active layer restart not exact");
 Epoch next(f,2e-12);fly.duration_s=next.e.dt_s;auto b2=prepare(f,next,{fly});ledger(f,b2);call(tfv4_copy_layers(f.h,1,layers.data(),1));check(b2.source_wall_J==0||layers[0].status==2,"continuation redebited birth heat");
 Epoch retry(restored,2e-12);auto b3=prepare(restored,retry,{fly});check(std::strcmp(f.info(1).state_sha256,restored.info(1).state_sha256)==0,"active-layer restart continuation differs");call(tfv4_commit(f.h,next.tx.c_str(),next.base.c_str()));call(tfv4_discard(restored.h,retry.tx.c_str()));
 // A separate zero-B handoff followed by a same-particle boundary collection
 // tests the exact collection API. It claims no intervening PIC trajectory.
 Fixture collected;Epoch born(collected,2e-12);born.e.B_T[0]=0;auto a=begin(91);auto move=a;move.kind=2;move.duration_s=born.e.dt_s;prepare(collected,born,{a,move});call(tfv4_copy_layers(collected.h,1,layers.data(),1));check(layers[0].status==1,"zero-B source failed handoff");auto K=layers[0].current_weighted_K_J;call(tfv4_commit(collected.h,born.tx.c_str(),born.base.c_str()));Epoch ce(collected,1e-12);a.kind=3;a.incoming_weighted_K_J=K;a.collection_accommodation_admitted=1;auto cb=prepare(collected,ce,{a});check(cb.source_wall_J>0&&cb.source_particle_K_change_J<0,"actual collection source heat inactive");ledger(collected,cb);call(tfv4_copy_layers(collected.h,1,layers.data(),1));check(layers[0].status==2&&layers[0].current_weighted_K_J==0,"same electron not removed on collection");
 // Skipping an active electron must reject the whole heat candidate.
 Epoch skipped(restored,1e-15);std::vector<double> zero(8);auto d=draw();tfv4_budget zz{};check(tfv4_prepare(restored.h,&skipped.e,zero.data(),8,nullptr,0,&d,1,&zz)!=0,"unadvanced active electron accepted");
}
}
int main(int argc,char** argv){try{check(argc==2,"report directory required");cpu_set_t mask;CPU_ZERO(&mask);sched_getaffinity(0,sizeof(mask),&mask);check(CPU_COUNT(&mask)==1&&CPU_ISSET(7,&mask)&&getpriority(PRIO_PROCESS,0)==15,"CPU7 nice15 required");std::filesystem::path out=argv[1];std::filesystem::create_directories(out);water_and_identity(out);source_layers(out);std::ofstream f(out/"native-v4-results.json");f<<"{\"passed\":true,\"checks\":"<<checks<<",\"component_tests_only\":true,\"sealed_source_calls_active\":true,\"water_transfers_active\":true,\"GPU_work\":false}\n";std::cout<<"passed "<<checks<<" checks\n";}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
