#include "phase_guard.hpp"
#include "magnetic_guard_assessment.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sched.h>
#include <sys/resource.h>
namespace se=surface_emission;namespace vg=original_v1_guard;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
std::vector<std::string> split(std::string s){std::istringstream i(s);std::string x;std::vector<std::string>v;while(std::getline(i,x,','))v.push_back(x);return v;}
int main(int argc,char**argv){try{
 cpu_set_t mask;CPU_ZERO(&mask);need(sched_getaffinity(0,sizeof(mask),&mask)==0&&CPU_COUNT(&mask)==1&&CPU_ISSET(7,&mask)&&getpriority(PRIO_PROCESS,0)>=15,"CPU7 nice15 required");
 need(argc==6,"INPUT OUTPREFIX FIRST COUNT ENERGY_PANELS");std::ifstream in(argv[1]);need(bool(in),"input missing");std::string line;std::getline(in,line);std::vector<std::vector<std::string>> rows;while(std::getline(in,line))rows.push_back(split(line));int first=std::stoi(argv[3]),count=std::stoi(argv[4]),panels=std::stoi(argv[5]);need(first>=0&&count>0&&first+count<=int(rows.size()),"slice invalid");
 std::string out=argv[2];std::ofstream csv(out+"-map.csv"),phase(out+"-phases.csv"),report(out+"-checks.json");csv<<std::setprecision(17);phase<<std::setprecision(17);
 csv<<"legacy_row,face,historical_T_K,F_V_m,B_normal_T,B_tangent_T,nominal_J_A_m2,log_nominal_J,nominal_Q_W_m2,SN_guard_applies,nominal_admitted_fraction,nominal_rejected_fraction,nominal_admitted_J_A_m2,log_nominal_admitted_J,nominal_rejected_J_A_m2,log_nominal_rejected_J,nominal_admitted_Q_W_m2,nominal_rejected_Q_W_m2,minimum_cutoff_eV,maximum_cutoff_eV,normal_nodes_with_admissible_transverse_phase,retained_normal_nodes,any_nonzero_rejected_tail,maximum_common_dt_gyro_005_s\n";
 phase<<"legacy_row,normal_uniform,normal_energy_eV,quantum_overlap_m,source_scale_eV,magnetic_room_eV,maximum_transverse_energy_eV,conditional_admitted_fraction,transverse_uniform,transverse_energy_eV,computed_magnetic_fraction,native_begin_pass,azimuth_uniform\n";
 std::size_t checks=0,native_pass=0,native_reject=0,boundary_checks=0;double max_partition=0,max_ratio=0;
 for(int k=first;k<first+count;++k){const auto&r=rows[k];se::Vector n{},E{},B{};for(int j=0;j<3;++j){n[j]=std::stod(r[4+j]);E[j]=std::stod(r[9+j]);B[j]=std::stod(r[12+j]);}const double T=std::stod(r[18]),F=se::extraction_field(E,n),Bn=vg::dot(B,n),Bnorm=std::sqrt(vg::dot(B,B)),Bt=std::sqrt(std::max(0.,vg::dot(B,B)-Bn*Bn));
 se::PreparedFlux f;auto&s=f.stage;s.accepted_state=s.field_accepted_state=s.thermal_accepted_state=s.mesh_revision=s.material_revision=s.patch_id=1;s.T_K=T;s.total_E=E;s.total_B_T=B;s.outward_normal=n;s.area_m2=1;s.dt=1e-15;f.spectrum=se::evaluate(T,F,se::Config{4.3,panels,panels==128?64:96,panels==128?48.:64.});auto&d=f.spectrum;
 const auto rej=magnetized_emission::assess_rejected_supply(d,B,n);long double accepted=0;double cmin=INFINITY,cmax=0;std::size_t positive=0;
 for(const auto&node:d.nodes){auto p=vg::boundary(d,B,n,node.normal_energy_eV);double a=vg::admitted_fraction(p,node.normal_energy_eV,T);accepted+=node.probability*a;if(a>0)++positive;if(!p.thermal_step){cmin=std::min(cmin,p.maximum_transverse_eV);cmax=std::max(cmax,p.maximum_transverse_eV);}}
 const double a=std::clamp(double(accepted),0.,1.),q=rej.SchottkyNordheim_guard_applies?rej.rejected_fraction:0;max_partition=std::max(max_partition,std::abs(a+q-1));need(std::abs(a+q-1)<1e-10,"admitted/rejected partition mismatch");++checks;
 csv<<k<<','<<r[0]<<','<<T<<','<<F<<','<<Bn<<','<<Bt<<','<<d.current_A_m2<<','<<d.log_current_A_m2<<','<<d.wall_heat_W_m2<<','<<rej.SchottkyNordheim_guard_applies<<','<<a<<','<<q<<','<<d.current_A_m2*a<<','<<(a>0?d.log_current_A_m2+std::log(a):-INFINITY)<<','<<rej.rejected_current_A_m2<<','<<(rej.has_log_rejected_current?rej.log_rejected_current_A_m2:-INFINITY)<<','<<d.wall_heat_W_m2-rej.rejected_complete_wall_heat_W_m2<<','<<rej.rejected_complete_wall_heat_W_m2<<','<<cmin<<','<<cmax<<','<<positive<<','<<d.nodes.size()<<','<<(F>0&&rej.has_log_rejected_fraction)<<','<<.05*se::me/(se::qe*Bnorm)<<'\n';
 se::LayerPolicy policy;policy.local_linear_vacuum_extent_m=.01;
 for(double un:{.001,.01,.1,.5,.9,.99,.999}){const auto it=std::lower_bound(d.nodes.begin(),d.nodes.end(),un,[](const se::Node&x,double u){return x.cdf<u;});need(it!=d.nodes.end(),"normal quantile missing");auto p=vg::boundary(d,B,n,it->normal_energy_eV);const double ac=vg::admitted_fraction(p,it->normal_energy_eV,T);std::vector<double>uts{.001,.1,.5,.9,.999};if(ac>1e-8&&ac<1-1e-8){uts.push_back(ac*(1-1e-5));uts.push_back(ac+(1-ac)*1e-5);boundary_checks+=2;}
 for(double ut:uts){double Et=se::transverse_quantile_eV(it->normal_energy_eV,T,ut),ratio=vg::ratio(p,Et);for(double az:{.13,.67}){bool pass=true;double actual=0;try{auto v=se::begin_layer(f,policy,1,un,ut,az,1);actual=v.state.quantum_magnetic_energy_fraction;}catch(const std::invalid_argument&e){need(std::string(e.what()).find("magnetic energy criterion")!=std::string::npos,"native fixture failed outside magnetic predicate");pass=false;}need(pass==(ratio<=.2),"phase formula/native v1 disagreement");++checks;if(pass){++native_pass;max_ratio=std::max(max_ratio,std::abs(actual-ratio));}else ++native_reject;phase<<k<<','<<un<<','<<it->normal_energy_eV<<','<<p.x_m<<','<<p.scale_eV<<','<<p.room_eV<<','<<p.maximum_transverse_eV<<','<<ac<<','<<ut<<','<<Et<<','<<ratio<<','<<pass<<','<<az<<'\n';}}
 }
 }
 report<<std::setprecision(17)<<"{\"passed\":true,\"first\":"<<first<<",\"count\":"<<count<<",\"energy_panels\":"<<panels<<",\"checks\":"<<checks<<",\"native_predicate_pass\":"<<native_pass<<",\"native_predicate_reject\":"<<native_reject<<",\"near_boundary_draw_pairs\":"<<boundary_checks/2<<",\"maximum_supply_partition_error\":"<<max_partition<<",\"maximum_magnetic_ratio_absolute_difference\":"<<max_ratio<<",\"native_production_draws\":false,\"finite_offset_diagnostic_only\":true,\"fixture_geometry_extent_not_a_physical_admission\":true,\"GPU_executed\":false}\n";std::cout<<"slice "<<first<<" count "<<count<<" checks "<<checks<<" pass "<<native_pass<<" reject "<<native_reject<<'\n';
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
