#include "surface_emission.hpp"
#include "magnetic_guard_assessment.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <sched.h>
#include <sys/resource.h>
namespace se=surface_emission;namespace mg=magnetized_emission;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
std::vector<std::string>split(const std::string&s){std::istringstream in(s);std::vector<std::string>r;std::string x;while(std::getline(in,x,','))r.push_back(x);return r;}
double dot(const se::Vector&a,const se::Vector&b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
int main(int argc,char**argv){try{
 cpu_set_t mask;CPU_ZERO(&mask);need(sched_getaffinity(0,sizeof(mask),&mask)==0&&CPU_COUNT(&mask)==1&&CPU_ISSET(7,&mask)&&getpriority(PRIO_PROCESS,0)>=15,"CPU7 nice15 required");need(argc==3,"input output prefix");
 std::ifstream in(argv[1]);need(bool(in),"input missing");std::string line;std::getline(in,line);std::ofstream out(std::string(argv[2])+".csv");out<<std::setprecision(17)<<"legacy_row,face,T_K,F_V_m,relative_window_half_width,uniform_SN_potential_error_bound_eV,bound_over_kT,maximum_real_turningpoint_deficit_eV,nominal_rejected_J_A_m2,nominal_rejected_Q_W_m2,deep_turningpoints_rejected_J,local_turningpoints_rejected_J,above_saddle_rejected_J,deep_turningpoints_rejected_Q,local_turningpoints_rejected_Q,above_saddle_rejected_Q\n";
 std::size_t checks=0,positive=0;int row=-1;double max_p=0,max_q=0;
 while(std::getline(in,line)){++row;auto r=split(line);se::Vector n{},E{},B{};for(int j=0;j<3;++j){n[j]=std::stod(r[4+j]);E[j]=std::stod(r[9+j]);B[j]=std::stod(r[12+j]);}double F=se::extraction_field(E,n),T=std::stod(r[18]);if(F<=0)continue;++positive;
 auto d=se::evaluate(T,F,se::Config{4.3});auto total=mg::assess_rejected_supply(d,B,n);const double c=se::qe/(16*3.14159265358979323846*se::eps0),xs=std::sqrt(c/F),Us=4.3-d.barrier_lowering_eV;
 for(double a:{.1,.25,.5}){const double deficit=F*xs*a*a/(1+a),bound=F*xs*a*a*a/(1-a);double fraction[3]{},J[3]{},Q[3]{};
 for(int bin=0;bin<3;++bin){auto part=d;part.nodes.clear();long double mass=0;for(const auto&node:d.nodes){const int b=node.normal_energy_eV<Us-deficit?0:(node.normal_energy_eV<Us?1:2);if(b==bin){part.nodes.push_back(node);mass+=node.probability;}}if(mass==0)continue;for(auto&node:part.nodes)node.probability/=double(mass);part.current_A_m2=d.current_A_m2*double(mass);part.log_current_A_m2=d.log_current_A_m2+std::log(double(mass));auto rejected=mg::assess_rejected_supply(part,B,n);fraction[bin]=double(mass)*rejected.rejected_fraction;J[bin]=rejected.rejected_current_A_m2;Q[bin]=rejected.rejected_complete_wall_heat_W_m2;}
 const double pe=std::abs(fraction[0]+fraction[1]+fraction[2]-total.rejected_fraction),qe=std::abs(Q[0]+Q[1]+Q[2]-total.rejected_complete_wall_heat_W_m2)/std::max(1e-300,std::abs(d.wall_heat_W_m2));max_p=std::max(max_p,pe);max_q=std::max(max_q,qe);need(pe<1e-10&&qe<1e-10,"nominal rejected current/energy region partition failed");checks+=2;
 out<<row<<','<<r[0]<<','<<T<<','<<F<<','<<a<<','<<bound<<','<<bound/(se::kb*T/se::qe)<<','<<deficit<<','<<total.rejected_current_A_m2<<','<<total.rejected_complete_wall_heat_W_m2;for(double x:J)out<<','<<x;for(double x:Q)out<<','<<x;out<<'\n';
 }
 }
 std::ofstream report(std::string(argv[2])+"-checks.json");report<<std::setprecision(17)<<"{\"passed\":true,\"checks\":"<<checks<<",\"extracting_probes\":"<<positive<<",\"relative_windows\":[0.1,0.25,0.5],\"maximum_current_fraction_partition_error\":"<<max_p<<",\"maximum_heat_partition_error_over_nominal_heat\":"<<max_q<<",\"physical_current_error_bound_qualified\":false,\"above_saddle_continuation_qualified\":false,\"GPU_executed\":false}\n";
 std::cout<<"rejected-supply window partition checks="<<checks<<" positive probes="<<positive<<'\n';
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
