#include "magnetic_guard_assessment.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sched.h>
#include <sys/resource.h>

namespace se=surface_emission;namespace mg=magnetized_emission;
namespace {
void need(bool b,const char* s){if(!b)throw std::runtime_error(s);}
double dot(const se::Vector& a,const se::Vector& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
se::PreparedFlux flux(double T,double F,const se::Vector& B) {
    se::Stage s;s.accepted_state=s.field_accepted_state=s.thermal_accepted_state=s.mesh_revision=s.material_revision=s.patch_id=1;
    s.T_K=T;s.total_E={0,-F,0};s.total_B_T=B;s.outward_normal={0,1,0};s.area_m2=1.25e-4;s.dt=1e-14;
    s.material_temperature_phase_admitted=s.locally_planar_barrier_admitted=s.local_equilibrium_electron_supply_admitted=true;
    return se::prepare(s,se::Config{4.3});
}
std::vector<std::string> split(const std::string& line){std::istringstream s(line);std::vector<std::string> r;std::string a;while(std::getline(s,a,','))r.push_back(a);return r;}
void qualify(const std::filesystem::path& out) {
    std::filesystem::create_directories(out);int checks=0;auto test=[&](bool b,const char* text){need(b,text);++checks;};
    const auto f=flux(1490,97724679.40458213,{30,0,0});
    const auto zero=mg::assess_rejected_supply(f.spectrum,{0,0,0},f.stage.outward_normal);
    test(zero.rejected_fraction==0 && zero.rejected_energy_moment_eV==0,"B0 guard rejected nonzero supply");
    const auto all=mg::assess_rejected_supply(f.spectrum,{0,10000,0},f.stage.outward_normal);
    test(std::abs(all.rejected_fraction-1)<1e-12,"Strong-normal-B predicate must reject all supply");
    test(std::abs(all.rejected_energy_moment_eV-f.spectrum.mean_metal_energy_eV)<1e-10,"Rejected Nottingham moment failed full-spectrum limit");
    const auto ret=mg::assess_rejected_supply(flux(1490,-1e8,{30,0,0}).spectrum,{30,0,0},f.stage.outward_normal);
    test(!ret.SchottkyNordheim_guard_applies,"Retarding thermal-step branch was relabeled as a validated magnetized law");
    double max_probability_error=0,max_moment_error=0;
    std::ofstream table(out/"analytic-guard-vs-native-predicate.csv");table<<std::setprecision(17)<<"T_K,F_V_m,B_tangent_T,analytic_rejected_fraction,native_grid_rejected_fraction,analytic_rejected_energy_moment_eV,native_grid_rejected_energy_moment_eV\n";
    for(const auto& input:std::vector<std::array<double,3>>{{1490,97724679.40458213,30},{773.15,1e7,30},{300,97724679.40458213,10},{0,1e9,30}}) {
        auto f=flux(input[0],input[1],{input[2],0,0});const auto a=mg::assess_rejected_supply(f.spectrum,f.stage.total_B_T,f.stage.outward_normal);
        int rejected=0;long double moment=0;constexpr int Nn=256,Nt=512;
        se::LayerPolicy p;p.local_linear_vacuum_extent_m=.001;
        for(int j=0;j<Nn;++j) {
            const double un=(j+.5)/Nn;const auto node=std::lower_bound(f.spectrum.nodes.begin(),f.spectrum.nodes.end(),un,[](const se::Node& n,double x){return n.cdf<x;});need(node!=f.spectrum.nodes.end(),"Bad test normal quantile");
            for(int k=0;k<Nt;++k) {
                const double ut=(k+.5)/Nt;bool pass=true;
                try{se::begin_layer(f,p,1,un,ut,.1,1);}catch(const std::invalid_argument& e){need(std::string(e.what()).find("magnetic energy criterion")!=std::string::npos,"Native predicate comparison hit a different rejection");pass=false;}
                if(!pass){++rejected;moment+=node->normal_energy_eV+se::transverse_quantile_eV(node->normal_energy_eV,input[0],ut);}
            }
        }
        const double numerical=double(rejected)/(Nn*Nt),numerical_moment=double(moment)/(Nn*Nt);
        max_probability_error=std::max(max_probability_error,std::abs(a.rejected_fraction-numerical));max_moment_error=std::max(max_moment_error,std::abs(a.rejected_energy_moment_eV-numerical_moment));
        test(std::abs(a.rejected_fraction-numerical)<.006,"Analytic conditional rejected rate disagrees with native predicate grid");
        test(std::abs(a.rejected_energy_moment_eV-numerical_moment)<.03,"Rejected phase energy moment disagrees with same-spectrum native draws");
        table<<input[0]<<','<<input[1]<<','<<input[2]<<','<<a.rejected_fraction<<','<<numerical<<','<<a.rejected_energy_moment_eV<<','<<numerical_moment<<'\n';
    }
    double maximum_refinement=0;
    for(double T:{300.,773.15,1490.})for(double F:{213175.33729099843,1e7,97724679.40458213}) {
        auto base=se::evaluate(T,F,se::Config{4.3});auto fine=se::evaluate(T,F,se::Config{4.3,256,96,64});
        const auto a=mg::assess_rejected_supply(base,{30,27,0},{0,1,0}),b=mg::assess_rejected_supply(fine,{30,27,0},{0,1,0});
        maximum_refinement=std::max(maximum_refinement,std::abs(a.rejected_fraction-b.rejected_fraction));
        test(std::abs(a.rejected_fraction-b.rejected_fraction)<2e-6,"Guard rate normal spectrum refinement failed");
        test(std::abs(a.rejected_energy_moment_eV-b.rejected_energy_moment_eV)<1e-5,"Rejected Nottingham refinement failed");
    }
    std::ofstream report(out/"guard-assessment-qualification.json");report<<std::setprecision(17)<<"{\"passed\":true,\"checks\":"<<checks<<",\"nominal_guard\":0.2,\"analytic_transverse_integration\":true,\"native_predicate_grid_max_probability_error\":"<<max_probability_error<<",\"native_predicate_grid_max_energy_moment_error_eV\":"<<max_moment_error<<",\"normal_spectrum_refinement_max_absolute_fraction_error\":"<<maximum_refinement<<",\"actual_magnetized_emission_qualified\":false,\"zero_unknown_contributions\":false}\n";
    std::cout<<"Guard assessment checks="<<checks<<'\n';
}
void map(const char* probes,const std::filesystem::path& out,int first,int count) {
    std::ifstream input(probes);need(bool(input),"Missing legacy probes");std::string line;std::getline(input,line);std::vector<std::vector<std::string>> rows;while(std::getline(input,line))rows.push_back(split(line));
    need(first>=0 && count>0 && first+count<=int(rows.size()),"Invalid bounded legacy map slice");std::filesystem::create_directories(out);
    std::ofstream csv(out/("guard-rate-map-"+std::to_string(first)+".csv"));csv<<std::setprecision(17)<<"legacy_row,face,T_K,signed_F_V_m,B_normal_T,B_tangent_T,nominal_MG_current_A_m2,nominal_log_current_A_m2,SN_guard_applies,nominal_rejected_fraction,has_log_rejected_current,log_nominal_rejected_current_A_m2,nominal_rejected_current_A_m2,nominal_rejected_complete_heat_W_m2,nominal_complete_heat_W_m2,geometry_class,actual_magnetized_current_known\n";
    for(int index=first;index<first+count;++index)for(double T:{300.,773.15,1490.}) {
        const auto& r=rows[index];se::Vector n{},B{},E{};for(int j=0;j<3;++j){n[j]=std::stod(r[4+j]);E[j]=std::stod(r[9+j]);B[j]=std::stod(r[12+j]);}
        const double F=se::extraction_field(E,n),Bn=dot(B,n),Bt=std::sqrt(std::max(0.,dot(B,B)-Bn*Bn));
        const auto d=se::evaluate(T,F,se::Config{4.3});const auto a=mg::assess_rejected_supply(d,B,n);
        csv<<index<<','<<r[0]<<','<<T<<','<<F<<','<<Bn<<','<<Bt<<','<<d.current_A_m2<<','<<d.log_current_A_m2<<','<<a.SchottkyNordheim_guard_applies<<',';
        if(a.SchottkyNordheim_guard_applies)csv<<a.rejected_fraction<<','<<a.has_log_rejected_current<<','<<a.log_rejected_current_A_m2<<','<<a.rejected_current_A_m2<<','<<a.rejected_complete_wall_heat_W_m2;
        else csv<<"NA,0,NA,NA,NA";
        csv<<','<<d.wall_heat_W_m2<<','<<(std::abs(Bn)<=1e-12*std::max(1.,std::sqrt(dot(B,B)))?"purely_tangential":"oblique_or_normal_no_tangential_replacement")<<",0\n";
    }
    std::cout<<"Mapped legacy rows "<<first<<".."<<first+count-1<<" at3temperatures; actual magnetized rates remain unknown\n";
}
}
int main(int argc,char** argv) {
    try {
        cpu_set_t mask;CPU_ZERO(&mask);need(sched_getaffinity(0,sizeof(mask),&mask)==0 && CPU_COUNT(&mask)==1 && CPU_ISSET(7,&mask) && getpriority(PRIO_PROCESS,0)>=15,"CPU7/nice15 required");
        need(argc>=3,"Expected qualify or map arguments");
        if(std::string(argv[1])=="qualify"){need(argc==3,"qualify OUT");qualify(argv[2]);}
        else {need(argc==6 && std::string(argv[1])=="map","map PROBES OUT FIRST COUNT");map(argv[2],argv[3],std::stoi(argv[4]),std::stoi(argv[5]));}
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
