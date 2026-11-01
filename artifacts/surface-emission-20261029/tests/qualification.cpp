#include "surface_emission.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include <sched.h>
#include <sys/resource.h>
namespace se=surface_emission;
namespace {
std::vector<std::string> passed;
void check(bool okay,const std::string& name){if(!okay)throw std::runtime_error(name);passed.push_back(name);}
template<class F>void rejects(F f,const std::string& name){bool failed=false;try{f();}catch(const std::exception&){failed=true;}check(failed,name);}
bool close(double a,double b,double rtol,double atol=0){return std::abs(a-b)<=atol+rtol*std::max(std::abs(a),std::abs(b));}
double radical(unsigned i,unsigned base){double f=1,out=0;while(i){f/=base;out+=f*(i%base);i/=base;}return out;}
se::Stage stage(double T,double F) {
    return {1,1,1,1,1,1,0,1e-12,0,0,0,{0,-F,0},{0,1,0},T,1e-8,true,true,true};
}
std::vector<std::string> split(const std::string& line){std::vector<std::string> out;std::istringstream in(line);std::string s;while(std::getline(in,s,','))out.push_back(s);return out;}
double quantile(std::vector<double> data,double fraction){std::sort(data.begin(),data.end());return data.at(std::size_t(fraction*(data.size()-1)));}
}
int main(int argc,char** argv) {
    std::filesystem::path output;
    try {
        if(argc<2)throw std::invalid_argument("Expected report directory and optional original probe CSV");
        output=argv[1];std::filesystem::create_directories(output);
        cpu_set_t mask;CPU_ZERO(&mask);
        check(sched_getaffinity(0,sizeof(mask),&mask)==0 && CPU_COUNT(&mask)==1 && CPU_ISSET(7,&mask) && getpriority(PRIO_PROCESS,0)>=15,"CPU7/nice15 bounded test policy");
        const se::Config config{4.3};const double k=se::kb/se::qe;
        check(close(se::ideal_richardson_A_m2_K2(),1.2017323e6,1e-7),"Ideal free-electron Richardson SI prefactor");
        for(double T:{300.,773.15,1000.,1500.,2000.}) {
            const auto d=se::evaluate(T,0,config);const double reference=se::ideal_richardson_A_m2_K2()*T*T*std::exp(-4.3/(k*T));
            check(close(d.current_A_m2,reference,1e-10),"Zero-field Richardson limit T="+std::to_string(T));
            check(close(d.mean_metal_energy_eV,4.3+2*k*T,1e-10),"Thermionic reservoir energy includes phi+2kT once T="+std::to_string(T));
            const auto retarding=se::evaluate(T,-1e8,config);
            check(retarding.current_A_m2==d.current_A_m2 && retarding.current_A_m2>0,"Retarding surface keeps outward thermal supply T="+std::to_string(T));
        }
        check(se::evaluate(0,0,config).current_A_m2==0 && se::evaluate(0,-1e8,config).nodes.empty(),"Cold retarding/zero field emits no electrons");
        const auto nearzero=se::evaluate(2000,1,config),zero=se::evaluate(2000,0,config);
        check(close(nearzero.current_A_m2,zero.current_A_m2,4e-4),"Small positive field approaches Richardson without sum double count");
        const auto schottky=se::evaluate(2000,1e6,config);
        const double RS=se::ideal_richardson_A_m2_K2()*4e6*std::exp(-(4.3-se::schottky_lowering_eV(1e6))/(k*2000));
        check(close(schottky.current_A_m2,RS,1e-3),"Low-field hot Richardson-Schottky limit");
        const double F=1e9,G=se::action(0,F,config),derivative=(se::action(1e-4,F,config)-se::action(-1e-4,F,config))/2e-4;
        const auto cold=se::evaluate(0,F,config);
        const double asymptotic=se::ideal_richardson_A_m2_K2()/(k*k)*std::exp(-G)/(derivative*derivative);
        check(close(cold.current_A_m2,asymptotic,.04),"Cold Murphy-Good linear-action asymptotic with image barrier");
        check(cold.wall_heat_W_m2>0 && cold.mean_metal_energy_eV<0,"Cold-field Nottingham heating from same spectrum");
        check(se::evaluate(1500,F,config).wall_heat_W_m2<0,"Thermal regime Nottingham cooling sign");
        const double top=4.3-se::schottky_lowering_eV(F);
        check(close(se::transmission(top,F,config),.5,1e-13),"Kemble saddle transmission one half");
        check(se::transmission(top-1e-5,F,config)<.5 && se::transmission(top+1e-5,F,config)>.5,"Continuous near-top tunneling/flyover branch");
        check(close(se::action(top-1e-5,F,config),-se::action(top+1e-5,F,config),.001),"Action slope matches inverted-parabola continuation");
        const double triangular=6.830890e9*std::pow(4.3,1.5)/1e5;
        check(close(se::action(0,1e5,config),triangular,5e-5),"Small scaled-field SN action approaches triangular analytic barrier");
        double max_convergence=0;
        for(const auto& condition:std::vector<std::pair<double,double>>{{0,1e9},{300,3e9},{1500,2e8},{2500,1e9},{3000,5e9},{300,1e6}}) {
            auto c=config;c.energy_panels=64;c.action_nodes=32;const auto coarse=se::evaluate(condition.first,condition.second,c);
            c.energy_panels=256;c.action_nodes=128;c.tail_widths=64;const auto fine=se::evaluate(condition.first,condition.second,c);
            const double error=std::abs(std::expm1(coarse.log_current_A_m2-fine.log_current_A_m2));max_convergence=std::max(max_convergence,error);
            check(error<2e-4,"Current quadrature/action/tail refinement T="+std::to_string(condition.first));
            check(std::abs(coarse.mean_metal_energy_eV-fine.mean_metal_energy_eV)<2e-4*4.3,"Energy moment quadrature/action/tail refinement T="+std::to_string(condition.first));
        }
        for(const auto& n:std::vector<se::Vector>{{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}}) {
            se::Vector E{-2e7*n[0],-2e7*n[1],-2e7*n[2]};check(se::extraction_field(E,n)==2e7,"Material-outward electron-force sign on face");
            for(double& v:E)v=-v;
            check(se::extraction_field(E,n)==-2e7,"Retarding sign on opposite face");
        }
        rejects([&]{se::extraction_field({1,2,3},{0,2,0});},"Nonunit surface normal rejected");
        rejects([&]{se::evaluate(1e6,1e7,config);},"Formal overheated original history rejected, never clamped");
        rejects([&]{se::evaluate(1000,6e9,config);},"Field outside numerical domain rejected");
        rejects([&]{se::evaluate(1000,4e9,se::Config{2});},"Image barrier collapse rejected");
        rejects([&]{se::evaluate(std::numeric_limits<double>::quiet_NaN(),1e6,config);},"NaN temperature rejected");
        auto s=stage(1500,1e8);const auto flux=se::prepare(s,config);
        check(close(flux.expected_physical_electrons,flux.spectrum.current_A_m2*s.area_m2*s.dt/se::qe,1e-15),"Area/time/depth counted exactly once");
        auto wrong=s;wrong.field_time=wrong.dt/2;rejects([&]{se::prepare(wrong,config);},"Mismatched solved field and evolved T stage refused");
        wrong=s;wrong.thermal_accepted_state=2;rejects([&]{se::prepare(wrong,config);},"Mixed accepted revisions refused");
        wrong=s;wrong.material_temperature_phase_admitted=false;rejects([&]{se::prepare(wrong,config);},"Unadmitted material phase/T refused");
        wrong=s;wrong.locally_planar_barrier_admitted=false;rejects([&]{se::prepare(wrong,config);},"Unadmitted microscopic planar geometry refused");
        se::MatchingPlane plane;rejects([&]{se::sample(flux,plane,1,.3,.4,.5,1);},"Missing microscopic matching/energy ownership fails closed");
        plane={se::PlaneKind::ExplicitSchottkyNordheimHandoff,se::required_SN_matching_distance_m(flux.spectrum)*1.01,true,true,true};
        std::vector<se::Packet> packets;packets.reserve(8192);long double metal=0,K=0;
        for(unsigned i=1;i<=8192;++i) {
            const auto p=se::sample(flux,plane,1,radical(i,2),radical(i,3),radical(i,5),i);
            const double scale=std::max({std::abs(p.ledger.wall_heat_J),p.ledger.outgoing_K_J,4.3*se::qe});
            if(std::abs(p.ledger.residual_J())>2e-13*scale)throw std::runtime_error("Packet boundary ledger closure");
            if(!close(se::kinetic_energy_J(p.proper_velocity_m_s),p.ledger.outgoing_K_J,1e-12))throw std::runtime_error("Proper-velocity particle energy equals ledger");
            if(p.proper_velocity_m_s[1]<0)throw std::runtime_error("Emission points into metal");
            metal+=p.normal_metal_energy_eV+p.transverse_energy_eV;K+=p.ledger.outgoing_K_J/se::qe;packets.push_back(p);
        }
        check(close(double(metal/8192),flux.spectrum.mean_metal_energy_eV,.002),"Sampled joint FD normal/transverse spectrum matches integrated Nottingham moment");
        check(std::abs(double(K/8192)-(flux.spectrum.mean_metal_energy_eV-4.3+1e8*plane.outward_distance_m+se::qe/(16*3.14159265358979323846*se::eps0*plane.outward_distance_m)))<.003,
              "Sampled launch energy includes declared reference/image/field matching work");
        const auto thermal_flux=se::prepare(stage(1500,0),config);
        const se::MatchingPlane thermal_surface{se::PlaneKind::ThermalVacuumSurfaceStep,0,true,true,true};
        const auto field_flux=se::prepare(stage(0,3e9),config);
        const se::MatchingPlane field_surface{se::PlaneKind::ExplicitSchottkyNordheimHandoff,se::required_SN_matching_distance_m(field_flux.spectrum)*1.01,true,true,true};
        long double cosine2=0,thermal_K=0,field_metal=0,field_azimuth_x=0,field_azimuth_z=0;
        for(unsigned i=1;i<=8192;++i) {
            const auto a=se::sample(thermal_flux,thermal_surface,1,radical(i,2),radical(i,3),radical(i,5),i);
            const double u2=std::inner_product(a.proper_velocity_m_s.begin(),a.proper_velocity_m_s.end(),a.proper_velocity_m_s.begin(),0.);
            cosine2+=a.proper_velocity_m_s[1]*a.proper_velocity_m_s[1]/u2;thermal_K+=a.ledger.outgoing_K_J/se::qe;
            const auto b=se::sample(field_flux,field_surface,1,radical(i,2),radical(i,3),radical(i,5),i);
            if(!(b.normal_metal_energy_eV<0 && b.ledger.outgoing_K_J>0 && b.proper_velocity_m_s[1]>0))throw std::runtime_error("Cold field spectrum/launch differs from zero-T half-Maxwellian");
            field_metal+=b.normal_metal_energy_eV+b.transverse_energy_eV;
            field_azimuth_x+=b.proper_velocity_m_s[0];field_azimuth_z+=b.proper_velocity_m_s[2];
        }
        check(std::abs(double(cosine2/8192)-.5)<.002,"Zero-field Boltzmann limit recovers cosine angular law");
        check(close(double(thermal_K/8192),2*k*1500,.002),"Only zero-field Boltzmann limit recovers flux half-Maxwellian mean2kT");
        check(std::abs(double(field_metal/8192)-field_flux.spectrum.mean_metal_energy_eV)<.002,"Cold-field phase-space spectrum samples negative Fermi-relative energies");
        check(std::abs(double(field_azimuth_x/8192))<500 && std::abs(double(field_azimuth_z/8192))<500,"Cold-field tangential azimuth remains symmetric");
        check(se::transverse_quantile_eV(-.2,0,.25)==.05 && close(se::transverse_quantile_eV(-.2,0,.75),.15,1e-15),"Cold conditional transverse energy is uniform on[0,-epsilon]");
        const auto negligible=se::evaluate(0,1e6,config);
        check(negligible.current_underflow && negligible.current_A_m2==0 && negligible.has_log_current && std::isfinite(negligible.log_current_A_m2),"Negligible cold-field current keeps finite logarithmic magnitude and explicit underflow flag");
        const auto total=se::sum_packets(packets);
        check(std::abs(total.residual_J())<1e-11*(std::abs(total.wall_heat_J)+total.outgoing_K_J+std::abs(total.electron_reference_change_J)),"Weighted packet batch energy closure");
        check(total.particle_charge_C+total.electrode_charge_C==0,"Particle/electrode signed charge closes");
        auto duplicate=packets;duplicate.push_back(packets.front());rejects([&]{se::sum_packets(duplicate);},"Duplicate source event refused");
        const auto retarding=se::prepare(stage(1500,-1e7),config);
        const se::MatchingPlane surface{se::PlaneKind::ThermalVacuumSurfaceStep,0,true,true,true};
        const auto emission=se::sample(retarding,surface,3.25,.31,.43,.61,99);
        const auto returned=se::collect_absorbed(retarding.spectrum,surface,emission.ledger.outgoing_K_J/emission.weight,emission.weight,true);
        check(std::abs(emission.ledger.wall_heat_J+returned.wall_heat_J)<1e-13*std::abs(emission.ledger.wall_heat_J),"Retarded returning thermal electron restores sampled wall energy, not mean2kT plus sampledK");
        check(std::abs(emission.ledger.electron_reference_change_J+returned.electron_reference_change_J)<1e-30,"Return reference exchange cancels once");
        const double vn=emission.proper_velocity_m_s[1]/std::sqrt(1+std::pow(emission.proper_velocity_m_s[0]/se::c_light,2)+std::pow(emission.proper_velocity_m_s[1]/se::c_light,2)+std::pow(emission.proper_velocity_m_s[2]/se::c_light,2));
        const double a=se::qe*1e7/se::me,treturn=2*vn/a,position=vn*treturn-.5*a*treturn*treturn;
        check(std::abs(position)<1e-15,"Manufactured retarding PIC flight returns to material surface without deleting source");
        auto double_count=emission.ledger;double_count.wall_heat_J-=4.3*se::qe*emission.weight;
        check(std::abs(double_count.residual_J())>1e-20,"Extra work-function heat debit is detected as double counting");
        rejects([&]{se::collect_absorbed(retarding.spectrum,surface,1e-19,1,false);},"Unadmitted collection/secondary/reflection ownership refused");
        se::MatchingPlane too_close=plane;too_close.outward_distance_m=1e-12;
        rejects([&]{se::sample(flux,too_close,1,.3,.4,.5,1);},"Spectrum truncating microscopic plane refused");
        if(argc>=3) {
            std::ifstream in(argv[2]);std::string line;std::getline(in,line);std::vector<double> signed_fields,positive;std::map<std::string,std::vector<double>> faces;
            while(std::getline(in,line)){const auto r=split(line);se::Vector E{},n{};for(int j=0;j<3;++j){n[j]=std::stod(r[4+j]);E[j]=std::stod(r[9+j]);}
                const double field=se::extraction_field(E,n);signed_fields.push_back(field);faces[r[0]].push_back(field);if(field>0)positive.push_back(field);}
            check(signed_fields.size()==1056 && !positive.empty(),"Original1056 field probes loaded for bounded magnitude sensitivity only");
            const std::vector<double> field_scan={*std::min_element(signed_fields.begin(),signed_fields.end()),0,quantile(positive,.5),quantile(positive,.95),*std::max_element(positive.begin(),positive.end())};
            std::ofstream csv(output/"original-field-sensitivity.csv");csv<<std::setprecision(17)<<"T_K,phi_eV,signed_F_V_m,J_A_m2,log_J_A_m2,wall_heat_W_m2,mean_metal_eV,Richardson_J_A_m2\n";
            for(double T:{295.,773.15,1500.,3000.})for(double phi:{4.1,4.3,4.5})for(double field:field_scan) {
                const auto d=se::evaluate(T,field,se::Config{phi});const double rd=se::ideal_richardson_A_m2_K2()*T*T*std::exp(-phi/(k*T));
                csv<<T<<','<<phi<<','<<field<<','<<d.current_A_m2<<','<<d.log_current_A_m2<<','<<d.wall_heat_W_m2<<','<<d.mean_metal_energy_eV<<','<<rd<<'\n';
            }
            std::ofstream field_json(output/"original-field-range.json");field_json<<std::setprecision(17)<<"{\"scope\":\"Legacy finite-offset probes; magnitude sensitivity only, NOT admitted one-sided solved surface field\",\"minimum_signed_extraction_V_m\":"<<field_scan[0]<<",\"positive_median_V_m\":"<<field_scan[2]<<",\"positive_p95_V_m\":"<<field_scan[3]<<",\"maximum_positive_V_m\":"<<field_scan[4]<<",\"faces\":{";
            bool first=true;for(const auto& [face,v]:faces){field_json<<(first?"":",")<<std::quoted(face)<<":["<<*std::min_element(v.begin(),v.end())<<','<<*std::max_element(v.begin(),v.end())<<']';first=false;}field_json<<"}}\n";
        }
        std::ofstream report(output/"qualification.json");report<<std::setprecision(17)<<"{\"passed\":true,\"checks\":"<<passed.size()<<",\"CPU\":7,\"nice\":15,\"GPU_executed\":false,\"physical_foil_model_admitted\":false,\"max_quadrature_refinement_relative_current\":"<<max_convergence<<",\"sampled_metal_mean_eV\":"<<double(metal/8192)<<",\"integrated_metal_mean_eV\":"<<flux.spectrum.mean_metal_energy_eV<<",\"cold_numeric_over_MG_asymptotic\":"<<cold.current_A_m2/asymptotic<<",\"passed_checks\":[";
        for(std::size_t i=0;i<passed.size();++i)report<<(i?",":"")<<std::quoted(passed[i]);
        report<<"]}\n";
        std::cout<<"{\"passed\":true,\"checks\":"<<passed.size()<<",\"physical_foil_model_admitted\":false}\n";
    }catch(const std::exception& error){
        if(!output.empty()){std::ofstream file(output/"qualification-failure.json");file<<"{\"passed\":false,\"error\":"<<std::quoted(error.what())<<",\"passed_checks_before_failure\":"<<passed.size()<<"}\n";}
        std::cerr<<error.what()<<'\n';return 1;
    }
}
