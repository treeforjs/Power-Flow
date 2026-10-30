#pragma once
#include "axial_system.hpp"
#include <iostream>
namespace axial {
struct Sample {State field,initial;std::vector<GapParticle>arena;double Tmax,heat,voltage;};
inline Sample time_sample(const Mesh&g,const CurrentTrace&trace,double dt,int steps){
    Thermal thermal(g,trace.start_s());Fields fields(g,thermal.conductivity(),trace.current(trace.start_s()));auto state=fields.start(trace.start_s());
    Sample result;result.initial=state;std::vector<GapParticle>arena;
    for(int sp=0;sp<2;++sp){GapParticle p;p.id=17+sp;p.position_m={.00013,-.00021,.003};p.proper_velocity_m_s={1e5,2e5,3e5};
        p.charge_C=sp?pic_e:-pic_e;p.mass_kg=sp?hall::proton_mass:pic_me;p.weight=1e7;arena.push_back(p);}
    for(int k=0;k<steps;++k){auto trial=prepare_joint(g,fields,state,arena,thermal,trace,dt);result.voltage=trial.receipt.terminal_V;
        hall::require(commit_joint(state,arena,thermal,std::move(trial)),"refinement common commit");}
    result.field=std::move(state);result.arena=std::move(arena);result.Tmax=thermal.max_temperature_K();result.heat=thermal.energy_J();return result;
}
inline Json run_refinements(const CurrentTrace&trace){
    Mesh g;auto coarse=time_sample(g,trace,1e-14,16),middle=time_sample(g,trace,5e-15,32),fine=time_sample(g,trace,2.5e-15,64);
    auto weighted=[&](const Vec&a,const Vec&b,const Vec&w){long double s=0;for(std::size_t i=0;i<a.size();++i)s+=(long double)w[i]*(a[i]-b[i])*(a[i]-b[i]);return std::sqrt(double(s));};
    double e01=weighted(coarse.field.E,middle.field.E,g.M),e12=weighted(middle.field.E,fine.field.E,g.M);
    double b01=weighted(coarse.field.B,middle.field.B,g.H),b12=weighted(middle.field.B,fine.field.B,g.H);
    double es=weighted(fine.field.E,fine.initial.E,g.M),bs=weighted(fine.field.B,fine.initial.B,g.H);
    Json result{{"time",{{"dt_s",{1e-14,5e-15,2.5e-15}},{"duration_s",1.6e-13},{"E_difference_over_actual_response",{e01/es,e12/es}},
        {"B_difference_over_actual_response",{b01/bs,b12/bs}},{"E_ratio",e01/e12},{"B_ratio",b01/b12},
        {"Tmax_K",{coarse.Tmax,middle.Tmax,fine.Tmax}},{"thermal_energy_J",{coarse.heat,middle.heat,fine.heat}},
        {"note","Maxwell/native orbit midpoint; accepted-temperature conductivity and thermal diffusion are lagged. No general second-order coupled thermal claim."}}}};
    hall::require(e12<e01&&b12<b01,"coupled temporal differences failed to decrease");
    std::cout<<result.dump()<<std::endl;
    Json space=Json::array();
    for(auto setting:std::vector<std::array<double,3>>{{45,4,.009},{90,4,.009},{180,4,.009},{45,8,.009},{45,16,.009},{85,4,.017},{125,4,.025}}){
      try {
        // Keep the physical axial air padding fixed at 6.25mm when refining z.
        Mesh mesh(int(setting[0]),int(setting[1]),setting[2],int(setting[1])/4);Thermal heat(mesh,trace.start_s());Fields fields(mesh,heat.conductivity(),trace.current(trace.start_s()));
        auto initial=fields.start(trace.start_s());Path sample(mesh,{.0000125,.0000125,0},{0,0,0});
        double magnetic=0;for(int b=0;b<mesh.nb;++b)magnetic+=.5*mesh.H[b]*initial.B[b]*initial.B[b];
        space.push_back({{"mesh",mesh.json()},{"center_B_T",sample.gather_B(initial.B)},{"center_E_V_m",sample.gather_E(initial.E)},
            {"predicted_DC_voltage_V",fields.initial.terminal_V},{"magnetic_energy_J",magnetic},
            {"loop_inductance_from_DC_energy_H",2*magnetic/std::pow(trace.current(trace.start_s()),2)},
            {"interior_Ampere_residual_A_m",fields.initial.magnetic_residual_A_m}});
        std::cout<<space.back().dump()<<std::endl;
        hall::require(std::abs(fields.initial.terminal_V-2*.025*trace.current(trace.start_s())*71.9e-8/(.005*.0002))<=1e-8,"refined DC voltage reference");
      }catch(const std::exception&e){throw std::runtime_error("sensitivity nx="+std::to_string(int(setting[0]))+" nz_foil="+std::to_string(int(setting[1]))+" side="+std::to_string(setting[2])+": "+e.what());}
    }
    result["space_and_exterior"]=space;
    result["spatial_qualification"]="sensitivity measurements only; no mesh/exterior/axial-instability acceptance. Production exterior must continue the independently measured wide graded domain, not revert to 9mm.";
    return result;
}
} // namespace axial
