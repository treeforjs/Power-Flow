#include "raikou/mfem_gap_pic.h"
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <filesystem>

using namespace raikou::mfem_backend;
namespace {
void check(bool pass,const char* why) { if (!pass) throw std::runtime_error(why); }
void close(double a,double b,double tolerance,const char* why) {
    check(std::isfinite(a) && std::abs(a-b)<=tolerance,why);
}
void vacuum_capacitor() {
    GapPicOptions o; o.nx=7; o.ny=11; o.lower_potential_V=-40; o.upper_potential_V=60;
    GapPic pic(o); pic.solve_fields();
    for (double x : {.01,.31,.99}) for (double y : {.01,.45,.99}) {
        auto E=pic.electric_field(x*o.width_m,y*o.gap_m);
        close(E[0],0,1e-5,"vacuum lateral E");
        close(E[1],-100/o.gap_m,1e-5,"both electrode potentials determine E");
    }
    auto d=pic.diagnostics();
    const double C=8.8541878188e-12*o.width_m*o.depth_m/o.gap_m;
    close(d.lower_field_charge_C,-C*100,1e-21,"lower electrode capacitor charge");
    close(d.upper_field_charge_C,C*100,1e-21,"upper electrode capacitor charge");
    close(d.field_energy_J,.5*C*10000,1e-19,"SI capacitor field energy");
    close(d.gauss_defect_C,0,1e-21,"vacuum global Gauss law");
    pic.set_empty_potentials(-80,120);
    close(pic.electric_field(.4e-3,.5e-3)[1],-200/o.gap_m,2e-5,"empty pulsed bias rescaling");
    pic.set_potentials(10,10); pic.solve_fields();
    close(pic.electric_field(.4e-3,.5e-3)[1],0,1e-5,"equal potential electrodes");
}
void space_charge_and_absorption() {
    GapPicOptions o; o.nx=8; o.ny=12;
    GapPic pic(o);
    GapParticle p; p.position_m={.4e-3,.5e-3,0}; p.weight=1e5;
    pic.add_particle(p); pic.solve_fields();
    auto d=pic.diagnostics();
    close(d.deposited_charge_C,p.charge_C*p.weight,1e-27,"nodal deposition preserves charge");
    close(d.gauss_defect_C,0,1e-24,"particle plus both electrodes satisfies Gauss law");
    check(pic.electric_field(.4e-3,.2e-3)[1]>0,"negative particle pulls field upward below it");
    check(pic.electric_field(.4e-3,.8e-3)[1]<0,"negative particle pulls field downward above it");
    GapPic walls(o);
    p.weight=1; p.position_m={.5e-3,1e-6,0}; p.proper_velocity_m_s={0,-1e6,0};
    walls.add_particle(p);
    p.position_m[1]=o.gap_m-1e-6; p.proper_velocity_m_s[1]=1e6; walls.add_particle(p);
    walls.advance(1e-11);
    d=walls.diagnostics();
    check(walls.particles().empty(),"both electrodes absorb exiting particles");
    close(d.lower_absorbed_charge_C,-pic_e,1e-33,"lower collected charge");
    close(d.upper_absorbed_charge_C,-pic_e,1e-33,"upper collected charge");
    close(d.particle_charge_defect_C,0,1e-32,"surface collection conserves charge");
    check(d.lower_absorbed_energy_J>0 && d.upper_absorbed_energy_J>0,"both surfaces own incident kinetic energy");
}
void magnetic_orbits_and_drift() {
    const double qm=-pic_e/pic_me, B=.01;
    PicVector u{1e7,0,0};
    const double g=std::sqrt(1+1e14/(pic_c*pic_c));
    const double period=2*std::numbers::pi*g/std::abs(qm*B);
    const int steps=2000;
    for (int n=0;n<steps;++n) u=pic_boris(u,qm,{0,0,0},{0,0,B},period/steps);
    close(std::hypot(u[0],u[1]),1e7,2e-6,"Boris magnetic rotation preserves momentum norm");
    check(std::abs(u[1])/1e7<6e-6,"relativistic gyroperiod converges");
    // E along x, B along z => E cross B drift toward negative y.
    const double E=100, omega=std::abs(qm*B), dt=2*std::numbers::pi/(omega*1000);
    u={0,0,0}; double y=0;
    for (int n=0;n<10000;++n) {
        const auto next=pic_boris(u,qm,{E,0,0},{0,0,B},dt);
        y+=.5*dt*(u[1]+next[1]); u=next;
    }
    close(y/(dt*10000),-E/B,.05,"E cross B drift direction and speed");
}
void thermal_boundary_flux() {
    GapPicOptions o; o.nx=2; o.ny=2;
    GapPic pic(o);
    const double T=2000,J=3,dt=1e-12;
    pic.inject_thermal_electrons(GapElectrode::Lower,T,J,dt,30000,123);
    pic.inject_thermal_electrons(GapElectrode::Upper,T,J,dt,30000,456);
    pic.solve_fields();
    const auto d=pic.diagnostics();
    const double Q=-J*o.width_m*o.depth_m*dt;
    close(d.lower_injected_charge_C,Q,1e-32,"lower thermal flux normalization");
    close(d.upper_injected_charge_C,Q,1e-32,"upper thermal flux normalization");
    close(d.particle_charge_defect_C,0,1e-29,"thermal source global charge ledger");
    double energy=0;
    for (const auto& p : pic.particles()) {
        energy+=pic_kinetic_energy_J(p);
        check(p.position_m[1]==0 ? p.proper_velocity_m_s[1]>0 : p.proper_velocity_m_s[1]<0,
              "both incoming distributions point into gap");
    }
    close(energy/pic.particles().size(),2*1.380649e-23*T,.015*2*1.380649e-23*T,
          "flux weighted thermal mean energy is 2 kT");
    close(pic_richardson_supply_A_m2(2000,4.3,1201732.289489034),70.2012,.0001,
          "temperature dependent unrestricted thermionic supply");
}
void surface_profile_source() {
    GapPicOptions o;o.nx=2;o.ny=2;GapPic pic(o);
    const double dt=1e-12;
    pic.inject_surface_profile(GapElectrode::Lower,GapSurface::Gap,{1000,3000},{0,8},dt,1000,1234);
    pic.solve_fields();
    const double expected=-4*o.width_m*o.depth_m*dt;
    close(pic.diagnostics().lower_injected_charge_C,expected,1e-29,"patch integrated source preserves charge");
    for(const auto& p:pic.particles())check(p.position_m[0]>=.5*o.width_m,"zero source patches emit no particles");
    bool rejected=false;
    try{pic.inject_surface_profile(GapElectrode::Lower,GapSurface::Gap,{1000},{-1},dt,1,42);}
    catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"negative patch source rejected");
}
void checkpoint_restart() {
    GapPicOptions o;o.nx=3;o.ny=5;GapPic continuous(o);
    continuous.inject_thermal_electrons(GapElectrode::Lower,2000,1,1e-12,100,45);
    continuous.advance(1e-13);
    const auto path=(std::filesystem::temp_directory_path()/"raikou-gap-pic-test.bin").string();
    continuous.save_checkpoint(path);GapPic resumed(o);resumed.load_checkpoint(path);
    close(resumed.diagnostics().time_s,continuous.diagnostics().time_s,0,"restart preserves exact clock");
    continuous.advance(1e-13);resumed.advance(1e-13);
    check(continuous.particles().size()==resumed.particles().size(),"restart preserves particle count");
    for(std::size_t i=0;i<continuous.particles().size();++i) {
        const auto& a=continuous.particles()[i];const auto& b=resumed.particles()[i];
        check(a.id==b.id && a.weight==b.weight,"restart preserves identities and weights");
        for(int k=0;k<3;++k){close(a.position_m[k],b.position_m[k],1e-20,"restart trajectory position");close(a.proper_velocity_m_s[k],b.proper_velocity_m_s[k],1e-8,"restart trajectory momentum");}
    }
    close(continuous.diagnostics().particle_charge_defect_C,resumed.diagnostics().particle_charge_defect_C,1e-30,"restart charge ledger");
    std::filesystem::remove(path);
    GapPic empty(o);empty.solve_fields();empty.advance_empty(1e-9);close(empty.diagnostics().time_s,1e-9,0,"empty vacuum clock advance");
    bool rejected=false;try{continuous.advance_empty(1e-9);}catch(const std::invalid_argument&){rejected=true;}check(rejected,"cannot skip populated particle evolution");
}
void extrusion_depth_scaling() {
    GapPicOptions a;a.nx=4;a.ny=6;a.depth_m=.01;auto b=a;b.depth_m=.025;
    GapPic first(a),second(b);
    first.inject_thermal_electrons(GapElectrode::Lower,2000,3,1e-12,100,123);
    second.inject_thermal_electrons(GapElectrode::Lower,2000,3,1e-12,100,123);
    first.advance(1e-13);second.advance(1e-13);
    close(second.diagnostics().lower_injected_charge_C,2.5*first.diagnostics().lower_injected_charge_C,1e-28,"25mm extrusion scales source charge");
    close(second.diagnostics().particle_kinetic_energy_J,2.5*first.diagnostics().particle_kinetic_energy_J,1e-26,"25mm extrusion scales kinetic energy");
    auto E=first.electric_field(.5*a.width_m,.5*a.gap_m),F=second.electric_field(.5*a.width_m,.5*a.gap_m);
    for(int k=0;k<3;++k)close(E[k],F[k],1e-8,"extrusion scaling preserves local field and volume density");
}
void finite_axial_short_and_feed() {
    GapPicOptions o;o.nx=o.ny=45;o.width_m=o.gap_m=.009;o.depth_m=.025;
    o.finite_foils=true;o.finite_axial_extent=true;o.lower_foil={.002,.007,.0018,.002};o.upper_foil={.002,.007,.007,.0072};
    GapParticle p;p.position_m={.0045,.0045,o.short_z_m+1e-6};p.proper_velocity_m_s={0,0,-1e6};p.weight=10000;
    GapPic shorted(o);shorted.add_particle(p);shorted.advance(2e-12);
    auto d=shorted.diagnostics();check(shorted.particles().empty(),"out-of-plane short absorbs incident particle");
    close(d.short_absorbed_charge_C,-pic_e*p.weight,1e-30,"short owns collected charge");
    check(d.short_absorbed_energy_J>0,"short owns collected kinetic energy");
    close(d.particle_charge_defect_C,0,1e-29,"short sink closes global charge ledger");
    close(d.maximum_nodal_continuity_defect_C,0,1e-28,"projected continuity includes axial sink");
    p.position_m[2]=o.feed_z_m-1e-6;p.proper_velocity_m_s[2]=1e6;
    GapPic feed(o);feed.add_particle(p);feed.advance(2e-12);d=feed.diagnostics();
    check(feed.particles().empty() && d.axial_escaped_charge_C<0,"feed end records upstream escape");
    close(d.particle_charge_defect_C,0,1e-29,"feed escape closes global charge ledger");
    p.position_m={.0005,.0045,o.short_z_m+1e-6};p.proper_velocity_m_s[2]=-1e6;
    GapPic outside(o);outside.add_particle(p);outside.advance(2e-12);d=outside.diagnostics();
    check(d.short_absorbed_charge_C==0 && d.axial_escaped_charge_C<0,"5x5mm short does not absorb outside its footprint");
    GapPic source(o);source.inject_thermal_electrons(GapElectrode::Upper,2000,1,1e-12,5000,123);
    double mean=0;for(const auto& q:source.particles()){check(q.position_m[2]>=o.short_z_m && q.position_m[2]<=o.feed_z_m,"surface source respects full25mm height");mean+=q.position_m[2];}
    check(std::abs(mean/source.particles().size())<.0004,"surface source spans both halves of observation plane");
    source.solve_fields();const auto path=(std::filesystem::temp_directory_path()/"raikou-axial-checkpoint-test.bin").string();source.save_checkpoint(path);GapPic resumed(o);resumed.load_checkpoint(path);
    check(resumed.particles().size()==source.particles().size(),"finite target checkpoint particle count");
    for(std::size_t i=0;i<source.particles().size();++i)close(resumed.particles()[i].position_m[2],source.particles()[i].position_m[2],0,"finite target checkpoint preserves axial coordinates");
    std::filesystem::remove(path);
    GapPic spatial(o);p.position_m={.0045,.0045,.0001};p.proper_velocity_m_s={0,0,0};p.weight=1;spatial.add_particle(p);
    spatial.set_external_electric_3d([](double,double,double z,double){return PicVector{0,0,1000+1e6*z};});spatial.advance(1e-12);
    close(spatial.particles()[0].proper_velocity_m_s[2],-pic_e/pic_me*1100*1e-12,1e-5,"3D prescribed field gathers actual axial position");
}
void deposited_poisson_refinement() {
    double previous=0;
    for(int n:{4,8,16}) {
        GapPicOptions o; o.nx=3; o.ny=n;
        GapPic pic(o);
        const double rho=-1e-8,volume=o.width_m*o.gap_m*o.depth_m;
        for(int j=0;j<n;++j) for(int i=0;i<o.nx;++i) {
            GapParticle p; p.position_m={(i+.5)*o.width_m/o.nx,(j+.5)*o.gap_m/n,0};
            p.weight=-rho*volume/(pic_e*o.nx*n); pic.add_particle(p);
        }
        pic.solve_fields();
        mfem::FunctionCoefficient exact([&](const mfem::Vector& x) {
            return rho*x[1]*(o.gap_m-x[1])/(2*8.8541878188e-12);
        });
        const double error=pic.potential().ComputeL2Error(exact);
        if(previous>0) check(previous/error>3.99,"deposited uniform-charge Poisson second-order potential convergence");
        previous=error;
    }
}
void electric_acceleration_and_guards() {
    GapPicOptions o; o.nx=3; o.ny=4; o.upper_potential_V=1;
    GapPic pic(o); GapParticle p; p.charge_C=0; p.position_m={.1e-3,.5e-3,0};
    p.proper_velocity_m_s={-1e6,0,0}; pic.add_particle(p); pic.advance(4e-10);
    close(pic.particles()[0].position_m[0],.3e-3,3e-9,"specular side reflection");
    auto u=pic_boris({0,0,0},-pic_e/pic_me,{0,-1000,0},{0,0,0},1e-12);
    close(u[1],pic_e/pic_me*1000e-12,1e-9,"electron accelerates toward higher potential");
    bool rejected=false; try { pic.advance(0); } catch(const std::invalid_argument&) { rejected=true; }
    check(rejected,"reject zero time step");
    rejected=false; try { pic.inject_thermal_electrons(GapElectrode::Lower,0,1,1,1,0); }
    catch(const std::invalid_argument&) { rejected=true; }
    check(rejected,"reject invalid surface temperature");
}
void path_continuity() {
    GapPicOptions o; o.nx=8; o.ny=11;
    GapPic pic(o);
    for(int i=0;i<15;++i) {
        GapParticle p;
        p.position_m={(.04+.061*i)*o.width_m,(.1+.05*i)*o.gap_m,0};
        p.proper_velocity_m_s={i%2==0 ? 2e7 : -2e7,i%3==0 ? -1e7 : 1e7,3e6};
        p.weight=1e-6;
        pic.add_particle(p);
    }
    pic.advance(9e-11);
    auto d=pic.diagnostics();
    close(d.maximum_nodal_continuity_defect_C,0,1e-37,
          "path currents close nodal continuity across cells, reflections and absorption");
    close(d.particle_charge_defect_C,0,1e-37,"multi-path particle balance");
    bool current=false;
    for(double q:pic.horizontal_transfers_C()) current=current || q!=0;
    check(current,"oriented path current is deposited");
}
void finite_foil_edges() {
    GapPicOptions o; o.nx=45;o.ny=45;o.width_m=.009;o.gap_m=.009;o.depth_m=.025;
    o.finite_foils=true;o.lower_foil={.002,.007,.0018,.002};o.upper_foil={.002,.007,.007,.0072};
    o.upper_potential_V=-20;
    GapPic pic(o);pic.solve_fields();
    auto left=pic.electric_field(.0015,.0065),right=pic.electric_field(.0075,.0065);
    check(std::abs(left[0])>10,"finite-width edge produces lateral fringing field");
    close(left[0],-right[0],1e-5,"mirror-symmetric foil edge fields");
    close(pic.diagnostics().gauss_defect_C,0,1e-20,"finite conductor reaction charges close Gauss law");
    GapParticle p;p.position_m={.001,.0065,0};p.proper_velocity_m_s={0,2e7,0};p.weight=1e-8;
    pic.add_particle(p);pic.advance(5e-11);
    check(pic.particles().size()==1 && pic.particles()[0].position_m[1]>.0072,
          "particle outside foil width passes around electrode edge");
    GapPic hit(o);p.position_m={.0019,.0071,0};p.proper_velocity_m_s={2e7,0,0};
    hit.add_particle(p);hit.advance(2e-11);
    check(hit.particles().empty() && hit.diagnostics().upper_absorbed_charge_C<0,
          "finite-thickness side face absorbs first ray intersection");
    GapPic escaping(o);p.position_m={.0089,.0045,0};p.proper_velocity_m_s={2e7,0,0};
    escaping.add_particle(p);escaping.advance(2e-11);
    check(escaping.particles().empty() && escaping.diagnostics().escaped_charge_C<0,
          "far vacuum boundary absorbs and records escape");
    close(escaping.diagnostics().particle_charge_defect_C,0,1e-38,"escape charge ledger closes");
    close(escaping.diagnostics().maximum_nodal_continuity_defect_C,0,1e-37,"edge and escape path continuity");
    GapPic emission(o);
    for(auto side:{GapSurface::Gap,GapSurface::Back,GapSurface::Left,GapSurface::Right})
        emission.inject_thermal_electrons(GapElectrode::Upper,2000,1,1e-12,100,123,side);
    emission.advance(1e-13);
    check(emission.particles().size()==400,"all four foil faces emit outward into vacuum");
}
}
int main() {
    try {
        vacuum_capacitor(); space_charge_and_absorption(); magnetic_orbits_and_drift();
        thermal_boundary_flux(); electric_acceleration_and_guards(); path_continuity();
        deposited_poisson_refinement();surface_profile_source();checkpoint_restart();extrusion_depth_scaling();finite_axial_short_and_feed();
        finite_foil_edges();
        std::cout<<"Two-electrode MFEM PIC foundation tests passed\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
