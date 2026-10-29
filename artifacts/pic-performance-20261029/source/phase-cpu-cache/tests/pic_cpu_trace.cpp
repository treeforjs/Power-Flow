#include "raikou/mfem_gap_pic.h"
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace raikou::mfem_backend;
int main() {
    GapPicOptions options;options.nx=13;options.ny=17;
    options.upper_potential_V=3.;options.magnetic_field_T={.1,.2,.3};
    GapPic pic(options);
    for(int n=0;n<120;++n) {
        GapParticle p;p.position_m={(.01+.008*n)*options.width_m,(.07+.007*n)*options.gap_m,0};
        p.proper_velocity_m_s={n%2?2e7:-2e7,n%3?1e7:-1e7,3e6};p.weight=1e-5;
        pic.add_particle(p);
    }
    std::cout<<std::hexfloat;
    for(int step=0;step<12;++step) {
        if(step%3==0)pic.inject_thermal_electrons(GapElectrode::Lower,2000,1,1e-12,32,100+step);
        pic.advance(5e-12);
        const auto d=pic.diagnostics();
        std::cout<<step<<' '<<d.time_s<<' '<<d.particle_charge_C<<' '<<d.deposited_charge_C<<' '
            <<d.gauss_defect_C<<' '<<d.particle_charge_defect_C<<' '<<d.field_energy_J<<' '
            <<d.particle_kinetic_energy_J<<' '<<d.maximum_nodal_continuity_defect_C<<'\n';
        for(const auto& p:pic.particles()) {
            std::cout<<p.id<<' '<<p.weight;
            for(double x:p.position_m)std::cout<<' '<<x;
            for(double u:p.proper_velocity_m_s)std::cout<<' '<<u;
            std::cout<<'\n';
        }
        for(int i=0;i<pic.potential().Size();++i)std::cout<<pic.potential()[i]<<' ';
        std::cout<<'\n';
        for(double q:pic.horizontal_transfers_C())std::cout<<q<<' ';
        std::cout<<'\n';
        for(double q:pic.vertical_transfers_C())std::cout<<q<<' ';
        std::cout<<'\n';
    }
}
