#include "raikou/mfem_gap_pic.h"
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <cstring>
using namespace raikou::mfem_backend;
void trace(GapPic& p) {
    const auto d=p.diagnostics();
    static_assert(sizeof(d)%sizeof(double)==0);
    for(std::size_t i=0;i<sizeof(d)/sizeof(double);++i){double v;std::memcpy(&v,reinterpret_cast<const char*>(&d)+i*sizeof(double),sizeof(v));std::cout<<v<<" ";}std::cout<<"\n";
    for(int i=0;i<p.potential().Size();++i)std::cout<<p.potential()[i]<<" ";std::cout<<"\n";
    const auto& o=p.options();for(double x:{0.,.1,.5,.9,1.})for(double y:{0.,.1,.5,.9,1.}){auto e=p.electric_field(x*o.width_m,y*o.gap_m);for(double v:e)std::cout<<v<<" ";}std::cout<<"\n";
    for(const auto& a:p.particles()){std::cout<<a.id<<" "<<a.charge_C<<" "<<a.mass_kg<<" "<<a.weight<<" ";for(double v:a.position_m)std::cout<<v<<" ";for(double v:a.proper_velocity_m_s)std::cout<<v<<" ";}std::cout<<"\n";
    for(double v:p.horizontal_transfers_C())std::cout<<v<<" ";std::cout<<"\n";for(double v:p.vertical_transfers_C())std::cout<<v<<" ";std::cout<<"\n";
}
int main(int argc,char** argv){
    if(argc!=2)return 2;std::cout<<std::hexfloat;
    GapPicOptions o;o.nx=9;o.ny=11;o.lower_potential_V=-2;o.upper_potential_V=3;
    GapPic p(o);p.solve_fields();trace(p);p.set_empty_potentials(4,-6);trace(p);p.release_field_solver();trace(p);p.solve_fields();trace(p);
    for(int i=0;i<20;++i){GapParticle a;a.position_m={(.03+.045*i)*o.width_m,(.05+.04*i)*o.gap_m,0};a.proper_velocity_m_s={i%2?2e7:-2e7,i%3?-1e7:1e7,3e6};a.weight=1e-4;p.add_particle(a);}
    for(int s=0;s<5;++s){p.set_potentials(s*.3,-s*.7);p.advance(1e-11);trace(p);}
    const auto checkpoint=(std::filesystem::path(argv[1])/"restart.bin").string();std::filesystem::create_directories(argv[1]);p.save_checkpoint(checkpoint);
    auto opts=p.options();GapPic resumed(opts);resumed.load_checkpoint(checkpoint);trace(resumed);p.advance(1e-12);resumed.advance(1e-12);trace(p);trace(resumed);
    GapPicOptions f;f.nx=f.ny=45;f.width_m=f.gap_m=.009;f.depth_m=.025;f.finite_foils=true;f.finite_axial_extent=true;f.lower_foil={.002,.007,.0018,.002};f.upper_foil={.002,.007,.007,.0072};f.upper_potential_V=-10;
    GapPic foil(f);foil.solve_fields();trace(foil);
    for(auto electrode:{GapElectrode::Lower,GapElectrode::Upper})for(auto side:{GapSurface::Gap,GapSurface::Back,GapSurface::Left,GapSurface::Right})foil.inject_thermal_electrons(electrode,2000,1,1e-12,3,123,side);
    foil.advance(1e-13);trace(foil);foil.set_potentials(1,-9);foil.advance(1e-13);trace(foil);
}
