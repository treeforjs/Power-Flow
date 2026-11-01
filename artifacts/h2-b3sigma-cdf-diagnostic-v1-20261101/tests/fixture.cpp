#include "component.hpp"
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace pfc;
void vec(V3 p) {std::cout<<'['<<p.x<<','<<p.y<<','<<p.z<<']';}
template<class F> void refuses(F f) {bool refused=false;try{f();}catch(const std::runtime_error&){refused=true;}if(!refused)throw std::runtime_error("negative control admitted");}
int main() {
    using namespace h2diag;
    std::cout<<std::setprecision(17);
    const ChemKey key{123,4,5,PFC_PURPOSE_FIXTURE_TEST};
    const V3 direction=v3(0,0,1),zero=v3(0,0,0);
    const auto mode=Mode::diagnostic_primary_cdf_iso;
    const auto frame=EnergyFrame::target_rest;
    refuses([&]{cross_section(8,frame);});
    refuses([&]{cross_section(8,EnergyFrame::laboratory,mode);});
    for(double x:{0.,-1.,4.47713,std::nextafter(energies[0],0.),std::nextafter(150.,151.),1000.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
        refuses([&]{cross_section(x,frame,mode);});
    refuses([&]{sample(8,frame,zero,zero,key,0,mode);});
    refuses([&]{sample(8,frame,v3(0,0,2),zero,key,0,mode);});
    refuses([&]{sample(8,frame,direction,v3(std::numeric_limits<double>::infinity(),0,0),key,0,mode);});
    auto bad=key;bad.purpose=PFC_PURPOSE_MCC;
    refuses([&]{sample(8,frame,direction,zero,bad,0,mode);});
    bad=key;bad.carrier=PFC_CARRIER_LIMIT;
    refuses([&]{sample(8,frame,direction,zero,bad,0,mode);});
    refuses([&]{sample(8,frame,direction,zero,key,PFC_EVENT_LIMIT,mode);});
    for(unsigned i=0;i<point_count;++i){
        std::cout<<"{\"type\":\"rate\",\"E\":"<<energies[i]<<",\"sigma\":"<<cross_section(energies[i],frame,mode)<<"}\n";
        if(i+1<point_count){const double e=(energies[i]+energies[i+1])/2;std::cout<<"{\"type\":\"rate\",\"E\":"<<e<<",\"sigma\":"<<cross_section(e,frame,mode)<<"}\n";}
    }
    // Integral-rate support is broader than the source-conditioned event support.
    for(double E:{4.5,5.,7.,threshold()-1e-8}) {
        if(!(cross_section(E,frame,mode)>0))throw std::runtime_error("rate must remain positive");
        refuses([&]{sample(E,frame,direction,zero,key,0,mode);});
    }
    for(int sign:{-1,1}) {
        const double t=threshold()+sign*1e-7;
        bool ok=true;try{sample(t,frame,direction,zero,key,0,mode);}catch(const std::runtime_error&){ok=false;}
        if(ok!=(sign==1))throw std::runtime_error("source support threshold control");
        std::cout<<"{\"type\":\"threshold\",\"E\":"<<t<<",\"accepted\":"<<ok<<"}\n";
    }
    for(unsigned i=0;i<4096;++i){
        const double E=(i%4==0)?7.2:(i%4==1?10.:(i%4==2?20.:150.));
        // Cold and moving-target tests; energy remains the target-rest operand.
        const V3 u=(i%5==0)?v3(12345,-6789,4321):zero;
        const ChemKey k{0x1234abcdefULL,17,1+i,PFC_PURPOSE_FIXTURE_TEST};
        const auto r=sample(E,frame,direction,u,k,0,mode);
        std::cout<<"{\"type\":\"event\",\"index\":"<<i<<",\"E\":"<<E<<",\"target_u\":";vec(u);
        std::cout<<",\"incoming\":[";vec(r.incoming[0]);std::cout<<',';vec(r.incoming[1]);std::cout<<"],\"outgoing\":[";
        for(unsigned j=0;j<3;++j){if(j)std::cout<<',';vec(r.outgoing[j]);}
        std::cout<<"],\"kinetic_in\":["<<r.kinetic_in[0]<<','<<r.kinetic_in[1]<<"],\"kinetic_out\":["<<r.kinetic_out[0]<<','<<r.kinetic_out[1]<<','<<r.kinetic_out[2]<<"],\"internal\":"<<r.internal_eV<<",\"quantile\":"<<r.quantile<<",\"quantile_hi\":"<<r.quantile_hi<<",\"KER\":"<<r.fragment_KER_eV<<",\"loss\":"<<r.kinematic_loss_eV<<"}\n";
    }
    std::cout<<"{\"type\":\"controls\",\"refusals\":20,\"events\":4096}\n";
}
