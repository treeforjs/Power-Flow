#include "component.hpp"
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace pfc;
void vec(V3 p) {std::cout<<'['<<p.x<<','<<p.y<<','<<p.z<<']';}
template<class F> void refuses(F f) {bool refused=false;try{f();}catch(const std::runtime_error&){refused=true;}if(!refused)throw std::runtime_error("negative control admitted");}
int main() {
    using namespace ohdiag;
    std::cout<<std::setprecision(17);
    const ChemKey key{123,4,5,PFC_PURPOSE_FIXTURE_TEST};
    const V3 direction=v3(0,0,1),zero=v3(0,0,0);
    const auto mode=Mode::diagnostic_vertical_iso;
    const auto frame=EnergyFrame::target_rest;
    refuses([&]{cross_section(8,frame);});
    refuses([&]{cross_section(8,EnergyFrame::laboratory,mode);});
    for(double x:{0.,-1.,7.3537,std::nextafter(energies[0],0.),std::nextafter(10.,11.),100.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
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
    // Isolate the actual solver's threshold from the more conservative source window.
    TableView table{};table.n_species=4;table.species=species;
    GroupRec group{};ProcessRec process{};
    process.kind=PFC_P_DISSOC;process.angular=PFC_ANG_ISO;process.n_out=3;
    process.out_species[0]=0;process.out_species[1]=2;process.out_species[2]=3;
    process.loss_eV=loss_eV;process.q_eV=reference_cost_eV;process.chain_q_eV=reference_cost_eV;process.kinetic_fraction=1;
    for(int sign:{-1,1}) {
        const double t=threshold()+sign*1e-8;V3 p[PFC_MAX_OUT];double masses[PFC_MAX_OUT],internal;
        const int ok=event_kinematics(table,process,group,0,key,0,species[0].mass_eV,species[1].mass_eV,
            v3(0,0,std::sqrt(t*(t+2*species[0].mass_eV))),t,zero,1,p,masses,&internal);
        if(ok!=(sign==1))throw std::runtime_error("solver threshold control");
        std::cout<<"{\"type\":\"threshold\",\"E\":"<<t<<",\"accepted\":"<<ok<<"}\n";
    }
    for(unsigned i=0;i<2048;++i){
        const double E=(i%3==0)?energies[0]:(i%3==1?8.5:10.);
        // Cold and moving-target tests; energy remains the target-rest operand.
        const V3 u=(i%4==0)?v3(12345,-6789,4321):zero;
        const ChemKey k{0x1234abcdefULL,17,1+i,PFC_PURPOSE_FIXTURE_TEST};
        const auto r=sample(E,frame,direction,u,k,0,mode);
        std::cout<<"{\"type\":\"event\",\"index\":"<<i<<",\"E\":"<<E<<",\"target_u\":";vec(u);
        std::cout<<",\"incoming\":[";vec(r.incoming[0]);std::cout<<',';vec(r.incoming[1]);std::cout<<"],\"outgoing\":[";
        for(unsigned j=0;j<3;++j){if(j)std::cout<<',';vec(r.outgoing[j]);}
        std::cout<<"],\"kinetic_in\":["<<r.kinetic_in[0]<<','<<r.kinetic_in[1]<<"],\"kinetic_out\":["<<r.kinetic_out[0]<<','<<r.kinetic_out[1]<<','<<r.kinetic_out[2]<<"],\"internal\":"<<r.internal_eV<<"}\n";
    }
    std::cout<<"{\"type\":\"controls\",\"refusals\":16,\"events\":2048}\n";
}
