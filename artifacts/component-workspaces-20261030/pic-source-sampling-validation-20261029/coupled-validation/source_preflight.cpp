#include "native_source_only.hpp"
#include "records.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
namespace {
using namespace cv;
std::vector<Event> source_events(const std::filesystem::path& path){
    std::ifstream f(path,std::ios::binary);std::uint64_t magic=0,n=0;read(f,magic);read(f,n);
    require(magic==0x5352434556543031ULL&&n>0&&n<=100000,"invalid source-only fixture");std::vector<Event> result;
    for(std::uint64_t k=0;k<n;++k){Event e;read(f,e.injection_time_s);read(f,e.span_s);read(f,e.native_seed);read(f,e.step);int face=0,lower=0;std::uint64_t size=0;read(f,face);read(f,lower);read(f,size);
        require(face>=0&&face<4&&(lower==0||lower==1),"invalid source face");
        e.electrode=lower?GapElectrode::Lower:GapElectrode::Upper;
        e.surface=face==0?GapSurface::Back:face==1?GapSurface::Gap:face==2?GapSurface::Left:GapSurface::Right;
        vector_read(f,e.T,size);vector_read(f,e.J,size);valid(e);result.push_back(std::move(e));}
    require(f.peek()==std::char_traits<char>::eof(),"trailing source fixture");return result;
}
}
int main(int argc,char** argv){try{
    require(argc==4,"Usage: source_preflight COPIED_INPUT EVENTS OWNED_OUTPUT");resource_check();require_owned_output(argv[3]);
    auto initial=read_input(argv[1]);auto events=source_events(argv[2]);std::uint64_t comparisons=0;
    long double target=0;for(const auto& e:events)target+=physical_electrons(e,initial.options);
    double maximum_event_error=0,maximum_total_error=0;std::array<std::array<long double,3>,3> totals{};
    for(int seed=0;seed<3;++seed)for(std::size_t index=0;index<events.size();++index){const auto& e=events[index];
        GapPic ref(initial.options);ref.next_particle_id_=initial.next_id;
        ref.inject_surface_profile(e.electrode,e.surface,e.T,e.J,e.span_s,64,e.native_seed+std::uint64_t(seed)*0x9e3779b97f4a7c15ULL);
        for(std::size_t m=0;m<counts.size();++m){int count=counts[m];GapPic sampler(initial.options);sampler.next_particle_id_=initial.next_id;
            sampler.inject_surface_profile(e.electrode,e.surface,e.T,e.J,e.span_s,count,e.native_seed+std::uint64_t(seed)*0x9e3779b97f4a7c15ULL);
            require(sampler.particles_.size()==std::size_t(count),"wrong macro count");long double sum=0;
            for(int n=0;n<count;++n){const auto& p=sampler.particles_[n];const auto& r=ref.particles_[n];
                require(p.position_m==r.position_m&&p.proper_velocity_m_s==r.proper_velocity_m_s&&p.weight*count==r.weight*64,"native sampling prefix differs");sum+=p.weight;++comparisons;}
            maximum_event_error=std::max(maximum_event_error,double(std::abs(sum-physical_electrons(e,initial.options))/physical_electrons(e,initial.options)));totals[seed][m]+=sum;}}
    for(const auto& seed:totals)for(auto n:seed)maximum_total_error=std::max(maximum_total_error,double(std::abs(n-target)/target));
    require(maximum_event_error<=1e-12&&maximum_total_error<=1e-12,"source physical charge changed");
    GapPic zero(initial.options);zero.inject_surface_profile(GapElectrode::Upper,GapSurface::Gap,{1500,2000},{0,0},1e-12,16,42);
    require(zero.particles_.empty(),"zero current source emitted macros");
    Trace fixture;fixture.initial_time_s=0;fixture.end_time_s=1e-14;fixture.seed_index=0;fixture.input_identity=std::string(64,'0');
    auto event=events.front();event.injection_time_s=0;event.step=initial.step;fixture.steps.push_back({initial.step,0,1e-14,{event}});
    std::filesystem::path out=argv[3];trace_write(out/"manufactured-trace.bin",fixture);auto readback=trace_read(out/"manufactured-trace.bin");
    require(!readback.coupled_verified&&same_event(readback.steps[0].events[0],event),"trace roundtrip/provenance lost");
    bool invalid=false;auto bad=fixture;bad.steps[0].events[0].injection_time_s=1e-14;
    try{valid(bad);}catch(const std::runtime_error&){invalid=true;}require(invalid,"post-advance injection was accepted");
    invalid=false;bad=fixture;bad.steps[0].dt_s=-1e-14;try{valid(bad);}catch(const std::runtime_error&){invalid=true;}require(invalid,"negative trace interval accepted");
    std::ofstream f(out/"source-preflight.json");f<<std::setprecision(17)<<"{\"passed\":true,\"source_only\":true,\"coupled_trace_qualified\":false,\"counts\":[64,32,16],\"seeds\":3,\"events\":"<<events.size()<<",\"prefix_comparisons\":"<<comparisons<<",\"physical_electrons\":"<<target<<",\"max_event_charge_relative_error\":"<<maximum_event_error<<",\"max_cumulative_charge_relative_error\":"<<maximum_total_error<<",\"trace_roundtrip_passed\":true,\"wrong_injection_time_rejected\":true,\"negative_dt_rejected\":true}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
