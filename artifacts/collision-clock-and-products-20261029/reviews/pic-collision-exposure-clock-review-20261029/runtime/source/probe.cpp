#define raikou legacy_raikou
#include "../snapshot/inputs/collision_scheduler.hpp"
#undef raikou
#include "../snapshot/source/collision_scheduler.hpp"
#include <bit>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <type_traits>
namespace now=raikou::pic_collisions;
namespace old=legacy_raikou::pic_collisions;
std::size_t checks=0,rejected=0,legacy_pairs=0,exposure_pairs=0;
void check(bool yes,const char* why){++checks;if(!yes)throw std::runtime_error(why);}
bool bits(double a,double b){return std::bit_cast<std::uint64_t>(a)==std::bit_cast<std::uint64_t>(b);}
template<class C>C candidate(std::string id,double rate){
    C c;c.projectile_id="projectile";c.target_id="target-"+id;c.channel_id=id;
    c.data_sha256=std::string(64,'d');c.kernel_sha256=std::string(64,'e');c.cell=1;
    c.event_weight=1;c.projectile_before=1;c.target_before=1;c.rate_s=rate;
    c.energy_eV=100;c.flux_m_s=1;c.sigma_m2=rate;c.volume_m3=1;return c;
}
template<class C>C clock_at(double begin,double dt,double elapsed=0){
    C c;c.begin_s=begin;c.end_s=begin+dt;c.frozen_physical_stage_s=c.end_s;
    c.collision_elapsed_s=elapsed;c.state_sha256=std::string(64,'a');
    c.mesh_sha256=std::string(64,'b');c.mechanism_sha256=std::string(64,'c');c.step=42;c.revision=23;return c;
}
template<class P,class Q>void equal_choice(const P& a,const Q& b){
    check(a.event().has_value()==b.event().has_value(),"event parity");
    check(bits(a.total_rate_s(),b.total_rate_s()),"rate parity");
    check(bits(a.event_subflow_elapsed_s(),b.event_subflow_elapsed_s()),"wait parity");
    if(a.event())check(a.event()->channel_id==b.event()->channel_id&&bits(a.event()->rate_s,b.event()->rate_s),"chosen identity/rate parity");
}
template<class F>void refuses(now::Scheduler& s,F f,const char* label){
    const auto before=s.checkpoint();bool threw=false;try{f();}catch(const now::Error&){threw=true;}
    check(threw,label);check(s.checkpoint()==before,"refusal mutated RNG/generation");++rejected;
}
int main(int argc,char**argv){try{
    if(argc!=2)throw std::runtime_error("output required");
    const std::vector<now::Candidate> nc={candidate<now::Candidate>("z",2e14),candidate<now::Candidate>("a",3e14),candidate<now::Candidate>("m",7e14)};
    const std::vector<old::Candidate> oc={candidate<old::Candidate>("z",2e14),candidate<old::Candidate>("a",3e14),candidate<old::Candidate>("m",7e14)};
    // A physically shifted explicit exposure must equal the untouched legacy
    // clock at local origin, including full RNG state after every publication.
    for(std::uint64_t seed=0;seed<128;++seed){
        now::Scheduler a(seed);old::Scheduler b(seed);
        for(double epoch:{0.,7.49788055492e-10,2.64784435e-8,1.,0x1p40}){
            const double dt=epoch<1?5e-15:epoch<2?3e-16:0x1p-10;
            for(double exposure:{dt,dt*.5,std::nextafter(dt*.5,0.)}){
                for(double fraction:{0.,.375,.875}){
                    const double elapsed=exposure*fraction;
                    const auto ca=clock_at<now::Clock>(epoch,dt,elapsed);
                    const auto cb=clock_at<old::Clock>(0.,exposure,elapsed);
                    const auto x=a.prepare_exposure({ca,dt,exposure},nc);const auto y=b.prepare(cb,oc);
                    equal_choice(x,y);check(bits(x.exposure_duration_s(),exposure),"exposure not preserved");
                    check(bits(x.clock().begin_s,epoch)&&bits(x.clock().end_s,ca.end_s)&&bits(x.frozen_physical_stage_s(),ca.frozen_physical_stage_s),"physical timestamps changed");
                    check(x.clock().step==42&&x.clock().revision==23&&x.clock().state_sha256==ca.state_sha256&&x.clock().collision_elapsed_s==elapsed,"clock provenance changed");
                    check(a.try_publish(x)&&b.try_publish(y),"publish parity");
                    check(a.checkpoint()==b.checkpoint(),"full RNG/generation parity");++exposure_pairs;
                }
            }
        }
    }
    // Direct old API comparison at rounded nonzero epochs. Numeric body must
    // retain its old end-begin semantics rather than adopt a new authority.
    for(std::uint64_t seed=0;seed<256;++seed){
        now::Scheduler a(seed);old::Scheduler b(seed);
        for(double epoch:{0.,1e-9,2.64784435e-8,1.}){
            const double dt=5e-15;
            const auto x=a.prepare(clock_at<now::Clock>(epoch,dt),nc);const auto y=b.prepare(clock_at<old::Clock>(epoch,dt),oc);
            equal_choice(x,y);check(a.try_publish(x)&&b.try_publish(y)&&a.checkpoint()==b.checkpoint(),"legacy API RNG parity");++legacy_pairs;
        }
    }
    now::Scheduler s(911);auto good=now::ExposureClock{clock_at<now::Clock>(1.,3e-16),3e-16,1.5e-16};
    const double inf=INFINITY,nan=NAN,tiny=std::numeric_limits<double>::denorm_min(),max=std::numeric_limits<double>::max();
    for(double bad:{0.,-0.,-tiny,-1.,inf,-inf,nan}){
        auto q=good;q.common_dt_s=bad;refuses(s,[&]{s.prepare_exposure(q,nc);},"bad common duration admitted");
        q=good;q.exposure_s=bad;refuses(s,[&]{s.prepare_exposure(q,nc);},"bad exposure admitted");
    }
    for(double bad:{-1.,inf,-inf,nan}){
        auto q=good;q.stage.begin_s=bad;refuses(s,[&]{s.prepare_exposure(q,nc);},"bad begin admitted");
        q=good;q.stage.end_s=bad;refuses(s,[&]{s.prepare_exposure(q,nc);},"bad end admitted");
        q=good;q.stage.frozen_physical_stage_s=bad;refuses(s,[&]{s.prepare_exposure(q,nc);},"bad physical stage admitted");
        q=good;q.stage.collision_elapsed_s=bad;refuses(s,[&]{s.prepare_exposure(q,nc);},"bad elapsed admitted");
    }
    for(double endpoint:{std::nextafter(good.stage.end_s,0.),std::nextafter(good.stage.end_s,inf)}){
        auto q=good;q.stage.end_s=endpoint;refuses(s,[&]{s.prepare_exposure(q,nc);},"rounded wrong endpoint admitted");
    }
    {auto q=good;q.stage=clock_at<now::Clock>(max,max);q.common_dt_s=max;q.exposure_s=1;refuses(s,[&]{s.prepare_exposure(q,nc);},"overflow common endpoint admitted");}
    {auto q=good;q.stage=clock_at<now::Clock>(1.,tiny);q.common_dt_s=tiny;q.exposure_s=tiny;refuses(s,[&]{s.prepare_exposure(q,nc);},"collapsed common interval admitted");}
    {auto q=good;q.exposure_s=std::nextafter(q.common_dt_s,inf);refuses(s,[&]{s.prepare_exposure(q,nc);},"overlong exposure admitted");}
    for(double elapsed:{good.exposure_s,std::nextafter(good.exposure_s,inf)}){
        auto q=good;q.stage.collision_elapsed_s=elapsed;refuses(s,[&]{s.prepare_exposure(q,nc);},"elapsed end admitted");}
    {auto q=good;q.stage.frozen_physical_stage_s=std::nextafter(q.stage.begin_s,0.);refuses(s,[&]{s.prepare_exposure(q,nc);},"before stage admitted");}
    {auto q=good;q.stage.frozen_physical_stage_s=std::nextafter(q.stage.end_s,inf);refuses(s,[&]{s.prepare_exposure(q,nc);},"after stage admitted");}
    {auto q=good;q.stage.state_sha256="bad";refuses(s,[&]{s.prepare_exposure(q,nc);},"unbound state admitted");}
    refuses(s,[&]{s.prepare_exposure(good,nc,2);},"candidate size guard bypassed");
    {auto bad=nc;bad.back().rate_s=nan;refuses(s,[&]{s.prepare_exposure(good,bad);},"invalid candidate bypassed");}
    // Tiny local exposure can remain valid at a representable common interval;
    // it is not an absolute physical timestamp or a source of clock advance.
    {auto q=good;q.exposure_s=tiny;const auto before=s.checkpoint();const auto p=s.prepare_exposure(q,nc);check(!p.event()&&bits(p.exposure_duration_s(),tiny),"subnormal exposure lost");check(s.checkpoint()==before,"tiny prepare mutation");check(s.try_publish(p),"tiny publish failed");}
    // Failed proposal construction cannot invalidate an earlier valid handle.
    const auto handle=s.prepare_exposure(good,nc);auto invalid=good;invalid.common_dt_s=nan;
    refuses(s,[&]{s.prepare_exposure(invalid,nc);},"bad proposal accepted");
    check(s.can_publish(handle)&&s.try_publish(handle)&&!s.try_publish(handle),"invalid prepare/stale handle semantics");
    auto restored=now::Scheduler::restore(s.checkpoint());
    for(int n=0;n<8;++n){const auto a=s.prepare_exposure(good,nc),b=restored.prepare_exposure(good,nc);equal_choice(a,b);check(s.try_publish(a)&&restored.try_publish(b)&&s.checkpoint()==restored.checkpoint(),"restart parity");}
    // A deterministic waiting-time threshold demonstrates both directions of
    // timestamp subtraction error without using the implementation's tests.
    std::size_t boundary_differences=0;
    for(double dt:{3e-16,4e-16}){
        const double begin=1.,rounded=(begin+dt)-begin,exposure=dt/2,wrong=rounded/2;
        std::mt19937_64 rng(147);const double u=(static_cast<double>(rng()>>12)+.5)*0x1p-52;
        const long double desired=(static_cast<long double>(exposure)+wrong)/2;
        const double rate=static_cast<double>(-std::log(u)/desired);
        const std::vector<now::Candidate> pair={candidate<now::Candidate>("boundary",rate)};
        now::Scheduler a(147),b(147);auto ca=clock_at<now::Clock>(begin,dt);
        const auto x=a.prepare_exposure({ca,dt,exposure},pair),y=b.prepare_exposure({ca,rounded,wrong},pair);
        check(bool(x.event())==(desired<exposure)&&bool(y.event())==(desired<wrong),"independent boundary decision");
        check(bool(x.event())!=bool(y.event()),"boundary durations indistinguishable");++boundary_differences;
    }
    std::ofstream out(argv[1]);out<<"{\"status\":\"passed\",\"assertions\":"<<checks<<",\"legacy_comparisons\":"<<legacy_pairs<<",\"shifted_exposure_comparisons\":"<<exposure_pairs<<",\"invalid_state_refusals\":"<<rejected<<",\"rounding_direction_boundary_pairs\":"<<boundary_differences<<",\"physical_time_advanced\":false}\n";
    std::cout<<checks<<" assertions passed\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
