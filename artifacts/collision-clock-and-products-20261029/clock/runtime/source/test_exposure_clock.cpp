#include "collision_scheduler.hpp"
#include <fstream>
#include <iostream>
using namespace raikou::pic_collisions;
std::size_t checks=0;
void check(bool yes,const char* name){++checks;if(!yes)throw std::runtime_error(name);}
template<class F>void reject(F&& f,const char* name){bool failed=false;try{f();}catch(const Error&){failed=true;}check(failed,name);}
Clock physical(double begin,double dt,double stage){return {5,9,begin,begin+dt,std::string(64,'a'),std::string(64,'b'),std::string(64,'c'),stage,0};}
std::vector<Candidate> one(double rate){
    Packet a{"e","electron",7,1,9.1093837139e-31,{1000000,0,0}};
    Packet b{"n","H2O",7,1,2.9915e-26,{0,0,0}};
    Channel c0{"parent","electron","H2O",std::string(64,'d'),std::string(64,'e'),Observable::ExclusiveEvent,[rate](double){return rate/1000000.;}};
    // Divide by the exact computed relativistic flux so the requested total
    // propensity remains the independently chosen boundary-test operand.
    const auto k=relative(a,b);c0.sigma_m2=[rate,k](double){return rate/k.moller_flux_m_s;};
    return candidates({a},{b},{{7,1.}},{c0});
}
int main(int argc,char**argv){try{
    if(argc!=2)throw std::runtime_error("output JSON required");
    const auto pair=one(1e14);
    // Same physical/exposure duration recovers every existing random decision.
    for(std::uint64_t seed=0;seed<256;++seed){
        Scheduler a(seed),b(seed);auto c=physical(0,5e-15,0);
        for(int n=0;n<3;++n){auto x=a.prepare(c,pair);auto y=b.prepare_exposure({c,5e-15,5e-15},pair);
            check(bool(x.event())==bool(y.event()),"legacy decision parity");
            check(x.event_subflow_elapsed_s()==y.event_subflow_elapsed_s(),"legacy waiting time parity");
            check(x.total_rate_s()==y.total_rate_s(),"legacy rate parity");
            check(x.exposure_duration_s()==y.exposure_duration_s(),"legacy duration parity");
            check(a.try_publish(x)&&b.try_publish(y)&&a.checkpoint()==b.checkpoint(),"legacy RNG parity");}
    }
    // A common interval far from zero deliberately loses interval-length bits.
    // The explicitly supplied local half-exposure must never be reconstructed.
    const double begin=1.,dt=3e-16,half=dt*.5;
    const auto c=physical(begin,dt,begin+dt);
    check(c.end_s-c.begin_s!=dt,"test actually distinguishes rounded absolute difference");
    std::mt19937_64 random(123);const double u=(static_cast<double>(random()>>12)+.5)*0x1p-52;
    const long double logarithm=-std::log(u);
    const double rounded_half=(c.end_s-c.begin_s)*.5;
    const double wait=(half+rounded_half)*.5;
    const auto boundary_pair=one(static_cast<double>(logarithm/wait));
    Scheduler exact(123),wrong_duration(123);
    auto actual=exact.prepare_exposure({c,dt,half},boundary_pair);
    auto rounded=wrong_duration.prepare_exposure({c,c.end_s-c.begin_s,rounded_half},boundary_pair);
    check(bool(actual.event())==(wait<half)&&bool(rounded.event())==(wait<rounded_half)&&bool(actual.event())!=bool(rounded.event()),"explicit duration changes boundary event correctly");
    check(actual.frozen_physical_stage_s()==c.end_s,"true endpoint physical stage retained");
    check(actual.exposure_duration_s()==half,"exact local exposure retained");
    // The same local problem is independent of the absolute physical origin.
    for(double epoch:{0.,7.49788055492e-10,2.647868472199812e-8,1.}){
        const double step=5e-15;const auto local=physical(0,step,step),global=physical(epoch,step,epoch+step);
        Scheduler a(456),b(456);
        for(int n=0;n<4;++n){auto x=a.prepare_exposure({local,step,step*.5},pair);
            auto y=b.prepare_exposure({global,step,step*.5},pair);
            check(bool(x.event())==bool(y.event())&&x.event_subflow_elapsed_s()==y.event_subflow_elapsed_s(),"exposure invariant under physical epoch");
            check(y.frozen_physical_stage_s()==epoch+step,"endpoint is not a local clock facade");
            check(a.try_publish(x)&&b.try_publish(y)&&a.checkpoint()==b.checkpoint(),"epoch invariant both random decisions");}
    }
    Scheduler s(789);auto before=s.checkpoint();auto good=physical(7.49788055492e-10,5e-15,7.49788055492e-10);
    ExposureClock q{good,5e-15,2.5e-15};
    for(double bad:{0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}){
        auto a=q;a.common_dt_s=bad;reject([&]{s.prepare_exposure(a,pair);},"invalid common duration");
        auto b=q;b.exposure_s=bad;reject([&]{s.prepare_exposure(b,pair);},"invalid exposure duration");}
    {auto a=q;a.exposure_s=1e-14;reject([&]{s.prepare_exposure(a,pair);},"exposure longer than common interval");}
    {auto a=q;a.stage.end_s=std::nextafter(a.stage.end_s,INFINITY);reject([&]{s.prepare_exposure(a,pair);},"endpoint arithmetic mismatch");}
    {auto a=q;a.stage.frozen_physical_stage_s=std::nextafter(a.stage.end_s,INFINITY);reject([&]{s.prepare_exposure(a,pair);},"physical stage outside interval");}
    {auto a=q;a.stage.collision_elapsed_s=a.exposure_s;reject([&]{s.prepare_exposure(a,pair);},"elapsed at exposure end");}
    {auto a=q;a.stage.collision_elapsed_s=-1;reject([&]{s.prepare_exposure(a,pair);},"negative elapsed");}
    {auto a=q;a.stage.end_s=a.stage.begin_s;a.common_dt_s=1e-300;a.exposure_s=1e-301;reject([&]{s.prepare_exposure(a,pair);},"physical clock cannot advance");}
    check(s.checkpoint()==before,"all rejected prepares retain random state");
    auto p=s.prepare_exposure(q,pair);check(s.checkpoint()==before,"successful prepare retains random state");
    check(s.try_publish(p),"accepted exposure publishes");check(!s.try_publish(p),"double publication refused");
    auto restored=Scheduler::restore(s.checkpoint());auto left=s.prepare_exposure(q,pair),right=restored.prepare_exposure(q,pair);
    check(bool(left.event())==bool(right.event())&&left.event_subflow_elapsed_s()==right.event_subflow_elapsed_s(),"exposure restart next event exact");
    check(s.try_publish(left)&&restored.try_publish(right)&&s.checkpoint()==restored.checkpoint(),"exposure restart next RNG exact");
    std::ofstream out(argv[1]);out<<"{\"status\":\"passed\",\"checks\":"<<checks<<",\"physical_trajectory_advanced\":false,\"explicit_duration_boundary_test\":true,\"legacy_rng_parity\":true}\n";
    std::cout<<"Passed "<<checks<<" explicit exposure-clock checks\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
