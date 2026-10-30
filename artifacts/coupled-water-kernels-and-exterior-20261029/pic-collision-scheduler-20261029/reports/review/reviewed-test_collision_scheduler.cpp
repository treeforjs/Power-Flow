#include "collision_scheduler.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace raikou::pic_collisions;
int checks=0;
void check(bool ok,const char* what){++checks;if(!ok)throw Error(what);}
template<class F>void rejects(F f,const char* what){bool failed=false;try{f();}catch(const Error&){failed=true;}check(failed,what);}
bool close(double a,double b,double tol=2e-13){return std::abs(a-b)<=tol*std::max({std::abs(a),std::abs(b),1e-300});}
int main(int argc,char**argv){try{
    require(argc==2,"results path");
    constexpr double me=9.1093837139e-31,mn=2.9915e-26;
    const std::string hash(64,'a'),kernel(64,'b');
    const double energy=8000,rest=me*c*c/qe;
    const double u=c*std::sqrt((energy/rest)*(energy/rest+2));
    Packet e{"e0","electron",3,2,me,{u,0,0}},n{"n0","H2O",3,8,mn,{0,0,0}};
    std::vector<Cell> cells{{3,1e-12}};
    Channel channel{"ionization","electron","H2O",hash,kernel,Observable::ExclusiveEvent,[](double E){if(E<0||E>10000)throw Error("unsupported energy");return 1e-20;}};
    const auto kin=relative(e,n);
    check(close(kin.target_frame_energy_eV,8000),"target-rest incident energy");
    check(close(kin.moller_flux_m_s,u/std::sqrt(1+u*u/(c*c))),"stationary target flux");
    Packet moving=n;moving.proper_velocity_m_s=e.proper_velocity_m_s;
    check(relative(e,moving).moller_flux_m_s==0,"comoving flux is zero");
    check(relative(e,moving).target_frame_energy_eV<1e-25,"comoving energy is zero");
    // Common Lorentz boost of a stationary target and initial electron.
    const double vb=.17*c,gb=1/std::sqrt(1-vb*vb/(c*c));
    Packet eb=e,nb=n;const double ge=std::sqrt(1+u*u/(c*c));
    eb.proper_velocity_m_s[0]=gb*(u+vb*ge);nb.proper_velocity_m_s[0]=gb*vb;
    check(close(relative(eb,nb).target_frame_energy_eV,energy),"target-frame energy invariant under common boost");
    auto rows=candidates({e},{n},cells,{channel});check(rows.size()==1,"one supported pair");
    check(rows[0].event_weight==2&&rows[0].projectile_remainder==0&&rows[0].target_remainder==6,"smaller packet event with target remainder");
    check(close(rows[0].rate_s*rows[0].event_weight,e.weight*n.weight*1e-20*kin.moller_flux_m_s/1e-12),"weighted physical reaction rate identity");
    auto large=e;large.weight=16;rows=candidates({large},{n},cells,{channel});
    check(rows[0].event_weight==8&&rows[0].projectile_remainder==8&&rows[0].target_remainder==0,"projectile remainder");
    check(close(rows[0].rate_s*rows[0].event_weight,large.weight*n.weight*1e-20*kin.moller_flux_m_s/1e-12),"reverse weight-rate identity");
    auto another=n;another.cell=4;auto two_cells=cells;two_cells.push_back({4,1e-12});
    check(candidates({e},{another},two_cells,{channel}).empty(),"different physical cells never collide");
    auto zero=channel;zero.sigma_m2=[](double){return 0.;};check(candidates({e},{n},cells,{zero}).empty(),"zero cross section");
    auto aggregate=channel;aggregate.observable=Observable::AggregateInelastic;
    rejects([&]{candidates({e},{n},cells,{aggregate});},"aggregate observable rejected");
    aggregate.observable=Observable::InclusiveProductYield;rejects([&]{candidates({e},{n},cells,{aggregate});},"inclusive yields rejected");
    auto absent=channel;absent.kernel_sha256.clear();rejects([&]{candidates({e},{n},cells,{absent});},"missing product law rejected");
    auto gap=channel;gap.sigma_m2=[](double)->double{throw Error("table support gap");};rejects([&]{candidates({e},{n},cells,{gap});},"unsupported energy not filled");
    auto negative=channel;negative.sigma_m2=[](double){return -1.;};rejects([&]{candidates({e},{n},cells,{negative});},"negative cross section rejected");
    auto duplicate=n;duplicate.id=e.id;rejects([&]{candidates({e},{duplicate},cells,{channel});},"duplicate packet identity rejected");
    rejects([&]{candidates({e},{n},cells,{channel},0);},"pair-work budget enforced");
    auto bad=e;bad.weight=0;rejects([&]{candidates({bad},{n},cells,{channel});},"zero weight rejected");
    Clock clock{1,0,0,100,hash,hash,hash};
    Scheduler owner(99301),foreign(99301);auto saved=owner.checkpoint();
    auto input=candidates({e},{n},cells,{channel});
    auto proposal=owner.prepare(clock,input);
    check(owner.checkpoint()==saved,"preparation does not mutate RNG");
    check(!foreign.can_publish(proposal),"foreign token rejected");
    check(owner.try_publish(proposal)&&!owner.try_publish(proposal),"publish exactly once");
    auto restarted=Scheduler::restore(owner.checkpoint());
    auto p1=owner.prepare(clock,input),p2=restarted.prepare(clock,input);
    check(p1.event_time_s()==p2.event_time_s()&&p1.event()->channel_id==p2.event()->channel_id,"restart exact event/time");
    check(owner.try_publish(p1)&&restarted.try_publish(p2)&&owner.checkpoint()==restarted.checkpoint(),"restart exact RNG continuation");
    rejects([&]{Scheduler::restore(saved+" garbage");},"checkpoint trailing data rejected");
    auto badclock=clock;badclock.mesh_sha256.clear();saved=owner.checkpoint();
    rejects([&]{owner.prepare(badclock,input);},"missing geometry binding rejected");
    check(owner.checkpoint()==saved,"failed proposal restores RNG");
    auto poisoned=input;poisoned[0].rate_s*=2;
    rejects([&]{owner.prepare(clock,poisoned);},"rate/weight mismatch rejected");
    const auto empty=owner.prepare(clock,{});check(!empty.event()&&owner.checkpoint()==saved,"empty population no event/no RNG change");
    // Controlled two-channel exponential clock; rates1/s and2/s exactly by construction.
    auto c1=channel,c2=channel;c1.id="one";c2.id="two";
    const double factor=1e-12/(n.weight*kin.moller_flux_m_s);
    c1.sigma_m2=[=](double){return factor;};c2.sigma_m2=[=](double){return 2*factor;};
    input=candidates({e},{n},cells,{c2,c1});
    Scheduler ensemble(71);long double mean=0,second=0;int count_one=0,cdf_count=0;
    constexpr int samples=50000;
    for(int i=0;i<samples;++i){const auto p=ensemble.prepare(clock,input);require(p.event().has_value(),"controlled window unexpected truncation");mean+=p.event_time_s();second+=p.event_time_s()*p.event_time_s();count_one+=p.event()->channel_id=="one";cdf_count+=p.event_time_s()<1./3;require(ensemble.try_publish(p),"ensemble commit");}
    const double observed=static_cast<double>(mean/samples),fraction=static_cast<double>(count_one)/samples,cdf=static_cast<double>(cdf_count)/samples;
    check(std::abs(observed-1./3)<6*(1./3)/std::sqrt(samples),"exponential mean six-sigma check");
    check(std::abs(static_cast<double>(second/samples)-2./9)<.012,"exponential second moment");
    check(std::abs(fraction-1./3)<6*std::sqrt((1./3)*(2./3)/samples),"branch frequencies six-sigma check");
    check(std::abs(cdf-(1-std::exp(-1.)))<6*std::sqrt(.25/samples),"waiting-time CDF");
    auto reverse=input;std::reverse(reverse.begin(),reverse.end());
    const auto x=ensemble.prepare(clock,input),y=ensemble.prepare(clock,reverse);
    check(x.event_time_s()==y.event_time_s()&&x.event()->channel_id==y.event()->channel_id,"input order does not alter sampled event");
    auto very_short=clock;very_short.end_s=1e-100;
    const auto noevent=ensemble.prepare(very_short,input);check(!noevent.event(),"finite interval can accept no event");
    std::ofstream out(argv[1]);out<<std::setprecision(17)<<"{\n  \"status\": \"passed\",\n  \"checks\": "<<checks<<",\n  \"samples\": "<<samples<<",\n  \"mean_wait_s\": "<<observed<<",\n  \"expected_mean_wait_s\": "<<1./3<<",\n  \"channel_one_fraction\": "<<fraction<<",\n  \"waiting_time_cdf_at_one_over_rate\": "<<cdf<<",\n  \"physical_H2O_event_data_bound\": false,\n  \"native_field_evolution\": false,\n  \"status_scope\": \"weighted finite-rate scheduler, no product or native mutation\"\n}\n";
    require(bool(out),"write result");std::cout<<checks<<" scheduler checks passed\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
