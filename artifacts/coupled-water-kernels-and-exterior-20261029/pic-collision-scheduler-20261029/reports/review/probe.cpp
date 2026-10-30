#include "reviewed-collision_scheduler.hpp"
#include <iomanip>
#include <iostream>
using namespace raikou::pic_collisions;
int main(){
 std::cout<<std::setprecision(17);constexpr double me=9.1093837139e-31,mn=2.9915e-26;const std::string H(64,'a'),K(64,'b');
 Packet e{"e","e",1,1,me,{1e7,0,0}},n{"n","H2O",1,1e20,mn,{0,0,0}};Cell cell{1,1e-12};Channel ch{"c","e","H2O",H,K,Observable::ExclusiveEvent,[](double){return 1e-20;}};
 auto a=candidates({e},{n},{cell},{ch});Clock clock{1,0,0,100,H,H,H};Scheduler s(1);auto p=s.prepare(clock,a);
 std::cout<<"large_remainder original="<<n.weight<<" remainder="<<a[0].target_remainder<<" actual_binary64_debit="<<n.weight-a[0].target_remainder<<" event_weight="<<a[0].event_weight<<" scheduler_admits="<<p.event().has_value()<<"\n";
 int queries=0;std::vector<Channel>channels;for(int i=0;i<128;++i){auto z=ch;z.sigma_m2=[&](double){++queries;return 1e-20;};z.id="ch"+std::to_string(i);channels.push_back(z);}auto rows=candidates({e},{n},{cell},channels,1);std::cout<<"one_pair_limit rows="<<rows.size()<<" channel_queries="<<queries<<"\n";
 for(double u:{1.,1e8,1e12,1e16,1e20,1e29,1e100,1e300}){auto b=e,t=n;b.proper_velocity_m_s=t.proper_velocity_m_s={u,u*.37,-u*.21};try{auto k=relative(b,t);std::cout<<"comoving u="<<u<<" E_eV="<<k.target_frame_energy_eV<<" flux="<<k.moller_flux_m_s<<"\n";}catch(const Error&x){std::cout<<"comoving u="<<u<<" REJECT="<<x.what()<<"\n";}}
 auto low=e;low.proper_velocity_m_s={1e-4,0,0};auto k=relative(low,n);double exact=me*1e-8/(std::sqrt(1+1e-8/(c*c))+1)/qe;std::cout<<"low_rest E="<<k.target_frame_energy_eV<<" reference="<<exact<<" relative_error="<<(k.target_frame_energy_eV/exact-1)<<"\n";
 auto rejected=[&](auto work){try{work();return false;}catch(const Error&){return true;}};
 auto huge=e;huge.weight=1e308;auto hn=n;hn.weight=1e308;auto hx=ch;hx.sigma_m2=[](double){return 1e308;};
 std::cout<<"propensity_overflow_rejected="<<rejected([&]{(void)candidates({huge},{hn},{Cell{1,std::numeric_limits<double>::denorm_min()}},{hx});})<<"\n";
 auto tiny=e;tiny.weight=std::numeric_limits<double>::denorm_min();auto tn=n;tn.weight=tiny.weight;auto tx=ch;tx.sigma_m2=[](double){return std::numeric_limits<double>::denorm_min();};
 std::cout<<"propensity_underflow_rejected="<<rejected([&]{(void)candidates({tiny},{tn},{Cell{1,1e308}},{tx});})<<"\n";
 auto late=clock;late.begin_s=1;late.end_s=std::nextafter(1.,2.);std::cout<<"unresolvable_event_time_rejected="<<rejected([&]{(void)s.prepare(late,a);})<<"\n";
 std::optional<Scheduler> life;life.emplace(99);auto stale=life->prepare(clock,a);life.reset();life.emplace(123);std::cout<<"reused_owner_storage_accepts_old_token="<<life->can_publish(stale)<<"\n";
 std::cout<<"pending_clock original_revision="<<p.clock().revision<<" generation_still_publishable="<<s.can_publish(p)<<" note=physical_owner_must_validate_current_clock\n";
}
