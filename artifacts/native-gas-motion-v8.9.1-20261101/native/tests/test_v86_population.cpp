// Actual native candidate inventory controls; diagnostic chemistry only.
#include "population_step_fixture.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
namespace n=raikou::mfem_backend;
using J=nlohmann::json;
namespace raikou::mfem_backend::full_time_detail {std::string candidate_population_digest(const FullTimeCandidatePopulation &);}
int checks=0;
void check(bool ok,const char*why){++checks;if(!ok)throw std::runtime_error(why);}
template<class F>void refuses(F f,const char*why){bool bad=false;try{f();}catch(const std::exception&){bad=true;}check(bad,why);}
n::FullTimeStepBinding binding(const chemfix::Arena&A){const auto b=*A.arena->full_time_boundary();const double dt=1e-12;return {b.accepted_step,b.accepted_arena_revision,b.time_s,b.time_s+dt,dt,b.accepted_arena_sha256,b.field_sha256,b.mesh_sha256,b.material_sha256,chemfix::sha("population-control-clock")};}
struct Attempt {
 n::FullTimeStepBinding b; n::PrivateFullTimeStep outer;n::CollisionScope scope;
 Attempt(chemfix::Arena&A):b(binding(A)),outer(A.arena->begin_full_time_step(b,120000,64)),scope(outer.first_exposure({b,n::FullTimeSlot::before_transport,b.begin_s,.5*b.dt_s,A.arena->full_time_chemistry().stages,b.accepted_field_sha256})){}
};
std::vector<n::GapParticle> rows(const n::FullTimeCandidatePopulation&v){std::vector<n::GapParticle> r;for(const auto&p:v.packets)r.insert(r.end(),p.population.carriers.begin(),p.population.carriers.end());std::sort(r.begin(),r.end(),[](const auto&a,const auto&b){return a.id<b.id;});return r;}
bool same(const n::GapParticle&a,const n::GapParticle&b){return a.id==b.id && a.weight==b.weight && a.mass_kg==b.mass_kg && a.charge_C==b.charge_C && a.position_m==b.position_m && a.proper_velocity_m_s==b.proper_velocity_m_s;}
void compare(const n::FullTimeCandidatePopulation&v,const n::ArenaSnapshot&s,std::span<const n::GapParticle> actual){
 auto r=rows(v);check(r.size()==actual.size() && std::equal(r.begin(),r.end(),actual.begin(),same),"actual candidate and published rows differ");
 check(v.packets.size()==s.packets.size(),"actual packet cardinality");
 for(const auto&p:v.packets){const auto i=s.packets.find(p.population.packet_id);check(i!=s.packets.end(),"actual packet identity");const auto&a=i->second;check(a.state==p.population.state && a.count==p.population.count && a.origin==p.origin && a.last_event==p.last_event,"actual packet count/state/provenance");std::vector<std::uint64_t> ids;for(const auto&r0:p.population.carriers)ids.push_back(r0.id);check(ids==a.native_ids,"actual packet carrier ownership");}
}
std::shared_ptr<const n::FullTimeCandidatePopulation> inspect(n::CollisionScope&s,n::PreparedArena&t){
 refuses([&]{(void)s.candidate_population(t,0,1u<<28);},"zero record cap refused before cache");
 refuses([&]{(void)s.candidate_population(t,120000,1);},"small payload budget refused before cache");
 auto v=s.candidate_population(t,120000,1u<<28);n::validate_full_time_candidate_population(*v);
 check(v->admitted_export_sha256==t.candidate_chemistry_stage()->admitted_export->sha256,"receipt export bound");
 check(v==s.candidate_population(t,120000,1u<<28),"same immutable view reused");
 check(v==s.candidate_population(t,std::max<std::uint64_t>(v->carriers,1),v->payload_bytes),"adequate changed cap retains view");
 refuses([&]{(void)s.candidate_population(t,120000,v->payload_bytes-1);},"payload guard on cache hit");
 if(v->carriers>1)refuses([&]{(void)s.candidate_population(t,v->carriers-1,1u<<28);},"record guard on cache hit");
 return v;
}
int main(int argc,char**argv){try{
 if(argc!=2)throw std::runtime_error("output path required");
 std::filesystem::path out=argv[1];std::filesystem::create_directories(out);
 auto o=chemfix::options();auto c=chemfix::chemistry_configuration(o);c.retain_stage_inputs=true;
 chemfix::Population p;p.carriers=3000;p.seed=11;p.clustered=true;p.water_weight=2e13;p.charged_weight=1e-10;p.electron_eV_min=5;p.electron_eV_max=300;
 auto packets=chemfix::population(p);auto A=chemfix::make_arena(o,c,packets,120000,2);
 const auto before=A.arena->state_prefix(),rng=A.arena->full_time_chemistry().state_sha256;
 std::shared_ptr<const n::FullTimeCandidatePopulation> saved;
 {Attempt a(A);auto t=a.scope.prepare_chemistry_v3();saved=inspect(a.scope,t);
  auto other=chemfix::make_arena(o,c,packets,120000,2);Attempt f(other);
  refuses([&]{(void)f.scope.candidate_population(t,120000,1u<<28);},"foreign equal owner token refused");
  auto moved=std::move(t);refuses([&]{(void)a.scope.candidate_population(t,120000,1u<<28);},"moved-from token refused");
  check(a.scope.candidate_population(moved,120000,1u<<28)==saved,"moved live token keeps snapshot");
  auto scope=std::move(a.scope);refuses([&]{(void)a.scope.candidate_population(moved,120000,1u<<28);},"moved-from scope refused");
  check(scope.candidate_population(moved,120000,1u<<28)==saved,"moved live scope preserves binding");
  check(scope.try_publish(moved),"private chemistry publish");
  refuses([&]{(void)scope.candidate_population(moved,120000,1u<<28);},"consumed token refused");
 }
 check(A.arena->state_prefix()==before && A.arena->full_time_chemistry().state_sha256==rng,"abandoned private stage preserves accepted state and RNG");
 {Attempt retry(A);auto t=retry.scope.prepare_chemistry_v3();auto v=inspect(retry.scope,t);check(v!=saved && v->sha256==saved->sha256,"same-content retry gets fresh immutable snapshot");}
 n::validate_full_time_candidate_population(*saved);
 {auto dead=chemfix::make_arena(o,c,packets,120000,2);Attempt a(dead);auto t=a.scope.prepare_chemistry_v3();auto v=inspect(a.scope,t);dead.arena.reset();refuses([&]{(void)a.scope.candidate_population(t,120000,1u<<28);},"destroyed native owner refused");n::validate_full_time_candidate_population(*v);}
 // Fully completed two-stage actual candidate snapshots, including compensated
 // target packet remainder carriers and release products.
 std::vector<std::shared_ptr<const n::FullTimeCandidatePopulation>> stages;
 auto R=chemfix::step_with_population(A,1e-12,0,[&](auto&s){auto t=s.prepare_chemistry_v3();stages.push_back(inspect(s,t));return t;},{},{},true,[&](auto&trial){compare(*stages.front(),*trial.snapshot(),trial.census_v7().carriers);});
 check(stages.size()==2 && R.stage1->admitted_events>0 && R.stage2->released_packets>0,"actual populated consumption and release controls");
 compare(*stages.back(),*A.arena->snapshot(),A.arena->full_time_population_view().rows);
 std::size_t multi=0;for(const auto&s:stages)for(const auto&x:s->packets)multi+=x.population.carriers.size()>1;
 check(multi>0,"actual nonprojected compensated multi-carrier count");
 for(const auto&s:stages)check(s->states.size()==A.states.size(),"complete state internal/reference catalogue");
 // Re-signed malformed facts cannot bypass independent structural validation.
 auto tamper=*stages.back();tamper.carriers=std::numeric_limits<std::uint64_t>::max();
 refuses([&]{n::validate_full_time_candidate_population(tamper);},"malformed cardinality refused before reserve");
 for(int k=0;k<4;++k){auto bad=*stages.back();auto&entry=bad.packets.front();
  if(k==0)entry.population.count=entry.population.count.credited(1);
  if(k==1)entry.population.carriers.front().weight*=2;
  if(k==2)entry.population.carriers.front().id=0;
  if(k==3)entry.population.carried_internal_energy_J=1e-19;
  bad.sha256=n::full_time_detail::candidate_population_digest(bad);
  refuses([&]{n::validate_full_time_candidate_population(bad);},"re-signed malformed packet refused");
 }
 // Empty stage and an actual authored birth with nonzero packet-owned energy.
 // Manufactured input exercises native ownership; not source physics qualification.
 auto empty=chemfix::make_arena(o,c,{},120000,2);stages.clear();
 auto E=chemfix::step_with_population(empty,1e-12,0,[&](auto&s){auto t=s.prepare_chemistry_v3();stages.push_back(inspect(s,t));return t;},{},{},true,{},2e-19);
 check(E.stage1->events==0 && E.stage2->events==0 && stages.front()->carriers==0,"actual no-event stages");
 check(stages.back()->packets.size()==1,"authored internal-energy birth persisted");const auto&ip=stages.back()->packets.front();
 check(ip.has_carried_internal_energy && ip.internal_birth_weight==2 && ip.population.carried_internal_energy_J==2e-19 && ip.population.state=="pfc-v3/H2O","actual carried internal owner and original operand retained");
 compare(*stages.back(),*empty.arena->snapshot(),empty.arena->full_time_population_view().rows);
 empty.arena->validate_full_time_history();A.arena->validate_full_time_history();
 J result={{"passed",true},{"checks",checks},{"source",n::PicArena::source_identity()},{"stage1_events",R.stage1->admitted_events},{"stage2_releases",R.stage2->released_packets},{"multi_carrier_packets",multi},{"production_qualified",false},{"device_view_qualified",false},{"internal_energy_fixture","manufactured authored native birth, not coupled source"}};
 std::ofstream(out/"result.json")<<result.dump(2)<<"\n";std::cout<<result.dump()<<"\n";return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
