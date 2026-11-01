#include "population_step_fixture.hpp"
#include "chem_v3_native.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
namespace n=raikou::mfem_backend;
using J=nlohmann::json;
int checks=0;
void check(bool b,const char*s){++checks;if(!b)throw std::runtime_error(s);}
template<class F>void refuses(F f,const char*s){bool b=false;try{f();}catch(const std::exception&){b=true;}check(b,s);}
void exact(const n::FullTimeExactValue&x,const pfc::ExactSum&actual){check(x.digits==n::chem8::exact_digits(actual),"untouched exact stored digits");check(x.value==n::chem8::exact_value(actual).value,"stored diagnostic rounded value");}
void compare(const n::FullTimeChemistryStockFrame&x,const n::ArenaSnapshot&snap,const pfc::Blob&blob){
 const auto&st=*snap.full_time_chemistry;check(x.chemistry_state_sha256==n::chem8::state_digest(st) && x.stages==st.stages && x.pending==st.pending,"actual stored chemistry root/phase");
 check(x.reservoirs.size()==st.reservoirs.size(),"actual reservoir key cardinality");std::size_t i=0;
 for(const auto&[key,r]:st.reservoirs){const auto&v=x.reservoirs[i++];check(v.atlas_cell==key.cell && v.species==key.species && v.state=="pfc-v3/"+blob.species_names[key.species],"actual reservoir key/species");
  check(v.rest_energy_eV==blob.view.species[key.species].mass_eV && v.reference_eV==blob.view.species[key.species].ref_eV,"pinned C3 conversion operands");
  exact(v.count,r.count);exact(v.kinetic_eV,r.kinetic);for(int a=0;a<3;++a)exact(v.momentum_eV_c[a],r.momentum[a]);
  check(v.parcels.size()==r.parcels.size(),"actual persistent parcels count");pfc::ExactSum C,K;std::array<pfc::ExactSum,3>P;
  for(std::size_t k=0;k<v.parcels.size();++k){const auto&p=v.parcels[k];const auto&q=r.parcels[k];check(p.weight==q.weight && std::equal(p.position_m.begin(),p.position_m.end(),q.x) && std::equal(p.proper_velocity_m_s.begin(),p.proper_velocity_m_s.end(),q.u),"actual unchanged parcel bytes");
   const auto u=pfc::v3(p.proper_velocity_m_s[0],p.proper_velocity_m_s[1],p.proper_velocity_m_s[2]);const auto mom=pfc::p_from_u(v.rest_energy_eV,u);C.add(p.weight);K.add_product(p.weight,pfc::ke_from_u(v.rest_energy_eV,u));P[0].add_product(p.weight,mom.x);P[1].add_product(p.weight,mom.y);P[2].add_product(p.weight,mom.z);
  }
  exact(v.count,C);exact(v.kinetic_eV,K);for(int a=0;a<3;++a)exact(v.momentum_eV_c[a],P[a]);
 }
 for(int k=0;k<6;++k)exact(x.cumulative_internal_sink_eV[k],st.cumulative.c3.internal[k]);
 check(x.cumulative_gas_transfer.size()==st.cumulative.gas_heating.size(),"actual cumulative transfer cardinality");i=0;
 for(const auto&[key,e]:st.cumulative.gas_heating){const auto&t=x.cumulative_gas_transfer[i++];check(t.atlas_cell==key.cell && t.species==key.species,"actual transfer key");exact(t.energy_eV,e);}
 if(snap.full_time_v81){const auto&f=*snap.full_time_v81;check(x.packet_internal.size()==f.internal.size(),"packet energy owner cardinality");i=0;for(const auto&[id,e]:f.internal){const auto&p=x.packet_internal[i++];check(p.packet_id==id && p.energy_J==e.energy_J && p.original_birth_weight==e.weight,"actual packet internal energy");}
  exact(x.cumulative_internal_born_J,f.internal_born_J);exact(x.cumulative_internal_sunk_J,f.internal_sunk_J);exact(x.cumulative_internal_released_J,f.internal_released_J);
 }else check(x.packet_internal.empty(),"no manufactured packet internal stock");
}
std::uint64_t parcels(const n::FullTimeChemistryStockFrame&f){std::uint64_t count=0;for(const auto&r:f.reservoirs)count+=r.parcels.size();return count;}
int main(int argc,char**argv){try{
 if(argc!=2)throw std::runtime_error("output path required");
 std::filesystem::path out=argv[1];std::filesystem::create_directories(out);const auto blob=pfc::load_blob(chemfix::blob_path(),chemfix::blob_sha());
 const auto o=chemfix::options();auto cfg=chemfix::chemistry_configuration(o);chemfix::Population p;p.carriers=3000;p.seed=11;p.clustered=true;p.water_weight=2e13;p.charged_weight=1e-10;p.electron_eV_min=5;p.electron_eV_max=300;
 auto A=chemfix::make_arena(o,cfg,chemfix::population(p),120000,2);std::vector<std::shared_ptr<const n::FullTimeCandidateStocks>> views;
 auto capture=[&](auto&scope){auto token=scope.prepare_chemistry_v3();
  refuses([&]{(void)scope.candidate_stocks(token,0,1u<<28);},"zero stock capacity refusal");
  refuses([&]{(void)scope.candidate_stocks(token,1000000,1);},"stock payload refusal before allocate");
  auto v=scope.candidate_stocks(token,1000000,1u<<28);auto population=scope.candidate_population(token,120000,1u<<28);
  check(v->admitted_export_sha256==population->admitted_export_sha256 && v->candidate_state_prefix==population->candidate_state_prefix && v->before_state_prefix==population->before_state_prefix,"same actual candidate stock/population binding");
  check(v==scope.candidate_stocks(token,1000000,1u<<28),"immutable stock cache reused");
  check(v==scope.candidate_stocks(token,std::max<std::uint64_t>(1,v->entries_and_parcels),v->payload_bound_bytes),"sufficient changed cache capacities");
  refuses([&]{(void)scope.candidate_stocks(token,1000000,v->payload_bound_bytes-1);},"cached stock payload refusal");
  if(v->entries_and_parcels>1)refuses([&]{(void)scope.candidate_stocks(token,v->entries_and_parcels-1,1u<<28);},"cached stock cardinality refusal");
  compare(v->before,*scope.snapshot(),blob);views.push_back(v);
  auto moved=std::move(token);refuses([&]{(void)scope.candidate_stocks(token,1000000,1u<<28);},"moved-from stock token refused");
  check(scope.candidate_stocks(moved,1000000,1u<<28)==v,"moved live stock token retains cache");return moved;};
 auto result=chemfix::step_with_population(A,1e-12,0,capture,{},{},true,[&](auto&trial){compare(views.front()->after,*trial.snapshot(),blob);});
 check(views.size()==2 && parcels(views.front()->after)>0,"actual stage1 pending reservoir population");
 check(parcels(views.back()->before)==parcels(views.front()->after) && parcels(views.back()->after)==0 && result.stage2->released_packets>0,"actual stage2 drains stored parcels");
 compare(views.back()->after,*A.arena->snapshot(),blob);
 std::uint64_t pending=parcels(views.front()->after);bool internal=false;for(const auto&v:views.front()->after.cumulative_internal_sink_eV)internal|=v.value!=0;
 // Actual private rollback leaves both the stored reservoir and RNG unchanged.
 const auto base=A.arena->state_prefix(),rng=A.arena->full_time_chemistry().state_sha256;views.clear();
 (void)chemfix::step_with_population(A,1e-12,1,capture,{},{},false);
 check(A.arena->state_prefix()==base && A.arena->full_time_chemistry().state_sha256==rng,"discarded stock export leaves accepted stock/RNG unchanged");
 const auto retry_base=views.front();views.clear();
 (void)chemfix::step_with_population(A,1e-12,1,capture,{},{},false);
 check(views.front()!=retry_base && views.front()->sha256==retry_base->sha256,"identical retry has separate immutable stock identity");
 // Authored native neutral birth carries stored U separately from cumulative sinks.
 auto empty=chemfix::make_arena(o,cfg,{},120000,2);views.clear();
 (void)chemfix::step_with_population(empty,1e-12,0,capture,{},{},true,{},2e-19);
 compare(views.back()->after,*empty.arena->snapshot(),blob);
 check(views.back()->after.packet_internal.size()==1 && views.back()->after.packet_internal.front().energy_J==2e-19,"actual packet U owner retained");
 check(views.back()->after.cumulative_internal_born_J.value==2e-19 && views.back()->after.cumulative_internal_released_J.value==0,"stored U separated from cumulative born/released ledgers");
 J r={{"passed",true},{"checks",checks},{"source",n::PicArena::source_identity()},{"actual_stage1_pending_parcels",pending},{"stage2_released_packets",result.stage2->released_packets},{"nonzero_cumulative_internal_fixture",internal},{"strict_energy_owner_gate_open",false},{"device_view_qualified",false}};
 std::ofstream(out/"result.json")<<r.dump(2)<<"\n";std::cout<<r.dump()<<"\n";return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
