#include "chem_v8_fixture.hpp"
#include "chem_v3_native.hpp"
#include "census_charge.hpp"
#include "raikou/json.h"
#include <fstream>
#include <iostream>
namespace n=raikou::mfem_backend;
namespace h=hall_corrective_v3;
using J=raikou::Json;
pfc::ExactSum total(std::span<const n::ArenaRateEntry> rows){pfc::ExactSum q;for(const auto&r:rows)q.add_product(r.native.weight,r.native.charge_C);return q;}
int main(int argc,char**argv){try{
 if(argc!=2)throw std::runtime_error("output path required");
 const auto options=chemfix::options();auto cfg=chemfix::chemistry_configuration(options);
 // This diagnostic explicitly enables unowned target-energy events solely to
 // examine the actual candidate's spatial charge. No candidate is published.
 std::vector<n::FullTimeAuthoredPacket> packets;
 for(unsigned k=0;k<16;++k){const bool electron=k%2==0;const auto name="diagnostic/"+std::to_string(k);
  packets.push_back({name,electron?"pfc-v3/e":"pfc-v3/H2O+","manufactured/local-charge-audit",n::ArenaCount(1e15,0),
    {electron?.00405:.00415,.0041,.0001},electron?n::PicVector{187553.7262105,0,0}:n::PicVector{0,0,0}});
 }
 auto A=chemfix::make_arena(options,cfg,packets,1000);const auto base=A.arena->state_prefix(),rng=A.arena->full_time_chemistry().state_sha256;
 const auto before=A.arena->rate_population(1000);J result;
 {
  const auto boundary=A.arena->full_time_boundary();const double dt=1e-12;
  n::FullTimeStepBinding common{0,boundary->accepted_arena_revision,boundary->time_s,boundary->time_s+dt,dt,boundary->accepted_arena_sha256,boundary->field_sha256,boundary->mesh_sha256,boundary->material_sha256,chemfix::sha("diagnostic/local-charge-step")};
  auto outer=A.arena->begin_full_time_step(common,1000,64);
  n::ExposureBinding exposure{common,n::FullTimeSlot::before_transport,common.begin_s,.5*dt,0,common.accepted_field_sha256};
  auto scope=outer.first_exposure(exposure);auto token=scope.prepare_chemistry_v3();
  const auto receipt=token.candidate_chemistry_stage();const auto population=scope.candidate_population(token,1000,1u<<24);
  std::vector<n::ArenaRateEntry> after;
  for(const auto&p:population->packets)for(const auto&carrier:p.population.carriers)after.push_back({p.population.packet_id,p.population.state,carrier});
  std::sort(after.begin(),after.end(),[](const auto&a,const auto&b){return a.native.id<b.native.id;});
  const auto export_=token.candidate_admitted_export(true);std::size_t charged_debits=0;
  for(const auto&d:export_->target_debits)if(A.reg.species.at(d.species).charge_e)++charged_debits;
  axial::Mesh mesh(45,4,.009,1);auto q0=h::charge_of_census(mesh,options,before,1000),q1=h::charge_of_census(mesh,options,after,1000);
  long double max=0,l1=0;std::size_t nodes=0;
  for(std::size_t k=0;k<q0.net_C.size();++k){const auto delta=std::abs((q1.net_C[k]-q0.net_C[k]).wide());max=std::max(max,delta);l1+=delta;nodes+=delta>0;}
  auto global=total(before);global.sub(total(after));
  const bool local=h::same_net_charge(q0,q1);
  if(receipt->admitted_events==0||charged_debits==0||!global.is_zero()||local)throw std::runtime_error("fixture did not reproduce charged-target local charge change");
  result={{"passed",true},{"scope","actual native diagnostic candidate only; no accepted chemistry or Hall step"},{"admitted_events",receipt->admitted_events},{"charged_target_debits",charged_debits},{"global_exact_charge_residual_zero",global.is_zero()},{"hall_same_net_charge_accepts",local},{"max_node_charge_change_C",double(max)},{"L1_node_charge_change_C",double(l1)},{"changed_nodes",nodes},{"before_carriers",before.size()},{"after_carriers",after.size()},{"candidate_population_sha256",population->sha256},{"admitted_export_sha256",export_->sha256},{"source_identity",n::PicArena::source_identity()}};
 }
 if(A.arena->state_prefix()!=base||A.arena->full_time_chemistry().state_sha256!=rng)throw std::runtime_error("diagnostic changed accepted state/RNG");
 result["accepted_native_and_rng_unchanged"]=true;
 std::ofstream(argv[1])<<result.dump(2)<<'\n';std::cout<<result.dump()<<'\n';return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
