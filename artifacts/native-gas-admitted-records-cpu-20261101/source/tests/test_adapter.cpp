// Actual native export / gas owner component check. Synthetic populated input;
// diagnostic chemistry policy, no field/heating feedback or production claim.
#include "chem_v8_fixture.hpp"
#include "native_gas_stage.hpp"
#include "staged_fixture.hpp"
#include <fstream>
#include <optional>
#include <iostream>
namespace n=raikou::mfem_backend;
namespace g=pf_gas_v1;
namespace a=native_gas_stage;
using J=nlohmann::json;
int checks=0;
void check(bool v,const char *why){++checks;if(!v)throw std::runtime_error(why);}
template<class F> void refuses(F f,const char *why){bool failed=false;try{f();}catch(const std::exception&){failed=true;}check(failed,why);}
g::Configuration config(const chemfix::Arena &A,const pfc::Blob &b){
  g::Configuration c;c.preset="diagnostic-native-admitted-records";
  auto grid=n::chem8::make_grid(A.o,A.cfg.chemistry.grid_counts);
  c.grid={grid.low,grid.high,grid.n,grid.volume};
  c.chemistry_blob_sha256=b.sha256;c.chemistry_registry_sha256=A.reg.registry_sha256;
  for(const auto &s:A.reg.species)c.species.push_back({s.state,s.mass_kg,s.charge_e,s.target});
  for(std::uint32_t k=0;k<b.n_proc;++k){const auto &p=b.view.procs[k];
    c.processes.push_back({b.view.groups[p.group].projectile,p.n_out,a::internal_carrier(p.kind,p.n_out),p.kind});}
  return c;
}
n::FullTimeStepBinding binding(const chemfix::Arena &A,double dt){
  const auto x=*A.arena->full_time_boundary();
  return {x.accepted_step,x.accepted_arena_revision,x.time_s,x.time_s+dt,dt,x.accepted_arena_sha256,
    x.field_sha256,x.mesh_sha256,x.material_sha256,chemfix::sha("adapter-component-clock")};
}
int main(int argc,char **argv){try{
  if(argc!=2)throw std::runtime_error("usage: test_adapter RESULT.json");
  auto blob=pfc::load_blob(chemfix::blob_path(),chemfix::blob_sha());
  auto options=chemfix::options();auto chemistry=chemfix::chemistry_configuration(options);
  chemfix::Population p;p.carriers=3000;p.seed=11;p.clustered=true;
  p.electron_fraction=.49;p.ion_fraction=.49;p.water_weight=2e13;p.charged_weight=1e-10;p.electron_eV_min=5;p.electron_eV_max=300;
  auto packets=chemfix::population(p);auto A=chemfix::make_arena(options,chemistry,packets,120000,2);
  auto cfg=config(A,blob);g::Owner gas(cfg);a::Binding bindings(gas,A.reg,blob,options);
  const auto accepted=gas.accepted_sha256(),native=A.arena->full_time_boundary()->accepted_arena_sha256;
  const auto rng=A.arena->full_time_chemistry().state_sha256;
  J counts;std::size_t nr=0,nf=0;
  {
    auto b=binding(A,1e-12);gas.begin_step(b.accepted_step);
    const auto unchanged=gas.candidate_sha256();
    auto outer=A.arena->begin_full_time_step(b,120000,64);
    n::ExposureBinding exposure{b,n::FullTimeSlot::before_transport,b.begin_s,.5*b.dt_s,
      A.arena->full_time_chemistry().stages,b.accepted_field_sha256};
    auto scope=outer.first_exposure(exposure);std::optional<a::Batch> retained;
    {
      auto token=scope.prepare_chemistry_v3();
      retained.emplace(a::prepare(scope,token,gas,bindings,120000));
      const auto &x=retained->source();nr=retained->recoils().size();nf=retained->flows().size();
      check(nr>0 && nf>0,"actual populated records and target debits");
      check(nr==x.events.size(),"only admitted event count");
      std::size_t gas_returns=0,consumed=0,partial=0,multi=0;
      for(std::size_t i=0;i<nr;++i){const auto &r=retained->recoils()[i];const auto &e=x.events[i];
        check(std::memcmp(&r.ev,e.event_bytes.data(),sizeof r.ev)==0,"admitted bytes preserved");
        check(r.grid_cell==e.atlas_cell && r.stage==1 && r.accepted_step==b.accepted_step,"event stage/atlas identity");
        bool returns=false;for(std::uint32_t k=0;k<r.ev.n_out;++k)returns|=r.ev.out_dest[k]==pfg::kDestGas;
        gas_returns+=returns;consumed+=e.consumes_tracked_target;
      }
      std::size_t idx=0;
      for(const auto &d:x.target_debits){partial+=!d.remainder.packet_id.empty();multi+=d.remainder.carriers.size()>1;
        for(const auto *snap:{&d.before,&d.remainder})for(const auto &row:snap->carriers){
          const auto &f=retained->flows()[idx++];
          check(f.native_id==row.id && a::same(f.weight,row.weight) && f.grid_cell==d.atlas_cell &&
            f.species==d.species && f.kind==pfg::kFlowDebit && f.negative==(snap==&d.before?1u:0u),"actual debit operands retained");
          check(std::memcmp(f.u,row.proper_velocity_m_s.data(),sizeof f.u)==0,"actual proper velocity retained");
        }
      }
      check(idx==nf && x.releases.empty(),"stage1 has no release flows");
      check(gas_returns>0 && consumed>0 && partial>0 && multi>0,"GAS return and partial consuming target cases exercised");
      retained->validate_current(scope,token,gas);
      refuses([&]{(void)a::prepare(scope,token,gas,bindings,nr+nf-1);},"combined record capacity enforced");
      check(gas.candidate_sha256()==unchanged,"adapter capacity refusal leaves gas unchanged");
      auto foreign=chemfix::make_arena(options,chemistry,packets,120000,2);
      auto foreign_outer=foreign.arena->begin_full_time_step(b,120000,64);
      auto foreign_scope=foreign_outer.first_exposure(exposure);
      refuses([&]{retained->validate_current(foreign_scope,token,gas);},"foreign native scope refused");
      g::Owner foreign_gas(cfg);foreign_gas.begin_step(b.accepted_step);
      check(foreign_gas.candidate_sha256()==gas.candidate_sha256(),"equal state foreign gas control");
      refuses([&]{retained->validate_current(scope,token,foreign_gas);},"equal-state foreign gas owner refused");
      gas.abandon();gas.begin_step(b.accepted_step);
      check(gas.candidate_sha256()==unchanged,"gas-only retry same digest control");
      refuses([&]{retained->validate_current(scope,token,gas);},"gas-only same-digest retry refused");
      retained.emplace(a::prepare(scope,token,gas,bindings,120000));
      gas.prepare_keys({}, {}, {});
      check(gas.candidate_sha256()==unchanged,"same-state key preparation control");
      refuses([&]{retained->validate_current(scope,token,gas);},"same-digest candidate operation invalidates binding");
      retained.emplace(a::prepare(scope,token,gas,bindings,120000));
      g::Owner wrong_step(cfg);wrong_step.begin_step(1);
      refuses([&]{(void)a::prepare(scope,token,wrong_step,bindings,120000);},"gas/native step mismatch refused before returning records");
      auto moved=std::move(token);
      refuses([&]{retained->validate_current(scope,token,gas);},"moved token refused");
      retained->validate_current(scope,moved,gas);
      counts={{"events",nr},{"flows",nf},{"gas_returns",gas_returns},{"consuming_events",consumed},
        {"partial_target_debits",partial},{"multicarrier_remainders",multi},{"export",x.sha256}};
    }
    auto retry=scope.prepare_chemistry_v3();
    check(retry.candidate_admitted_export()->sha256==retained->source().sha256,"deterministic retry same digest");
    refuses([&]{retained->validate_current(scope,retry,gas);},"stale same-digest native attempt refused");
    auto fresh=a::prepare(scope,retry,gas,bindings,120000);
    fresh.validate_current(scope,retry,gas);
    const auto taken=gas.take(b.accepted_step,1,fresh.recoils(),fresh.flows());
    check(taken!=unchanged,"actual gas stage consumes adapter records");
    refuses([&]{fresh.validate_current(scope,retry,gas);},"gas candidate mutation invalidates batch");
    gas.abandon();
    refuses([&]{fresh.validate_current(scope,retry,gas);},"closed gas candidate refused");
    // The native outer step is deliberately abandoned, never committed alone.
  }
  check(gas.accepted_sha256()==accepted && !gas.step_open(),"gas accepted state retained after abandon");
  check(A.arena->full_time_boundary()->accepted_arena_sha256==native && A.arena->full_time_chemistry().state_sha256==rng,
    "native accepted population/RNG retained after abandon");
  for(int k=0;k<7;++k){auto bad=cfg;
    if(k==0)bad.chemistry_blob_sha256=chemfix::sha("different-blob");
    if(k==1)bad.chemistry_registry_sha256=chemfix::sha("different-registry");
    if(k==2)bad.species[0].mass_kg=std::nextafter(bad.species[0].mass_kg,INFINITY);
    if(k==3)bad.species[0].charge=0;
    if(k==4)bad.processes[0].projectile=(bad.processes[0].projectile+1)%bad.species.size();
    if(k==5)bad.processes[0].internal_carrier=1;
    if(k==6)bad.processes.pop_back();
    refuses([&]{a::configuration(bad,A.reg,blob);},"different blob/process/species operands refused");
  }
  auto shifted=options;shifted.lower_foil.top_m+=1e-4;
  refuses([&]{a::Binding wrong(gas,A.reg,blob,shifted);},"different solid collision geometry refused");
  // Run both real native stages through a staged zero-field motion fixture.
  // A small diagnostic release weight exercises actual reservoir releases.
  auto release_chem=chemistry;release_chem.release_weight=1;
  auto release_pop=p;release_pop.charged_weight=1;
  auto S=chemfix::make_arena(options,release_chem,chemfix::population(release_pop),120000,2);
  g::Owner sg(config(S,blob));a::Binding sb(sg,S.reg,blob,options);
  const auto saccepted=sg.accepted_sha256(), snative=S.arena->full_time_boundary()->accepted_arena_sha256;
  const auto srng=S.arena->full_time_chemistry().state_sha256;
  sg.begin_step(S.arena->full_time_boundary()->accepted_step);
  std::size_t release_packets=0, release_rows=0, release_internal=0, stage2_admitted=0, refused_stage=0;
  adapter_fixture::staged_step(S,1e-12,0,[&](n::CollisionScope &s,n::PreparedArena &t){
    auto x=t.candidate_admitted_export(false);bool unsupported=false;
    for(const auto &d:x->target_debits)unsupported|=d.also_projectile||d.internal_energy_deferred;
    for(const auto &r:x->releases){++release_packets;release_rows+=r.carriers.size();
      release_internal+=r.carried_internal_energy_J!=0;unsupported|=r.carried_internal_energy_J!=0;}
    const auto unchanged=sg.candidate_sha256();
    if(unsupported){refuses([&]{(void)a::prepare(s,t,sg,sb,120000);},"actual unsupported internal ownership refused");
      check(sg.candidate_sha256()==unchanged && s.can_publish(t),"unsupported export does not mutate either candidate");++refused_stage;return;}
    auto batch=a::prepare(s,t,sg,sb,120000);batch.validate_current(s,t,sg);
    if(x->stage==2){++stage2_admitted;std::size_t count=0;
      for(const auto &f:batch.flows())if(f.kind==pfg::kFlowRelease){++count;check(!f.negative,"release flow positive");}
      check(count==release_rows,"all actual release carrier operands exported");
      std::vector<const pfg::FlowRecord *> actual;
      for(const auto &f:batch.flows())if(f.kind==pfg::kFlowRelease)actual.push_back(&f);
      std::size_t at=0;
      for(const auto &release:x->releases)for(const auto &row:release.carriers){
        const auto &f=*actual.at(at++);
        check(f.native_id==row.id && a::same(f.weight,row.weight) && f.species==sb.species().at(release.state) &&
          f.grid_cell==sb.grid().locate(row.position_m.data()),"actual release count/species/cell/ID preserved");
        check(std::memcmp(f.u,row.proper_velocity_m_s.data(),sizeof f.u)==0,"actual release proper velocity preserved");
      }
    }
    sg.take(x->accepted_step,x->stage,batch.recoils(),batch.flows());
  });
  check(release_packets>0 && release_rows>0,"actual stage2 native releases exercised");
  sg.abandon();check(sg.accepted_sha256()==saccepted && S.arena->full_time_boundary()->accepted_arena_sha256==snative &&
    S.arena->full_time_chemistry().state_sha256==srng,"both-stage fixture abandons all accepted changes");
  counts["stage2"]={{"release_packets",release_packets},{"release_carriers",release_rows},
    {"nonzero_internal_release_packets",release_internal},{"admitted_stages",stage2_admitted},{"refused_stages",refused_stage}};
  auto empty=chemfix::make_arena(options,chemistry,{},1000,2);g::Owner eg(config(empty,blob));
  a::Binding eb(eg,empty.reg,blob,options);auto b=binding(empty,1e-12);eg.begin_step(b.accepted_step);
  auto outer=empty.arena->begin_full_time_step(b,1000,64);
  n::ExposureBinding e{b,n::FullTimeSlot::before_transport,b.begin_s,.5*b.dt_s,
    empty.arena->full_time_chemistry().stages,b.accepted_field_sha256};
  auto scope=outer.first_exposure(e);auto t=scope.prepare_chemistry_v3();auto z=a::prepare(scope,t,eg,eb,1);
  check(z.recoils().empty()&&z.flows().empty(),"actual empty stage exports no records");
  J result={{"passed",true},{"checks",checks},{"populated",counts},{"diagnostic_inputs",true},
    {"joint_publication_qualified",false},{"feedback_injection_qualified",false},{"production_qualified",false}};
  std::ofstream(argv[1])<<result.dump(2)<<'\n';std::cout<<result.dump(2)<<'\n';return 0;
}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}}
