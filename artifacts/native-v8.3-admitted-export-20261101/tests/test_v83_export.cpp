// Diagnostic component controls, never a physical breakdown qualification.
#include "chem_v8_fixture.hpp"
#include "chem_v3_native.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

namespace m = raikou::mfem_backend;
namespace c = m::chem8;
using J = nlohmann::json;
// Test-only construction of self-consistent diagnostic facts. These are not
// PreparedArena capabilities and cannot authorize an adapter transaction.
namespace raikou::mfem_backend::full_time_detail {
std::string export_admitted_root(const FullTimeChemistryStageExport &);
std::string export_debits_root(const FullTimeChemistryStageExport &);
std::string export_releases_root(const FullTimeChemistryStageExport &);
std::string export_digest(const FullTimeChemistryStageExport &);
}
int checks = 0;
void check(bool ok, const char *why) { ++checks; if (!ok) throw std::runtime_error(why); }
template<class F> void refuses(F f, const char *why) {
  bool caught = false; try { f(); } catch (const std::exception &) { caught = true; }
  check(caught, why);
}
void seal(m::FullTimeChemistryStageExport &x) {
  namespace d = m::full_time_detail;
  x.admitted_root=d::export_admitted_root(x); x.target_debits_root=d::export_debits_root(x);
  x.releases_root=d::export_releases_root(x); x.sha256=d::export_digest(x);
}
pfc::EventRec event(const m::FullTimeChemistryAdmittedEvent &r) {
  pfc::EventRec e{}; check(r.event_bytes.size()==sizeof e,"record ABI");
  std::memcpy(&e,r.event_bytes.data(),sizeof e); return e;
}
bool gas_return(const pfc::EventRec &e) {
  for (pfc_u32 k=0;k<e.n_out;++k) if(e.out_dest[k]==PFC_DEST_GAS) return true;
  return false;
}
int main(int argc,char **argv) {
  try {
    if(argc!=2) throw std::runtime_error("usage: test_v83_export OUT");
    const std::filesystem::path out=argv[1];std::filesystem::create_directories(out);
    auto o=chemfix::options();auto cfg=chemfix::chemistry_configuration(o);
    chemfix::Population pop;pop.carriers=3000;pop.seed=11;pop.clustered=true;
    pop.water_weight=2e13;pop.electron_eV_min=5;pop.electron_eV_max=300;
    auto packets=chemfix::population(pop);
    auto A=chemfix::make_arena(o,cfg,packets,120000,2);
    std::vector<std::shared_ptr<const m::FullTimeChemistryStageExport>> exports;
    std::optional<pfc::EventRec> consume,gas;
    std::size_t actual_whole=0,actual_partial=0,actual_overlap=0;
    J trace=J::array();
    for(int k=0;k<2;++k) {
      auto r=chemfix::step(A,1e-12,k);
      for(const auto &receipt:{r.stage1,r.stage2}) {
        check(bool(receipt->admitted_export),"attached immutable export");
        const auto &x=*receipt->admitted_export;
        m::validate_full_time_chemistry_stage_export(x,false);
        check(x.raw_records_root==receipt->records_root && x.events.size()==receipt->admitted_events,
              "raw/admitted receipt binding");
        check(x.target_debits.size()==receipt->debited_packets && x.releases.size()==receipt->released_packets,
              "actual flow counts");
        check(x.stage==(receipt==r.stage1?1u:2u),"both stage identities");
        for(const auto &v:x.events) {
          auto e=event(v);
          if(v.consumes_tracked_target && !consume) consume=e;
          if(gas_return(e) && !gas) gas=e;
        }
        for(const auto &d:x.target_debits) {
          actual_whole+=d.remainder.packet_id.empty();actual_partial+=!d.remainder.packet_id.empty();
          actual_overlap+=d.also_projectile;
        }
        if(std::any_of(x.target_debits.begin(),x.target_debits.end(),[](const auto &d){return d.also_projectile||d.internal_energy_deferred;}))
          refuses([&]{m::validate_full_time_chemistry_stage_export(x);},"actual ambiguous target refuses consuming export");
        exports.push_back(receipt->admitted_export);
        trace.push_back({{"export",x.sha256},{"raw",x.raw_records_root},{"admitted",x.admitted_root},
                         {"debits",x.target_debits_root},{"releases",x.releases_root}});
      }
    }
    check(bool(consume)&&bool(gas),"populated consuming and GAS-return events");
    A.arena->validate_full_time_history();
    // Live capability controls: deterministic retries reproduce facts but get a
    // fresh immutable object. Data lifetime never extends publication authority.
    auto C=chemfix::make_arena(o,cfg,packets,120000,2);
    auto boundary=*C.arena->full_time_boundary();
    const auto accepted_rng=C.arena->full_time_chemistry().state_sha256;
    std::shared_ptr<const m::FullTimeChemistryStageExport> retained;
    {
      m::FullTimeStepBinding b{boundary.accepted_step,boundary.accepted_arena_revision,boundary.time_s,
        boundary.time_s+1e-12,1e-12,boundary.accepted_arena_sha256,boundary.field_sha256,
        boundary.mesh_sha256,boundary.material_sha256,chemfix::sha("live-export-clock")};
      auto outer=C.arena->begin_full_time_step(b,120000,64);
      m::ExposureBinding exposure{b,m::FullTimeSlot::before_transport,b.begin_s,.5*b.dt_s,
        C.arena->full_time_chemistry().stages,b.accepted_field_sha256};
      auto scope=outer.first_exposure(exposure);
      {
        auto token=scope.prepare_chemistry_v3();
        check(scope.can_publish(token),"live candidate publication preflight");
        retained=token.candidate_admitted_export(false);
        check(retained==token.candidate_chemistry_stage()->admitted_export,"getter returns same immutable object");
        check(scope.can_publish(token),"export extraction preserves live capability");
        auto moved=std::move(token);
        check(!scope.can_publish(token) && scope.can_publish(moved),"moved token authority transfers");
        refuses([&]{(void)token.candidate_admitted_export(false);},"moved token cannot extract facts");
      }
      auto retry=scope.prepare_chemistry_v3();
      auto again=retry.candidate_admitted_export(false);
      check(again.get()!=retained.get() && again->sha256==retained->sha256,
            "deterministic retry has fresh immutable object with identical facts");
      auto foreign=chemfix::make_arena(o,cfg,packets,120000,2);
      auto other_outer=foreign.arena->begin_full_time_step(b,120000,64);
      auto other_scope=other_outer.first_exposure(exposure);
      check(!other_scope.can_publish(retry),"identical foreign scope cannot authorize export token");
      check(scope.try_publish(retry) && !scope.can_publish(retry),"exposure publishes exactly once");
      check(retry.candidate_admitted_export(false)==again,"published token retains read-only facts");
      // Outer step intentionally abandoned; no accepted population/RNG changes.
    }
    check(C.arena->full_time_boundary()->accepted_arena_sha256==boundary.accepted_arena_sha256 &&
          C.arena->full_time_chemistry().state_sha256==accepted_rng,"abandoned outer preserves accepted state and RNG");
    m::validate_full_time_chemistry_stage_export(*retained,false);
    // Empty stages still have a bound export, with zero actual flows.
    auto empty=chemfix::make_arena(o,cfg,{},1000,2);
    auto er=chemfix::step(empty,1e-12,0);
    for(const auto &r:{er.stage1,er.stage2}) {
      m::validate_full_time_chemistry_stage_export(*r->admitted_export);
      check(r->admitted_export->events.empty() && r->admitted_export->target_debits.empty() &&
            r->admitted_export->releases.empty(),"empty stage has zero admitted flows");
    }
    // Source/count/clock tampering retains the old commitment and must refuse.
    for(int k=0;k<8;++k) {
      auto x=*exports.front();
      if(k==0)x.raw_records_root=chemfix::sha("wrong raw root");
      if(k==1)++x.accepted_step;
      if(k==2)x.common_clock_sha256=chemfix::sha("wrong clock");
      if(k==3)x.frame_snapshot_sha256=chemfix::sha("wrong frame");
      if(k==4)x.chemistry_before_sha256=chemfix::sha("wrong chemistry");
      if(k==5)x.events.at(0).weight*=2;
      if(k==6)x.admissions.at(0).raw_begin++;
      if(k==7)x.target_debits.at(0).amount=x.target_debits.at(0).amount.credited(1);
      refuses([&]{m::validate_full_time_chemistry_stage_export(x,false);},"tampered export refused");
    }
    // Native resolution prefix: first consuming event fits, second does not;
    // the entire remaining suffix is refused after that inventory conflict.
    auto setup=c::make_setup(cfg,o);auto before=c::initial_state(*setup);
    auto active_for=[&](pfc::EventRec e,unsigned n) {
      c::StageActive a;a.cell_atlas={0};
      c::ActiveCarrier ac{};ac.id=1;ac.row=0;ac.cell=0;ac.begin=0;ac.end=n;ac.weight=1;
      ac.species=setup->reg->blob->view.groups[setup->reg->blob->view.procs[e.process].group].projectile;
      ac.x[0]=ac.x[1]=.004;
      a.carriers.push_back(ac);a.totals.carriers=1;a.totals.candidates=n;
      a.totals.status_counts[PFC_ST_OK]=1;
      a.totals.exposure.resize(setup->reg->blob->n_regions);
      a.totals.bound.resize(setup->reg->blob->n_regions);
      a.totals.support.resize(setup->reg->blob->n_regions);
      for(unsigned j=0;j<n;++j){e.carrier=1;e.cell=0;e.candidate=j;e.weight=1;a.events.push_back(e);}
      return a;
    };
    auto active=active_for(*consume,3);
    c::DebitAccess da;
    da.inventory=[](auto,auto){c::ESum n;n.add(1.5);return n;};
    da.packets=[](auto,auto){c::ESum n;n.add(1.5);return std::vector<std::pair<std::uint32_t,c::ESum>>{{7,n}};};
    const auto raw=c::sha_hex(std::string(reinterpret_cast<const char *>(active.events.data()),active.events.size()*sizeof(pfc::EventRec)));
    auto res=c::resolve(*setup,active,{0,1,1e-12,10000},*before,da);
    check(res.admissions.size()==1 && res.admissions[0].admitted_count==1 &&
          res.admissions[0].raw_count==3 && res.admitted_records.size()==1 && res.debits.refused==1,
          "prefix rejects entire raw suffix");
    check(res.plan.debits.size()==1 && !res.plan.debits[0].whole,"actual partial target debit");
    da.inventory=[](auto,auto){c::ESum n;n.add(.5);return n;};
    auto rejected=c::resolve(*setup,active,{0,1,1e-12,10000},*before,da);
    check(rejected.admitted_records.empty() && rejected.plan.debits.empty(),"shorter admitted prefix");
    check(raw==c::sha_hex(std::string(reinterpret_cast<const char *>(active.events.data()),active.events.size()*sizeof(pfc::EventRec))),
          "admission changes do not rewrite raw records");
    active=active_for(*gas,1);
    da.inventory=[](auto,auto)->c::ESum{throw std::runtime_error("GAS event requested inventory");};
    da.packets=[](auto,auto)->std::vector<std::pair<std::uint32_t,c::ESum>>{throw std::runtime_error("GAS event requested packets");};
    auto returned=c::resolve(*setup,active,{0,1,1e-12,10000},*before,da);
    check(returned.admitted_records.size()==1 && returned.plan.debits.empty() &&
          !returned.admitted_records.front().consumes_tracked_target,"GAS return has no count debit");
    // Public contract cases with exact count and moment operands, including a
    // remainder needing multiple positive carriers. These are diagnostic facts,
    // explicitly not live candidate publication capabilities.
    auto x=*exports.front();x.admissions.clear();x.events.clear();x.target_debits.clear();x.releases.clear();
    x.stage=1;x.purpose=PFC_PURPOSE_MCC;x.frozen_physical_stage_s=x.begin_s;x.next_native_id=x.base_native_id+2;
    auto e=*consume;e.carrier=1;e.cell=0;e.candidate=0;e.weight=1e-20;
    x.admissions.push_back({1,0,1,1,0,0});
    x.events.push_back({0,0,e.target_species,true,e.weight,std::string(reinterpret_cast<const char *>(&e),sizeof e)});
    m::FullTimeChemistryTargetDebit d;d.species=e.target_species;d.amount=m::ArenaCount(e.weight);
    d.before.packet_id="manufactured/target";d.before.state=setup->reg->states.at(e.target_species);d.before.count=m::ArenaCount(1);
    m::GapParticle g;g.id=2;g.position_m={.004,.004,0};g.proper_velocity_m_s={1,2,3};g.weight=1;
    g.mass_kg=setup->reg->mass_kg.at(e.target_species);g.charge_C=setup->reg->charge.at(e.target_species)*m::pic_e;
    d.before.carriers={g};d.remainder=d.before;d.remainder.count=d.before.count.debit(e.weight);d.remainder.carriers.clear();
    const auto weights=d.remainder.count.positive_carriers();
    for(std::size_t k=0;k<weights.size();++k){auto r=g;r.id=k?x.base_native_id+k-1:2;r.weight=weights[k];d.remainder.carriers.push_back(r);}
    x.target_debits={d};seal(x);m::validate_full_time_chemistry_stage_export(x);
    check(weights.size()>1,"nonbinary remainder preserves exact count and individual weights");
    auto overlap=x;overlap.target_debits[0].also_projectile=true;seal(overlap);
    m::validate_full_time_chemistry_stage_export(overlap,false);
    refuses([&]{m::validate_full_time_chemistry_stage_export(overlap);},"ambiguous projectile/target debit refused");
    auto internal=x;internal.target_debits[0].before.carried_internal_energy_J=1;
    internal.target_debits[0].remainder.carried_internal_energy_J=1;
    internal.target_debits[0].internal_energy_deferred=true;seal(internal);
    m::validate_full_time_chemistry_stage_export(internal,false);
    refuses([&]{m::validate_full_time_chemistry_stage_export(internal);},"deferred internal ownership refused");
    // Real whole-packet target consumption: equal charged/neutral packet weights.
    // Only stage one is prepared; no enormous neutral release workload is run.
    auto whole_pop=pop;whole_pop.charged_weight=whole_pop.water_weight;
    auto W=chemfix::make_arena(o,cfg,chemfix::population(whole_pop),120000,2);
    {
      auto wb=*W.arena->full_time_boundary();
      m::FullTimeStepBinding b{wb.accepted_step,wb.accepted_arena_revision,wb.time_s,wb.time_s+1e-12,
        1e-12,wb.accepted_arena_sha256,wb.field_sha256,wb.mesh_sha256,wb.material_sha256,chemfix::sha("whole-target-clock")};
      auto outer=W.arena->begin_full_time_step(b,120000,64);
      m::ExposureBinding exposure{b,m::FullTimeSlot::before_transport,b.begin_s,.5*b.dt_s,
        W.arena->full_time_chemistry().stages,b.accepted_field_sha256};
      auto scope=outer.first_exposure(exposure);auto token=scope.prepare_chemistry_v3();
      auto facts=token.candidate_admitted_export(false);
      check(std::any_of(facts->target_debits.begin(),facts->target_debits.end(),[](const auto &d){return d.remainder.packet_id.empty();}),
            "actual whole-packet target debit exported with absent remainder");
      check(scope.can_publish(token),"whole-target diagnostic token remains live");
    }
    // Genuine native path with nonbinary target remainders and explicit host fallback.
    pop.charged_weight=1e-10;auto B=chemfix::make_arena(o,cfg,chemfix::population(pop),120000);
    auto fractional=chemfix::step(B,1e-12,0);bool multi=false;
    for(const auto &r:{fractional.stage1,fractional.stage2})for(const auto &f:r->admitted_export->target_debits)
      multi|=f.remainder.carriers.size()>1;
    check(multi,"native exact target remainder exported as actual positive carriers");
    std::ofstream(out/"export-controls.json")<<J{{"checks",checks},{"passed",true},{"trace",trace},
      {"actual_whole_debits",actual_whole},{"actual_partial_debits",actual_partial},{"actual_projectile_target_overlaps",actual_overlap},
      {"note","Diagnostic component facts and contract controls; no gas feedback or energy closure."}}.dump(2)<<"\n";
    std::cout<<checks<<" admitted-export controls passed\n";
    return 0;
  }catch(const std::exception &e){std::cerr<<"error: "<<e.what()<<"\n";return 1;}
}
