#include "common_step.hpp"
#include "request_fixture.hpp"
#include "positive_source_fixture.hpp"
#include <cstdlib>
#include <new>
namespace allocation_probe { bool enabled=false;std::size_t count=0; }
void* operator new(std::size_t n){if(allocation_probe::enabled)++allocation_probe::count;if(auto p=std::malloc(n?n:1))return p;throw std::bad_alloc();}
void* operator new[](std::size_t n){return ::operator new(n);}
void operator delete(void*p)noexcept{std::free(p);}void operator delete[](void*p)noexcept{std::free(p);}
void operator delete(void*p,std::size_t)noexcept{std::free(p);}void operator delete[](void*p,std::size_t)noexcept{std::free(p);}
void* operator new(std::size_t n,std::align_val_t a){if(allocation_probe::enabled)++allocation_probe::count;void*p=nullptr;if(posix_memalign(&p,std::size_t(a),n?n:1))throw std::bad_alloc();return p;}
void* operator new[](std::size_t n,std::align_val_t a){return ::operator new(n,a);}
void operator delete(void*p,std::align_val_t)noexcept{std::free(p);}void operator delete[](void*p,std::align_val_t)noexcept{std::free(p);}
void operator delete(void*p,std::size_t,std::align_val_t)noexcept{std::free(p);}void operator delete[](void*p,std::size_t,std::align_val_t)noexcept{std::free(p);}
namespace powerflow::common_step_v9 {
struct StepTestAccess {static void publish(Step&s)noexcept{s.publish_prevalidated();}};
}
namespace cs=powerflow::common_step_v9;
J roots(const std::shared_ptr<cs::Run>&run){
 const auto&o=run->owners();J result={{"native",o.native->state_prefix()},{"gas",o.gas->accepted_sha256()},
 {"field",o.hall->field_sha256()},{"heat",o.hall->heat_sha256()},{"source",o.source_sha256},
 {"source_metadata",o.source_states},{"arrival",o.sources.arrival->state_sha256()},
 {"arrival_metadata",o.arrival_state_sha256},{"common",run->receipt_root()},
 {"configuration",run->configuration_sha256()},{"conservation",run->ledger()?run->ledger()->sha256:std::string{}},{"sources",J::array()}};
 for(const auto&p:o.sources.ports)result["sources"].push_back(p.source()->state_sha256());return result;
}
cg::Request positive_request(const char*trace){
 auto r=request(trace);auto arrivals=r.arrival_policy->configuration();
 powerflow::positive_source_fixture_v1::Spec spec;spec.water_enabled=true;spec.work_function_eV=2.;spec.electron_weight=1e-23;
 powerflow::positive_source_fixture_v1::apply(r.ports,arrivals,r.options,r.native_configuration.chemistry,spec);
 r.arrival_policy=std::make_shared<const ss::pa::Policy>(arrivals);
 std::vector<ss::hb::SliceSource>slices;for(const auto&p:r.ports)slices.push_back(p.slice);
 h::MotionForms forms;forms.native_v7=forms.native_v81_bridges=true;forms.neutral_stream_dt_max_s=1e-11;
 r.model=std::make_shared<h::Model>(h::MeshSpec{},trace,h::make_cpu_pcg(),h::CarrierLimits{},ss::hb::declaration(slices),1e-14,forms);
 r.seed=r.model->dc_seed();return r;
}
J receipt(const cs::Receipt&r){
 J out={{"step",r.step},{"begin_s",r.begin_s},{"end_s",r.end_s},{"attempts",r.attempts},
 {"final_carriers",r.final_carriers},{"stage1_events",r.stage1_events},{"stage2_events",r.stage2_events},
 {"parent",r.parent},{"sha256",r.sha256},{"final_census",r.final_census},
 {"conservation_sha256",r.conservation_sha256},{"raw_energy_residual_J",r.raw_energy_residual_J},
 {"accounted_energy_residual_J",r.accounted_energy_residual_J},{"energy_bound_J",r.energy_bound_J},
 {"phases",J::array()},{"iterations",J::array()}};
 for(const auto&p:r.phases)out["phases"].push_back({{"name",p.phase},{"wall_s",p.time.wall_s},{"cpu_s",p.time.cpu_s}});
 for(const auto&i:r.iterations){const auto&t=i.hall;out["iterations"].push_back({{"fork",i.fork},{"carriers",i.carriers},{"electrons",i.electron_births},
 {"neutrals",i.neutral_births},{"bridges",i.bridges},{"identical_midpoint",i.identical_midpoint},
 {"source_wall_s",i.source.wall_s},{"source_cpu_s",i.source.cpu_s},{"field_wall_s",i.field.wall_s},{"field_cpu_s",i.field.cpu_s},
 {"hall",{{"admit_census_s",t.admit_census_s},{"motion_s",t.motion_s},{"support_s",t.support_s},{"field_prepare_s",t.field_prepare_s},
 {"solve_subset_s",t.solve_s},{"residual_check_subset_s",t.residual_check_s},{"seal_s",t.seal_s},{"heat_s",t.heat_s},
 {"motion_record_seal_s",t.motion_v2_seal_s},{"iterations",t.iterations},{"total_s",t.total_s}}}});}
 return out;
}
int main(int argc,char**argv){
 if(argc!=3)return 2;J results={{"scope","manufactured neutral desorption common integration with independent whole-system conservation and cumulative ledger plus authentic magnetic electron-source refusal; no conductor sinks; positive emission, full plasma or common restart unqualified"},
 {"cases",J::array()},{"steps",J::array()},{"passed",true}};int failures=0;
 const auto out=std::filesystem::path(argv[2]);std::filesystem::create_directories(out);
 auto check=[&](const char*name,bool ok){results["cases"].push_back({{"name",name},{"passed",ok}});failures+=!ok;std::cout<<(ok?"PASS ":"FAIL ")<<name<<std::endl;};
 try{
  ss::s7::BirthPolicy policy;policy.subcell_extent_m=20e-6;policy.return_extent_m=0;
  for(int bad=0;bad<2;++bad){auto config=request(argv[1]);auto&spec=config.ports.front().slice.settings.front().spectrum;
   if(bad==0)spec.electron_mass_kg=std::nextafter(spec.electron_mass_kg,INFINITY);else spec.electron_mass_identity="forged";
   bool refused=false;try{auto unused=cs::Run::create(config);}catch(const std::exception&){refused=true;}
   check("authored base mass or provenance error refuses before rebind",refused);
  }
  {
   auto run=cs::Run::create(positive_request(argv[1]));const auto before=roots(run);bool refused=false;std::string reason;
   try{auto unused=cs::Step::prepare(run,1e-14,policy);}catch(const std::exception&e){refused=true;reason=e.what();}
   check("positive source refuses unsupported straight magnetic projection",refused&&reason.find("Straight matching projection fails local B gyro/deflection criterion")!=std::string::npos);
   check("positive source refusal preserves every accepted owner",roots(run)==before&&!run->owners().gas->step_open());
   results["positive_source_refusal"]={{"reason",reason},{"accepted_unchanged",roots(run)==before},{"scope","not a passing electron-source evolution; time-resolved image/Lorentz layer required"}};
  }
  for(bool positive:{false}){
   auto config=positive?positive_request(argv[1]):request(argv[1]);
   auto run=cs::Run::create(config);auto reference=cs::Run::create(config);
   check("root conservation contract is distinct and reproducible",run->configuration_sha256()!=run->owners().configuration_sha256&&run->configuration_sha256()==reference->configuration_sha256());
   for(unsigned step=0;step<2;++step){
    const auto before=roots(run);const auto before_ledger=run->ledger();std::string abandoned;
    {auto candidate=cs::Step::prepare(run,1e-14,policy);check("candidate holds every accepted owner unchanged",roots(run)==before);
     check("prepared candidate preflights",candidate.preflight());abandoned=candidate.receipt().sha256;}
    check("abandon rolls back native gas Hall sources arrival",roots(run)==before&&!run->owners().gas->step_open());
    auto candidate=cs::Step::prepare(run,1e-14,policy);
    check("abandon retry reproduces exact receipt",candidate.receipt().sha256==abandoned);
    check("retry all-owner preflight",candidate.preflight());
    const auto expected=candidate.receipt().sha256;
    allocation_probe::count=0;allocation_probe::enabled=true;
    cs::StepTestAccess::publish(candidate);
    allocation_probe::enabled=false;
    check("joint accepted publication allocates zero times",allocation_probe::count==0);
    check("receipt preserved by publication",candidate.receipt().sha256==expected&&run->receipt_root()==expected);
    const auto ledger=run->ledger();
    check("accepted immutable conservation ledger matches common receipt",ledger&&ledger->sha256==candidate.receipt().conservation_sha256&&cs::ledger_hash(*ledger)==ledger->sha256&&ledger->accepted_step==step+1&&ledger->time_s==candidate.receipt().end_s);
    check("independent per-step energy balance meets derived bound",std::isfinite(candidate.receipt().energy_bound_J)&&std::abs(candidate.receipt().accounted_energy_residual_J)<=candidate.receipt().energy_bound_J);
    auto lower=ledger->accounted_residual_J,upper=lower;lower.merge(ledger->admitted_bound_sum_J);upper.merge(ledger->admitted_bound_sum_J,1);
    check("exact cumulative energy balance meets cumulative derived bound",lower.normalised().back()>=0&&(upper.normalised().back()<0||upper.zero()));
    auto H=cs::cm::diff(ledger->accepted_material.hydrogen,ledger->initial_material.hydrogen);H.merge(ledger->boundary.conductor.hydrogen);H.merge(ledger->boundary.exterior.hydrogen);
    auto O=cs::cm::diff(ledger->accepted_material.oxygen,ledger->initial_material.oxygen);O.merge(ledger->boundary.conductor.oxygen);O.merge(ledger->boundary.exterior.oxygen);
    auto Q=cs::cm::diff(ledger->accepted_material.charge_e,ledger->initial_material.charge_e);Q.merge(ledger->boundary.conductor.charge_e);Q.merge(ledger->boundary.exterior.charge_e);Q.merge(ledger->emitted_electrode_charge_e);
    check("exact whole-history hydrogen oxygen and charge close",H.zero()&&O.zero()&&Q.zero());
    check("conservation ledger preserves immutable accepted anchor",!before_ledger||(ledger->initial_energy_J==before_ledger->initial_energy_J&&cs::same_inventory(ledger->initial_material,before_ledger->initial_material)&&ledger->begin_s==before_ledger->begin_s&&before_ledger->accepted_step==step));
    check("published token cannot publish twice",!candidate.preflight()&&!candidate.try_publish());
    const auto&owners=run->owners();const auto boundary=owners.native->full_time_boundary();
    check("all common clocks advance together",boundary->time_s==owners.hall->time_s()&&boundary->accepted_step==step+1&&owners.hall->revision()==step+1&&!owners.gas->step_open());
    check("actual source hashes equal accepted metadata",owners.sources.arrival->state_sha256()==owners.arrival_state_sha256&&
     [&]{for(std::size_t i=0;i<owners.sources.ports.size();++i)if(owners.sources.ports[i].source()->state_sha256()!=owners.source_states[i])return false;return true;}());
    auto reference_candidate=cs::Step::prepare(reference,1e-14,policy);check("independent continuous run accepts",reference_candidate.try_publish());
    check("abandon retry matches independent continuous roots",roots(run)==roots(reference));
    auto row=receipt(candidate.receipt());row["positive_electron_fixture"]=positive;
    row["cumulative"]={{"raw_energy_residual_J",ledger->raw_residual_J.value()},{"accounted_energy_residual_J",ledger->accounted_residual_J.value()},{"energy_bound_sum_J",ledger->admitted_bound_sum_J.value()},{"events",ledger->events}};
    results["steps"].push_back(std::move(row));
    std::ofstream(out/"common.json")<<results.dump(2)<<'\n';
   }
  }
 }catch(const std::exception&e){allocation_probe::enabled=false;results["error"]=e.what();check("unexpected refusal",false);std::cerr<<e.what()<<std::endl;}
 results["passed"]=failures==0;results["failures"]=failures;std::ofstream(out/"common.json")<<results.dump(2)<<'\n';return failures?1:0;
}
