#include "chem_v8_fixture.hpp"
#include "audit_cold_support.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
namespace m = raikou::mfem_backend;
void need(bool ok, const char* msg) { if (!ok) throw std::runtime_error(msg); }
int main(int argc, char**argv) {
  try {
    const std::filesystem::path out = argc > 1 ? argv[1] : ".";
    std::filesystem::create_directories(out);
    auto o = chemfix::options(45);
    auto diagnostic = chemfix::chemistry_configuration(o);
    need(std::any_of(diagnostic.policy.admitted_regions.begin(), diagnostic.policy.admitted_regions.end(),
                    [](const std::string &r){return r.size()>=1024;}), "long valid region fixture");
    chemfix::Population p;
    p.carriers=1500; p.seed=71; p.water_weight=2e13; p.clustered=true;
    p.electron_eV_min=5; p.electron_eV_max=300;
    auto A=chemfix::make_arena(o,diagnostic,chemfix::population(p),40000,2);
    const auto r=chemfix::step(A,1e-12,0);
    need(r.stage1->admitted_events+r.stage2->admitted_events>0,"fixture must exercise events");
    need(r.stage1->diagnostic_field_target_approximation && !r.stage1->full_system_energy_closure,
         "receipt must explicitly exclude full-system closure");
    const auto status=A.arena->full_time_chemistry();
    need(status.diagnostic_field_target_approximation && !status.full_system_energy_closure,"diagnostic status");
    auto strict=diagnostic; strict.diagnostic_field_target_approximation=false;
    auto B=chemfix::make_arena(o,strict,chemfix::population(p),40000,2);
    const auto before=*B.arena->full_time_boundary();
    bool refused=false;
    try { (void)chemfix::step(B,1e-12,0); }
    catch(const m::FullTimeChemistryRefusal&e) { refused=e.kind()==m::FullTimeChemistryRefusal::Kind::target_energy_owner; }
    need(refused,"strict chemistry must refuse unsupported target energy closure");
    need(B.arena->full_time_boundary()->accepted_arena_sha256==before.accepted_arena_sha256,"refusal mutated state");
    const auto boundary=*A.arena->full_time_boundary();
    A.arena->save((out/"diagnostic.json").string());
    auto restored=m::PicArena::restore_full_time((out/"diagnostic.json").string(),o,A.states,A.channels,boundary);
    need(restored->full_time_chemistry().diagnostic_field_target_approximation &&
         !restored->full_time_chemistry().full_system_energy_closure,"restore lost diagnostic scope");
    nlohmann::json saved;
    {std::ifstream f(out/"diagnostic.json"); f>>saved;}
    auto payload=saved.at("payload");
    payload["configuration"]["v8"]["chemistry"]["policy"]["admitted_regions"][0]="unknown/"+std::string(2048,'x');
    {std::ofstream f(out/"unknown-region.json");
     f<<nlohmann::json{{"payload",payload},{"sha256",chemfix::sha(payload.dump())}}.dump();}
    std::string message;
    try{(void)m::PicArena::restore_full_time((out/"unknown-region.json").string(),o,A.states,A.channels,boundary);}
    catch(const std::exception&e){message=e.what();}
    need(message.find("unknown admitted region")!=std::string::npos,"unknown long region must refuse exact registry validation");
    std::cout<<"audit target-energy guard: diagnostic scope, strict refusal, no mutation, long valid region restore and unknown region refusal passed\n";
    std::cout << audit::cold_support(out/"cold-host") << " cold reachable-support controls passed\n";
    return 0;
  } catch(const std::exception&e) {std::cerr<<e.what()<<'\n';return 1;}
}
