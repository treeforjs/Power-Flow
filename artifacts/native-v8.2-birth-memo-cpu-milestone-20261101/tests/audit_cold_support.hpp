#pragma once
#include "chem_v8_fixture.hpp"
#include "raikou/json.hpp"
#include <filesystem>
#include <functional>
#include <stdexcept>
namespace audit {
namespace m=raikou::mfem_backend;
inline int cold_support(const std::filesystem::path &out,
                        m::FullTimeChemistryExecution execution=m::FullTimeChemistryExecution::host,
                        const std::function<void(m::PicArena&)> &attach={}) {
  auto need=[](bool ok,const char *s){if(!ok)throw std::runtime_error(s);};
  const auto o=chemfix::options(45);
  auto c=chemfix::chemistry_configuration(o);
  c.execution=execution;
  const m::PicVector x{.00215,.00215,0.0}, u{0,0,0};
  const std::vector<m::FullTimeAuthoredPacket> p{
      {"cold/e","pfc-v3/e","manufactured/cold-support",m::ArenaCount(1,0),x,u},
      {"cold/OH+","pfc-v3/OH+","manufactured/cold-support",m::ArenaCount(1,0),x,u}};
  for(const double dt:{1e-12,1e-15}) {
    auto strict=c; strict.policy.admitted_regions.clear();
    auto A=chemfix::make_arena(o,strict,p,1000);
    if(attach) attach(*A.arena);
    const auto before=*A.arena->full_time_boundary();
    const auto status=A.arena->full_time_chemistry();
    bool refused=false;
    try{(void)chemfix::step(A,dt,0);}
    catch(const m::FullTimeChemistryRefusal&e){refused=e.kind()==m::FullTimeChemistryRefusal::Kind::coverage_policy;}
    need(refused,"cold unknown support must refuse regardless of timestep");
    need(A.arena->full_time_boundary()->accepted_arena_sha256==before.accepted_arena_sha256,
         "cold refusal mutated accepted state/RNG/stock");
    need(A.arena->full_time_chemistry().coverage_json==status.coverage_json,
         "cold refusal mutated coverage ledger");
  }
  auto A=chemfix::make_arena(o,c,p,1000,1);
  if(attach)attach(*A.arena);
  const auto r=chemfix::step(A,1e-12,0);
  need(r.stage1->admitted_events+r.stage2->admitted_events==0,"cold diagnostic unexpectedly had events");
  const auto status=A.arena->full_time_chemistry();
  const auto coverage=nlohmann::json::parse(status.coverage_json);
  std::uint64_t support=0;
  for(const auto&row:coverage.at("cumulative")) {
    support+=row.at("reachable_support").get<std::uint64_t>();
    need(row.at("exposure").get<double>()==0 && row.at("bound").get<double>()==0,
         "cold support fabricated a physical rate");
  }
  need(support>0,"cold diagnostic lost reachable support");
  std::filesystem::create_directories(out);
  A.arena->save((out/"cold-support.json").string());
  auto restored=m::PicArena::restore_full_time((out/"cold-support.json").string(),o,A.states,A.channels,
                                             *A.arena->full_time_boundary());
  need(restored->full_time_chemistry().coverage_json==status.coverage_json,"cold support restore mismatch");
  return 7;
}
}
