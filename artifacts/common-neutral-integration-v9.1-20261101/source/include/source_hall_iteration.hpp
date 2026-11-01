#pragma once
#include "source_stage_v8.hpp"
#include "native_record_v81_view.hpp"
#include <chrono>
#include <ctime>
#include <memory>
#include <optional>

namespace powerflow::source_hall_iteration_v1 {
namespace s=::powerflow::source_stage_v8;
namespace h=::hall_corrective_v3;
namespace n=::raikou::mfem_backend;
inline void require(bool ok,const char *why){if(!ok)throw std::invalid_argument(std::string("source/Hall iteration: ")+why);}
struct Time {double wall_s=0,cpu_s=0;};
struct Clock {
  std::chrono::steady_clock::time_point wall=std::chrono::steady_clock::now();
  timespec cpu{};
  Clock(){if(clock_gettime(CLOCK_THREAD_CPUTIME_ID,&cpu))throw std::runtime_error("thread clock unavailable");}
  Time elapsed()const{timespec now{};if(clock_gettime(CLOCK_THREAD_CPUTIME_ID,&now))throw std::runtime_error("thread clock unavailable");return {std::chrono::duration<double>(std::chrono::steady_clock::now()-wall).count(),double(now.tv_sec-cpu.tv_sec)+1e-9*double(now.tv_nsec-cpu.tv_nsec)};}
};
struct IterationReceipt {
  unsigned fork=0;
  std::size_t carriers=0,electron_births=0,neutral_births=0,bridges=0;
  std::string source_field_sha256,hall_field_sha256;
  bool identical_midpoint=false;
  Time source,field;
  h::ProposalTimings hall;
};
// A genuine prior Hall proposal may seed a new native outer attempt after the
// four-fork cap. It retains no authority to publish a native or gas candidate.
// SourceStage revalidates its exact Hall owner, census and common clock.
struct RetrySeed {
  h::JointProposal proposal;
  std::vector<n::ArenaRateEntry> census;
  RetrySeed(h::JointProposal p,std::vector<n::ArenaRateEntry> c):proposal(std::move(p)),census(std::move(c)){}
};
struct Ready {
  s::PreparedSourceStage source;
  n::PrivateTransportTrial transport;
  std::vector<n::ArenaRateEntry> census;
  h::JointProposal proposal;
  s::PublicationState future;
  Ready(s::PreparedSourceStage src,n::PrivateTransportTrial trial,std::vector<n::ArenaRateEntry> c,
        h::JointProposal p,s::PublicationState f):source(std::move(src)),transport(std::move(trial)),
        census(std::move(c)),proposal(std::move(p)),future(std::move(f)){}
  // Native preparation uses the exact admitted Hall motion, including B3.
  // Acceptance must still go through the native/gas coupled motion capability.
  n::PreparedPrivateMotion prepare_motion(const n::FullTimeStepBinding &common){
    n::FullTimeMotionBinding binding{common,proposal.motion().record_v7_sha256(),common.accepted_field_sha256,
      proposal.candidate_field_sha256(),proposal.candidate_heat_sha256(),source.source_after_sha256()};
    return transport.prepare_total_field_motion(binding,h::matching_motion_admission(proposal));
  }
};
struct Result {
  std::unique_ptr<Ready> ready;
  std::unique_ptr<RetrySeed> retry;
  std::vector<IterationReceipt> iterations;
};
inline Result prepare(n::PrivateFullTimeStep &outer,const s::OwnerSet &owners,const s::AcceptedSourceState &accepted,
                      const n::FullTimeStepBinding &common,const n::GapPicOptions &options,
                      const n::FullTimeChemistryRegistry &registry,const s::s7::BirthPolicy &birth_policy,
                      std::size_t carrier_cap,unsigned maximum_forks=4,const RetrySeed *seed=nullptr,
                      std::span<const s::s7::IntervalEvidence> intervals={}) {
  require(!owners.ports.empty()&&carrier_cap>0,"source owners and finite carrier cap required");
  require(maximum_forks>=1&&maximum_forks<=4,"native outer admits at most four transport forks");
  auto hall=owners.ports.front().hall();require(bool(hall),"missing actual Hall owner");
  auto scope=outer.birth_scope();Result result;
  Clock source_clock;
  auto source=std::make_unique<s::PreparedSourceStage>(seed?
    s::PreparedSourceStage::prepare_against_proposal(owners,accepted,common,scope,seed->proposal,seed->census,birth_policy,intervals):
    s::PreparedSourceStage::prepare(owners,accepted,common,scope,birth_policy,intervals));
  Time source_time=source_clock.elapsed();
  for(unsigned k=0;k<maximum_forks;++k){
    auto current=std::move(source); // outlives this iteration's native trial and Hall proposal
    require(current->preflight_sources(owners),"accepted source/arrival owners changed");
    auto transport=outer.fork_transport_trial(current->native_births());
    auto census=transport.native_census(carrier_cap);const auto &native_census=transport.census_v7();
    std::vector<h::Newborn> births;for(const auto &b:transport.birth_flights())births.push_back({b.native_id,b.birth_time_s});
    hall_driver::MotionBinding motion;
    motion.phase=hall_driver::ParticlePhase::full_time_position_and_proper_velocity;
    motion.accepted_step=common.accepted_step;motion.staged_arena_revision=transport.staged_arena_revision();
    motion.staged_arena_sha256=transport.state_prefix();motion.trial_stage=1;
    motion.begin_s=common.begin_s;motion.end_s=common.end_s;motion.dt_s=common.dt_s;
    motion.mesh_sha256=common.mesh_sha256;motion.material_sha256=common.material_sha256;
    motion.species_sha256=registry.registry_sha256;motion.common_clock_sha256=common.common_clock_sha256;
    motion.kinetic_policy_sha256=hall->motion_context().policy_sha256;
    h::CpuMotionProvider provider;
    h::StepInputs inputs{births,current->loads(),true,native_census.states,&current->native_births()};
    Clock field_clock;
    auto proposal=h::prepare_closed_proposal(*hall,options,census,motion,current->source_after_sha256(),provider,inputs);
    const auto field_time=field_clock.elapsed();
    const auto &field=current->source_field();
    const bool identical=field->matches({field->frame(),hall->mesh(),proposal.trial_E(),proposal.trial_B()});
    result.iterations.push_back({k+1,census.size(),current->reference().electron_births(),current->reference().neutral_births(),
      proposal.motion().bridges().size(),field->frame().field_hash,proposal.midpoint().total_field_sha256(),identical,source_time,field_time});
    result.iterations.back().hall=proposal.timings();
    if(identical){
      require(hall->preflight(proposal,census)&&current->preflight_sources(owners),"actual source/Hall proposal preflight");
      // This also enforces the original matching-field bridge-work tolerance;
      // matching arrays alone cannot replace the actual work comparison.
      auto future=current->prepare_publication(owners,proposal);
      result.ready=std::make_unique<Ready>(std::move(*current),std::move(transport),std::move(census),std::move(proposal),std::move(future));
      return result;
    }
    if(k+1==maximum_forks){
      result.retry=std::make_unique<RetrySeed>(std::move(proposal),std::move(census));
      return result;
    }
    Clock retry_clock;
    auto next=std::make_unique<s::PreparedSourceStage>(s::PreparedSourceStage::prepare_against_proposal(
      owners,accepted,common,scope,proposal,census,birth_policy,intervals));
    source_time=retry_clock.elapsed();
    // The old transport/proposal are abandoned at the end of this iteration;
    // no source inventory, arrival marks, Hall state or native state publishes.
    source=std::move(next);
  }
  throw std::logic_error("unreachable bounded source iteration");
}
} // namespace powerflow::source_hall_iteration_v1
