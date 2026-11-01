// Derived from frozen chem_v8_fixture.hpp; only adds stage preparer injection.
#pragma once
#include "chem_v8_fixture.hpp"
#include "coupled_trial.hpp"
#include "allocation_probe.hpp"
namespace chemfix {
inline long measured_motion_allocations=-1;
inline StepReport coupled_step(Arena &A,double dt,int k,std::shared_ptr<pf_gas_v1::Owner> gas,
 const std::function<void(native_gas_common::CoupledTrial&,m::CollisionScope&,unsigned)>&before_stage,
 const std::function<void(native_gas_common::CoupledTrial&)>&before_seal,
 const std::function<void(native_gas_common::CoupledTrial&)>&before_publish,
 nlohmann::json &receipt,std::shared_ptr<const native_gas_common::CoupledTrial::Continuation>&continuation,bool publish=true,const std::string&motion_control={}) {
  StepReport R;
  Timer T;
  auto &arena = *A.arena;
  const auto &o = A.o;
  const auto bnd = arena.full_time_boundary();
  m::FullTimeStepBinding b{bnd->accepted_step, bnd->accepted_arena_revision, bnd->time_s, bnd->time_s + dt, dt,
                           bnd->accepted_arena_sha256, bnd->field_sha256, bnd->mesh_sha256, bnd->material_sha256,
                           sha("clock/" + std::to_string(k))};
  const auto stages = arena.full_time_chemistry().stages;
  T.start();
  auto raw_outer = arena.begin_full_time_step(b, A.cfg.limits.maximum_carriers, 64);
  native_gas_common::CoupledTrial coupled(std::move(raw_outer),gas,A.cfg.chemistry,A.o,b.accepted_step,A.cfg.limits.maximum_carriers,continuation);
  auto &outer=coupled.outer();
  T.stop(R.phases["begin"]);
  m::ExposureBinding e1{b, m::FullTimeSlot::before_transport, b.begin_s, .5 * dt, stages, b.accepted_field_sha256};
  T.start();
  auto s1 = outer.first_exposure(e1);
  if(before_stage)before_stage(coupled,s1,1);
  R.stage1=coupled.run_stage(s1);
  T.stop(R.phases["stage1"]);

  T.start();
  std::optional<m::PreparedBirths> born;
  if(motion_control=="motion-flows" || motion_control=="motion-empty-flows" || motion_control.starts_with("allocation-motion-")) {
    auto bs=outer.birth_scope();m::FullTimeBirthCredit credit;
    credit.arrival_step=b.accepted_step;credit.base_native_id=bs.base_native_id();
    credit.plan_sha256=sha("manufactured/motion-plan");credit.staged_token_sha256=sha("manufactured/motion-token");
    credit.source_before_sha256=bs.source_before_sha256();credit.source_after_sha256=sha("manufactured/motion-source-after");
    credit.streams.push_back({0,"pfc-v3/H2O",m::ArenaCount(8,0)});
    for(unsigned i=0;i<4;++i) {
      m::FullTimeBirth v;v.step=b.accepted_step;v.index=i;v.native_id=bs.base_native_id()+i;
      v.packet_id=m::full_time_birth_packet_id(b.accepted_step,0,i);v.state="pfc-v3/H2O";
      v.origin="manufactured/authored-neutral-motion-fixture";v.weight=2;
      v.birth_time_s=b.begin_s;v.remaining_flight_s=b.end_s-v.birth_time_s;
      v.position_m={.006,.004,0};v.internal_energy_J=i==0?1e-19:0;
      if(i==1){v.position_m={.004,.002+1e-10,0};v.proper_velocity_m_s={0,-650,0};}
      if(i==2){v.position_m={o.width_m-1e-10,.004,0};v.proper_velocity_m_s={650,0,0};}
      if(i==3){v.position_m={.0042-1e-10,.004,0};v.proper_velocity_m_s={650,0,0};}
      credit.births.push_back(v);
    }
    born.emplace(bs.prepare(credit));
  }
  if(motion_control=="empty-stream-caps") {
    auto bs=outer.birth_scope();m::FullTimeBirthCredit credit;
    credit.arrival_step=b.accepted_step;credit.base_native_id=bs.base_native_id();
    credit.plan_sha256=sha("manufactured/empty-stream-plan");credit.staged_token_sha256=sha("manufactured/empty-stream-token");
    credit.source_before_sha256=bs.source_before_sha256();credit.source_after_sha256=sha("manufactured/empty-stream-source");
    for(unsigned i=0;i<20;++i)credit.streams.push_back({i,"pfc-v3/H2O",m::ArenaCount(0,0)});
    born.emplace(bs.prepare(credit));
  }
  auto trial=born?outer.fork_transport_trial(*born):outer.fork_transport_trial();
  auto flights=trial.birth_flights();
  R.census_sha256 = trial.census_sha256();
  T.stop(R.phases["fork_census"]);
  // ---- producer (untimed): straight drift of the trial census
  const auto &census = trial.census_v7();
  std::vector<m::GapParticle> rows(census.carriers.begin(), census.carriers.end());
  R.census = rows.size();
  const m::PicVector shift{-.5 * o.width_m, -.5 * o.gap_m, 0};
  std::vector<m::GapParticle> ends;
  std::vector<m::FullTimeMotionPath> paths;
  std::vector<m::FullTimeSinkRecord> sinks;
  std::vector<pf_v7::Sink> form_sinks;
  for (const auto&p:rows) {
    double F=dt;for(const auto&flight:flights)if(flight.native_id==p.id)F=flight.flight_s;
    double u2=0;for(double v:p.proper_velocity_m_s)u2+=v*v;
    const double gamma=std::sqrt(1+u2/(m::pic_c*m::pic_c));
    m::PicVector D{};for(int d=0;d<3;++d)D[d]=F*p.proper_velocity_m_s[d]/gamma;
    const auto hit=m::full_time_first_kinetic_hit(o,p.position_m,D);
    m::FullTimeMotionPath pa;pa.native_id=p.id;
    for(int d=0;d<3;++d)pa.field_start_m[d]=p.position_m[d]+shift[d];
    if(hit.surface>=0) {
      m::FullTimeSinkRecord z;z.native_id=p.id;z.surface=static_cast<m::FullTimeSinkSurface>(hit.surface);
      z.kind=(hit.surface==0 || hit.surface==1 || hit.surface==3)?m::FullTimeSinkKind::collected_at_conductor:m::FullTimeSinkKind::escaped_exterior;
      z.axis=hit.axis;z.plane_m=hit.plane_m;z.fraction=hit.fraction;z.segment_m=D;
      for(int d=0;d<3;++d)z.hit_point_m[d]=p.position_m[d]+hit.fraction*D[d];
      z.hit_point_m[z.axis]=z.plane_m;z.incident_proper_velocity_m_s=p.proper_velocity_m_s;
      z.flight_s=F*hit.fraction;z.discarded_flight_s=F-z.flight_s;
      for(int d=0;d<3;++d)pa.field_displacement_m[d]=(z.hit_point_m[d]+shift[d])-pa.field_start_m[d];
      pa.flight_s=z.flight_s;sinks.push_back(z);
      pf_v7::Sink f{};f.id=z.native_id;f.kind=unsigned(z.kind);f.surface=unsigned(z.surface);f.axis=z.axis;
      f.plane_m=z.plane_m;f.fraction=z.fraction;f.flight_s=z.flight_s;f.discarded_flight_s=z.discarded_flight_s;
      for(int d=0;d<3;++d){f.segment_m[d]=z.segment_m[d];f.hit_point_m[d]=z.hit_point_m[d];f.incident_proper_velocity_m_s[d]=z.incident_proper_velocity_m_s[d];}
      form_sinks.push_back(f);
    }else {
      auto e=p;for(int d=0;d<3;++d){e.position_m[d]+=D[d];pa.field_displacement_m[d]=(e.position_m[d]+shift[d])-pa.field_start_m[d];}
      pa.flight_s=F;ends.push_back(e);
    }
    paths.push_back(pa);
  }
  const pf_v7::HallAxialMesh hm{std::uint64_t(o.nx), 6};
  static std::map<std::uint64_t, std::vector<double>> Es, Bs, Js;
  auto &E = Es[hm.n], &B = Bs[hm.n], &Jc = Js[hm.n];
  if (E.empty()) {
    E.assign(hm.ne(), 1e3);
    B.assign(hm.nb(), 0.0);
    Jc.assign(hm.ne(), 0.0);
  }
  m::FullTimeMotionRecordV7 x;
  x.accepted_step = b.accepted_step;
  x.staged_arena_revision = trial.staged_arena_revision();
  x.begin_s = b.begin_s;
  x.end_s = b.end_s;
  x.dt_s = dt;
  x.staged_arena_sha256 = trial.state_prefix();
  x.trial_field_sha256 = pf_v7::trial_field_blocks_sha256(E, B);
  x.mesh_sha256 = b.mesh_sha256;
  x.material_sha256 = b.material_sha256;
  x.species_sha256 = sha("species");
  x.common_clock_sha256 = b.common_clock_sha256;
  x.kinetic_policy_sha256 = hall_driver::kinetic_policy_sha256(hall_driver::declared_kinetic_policy(o));
  const auto K = pf_v7::kinetic_change_J({reinterpret_cast<const pf_v7::Carrier *>(rows.data()), rows.size()},
                                         {reinterpret_cast<const pf_v7::Carrier *>(ends.data()), ends.size()}, form_sinks);
  const auto W = pf_v7::field_work_J(dt, E, Jc);
  x.budget = {K.value, W.value, K.value - W.value, 0, 0, 0, 0};
  x.census_sha256 = R.census_sha256;
  x.current_A_m_sha256 = pf_v7::current_blocks_sha256(Jc);
  x.endpoint_charge_C_sha256 = sha("endpoint-charge");
  x.initial_charge_C_sha256 = sha("initial-charge");
  x.kinetic_policy = {o.width_m, o.gap_m, o.short_z_m, o.feed_z_m, o.lower_foil.left_m, o.lower_foil.right_m,
                      o.lower_foil.bottom_m, o.lower_foil.top_m, o.upper_foil.left_m, o.upper_foil.right_m,
                      o.upper_foil.bottom_m, o.upper_foil.top_m};
  x.native_to_field_m = shift;
  x.neutral_stream_dt_max_s = A.cfg.motion.neutral_stream_dt_max_s;
  x.endpoints = ends;
  x.paths = paths;
  x.sinks=sinks;
  x.births=flights;
  if(born)x.schema=m::motion_v81_schema;
  x.record_sha256 = born?m::full_time_motion_record_v81_sha256(x):m::full_time_motion_record_v7_sha256(x);
  x.trial_E = E;
  x.trial_B = B;
  x.current_A_m = Jc;
  m::FullTimeMotionBinding mb{b, x.record_sha256, b.accepted_field_sha256, sha("f/" + std::to_string(k)),
                              sha("h/" + std::to_string(k)), trial.source_after_sha256()};
  T.start();
  auto moved = trial.prepare_total_field_motion(mb, m::bind_total_field_motion_v7(x));
  if(motion_control=="snapshot-controls" || motion_control=="empty-stream-caps") {
  auto check_small_cap=[&](){bool small=false;try{(void)coupled.motion_snapshot(moved,1u<<28,1);}catch(const m::FullTimeResourceRefusal&e){small=e.kind()==m::FullTimeResourceRefusal::Kind::carrier_cap;}
    if(!small)throw std::runtime_error("motion receipt cap before/after cache must refuse");};
  if(motion_control=="empty-stream-caps")check_small_cap();
  const auto view=coupled.motion_snapshot(moved);
  if(motion_control=="empty-stream-caps")check_small_cap();
  if(!coupled.can_inspect_motion(moved,*view) || coupled.motion_snapshot(moved).get()!=view.get())throw std::runtime_error("motion snapshot cache/live preflight");
  auto copied=*view;if(coupled.can_inspect_motion(moved,copied))throw std::runtime_error("copied snapshot gained authority");
  bool capped=false;try{(void)coupled.motion_snapshot(moved,1);}catch(const m::FullTimeResourceRefusal&){capped=true;}
  if(!capped || !coupled.can_inspect_motion(moved,*view))throw std::runtime_error("snapshot capacity refusal changed live capability");
  auto other=trial.prepare_total_field_motion(mb,m::bind_total_field_motion_v7(x));
  const auto other_view=coupled.motion_snapshot(other);
  if(view==other_view || coupled.can_inspect_motion(other,*view) || coupled.can_inspect_motion(moved,*other_view))throw std::runtime_error("equal motion retry snapshot crossed attempts");
  auto relocated=std::move(moved);
  if(coupled.can_inspect_motion(moved,*view) || !coupled.can_inspect_motion(relocated,*view))throw std::runtime_error("moved motion capability");
  if(motion_control=="same-digest-motion"){const auto prior=gas->candidate_sha256();gas->prepare_keys({},{},{});if(prior!=gas->candidate_sha256())throw std::runtime_error("motion mutation fixture changed digest");}
  coupled.accept_transport(std::move(relocated));
  if(coupled.can_inspect_motion(relocated,*view) || coupled.can_inspect_motion(other,*other_view))throw std::runtime_error("consumed motion retained authority");
  } else {
    if(motion_control=="same-digest-motion"){const auto prior=gas->candidate_sha256();gas->prepare_keys({},{},{});if(prior!=gas->candidate_sha256())throw std::runtime_error("motion mutation fixture changed digest");}
    if(motion_control.starts_with("allocation-motion-"))allocation_probe::fail_after=motion_control=="allocation-motion-count"?std::numeric_limits<long>::max():std::stol(motion_control.substr(18));
    coupled.accept_transport(std::move(moved));
    if(motion_control=="allocation-motion-count"){measured_motion_allocations=std::numeric_limits<long>::max()-allocation_probe::fail_after;allocation_probe::fail_after=-1;}
  }
  for(const auto&[phase,seconds]:coupled.motion_phase_seconds())R.phases["motion_"+phase].wall=seconds;
  T.stop(R.phases["admission_accept"]);
  m::ExposureBinding e2 = e1;
  e2.slot = m::FullTimeSlot::after_transport;
  e2.frozen_physical_stage_s = b.end_s;
  e2.owner_stage_sequence = stages + 1;
  e2.frozen_field_sha256 = mb.field_after_sha256;
  T.start();
  auto s2 = outer.second_exposure(e2);
  if(before_stage)before_stage(coupled,s2,2);
  R.stage2=coupled.run_stage(s2);
  T.stop(R.phases["stage2"]);
  native_gas_common::clone_proof_json(receipt,coupled.stage_receipts());
  if (!publish) return R; // deliberately abandoned common candidate
  if (before_seal)
    before_seal(coupled);
  T.start();
  coupled.seal({mb.field_after_sha256, mb.heat_after_sha256, mb.source_after_sha256, R.stage2->state_after_sha256});
  native_gas_common::clone_proof_json(receipt,coupled.receipt());
  const auto final=coupled.final_candidate_population();
  const auto boundary=coupled.candidate_boundary();
  if(!final || !boundary || boundary->accepted_step!=b.accepted_step+1)throw std::runtime_error("genuine final candidate census/boundary absent");
  if(before_publish)before_publish(coupled);
  allocation_probe::fail_after=0;
  const bool ok=coupled.can_publish() && coupled.try_publish();
  const auto remaining=allocation_probe::fail_after;allocation_probe::fail_after=-1;
  if(!ok || remaining!=0)throw std::runtime_error("common publication refused or allocated");
  T.stop(R.phases["seal_publish"]);
  continuation=coupled.continuation();
  R.boundary_sha256 = arena.full_time_boundary()->accepted_arena_sha256;
  return R;
}
}
