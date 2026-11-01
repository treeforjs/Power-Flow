// Actual sealed instantaneous source tokens; native bridge/first-hit qualification
// belongs to the integration adapter and is not implied by these component tests.
void surface_birth_group(Results&R,const std::filesystem::path&out){
 axial::Mesh g(45,4,.009,1);const char*file=std::getenv("PF_SOURCE_MASS_CONFIG");if(!file)throw std::runtime_error("actual source mass config missing");std::ifstream input(file);J bound;input>>bound;
 se::Config cfg{4.3};cfg.electron_mass_kg=bound.at("mass_kg").get<double>();cfg.electron_mass_identity=bound.at("identity").get<std::string>();
 s7::BirthPolicy policy;policy.subcell_extent_m=20e-6;policy.return_extent_m=0;
 for(double F:{0.,-1e6,3e9}){
  auto hv=thermal_view(patch(),1500);auto o=owner<VB>(patch(),1500,*hv,cfg);auto pf=field<VB>(g,*o,1,F);auto actual=view<VB>(pf,g);auto p=prepared<VB>(*o,pf.f);auto e=o->make_epoch(1e-14,"trial-midpoint-field");
  const auto base=o->state_sha256();const double t=e.t+.25*e.dt;
  auto b=o->prepare_surface_birth(p,2.,.6,.35,.1,42,t,event(1),policy);const auto raw=o->sample_birth(p,2.,.6,.35,.1,42,policy);const auto&packet=b.phase().packet;
  bool phase=packet.particle.electron_mass_kg==cfg.electron_mass_kg&&packet.particle.proper_velocity_m_s==raw.packet.particle.proper_velocity_m_s&&s7::detail::ledger_values(packet.particle.ledger)==s7::detail::ledger_values(raw.packet.particle.ledger)&&b.origin_m()==origin&&b.patch_index()==0&&b.birth_time_s()==t&&b.canonical_id()==42&&b.phase().outcome==s7::BirthOutcome::Handoff;
  auto tr=o->prepare_with_surface_births(e,pf.f,actual,std::span<const s7::SurfaceBirth>(&b,1),{{}});const auto&l=packet.particle.ledger;
  bool heat=tr.loads().size()==1&&tr.loads()[0].certified_source_exchange&&tr.loads()[0].source_wall_J==l.wall_heat_J&&tr.loads()[0].source_field_work_J==l.resolved_matching_field_work_J&&tr.loads()[0].source_electron_reference_change_J==l.electron_reference_change_J&&tr.loads()[0].authored_heat_J==0;
  R.record("surface-birth","sealed exact-origin mass/phase-space and single certified heat receipt, F="+std::to_string(F),phase&&heat&&o->state_sha256()==base&&o->preflight(b,pf.f,actual)&&o->preflight(tr,actual)&&tr.surface_births().size()==1,{{"F_V_m",F},{"matching_m",packet.matching_distance_m},{"wall_J",l.wall_heat_J},{"matching_work_J",l.resolved_matching_field_work_J},{"reference_J",l.electron_reference_change_J}});
  double work=0;for(int k=0;k<3;++k)work+=p.flux().stage.total_E[k]*packet.boundary_charge_displacement_C_m[k];
  const double tolerance=256*std::numeric_limits<double>::epsilon()*(std::abs(work)+std::abs(l.resolved_matching_field_work_J));
  R.record("surface-birth","bridge displacement and source matching work conjugate; no source preview current, F="+std::to_string(F),std::abs(work-l.resolved_matching_field_work_J)<=tolerance&&((F>0)==(packet.matching_distance_m>0))&&tr.source_candidates().empty()&&tr.source_handoffs().empty());
  J refused=J::object();auto bad=policy;bad.return_extent_m=20e-6;refused["unqualified_return"]=attempt([&]{o->prepare_surface_birth(p,2.,.6,.35,.1,42,t,event(1),bad);}).threw;
  refused["end_time"]=attempt([&]{o->prepare_surface_birth(p,2.,.6,.35,.1,42,e.t+e.dt,event(1),policy);}).threw;
  refused["earlier_time"]=attempt([&]{o->prepare_surface_birth(p,2.,.6,.35,.1,42,std::nextafter(e.t,-INFINITY),event(1),policy);}).threw;
  refused["NaN_time"]=attempt([&]{o->prepare_surface_birth(p,2.,.6,.35,.1,42,NAN,event(1),policy);}).threw;
  refused["zero_id"]=attempt([&]{o->prepare_surface_birth(p,2.,.6,.35,.1,0,t,event(1),policy);}).threw;
  refused["zero_weight"]=attempt([&]{o->prepare_surface_birth(p,0.,.6,.35,.1,42,t,event(1),policy);}).threw;
  refused["closed_mark"]=attempt([&]{o->prepare_surface_birth(p,2.,0.,.35,.1,42,t,event(1),policy);}).threw;
  auto ev=event(1);ev.energy_owner.clear();refused["empty_energy_owner"]=attempt([&]{o->prepare_surface_birth(p,2.,.6,.35,.1,42,t,ev,policy);}).threw;
  auto other=owner<VB>(patch(),1500,*hv,cfg);refused["foreign_owner"]=!other->preflight(b,pf.f,actual)&&attempt([&]{other->prepare_with_surface_births(e,pf.f,actual,std::span<const s7::SurfaceBirth>(&b,1),{{}});}).threw;
  auto different=field<VB>(g,*o,2,F);refused["foreign_field_token"]=!o->preflight(b,different.f,view<VB>(different,g));
  std::array<s7::SurfaceBirth,2>duplicates{b,b};refused["duplicate_id"]=attempt([&]{o->prepare_with_surface_births(e,pf.f,actual,duplicates,{{}});}).threw;
  refused["duplicate_heat"]=attempt([&]{o->prepare_with_surface_births(e,pf.f,actual,std::span<const s7::SurfaceBirth>(&b,1),{{}},{b.heat()});}).threw;
  bool ok=true;for(const auto&x:refused)ok=ok&&x.get<bool>();R.record("surface-birth","time/lifetime/count/identity/duplicate refusals leave source unchanged, F="+std::to_string(F),ok&&o->state_sha256()==base,refused);
  auto invalid_arrays=pf.a;invalid_arrays.E[0]=std::nextafter(invalid_arrays.E[0],INFINITY);s7::FinalFieldView changed{pf.f->frame(),g,invalid_arrays.E,invalid_arrays.B};R.record("surface-birth","changed final field refuses both surface token and owning trial",!o->preflight(b,pf.f,changed)&&!o->preflight(tr,changed));
  const s7::CandidateThermalClock clock{o->bound_thermal_configuration_signature(),e.t+e.dt,e.revision+1};auto predicted=o->candidate_state_sha256(tr,actual,clock);mfem::ConstantCoefficient zero(0);auto heat_trial=hv->prepare(e,zero,tr.loads(),{{}});bool ready=hv->preflight(heat_trial)&&o->preflight(tr,actual);auto moved=std::move(tr);bool empty=tr.surface_births().empty();bool hot=hv->commit(std::move(heat_trial));o->publish_prevalidated(std::move(moved));
  auto tag="surface-birth-"+std::to_string(static_cast<long long>(F));auto hp=out/(tag+"-thermal.bin"),sp=out/(tag+"-source.bin");hv->save(hp.string());o->save(sp.string());auto restored_heat=thermal_view(patch(),1500);restored_heat->load(hp.string());auto restored=owner<VB>(patch(),1500,*restored_heat,cfg);restored->load(sp.string());
  R.record("surface-birth","projection/publication/restore and stale surface-token rejection, F="+std::to_string(F),ready&&empty&&hot&&o->state_sha256()==predicted&&restored->state_sha256()==predicted&&!o->preflight(b,pf.f,actual)&&!restored->preflight(b,pf.f,actual));
 }
 // Physical guards of the existing conditional model remain executable.
 auto hv=thermal_view(patch(),1500);auto o=owner<VB>(patch(),1500,*hv,cfg);auto pf=field<VB>(g,*o,1,3e9);auto p=prepared<VB>(*o,pf.f);const auto t=pf.f->frame().begin_s;
 auto tiny=policy;tiny.subcell_extent_m=1e-12;auto fast=policy;fast.source_log_current_timescale_s=1e-18;auto gyro=field<VB>(g,*o,2,3e9,25.);auto gp=prepared<VB>(*o,gyro.f);
 R.record("surface-birth","actual linear-extent, source-variation and matching gyro guards remain closed",attempt([&]{o->prepare_surface_birth(p,1.,.6,.35,.1,42,t,event(1),tiny);}).threw&&attempt([&]{o->prepare_surface_birth(p,1.,.6,.35,.1,42,t,event(1),fast);}).threw&&attempt([&]{o->prepare_surface_birth(gp,1.,.6,.35,.1,42,t,event(1),policy);}).threw);
}
