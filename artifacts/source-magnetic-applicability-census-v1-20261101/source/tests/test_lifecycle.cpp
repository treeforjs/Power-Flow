#include "fixture.hpp"
#include "magnetic_diagnostic.hpp"
namespace md=powerflow::magnetic_source_census_v1;
struct Mark {std::size_t port;double un,ut,ua,weight,time;std::uint64_t id;};
namespace powerflow::source_stage_v8 {
struct SourceStageTestAccess {
 static auto origin(const eh::Coupler&h,std::size_t p,double z){return PreparedSourceStage::patch_origin(h,p,z);}
 static auto binding(const eh::Coupler&h,std::size_t p,double z){return PreparedSourceStage::patch_binding(h,p,z);}
};
}
J lifecycle(const char*trace,bool water,double weight,double phi,std::vector<Mark>*captured=nullptr,const std::vector<Mark>*probes=nullptr){
 Fixture f(trace,water,weight,phi);J report={{"water_enabled",water},{"electron_weight",weight},{"electron_weight_bits",f.declared.electron_weight_bits},{"work_function_eV",phi},{"start_s",f.common.begin_s},{"dt_s",f.dt},{"attempts",J::array()},{"source_spectra",J::array()},{"reference_probes",J::array()},{"passed",false}};
 const auto&g=f.hall->mesh();const auto actual=f.hall->model().at(f.hall->field());
 ss::s7::v3::TrialFrame frame;frame.role=ss::s7::v3::FieldRole::TotalMidpoint;frame.accepted_revision=f.hall->revision();frame.trial_stage=1;frame.mesh_revision=1;frame.material_revision=2;
 frame.begin_s=f.common.begin_s;frame.end_s=f.common.end_s;frame.dt_s=f.dt;frame.field_stage_s=frame.begin_s+.5*f.dt;
 frame.arena_hash=f.common.accepted_arena_sha256;frame.field_hash="source-lifecycle/actual-accepted-field/"+actual.total_field_sha256();frame.mesh_hash=f.common.mesh_sha256;frame.material_hash=f.common.material_sha256;frame.species_hash=f.registry.registry_sha256;
 const auto mesh=ss::s7::MeshIdentity::of(g);const auto field=std::make_shared<const ss::s7::Field>(mesh,frame,actual.E(),actual.B(),f.hall->model().dt_max());
 const ss::s7::FinalFieldView view{frame,g,actual.E(),actual.B()};
 std::vector<ss::s7::PreparedPatch>patches;std::map<std::uint32_t,std::size_t>streams;ss::pa::StepInput input;input.step=f.common.accepted_step;
 for(std::size_t i=0;i<f.owners.ports.size();++i){const auto&p=f.owners.ports[i];const auto&r=p.registration().front();const auto&heat=f.hall->thermal().slice(p.slice().body,p.slice().slice);const auto epoch=p.source()->make_epoch(f.dt,frame.field_hash);const auto origin=ss::SourceStageTestAccess::origin(heat,0,r.z_center_m);
  patches.push_back(p.source()->prepare_patch(epoch,0,ss::SourceStageTestAccess::binding(heat,0,r.z_center_m),field,origin,{frame.begin_s,f.dt,1,"source-stage-v8/fixed-patch-centre-quadrature"}));const auto&ex=patches.back();auto spec=md::spectrum(ex.flux(),p.slice().settings.at(0).layer);spec["source_port"]=i;spec["body"]=p.slice().body;spec["slice"]=p.slice().slice;spec["origin_m"]=origin;spec["field_sampling"]="actual field interpolation at patch origin";spec["thermal_sampling"]="actual accepted patch area mean (not pointwise T at origin)";spec["z_center_m"]=r.z_center_m;spec["exposure_key"]=ex.exposure_key();spec["expected_electrons"]=ex.expected_number();spec["log_expected_electrons"]=ex.log_expected_number();report["source_spectra"].push_back(std::move(spec));streams.emplace(*r.electron_stream,i);input.electrons.push_back({*r.electron_stream,ex.expected_number(),ex.log_expected_number(),std::isfinite(ex.log_expected_number()),frame.begin_s,f.dt,frame.end_s});
  if(r.neutral_stream){const auto q=p.source()->release_quotes(epoch).front();input.neutrals.push_back({*r.neutral_stream,q.released,q.log_released,p.source()->inventory(0),frame.begin_s,frame.end_s});}
 }
 auto outer=f.arena->begin_full_time_step(f.common,1000,128);n::ExposureBinding ex{f.common,n::FullTimeSlot::before_transport,f.common.begin_s,.5*f.dt,0,f.common.accepted_field_sha256};auto first_scope=outer.first_exposure(ex);auto first=first_scope.prepare_chemistry_v3();if(!first_scope.try_publish(first))throw std::runtime_error("first diagnostic private chemistry refused");auto scope=outer.birth_scope();auto staged=f.owners.arrival->stage(f.owners.arrival->plan(input),scope.base_native_id());
 unsigned active=0,returned=0,handed=0,refused=0;std::size_t paths=0;
 for(const auto&b:staged.births())if(b.kind==ss::pa::Kind::electron){const auto i=streams.at(b.stream_id);const auto&p=f.owners.ports[i];const auto&prepared=patches[i];const auto&flux=prepared.flux();J row={{"source_port",i},{"body",p.slice().body},{"attempt_id",b.native_id},{"packet",b.packet_id()},{"time_s",b.time_s},{"weight",b.weight},{"F_V_m",flux.spectrum.extraction_field_V_m},{"B_T",flux.stage.total_B_T},{"expected_electrons",prepared.expected_number()},{"emitted",false}};
  if(report["attempts"].size()>=128)throw std::runtime_error("attempt census cap exceeded");
  const auto diagnostic=md::operands(flux,p.slice().settings.at(0).layer,b.marks[0],b.marks[1],b.marks[2]);row["diagnostic"]=diagnostic;
  if(captured)captured->push_back({i,b.marks[0],b.marks[1],b.marks[2],b.weight,b.time_s,b.native_id});
  try{ss::s7::v3::BirthSample sample;sample.physical_weight=b.weight;sample.normal_uniform=b.marks[0];sample.transverse_uniform=b.marks[1];sample.azimuth_uniform=b.marks[2];sample.event_id=b.ordinal+1;sample.canonical_particle_id=b.native_id;sample.count_owner_revision=frame.accepted_revision;sample.count_owner_base_hash=staged.owner_base_sha256();sample.count_owner_token=staged.token();sample.exposure_key=prepared.exposure_key();sample.time_policy=ss::s7::v3::BirthTime::SpecifiedPhysicalTime;sample.specified_birth_s=b.time_s;sample.weight_provenance=ss::s7::v3::WeightProvenance::StagedCanonicalOwner;
   ss::eh::SourceEventIdentity event{"source-lifecycle/"+b.packet_id(),"pfc-v3/e","electrode-electron-reservoir","time-resolved-Lorentz-image",b.packet_id(),b.ordinal+1,b.ordinal+1};auto candidate=p.source()->begin_source(prepared,sample,event);md::compare_begin(diagnostic,candidate.state().native);row["emitted"]=true;row["begin_status"]=int(candidate.state().native.status);row["quantum_magnetic_fraction"]=candidate.state().native.quantum_magnetic_energy_fraction;
   if(candidate.state().native.status==surface_emission::LayerStatus::Active){event.producer+="/advance";++event.sequence;++event.owner_sequence;candidate=p.source()->advance_birth(prepared,candidate,event);}
   if(!p.source()->preflight(candidate,field,view))throw std::runtime_error("actual source candidate lost preflight");const auto&s=candidate.state().native;const auto&e=candidate.exchange();row["end_status"]=int(s.status);row["end_time_s"]=s.current_time_s;row["position_from_surface_m"]=s.position_from_surface_m;row["proper_velocity_m_s"]=s.particle.proper_velocity_m_s;row["paths"]=candidate.paths().size();row["wall_heat_J"]=e.wall_heat_J;row["particle_K_change_J"]=e.particle_K_change_J;row["electron_reference_change_J"]=e.electron_reference_change_J;row["image_energy_change_J"]=e.image_energy_change_J;row["field_work_J"]=e.field_work_J;row["energy_residual_J"]=e.residual_J();paths+=candidate.paths().size();
   if(s.status==surface_emission::LayerStatus::Active)++active;else if(s.status==surface_emission::LayerStatus::Recollected)++returned;else if(s.status==surface_emission::LayerStatus::PICHandedOff)++handed;else throw std::runtime_error("unknown layer status");
  }catch(const std::exception&e){const std::string reason=e.what();if(reason.find("diagnostic")!=std::string::npos)throw;if(reason.find("Unmagnetized MG")!=std::string::npos&&!diagnostic.at("magnetic_guard_would_refuse").get<bool>())throw std::runtime_error("magnetic refusal not explained by preserved operands");++refused;row["refusal"]=reason;}
  report["attempts"].push_back(std::move(row));
 }
 if(probes){
  if(probes->size()>128)throw std::runtime_error("probe census cap exceeded");
  for(const auto&m:*probes){const auto&flux=patches.at(m.port).flux();const auto&policy=f.owners.ports.at(m.port).slice().settings.at(0).layer;
   const auto d=md::operands(flux,policy,m.un,m.ut,m.ua);J row={{"source_port",m.port},{"source_mark_id",m.id},{"source_mark_time_s",m.time},{"diagnostic",d},{"begin_admitted",false},{"scope","same quantiles from manufactured phi2 arrivals; nominal phi4.3 spectrum probe only, NOT a physical arrival"}};
   try{const auto trial=surface_emission::begin_layer(flux,policy,m.weight,m.un,m.ut,m.ua,m.id);md::compare_begin(d,trial.state);row["begin_admitted"]=true;row["begin_status"]=int(trial.state.status);}
   catch(const std::exception&e){const std::string reason=e.what();if(reason.find("diagnostic")!=std::string::npos)throw;if(reason.find("Unmagnetized MG")!=std::string::npos&&!d.at("magnetic_guard_would_refuse").get<bool>())throw std::runtime_error("probe refusal lacks operand proof");row["refusal"]=reason;}
   report["reference_probes"].push_back(std::move(row));
  }
 }
 report["active"]=active;report["returned"]=returned;report["handed_off"]=handed;report["refused"]=refused;report["path_count"]=paths;report["accepted_unchanged"]=f.unchanged();report["all_actual_sample_launches_admitted"]=refused==0;report["has_actual_arrival_samples"]=!report["attempts"].empty();report["common_step_published"]=false;report["passed"]=f.unchanged()&&(phi==2.?report["attempts"].size()==47&&refused==21&&active==6&&handed==20:report["attempts"].empty()&&probes&&report["reference_probes"].size()==47);return report;
}
int main(int argc,char**argv){try{
 if(argc!=4)throw std::runtime_error("TRACE OUT WEIGHT");
 {
  surface_emission::PreparedFlux empty;surface_emission::LayerPolicy policy;bool refused=false;
  try{(void)md::operands(empty,policy,0.,.5,.5);}catch(const std::invalid_argument&){refused=true;}
  if(!refused)throw std::runtime_error("closed mark admitted by bounded census");
  refused=false;try{(void)md::operands(empty,policy,.5,.5,.5);}catch(const std::invalid_argument&){refused=true;}
  if(!refused)throw std::runtime_error("empty spectrum admitted by census");
  empty.spectrum.nodes.resize(8193);refused=false;
  try{(void)md::spectrum(empty,policy);}catch(const std::runtime_error&){refused=true;}
  if(!refused)throw std::runtime_error("spectrum census cap not enforced");
 }
 std::filesystem::create_directories(argv[2]);J result={{"schema","actual-source-applicability-census-v1"},{"scope","read-only original arrival replay and separately labelled nominal-reference probes; no publication or material calibration"},{"cases",J::array()},{"passed",true},{"production_qualified",false}};
 std::vector<Mark> marks;for(bool water:{false,true}){auto r=lifecycle(argv[1],water,std::stod(argv[3]),2.,water?nullptr:&marks);result["passed"]=result["passed"].get<bool>()&&r["passed"].get<bool>();result["cases"].push_back(std::move(r));}
 auto r=lifecycle(argv[1],false,std::stod(argv[3]),4.3,nullptr,&marks);result["passed"]=result["passed"].get<bool>()&&r["passed"].get<bool>();result["cases"].push_back(std::move(r));
 std::ofstream(std::filesystem::path(argv[2])/"lifecycle.json")<<result.dump(2)<<'\n';std::cout<<"census_passed="<<result["passed"]<<" cases="<<result["cases"].size()<<'\n';return result["passed"].get<bool>()?0:1;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
