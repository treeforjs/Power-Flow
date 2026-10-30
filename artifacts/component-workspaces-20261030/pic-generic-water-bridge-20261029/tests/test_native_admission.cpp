#include "generic_water_native_adapter.hpp"
#include <fstream>
#include <iostream>
namespace g=raikou::generic_water;namespace m=raikou::mfem_backend;using J=g::Json;
namespace raikou::mfem_backend {class PicReactionBatch {public:
 static std::vector<double> nodal(GapPic& pic,const std::vector<GapParticle>& particles){mfem::LinearForm rhs(&pic.space_);rhs=0;mfem::Array<int>dofs;mfem::Vector shape;for(const auto&p:particles){mfem::IntegrationPoint point;const int element=pic.locate(p.position_m[0],p.position_m[1],point);const auto*fe=pic.space_.GetFE(element);shape.SetSize(fe->GetDof());fe->CalcShape(point,shape);pic.space_.GetElementDofs(element,dofs);rhs.AddElementVector(dofs,shape*=p.charge_C*p.weight/pic.options_.depth_m);}std::vector<double> out(rhs.Size());for(int i=0;i<rhs.Size();++i)out[i]=rhs[i]*pic.options_.depth_m;return out;}
};}
std::size_t checks=0;void check(bool yes,const char*s){++checks;if(!yes)throw std::runtime_error(s);}
constexpr double amu=1.66053906660e-27,M=18.01056468403*amu,H=1.00782503223*amu;
m::GapPicOptions options(){m::GapPicOptions o;o.nx=o.ny=45;o.width_m=o.gap_m=.009;o.depth_m=.025;o.finite_foils=true;o.finite_axial_extent=true;o.short_z_m=0;o.feed_z_m=.025;o.lower_foil={.002,.007,.0018,.002};o.upper_foil={.002,.007,.007,.0072};return o;}
int main(int argc,char**argv){try{
 if(argc!=4)return 2;
 auto data=std::make_shared<const raikou::kinetics::StrictCollisionDataset>(raikou::kinetics::StrictCollisionDataset::read(argv[1]));J binding;{std::ifstream f(argv[2]);f>>binding;}
 g::Bridge bridge(data,binding.at("dataset_sha256"),binding.at("kernel_manifest_sha256"),M,H,g::SitePolicy::PreserveReactantSites);const auto catalog=g::native_catalog(bridge);const std::string out=argv[3];J rows=J::array();
 std::size_t model=0;for(const auto&s:bridge.specs()){
  const double E=model==0?100:model==1||model==2?6.5:model==3?20:8084.470604111418;
  std::vector<std::vector<double>> draws{s.variates==2?std::vector<double>{.37,.61}:s.variates==3?std::vector<double>{.5,.37,.61}:std::vector<double>{.37,.61,.42,.31}};
  if(model==1){draws.clear();auto w=data->kernel_status(s.channel,s.model).at("model").at("OH_vibrational_relative_intensities").get<std::vector<double>>();double sum=0,prefix=0;for(double x:w)sum+=x;for(double x:w){draws.push_back({(prefix+.5*x)/sum,.37,.61});prefix+=x;}}
  for(std::size_t branch=0;branch<draws.size();++branch)for(double electron_weight:{1.,2.})for(int placement=0;placement<(s.captures_primary?3:1);++placement){
   // Sample attachment only at the physically admitted co-located sites. The
   // two changed native placements deliberately test native guard independence.
   const double x=E*g::qe/(g::me*g::c*g::c),u=g::c*std::sqrt(x*(x+2));
   g::PhysicalPacket e{{"primary","electron",7,electron_weight,g::me,{.6*u,0,.8*u}},{.004,.004,.0125}};
   g::PhysicalPacket n{{"water","H2O",7,1,M,{}},s.captures_primary?e.position_m:g::Vec{.00403,.004,.0125}};
   auto batch=bridge.candidates({e},{n},{{7,1e-12}},g::Scope::CompareSupportedSubset);const auto it=std::find_if(batch.rows.begin(),batch.rows.end(),[&](const auto&a){return a.channel_id==s.channel;});check(it!=batch.rows.end(),"native fixture candidate");
   auto result=bridge.sample(*it,e,n,draws[branch]);auto channel=g::fixed_native_comparison_channel(result,binding.at("kernel_binding_file_sha256"));m::GapPic pic(options());
   m::GapParticle pe,pn;pe.position_m=e.position_m;pe.proper_velocity_m_s=e.packet.proper_velocity_m_s;pe.weight=e.packet.weight;pe.mass_kg=e.packet.mass_kg;
   pn.position_m=n.position_m;pn.proper_velocity_m_s=n.packet.proper_velocity_m_s;pn.weight=n.packet.weight;pn.mass_kg=n.packet.mass_kg;pn.charge_C=0;
   if(placement==1)pn.position_m[0]+=.00003;
   if(placement==2)pn.position_m[2]+=.001;
   const auto ids=pic.add_authored_batch({pe,pn});pic.solve_fields();m::PicArena arena(pic,catalog,{channel},{{"primary","electron","generic-comparison","",{electron_weight,0},{ids[0]}},{"water","H2O","generic-comparison","",{1,0},{ids[1]}}},"generic-water-forced-native-proof-v1");
   const auto before=arena.state_prefix(),native_before=m::PicArena::fingerprint(pic);const auto q_before=m::PicReactionBatch::nodal(pic,pic.particles());
   auto abs_particles=pic.particles();for(auto&particle:abs_particles)particle.charge_C=std::abs(particle.charge_C);const auto absolute_before=m::PicReactionBatch::nodal(pic,abs_particles);const auto count_before=pic.particles().size();
   m::ArenaPhase phase;phase.state_prefix=before;phase.mesh_prefix=std::string(64,'a');phase.clock_prefix=std::string(64,'b');
   m::ArenaProvenance provenance{m::ArenaProvenance::Kind::ConditionalPhysicalData,result.data_sha256,result.kernel_sha256,binding.at("kernel_binding_file_sha256"),result.model,"forced_kernel_integration_only",before,0};
   auto event=g::native_event(result,"comparison/"+std::to_string(model)+"/"+std::to_string(branch),"primary","water",{electron_weight,0},phase,provenance);
   J row={{"model",s.model},{"OH_level",result.OH_level},{"products",result.products.size()},{"capture",s.captures_primary},{"placement",placement},{"electron_weight",electron_weight},{"native_scope","forced_single_fixed_product_registry"}};
   bool accepted=false;try{auto token=arena.prepare(event);check(arena.state_prefix()==before&&m::PicArena::fingerprint(pic)==native_before,"native prepare mutation");check(arena.try_publish(token),"native publish");accepted=true;
    // Exact physical charge closure is per reactant site; changed assembly
    // ordering need only satisfy its own finite-operation roundoff bound.
    int qp=0,qt=0;for(const auto&p:result.products)(p.at_primary?qp:qt)+=p.state.q;
    check(s.captures_primary?(pe.position_m==pn.position_m&&qp+qt==-1):(qp==-1&&qt==0),"exact per-site reaction charge balance");
    const auto after_q=m::PicReactionBatch::nodal(pic,pic.particles());auto abs_after_particles=pic.particles();for(auto&particle:abs_after_particles)particle.charge_C=std::abs(particle.charge_C);const auto absolute_after=m::PicReactionBatch::nodal(pic,abs_after_particles);
    const long double operations=count_before+pic.particles().size()+12,epsilon=std::numeric_limits<double>::epsilon();const long double gamma=operations*epsilon/(1-operations*epsilon);long double maximum=0,maximum_bound=0,ratio=0;std::size_t nonzero=0;
    for(std::size_t i=0;i<q_before.size();++i){const long double delta=std::abs(static_cast<long double>(after_q[i])-q_before[i]);const long double bound=gamma*(static_cast<long double>(absolute_before[i])+absolute_after[i])/(1-gamma);check(delta<=bound,"native Q1 assembly exceeds operation roundoff bound");maximum=std::max(maximum,delta);maximum_bound=std::max(maximum_bound,bound);if(delta>0){++nonzero;ratio=std::max(ratio,delta/bound);}}
    row["native_node_delta_Linf_C"]=static_cast<double>(maximum);row["native_node_roundoff_bound_Linf_C"]=static_cast<double>(maximum_bound);row["native_node_worst_bound_fraction"]=static_cast<double>(ratio);row["native_nonzero_roundoff_nodes"]=nonzero;
    if(s.captures_primary)check(!arena.snapshot()->packets.contains("primary")&&(electron_weight>1?arena.snapshot()->packets.contains(event.untouched_primary_id):event.untouched_primary_id.empty()),"capture with full retirement or positive untouched remainder");
    const auto path=out+"/restart-"+std::to_string(model)+"-"+std::to_string(branch)+"-w"+std::to_string(int(electron_weight))+".json";arena.save(path);m::GapPic fresh(options());auto restored=m::PicArena::restore(fresh,path,catalog,{channel},"generic-water-forced-native-proof-v1");check(restored->state_prefix()==arena.state_prefix()&&m::PicArena::fingerprint(fresh)==m::PicArena::fingerprint(pic),"variable product restart exact");
    check(arena.try_rollback(token)&&arena.state_prefix()==before&&m::PicArena::fingerprint(pic)==native_before,"native rollback exact");row["restart_exact"]=true;row["rollback_exact"]=true;
   }catch(const std::exception&error){if(accepted)throw std::runtime_error(s.model+" native accepted check: "+error.what());row["reason"]=error.what();check(arena.state_prefix()==before&&m::PicArena::fingerprint(pic)==native_before,"native refusal mutation");}
   check(accepted==(placement==0),"native admission policy mismatch");row["accepted"]=accepted;rows.push_back(row);
  }
  ++model;
 }
 std::ofstream(out+"/native-admission.json")<<J{{"status","passed"},{"assertions",checks},{"records",rows},{"native_source_identity",m::PicArena::source_identity()},{"competing_native_completion_integrated",false},{"physical_time_advanced",false}}.dump(2)<<'\n';
 std::ofstream(out+"/mechanism-contract.json")<<g::mechanism_contract(bridge,binding.at("kernel_binding_file_sha256")).dump(2)<<'\n';std::cout<<checks<<" native admission assertions across "<<rows.size()<<" records\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
