#include "metal_transfer.hpp"
#include <fstream>
#include <iostream>
using namespace metal_transfer;
J checks=J::array(),proofs=J::array(),boundaries=J::array();
void check(bool v,const std::string& name){checks.push_back({{"name",name},{"passed",v}});if(!v)throw std::runtime_error(name);}
template<class F>void rejects(F fn,const std::string& name){bool failed=false;try{fn();}catch(const std::exception&){failed=true;}check(failed,name);}
m::GapPicOptions options(){m::GapPicOptions o;o.nx=o.ny=45;o.width_m=o.gap_m=.009;o.depth_m=.025;o.finite_foils=true;o.finite_axial_extent=true;o.lower_foil={.002,.007,.0018,.002};o.upper_foil={.002,.007,.007,.0072};return o;}
Config config(double internal=0){Config c;c.material_seal="1a1845e04f1f6fda33338d9bd284f5e1ead9cbc4b7747d0c6355e32cb5661c2d";c.initial_prefix="manufactured-native-metal-integration";c.geometry.nx=c.geometry.ny=1;c.geometry.x0=.002;c.geometry.y0=.007;c.geometry.revision="fixed-upper-foil-uniform";c.geometry.electrode="upper";c.internal_per_atom_J=internal;c.manufactured_internal_state=internal!=0;return c;}
struct Case {
 m::GapPic pic{options()};std::unique_ptr<m::PicArena> arena;std::unique_ptr<Adapter> adapter;Config cfg;h::Count initial;
 Case(double T=2500,bool incoming=false,double internal=0):cfg(config(internal)){
  auto spec=Adapter::required_registry_state(internal);initial={8000*.005*.0002*.025/spec.mass_kg,0};std::vector<m::ArenaPacket> packets;
  if(incoming){m::GapParticle p;p.mass_kg=spec.mass_kg;p.charge_C=0;p.weight=1e15;p.position_m={.0045,std::nextafter(.0072,INFINITY),0};p.proper_velocity_m_s={100,-2000,300};auto ids=pic.add_authored_batch({p});packets.push_back({"incident",spec.id,"manufactured-incident", "",{1e15,0},ids});}
  pic.solve_fields();arena=std::make_unique<m::PicArena>(pic,std::vector<m::ArenaStateSpec>{spec},std::vector<m::ArenaChannel>{},packets,cfg.initial_prefix);adapter=std::make_unique<Adapter>(cfg,*arena,initial,T);
 }
 Event event(std::string id,double count0=1e16){Event e;e.id=id;e.state_prefix=adapter->prefix();e.base_revision=adapter->snapshot()->revision;
  m::MetalSurfacePacket p;p.id=id+"/packet";p.count={count0,0};p.pose.position_m={.0045,std::nextafter(.0072,INFINITY),0};p.pose.proper_velocity_m_s={1000,2000,500};e.packets.push_back(p);return e;
 }
};
int main(int argc,char**argv){try{need(argc==2,"report directory");std::string out=argv[1];
 {
  Case c;auto before=c.adapter->checkpoint();auto prefix=c.adapter->prefix();auto native=m::PicArena::fingerprint(c.pic);auto e=c.event("emit");auto p=c.adapter->prepare(e);
  check(c.adapter->checkpoint()==before&&c.adapter->prefix()==prefix&&m::PicArena::fingerprint(c.pic)==native,"prepare leaves canonical metal and native PIC intact");
  auto audit=p.audit();auto eb=D::read(audit["solid_before"]["energy"]),ea=D::read(audit["solid_after"]["energy"]),ep=D::read(audit["outward_packet_energy_J"]);
  check((eb-ea-ep).compare(D(0))==0,"exact heat+native particle reference internal kinetic energy closure");
  check(D::read(audit["actual_state_internal_J"]).compare(D(0))==0,"ground-state internal energy not invented");
  for(int d=0;d<3;++d)check((D::read(audit["solid_after"]["support_impulse"][d])+D::read(audit["gas_momentum"][d])).compare(D(0))==0,"fixed support receives exact proper momentum component"+std::to_string(d));
  auto stale=c.adapter->prepare(c.event("stale"));const auto diag=c.pic.diagnostics();std::vector<double>phi(c.pic.potential().GetData(),c.pic.potential().GetData()+c.pic.potential().Size());
  check(c.adapter->try_publish(p),"joint gas and canonical thermal publish");
  check(c.arena->snapshot()->packets.size()==1&&c.pic.particles().size()==1,"one canonical native gas owner");
  check(c.pic.diagnostics().authored_kinetic_energy_J==diag.authored_kinetic_energy_J,"no authored-reservoir duplicate source debit");
  check(std::equal(phi.begin(),phi.end(),c.pic.potential().GetData()),"neutral source preserves actual solved field bitwise");
  check(!c.adapter->try_publish(p)&&!c.adapter->try_publish(stale),"published and stale trials rejected");
  rejects([&]{(void)c.adapter->prepare(e);},"replayed event/prefix rejected");
  check(c.adapter->try_rollback(p)&&c.adapter->checkpoint()==before&&m::PicArena::fingerprint(c.pic)==native,"rollback restores exact heat inventory native particles fields IDs and ledger");
  check(!c.adapter->try_publish(stale)&&!c.adapter->try_rollback(p),"generation prevents stale resurrection and double rollback");proofs.push_back(audit);
 }
 for(int fault=1;fault<=4;++fault){Case c;auto before=c.adapter->checkpoint();auto native=m::PicArena::fingerprint(c.pic);rejects([&]{(void)c.adapter->prepare(c.event("fault"),fault);},"fault injection"+std::to_string(fault));check(before==c.adapter->checkpoint()&&native==m::PicArena::fingerprint(c.pic),"native and heat intact after fault"+std::to_string(fault));}
 {
  Case c;auto before=c.adapter->checkpoint();auto native=m::PicArena::fingerprint(c.pic);auto p=c.adapter->prepare(c.event("publish-fault"));check(!c.adapter->try_publish(p,true)&&before==c.adapter->checkpoint()&&native==m::PicArena::fingerprint(c.pic),"failure after gas publish rolls back before exposing thermal acceptance");
 }
 {
  Case c(2500,true,1e-20);auto e=c.event("collect");e.packets.clear();m::MetalSurfacePacket packet;packet.emit=false;packet.id="incident";packet.count={1e15,0};e.packets.push_back(packet);auto p=c.adapter->prepare(e);auto a=p.audit();check(D::read(a["actual_native_kinetic_J"]).compare(D(0))<0&&D::read(a["actual_state_internal_J"]).compare(D(0))<0&&D::read(a["actual_state_reference_J"]).compare(D(0))<0,"actual incident kinetic internal and reference returned to metal");check(c.adapter->try_publish(p)&&c.pic.particles().empty()&&c.arena->snapshot()->packets.empty(),"condensation removes actual canonical particle exactly once");check(count(c.adapter->snapshot()->solid_atoms).compare(count(c.initial)+D(1e15))==0,"condensation credits one canonical material inventory");proofs.push_back(a);
 }
 {
  Case c;auto old_count=c.adapter->snapshot()->solid_atoms;auto old_energy=c.adapter->snapshot()->energy;auto p=c.adapter->prepare(c.event("one-atom",1));auto a=p.audit();check(c.adapter->try_publish(p),"sub-ULP single atom source accepted");auto now=c.adapter->snapshot();check(now->solid_atoms.hi==old_count.hi&&now->solid_atoms.lo==-1&&count(now->solid_atoms).compare(count(old_count)-D(1))==0,"compensated canonical inventory retains single atom below high-word ULP");check((old_energy-now->energy).compare(D::read(a["outward_packet_energy_J"]))==0,"exact energy retains atom debit below total-energy ULP");proofs.push_back(a);
  auto over=c.event("overdraw",old_count.hi);auto before=c.adapter->checkpoint();rejects([&]{(void)c.adapter->prepare(over);},"high-word rounded inventory cannot overdraw exact count");check(c.adapter->checkpoint()==before,"exact overdraw leaves owner unchanged");
 }
 for(auto [T,heat]:std::vector<std::pair<double,double>>{{3000,std::numeric_limits<double>::denorm_min()},{298.2,-std::numeric_limits<double>::denorm_min()}}){Case c(T);auto e=c.event("exact-boundary");e.packets.clear();e.other_heat_J=heat;e.heat_owner="manufactured-boundary-energy";auto before=c.adapter->checkpoint();rejects([&]{(void)c.adapter->prepare(e);},"exact boundary rejects heat invisible to rounded total"+std::to_string(T));check(before==c.adapter->checkpoint(),"failed exact boundary is immutable");boundaries.push_back({{"temperature",T},{"heat_J",heat},{"before",before}});}
 {
  Case c;auto e=c.event("tiny-heat");e.packets.clear();e.other_heat_J=1e-300;e.heat_owner="manufactured-tiny-solid-heat";auto old=c.adapter->snapshot()->energy;auto p=c.adapter->prepare(e);check(c.adapter->try_publish(p)&&(c.adapter->snapshot()->energy-old).compare(D(1e-300))==0,"subresolution heat retained exactly without temperature clamp");proofs.push_back(p.audit());
  auto bad=c.event("weighted-underflow",1e-308);auto before=c.adapter->checkpoint();rejects([&]{(void)c.adapter->prepare(bad);},"unrepresentable actual particle weighted energy rejected without source floor");check(before==c.adapter->checkpoint(),"unrepresentable particle refuses whole transfer");
 }
 {
  Case c;auto e=c.event("restart");auto p=c.adapter->prepare(e);check(c.adapter->try_publish(p),"restart fixture publish");c.arena->save(out+"/arena.restart.json");auto record=c.adapter->checkpoint();std::ofstream(out+"/metal.restart.json")<<record.dump(2);
  m::GapPic restored_pic(options());restored_pic.solve_fields();auto restored=m::PicArena::restore(restored_pic,out+"/arena.restart.json",{Adapter::required_registry_state()},{},c.cfg.initial_prefix);Adapter a(c.cfg,*restored,c.initial,2500);a.restore(record);check(a.checkpoint()==record&&a.prefix()==c.adapter->prefix(),"restart binds and reproduces actual native arena and exact metal state");
  auto bad=record;bad["body"]["state"]["energy"]["exponent"]=4;rejects([&]{Adapter fresh(c.cfg,*restored,c.initial,2500);fresh.restore(bad);},"restart energy tamper rejected before acceptance");
  bad=record;bad["body"]["state"]["energy"]=(D::read(bad["body"]["state"]["energy"])+D(.1)).json();bad["sha256"]=hash(bad["body"].dump());rejects([&]{Adapter fresh(c.cfg,*restored,c.initial,2500);fresh.restore(bad);},"rehash cannot hide inconsistent canonical energy journal");
  bad=record;bad["body"]["state"]["support_impulse"][0]=D(1).json();bad["sha256"]=hash(bad["body"].dump());rejects([&]{Adapter fresh(c.cfg,*restored,c.initial,2500);fresh.restore(bad);},"rehash cannot hide support momentum mismatch");
  auto changed=c.cfg;changed.material_seal=std::string(64,'a');rejects([&]{Adapter fresh(changed,*restored,c.initial,2500);},"material identity change rejected");
  auto fresh=c.event("continued");auto p1=c.adapter->prepare(fresh);Event e2=fresh;e2.state_prefix=a.prefix();auto p2=a.prepare(e2);check(p1.audit()["solid_after"]==p2.audit()["solid_after"],"restart next prepared canonical state identical");
 }
 {
  Case c(2500,true);auto event=c.event("actual-native-collection");event.packets.clear();event.dt_s=1e-9;auto transport=c.arena->prepare_transport(event.dt_s);check(transport.audit()["species_sink_receipts"].size()==1,"actual native drift produces owned wall sink receipt");auto before=c.adapter->checkpoint();auto native=m::PicArena::fingerprint(c.pic);auto p=c.adapter->prepare_native_collection(event,std::move(transport));check(c.adapter->checkpoint()==before&&m::PicArena::fingerprint(c.pic)==native,"native sink and thermal preparation remain unpublished");check(c.adapter->try_publish(p)&&c.pic.particles().empty()&&c.adapter->snapshot()->consumed_sinks.size()==1,"real drift absorption and material condensation publish once");check(c.pic.diagnostics().time_s==1e-9,"actual native transport clock accepted with collection");proofs.push_back(p.audit());
  c.arena->save(out+"/collected.arena.json");auto record=c.adapter->checkpoint();m::GapPic pic(options());pic.solve_fields();auto arena=m::PicArena::restore(pic,out+"/collected.arena.json",{Adapter::required_registry_state()},{},c.cfg.initial_prefix);Adapter restored(c.cfg,*arena,c.initial,2500);restored.restore(record);check(restored.checkpoint()==record,"actual consumed sink and exact metal state survive restart");
  auto corrupt=record;corrupt["body"]["state"]["consumed_sinks"]=J::array({"not-a-native-receipt"});corrupt["sha256"]=hash(corrupt["body"].dump());rejects([&]{Adapter fresh(c.cfg,*arena,c.initial,2500);fresh.restore(corrupt);},"rehash cannot replace actual consumed native sink identity");
  check(c.adapter->try_rollback(p)&&c.adapter->checkpoint()==before&&m::PicArena::fingerprint(c.pic)==native,"actual drift/collection rollback restores clock fields particle and material");
 }
 {
  Case c(300);auto e=c.event("cold-real-heat");e.packets.clear();ss304v2::Model model(ss304v2::Qualification::source_conditioned_liquid_3000,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,8000);e.other_heat_J=(count(c.initial)*D(Adapter::required_registry_state().mass_kg)*D(model.evaluate(299).u_model_j_kg)).value();e.heat_owner="manufactured-cooling-check";auto p=c.adapter->prepare(e);check(c.adapter->try_publish(p)&&std::abs(c.adapter->temperature()-299)<1e-9,"actual native heat cools below300K without clamp");proofs.push_back(p.audit());
 }
 {
  Case c;auto e=c.event("two-word-gas");e.packets[0].count={1e16,-1};auto p=c.adapter->prepare(e);check(c.adapter->try_publish(p),"two-word gas source accepted");auto carriers=c.arena->rate_carriers(e.packets[0].id);D total;for(auto& x:carriers)total=total+D(x.weight);check(carriers.size()>1&&total.compare(count(e.packets[0].count))==0,"actual positive native carriers preserve exact two-word gas count");proofs.push_back(p.audit());
 }
 {
  for(auto [start,target]:std::vector<std::pair<double,double>>{{1600,1750},{1750,1600}}){Case c(start);ss304v2::Model model(ss304v2::Qualification::source_conditioned_liquid_3000,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,8000);auto e=c.event("native-phase-crossing");e.packets.clear();D kg=count(c.initial)*D(Adapter::required_registry_state().mass_kg);e.other_heat_J=(kg*D(model.evaluate(target).u_model_j_kg)-c.adapter->snapshot()->energy).value();e.heat_owner="manufactured-phase-integration";auto p=c.adapter->prepare(e);check(c.adapter->try_publish(p)&&std::abs(c.adapter->temperature()-target)<1e-8,"actual native heat crosses fusion branch once "+std::to_string(start));proofs.push_back(p.audit());}
  Case c(1600);rejects([&]{(void)c.adapter->prepare(c.event("cold-vapor-emission"));},"cold material cannot invoke unsupported vapor-emission branch");
 }
 {
  Case c;h::Packet sampled;auto e=c.event("caller-sampled");e.packets={c.adapter->sample("sampled/atom",{1e15,0},{.4,.35,.7,.5,.5},&sampled)};auto p=c.adapter->prepare(e);auto a=p.audit();auto native=e.packets[0].pose;
  check(sampled.kinetic_J==m::pic_kinetic_energy_J(native),"actual heat sampler proper velocity kinetic equals native particle kinetic");
  check(sampled.kinetic_J!=2*h::kb*c.adapter->temperature(),"sampled kinetic energy is not replaced by mean2kT");
  check(D::read(a["actual_native_kinetic_J"]).compare(D(1e15)*D(sampled.kinetic_J))==0,"actual sampled kinetic debited once");
  rejects([&]{(void)c.adapter->sample("invalid-draw",{1,0},{0,.5,.5,.5,.5});},"closed or hidden sampling draws rejected");
  check(c.adapter->try_publish(p),"caller-owned actual sample publishes through native arena");a["sampled_classical_per_atom_J"]=sampled.classical_sample_K_J;a["sampled_relativistic_conversion_J"]=sampled.relativistic_conversion_J;proofs.push_back(a);
 }
 {
  Case c;auto e=c.event("wrong-stage");e.physical_stage_s=1;rejects([&]{(void)c.adapter->prepare(e);},"wrong native physical stage rejected");e=c.event("wrong-source-face");e.packets[0].pose.position_m[1]=.004;rejects([&]{(void)c.adapter->prepare(e);},"interior gas is not an authored wall impact");e=c.event("duplicate-field");e.field_work_J=1;rejects([&]{(void)c.adapter->prepare(e);},"PIC field work cannot be deposited twice");e=c.event("duplicate-latent");e.extra_latent_J=1;rejects([&]{(void)c.adapter->prepare(e);},"latent phase energy cannot be debited twice");auto cfg=c.cfg;cfg.fixed_support=false;rejects([&]{Adapter a(cfg,*c.arena,c.initial,2500);},"moving recoil requires separately admitted work closure");
 }
 {
  auto opts=options();opts.short_z_m=0;opts.feed_z_m=.025;m::GapPic pic(opts);pic.solve_fields();auto cfg=config();cfg.axial_center_m=.0125;auto spec=Adapter::required_registry_state();m::PicArena arena(pic,{spec},{},{},cfg.initial_prefix);h::Count initial{8000*.005*.0002*.025/spec.mass_kg,0};Adapter a(cfg,arena,initial,2500);Event e;e.id="explicit-shifted-z-frame";e.state_prefix=a.prefix();e.packets={a.sample("shifted-z/packet",{1e15,0},{.4,.35,.7,.5,.5})};auto p=a.prepare(e);check(a.try_publish(p)&&pic.particles().front().position_m[2]==.0125,"zero-based native z frame bound explicitly without hidden translation");
 }
 {
  Case c;auto cfg=c.cfg;cfg.geometry.width=.0045;Adapter a(cfg,*c.arena,c.initial,2500);Event e;e.id="wrong-geometry";e.state_prefix=a.prefix();e.base_revision=0;rejects([&]{(void)a.prepare(e);},"thermal geometry must match actual native foil");
 }
 {
  m::GapPic pic(options());auto metal=Adapter::required_registry_state();m::ArenaStateSpec water{"H2O",2,1,0,2.9915e-26,0,0};m::GapParticle p;p.mass_kg=water.mass_kg;p.charge_C=0;p.weight=1e10;p.position_m={.0045,std::nextafter(.0072,INFINITY),0};p.proper_velocity_m_s={0,-1000,0};auto ids=pic.add_authored_batch({p});pic.solve_fields();auto cfg=config();m::PicArena arena(pic,{metal,water},{},{{"water",water.id,"manufactured-fixture","",{p.weight,0},ids}},cfg.initial_prefix);h::Count n{8000*.005*.0002*.025/metal.mass_kg,0};Adapter a(cfg,arena,n,2500);Event e;e.id="foreign-native-sink";e.state_prefix=a.prefix();e.dt_s=1e-9;auto transport=arena.prepare_transport(e.dt_s);rejects([&]{(void)a.prepare_native_collection(e,std::move(transport));},"actual H2O sink cannot become surrogate metal mass or heat");check(arena.can_publish(transport),"refused foreign sink preserves caller transport token");
 }
 {
  Case c;Adapter huge(c.cfg,*c.arena,{1e308,0},2500);Event e;e.id="native-activity-overflow";e.state_prefix=huge.prefix();auto before=huge.checkpoint();rejects([&]{(void)huge.prepare(e);},"nonfinite native heat norm cannot create an infinite acceptance tolerance");check(huge.checkpoint()==before,"unrepresentable heat norm refuses unchanged state");
 }
 {
  m::GapPic pic(options());pic.solve_fields();auto spec=Adapter::required_registry_state();m::ArenaStateSpec electron{"electron",0,0,-1,m::pic_me,0,0};m::ArenaChannel invalid{"manufactured/metal-pretending-water","electron",spec.id,{"electron",spec.id},"","","",""};m::PicArena arena(pic,{spec,electron},{invalid},{},"invalid-metal-channel");auto cfg=config();h::Count n{1e20,0};Adapter a(cfg,arena,n,2500);Event e;e.id="invalid-channel-source";e.state_prefix=a.prefix();e.base_revision=0;rejects([&]{(void)a.prepare(e);},"metal surrogate cannot inherit any registered water or collision channel");
 }
 std::ofstream(out+"/native-tests.json")<<J{{"status","passed"},{"checks",checks.size()},{"cases",checks},{"proofs",proofs},{"boundaries",boundaries},{"scope","uniform1x1 native thermal block; actual canonical native particle arena; no coordinator or full60ns qualification"}}.dump(2);std::cout<<"native metal transfer passed "<<checks.size()<<" checks\n";
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
