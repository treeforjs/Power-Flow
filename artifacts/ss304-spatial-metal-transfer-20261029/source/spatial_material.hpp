#pragma once
#include "dyadic.hpp"
#include "material/ss304_thermal_v2.hpp"
#include <openssl/sha.h>
#include <array>
#include <algorithm>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <vector>

// A material/thermal owner only. There is no gas count debit, native ID allocator,
// collision law, accepted temperature cache, or second thermal reservoir here.
namespace spatial_metal {
using metal_transfer::D;
using metal_transfer::need;
using J=nlohmann::json;
inline constexpr const char* material_seal="1a1845e04f1f6fda33338d9bd284f5e1ead9cbc4b7747d0c6355e32cb5661c2d";
inline std::string hash(const std::string& s){unsigned char out[SHA256_DIGEST_LENGTH];SHA256(reinterpret_cast<const unsigned char*>(s.data()),s.size(),out);std::ostringstream o;for(auto x:out)o<<std::hex<<std::setfill('0')<<std::setw(2)<<unsigned(x);return o.str();}
inline bool seal(const std::string& s){return s.size()==64&&s.find_first_not_of("0123456789abcdef")==std::string::npos;}
inline D zero(){return D(0.);}
inline bool eq(const D&a,const D&b){return a.compare(b)==0;}
inline bool finite(double x){return std::isfinite(x);}
enum class Face {Left=0,Right=1,Bottom=2,Top=3};
struct Patch {std::string id;Face face;double lo=0,hi=0,zlo=0,zhi=0;};
struct Config {
 std::vector<double> x,y;
 double zlo=0,zhi=0,rho_reference=8000,max_fraction_loss=0,max_fraction_gain=0;
 std::vector<Patch> patches;
 std::string material_sha256=material_seal,mesh_sha256,initial_common_sha256;
 bool fixed_support=false,conditional_congruent_neutral=false;
 // Explicit fixture gate. A production port must certify actual typed births.
 bool permit_native_genesis_fixture=false;
 std::size_t maximum_cells=4096,maximum_substeps=10000,maximum_records=4096;
};
struct Cell {D mass_kg,energy_J;std::array<D,3>support_impulse;};
struct State {
 std::vector<Cell> cells;
 std::set<std::string> consumed_receipts;
 std::map<std::string,std::uint64_t> heat_watermarks;
 std::uint64_t revision=0,accepted_step=0;
 double time_s=0;
 J journal=J::array();
};
// This binding is owned by the outer coordinator. dt is authoritative; it is
// never reconstructed by subtracting the two absolute endpoint doubles.
struct Binding {
 std::uint64_t accepted_step=0;
 double begin_s=0,end_s=0,dt_s=0;
 std::string material_prefix,common_sha256,field_sha256,mesh_sha256;
 std::string gas_base_sha256,gas_trial_sha256,source_sha256,native_source_sha256;
};
struct Heat {std::string id,owner;std::uint64_t sequence=0;std::size_t cell=0;D energy_J;};
struct Step {Binding binding;std::vector<Heat> heat;double field_work_J=0,extra_latent_J=0;};
// GasPort is an external native-owner attestation boundary, not a parser that
// promotes arbitrary JSON to an actual source/sink. Root must implement it from
// its typed private outer trial and validate full native IDs/counts/kinematics.
// The supplied v6.2 test port is deliberately fixture-only for source births.
class GasPort {
 public:
 virtual ~GasPort()=default;
 virtual const J& certified_batch(const Binding&)const=0;
};
class Owner;
class Prepared {
 friend class Owner;
 Owner* owner_=nullptr;
 std::uint64_t generation_=0;
 std::shared_ptr<const State> before_,after_;
 J audit_;
 enum class Status {Staged,Published,RolledBack};Status status_=Status::Staged;
 public:
 Prepared()=default;Prepared(Prepared&&)noexcept=default;Prepared&operator=(Prepared&&)noexcept=default;
 Prepared(const Prepared&)=delete;Prepared&operator=(const Prepared&)=delete;
 const J& audit()const{return audit_;}
};
class Owner {
 Config cfg_;
 ss304v2::Model law_;
 ss304v2::NeutralVapor vapor_;
 std::vector<D> initial_mass_;
 std::vector<double> initial_T_;
 double initial_time_=0,atom_mass_=0,reference_=0;
 std::string signature_;
 std::shared_ptr<const State> accepted_;
 std::uint64_t generation_=0;
 struct Edge {std::size_t a,b;double area_over_distance;};
 std::vector<Edge> edges_;
 std::size_t nx()const{return cfg_.x.size()-1;}
 std::size_t ny()const{return cfg_.y.size()-1;}
 const Patch& patch(const std::string&id)const{for(auto& p:cfg_.patches)if(p.id==id)return p;throw std::invalid_argument("unknown physical surface patch");}
 static J binding_json(const Binding&b){return {{"accepted_step",b.accepted_step},{"begin_s",b.begin_s},{"end_s",b.end_s},{"dt_s",b.dt_s},{"material_prefix",b.material_prefix},{"common_sha256",b.common_sha256},{"field_sha256",b.field_sha256},{"mesh_sha256",b.mesh_sha256},{"gas_base_sha256",b.gas_base_sha256},{"gas_trial_sha256",b.gas_trial_sha256},{"source_sha256",b.source_sha256},{"native_source_sha256",b.native_source_sha256}};}
 static Binding read_binding(const J&j){need(j.size()==12,"binding schema");Binding b;b.accepted_step=j.at("accepted_step");b.begin_s=j.at("begin_s");b.end_s=j.at("end_s");b.dt_s=j.at("dt_s");b.material_prefix=j.at("material_prefix");b.common_sha256=j.at("common_sha256");b.field_sha256=j.at("field_sha256");b.mesh_sha256=j.at("mesh_sha256");b.gas_base_sha256=j.at("gas_base_sha256");b.gas_trial_sha256=j.at("gas_trial_sha256");b.source_sha256=j.at("source_sha256");b.native_source_sha256=j.at("native_source_sha256");return b;}
 static J heat_json(const Heat&h){return {{"id",h.id},{"owner",h.owner},{"sequence",h.sequence},{"cell",h.cell},{"energy_J",h.energy_J.json()}};}
 static Heat read_heat(const J&j){need(j.size()==5,"heat schema");Heat h;h.id=j.at("id");h.owner=j.at("owner");h.sequence=j.at("sequence");h.cell=j.at("cell");h.energy_J=D::read(j.at("energy_J"));return h;}
 J state_json(const State&s,bool journal=true)const{J cells=J::array();for(auto&c:s.cells){J p=J::array();for(auto&v:c.support_impulse)p.push_back(v.json());cells.push_back({{"mass_kg",c.mass_kg.json()},{"energy_J",c.energy_J.json()},{"support_impulse",p}});}J j={{"cells",cells},{"consumed_receipts",s.consumed_receipts},{"heat_watermarks",s.heat_watermarks},{"revision",s.revision},{"accepted_step",s.accepted_step},{"time_s",s.time_s}};if(journal)j["journal"]=s.journal;return j;}
 double temperature(const Cell&c)const{
  need(c.mass_kg.compare(zero())>0,"nonpositive material stock needs geometry handoff");
  need(c.energy_J.compare(c.mass_kg*D(law_.min_u()))>=0&&c.energy_J.compare(c.mass_kg*D(law_.max_u()))<=0,"local exact energy outside298.2..3000K");
  return law_.temperature_from_u(static_cast<double>(c.energy_J.ld()/c.mass_kg.ld()));
 }
 void admit_cell(const Cell&c,std::size_t i)const{
  need(c.mass_kg.compare(initial_mass_[i]*(D(1.)-D(cfg_.max_fraction_loss)))>=0,"local recession budget exhausted; geometry handoff required");
  need(c.mass_kg.compare(initial_mass_[i]*(D(1.)+D(cfg_.max_fraction_gain)))<=0,"local deposition budget exhausted; geometry handoff required");
  (void)temperature(c);
 }
 static std::size_t interval(const std::vector<double>& edges,double q){need(finite(q)&&q>=edges.front()&&q<=edges.back(),"surface tangent outside solid");return std::min(std::size_t(std::upper_bound(edges.begin(),edges.end(),q)-edges.begin()-1),edges.size()-2);}
 std::size_t cell_on_patch(const Patch&p,const std::array<double,3>&x)const{
  for(double v:x)need(finite(v),"nonfinite surface point");
  need(x[2]>=p.zlo&&x[2]<=p.zhi,"surface axial point outside patch");
  const bool vertical=p.face==Face::Left||p.face==Face::Right;
  const double tangent=vertical?x[1]:x[0],normal=vertical?x[0]:x[1];
  const double boundary=p.face==Face::Left?cfg_.x.front():p.face==Face::Right?cfg_.x.back():p.face==Face::Bottom?cfg_.y.front():cfg_.y.back();
  // Native sources can start one nextafter outside a face. Collection is at
  // the actual first-hit point. No broad coordinate-scale epsilon admission.
  need(normal==boundary||normal==std::nextafter(boundary,INFINITY)||normal==std::nextafter(boundary,-INFINITY),"receipt is not on actual material face");
  need(tangent>=p.lo&&tangent<=p.hi,"surface tangent outside selected physical patch");
  auto ix=vertical?(p.face==Face::Left?0:nx()-1):interval(cfg_.x,x[0]);
  auto iy=vertical?interval(cfg_.y,x[1]):(p.face==Face::Bottom?0:ny()-1);
  return iy*nx()+ix;
 }
 std::array<double,3> normal(Face f)const{return f==Face::Left?std::array<double,3>{-1,0,0}:f==Face::Right?std::array<double,3>{1,0,0}:f==Face::Bottom?std::array<double,3>{0,-1,0}:std::array<double,3>{0,1,0};}
 J conduct(State&next,double dt)const{
  // Cell-centered finite-volume conduction. Each edge flux is one exact binary
  // currency subtracted from one cell and added to its neighbor. Temperature is
  // a derived projection only. Stability uses NEW material mass after transfer.
  // Conservative global bounds: cp>=400 J/kg/K and k<=100 W/m/K over sealed v2.
  std::vector<double> gmax(next.cells.size());for(auto&e:edges_){double v=100*e.area_over_distance;gmax[e.a]+=v;gmax[e.b]+=v;}
  double maxratio=0;for(std::size_t i=0;i<next.cells.size();++i){double c=next.cells[i].mass_kg.value()*400;need(finite(c)&&c>0,"unrepresentable conduction capacity");maxratio=std::max(maxratio,dt*gmax[i]/(.2*c));}
  need(finite(maxratio)&&maxratio<=double(cfg_.maximum_substeps),"conduction substep resource/stability cap");
  const std::size_t n=std::max(std::size_t(1),std::size_t(std::ceil(maxratio)));
  need(edges_.empty()||n<=16384/edges_.size(),"total spatial flux-record resource cap");
  D elapsed;const D usual=D(dt)*D(1./double(n));J rows=J::array();
  for(std::size_t k=0;k<n;++k){D duration=k+1==n?D(dt)-elapsed:usual;elapsed=elapsed+duration;need(duration.compare(zero())>0,"unresolved conduction duration");std::vector<double>T;for(auto&c:next.cells)T.push_back(temperature(c));std::vector<D> delta(next.cells.size());J fluxes=J::array();
   for(auto&e:edges_){auto a=law_.evaluate(T[e.a]),b=law_.evaluate(T[e.b]);need(a.apparent_heat_capacity_j_kg_k>=400&&b.apparent_heat_capacity_j_kg_k>=400&&a.k_w_m_k<=100&&b.k_w_m_k<=100,"sealed material conduction bound violated");double g=(2*a.k_w_m_k*b.k_w_m_k/(a.k_w_m_k+b.k_w_m_k))*e.area_over_distance;need(finite(g)&&g>0,"conduction conductance");D flux=duration*D(g)*D(T[e.a]-T[e.b]);delta[e.a]=delta[e.a]-flux;delta[e.b]=delta[e.b]+flux;fluxes.push_back({{"a",e.a},{"b",e.b},{"conductance_W_K",g},{"temperature_difference_K",T[e.a]-T[e.b]},{"energy_J",flux.json()}});}
   for(std::size_t i=0;i<next.cells.size();++i){next.cells[i].energy_J=next.cells[i].energy_J+delta[i];admit_cell(next.cells[i],i);}
   rows.push_back({{"duration_s",duration.json()},{"fluxes",fluxes}});
  }
  need(eq(elapsed,D(dt)),"conduction does not span authoritative dt");return rows;
 }
 Prepared prepare_batch(const Step&s,const J&batch){
  const auto&b=s.binding;need(b.material_prefix==prefix()&&b.accepted_step==accepted_->accepted_step&&b.begin_s==accepted_->time_s&&finite(b.dt_s)&&b.dt_s>0&&finite(b.end_s)&&b.end_s==b.begin_s+b.dt_s&&b.end_s>b.begin_s,"stale material/common clock binding");
  for(auto*x:{&b.common_sha256,&b.field_sha256,&b.mesh_sha256,&b.gas_base_sha256,&b.gas_trial_sha256,&b.source_sha256,&b.native_source_sha256})need(seal(*x),"missing common owner identity");
  need(b.mesh_sha256==cfg_.mesh_sha256&&s.field_work_J==0&&s.extra_latent_J==0,"mesh mismatch or duplicate field/latent debit");
  need(generation_<UINT64_MAX&&accepted_->revision<UINT64_MAX&&accepted_->accepted_step<UINT64_MAX,"material generation exhaustion");
  need(accepted_->journal.size()<cfg_.maximum_records,"material journal entry cap");
  need(batch.at("binding")==binding_json(b)&&batch.at("material_sha256")==cfg_.material_sha256,"native receipt batch identity");
  std::string scope=batch.at("scope");need(scope=="full-time-outer-native-owner"||(cfg_.permit_native_genesis_fixture&&scope=="v6.2-native-genesis-source-and-first-hit-fixture"),"unqualified native source receipt scope");
  need(batch.at("packets").is_array()&&batch.at("packets").size()+s.heat.size()<=cfg_.maximum_records,"transfer resource cap");
  auto next=std::make_shared<State>(*accepted_);std::vector<D> outward_mass(next->cells.size()),outward_energy(next->cells.size()),heat(next->cells.size());J mapped=J::array();std::set<std::uint64_t>ids;std::set<std::string>packet_ids;
  for(const auto&r:batch.at("packets")){
   const std::string id=r.at("receipt_id"),pid=r.at("packet_id"),patchid=r.at("patch");bool emit=r.at("emit");
   need(!id.empty()&&id.size()<=256&&!pid.empty()&&pid.size()<=256&&patchid.size()<=256&&packet_ids.insert(pid).second&&next->consumed_receipts.insert("gas/"+id).second,"duplicate/stale native surface receipt");
   need(r.at("species")=="SS304_CONGRUENT_NEUTRAL"&&r.at("mass_per_atom_kg")==atom_mass_&&r.at("internal_per_atom_J")==0&&r.at("reference_per_atom_J")==reference_,"unsupported metal chemistry/reference");
   need(r.at("full_packet")==true&&r.at("native_removed")==!emit,"source birth/first-hit ownership mismatch");
   const auto x=r.at("position_m").get<std::array<double,3>>(),u=r.at("proper_velocity_m_s").get<std::array<double,3>>();const auto&p=patch(patchid);std::size_t cell=cell_on_patch(p,x);auto norm=normal(p.face);double un=0;for(int d=0;d<3;++d){need(finite(u[d]),"nonfinite native momentum");un+=u[d]*norm[d];}need(finite(un)&&(emit?un>=0:un<=0),"wrong native surface orientation");
   double when=r.at("time_s");need(finite(when)&&when>=b.begin_s&&when<=b.end_s,"surface event outside common interval");
   if(emit){double t=temperature(accepted_->cells[cell]);need(t>=1800&&r.at("sampled_temperature_K")==t,"source does not use local accepted material temperature");}
   else need(r.at("first_hit")==true,"collection lacks actual native first-hit receipt");
   D atoms,energy;std::array<D,3>momentum;need(!r.at("carriers").empty(),"empty native surface packet");
   for(const auto&c:r.at("carriers")){std::uint64_t nid=c.at("native_id");double w=c.at("weight"),K=c.at("kinetic_per_atom_J");need(nid>0&&ids.insert(nid).second&&finite(w)&&w>0&&finite(K)&&K>=0,"invalid/duplicate native carrier");atoms=atoms+D(w);energy=energy+D(w)*(D(K)+D(reference_));for(int d=0;d<3;++d)momentum[d]=momentum[d]+D(w)*D(atom_mass_)*D(u[d]);}
   need(eq(atoms,D::read(r.at("exact_atoms"))),"native positive-carrier sum differs from canonical gas count");D sign(emit?1.:-1.);D mass=sign*atoms*D(atom_mass_);energy=sign*energy;
   outward_mass[cell]=outward_mass[cell]+mass;outward_energy[cell]=outward_energy[cell]+energy;for(int d=0;d<3;++d)next->cells[cell].support_impulse[d]=next->cells[cell].support_impulse[d]-sign*momentum[d];
   mapped.push_back({{"receipt_id",id},{"cell",cell},{"outward_mass_kg",mass.json()},{"outward_energy_J",energy.json()}});
  }
  J heats=J::array();for(auto&h:s.heat){need(h.cell<next->cells.size()&&!h.id.empty()&&h.id.size()<=256&&!h.owner.empty()&&h.owner.size()<=256&&h.owner!="metal-packet"&&h.sequence>next->heat_watermarks[h.owner]&&next->consumed_receipts.insert("heat/"+h.id).second,"invalid/duplicate canonical heat receipt");next->heat_watermarks[h.owner]=h.sequence;heat[h.cell]=heat[h.cell]+h.energy_J;heats.push_back(heat_json(h));}
  J before=diagnostics(*accepted_);for(std::size_t i=0;i<next->cells.size();++i){auto&c=next->cells[i];c.mass_kg=c.mass_kg-outward_mass[i];c.energy_J=c.energy_J+heat[i]-outward_energy[i];admit_cell(c,i);}
  J after_transfer=diagnostics(*next);J fluxes=conduct(*next,b.dt_s);
  next->time_s=b.end_s;++next->accepted_step;++next->revision;
  D dmass,denergy,packetmass,packetenergy,Q;for(std::size_t i=0;i<next->cells.size();++i){dmass=dmass+next->cells[i].mass_kg-accepted_->cells[i].mass_kg;denergy=denergy+next->cells[i].energy_J-accepted_->cells[i].energy_J;packetmass=packetmass+outward_mass[i];packetenergy=packetenergy+outward_energy[i];Q=Q+heat[i];}
  need(eq(dmass+packetmass,zero())&&eq(denergy+packetenergy-Q,zero()),"exact global material transfer/conduction closure");
  J row={{"binding",binding_json(b)},{"certified_batch_sha256",hash(batch.dump())},{"heat",heats},{"mapping",mapped},{"conduction",fluxes},{"before",state_json(*accepted_,false)},{"after",state_json(*next,false)}};constexpr std::size_t max_journal_bytes=8*1024*1024;auto previous_bytes=accepted_->journal.dump().size(),row_bytes=row.dump().size();need(previous_bytes<=max_journal_bytes&&row_bytes<=max_journal_bytes-previous_bytes,"bounded material journal bytes; anchor required");next->journal.push_back(row);
  Prepared p;p.owner_=this;p.generation_=generation_;p.before_=accepted_;p.after_=next;p.audit_={{"signature",signature_},{"row",row},{"before_diagnostics",before},{"after_transfer_diagnostics",after_transfer},{"after_diagnostics",diagnostics(*next)},{"outward_mass_kg",packetmass.json()},{"outward_energy_J",packetenergy.json()},{"external_heat_J",Q.json()},{"exact_mass_residual",(dmass+packetmass).json()},{"exact_energy_residual",(denergy+packetenergy-Q).json()},{"fixed_support_work_J",0},{"extra_latent_J",0},{"field_work_J",0}};return p;
 }
 public:
 Owner(Config cfg,std::vector<double> initial_T,double initial_time=0):cfg_(std::move(cfg)),law_(ss304v2::Qualification::source_conditioned_liquid_3000,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,cfg_.rho_reference),vapor_(law_),initial_T_(std::move(initial_T)),initial_time_(initial_time){
  need(cfg_.x.size()>=2&&cfg_.y.size()>=2&&cfg_.x.size()<=4097&&cfg_.y.size()<=4097,"bounded spatial mesh required");need(nx()*ny()<=cfg_.maximum_cells&&cfg_.maximum_cells<=4096&&cfg_.maximum_records>0&&cfg_.maximum_records<=4096&&cfg_.maximum_substeps>0&&cfg_.maximum_substeps<=10000,"bounded spatial resources");
  for(auto*edges:{&cfg_.x,&cfg_.y})for(std::size_t i=0;i<edges->size();++i)need(finite((*edges)[i])&&(i==0||(*edges)[i]>(*edges)[i-1]),"strict finite cell edges");
  need(finite(cfg_.zlo)&&finite(cfg_.zhi)&&cfg_.zhi>cfg_.zlo&&finite(initial_time)&&finite(cfg_.rho_reference)&&cfg_.rho_reference>0,"spatial material geometry/time");
  need(cfg_.max_fraction_loss>=0&&cfg_.max_fraction_loss<1&&cfg_.max_fraction_gain>=0&&cfg_.max_fraction_gain<1&&cfg_.fixed_support&&cfg_.conditional_congruent_neutral&&cfg_.material_sha256==material_seal&&seal(cfg_.mesh_sha256)&&seal(cfg_.initial_common_sha256),"explicit sealed material/recession/fixed-support scenario required");
  need(initial_T_.size()==nx()*ny(),"one initial temperature per thermal cell required");
  atom_mass_=ss304v2::NeutralVapor::molar_mass/ss304v2::NeutralVapor::avogadro;reference_=atom_mass_*vapor_.atomic_u0();auto s=std::make_shared<State>();s->time_s=initial_time;
  const double depth=cfg_.zhi-cfg_.zlo;for(std::size_t y=0;y<ny();++y)for(std::size_t x=0;x<nx();++x){D m=D(cfg_.rho_reference)*D(cfg_.x[x+1]-cfg_.x[x])*D(cfg_.y[y+1]-cfg_.y[y])*D(depth);need(m.value()>0,"unrepresentable initial cell mass");initial_mass_.push_back(m);Cell c;c.mass_kg=m;c.energy_J=m*D(law_.evaluate(initial_T_[y*nx()+x]).u_model_j_kg);s->cells.push_back(c);if(x+1<nx())edges_.push_back({y*nx()+x,y*nx()+x+1,(cfg_.y[y+1]-cfg_.y[y])*depth/((cfg_.x[x+2]-cfg_.x[x])/2)});if(y+1<ny())edges_.push_back({y*nx()+x,(y+1)*nx()+x,(cfg_.x[x+1]-cfg_.x[x])*depth/((cfg_.y[y+2]-cfg_.y[y])/2)});}
  std::set<std::string>pids;J patches=J::array();for(auto&p:cfg_.patches){need(!p.id.empty()&&p.id.size()<=256&&pids.insert(p.id).second&&int(p.face)>=0&&int(p.face)<=3&&finite(p.lo)&&finite(p.hi)&&finite(p.zlo)&&finite(p.zhi)&&p.hi>p.lo&&p.zhi>p.zlo&&p.zlo>=cfg_.zlo&&p.zhi<=cfg_.zhi,"finite distinct surface patches");const auto&e=(p.face==Face::Left||p.face==Face::Right)?cfg_.y:cfg_.x;need(p.lo>=e.front()&&p.hi<=e.back(),"patch outside thermal face");patches.push_back({p.id,int(p.face),p.lo,p.hi,p.zlo,p.zhi});}
  J config={{"schema","spatial-material-v1"},{"compiled_source",SPATIAL_MATERIAL_SOURCE_SHA256},{"native_receipt_contract","full-time-outer-r1;v6.2-fixture-explicit"},{"material",cfg_.material_sha256},{"mesh",cfg_.mesh_sha256},{"initial_common",cfg_.initial_common_sha256},{"x",cfg_.x},{"y",cfg_.y},{"z",{cfg_.zlo,cfg_.zhi}},{"rho",cfg_.rho_reference},{"recession",{cfg_.max_fraction_loss,cfg_.max_fraction_gain}},{"patches",patches},{"initial_T",initial_T_},{"initial_time",initial_time_},{"fixture",cfg_.permit_native_genesis_fixture},{"resources",{cfg_.maximum_cells,cfg_.maximum_substeps,cfg_.maximum_records}}};signature_=hash(config.dump());accepted_=s;
 }
 Owner(const Owner&)=delete;Owner&operator=(const Owner&)=delete;
 const std::shared_ptr<const State>& snapshot()const{return accepted_;}
 std::string prefix()const{return hash(signature_+state_json(*accepted_).dump());}
 std::string signature()const{return signature_;}
 double atom_mass()const{return atom_mass_;}double reference_per_atom()const{return reference_;}
 double cell_temperature(std::size_t i)const{return temperature(accepted_->cells.at(i));}
 J diagnostics(const State&s)const{J rows=J::array();for(std::size_t i=0;i<s.cells.size();++i){const auto&c=s.cells[i];double T=temperature(c);auto v=law_.evaluate(T);rows.push_back({{"cell",i},{"T_K",T},{"mass_kg",c.mass_kg.value()},{"energy_J",c.energy_J.json()},{"capacity_J_K",(c.mass_kg*D(v.apparent_heat_capacity_j_kg_k)).value()},{"k_W_m_K",v.k_w_m_k},{"sigma_S_m",v.sigma_s_m},{"phase",v.phase},{"hot_headroom_J",(c.mass_kg*D(law_.max_u())-c.energy_J).json()},{"cold_headroom_J",(c.energy_J-c.mass_kg*D(law_.min_u())).json()}});}return rows;}
 // Local half-Maxwell proposal, five caller-owned open uniforms. Native owner
 // must create actual IDs and certify its per-atom K before material preparation.
 J sample(const std::string&patchid,const std::array<double,5>&u)const{
  for(double v:u)need(finite(v)&&v>0&&v<1,"open source draws required");auto&p=patch(patchid);auto n=normal(p.face);bool vertical=p.face==Face::Left||p.face==Face::Right;std::array<double,3>x{};double t=p.lo+u[3]*(p.hi-p.lo);x[0]=vertical?(p.face==Face::Left?cfg_.x.front():cfg_.x.back()):t;x[1]=vertical?t:(p.face==Face::Bottom?cfg_.y.front():cfg_.y.back());x[2]=p.zlo+u[4]*(p.zhi-p.zlo);auto i=cell_on_patch(p,x);double T=cell_temperature(i);need(T>=1800,"local source below vapor scenario range");double vn=std::sqrt(-2*1.380649e-23*T/atom_mass_*std::log(u[0])),vt=std::sqrt(-2*1.380649e-23*T/atom_mass_*std::log(u[1])),a=2*std::acos(-1.)*u[2];std::array<double,3>v{n[0]*vn,n[1]*vn,vt*std::sin(a)};v[vertical?1:0]=vt*std::cos(a);double v2=0;for(double w:v)v2+=w*w;need(finite(v2)&&v2<299792458.*299792458.,"sample physical velocity outside native domain");double gamma=1/std::sqrt(1-v2/(299792458.*299792458.));std::array<double,3>proper;for(int d=0;d<3;++d)proper[d]=v[d]*gamma;return {{"patch",patchid},{"cell",i},{"sampled_temperature_K",T},{"position_m",x},{"proper_velocity_m_s",proper}};
 }
 // An integrated physical patch heat is partitioned by actual overlap area.
 // Last-cell exact remainder guarantees identical total J, even for non-dyadic
 // geometric fractions. Point particle transfers instead map to ONE FV cell.
 std::vector<Heat> patch_heat(const std::string&patchid,const D&Q,const std::string&id,const std::string&owner,std::uint64_t first_sequence)const{
  const auto&p=patch(patchid);bool vertical=p.face==Face::Left||p.face==Face::Right;const auto&e=vertical?cfg_.y:cfg_.x;std::vector<std::pair<std::size_t,double>>parts;for(std::size_t i=0;i+1<e.size();++i){double overlap=std::min(e[i+1],p.hi)-std::max(e[i],p.lo);if(overlap>0){std::size_t cell=vertical?i*nx()+(p.face==Face::Left?0:nx()-1):(p.face==Face::Bottom?0:(ny()-1)*nx())+i;parts.push_back({cell,overlap/(p.hi-p.lo)});}}
  need(!parts.empty()&&first_sequence>0&&first_sequence<=UINT64_MAX-parts.size(),"patch partition sequence/resource");std::vector<Heat>result;D used;for(std::size_t k=0;k<parts.size();++k){D q=k+1==parts.size()?Q-used:Q*D(parts[k].second);need(q.compare(zero())==0||q.compare(zero())==Q.compare(zero()),"negative patch partition remainder");used=used+q;result.push_back({id+"/"+std::to_string(k),owner,first_sequence+k,parts[k].first,std::move(q)});}need(eq(used,Q),"patch heat is not exactly conservative");return result;
 }
 Prepared prepare(const Step&s,const GasPort&port){return prepare_batch(s,port.certified_batch(s.binding));}
 // Root preflights native/field/source/RNG owners separately against the SAME
 // binding. This local preflight performs no callbacks or allocations.
 bool can_publish(const Prepared&p)const noexcept{return p.owner_==this&&p.status_==Prepared::Status::Staged&&p.generation_==generation_&&p.before_==accepted_&&generation_<UINT64_MAX;}
 bool publish_prevalidated(Prepared&p)noexcept{if(!can_publish(p))return false;accepted_=p.after_;p.status_=Prepared::Status::Published;++generation_;return true;}
 bool rollback(Prepared&p)noexcept{if(p.owner_!=this||p.status_!=Prepared::Status::Published||accepted_!=p.after_)return false;accepted_=p.before_;p.status_=Prepared::Status::RolledBack;if(generation_<UINT64_MAX)++generation_;return true;}
 J checkpoint(const J&outer_boundary)const{need(outer_boundary.at("accepted_boundary")==true&&outer_boundary.at("material_prefix")==prefix()&&outer_boundary.at("time_s")==accepted_->time_s&&outer_boundary.at("accepted_step")==accepted_->accepted_step,"partial/mismatched common checkpoint boundary");for(auto*k:{"gas_sha256","field_sha256","source_sha256","rng_sha256"})need(seal(outer_boundary.at(k)),"missing common checkpoint owner seal");J body={{"schema","spatial-material-v1"},{"signature",signature_},{"outer_boundary",outer_boundary},{"state",state_json(*accepted_)}};return {{"body",body},{"sha256",hash(body.dump())}};}
 // External history verifier must bind each archived certified batch to the
 // actual native journal. Rehashed material-only state cannot bypass replay.
 void restore(const J&record,const J&expected_outer_boundary,const GasPort&history){
  need(accepted_->revision==0&&generation_==0&&record.size()==2&&record.at("sha256")==hash(record.at("body").dump()),"restart digest/freshness");const auto&b=record.at("body");need(b.size()==4&&b.at("schema")=="spatial-material-v1"&&b.at("signature")==signature_&&b.at("outer_boundary")==expected_outer_boundary,"restart source/outer boundary identity");
  const auto&state=b.at("state");need(state.size()==7&&state.at("journal").is_array()&&state.at("journal").size()<=cfg_.maximum_records,"restart journal schema/resource");Owner replay(cfg_,initial_T_,initial_time_);
  for(const auto&row:state.at("journal")){Step s;s.binding=read_binding(row.at("binding"));for(const auto&h:row.at("heat"))s.heat.push_back(read_heat(h));auto p=replay.prepare(s,history);need(p.audit().at("row")==row,"restart native receipt/conduction replay mismatch");need(replay.publish_prevalidated(p),"restart local publication");}
  need(replay.state_json(*replay.accepted_)==state,"restart canonical material state mismatch");need(expected_outer_boundary.at("accepted_boundary")==true&&expected_outer_boundary.at("material_prefix")==replay.prefix()&&expected_outer_boundary.at("time_s")==replay.accepted_->time_s&&expected_outer_boundary.at("accepted_step")==replay.accepted_->accepted_step,"restart partial or mismatched common boundary");for(auto*k:{"gas_sha256","field_sha256","source_sha256","rng_sha256"})need(seal(expected_outer_boundary.at(k)),"restart missing owner identity");accepted_=replay.accepted_;++generation_;
 }
 static J describe_binding(const Binding&b){return binding_json(b);}
};
} // namespace spatial_metal
