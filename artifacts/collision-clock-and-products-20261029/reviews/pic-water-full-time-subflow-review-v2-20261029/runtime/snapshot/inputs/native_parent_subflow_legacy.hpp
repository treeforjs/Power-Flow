#pragma once
#include "raikou/mfem_pic_arena.h"
#include "water_parent_bridge.hpp"
#include <openssl/sha.h>
#include <iomanip>

// Finite-rate, parent-ionization-only integration. This is an operator subflow
// at a frozen particle stage, not a replacement for the common physical clock.
namespace raikou::native_water {
namespace m = mfem_backend;
namespace pc = pic_collisions;
namespace w = pic_water;
using Json = m::ArenaJson;
inline std::string sha(const std::string& bytes) {
  unsigned char out[SHA256_DIGEST_LENGTH];
  SHA256(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size(),out);
  std::ostringstream s;
  for(unsigned char c:out)s<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(c);
  return s.str();
}
struct Limits {
  std::size_t reactions=128, carriers=100000, cells=1000000,
              pairs=1000000, queries=1000000, candidates=100000;
};
struct Binding {
  std::string dataset_sha256,kernel_sha256,kernel_binding_file_sha256;
};
struct Interval {
  std::uint64_t sequence=0;
  double begin_s=0,end_s=0;
  m::ArenaPhase phase;
  std::string common_clock_sha256;
};
struct Grid {
  std::array<double,3> low{},high{};
  std::array<std::uint32_t,3> counts{};
  std::vector<pc::Cell> cells;
  std::vector<unsigned char> vacuum;
  std::vector<m::GapFoil> foils;
  std::string identity;
  Grid(const m::GapPicOptions& o,std::array<std::uint32_t,3> n,std::size_t cap)
      :low{0,0,o.short_z_m},high{o.width_m,o.gap_m,o.feed_z_m},counts(n) {
    pc::require(o.finite_axial_extent,"collision grid requires declared physical axial bounds");
    std::size_t total=1;
    long double volume=1;
    for(int d=0;d<3;++d){
      pc::require(n[d]>0&&total<=cap/n[d],"collision grid cell cap");total*=n[d];
      pc::require(std::isfinite(low[d])&&std::isfinite(high[d])&&high[d]>low[d],"collision grid extent");
      volume*=(static_cast<long double>(high[d])-low[d])/n[d];
    }
    pc::require(o.depth_m==o.feed_z_m-o.short_z_m,"collision depth differs from physical height");
    pc::positive(static_cast<double>(volume),"unresolved collision volume");
    std::vector<std::array<std::uint64_t,4>> rectangles;
    auto aligned=[&](double x,int d){
      const long double q=(static_cast<long double>(x)-low[d])*n[d]/(static_cast<long double>(high[d])-low[d]);
      const long double nearest=std::round(q);
      pc::require(nearest>=0&&nearest<=n[d]&&std::abs(q-nearest)<=
        64*std::numeric_limits<double>::epsilon()*std::max(1.L,std::abs(q)),
        "collision grid cuts a foil; exact vacuum clipping is not implemented");
      return static_cast<std::uint64_t>(nearest);
    };
    if(o.finite_foils)for(const auto& f:{o.lower_foil,o.upper_foil}){
      rectangles.push_back({aligned(f.left_m,0),aligned(f.right_m,0),aligned(f.bottom_m,1),aligned(f.top_m,1)});
      foils.push_back(f);
    }
    cells.reserve(total);vacuum.resize(total,1);
    for(std::size_t k=0;k<n[2];++k)for(std::size_t j=0;j<n[1];++j)for(std::size_t i=0;i<n[0];++i){
      const auto id=i+n[0]*(j+n[1]*k);
      for(const auto& r:rectangles)if(i>=r[0]&&i<r[1]&&j>=r[2]&&j<r[3])vacuum[id]=0;
      if(vacuum[id])cells.push_back({id,static_cast<double>(volume)});
    }
    identity=sha(Json{{"schema","physical-cartesian-collision-grid-v1"},{"low",low},
      {"high",high},{"counts",counts},{"solid_xy_index_ranges",rectangles},
      {"volume_m3",static_cast<double>(volume)}}.dump());
  }
  std::uint64_t locate(const m::PicVector& x)const {
    for(const auto& f:foils){
      const bool vertical=(x[0]==f.left_m||x[0]==f.right_m)&&x[1]>=f.bottom_m&&x[1]<=f.top_m;
      const bool horizontal=(x[1]==f.bottom_m||x[1]==f.top_m)&&x[0]>=f.left_m&&x[0]<=f.right_m;
      pc::require(!vertical&&!horizontal,"collision on material surface requires side-resolved transport stage");
    }
    std::array<std::uint64_t,3> i{};
    for(int d=0;d<3;++d){
      pc::require(std::isfinite(x[d])&&x[d]>=low[d]&&x[d]<=high[d],"carrier outside collision domain");
      long double f=(static_cast<long double>(x[d])-low[d])/(static_cast<long double>(high[d])-low[d]);
      i[d]=std::min<std::uint64_t>(counts[d]-1,static_cast<std::uint64_t>(f*counts[d]));
    }
    const auto id=i[0]+counts[0]*(i[1]+counts[1]*i[2]);
    pc::require(vacuum.at(id),"collision carrier is in a solid cell");return id;
  }
};
struct State {
  std::string scheduler,product_rng;
  std::uint64_t next_event=1,next_interval=0;
};
class ParentSubflow;
class Prepared {
  friend class ParentSubflow;
  std::shared_ptr<const char> lifetime;
  std::uint64_t generation=0;
  std::unique_ptr<m::PreparedArena> arena;
  std::unique_ptr<State> after;
  Json audit_;
public:
  const Json& audit()const{return audit_;}
};
class ParentSubflow {
  m::PicArena& arena_;
  const Grid grid_;
  const Binding binding_;
  const Limits limits_;
  const w::ParentBridge bridge_;
  std::unique_ptr<State> accepted_;
  const std::shared_ptr<const char> lifetime_=std::make_shared<const char>(0);
  std::uint64_t generation_=0;
  static constexpr const char* owner_id="native_water_parent_subflow/v1";
  static Json payload(const State& s){return {{"scheduler",s.scheduler},{"product_rng",s.product_rng},
    {"next_event",s.next_event},{"next_interval",s.next_interval}};}
  void check_completion(const State& s)const {
    auto prior=arena_.stochastic_completion(owner_id);
    if(prior){
      pc::require(prior->stage_sequence!=UINT64_MAX&&prior->stage_sequence+1==s.next_interval&&
        prior->after_rng_clock_sha256==sha(payload(s).dump()),"RNG/clock state does not match accepted native exposure");
    }else pc::require(s.next_interval==0,"missing accepted native exposure history");
  }
  static std::string encode(const std::mt19937_64& rng){std::ostringstream s;s<<rng;return s.str();}
  static std::mt19937_64 decode(const std::string& text){
    pc::require(text.size()<10000,"product RNG checkpoint cap");
    std::istringstream s(text);std::mt19937_64 r;
    pc::require(bool(s>>r),"product RNG checkpoint");s>>std::ws;
    pc::require(s.eof(),"product RNG trailing bytes");return r;
  }
  static double uniform(std::mt19937_64& rng){return (static_cast<double>(rng()>>12)+.5)*0x1p-52;}
  static pc::Vec direction(std::mt19937_64& rng){
    const double z=2*uniform(rng)-1,phi=6.283185307179586476925286766559*uniform(rng);
    const double r=std::sqrt((1-z)*(1+z));return {r*std::cos(phi),r*std::sin(phi),z};
  }
  struct Population {
    std::vector<pc::Packet> electrons,water;
    std::map<std::string,std::string> parent;
    std::map<std::string,pc::Packet> by_id;
  };
  Population population(const m::PicArena& a)const {
    Population p;
    // One bounded census validates all species and exact counts once. Calling
    // rate_carriers separately for every packet repeated whole-population work.
    for(const auto& entry:a.rate_population(limits_.carriers)){
      if(entry.state!="electron"&&entry.state!="H2O")continue;
      const auto& carrier=entry.native;
        const auto name="native-carrier/"+std::to_string(carrier.id);
        pc::Packet q{name,entry.state,grid_.locate(carrier.position_m),carrier.weight,
                     carrier.mass_kg,carrier.proper_velocity_m_s};
        pc::require(p.parent.emplace(name,entry.packet_id).second,"duplicate native rate carrier");
        p.by_id.emplace(name,q);
        (entry.state=="electron"?p.electrons:p.water).push_back(q);
    }
    return p;
  }
public:
  ParentSubflow(m::PicArena& arena,std::array<std::uint32_t,3> grid_counts,Binding binding,Limits limits,
                const w::ParentBridge& bridge,std::uint64_t clock_seed,std::uint64_t product_seed)
      :arena_(arena),grid_(arena.native_options(),grid_counts,limits.cells),binding_(std::move(binding)),limits_(limits),bridge_(bridge),accepted_(std::make_unique<State>()) {
    for(const auto& s:{binding_.dataset_sha256,binding_.kernel_sha256,binding_.kernel_binding_file_sha256})
      pc::require(pc::digest(s),"unbound native subflow inputs");
    const auto c=bridge_.channel();
    pc::require(c.data_sha256==binding_.dataset_sha256&&c.kernel_sha256==binding_.kernel_sha256,"bridge/binding mismatch");
    const auto& catalog=arena_.channel_binding(w::parent_channel);
    pc::require(catalog.dataset_sha256==binding_.dataset_sha256&&catalog.kernel_manifest_sha256==binding_.kernel_sha256&&
      catalog.kernel_binding_file_sha256==binding_.kernel_binding_file_sha256&&catalog.model==w::parent_model,
      "native registered channel differs from subflow binding");
    pc::require(limits_.reactions>0&&limits_.carriers>0&&limits_.cells>=grid_.cells.size(),"native subflow limits");
    accepted_->scheduler=pc::Scheduler(clock_seed).checkpoint();
    accepted_->product_rng=encode(std::mt19937_64(product_seed));
  }
  ParentSubflow(const ParentSubflow&)=delete;
  ParentSubflow& operator=(const ParentSubflow&)=delete;
  Json checkpoint()const {
    check_completion(*accepted_);
    return {{"schema","native-water-parent-subflow-v2"},{"grid",grid_.identity},
      {"dataset",binding_.dataset_sha256},{"kernel",binding_.kernel_sha256},
      {"binding_file",binding_.kernel_binding_file_sha256},{"arena",arena_.state_prefix()},
      {"scheduler",accepted_->scheduler},{"product_rng",accepted_->product_rng},
      {"next_event",accepted_->next_event},{"next_interval",accepted_->next_interval}};
  }
  // Called only beside a matching restored native arena. Validation precedes swap.
  void restore(const Json& j){
    pc::require(generation_!=UINT64_MAX,"subflow revision exhausted");
    pc::require(j.at("schema")=="native-water-parent-subflow-v2"&&j.at("grid")==grid_.identity&&
      j.at("dataset")==binding_.dataset_sha256&&j.at("kernel")==binding_.kernel_sha256&&
      j.at("binding_file")==binding_.kernel_binding_file_sha256&&j.at("arena")==arena_.state_prefix(),"subflow restart binding mismatch");
    auto s=std::make_unique<State>();s->scheduler=j.at("scheduler");s->product_rng=j.at("product_rng");
    pc::require(s->scheduler.size()<10000,"scheduler restart size cap");
    (void)pc::Scheduler::restore(s->scheduler);(void)decode(s->product_rng);
    pc::require(j.at("next_event").is_number_unsigned()&&j.at("next_interval").is_number_unsigned(),"restart sequence type");
    s->next_event=j.at("next_event");s->next_interval=j.at("next_interval");
    pc::require(s->next_event>0&&s->next_event!=UINT64_MAX&&s->next_interval!=UINT64_MAX,"restart event/interval sequence exhausted");
    check_completion(*s);accepted_.swap(s);++generation_;
  }
  Prepared prepare(const Interval& interval){
    check_completion(*accepted_);
    pc::require(generation_!=UINT64_MAX&&interval.sequence==accepted_->next_interval&&interval.sequence!=UINT64_MAX,"duplicate/stale subflow interval");
    pc::require(pc::digest(interval.common_clock_sha256)&&interval.phase.mesh_prefix==grid_.identity&&
      interval.phase.clock_prefix==interval.common_clock_sha256&&interval.phase.state_prefix==arena_.state_prefix(),"native subflow common stage binding");
    auto clock=std::unique_ptr<pc::Scheduler>(new pc::Scheduler(pc::Scheduler::restore(accepted_->scheduler)));
    auto rng=decode(accepted_->product_rng);
    auto next=std::make_unique<State>(*accepted_);
    Json events=Json::array();std::size_t pairs=0,queries=0;
    double elapsed=0;bool terminal=false;
    const auto before_hash=sha(payload(*accepted_).dump());
    auto token=arena_.prepare_compound_with_completion(limits_.reactions,[&](m::PicArena& scoped){
      while(!terminal){
        const auto p=population(scoped);
        pc::CandidateInputLimits caps;caps.projectiles=caps.targets=limits_.carriers;caps.cells=limits_.cells;
        pc::CandidateWork work;
        auto candidates=pc::candidates(p.electrons,p.water,grid_.cells,{bridge_.channel()},
          limits_.pairs-pairs,limits_.queries-queries,limits_.candidates,caps,&work);
        pairs+=work.pair_visits;queries+=work.channel_queries;
        const auto prefix=scoped.state_prefix();
        pc::Clock stage{interval.phase.accepted_step,scoped.snapshot()->revision,interval.begin_s,
          interval.end_s,prefix,grid_.identity,binding_.kernel_sha256,interval.phase.physical_stage_s,elapsed};
        const auto proposal=clock->prepare(stage,candidates,limits_.candidates);
        if(!proposal.event()){
          pc::require(clock->try_publish(proposal),"terminal scheduler publication");terminal=true;break;
        }
        pc::require(events.size()<limits_.reactions&&next->next_event!=UINT64_MAX,"native reaction/event resource limit");
        const auto& c=*proposal.event();
        const auto products=bridge_.sample(c,p.by_id.at(c.projectile_id),p.by_id.at(c.target_id),direction(rng));
        const std::string id="conditional/native-parent/"+std::to_string(next->next_event++);
        m::ArenaEvent e;e.id=id;e.channel=c.channel_id;e.primary=p.parent.at(c.projectile_id);
        e.target=p.parent.at(c.target_id);e.event_weight=c.event_weight;e.scattered_primary_id=e.primary;
        const auto remaining=m::PicArena::debit(scoped.snapshot()->packets.at(e.primary).count,e.event_weight);
        if(remaining.hi>0||remaining.lo>0)e.untouched_primary_id=id+"/unscattered-primary";
        e.phase=interval.phase;e.phase.state_prefix=prefix;
        e.provenance={m::ArenaProvenance::Kind::ConditionalPhysicalData,binding_.dataset_sha256,
          binding_.kernel_sha256,binding_.kernel_binding_file_sha256,w::parent_model,
          "ssa_weighted_pair_clock_v1",prefix,proposal.event_subflow_elapsed_s()};
        for(int i=0;i<3;++i)e.products.push_back({i==0?e.primary:id+(i==1?"/secondary":"/ion"),
          products.particles[i].state,products.particles[i].proper_velocity_m_s,
          i==0?m::ArenaProduct::At::Primary:m::ArenaProduct::At::Target});
        auto reaction=scoped.prepare(e);
        pc::require(clock->can_publish(proposal)&&scoped.can_publish(reaction),"reaction common preflight");
        pc::require(scoped.try_publish(reaction),"native reaction publication");
        pc::require(clock->try_publish(proposal),"reaction RNG publication");
        elapsed=proposal.event_subflow_elapsed_s();
        events.push_back({{"event",id},{"cell",c.cell},{"physical_volume_m3",c.volume_m3},
          {"weight",c.event_weight},{"rate_s",c.rate_s},{"total_rate_s",proposal.total_rate_s()},
          {"energy_eV",c.energy_eV},{"physical_stage_s",stage.frozen_physical_stage_s},
          {"collision_elapsed_s",elapsed},{"native",reaction.audit()}});
      }
      next->scheduler=clock->checkpoint();next->product_rng=encode(rng);++next->next_interval;
      m::ArenaStochasticCompletion completion;
      completion.integration_owner_id=owner_id;
      completion.before_rng_clock_sha256=before_hash;
      completion.after_rng_clock_sha256=sha(payload(*next).dump());
      completion.channel_id=w::parent_channel;
      completion.dataset_sha256=binding_.dataset_sha256;
      completion.kernel_manifest_sha256=binding_.kernel_sha256;
      completion.kernel_binding_file_sha256=binding_.kernel_binding_file_sha256;
      completion.model=w::parent_model;
      completion.common_clock_sha256=interval.common_clock_sha256;
      completion.stage_sequence=interval.sequence;
      completion.frozen_physical_stage_s=interval.phase.physical_stage_s;
      completion.common_interval_begin_s=interval.begin_s;
      completion.common_interval_end_s=interval.end_s;
      return completion;
    });
    Prepared out;out.lifetime=lifetime_;out.generation=generation_;
    out.arena=std::make_unique<m::PreparedArena>(std::move(token));out.after=std::move(next);
    out.audit_={{"schema","native-parent-subflow-receipt-v1"},{"conditional_parent_only",true},
      {"full_mechanism",false},{"physical_time_advanced",false},{"interval_sequence",interval.sequence},
      {"interval_begin_s",interval.begin_s},{"interval_end_s",interval.end_s},
      {"events",std::move(events)},{"pairs",pairs},{"queries",queries},{"terminal_clock_sampled",terminal},
      {"compound",out.arena->audit()}};
    return out;
  }
  bool can_publish(const Prepared& p)const noexcept {
    return p.lifetime==lifetime_&&p.generation==generation_&&generation_!=UINT64_MAX&&p.after&&p.arena&&arena_.can_publish(*p.arena);
  }
  bool try_publish(Prepared& p)noexcept {
    if(!can_publish(p))return false;
    if(!arena_.try_publish(*p.arena))return false;
    accepted_.swap(p.after);++generation_;return true;
  }
  const std::string& grid_identity()const{return grid_.identity;}
};
} // namespace raikou::native_water
