#pragma once
#include "native_parent_subflow.hpp"
#include "root_api.hpp"
#include <bit>
#include <map>
#include <span>

// Species and frame identities derived from live owner objects. A digest here
// is a recomputed comparison value, never an authority token by itself.
// Coordinator v5: unchanged from v4 except that the Hall identities come from
// Hall corrective v2 (cached, no per-call JSON). Record tags whose encoding did
// not change keep their v4 text, so their digests change only through inputs.
namespace raikou::native_water::v5 {
namespace h = ::hall_driver;
namespace hc = ::hall_corrective_v2;

// Length-prefixed little-endian binary encoding (same convention as the Hall
// motion record and the native census encoding).
struct Encoding {
  std::string bytes;
  void u64(std::uint64_t v){for(int i=0;i<8;++i)bytes.push_back(char((v>>(8*i))&255));}
  void i64(std::int64_t v){u64(static_cast<std::uint64_t>(v));}
  void f64(double v){u64(std::bit_cast<std::uint64_t>(v));}
  void text(const std::string& s){u64(s.size());bytes.append(s);}
  std::string finish()const{return sha(bytes);}
};

inline std::map<std::string,const m::ArenaStateSpec*> state_index(const std::vector<m::ArenaStateSpec>& states){
  std::map<std::string,const m::ArenaStateSpec*> out;
  for(const auto& s:states)
    pc::require(!s.id.empty()&&out.emplace(s.id,&s).second,"species catalog requires unique nonempty identities");
  pc::require(!out.empty()&&out.size()<=256,"bounded nonempty species catalog");
  return out;
}

inline bool same_channel(const m::ArenaChannel& a,const m::ArenaChannel& b){
  return a.id==b.id&&a.primary_state==b.primary_state&&a.target_state==b.target_state&&
    a.product_states==b.product_states&&a.dataset_sha256==b.dataset_sha256&&
    a.kernel_manifest_sha256==b.kernel_manifest_sha256&&
    a.kernel_binding_file_sha256==b.kernel_binding_file_sha256&&a.model==b.model;
}

// Canonical identity of the complete species/channel catalog, including the
// exact binary64 charge every carrier of each species must carry.
inline std::string species_sha256(const std::vector<m::ArenaStateSpec>& states,
                                  const std::vector<m::ArenaChannel>& channels){
  const auto index=state_index(states);
  std::map<std::string,const m::ArenaChannel*> catalog;
  for(const auto& c:channels){
    pc::require(!c.id.empty()&&catalog.emplace(c.id,&c).second,"channel catalog requires unique nonempty identities");
    pc::require(index.count(c.primary_state)&&index.count(c.target_state),"channel reactant species undeclared");
    for(const auto& p:c.product_states)pc::require(index.count(p)>0,"channel product species undeclared");
  }
  Encoding e;e.text("closed-parent-species-catalog-v4");e.u64(index.size());
  for(const auto& [id,s]:index){
    e.text(id);e.i64(s->hydrogen);e.i64(s->oxygen);e.i64(s->charge_e);
    e.f64(s->mass_kg);e.f64(s->internal_J);e.f64(s->reference_J);e.f64(s->charge_e*m::pic_e);
  }
  e.u64(catalog.size());
  for(const auto& [id,c]:catalog){
    e.text(id);e.text(c->primary_state);e.text(c->target_state);e.u64(c->product_states.size());
    for(const auto& p:c->product_states)e.text(p);
    e.text(c->dataset_sha256);e.text(c->kernel_manifest_sha256);
    e.text(c->kernel_binding_file_sha256);e.text(c->model);
  }
  return e.finish();
}

// Every actual carrier must carry exactly the charge and mass of its species.
inline void require_census_species(std::span<const m::ArenaRateEntry> census,
                                   const std::vector<m::ArenaStateSpec>& states){
  const auto index=state_index(states);
  for(const auto& row:census){
    const auto it=index.find(row.state);
    pc::require(it!=index.end(),"census carrier has an undeclared species");
    pc::require(row.native.charge_C==it->second->charge_e*m::pic_e&&row.native.mass_kg==it->second->mass_kg,
      "census carrier charge/mass differs from its declared species");
  }
}

// The declared catalog is compared with the registered native channels, every
// accepted packet species and every accepted native carrier.
inline void require_arena_species(const m::PicArena& arena,const std::vector<m::ArenaStateSpec>& states,
                                  const std::vector<m::ArenaChannel>& channels,std::size_t cap){
  const auto index=state_index(states);
  for(const auto& c:channels)
    pc::require(same_channel(arena.channel_binding(c.id),c),"registered native channel differs from the declared catalog");
  for(const auto& [id,p]:arena.snapshot()->packets){
    (void)id;pc::require(index.count(p.state)>0,"accepted native packet has an undeclared species");
  }
  require_census_species(arena.rate_population(cap),states);
}

// Exactly the translation native admission compares (native -> Hall field).
inline m::PicVector native_to_field_shift(const m::GapPicOptions& o){
  return {-.5*o.width_m,-.5*o.gap_m,0.};
}

inline std::string frame_sha256(const m::GapPicOptions& o,const hc::Model& model){
  Encoding e;e.text("native-to-hall-translation-frame-v4");
  for(double x:native_to_field_shift(o))e.f64(x);
  e.text(model.motion_mesh_sha256(o));
  e.text(h::kinetic_policy_sha256(h::declared_kinetic_policy(o)));
  return e.finish();
}
} // namespace raikou::native_water::v5
