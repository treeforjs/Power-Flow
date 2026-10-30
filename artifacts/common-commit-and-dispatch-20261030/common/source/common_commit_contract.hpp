#pragma once
#include "raikou/mfem_pic_full_time.h"
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>

namespace raikou::native_water::common_commit {
using Boundary = ::raikou::mfem_backend::FullTimeAcceptedBoundary;

// Observations are derived by the coordinator from real owner objects. This
// predicate checks cross-owner agreement; these public value types are NOT a
// factory for a native, Hall, source or RNG authority token.
struct ObservedOwners {
  Boundary native;
  std::string actual_native_prefix, field, heat, source, rng;
  std::string mesh, material, frame, configuration;
  std::uint64_t actual_native_revision=0, actual_native_step=0;
  std::uint64_t field_revision=0, heat_revision=0;
  double field_time_s=0, heat_time_s=0;
};

struct ExpectedStep {
  Boundary before;
  double dt_s=0, end_s=0;
  std::string common_clock;
  std::string configuration;
};

struct ProposedOwners {
  ObservedOwners after;
  std::string motion_midpoint, field_midpoint;
};

inline bool digest(const std::string& s) noexcept {
  if(s.size()!=64)return false;
  for(unsigned char c:s)
    if(!((c>='0'&&c<='9')||(c>='a'&&c<='f')))return false;
  return true;
}

inline bool valid_boundary(const Boundary& b) noexcept {
  return std::isfinite(b.time_s)&&b.time_s>=0&&
    digest(b.accepted_arena_sha256)&&digest(b.field_sha256)&&
    digest(b.heat_sha256)&&digest(b.source_sha256)&&digest(b.rng_sha256)&&
    digest(b.mesh_sha256)&&digest(b.material_sha256)&&
    digest(b.frame_sha256)&&digest(b.common_clock_sha256);
}

inline bool same_boundary(const Boundary& a,const Boundary& b) noexcept {
  return a.accepted_step==b.accepted_step&&
    a.accepted_arena_revision==b.accepted_arena_revision&&a.time_s==b.time_s&&
    a.accepted_arena_sha256==b.accepted_arena_sha256&&
    a.field_sha256==b.field_sha256&&a.heat_sha256==b.heat_sha256&&
    a.source_sha256==b.source_sha256&&a.rng_sha256==b.rng_sha256&&
    a.mesh_sha256==b.mesh_sha256&&a.material_sha256==b.material_sha256&&
    a.frame_sha256==b.frame_sha256&&a.common_clock_sha256==b.common_clock_sha256&&
    a.authored_full_time_genesis==b.authored_full_time_genesis;
}

inline bool owner_values_match(const ObservedOwners& o) noexcept {
  const auto& b=o.native;
  return valid_boundary(b)&&o.actual_native_prefix==b.accepted_arena_sha256&&
    o.actual_native_revision==b.accepted_arena_revision&&
    o.actual_native_step==b.accepted_step&&
    o.field==b.field_sha256&&o.heat==b.heat_sha256&&
    o.source==b.source_sha256&&o.rng==b.rng_sha256&&
    o.mesh==b.mesh_sha256&&o.material==b.material_sha256&&
    o.frame==b.frame_sha256&&digest(o.configuration)&&
    o.field_revision==b.accepted_step&&o.heat_revision==b.accepted_step&&
    o.field_time_s==b.time_s&&o.heat_time_s==b.time_s;
}

// Closed resident scope only: source and geometry do not change. A future
// source/geometry-changing common step needs a different explicit contract.
// Call this after obtaining fresh observations and before ANY accepted swap.
inline bool accepts_closed_step(const ExpectedStep& e,
    const ObservedOwners& current,const ProposedOwners& candidate) noexcept {
  const auto& b=e.before;const auto& a=candidate.after.native;
  if(!valid_boundary(b)||!owner_values_match(current)||
     !same_boundary(b,current.native)||!owner_values_match(candidate.after))return false;
  if(!std::isfinite(e.dt_s)||e.dt_s<=0||e.dt_s>1e-14||.5*e.dt_s<=0||!std::isfinite(e.end_s)||
     e.end_s!=b.time_s+e.dt_s||e.end_s<=b.time_s||
     !digest(e.common_clock)||!digest(e.configuration)||
     current.configuration!=e.configuration||
     candidate.after.configuration!=e.configuration)return false;
  constexpr auto max=std::numeric_limits<std::uint64_t>::max();
  if(b.accepted_step>=max-1||b.accepted_arena_revision>=max-1||
     a.accepted_step!=b.accepted_step+1||
     a.accepted_arena_revision<=b.accepted_arena_revision||
     a.accepted_arena_revision==max||a.time_s!=e.end_s||
     a.authored_full_time_genesis)return false;
  if(a.common_clock_sha256!=e.common_clock||
     a.mesh_sha256!=b.mesh_sha256||a.material_sha256!=b.material_sha256||
     a.frame_sha256!=b.frame_sha256||a.source_sha256!=b.source_sha256)return false;
  return digest(candidate.motion_midpoint)&&
    candidate.motion_midpoint==candidate.field_midpoint;
}
} // namespace raikou::native_water::common_commit
