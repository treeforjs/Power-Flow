#pragma once
// Coordinator v8: accepted rate samples of the source v7 slice owners and the source-owner interval
// bound (source v7 API-v7.txt 3.2, INTEGRATION.txt D5).
//   After genesis (or a restore) and after every publication, each slice owner evaluates
//   Owner::rate_sample(epoch, field, origins) on the accepted state: the epoch of the accepted thermal
//   state with length rate_sample_dt_s (the controller cap; it sets the per-step count of the policy's
//   count floor), the source Field of the accepted Hall E and B with the configured trial bound, one
//   origin per patch (the face centre, as the emission exposure). The sample of the state before is kept.
//   The bound of a step is min over slice owners of Owner::interval_bound({previous, current}) (the owner's
//   own function, recomputed by the owner in prepare from the same evidence); the controller applies it as
//   the criterion "source_interval". While only one accepted sample exists (start-up) the bound is the
//   start-up step (the particle cap) and the evidence's previous sample is the current one re-timed
//   (joint_step_v8.hpp prepare_v8).
// Only runs whose source owners have an interval policy (nonzero error targets) sample rates.
#include "joint_step_v8.hpp"

namespace raikou::native_water::v8 {
struct SourceBound {
  double bound_s = INFINITY;
  bool startup = false;
  std::string limiting = "none";
  double rate_per_s = 0;
  std::size_t slice = 0, index = 0;
};

// Source trial frame of the accepted state with step length dt (the frame of prepare_v8 for that dt).
template <class H> v3k::TrialFrame accepted_source_frame(const OwnersV8<H>& o, double dt_s, const std::string& start_field_digest) {
  v3k::TrialFrame f;
  f.role = v3k::FieldRole::TotalMidpoint;
  f.accepted_revision = o.hall->revision();
  f.trial_stage = 1;
  f.mesh_revision = 1;
  f.material_revision = 2;
  f.begin_s = o.hall->time_s();
  f.end_s = f.begin_s + dt_s;
  f.dt_s = dt_s;
  f.field_stage_s = f.begin_s + .5 * dt_s;
  f.arena_hash = o.arena->state_prefix();
  f.field_hash = "hall-v3-accepted-start-field/" + start_field_digest;
  f.mesh_hash = immutable_digests(o).mesh;
  f.material_hash = o.hall->model().material_sha256();
  f.species_hash = o.inputs.species_sha256;
  return f;
}

// Rate samples of the accepted state (shifts current -> previous). No-op without an interval policy.
template <class H> void sample_source_rates(OwnersV8<H>& o) {
  if (!o.interval_policy()) return;
  const auto acc = o.accepted();
  const auto digest = h::total_field_sha256(*acc.E, *acc.B);
  const auto& mesh = o.hall->mesh();
  o.mesh_identity = s7::MeshIdentity::of(mesh);
  auto field = std::make_shared<const s7::Field>(o.mesh_identity, accepted_source_frame(o, o.rate_sample_dt_s, digest), *acc.E, *acc.B,
                                                 o.source_options.interval.cap_s);
  std::vector<sv7::RateSample> next;
  for (std::size_t i = 0; i < o.sources.size(); ++i) {
    std::vector<s7::V3> origins(o.slices[i].patches.size());
    for (const auto& r : o.patches)
      if (r.slice == i) origins.at(r.local) = patch_origin(mesh, r);
    const auto e = o.sources[i]->make_epoch(o.rate_sample_dt_s, "coordinator-v8/rate-sample/" + digest);
    next.push_back(o.sources[i]->rate_sample(e, field, origins));
  }
  if (o.rate_current.size() == next.size()) o.rate_previous = std::move(o.rate_current);
  else o.rate_previous.clear();
  o.rate_current = std::move(next);
}

// Bound of the next step from the accepted samples (the owners' interval_bound).
template <class H> SourceBound source_bound(const OwnersV8<H>& o) {
  SourceBound b;
  if (!o.interval_policy()) return b;
  need(o.rate_current.size() == o.sources.size(), "source interval policy without accepted rate samples");
  if (o.rate_previous.size() != o.sources.size()) {
    b.bound_s = o.interval_startup_dt_s;
    b.startup = true;
    b.limiting = "start-up";
    return b;
  }
  for (std::size_t i = 0; i < o.sources.size(); ++i) {
    const auto d = o.sources[i]->interval_bound({o.rate_previous[i], o.rate_current[i]});
    if (d.bound_s < b.bound_s) {
      b.bound_s = d.bound_s;
      b.limiting = d.limiting;
      b.rate_per_s = d.rate_per_s;
      b.slice = i;
      b.index = d.index;
    }
  }
  return b;
}
}  // namespace raikou::native_water::v8
