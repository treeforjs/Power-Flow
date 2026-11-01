#pragma once
// Interface of the sealed heat/source owner v5 harness (source_v5_harness.cpp,
// compiled in its own translation unit with the source v5 include paths). Only
// standard types and the arrival owner v1 types cross this boundary, so the
// native v6 test translation unit never sees the source v5 headers.
#include "pf_arrival/arrival_owner.hpp"
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace source_harness {
using V3d = std::array<double, 3>;
// Source-frame (field-frame) kinematics of one birth.
struct Birth {
  double weight = 0, time = 0;
  V3d x{}, u{}, outward_normal{};
  std::string event;
  bool active_layer = false;
};
struct Step {
  bool admitted = false;
  std::string outcome, refusal;
  pf_arrival::Count debit, inventory_before, inventory_after;
  std::vector<Birth> births;
  std::string source_before, source_after;
};
// Neutral desorption patch (lower face of the upper foil), quote -> leaves.
class Neutral {
public:
  // start_s: the common clock at genesis (equal to the native genesis time).
  Neutral(double N0, double rate, const std::string &sampling_identity, double start_s);
  ~Neutral();
  std::uint64_t revision() const;
  std::string state_sha256() const;
  // Release quote of the next step for an arrival neutral stream.
  pf_arrival::NeutralRelease quote(const pf_arrival::Policy &, std::uint32_t stream,
                                   double dt);
  // Owner::prepare with the staged leaves; binds the source release times.
  Step prepare(pf_arrival::Staged &, std::uint32_t stream);
  bool commit();

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
// Electron emission patch (same surface), exposure -> births via begin_source.
class Electron {
public:
  Electron(const std::string &sampling_identity, double field_V_m, double start_s);
  ~Electron();
  std::uint64_t revision() const;
  std::string state_sha256() const;
  // Returns false (with the refusal text) when source v5 refuses the exposure.
  bool expose(std::uint32_t stream, double area_fraction, double dt,
              pf_arrival::ElectronExposure &out, std::string &refusal);
  Step births(const pf_arrival::Staged &);

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
double probe_neutral_debit(double N0, double rate, double dt);
double probe_electron_count(double field_V_m, double area_fraction, double dt);
// Source frame geometry (eh::Geometry defaults): the patch is the face y = y0.
inline constexpr double patch_face_y_m = .0025;
} // namespace source_harness
