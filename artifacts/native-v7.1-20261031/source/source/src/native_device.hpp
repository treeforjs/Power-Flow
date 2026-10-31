// native_device.hpp -- internal interface of the native v7.1 device lane
// (implementation native_device.cpp; public API mfem_pic_full_time_device.h).
#pragma once
#include "native_rows.hpp"
#include "pf_v7/forms.hpp"
#include "raikou/mfem_pic_full_time_device.h"
#include <array>
#include <memory>
#include <string>

namespace raikou::mfem_backend {
struct DevicePool;
// One pooled device allocation (returned to the pool on destruction).
struct DeviceBlock {
  std::shared_ptr<DevicePool> pool;
  std::uint64_t address = 0;
  std::size_t bytes = 0;
  int slot = -1;
  DeviceBlock() = default;
  DeviceBlock(const DeviceBlock &) = delete;
  DeviceBlock &operator=(const DeviceBlock &) = delete;
  ~DeviceBlock();
};
// Census columns of device rows: packet key and state index per row.
struct DeviceColumns {
  std::shared_ptr<DeviceBlock> keys, states;
  std::size_t n = 0;
};
struct DeviceRows {
  std::shared_ptr<FullTimeDevice::Impl> dev;
  std::shared_ptr<DeviceBlock> rows; // n x 80 B (GapParticle)
  std::shared_ptr<const DeviceColumns> columns; // may be null
  std::size_t n = 0;
  std::string carrier_root; // v5 carrier root (hex) when known
};
// Rows (+ optional census columns) uploaded from the host.
std::shared_ptr<const DeviceRows> device_upload_rows(FullTimeDevice &, const std::vector<GapParticle> &,
                                                     const std::vector<std::uint64_t> *keys,
                                                     const std::vector<std::uint64_t> *states);
// Census v2 of device rows with columns (packet keys, states).
FullTimeCensusDevice device_census(FullTimeDevice &, const DeviceRows &, const std::string &states_sha256);
// Device admission of a motion record without births and sinks.
struct DeviceMotionInput {
  const DeviceRows *base = nullptr;
  const FullTimeMotionRecordV7 *record = nullptr;
  double width = 0, gap = 0, short_z = 0, feed_z = 0;
  std::array<double, 4> lower{}, upper{};
  std::array<double, 3> shift{};
  double dt = 0, neutral_bound = 0;
  bool finite_foils = false, finite_axial = false;
  std::string states_sha256;
  // Census v2 root of the base rows when already computed for the trial
  // (census_v7_device / census_sha256); empty: computed here.
  std::string census_root;
};
struct DeviceMotionOutput {
  bool carriers_pass = false;  // every per-carrier check passed
  std::uint64_t first_fail = ~std::uint64_t(0);
  std::uint32_t fail_code = 0;
  pf_v7::SumValue kinetic, work;   // K; W (value = RN(dt x sum))
  pf_v7::WorkScales scales;        // motion-work-gross-v1 scales and tolerance
  std::string census_root, census_sha256;
  pf_v7::MotionSectionsV7 sections;
  std::string E_root, B_root, J_root;
  bool mesh_ok = false;
  // record.trial_field_sha256 / current_A_m_sha256 equal the block form or,
  // checked only when it does not, the v7 streaming form (hashed on the host).
  bool field_digest_ok = false, current_digest_ok = false;
  std::shared_ptr<const DeviceRows> after; // endpoints copied (base columns)
  std::string after_carrier_root, paths_v5_root;
};
DeviceMotionOutput device_admit_motion(FullTimeDevice &, const DeviceMotionInput &);
FullTimeDeviceCounters &device_counters(FullTimeDevice &);
// Reads bytes at a producer device address (copy into a pooled block, then a
// download).
void device_read(FullTimeDevice &, std::uint64_t address, void *out, std::size_t bytes);
} // namespace raikou::mfem_backend
