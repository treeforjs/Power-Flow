// native_device.hpp -- internal interface of the native v7.1 device lane
// (implementation native_device.cpp; public API mfem_pic_full_time_device.h).
#pragma once
#include "native_rows.hpp"
#include "chem_v3_native.hpp"
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
// Device admission of a motion record (v7.1: without births and sinks; v8:
// with first-hit sinks and newborns, whose O(sinks + births) checks run on the
// host from the gathered rows below).
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
  // v8: sinks and newborns in the canonical layouts (record.device.sinks /
  // births when given there, else these host copies are uploaded).
  std::span<const pf_v7::Sink> sinks;
  std::span<const pf_v7::Birth> births;
  // v8.1: extra current (record.device.extra_J when given there, else this
  // host copy is uploaded; empty and no device address: none).
  std::span<const double> extra_J;
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
  std::shared_ptr<const DeviceRows> after; // endpoints copied (base columns; compacted with sinks)
  std::string after_carrier_root, paths_v5_root;
  // v8: per sink (ascending ID) the census row, path and packet key; per
  // newborn (ascending ID) the census row and its endpoint (zero row when sunk).
  std::vector<GapParticle> sunk_rows, born_rows, born_ends;
  std::vector<FullTimeMotionPath> sunk_paths;
  std::vector<std::uint64_t> sunk_keys, born_end_index;
  // v8.1: extra current given, its digest (block or streaming form) and the
  // exact extra work (value = dt x RN(sum E J_extra)).
  bool extra_given = false, extra_digest_ok = false;
  pf_v7::SumValue extra_work;
};
DeviceMotionOutput device_admit_motion(FullTimeDevice &, const DeviceMotionInput &);
FullTimeDeviceCounters &device_counters(FullTimeDevice &);
// Reads bytes at a producer device address (copy into a pooled block, then a
// download).
void device_read(FullTimeDevice &, std::uint64_t address, void *out, std::size_t bytes);
// ---- v8 chemistry v3 device stage (API-008.txt 3.7)
struct ChemDeviceRequest {
  const DeviceRows *rows = nullptr;       // frame rows with census columns
  const chem8::Setup *setup = nullptr;
  std::vector<std::uint32_t> state_species; // census state index -> registry species (none)
  pfc::StepConfig cfg{};
  std::uint32_t partitions = 0; // launches over the carriers (0: memory budget)
};
struct ChemDeviceStage {
  int error = 0;                 // 1 grid, 2 projectile packet, 3 ID domain, 4 cell targets
  std::uint64_t error_index = 0;
  chem8::StageActive active;
  std::string inputs_sha256, records_root;
  std::uint64_t projectiles = 0, retries = 0;
  double device_wall_ms = 0;
  std::shared_ptr<void> keep;    // device scratch kept for the debit callbacks
};
ChemDeviceStage device_chem_stage(FullTimeDevice &, const ChemDeviceRequest &);
// Exact stage-start count of (compact cell, species) (the inventory).
chem8::ESum device_chem_inventory(FullTimeDevice &, const ChemDeviceStage &, std::uint32_t cell, std::uint32_t species);
// Rows of (compact cell, species): (row, native ID, key, flags, weight) ascending native ID.
struct ChemSegmentRow {
  std::uint32_t row = 0, flags = 0;
  std::uint64_t id = 0, key = 0;
  double weight = 0;
  std::array<double, 3> position{}, proper_velocity{};
};
std::vector<ChemSegmentRow> device_chem_segment(FullTimeDevice &, const ChemDeviceStage &, std::uint32_t cell,
                                                std::uint32_t species);
struct ChemDeviceEdits {
  std::vector<std::pair<std::uint32_t, std::array<double, 3>>> velocity;
  std::vector<std::pair<std::uint32_t, double>> weight;
  std::vector<std::uint32_t> drop;
  std::vector<GapParticle> append; // ascending new native ID
  std::vector<std::uint64_t> append_keys, append_states;
};
// New device rows (with census columns and the v5 carrier root) after the edits.
std::shared_ptr<const DeviceRows> device_chem_apply(FullTimeDevice &, const DeviceRows &, const ChemDeviceEdits &);
} // namespace raikou::mfem_backend
