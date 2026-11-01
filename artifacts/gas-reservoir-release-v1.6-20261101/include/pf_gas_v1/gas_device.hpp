// gas_device.hpp -- gas target owner v1 device lane: GPU-resident candidate
// and accepted state (entry digits, rounded entry values, global digits) on
// one hall_gpu::Session stream; kernels in device/gas_kernels.cu compiled by
// NVRTC from the gas_core.h and forms_core.h bytes of the host build.
#pragma once
#include "pf_gas_v1/gas_owner.hpp"

#include <array>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace pf_gas_v1 {

using Leaf = std::array<pf7_u32, 8>;

class DeviceLane {
public:
  DeviceLane(const Configuration &, const std::vector<pfg::SpeciesInfo> &, const std::vector<pfg::ProcessInfo> &,
             pfg::Layout);
  ~DeviceLane();
  DeviceLane(const DeviceLane &) = delete;
  DeviceLane &operator=(const DeviceLane &) = delete;

  void begin_step(std::uint64_t step);   // candidate = accepted
  // Make the device key table a superset of keys (both candidate and accepted
  // arrays are re-keyed; zero entries are not part of the state).
  void ensure_keys(const std::vector<std::uint64_t> &keys);
  // Applies the inputs into a zeroed delta, refuses (Refusal::input, nothing
  // merged) on status bits, else merges the delta into the candidate and
  // normalises it. Returns the delta counters (records/flows/movers are counted
  // by the caller).
  std::array<std::uint64_t, pfg::C_COUNTERS> apply(std::span<const RecoilRecord>, std::span<const FlowRecord>,
                                                   std::span<const MoverRecord>);
  FeedbackReceipt feedback(std::span<const std::uint32_t> atlas, std::vector<CellTargetsAbi> &cells, double volume,
                           bool candidate);
  void commit() noexcept;
  void abandon() noexcept;
  void leaves(bool candidate, std::vector<Leaf> &entry, std::vector<Leaf> &global);
  // Normalised state download (nonzero entries only; all globals).
  void download(bool candidate, std::vector<std::uint64_t> &keys, std::vector<std::int64_t> &entry,
                std::vector<std::int64_t> &global);
  void load(const std::vector<std::uint64_t> &keys, const std::vector<std::int64_t> &entry,
            const std::vector<std::int64_t> &global);
  std::string session_json() const;
  double kernel_ms_total() const noexcept;

private:
  struct Impl;
  std::unique_ptr<Impl> p_;
};

}  // namespace pf_gas_v1
