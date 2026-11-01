#pragma once
// Native v7.1 device lane (API-071.txt section 2). A FullTimeDevice runs on
// the producer's hall-gpu-resident-v1 Session (same device, same stream) and
// its pf_v7 forms module. Attached to a full-time arena it keeps the accepted
// native carrier rows resident on the device and admits motion records on the
// device: per-carrier checks, exact sums, motion-work-gross-v1 scales, section
// roots, census v2, field/current digests, the after-state v5 carrier root and
// the motion row path root. The accept/refuse decision and every digest are
// those of the host path (v7 code); a refusal is reported with the host
// path's text (the host path is re-run for it). Execution modes:
//   host            v7 host path only (reference)
//   device          device lane where eligible, host path otherwise
//   device_verified both; outputs compared, a difference is an internal error
#include "raikou/mfem_pic_full_time.h"
#include <cstdint>
#include <memory>

namespace hall_gpu {
class Session;
}
namespace pf_v7::gpu {
class Forms;
}

namespace raikou::mfem_backend {
enum class FullTimeExecution { host = 0, device = 1, device_verified = 2 };
struct FullTimeDeviceCounters {
  std::uint64_t device_admissions = 0;   // motion records admitted on the device lane
  std::uint64_t host_admissions = 0;     // host path (ineligible step or host mode)
  std::uint64_t host_refusal_texts = 0;  // device refusal re-run on the host for its text
  std::uint64_t verified_comparisons = 0;
  std::uint64_t row_uploads = 0, row_upload_bytes = 0;
  std::uint64_t row_downloads = 0, row_download_bytes = 0;
  std::uint64_t array_uploads = 0, array_upload_bytes = 0; // record arrays given only on the host
  std::uint64_t array_downloads = 0, array_download_bytes = 0; // device field arrays read back (streaming forms)
  std::uint64_t streaming_field_digests = 0; // v7 streaming digest forms (hashed on the host)
};
class FullTimeDevice {
public:
  FullTimeDevice(hall_gpu::Session &, pf_v7::gpu::Forms &);
  ~FullTimeDevice();
  FullTimeDevice(const FullTimeDevice &) = delete;
  FullTimeDevice &operator=(const FullTimeDevice &) = delete;
  FullTimeDeviceCounters counters() const;
  // Host row materialisations of every arena using this device (NativeRows).
  static std::uint64_t row_materialisations() noexcept;
  struct Impl;
  std::shared_ptr<Impl> impl_;
};
// Census v2 of a trial on the device (zero copy for the producer): rows (80 B
// GapParticle records), packet keys and state indices (u64 each), ascending
// native ID; addresses on the shared session; valid while the trial lives.
struct FullTimeCensusDevice {
  std::uint64_t rows = 0, packet_keys = 0, state_indices = 0, carriers = 0;
  std::string states_sha256, carriers_root, sha256;
};
} // namespace raikou::mfem_backend
