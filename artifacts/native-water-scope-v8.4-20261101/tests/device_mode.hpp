#pragma once
// Native v7.1 qualification: runs an unchanged v7 test suite in a device
// execution mode. Environment:
//   PF_V71_EXECUTION  host (default) | device | device_verified
//   PF_V71_NVRTC      absolute path of libnvrtc.so.12
//   PF_V71_UUID       device UUID (GPU-...)
// install_from_environment() creates one session, forms module and device
// lane for the process and makes every full-time arena created or restored
// afterwards attach it (PicArena::set_process_full_time_device).
#include "forms_gpu.hpp"
#include "hall_gpu_common.hpp"
#include "raikou/mfem_pic_arena.h"
#include "raikou/mfem_pic_full_time_device.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

namespace v71 {
namespace mm = raikou::mfem_backend;
inline mm::FullTimeExecution env_mode() {
  const char *v = std::getenv("PF_V71_EXECUTION");
  const std::string s = v ? v : "host";
  if (s == "host")
    return mm::FullTimeExecution::host;
  if (s == "device")
    return mm::FullTimeExecution::device;
  if (s == "device_verified")
    return mm::FullTimeExecution::device_verified;
  throw std::runtime_error("PF_V71_EXECUTION must be host, device or device_verified");
}
inline const char *mode_name(mm::FullTimeExecution m) {
  return m == mm::FullTimeExecution::host ? "host" : m == mm::FullTimeExecution::device ? "device" : "device_verified";
}
struct DeviceEnv {
  std::unique_ptr<hall_gpu::Session> session;
  std::unique_ptr<pf_v7::gpu::Forms> forms;
  std::shared_ptr<mm::FullTimeDevice> device;
  mm::FullTimeExecution mode = mm::FullTimeExecution::host;
  ~DeviceEnv() {
    mm::PicArena::set_process_full_time_device(nullptr, mm::FullTimeExecution::host);
    device.reset();
    forms.reset();
    if (session)
      (void)session->close();
  }
};
inline DeviceEnv &env() {
  static DeviceEnv e;
  return e;
}
inline void install_from_environment() {
  const auto mode = env_mode();
  if (mode == mm::FullTimeExecution::host)
    return;
  auto &e = env();
  hall_gpu::SessionOptions so;
  const char *uuid = std::getenv("PF_V71_UUID"), *nvrtc = std::getenv("PF_V71_NVRTC");
  if (!uuid || !nvrtc)
    throw std::runtime_error("PF_V71_UUID and PF_V71_NVRTC are required in device modes");
  so.require_uuid = uuid;
  so.nvrtc_library = nvrtc;
  e.session = std::make_unique<hall_gpu::Session>(so);
  e.forms = std::make_unique<pf_v7::gpu::Forms>(*e.session);
  e.device = std::make_shared<mm::FullTimeDevice>(*e.session, *e.forms);
  e.mode = mode;
  mm::PicArena::set_process_full_time_device(e.device, mode);
}
inline void write_report_only(const std::string &path);
// Releases the device lane and closes the session (before static destruction).
inline void shutdown() {
  auto &e = env();
  mm::PicArena::set_process_full_time_device(nullptr, mm::FullTimeExecution::host);
  e.device.reset();
  e.forms.reset();
  if (e.session) {
    (void)e.session->close();
    e.session.reset();
  }
}
// Execution report (separate from the suite's own output, which must be equal
// in every mode); then shutdown().
inline void write_report(const std::string &path) {
  write_report_only(path);
  shutdown();
}
inline void write_report_only(const std::string &path) {
  const auto &e = env();
  std::ofstream f(path);
  f << "{\"execution\":\"" << mode_name(e.mode) << "\"";
  if (e.device) {
    const auto c = e.device->counters();
    f << ",\"device_admissions\":" << c.device_admissions << ",\"host_admissions\":" << c.host_admissions
      << ",\"host_refusal_texts\":" << c.host_refusal_texts << ",\"verified_comparisons\":" << c.verified_comparisons
      << ",\"row_uploads\":" << c.row_uploads << ",\"row_upload_bytes\":" << c.row_upload_bytes
      << ",\"row_downloads\":" << c.row_downloads << ",\"row_download_bytes\":" << c.row_download_bytes
      << ",\"array_uploads\":" << c.array_uploads << ",\"array_upload_bytes\":" << c.array_upload_bytes
      << ",\"array_downloads\":" << c.array_downloads << ",\"array_download_bytes\":" << c.array_download_bytes
      << ",\"streaming_field_digests\":" << c.streaming_field_digests
      << ",\"host_row_materialisations\":" << mm::FullTimeDevice::row_materialisations();
  }
  f << "}\n";
}
} // namespace v71
