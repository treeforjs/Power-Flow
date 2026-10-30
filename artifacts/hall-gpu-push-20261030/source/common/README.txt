common/ -- device session, NVRTC compilation and deterministic device primitives
for the Hall GPU-resident step (shared by push/ and field/). Stable interface;
changes are additive.

Files
  hall_gpu_common.hpp/.cpp     Session, buffers, NVRTC compile, timing, teardown, JSON helpers
  hall_gpu_primitives.hpp/.cpp host launchers of the deterministic primitives
  device/deterministic_primitives.cuh   NVRTC device code: u32 exclusive scan, stable
                               LSD radix sort (u32 key + u32 payload), segment bounds,
                               exact max(f64 >= 0) / min(u32) reductions

Build
  c++ -std=c++20 -O2 -ffp-contract=off -fno-fast-math
      -isystem <cuda_runtime 12.8 include (cuda.h)> -isystem <cuda_nvrtc 12.8 include (nvrtc.h)>
      hall_gpu_common.cpp hall_gpu_primitives.cpp ... -ldl -lcrypto
  No link against libcuda or NVRTC: libcuda.so.1 is dlopen'ed (path and sha256 recorded),
  NVRTC is dlopen'ed from the absolute path in SessionOptions::nvrtc_library.
  hall_gpu_primitives.cpp embeds the .cuh with #embed (GCC 16; -Wno-c++26-extensions).

Session (one device, one primary context, one stream)
  SessionOptions{require_uuid, ordinal, nvrtc_library, blocking_sync=true,
                 check_each_launch=false, event_pool=96}
  Session s(opts);                 selects the device by UUID, sets CU_CTX_SCHED_BLOCKING_SYNC
                                   on the primary context (refuses if it is already active
                                   with another policy), retains + pushes it, creates ONE
                                   CU_STREAM_NON_BLOCKING stream.
  s.device()                       DeviceInfo: name, uuid, cc, SMs, driver API version,
                                   /proc/driver/nvidia/version line, libcuda path + sha256.
  s.alloc(bytes, tag) -> DeviceBuffer      tracked; own_peak_bytes recorded
  s.alloc_pinned(bytes, tag) -> HostBuffer page-locked staging
  s.release(buffer)                drains the stream, then frees
  s.upload(dst, off, host, off, n) / s.download(...) / s.zero(...)   async, on the stream
  s.launch(kernel, grid, block, shmem, args)                          on the stream
  s.synchronize(where)             stream sync; surfaces asynchronous faults
  s.mark(label) / s.wait_marker(m) / s.elapsed_ms(a, b) / s.clear_markers()
  s.compile_and_load(CompileRequest{main, headers, extra_options}) -> Module
  s.kernel(module, "name") -> Kernel       (kernels are extern "C")
  s.memory() -> MemoryReport       own in-use/peak/pinned bytes, device free/total/min-free
  s.close() -> TeardownReport      also run by the destructor

Ordering rules (closing the two earlier review findings)
  1. Static data is never uploaded with synchronous/pageable copies or on the legacy
     default stream. Every upload is cuMemcpyHtoDAsync from pinned memory on the
     session stream, so any later launch on that stream sees it; there is no implicit
     cross-stream ordering anywhere.
  2. The host writes a HostBuffer only after a synchronize() (or wait_marker() on a
     marker recorded after the last operation reading it). Downloads are readable only
     after the synchronize() that follows them.
  3. Teardown drains the stream (falling back to a context synchronize) before
     destroying events, modules, device or pinned memory. If completion cannot be
     established, nothing is freed (pinned pages may still be DMA targets), the report
     says quarantined=true and the process must exit non-zero. A failed driver call
     marks the session failed; all later calls refuse; buffers held by consumers are
     then released only by the session teardown.

NVRTC
  compile_and_load() picks --gpu-architecture=sm_<major><minor> from the device when
  NVRTC supports it (cubin, no driver JIT); otherwise compute_<highest supported <= device>
  PTX for driver JIT (image_kind "ptx-jit"). Always adds
    --fmad=false --prec-div=true --prec-sqrt=true --ftz=false --std=c++17
  and refuses extra options that touch architecture or floating-point behaviour.
  CompileRecord keeps NVRTC version, library and builtins sha256, supported archs,
  device/target cc, full option list, log, source and image sha256, compile time.

Deterministic primitives (device/deterministic_primitives.cuh, blockDim 256)
  Primitives p(session, module);   module must be compiled with
                                   Primitives::device_source() among its headers
  p.exclusive_scan(in, n, out, scratch, total)        u32, exact
  p.radix_sort(key_a, pay_a, key_b, pay_b, n, bits, scratch) -> result in b?
                                   stable; equal keys keep input order (warp match +
                                   popc ranks, no atomics in the scatter)
  p.segment_bounds(sorted_keys, n, begin, end)        begin/end zeroed by the caller
  p.reduce_max_f64 / p.reduce_min_u32                 exact, order independent
  Scratch sizes: Primitives::scan_scratch / sort_scratch / reduce_scratch.
  Floating-point sums are never formed by these primitives; a deterministic sum is
  obtained by sorting (index, record order) and summing each segment sequentially
  in the required order (see push/push_kernels.cu), or by a fixed reduction tree
  written for that purpose.

Execution under the frozen launcher
  arch-base-20261028/root contains usr/lib/libcuda.so.615.78.08 byte-identical to the
  host driver library (sha256 f73320931be3c4ea..., sealed in root-files.json) and
  run-frozen binds /dev, so programs using this layer run inside run-frozen without
  extra binds. Do not apply prlimit --as to GPU processes (the driver reserves large
  virtual address ranges).
