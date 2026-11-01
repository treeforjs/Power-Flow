// forms_gpu.hpp -- device producer of the native v7 forms (NVRTC kernels in
// forms_kernels.cu) on a hall-gpu-resident-v1 Session (common/). Every call
// is enqueued on the session stream; addresses are raw device addresses
// (DeviceBuffer::address + byte offset). Results are readable after the
// caller's Session::synchronize(). Bit identity with the CPU reference
// (pf_v7/forms.hpp) is the acceptance test (forms/tests/test_forms_gpu.cpp).
#pragma once
#include "hall_gpu_common.hpp"
#include "hall_gpu_primitives.hpp"
#include "pf_v7/forms.hpp"
#include <memory>

namespace pf_v7::gpu {
// Device result of an exact sum (layout of struct SumOut in forms_kernels.cu).
struct SumOut {
  double value = 0;  // RN(exact sum)
  double scaled = 0; // RN(scale x value); field work: scale = dt
  std::uint32_t status = 0, parts = 0;
  std::uint64_t terms = 0;
  std::int64_t digit[pf7::acc_digits]{}; // normalised digits (interchange form)
};
static_assert(sizeof(SumOut) == 32 + 8 * pf7::acc_digits);

class Forms {
public:
  static constexpr unsigned block = 256;
  // Compiles forms_kernels.cu (embedded) with forms_core.h and the
  // primitives on the session device.
  explicit Forms(hall_gpu::Session &);
  ~Forms();
  const hall_gpu::CompileRecord &compile_record() const;
  hall_gpu::Primitives &primitives();
  static const std::string &kernel_source_sha256(); // sha256 of the embedded sources (kernels + core)

  // ---- (a) exact sums. scratch: sum_scratch_bytes(blocks); out: one SumOut.
  // blocks: any count >= 1 (the result does not depend on it).
  static std::size_t sum_scratch_bytes(unsigned blocks);
  void sum_products(std::uint64_t n, std::uint64_t a, std::uint64_t a_stride, std::uint64_t b,
                    std::uint64_t b_stride, double scale, unsigned blocks, std::uint64_t scratch,
                    std::uint64_t out);
  void kinetic_aligned(std::uint64_t n, std::uint64_t w, std::uint64_t w_stride, std::uint64_t m,
                       std::uint64_t m_stride, std::uint64_t u0, std::uint64_t u0_stride, std::uint64_t u1,
                       std::uint64_t u1_stride, unsigned blocks, std::uint64_t scratch, std::uint64_t out);
  void kinetic_split(std::uint64_t ne, std::uint64_t endpoints, std::uint64_t ns, std::uint64_t sinks,
                     std::uint64_t nc, std::uint64_t census, std::uint64_t census_words, unsigned blocks,
                     std::uint64_t scratch, std::uint64_t out);

  // ---- (b) block lists. scratch: block_list_scratch_bytes(n); root_out: 8
  // u32 state words (bytes via pf_v7::digest_of_words). The leaf digests stay
  // in scratch at leaf_offset(n) (8 u32 each) for host verification.
  static std::size_t block_list_scratch_bytes(std::uint64_t n);
  static std::size_t leaf_offset(std::uint64_t n);
  void block_list(std::uint32_t section, std::uint64_t n, std::uint64_t records, std::uint64_t scratch,
                  std::uint64_t root_out);

  // ---- (c) census. Records: 12 u64 words (pf_v7::CensusRecord).
  void census_pack(std::uint64_t n, std::uint64_t rows, std::uint64_t row_words, std::uint64_t keys,
                   std::uint64_t states, std::uint64_t out);
  static std::size_t census_check_scratch_bytes(std::uint64_t n);
  // out: 2 u64 = first failing index (~0 when none), pf7::CensusCode.
  void census_check(std::uint64_t n, std::uint64_t records, std::uint64_t states, std::uint64_t scratch,
                    std::uint64_t out);
  static std::size_t census_order_scratch_bytes(std::uint64_t n);
  // Stable ascending-ID order of n records of rec_words words; writes the
  // ordered records to out (n < 2^32).
  void census_order(std::uint64_t n, std::uint64_t records, std::uint64_t rec_words, std::uint64_t scratch,
                    std::uint64_t out);

  // ---- native v7.1 (additive; API-071.txt)
  // Fused block lists: up to 8 sections, one launch per level (chunks,
  // leaves, trees); same digests as block_list().
  struct BlockListJob {
    std::uint32_t section = 0;
    std::uint64_t n = 0, records = 0, root_out = 0;
  };
  static std::size_t block_lists_scratch_bytes(std::span<const BlockListJob>);
  void block_lists(std::span<const BlockListJob>, std::uint64_t scratch);
  // v5 ID-block Merkle root (mfem_pic_canonical carrier/path roots) of
  // ascending-ID records of record_words words. kind 0: native-v5/carrier-*,
  // kind 1: native-v5/path-*. root_out: 8 state words; leaf count at
  // scratch + id_root_count_offset(n) (u32).
  static std::size_t id_root_scratch_bytes(std::uint64_t n);
  static std::size_t id_root_count_offset(std::uint64_t n);
  void id_root(std::uint32_t kind, std::uint64_t n, std::uint64_t records, std::uint64_t record_words,
               std::uint64_t scratch, std::uint64_t root_out);
  // Both kinds in one leaf launch: kind 0 over A, kind 1 over B, where A and B
  // carry the same ascending ID sequence (word 0 of each record; the leaf
  // partition is taken from A). Same roots as id_root(0, A) and id_root(1, B).
  static std::size_t id_root_pair_scratch_bytes(std::uint64_t n);
  void id_root_pair(std::uint64_t n, std::uint64_t recA, std::uint64_t wordsA, std::uint64_t recB,
                    std::uint64_t wordsB, std::uint64_t scratch, std::uint64_t rootA, std::uint64_t rootB);
  // Exact sum of |a_i| |b_i| (A1 G_W): as sum_products.
  void sum_abs_products(std::uint64_t n, std::uint64_t a, std::uint64_t a_stride, std::uint64_t b,
                        std::uint64_t b_stride, double scale, unsigned blocks, std::uint64_t scratch,
                        std::uint64_t out);
  // A1 E^: exact maximum of |E_e| over edges with an endpoint node on an edge
  // carrying J != 0 (Hall axial mesh n, nz); out = {f64 max, u32 nonfinite}.
  static std::size_t field_bound_scratch_bytes(std::uint64_t ne, std::uint64_t nn);
  void field_bound(std::uint64_t ne, std::uint64_t E, std::uint64_t J, std::uint64_t n, std::uint64_t nz,
                   std::uint64_t scratch, std::uint64_t out);
  // v7 streaming digests on one thread (kind 0: total-midpoint-E-B-binary64-v1
  // over E then B; kind 1: motion-array-binary64-v1 over E). out: 8 words.
  void stream_digest(std::uint32_t kind, std::uint64_t ne, std::uint64_t E, std::uint64_t nb, std::uint64_t B,
                     std::uint64_t out);

private:
  struct Impl;
  std::unique_ptr<Impl> p_;
};

// Blocking helpers over host arrays (verification, tests, small producers):
// allocate, upload through pinned staging, run, download, release.
class HostRunner {
public:
  explicit HostRunner(hall_gpu::Session &s, Forms &f) : s_(s), f_(f) {}
  SumOut sum_products(std::span<const double> a, std::span<const double> b, double scale, unsigned blocks);
  SumOut kinetic_aligned(std::span<const Carrier> census, std::span<const double> u_end, unsigned blocks);
  SumOut kinetic_split(std::span<const Carrier> census, std::span<const Carrier> endpoints,
                       std::span<const Sink> sinks, unsigned blocks);
  struct BlockListOut {
    Digest root{};
    std::vector<Digest> leaves;
  };
  BlockListOut block_list(std::uint32_t section, const void *records, std::uint64_t n);
  CensusCheck census_check(std::span<const CensusRecord>, std::uint64_t states);
  std::vector<CensusRecord> census_order(std::span<const CensusRecord>);
  double last_device_ms = 0; // kernel time of the last call (markers)

private:
  hall_gpu::Session &s_;
  Forms &f_;
};
} // namespace pf_v7::gpu
