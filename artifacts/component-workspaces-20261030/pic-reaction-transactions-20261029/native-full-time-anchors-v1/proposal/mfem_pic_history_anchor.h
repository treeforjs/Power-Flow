#pragma once
// Draft API ONLY. No implementation/qualification is implied by these declarations.
#include "raikou/mfem_pic_full_time.h"
namespace raikou::mfem_backend {
struct NativeMonotoneIdPolicy {
  std::string integration_owner_id;
  // Immutable closed single-owner namespace v1; counters start at zero.
};
enum class NativeAnchorRetention { FullReplaySegmentV1, CompactCommitmentsV1 };
struct NativeAnchorIdentity {
  std::string schema, source_sha256, policy_sha256, native_config_sha256;
  std::string registry_sha256, original_provenance_sha256;
  std::uint64_t anchor_sequence=0,next_event_sequence=0,next_native_id=0;
  FullTimeAcceptedBoundary accepted;
  std::string canonical_state_sha256,native_ledger_sha256;
  std::string previous_anchor_sha256,base_history_chain,tail_history_chain;
  std::string segment_sha256,payload_sha256;
  NativeAnchorRetention retention=NativeAnchorRetention::FullReplaySegmentV1;
};
// MUST come from a separately pinned trusted coordinator manifest/opaque native
// export. Do not fill this from the same untrusted file being restored.
struct ExpectedNativeAnchorIdentity { NativeAnchorIdentity independently_pinned; };
class PreparedNativeAnchor {
public:
  PreparedNativeAnchor(PreparedNativeAnchor&&) noexcept;
  PreparedNativeAnchor& operator=(PreparedNativeAnchor&&) noexcept;
  ~PreparedNativeAnchor();
  PreparedNativeAnchor(const PreparedNativeAnchor&)=delete;
  PreparedNativeAnchor& operator=(const PreparedNativeAnchor&)=delete;
  const NativeAnchorIdentity& identity() const;
  std::size_t anchor_bytes() const;
  std::size_t segment_bytes() const;
  // File IO/fsync happens before publication; root owns all-owner durable commit.
  void write_files(const std::string& anchor_path,const std::string& segment_path);
private:
  friend class PicArena;
  struct Impl;
  explicit PreparedNativeAnchor(std::unique_ptr<Impl>);
  std::unique_ptr<Impl> impl_;
};
// Proposed PicArena member additions (shown separately to avoid declaring a
// competing class while the coordinator reviews the exact API):
// static unique_ptr<PicArena> create_full_time_anchored(options,states,channels,
//     authored,genesis,NativeMonotoneIdPolicy);
// const NativeMonotoneIdPolicy& native_id_policy() const;
// uint64_t next_event_sequence() const;
// PreparedNativeAnchor prepare_history_anchor(NativeAnchorRetention);
// bool can_publish(const PreparedNativeAnchor&) const noexcept;
// bool try_publish(PreparedNativeAnchor&) noexcept;
// static unique_ptr<PicArena> restore_history_anchor(path,options,states,channels,
//     const ExpectedNativeAnchorIdentity&);
// void save_anchored_checkpoint(path) const;
// static unique_ptr<PicArena> restore_anchored_checkpoint(path,options,states,
//     channels,const ExpectedNativeAnchorIdentity& expected_base,
//     const FullTimeAcceptedBoundary& expected_current);
// Pure naming helpers; they confer no mutation authority:
std::string native_monotone_event_id(const NativeMonotoneIdPolicy&,std::uint64_t);
std::string native_monotone_product_id(const NativeMonotoneIdPolicy&,
                                      std::uint64_t,std::size_t product_slot);
std::string native_monotone_primary_remainder_id(const NativeMonotoneIdPolicy&,
                                                std::uint64_t);
std::string native_monotone_genesis_id(std::size_t genesis_slot);
}
