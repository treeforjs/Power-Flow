#pragma once
#include "raikou/mfem_gap_pic.h"
#include "reaction_transactions.hpp"
#include <map>
#include <memory>

namespace raikou::mfem_backend {
using ReactionPacketBindings = std::map<std::string, std::uint64_t>;
struct ReactionPhase {
    enum class Kind { InitialSynchronizedFieldProof, MidXEqualStepBoundary };
    Kind kind = Kind::InitialSynchronizedFieldProof;
    // Frozen physical operator placement. SSA subflow event-clock time is not
    // an instruction to drift the particles or advance this native clock.
    double expected_time_s = 0, previous_dt_s = 0, next_dt_s = 0;
    std::uint64_t accepted_step = 0;
};
struct PrimaryWeightSlice {
    double event_weight = 0;
    std::uint64_t untouched_primary_native_id = 0;
    // v1 admits only an entire primary. Future partial splitting must stage the
    // exact untouched remainder with this fresh reserved ID and unchanged u/x.
};
enum class NativeReactionFault { None, AfterPopulationStage, AfterFieldSolve, BeforePreparedReturn };
class PreparedPicReaction {
  public:
    PreparedPicReaction(PreparedPicReaction &&) noexcept;
    PreparedPicReaction &operator=(PreparedPicReaction &&) noexcept;
    ~PreparedPicReaction();
    PreparedPicReaction(const PreparedPicReaction &) = delete;
    PreparedPicReaction &operator=(const PreparedPicReaction &) = delete;
    const ReactionPacketBindings &output_bindings() const;
    const raikou::pic_reactions::Json &audit() const;

  private:
    friend class PicReactionBatch;
    struct Impl;
    explicit PreparedPicReaction(std::unique_ptr<Impl>);
    std::unique_ptr<Impl> impl_;
};
struct RestoredPicReaction {
    std::unique_ptr<GapPic> pic;
    std::unique_ptr<raikou::pic_reactions::Transactions> reactions;
    ReactionPacketBindings bindings;
};
// Single-owner exclusive-barrier API. No authored source is added. Field work is
// staged against a complete population and restored before prepare returns.
class PicReactionBatch {
  public:
    struct NativeState;
    static PreparedPicReaction prepare(GapPic &, raikou::pic_reactions::Transactions &,
                                       const raikou::pic_reactions::PreparedReaction &,
                                       ReactionPacketBindings &, const ReactionPhase &,
                                       const PrimaryWeightSlice &,
                                       NativeReactionFault = NativeReactionFault::None);
    static bool can_publish(const GapPic &, const raikou::pic_reactions::Transactions &,
                            const ReactionPacketBindings &, const PreparedPicReaction &) noexcept;
    static bool try_publish(GapPic &, raikou::pic_reactions::Transactions &, ReactionPacketBindings &,
                            PreparedPicReaction &) noexcept;
    static bool try_rollback(GapPic &, raikou::pic_reactions::Transactions &, ReactionPacketBindings &,
                             PreparedPicReaction &) noexcept;
    static std::string fingerprint(const GapPic &);
    static const char *source_identity() noexcept;
    static void save_bundle(const GapPic &, const raikou::pic_reactions::Transactions &,
                            const ReactionPacketBindings &, const std::string &path);
    static RestoredPicReaction load_bundle(const raikou::pic_reactions::Catalog &,
                                           const GapPicOptions &expected_options,
                                           const std::string &expected_coupling_prefix,
                                           const std::string &path);

  private:
    static NativeState capture(const GapPic &);
    static void swap_state(GapPic &, NativeState &) noexcept;
    static bool matches(const GapPic &, const NativeState &) noexcept;
};
} // namespace raikou::mfem_backend
