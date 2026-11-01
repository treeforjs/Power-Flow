#pragma once

#include <pf_gas_v1/gas_owner.hpp>

#include <memory>
#include <stdexcept>
#include <utility>

namespace pf_coupled {

// One coordinator trial's gas rollback scope. It keeps the owner alive and
// abandons only the attempt it began, even after mutation or a failed stage.
// The coordinator must still jointly preflight all owners before publication.
// Serial use only; this is not a thread synchronization mechanism.
class GasStepScope {
public:
  GasStepScope() = default;
  GasStepScope(const GasStepScope &) = delete;
  GasStepScope &operator=(const GasStepScope &) = delete;
  GasStepScope(GasStepScope &&) noexcept = default;

  GasStepScope &operator=(GasStepScope &&other) noexcept {
    if (this != &other) {
      reset();
      owner_ = std::move(other.owner_);
      attempt_ = std::move(other.attempt_);
    }
    return *this;
  }

  ~GasStepScope() noexcept { reset(); }

  static GasStepScope begin(std::shared_ptr<pf_gas_v1::Owner> owner,
                            std::uint64_t step) {
    if (!owner) throw std::invalid_argument("gas step scope requires an owner");
    auto attempt = owner->begin_attempt(step);
    return GasStepScope(std::move(owner), std::move(attempt));
  }

  bool active() const noexcept {
    return owner_ && owner_->matches_attempt(attempt_);
  }

  void reset() noexcept {
    if (owner_) owner_->abandon_attempt(attempt_);
    owner_.reset();
  }

private:
  GasStepScope(std::shared_ptr<pf_gas_v1::Owner> owner,
               pf_gas_v1::Owner::AttemptBinding attempt) noexcept
      : owner_(std::move(owner)), attempt_(std::move(attempt)) {}

  std::shared_ptr<pf_gas_v1::Owner> owner_;
  pf_gas_v1::Owner::AttemptBinding attempt_;
};

}  // namespace pf_coupled
