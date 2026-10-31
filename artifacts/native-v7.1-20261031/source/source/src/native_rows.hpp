// native_rows.hpp -- native v7.1 carrier rows: immutable shared content with
// copy-on-write, held on the host, on the device, or both (API-071.txt 2).
//
// Copying a NativeRows copies two pointers (frames, tokens and trials share
// the rows of their parent until one of them changes them). Read access
// materialises the host copy from the device copy when only the device copy
// exists (counted: a device-mode step without host access performs none).
// Mutation makes the host copy unique and drops the device copy.
#pragma once
#include "raikou/mfem_gap_pic.h"
#include <cstddef>
#include <memory>
#include <vector>

namespace raikou::mfem_backend {
struct DeviceRows; // defined by the device lane (native_device.cpp)
// Device lane hooks (native_device.cpp).
std::vector<GapParticle> device_rows_download(const DeviceRows &);

class NativeRows {
public:
  using value_type = GapParticle;
  using size_type = std::size_t;
  using const_iterator = const GapParticle *;
  using iterator = GapParticle *;
  NativeRows() = default;
  NativeRows(std::vector<GapParticle> v)
      : host_(std::make_shared<std::vector<GapParticle>>(std::move(v))) {
    n_ = host_->size();
  }
  NativeRows &operator=(std::vector<GapParticle> v) {
    host_ = std::make_shared<std::vector<GapParticle>>(std::move(v));
    n_ = host_->size();
    device_.reset();
    return *this;
  }
  // ---- read access (materialises the host copy when needed)
  size_type size() const noexcept { return host_ ? host_->size() : n_; }
  bool empty() const noexcept { return size() == 0; }
  const std::vector<GapParticle> &host() const {
    if (!host_) {
      if (device_) {
        host_ = std::make_shared<std::vector<GapParticle>>(device_rows_download(*device_));
        ++materialisations_;
      } else
        host_ = std::make_shared<std::vector<GapParticle>>();
    }
    return *host_;
  }
  operator const std::vector<GapParticle> &() const { return host(); }
  const GapParticle *data() const { return host().data(); }
  const_iterator begin() const { return host().data(); }
  const_iterator end() const { return host().data() + host().size(); }
  const_iterator cbegin() const { return begin(); }
  const_iterator cend() const { return end(); }
  const GapParticle &operator[](size_type i) const { return host()[i]; }
  const GapParticle &front() const { return host().front(); }
  const GapParticle &back() const { return host().back(); }
  // ---- mutation (copy-on-write)
  std::vector<GapParticle> &mutate() {
    (void)host();
    if (host_.use_count() != 1)
      host_ = std::make_shared<std::vector<GapParticle>>(*host_);
    device_.reset();
    return const_cast<std::vector<GapParticle> &>(*host_);
  }
  iterator begin() { return mutate().data(); }
  iterator end() {
    auto &v = mutate();
    return v.data() + v.size();
  }
  GapParticle &operator[](size_type i) { return mutate()[i]; }
  void reserve(size_type n) { mutate().reserve(n); }
  void push_back(const GapParticle &p) {
    mutate().push_back(p);
    ++n_;
  }
  void erase(iterator first, iterator last) {
    auto &v = mutate();
    v.erase(v.begin() + (first - v.data()), v.begin() + (last - v.data()));
    n_ = v.size();
  }
  void clear() { *this = NativeRows(); }
  void swap(NativeRows &o) noexcept {
    host_.swap(o.host_);
    device_.swap(o.device_);
    std::swap(n_, o.n_);
  }
  // Exchanges the content with a plain vector (GapPic rows).
  void swap(std::vector<GapParticle> &v) {
    auto &mine = mutate();
    mine.swap(v);
    n_ = mine.size();
  }
  // ---- sharing and device copies
  bool same_as(const NativeRows &o) const noexcept {
    return size() == o.size() && ((host_ && host_ == o.host_) || (device_ && device_ == o.device_));
  }
  bool host_present() const noexcept { return bool(host_); }
  // Shared host copy (materialised); keeps the content alive (memo keys).
  std::shared_ptr<const std::vector<GapParticle>> host_shared() const {
    (void)host();
    return host_;
  }
  const std::shared_ptr<const DeviceRows> &device() const noexcept { return device_; }
  // Attaches a device copy of the same content (host copy kept if present).
  void attach_device(std::shared_ptr<const DeviceRows> d) const { device_ = std::move(d); }
  // Device-only content of n rows (no host copy until read).
  static NativeRows from_device(std::shared_ptr<const DeviceRows> d, size_type n) {
    NativeRows r;
    r.device_ = std::move(d);
    r.n_ = n;
    return r;
  }
  static std::uint64_t materialisations() noexcept { return materialisations_; }

private:
  mutable std::shared_ptr<const std::vector<GapParticle>> host_;
  mutable std::shared_ptr<const DeviceRows> device_;
  size_type n_ = 0;
  static inline std::uint64_t materialisations_ = 0;
};
} // namespace raikou::mfem_backend
