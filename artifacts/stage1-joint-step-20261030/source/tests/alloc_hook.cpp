// Global allocation counter and one-shot countdown failure for the Stage 1
// tests. Single-threaded process only. A countdown of k makes the (k+1)-th
// subsequent operator new throw std::bad_alloc once; -1 disables injection.
#include <cstdlib>
#include <new>

namespace {
long long allocations = 0;
long long countdown = -1;
long long injected = 0;

void *allocate(std::size_t n) {
  ++allocations;
  if (countdown == 0) {
    countdown = -1;
    ++injected;
    throw std::bad_alloc();
  }
  if (countdown > 0)
    --countdown;
  if (void *p = std::malloc(n ? n : 1))
    return p;
  throw std::bad_alloc();
}

void *allocate_aligned(std::size_t n, std::align_val_t a) {
  ++allocations;
  if (countdown == 0) {
    countdown = -1;
    ++injected;
    throw std::bad_alloc();
  }
  if (countdown > 0)
    --countdown;
  const std::size_t align = static_cast<std::size_t>(a);
  std::size_t size = n ? n : 1;
  size = (size + align - 1) / align * align;
  if (void *p = std::aligned_alloc(align, size))
    return p;
  throw std::bad_alloc();
}
} // namespace

extern "C" long long stage1_allocation_count() noexcept { return allocations; }
extern "C" void stage1_set_allocation_countdown(long long k) noexcept { countdown = k; }
extern "C" long long stage1_allocation_countdown() noexcept { return countdown; }
extern "C" long long stage1_injected_allocation_failures() noexcept { return injected; }

void *operator new(std::size_t n) { return allocate(n); }
void *operator new[](std::size_t n) { return allocate(n); }
void *operator new(std::size_t n, const std::nothrow_t &) noexcept {
  try {
    return allocate(n);
  } catch (...) {
    return nullptr;
  }
}
void *operator new[](std::size_t n, const std::nothrow_t &) noexcept {
  try {
    return allocate(n);
  } catch (...) {
    return nullptr;
  }
}
void *operator new(std::size_t n, std::align_val_t a) { return allocate_aligned(n, a); }
void *operator new[](std::size_t n, std::align_val_t a) { return allocate_aligned(n, a); }
void *operator new(std::size_t n, std::align_val_t a, const std::nothrow_t &) noexcept {
  try {
    return allocate_aligned(n, a);
  } catch (...) {
    return nullptr;
  }
}
void *operator new[](std::size_t n, std::align_val_t a, const std::nothrow_t &) noexcept {
  try {
    return allocate_aligned(n, a);
  } catch (...) {
    return nullptr;
  }
}
void operator delete(void *p) noexcept { std::free(p); }
void operator delete[](void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }
void operator delete[](void *p, std::size_t) noexcept { std::free(p); }
void operator delete(void *p, std::align_val_t) noexcept { std::free(p); }
void operator delete[](void *p, std::align_val_t) noexcept { std::free(p); }
void operator delete(void *p, std::size_t, std::align_val_t) noexcept { std::free(p); }
void operator delete[](void *p, std::size_t, std::align_val_t) noexcept { std::free(p); }
void operator delete(void *p, const std::nothrow_t &) noexcept { std::free(p); }
void operator delete[](void *p, const std::nothrow_t &) noexcept { std::free(p); }
void operator delete(void *p, std::align_val_t, const std::nothrow_t &) noexcept { std::free(p); }
void operator delete[](void *p, std::align_val_t, const std::nothrow_t &) noexcept { std::free(p); }
