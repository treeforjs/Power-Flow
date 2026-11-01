#pragma once
#include <cstdlib>
#include <new>
#include <cstddef>
namespace allocation_probe {
inline thread_local long fail_after=-1;
inline thread_local std::size_t allocations=0;
inline void tick() {
  if(fail_after==0) throw std::bad_alloc();
  if(fail_after>0) --fail_after;
  ++allocations;
}
}
