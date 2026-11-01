#include "allocation_probe.hpp"
void *operator new(std::size_t size) {
  allocation_probe::tick();
  if(void *p=std::malloc(size?size:1))return p;
  throw std::bad_alloc();
}
void *operator new[](std::size_t size){return ::operator new(size);}
void operator delete(void *p) noexcept {std::free(p);}
void operator delete[](void *p) noexcept {std::free(p);}
void operator delete(void *p,std::size_t) noexcept {std::free(p);}
void operator delete[](void *p,std::size_t) noexcept {std::free(p);}
void *operator new(std::size_t size,std::align_val_t alignment) {
  allocation_probe::tick();void *p=nullptr;
  if(posix_memalign(&p,std::size_t(alignment),size?size:1)!=0)throw std::bad_alloc();
  return p;
}
void *operator new[](std::size_t size,std::align_val_t a){return ::operator new(size,a);}
void operator delete(void *p,std::align_val_t) noexcept {std::free(p);}
void operator delete[](void *p,std::align_val_t) noexcept {std::free(p);}
void operator delete(void *p,std::size_t,std::align_val_t) noexcept {std::free(p);}
void operator delete[](void *p,std::size_t,std::align_val_t) noexcept {std::free(p);}
